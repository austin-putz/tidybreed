# ---------------------------------------------------------------------------
# Genome-effects schema helpers: view SQL, family signature, origin-predicate
# containment, and the cross-row validator.
#
# Storage is three tables — one coefficient (genome_effects) over one or more
# loci (genome_effect_members), each optionally scoped to allele copies of a
# given line and/or parent of origin (genome_effect_member_origins). Row-local
# invariants are declared SQL constraints in define_genome(); everything here is
# a cross-row rule that SQL cannot express.
#
# Design and worked cases: plans/update_genome_effects_v4.md (v4.9).
# The fixtures these rules were derived from, with hand-computed values:
# tests/testthat/helper-genome-effects.R and plans/update_genome_effects_phase_A.md.
# ---------------------------------------------------------------------------


# ── View SQL ────────────────────────────────────────────────────────────────

#' SQL for the `genome_effect_terms` view (one row per term)
#'
#' `effect_order` (member count) and `family_key` (the family signature) are
#' **derived here and never stored**. Terms sharing a `family_key` are scope
#' variants of one mathematical term and compete under specificity fallback;
#' terms with different keys sum. That sentence is the whole mental model, and
#' this column is what lets a user check it without reading the design notes.
#'
#' @keywords internal
#' @noRd
.genome_effect_terms_view_sql <- function() {
  paste0(
    "CREATE VIEW genome_effect_terms AS ",
    "SELECT e.trait_name, e.effect_owner, e.effect_name, e.id_genome_effect, ",
    "       m.effect_order, m.contrast_signature, ",
    "       e.trait_name || '|' || e.effect_owner || '|' || m.member_signature ",
    "         AS family_key, ",
    "       COALESCE(o.scope_description, 'common') AS scope_description, ",
    "       e.genome_value ",
    "FROM genome_effects e ",
    "JOIN ( ",
    "  SELECT id_genome_effect, ",
    "         COUNT(*) AS effect_order, ",
    "         string_agg(contrast_name, '+' ORDER BY member_slot) ",
    "           AS contrast_signature, ",
    "         string_agg(locus_id || ':' || contrast_name || ':' || ",
    "                    COALESCE(CAST(copy_count_value AS VARCHAR), '-') || ':' || ",
    "                    COALESCE(CAST(dosage_value AS VARCHAR), '-'), ",
    "                    '&' ORDER BY locus_id) AS member_signature ",
    "  FROM genome_effect_members GROUP BY id_genome_effect ",
    ") m ON m.id_genome_effect = e.id_genome_effect ",
    "LEFT JOIN ( ",
    "  SELECT id_genome_effect, ",
    "         string_agg(member_slot || ':' || line_match_type || ",
    "                    COALESCE('(' || line_name || ')', '') || ",
    "                    COALESCE('@p' || parent_origin, '') || 'x' || copy_count, ",
    "                    ', ' ORDER BY member_slot, origin_slot) ",
    "           AS scope_description ",
    "  FROM genome_effect_member_origins GROUP BY id_genome_effect ",
    ") o ON o.id_genome_effect = e.id_genome_effect"
  )
}

#' SQL for the `genome_effect_loci` view (one row per term x locus)
#'
#' `locus_name` lives only in this view: it was removed from
#' `genome_effect_members` so there is no id/name agreement invariant to keep.
#'
#' @keywords internal
#' @noRd
.genome_effect_loci_view_sql <- function() {
  paste0(
    "CREATE VIEW genome_effect_loci AS ",
    "SELECT e.trait_name, e.effect_owner, m.id_genome_effect, m.member_slot, ",
    "       m.locus_id, g.locus_name, m.contrast_name, e.genome_value ",
    "FROM genome_effect_members m ",
    "JOIN genome_effects e ON e.id_genome_effect = m.id_genome_effect ",
    "JOIN genome_meta     g ON g.locus_id         = m.locus_id"
  )
}

#' SQL for the `ind_tgv_total` view (one row per individual x trait)
#'
#' The total is derived, never stored: a stored `'total'` row would make every
#' `SUM(tgv_value)` double-count.
#'
#' Phenotypes read this view, so it must be bit-identical whatever DuckDB's
#' thread count (CLAUDE.md). With up to four components per individual x trait
#' a parallel floating `SUM()` is not, so the components are added in a fixed
#' order (`component_name`) with `list_sum(list(... ORDER BY ...))`: the result
#' is a function of the stored rows alone. Not `GEV_ACC_TYPE`, which the
#' evaluator uses for its many-term sum: a `DOUBLE -> DECIMAL(38, 18) ->
#' DOUBLE` round trip is not exact (it moved ~10% of single-component totals
#' by one ulp), while an ordered floating sum of one row returns that row bit
#' for bit, so an additive-only trait's total *is* its breeding value.
#'
#' @keywords internal
#' @noRd
.ind_tgv_total_view_sql <- function() {
  paste0(
    "CREATE VIEW ind_tgv_total AS ",
    "SELECT id_ind, trait_name, ",
    "list_sum(list(tgv_value ORDER BY component_name)) AS tgv_total ",
    "FROM ind_tgv GROUP BY id_ind, trait_name"
  )
}


# ── Family signature ────────────────────────────────────────────────────────

#' Family signatures, one per term
#'
#' Must agree with the `family_key` column of `genome_effect_terms`.
#' `center_value` is deliberately excluded (a line-specific variant legitimately
#' carries a different frequency from its common fallback), as are origin rows
#' and `genome_value` (they distinguish variants *within* a family).
#' `effect_owner` is included: owners always sum, so two owners defining the
#' same common term must both fire rather than collide as a tie.
#'
#' Vectorized over terms rather than called per term: every caller needs the key
#' for a whole model at once, and a per-term filter makes the validator and the
#' evaluator quadratic in a model with a QTL per locus.
#'
#' @param terms Data frame with `id_genome_effect`, `trait_name`, `effect_owner`.
#' @param members Data frame of member rows for those terms.
#' @return Character vector, one key per row of `terms`.
#' @keywords internal
#' @noRd
.ge_family_keys <- function(terms, members) {
  if (nrow(terms) == 0L) return(character(0))
  members <- members[order(members$id_genome_effect, members$locus_id), ,
                     drop = FALSE]
  piece <- paste(members$locus_id, members$contrast_name,
                 ifelse(is.na(members$copy_count_value), "-",
                        members$copy_count_value),
                 ifelse(is.na(members$dosage_value), "-", members$dosage_value),
                 sep = ":")
  sig <- tapply(piece, as.character(members$id_genome_effect),
                function(v) paste(v, collapse = "&"))
  paste0(terms$trait_name, "|", terms$effect_owner, "|",
         as.character(sig[as.character(terms$id_genome_effect)]))
}


# ── Origin predicates and containment ───────────────────────────────────────

#' Canonical line token for an origin row
#'
#' `'exact'` rows are identified by their line name, `'unknown'` and `'any'` by
#' the match type itself. Prefixed so a line literally named "unknown" cannot
#' masquerade as the unknown-line scope. Used on both dimensions of the lattice,
#' so an `'unknown'` demand inside a genotype multiset compares like any other.
#'
#' @keywords internal
#' @noRd
.ge_line_token <- function(line_match_type, line_name) {
  ifelse(line_match_type == "exact", paste0("exact:", line_name), line_match_type)
}

#' Build a variant's origin predicate: one per-member predicate, in slot order
#'
#' @param members,origins Data frames restricted to one `id_genome_effect`.
#' @return List of per-member predicates.
#' @keywords internal
#' @noRd
.ge_predicate <- function(members, origins) {
  members <- members[order(members$member_slot), , drop = FALSE]
  lapply(seq_len(nrow(members)), function(i) {
    slot <- members$member_slot[i]
    kind <- if (members$contrast_name[i] == "additive") "additive" else "genotype"
    oo <- origins[origins$member_slot == slot, , drop = FALSE]
    if (nrow(oo) == 0L) return(list(kind = kind, scope = "any"))
    oo <- oo[order(oo$origin_slot), , drop = FALSE]
    if (kind == "additive") {
      list(kind   = kind,
           scope  = "scoped",
           line   = .ge_line_token(oo$line_match_type[1], oo$line_name[1]),
           parent = oo$parent_origin[1])
    } else {
      items <- do.call(rbind, lapply(seq_len(nrow(oo)), function(k) {
        data.frame(line   = rep(.ge_line_token(oo$line_match_type[k],
                                               oo$line_name[k]),
                                oo$copy_count[k]),
                   parent = rep(oo$parent_origin[k], oo$copy_count[k]),
                   stringsAsFactors = FALSE)
      }))
      list(kind = kind, scope = "scoped", items = items)
    }
  })
}

#' Human-readable label for a variant's origin predicate
#'
#' Used in warnings and errors, where naming an `id_genome_effect` would tell a
#' user nothing about which of their calls produced the row.
#'
#' @keywords internal
#' @noRd
.ge_scope_label <- function(P) {
  paste(vapply(P, function(p) {
    if (identical(p$scope, "any")) return("common")
    if (p$kind == "additive") {
      paste0(p$line, if (is.na(p$parent)) "" else paste0("@p", p$parent))
    } else {
      paste(paste0(p$items$line,
                   ifelse(is.na(p$items$parent), "",
                          paste0("@p", p$items$parent))), collapse = "+")
    }
  }, character(1)), collapse = " x ")
}

#' Is predicate `P` contained in predicate `Q`? (componentwise, per member)
#'
#' `P <= Q` reads "P is at least as specific as Q", so the variant selected for
#' a tuple is the minimum among those that match it.
#'
#' @keywords internal
#' @noRd
.ge_pred_leq <- function(P, Q) {
  if (length(P) != length(Q)) return(FALSE)
  all(vapply(seq_along(P), function(i) .ge_member_leq(P[[i]], Q[[i]]), logical(1)))
}

#' @keywords internal
#' @noRd
.ge_member_leq <- function(p, q) {
  if (identical(q$scope, "any")) return(TRUE)
  if (identical(p$scope, "any")) return(FALSE)
  if (p$kind == "additive") {
    line_ok   <- identical(q$line, "any") || identical(p$line, q$line)
    parent_ok <- is.na(q$parent) || (!is.na(p$parent) && p$parent == q$parent)
    return(line_ok && parent_ok)
  }
  if (nrow(p$items) != nrow(q$items)) return(FALSE)
  if (!identical(sort(p$items$line), sort(q$items$line))) return(FALSE)
  .ge_bijection(p$items, q$items, mode = "refines")
}

#' Could some reachable copy label satisfy both predicates?
#'
#' Decided analytically, not by enumerating a population: a write must be
#' validated before any individual exists.
#'
#' @keywords internal
#' @noRd
.ge_pred_overlap <- function(P, Q) {
  if (length(P) != length(Q)) return(FALSE)
  all(vapply(seq_along(P), function(i) .ge_member_overlap(P[[i]], Q[[i]]),
             logical(1)))
}

#' @keywords internal
#' @noRd
.ge_member_overlap <- function(p, q) {
  if (identical(p$scope, "any") || identical(q$scope, "any")) return(TRUE)
  if (p$kind == "additive") {
    line_ok <- identical(p$line, "any") || identical(q$line, "any") ||
      identical(p$line, q$line)
    parent_ok <- is.na(p$parent) || is.na(q$parent) || p$parent == q$parent
    return(line_ok && parent_ok)
  }
  if (nrow(p$items) != nrow(q$items)) return(FALSE)
  if (!identical(sort(p$items$line), sort(q$items$line))) return(FALSE)
  .ge_bijection(p$items, q$items, mode = "compatible")
}

#' Exhaustive bijection between two item multisets
#'
#' Greedy consumption is **unsound** here: a demand set mixing a
#' parent-qualified row with an ANY-parent row (`{A:1@p1, A:1}`) can be
#' satisfiable while a greedy pass fails it, because the ANY row eats the copy
#' the qualified row needed. Established by the Phase A fixtures. At diploidy
#' the item lists are at most 2 long, so the exhaustive search is free.
#'
#' `mode = "refines"`: every `q` is ANY or equal to its partner (containment).
#' `mode = "compatible"`: either side may be ANY (overlap).
#'
#' @keywords internal
#' @noRd
.ge_bijection <- function(p_items, q_items, mode = c("refines", "compatible")) {
  mode <- match.arg(mode)
  if (nrow(p_items) == 0L) return(TRUE)
  d <- p_items[1, ]
  for (i in seq_len(nrow(q_items))) {
    q <- q_items[i, ]
    line_ok <- !is.na(q$line) && !is.na(d$line) && q$line == d$line
    parent_ok <- if (mode == "refines") {
      is.na(q$parent) || (!is.na(d$parent) && d$parent == q$parent)
    } else {
      is.na(q$parent) || is.na(d$parent) || d$parent == q$parent
    }
    if (line_ok && parent_ok &&
        .ge_bijection(p_items[-1, , drop = FALSE], q_items[-i, , drop = FALSE],
                      mode)) {
      return(TRUE)
    }
  }
  FALSE
}


# ── Validator ───────────────────────────────────────────────────────────────

#' Cross-row validation of the genome-effect tables
#'
#' Row-local invariants are declared SQL constraints in [define_genome()]. This
#' checks what SQL cannot: the additive-member origin invariant, `'any'`
#' restricted to additive members, exact-multiset satisfiability, and — within
#' each fallback family — no duplicate scope identity and no overlapping but
#' incomparable scopes.
#'
#' Called inside every genome-effect write transaction, before `COMMIT`.
#'
#' @param conn A DBI connection.
#' @param labels Optional `id_genome_effect` -> user label map; see
#'   `.ge_validate_frames()`.
#' @return `invisible(NULL)`; errors listing every violation found.
#' @keywords internal
#' @noRd
validate_genome_effects <- function(conn, labels = NULL) {
  terms   <- DBI::dbGetQuery(conn, "SELECT * FROM genome_effects")
  members <- DBI::dbGetQuery(conn, "SELECT * FROM genome_effect_members")
  origins <- DBI::dbGetQuery(conn, "SELECT * FROM genome_effect_member_origins")
  v <- c(.ge_validate_frames(terms, members, origins, labels = labels),
         .ge_validate_dominance_ploidy(conn, members, origins, labels = labels))
  if (length(v) > 0L) {
    stop("Invalid genome effects:\n  - ", paste(v, collapse = "\n  - "),
         call. = FALSE)
  }
  invisible(NULL)
}

#' Refuse a `dominance` member at a locus that is not reliably diploid
#'
#' Cockerham dominance coding is defined over the three diploid genotype states,
#' so a locus whose resolved `chr_inheritance` is not `1,1` for every applicable
#' offspring sex would have no value to give at evaluation time. The one escape
#' is an exact origin multiset **on that same member** demanding two copies —
#' a diploid-proving scope on a *different* member of the term proves nothing
#' about this locus.
#'
#' Needs the database (`chr_inheritance`, `genome_meta`), so it is separate from
#' the frame-level rules.
#'
#' @param conn A DBI connection.
#' @param members,origins Data frames of the effect member and origin rows.
#' @return Character vector of violations.
#' @keywords internal
#' @noRd
.ge_validate_dominance_ploidy <- function(conn, members, origins,
                                         labels = NULL) {
  tag_of <- function(id) {
    lab <- if (is.null(labels)) NA_character_ else
      unname(labels[match(as.character(id), names(labels))])
    if (length(lab) != 1L || is.na(lab)) paste0("term ", id)
    else paste0("term_id '", lab, "'")
  }
  dom <- members[members$contrast_name == "dominance", , drop = FALSE]
  if (nrow(dom) == 0L) return(character(0))

  chr_of <- DBI::dbGetQuery(conn, "SELECT locus_id, chr_name FROM genome_meta")
  inh <- lapply(c("M", "F"), function(sx) resolve_chr_inheritance(conn, sx))
  diploid <- Reduce(intersect, lapply(inh, function(d) {
    d$chr_name[d$from_parent_1 == 1L & d$from_parent_2 == 1L]
  }))

  v <- character(0)
  # Only members off a diploid chromosome need their origin rows.
  chrs <- chr_of$chr_name[match(dom$locus_id, chr_of$locus_id)]
  for (i in which(is.na(chrs) | !chrs %in% diploid)) {
    id   <- dom$id_genome_effect[i]
    slot <- dom$member_slot[i]
    chr  <- chrs[i]
    so <- origins[origins$id_genome_effect == id & origins$member_slot == slot,
                  , drop = FALSE]
    if (nrow(so) > 0L && sum(so$copy_count) == 2L) next
    v <- c(v, paste0(
      tag_of(id), " member ", slot, ": 'dominance' needs a diploid locus, but ",
      "chromosome '", chr, "' does not resolve to one copy from each parent for ",
      "both offspring sexes. Give this member an exact origin multiset ",
      "demanding two copies, or use 'indicator' states instead"))
  }
  v
}

#' Guidance for a duplicate caused only by a different centring
#'
#' `center_value` is deliberately outside the family signature -- a line-A
#' variant legitimately carries a different frequency from its common fallback
#' -- so a functional `additive`@0.5 and a Cockerham `additive`@p at one locus
#' under one owner collide as a duplicate. They are genuinely combinable:
#' `a1(g - c1) + a2(g - c2) = (a1 + a2)(g - c')` with
#' `c' = (a1 c1 + a2 c2) / (a1 + a2)`. Rejecting keeps the duplicate-term guard;
#' naming the single equivalent term keeps the rejection actionable.
#'
#' Order-one `additive` only: the Cockerham dominance contrast is not linear in
#' its centre, so no such combination exists for it.
#'
#' @keywords internal
#' @noRd
.ge_duplicate_hint <- function(terms, members, id1, id2) {
  m1 <- members[members$id_genome_effect == id1, , drop = FALSE]
  m2 <- members[members$id_genome_effect == id2, , drop = FALSE]
  if (nrow(m1) != 1L || nrow(m2) != 1L) return("")
  if (!identical(m1$contrast_name, "additive")) return("")
  c1 <- m1$center_value; c2 <- m2$center_value
  if (isTRUE(all.equal(c1, c2))) return("")
  a1 <- terms$genome_value[terms$id_genome_effect == id1]
  a2 <- terms$genome_value[terms$id_genome_effect == id2]
  if (isTRUE(all.equal(a1 + a2, 0))) {
    return(paste0(". They differ only in center_value (", c1, " vs ", c2,
                  "), but their coefficients cancel, so the combined term is ",
                  "the constant ", signif(a2 * (c2 - c1), 6),
                  " -- an intercept, which this model has no place for"))
  }
  paste0(". They differ only in center_value (", c1, " vs ", c2,
         "), which is outside the family signature on purpose. Write the one ",
         "equivalent term instead: genome_value = ", signif(a1 + a2, 6),
         ", center_value = ", signif((a1 * c1 + a2 * c2) / (a1 + a2), 6))
}

#' Validate genome-effect rows held as data frames
#'
#' Split from `validate_genome_effects()` so a writer can validate a candidate
#' set before issuing any SQL, and so the rules are testable without a database.
#'
#' @param terms,members,origins Data frames matching the three tables.
#' @param labels Optional named character vector mapping `id_genome_effect` to
#'   the label to print for it. A writer passes the `term_id` the user typed, so
#'   a rejected call never names a surrogate id the user has not seen.
#' @return Character vector of violations; empty when the rows are writable.
#' @keywords internal
#' @noRd
.ge_validate_frames <- function(terms, members, origins, labels = NULL) {
  tag_of <- function(id) {
    lab <- if (is.null(labels)) rep(NA_character_, length(id)) else
      unname(labels[match(as.character(id), names(labels))])
    ifelse(is.na(lab), paste0("term ", id), paste0("term_id '", lab, "'"))
  }
  v <- character(0)
  # Only a genuinely empty model is trivially valid. Members with no term at all
  # must still be reported, not skipped.
  if (nrow(terms) == 0L && nrow(members) == 0L && nrow(origins) == 0L) return(v)

  # Every member must belong to a term, and every origin row to a member.
  orphan_m <- setdiff(members$id_genome_effect, terms$id_genome_effect)
  if (length(orphan_m) > 0L) {
    v <- c(v, paste0("member rows with no term: ",
                     paste(sort(unique(orphan_m)), collapse = ", ")))
  }
  m_key <- paste(members$id_genome_effect, members$member_slot)
  o_key <- paste(origins$id_genome_effect, origins$member_slot)
  if (any(!o_key %in% m_key)) {
    v <- c(v, paste0("origin rows with no member: ",
                     paste(sort(unique(o_key[!o_key %in% m_key])), collapse = "; ")))
  }

  # Per-term rules, decided for every term at once: a per-term subset of the
  # member and origin frames is quadratic in the model size, and this runs on
  # the whole stored model before every COMMIT. The messages are then put in
  # the order a loop over terms (in row order) would produce them: term, then
  # rule, then member (in member row order), then the member's rules.
  v <- c(v, .ge_frame_term_messages(terms, members, origins, tag_of))

  # Fallback families: precedence operates only within a family.
  if (nrow(terms) == 0L) return(unique(v))
  keys  <- .ge_family_keys(terms, members)
  t_pos <- seq_len(nrow(terms))
  m_rows <- split(seq_len(nrow(members)),
                  factor(match(members$id_genome_effect, terms$id_genome_effect),
                         levels = t_pos))
  o_rows <- split(seq_len(nrow(origins)),
                  factor(match(origins$id_genome_effect, terms$id_genome_effect),
                         levels = t_pos))
  for (fam in split(t_pos, keys)) {
    if (length(fam) < 2L) next
    preds <- lapply(fam, function(i) {
      .ge_predicate(members[m_rows[[i]], , drop = FALSE],
                    origins[o_rows[[i]], , drop = FALSE])
    })
    ids <- terms$id_genome_effect[fam]
    for (i in seq_along(fam)) for (j in seq_along(fam)) {
      if (j <= i) next
      le <- .ge_pred_leq(preds[[i]], preds[[j]])
      ge <- .ge_pred_leq(preds[[j]], preds[[i]])
      if (le && ge) {
        v <- c(v, paste0(tag_of(ids[i]), " and ", tag_of(ids[j]),
                         " are the same logical term at the same scope",
                         " (duplicate family + scope identity)",
                         .ge_duplicate_hint(terms, members, ids[i], ids[j])))
      } else if (!le && !ge && .ge_pred_overlap(preds[[i]], preds[[j]])) {
        v <- c(v, paste0(tag_of(ids[i]), " and ", tag_of(ids[j]),
                         " have overlapping but incomparable scopes in one",
                         " family: neither is more specific, so no variant can",
                         " be selected"))
      }
    }
  }
  unique(v)
}

#' The per-term structural messages of `.ge_validate_frames()`, in loop order
#'
#' For each term (row order; a repeated id is reported at its first row):
#' "has no members" alone, or else a repeated locus, then the `member_slot`
#' rule (or, when the slots are valid, the locus order), then each member's
#' origin rules in member row order.
#'
#' @keywords internal
#' @noRd
.ge_frame_term_messages <- function(terms, members, origins, tag_of) {
  ids <- unique(terms$id_genome_effect)
  K   <- length(ids)
  if (K == 0L) return(character(0))
  tpos <- match(members$id_genome_effect, ids)
  keep <- which(!is.na(tpos))
  tp   <- tpos[keep]
  n_m  <- tabulate(tp, K)
  tag  <- tag_of(ids)
  out_tp <- integer(0); out_rule <- integer(0); out_m <- integer(0)
  out_msg <- character(0)
  # paste() recycles a zero-length argument to "", so an empty selection
  # still yields one string: append only when something was selected.
  add <- function(t, rule, m, msg) {
    if (length(t) == 0L) return(invisible(NULL))
    out_tp   <<- c(out_tp, t);    out_rule <<- c(out_rule, rule)
    out_m    <<- c(out_m, m);     out_msg  <<- c(out_msg, msg)
  }

  none <- which(n_m == 0L)
  add(none, rep(1L, length(none)), rep(0L, length(none)),
      paste(tag[none], "has no members"))

  locus <- members$locus_id[keep]
  dup_locus <- tabulate(tp[duplicated(.ge_pair_key(tp, locus))], K) > 0L
  w <- which(dup_locus)
  add(w, rep(2L, length(w)), rep(0L, length(w)),
      paste(tag[w], "names the same locus more than once"))

  # identical(sort(slots), seq_len(n)): integer slots that are exactly 1..n.
  slot <- members$member_slot[keep]
  bad_slot_row <- is.na(slot) | slot < 1L | slot > n_m[tp] |
    duplicated(.ge_pair_key(tp, slot))
  slot_bad <- if (!is.integer(members$member_slot)) n_m > 0L else
    tabulate(tp[bad_slot_row], K) > 0L
  w <- which(slot_bad)
  add(w, rep(3L, length(w)), rep(0L, length(w)),
      paste(tag[w], "member_slot must be 1..n in ascending locus_id order"))
  # Valid slots: are the loci ascending in slot order?
  o  <- order(tp, slot)
  so <- tp[o]; lo <- locus[o]
  down <- c(FALSE, so[-1] == so[-length(so)] & lo[-1] < lo[-length(lo)])
  unsorted <- tabulate(so[down], K) > 0L & !slot_bad
  w <- which(unsorted)
  add(w, rep(3L, length(w)), rep(0L, length(w)),
      paste(tag[w], "members are not canonicalized by ascending locus_id"))

  # Member origin rules. Each member reads the origin rows of its own
  # (term, slot), as a per-term subset would.
  if (nrow(origins) > 0L && length(keep) > 0L) {
    okey <- paste(origins$id_genome_effect, origins$member_slot)
    ukey <- unique(okey)
    g    <- match(okey, ukey)
    G    <- length(ukey)
    n_o  <- tabulate(g, G)
    os   <- origins$origin_slot
    os_bad <- if (!is.integer(os)) n_o > 0L else
      tabulate(g[is.na(os) | os < 1L | os > n_o[g] |
                   duplicated(.ge_pair_key(g, os))], G) > 0L
    cc <- origins$copy_count
    cc_ne1  <- tabulate(g[is.na(cc) | cc != 1L], G) > 0L
    any_any <- tabulate(g[origins$line_match_type %in% "any"], G) > 0L
    sum_cc  <- as.numeric(rowsum(as.numeric(cc), factor(g, levels = seq_len(G)),
                                 reorder = TRUE))

    mg   <- match(paste(members$id_genome_effect[keep], slot), ukey)
    has  <- which(!is.na(mg))
    if (length(has) > 0L) {
      r    <- keep[has]                       # member rows (original order)
      gm   <- mg[has]
      kind <- members$contrast_name[r]
      mtag <- paste0(tag[tp[has]], " member ", slot[has])
      sub1 <- ifelse(os_bad[gm], paste(mtag, "origin_slot must be 1..n"),
                     NA_character_)
      add_k <- kind == "additive"
      want <- ifelse(kind == "dominance", 2L, members$copy_count_value[r])
      sub2 <- ifelse(add_k,
        ifelse(n_o[gm] > 1L,
               paste(mtag, "additive members take at most one origin row",
                     "(expand alternatives into separate variants)"),
               NA_character_),
        ifelse(any_any[gm],
               paste(mtag, "'any' is permitted only on additive members:",
                     "a genotype scope must name lines and sum to the",
                     "realized copy count, so an 'any' entry constrains",
                     "nothing"),
               NA_character_))
      sub3 <- ifelse(add_k,
        ifelse(cc_ne1[gm],
               paste(mtag, "additive origin requires copy_count = 1",
                     "(matching is per allele copy)"),
               NA_character_),
        ifelse(!is.na(want) & sum_cc[gm] != want,
               paste0(mtag, ": origin multiset demands ", sum_cc[gm],
                      " copies but the member's state is defined over ",
                      want, ifelse(kind == "dominance",
                                   " (Cockerham dominance is diploid)",
                                   " (copy_count_value)")),
               NA_character_))
      msg <- c(rbind(sub1, sub2, sub3))
      m_pos <- rep(r, each = 3L)
      t_of  <- rep(tp[has], each = 3L)
      ok <- !is.na(msg)
      add(t_of[ok], rep(4L, sum(ok)), m_pos[ok], msg[ok])
    }
  }

  ord <- order(out_tp, out_rule, out_m, seq_along(out_tp))
  out_msg[ord]
}

# ── Dosage collection ───────────────────────────────────────────────────────

#' Refuse an `n x m` dosage matrix above `QTL_REALISED_MAX_CELLS`
#'
#' Separate from the collector so a caller can check a knowable size before
#' any expensive work (`extract_genetic_variance()` checks before evaluating).
#'
#' @keywords internal
#' @noRd
.dosage_guard <- function(n, m, who, size_fix) {
  if (as.numeric(n) * m > QTL_REALISED_MAX_CELLS) {
    stop(who, " would collect a ", n, " x ", m, " genotype matrix (",
         format(as.numeric(n) * m, big.mark = ","), " cells), above the ",
         "limit of ", format(QTL_REALISED_MAX_CELLS, big.mark = ",",
                             scientific = FALSE), ". ", size_fix,
         call. = FALSE)
  }
  invisible(NULL)
}

#' Collect dosages, `id_ind` then `locus_id`, for a set of individuals
#'
#' Integer sums only (exact). One implementation for every caller that needs
#' an in-memory genotype matrix (the realised anchor of
#' `define_additive_effects()`, `extract_genetic_variance()`); each caller
#' passes its own wording and fix. Refuses fewer than 2 individuals, a matrix
#' above `QTL_REALISED_MAX_CELLS`, and any individual without both allele
#' copies at every requested locus.
#'
#' @param conn A DBI connection.
#' @param ids_sql SQL returning one column `id_ind` (distinct individuals).
#'   Individual ids never appear in it as literals.
#' @param locus_ids Integer locus ids.
#' @param who,where,size_fix,incomplete Message pieces: the caller, what
#'   selected the individuals, the fix for a too-large matrix, and the whole
#'   explanation for incomplete genotypes.
#' @return list(X = n x m matrix with `id_ind` row names, id_ind, locus_id),
#'   rows in `id_ind` order and columns in `locus_id` order.
#' @keywords internal
#' @noRd
.collect_dosages <- function(conn, ids_sql, locus_ids, who, where, size_fix,
                             incomplete) {
  n <- DBI::dbGetQuery(conn, paste0("SELECT COUNT(*) AS n FROM (", ids_sql,
                                    ")"))$n
  m <- length(locus_ids)
  if (n < 2L) {
    stop(who, " needs at least 2 individuals in ", where, "; it selects ", n,
         ".", call. = FALSE)
  }
  .dosage_guard(n, m, who, size_fix)
  locus_ids <- sort(as.integer(locus_ids))
  lst <- paste(locus_ids, collapse = ", ")
  d <- DBI::dbGetQuery(conn, paste0(
    "SELECT h.id_ind, h.locus_id, CAST(SUM(h.allele) AS INTEGER) AS dosage, ",
    "COUNT(*) AS n_copies FROM ind_haplotype h ",
    "JOIN (", ids_sql, ") ids USING (id_ind) ",
    "WHERE h.locus_id IN (", lst, ") ",
    "GROUP BY h.id_ind, h.locus_id ORDER BY h.id_ind, h.locus_id"))
  per_ind <- table(factor(d$id_ind))
  # The SQL order, not R's sort(): row names must name the rows they label,
  # whatever the session's collation.
  ids <- unique(d$id_ind)
  if (length(ids) < n || any(per_ind != m) || any(d$n_copies != 2L)) {
    stop(who, ": ", incomplete, call. = FALSE)
  }
  X <- matrix(as.numeric(d$dosage), nrow = length(ids), ncol = m, byrow = TRUE,
              dimnames = list(ids, NULL))
  list(X = X, id_ind = ids, locus_id = locus_ids)
}
