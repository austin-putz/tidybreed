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
    "       m.locus_id, g.locus_name, m.contrast_name ",
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
#' @keywords internal
#' @noRd
.ind_tgv_total_view_sql <- function() {
  paste0(
    "CREATE VIEW ind_tgv_total AS ",
    "SELECT id_ind, trait_name, SUM(tgv_value) AS tgv_total ",
    "FROM ind_tgv GROUP BY id_ind, trait_name"
  )
}


# ── Family signature ────────────────────────────────────────────────────────

#' Family signature for one term, from its member rows
#'
#' Must agree with the `family_key` column of `genome_effect_terms`.
#' `center_value` is deliberately excluded (a line-specific variant legitimately
#' carries a different frequency from its common fallback), as are origin rows
#' and `genome_value` (they distinguish variants *within* a family).
#' `effect_owner` is included: owners always sum, so two owners defining the
#' same common term must both fire rather than collide as a tie.
#'
#' @param trait_name,effect_owner Character scalars.
#' @param members Data frame of that term's member rows.
#' @return Character scalar.
#' @keywords internal
#' @noRd
.ge_family_key <- function(trait_name, effect_owner, members) {
  members <- members[order(members$locus_id), , drop = FALSE]
  paste0(
    trait_name, "|", effect_owner, "|",
    paste(members$locus_id, members$contrast_name,
          ifelse(is.na(members$copy_count_value), "-", members$copy_count_value),
          ifelse(is.na(members$dosage_value), "-", members$dosage_value),
          sep = ":", collapse = "&")
  )
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
#' @return `invisible(NULL)`; errors listing every violation found.
#' @keywords internal
#' @noRd
validate_genome_effects <- function(conn) {
  terms   <- DBI::dbGetQuery(conn, "SELECT * FROM genome_effects")
  members <- DBI::dbGetQuery(conn, "SELECT * FROM genome_effect_members")
  origins <- DBI::dbGetQuery(conn, "SELECT * FROM genome_effect_member_origins")
  v <- c(.ge_validate_frames(terms, members, origins),
         .ge_validate_dominance_ploidy(conn, members, origins))
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
.ge_validate_dominance_ploidy <- function(conn, members, origins) {
  dom <- members[members$contrast_name == "dominance", , drop = FALSE]
  if (nrow(dom) == 0L) return(character(0))

  chr_of <- DBI::dbGetQuery(conn, "SELECT locus_id, chr_name FROM genome_meta")
  inh <- lapply(c("M", "F"), function(sx) resolve_chr_inheritance(conn, sx))
  diploid <- Reduce(intersect, lapply(inh, function(d) {
    d$chr_name[d$from_parent_1 == 1L & d$from_parent_2 == 1L]
  }))

  v <- character(0)
  for (i in seq_len(nrow(dom))) {
    id   <- dom$id_genome_effect[i]
    slot <- dom$member_slot[i]
    chr  <- chr_of$chr_name[match(dom$locus_id[i], chr_of$locus_id)]
    if (!is.na(chr) && chr %in% diploid) next
    so <- origins[origins$id_genome_effect == id & origins$member_slot == slot,
                  , drop = FALSE]
    if (nrow(so) > 0L && sum(so$copy_count) == 2L) next
    v <- c(v, paste0(
      "term ", id, " member ", slot, ": 'dominance' needs a diploid locus, but ",
      "chromosome '", chr, "' does not resolve to one copy from each parent for ",
      "both offspring sexes. Give this member an exact origin multiset ",
      "demanding two copies, or use 'indicator' states instead"))
  }
  v
}

#' Validate genome-effect rows held as data frames
#'
#' Split from [validate_genome_effects()] so a writer can validate a candidate
#' set before issuing any SQL, and so the rules are testable without a database.
#'
#' @param terms,members,origins Data frames matching the three tables.
#' @return Character vector of violations; empty when the rows are writable.
#' @keywords internal
#' @noRd
.ge_validate_frames <- function(terms, members, origins) {
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

  for (id in terms$id_genome_effect) {
    mm <- members[members$id_genome_effect == id, , drop = FALSE]
    oo <- origins[origins$id_genome_effect == id, , drop = FALSE]
    tag <- paste0("term ", id)

    if (nrow(mm) == 0L) {
      v <- c(v, paste(tag, "has no members"))
      next
    }
    if (anyDuplicated(mm$locus_id)) {
      v <- c(v, paste(tag, "names the same locus more than once"))
    }
    if (!identical(sort(mm$member_slot), seq_len(nrow(mm)))) {
      v <- c(v, paste(tag, "member_slot must be 1..n in ascending locus_id order"))
    } else if (is.unsorted(mm$locus_id[order(mm$member_slot)])) {
      v <- c(v, paste(tag, "members are not canonicalized by ascending locus_id"))
    }

    for (i in seq_len(nrow(mm))) {
      slot <- mm$member_slot[i]
      kind <- mm$contrast_name[i]
      so   <- oo[oo$member_slot == slot, , drop = FALSE]
      mtag <- paste0(tag, " member ", slot)
      if (nrow(so) == 0L) next
      if (!identical(sort(so$origin_slot), seq_len(nrow(so)))) {
        v <- c(v, paste(mtag, "origin_slot must be 1..n"))
      }
      if (kind == "additive") {
        # Two origin rows on one additive member are incoherent: under OR the
        # scope matches more copies than either alone while looking narrower;
        # under AND no single copy can be both. "A or B" is expanded by the
        # writer into separate variants.
        if (nrow(so) > 1L) {
          v <- c(v, paste(mtag, "additive members take at most one origin row",
                          "(expand alternatives into separate variants)"))
        }
        if (any(so$copy_count != 1L)) {
          v <- c(v, paste(mtag, "additive origin requires copy_count = 1",
                          "(matching is per allele copy)"))
        }
      } else {
        if (any(so$line_match_type == "any")) {
          v <- c(v, paste(mtag, "'any' is permitted only on additive members:",
                          "a genotype scope must name lines and sum to the",
                          "realized copy count, so an 'any' entry constrains",
                          "nothing"))
        }
        # A genotype unit is the whole locus state, so the multiset must account
        # for every copy the state is defined over: 2 for Cockerham dominance,
        # and the declared copy_count_value for an indicator.
        want <- if (kind == "dominance") 2L else mm$copy_count_value[i]
        if (!is.na(want) && sum(so$copy_count) != want) {
          v <- c(v, paste0(mtag, ": origin multiset demands ", sum(so$copy_count),
                           " copies but the member's state is defined over ",
                           want, if (kind == "dominance")
                             " (Cockerham dominance is diploid)" else
                             " (copy_count_value)"))
        }
      }
    }
  }

  # Fallback families: precedence operates only within a family.
  if (nrow(terms) == 0L) return(unique(v))
  keys <- vapply(terms$id_genome_effect, function(id) {
    .ge_family_key(terms$trait_name[terms$id_genome_effect == id],
                   terms$effect_owner[terms$id_genome_effect == id],
                   members[members$id_genome_effect == id, , drop = FALSE])
  }, character(1))
  for (fam in split(terms$id_genome_effect, keys)) {
    if (length(fam) < 2L) next
    preds <- lapply(fam, function(id) {
      .ge_predicate(members[members$id_genome_effect == id, , drop = FALSE],
                    origins[origins$id_genome_effect == id, , drop = FALSE])
    })
    for (i in seq_along(fam)) for (j in seq_along(fam)) {
      if (j <= i) next
      le <- .ge_pred_leq(preds[[i]], preds[[j]])
      ge <- .ge_pred_leq(preds[[j]], preds[[i]])
      if (le && ge) {
        v <- c(v, paste0("terms ", fam[i], " and ", fam[j],
                         " are the same logical term at the same scope",
                         " (duplicate family + scope identity)"))
      } else if (!le && !ge && .ge_pred_overlap(preds[[i]], preds[[j]])) {
        v <- c(v, paste0("terms ", fam[i], " and ", fam[j],
                         " have overlapping but incomparable scopes in one",
                         " family: neither is more specific, so no variant can",
                         " be selected"))
      }
    }
  }
  unique(v)
}
