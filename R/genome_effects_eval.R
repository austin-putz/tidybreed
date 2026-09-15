# ---------------------------------------------------------------------------
# The genome-effect evaluator.
#
# One evaluator serves every consumer: add_tgv() writes all of its components to
# ind_tgv, and add_tbv() is a filtered call into it (reserved owner, order-one
# `additive` terms only). There is deliberately no second implementation of the
# effect math anywhere in the package.
#
# Built to the "Evaluation strategy" section of plans/update_genome_effects_v4.md
# (v4.9). The semantic definition is per evaluation tuple -- one unit from each
# member, one variant selected per tuple by predicate containment -- but it is
# never executed that way. Two facts collapse it into set-based SQL:
#
#   1. An origin predicate reads a copy's (line_origin, parent_origin) *label*
#      and nothing else, so the winning variant is a function of the label, not
#      of who carries it.
#   2. Tuples therefore group by label-vector, and the inner sum factors inside
#      each group.
#
# So resolution runs once over a table of size #families x prod(|labels_m|) and
# never during evaluation, and the evaluation itself is a fixed number of SQL
# statements whatever the number of individuals (gates 51-52).
#
# The reference implementations this must agree with are the two independent
# evaluators in tests/testthat/helper-genome-effects.R, written in Phase A
# before any of this existed.
# ---------------------------------------------------------------------------


#' Model-structure component for one term
#'
#' These are **declared model structure, not variance components**: a functional
#' A x A term contributes to A, D *and* I in the statistical sense, so the names
#' carry the order and can never be read as `V_A` / `V_D` / `V_I`.
#'
#' @keywords internal
#' @noRd
.gev_component <- function(contrast_name, effect_order) {
  unname(ifelse(effect_order > 1L, "interaction",
                c(additive  = "order1_additive",
                  dominance = "order1_dominance",
                  indicator = "order1_other")[contrast_name]))
}

#' Evaluation-unit kind for a member contrast
#'
#' `additive` matches per allele copy; `dominance` and `indicator` match the
#' locus local state, which is one unit whether or not any copy exists.
#'
#' @keywords internal
#' @noRd
.gev_member_kind <- function(contrast_name) {
  ifelse(contrast_name == "additive", "additive", "genotype")
}


# -- Reading the model ------------------------------------------------------

#' Read the stored effect model for a set of traits
#'
#' @param conn A DBI connection.
#' @param trait_names Character vector.
#' @param effect_owner Optional owner filter (a single reserved or custom name).
#' @return List of three data frames, `terms`, `members`, `origins`.
#' @keywords internal
#' @noRd
.gev_read_model <- function(conn, trait_names, effect_owner = NULL) {
  terms <- DBI::dbGetQuery(conn, paste0(
    "SELECT id_genome_effect, trait_name, effect_owner, effect_name, ",
    "       genome_value FROM genome_effects ",
    "WHERE trait_name IN (", sql_in_list(trait_names, "trait name"), ")",
    if (is.null(effect_owner)) "" else
      paste0(" AND effect_owner = ", sql_in_list(effect_owner, "effect owner")),
    " ORDER BY id_genome_effect"))

  if (nrow(terms) == 0L) {
    terms$effect_order   <- integer(0)
    terms$component_name <- character(0)
    terms$family_key     <- character(0)
    return(list(terms = terms,
                members = .gev_empty_members(),
                origins = .gev_empty_origins()))
  }

  ids <- terms$id_genome_effect
  members <- DBI::dbGetQuery(conn, paste0(
    "SELECT id_genome_effect, member_slot, locus_id, contrast_name, ",
    "       copy_count_value, dosage_value, center_value ",
    "FROM genome_effect_members WHERE id_genome_effect IN (",
    paste(ids, collapse = ", "), ") ORDER BY id_genome_effect, member_slot"))
  origins <- DBI::dbGetQuery(conn, paste0(
    "SELECT id_genome_effect, member_slot, origin_slot, line_match_type, ",
    "       line_name, parent_origin, copy_count ",
    "FROM genome_effect_member_origins WHERE id_genome_effect IN (",
    paste(ids, collapse = ", "), ") ORDER BY id_genome_effect, member_slot, origin_slot"))

  ord <- tabulate(match(members$id_genome_effect, ids), nbins = length(ids))
  terms$effect_order <- ord

  m1 <- members[!duplicated(members$id_genome_effect), , drop = FALSE]
  terms$component_name <- .gev_component(
    m1$contrast_name[match(terms$id_genome_effect, m1$id_genome_effect)],
    terms$effect_order)

  terms$family_key <- .ge_family_keys(terms, members)

  list(terms = terms, members = members, origins = origins)
}

#' @keywords internal
#' @noRd
.gev_empty_members <- function() {
  data.frame(id_genome_effect = integer(0), member_slot = integer(0),
             locus_id = integer(0), contrast_name = character(0),
             copy_count_value = integer(0), dosage_value = integer(0),
             center_value = numeric(0), stringsAsFactors = FALSE)
}

#' @keywords internal
#' @noRd
.gev_empty_origins <- function() {
  data.frame(id_genome_effect = integer(0), member_slot = integer(0),
             origin_slot = integer(0), line_match_type = character(0),
             line_name = character(0), parent_origin = integer(0),
             copy_count = integer(0), stringsAsFactors = FALSE)
}


#' The subset of a model that `add_tbv()` reads
#'
#' Order-one terms whose single member's contrast is `additive`, under the
#' reserved owner. One definition, because "what is a breeding-value
#' coefficient" is a contract rather than a convenience filter: an additive
#' member sitting inside an interaction is not one, and neither is a coefficient
#' a user hand-wrote under their own owner.
#'
#' @keywords internal
#' @noRd
.gev_reserved_additive <- function(model) {
  keep <- model$terms$effect_order == 1L &
    model$terms$effect_owner == GE_ADDITIVE_OWNER &
    model$terms$id_genome_effect %in%
      model$members$id_genome_effect[model$members$contrast_name == "additive"]
  terms <- model$terms[keep, , drop = FALSE]
  list(terms   = terms,
       members = model$members[model$members$id_genome_effect %in%
                                 terms$id_genome_effect, , drop = FALSE],
       origins = model$origins[model$origins$id_genome_effect %in%
                                 terms$id_genome_effect, , drop = FALSE])
}


# -- Label alphabets --------------------------------------------------------

# A copy's label is `line_origin@parent_origin`, with the empty string standing
# for an unknown line -- which is unambiguous because `add_founders()` requires
# `^[a-zA-Z][a-zA-Z0-9_]*$` of every line name. A genotype label is the sorted
# multiset of the copy labels at that locus, comma-joined; the empty string is
# the zero-copy state, which has no haplotype row and must be materialized
# rather than joined away.

#' Copy-label alphabet: `SELECT DISTINCT line_origin, parent_origin`
#'
#' @param conn A DBI connection.
#' @param ind_tmp Name of a registered one-column frame of `id_ind`.
#' @return Data frame with `line` (may be `NA`), `parent`, `label`.
#' @keywords internal
#' @noRd
.gev_additive_alphabet <- function(conn, ind_tmp) {
  lab <- DBI::dbGetQuery(conn, paste0(
    "SELECT DISTINCT h.line_origin AS line, h.parent_origin AS parent ",
    "FROM ind_haplotype h JOIN ", ind_tmp, " i ON i.id_ind = h.id_ind ",
    "ORDER BY 1, 2"))
  if (nrow(lab) == 0L) {
    return(data.frame(line = character(0), parent = integer(0),
                      label = character(0), stringsAsFactors = FALSE))
  }
  bad <- lab$line[!is.na(lab$line) & grepl("[,@]", lab$line)]
  if (length(bad) > 0L) {
    stop("Line name(s) ", paste(unique(bad), collapse = ", "),
         " contain ',' or '@', which the evaluator uses to encode origin ",
         "labels. Line names must match ^[a-zA-Z][a-zA-Z0-9_]*$.", call. = FALSE)
  }
  lab$label <- paste0(ifelse(is.na(lab$line), "", lab$line), "@", lab$parent)
  lab
}

#' Genotype-label alphabet: the distinct copy multisets that actually occur
#'
#' Computed over the loci that genotype members name, plus the zero-copy state.
#' The left join is against every such locus rather than against what
#' `chr_inheritance` says an individual should carry: a `(copy_count 0,
#' dosage 0)` indicator is a "carries no copy here" effect and must match every
#' individual with no row, whatever the reason.
#'
#' @return Data frame with `label`; the items are recovered by
#'   `.gev_parse_label()`.
#' @keywords internal
#' @noRd
.gev_genotype_alphabet <- function(conn, ind_tmp, locus_ids) {
  if (length(locus_ids) == 0L) return(data.frame(label = character(0)))
  labs <- DBI::dbGetQuery(conn, paste0(
    "WITH st AS ( ",
    "  SELECT i.id_ind, g.locus_id, ",
    "         COALESCE(string_agg(COALESCE(h.line_origin, '') || '@' || ",
    "                             h.parent_origin, ',' ",
    "                  ORDER BY COALESCE(h.line_origin, ''), h.parent_origin), ",
    "                  '') AS label ",
    "  FROM ", ind_tmp, " i ",
    "  CROSS JOIN (SELECT UNNEST([", paste(locus_ids, collapse = ", "),
    "]) AS locus_id) g ",
    "  LEFT JOIN ind_haplotype h ",
    "    ON h.id_ind = i.id_ind AND h.locus_id = g.locus_id ",
    "  GROUP BY 1, 2 ) ",
    "SELECT DISTINCT label FROM st ORDER BY label"))
  labs
}

#' Recover a label's items: one `(line, parent)` per copy
#'
#' @keywords internal
#' @noRd
.gev_parse_label <- function(label) {
  if (!nzchar(label)) {
    return(data.frame(line = character(0), parent = integer(0),
                      stringsAsFactors = FALSE))
  }
  parts <- strsplit(strsplit(label, ",", fixed = TRUE)[[1]], "@", fixed = TRUE)
  data.frame(
    line   = vapply(parts, function(p) if (nzchar(p[1])) p[1] else NA_character_,
                    character(1)),
    parent = vapply(parts, function(p) as.integer(p[2]), integer(1)),
    stringsAsFactors = FALSE)
}


# -- Predicate matching -----------------------------------------------------

#' Does a per-member predicate match a unit's label?
#'
#' The containment lattice itself lives in `.ge_pred_leq()`; this is the other
#' half, asking whether a predicate *fires* on a concrete label rather than how
#' two predicates order.
#'
#' @keywords internal
#' @noRd
.gev_member_match <- function(p, unit) {
  if (identical(p$scope, "any")) return(TRUE)
  if (p$kind == "additive") {
    line_ok <- if (identical(p$line, "any")) TRUE
               else if (identical(p$line, "unknown")) is.na(unit$line)
               else !is.na(unit$line) && identical(p$line, paste0("exact:", unit$line))
    parent_ok <- is.na(p$parent) || unit$parent == p$parent
    return(line_ok && parent_ok)
  }
  .gev_match_items(p$items, unit$items)
}

#' Exact-multiset satisfaction: every demand consumed by a distinct copy
#'
#' Exhaustive rather than greedy, for the reason `.ge_bijection()` records: a
#' demand set mixing an ANY-parent row with a parent-qualified one can be
#' satisfiable while a greedy pass fails it.
#'
#' @keywords internal
#' @noRd
.gev_match_items <- function(demands, copies) {
  if (nrow(demands) != nrow(copies)) return(FALSE)
  if (nrow(demands) == 0L) return(TRUE)
  d <- demands[1, ]
  for (i in seq_len(nrow(copies))) {
    cc <- copies[i, ]
    line_ok <- if (identical(d$line, "unknown")) is.na(cc$line)
               else !is.na(cc$line) && identical(d$line, paste0("exact:", cc$line))
    parent_ok <- is.na(d$parent) || cc$parent == d$parent
    if (line_ok && parent_ok &&
        .gev_match_items(demands[-1, , drop = FALSE],
                         copies[-i, , drop = FALSE])) {
      return(TRUE)
    }
  }
  FALSE
}


# -- Which members need a label at all --------------------------------------

#' Per-family, which member slots no variant scopes
#'
#' A member that **no** variant in its family scopes has predicate `any` for
#' every variant, so the selected variant cannot depend on that member's label.
#' The inner sum over its units then runs over all of them at once:
#'
#' ```text
#'   sum over L_m sum over units carrying L_m  =  sum over all units
#' ```
#'
#' which is the plan's fast path stated per member rather than per family. It
#' matters: without it, a 50-locus unscoped dominance term would enumerate
#' `|labels|^50` label-vectors and trip the resource guard, even though
#' evaluating it is one `GROUP BY`. Such a member carries the sentinel label
#' `"*"`, which no real label can collide with because every real label
#' contains an `@`.
#'
#' @return Logical vector, one element per row of `model$members`.
#' @keywords internal
#' @noRd
.gev_slot_freedom <- function(model) {
  scoped <- paste(model$origins$id_genome_effect, model$origins$member_slot)
  fam_of <- model$terms$family_key[match(model$members$id_genome_effect,
                                         model$terms$id_genome_effect)]
  key <- paste(fam_of, model$members$member_slot)
  any_scoped <- paste(model$members$id_genome_effect,
                      model$members$member_slot) %in% scoped
  # A slot is free when no variant of its family scopes that slot.
  !(key %in% unique(key[any_scoped]))
}


# -- The resolved variant map -----------------------------------------------

#' Canonical signature of a family's resolution problem
#'
#' Resolution depends on the member kinds and the variants' origin rows —
#' **never on which loci the members name.** A 500-QTL model with a common
#' variant and two line-specific ones is therefore 500 copies of one problem,
#' and solving it once and reusing the answer is the difference between a fixed
#' cost of seconds and of milliseconds. Locus ids are deliberately absent from
#' the key, and the key is built from the raw origin rows rather than from
#' built predicates so that a cache hit costs no predicate construction either.
#'
#' @param kinds,use_star Per-member vectors, in slot order.
#' @param org_sig Per-variant origin signature, in the family's variant order.
#' @keywords internal
#' @noRd
.gev_family_sig <- function(kinds, use_star, org_sig) {
  paste0(paste(kinds, collapse = ","), "|",
         paste(as.integer(use_star), collapse = ","), "|",
         paste(org_sig, collapse = ";"))
}

#' Per-term signature of its origin rows, order-insensitive
#'
#' @keywords internal
#' @noRd
.gev_origin_sigs <- function(origins, ids) {
  out <- stats::setNames(rep("", length(ids)), as.character(ids))
  if (nrow(origins) == 0L) return(out)
  row <- paste(origins$member_slot, origins$line_match_type, origins$line_name,
               origins$parent_origin, origins$copy_count, sep = ":")
  agg <- tapply(row, as.character(origins$id_genome_effect),
                function(v) paste(sort(v), collapse = "|"))
  out[names(agg)] <- as.character(agg)
  out
}

#' Resolve one variant per (family, label-vector), once, before evaluation
#'
#' Containment search runs **only** here, over a table whose size is
#' `#families x prod(|labels_m|)`. The map is total and unambiguous by
#' construction, because overlapping-but-incomparable predicates are refused at
#' write time -- so no tie can reach this point, and one arriving is a bug in
#' the writer rather than a user error.
#'
#' @return Long map: one row per `(map_id, member_slot)`, carrying the winning
#'   `id_genome_effect`, that slot's `label`, and the term's member count.
#' @keywords internal
#' @noRd
.gev_variant_map <- function(model, alphabets, free) {
  terms   <- model$terms
  members <- model$members
  ids_by  <- split(terms$id_genome_effect, terms$family_key)
  key_of  <- function(id) as.character(id)

  # Column-wise splits, not data-frame subsetting: the loop below runs once per
  # family and a 500-QTL model has 500 of them.
  mid      <- as.character(members$id_genome_effect)
  slot_by  <- split(members$member_slot, mid)
  kind_by  <- split(.gev_member_kind(members$contrast_name), mid)
  star_by  <- split(free, mid)
  org_sig  <- .gev_origin_sigs(model$origins, terms$id_genome_effect)
  mem_by   <- NULL   # built lazily: only a cache miss needs the full rows
  org_by   <- NULL

  acc_map <- acc_id <- acc_slot <- acc_lab <- acc_n <- list()
  n_out <- 0L
  map_id <- 0L
  cache <- new.env(parent = emptyenv())

  for (key in names(ids_by)) {
    fam_ids <- ids_by[[key]]
    k1      <- key_of(fam_ids[1])
    slots   <- slot_by[[k1]]
    kinds   <- kind_by[[k1]]
    stars   <- star_by[[k1]]

    labs <- lapply(seq_along(kinds), function(m) {
      if (stars[m]) "*" else alphabets[[kinds[m]]]
    })
    if (any(lengths(labs) == 0L)) next

    sigs <- org_sig[key_of(fam_ids)]
    if (length(fam_ids) == 1L && !nzchar(sigs[1])) {
      # The common scope matches every label-vector: no containment search.
      # This is the overwhelmingly common family shape.
      picks <- rep(fam_ids, prod(lengths(labs)))
      vecs  <- .gev_label_grid(labs)
    } else {
      vecs <- .gev_label_grid(labs)
      sig  <- .gev_family_sig(kinds, stars, sigs)
      idx  <- cache[[sig]]
      if (is.null(idx)) {
        if (is.null(mem_by)) {
          mem_by <- split(members, mid)
          org_by <- split(model$origins, as.character(model$origins$id_genome_effect))
        }
        preds <- lapply(key_of(fam_ids), function(id) {
          oo <- org_by[[id]]
          .ge_predicate(mem_by[[id]], if (is.null(oo)) .gev_empty_origins() else oo)
        })
        idx <- vapply(seq_len(nrow(vecs)), function(r) {
          units <- lapply(seq_along(labs), function(m) {
            if (identical(vecs[r, m], "*")) return(list(line = NA, parent = NA))
            it <- .gev_parse_label(vecs[r, m])
            if (kinds[m] == "additive") list(line = it$line[1], parent = it$parent[1])
            else list(items = it)
          })
          .gev_select_variant(preds, units, fam_ids, key)
        }, integer(1))
        cache[[sig]] <- idx
      }
      picks <- ifelse(is.na(idx), NA_integer_, fam_ids[idx])
    }

    keep <- which(!is.na(picks))
    if (length(keep) == 0L) next
    n_mem <- length(slots)
    n_out <- n_out + 1L
    acc_map[[n_out]]  <- rep(map_id + seq_along(keep), each = n_mem)
    acc_id[[n_out]]   <- rep(picks[keep], each = n_mem)
    acc_slot[[n_out]] <- rep(slots, times = length(keep))
    acc_lab[[n_out]]  <- as.vector(t(vecs[keep, , drop = FALSE]))
    acc_n[[n_out]]    <- rep(n_mem, n_mem * length(keep))
    map_id <- map_id + length(keep)
  }

  if (n_out == 0L) {
    return(data.frame(map_id = integer(0), id_genome_effect = integer(0),
                      member_slot = integer(0), label = character(0),
                      n_members = integer(0), stringsAsFactors = FALSE))
  }
  data.frame(map_id           = unlist(acc_map, use.names = FALSE),
             id_genome_effect = unlist(acc_id, use.names = FALSE),
             member_slot      = unlist(acc_slot, use.names = FALSE),
             label            = unlist(acc_lab, use.names = FALSE),
             n_members        = unlist(acc_n, use.names = FALSE),
             stringsAsFactors = FALSE)
}

#' Every label-vector for a family, one row per vector
#'
#' @keywords internal
#' @noRd
.gev_label_grid <- function(labs) {
  grid <- as.matrix(expand.grid(lapply(labs, seq_along), KEEP.OUT.ATTRS = FALSE))
  vecs <- matrix(NA_character_, nrow(grid), ncol(grid))
  for (m in seq_along(labs)) vecs[, m] <- labs[[m]][grid[, m]]
  vecs
}

#' Index of the most specific matching variant, or `NA` when nothing matches
#'
#' Returns a **position** in `preds`, not an `id_genome_effect`: the caller
#' caches the answer across families that share a resolution problem, and only
#' positions are portable between them. `ids` is used for the error message.
#'
#' @keywords internal
#' @noRd
.gev_select_variant <- function(preds, units, ids, key) {
  hit <- which(vapply(preds, function(P) {
    all(vapply(seq_along(units),
               function(m) .gev_member_match(P[[m]], units[[m]]), logical(1)))
  }, logical(1)))
  if (length(hit) == 0L) return(NA_integer_)
  minimal <- hit[vapply(hit, function(i) {
    !any(vapply(setdiff(hit, i),
                function(j) .ge_pred_leq(preds[[j]], preds[[i]]), logical(1)))
  }, logical(1))]
  if (length(minimal) != 1L) {
    stop("Ambiguous specificity in family '", key, "' among genome effect(s) ",
         paste(ids[minimal], collapse = ", "),
         ". Overlapping but incomparable scopes are refused at write time, so ",
         "this indicates the stored model was not written through ",
         "define_genome_effects() or define_additive_effects().", call. = FALSE)
  }
  minimal
}


# -- Label-vector preflight -------------------------------------------------

#' Warn or stop before building a map that is too large
#'
#' The subject is the resolved map -- `prod(|labels_m|)` per family -- which is
#' the object that actually grows with effect order. An earlier per-individual
#' framing made the cap either meaningless or unreachable.
#'
#' @keywords internal
#' @noRd
.gev_preflight <- function(model, alphabets, free) {
  warn_at <- getOption("tidybreed.label_vector_warn", 1e4)
  stop_at <- getOption("tidybreed.label_vector_max",  1e6)
  terms   <- model$terms
  mid     <- as.character(model$members$id_genome_effect)
  kind_by <- split(.gev_member_kind(model$members$contrast_name), mid)
  star_by <- split(free, mid)
  ids_by  <- split(terms$id_genome_effect, terms$family_key)

  worst <- 0
  worst_id <- NA_integer_
  for (key in names(ids_by)) {
    k1    <- as.character(ids_by[[key]][1])
    sizes <- ifelse(star_by[[k1]], 1L, lengths(alphabets[kind_by[[k1]]]))
    n <- prod(sizes)
    if (n > worst) { worst <- n; worst_id <- ids_by[[key]][1] }
  }
  # A family key can be fifty locus ids long, which is unreadable in a message.
  # Name the trait and one id the user can look up in genome_effect_terms.
  where <- if (is.na(worst_id)) "" else paste0(
    "trait '", terms$trait_name[terms$id_genome_effect == worst_id],
    "', genome effect ", worst_id, " (", terms$effect_order[
      terms$id_genome_effect == worst_id], " members)")
  if (worst > stop_at) {
    stop("Label-vector enumeration for ", where, " would be ",
         format(worst, scientific = FALSE), " rows, above the hard cap of ",
         format(stop_at, scientific = FALSE),
         ". Reduce the term's order or the number of distinct ",
         "(line_origin, parent_origin) labels, or raise ",
         "options(tidybreed.label_vector_max = ).", call. = FALSE)
  }
  if (worst > warn_at) {
    warning("Label-vector enumeration for ", where, " is ",
            format(worst, scientific = FALSE), " rows (warning threshold ",
            format(warn_at, scientific = FALSE),
            "). Set options(tidybreed.label_vector_warn = ) to change.",
            call. = FALSE)
  }
  invisible(worst)
}


# -- Evaluation -------------------------------------------------------------

#' Evaluate the stored effect model for a set of individuals and traits
#'
#' Member reduction, tuple grouping and summation are one SQL statement, so the
#' statement count is a function of the model, never of the number of
#' individuals.
#'
#' @param conn A DBI connection.
#' @param id_ind Character vector of individuals.
#' @param trait_names Character vector of traits.
#' @param model The model to evaluate. `NULL` reads the whole stored model for
#'   `trait_names`; `add_tbv()` passes `.gev_reserved_additive()` of one.
#' @return Data frame `id_ind`, `trait_name`, `component_name`, `tgv_value`.
#'   Individuals contributing nothing are absent rather than zero -- the caller
#'   decides whether that is an error.
#' @keywords internal
#' @noRd
.gev_evaluate <- function(conn, id_ind, trait_names, model = NULL) {

  if (is.null(model)) model <- .gev_read_model(conn, trait_names)
  empty <- data.frame(id_ind = character(0), trait_name = character(0),
                      component_name = character(0), tgv_value = numeric(0),
                      stringsAsFactors = FALSE)
  if (nrow(model$terms) == 0L) return(empty)

  stamp    <- as.character(round(as.numeric(Sys.time()) * 1000))
  ind_tmp  <- paste0("_gev_ind_",  stamp)
  mem_tmp  <- paste0("_gev_mem_",  stamp)
  map_tmp  <- paste0("_gev_map_",  stamp)
  term_tmp <- paste0("_gev_term_", stamp)
  on.exit({
    for (nm in c(ind_tmp, mem_tmp, map_tmp, term_tmp)) {
      try(duckdb::duckdb_unregister(conn, nm), silent = TRUE)
    }
  }, add = TRUE)

  duckdb::duckdb_register(conn, ind_tmp,
                          data.frame(id_ind = unique(id_ind),
                                     stringsAsFactors = FALSE))

  members <- model$members
  members$member_kind <- .gev_member_kind(members$contrast_name)
  free <- .gev_slot_freedom(model)
  members$use_star <- free

  alphabets <- list(
    additive = .gev_additive_alphabet(conn, ind_tmp)$label,
    genotype = .gev_genotype_alphabet(
      conn, ind_tmp,
      sort(unique(members$locus_id[members$member_kind == "genotype"])))$label)

  .gev_preflight(model, alphabets, free)
  map <- .gev_variant_map(model, alphabets, free)
  if (nrow(map) == 0L) return(empty)

  duckdb::duckdb_register(conn, mem_tmp, members)
  duckdb::duckdb_register(conn, map_tmp, map)
  duckdb::duckdb_register(conn, term_tmp,
                          model$terms[, c("id_genome_effect", "trait_name",
                                          "genome_value", "component_name")])

  res <- DBI::dbGetQuery(conn, .gev_sql(ind_tmp, mem_tmp, map_tmp, term_tmp))

  bad <- res[res$n_bad > 0L, , drop = FALSE]
  if (nrow(bad) > 0L) {
    ids <- sort(unique(bad$id_ind))
    stop("A 'dominance' contrast was selected for individual(s) ",
         paste(utils::head(ids, 5), collapse = ", "),
         if (length(ids) > 5) ", ..." else "",
         " at a locus where they carry a number of copies other than 2. ",
         "define_genome_effects() refuses a dominance member at a locus whose ",
         "chr_inheritance is not 1,1, so this means those individuals' ",
         "ind_haplotype rows disagree with chr_inheritance.", call. = FALSE)
  }
  res$n_bad <- NULL
  res
}

#' The single evaluation statement
#'
#' `add_red` and `geno_red` are the member reduction: one value per
#' `(id_ind, id_genome_effect, member_slot, label)`, where a member no variant
#' scopes (`use_star`) reduces straight to the sentinel label `"*"` — one row
#' summed over every unit instead of one row per label. Each member reduces to one
#' value per label **before** members are combined, so no row ever carries a
#' Cartesian pairing of raw haplotype rows. `hit` then multiplies one reduced
#' value per member -- the label-vector factorization -- and `HAVING` every slot
#' matched is what makes a partially matched label-vector contribute nothing
#' rather than a truncated product.
#'
#' @keywords internal
#' @noRd
.gev_sql <- function(ind_tmp, mem_tmp, map_tmp, term_tmp) {
  paste0(
    "WITH add_red AS ( ",
    "  SELECT h.id_ind, m.id_genome_effect, m.member_slot, ",
    "         CASE WHEN m.use_star THEN '*' ELSE ",
    "           COALESCE(h.line_origin, '') || '@' || h.parent_origin END AS label, ",
    "         SUM(h.allele - m.center_value) AS x ",
    "  FROM ind_haplotype h ",
    "  JOIN ", ind_tmp, " i ON i.id_ind = h.id_ind ",
    "  JOIN ", mem_tmp, " m ON m.locus_id = h.locus_id ",
    "                      AND m.member_kind = 'additive' ",
    "  GROUP BY 1, 2, 3, 4 ), ",
    "geno_loci AS ( ",
    "  SELECT DISTINCT locus_id FROM ", mem_tmp,
    "  WHERE member_kind = 'genotype' ), ",
    "state AS ( ",
    "  SELECT i.id_ind, g.locus_id, ",
    "         COUNT(h.locus_id) AS cc, ",
    "         COALESCE(SUM(h.allele), 0) AS dos, ",
    "         COALESCE(string_agg(COALESCE(h.line_origin, '') || '@' || ",
    "                             h.parent_origin, ',' ",
    "                  ORDER BY COALESCE(h.line_origin, ''), h.parent_origin), ",
    "                  '') AS label ",
    "  FROM ", ind_tmp, " i CROSS JOIN geno_loci g ",
    "  LEFT JOIN ind_haplotype h ",
    "    ON h.id_ind = i.id_ind AND h.locus_id = g.locus_id ",
    "  GROUP BY 1, 2 ), ",
    "geno_red AS ( ",
    "  SELECT s.id_ind, m.id_genome_effect, m.member_slot, ",
    "         CASE WHEN m.use_star THEN '*' ELSE s.label END AS label, ",
    "         CASE WHEN m.contrast_name = 'indicator' THEN ",
    "                CASE WHEN s.cc = m.copy_count_value ",
    "                      AND s.dos = m.dosage_value THEN 1.0 ELSE 0.0 END ",
    "              WHEN s.cc <> 2 THEN NULL ",
    "              WHEN s.dos = 0 THEN -2 * m.center_value * m.center_value ",
    "              WHEN s.dos = 1 THEN  2 * m.center_value * (1 - m.center_value) ",
    "              ELSE -2 * (1 - m.center_value) * (1 - m.center_value) ",
    "         END AS x ",
    "  FROM state s ",
    "  JOIN ", mem_tmp, " m ON m.locus_id = s.locus_id ",
    "                      AND m.member_kind = 'genotype' ), ",
    "red AS (SELECT * FROM add_red UNION ALL SELECT * FROM geno_red), ",
    "hit AS ( ",
    "  SELECT r.id_ind, mp.map_id, mp.id_genome_effect, ",
    "         product(r.x) AS prod, ",
    "         SUM(CASE WHEN r.x IS NULL THEN 1 ELSE 0 END) AS n_bad ",
    "  FROM ", map_tmp, " mp ",
    "  JOIN red r ON r.id_genome_effect = mp.id_genome_effect ",
    "            AND r.member_slot      = mp.member_slot ",
    "            AND r.label            = mp.label ",
    "  GROUP BY 1, 2, 3, mp.n_members ",
    "  HAVING COUNT(*) = mp.n_members ) ",
    "SELECT h.id_ind, t.trait_name, t.component_name, ",
    "       SUM(t.genome_value * h.prod) AS tgv_value, ",
    "       CAST(SUM(h.n_bad) AS INTEGER) AS n_bad ",
    "FROM hit h JOIN ", term_tmp, " t ",
    "  ON t.id_genome_effect = h.id_genome_effect ",
    "GROUP BY 1, 2, 3"
  )
}


# -- Shared entry-point plumbing --------------------------------------------

#' Individuals selected by a `tidybreed_table`, as a character vector
#'
#' @keywords internal
#' @noRd
.gev_subset_ids <- function(tbl, what) {
  pop <- tbl$pop
  if (length(tbl$pending_filter) == 0L) {
    ids <- dplyr::collect(get_table(pop, "ind_meta"))$id_ind
  } else {
    collected <- dplyr::collect(tbl)
    if (!"id_ind" %in% names(collected)) {
      stop("Filtered table '", tbl$table_name, "' must contain 'id_ind' to ",
           "subset individuals for ", what, " computation.", call. = FALSE)
    }
    sel <- unique(collected[["id_ind"]])
    ids <- get_table(pop, "ind_meta") |>
      dplyr::filter(.data$id_ind %in% !!sel) |>
      dplyr::collect() |>
      dplyr::pull("id_ind")
  }
  unique(ids)
}

#' Resolve and check a `trait_name` argument against `trait_meta`
#'
#' @keywords internal
#' @noRd
.gev_resolve_traits <- function(conn, trait_name) {
  if (is.null(trait_name)) {
    trait_name <- DBI::dbGetQuery(
      conn, "SELECT trait_name FROM trait_meta ORDER BY id_trait")$trait_name
    if (length(trait_name) == 0L) {
      stop("No traits found in trait_meta. Define traits with define_trait() ",
           "first.", call. = FALSE)
    }
    return(trait_name)
  }
  stopifnot(is.character(trait_name), length(trait_name) >= 1)
  lapply(trait_name, validate_sql_identifier, what = "trait name")
  found <- DBI::dbGetQuery(conn, paste0(
    "SELECT trait_name FROM trait_meta WHERE trait_name IN (",
    sql_in_list(trait_name, what = "trait name"), ")"))$trait_name
  missing_t <- setdiff(trait_name, found)
  if (length(missing_t) > 0L) {
    stop("Traits not found: ", paste(missing_t, collapse = ", "), call. = FALSE)
  }
  trait_name
}

#' Stop when a requested trait has no terms the caller can evaluate
#'
#' Takes the model the caller already read rather than issuing another query:
#' three more reads per trait would be pure waste, and the check is a row count.
#'
#' @keywords internal
#' @noRd
.gev_require_terms <- function(model, trait, tbv = FALSE) {
  if (sum(model$terms$trait_name == trait) > 0L) return(invisible(NULL))
  if (tbv) {
    stop("No order-one additive effects found for trait '", trait,
         "' under effect owner '", GE_ADDITIVE_OWNER, "'. ",
         "Call define_additive_effects() first. Terms written through ",
         "define_genome_effects() contribute to ind_tgv but never redefine ",
         "the breeding value.", call. = FALSE)
  }
  stop("No genome effects found for trait '", trait,
       "'. Call define_additive_effects() or define_genome_effects() first.",
       call. = FALSE)
}

#' Warn when `add_tbv()`'s coefficients have stopped being average effects
#'
#' `add_tbv()` reads only the reserved owner's order-one `additive` terms. That
#' is not a partial answer — a breeding value is the additive component by
#' definition — but it stops being *the model's* breeding value as soon as some
#' other term either contributes to it or shifts the coefficients it reads:
#'
#' - a non-reserved order-one `additive` term contributes to A and is skipped;
#' - an `indicator` surface is raw functional coding, so its additive projection
#'   is real and, at a locus that also carries a generated additive term, the
#'   stored `a` is no longer the average effect: `alpha = a + d(q - p)`;
#' - an interaction has an additive projection that depends on the other loci
#'   and on LD, so there is no local correction at all;
#' - an order-one `dominance` term is the **exception**. Cockerham coding is
#'   HWE-orthogonal, so it contributes nothing to A and leaves the co-located
#'   additive coefficient alone — provided it is centred at the same frequency.
#'   Warning on those would cry wolf on the common case, so it stays silent.
#'
#' @param full The whole stored model for one trait, every owner.
#' @keywords internal
#' @noRd
.gev_warn_tbv_stale <- function(full, trait) {
  terms <- full$terms[full$terms$trait_name == trait, , drop = FALSE]
  if (nrow(terms) == 0L) return(invisible(NULL))
  mem <- full$members[full$members$id_genome_effect %in% terms$id_genome_effect, ,
                      drop = FALSE]

  reserved <- terms$id_genome_effect[terms$effect_owner == GE_ADDITIVE_OWNER &
                                       terms$effect_order == 1L]
  res_mem <- mem[mem$id_genome_effect %in% reserved &
                   mem$contrast_name == "additive", , drop = FALSE]
  if (nrow(res_mem) == 0L) return(invisible(NULL))   # add_tbv() errors anyway
  others <- terms[!terms$id_genome_effect %in% res_mem$id_genome_effect, ,
                  drop = FALSE]
  if (nrow(others) == 0L) return(invisible(NULL))

  # A dominance member is orthogonal only if it is centred where the generated
  # additive term at that same locus is centred.
  centred_ok <- function(locus, centre) {
    have <- res_mem$center_value[res_mem$locus_id == locus]
    length(have) > 0L && any(abs(have - centre) < 1e-9)
  }

  kind_of <- vapply(others$id_genome_effect, function(id) {
    mm <- mem[mem$id_genome_effect == id, , drop = FALSE]
    if (nrow(mm) > 1L) return("interaction")
    if (mm$contrast_name == "dominance" &&
        centred_ok(mm$locus_id, mm$center_value)) return("orthogonal")
    mm$contrast_name
  }, character(1))
  stale <- kind_of != "orthogonal"
  if (!any(stale)) return(invisible(NULL))

  hit <- others[stale, , drop = FALSE]
  warning(
    "Trait '", trait, "' has ", nrow(hit), " term(s) add_tbv() does not read ",
    "(owner(s) ", paste0("'", sort(unique(hit$effect_owner)), "'", collapse = ", "),
    "; ", paste(sort(unique(kind_of[stale])), collapse = ", "), ") that ",
    "contribute to the additive component. tbv_value is the sum of the stored ",
    "'", GE_ADDITIVE_OWNER, "' additive coefficients, which are no longer ",
    "average effects: under functional coding alpha = a + d(q - p), and under ",
    "epistasis the average effect depends on other loci and on LD. ",
    "Use add_tgv() for the full genetic value; deriving average effects from a ",
    "general non-additive model is a separate calculation. A Cockerham ",
    "'dominance' term centred at the same frequency does not trigger this and ",
    "leaves tbv_value exact.",
    call. = FALSE)
  invisible(NULL)
}

#' Stop when an individual's every term missed
#'
#' A tuple whose predicate is unmatched contributes 0; an individual for whom
#' **every** term of a trait matches nothing is an error, because that is the
#' shape of a silently-zero trait -- a paternally qualified X-linked effect
#' scoring every male at 0, say.
#'
#' @keywords internal
#' @noRd
.gev_require_contribution <- function(ids, got, trait) {
  missing_ids <- setdiff(ids, got)
  if (length(missing_ids) == 0L) return(invisible(NULL))
  stop(
    "No allele copy matched any effect for trait '", trait,
    "' for individual(s): ",
    paste(utils::head(missing_ids, 5), collapse = ", "),
    if (length(missing_ids) > 5) ", ..." else "",
    ". Usual causes: (a) every locus with an effect for the trait sits on a ",
    "chromosome the individual does not inherit (chr_inheritance ",
    "from_parent_1 = 0 and from_parent_2 = 0, e.g. Y in females) -- see ",
    "define_chromosome(); or (b) every term for the trait is scoped to a ",
    "line or parent_origin the individual has no copies of (a male's X is ",
    "from_parent_1 = 0). Place the trait's loci on a chromosome these ",
    "individuals carry, widen the scope, or exclude them from the subset.",
    call. = FALSE)
}
