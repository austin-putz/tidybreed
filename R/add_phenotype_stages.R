#' The three stages of `add_phenotype()`
#'
#' @description
#' `add_phenotype()` runs in three stages with a strict boundary between them:
#'
#' * **Stage 1 — PLAN** (`.ap_plan()`): decide the final record list for every
#'   phenotype. Sex expression, the repeatable guard, fixed-effect
#'   contributions and `null_class_action = "skip"`, formula/composite TBV
#'   evaluation and `missing_component_action`, path classification,
#'   `pheno_number` assignment, the residual condition value of every record,
#'   and the random-effect level every record touches. **No random number is
#'   drawn and nothing is written** (the one prerequisite write is
#'   [add_tbv()], which materializes the TBVs the plan reads and is
#'   RNG-neutral).
#' * **Stage 2 — RESOLVE** (`.ap_resolve()`): every random draw of the call,
#'   in a fixed order, and the liability / type conversion — all in memory.
#'   Nothing is written. Two adapters over [find_covariance_blocks()] and
#'   [resolve_correlated_draws()] draw everything: the named-effect adapter
#'   (`.ap_resolve_named_effects()`, entity `(effect_name, level)`, stored
#'   in `phenotype_random_effects`, one draw per level reused forever), then
#'   the residual adapter (`.ap_resolve_residuals()`, entity
#'   `(id_ind, pheno_number)`, stored in `ind_phenotype`, one draw per
#'   record). In both, an entity draws its planned coordinates from the
#'   block's Gaussian conditional on the coordinates it has already
#'   realized — stored on disk from an earlier call, or (residuals only)
#'   fixed by `user_residual`.
#' * **Stage 3 — COMMIT** (`.ap_commit()`): one transaction that inserts the
#'   new `phenotype_random_effects` rows and the `ind_phenotype` records via
#'   `duckdb_register()` + `INSERT`. No random number is drawn, so the stream
#'   a call consumes is a function of the model and the plan only.
#'
#' The consequence is that a record that will not exist — an individual
#' skipped by `null_class_action`, excluded by a missing component, or
#' refused by the repeatable guard — never consumes RNG and never leaves
#' stochastic state behind. Planned records are ordered by `id_ind` within a
#' phenotype; that order is what `user_values` and `user_residual` match
#' positionally, and it does not depend on physical row order in the
#' database.
#'
#' See `plans/sample_correlated_effects.md` §5.5.
#'
#' @name add_phenotype_stages
#' @keywords internal
NULL


# ── Stage 1: PLAN ─────────────────────────────────────────────────────────────

#' Stage 1: plan every record of an `add_phenotype()` call
#'
#' @param tbl The `tidybreed_table` passed to `add_phenotype()`.
#' @param phenos Character vector of validated phenotype names.
#' @param user_values The `user_values` argument, or `NULL`.
#' @return `NULL` when no individual matched (a warning is issued), otherwise
#'   a list with `pop` (after `add_tbv()`), `phenos` (in evaluation order),
#'   `pheno_meta` (rows in that order) and `entries`: one list per phenotype,
#'   in the same order, each with
#'   \describe{
#'     \item{`phenotype_name`}{}
#'     \item{`path`}{`"model"`, `"derived_formula"` or `"user_values"`.}
#'     \item{`id_ind`, `pheno_number`}{The planned records, `id_ind`-ordered.
#'       `pheno_number` is the value Stage 3 writes.}
#'     \item{`tbv`}{Numeric per record (`"model"` path only).}
#'     \item{`fixed`}{Fixed-effect contribution per record (`"model"` path;
#'       zeros for `"derived_formula"`, which has no model terms).}
#'     \item{`random`}{Random-effect terms: a list of
#'       `list(effect_name, distribution, level)` where `level` is the
#'       character level of the grouping column per record (`NA` = none).}
#'     \item{`condition_table`, `condition_column`, `condition_value`}{The
#'       residual stratum lookup: `NULL` when the phenotype's residual block
#'       is unconditional or absent, otherwise the raw condition value per
#'       record as character (`NA` = `NULL` in the table).}
#'     \item{`user_values`}{Numeric per record (`"user_values"` path).}
#'     \item{`formula`}{The derived formula string (`"derived_formula"` path).}
#'   }
#'   and `residual_blocks`: the residual covariance blocks touching the
#'   call's phenotypes, from [find_covariance_blocks()].
#' @keywords internal
.ap_plan <- function(tbl, phenos, user_values = NULL) {
  pop  <- tbl$pop
  conn <- pop$db_conn

  # ── Subset ──────────────────────────────────────────────────────────────
  subset_ids <- resolve_subset_ids(tbl, "phenotyping")
  if (!is.null(subset_ids)) {
    ind_meta_subset <- get_table(pop, "ind_meta") |>
      dplyr::filter(.data$id_ind %in% !!subset_ids) |>
      dplyr::collect()
  } else {
    ind_meta_subset <- dplyr::collect(get_table(pop, "ind_meta"))
  }
  if (nrow(ind_meta_subset) == 0) {
    warning("No individuals matched the filter; no phenotypes generated.",
            call. = FALSE)
    return(NULL)
  }
  # Stable record order: sorted id_ind (byte order, so the same on every
  # platform and locale), never physical row order.
  ind_meta_subset <- ind_meta_subset[order(ind_meta_subset$id_ind,
                                           method = "radix"), , drop = FALSE]

  # ── Phenotype metadata ──────────────────────────────────────────────────
  phenos_in  <- .pvc_in_list(conn, phenos)
  pheno_meta <- DBI::dbGetQuery(conn, paste0(
    "SELECT * FROM phenotype_meta WHERE phenotype_name IN (", phenos_in, ")"))
  missing_p <- setdiff(phenos, pheno_meta$phenotype_name)
  if (length(missing_p) > 0) {
    stop("Phenotypes not found in phenotype_meta: ",
         paste(missing_p, collapse = ", "),
         ". Call define_phenotype() first.", call. = FALSE)
  }
  pheno_meta <- pheno_meta[match(phenos, pheno_meta$phenotype_name), , drop = FALSE]

  components_by_pheno <- lapply(phenos, function(t) DBI::dbGetQuery(conn, paste0(
    "SELECT * FROM phenotype_components WHERE phenotype_name = ",
    DBI::dbQuoteLiteral(conn, t))))
  names(components_by_pheno) <- phenos
  has_components <- vapply(components_by_pheno, nrow, integer(1)) > 0L

  .blank_to_na <- function(v) {
    v[is.na(v) | !nzchar(v)] <- NA_character_
    stats::setNames(as.character(v), phenos)
  }
  formula_tbv_str <- .blank_to_na(pheno_meta$formula_tbv)
  formula_str     <- .blank_to_na(pheno_meta$formula)
  has_formula_tbv <- !is.na(formula_tbv_str)
  has_formula     <- !is.na(formula_str)

  # Simple phenotypes need a term add_tbv() can read: order one, contrast
  # 'additive', under the reserved owner. A line-specific-only model
  # qualifies.
  for (t in phenos[!has_components & !has_formula_tbv & !has_formula]) {
    n_eff <- nrow(.gev_reserved_additive(
      .gev_read_model(conn, t, GE_ADDITIVE_OWNER))$terms)
    if (n_eff == 0L) {
      stop(
        "No additive effects found for phenotype '", t, "' in genome_effects. ",
        "For simple phenotypes call define_additive_effects() first. ",
        "For composite phenotypes supply 'components' or 'formula_tbv' in define_phenotype(). ",
        "For derived phenotypes (no genetic architecture) ",
        "use type = 'derived_formula'.",
        call. = FALSE
      )
    }
  }

  # Derived formulas read the phenotypes they reference, so those go first.
  if (any(has_formula)) {
    phenos              <- .topo_sort_phenotypes(pheno_meta)
    pheno_meta          <- pheno_meta[match(phenos, pheno_meta$phenotype_name), , drop = FALSE]
    has_components      <- has_components[phenos]
    has_formula_tbv     <- has_formula_tbv[phenos]
    has_formula         <- has_formula[phenos]
    formula_tbv_str     <- formula_tbv_str[phenos]
    formula_str         <- formula_str[phenos]
    components_by_pheno <- components_by_pheno[phenos]
  }

  # ── Sex expression and the repeatable guard ─────────────────────────────
  subset_by_pheno <- lapply(seq_along(phenos), function(i) {
    ex_sex <- pheno_meta$expressed_sex[i]
    if (is.null(ex_sex) || is.na(ex_sex) || ex_sex == "both") return(ind_meta_subset)
    ind_meta_subset[ind_meta_subset$sex == ex_sex, , drop = FALSE]
  })
  names(subset_by_pheno) <- phenos

  for (i in seq_along(phenos)) {
    t     <- phenos[i]
    ids_t <- subset_by_pheno[[t]]$id_ind
    if (length(ids_t) == 0 || isTRUE(pheno_meta$repeatable[i])) next
    already_done <- .ap_phenotyped_ids(conn, t, ids_t)
    n_rejected   <- length(already_done)
    if (n_rejected > 0) {
      keep <- ids_t[!ids_t %in% already_done]
      warning(
        "Phenotype '", t, "' is not repeatable: ",
        n_rejected, " individual(s) already phenotyped were skipped; ",
        length(keep), " individual(s) will receive a new phenotype record.",
        call. = FALSE
      )
      subset_by_pheno[[t]] <- subset_by_pheno[[t]][
        subset_by_pheno[[t]]$id_ind %in% keep, , drop = FALSE]
    }
  }

  # ── TBVs: the one write before Stage 3 (idempotent, RNG-neutral) ────────
  pop <- .ap_materialize_tbvs(pop, tbl, phenos, has_components, has_formula_tbv,
                              has_formula, components_by_pheno, formula_tbv_str,
                              subset_by_pheno)

  # ── Per-phenotype record planning ───────────────────────────────────────
  entries <- vector("list", length(phenos))
  names(entries) <- phenos
  for (i in seq_along(phenos)) {
    t <- phenos[i]
    path <- if (!is.null(user_values)) "user_values"
            else if (has_formula[t])   "derived_formula"
            else                       "model"
    entries[[t]] <- .ap_plan_phenotype(
      pop, t, pheno_meta[i, , drop = FALSE], subset_by_pheno[[t]], path,
      tbv_kind = if (has_formula_tbv[t]) "formula_tbv"
                 else if (has_components[t]) "components" else "simple",
      formula_tbv = formula_tbv_str[[t]], formula = formula_str[[t]],
      comp_rows = components_by_pheno[[t]], user_values = user_values,
      n_phenos = length(phenos))
  }

  # ── Residual stratum of every planned record ────────────────────────────
  # The blocks are loaded once, here; Stage 2 draws through the same list.
  residual_blocks <- find_covariance_blocks(conn, "residual", phenos)
  for (b in residual_blocks) {
    if (is.null(b$condition_column)) next
    for (t in intersect(b$phenotypes, phenos)) {
      e <- entries[[t]]
      if (e$path != "model" || length(e$id_ind) == 0L) next
      e$condition_table  <- b$condition_table
      e$condition_column <- b$condition_column
      e$condition_value  <- .ap_condition_values(conn, b$condition_table,
                                                 b$condition_column, e$id_ind)
      entries[[t]] <- e
    }
  }

  # ── Named-effect blocks: one loader call per effect, in sorted order ────
  # Only phenotypes that carry a random term for the effect are targets; a
  # block member without such a term is a latent coordinate.
  named_targets <- .ap_named_effect_targets(entries)
  named_blocks  <- lapply(names(named_targets), function(eff)
    find_covariance_blocks(conn, eff, named_targets[[eff]]))
  names(named_blocks) <- names(named_targets)

  list(pop = pop, phenos = phenos, pheno_meta = pheno_meta, entries = entries,
       residual_blocks = residual_blocks, named_targets = named_targets,
       named_blocks = named_blocks)
}


#' The model-path phenotypes carrying each random effect, by effect name
#'
#' @return A list named by `effect_name` (byte-sorted) of the phenotypes
#'   whose planned records touch that effect, in plan order.
#' @keywords internal
.ap_named_effect_targets <- function(entries) {
  pairs <- do.call(rbind, lapply(names(entries), function(t) {
    e <- entries[[t]]
    if (e$path != "model" || length(e$random) == 0L) return(NULL)
    data.frame(effect_name = vapply(e$random, `[[`, character(1), "effect_name"),
               phenotype_name = t, stringsAsFactors = FALSE)
  }))
  if (is.null(pairs)) return(list())
  effs <- sort(unique(pairs$effect_name), method = "radix")
  out <- lapply(effs, function(eff)
    pairs$phenotype_name[pairs$effect_name == eff])
  names(out) <- effs
  out
}


#' Materialize the TBVs Stage 1 reads (simple, composite, formula_tbv)
#' @keywords internal
.ap_materialize_tbvs <- function(pop, tbl, phenos, has_components,
                                 has_formula_tbv, has_formula,
                                 components_by_pheno, formula_tbv_str,
                                 subset_by_pheno) {
  conn <- pop$db_conn
  simple_phenos      <- phenos[!has_components & !has_formula_tbv & !has_formula]
  composite_phenos   <- phenos[has_components]
  formula_tbv_phenos <- phenos[has_formula_tbv]

  # Simple: phenotype_name == trait_name in trait_meta
  if (length(simple_phenos) > 0) {
    pop <- add_tbv(tbl, trait_name = simple_phenos)
  }

  parent_ids <- function(x) { x <- as.character(x); x[!is.na(x)] }

  # Composite — gather all contributor IDs + source traits, then add_tbv once
  if (length(composite_phenos) > 0) {
    all_source_traits   <- character(0)
    all_contributor_ids <- character(0)
    for (t in composite_phenos) {
      comp_rows <- components_by_pheno[[t]]
      all_source_traits <- unique(c(all_source_traits,
                                    as.character(comp_rows$source_trait_name)))
      subset_df <- subset_by_pheno[[t]]
      for (ct in unique(as.character(comp_rows$contributor_type))) {
        ids <- switch(ct,
          self = subset_df$id_ind,
          dam  = parent_ids(subset_df$id_parent_2),
          sire = parent_ids(subset_df$id_parent_1),
          group = {
            grp_rows <- comp_rows[comp_rows$contributor_type == "group", , drop = FALSE]
            unlist(lapply(seq_len(nrow(grp_rows)), function(gi) {
              g <- grp_rows[gi, , drop = FALSE]
              .group_members(conn, subset_df$id_ind, g$group_column, g$group_table,
                             what = paste0("Phenotype '", t, "', component '",
                                           g$source_trait_name, "' (group)"))
            }))
          })
        all_contributor_ids <- unique(c(all_contributor_ids, ids))
      }
    }
    if (length(all_contributor_ids) > 0 && length(all_source_traits) > 0) {
      contrib_tbl <- get_table(pop, "ind_meta") |>
        dplyr::filter(.data$id_ind %in% !!all_contributor_ids)
      pop <- add_tbv(contrib_tbl, trait_name = all_source_traits)
    }
  }

  # formula_tbv — gather all contributor IDs + source traits via AST walk
  if (length(formula_tbv_phenos) > 0) {
    all_source_traits   <- character(0)
    all_contributor_ids <- character(0)
    for (t in formula_tbv_phenos) {
      expr      <- parse(text = formula_tbv_str[[t]], keep.source = FALSE)[[1]]
      walk_res  <- .walk_formula_tbv_ast(expr)
      subset_df <- subset_by_pheno[[t]]
      all_source_traits <- unique(c(all_source_traits,
        vapply(walk_res$trait_refs, `[[`, character(1), "trait")))
      for (ref in walk_res$trait_refs) {
        ids <- switch(ref$type,
          self = as.character(subset_df$id_ind),
          dam  = parent_ids(subset_df$id_parent_2),
          sire = parent_ids(subset_df$id_parent_1),
          group_sum = , group_mean =
            .group_members(conn, subset_df$id_ind, ref$col, ref$table,
                           what = paste0("formula_tbv for phenotype '", t, "'")))
        all_contributor_ids <- unique(c(all_contributor_ids, ids))
      }
    }
    if (length(all_contributor_ids) > 0 && length(all_source_traits) > 0) {
      contrib_tbl <- get_table(pop, "ind_meta") |>
        dplyr::filter(.data$id_ind %in% !!all_contributor_ids)
      pop <- add_tbv(contrib_tbl, trait_name = all_source_traits)
    }
  }

  pop
}


#' Plan the records of one phenotype
#'
#' Applies the covariate skip and the TBV exclusions, reads the TBV, and
#' assigns `pheno_number`. Emits the same warnings and messages the
#' exclusions always have. No RNG, no writes.
#'
#' @keywords internal
.ap_plan_phenotype <- function(pop, t, m, subset_df, path, tbv_kind,
                               formula_tbv, formula, comp_rows, user_values,
                               n_phenos) {
  conn  <- pop$db_conn
  empty <- list(phenotype_name = t, path = path,
                id_ind = character(0), pheno_number = integer(0),
                tbv = numeric(0), fixed = numeric(0), random = list(),
                condition_table = NULL, condition_column = NULL,
                condition_value = NULL, user_values = NULL, formula = formula)

  # ── user_values: the model is not evaluated at all ──────────────────────
  if (path == "user_values") {
    vals <- if (n_phenos == 1 && !is.list(user_values)) user_values
            else user_values[[t]]
    if (is.null(vals)) {
      stop("user_values missing entry for phenotype '", t, "'.", call. = FALSE)
    }
    ids_t <- subset_df$id_ind
    if (!is.null(names(vals))) {
      unknown <- setdiff(names(vals), ids_t)
      if (length(unknown) > 0 || anyDuplicated(names(vals))) {
        stop("user_values for '", t, "': names must be individuals in the ",
             "subset for this phenotype (after sex expression and the ",
             "repeatable guard), each once. ",
             if (length(unknown) > 0) paste0(
               length(unknown), " unknown (e.g. ",
               paste(head(unknown, 5), collapse = ", "), "). "),
             if (anyDuplicated(names(vals))) "Duplicate names present.",
             call. = FALSE)
      }
      # Planned order is sorted id_ind, whatever order the names came in
      ids_t <- ids_t[ids_t %in% names(vals)]
      vals  <- unname(vals[ids_t])
    } else if (length(vals) != length(ids_t)) {
      stop("user_values for '", t, "' must have length equal to subset (",
           length(ids_t), ") or be a named vector.", call. = FALSE)
    }
    if (length(vals) == 0) return(empty)
    empty$id_ind       <- ids_t
    empty$user_values  <- as.numeric(vals)
    empty$pheno_number <- next_pheno_numbers(pop, t, ids_t)
    return(empty)
  }

  ids_t <- subset_df$id_ind
  if (length(ids_t) == 0) return(empty)

  # ── derived_formula: arithmetic over other records, no model terms ──────
  if (path == "derived_formula") {
    empty$id_ind       <- ids_t
    empty$pheno_number <- next_pheno_numbers(pop, t, ids_t)
    empty$tbv          <- rep(0, length(ids_t))
    empty$fixed        <- rep(0, length(ids_t))
    return(empty)
  }

  # ── Covariate terms; null_class_action = "skip" ─────────────────────────
  terms <- .ap_covariate_terms(pop, t, subset_df)
  skip_mask <- is.na(terms$fixed)
  if (any(skip_mask)) {
    warning("Phenotype '", t, "': ", sum(skip_mask),
            " individual(s) excluded due to null_class_action = 'skip'.",
            call. = FALSE)
    keep      <- !skip_mask
    ids_t     <- ids_t[keep]
    subset_df <- subset_df[keep, , drop = FALSE]
    terms     <- .ap_subset_terms(terms, keep)
    if (length(ids_t) == 0) {
      message("Phenotype '", t, "': all individuals skipped; no records written.")
      return(empty)
    }
  }

  # ── TBV ─────────────────────────────────────────────────────────────────
  if (tbv_kind == "formula_tbv") {
    mca     <- .ap_missing_action(m)
    raw_tbv <- .eval_formula_tbv(pop, formula_tbv, subset_df, t)
    tbv     <- unname(raw_tbv[ids_t])
    excl    <- is.na(tbv)
    if (any(excl)) {
      n_excl   <- sum(excl)
      excl_ids <- ids_t[excl]
      msg <- paste0(
        n_excl, " individual(s) had one or more missing components for phenotype '",
        t, "' (formula_tbv: ", formula_tbv, ") and were excluded. ",
        "(IDs: ", paste(head(excl_ids, 5), collapse = ", "),
        if (n_excl > 5) paste0(" ... +", n_excl - 5L, " more") else "", ")"
      )
      if (mca == "error") stop(msg, call. = FALSE) else warning(msg, call. = FALSE)
      ids_t <- ids_t[!excl]; tbv <- tbv[!excl]; terms <- .ap_subset_terms(terms, !excl)
      if (length(ids_t) == 0) {
        message("Phenotype '", t, "': all individuals excluded (no formula_tbv result).")
        return(empty)
      }
    }
  } else if (tbv_kind == "components") {
    mca <- .ap_missing_action(m)
    tbv  <- unname(.assemble_composite_tbv(pop, t, comp_rows, subset_df,
                                           missing_component_action = mca)[ids_t])
    excl <- is.na(tbv)
    if (any(excl)) {
      ids_t <- ids_t[!excl]; tbv <- tbv[!excl]; terms <- .ap_subset_terms(terms, !excl)
      if (length(ids_t) == 0) {
        message("Phenotype '", t, "': all individuals excluded (no composite TBV).")
        return(empty)
      }
    }
  } else {
    tbv_rows <- .ap_read_by_id(conn, "ind_tbv", ids_t, "tbv_value",
                               where = paste0("t.trait_name = ",
                                              DBI::dbQuoteLiteral(conn, t)))
    tbv <- unname(stats::setNames(tbv_rows$tbv_value, tbv_rows$id_ind)[ids_t])
  }

  empty$id_ind       <- ids_t
  empty$pheno_number <- next_pheno_numbers(pop, t, ids_t)
  empty$tbv          <- as.numeric(tbv)
  empty$fixed        <- terms$fixed
  empty$random       <- terms$random
  empty
}


.ap_missing_action <- function(m) {
  v <- m$missing_component_action
  if (is.na(v) || !nzchar(v)) "skip" else v
}

.ap_subset_terms <- function(terms, keep) {
  terms$fixed  <- terms$fixed[keep]
  terms$random <- lapply(terms$random, function(r) { r$level <- r$level[keep]; r })
  terms
}


#' Rows of `table` for a set of individuals, by registered-view join
#'
#' The ids never enter the SQL text. Returns whatever rows match — zero,
#' one or several per id — so the caller applies its own contract.
#'
#' @param columns Character vector of columns to return besides `id_ind`.
#' @param where Optional extra predicate on `t` (already quoted SQL).
#' @keywords internal
.ap_read_by_id <- function(conn, table, ids, columns, where = NULL) {
  tmp <- "__ap_ids"
  duckdb::duckdb_register(conn, tmp, data.frame(id_ind = unique(ids),
                                                stringsAsFactors = FALSE))
  on.exit(try(duckdb::duckdb_unregister(conn, tmp), silent = TRUE), add = TRUE)
  columns <- setdiff(columns, "id_ind")
  cols <- if (length(columns) == 0) "" else
    paste0(", ", paste0("t.\"", columns, "\"", collapse = ", "))
  DBI::dbGetQuery(conn, paste0(
    "SELECT t.id_ind", cols, " FROM ", table, " AS t ",
    "JOIN ", tmp, " AS f USING (id_ind)",
    if (!is.null(where)) paste0(" WHERE ", where) else ""))
}


#' Distinct ids among `ids` that already have a record for `phenotype_name`
#' @keywords internal
.ap_phenotyped_ids <- function(conn, phenotype_name, ids) {
  unique(.ap_read_by_id(conn, "ind_phenotype", ids, "phenotype_name",
    where = paste0("t.phenotype_name = ",
                   DBI::dbQuoteLiteral(conn, phenotype_name)))$id_ind)
}


#' The residual condition value of every planned record
#'
#' `.read_one_per_id()` on the condition table (exactly one row per planned
#' id), returned as character (`NA` for `NULL`) to match how
#' `phenotype_var_comp.condition_level` is stored.
#' @keywords internal
.ap_condition_values <- function(conn, condition_table, condition_column, ids) {
  v <- .read_one_per_id(conn, condition_table, condition_column, ids,
                        what = "Residual condition lookup")
  out <- as.character(v)
  out[is.na(v)] <- NA_character_
  out
}


# ── Stage 2: RESOLVE ──────────────────────────────────────────────────────────

#' Stage 2: every random draw and the in-memory record assembly
#'
#' Order of RNG consumption, fixed regardless of database row order: the
#' named-effect adapter first (`.ap_resolve_named_effects()`: effects in
#' byte-sorted `effect_name` order, blocks in `find_covariance_blocks()`
#' order, one resolver call per sample-set group in sorted group order),
#' then the residual adapter (`.ap_resolve_residuals()`: one block at a
#' time, in loader order; within a block one resolver call per
#' `(stratum, sample set)` group in sorted group order). Record assembly
#' (derived formulas, liability, type conversion) follows in plan order and
#' draws nothing. Nothing is written.
#'
#' @param plan The Stage-1 plan.
#' @param user_residual The `user_residual` argument, or `NULL`.
#' @return A list with `records` (one tibble per phenotype, in plan order,
#'   possibly empty) and `random_effects` (a data frame of new
#'   `phenotype_random_effects` rows, possibly empty).
#' @keywords internal
.ap_resolve <- function(plan, user_residual = NULL) {
  pop     <- plan$pop
  phenos  <- plan$phenos
  entries <- plan$entries

  # Only model-path phenotypes with planned records draw anything; a
  # user_values or derived_formula call is RNG-neutral.
  is_model <- vapply(entries, function(e) e$path == "model" &&
                                             length(e$id_ind) > 0L, logical(1))

  # user_residual is checked against the plan whether or not anything draws
  fixed <- .ap_fixed_residuals(entries, user_residual)

  # ── Draws: named effects, then residuals ─────────────────────────────────
  pending_re <- .ap_empty_random_effects()
  random_contrib <- list()
  residuals <- list()
  if (any(is_model)) {
    ne <- .ap_resolve_named_effects(plan)
    pending_re     <- ne$pending
    random_contrib <- ne$contribution
    residuals      <- .ap_resolve_residuals(plan, fixed)
  }

  # ── Record assembly, in plan order (no RNG) ──────────────────────────────
  records <- vector("list", length(phenos))
  names(records) <- phenos
  for (i in seq_along(phenos)) {
    t <- phenos[i]
    e <- entries[[t]]
    m <- plan$pheno_meta[i, , drop = FALSE]
    n <- length(e$id_ind)

    # (An element left NULL stays in the list; `records[[t]] <- NULL` would
    # drop it and shift the plan-order positions.)
    if (n == 0) next
    if (e$path == "user_values") {
      records[[t]] <- tibble::tibble(
        id_ind = e$id_ind, phenotype_name = t,
        pheno_value = e$user_values, pheno_number = e$pheno_number)
      next
    }

    if (e$path == "derived_formula") {
      derived_vals <- .eval_derived_formula(
        pop, e$formula, e$id_ind, t, pending = dplyr::bind_rows(records))
      records[[t]] <- tibble::tibble(
        id_ind = e$id_ind, phenotype_name = t,
        pheno_value = as.numeric(derived_vals), pheno_number = e$pheno_number)
      next
    }

    r <- residuals[[t]]
    pheno_mean <- if (is.na(m$mean)) 0 else m$mean
    liability  <- pheno_mean + e$fixed + random_contrib[[t]] + e$tbv + r$value
    records[[t]] <- .ap_liability_records(pop, t, m, e, liability, r)
  }

  list(records = records, random_effects = pending_re)
}


# ── Stage 2: the named-effect adapter ─────────────────────────────────────────

#' The named-effect adapter: every random-effect draw of the call
#'
#' Implements `plans/sample_correlated_effects.md` §5.6 and §5.8 for every
#' `effect_name` other than `'residual'`. The entity is the *level* — a pen,
#' a herd, an `id_ind` for a permanent-environment effect — and a level's
#' draw is realized once and reused by every record that ever touches it,
#' in this call or any later one. Effects are processed in byte-sorted
#' order, each through its blocks (from the plan, in
#' [find_covariance_blocks()] order): a level draws its planned coordinates
#' conditional on the coordinates already stored in
#' `phenotype_random_effects` for the block's other phenotypes. Every
#' model-path phenotype with a random term for the effect must be in a
#' block (else "No variance stored"); the §5.6 checks are re-run here as
#' the backstop; a 1 x 1 block whose effect is `gamma` or `uniform` keeps
#' its marginal sampler.
#'
#' @param plan The Stage-1 plan.
#' @return A list with `contribution` (named by model-path phenotype with
#'   planned records: the summed random-effect value per record, `0` for a
#'   record whose level is `NULL`) and `pending` (the new
#'   `phenotype_random_effects` rows for Stage 3).
#' @keywords internal
.ap_resolve_named_effects <- function(plan) {
  entries <- plan$entries
  targets_all <- names(entries)[vapply(entries, function(e)
    e$path == "model" && length(e$id_ind) > 0L, logical(1))]
  contribution <- lapply(targets_all, function(t)
    rep(0, length(entries[[t]]$id_ind)))
  names(contribution) <- targets_all
  pending <- .ap_empty_random_effects()

  for (eff in names(plan$named_targets)) {
    targets <- plan$named_targets[[eff]]
    blocks  <- Filter(function(b) any(b$phenotypes %in% targets),
                      plan$named_blocks[[eff]])
    covered <- unlist(lapply(blocks, `[[`, "phenotypes"))
    missing <- setdiff(targets, covered)
    if (length(missing) > 0L) {
      stop("No variance stored for random effect '", eff, "' / phenotype '",
           missing[[1L]], "'. Specify via define_effect_random(variance = ...) ",
           "or define_effect_cov_matrix().", call. = FALSE)
    }
    for (b in blocks) {
      res <- .ap_named_effect_block(plan, b, targets)
      pending <- rbind(pending, res$pending)
      for (t in names(res$contribution)) {
        contribution[[t]] <- contribution[[t]] + res$contribution[[t]]
      }
    }
  }
  list(contribution = contribution, pending = pending)
}


.ap_empty_random_effects <- function() {
  data.frame(phenotype_name = character(0), effect_name = character(0),
             level = character(0), draw_value = numeric(0),
             stringsAsFactors = FALSE)
}


#' Resolve the draws of one named-effect covariance block
#'
#' @param plan The Stage-1 plan.
#' @param b One block from [find_covariance_blocks()] for a named effect.
#' @param targets The model-path phenotypes of the call with planned
#'   records and a random term for `b$effect_name`.
#' @return A list with `contribution` (named by the block's in-call
#'   phenotypes, one value per planned record) and `pending` (new
#'   `phenotype_random_effects` rows: every level drawn here, by coordinate
#'   then level).
#' @keywords internal
.ap_named_effect_block <- function(plan, b, targets) {
  conn    <- plan$pop$db_conn
  entries <- plan$entries
  eff     <- b$effect_name
  coords  <- b$phenotypes
  in_call <- intersect(coords, targets)      # sorted (block order)

  # A named-effect block has exactly one, unconditional stratum: the writers
  # never store a condition on it, so anything else was edited by hand.
  if (is.null(b$unconditional) || length(b$conditional) > 0L) {
    stop("The '", eff, "' covariance block ", .pvc_set(coords),
         " has conditional strata; condition_column is residual-only. ",
         "Redeclare the block with define_effect_cov_matrix().", call. = FALSE)
  }
  R <- b$unconditional

  # §5.6 backstop: normal, random, one (source_column, source_table)
  validate_named_effect_block(conn, eff, coords, caller = "add_phenotype()")

  # ── The block's random term of each in-call phenotype ───────────────────
  terms <- lapply(in_call, function(t) {
    r <- entries[[t]]$random
    r[[which(vapply(r, `[[`, character(1), "effect_name") == eff)]]
  })
  names(terms) <- in_call

  # ── Entities: the planned levels, sorted ─────────────────────────────────
  planned <- do.call(rbind, lapply(in_call, function(t) {
    lv <- terms[[t]]$level
    data.frame(level = unique(lv[!is.na(lv)]), phenotype_name = t,
               stringsAsFactors = FALSE)
  }))
  levels <- sort(unique(planned$level), method = "radix")
  n_ent  <- length(levels)
  contribution <- lapply(in_call, function(t)
    rep(0, length(entries[[t]]$id_ind)))
  names(contribution) <- in_call
  pending <- .ap_empty_random_effects()
  if (n_ent == 0L) return(list(contribution = contribution, pending = pending))

  # ── Stored coordinates of every block member at the planned levels ──────
  tmp <- "__ap_levels"
  duckdb::duckdb_register(conn, tmp, data.frame(level = levels,
                                                stringsAsFactors = FALSE))
  on.exit(try(duckdb::duckdb_unregister(conn, tmp), silent = TRUE), add = TRUE)
  stored <- DBI::dbGetQuery(conn, paste0(
    "SELECT r.phenotype_name, r.level, r.draw_value ",
    "FROM phenotype_random_effects AS r JOIN ", tmp, " AS l USING (level) ",
    "WHERE r.effect_name = ", DBI::dbQuoteLiteral(conn, eff), " ",
    "AND r.phenotype_name IN (", .pvc_in_list(conn, coords), ")"))

  value <- matrix(NA_real_, n_ent, length(coords), dimnames = list(NULL, coords))
  value[cbind(match(stored$level, levels),
              match(stored$phenotype_name, coords))] <- stored$draw_value

  # ── Sample set per entity: planned and not yet stored ───────────────────
  sample <- matrix(FALSE, n_ent, length(in_call), dimnames = list(NULL, in_call))
  for (t in in_call) {
    rows <- match(planned$level[planned$phenotype_name == t], levels)
    sample[rows, t] <- is.na(value[rows, t])
  }

  # A 1 x 1 block is the one place a non-normal distribution is legal
  # (§5.6); its marginal sampler is kept. Everything else is the resolver.
  dist <- if (length(coords) == 1L) terms[[1L]]$distribution else "normal"
  if (!is.na(dist) && dist %in% c("gamma", "uniform")) {
    new <- which(sample[, 1L])
    v   <- R[1L, 1L]
    value[new, 1L] <- switch(
      dist,
      gamma   = stats::rgamma(length(new), shape = 1, rate = 1 / sqrt(v)),
      uniform = stats::runif(length(new), min = -sqrt(3 * v),
                             max = sqrt(3 * v)))
  } else {
    # ── One resolver call per sample set, in sorted order ─────────────────
    sample_key <- rep("", n_ent)
    for (t in in_call) {
      sample_key <- ifelse(sample[, t], paste(sample_key, t, sep = ","),
                           sample_key)
    }
    for (g in sort(unique(sample_key[nzchar(sample_key)]), method = "radix")) {
      rows     <- which(sample_key == g)
      S        <- in_call[sample[rows[[1L]], ]]
      obs_cols <- setdiff(coords, S)
      obs <- if (length(obs_cols) > 0L)
        value[rows, obs_cols, drop = FALSE] else NULL
      value[rows, S] <- resolve_correlated_draws(
        R, S, entity_keys = data.frame(level = levels[rows],
                                       stringsAsFactors = FALSE),
        observed = obs)
    }
  }

  # ── New rows for Stage 3, and the per-record contribution ───────────────
  for (t in in_call) {
    new <- which(sample[, t])
    if (length(new) > 0L) {
      pending <- rbind(pending, data.frame(
        phenotype_name = t, effect_name = eff, level = levels[new],
        draw_value = unname(value[new, t]), stringsAsFactors = FALSE))
    }
    lv <- terms[[t]]$level
    per_record <- unname(value[match(lv, levels), t])
    per_record[is.na(lv)] <- 0
    contribution[[t]] <- per_record
  }
  list(contribution = contribution, pending = pending)
}


# ── Stage 2: the residual adapter ─────────────────────────────────────────────

#' The residual adapter: every model-path residual of the call
#'
#' Implements `plans/sample_correlated_effects.md` §5.3–§5.5 for
#' `effect_name = 'residual'`. Every model-path phenotype with at least one
#' planned record must belong to a residual covariance block (else "No
#' residual variance found") unless its residuals are all fixed by
#' `user_residual`, in which case nothing is drawn or conditioned for it.
#' Blocks are the plan's, processed in
#' [find_covariance_blocks()] order; within a block the entity is
#' `(id_ind, pheno_number)`, the coordinates are the block's phenotypes, and
#' each entity's residual is drawn from the stratum its condition value
#' selects, conditional on what it has already realized — stored residuals
#' of *any* block member at the same `pheno_number` (subject to D2) plus the
#' `user_residual` values fixed in this call.
#'
#' @param plan The Stage-1 plan.
#' @param fixed The validated `user_residual` list from
#'   `.ap_fixed_residuals()`.
#' @return A list named by model-path phenotype (those with planned
#'   records). Each element has `value` (numeric per planned record),
#'   `level` (the `residual_condition_level` per record: the selected
#'   stratum, `NA` for the unconditional `R`) and `var_unconditional` (the
#'   phenotype's unconditional residual variance, `NA` if no unconditional
#'   stratum is stored).
#' @keywords internal
.ap_resolve_residuals <- function(plan, fixed = list()) {
  entries <- plan$entries
  targets <- names(entries)[vapply(entries, function(e) e$path == "model" &&
                                                length(e$id_ind) > 0L,
                                   logical(1))]
  out <- list()
  if (length(targets) == 0L) return(out)

  blocks  <- Filter(function(b) any(b$phenotypes %in% targets),
                    plan$residual_blocks)
  covered <- unlist(lapply(blocks, `[[`, "phenotypes"))
  missing <- setdiff(targets, covered)
  needs_draw <- setdiff(missing, names(fixed))
  if (length(needs_draw) > 0L) {
    stop("No residual variance found for phenotype '", needs_draw[[1L]], "'. ",
         "Specify via define_phenotype(residual_var = ...) or ",
         "define_residual_cov().", call. = FALSE)
  }
  for (t in missing) {                      # fixed by the caller, no block
    out[[t]] <- list(value = fixed[[t]],
                     level = rep(NA_character_, length(fixed[[t]])),
                     var_unconditional = NA_real_)
  }

  for (b in blocks) {
    res <- .ap_residual_block(plan, b, targets, fixed)
    for (t in names(res)) {
      res[[t]]$var_unconditional <-
        if (is.null(b$unconditional)) NA_real_ else b$unconditional[t, t]
    }
    out[names(res)] <- res
  }
  out
}


#' Validate `user_residual` against the plan
#'
#' A single model-path phenotype takes a numeric vector; otherwise a named
#' list keyed by `phenotype_name` that may name any subset of the model-path
#' phenotypes. Each vector is positional over that phenotype's planned
#' records.
#'
#' @param entries The plan's entries.
#' @param user_residual The `user_residual` argument, or `NULL`.
#' @return A list named by phenotype of finite numeric vectors (possibly
#'   empty).
#' @keywords internal
.ap_fixed_residuals <- function(entries, user_residual) {
  if (is.null(user_residual)) return(list())
  model <- names(entries)[vapply(entries, function(e) e$path == "model",
                                 logical(1))]
  if (is.list(user_residual)) {
    nm <- names(user_residual)
    if (length(user_residual) > 0L &&
        (is.null(nm) || anyNA(nm) || !all(nzchar(nm)) || anyDuplicated(nm))) {
      stop("user_residual must be a named list keyed by phenotype_name, ",
           "each name once.", call. = FALSE)
    }
    bad <- setdiff(nm, model)
    if (length(bad) > 0L) {
      stop("user_residual names phenotype(s) that are not generated from the ",
           "model in this call: ", .pvc_set(bad), ". Model-path phenotypes: ",
           .pvc_set(model), ".", call. = FALSE)
    }
    vals <- user_residual
  } else {
    if (length(model) == 0L) {
      stop("user_residual was supplied but no phenotype in the call is ",
           "generated from the model (derived_formula phenotypes have no ",
           "residual).", call. = FALSE)
    }
    if (length(model) > 1L) {
      stop("user_residual must be a named list keyed by phenotype_name when ",
           "more than one phenotype is generated from the model in the call ",
           "(", .pvc_set(model), ").", call. = FALSE)
    }
    vals <- stats::setNames(list(user_residual), model)
  }
  for (t in names(vals)) {
    v <- vals[[t]]
    n <- length(entries[[t]]$id_ind)
    if (!is.null(names(v))) {
      stop("user_residual for phenotype '", t, "' must be an unnamed numeric ",
           "vector matched by position to the planned records (id_ind ",
           "order); per-id_ind names are not supported.", call. = FALSE)
    }
    if (!is.numeric(v) || length(v) != n) {
      stop("user_residual length for phenotype '", t, "' must equal ", n,
           " (the planned records, in id_ind order); got a ",
           if (is.numeric(v)) paste0("length-", length(v), " vector")
           else class(v)[[1L]], ".", call. = FALSE)
    }
    if (!all(is.finite(v))) {
      stop("user_residual for phenotype '", t, "' must be finite.",
           call. = FALSE)
    }
    vals[[t]] <- as.numeric(v)
  }
  vals
}


#' Resolve the residuals of one covariance block
#'
#' @param plan The Stage-1 plan.
#' @param b One block from [find_covariance_blocks()].
#' @param targets The model-path phenotypes of the call with planned records.
#' @param fixed The validated `user_residual` list.
#' @return A list named by the block's in-call phenotypes, each
#'   `list(value, level)` per planned record.
#' @keywords internal
.ap_residual_block <- function(plan, b, targets, fixed) {
  conn    <- plan$pop$db_conn
  entries <- plan$entries
  in_call <- intersect(b$phenotypes, targets)      # sorted (block order)
  coords  <- b$phenotypes

  # ── Planned coordinates: one row per (entity, phenotype) ─────────────────
  planned <- do.call(rbind, lapply(in_call, function(t) {
    e <- entries[[t]]
    data.frame(id_ind = e$id_ind, pheno_number = e$pheno_number,
               phenotype_name = t,
               condition_value = if (is.null(e$condition_value))
                 rep(NA_character_, length(e$id_ind)) else e$condition_value,
               stringsAsFactors = FALSE)
  }))

  # ── Entities, sorted ─────────────────────────────────────────────────────
  ents <- unique(planned[, c("id_ind", "pheno_number")])
  ents <- ents[order(ents$id_ind, ents$pheno_number, method = "radix"), ,
               drop = FALSE]
  rownames(ents) <- NULL
  n_ent <- nrow(ents)
  ent_key <- function(id, pn) paste(id, pn, sep = "\r")
  ents$key    <- ent_key(ents$id_ind, ents$pheno_number)
  planned$row <- match(ent_key(planned$id_ind, planned$pheno_number), ents$key)

  # ── Stratum per entity (D2 selection and fallback) ───────────────────────
  cond <- planned$condition_value[match(ents$key,
                                        ent_key(planned$id_ind,
                                                planned$pheno_number))]
  stratum <- ifelse(!is.na(cond) & cond %in% names(b$conditional),
                    cond, NA_character_)
  no_match <- is.na(stratum)
  if (any(no_match) && is.null(b$unconditional)) {
    lv <- cond[no_match]
    lv <- ifelse(is.na(lv), "NULL", paste0("'", lv, "'"))
    stop("add_phenotype(): ", sum(no_match), " planned record(s) of block ",
         .pvc_set(coords), " resolve to no residual stratum (",
         b$condition_table, ".", b$condition_column, " = ",
         paste(unique(lv), collapse = ", "), "; strata stored: ",
         .pvc_set(names(b$conditional)),
         ") and the block has no unconditional stratum to fall back on ",
         "(e.g. ", paste(utils::head(ents$id_ind[no_match], 5L), collapse = ", "),
         ").", call. = FALSE)
  }
  unmatched <- no_match & !is.na(cond)
  if (any(unmatched)) {
    warning("add_phenotype(): ", sum(unmatched), " planned record(s) of block ",
            .pvc_set(coords), " have a ", b$condition_column,
            " value matching no residual stratum (",
            paste0("'", unique(cond[unmatched]), "'", collapse = ", "),
            "); drawn from the unconditional R (e.g. ",
            paste(utils::head(ents$id_ind[unmatched], 5L), collapse = ", "),
            ").", call. = FALSE)
  }

  # ── Stored coordinates at the same pheno_number (any block member) ──────
  tmp <- "__ap_entities"
  duckdb::duckdb_register(conn, tmp, ents[, c("id_ind", "pheno_number")])
  on.exit(try(duckdb::duckdb_unregister(conn, tmp), silent = TRUE), add = TRUE)
  stored <- DBI::dbGetQuery(conn, paste0(
    "SELECT p.id_ind, p.pheno_number, p.phenotype_name, p.residual_value, ",
    "p.residual_condition_level ",
    "FROM ind_phenotype AS p JOIN ", tmp, " AS e USING (id_ind, pheno_number) ",
    "WHERE p.phenotype_name IN (", .pvc_in_list(conn, coords), ") ",
    "AND p.residual_value IS NOT NULL"))
  if (anyDuplicated(stored[, c("id_ind", "pheno_number", "phenotype_name")])) {
    stop("ind_phenotype has more than one record for the same (id_ind, ",
         "phenotype_name, pheno_number) with a residual_value; the residual ",
         "conditioning of block ", .pvc_set(coords), " is undefined.",
         call. = FALSE)
  }
  stored$row <- match(ent_key(stored$id_ind, stored$pheno_number), ents$key)

  # ── D6 (always), then D2 on the stored coordinates ──────────────────────
  .check_condition_change_agreement(conn, coords, caller = "add_phenotype()")
  if (nrow(stored) > 0L) {
    now  <- stratum[stored$row]
    was  <- stored$residual_condition_level
    same <- (is.na(now) & is.na(was)) | (!is.na(now) & !is.na(was) & now == was)
    if (!all(same)) {
      m <- plan$pheno_meta[match(in_call[[1L]], plan$pheno_meta$phenotype_name),
                           "condition_change_action"]
      action <- if (is.na(m)) "error" else m
      bad <- stored[!same, , drop = FALSE]
      show <- function(x) ifelse(is.na(x), "NULL", paste0("'", x, "'"))
      examples <- paste0(bad$id_ind, " (", bad$phenotype_name, ": stored under ",
                         show(bad$residual_condition_level), ", now ",
                         show(now[!same]), ")")
      examples <- paste(utils::head(unique(examples), 5L), collapse = "; ")
      n_ent_bad <- length(unique(bad$row))
      if (action == "independent") {
        warning("add_phenotype(): ", n_ent_bad, " record(s) in block ",
                .pvc_set(coords), " have a stored residual drawn under a ",
                "different stratum; those coordinates were dropped from the ",
                "conditioning set (condition_change_action = 'independent'): ",
                "dropped ", .pvc_set(sort(unique(bad$phenotype_name))),
                "; e.g. ", examples, ".", call. = FALSE)
        stored <- stored[same, , drop = FALSE]
      } else {
        stop("add_phenotype(): ", n_ent_bad, " record(s) in block ",
             .pvc_set(coords), " have a stored residual drawn under a ",
             "different residual stratum than the current record resolves ",
             "to (e.g. ", examples, "). Nothing defines the covariance ",
             "across strata. To draw independently of the mismatched values, ",
             "set condition_change_action = 'independent' on every phenotype ",
             "in the block with define_phenotype(..., overwrite = TRUE).",
             call. = FALSE)
      }
    }
  }

  # ── Observed matrix (stored + fixed) and sample set per entity ──────────
  observed <- matrix(NA_real_, n_ent, length(coords),
                     dimnames = list(NULL, coords))
  if (nrow(stored) > 0L) {
    observed[cbind(stored$row, match(stored$phenotype_name, coords))] <-
      stored$residual_value
  }
  value  <- observed * NA_real_
  sample <- matrix(FALSE, n_ent, length(in_call), dimnames = list(NULL, in_call))
  for (t in in_call) {
    rows <- planned$row[planned$phenotype_name == t]
    if (is.null(fixed[[t]])) {
      sample[rows, t] <- TRUE
    } else {
      observed[rows, t] <- fixed[[t]]
      value[rows, t]    <- fixed[[t]]
    }
  }

  # ── One resolver call per (stratum, sample set), in sorted order ────────
  sample_key <- rep("", n_ent)
  for (t in in_call) {
    sample_key <- ifelse(sample[, t], paste(sample_key, t, sep = ","), sample_key)
  }
  group_key <- paste(ifelse(is.na(stratum), 0L, 1L),
                     ifelse(is.na(stratum), "", stratum), sample_key,
                     sep = "\r")
  for (g in sort(unique(group_key), method = "radix")) {
    rows <- which(group_key == g)
    S    <- in_call[sample[rows[[1L]], ]]
    if (length(S) == 0L) next                       # every coordinate fixed
    lv   <- stratum[[rows[[1L]]]]
    R    <- if (is.na(lv)) b$unconditional else b$conditional[[lv]]
    obs_cols <- setdiff(coords, S)
    obs <- if (length(obs_cols) > 0L)
      observed[rows, obs_cols, drop = FALSE] else NULL
    value[rows, S] <- resolve_correlated_draws(
      R, S, entity_keys = ents[rows, c("id_ind", "pheno_number")],
      observed = obs)
  }

  # ── Back to planned-record order per phenotype ──────────────────────────
  out <- list()
  for (t in in_call) {
    rows <- planned$row[planned$phenotype_name == t]
    out[[t]] <- list(value = unname(value[rows, t]),
                     level = unname(stratum[rows]))
  }
  out
}


#' Liability to phenotype records for one phenotype (in memory)
#'
#' @param r The phenotype's element of the residual adapter's result:
#'   `value`, `level`, `var_unconditional`.
#' @keywords internal
.ap_liability_records <- function(pop, t, m, e, liability, r) {
  pheno_type <- if (is.na(m$type)) "continuous" else m$type
  cat_idx    <- NULL

  if (pheno_type == "categorical") {
    has_thresh <- !is.na(m$thresholds) && nzchar(m$thresholds)
    if (has_thresh) {
      thresh_vec <- as.numeric(strsplit(m$thresholds, ",", fixed = TRUE)[[1]])
    } else {
      # Prevalence threshold on the liability scale: mean + z * sqrt(Va + Ve),
      # with Ve the unconditional (marginal) residual variance.
      if (is.na(r$var_unconditional)) {
        stop("Phenotype '", t, "': the prevalence threshold needs an ",
             "unconditional residual variance, but none is stored for it ",
             "(only conditional strata, or no residual block at all). Add ",
             "an unconditional stratum (define_phenotype(residual_var = ) or ",
             "define_residual_cov() without condition_column) or give ",
             "explicit thresholds (define_phenotype(thresholds = )).",
             call. = FALSE)
      }
      pheno_mean <- if (is.na(m$mean)) 0 else m$mean
      va <- get_trait_var(pop, "gen_add", t)
      va <- if (is.na(va)) 0 else va
      thresh_vec <- pheno_mean +
        stats::qnorm(1 - m$prevalence) * sqrt(va + r$var_unconditional)
    }
    cat_idx <- liability_to_categorical(liability, thresh_vec)
    has_cv  <- !is.na(m$cat_values) && nzchar(m$cat_values)
    value <- if (has_cv) {
      cv <- as.numeric(strsplit(m$cat_values, ",", fixed = TRUE)[[1]])
      as.numeric(cv[cat_idx])
    } else {
      as.numeric(cat_idx)
    }
  } else {
    value <- switch(
      pheno_type,
      continuous = liability,
      count      = as.numeric(clip_count(liability, m$min_value, m$max_value)),
      liability
    )
  }

  records <- tibble::tibble(
    id_ind                   = e$id_ind,
    phenotype_name           = t,
    pheno_value              = as.numeric(value),
    pheno_number             = e$pheno_number,
    residual_value           = as.numeric(r$value),
    residual_condition_level = as.character(r$level)
  )
  if (isTRUE(m$store_liability) && !is.null(cat_idx)) {
    records$liability_value <- as.numeric(liability)
  }
  has_cn <- !is.null(cat_idx) && !is.na(m$cat_names) && nzchar(m$cat_names)
  if (has_cn) {
    cn <- strsplit(m$cat_names, ",", fixed = TRUE)[[1]]
    records$cat_name <- cn[cat_idx]
  }
  records
}


# ── Stage 3: COMMIT ───────────────────────────────────────────────────────────

#' Stage 3: write the resolved call in one transaction
#'
#' `duckdb_register()` + `INSERT` only — never `dbWriteTable()`, which
#' advances R's RNG. Any failure rolls back both tables. The per-phenotype
#' "Wrote ..." messages are emitted after the commit.
#'
#' @param plan The Stage-1 plan (for the path of each phenotype).
#' @param resolved The Stage-2 result.
#' @param extra_cols Scalar extra columns for `ind_phenotype`.
#' @keywords internal
.ap_commit <- function(pop, plan, resolved, extra_cols = list()) {
  conn    <- pop$db_conn
  records <- resolved$records
  re_rows <- resolved$random_effects

  all_rec <- dplyr::bind_rows(records)
  if (nrow(all_rec) > 0) {
    all_rec$id_phenotype <- next_phenotype_ids(pop, nrow(all_rec))
    all_rec <- all_rec[, c("id_phenotype", setdiff(names(all_rec), "id_phenotype"))]
  }

  tmp_re  <- "__ap_random_effects"
  tmp_rec <- "__ap_records"
  DBI::dbExecute(conn, "BEGIN TRANSACTION")
  committed <- FALSE
  on.exit({
    try(duckdb::duckdb_unregister(conn, tmp_re),  silent = TRUE)
    try(duckdb::duckdb_unregister(conn, tmp_rec), silent = TRUE)
    if (!committed) try(DBI::dbExecute(conn, "ROLLBACK"), silent = TRUE)
  }, add = TRUE)

  if (nrow(re_rows) > 0) {
    re_rows$date_sampled <- Sys.Date()
    duckdb::duckdb_register(conn, tmp_re, re_rows)
    DBI::dbExecute(conn, paste0(
      "INSERT INTO phenotype_random_effects ",
      "(phenotype_name, effect_name, level, draw_value, date_sampled) ",
      "SELECT phenotype_name, effect_name, level, draw_value, date_sampled ",
      "FROM ", tmp_re))
  }

  if (nrow(all_rec) > 0) {
    if (length(extra_cols) > 0) {
      prepped <- prepare_extra_cols(extra_cols, nrow(all_rec), "ind_phenotype", conn)
      for (nm in names(prepped)) all_rec[[nm]] <- prepped[[nm]]
    }
    cols <- paste0("\"", names(all_rec), "\"", collapse = ", ")
    duckdb::duckdb_register(conn, tmp_rec, as.data.frame(all_rec))
    DBI::dbExecute(conn, paste0(
      "INSERT INTO ind_phenotype (", cols, ") SELECT ", cols, " FROM ", tmp_rec))
  }

  DBI::dbExecute(conn, "COMMIT")
  committed <- TRUE

  for (t in plan$phenos) {
    n <- if (is.null(records[[t]])) 0L else nrow(records[[t]])
    if (n == 0L) next
    what <- switch(plan$entries[[t]]$path,
                   derived_formula = " derived phenotype records for '",
                   user_values     = " user-supplied phenotype records for '",
                   " phenotype records for '")
    message("Wrote ", n, what, t, "'.")
  }
  invisible(NULL)
}
