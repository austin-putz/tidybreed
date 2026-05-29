#' Generate phenotype records for a subset of individuals
#'
#' @description
#' Simulates phenotype values for one or more phenotypes and writes them to
#' `ind_phenotype`. Also computes and stores the underlying true breeding
#' value (TBV) per individual per trait in `ind_tbv`.
#'
#' **Model** (per phenotype, on the liability / continuous scale):
#' \preformatted{
#'   y_i = mean + sum(fixed_shifts) + sum(random_shifts) + TBV_i + e_i
#' }
#'
#' * `mean` comes from `phenotype_meta.mean`.
#' * Fixed and random shifts come from `trait_effects` rows.
#' * For **simple** phenotypes (`phenotype_name == trait_name`), `TBV_i` is the
#'   standard additive TBV from `genome_effects`.
#' * For **composite** phenotypes (rows in `phenotype_components`), `TBV_i` is
#'   the weighted sum of contributor TBVs (self, dam, sire, or group).
#' * `e_i` is residual: drawn from `MVN(0, R)` across phenotypes when a
#'   covariance matrix is stored in `phenotype_residual_cov` and multiple
#'   phenotypes share the same subset; otherwise drawn independently.
#'
#' **Subset selection**: pipe a `tidybreed_table` (from [get_table()] and
#' optionally [filter()]) as the first argument.
#'
#' **Escape hatches**:
#' * `user_values`: skip model computation and write these values as phenotype
#'   records for the subset.
#' * `user_residual`: supply a numeric vector (or named list) to override the
#'   residual draw.
#'
#' @param tbl A `tidybreed_table` object from [get_table()] (optionally piped
#'   through [filter()]). The table must contain an `id_ind` column.
#' @param trait_name Character vector of phenotype name(s) (accepts phenotype
#'   names — the argument is named `trait_name` for backward compatibility).
#'   When `NULL` (default), all phenotypes in `phenotype_meta` are used in
#'   `id_phenotype_meta` order.
#' @param user_residual Optional override for residual draws. Numeric vector
#'   of length `n_subset` for single phenotype, or a named list keyed by
#'   phenotype name for multi-phenotype.
#' @param user_values Optional override for the full phenotype value. If
#'   supplied, the model is not evaluated (TBVs are still computed and stored).
#' @param seed Optional integer for reproducibility.
#' @param ... Optional scalar extra columns written to `ind_phenotype`
#'   (broadcast to all records). Supply per-record vectors with
#'   [mutate_table()] after the call.
#'
#' @return The modified `tidybreed_pop` (invisibly).
#'
#' @seealso [define_phenotype()], [define_trait()], [define_additive_effects()],
#'   [define_residual_cov()], [define_effect_cov_matrix()],
#'   [define_effect_fixed_class()], [define_effect_fixed_cov()],
#'   [define_effect_random()], [add_tbv()]
#'
#' @examples
#' \dontrun{
#' # All individuals — all phenotypes
#' pop <- pop |> get_table("ind_meta") |> add_phenotype()
#'
#' # Named phenotype, filtered subset
#' pop <- pop |>
#'   get_table("ind_meta") |>
#'   dplyr::filter(sex == "F", gen == 1L) |>
#'   add_phenotype("ADG")
#' }
#' @export
add_phenotype <- function(tbl,
                          trait_name    = NULL,
                          user_residual = NULL,
                          user_values   = NULL,
                          seed          = NULL,
                          ...) {

  stopifnot(inherits(tbl, "tidybreed_table"))
  pop <- tbl$pop
  validate_tidybreed_pop(pop)

  # Resolve phenotype names from phenotype_meta
  if (is.null(trait_name)) {
    trait_name <- DBI::dbGetQuery(
      pop$db_conn,
      "SELECT phenotype_name FROM phenotype_meta ORDER BY id_phenotype_meta"
    )$phenotype_name
    if (length(trait_name) == 0L)
      stop("No phenotypes found in phenotype_meta. ",
           "Call define_phenotype() first.", call. = FALSE)
  }

  stopifnot(is.character(trait_name), length(trait_name) >= 1)
  lapply(trait_name, validate_sql_identifier, what = "phenotype name")
  phenos <- trait_name

  extra_cols <- list(...)
  if (length(extra_cols) > 0) {
    for (nm in names(extra_cols)) {
      if (length(extra_cols[[nm]]) != 1L) {
        stop("Custom field '", nm, "' in add_phenotype() must be a scalar ",
             "(broadcast to all records). Supply per-record vectors with ",
             "mutate_table() after the call.", call. = FALSE)
      }
    }
  }

  if (!is.null(seed)) set.seed(seed)

  # ── 1. Resolve subset ──────────────────────────────────────────────────────

  if (length(tbl$pending_filter) == 0) {
    subset_ids <- NULL
  } else {
    collected <- dplyr::collect(tbl)
    if (!"id_ind" %in% names(collected)) {
      stop("Filtered table '", tbl$table_name,
           "' must contain 'id_ind' to subset individuals for phenotyping.",
           call. = FALSE)
    }
    subset_ids <- unique(collected[["id_ind"]])
  }

  # ── 2. Pull candidate ind_meta rows ───────────────────────────────────────

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
    return(invisible(pop))
  }

  # ── 3. Validate phenotypes and collect metadata from phenotype_meta ────────

  phenos_in   <- paste0("'", gsub("'", "''", phenos), "'", collapse = ", ")
  pheno_meta  <- DBI::dbGetQuery(
    pop$db_conn,
    paste0("SELECT * FROM phenotype_meta WHERE phenotype_name IN (", phenos_in, ")")
  )
  missing_p <- setdiff(phenos, pheno_meta$phenotype_name)
  if (length(missing_p) > 0) {
    stop("Phenotypes not found in phenotype_meta: ",
         paste(missing_p, collapse = ", "),
         ". Call define_phenotype() first.", call. = FALSE)
  }
  pheno_meta <- pheno_meta[match(phenos, pheno_meta$phenotype_name), , drop = FALSE]

  # ── 4 (G1): Check which phenotypes have components ────────────────────────

  has_components       <- stats::setNames(logical(length(phenos)), phenos)
  components_by_pheno  <- stats::setNames(vector("list", length(phenos)), phenos)

  if ("phenotype_components" %in% DBI::dbListTables(pop$db_conn)) {
    for (t in phenos) {
      comp_rows <- DBI::dbGetQuery(
        pop$db_conn,
        paste0("SELECT * FROM phenotype_components WHERE phenotype_name = '",
               gsub("'", "''", t), "'")
      )
      has_components[t]      <- nrow(comp_rows) > 0
      components_by_pheno[[t]] <- comp_rows
    }
  }

  # ── 4b. Classify formula types for each phenotype ─────────────────────────

  has_formula_tbv <- stats::setNames(logical(length(phenos)), phenos)
  has_formula     <- stats::setNames(logical(length(phenos)), phenos)
  formula_tbv_str <- stats::setNames(character(length(phenos)), phenos)
  formula_str     <- stats::setNames(character(length(phenos)), phenos)

  for (t in phenos) {
    m_t  <- pheno_meta[pheno_meta$phenotype_name == t, ]
    ftbv <- if ("formula_tbv" %in% names(m_t) && !is.na(m_t$formula_tbv) &&
                 nzchar(m_t$formula_tbv)) m_t$formula_tbv else NA_character_
    fder <- if ("formula" %in% names(m_t) && !is.na(m_t$formula) &&
                 nzchar(m_t$formula)) m_t$formula else NA_character_
    has_formula_tbv[t] <- !is.na(ftbv)
    has_formula[t]     <- !is.na(fder)
    formula_tbv_str[t] <- if (is.na(ftbv)) "" else ftbv
    formula_str[t]     <- if (is.na(fder)) "" else fder
  }

  # G1: Validate genome_effects only for simple (non-composite) phenotypes
  for (t in phenos[!has_components & !has_formula_tbv & !has_formula]) {
    n_eff <- DBI::dbGetQuery(
      pop$db_conn,
      paste0("SELECT COUNT(*) AS n FROM genome_effects ",
             "WHERE trait_name = '", gsub("'", "''", t), "' ",
             "AND genome_effect_type = 'additive' AND line_name IS NULL")
    )$n
    if (n_eff == 0L) {
      stop(
        "No additive effects found for phenotype '", t, "' in genome_effects. ",
        "For simple phenotypes call define_additive_effects() first. ",
        "For composite phenotypes supply 'components' or 'formula_tbv' in define_phenotype(). ",
        "For derived phenotypes (no genetic architecture) use trait_type = 'derived_formula'.",
        call. = FALSE
      )
    }
  }

  # ── 4c. Topological sort when derived_formula phenotypes are present ────────

  if (any(has_formula)) {
    ordered_phenos <- .topo_sort_phenotypes(pheno_meta)
    phenos             <- ordered_phenos
    has_components     <- has_components[phenos]
    has_formula_tbv    <- has_formula_tbv[phenos]
    has_formula        <- has_formula[phenos]
    formula_tbv_str    <- formula_tbv_str[phenos]
    formula_str        <- formula_str[phenos]
    components_by_pheno <- components_by_pheno[phenos]
    pheno_meta         <- pheno_meta[match(phenos, pheno_meta$phenotype_name), ]
  }

  # ── 5. Sex filter using phenotype_meta.expressed_sex ──────────────────────

  subset_by_pheno <- lapply(seq_len(nrow(pheno_meta)), function(i) {
    ex_sex <- pheno_meta$expressed_sex[i]
    if (is.null(ex_sex) || is.na(ex_sex) || ex_sex == "both") return(ind_meta_subset)
    ind_meta_subset[ind_meta_subset$sex == ex_sex, , drop = FALSE]
  })
  names(subset_by_pheno) <- phenos

  # ── 5.5. Repeatable guard ─────────────────────────────────────────────────

  for (t in phenos) {
    m_t   <- pheno_meta[pheno_meta$phenotype_name == t, ]
    ids_t <- subset_by_pheno[[t]]$id_ind
    if (length(ids_t) == 0 || isTRUE(m_t$repeatable)) next

    pn_safe  <- gsub("'", "''", t)
    ids_sql  <- paste0("'", ids_t, "'", collapse = ", ")
    already_done <- DBI::dbGetQuery(
      pop$db_conn,
      paste0("SELECT DISTINCT id_ind FROM ind_phenotype ",
             "WHERE phenotype_name = '", pn_safe, "' ",
             "AND id_ind IN (", ids_sql, ")")
    )$id_ind
    n_rejected <- length(already_done)
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

  # ── 6. Compute TBVs ───────────────────────────────────────────────────────

  # Derived formula phenotypes have no TBV; exclude from TBV computation
  simple_phenos        <- phenos[!has_components & !has_formula_tbv & !has_formula]
  composite_phenos     <- phenos[has_components]
  formula_tbv_phenos   <- phenos[has_formula_tbv]

  # Simple: phenotype_name == trait_name in trait_meta
  if (length(simple_phenos) > 0) {
    pop <- add_tbv(tbl, trait_name = simple_phenos)
  }

  # G2: Composite — gather all contributor IDs + source traits, then add_tbv once
  if (length(composite_phenos) > 0) {
    all_source_traits   <- character(0)
    all_contributor_ids <- character(0)

    for (t in composite_phenos) {
      comp_rows <- components_by_pheno[[t]]
      all_source_traits <- unique(c(all_source_traits,
                                    as.character(comp_rows$source_trait_name)))
      subset_df <- subset_by_pheno[[t]]

      for (ct in unique(as.character(comp_rows$contributor_type))) {
        if (ct == "self") {
          all_contributor_ids <- unique(c(all_contributor_ids, subset_df$id_ind))
        } else if (ct == "dam") {
          dids <- as.character(subset_df$id_parent_2)
          dids <- dids[!is.na(dids) & dids != "NA" & nzchar(dids)]
          all_contributor_ids <- unique(c(all_contributor_ids, dids))
        } else if (ct == "sire") {
          sids <- as.character(subset_df$id_parent_1)
          sids <- sids[!is.na(sids) & sids != "NA" & nzchar(sids)]
          all_contributor_ids <- unique(c(all_contributor_ids, sids))
        } else if (ct == "group") {
          # Resolve all group-members so their TBVs can be pre-computed
          grp_comp_rows <- comp_rows[as.character(comp_rows$contributor_type) == "group",
                                     , drop = FALSE]
          for (gi in seq_len(nrow(grp_comp_rows))) {
            grp_row <- grp_comp_rows[gi, ]
            grp_col <- as.character(grp_row$group_column)
            grp_tbl <- if (is.na(grp_row$group_table) || !nzchar(grp_row$group_table))
                         "ind_meta" else as.character(grp_row$group_table)
            focal_sql <- paste0("'", subset_df$id_ind, "'", collapse = ", ")
            grp_vals  <- DBI::dbGetQuery(pop$db_conn, paste0(
              "SELECT DISTINCT \"", grp_col, "\" AS gv FROM ", grp_tbl,
              " WHERE id_ind IN (", focal_sql,
              ") AND \"", grp_col, "\" IS NOT NULL"
            ))$gv
            if (length(grp_vals) > 0) {
              gv_sql     <- paste0("'", grp_vals, "'", collapse = ", ")
              member_ids <- DBI::dbGetQuery(pop$db_conn, paste0(
                "SELECT id_ind FROM ", grp_tbl,
                " WHERE \"", grp_col, "\" IN (", gv_sql, ")"
              ))$id_ind
              all_contributor_ids <- unique(c(all_contributor_ids, member_ids))
            }
          }
        }
      }
    }

    if (length(all_contributor_ids) > 0 && length(all_source_traits) > 0) {
      contrib_tbl <- get_table(pop, "ind_meta") |>
        dplyr::filter(.data$id_ind %in% !!all_contributor_ids)
      pop <- add_tbv(contrib_tbl, trait_name = all_source_traits)
    }
  }

  # G2b: formula_tbv — gather all contributor IDs + source traits via AST walk
  if (length(formula_tbv_phenos) > 0) {
    all_source_traits   <- character(0)
    all_contributor_ids <- character(0)

    for (t in formula_tbv_phenos) {
      ftbv      <- formula_tbv_str[t]
      expr      <- parse(text = ftbv, keep.source = FALSE)[[1]]
      walk_res  <- .walk_formula_tbv_ast(expr)
      subset_df <- subset_by_pheno[[t]]

      all_source_traits <- unique(c(all_source_traits,
        vapply(walk_res$trait_refs, `[[`, character(1), "trait")
      ))

      for (ref in walk_res$trait_refs) {
        if (ref$type == "self") {
          all_contributor_ids <- unique(c(all_contributor_ids,
                                          as.character(subset_df$id_ind)))
        } else if (ref$type == "dam") {
          dids <- as.character(subset_df$id_parent_2)
          dids <- dids[!is.na(dids) & dids != "NA" & nzchar(dids)]
          all_contributor_ids <- unique(c(all_contributor_ids, dids))
        } else if (ref$type == "sire") {
          sids <- as.character(subset_df$id_parent_1)
          sids <- sids[!is.na(sids) & sids != "NA" & nzchar(sids)]
          all_contributor_ids <- unique(c(all_contributor_ids, sids))
        } else if (ref$type %in% c("group_sum", "group_mean")) {
          grp_col   <- ref$col
          grp_tbl   <- ref$table
          focal_sql <- paste0("'", subset_df$id_ind, "'", collapse = ", ")
          grp_vals  <- tryCatch(
            DBI::dbGetQuery(pop$db_conn, paste0(
              "SELECT DISTINCT \"", grp_col, "\" AS gv FROM ", grp_tbl,
              " WHERE id_ind IN (", focal_sql,
              ") AND \"", grp_col, "\" IS NOT NULL"
            ))$gv,
            error = function(e) character(0)  # column/table errors caught later
          )
          if (length(grp_vals) > 0) {
            gv_sql     <- paste0("'", grp_vals, "'", collapse = ", ")
            member_ids <- DBI::dbGetQuery(pop$db_conn, paste0(
              "SELECT id_ind FROM ", grp_tbl,
              " WHERE \"", grp_col, "\" IN (", gv_sql, ")"
            ))$id_ind
            all_contributor_ids <- unique(c(all_contributor_ids, member_ids))
          }
        }
      }
    }

    if (length(all_contributor_ids) > 0 && length(all_source_traits) > 0) {
      contrib_tbl <- get_table(pop, "ind_meta") |>
        dplyr::filter(.data$id_ind %in% !!all_contributor_ids)
      pop <- add_tbv(contrib_tbl, trait_name = all_source_traits)
    }
  }

  # ── 7. user_values short-circuit ──────────────────────────────────────────

  if (!is.null(user_values)) {
    write_user_phenotype_values(pop, phenos, subset_by_pheno, user_values,
                                extra_cols)
    return(invisible(pop))
  }

  # ── 7.5. Pre-draw correlated random effects (joint MVN across phenotypes) ──

  if ("trait_effect_cov" %in% pop$tables && length(phenos) >= 2) {
    phenos_sql  <- paste0("'", gsub("'", "''", phenos), "'", collapse = ", ")
    cov_effects <- DBI::dbGetQuery(
      pop$db_conn,
      paste0("SELECT DISTINCT effect_name FROM trait_effect_cov ",
             "WHERE effect_name NOT IN ('gen_add', 'residual') ",
             "AND trait_name_1 IN (", phenos_sql, ") ",
             "AND trait_name_2 IN (", phenos_sql, ")")
    )$effect_name

    for (eff in cov_effects) {
      eff_safe    <- gsub("'", "''", eff)
      eff_phenos_q <- DBI::dbGetQuery(
        pop$db_conn,
        paste0("SELECT DISTINCT trait_name_1 AS p FROM trait_effect_cov ",
               "WHERE effect_name = '", eff_safe, "' ",
               "AND trait_name_1 IN (", phenos_sql, ") ",
               "AND trait_name_2 IN (", phenos_sql, ")")
      )$p
      eff_phenos <- intersect(phenos, eff_phenos_q)
      if (length(eff_phenos) < 2) next

      R_eff <- load_effect_cov(pop, eff, eff_phenos)
      if (is.null(R_eff)) next

      eff_rows <- DBI::dbGetQuery(
        pop$db_conn,
        paste0("SELECT phenotype_name, source_column, source_table ",
               "FROM trait_effects WHERE effect_name = '", eff_safe, "' ",
               "AND phenotype_name IN (", phenos_sql, ")")
      )

      all_levels <- character(0)
      for (et in eff_phenos) {
        er <- eff_rows[eff_rows$phenotype_name == et, , drop = FALSE]
        if (nrow(er) == 0) next
        src_tbl <- if (is.na(er$source_table[1]) ||
                       !nzchar(er$source_table[1])) "ind_meta" else er$source_table[1]
        src_col <- er$source_column[1]
        ids_t   <- subset_by_pheno[[et]]$id_ind
        if (length(ids_t) == 0) next
        if (src_tbl == "ind_meta") {
          grp_vals <- subset_by_pheno[[et]][[src_col]]
        } else {
          ids_sql_et <- paste0("'", ids_t, "'", collapse = ", ")
          grp_df <- DBI::dbGetQuery(
            pop$db_conn,
            paste0("SELECT ", src_col, " FROM ", src_tbl,
                   " WHERE id_ind IN (", ids_sql_et, ")")
          )
          grp_vals <- grp_df[[src_col]]
        }
        all_levels <- union(all_levels,
                            unique(as.character(grp_vals[!is.na(grp_vals)])))
      }
      if (length(all_levels) == 0) next

      new_levels <- character(0)
      for (et in eff_phenos) {
        pn_safe2 <- gsub("'", "''", et)
        existing_lvls <- DBI::dbGetQuery(
          pop$db_conn,
          paste0("SELECT level FROM trait_random_effects ",
                 "WHERE phenotype_name = '", pn_safe2,
                 "' AND effect_name = '", eff_safe, "'")
        )$level
        new_levels <- union(new_levels, setdiff(all_levels, existing_lvls))
      }
      if (length(new_levels) == 0) next

      if (!requireNamespace("MASS", quietly = TRUE)) {
        stop("Package 'MASS' is required for correlated random effect sampling.",
             call. = FALSE)
      }
      draws_mat <- MASS::mvrnorm(
        n     = length(new_levels),
        mu    = rep(0, length(eff_phenos)),
        Sigma = R_eff
      )
      if (!is.matrix(draws_mat)) draws_mat <- matrix(draws_mat, nrow = 1)
      colnames(draws_mat) <- eff_phenos
      rownames(draws_mat) <- new_levels

      for (et in eff_phenos) {
        pn_safe2 <- gsub("'", "''", et)
        existing_lvls <- DBI::dbGetQuery(
          pop$db_conn,
          paste0("SELECT level FROM trait_random_effects ",
                 "WHERE phenotype_name = '", pn_safe2,
                 "' AND effect_name = '", eff_safe, "'")
        )$level
        new_for_pheno <- setdiff(new_levels, existing_lvls)
        if (length(new_for_pheno) == 0) next
        new_df <- tibble::tibble(
          phenotype_name = et,
          effect_name    = eff,
          level          = new_for_pheno,
          draw_value     = draws_mat[new_for_pheno, et],
          date_sampled   = Sys.Date()
        )
        DBI::dbWriteTable(pop$db_conn, "trait_random_effects", new_df,
                          append = TRUE)
      }
    }
  }

  # ── 8. Residual covariance info (G3) ──────────────────────────────────────

  resid_info <- get_residual_cov(pop, phenos, ind_meta_subset)

  # ── 8.5. Joint residual draw for multi-phenotype case (G3 + G5) ───────────

  joint_resid <- NULL

  if (length(phenos) >= 2 && is.null(user_residual)) {
    subset_ids_list <- lapply(subset_by_pheno, function(df) sort(df$id_ind))
    id_strs   <- vapply(subset_ids_list, paste, character(1), collapse = "|")
    all_equal <- length(unique(id_strs)) == 1

    if (all_equal && !is.null(resid_info$R_unconditional)) {
      n_common <- length(subset_ids_list[[1]])

      if (is.null(resid_info$condition_column)) {
        # Fast path: single unconditional MVN draw
        var_vec <- stats::setNames(diag(resid_info$R_unconditional), phenos)
        draws   <- sample_residuals(n_common, var_vec,
                                    R = resid_info$R_unconditional)
        rownames(draws) <- subset_ids_list[[1]]
        joint_resid <- draws

      } else if (!is.null(resid_info$R_by_level)) {
        # G5: Heterogeneous residuals — draw per condition group
        cond_col <- resid_info$condition_column
        cond_tbl <- resid_info$condition_table %||% "ind_meta"
        common_ids <- subset_ids_list[[1]]

        if (cond_tbl == "ind_meta") {
          ref_df    <- subset_by_pheno[[phenos[1]]]
          cond_vals <- ref_df[[cond_col]][match(common_ids, ref_df$id_ind)]
        } else {
          ids_sql <- paste0("'", common_ids, "'", collapse = ", ")
          cv_df   <- DBI::dbGetQuery(
            pop$db_conn,
            paste0("SELECT id_ind, ", cond_col, " FROM ", cond_tbl,
                   " WHERE id_ind IN (", ids_sql, ")")
          )
          cond_vals <- cv_df[[cond_col]][match(common_ids, cv_df$id_ind)]
        }

        draws_mat <- matrix(NA_real_, n_common, length(phenos),
                            dimnames = list(common_ids, phenos))

        for (lvl in names(resid_info$R_by_level)) {
          grp_mask <- !is.na(cond_vals) & as.character(cond_vals) == lvl
          if (!any(grp_mask)) next
          grp_ids <- common_ids[grp_mask]
          R_grp   <- resid_info$R_by_level[[lvl]]
          var_vec <- stats::setNames(diag(R_grp), phenos)
          grp_draws <- sample_residuals(length(grp_ids), var_vec, R = R_grp)
          draws_mat[grp_ids, ] <- grp_draws
        }

        # Individuals without a level match: fall back to unconditional R
        no_match <- which(rowSums(is.na(draws_mat)) > 0)
        if (length(no_match) > 0) {
          if (!is.null(resid_info$R_unconditional)) {
            warning(length(no_match),
                    " individual(s) have no matching condition level in ",
                    "phenotype_residual_cov; using unconditional R.",
                    call. = FALSE)
            var_vec <- stats::setNames(diag(resid_info$R_unconditional), phenos)
            fb_draws <- sample_residuals(length(no_match), var_vec,
                                         R = resid_info$R_unconditional)
            draws_mat[no_match, ] <- fb_draws
          } else {
            warning(length(no_match),
                    " individual(s) have no matching condition level; residuals ",
                    "set to 0 (no unconditional fallback).", call. = FALSE)
            draws_mat[no_match, ] <- 0
          }
        }

        joint_resid <- draws_mat
      }
    }
  }

  # ── 9. Per-phenotype generation loop ──────────────────────────────────────

  for (t_idx in seq_along(phenos)) {
    t         <- phenos[t_idx]
    m         <- pheno_meta[pheno_meta$phenotype_name == t, ]
    subset_df <- subset_by_pheno[[t]]
    ids_t     <- subset_df$id_ind
    n_ind     <- length(ids_t)
    if (n_ind == 0) next

    # 9a. Covariate contribution (fixed + random effects)
    cov_result        <- compute_covariate_contribution(pop, t, subset_df)
    covariate_contrib <- cov_result$contribution

    # G4: Remove individuals skipped by null_class_action = "skip"
    skip_mask <- is.na(covariate_contrib)
    if (any(skip_mask)) {
      n_skip <- sum(skip_mask)
      warning("Phenotype '", t, "': ", n_skip,
              " individual(s) excluded due to null_class_action = 'skip'.",
              call. = FALSE)
      keep_mask         <- !skip_mask
      ids_t             <- ids_t[keep_mask]
      n_ind             <- length(ids_t)
      subset_df         <- subset_df[keep_mask, , drop = FALSE]
      covariate_contrib <- covariate_contrib[keep_mask]
      if (n_ind == 0) {
        message("Phenotype '", t, "': all individuals skipped; no records written.")
        next
      }
    }

    # 9b. TBV: four-way dispatch
    if (has_formula[t]) {
      # ── derived_formula: evaluate arithmetic over ind_phenotype records ──
      derived_vals <- .eval_derived_formula(pop, formula_str[t], ids_t, t)

      # Write records directly — no residual, no liability conversion
      records <- tibble::tibble(
        id_phenotype   = next_phenotype_ids(pop, n_ind),
        id_ind         = ids_t,
        phenotype_name = t,
        pheno_value    = as.numeric(derived_vals),
        pheno_number   = next_pheno_numbers(pop, t, ids_t)
      )
      if (length(extra_cols) > 0) {
        prepped <- prepare_extra_cols(extra_cols, nrow(records), "ind_phenotype",
                                     pop$db_conn)
        for (nm in names(prepped)) records[[nm]] <- prepped[[nm]]
      }
      DBI::dbWriteTable(pop$db_conn, "ind_phenotype", records, append = TRUE)
      message("Wrote ", n_ind, " derived phenotype records for '", t, "'.")
      next  # skip residual draw, liability conversion, etc.

    } else if (has_formula_tbv[t]) {
      # ── formula_tbv composite: evaluate DSL formula ──
      mca <- if ("missing_component_action" %in% names(m) &&
                 !is.na(m$missing_component_action) &&
                 nzchar(m$missing_component_action)) m$missing_component_action else "skip"

      raw_tbv <- .eval_formula_tbv(pop, formula_tbv_str[t], subset_df, t)
      tbv     <- raw_tbv[ids_t]

      excl_mask <- is.na(tbv)
      if (any(excl_mask)) {
        n_excl    <- sum(excl_mask)
        excl_ids  <- ids_t[excl_mask]
        msg <- paste0(
          n_excl, " individual(s) had one or more missing components for phenotype '",
          t, "' (formula_tbv: ", formula_tbv_str[t], ") and were excluded. ",
          "(IDs: ", paste(head(excl_ids, 5), collapse = ", "),
          if (n_excl > 5) paste0(" ... +", n_excl - 5L, " more") else "", ")"
        )
        if (mca == "error") stop(msg, call. = FALSE)
        else warning(msg, call. = FALSE)

        ids_t             <- ids_t[!excl_mask]
        n_ind             <- length(ids_t)
        subset_df         <- subset_df[!excl_mask, , drop = FALSE]
        covariate_contrib <- covariate_contrib[!excl_mask]
        tbv               <- tbv[ids_t]
        if (n_ind == 0) {
          message("Phenotype '", t, "': all individuals excluded (no formula_tbv result).")
          next
        }
      }
      # Fall through to 9c (residual draw) — same as simple/composite

    } else if (has_components[t]) {
      # ── components data frame path (unchanged) ──
      mca <- if ("missing_component_action" %in% names(m) &&
                 !is.na(m$missing_component_action) &&
                 nzchar(m$missing_component_action)) m$missing_component_action else "skip"
      comp_result <- .assemble_composite_tbv(pop, t, components_by_pheno[[t]],
                                             subset_df,
                                             missing_component_action = mca)
      tbv <- comp_result$composite_tbv[ids_t]

      excl_mask <- is.na(tbv)
      if (any(excl_mask)) {
        ids_t             <- ids_t[!excl_mask]
        n_ind             <- length(ids_t)
        subset_df         <- subset_df[!excl_mask, , drop = FALSE]
        covariate_contrib <- covariate_contrib[!excl_mask]
        tbv               <- tbv[ids_t]
        if (n_ind == 0) {
          message("Phenotype '", t, "': all individuals excluded (no composite TBV).")
          next
        }
      }
    } else {
      # ── Simple phenotype: read TBV from ind_tbv ──
      ids_sql  <- paste0("'", ids_t, "'", collapse = ", ")
      t_safe   <- gsub("'", "''", t)
      tbv_rows <- DBI::dbGetQuery(
        pop$db_conn,
        paste0("SELECT id_ind, tbv_value FROM ind_tbv ",
               "WHERE trait_name = '", t_safe, "' ",
               "AND id_ind IN (", ids_sql, ")")
      )
      tbv <- stats::setNames(tbv_rows$tbv_value, tbv_rows$id_ind)[ids_t]
    }

    # 9c. Residual draw
    if (!is.null(user_residual)) {
      resid <- if (is.list(user_residual)) user_residual[[t]] else user_residual
      if (length(resid) != n_ind)
        stop("user_residual length for phenotype '", t, "' must equal ",
             n_ind, ".", call. = FALSE)

    } else if (!is.null(joint_resid) && t %in% colnames(joint_resid)) {
      resid <- joint_resid[ids_t, t]

    } else {
      # Independent per-phenotype draw
      resid_var <- resid_info$residual_var_unconditional[t]
      if (is.na(resid_var)) {
        # Backward-compat fallback for databases without phenotype_residual_cov
        resid_var <- get_effect_var(pop, "residual", t)
      }
      if (is.na(resid_var)) {
        stop("No residual variance found for phenotype '", t, "'. ",
             "Specify via define_phenotype(residual_var = ...) or ",
             "define_residual_cov().", call. = FALSE)
      }
      resid <- stats::rnorm(n_ind, sd = sqrt(resid_var))
    }

    # 9d. Liability
    pheno_mean <- if (!is.null(m$mean) && !is.na(m$mean)) m$mean else 0
    liability  <- pheno_mean + covariate_contrib + as.numeric(tbv) + resid

    # 9e. Phenotype-type conversion
    trait_type <- if (is.null(m$trait_type) || is.na(m$trait_type)) {
      "continuous"
    } else {
      m$trait_type
    }

    cat_idx <- NULL

    if (trait_type == "categorical") {
      has_thresh <- !is.null(m$thresholds) && !is.na(m$thresholds) &&
                   nzchar(m$thresholds)
      if (has_thresh) {
        thresh_vec <- as.numeric(strsplit(m$thresholds, ",", fixed = TRUE)[[1]])
      } else {
        va <- get_effect_var(pop, "gen_add", t)
        vr <- resid_info$residual_var_unconditional[t]
        if (is.na(vr)) vr <- get_effect_var(pop, "residual", t)
        va <- if (is.na(va)) 0 else va
        vr <- if (is.na(vr)) 0 else vr
        thresh_vec <- pheno_mean +
                      stats::qnorm(1 - m$prevalence) * sqrt(va + vr)
      }
      cat_idx <- liability_to_categorical(liability, thresh_vec)
      has_cv  <- !is.null(m$cat_values) && !is.na(m$cat_values) &&
                 nzchar(m$cat_values)
      value <- if (has_cv) {
        cv <- as.numeric(strsplit(m$cat_values, ",", fixed = TRUE)[[1]])
        as.numeric(cv[cat_idx])
      } else {
        as.numeric(cat_idx)
      }
    } else {
      value <- switch(
        trait_type,
        continuous = liability,
        count      = as.numeric(clip_count(liability, m$min_value, m$max_value)),
        liability
      )
    }

    # 9f. Build and write records
    records <- tibble::tibble(
      id_phenotype   = next_phenotype_ids(pop, n_ind),
      id_ind         = ids_t,
      phenotype_name = t,
      pheno_value    = as.numeric(value),
      pheno_number   = next_pheno_numbers(pop, t, ids_t)
    )
    if (length(extra_cols) > 0) {
      prepped <- prepare_extra_cols(extra_cols, nrow(records), "ind_phenotype",
                                   pop$db_conn)
      for (nm in names(prepped)) records[[nm]] <- prepped[[nm]]
    }

    # Optional: store raw liability for categorical traits
    if (isTRUE(m$store_liability) && !is.null(cat_idx)) {
      pheno_cols <- DBI::dbListFields(pop$db_conn, "ind_phenotype")
      if (!"liability_value" %in% pheno_cols)
        DBI::dbExecute(pop$db_conn,
          "ALTER TABLE ind_phenotype ADD COLUMN liability_value DOUBLE")
      records$liability_value <- as.numeric(liability)
    }

    # Optional: store category label for categorical traits
    has_cn <- !is.null(cat_idx) &&
              !is.null(m$cat_names) && !is.na(m$cat_names) && nzchar(m$cat_names)
    if (has_cn) {
      cn <- strsplit(m$cat_names, ",", fixed = TRUE)[[1]]
      pheno_cols <- DBI::dbListFields(pop$db_conn, "ind_phenotype")
      if (!"cat_name" %in% pheno_cols)
        DBI::dbExecute(pop$db_conn,
          "ALTER TABLE ind_phenotype ADD COLUMN cat_name VARCHAR")
      records$cat_name <- cn[cat_idx]
    }

    DBI::dbWriteTable(pop$db_conn, "ind_phenotype", records, append = TRUE)
    message("Wrote ", n_ind, " phenotype records for '", t, "'.")
  }

  invisible(pop)
}


# ── Internal helpers ───────────────────────────────────────────────────────────

#' Assemble composite TBV from phenotype_components for a single phenotype
#'
#' Reads contributor TBVs from `ind_tbv` (which must already be populated by
#' `add_tbv()` for all relevant source traits and contributor IDs) and
#' multiplies by component weights, summing across all components.
#'
#' @param pop A `tidybreed_pop` object.
#' @param phenotype_name Character. The composite phenotype name.
#' @param comp_rows Data frame. Rows from `phenotype_components` for this phenotype.
#' @param subset_df Data frame. Sex-filtered (and skip-masked) `ind_meta` rows.
#' @return A list with `composite_tbv`: named numeric vector (NA = excluded).
#' @keywords internal
.assemble_composite_tbv <- function(pop, phenotype_name, comp_rows, subset_df,
                                    missing_component_action = "skip") {
  focal_ids <- subset_df$id_ind
  n         <- length(focal_ids)
  composite <- stats::setNames(rep(0, n), focal_ids)

  for (i in seq_len(nrow(comp_rows))) {
    comp           <- comp_rows[i, ]
    source_trait   <- as.character(comp$source_trait_name)
    contr_type     <- as.character(comp$contributor_type)
    weight_val     <- if (is.null(comp$weight)       || is.na(comp$weight))       1.0     else as.numeric(comp$weight)
    weight_type_c  <- if (is.null(comp$weight_type)  || is.na(comp$weight_type))  "fixed" else as.character(comp$weight_type)

    # Compute effective per-individual weight (used by all contributor types)
    if (weight_type_c == "covariate" &&
        !is.null(comp$covariate_name) && !is.na(comp$covariate_name) &&
        nzchar(comp$covariate_name)) {
      cov_col <- as.character(comp$covariate_name)
      cov_tbl <- if (is.null(comp$covariate_table) || is.na(comp$covariate_table) ||
                     !nzchar(comp$covariate_table)) "ind_meta"
                 else as.character(comp$covariate_table)
      if (cov_tbl == "ind_meta") {
        cov_vals <- as.numeric(subset_df[[cov_col]])
      } else {
        ids_sql2 <- paste0("'", focal_ids, "'", collapse = ", ")
        cv_df    <- DBI::dbGetQuery(pop$db_conn, paste0(
          "SELECT id_ind, ", cov_col, " FROM ", cov_tbl,
          " WHERE id_ind IN (", ids_sql2, ")"
        ))
        cov_vals <- as.numeric(cv_df[[cov_col]][match(focal_ids, cv_df$id_ind)])
      }
      eff_weight <- cov_vals * weight_val
    } else {
      eff_weight <- rep(weight_val, n)
    }

    # Resolve contributor IDs (parallel to focal_ids, NA where missing)
    if (contr_type == "self") {
      contr_ids <- focal_ids
    } else if (contr_type == "dam") {
      contr_ids <- as.character(subset_df$id_parent_2)
    } else if (contr_type == "sire") {
      contr_ids <- as.character(subset_df$id_parent_1)
    } else if (contr_type == "group") {
      # ── Group contributor (SGE / Bijma model) ───────────────────────────────
      grp_col    <- as.character(comp$group_column)
      grp_tbl    <- if (is.na(comp$group_table) || !nzchar(comp$group_table))
                      "ind_meta" else as.character(comp$group_table)
      agg_method <- if (is.na(comp$aggregation) || !nzchar(comp$aggregation))
                      "sum" else as.character(comp$aggregation)
      st_safe    <- gsub("'", "''", source_trait)

      focal_sql    <- paste0("'", focal_ids, "'", collapse = ", ")
      focal_grp_df <- DBI::dbGetQuery(pop$db_conn, paste0(
        "SELECT id_ind, \"", grp_col, "\" AS group_val FROM ", grp_tbl,
        " WHERE id_ind IN (", focal_sql, ")"
      ))
      focal_grp_map <- stats::setNames(
        as.character(focal_grp_df$group_val),
        focal_grp_df$id_ind
      )

      non_na_grps <- unique(focal_grp_map[
        !is.na(focal_grp_map) & nzchar(focal_grp_map)
      ])

      all_members_df <- if (length(non_na_grps) > 0) {
        gv_sql <- paste0("'", non_na_grps, "'", collapse = ", ")
        DBI::dbGetQuery(pop$db_conn, paste0(
          "SELECT id_ind, \"", grp_col, "\" AS group_val FROM ", grp_tbl,
          " WHERE \"", grp_col, "\" IN (", gv_sql, ")"
        ))
      } else {
        data.frame(id_ind = character(0), group_val = character(0),
                   stringsAsFactors = FALSE)
      }
      all_members_df$group_val <- as.character(all_members_df$group_val)

      all_member_ids <- unique(all_members_df$id_ind)
      tbv_map_grp    <- stats::setNames(numeric(0), character(0))
      if (length(all_member_ids) > 0) {
        mem_sql  <- paste0("'", all_member_ids, "'", collapse = ", ")
        tbv_rows <- DBI::dbGetQuery(pop$db_conn, paste0(
          "SELECT id_ind, tbv_value FROM ind_tbv WHERE trait_name = '",
          st_safe, "' AND id_ind IN (", mem_sql, ")"
        ))
        if (nrow(tbv_rows) > 0)
          tbv_map_grp <- stats::setNames(tbv_rows$tbv_value, tbv_rows$id_ind)
      }

      for (j in seq_len(n)) {
        fid <- focal_ids[j]
        if (is.na(composite[fid])) next  # already excluded

        grp_val <- focal_grp_map[fid]
        if (is.na(grp_val) || !nzchar(grp_val)) {
          composite[fid] <- NA_real_   # no group assignment
          next
        }

        mate_ids   <- all_members_df$id_ind[
          all_members_df$group_val == grp_val & all_members_df$id_ind != fid
        ]
        valid_tbvs <- tbv_map_grp[mate_ids]
        valid_tbvs <- valid_tbvs[!is.na(valid_tbvs)]

        social_val <- if (length(valid_tbvs) == 0) 0
                      else if (agg_method == "mean") mean(valid_tbvs)
                      else sum(valid_tbvs)

        composite[fid] <- composite[fid] + eff_weight[j] * social_val
      }
      next   # skip the generic contr_ids path below
    } else {
      warning("contributor_type '", contr_type, "' is not yet implemented; ",
              "skipping component '", source_trait, "'.", call. = FALSE)
      next
    }

    missing_contr <- is.na(contr_ids) | contr_ids == "NA" | !nzchar(contr_ids)

    # Fetch TBVs for non-missing contributors
    notna_ids <- unique(contr_ids[!missing_contr])
    tbv_map   <- stats::setNames(numeric(0), character(0))

    if (length(notna_ids) > 0) {
      st_safe  <- gsub("'", "''", source_trait)
      ids_sql  <- paste0("'", notna_ids, "'", collapse = ", ")
      tbv_rows <- DBI::dbGetQuery(pop$db_conn, paste0(
        "SELECT id_ind, tbv_value FROM ind_tbv ",
        "WHERE trait_name = '", st_safe, "' AND id_ind IN (", ids_sql, ")"
      ))
      if (nrow(tbv_rows) > 0) {
        tbv_map <- stats::setNames(tbv_rows$tbv_value, tbv_rows$id_ind)
      }
    }

    # Add contribution per focal individual
    for (j in seq_len(n)) {
      fid <- focal_ids[j]
      if (is.na(composite[fid])) next  # already excluded by a prior component

      if (missing_contr[j]) {
        composite[fid] <- NA_real_
        next
      }
      cid     <- contr_ids[j]
      tbv_val <- tbv_map[cid]
      if (is.na(tbv_val)) {
        composite[fid] <- NA_real_
        next
      }
      composite[fid] <- composite[fid] + eff_weight[j] * tbv_val
    }
  }

  # Generic missing-component handler (covers group, dam, sire, any future type)
  n_missing <- sum(is.na(composite))
  if (n_missing > 0) {
    missing_ids <- focal_ids[is.na(composite)]
    msg <- paste0(
      n_missing, " individual(s) had one or more missing components for phenotype '",
      phenotype_name, "' and were excluded.",
      " (IDs: ", paste(head(missing_ids, 5), collapse = ", "),
      if (n_missing > 5) paste0(" ... +", n_missing - 5L, " more") else "", ")"
    )
    if (missing_component_action == "error") stop(msg, call. = FALSE)
    else warning(msg, call. = FALSE)
  }

  list(composite_tbv = composite)
}


#' Write user-supplied phenotype values verbatim
#' @keywords internal
write_user_phenotype_values <- function(pop, phenos, subset_by_pheno,
                                        user_values, extra_cols = list()) {
  if (length(phenos) == 1 && !is.list(user_values)) {
    user_values <- stats::setNames(list(user_values), phenos)
  }
  for (t in phenos) {
    vals <- user_values[[t]]
    if (is.null(vals)) {
      stop("user_values missing entry for phenotype '", t, "'.", call. = FALSE)
    }
    ids_t <- subset_by_pheno[[t]]$id_ind

    if (!is.null(names(vals))) {
      ids_t <- names(vals)
      vals  <- unname(vals)
    } else if (length(vals) != length(ids_t)) {
      stop("user_values for '", t, "' must have length equal to subset (",
           length(ids_t), ") or be a named vector.", call. = FALSE)
    }
    if (length(vals) == 0) next

    records <- tibble::tibble(
      id_phenotype   = next_phenotype_ids(pop, length(vals)),
      id_ind         = ids_t,
      phenotype_name = t,
      pheno_value    = as.numeric(vals),
      pheno_number   = next_pheno_numbers(pop, t, ids_t)
    )
    if (length(extra_cols) > 0) {
      prepped <- prepare_extra_cols(extra_cols, nrow(records), "ind_phenotype",
                                   pop$db_conn)
      for (nm in names(prepped)) records[[nm]] <- prepped[[nm]]
    }
    DBI::dbWriteTable(pop$db_conn, "ind_phenotype", records, append = TRUE)
    message("Wrote ", nrow(records), " user-supplied phenotype records for '",
            t, "'.")
  }
  invisible(NULL)
}


#' Upsert TBV rows — column-named INSERT with ON CONFLICT DO UPDATE SET
#'
#' Only the columns present in tbv_df are updated on conflict; all other
#' columns (user-defined extras) are left untouched in existing rows.
#'
#' @keywords internal
upsert_ind_tbv <- function(pop, tbv_df) {
  if (nrow(tbv_df) == 0) return(invisible(NULL))

  start  <- next_int_id(pop$db_conn, "ind_tbv", "id_tbv")
  tbv_df <- tibble::add_column(tbv_df,
                                id_tbv = seq.int(start, start + nrow(tbv_df) - 1L),
                                .before = 1)

  tmp <- paste0("_tbv_tmp_", as.character(round(as.numeric(Sys.time()) * 1000)))
  duckdb::duckdb_register(pop$db_conn, tmp, as.data.frame(tbv_df))
  on.exit(duckdb::duckdb_unregister(pop$db_conn, tmp), add = TRUE)

  cols        <- names(tbv_df)
  key_cols    <- c("id_tbv", "id_ind", "trait_name")
  update_cols <- setdiff(cols, key_cols)

  col_list   <- paste(cols, collapse = ", ")
  update_set <- paste(
    paste0(update_cols, " = EXCLUDED.", update_cols),
    collapse = ", "
  )
  DBI::dbExecute(pop$db_conn, paste0(
    "INSERT INTO ind_tbv (", col_list, ") ",
    "SELECT ", col_list, " FROM ", tmp, " ",
    "ON CONFLICT (id_ind, trait_name) DO UPDATE SET ", update_set
  ))
  invisible(NULL)
}
