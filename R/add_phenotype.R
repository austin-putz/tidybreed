#' Generate phenotype records for a subset of individuals
#'
#' @description
#' Simulates phenotype values for one or more traits and writes them to
#' `ind_phenotype`. Also computes and stores the underlying true breeding
#' value (TBV) per individual per trait in `ind_tbv`.
#'
#' **Model** (per trait, on the liability / continuous scale):
#' \preformatted{
#'   y_i = target_add_mean + sum(fixed_shifts) + sum(random_shifts) + TBV_i + e_i
#' }
#'
#' * `target_add_mean` comes from `trait_meta`.
#' * Fixed and random shifts come from `trait_effects` rows.
#' * `TBV_i` = sum over QTL of `add_{trait}` * genotype dose (or haplotype
#'   dose for imprinted traits).
#' * `e_i` is residual: drawn from `MVN(0, R)` across traits when a residual
#'   covariance matrix is stored (see [add_effect_cov_matrix()]) and multiple
#'   traits share the same subset; otherwise drawn independently.
#'
#' Trait-type specific output:
#' * `"continuous"`: liability written verbatim.
#' * `"count"`: liability rounded and clipped to `[min_value, max_value]`.
#' * `"binary"`: 0/1 via quantile threshold at `1 - prevalence`.
#' * `"categorical"`: integer level via `thresholds` cutpoints.
#'
#' **Subset selection**: pipe a `tidybreed_table` (from [get_table()] and
#' optionally [filter()]) as the first argument. The unique `id_ind` values in
#' the collected table are the candidate individuals. The candidate set is then
#' intersected with `trait_meta$expressed_sex` (e.g. `"F"`-only traits skip
#' males).
#'
#' **Escape hatches**:
#' * `user_values`: skip model computation and write these values as
#'   phenotype records for the subset. For multi-trait calls, supply a
#'   named list keyed by trait.
#' * `user_residual`: supply a numeric vector (or named list for multi-trait)
#'   to override the residual draw.
#'
#' @param tbl A `tidybreed_table` object from [get_table()] (optionally piped
#'   through [filter()]). The table must contain an `id_ind` column; unique
#'   values in that column determine which individuals are phenotyped.
#' @param trait Character vector of trait name(s).
#' @param user_residual Optional override for residual draws. Numeric vector
#'   of length `n_subset` for single trait, or a named list keyed by trait
#'   for multi-trait.
#' @param user_values Optional override for the full phenotype value. If
#'   supplied, the model is not evaluated — these values are written
#'   directly (but TBVs are still computed and stored).
#' @param seed Optional integer for reproducibility.
#'
#' @return The modified `tidybreed_pop` (invisibly).
#'
#' @seealso [add_trait()], [define_qtl()], [set_qtl_effects()],
#'   [add_effect_cov_matrix()], [add_trait_covariate()], [add_tbv()]
#'
#' @examples
#' \dontrun{
#' # All individuals
#' pop <- pop |> get_table("ind_meta") |> add_phenotype("ADG")
#'
#' # Filtered subset
#' pop <- pop |>
#'   get_table("ind_meta") |>
#'   dplyr::filter(sex == "F", gen == 1L) |>
#'   add_phenotype("ADG")
#'
#' # Pre-select by phenotype record value from another table
#' pop <- pop |>
#'   get_table("ind_phenotype") |>
#'   dplyr::filter(value > 500) |>
#'   add_phenotype("ADG2")
#' }
#' @export
add_phenotype <- function(tbl,
                          trait,
                          user_residual = NULL,
                          user_values   = NULL,
                          seed          = NULL) {

  stopifnot(inherits(tbl, "tidybreed_table"))
  pop <- tbl$pop
  validate_tidybreed_pop(pop)
  stopifnot(is.character(trait), length(trait) >= 1)
  lapply(trait, validate_sql_identifier, what = "trait name")

  if (!is.null(seed)) set.seed(seed)

  # 1. Resolve subset: collect filtered table and extract unique id_ind
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

  # 2. Pull candidate ind_meta rows (push filter down when possible)
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

  # 3. Validate each trait and collect metadata
  meta_rows <- DBI::dbGetQuery(
    pop$db_conn,
    paste0("SELECT * FROM trait_meta WHERE trait_name IN (",
           paste0("'", trait, "'", collapse = ", "), ")")
  )
  missing_t <- setdiff(trait, meta_rows$trait_name)
  if (length(missing_t) > 0) {
    stop("Traits not found: ", paste(missing_t, collapse = ", "),
         call. = FALSE)
  }
  meta_rows <- meta_rows[match(trait, meta_rows$trait_name), , drop = FALSE]

  genome_cols <- DBI::dbListFields(pop$db_conn, "genome_meta")
  for (t in trait) {
    if (!paste0("is_QTL_", t) %in% genome_cols) {
      stop("QTL column 'is_QTL_", t, "' not found. Call define_qtl('", t,
           "', ...) first.", call. = FALSE)
    }
    if (!paste0("add_", t) %in% genome_cols) {
      stop("Additive-effect column 'add_", t, "' not found. Call ",
           "set_qtl_effects('", t, "', ...) first.", call. = FALSE)
    }
  }

  # 4. Apply sex filter (intersect subset with expressed_sex)
  # Record per-trait subset (may differ because of sex-limited traits).
  subset_by_trait <- lapply(seq_len(nrow(meta_rows)), function(i) {
    ex_sex <- meta_rows$expressed_sex[i]
    if (ex_sex == "both") return(ind_meta_subset)
    ind_meta_subset[ind_meta_subset$sex == ex_sex, , drop = FALSE]
  })
  names(subset_by_trait) <- trait

  # 5. Pull genome data needed for TBV. Only push a subset filter to DuckDB
  #    when the user actually filtered — a 1-million-wide IN() list for the
  #    full population is slower than a plain SELECT *.
  genome <- dplyr::collect(get_table(pop, "genome_meta"))
  geno_mat_full <- get_genotype_matrix(pop, subset_ids = subset_ids)

  # 6. Compute TBV per trait for each trait's subset + store in ind_tbv
  tbv_by_trait <- list()
  for (t in trait) {
    m <- meta_rows[meta_rows$trait_name == t, ]
    ids_t <- subset_by_trait[[t]]$id_ind
    if (length(ids_t) == 0) {
      tbv_by_trait[[t]] <- numeric(0)
      next
    }
    a <- genome[[paste0("add_", t)]]
    a[is.na(a)] <- 0

    p_base_col <- paste0("base_allele_freq_", t)
    p_base <- if (p_base_col %in% names(genome)) {
      pv <- as.numeric(genome[[p_base_col]])
      pv[is.na(pv)] <- 0
      pv
    } else {
      rep(0, length(a))
    }

    if (m$expressed_parent == "both") {
      rows_idx <- match(ids_t, rownames(geno_mat_full))
      tbv <- as.numeric(geno_mat_full[rows_idx, , drop = FALSE] %*% a) -
             2 * sum(p_base * a)
    } else {
      parent_origin <- if (m$expressed_parent == "parent_1") 1L else 2L
      hap_mat <- get_haplotype_matrix(pop, parent_origin, ids_t)
      tbv <- as.numeric(hap_mat %*% a) - sum(p_base * a)
    }
    tbv_by_trait[[t]] <- stats::setNames(tbv, ids_t)

    tbv_df <- tibble::tibble(
      id_ind     = ids_t,
      trait_name = t,
      tbv        = tbv,
      date_calc  = Sys.Date()
    )
    upsert_ind_tbv(pop, tbv_df)
  }

  # 7. If user_values supplied, short-circuit the model
  if (!is.null(user_values)) {
    write_user_phenotype_values(pop, trait, subset_by_trait, user_values)
    return(invisible(pop))
  }

  # 7.5. Pre-draw correlated random effects (joint MVN) for effects that have
  #      a covariance matrix stored in trait_effect_cov. The draws are written
  #      to trait_random_effects so compute_covariate_contribution() reuses them.
  if ("trait_effect_cov" %in% pop$tables && length(trait) >= 2) {
    traits_sql <- paste0("'", trait, "'", collapse = ", ")
    cov_effects <- DBI::dbGetQuery(
      pop$db_conn,
      paste0("SELECT DISTINCT effect_name FROM trait_effect_cov ",
             "WHERE effect_name NOT IN ('gen_add', 'residual') ",
             "AND trait_1 IN (", traits_sql, ") ",
             "AND trait_2 IN (", traits_sql, ")")
    )$effect_name

    for (eff in cov_effects) {
      eff_traits_q <- DBI::dbGetQuery(
        pop$db_conn,
        paste0("SELECT DISTINCT trait_1 AS t FROM trait_effect_cov ",
               "WHERE effect_name = '", eff, "' ",
               "AND trait_1 IN (", traits_sql, ") ",
               "AND trait_2 IN (", traits_sql, ")")
      )$t
      eff_traits <- intersect(trait, eff_traits_q)
      if (length(eff_traits) < 2) next

      R_eff <- load_effect_cov(pop, eff, eff_traits)
      if (is.null(R_eff)) next

      # Gather all unique levels across all trait subsets for this effect
      eff_rows <- DBI::dbGetQuery(
        pop$db_conn,
        paste0("SELECT trait_name, source_column, source_table ",
               "FROM trait_effects WHERE effect_name = '", eff, "' ",
               "AND trait_name IN (", traits_sql, ")")
      )

      all_levels <- character(0)
      for (et in eff_traits) {
        er <- eff_rows[eff_rows$trait_name == et, , drop = FALSE]
        if (nrow(er) == 0) next
        src_tbl <- if (is.na(er$source_table[1]) ||
                       !nzchar(er$source_table[1])) "ind_meta" else er$source_table[1]
        src_col <- er$source_column[1]
        ids_t   <- subset_by_trait[[et]]$id_ind
        if (length(ids_t) == 0) next
        if (src_tbl == "ind_meta") {
          grp_vals <- subset_by_trait[[et]][[src_col]]
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

      # Find levels that are new (missing) for any trait
      new_levels <- character(0)
      for (et in eff_traits) {
        existing_lvls <- DBI::dbGetQuery(
          pop$db_conn,
          paste0("SELECT level FROM trait_random_effects ",
                 "WHERE trait_name = '", et, "' AND effect_name = '", eff, "'")
        )$level
        new_levels <- union(new_levels, setdiff(all_levels, existing_lvls))
      }
      if (length(new_levels) == 0) next

      # Draw correlated values for new levels
      if (!requireNamespace("MASS", quietly = TRUE)) {
        stop("Package 'MASS' is required for correlated random effect sampling.",
             call. = FALSE)
      }
      draws_mat <- MASS::mvrnorm(
        n     = length(new_levels),
        mu    = rep(0, length(eff_traits)),
        Sigma = R_eff
      )
      if (!is.matrix(draws_mat)) {
        draws_mat <- matrix(draws_mat, nrow = 1)
      }
      colnames(draws_mat) <- eff_traits
      rownames(draws_mat) <- new_levels

      # Upsert into trait_random_effects per trait
      for (et in eff_traits) {
        existing_lvls <- DBI::dbGetQuery(
          pop$db_conn,
          paste0("SELECT level FROM trait_random_effects ",
                 "WHERE trait_name = '", et, "' AND effect_name = '", eff, "'")
        )$level
        new_for_trait <- setdiff(new_levels, existing_lvls)
        if (length(new_for_trait) == 0) next
        new_df <- tibble::tibble(
          trait_name   = et,
          effect_name  = eff,
          level        = new_for_trait,
          draw_value   = draws_mat[new_for_trait, et],
          date_sampled = Sys.Date()
        )
        DBI::dbWriteTable(pop$db_conn, "trait_random_effects", new_df,
                          append = TRUE)
      }
    }
  }

  # 8. Joint residuals if all traits share the same subset + R is stored
  joint_resid <- NULL
  if (length(trait) >= 2 && is.null(user_residual)) {
    subset_ids_list <- lapply(subset_by_trait, function(df) sort(df$id_ind))
    all_equal <- length(unique(subset_ids_list)) == 1
    if (all_equal) {
      R_mat <- load_effect_cov(pop, "residual", trait)
      if (!is.null(R_mat)) {
        n_common <- length(subset_ids_list[[1]])
        var_vec <- stats::setNames(
          vapply(meta_rows$trait_name,
                 function(t) get_effect_var(pop, "residual", t),
                 numeric(1)),
          meta_rows$trait_name
        )
        draws <- sample_residuals(n_common, var_vec, R = R_mat)
        rownames(draws) <- subset_ids_list[[1]]
        joint_resid <- draws
      }
    }
  }

  # 9. Per-trait phenotype generation
  for (t_idx in seq_along(trait)) {
    t <- trait[t_idx]
    m <- meta_rows[meta_rows$trait_name == t, ]
    subset_df <- subset_by_trait[[t]]
    ids_t <- subset_df$id_ind
    n_ind <- length(ids_t)
    if (n_ind == 0) next

    # Covariate (fixed + random) contribution
    covariate_contrib <- compute_covariate_contribution(pop, t, subset_df)

    # Residual
    if (!is.null(user_residual)) {
      resid <- if (is.list(user_residual)) user_residual[[t]] else user_residual
      if (length(resid) != n_ind) {
        stop("user_residual length for trait '", t, "' must equal ",
             n_ind, ".", call. = FALSE)
      }
    } else if (!is.null(joint_resid)) {
      resid <- joint_resid[ids_t, t]
    } else {
      resid_var <- get_effect_var(pop, "residual", t)
      if (is.na(resid_var)) {
        stop("No residual variance found for trait '", t, "'. ",
             "Specify via add_effect_cov_matrix(pop, 'residual', ...) or ",
             "add_trait(pop, '", t, "', residual_var = ...).",
             call. = FALSE)
      }
      resid <- stats::rnorm(n_ind, sd = sqrt(resid_var))
    }

    tbv <- tbv_by_trait[[t]]
    liability <- m$target_add_mean + covariate_contrib + as.numeric(tbv) + resid

    value <- switch(
      m$trait_type,
      continuous  = liability,
      count       = as.numeric(clip_count(liability, m$min_value, m$max_value)),
      binary      = as.numeric(liability_to_binary(liability, m$prevalence)),
      categorical = as.numeric(liability_to_categorical(
                       liability,
                       as.numeric(strsplit(m$thresholds, ",", fixed = TRUE)[[1]]))),
      liability
    )

    records <- tibble::tibble(
      id_record    = next_record_ids(pop, t, n_ind),
      id_ind       = ids_t,
      trait_name   = t,
      value        = as.numeric(value),
      pheno_number = next_pheno_numbers(pop, t, ids_t)
    )
    DBI::dbWriteTable(pop$db_conn, "ind_phenotype", records, append = TRUE)
    message("Wrote ", n_ind, " phenotype records for trait '", t, "'.")
  }

  invisible(pop)
}


#' Write user-supplied phenotype values verbatim
#'
#' @keywords internal
write_user_phenotype_values <- function(pop, trait, subset_by_trait,
                                        user_values) {
  if (length(trait) == 1 && !is.list(user_values)) {
    user_values <- stats::setNames(list(user_values), trait)
  }
  for (t in trait) {
    vals <- user_values[[t]]
    if (is.null(vals)) {
      stop("user_values missing entry for trait '", t, "'.", call. = FALSE)
    }
    ids_t <- subset_by_trait[[t]]$id_ind

    # If the user supplies a named vector, honour those ids
    if (!is.null(names(vals))) {
      ids_t <- names(vals)
      vals <- unname(vals)
    } else if (length(vals) != length(ids_t)) {
      stop("user_values for '", t, "' must have length equal to subset (",
           length(ids_t), ") or be a named vector.", call. = FALSE)
    }
    if (length(vals) == 0) next

    records <- tibble::tibble(
      id_record    = next_record_ids(pop, t, length(vals)),
      id_ind       = ids_t,
      trait_name   = t,
      value        = as.numeric(vals),
      pheno_number = next_pheno_numbers(pop, t, ids_t)
    )
    DBI::dbWriteTable(pop$db_conn, "ind_phenotype", records, append = TRUE)
    message("Wrote ", nrow(records), " user-supplied phenotype records for '",
            t, "'.")
  }
  invisible(NULL)
}


#' Upsert TBV rows (delete existing by id_ind + trait, then insert)
#'
#' @keywords internal
upsert_ind_tbv <- function(pop, tbv_df) {
  if (nrow(tbv_df) == 0) return(invisible(NULL))
  trait <- unique(tbv_df$trait_name)
  stopifnot(length(trait) == 1)
  ids_sql <- paste0("'", tbv_df$id_ind, "'", collapse = ", ")
  DBI::dbExecute(
    pop$db_conn,
    paste0("DELETE FROM ind_tbv WHERE trait_name = '", trait,
           "' AND id_ind IN (", ids_sql, ")")
  )
  DBI::dbWriteTable(pop$db_conn, "ind_tbv", tbv_df, append = TRUE)
  invisible(NULL)
}
