#' Encode a named numeric vector as a small JSON string
#'
#' @param x Named numeric vector.
#' @return A single character string, or `NA_character_` if `x` is empty.
#' @keywords internal
encode_levels_json <- function(x) {
  if (is.null(x) || length(x) == 0) return(NA_character_)
  stopifnot(!is.null(names(x)))
  keys <- names(x)
  bad <- grepl("[\",:]", keys)
  if (any(bad)) {
    stop("Level names must not contain \" , or : (offending: ",
         paste(keys[bad], collapse = ", "), ")", call. = FALSE)
  }
  vals <- format(as.numeric(x), scientific = FALSE, trim = TRUE)
  paste0("{", paste0("\"", keys, "\":", vals, collapse = ","), "}")
}


#' Decode a level-map JSON string produced by `encode_levels_json()`
#'
#' @param s Character scalar, or `NA`.
#' @return Named numeric vector (possibly length 0), or `NULL` if missing/empty.
#' @keywords internal
decode_levels_json <- function(s) {
  if (is.na(s) || !nzchar(s)) return(NULL)
  inner <- sub("^\\{", "", sub("\\}$", "", s))
  if (!nzchar(inner)) return(stats::setNames(numeric(0), character(0)))
  pairs <- strsplit(inner, ",", fixed = TRUE)[[1]]
  split_on_colon <- strsplit(pairs, ":", fixed = TRUE)
  keys <- vapply(split_on_colon, function(p) gsub("^\"|\"$", "", p[1]), character(1))
  vals <- vapply(split_on_colon, function(p) as.numeric(p[2]), numeric(1))
  stats::setNames(vals, keys)
}


#' Sample residuals for one or more phenotypes
#'
#' @param n_ind Integer. Number of individuals.
#' @param residual_var Named numeric vector of per-phenotype residual variances.
#' @param R Optional covariance matrix with row/col names matching
#'   `names(residual_var)`.
#' @return Numeric matrix of shape `n_ind × length(residual_var)`, columns
#'   named by phenotype.
#' @keywords internal
sample_residuals <- function(n_ind, residual_var, R = NULL) {
  phenos <- names(residual_var)
  if (is.null(R)) {
    out <- vapply(
      phenos,
      function(p) stats::rnorm(n_ind, sd = sqrt(residual_var[[p]])),
      numeric(n_ind)
    )
    if (!is.matrix(out)) out <- matrix(out, nrow = n_ind, ncol = length(phenos))
    colnames(out) <- phenos
    return(out)
  }
  if (!requireNamespace("MASS", quietly = TRUE)) {
    stop("Package 'MASS' is required for joint residual sampling.", call. = FALSE)
  }
  R <- R[phenos, phenos, drop = FALSE]
  draws <- MASS::mvrnorm(n = n_ind, mu = rep(0, length(phenos)), Sigma = R)
  if (is.null(dim(draws))) draws <- matrix(draws, nrow = 1)
  colnames(draws) <- phenos
  draws
}


#' Convert a liability vector to ordered integer categories
#'
#' @param liability Numeric vector.
#' @param thresholds Numeric vector of cutpoints (ascending).
#' @return Integer vector of category indices (1-based).
#' @keywords internal
liability_to_categorical <- function(liability, thresholds) {
  thresholds <- sort(thresholds)
  findInterval(liability, thresholds) + 1L
}


#' Clip a count vector to [min_value, max_value]
#'
#' @param x Numeric vector.
#' @param min_value,max_value Numeric scalars or `NA`.
#' @return Integer vector (rounded) clipped to bounds.
#' @keywords internal
clip_count <- function(x, min_value = NA_real_, max_value = NA_real_) {
  y <- round(x)
  if (!is.na(min_value)) y <- pmax(y, min_value)
  if (!is.na(max_value)) y <- pmin(y, max_value)
  as.integer(y)
}


#' Compute the covariate contribution for each individual for a single phenotype
#'
#' @param pop A `tidybreed_pop` object.
#' @param phenotype_name Character. Phenotype name to look up in `trait_effects`.
#' @param subset_df Data frame: the per-phenotype subset of `ind_meta` (already
#'   sex-filtered). Must contain `id_ind` and any `ind_meta` columns referenced
#'   by effects.
#' @return A list:
#'   - `contribution`: numeric vector of length `nrow(subset_df)`. `NA` for
#'     individuals excluded by `null_class_action = "skip"`.
#'   - `n_skipped`: integer count of NAs introduced.
#' @keywords internal
compute_covariate_contribution <- function(pop, phenotype_name, subset_df) {
  pn_safe <- gsub("'", "''", phenotype_name)
  effects  <- DBI::dbGetQuery(
    pop$db_conn,
    paste0("SELECT * FROM trait_effects WHERE phenotype_name = '", pn_safe, "'")
  )
  n_ind <- nrow(subset_df)
  if (nrow(effects) == 0) return(list(contribution = rep(0, n_ind), n_skipped = 0L))

  ids_t   <- subset_df$id_ind
  ids_sql <- paste0("'", ids_t, "'", collapse = ", ")
  total   <- rep(0, n_ind)
  skip_mask <- rep(FALSE, n_ind)  # tracks which individuals to skip

  for (i in seq_len(nrow(effects))) {
    e <- effects[i, ]

    src_tbl <- if (is.na(e$source_table) || !nzchar(e$source_table)) {
      "ind_meta"
    } else {
      e$source_table
    }

    if (src_tbl == "ind_meta") {
      src_df <- subset_df
    } else {
      if (!src_tbl %in% DBI::dbListTables(pop$db_conn)) {
        stop("Effect '", e$effect_name, "': source table '", src_tbl,
             "' does not exist.", call. = FALSE)
      }
      src_df <- DBI::dbGetQuery(
        pop$db_conn,
        paste0("SELECT * FROM ", src_tbl, " WHERE id_ind IN (", ids_sql, ")")
      )
      if (any(duplicated(src_df$id_ind))) {
        stop("Effect '", e$effect_name, "': source table '", src_tbl,
             "' has duplicate id_ind rows for the current subset.", call. = FALSE)
      }
      src_df <- src_df[match(ids_t, src_df$id_ind), , drop = FALSE]
    }

    if (!e$source_column %in% names(src_df)) {
      stop("Effect '", e$effect_name, "': column '", e$source_column,
           "' not found in table '", src_tbl, "'.", call. = FALSE)
    }

    group <- src_df[[e$source_column]]
    ec    <- e$effect_class

    if (ec %in% c("fixed", "fixed_class")) {
      # Resolve null_class_action (default to "skip" for old rows without column)
      nca <- if (!is.null(e$null_class_action) && !is.na(e$null_class_action)) {
        e$null_class_action
      } else {
        "skip"
      }

      null_mask <- is.na(group) | (as.character(group) == "NA") |
                   (as.character(group) == "")
      if (any(null_mask)) {
        if (nca == "error") {
          stop("Effect '", e$effect_name, "': ", sum(null_mask),
               " individual(s) have NULL in '", e$source_column, "'. ",
               "Set null_class_action = 'skip' or 'zero', or fix the data.",
               call. = FALSE)
        } else if (nca == "skip") {
          skip_mask <- skip_mask | null_mask
          group[null_mask] <- NA_character_
        } else {
          # "zero": treat NULL as reference level (zero contribution)
          group[null_mask] <- NA_character_
        }
      }

      level_map <- decode_levels_json(e$levels_json)
      group_chr <- as.character(group)
      shifts    <- unname(level_map[group_chr])

      # For null_class_action = "zero", NAs get 0
      if (nca == "zero") shifts[is.na(shifts) & null_mask] <- 0

      missing_lvls <- unique(group_chr[!is.na(group_chr) & is.na(shifts)])
      if (length(missing_lvls) > 0) {
        stop("Effect '", e$effect_name, "': levels have no shift defined: ",
             paste(missing_lvls, collapse = ", "),
             ". Update define_effect_fixed_class() to include all levels.",
             call. = FALSE)
      }
      shifts[is.na(shifts)] <- 0
      total <- total + shifts

    } else if (ec == "fixed_cov") {
      vals       <- as.numeric(src_df[[e$source_column]])
      center_val <- if (is.na(e$center)) 0 else e$center
      slope_val  <- if (is.na(e$slope))  0 else e$slope
      p_ord      <- if (!is.null(e$poly_order) && !is.na(e$poly_order) &&
                        e$poly_order >= 1L)
                      as.integer(e$poly_order) else 1L
      total <- total + slope_val * (vals - center_val)^p_ord

    } else {
      # Random effect — store/reuse draws in trait_random_effects
      unique_lvls <- unique(as.character(group[!is.na(group)]))

      existing_draws <- DBI::dbGetQuery(
        pop$db_conn,
        paste0("SELECT level, draw_value FROM trait_random_effects ",
               "WHERE phenotype_name = '", pn_safe,
               "' AND effect_name = '", gsub("'", "''", e$effect_name), "'")
      )
      existing_map <- if (nrow(existing_draws) > 0) {
        stats::setNames(existing_draws$draw_value, existing_draws$level)
      } else {
        stats::setNames(numeric(0), character(0))
      }

      new_lvls <- setdiff(unique_lvls, names(existing_map))
      if (length(new_lvls) > 0) {
        effect_var <- get_phenotype_var(pop, e$effect_name, phenotype_name)
        if (is.na(effect_var)) {
          stop("No variance stored for random effect '", e$effect_name,
               "' / phenotype '", phenotype_name, "'.", call. = FALSE)
        }
        new_draws <- switch(
          e$distribution %||% "normal",
          normal  = stats::rnorm(length(new_lvls), sd = sqrt(effect_var)),
          gamma   = stats::rgamma(length(new_lvls), shape = 1,
                                  rate = 1 / sqrt(effect_var)),
          uniform = stats::runif(length(new_lvls),
                                 min = -sqrt(3 * effect_var),
                                 max =  sqrt(3 * effect_var)),
          stats::rnorm(length(new_lvls), sd = sqrt(effect_var))
        )
        new_df <- tibble::tibble(
          phenotype_name = phenotype_name,
          effect_name    = e$effect_name,
          level          = new_lvls,
          draw_value     = new_draws,
          date_sampled   = as.Date(Sys.Date())
        )
        DBI::dbWriteTable(pop$db_conn, "trait_random_effects", new_df,
                          append = TRUE)
        existing_map <- c(existing_map, stats::setNames(new_draws, new_lvls))
      }

      per_ind <- unname(existing_map[as.character(group)])
      per_ind[is.na(per_ind)] <- 0
      total <- total + per_ind
    }
  }

  # Apply skip_mask: set contribution to NA for skipped individuals
  total[skip_mask] <- NA_real_
  list(contribution = total, n_skipped = sum(skip_mask))
}


#' Look up residual covariance information from phenotype_var_comp
#'
#' @param pop A `tidybreed_pop` object.
#' @param phenotype_names Character vector of phenotype names.
#' @param subset_df Data frame with id_ind column (used to look up condition
#'   group values when heterogeneous residuals are defined).
#' @return A list:
#'   - `R_unconditional`: numeric matrix (the unconditional residual R) or NULL.
#'   - `condition_column`: character scalar or NULL.
#'   - `condition_table`: character scalar or NULL.
#'   - `R_by_level`: named list of matrices keyed by condition_level, or NULL.
#'   - `residual_var_unconditional`: named numeric vector of diagonal variances.
#' @keywords internal
get_residual_cov <- function(pop, phenotype_names, subset_df = NULL) {
  n <- length(phenotype_names)
  pn_in <- paste0("'", gsub("'", "''", phenotype_names), "'", collapse = ", ")

  if (!"phenotype_var_comp" %in% DBI::dbListTables(pop$db_conn)) {
    return(list(R_unconditional = NULL, condition_column = NULL,
                condition_table = NULL, R_by_level = NULL,
                residual_var_unconditional = stats::setNames(rep(NA_real_, n),
                                                             phenotype_names)))
  }

  all_rows <- DBI::dbGetQuery(pop$db_conn, paste0(
    "SELECT * FROM phenotype_var_comp ",
    "WHERE effect_name = 'residual' ",
    "  AND phenotype_name_1 IN (", pn_in, ") ",
    "  AND phenotype_name_2 IN (", pn_in, ")"
  ))

  if (nrow(all_rows) == 0) {
    return(list(R_unconditional = NULL, condition_column = NULL,
                condition_table = NULL, R_by_level = NULL,
                residual_var_unconditional = stats::setNames(rep(NA_real_, n),
                                                             phenotype_names)))
  }

  # Separate unconditional rows from conditional rows
  uncond_rows <- all_rows[is.na(all_rows$condition_column) |
                          all_rows$condition_column == "", , drop = FALSE]
  cond_rows   <- all_rows[!is.na(all_rows$condition_column) &
                          all_rows$condition_column != "", , drop = FALSE]

  # Build unconditional R matrix
  R_unconditional <- NULL
  resid_var       <- stats::setNames(rep(NA_real_, n), phenotype_names)

  if (nrow(uncond_rows) > 0) {
    R_u <- matrix(NA_real_, n, n, dimnames = list(phenotype_names, phenotype_names))
    for (k in seq_len(nrow(uncond_rows))) {
      r <- uncond_rows[k, ]
      if (r$phenotype_name_1 %in% phenotype_names &&
          r$phenotype_name_2 %in% phenotype_names) {
        R_u[r$phenotype_name_1, r$phenotype_name_2] <- r$cov_value
      }
    }
    # Extract diagonal variances for phenotypes that have them (even if others don't)
    diag_vals <- diag(R_u)
    resid_var[!is.na(diag_vals)] <- diag_vals[!is.na(diag_vals)]
    # Full matrix only usable for joint MVN draw when entirely non-NA
    if (!any(is.na(R_u))) R_unconditional <- R_u
  }

  # Conditional rows
  R_by_level      <- NULL
  condition_column <- NULL
  condition_table  <- NULL

  if (nrow(cond_rows) > 0) {
    # Use the first condition_column found (assuming one per phenotype set)
    condition_column <- cond_rows$condition_column[1]
    condition_table  <- cond_rows$condition_table[1]
    if (is.na(condition_table) || !nzchar(condition_table)) {
      condition_table <- "ind_meta"
    }

    levels_found <- unique(cond_rows$condition_level[!is.na(cond_rows$condition_level)])
    R_by_level <- list()
    for (lvl in levels_found) {
      lvl_rows <- cond_rows[!is.na(cond_rows$condition_level) &
                            cond_rows$condition_level == lvl, , drop = FALSE]
      R_l <- matrix(NA_real_, n, n,
                    dimnames = list(phenotype_names, phenotype_names))
      for (k in seq_len(nrow(lvl_rows))) {
        r <- lvl_rows[k, ]
        if (r$phenotype_name_1 %in% phenotype_names &&
            r$phenotype_name_2 %in% phenotype_names) {
          R_l[r$phenotype_name_1, r$phenotype_name_2] <- r$cov_value
        }
      }
      if (!any(is.na(R_l))) R_by_level[[lvl]] <- R_l
    }
    if (length(R_by_level) == 0) R_by_level <- NULL
  }

  list(
    R_unconditional              = R_unconditional,
    condition_column             = condition_column,
    condition_table              = condition_table,
    R_by_level                   = R_by_level,
    residual_var_unconditional   = resid_var
  )
}


#' Null-coalescing operator for internal use
#' @keywords internal
`%||%` <- function(x, y) if (is.null(x) || (length(x) == 1 && is.na(x))) y else x


#' Compute the next pheno_number for each individual for a given phenotype
#'
#' @param pop A `tidybreed_pop` object.
#' @param phenotype_name Character. Phenotype name.
#' @param ids Character vector of individual IDs.
#' @return Integer vector, same length and order as `ids`.
#' @keywords internal
next_pheno_numbers <- function(pop, phenotype_name, ids) {
  if (length(ids) == 0) return(integer(0))
  ids_df <- data.frame(id_ind = ids, stringsAsFactors = FALSE)
  DBI::dbWriteTable(pop$db_conn, "_tb_pn_ids", ids_df, overwrite = TRUE)
  pn_safe <- gsub("'", "''", phenotype_name)
  res <- DBI::dbGetQuery(pop$db_conn, paste0(
    "SELECT f.id_ind, COALESCE(MAX(p.pheno_number), 0) + 1 AS next_pn ",
    "FROM _tb_pn_ids AS f ",
    "LEFT JOIN ind_phenotype AS p ",
    "  ON p.id_ind = f.id_ind AND p.phenotype_name = '", pn_safe, "' ",
    "GROUP BY f.id_ind"
  ))
  DBI::dbExecute(pop$db_conn, "DROP TABLE _tb_pn_ids")
  as.integer(stats::setNames(res$next_pn, res$id_ind)[ids])
}


#' Generate next global integer IDs for ind_phenotype
#'
#' @param pop A `tidybreed_pop` object.
#' @param n Number of IDs to generate.
#' @return Integer vector of length `n`.
#' @keywords internal
next_phenotype_ids <- function(pop, n) {
  if (n == 0L) return(integer(0))
  start <- next_int_id(pop$db_conn, "ind_phenotype", "id_phenotype")
  seq.int(start, start + n - 1L)
}
