#' Encode a named numeric vector as a small JSON string
#'
#' @description
#' Serialises a named numeric vector (e.g. `c(M = 30, F = 0)`) into a
#' JSON-style string `{"M":30,"F":0}` that fits a VARCHAR column. Used by
#' [add_trait_covariate()] to persist fixed-effect level maps.
#'
#' Keys containing `"`, `,`, or `:` are not supported (fixed-effect levels in
#' breeding programs are typically short tokens like `"M"`, `"A"`, `"gen1"`).
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
  paste0("{",
         paste0("\"", keys, "\":", vals, collapse = ","),
         "}")
}


#' Decode a level-map JSON string produced by `encode_levels_json()`
#'
#' @param s Character scalar (one of the encoded strings, or `NA`).
#' @return Named numeric vector (possibly length 0), or `NULL` if input is
#'   missing / empty.
#' @keywords internal
decode_levels_json <- function(s) {
  if (is.na(s) || !nzchar(s)) return(NULL)
  inner <- sub("^\\{", "", sub("\\}$", "", s))
  if (!nzchar(inner)) return(stats::setNames(numeric(0), character(0)))
  pairs <- strsplit(inner, ",", fixed = TRUE)[[1]]
  split_on_colon <- strsplit(pairs, ":", fixed = TRUE)
  keys <- vapply(split_on_colon,
                 function(p) gsub("^\"|\"$", "", p[1]),
                 character(1))
  vals <- vapply(split_on_colon,
                 function(p) as.numeric(p[2]),
                 numeric(1))
  stats::setNames(vals, keys)
}


#' Sample residuals for one or more traits
#'
#' @description
#' Draws a length-`n_ind` residual vector per trait. When `R` is supplied,
#' draws from `MVN(0, R)`; otherwise draws independently from
#' `N(0, residual_var[k])` for each trait.
#'
#' @param n_ind Integer, number of individuals.
#' @param residual_var Numeric vector of per-trait residual variances, named
#'   by trait.
#' @param R Optional covariance matrix with row/col names matching
#'   `names(residual_var)`.
#' @return Numeric matrix of shape `n_ind x length(residual_var)`, columns
#'   named by trait.
#' @keywords internal
sample_residuals <- function(n_ind, residual_var, R = NULL) {
  traits <- names(residual_var)
  if (is.null(R)) {
    out <- vapply(
      traits,
      function(t) stats::rnorm(n_ind, sd = sqrt(residual_var[[t]])),
      numeric(n_ind)
    )
    if (!is.matrix(out)) out <- matrix(out, nrow = n_ind, ncol = length(traits))
    colnames(out) <- traits
    return(out)
  }
  if (!requireNamespace("MASS", quietly = TRUE)) {
    stop("Package 'MASS' is required for joint residual sampling.",
         call. = FALSE)
  }
  R <- R[traits, traits, drop = FALSE]
  draws <- MASS::mvrnorm(n = n_ind, mu = rep(0, length(traits)), Sigma = R)
  if (is.null(dim(draws))) draws <- matrix(draws, nrow = 1)
  colnames(draws) <- traits
  draws
}


#' Convert a liability vector to binary 0/1 using a prevalence threshold
#'
#' @param liability Numeric vector on the continuous scale.
#' @param prevalence Numeric between 0 and 1.
#' @return Integer vector of 0/1.
#' @keywords internal
liability_to_binary <- function(liability, prevalence) {
  # Threshold at the empirical (1 - prevalence) quantile so the realised rate
  # matches the target even when variance components differ from 1.
  thresh <- stats::quantile(liability, probs = 1 - prevalence, na.rm = TRUE,
                            names = FALSE)
  as.integer(liability > thresh)
}


#' Convert a liability vector to ordered integer categories
#'
#' @param liability Numeric vector on the continuous scale.
#' @param thresholds Numeric vector of cutpoints (length = n_categories - 1),
#'   in ascending order.
#' @return Integer vector: 1 for liability <= thresholds[1], 2 for next,
#'   ..., `length(thresholds) + 1` for above all.
#' @keywords internal
liability_to_categorical <- function(liability, thresholds) {
  thresholds <- sort(thresholds)
  findInterval(liability, thresholds) + 1L
}


#' Clip a count vector to [min_value, max_value] (NA = no bound)
#'
#' @param x Numeric vector.
#' @param min_value,max_value Numeric scalars or `NA`.
#' @return Integer vector (rounded from `x`) clipped to bounds.
#' @keywords internal
clip_count <- function(x, min_value = NA_real_, max_value = NA_real_) {
  y <- round(x)
  if (!is.na(min_value)) y <- pmax(y, min_value)
  if (!is.na(max_value)) y <- pmin(y, max_value)
  as.integer(y)
}


#' Compute the covariate contribution (fixed + random shifts) for each
#' individual in a subset, for a single trait
#'
#' @param pop A `tidybreed_pop` object (connection used for queries).
#' @param trait Trait name.
#' @param subset_df Data frame: the per-trait subset of `ind_meta` (already
#'   sex-filtered). Must contain `id_ind` and any `ind_meta` columns referenced
#'   by effects. Effects from other tables are fetched from the DB directly.
#' @return Numeric vector of length `nrow(subset_df)` giving the summed
#'   covariate contribution per individual.
#' @keywords internal
compute_covariate_contribution <- function(pop, trait, subset_df) {

  effects <- DBI::dbGetQuery(
    pop$db_conn,
    paste0("SELECT * FROM trait_effects WHERE trait_name = '", trait, "'")
  )
  n_ind <- nrow(subset_df)
  if (nrow(effects) == 0) return(rep(0, n_ind))

  ids_t   <- subset_df$id_ind
  ids_sql <- paste0("'", ids_t, "'", collapse = ", ")
  total   <- rep(0, n_ind)

  for (i in seq_len(nrow(effects))) {
    e <- effects[i, ]

    # Resolve source table (NA or empty means ind_meta for old rows)
    src_tbl <- if (is.na(e$source_table) || !nzchar(e$source_table)) {
      "ind_meta"
    } else {
      e$source_table
    }

    # Fetch source data frame (use subset_df for ind_meta; query DB otherwise)
    if (src_tbl == "ind_meta") {
      src_df <- subset_df
    } else {
      if (!src_tbl %in% DBI::dbListTables(pop$db_conn)) {
        stop("Effect '", e$effect_name, "': source table '", src_tbl,
             "' does not exist in the database.", call. = FALSE)
      }
      src_df <- DBI::dbGetQuery(
        pop$db_conn,
        paste0("SELECT * FROM ", src_tbl, " WHERE id_ind IN (", ids_sql, ")")
      )
      if (any(duplicated(src_df$id_ind))) {
        stop("Effect '", e$effect_name, "': source table '", src_tbl,
             "' has duplicate id_ind rows for the current subset. ",
             "Each individual must appear at most once.", call. = FALSE)
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
      # Discrete fixed effect — error on unknown level
      level_map <- decode_levels_json(e$levels_json)
      shifts <- unname(level_map[as.character(group)])
      missing_lvls <- unique(as.character(group[is.na(shifts)]))
      if (length(missing_lvls) > 0) {
        stop("Effect '", e$effect_name, "': the following levels in column '",
             e$source_column, "' have no shift defined: ",
             paste(missing_lvls, collapse = ", "),
             ". Update add_effect_fixed_class() to include all levels.",
             call. = FALSE)
      }
      total <- total + shifts

    } else if (ec == "fixed_cov") {
      # Continuous covariate
      vals       <- as.numeric(src_df[[e$source_column]])
      center_val <- if (is.na(e$center)) 0 else e$center
      slope_val  <- if (is.na(e$slope))  0 else e$slope
      total <- total + slope_val * (vals - center_val)

    } else {
      # Random effect — store/reuse draws in trait_random_effects
      unique_lvls <- unique(as.character(group[!is.na(group)]))

      existing_draws <- DBI::dbGetQuery(
        pop$db_conn,
        paste0("SELECT level, draw_value FROM trait_random_effects ",
               "WHERE trait_name = '", trait,
               "' AND effect_name = '", e$effect_name, "'")
      )
      existing_map <- if (nrow(existing_draws) > 0) {
        stats::setNames(existing_draws$draw_value, existing_draws$level)
      } else {
        stats::setNames(numeric(0), character(0))
      }

      new_lvls <- setdiff(unique_lvls, names(existing_map))
      if (length(new_lvls) > 0) {
        effect_var <- get_effect_var(pop, e$effect_name, trait)
        if (is.na(effect_var)) {
          stop("No variance stored for random effect '", e$effect_name,
               "' / trait '", trait, "'. ",
               "Call add_effect_cov_matrix(pop, '", e$effect_name,
               "', ...) or add_effect_random(..., variance = ...) first.",
               call. = FALSE)
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
          trait_name   = trait,
          effect_name  = e$effect_name,
          level        = new_lvls,
          draw_value   = new_draws,
          date_sampled = as.Date(Sys.Date())
        )
        DBI::dbWriteTable(pop$db_conn, "trait_random_effects", new_df,
                          append = TRUE)
        existing_map <- c(existing_map,
                          stats::setNames(new_draws, new_lvls))
      }

      per_ind <- unname(existing_map[as.character(group)])
      per_ind[is.na(per_ind)] <- 0
      total <- total + per_ind
    }
  }
  total
}


#' Null-coalescing operator for internal use
#' @keywords internal
`%||%` <- function(x, y) if (is.null(x) || (length(x) == 1 && is.na(x))) y else x


#' Generate per-trait next-record IDs for ind_phenotype
#'
#' @param pop A `tidybreed_pop` object.
#' @param trait Trait name.
#' @param n Number of IDs to generate.
#' @return Character vector of length `n` with unique ids like `"ADG-42"`.
#' @keywords internal
next_record_ids <- function(pop, trait, n) {
  if (n == 0) return(character(0))
  q <- paste0(
    "SELECT COALESCE(",
    "MAX(CAST(SUBSTRING(id_record FROM POSITION('-' IN id_record) + 1) AS INTEGER)),",
    " 0) AS mx FROM ind_phenotype WHERE id_record LIKE '", trait, "-%'"
  )
  start <- DBI::dbGetQuery(pop$db_conn, q)$mx + 1L
  paste0(trait, "-", seq.int(start, start + n - 1L))
}


#' Pull the haplotype for a single parent origin, as individuals x loci matrix
#'
#' @param pop A `tidybreed_pop` object.
#' @param parent_origin Integer 1 (paternal) or 2 (maternal).
#' @param id_ind Character vector of individuals to pull.
#' @return Numeric matrix, rows in the order of `id_ind`, columns in
#'   `locus_id` order.
#' @keywords internal
get_haplotype_matrix <- function(pop, parent_origin, id_ind) {
  ids_sql <- paste0("'", id_ind, "'", collapse = ", ")
  df <- DBI::dbGetQuery(
    pop$db_conn,
    paste0("SELECT * FROM genome_haplotype WHERE parent_origin = ",
           parent_origin, " AND id_ind IN (", ids_sql, ")")
  )
  if (nrow(df) == 0) {
    stop("No haplotype rows found for parent_origin = ", parent_origin,
         call. = FALSE)
  }
  locus_cols <- setdiff(names(df), c("id_ind", "parent_origin"))
  locus_order <- order(as.integer(sub("^locus_", "", locus_cols)))
  mat <- as.matrix(df[match(id_ind, df$id_ind), locus_cols[locus_order],
                      drop = FALSE])
  rownames(mat) <- id_ind
  mat
}
