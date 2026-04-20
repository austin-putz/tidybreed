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
#' @param ind_meta_df Data frame containing at least `id_ind` and any columns
#'   referenced by trait_effects rows for this trait.
#' @return Numeric vector of length `nrow(ind_meta_df)` giving the summed
#'   covariate contribution per individual.
#' @keywords internal
compute_covariate_contribution <- function(pop, trait, ind_meta_df) {

  effects <- DBI::dbGetQuery(
    pop$db_conn,
    paste0("SELECT * FROM trait_effects WHERE trait_name = '", trait, "'")
  )
  n_ind <- nrow(ind_meta_df)
  if (nrow(effects) == 0) return(rep(0, n_ind))

  total <- rep(0, n_ind)
  for (i in seq_len(nrow(effects))) {
    e <- effects[i, ]
    if (!e$source_column %in% names(ind_meta_df)) {
      warning("Covariate '", e$effect_name, "' references missing column '",
              e$source_column, "' in ind_meta; skipping.", call. = FALSE)
      next
    }
    group <- ind_meta_df[[e$source_column]]
    if (e$effect_class == "fixed") {
      level_map <- decode_levels_json(e$levels_json)
      shifts <- unname(level_map[as.character(group)])
      shifts[is.na(shifts)] <- 0
      total <- total + shifts
    } else {  # random
      unique_levels <- unique(group[!is.na(group)])
      draws <- switch(
        e$distribution,
        normal  = stats::rnorm(length(unique_levels),
                               sd = sqrt(e$variance)),
        gamma   = stats::rgamma(length(unique_levels),
                                shape = 1, rate = 1 / sqrt(e$variance)),
        uniform = stats::runif(length(unique_levels),
                               min = -sqrt(3 * e$variance),
                               max =  sqrt(3 * e$variance)),
        stats::rnorm(length(unique_levels), sd = sqrt(e$variance))
      )
      level_draw <- stats::setNames(draws, as.character(unique_levels))
      per_ind <- unname(level_draw[as.character(group)])
      per_ind[is.na(per_ind)] <- 0
      total <- total + per_ind
    }
  }
  total
}


#' Load the residual covariance matrix for a set of traits
#'
#' @param pop A `tidybreed_pop` object.
#' @param traits Character vector of trait names.
#' @return Numeric matrix with rows/columns named by `traits`, or `NULL` if
#'   not all pairs are present in `trait_residual_cov`.
#' @keywords internal
load_residual_cov <- function(pop, traits) {
  if (!"trait_residual_cov" %in% pop$tables) return(NULL)
  n <- length(traits)
  R <- matrix(NA_real_, nrow = n, ncol = n, dimnames = list(traits, traits))
  rows <- DBI::dbGetQuery(
    pop$db_conn,
    paste0("SELECT trait_1, trait_2, cov FROM trait_residual_cov ",
           "WHERE trait_1 IN (",
           paste0("'", traits, "'", collapse = ", "), ") ",
           "AND trait_2 IN (",
           paste0("'", traits, "'", collapse = ", "), ")")
  )
  if (nrow(rows) == 0) return(NULL)
  for (i in seq_len(nrow(rows))) {
    R[rows$trait_1[i], rows$trait_2[i]] <- rows$cov[i]
  }
  if (any(is.na(R))) return(NULL)
  R
}


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
