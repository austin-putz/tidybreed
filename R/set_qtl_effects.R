#' Sample or assign additive QTL effects for a trait
#'
#' @description
#' Writes the `add_{trait}` column in `genome_meta`, giving the additive
#' substitution effect for each locus (with `NA` at non-QTL loci). Two modes:
#'
#' * **Manual**: pass `effects`, a numeric vector of length equal to the
#'   number of QTL for this trait (loci where `is_QTL_{trait}` is `TRUE`), in
#'   ascending `locus_id` order.
#' * **Sampled**: draw effects from `distribution` (`"normal"` or
#'   `"gamma"`). If `scale_to_target = TRUE`, effects are then rescaled so
#'   that the realised `var(X %*% a)` across the currently genotyped
#'   individuals matches the trait's `target_add_var`.
#'
#' Requires [define_qtl()] has been called for `trait`.
#'
#' @param pop A `tidybreed_pop` object.
#' @param trait_name Character. Name of an existing trait.
#' @param effects Optional numeric vector of length `n_qtl` (manual mode).
#' @param distribution Character. `"normal"` (default) or `"gamma"`, used
#'   when `effects` is `NULL`.
#' @param scale_to_target Logical. If `TRUE`, rescale sampled effects so the
#'   realised additive variance on current genotyped individuals equals
#'   `trait_meta$target_add_var`.
#' @param seed Optional integer. Passed to [set.seed()] for reproducibility.
#'
#' @return The modified `tidybreed_pop` (invisibly).
#'
#' @seealso [set_qtl_effects_multi()] for correlated multi-trait effects.
#'
#' @examples
#' \dontrun{
#' pop <- pop |>
#'   add_trait("ADG", target_add_var = 0.25, residual_var = 0.75) |>
#'   define_qtl("ADG", n = 100) |>
#'   set_qtl_effects("ADG", distribution = "normal")
#' }
#' @export
set_qtl_effects <- function(pop,
                            trait_name,
                            effects         = NULL,
                            distribution    = c("normal", "gamma"),
                            scale_to_target = TRUE,
                            seed            = NULL) {

  stopifnot(inherits(pop, "tidybreed_pop"))
  validate_tidybreed_pop(pop)
  validate_sql_identifier(trait_name, what = "trait name")
  distribution <- match.arg(distribution)

  if (!is.null(seed)) set.seed(seed)

  qtl_col <- paste0("is_QTL_", trait_name)
  genome_cols <- DBI::dbListFields(pop$db_conn, "genome_meta")
  if (!qtl_col %in% genome_cols) {
    stop("QTL column '", qtl_col, "' not found. Call define_qtl('", trait_name,
         "', ...) first.", call. = FALSE)
  }

  trait_row <- DBI::dbGetQuery(
    pop$db_conn,
    paste0("SELECT target_add_var FROM trait_meta WHERE trait_name = '",
           trait_name, "'")
  )
  if (nrow(trait_row) == 0) {
    stop("Trait '", trait_name, "' not found in trait_meta.", call. = FALSE)
  }
  target_add_var <- trait_row$target_add_var

  genome <- dplyr::collect(get_table(pop, "genome_meta"))
  qtl_tf <- as.logical(genome[[qtl_col]])
  qtl_tf[is.na(qtl_tf)] <- FALSE
  n_qtl <- sum(qtl_tf)

  if (n_qtl == 0) {
    stop("No QTL defined for trait '", trait_name, "'.", call. = FALSE)
  }

  if (!is.null(effects)) {
    if (!is.numeric(effects)) {
      stop("`effects` must be numeric.", call. = FALSE)
    }
    if (length(effects) != n_qtl) {
      stop("`effects` length (", length(effects), ") must equal number of QTL (",
           n_qtl, ").", call. = FALSE)
    }
    qtl_effects <- as.numeric(effects)
  } else {
    qtl_effects <- switch(
      distribution,
      normal = stats::rnorm(n_qtl),
      gamma  = stats::rgamma(n_qtl, shape = 0.4, rate = 1.66) *
                 sample(c(-1, 1), n_qtl, replace = TRUE)
    )
    if (scale_to_target) {
      qtl_effects <- rescale_effects_to_target(
        pop, trait_name, qtl_tf, qtl_effects, target_add_var
      )
    }
  }

  full_effects <- rep(NA_real_, nrow(genome))
  full_effects[qtl_tf] <- qtl_effects

  args   <- setNames(list(full_effects), paste0("add_", trait_name))
  result <- do.call(mutate_table, c(list(tbl_obj = get_table(pop, "genome_meta")), args))

  message("Set additive effects for ", n_qtl, " QTL on trait '", trait_name, "'.")
  invisible(result)
}


#' Sample correlated additive QTL effects across multiple traits
#'
#' @description
#' Draws additive effects for multiple traits from a multivariate normal
#' distribution keyed by the additive-genetic covariance matrix `G`. Two
#' selection methods:
#'
#' * `method = "shared"` — effects are drawn only at loci that are QTL for
#'   **all** listed traits (pure pleiotropy). Loci that are QTL for only a
#'   subset of traits receive independent draws (with the diagonal variance
#'   of `G` for that trait).
#' * `method = "union"` — effects are drawn at the union of QTL across traits
#'   using `MASS::mvrnorm`; loci that are not QTL for a particular trait have
#'   their effect set to `NA` for that trait.
#'
#' If `scale_to_target = TRUE`, each trait's vector of effects is then
#'   rescaled independently so the realised `var(X %*% a_k)` matches its
#'   `target_add_var` — this slightly shifts the realised off-diagonal of
#'   `G` but matches conventional simulation practice.
#'
#' @param pop A `tidybreed_pop` object.
#' @param trait_names Character vector of trait names (length >= 2). All must
#'   exist in `trait_meta` and have `is_QTL_{trait_name}` already defined.
#' @param G Numeric matrix of additive-genetic (co)variances. Must be square
#'   with side length `length(trait_names)` and symmetric positive semi-definite.
#' @param method Character. `"shared"` (default) or `"union"`.
#' @param scale_to_target Logical. Rescale each trait's effects to its
#'   `target_add_var`.
#' @param seed Optional integer for reproducibility.
#'
#' @return The modified `tidybreed_pop` (invisibly).
#'
#' @examples
#' \dontrun{
#' G <- matrix(c(0.25, 0.10, 0.10, 0.30), 2, 2)
#' pop <- pop |>
#'   add_trait("ADG", target_add_var = 0.25, residual_var = 0.75) |>
#'   add_trait("BW",  target_add_var = 0.30, residual_var = 0.70) |>
#'   define_qtl("ADG", n = 200) |>
#'   define_qtl("BW",  n = 200) |>
#'   set_qtl_effects_multi(trait_names = c("ADG", "BW"), G = G)
#' }
#' @export
set_qtl_effects_multi <- function(pop,
                                  trait_names,
                                  G,
                                  method          = c("shared", "union"),
                                  scale_to_target = TRUE,
                                  seed            = NULL) {

  stopifnot(inherits(pop, "tidybreed_pop"))
  validate_tidybreed_pop(pop)
  stopifnot(is.character(trait_names), length(trait_names) >= 2)
  lapply(trait_names, validate_sql_identifier, what = "trait name")
  method <- match.arg(method)

  if (!requireNamespace("MASS", quietly = TRUE)) {
    stop("Package 'MASS' is required for set_qtl_effects_multi(). ",
         "Install with install.packages('MASS').", call. = FALSE)
  }

  if (!is.matrix(G) || nrow(G) != length(trait_names) || ncol(G) != length(trait_names)) {
    stop("`G` must be a square matrix with side = length(trait_names).",
         call. = FALSE)
  }
  if (!isSymmetric(unname(G))) {
    stop("`G` must be symmetric.", call. = FALSE)
  }

  if (!is.null(seed)) set.seed(seed)

  genome <- dplyr::collect(get_table(pop, "genome_meta"))
  n_loci <- nrow(genome)

  # Validate trait_names and QTL columns
  trait_meta_rows <- DBI::dbGetQuery(
    pop$db_conn,
    paste0("SELECT trait_name, target_add_var FROM trait_meta WHERE trait_name IN (",
           paste0("'", trait_names, "'", collapse = ", "), ")")
  )
  missing_traits <- setdiff(trait_names, trait_meta_rows$trait_name)
  if (length(missing_traits) > 0) {
    stop("Traits not found in trait_meta: ",
         paste(missing_traits, collapse = ", "), call. = FALSE)
  }
  target_var <- stats::setNames(
    trait_meta_rows$target_add_var[match(trait_names, trait_meta_rows$trait_name)],
    trait_names
  )

  qtl_tf_mat <- matrix(FALSE, nrow = n_loci, ncol = length(trait_names),
                       dimnames = list(NULL, trait_names))
  for (t in trait_names) {
    col <- paste0("is_QTL_", t)
    if (!col %in% colnames(genome)) {
      stop("QTL column '", col, "' not found. Call define_qtl('", t,
           "', ...) first.", call. = FALSE)
    }
    v <- as.logical(genome[[col]])
    v[is.na(v)] <- FALSE
    qtl_tf_mat[, t] <- v
  }

  effects_mat <- matrix(NA_real_, nrow = n_loci, ncol = length(trait_names),
                        dimnames = list(NULL, trait_names))

  if (method == "shared") {
    shared <- apply(qtl_tf_mat, 1, all)
    n_shared <- sum(shared)
    if (n_shared == 0) {
      warning("No loci are QTL for all traits; using method='union' fallback.",
              call. = FALSE)
    } else {
      draws <- MASS::mvrnorm(n = n_shared, mu = rep(0, length(trait_names)),
                             Sigma = G)
      if (is.null(dim(draws))) draws <- matrix(draws, nrow = 1)
      effects_mat[shared, ] <- draws
    }

    # Independent draws for trait-specific (non-shared) QTL
    for (k in seq_along(trait_names)) {
      t <- trait_names[k]
      solo <- qtl_tf_mat[, t] & !shared
      if (sum(solo) > 0) {
        effects_mat[solo, t] <- stats::rnorm(sum(solo), sd = sqrt(G[k, k]))
      }
    }
  } else {  # union
    any_qtl <- apply(qtl_tf_mat, 1, any)
    n_any <- sum(any_qtl)
    draws <- MASS::mvrnorm(n = n_any, mu = rep(0, length(trait_names)),
                           Sigma = G)
    if (is.null(dim(draws))) draws <- matrix(draws, nrow = 1)
    # Zero out where a locus is not QTL for that trait
    block_mask <- qtl_tf_mat[any_qtl, , drop = FALSE]
    draws[!block_mask] <- NA_real_
    effects_mat[any_qtl, ] <- draws
  }

  # Rescale per-trait to hit target_add_var (treating NAs as 0 for the calc)
  if (scale_to_target) {
    geno_mat <- get_genotype_matrix(pop)
    for (k in seq_along(trait_names)) {
      t <- trait_names[k]
      a <- effects_mat[, t]
      a_filled <- ifelse(is.na(a), 0, a)
      bv <- as.numeric(geno_mat %*% a_filled)
      realised <- stats::var(bv)
      if (realised > 0 && is.finite(realised)) {
        scale <- sqrt(target_var[[t]] / realised)
        effects_mat[, t] <- a * scale
      }
    }
  }

  # Write one column per trait
  for (t in trait_names) {
    args <- setNames(list(effects_mat[, t]), paste0("add_", t))
    pop  <- do.call(mutate_table, c(list(tbl_obj = get_table(pop, "genome_meta")), args))
  }

  message("Set correlated additive effects for traits: ",
          paste(trait_names, collapse = ", "), " (method: ", method, ")")
  invisible(pop)
}


#' Rescale an effect vector to hit a target additive variance
#'
#' @param pop A `tidybreed_pop` object.
#' @param trait Trait name (unused except for messaging).
#' @param qtl_tf Logical mask of QTL loci, length `n_loci`.
#' @param qtl_effects Numeric effects at QTL loci, length `sum(qtl_tf)`.
#' @param target_add_var Target additive variance.
#' @return Rescaled `qtl_effects` vector.
#' @keywords internal
rescale_effects_to_target <- function(pop, trait, qtl_tf, qtl_effects,
                                      target_add_var) {

  geno_mat <- get_genotype_matrix(pop)
  full_a <- rep(0, ncol(geno_mat))
  full_a[qtl_tf] <- qtl_effects
  bv <- as.numeric(geno_mat %*% full_a)
  realised <- stats::var(bv)
  if (realised <= 0 || !is.finite(realised)) {
    warning("Realised additive variance is zero; cannot rescale. ",
            "Effects returned unchanged.", call. = FALSE)
    return(qtl_effects)
  }
  scale <- sqrt(target_add_var / realised)
  qtl_effects * scale
}


#' Pull the genotype matrix (individuals x loci) into memory
#'
#' @param pop A `tidybreed_pop` object.
#' @param subset_ids Optional character vector of `id_ind`. When supplied,
#'   only those rows are pulled from DuckDB (push-down filter) — important
#'   for large genomes where the full genotype matrix wouldn't fit in R.
#' @return Numeric matrix with rows in `id_ind` order and columns in
#'   `locus_id` order.
#' @keywords internal
get_genotype_matrix <- function(pop, subset_ids = NULL) {
  if (!"genome_genotype" %in% pop$tables) {
    stop("genome_genotype table does not exist. Cannot rescale effects.",
         call. = FALSE)
  }
  if (is.null(subset_ids)) {
    df <- DBI::dbGetQuery(pop$db_conn, "SELECT * FROM genome_genotype")
  } else {
    if (length(subset_ids) == 0) {
      stop("subset_ids is empty.", call. = FALSE)
    }
    ids_sql <- paste0("'", subset_ids, "'", collapse = ", ")
    df <- DBI::dbGetQuery(
      pop$db_conn,
      paste0("SELECT * FROM genome_genotype WHERE id_ind IN (", ids_sql, ")")
    )
  }
  if (nrow(df) == 0) {
    stop("No genotype rows available for the requested subset.",
         call. = FALSE)
  }
  locus_cols <- setdiff(names(df), "id_ind")
  locus_order <- order(as.integer(sub("^locus_", "", locus_cols)))
  mat <- as.matrix(df[, locus_cols[locus_order], drop = FALSE])
  rownames(mat) <- df$id_ind
  mat
}
