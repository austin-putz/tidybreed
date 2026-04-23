#' Sample or assign additive QTL effects for a trait
#'
#' @description
#' Writes the `add_{trait}` and `base_allele_freq_{trait}` columns in
#' `genome_meta`. Two modes:
#'
#' * **Manual**: pass `effects`, a numeric vector of length equal to the
#'   number of QTL for this trait (loci where `is_QTL_{trait}` is `TRUE`), in
#'   ascending `locus_id` order.
#' * **Sampled**: draw effects from `distribution` (`"normal"` or
#'   `"gamma"`). If `scale_to_target = TRUE`, effects are rescaled using the
#'   Falconer formula so that the expected additive variance in the base
#'   population equals `trait_meta$target_add_var`.
#'
#' The `base` argument controls which allele frequencies are used for effect
#' scaling and TBV centering:
#'
#' * `"founder_haplotypes"` (default) — uses `founder_allele_freq` from
#'   `genome_meta` (requires `initialize_genome()` was called with
#'   `n_haplotypes`).
#' * `"current_pop"` — computes allele frequencies from the current
#'   `genome_haplotype` table. Pass a filtered `tidybreed_table` as the first
#'   argument to restrict which individuals define the base population.
#'
#' Requires [define_qtl()] has been called for `trait`.
#'
#' @param x A `tidybreed_pop` object, or a `tidybreed_table` (from
#'   [get_table()] plus optional [filter()]) when `base = "current_pop"` to
#'   specify which individuals define the base allele frequencies.
#' @param trait_name Character. Name of an existing trait.
#' @param effects Optional numeric vector of length `n_qtl` (manual mode).
#' @param distribution Character. `"normal"` (default) or `"gamma"`, used
#'   when `effects` is `NULL`.
#' @param base Character. `"founder_haplotypes"` (default) or `"current_pop"`.
#'   Determines which allele frequencies are used for Falconer variance scaling
#'   and TBV centering.
#' @param scale_to_target Logical. If `TRUE`, rescale sampled effects using the
#'   Falconer formula so the expected additive variance equals
#'   `trait_meta$target_add_var`.
#' @param seed Optional integer. Passed to [set.seed()] for reproducibility.
#'
#' @return The modified `tidybreed_pop` (invisibly).
#'
#' @seealso [set_qtl_effects_multi()] for correlated multi-trait effects.
#'
#' @examples
#' \dontrun{
#' # Default: scale using founder haplotype allele frequencies
#' pop <- pop |>
#'   add_trait("ADG", target_add_var = 0.25, residual_var = 0.75) |>
#'   define_qtl("ADG", n = 100) |>
#'   set_qtl_effects("ADG", distribution = "normal")
#'
#' # current_pop: use only generation-0 animals to define base frequencies
#' pop <- pop |>
#'   get_table("ind_meta") |>
#'   dplyr::filter(gen == 0L) |>
#'   set_qtl_effects("ADG", base = "current_pop", distribution = "normal")
#' }
#' @export
set_qtl_effects <- function(x,
                            trait_name,
                            effects         = NULL,
                            distribution    = c("normal", "gamma"),
                            base            = c("founder_haplotypes", "current_pop"),
                            scale_to_target = TRUE,
                            seed            = NULL) {

  if (inherits(x, "tidybreed_table")) {
    pop      <- x$pop
    base_ids <- unique(dplyr::collect(x$tbl)[["id_ind"]])
  } else if (inherits(x, "tidybreed_pop")) {
    pop      <- x
    base_ids <- NULL
  } else {
    stop("x must be a tidybreed_pop or tidybreed_table.", call. = FALSE)
  }

  validate_tidybreed_pop(pop)
  validate_sql_identifier(trait_name, what = "trait name")
  distribution <- match.arg(distribution)
  base         <- match.arg(base)

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

  # Compute base allele frequencies (always needed for TBV centering)
  p_base <- compute_base_allele_freq(pop, base, base_ids)

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
      qtl_effects <- rescale_effects_to_target(qtl_tf, qtl_effects, target_add_var, p_base)
    }
  }

  full_effects <- rep(NA_real_, nrow(genome))
  full_effects[qtl_tf] <- qtl_effects

  args <- c(
    setNames(list(full_effects), paste0("add_", trait_name)),
    setNames(list(p_base),       paste0("base_allele_freq_", trait_name))
  )
  result <- do.call(mutate_table, c(list(tbl_obj = get_table(pop, "genome_meta")), args))

  message("Set additive effects for ", n_qtl, " QTL on trait '", trait_name,
          "' (base: ", base, ").")
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
#' If `scale_to_target = TRUE`, each trait's effects are rescaled independently
#' using the Falconer formula so the expected additive variance in the base
#' population matches its `target_add_var`. See [set_qtl_effects()] for details
#' on the `base` argument.
#'
#' @param x A `tidybreed_pop` object, or a `tidybreed_table` when
#'   `base = "current_pop"` to specify which individuals define the base allele
#'   frequencies.
#' @param trait_names Character vector of trait names (length >= 2). All must
#'   exist in `trait_meta` and have `is_QTL_{trait_name}` already defined.
#' @param G Numeric matrix of additive-genetic (co)variances. Must be square
#'   with side length `length(trait_names)` and symmetric positive semi-definite.
#' @param method Character. `"shared"` (default) or `"union"`.
#' @param base Character. `"founder_haplotypes"` (default) or `"current_pop"`.
#' @param scale_to_target Logical. Rescale each trait's effects to its
#'   `target_add_var` using the Falconer formula.
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
set_qtl_effects_multi <- function(x,
                                  trait_names,
                                  G,
                                  method          = c("shared", "union"),
                                  base            = c("founder_haplotypes", "current_pop"),
                                  scale_to_target = TRUE,
                                  seed            = NULL) {

  if (inherits(x, "tidybreed_table")) {
    pop      <- x$pop
    base_ids <- unique(dplyr::collect(x$tbl)[["id_ind"]])
  } else if (inherits(x, "tidybreed_pop")) {
    pop      <- x
    base_ids <- NULL
  } else {
    stop("x must be a tidybreed_pop or tidybreed_table.", call. = FALSE)
  }

  validate_tidybreed_pop(pop)
  base <- match.arg(base)
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

  # Compute base allele frequencies (always needed for TBV centering)
  p_base <- compute_base_allele_freq(pop, base, base_ids)

  # Rescale per-trait using Falconer formula
  if (scale_to_target) {
    for (k in seq_along(trait_names)) {
      t      <- trait_names[k]
      qtl_tf_k <- qtl_tf_mat[, t]
      a_qtl  <- effects_mat[qtl_tf_k, t]
      if (any(!is.na(a_qtl))) {
        effects_mat[qtl_tf_k, t] <- rescale_effects_to_target(
          qtl_tf_k, a_qtl, target_var[[t]], p_base
        )
      }
    }
  }

  # Write one column per trait plus base allele frequencies
  for (t in trait_names) {
    args <- c(
      setNames(list(effects_mat[, t]), paste0("add_", t)),
      setNames(list(p_base),           paste0("base_allele_freq_", t))
    )
    pop <- do.call(mutate_table, c(list(tbl_obj = get_table(pop, "genome_meta")), args))
  }

  message("Set correlated additive effects for traits: ",
          paste(trait_names, collapse = ", "), " (method: ", method, ")")
  invisible(pop)
}


#' Rescale QTL effects to hit a target additive variance using Falconer formula
#'
#' @param qtl_tf Logical mask of QTL loci, length `n_loci`.
#' @param qtl_effects Numeric effects at QTL loci, length `sum(qtl_tf)`.
#' @param target_add_var Target additive variance.
#' @param p_base Numeric vector of base allele frequencies, length `n_loci`.
#' @return Rescaled `qtl_effects` vector.
#' @keywords internal
rescale_effects_to_target <- function(qtl_tf, qtl_effects, target_add_var, p_base) {
  p_qtl <- p_base[qtl_tf]
  V_A   <- sum(2 * p_qtl * (1 - p_qtl) * qtl_effects^2)
  if (V_A <= 0 || !is.finite(V_A)) {
    warning("Falconer V_A is zero or infinite; effects returned unchanged.", call. = FALSE)
    return(qtl_effects)
  }
  qtl_effects * sqrt(target_add_var / V_A)
}


#' Compute per-locus allele frequencies from the base population
#'
#' @param pop A `tidybreed_pop` object.
#' @param base Character. `"founder_haplotypes"` or `"current_pop"`.
#' @param base_ids Optional character vector of `id_ind` for `"current_pop"`.
#' @return Numeric vector of allele frequencies, length `n_loci`, in
#'   `locus_id` order.
#' @keywords internal
compute_base_allele_freq <- function(pop, base, base_ids = NULL) {
  if (base == "founder_haplotypes") {
    if (!"founder_haplotypes" %in% pop$tables) {
      stop(
        "founder_haplotypes table not found. ",
        "Did you call initialize_genome() with n_haplotypes? ",
        "Use base = 'current_pop' instead.",
        call. = FALSE
      )
    }
    df <- DBI::dbGetQuery(pop$db_conn, "SELECT * FROM founder_haplotypes")
    if (nrow(df) == 0) {
      stop("founder_haplotypes table is empty.", call. = FALSE)
    }
    locus_cols  <- setdiff(names(df), "hap_id")
    locus_order <- order(as.integer(sub("^locus_", "", locus_cols)))
    hap_mat     <- as.matrix(df[, locus_cols[locus_order], drop = FALSE])
    return(colMeans(hap_mat))
  }

  # current_pop: compute from genome_haplotype rows
  if (!"genome_haplotype" %in% pop$tables) {
    stop("genome_haplotype table does not exist.", call. = FALSE)
  }
  if (is.null(base_ids)) {
    df <- DBI::dbGetQuery(pop$db_conn, "SELECT * FROM genome_haplotype")
  } else {
    if (length(base_ids) == 0) stop("base_ids is empty.", call. = FALSE)
    ids_sql <- paste0("'", base_ids, "'", collapse = ", ")
    df <- DBI::dbGetQuery(
      pop$db_conn,
      paste0("SELECT * FROM genome_haplotype WHERE id_ind IN (", ids_sql, ")")
    )
  }
  if (nrow(df) == 0) {
    stop("No haplotype rows found for the base population.", call. = FALSE)
  }
  locus_cols  <- setdiff(names(df), c("id_ind", "parent_origin"))
  locus_order <- order(as.integer(sub("^locus_", "", locus_cols)))
  hap_mat     <- as.matrix(df[, locus_cols[locus_order], drop = FALSE])
  colMeans(hap_mat)
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
