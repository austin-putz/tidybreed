#' Define a variance-covariance matrix for any named effect
#'
#' @description
#' Single entry point for storing all variance and covariance data in tidybreed.
#' Routes to `trait_var_comp` for genetic effects and to `phenotype_var_comp` for
#' phenotype-level effects.
#'
#' Common `effect_name` values:
#'
#' * `"gen_add"` — additive genetic (co)variances (G matrix). Written to
#'   `trait_var_comp`. Used by [define_additive_effects()] when rescaling to
#'   target variance and as the sampling distribution for multi-trait draws.
#' * `"dominance"`, `"epistasis"` — future genetic effects. Written to
#'   `trait_var_comp`. Row/column names are trait names.
#' * `"residual"` — residual (co)variances (R matrix). Routed to
#'   `phenotype_var_comp` with `effect_name = "residual"`. Row/column names are
#'   phenotype names. Equivalent to calling [define_residual_cov()] with
#'   `condition_column = NULL`. Use this for a multi-phenotype correlated
#'   residual matrix; for a single scalar residual use `residual_var` in
#'   [define_phenotype()] instead.
#' * Any named random effect (`"hys"`, `"litter"`, `"pen"`, …) — written to
#'   `phenotype_var_comp`. Must match the `effect_name` used in
#'   [define_effect_random()]. Row/column names are phenotype names.
#'
#' `define_effect_cov_matrix()` can be called **before** [define_trait()] or
#' [define_effect_random()] — no prior setup is required.
#'
#' All n² pairs are stored. For `phenotype_var_comp` effects the names form a
#' *covariance block* that is declared in one call, as a complete matrix: a
#' call that names a fragment or a strict subset of an existing block is an
#' error, the matrix must be positive semi-definite, and a block cannot be
#' redefined once draws exist under it (the error gives the [remove_rows()]
#' call that clears them). In a block of two or more phenotypes every
#' [define_effect_random()] row for the effect must use
#' `distribution = "normal"` and read the same `(source_column, source_table)`.
#' See [define_residual_cov()] for the full rules. A rejected call changes
#' nothing.
#'
#' @param pop A `tidybreed_pop` object.
#' @param effect_name Character. Label for the variance component, e.g.
#'   `"gen_add"`, `"residual"`, `"hys"`.
#' @param cov_matrix A numeric square matrix. Must be symmetric within `tol`.
#'   Row and column names are used as trait/phenotype names when `trait_names`
#'   is not supplied.
#' @param trait_names Optional character vector of trait/phenotype names (length
#'   == `nrow(cov_matrix)`). Overrides the matrix's `rownames` / `colnames`.
#' @param tol Numeric. Tolerance for symmetry check (default `1e-9`).
#'
#' @return The modified `tidybreed_pop` (invisibly).
#'
#' @seealso [define_trait()], [define_effect_random()], [define_additive_effects()],
#'   [add_phenotype()]
#'
#' @examples
#' \dontrun{
#' # Additive genetic covariance matrix → trait_var_comp
#' G <- matrix(c(100, -20, -20, 50), 2, 2,
#'             dimnames = list(c("ADG", "BF"), c("ADG", "BF")))
#' pop <- pop |>
#'   define_effect_cov_matrix("gen_add", G)
#'
#' # Residual → phenotype_var_comp (effect_name = "residual")
#' R <- matrix(c(30, 5, 5, 10), 2, 2,
#'             dimnames = list(c("ADG", "BF"), c("ADG", "BF")))
#' pop <- pop |>
#'   define_effect_cov_matrix("residual", R)
#'
#' # Multi-phenotype HYS covariance → phenotype_var_comp (effect_name = "hys")
#' R_hys <- matrix(c(0.2, 0.05, 0.05, 0.3), 2, 2,
#'                 dimnames = list(c("ADG", "BF"), c("ADG", "BF")))
#' pop <- pop |>
#'   define_effect_cov_matrix("hys", R_hys)
#' }
#' @export
define_effect_cov_matrix <- function(pop,
                                   effect_name,
                                   cov_matrix,
                                   trait_names = NULL,
                                   tol         = 1e-9) {
  stopifnot(inherits(pop, "tidybreed_pop"))
  validate_tidybreed_pop(pop)
  validate_sql_identifier(effect_name, what = "effect name")

  if (!is.matrix(cov_matrix) || !is.numeric(cov_matrix)) {
    stop("`cov_matrix` must be a numeric matrix.", call. = FALSE)
  }
  n <- nrow(cov_matrix)
  if (ncol(cov_matrix) != n) {
    stop("`cov_matrix` must be square.", call. = FALSE)
  }

  if (!is.null(trait_names)) {
    if (length(trait_names) != n) {
      stop("`trait_names` length (", length(trait_names),
           ") must equal matrix dimension (", n, ").", call. = FALSE)
    }
  } else {
    trait_names <- rownames(cov_matrix)
    if (is.null(trait_names) || any(!nzchar(trait_names))) {
      stop("`cov_matrix` must have row names, or supply `trait_names`.",
           call. = FALSE)
    }
  }
  lapply(trait_names, validate_sql_identifier, what = "trait name")
  if (anyDuplicated(trait_names)) {
    stop("`trait_names` must not contain duplicates.", call. = FALSE)
  }
  dimnames(cov_matrix) <- list(trait_names, trait_names)

  if (!isSymmetric(unname(cov_matrix), tol = tol)) {
    stop("`cov_matrix` must be symmetric (max discrepancy: ",
         max(abs(cov_matrix - t(cov_matrix))), ").", call. = FALSE)
  }
  if (any(diag(cov_matrix) < 0)) {
    stop("Diagonal entries (variances) must be non-negative.", call. = FALSE)
  }

  # Genetic effects -> trait_var_comp; phenotype-level effects -> phenotype_var_comp
  genetic_effects <- c("gen_add", "dominance", "epistasis")

  if (effect_name %in% genetic_effects) {
    quoted_names <- paste0("'", trait_names, "'", collapse = ", ")
    DBI::dbExecute(
      pop$db_conn,
      paste0("DELETE FROM trait_var_comp WHERE effect_name = '", effect_name,
             "' AND trait_name_1 IN (", quoted_names,
             ") AND trait_name_2 IN (", quoted_names, ")")
    )
    # Build multi-row INSERT to avoid DBI::dbWriteTable() which consumes R's RNG
    n_rows    <- n * n
    start_id  <- next_int_id(pop$db_conn, "trait_var_comp", "id_trait_var_comp")
    value_rows <- character(n_rows)
    k <- 1L
    for (i in seq_len(n)) {
      for (j in seq_len(n)) {
        value_rows[k] <- paste0("(", start_id + k - 1L, ", '", effect_name, "', '",
                                trait_names[i], "', '", trait_names[j], "', ",
                                format(cov_matrix[i, j], scientific = FALSE), ")")
        k <- k + 1L
      }
    }
    DBI::dbExecute(
      pop$db_conn,
      paste0("INSERT INTO trait_var_comp ",
             "(id_trait_var_comp, effect_name, trait_name_1, trait_name_2, cov_value) VALUES ",
             paste(value_rows, collapse = ", "))
    )
  } else if (identical(effect_name, "residual")) {
    define_residual_cov(pop,
      phenotype_names  = trait_names,
      cov_matrix       = cov_matrix,
      condition_column = NULL)
  } else {
    # Named random effect: one complete, unconditional block (D1/D3/section 5.6)
    write_phenotype_cov_block(
      pop$db_conn, effect_name, trait_names, cov_matrix,
      caller = "define_effect_cov_matrix()", tol = tol)
  }

  message("Stored '", effect_name, "' covariance matrix for: ",
          paste(trait_names, collapse = ", "), ".")
  invisible(pop)
}


#' Get the variance (diagonal) for one trait from trait_var_comp
#'
#' @param pop A `tidybreed_pop` object.
#' @param effect_name Character.
#' @param trait_name Character.
#' @return Numeric scalar, or `NA_real_` if not found.
#' @keywords internal
get_trait_var <- function(pop, effect_name, trait_name) {
  row <- DBI::dbGetQuery(
    pop$db_conn,
    paste0("SELECT cov_value FROM trait_var_comp ",
           "WHERE effect_name = '", effect_name, "' ",
           "AND trait_name_1 = '", trait_name, "' ",
           "AND trait_name_2 = '", trait_name, "'")
  )
  if (nrow(row) == 0L) NA_real_ else row$cov_value[[1L]]
}


#' Get the variance (diagonal) for one phenotype from phenotype_var_comp
#'
#' @param pop A `tidybreed_pop` object.
#' @param effect_name Character.
#' @param phenotype_name Character.
#' @return Numeric scalar, or `NA_real_` if not found.
#' @keywords internal
get_phenotype_var <- function(pop, effect_name, phenotype_name) {
  eff_safe <- gsub("'", "''", effect_name)
  pn_safe  <- gsub("'", "''", phenotype_name)
  row <- DBI::dbGetQuery(
    pop$db_conn,
    paste0("SELECT cov_value FROM phenotype_var_comp ",
           "WHERE effect_name = '", eff_safe, "' ",
           "AND phenotype_name_1 = '", pn_safe, "' ",
           "AND phenotype_name_2 = '", pn_safe, "' ",
           "AND condition_column IS NULL")
  )
  if (nrow(row) == 0L) NA_real_ else row$cov_value[[1L]]
}


#' Load a full covariance matrix from trait_var_comp
#'
#' @param pop A `tidybreed_pop` object.
#' @param effect_name Character.
#' @param trait_names Character vector of trait names.
#' @return Named numeric matrix, or `NULL` if any entry is missing.
#' @keywords internal
load_trait_cov <- function(pop, effect_name, trait_names) {
  n <- length(trait_names)
  R <- matrix(NA_real_, nrow = n, ncol = n, dimnames = list(trait_names, trait_names))
  rows <- DBI::dbGetQuery(
    pop$db_conn,
    paste0("SELECT trait_name_1, trait_name_2, cov_value FROM trait_var_comp ",
           "WHERE effect_name = '", effect_name, "' ",
           "AND trait_name_1 IN (", paste0("'", trait_names, "'", collapse = ", "), ") ",
           "AND trait_name_2 IN (", paste0("'", trait_names, "'", collapse = ", "), ")")
  )
  if (nrow(rows) == 0L) return(NULL)
  for (i in seq_len(nrow(rows))) {
    R[rows$trait_name_1[i], rows$trait_name_2[i]] <- rows$cov_value[i]
  }
  if (any(is.na(R))) return(NULL)
  R
}


#' Load a full covariance matrix from phenotype_var_comp for a named random effect
#'
#' @param pop A `tidybreed_pop` object.
#' @param effect_name Character. The random effect name (not "residual").
#' @param phenotype_names Character vector of phenotype names.
#' @return Named numeric matrix, or `NULL` if any entry is missing.
#' @keywords internal
load_phenotype_cov <- function(pop, effect_name, phenotype_names) {
  n <- length(phenotype_names)
  R <- matrix(NA_real_, nrow = n, ncol = n,
              dimnames = list(phenotype_names, phenotype_names))
  eff_safe <- gsub("'", "''", effect_name)
  rows <- DBI::dbGetQuery(
    pop$db_conn,
    paste0("SELECT phenotype_name_1, phenotype_name_2, cov_value FROM phenotype_var_comp ",
           "WHERE effect_name = '", eff_safe, "' ",
           "AND phenotype_name_1 IN (", paste0("'", phenotype_names, "'", collapse = ", "), ") ",
           "AND phenotype_name_2 IN (", paste0("'", phenotype_names, "'", collapse = ", "), ")")
  )
  if (nrow(rows) == 0L) return(NULL)
  for (i in seq_len(nrow(rows))) {
    R[rows$phenotype_name_1[i], rows$phenotype_name_2[i]] <- rows$cov_value[i]
  }
  if (any(is.na(R))) return(NULL)
  R
}


#' Write a single diagonal entry to trait_var_comp
#'
#' Used internally by define_trait() to write a per-trait variance as a 1x1
#' diagonal entry. Uses dbExecute() to avoid consuming R's RNG.
#'
#' @param pop A `tidybreed_pop` object.
#' @param effect_name Character.
#' @param trait_name Character.
#' @param variance Numeric scalar.
#' @return The modified `tidybreed_pop` (invisibly).
#' @keywords internal
write_trait_var_diag <- function(pop, effect_name, trait_name, variance) {
  DBI::dbExecute(
    pop$db_conn,
    paste0("DELETE FROM trait_var_comp WHERE effect_name = '", effect_name,
           "' AND trait_name_1 = '", trait_name, "' AND trait_name_2 = '", trait_name, "'")
  )
  new_id <- next_int_id(pop$db_conn, "trait_var_comp", "id_trait_var_comp")
  DBI::dbExecute(
    pop$db_conn,
    paste0("INSERT INTO trait_var_comp ",
           "(id_trait_var_comp, effect_name, trait_name_1, trait_name_2, cov_value) VALUES (",
           new_id, ", '", effect_name, "', '", trait_name, "', '", trait_name, "', ",
           format(as.numeric(variance), scientific = FALSE), ")")
  )
  invisible(pop)
}
