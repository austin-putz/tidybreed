#' Store a variance-covariance matrix for any named effect
#'
#' @description
#' Writes a symmetric variance-covariance matrix into the unified
#' `trait_effect_cov` table under a user-supplied `effect_name`. This is the
#' single entry point for storing all variance and covariance data in
#' tidybreed.
#'
#' Common `effect_name` values:
#'
#' * `"gen_add"` — additive genetic (co)variances (G matrix). Used by
#'   [set_qtl_effects()] when rescaling to target variance and by
#'   [set_qtl_effects_multi()] as the sampling distribution.
#' * `"residual"` — residual (co)variances (R matrix). Used by [add_phenotype()]
#'   when drawing residuals.
#' * Any named random effect (`"litter"`, `"herd"`, `"dam"`, …) — must match
#'   the `effect_name` used in [add_effect_random()].
#'
#' `add_effect_cov_matrix()` can be called **before** [add_trait()] or
#' [add_effect_random()] — no prior setup is required. Intended to be called
#' early in a simulation so that downstream functions can look up stored
#' variances automatically.
#'
#' All n² pairs (diagonal and both off-diagonal directions) are stored to
#' support lookups in either (trait_1, trait_2) order. Previous entries for
#' this `effect_name` × traits combination are replaced.
#'
#' @param pop A `tidybreed_pop` object.
#' @param effect_name Character. Label for the variance component, e.g.
#'   `"gen_add"`, `"residual"`, `"litter"`.
#' @param cov_matrix A numeric square matrix. Must be symmetric within `tol`.
#'   Row and column names are used as trait names when `trait_names` is not
#'   supplied.
#' @param trait_names Optional character vector of trait names (length ==
#'   `nrow(cov_matrix)`). Overrides the matrix's `rownames` / `colnames`.
#' @param tol Numeric. Tolerance for symmetry check (default `1e-9`).
#'
#' @return The modified `tidybreed_pop` (invisibly).
#'
#' @seealso [add_trait()], [add_effect_random()], [set_qtl_effects_multi()],
#'   [add_phenotype()]
#'
#' @examples
#' \dontrun{
#' # Define additive genetic covariance matrix early in simulation setup
#' G <- matrix(c(100, -20, -20, 50), 2, 2,
#'             dimnames = list(c("ADG", "BF"), c("ADG", "BF")))
#' pop <- pop |>
#'   add_effect_cov_matrix("gen_add", G)
#'
#' R <- matrix(c(30, 5, 5, 10), 2, 2,
#'             dimnames = list(c("ADG", "BF"), c("ADG", "BF")))
#' pop <- pop |>
#'   add_effect_cov_matrix("residual", R)
#'
#' # Single-trait variance (1x1 matrix)
#' pop <- pop |>
#'   add_effect_cov_matrix("gen_add",
#'     matrix(100, 1, 1, dimnames = list("ADG", "ADG")))
#' }
#' @export
add_effect_cov_matrix <- function(pop,
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
    lapply(trait_names, validate_sql_identifier, what = "trait name")
  } else {
    trait_names <- rownames(cov_matrix)
    if (is.null(trait_names) || any(!nzchar(trait_names))) {
      stop("`cov_matrix` must have row names, or supply `trait_names`.",
           call. = FALSE)
    }
    lapply(trait_names, validate_sql_identifier, what = "trait name")
  }

  if (!isSymmetric(unname(cov_matrix), tol = tol)) {
    stop("`cov_matrix` must be symmetric (max discrepancy: ",
         max(abs(cov_matrix - t(cov_matrix))), ").", call. = FALSE)
  }
  if (any(diag(cov_matrix) < 0)) {
    stop("Diagonal entries (variances) must be non-negative.", call. = FALSE)
  }

  pop <- ensure_effect_cov_table(pop)

  quoted_traits <- paste0("'", trait_names, "'", collapse = ", ")
  DBI::dbExecute(
    pop$db_conn,
    paste0("DELETE FROM trait_effect_cov WHERE effect_name = '", effect_name,
           "' AND trait_1 IN (", quoted_traits,
           ") AND trait_2 IN (", quoted_traits, ")")
  )

  # Build multi-row INSERT to avoid DBI::dbWriteTable() which consumes R's RNG
  value_rows <- character(n * n)
  k <- 1L
  for (i in seq_len(n)) {
    for (j in seq_len(n)) {
      value_rows[k] <- paste0("('", effect_name, "', '", trait_names[i], "', '",
                              trait_names[j], "', ",
                              format(cov_matrix[i, j], scientific = FALSE), ")")
      k <- k + 1L
    }
  }
  DBI::dbExecute(
    pop$db_conn,
    paste0("INSERT INTO trait_effect_cov (effect_name, trait_1, trait_2, cov) VALUES ",
           paste(value_rows, collapse = ", "))
  )

  message("Stored '", effect_name, "' covariance matrix for trait(s): ",
          paste(trait_names, collapse = ", "), ".")
  invisible(pop)
}


#' Ensure the trait_effect_cov table exists in the database
#'
#' @param pop A `tidybreed_pop` object.
#' @return The `tidybreed_pop` with `$tables` updated.
#' @keywords internal
ensure_effect_cov_table <- function(pop) {
  if (!"trait_effect_cov" %in% DBI::dbListTables(pop$db_conn)) {
    DBI::dbExecute(pop$db_conn, "
      CREATE TABLE trait_effect_cov (
        effect_name VARCHAR,
        trait_1     VARCHAR,
        trait_2     VARCHAR,
        cov         DOUBLE,
        PRIMARY KEY (effect_name, trait_1, trait_2)
      )
    ")
  }
  pop$tables <- unique(c(pop$tables, "trait_effect_cov"))
  pop
}


#' Get the variance (diagonal) for one trait from trait_effect_cov
#'
#' @param pop A `tidybreed_pop` object.
#' @param effect_name Character.
#' @param trait_name Character.
#' @return Numeric scalar, or `NA_real_` if not found.
#' @keywords internal
get_effect_var <- function(pop, effect_name, trait_name) {
  if (!"trait_effect_cov" %in% DBI::dbListTables(pop$db_conn)) {
    return(NA_real_)
  }
  row <- DBI::dbGetQuery(
    pop$db_conn,
    paste0("SELECT cov FROM trait_effect_cov ",
           "WHERE effect_name = '", effect_name, "' ",
           "AND trait_1 = '", trait_name, "' ",
           "AND trait_2 = '", trait_name, "'")
  )
  if (nrow(row) == 0L) NA_real_ else row$cov[[1L]]
}


#' Load a full covariance matrix from trait_effect_cov
#'
#' @param pop A `tidybreed_pop` object.
#' @param effect_name Character.
#' @param traits Character vector of trait names.
#' @return Named numeric matrix, or `NULL` if any entry is missing.
#' @keywords internal
load_effect_cov <- function(pop, effect_name, traits) {
  if (!"trait_effect_cov" %in% DBI::dbListTables(pop$db_conn)) return(NULL)
  n <- length(traits)
  R <- matrix(NA_real_, nrow = n, ncol = n, dimnames = list(traits, traits))
  rows <- DBI::dbGetQuery(
    pop$db_conn,
    paste0("SELECT trait_1, trait_2, cov FROM trait_effect_cov ",
           "WHERE effect_name = '", effect_name, "' ",
           "AND trait_1 IN (", paste0("'", traits, "'", collapse = ", "), ") ",
           "AND trait_2 IN (", paste0("'", traits, "'", collapse = ", "), ")")
  )
  if (nrow(rows) == 0L) return(NULL)
  for (i in seq_len(nrow(rows))) {
    R[rows$trait_1[i], rows$trait_2[i]] <- rows$cov[i]
  }
  if (any(is.na(R))) return(NULL)
  R
}


#' Write a single diagonal entry to trait_effect_cov
#'
#' Used internally by add_trait() and add_effect_random() to write a
#' per-trait variance as a 1x1 diagonal entry. Uses dbExecute() instead of
#' dbWriteTable() to avoid consuming R's RNG (DuckDB's dbWriteTable touches it).
#'
#' @param pop A `tidybreed_pop` object.
#' @param effect_name Character.
#' @param trait_name Character.
#' @param variance Numeric scalar.
#' @return The modified `tidybreed_pop` (invisibly).
#' @keywords internal
write_effect_cov_diagonal <- function(pop, effect_name, trait_name, variance) {
  pop <- ensure_effect_cov_table(pop)
  DBI::dbExecute(
    pop$db_conn,
    paste0("DELETE FROM trait_effect_cov WHERE effect_name = '", effect_name,
           "' AND trait_1 = '", trait_name, "' AND trait_2 = '", trait_name, "'")
  )
  DBI::dbExecute(
    pop$db_conn,
    paste0("INSERT INTO trait_effect_cov (effect_name, trait_1, trait_2, cov) VALUES ('",
           effect_name, "', '", trait_name, "', '", trait_name, "', ",
           format(as.numeric(variance), scientific = FALSE), ")")
  )
  invisible(pop)
}
