#' Store a residual covariance matrix across traits
#'
#' @description
#' Writes the off-diagonal and diagonal entries of a residual covariance
#' matrix `R` to the `trait_residual_cov` table. Consumed by [add_phenotype()]
#' when multiple traits are phenotyped on the same filtered subset, causing
#' residuals to be drawn jointly from MVN(0, R) instead of independently.
#'
#' Diagonal entries of `R` must match each trait's stored `residual_var`
#' (within a small tolerance); otherwise an error is raised. This keeps
#' `trait_meta` and `trait_residual_cov` consistent.
#'
#' Previous entries involving any of the listed traits are overwritten.
#'
#' @param pop A `tidybreed_pop` object.
#' @param traits Character vector of trait names (length >= 2).
#' @param R Numeric matrix, square with side = `length(traits)`, symmetric
#'   positive semi-definite.
#' @param tol Numeric tolerance for diagonal match (default `1e-6`).
#'
#' @return The modified `tidybreed_pop` (invisibly).
#'
#' @seealso [set_qtl_effects_multi()], [add_phenotype()]
#'
#' @examples
#' \dontrun{
#' R <- matrix(c(0.75, 0.20, 0.20, 0.70), 2, 2)
#' pop <- pop |>
#'   set_residual_cov(traits = c("ADG", "BW"), R = R)
#' }
#' @export
set_residual_cov <- function(pop, traits, R, tol = 1e-6) {

  stopifnot(inherits(pop, "tidybreed_pop"))
  validate_tidybreed_pop(pop)
  stopifnot(is.character(traits), length(traits) >= 2)
  lapply(traits, validate_sql_identifier, what = "trait name")

  if (!is.matrix(R) || nrow(R) != length(traits) || ncol(R) != length(traits)) {
    stop("`R` must be a square matrix with side = length(traits).",
         call. = FALSE)
  }
  if (!isSymmetric(unname(R))) {
    stop("`R` must be symmetric.", call. = FALSE)
  }

  if (!"trait_meta" %in% pop$tables) {
    stop("No traits defined yet. Call add_trait() first.", call. = FALSE)
  }

  rows <- DBI::dbGetQuery(
    pop$db_conn,
    paste0("SELECT trait_name, residual_var FROM trait_meta WHERE trait_name IN (",
           paste0("'", traits, "'", collapse = ", "), ")")
  )
  missing_traits <- setdiff(traits, rows$trait_name)
  if (length(missing_traits) > 0) {
    stop("Traits not found in trait_meta: ",
         paste(missing_traits, collapse = ", "), call. = FALSE)
  }

  expected_diag <- rows$residual_var[match(traits, rows$trait_name)]
  observed_diag <- diag(R)
  if (any(abs(expected_diag - observed_diag) > tol)) {
    mismatch <- traits[abs(expected_diag - observed_diag) > tol]
    stop("Diagonal of R does not match trait_meta.residual_var for: ",
         paste(mismatch, collapse = ", "),
         ". Update trait_meta first (e.g. by re-adding the trait).",
         call. = FALSE)
  }

  # Clear old entries involving any listed trait
  quoted <- paste0("'", traits, "'", collapse = ", ")
  DBI::dbExecute(
    pop$db_conn,
    paste0("DELETE FROM trait_residual_cov WHERE trait_1 IN (", quoted,
           ") OR trait_2 IN (", quoted, ")")
  )

  # Write all unordered pairs + diagonal
  n <- length(traits)
  entries <- list()
  for (i in seq_len(n)) {
    for (j in seq_len(n)) {
      entries[[length(entries) + 1]] <- tibble::tibble(
        trait_1 = traits[i],
        trait_2 = traits[j],
        cov     = R[i, j]
      )
    }
  }
  df <- dplyr::bind_rows(entries)
  DBI::dbWriteTable(pop$db_conn, "trait_residual_cov", df, append = TRUE)

  message("Stored residual covariance for traits: ",
          paste(traits, collapse = ", "))
  invisible(pop)
}
