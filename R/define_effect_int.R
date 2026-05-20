#' Set the intercept (population mean) for a phenotype
#'
#' @description
#' Updates `mean` in `phenotype_meta` for the named phenotype. This value is
#' the overall population mean added to every individual's liability before
#' residuals and fixed/random shifts. It can also be set at definition time via
#' the `mean` argument in [define_phenotype()]; `define_effect_int()` lets you
#' update it without redefining the whole phenotype.
#'
#' @param pop A `tidybreed_pop` object.
#' @param phenotype_name Character. Name of an existing phenotype in
#'   `phenotype_meta`.
#' @param mean Numeric scalar. The new population mean / intercept value.
#'
#' @return The modified `tidybreed_pop` (invisibly).
#'
#' @seealso [define_phenotype()], [define_effect_fixed_class()],
#'   [define_effect_fixed_cov()], [define_effect_random()]
#'
#' @examples
#' \dontrun{
#' pop <- pop |> define_effect_int("ADG", mean = 850)
#' }
#' @export
define_effect_int <- function(pop, phenotype_name, mean) {

  stopifnot(inherits(pop, "tidybreed_pop"))
  validate_tidybreed_pop(pop)
  validate_sql_identifier(phenotype_name, what = "phenotype name")
  stopifnot(is.numeric(mean), length(mean) == 1, !is.na(mean))

  if (!"phenotype_meta" %in% pop$tables) {
    stop("No phenotype tables exist. Call define_phenotype() first.", call. = FALSE)
  }

  pn_safe <- gsub("'", "''", phenotype_name)
  existing <- DBI::dbGetQuery(
    pop$db_conn,
    paste0("SELECT COUNT(*) AS n FROM phenotype_meta WHERE phenotype_name = '",
           pn_safe, "'")
  )$n
  if (existing == 0) {
    stop("Phenotype '", phenotype_name, "' not found in phenotype_meta. ",
         "Call define_phenotype() first.", call. = FALSE)
  }

  DBI::dbExecute(
    pop$db_conn,
    paste0("UPDATE phenotype_meta SET mean = ", mean,
           " WHERE phenotype_name = '", pn_safe, "'")
  )

  message("Set mean for phenotype '", phenotype_name, "' to ", mean, ".")
  invisible(pop)
}
