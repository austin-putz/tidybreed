#' Set the intercept (overall mean) for a trait
#'
#' @description
#' Updates `target_add_mean` in `trait_meta` for the named trait. This value
#' is the intercept added to every individual's liability before residuals and
#' fixed/random shifts. It can also be set at trait-definition time via
#' [add_trait()]; `add_effect_int()` lets you update it without redefining
#' the whole trait.
#'
#' @param pop A `tidybreed_pop` object.
#' @param trait_name Character. Name of an existing trait.
#' @param mean Numeric scalar. The new intercept value.
#'
#' @return The modified `tidybreed_pop` (invisibly).
#'
#' @seealso [add_trait()], [add_effect_fixed_class()], [add_effect_fixed_cov()],
#'   [add_effect_random()]
#'
#' @examples
#' \dontrun{
#' pop <- pop |> add_effect_int("ADG", mean = 850)
#' }
#' @export
add_effect_int <- function(pop, trait_name, mean) {

  stopifnot(inherits(pop, "tidybreed_pop"))
  validate_tidybreed_pop(pop)
  validate_sql_identifier(trait_name, what = "trait name")
  stopifnot(is.numeric(mean), length(mean) == 1, !is.na(mean))

  if (!"trait_meta" %in% pop$tables) {
    stop("No trait tables exist. Call add_trait() first.", call. = FALSE)
  }

  existing <- DBI::dbGetQuery(
    pop$db_conn,
    paste0("SELECT COUNT(*) AS n FROM trait_meta WHERE trait_name = '",
           trait_name, "'")
  )$n
  if (existing == 0) {
    stop("Trait '", trait_name, "' not found.", call. = FALSE)
  }

  DBI::dbExecute(
    pop$db_conn,
    paste0("UPDATE trait_meta SET target_add_mean = ", mean,
           " WHERE trait_name = '", trait_name, "'")
  )

  message("Set intercept for trait '", trait_name, "' to ", mean, ".")
  invisible(pop)
}
