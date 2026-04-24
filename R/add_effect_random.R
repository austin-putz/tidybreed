#' Add a random group effect to a trait's phenotype model
#'
#' @description
#' Inserts a row into `trait_effects` for a random effect. One value is drawn
#' per distinct level of `source_column`; all individuals sharing that level
#' receive the same shift. Drawn values are stored in `trait_random_effects`
#' so they are reproducible across repeated calls to [add_phenotype()] without
#' requiring a fixed `seed`.
#'
#' To correlate this effect across multiple traits (e.g. the same herd affects
#' both ADG and BW), call [set_random_effect_cov()] after adding the effect to
#' each trait.
#'
#' @param pop A `tidybreed_pop` object.
#' @param trait_name Character. Name of an existing trait.
#' @param effect_name Character. Unique label for this effect within the trait.
#' @param source_column Character. Column in `source_table` whose distinct
#'   values define the groups (e.g. `"herd_id"`, `"litter"`).
#' @param variance Numeric. Variance of the random effect distribution.
#' @param distribution Character. Sampling distribution: `"normal"` (default),
#'   `"gamma"`, or `"uniform"`.
#' @param source_table Character. Database table that contains `source_column`.
#'   Defaults to `"ind_meta"`. Must contain an `id_ind` column.
#' @param overwrite Logical. Replace an existing effect with the same name.
#'
#' @return The modified `tidybreed_pop` (invisibly).
#'
#' @seealso [set_random_effect_cov()], [add_effect_fixed_class()],
#'   [add_effect_fixed_cov()], [add_phenotype()]
#'
#' @examples
#' \dontrun{
#' pop <- pop |>
#'   add_effect_random("ADG", "herd",
#'     source_column = "herd_id",
#'     variance = 150)
#' }
#' @export
add_effect_random <- function(pop,
                              trait_name,
                              effect_name,
                              source_column,
                              variance,
                              distribution = c("normal", "gamma", "uniform"),
                              source_table = "ind_meta",
                              overwrite    = FALSE) {

  stopifnot(inherits(pop, "tidybreed_pop"))
  validate_tidybreed_pop(pop)
  validate_sql_identifier(trait_name,    what = "trait name")
  validate_sql_identifier(effect_name,   what = "effect name")
  validate_sql_identifier(source_column, what = "source_column")
  stopifnot(is.character(source_table), nzchar(source_table))
  distribution <- match.arg(distribution)

  if (is.na(variance) || !is.numeric(variance) || variance < 0) {
    stop("`variance` must be a non-negative number.", call. = FALSE)
  }

  if (!"trait_effects" %in% pop$tables) {
    stop("No trait tables exist. Call add_trait() first.", call. = FALSE)
  }

  .check_trait_exists(pop, trait_name)
  .handle_effect_overwrite(pop, trait_name, effect_name, overwrite)

  row <- tibble::tibble(
    trait_name    = trait_name,
    effect_name   = effect_name,
    effect_class  = "random",
    source_column = source_column,
    source_table  = source_table,
    distribution  = distribution,
    variance      = as.numeric(variance),
    levels_json   = NA_character_,
    slope         = NA_real_,
    center        = NA_real_,
    value         = NA_real_
  )
  DBI::dbWriteTable(pop$db_conn, "trait_effects", row, append = TRUE)

  message("Added random effect '", effect_name, "' to trait '", trait_name,
          "' (variance = ", variance, ", distribution = ", distribution, ").")
  invisible(pop)
}
