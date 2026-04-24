#' Add a discrete fixed-class effect to a trait's phenotype model
#'
#' @description
#' Inserts a row into `trait_effects` representing a discrete fixed effect:
#' each level of `source_column` in `source_table` maps to a numeric shift on
#' the trait scale. At phenotyping time, an individual whose level is not in
#' `levels` causes an error (the user must supply a complete map).
#'
#' @param pop A `tidybreed_pop` object.
#' @param trait_name Character. Name of an existing trait.
#' @param effect_name Character. Unique label for this effect within the trait.
#'   Must be a valid SQL identifier.
#' @param source_column Character. Column in `source_table` that holds the
#'   group membership for each individual. Must be a valid SQL identifier.
#' @param levels Named numeric vector mapping level strings to shifts, e.g.
#'   `c(M = 30, F = 0)`. Every level present in the data at phenotyping time
#'   must appear here.
#' @param source_table Character. Database table that contains `source_column`.
#'   Defaults to `"ind_meta"`. Must contain an `id_ind` column.
#' @param overwrite Logical. Replace an existing effect with the same name.
#'
#' @return The modified `tidybreed_pop` (invisibly).
#'
#' @seealso [add_effect_int()], [add_effect_fixed_cov()], [add_effect_random()],
#'   [add_phenotype()]
#'
#' @examples
#' \dontrun{
#' pop <- pop |>
#'   add_effect_fixed_class("ADG", "sex",
#'     source_column = "sex",
#'     levels = c(M = 30, F = 0))
#' }
#' @export
add_effect_fixed_class <- function(pop,
                                   trait_name,
                                   effect_name,
                                   source_column,
                                   levels,
                                   source_table = "ind_meta",
                                   overwrite    = FALSE) {

  stopifnot(inherits(pop, "tidybreed_pop"))
  validate_tidybreed_pop(pop)
  validate_sql_identifier(trait_name,    what = "trait name")
  validate_sql_identifier(effect_name,   what = "effect name")
  validate_sql_identifier(source_column, what = "source_column")
  stopifnot(is.character(source_table), nzchar(source_table))

  if (!"trait_effects" %in% pop$tables) {
    stop("No trait tables exist. Call add_trait() first.", call. = FALSE)
  }

  .check_trait_exists(pop, trait_name)

  if (is.null(levels) || !is.numeric(levels) || is.null(names(levels))) {
    stop("`levels` must be a named numeric vector, e.g. c(M = 30, F = 0).",
         call. = FALSE)
  }
  if (any(!nzchar(names(levels)))) {
    stop("All entries in `levels` must have non-empty names.", call. = FALSE)
  }

  .handle_effect_overwrite(pop, trait_name, effect_name, overwrite)

  row <- tibble::tibble(
    trait_name    = trait_name,
    effect_name   = effect_name,
    effect_class  = "fixed_class",
    source_column = source_column,
    source_table  = source_table,
    distribution  = NA_character_,
    variance      = NA_real_,
    levels_json   = encode_levels_json(levels),
    slope         = NA_real_,
    center        = NA_real_,
    value         = NA_real_
  )
  DBI::dbWriteTable(pop$db_conn, "trait_effects", row, append = TRUE)

  message("Added fixed-class effect '", effect_name, "' to trait '",
          trait_name, "' (", length(levels), " levels).")
  invisible(pop)
}
