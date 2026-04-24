#' Add a continuous covariate (regression) effect to a trait's phenotype model
#'
#' @description
#' Inserts a row into `trait_effects` representing a linear regression term.
#' At phenotyping time the contribution per individual is
#' `slope * (source_column_value - center)`. When `center` is `NULL` the mean
#' of `source_column` across all rows currently in `source_table` is computed
#' and stored automatically.
#'
#' @param pop A `tidybreed_pop` object.
#' @param trait_name Character. Name of an existing trait.
#' @param effect_name Character. Unique label for this effect within the trait.
#' @param source_column Character. Column in `source_table` holding the
#'   continuous predictor values.
#' @param slope Numeric scalar. Regression coefficient (change in trait per
#'   unit of the predictor).
#' @param center Numeric scalar or `NULL`. The value subtracted from
#'   `source_column` before multiplying by `slope`. If `NULL`, the column mean
#'   is computed from `source_table` at definition time and stored.
#' @param source_table Character. Database table that contains `source_column`.
#'   Defaults to `"ind_meta"`. Must contain an `id_ind` column.
#' @param overwrite Logical. Replace an existing effect with the same name.
#'
#' @return The modified `tidybreed_pop` (invisibly).
#'
#' @seealso [add_effect_int()], [add_effect_fixed_class()], [add_effect_random()],
#'   [add_phenotype()]
#'
#' @examples
#' \dontrun{
#' # center auto-computed from current ind_meta$age_days
#' pop <- pop |>
#'   add_effect_fixed_cov("ADG", "age",
#'     source_column = "age_days",
#'     slope = 2.5)
#' }
#' @export
add_effect_fixed_cov <- function(pop,
                                 trait_name,
                                 effect_name,
                                 source_column,
                                 slope,
                                 center       = NULL,
                                 source_table = "ind_meta",
                                 overwrite    = FALSE) {

  stopifnot(inherits(pop, "tidybreed_pop"))
  validate_tidybreed_pop(pop)
  validate_sql_identifier(trait_name,    what = "trait name")
  validate_sql_identifier(effect_name,   what = "effect name")
  validate_sql_identifier(source_column, what = "source_column")
  stopifnot(is.character(source_table), nzchar(source_table))
  stopifnot(is.numeric(slope), length(slope) == 1, !is.na(slope))

  if (!"trait_effects" %in% pop$tables) {
    stop("No trait tables exist. Call add_trait() first.", call. = FALSE)
  }

  .check_trait_exists(pop, trait_name)

  if (!source_table %in% DBI::dbListTables(pop$db_conn)) {
    stop("Source table '", source_table, "' does not exist in the database.",
         call. = FALSE)
  }

  if (is.null(center)) {
    res <- DBI::dbGetQuery(
      pop$db_conn,
      paste0("SELECT AVG(CAST(", source_column, " AS DOUBLE)) AS center FROM ",
             source_table)
    )
    center <- res$center
    if (is.na(center)) {
      stop("Could not auto-compute center: '", source_column, "' in '",
           source_table, "' has all NULL values.", call. = FALSE)
    }
    message("Auto-computed center for '", effect_name, "': ", round(center, 6),
            " (mean of ", source_column, " in ", source_table, ")")
  } else {
    stopifnot(is.numeric(center), length(center) == 1, !is.na(center))
  }

  .handle_effect_overwrite(pop, trait_name, effect_name, overwrite)

  row <- tibble::tibble(
    trait_name    = trait_name,
    effect_name   = effect_name,
    effect_class  = "fixed_cov",
    source_column = source_column,
    source_table  = source_table,
    distribution  = NA_character_,
    variance      = NA_real_,
    levels_json   = NA_character_,
    slope         = as.numeric(slope),
    center        = as.numeric(center),
    value         = NA_real_
  )
  DBI::dbWriteTable(pop$db_conn, "trait_effects", row, append = TRUE)

  message("Added fixed-covariate effect '", effect_name, "' to trait '",
          trait_name, "' (slope = ", slope, ", center = ", round(center, 4), ").")
  invisible(pop)
}
