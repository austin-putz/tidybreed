#' Define a continuous covariate (regression) effect in a phenotype model
#'
#' @description
#' Inserts a row into `trait_effects` representing a polynomial regression term.
#' At phenotyping time the contribution per individual is
#' `slope * (source_column_value - center)^poly_order`. When `center` is `NULL`
#' the mean of `source_column` is computed and stored automatically.
#'
#' @param pop A `tidybreed_pop` object.
#' @param phenotype_name Character. Name of an existing phenotype in
#'   `phenotype_meta`.
#' @param effect_name Character. Unique label for this effect within the
#'   phenotype.
#' @param source_column Character. Column in `source_table` holding the
#'   continuous predictor values.
#' @param slope Numeric scalar. Regression coefficient.
#' @param center Numeric scalar or `NULL`. Value subtracted before raising to
#'   `poly_order`. If `NULL`, auto-computed as the column mean.
#' @param poly_order Integer ≥ 1. Polynomial order. Default `1` gives the
#'   standard linear `slope * (x - center)` term. Use `2` for quadratic, etc.
#' @param source_table Character. Database table containing `source_column`.
#'   Default `"ind_meta"`.
#' @param overwrite Logical. Replace an existing effect with the same name.
#'
#' @return The modified `tidybreed_pop` (invisibly).
#'
#' @seealso [define_effect_intercept()], [define_effect_fixed_class()],
#'   [define_effect_random()], [add_phenotype()]
#'
#' @examples
#' \dontrun{
#' # Linear regression on age (center auto-computed)
#' pop <- pop |>
#'   define_effect_fixed_cov("ADG", "age",
#'     source_column = "age_days",
#'     slope = 2.5)
#'
#' # Quadratic regression on DIM (dairy)
#' pop <- pop |>
#'   define_effect_fixed_cov("milk_td", "dim_sq",
#'     source_column = "dim",
#'     slope  = -0.02,
#'     center = 150,
#'     poly_order = 2L)
#' }
#' @export
define_effect_fixed_cov <- function(pop,
                                    phenotype_name,
                                    effect_name,
                                    source_column,
                                    slope,
                                    center       = NULL,
                                    poly_order   = 1L,
                                    source_table = "ind_meta",
                                    overwrite    = FALSE) {

  stopifnot(inherits(pop, "tidybreed_pop"))
  validate_tidybreed_pop(pop)
  validate_sql_identifier(phenotype_name, what = "phenotype name")
  validate_sql_identifier(effect_name,    what = "effect name")
  validate_sql_identifier(source_column,  what = "source_column")
  stopifnot(is.character(source_table), nzchar(source_table))
  stopifnot(is.numeric(slope), length(slope) == 1, !is.na(slope))

  poly_order <- as.integer(poly_order)
  if (is.na(poly_order) || poly_order < 1L) {
    stop("`poly_order` must be an integer >= 1.", call. = FALSE)
  }

  .check_phenotype_exists(pop, phenotype_name)

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

  .handle_effect_overwrite(pop, phenotype_name, effect_name, overwrite)

  row <- tibble::tibble(
    phenotype_name    = phenotype_name,
    effect_name       = effect_name,
    effect_class      = "fixed_cov",
    source_column     = source_column,
    source_table      = source_table,
    distribution      = NA_character_,
    levels_json       = NA_character_,
    slope             = as.numeric(slope),
    center            = as.numeric(center),
    value             = NA_real_,
    poly_order        = poly_order,
    null_class_action = NA_character_
  )
  DBI::dbWriteTable(pop$db_conn, "trait_effects", row, append = TRUE)

  message("Added fixed-covariate effect '", effect_name, "' to phenotype '",
          phenotype_name, "' (slope = ", slope,
          ", center = ", round(center, 4),
          ", poly_order = ", poly_order, ").")
  invisible(pop)
}
