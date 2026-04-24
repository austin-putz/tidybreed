#' Add a fixed or random covariate to a trait's phenotype model
#'
#' @description
#' Inserts one row into `trait_effects` describing a non-genetic, non-residual
#' term in the phenotype equation for a trait. Two classes:
#'
#' * **Fixed effect**: a lookup table of discrete shifts per level of an
#'   `ind_meta` column. Pass `levels` as a named numeric vector, e.g.
#'   `c(M = 30, F = 0)`.
#' * **Random effect**: a variance component tied to a grouping column. One
#'   draw per distinct level of `source_column`, shared across all
#'   individuals in that level.
#'
#' The user is responsible for ensuring `source_column` exists in `ind_meta`
#' (via [mutate_ind_meta()] or as an extra column in [add_offspring()]).
#'
#' @param pop A `tidybreed_pop` object.
#' @param trait_name Character. Name of an existing trait.
#' @param effect_name Character. Label for this effect (unique within a
#'   trait). Valid SQL identifier.
#' @param effect_class Character. `"fixed"` or `"random"`.
#' @param source_column Character. Column in `ind_meta` used as the grouping
#'   variable. Required for both fixed and random effects.
#' @param levels For fixed effects: a **named numeric vector** mapping levels
#'   (e.g. `"M"`, `"F"`) to shifts on the trait scale. Ignored for random.
#' @param distribution Character. For random effects: `"normal"` (default),
#'   `"gamma"`, or `"uniform"`.
#' @param variance Numeric. Variance component for random effects.
#' @param value Numeric. A scalar constant (rarely used; normally the shift
#'   lives in `levels`).
#' @param overwrite Logical. If `TRUE` and a row with the same (trait,
#'   effect_name) exists, replace it.
#'
#' @return The modified `tidybreed_pop` (invisibly).
#'
#' @examples
#' \dontrun{
#' pop <- pop |>
#'   add_trait_covariate("ADG", "sex",
#'     effect_class = "fixed", source_column = "sex",
#'     levels = c(M = 30, F = 0)) |>
#'   add_trait_covariate("ADG", "litter",
#'     effect_class = "random", source_column = "id_litter",
#'     distribution = "normal", variance = 0.5)
#' }
#' @export
add_trait_covariate <- function(pop,  # nolint: object_name_linter
                                trait_name,
                                effect_name,
                                effect_class  = c("fixed", "random"),
                                source_column = NULL,
                                levels        = NULL,
                                distribution  = c("normal", "gamma", "uniform"),
                                variance      = NA_real_,
                                value         = NA_real_,
                                overwrite     = FALSE) {

  .Deprecated(
    msg = paste0(
      "add_trait_covariate() is deprecated. ",
      "Use add_effect_fixed_class(), add_effect_fixed_cov(), or ",
      "add_effect_random() instead."
    )
  )
  stopifnot(inherits(pop, "tidybreed_pop"))
  validate_tidybreed_pop(pop)
  validate_sql_identifier(trait_name, what = "trait name")
  validate_sql_identifier(effect_name, what = "effect name")
  effect_class <- match.arg(effect_class)
  distribution <- match.arg(distribution)
  if (!is.null(source_column)) {
    validate_sql_identifier(source_column, what = "source_column")
  }

  if (!"trait_effects" %in% pop$tables) {
    stop("No trait tables exist. Call add_trait() first.", call. = FALSE)
  }

  trait_exists <- DBI::dbGetQuery(
    pop$db_conn,
    paste0("SELECT COUNT(*) AS n FROM trait_meta WHERE trait_name = '",
           trait_name, "'")
  )$n
  if (trait_exists == 0) {
    stop("Trait '", trait_name, "' not found.", call. = FALSE)
  }

  if (is.null(source_column) || !nzchar(source_column)) {
    stop("`source_column` is required.", call. = FALSE)
  }

  if (effect_class == "fixed") {
    if (is.null(levels) || !is.numeric(levels) || is.null(names(levels))) {
      stop("Fixed effects require `levels` as a named numeric vector.",
           call. = FALSE)
    }
    if (any(!nzchar(names(levels)))) {
      stop("All `levels` entries must have non-empty names.", call. = FALSE)
    }
    levels_json <- encode_levels_json(levels)
    variance <- NA_real_
  } else {
    if (is.na(variance) || variance < 0) {
      stop("Random effects require a non-negative `variance`.", call. = FALSE)
    }
    levels_json <- NA_character_
    pop <- write_effect_cov_diagonal(pop, effect_name, trait_name, as.numeric(variance))
  }

  # Handle existing row
  existing <- DBI::dbGetQuery(
    pop$db_conn,
    paste0("SELECT COUNT(*) AS n FROM trait_effects ",
           "WHERE trait_name = '", trait_name,
           "' AND effect_name = '", effect_name, "'")
  )$n
  if (existing > 0 && !overwrite) {
    stop("Effect '", effect_name, "' for trait '", trait_name,
         "' already exists. Use overwrite = TRUE.", call. = FALSE)
  }
  if (existing > 0 && overwrite) {
    DBI::dbExecute(
      pop$db_conn,
      paste0("DELETE FROM trait_effects WHERE trait_name = '", trait_name,
             "' AND effect_name = '", effect_name, "'")
    )
  }

  DBI::dbExecute(
    pop$db_conn,
    paste0(
      "INSERT INTO trait_effects ",
      "(trait_name, effect_name, effect_class, source_column, source_table, ",
      "distribution, levels_json, slope, center, value) VALUES (",
      "'", trait_name, "', '", effect_name, "', '", effect_class, "', ",
      "'", source_column, "', 'ind_meta', ",
      if (effect_class == "random") paste0("'", distribution, "'") else "NULL", ", ",
      if (is.na(levels_json)) "NULL" else paste0("'", levels_json, "'"), ", ",
      "NULL, NULL, ",
      if (is.na(value)) "NULL" else as.numeric(value),
      ")"
    )
  )

  message("Added ", effect_class, " effect '", effect_name,
          "' to trait '", trait_name, "'.")
  invisible(pop)
}
