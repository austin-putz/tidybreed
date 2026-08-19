#' Define a random group effect in a phenotype model
#'
#' @description
#' Inserts a row into `trait_effects` for a random effect. One value is drawn
#' per distinct level of `source_column`; all individuals sharing that level
#' receive the same shift. Drawn values are stored in `trait_random_effects`
#' so they are reproducible across repeated calls to [add_phenotype()] without
#' requiring a fixed `seed`.
#'
#' To correlate this effect across multiple phenotypes (e.g. the same herd
#' affects both ADG and BW), call [define_effect_cov_matrix()] with the
#' appropriate `effect_name` — either before or after this call.
#'
#' @param pop A `tidybreed_pop` object.
#' @param phenotype_name Character. Name of an existing phenotype in
#'   `phenotype_meta`.
#' @param effect_name Character. Unique label for this effect within the
#'   phenotype.
#' @param source_column Character. Column in `source_table` whose distinct
#'   values define the groups (e.g. `"herd_id"`, `"litter"`, `"id_ind"` for PE).
#' @param variance Numeric scalar. Variance of the random effect. Optional if
#'   already stored via [define_effect_cov_matrix()]; required otherwise.
#' @param distribution Character. Sampling distribution: `"normal"` (default),
#'   `"gamma"`, or `"uniform"`.
#' @param source_table Character. Table containing `source_column`. Default
#'   `"ind_meta"`.
#' @param overwrite Logical. Replace an existing effect with the same name.
#'
#' @return The modified `tidybreed_pop` (invisibly).
#'
#' @seealso [define_effect_cov_matrix()], [define_effect_fixed_class()],
#'   [define_effect_fixed_cov()], [add_phenotype()]
#'
#' @examples
#' \dontrun{
#' # Herd random effect
#' pop <- pop |>
#'   define_effect_random("ADG", "herd",
#'     source_column = "herd_id",
#'     variance = 150)
#'
#' # Permanent environment (PE) for repeatability — one draw per animal
#' pop <- pop |>
#'   define_effect_random("litter_size", "pe",
#'     source_column = "id_ind",
#'     variance = 0.3)
#' }
#' @export
define_effect_random <- function(pop,
                                 phenotype_name,
                                 effect_name,
                                 source_column,
                                 variance     = NULL,
                                 distribution = c("normal", "gamma", "uniform"),
                                 source_table = "ind_meta",
                                 overwrite    = FALSE) {

  stopifnot(inherits(pop, "tidybreed_pop"))
  validate_tidybreed_pop(pop)
  validate_sql_identifier(phenotype_name, what = "phenotype name")
  validate_sql_identifier(effect_name,    what = "effect name")
  validate_sql_identifier(source_column,  what = "source_column")
  stopifnot(is.character(source_table), nzchar(source_table))
  distribution <- match.arg(distribution)

  # Resolve variance: check phenotype_var_comp first
  stored_var <- get_phenotype_var(pop, effect_name, phenotype_name)
  if (!is.na(stored_var)) {
    variance <- stored_var
  } else {
    if (is.null(variance)) {
      stop("No variance found in phenotype_var_comp for effect '", effect_name,
           "' / phenotype '", phenotype_name, "'. ",
           "Either call define_effect_cov_matrix() first or supply `variance`.",
           call. = FALSE)
    }
    if (!is.numeric(variance) || length(variance) != 1 ||
        is.na(variance) || variance < 0) {
      stop("`variance` must be a non-negative number.", call. = FALSE)
    }
    pop <- write_phenotype_var_diag(pop, effect_name, phenotype_name,
                                    as.numeric(variance))
  }

  .check_phenotype_exists(pop, phenotype_name)
  .handle_effect_overwrite(pop, phenotype_name, effect_name, overwrite)

  row <- tibble::tibble(
    phenotype_name    = phenotype_name,
    effect_name       = effect_name,
    effect_class      = "random",
    source_column     = source_column,
    source_table      = source_table,
    distribution      = distribution,
    levels_json       = NA_character_,
    slope             = NA_real_,
    center            = NA_real_,
    value             = NA_real_,
    poly_order        = NA_integer_,
    null_class_action = NA_character_
  )
  DBI::dbWriteTable(pop$db_conn, "trait_effects", row, append = TRUE)

  message("Added random effect '", effect_name, "' to phenotype '", phenotype_name,
          "' (variance = ", variance, ", distribution = ", distribution, ").")
  invisible(pop)
}
