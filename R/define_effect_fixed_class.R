#' Define a discrete fixed-class effect in a phenotype model
#'
#' @description
#' Inserts a row into `trait_effects` representing a discrete fixed effect:
#' each level of `source_column` maps to a numeric shift on the phenotype
#' scale. The `null_class_action` argument controls what happens when an
#' individual has a `NULL` value for `source_column` at phenotyping time.
#'
#' @param pop A `tidybreed_pop` object.
#' @param phenotype_name Character. Name of an existing phenotype in
#'   `phenotype_meta`.
#' @param effect_name Character. Unique label for this effect within the
#'   phenotype. Must be a valid SQL identifier.
#' @param source_column Character. Column in `source_table` holding the group
#'   membership. Must be a valid SQL identifier.
#' @param levels Named numeric vector mapping level strings to shifts, e.g.
#'   `c(M = 30, F = 0)`. Every level present in the data at phenotyping time
#'   must appear here.
#' @param source_table Character. Database table containing `source_column`.
#'   Default `"ind_meta"`.
#' @param null_class_action Character. Behaviour when `source_column` is `NULL`
#'   for an individual:
#'   * `"skip"` (default) — individual is excluded from `ind_phenotype` for
#'     this phenotype; a warning lists the count.
#'   * `"error"` — hard stop when any `NULL` is encountered. Use for effects
#'     like `sex` that should never be missing.
#'   * `"zero"` — treat `NULL` as a zero-shift contribution (reference level).
#' @param overwrite Logical. Replace an existing effect with the same name.
#'
#' @return The modified `tidybreed_pop` (invisibly).
#'
#' @seealso [define_effect_int()], [define_effect_fixed_cov()],
#'   [define_effect_random()], [add_phenotype()]
#'
#' @examples
#' \dontrun{
#' pop <- pop |>
#'   define_effect_fixed_class("ADG", "sex",
#'     source_column = "sex",
#'     levels = c(M = 30, F = 0))
#' }
#' @export
define_effect_fixed_class <- function(pop,
                                      phenotype_name,
                                      effect_name,
                                      source_column,
                                      levels,
                                      source_table      = "ind_meta",
                                      null_class_action = c("skip", "error", "zero"),
                                      overwrite         = FALSE) {

  stopifnot(inherits(pop, "tidybreed_pop"))
  validate_tidybreed_pop(pop)
  validate_sql_identifier(phenotype_name, what = "phenotype name")
  validate_sql_identifier(effect_name,    what = "effect name")
  validate_sql_identifier(source_column,  what = "source_column")
  stopifnot(is.character(source_table), nzchar(source_table))
  null_class_action <- match.arg(null_class_action)

  if (!"trait_effects" %in% pop$tables) {
    stop("No trait tables exist. Call define_trait() / define_phenotype() first.",
         call. = FALSE)
  }

  .check_phenotype_exists(pop, phenotype_name)

  if (is.null(levels) || !is.numeric(levels) || is.null(names(levels))) {
    stop("`levels` must be a named numeric vector, e.g. c(M = 30, F = 0).",
         call. = FALSE)
  }
  if (any(!nzchar(names(levels)))) {
    stop("All entries in `levels` must have non-empty names.", call. = FALSE)
  }

  .handle_effect_overwrite(pop, phenotype_name, effect_name, overwrite)

  row <- tibble::tibble(
    phenotype_name    = phenotype_name,
    effect_name       = effect_name,
    effect_class      = "fixed_class",
    source_column     = source_column,
    source_table      = source_table,
    distribution      = NA_character_,
    levels_json       = encode_levels_json(levels),
    slope             = NA_real_,
    center            = NA_real_,
    value             = NA_real_,
    poly_order        = NA_integer_,
    null_class_action = null_class_action
  )
  DBI::dbWriteTable(pop$db_conn, "trait_effects", row, append = TRUE)

  message("Added fixed-class effect '", effect_name, "' to phenotype '",
          phenotype_name, "' (", length(levels), " levels).")
  invisible(pop)
}
