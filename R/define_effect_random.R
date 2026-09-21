#' Define a random group effect in a phenotype model
#'
#' @description
#' Inserts a row into `phenotype_effects` for a random effect. One value is drawn
#' per distinct level of `source_column`; all individuals sharing that level
#' receive the same shift. Drawn values are stored in `phenotype_random_effects`
#' so they are reproducible across repeated calls to [add_phenotype()] without
#' requiring a fixed `seed`.
#'
#' To correlate this effect across multiple phenotypes (e.g. the same herd
#' affects both ADG and BW), call [define_effect_cov_matrix()] with the
#' appropriate `effect_name` — either before or after this call. Once the
#' phenotype belongs to a block of two or more phenotypes for `effect_name`,
#' this call must use `distribution = "normal"` and the same
#' `(source_column, source_table)` as the block's other members, and
#' `variance` can no longer be set here — the block is redeclared as a whole
#' with [define_effect_cov_matrix()].
#'
#' @param pop A `tidybreed_pop` object.
#' @param phenotype_name Character. Name of an existing phenotype in
#'   `phenotype_meta`.
#' @param effect_name Character. Unique label for this effect within the
#'   phenotype.
#' @param source_column Character. Column in `source_table` whose distinct
#'   values define the groups (e.g. `"herd_id"`, `"litter"`, `"id_ind"` for PE).
#' @param variance Numeric scalar or `NULL`. Variance of the random effect.
#'   `NULL` (default) uses the value already stored in `phenotype_var_comp` via
#'   [define_effect_cov_matrix()] and errors if there is none. A number writes
#'   (or overwrites) a 1 × 1 block for this phenotype; it is an error when the
#'   phenotype is already in a multi-phenotype block for `effect_name`.
#' @param distribution Character. Sampling distribution: `"normal"` (default),
#'   `"gamma"`, or `"uniform"`.
#' @param source_table Character. Table containing `source_column`. Default
#'   `"ind_meta"`.
#' @param overwrite Logical. Replace an existing effect with the same name.
#'   The stored draws of that effect for this phenotype in
#'   `phenotype_random_effects` are discarded with it.
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
  if (identical(effect_name, "residual")) {
    stop("'residual' is reserved for the residual; use define_residual_cov().",
         call. = FALSE)
  }
  if (!is.null(variance) &&
      (!is.numeric(variance) || length(variance) != 1 ||
       is.na(variance) || variance < 0)) {
    stop("`variance` must be a non-negative number.", call. = FALSE)
  }

  .check_phenotype_exists(pop, phenotype_name)
  conn <- pop$db_conn

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

  # One transaction: drop the effect being overwritten (and its draws), write
  # the variance through the block writer, check the block, insert the row.
  DBI::dbExecute(conn, "BEGIN TRANSACTION")
  committed <- FALSE
  on.exit(if (!committed) try(DBI::dbExecute(conn, "ROLLBACK"), silent = TRUE),
          add = TRUE)

  .handle_effect_overwrite(pop, phenotype_name, effect_name, overwrite)

  if (!is.null(variance)) {
    .pvc_write_block(
      conn, effect_name, phenotype_name,
      matrix(as.numeric(variance), 1L, 1L,
             dimnames = list(phenotype_name, phenotype_name)),
      caller = "define_effect_random(variance = )")
  } else {
    variance <- get_phenotype_var(pop, effect_name, phenotype_name)
    if (is.na(variance)) {
      stop("No variance found in phenotype_var_comp for effect '", effect_name,
           "' / phenotype '", phenotype_name, "'. ",
           "Either call define_effect_cov_matrix() first or supply `variance`.",
           call. = FALSE)
    }
  }

  block <- .pvc_block_members(conn, effect_name, phenotype_name)
  validate_named_effect_block(conn, effect_name, block,
                              pending = as.data.frame(row),
                              caller  = "define_effect_random()")

  tmp <- "__define_effect_random_tmp"
  duckdb::duckdb_register(conn, tmp, as.data.frame(row))
  on.exit(try(duckdb::duckdb_unregister(conn, tmp), silent = TRUE), add = TRUE)
  DBI::dbExecute(conn, sprintf(
    "INSERT INTO phenotype_effects (%s) SELECT %s FROM %s",
    paste(names(row), collapse = ", "), paste(names(row), collapse = ", "), tmp))
  DBI::dbExecute(conn, "COMMIT")
  committed <- TRUE

  message("Added random effect '", effect_name, "' to phenotype '", phenotype_name,
          "' (variance = ", variance, ", distribution = ", distribution, ").")
  invisible(pop)
}
