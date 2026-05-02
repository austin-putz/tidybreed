#' Define a selection index
#'
#' @description
#' Registers a named selection index in the `index_meta` table by specifying
#' which traits are included and their weighting coefficients. Call this once
#' per index (before calling [add_index()]). Re-calling with the same
#' `index_name` and `trait_name` pairs updates the weights in place.
#'
#' Not all traits in the population need to appear in an index — you may have
#' five traits but index only three. Traits not listed in `trait_names` are
#' ignored when [add_index()] computes the index value.
#'
#' @param pop A `tidybreed_pop` object.
#' @param index_name Character scalar. Name of the index (e.g. `"terminal"`,
#'   `"maternal"`). Must be a valid SQL identifier.
#' @param trait_names Character vector. Trait names to include in the index.
#'   Each must be a valid SQL identifier. Order must match `index_wts`.
#' @param index_wts Numeric vector. Selection index weights, one per trait in
#'   `trait_names`. Positive weights favour higher trait values; negative
#'   weights favour lower values.
#' @param ... Optional extra columns to add to `index_meta`. Scalars are
#'   broadcast to all rows (one per trait). Vectors must have length equal to
#'   `length(trait_names)`. Types are inferred via [infer_duckdb_type()].
#'
#' @return The `tidybreed_pop` (invisibly). Assign the result back.
#'
#' @seealso [add_index()], [add_trait()]
#'
#' @examples
#' \dontrun{
#' pop <- define_index(pop, "terminal",
#'                     trait_names = c("ADG", "FCR"),
#'                     index_wts   = c(1.2, -0.8))
#'
#' # With an extra metadata column
#' pop <- define_index(pop, "maternal",
#'                     trait_names = c("LS", "MW"),
#'                     index_wts   = c(0.5, -0.3),
#'                     source = "genetic_team_v2")
#' }
#' @export
define_index <- function(pop, index_name, trait_names, index_wts, ...) {

  # ---- Input validation ----
  validate_tidybreed_pop(pop)

  stopifnot(is.character(index_name), length(index_name) == 1)
  validate_sql_identifier(index_name, what = "index name")

  stopifnot(is.character(trait_names), length(trait_names) >= 1)
  lapply(trait_names, validate_sql_identifier, what = "trait name")

  stopifnot(is.numeric(index_wts), length(index_wts) == length(trait_names))

  extra_cols <- list(...)
  n_rows <- length(trait_names)

  # Scalars are fine; vectors must match trait count
  for (nm in names(extra_cols)) {
    v <- extra_cols[[nm]]
    if (length(v) != 1L && length(v) != n_rows) {
      stop(
        "Vector length (", length(v), ") for field '", nm,
        "' must be 1 (scalar) or equal to length(trait_names) (", n_rows, ").",
        call. = FALSE
      )
    }
  }

  # ---- Build insertion data frame ----
  df <- data.frame(
    index_name = index_name,
    trait_name = trait_names,
    index_wt   = as.double(index_wts),
    stringsAsFactors = FALSE
  )

  # ---- Attach extra columns ----
  if (length(extra_cols) > 0) {
    prepped <- prepare_extra_cols(extra_cols, n_rows, "index_meta", pop$db_conn)
    for (nm in names(prepped)) df[[nm]] <- prepped[[nm]]
  }

  # ---- Upsert into index_meta ----
  tmp_tbl <- paste0("__define_index_tmp_", as.integer(Sys.time()))
  DBI::dbWriteTable(pop$db_conn, tmp_tbl, df, temporary = TRUE, overwrite = TRUE)

  non_key_cols <- setdiff(names(df), c("index_name", "trait_name"))
  update_clause <- paste(
    vapply(non_key_cols, function(col) {
      paste0(col, " = EXCLUDED.", col)
    }, character(1)),
    collapse = ", "
  )

  sql <- paste0(
    "INSERT INTO index_meta SELECT * FROM ", tmp_tbl,
    " ON CONFLICT (index_name, trait_name) DO UPDATE SET ", update_clause
  )
  DBI::dbExecute(pop$db_conn, sql)
  DBI::dbExecute(pop$db_conn, paste0("DROP TABLE IF EXISTS ", tmp_tbl))

  # ---- User message ----
  wt_str <- paste(
    vapply(seq_along(trait_names), function(i) {
      sprintf("%s (%.4g)", trait_names[i], index_wts[i])
    }, character(1)),
    collapse = ", "
  )
  message("Defined index '", index_name, "': ", wt_str)

  invisible(pop)
}
