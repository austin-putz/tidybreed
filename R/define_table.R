#' Define a custom table in the population database
#'
#' @description
#' Creates a new user-defined table inside the tidybreed DuckDB database and
#' registers it in `pop$tables` so that [get_table()] and [mutate_table()] work
#' on it immediately.
#'
#' Column names and types are specified as named `...` arguments using the same
#' typed-NA convention as [mutate_table()] schema pre-declaration:
#'
#' ```r
#' pop |> define_table(
#'   "sim_timing",
#'   run_id       = NA_integer_,
#'   duration_sec = NA_real_,
#'   label        = NA_character_
#' )
#' ```
#'
#' Types are inferred by [infer_duckdb_type()]. Use R's typed NA values
#' (`NA_integer_`, `NA_real_`, `NA_character_`, `as.Date(NA)`, etc.) to get the
#' intended DuckDB column type. There is no typed logical NA in R, so for a
#' `BOOLEAN` column pass a concrete placeholder (`FALSE` or `TRUE`) instead.
#' Bare `NA` defaults to `VARCHAR` with a warning.
#'
#' After creation, add rows via `DBI::dbAppendTable(pop$db_conn, table_name,
#' my_tibble)` and query via `get_table(pop, table_name) |> dplyr::collect()`.
#'
#' @param pop A `tidybreed_pop` object.
#' @param table_name Character scalar. Name for the new table. Must be a valid
#'   SQL identifier and must not conflict with a system-managed tidybreed table.
#' @param ... Named column definitions. Each name becomes a column; the value's
#'   R type determines the DuckDB column type. Use typed NAs to declare columns
#'   without supplying real data (e.g. `run_id = NA_integer_`).
#' @param primary_key Character scalar (optional). Name of one of the `...`
#'   columns to declare as the `PRIMARY KEY`. When supplied, `mutate_table()`
#'   filtered-row updates will work on this table (DuckDB enforces uniqueness).
#'   Defaults to `NULL` (no primary key).
#' @param overwrite Logical. If `FALSE` (default), an error is raised when
#'   `table_name` already exists. If `TRUE`, the existing table is dropped and
#'   recreated (all data in it will be lost).
#'
#' @return The `tidybreed_pop` (invisibly). Assign the result back.
#'
#' @seealso [get_table()], [mutate_table()]
#'
#' @examples
#' \dontrun{
#' pop <- open_pop(pop_name = "demo", db_name = ":memory:") |>
#'   define_genome(n_loci = 200, n_chr = 2, chr_len_Mb = 100)
#'
#' # Create a custom table for tracking simulation run metadata
#' pop <- pop |> define_table(
#'   "sim_timing",
#'   run_id       = NA_integer_,
#'   duration_sec = NA_real_,
#'   label        = NA_character_,
#'   primary_key  = "run_id"
#' )
#'
#' # Insert rows directly via DBI
#' DBI::dbAppendTable(
#'   pop$db_conn, "sim_timing",
#'   tibble::tibble(run_id = 1L, duration_sec = 42.3, label = "baseline")
#' )
#'
#' # Query via the tidy interface
#' get_table(pop, "sim_timing") |> dplyr::collect()
#'
#' # Add more columns later
#' pop <- pop |> get_table("sim_timing") |> mutate_table(notes = NA_character_)
#'
#' close_pop(pop)
#' }
#' @export
define_table <- function(pop, table_name, ..., primary_key = NULL,
                         overwrite = FALSE) {

  # ---- Validate pop ----
  validate_tidybreed_pop(pop)

  # ---- Validate table_name ----
  validate_sql_identifier(table_name, what = "table name")

  if (table_name %in% SYSTEM_TABLES) {
    stop(
      "'", table_name, "' is a system-managed tidybreed table and cannot be ",
      "redefined with define_table(). System tables: ",
      paste(SYSTEM_TABLES, collapse = ", "),
      call. = FALSE
    )
  }

  # ---- Collect and validate column definitions ----
  extra_args <- list(...)

  if (length(extra_args) == 0L) {
    stop(
      "define_table() requires at least one column definition in `...`. ",
      "Use named typed-NA arguments, e.g. run_id = NA_integer_.",
      call. = FALSE
    )
  }

  col_names <- names(extra_args)
  if (is.null(col_names) || any(col_names == "")) {
    stop("All `...` arguments must be named (each name becomes a column).",
         call. = FALSE)
  }

  for (nm in col_names) {
    validate_sql_identifier(nm, what = "column name")
  }

  # ---- Infer DuckDB types ----
  col_types <- vapply(extra_args, infer_duckdb_type, character(1L))

  # ---- Validate primary_key ----
  if (!is.null(primary_key)) {
    if (!is.character(primary_key) || length(primary_key) != 1L) {
      stop("`primary_key` must be a single character string.", call. = FALSE)
    }
    if (!primary_key %in% col_names) {
      stop(
        "`primary_key` '", primary_key, "' is not among the defined columns: ",
        paste(col_names, collapse = ", "),
        call. = FALSE
      )
    }
  }

  # ---- Handle existing table ----
  already_exists <- table_name %in% pop$tables ||
    DBI::dbExistsTable(pop$db_conn, table_name)

  if (already_exists) {
    if (!overwrite) {
      stop(
        "Table '", table_name, "' already exists. ",
        "Use overwrite = TRUE to drop and recreate it.",
        call. = FALSE
      )
    }
    DBI::dbExecute(pop$db_conn,
                   paste0("DROP TABLE IF EXISTS ", table_name))
    pop$tables <- setdiff(pop$tables, table_name)
  }

  # ---- Build and execute CREATE TABLE ----
  col_defs <- mapply(
    function(nm, tp) {
      pk_clause <- if (!is.null(primary_key) && nm == primary_key) {
        " PRIMARY KEY"
      } else {
        ""
      }
      paste0(nm, " ", tp, pk_clause)
    },
    col_names, col_types,
    SIMPLIFY = TRUE, USE.NAMES = FALSE
  )

  create_sql <- paste0(
    "CREATE TABLE ", table_name,
    " (", paste(col_defs, collapse = ", "), ")"
  )
  DBI::dbExecute(pop$db_conn, create_sql)

  # ---- Register in pop ----
  pop$tables <- unique(c(pop$tables, table_name))

  invisible(pop)
}
