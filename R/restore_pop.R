#' Reconnect to an existing tidybreed database
#'
#' @description
#' Re-opens a `.duckdb` file created by `initialize_genome()` and reconstructs
#' the `tidybreed_pop` R object so that simulation can continue. No metadata is
#' stored on the object; all genome statistics are derived from live DB queries
#' when needed by individual functions.
#'
#' This is the intended way to resume a simulation after `close_pop()`, after a
#' crash, or in a new R session. A common pattern is to call
#' `initialize_genome()` once, configure all traits and covariances, then start
#' each replicate with `restore_pop()` followed by `add_founders()`.
#'
#' @param db_path Character. Path to an existing `.duckdb` file created by
#'   `initialize_genome()`.
#' @param pop_name Character or `NULL`. Population name to assign to the
#'   restored object. When `NULL` (default) the name is inferred from the
#'   filename by stripping a trailing `_tidybreed.duckdb` suffix; if the file
#'   does not follow that convention the base filename (without extension) is
#'   used.
#'
#' @return A fully operational `tidybreed_pop` object.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # One-time setup
#' pop <- initialize_genome(
#'   pop_name     = "cattle",
#'   n_loci       = 50000,
#'   n_chr        = 29,
#'   chr_len_Mb   = 100,
#'   n_haplotypes = 200,
#'   db_path      = "cattle_sim.duckdb"
#' )
#' pop <- define_trait(pop, trait_name = "milk", ...)
#' close_pop(pop)
#'
#' # Each replicate
#' pop <- restore_pop("cattle_sim.duckdb")
#' pop <- add_founders(pop, n_males = 50, n_females = 200, line_name = "A")
#' }
restore_pop <- function(db_path, pop_name = NULL) {

  stopifnot(is.character(db_path), length(db_path) == 1)

  if (db_path == ":memory:") {
    stop(
      "restore_pop() cannot restore in-memory databases (':memory:'). ",
      "In-memory databases are not persisted to disk.",
      call. = FALSE
    )
  }

  if (!file.exists(db_path)) {
    stop("Database file '", db_path, "' does not exist.", call. = FALSE)
  }

  db_conn <- tryCatch(
    DBI::dbConnect(duckdb::duckdb(), dbdir = db_path),
    error = function(e) {
      stop(
        "Failed to open '", db_path, "' as a DuckDB database. ",
        "Is the file a valid DuckDB file and not open in another process?\n",
        "Original error: ", conditionMessage(e),
        call. = FALSE
      )
    }
  )

  existing_tables <- DBI::dbListTables(db_conn)

  if (!"genome_meta" %in% existing_tables) {
    DBI::dbDisconnect(db_conn, shutdown = TRUE)
    stop(
      "Table 'genome_meta' not found in '", db_path, "'. ",
      "Is this a valid tidybreed database created by initialize_genome()?",
      call. = FALSE
    )
  }

  n_loci <- DBI::dbGetQuery(db_conn, "SELECT COUNT(*) AS n FROM genome_meta")$n

  if (n_loci == 0L) {
    DBI::dbDisconnect(db_conn, shutdown = TRUE)
    stop("'genome_meta' is empty in '", db_path, "'.", call. = FALSE)
  }

  # Infer pop_name from filename when not supplied
  if (is.null(pop_name)) {
    base <- basename(db_path)
    stripped <- sub("_tidybreed\\.duckdb$", "", base, ignore.case = TRUE)
    pop_name <- if (stripped != base) stripped else tools::file_path_sans_ext(base)
  }

  pop <- new_tidybreed_pop(
    db_conn  = db_conn,
    pop_name = pop_name,
    db_path  = db_path,
    tables   = existing_tables
  )

  validate_tidybreed_pop(pop)

  if (!"_schema_meta" %in% existing_tables) {
    warning(
      "This database predates schema metadata (pre-v0.36.0). ",
      "Run `migrate_schema_meta(pop)` to add descriptions.",
      call. = FALSE
    )
  }

  n_ind <- DBI::dbGetQuery(db_conn, "SELECT COUNT(*) AS n FROM ind_meta")$n

  message("Restored population '", pop_name, "' from '", db_path, "'")
  message("  Loci: ", n_loci, ", Individuals: ", n_ind)

  pop
}
