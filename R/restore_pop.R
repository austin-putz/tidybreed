#' Reconnect to an existing tidybreed database
#'
#' @description
#' Re-opens a `.duckdb` file and reconstructs the `tidybreed_pop` R object so
#' that simulation can continue. No metadata is stored on the object; all
#' genome statistics are derived from live DB queries when needed.
#'
#' This is the intended way to resume a simulation after [close_pop()], after a
#' crash, or in a new R session. A common pattern is to call [open_pop()] once,
#' configure all traits and covariances, then start each replicate with
#' `restore_pop()` followed by [add_founders()].
#'
#' `run_dirs` are reconstructed from the database location and the `tools`
#' argument. Set `options(tidybreed.tools = ...)` before calling `restore_pop()`
#' to get the same tool dirs back automatically.
#'
#' @param db_path Character. Path to an existing `.duckdb` file.
#' @param pop_name Character or `NULL`. Population name to assign to the
#'   restored object. When `NULL` (default) the name is inferred from the
#'   filename by stripping a trailing `_tidybreed.duckdb` or `.duckdb` suffix.
#' @param tools Character vector or `NULL`. Tool subdirectory names to
#'   reconstruct in `pop$run_dirs`. Defaults to
#'   `getOption("tidybreed.tools", NULL)`.
#'
#' @return A fully operational `tidybreed_pop` object.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # One-time setup. db_name is placed inside open_pop()'s layer-2/3 folder
#' # structure (see ?open_pop), so capture pop$db_path rather than assuming
#' # the file lives at "cattle_sim.duckdb" in the working directory.
#' pop <- open_pop(pop_name = "cattle", db_name = "cattle_sim.duckdb") |>
#'   define_genome(n_loci = 50000, n_chr = 29, chr_len_Mb = 100) |>
#'   define_founder_haplotypes(n_haplotypes = 100, line_name = "A")
#' pop <- define_trait(pop, trait_name = "milk", target_add_var = 10)
#' db_path <- pop$db_path
#' close_pop(pop)
#'
#' # Each replicate
#' pop <- restore_pop(db_path)
#' pop <- pop |>
#'   get_table("founder_haplotypes") |>
#'   add_founders(n_males = 50, n_females = 200, line_name = "A")
#' }
restore_pop <- function(db_path,
                        pop_name = NULL,
                        tools    = getOption("tidybreed.tools", NULL)) {

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

  # Migrate old table names before reading table list
  .migrate_var_comp_tables(db_conn)

  existing_tables <- DBI::dbListTables(db_conn)

  # ind_meta is always created by open_pop() — use it as the minimum check
  if (!"ind_meta" %in% existing_tables) {
    DBI::dbDisconnect(db_conn, shutdown = TRUE)
    stop(
      "Table 'ind_meta' not found in '", db_path, "'. ",
      "Is this a valid tidybreed database created by open_pop()?",
      call. = FALSE
    )
  }

  # genome_meta is optional (pop may not have called define_genome() yet)
  has_genome <- "genome_meta" %in% existing_tables
  if (has_genome) {
    n_loci <- DBI::dbGetQuery(db_conn, "SELECT COUNT(*) AS n FROM genome_meta")$n
    if (n_loci == 0L) {
      warning("'genome_meta' is empty in '", db_path, "'.", call. = FALSE)
      n_loci <- 0L
    }
  } else {
    n_loci <- NA_integer_
    warning(
      "No 'genome_meta' table found in '", db_path, "'. ",
      "Call define_genome() to add genome structure.",
      call. = FALSE
    )
  }

  # Infer pop_name from filename when not supplied
  if (is.null(pop_name)) {
    base     <- basename(db_path)
    # Strip _tidybreed.duckdb or .duckdb suffix
    stripped <- sub("_tidybreed\\.duckdb$", "", base, ignore.case = TRUE)
    stripped <- sub("\\.duckdb$", "", stripped, ignore.case = TRUE)
    pop_name <- if (stripped != base) stripped else tools::file_path_sans_ext(base)
  }

  # Reconstruct run_dirs from db location + tools
  if (!is.null(tools) && length(tools) > 0) {
    db_dir   <- dirname(normalizePath(db_path))
    run_dirs <- setNames(file.path(db_dir, tools), tools)
    run_dirs <- c(base = db_dir, run_dirs)
  } else {
    run_dirs <- character(0)
  }

  pop <- new_tidybreed_pop(
    db_conn  = db_conn,
    pop_name = pop_name,
    db_path  = db_path,
    tables   = existing_tables,
    run_dirs = run_dirs
  )

  validate_tidybreed_pop(pop)

  if (!"_schema_meta" %in% existing_tables) {
    warning(
      "This database predates schema metadata (pre-v0.36.0); ",
      "table/column descriptions are unavailable.",
      call. = FALSE
    )
  }

  n_ind <- DBI::dbGetQuery(db_conn, "SELECT COUNT(*) AS n FROM ind_meta")$n

  message("Restored population '", pop_name, "' from '", db_path, "'")
  if (!is.na(n_loci))
    message("  Loci: ", n_loci, ", Individuals: ", n_ind)
  else
    message("  Individuals: ", n_ind)

  pop
}
