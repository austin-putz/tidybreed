#' Create a new tidybreed population object
#'
#' @description
#' Constructor for the `tidybreed_pop` S3 class. This class wraps a DuckDB
#' connection and provides a tidy interface for breeding program simulation.
#'
#' @param db_conn A DuckDB connection object
#' @param pop_name Character string naming the population
#' @param db_path Path to the DuckDB database file
#' @param tables Character vector of table names in the database
#' @param run_dirs Named character vector mapping tool names to their pre-created
#'   layer-4 output directories. `character(0)` for in-memory databases or when
#'   no tools were declared. Always includes a `"base"` entry pointing to the
#'   layer-3 scenario directory when non-empty.
#'
#' @return A `tidybreed_pop` object
#' @keywords internal
new_tidybreed_pop <- function(db_conn,
                               pop_name,
                               db_path,
                               tables   = character(),
                               run_dirs = character(0)) {

  stopifnot(inherits(db_conn, "duckdb_connection"))
  stopifnot(is.character(pop_name))
  stopifnot(is.character(db_path))
  stopifnot(is.character(run_dirs))

  structure(
    list(
      db_conn  = db_conn,
      pop_name = pop_name,
      db_path  = db_path,
      tables   = tables,
      run_dirs = run_dirs
    ),
    class = c("tidybreed_pop", "list")
  )
}


#' Validate tidybreed population object
#'
#' @param x A tidybreed_pop object to validate
#' @return The object (invisibly) if valid, error otherwise
#' @keywords internal
validate_tidybreed_pop <- function(x) {

  # Check required components
  if (!inherits(x$db_conn, "duckdb_connection")) {
    stop("db_conn must be a DuckDB connection", call. = FALSE)
  }

  # Check connection is still valid
  if (!DBI::dbIsValid(x$db_conn)) {
    stop("Database connection is no longer valid", call. = FALSE)
  }

  # Check pop_name
  if (!is.character(x$pop_name) || length(x$pop_name) != 1) {
    stop("pop_name must be a single character string", call. = FALSE)
  }

  # Check tables exist in database
  existing_tables <- DBI::dbListTables(x$db_conn)
  missing_tables <- setdiff(x$tables, existing_tables)
  if (length(missing_tables) > 0) {
    warning(
      "Tables listed in object not found in database: ",
      paste(missing_tables, collapse = ", "),
      call. = FALSE
    )
  }

  invisible(x)
}


#' Print method for tidybreed_pop
#'
#' @param x A tidybreed_pop object
#' @param ... Additional arguments (not used)
#' @return `x`, invisibly.
#' @export
print.tidybreed_pop <- function(x, ...) {
  cat("<tidybreed population>\n")
  cat("Population:", x$pop_name, "\n")
  cat("Database:", x$db_path, "\n")
  user_tables <- setdiff(x$tables, "_schema_meta")
  cat("Tables:", paste(user_tables, collapse = ", "), "\n")

  # Show number of individuals if ind_meta exists
  if ("ind_meta" %in% x$tables) {
    n_ind <- DBI::dbGetQuery(
      x$db_conn,
      "SELECT COUNT(*) as n FROM ind_meta"
    )$n
    cat("Individuals:", n_ind, "\n")
  }

  # Show number of loci if genome_meta exists
  if ("genome_meta" %in% x$tables) {
    n_loci <- DBI::dbGetQuery(
      x$db_conn,
      "SELECT COUNT(*) as n FROM genome_meta"
    )$n
    cat("Loci:", n_loci, "\n")
  }

  # Show number of traits if trait_meta exists
  if ("trait_meta" %in% x$tables) {
    n_traits <- DBI::dbGetQuery(
      x$db_conn,
      "SELECT COUNT(*) as n FROM trait_meta"
    )$n
    cat("Traits:", n_traits, "\n")
  }

  # Show tool dirs if any are registered
  tool_dirs <- x$run_dirs[names(x$run_dirs) != "base"]
  if (length(tool_dirs) > 0)
    cat("Tool dirs:", paste(names(tool_dirs), collapse = ", "), "\n")

  cat("Use schema(pop) or describe_table(pop, \"name\") to explore table structure.\n")

  invisible(x)
}


#' Get a table reference from a tidybreed population
#'
#' @description
#' Returns a `tidybreed_table` object that can be piped into [dplyr::filter()] and
#' [mutate_table()], or collected with [dplyr::collect()].  All existing
#' `get_table(pop, x) |> dplyr::filter(...) |> dplyr::collect()` patterns
#' continue to work unchanged.
#'
#' @param pop A `tidybreed_pop` object.
#' @param table_name Name of the table to retrieve.
#'
#' @return A `tidybreed_table` object.
#' @export
#'
#' @examples
#' \dontrun{
#' pop <- open_pop(pop_name = "A", db_name = ":memory:") |>
#'   define_genome(n_loci = 100, n_chr = 2, chr_len_Mb = 100)
#'
#' # Read-only query
#' get_table(pop, "genome_meta") |> dplyr::collect()
#'
#' # Mutate all rows
#' pop <- pop |> get_table("ind_meta") |> mutate_table(gen = 1L)
#'
#' # Mutate a filtered subset
#' pop <- pop |>
#'   get_table("ind_meta") |>
#'   filter(sex == "M") |>
#'   mutate_table(gen = 2L)
#' }
get_table <- function(pop, table_name) {

  stopifnot(inherits(pop, "tidybreed_pop"))

  if (!table_name %in% pop$tables) {
    stop(
      "Table '", table_name, "' does not exist in this population. ",
      "Available tables: ", paste(pop$tables, collapse = ", "),
      call. = FALSE
    )
  }

  new_tidybreed_table(pop, table_name)
}


#' Close tidybreed population database connection
#'
#' @param pop A tidybreed_pop object
#' @param results_dir Optional path to a directory. When provided, the `.duckdb`
#'   file is moved there after the connection is closed. The directory is created
#'   if it does not yet exist. Errors if the destination file already exists.
#'   Ignored silently for in-memory databases.
#' @return `NULL`, invisibly.
#' @export
close_pop <- function(pop, results_dir = NULL) {

  stopifnot(inherits(pop, "tidybreed_pop"))

  db_path <- pop$db_path

  DBI::dbDisconnect(pop$db_conn, shutdown = TRUE)

  if (!is.null(results_dir)) {
    if (db_path == ":memory:") {
      warning("Cannot move an in-memory database; results_dir ignored.", call. = FALSE)
      return(invisible(NULL))
    }

    if (!file.exists(db_path)) {
      warning("Database file '", db_path, "' not found; nothing moved.", call. = FALSE)
      return(invisible(NULL))
    }

    if (!dir.exists(results_dir)) {
      dir.create(results_dir, recursive = TRUE)
    }

    dest_path <- file.path(results_dir, basename(db_path))

    if (file.exists(dest_path)) {
      stop("Destination already exists: ", dest_path, call. = FALSE)
    }

    moved <- file.rename(db_path, dest_path)
    if (!moved) {
      if (!file.copy(db_path, dest_path)) {
        stop("Failed to move database to '", dest_path, "'", call. = FALSE)
      }
      file.remove(db_path)
    }
  }

  invisible(NULL)
}

