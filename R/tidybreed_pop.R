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
#' @param metadata List of additional metadata
#'
#' @return A `tidybreed_pop` object
#' @keywords internal
new_tidybreed_pop <- function(db_conn,
                               pop_name,
                               db_path,
                               tables = character(),
                               metadata = list(),
                               pending_filter = list()) {

  stopifnot(inherits(db_conn, "duckdb_connection"))
  stopifnot(is.character(pop_name))
  stopifnot(is.character(db_path))

  structure(
    list(
      db_conn = db_conn,
      pop_name = pop_name,
      db_path = db_path,
      tables = tables,
      metadata = metadata,
      pending_filter = pending_filter
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
#' @export
print.tidybreed_pop <- function(x, ...) {
  cat("<tidybreed population>\n")
  cat("Population:", x$pop_name, "\n")
  cat("Database:", x$db_path, "\n")
  cat("Tables:", paste(x$tables, collapse = ", "), "\n")

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

  # Show pending filter if set
  if (length(x$pending_filter) > 0) {
    exprs <- vapply(x$pending_filter, rlang::as_label, character(1))
    cat("Pending filter:", paste(exprs, collapse = " & "), "\n")
  }

  invisible(x)
}


#' Get table from tidybreed population
#'
#' @description
#' Retrieve a table from the population database as a lazy dplyr tibble.
#' This allows for efficient querying without loading the entire table into memory.
#'
#' @param pop A tidybreed_pop object
#' @param table_name Name of the table to retrieve
#'
#' @return A lazy tibble (tbl_duckdb_connection)
#' @export
#'
#' @examples
#' \dontrun{
#' pop <- initialize_genome(pop_name = "A", n_loci = 100, n_chr = 2)
#' genome <- get_table(pop, "genome_meta")
#' genome %>% dplyr::filter(chr == 1)
#' }
get_table <- function(pop, table_name) {

  stopifnot(inherits(pop, "tidybreed_pop"))

  if (!table_name %in% pop$tables) {
    stop(
      "Table '", table_name, "' not found. ",
      "Available tables: ", paste(pop$tables, collapse = ", "),
      call. = FALSE
    )
  }

  dplyr::tbl(pop$db_conn, table_name)
}


#' Close tidybreed population database connection
#'
#' @param pop A tidybreed_pop object
#' @export
close_pop <- function(pop) {
  stopifnot(inherits(pop, "tidybreed_pop"))
  DBI::dbDisconnect(pop$db_conn, shutdown = TRUE)
  invisible(NULL)
}


#' Stash a dplyr filter predicate on a tidybreed population
#'
#' @description
#' S3 method for [dplyr::filter()] on `tidybreed_pop` objects. Captures the
#' predicate expressions and stores them on the population in `$pending_filter`,
#' to be applied by the next consuming operation (currently [add_phenotype()]
#' and [add_tbv()]). Multiple `filter()` calls stack with AND semantics.
#'
#' Non-consuming operations ignore the stashed filter. If an operation that
#' does not consume the filter runs while one is set, the filter remains on
#' the object and will be applied by the next consumer.
#'
#' @param .data A `tidybreed_pop` object.
#' @param ... Unquoted predicate expressions evaluated against `ind_meta`.
#' @param .preserve Present for S3 signature compatibility; not used.
#'
#' @return The `tidybreed_pop` object with the expressions appended to
#'   `$pending_filter`.
#'
#' @examples
#' \dontrun{
#' pop |>
#'   dplyr::filter(sex == "F", gen == 1L) |>
#'   add_phenotype("ADG")
#' }
#' @export
filter.tidybreed_pop <- function(.data, ..., .preserve = FALSE) {
  new_quosures <- rlang::enquos(...)
  .data$pending_filter <- c(.data$pending_filter, new_quosures)
  .data
}


#' Resolve and clear pending filter on a population
#'
#' @description
#' Applies the stacked predicates in `pop$pending_filter` to `ind_meta`, pulls
#' the matching `id_ind` values, and returns the ids together with a copy of
#' `pop` with the pending filter cleared. Used internally by [add_phenotype()]
#' and [add_tbv()].
#'
#' @param pop A `tidybreed_pop` object.
#' @return A list with components `ids` (character vector of matching `id_ind`)
#'   and `pop` (the population with `pending_filter` cleared). If no filter is
#'   set, `ids` is `NULL` and `pop` is returned unchanged.
#' @keywords internal
resolve_pending_filter <- function(pop) {

  if (length(pop$pending_filter) == 0) {
    return(list(ids = NULL, pop = pop))
  }

  if (!"ind_meta" %in% pop$tables) {
    stop(
      "Pending filter cannot be applied: ind_meta does not exist yet.",
      call. = FALSE
    )
  }

  filtered <- dplyr::tbl(pop$db_conn, "ind_meta")
  filtered <- dplyr::filter(filtered, !!!pop$pending_filter)
  ids <- dplyr::pull(dplyr::collect(dplyr::select(filtered, "id_ind")), "id_ind")

  pop$pending_filter <- list()
  list(ids = ids, pop = pop)
}
