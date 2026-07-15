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
#' @description
#' Prints a grouped, human-readable summary of the population: a header rule
#' with the population name, the database location and connection status, then
#' one section each for the **genome** (loci, chromosomes, physical/genetic
#' length, founder haplotype pool), the **model** (genetic-component traits,
#' observed phenotypes, selection indices, QTL), the **individuals** (total,
#' broken down by sex and by line), and the **records** written so far
#' (phenotypes, TBVs, EBVs, index values).
#'
#' Sections are shown only when their underlying data exists, so a freshly
#' opened population collapses to just the header, the database line, and a
#' pointer to [schema()] / [describe_table()]. Counts are queried live from the
#' DuckDB connection; a closed connection prints as `[disconnected]` without
#' running any queries.
#'
#' @param x A tidybreed_pop object
#' @param ... Additional arguments (not used)
#' @return `x`, invisibly.
#' @export
print.tidybreed_pop <- function(x, ...) {
  width <- getOption("width", 80L)
  .tb_ruler(paste0("tidybreed population: ", x$pop_name), width)

  db_label <- if (identical(x$db_path, ":memory:")) "in-memory" else basename(x$db_path)

  # A dead connection: report and stop before touching the database.
  if (!DBI::dbIsValid(x$db_conn)) {
    cat("  Database   ", db_label, "  [disconnected]\n", sep = "")
    cat(strrep("─", width), "\n", sep = "")
    return(invisible(x))
  }

  cat("  Database   ", db_label, "  [connected]\n", sep = "")

  conn   <- x$db_conn
  tables <- x$tables
  # Guarded scalar COUNT/aggregate: a partial or unusual DB never errors print().
  n_of <- function(sql) {
    tryCatch(as.numeric(DBI::dbGetQuery(conn, sql)[[1L]]), error = function(e) 0)
  }
  fmt <- function(n) format(as.integer(round(n)), big.mark = ",")

  # ── Genome ────────────────────────────────────────────────────────────────
  n_loci <- if ("genome_meta" %in% tables) n_of("SELECT COUNT(*) FROM genome_meta") else 0
  if (n_loci > 0) {
    n_chr  <- n_of("SELECT COUNT(DISTINCT chr) FROM genome_meta")
    len_mb <- n_of("SELECT SUM(mx) FROM (SELECT MAX(pos_bp) mx FROM genome_meta GROUP BY chr)") / 1e6
    parts  <- c(paste(fmt(n_loci), "loci"), paste(fmt(n_chr), "chr"),
                paste0(fmt(len_mb), " Mb"))
    if ("genome_map" %in% tables) {
      # Genetic length = sum of per-chromosome cM spans on the default map.
      len_cm <- n_of(paste0(
        "SELECT SUM(span) FROM (",
        "SELECT MAX(pos_cM) - MIN(pos_cM) AS span FROM genome_map gm ",
        "JOIN genome_meta g ON g.locus_id = gm.locus_id ",
        "WHERE gm.map_name = 'default' GROUP BY g.chr)"
      ))
      if (len_cm > 0) parts <- c(parts, paste0(fmt(len_cm), " cM"))
    }
    cat("\n  Genome     ", paste(parts, collapse = " · "), "\n", sep = "")

    if ("founder_haplotypes" %in% tables) {
      # haplotype_id is unique only within a line; count (line_name, haplotype_id) pairs.
      n_hap <- n_of(paste0(
        "SELECT COUNT(*) FROM ",
        "(SELECT DISTINCT line_name, haplotype_id FROM founder_haplotypes)"
      ))
      if (n_hap > 0)
        cat("             founder pool: ", fmt(n_hap), " haplotypes\n", sep = "")
    }
  }

  # ── Model ───────────────────────────────────────────────────────────────────
  n_traits <- if ("trait_meta" %in% tables) n_of("SELECT COUNT(*) FROM trait_meta") else 0
  n_pheno  <- if ("phenotype_meta" %in% tables) n_of("SELECT COUNT(*) FROM phenotype_meta") else 0
  n_index  <- if ("index_meta" %in% tables)
    n_of("SELECT COUNT(DISTINCT index_name) FROM index_meta WHERE index_name IS NOT NULL") else 0
  n_qtl    <- if ("genome_effects" %in% tables)
    n_of("SELECT COUNT(DISTINCT locus_name) FROM genome_effects WHERE genome_effect_type = 'additive'") else 0
  model_parts <- c(
    if (n_traits > 0) paste(fmt(n_traits), if (n_traits == 1) "trait" else "traits"),
    if (n_pheno  > 0) paste(fmt(n_pheno),  if (n_pheno  == 1) "phenotype" else "phenotypes"),
    if (n_index  > 0) paste(fmt(n_index),  if (n_index  == 1) "index" else "indices"),
    if (n_qtl    > 0) paste(fmt(n_qtl), "QTL")
  )
  if (length(model_parts) > 0)
    cat("  Model      ", paste(model_parts, collapse = " · "), "\n", sep = "")

  # ── Individuals ─────────────────────────────────────────────────────────────
  n_ind <- if ("ind_meta" %in% tables) n_of("SELECT COUNT(*) FROM ind_meta") else 0
  if (n_ind > 0) {
    cat("\n  Individuals  ", fmt(n_ind), "\n", sep = "")

    by_sex <- tryCatch(
      DBI::dbGetQuery(conn,
        "SELECT sex, COUNT(*) AS n FROM ind_meta GROUP BY sex ORDER BY sex"),
      error = function(e) NULL)
    if (!is.null(by_sex) && nrow(by_sex) > 0)
      cat("    by sex     ",
          paste(paste(vapply(by_sex$n, fmt, character(1)), by_sex$sex), collapse = " · "),
          "\n", sep = "")

    by_line <- tryCatch(
      DBI::dbGetQuery(conn,
        "SELECT line_name, COUNT(*) AS n FROM ind_meta GROUP BY line_name ORDER BY line_name"),
      error = function(e) NULL)
    if (!is.null(by_line) && nrow(by_line) > 1)
      cat("    by line    ",
          paste(paste(by_line$line_name, vapply(by_line$n, fmt, character(1))), collapse = " · "),
          "\n", sep = "")
  }

  # ── Records ─────────────────────────────────────────────────────────────────
  rec_parts <- c(
    .tb_rec(conn, tables, "ind_phenotype", "phenotypes", fmt),
    .tb_rec(conn, tables, "ind_tbv",       "TBV",        fmt),
    .tb_rec(conn, tables, "ind_ebv",       "EBV",        fmt),
    .tb_rec(conn, tables, "ind_index",     "index",      fmt)
  )
  if (length(rec_parts) > 0)
    cat("\n  Records    ", paste(rec_parts, collapse = " · "), "\n", sep = "")

  # ── Tool dirs ───────────────────────────────────────────────────────────────
  tool_dirs <- x$run_dirs[names(x$run_dirs) != "base"]
  if (length(tool_dirs) > 0)
    cat("  Tool dirs  ", paste(names(tool_dirs), collapse = ", "), "\n", sep = "")

  cat("\n  schema(pop) · describe_table(pop, \"name\")\n")
  cat(strrep("─", width), "\n", sep = "")

  invisible(x)
}


# One "label N" record fragment, or NULL when the table is absent/empty.
.tb_rec <- function(conn, tables, tbl, label, fmt) {
  if (!tbl %in% tables) return(NULL)
  n <- tryCatch(as.numeric(DBI::dbGetQuery(conn, paste0("SELECT COUNT(*) FROM ", tbl))[[1L]]),
                error = function(e) 0)
  if (n > 0) paste(label, fmt(n)) else NULL
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

