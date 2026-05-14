#' Reconnect to an existing tidybreed database
#'
#' @description
#' Re-opens a `.duckdb` file created by `initialize_genome()` and reconstructs
#' the `tidybreed_pop` R object so that simulation can continue. All runtime
#' metadata (`n_loci`, `n_chr`, `chr_len_Mb`, `chr_names`, `n_individuals`) is
#' derived from the current state of `genome_meta` and `ind_meta` — no
#' additional tables are required and no metadata is stored in the database.
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
#' @details
#' **What is derived and how:**
#' - `n_loci` — `COUNT(*)` of `genome_meta`
#' - `n_chr` — count of distinct `chr` values in `genome_meta`
#' - `chr_names` — ordered distinct `chr_name` values from `genome_meta`
#' - `chr_len_Mb` — `MAX(pos_Mb)` per chromosome from `genome_meta`. This is
#'   functionally identical to the original `chr_len_Mb` parameter for
#'   recombination simulation: crossovers beyond the last locus position cannot
#'   affect any observed allele.
#' - `n_individuals` — `COUNT(*)` of `ind_meta` (0 if no founders added yet)
#'
#' Founder-haplotype generation parameters (`n_haplotypes`, `allele_freq_dist`,
#' etc.) are not restored — they are only used at initialisation time.
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
#' pop <- add_trait(pop, trait_name = "milk", ...)
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

  # Derive per-chromosome lengths and names from current genome_meta state
  gm <- DBI::dbGetQuery(
    db_conn,
    "SELECT chr, chr_name, MAX(pos_Mb) AS chr_len_Mb
     FROM genome_meta
     GROUP BY chr, chr_name
     ORDER BY chr"
  )

  if (nrow(gm) == 0) {
    DBI::dbDisconnect(db_conn, shutdown = TRUE)
    stop("'genome_meta' is empty in '", db_path, "'.", call. = FALSE)
  }

  n_loci     <- as.integer(
    DBI::dbGetQuery(db_conn, "SELECT COUNT(*) AS n FROM genome_meta")$n
  )
  n_chr      <- nrow(gm)
  chr_names  <- gm$chr_name
  chr_len_Mb <- gm$chr_len_Mb

  # n_individuals: always an integer (0L if no founders yet). Both add_founders()
  # and add_offspring() use is.null() to detect first-use — 0L hits the else
  # branch and computes 0L + n = n, which is correct.
  n_individuals <- 0L
  if ("ind_meta" %in% existing_tables) {
    n_individuals <- as.integer(
      DBI::dbGetQuery(db_conn, "SELECT COUNT(*) AS n FROM ind_meta")$n
    )
  }

  # Infer pop_name from filename when not supplied
  if (is.null(pop_name)) {
    base <- basename(db_path)
    stripped <- sub("_tidybreed\\.duckdb$", "", base, ignore.case = TRUE)
    pop_name <- if (stripped != base) stripped else tools::file_path_sans_ext(base)
  }

  metadata <- list(
    n_loci        = n_loci,
    n_chr         = as.integer(n_chr),
    chr_len_Mb    = chr_len_Mb,
    chr_names     = chr_names,
    n_individuals = n_individuals
  )

  pop <- new_tidybreed_pop(
    db_conn  = db_conn,
    pop_name = pop_name,
    db_path  = db_path,
    tables   = existing_tables,
    metadata = metadata
  )

  validate_tidybreed_pop(pop)

  message("Restored population '", pop_name, "' from '", db_path, "'")
  message(
    "  Loci: ", n_loci,
    ", Chromosomes: ", n_chr,
    ", Individuals: ", n_individuals
  )

  pop
}
