#' Archive a completed replicate and reset the working database
#'
#' @description
#' Stamps per-replicate tables with a `replicate` column and appends them to a
#' separate archive DuckDB file, copies static metadata tables to the archive on
#' the **first call only** (no `replicate` column), then deletes the working-DB
#' rows so the next replicate starts clean. This is the only function that
#' writes to an archive DB; the working DB stays lean throughout the loop.
#'
#' @param pop A `tidybreed_pop` object with an open DuckDB connection.
#' @param replicate Integer scalar. Replicate number stamped on every archived
#'   row in `store_and_reset` tables. Defaults to
#'   `getOption("tidybreed.replicate")` (starts at `1L`). Always increments the
#'   option by 1 after a successful archive so the next loop iteration uses the
#'   next number automatically.
#' @param archive_path Character or `NULL`. Full path to the archive DuckDB
#'   file. When `NULL` (default) the path is resolved from global options — see
#'   Details. Passing `NULL` and leaving both archive options unset means only
#'   the reset phases run (no archive is written).
#' @param store_and_reset Character vector. Tables to archive with a `replicate`
#'   column and then DELETE from the working DB after archiving succeeds.
#' @param store_once Character vector. Tables copied to the archive on the first
#'   call only (detected by whether the table already exists in the archive
#'   schema). Working-DB rows are left intact. Use for static configuration that
#'   does not change between replicates (trait definitions, genome structure,
#'   variance components, index weights).
#' @param reset_only Character vector. Tables whose rows are deleted from the
#'   working DB without being archived. Intended for large genomic tables
#'   (`genome_haplotype`, `genome_genotype`) that are unnecessary for most
#'   downstream analyses.
#'
#' @details
#' **Archive path resolution** — first non-`NULL` result wins:
#'
#' 1. Explicit `archive_path` argument (full file path)
#' 2. `file.path(tidybreed.archive_path, tidybreed.db_name_archive)` when both
#'    options are set
#' 3. `file.path(dirname(pop$db_path), tidybreed.db_name_archive)` when only
#'    `tidybreed.db_name_archive` is set (archive lands next to the working DB)
#' 4. `tidybreed.db_name_archive = NULL` → `archive_path` resolves to `NULL` →
#'    skip all archive writes; only the reset phases run
#'
#' **Individual-ID invariant**: in sequential loops, sequence counters keep
#' accumulating so `id_ind` stays globally unique across replicates. In HPC
#' parallel jobs each job reuses `id_ind` values — use `(replicate, id_ind)` as
#' the composite key in all archive queries.
#'
#' **HPC note**: concurrent writes from multiple jobs to one shared archive file
#' are **not supported**. Use per-job archive files
#' (`scenario_rep_007.duckdb`) and merge in a post-processing step.
#'
#' @section Usage patterns:
#'
#' **In-process loop**:
#' ```r
#' options(tidybreed.replicate = 1L)
#' options(tidybreed.db_name_archive = "scenario_1_all_reps.duckdb")
#'
#' pop <- open_pop(...) |> define_genome(...) |>
#'   define_trait(...) |> define_phenotype(...)
#'
#' for (i in seq_len(100)) {
#'   pop <- pop |>
#'     get_table("founder_haplotypes") |>
#'     add_founders(n_males = 10, n_females = 90, line_name = "A") |>
#'     get_table("ind_meta") |>
#'     add_phenotype("ADG") |>
#'     archive_replicate()
#' }
#'
#' # post-hoc analysis
#' archive <- restore_pop("scenario_1_all_reps.duckdb")
#' get_table(archive, "ind_phenotype") |>
#'   dplyr::filter(replicate == 5L) |>
#'   dplyr::collect()
#' ```
#'
#' **HPC / SLURM array** (one job per replicate):
#' ```r
#' rep_id <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID"))
#' options(tidybreed.replicate = rep_id)
#' options(tidybreed.db_name_archive = paste0("scenario_1_rep_", rep_id, ".duckdb"))
#' pop <- open_pop(...) |> ... |> archive_replicate()
#' ```
#'
#' @return `pop` invisibly (pipeline-compatible).
#' @export
archive_replicate <- function(
    pop,
    replicate       = getOption("tidybreed.replicate"),
    archive_path    = NULL,
    store_and_reset = c("ind_meta", "ind_phenotype", "ind_tbv",
                        "ind_ebv", "ind_index", "ind_true_index"),
    store_once      = c("genome_meta", "genome_effects",
                        "trait_meta", "trait_effects", "trait_var_comp",
                        "phenotype_meta", "phenotype_components",
                        "phenotype_var_comp", "index_meta"),
    reset_only      = c("genome_haplotype", "genome_genotype")
) {

  # ── Phase 0: Validation ───────────────────────────────────────────────────

  if (!inherits(pop, "tidybreed_pop"))
    stop("'pop' must be a tidybreed_pop object.", call. = FALSE)

  if (!DBI::dbIsValid(pop$db_conn))
    stop("The DuckDB connection in 'pop' is closed. ",
         "Call restore_pop() to reopen.", call. = FALSE)

  # Coerce replicate to integer; accept e.g. 1.0 but reject NA or < 1
  replicate <- suppressWarnings(as.integer(replicate[[1L]]))
  if (is.na(replicate) || replicate < 1L)
    stop("'replicate' must be a single positive integer.", call. = FALSE)

  # No table may appear in more than one category
  all_tbls <- c(store_and_reset, store_once, reset_only)
  dup_tbls <- unique(all_tbls[duplicated(all_tbls)])
  if (length(dup_tbls) > 0L)
    stop("Table(s) appear in more than one category: ",
         paste(dup_tbls, collapse = ", "), call. = FALSE)

  # Resolve archive path from options / argument
  resolved_path <- .resolve_archive_path(pop, archive_path)

  # Query the live table list before any ATTACH (avoids cross-DB confusion)
  conn      <- pop$db_conn
  db_tables <- DBI::dbListTables(conn)

  # Filter table lists to those that actually exist; warn about skipped ones
  store_and_reset <- .filter_archive_tables(store_and_reset, db_tables, "store_and_reset")
  store_once      <- .filter_archive_tables(store_once,      db_tables, "store_once")
  reset_only      <- .filter_archive_tables(reset_only,      db_tables, "reset_only")

  # Collision guard: 'replicate' column must not pre-exist in any table we'll stamp
  if (!is.null(resolved_path) && length(store_and_reset) > 0L) {
    for (tbl in store_and_reset) {
      info <- DBI::dbGetQuery(conn, paste0('DESCRIBE "', tbl, '"'))
      if ("replicate" %in% info$column_name)
        stop("Table '", tbl, "' already contains a column named 'replicate'. ",
             "Rename it before calling archive_replicate().", call. = FALSE)
    }
  }

  # Ensure archive directory exists
  if (!is.null(resolved_path)) {
    arc_dir <- dirname(resolved_path)
    if (!dir.exists(arc_dir))
      dir.create(arc_dir, recursive = TRUE, showWarnings = FALSE)
    if (file.access(arc_dir, 2L) != 0L)
      stop("Archive directory is not writable: ", arc_dir, call. = FALSE)
  }

  # ── Phase 1: Archive writes ───────────────────────────────────────────────

  if (!is.null(resolved_path)) {

    # Detach any stale 'archive' DB left over from a prior failed call
    tryCatch(DBI::dbExecute(conn, "DETACH archive"), error = function(e) NULL)

    # ATTACH creates the archive DB file if it doesn't exist yet
    DBI::dbExecute(conn, paste0("ATTACH '", resolved_path, "' AS archive"))
    on.exit(
      tryCatch(DBI::dbExecute(conn, "DETACH archive"), error = function(e) NULL),
      add = TRUE
    )

    # store_and_reset: append rows stamped with replicate number
    sar_counts <- integer(length(store_and_reset))
    names(sar_counts) <- store_and_reset
    for (tbl in store_and_reset) {
      src_cols <- .ensure_archive_table(conn, tbl, add_replicate = TRUE)
      qcols    <- paste(paste0('"', src_cols, '"'), collapse = ", ")
      sar_counts[[tbl]] <- DBI::dbExecute(conn, paste0(
        'INSERT INTO archive."', tbl, '" (', qcols, ', "replicate") ',
        'SELECT ', qcols, ', ', replicate, ' AS replicate ',
        'FROM "', tbl, '"'
      ))
    }

    # store_once: copy to archive only if table is absent (first replicate)
    arc_existing  <- .list_archive_tables(conn)
    store_once_copied  <- character(0L)
    store_once_skipped <- character(0L)
    for (tbl in store_once) {
      if (tbl %in% arc_existing) {
        store_once_skipped <- c(store_once_skipped, tbl)
        next
      }
      src_cols <- .ensure_archive_table(conn, tbl, add_replicate = FALSE)
      qcols    <- paste(paste0('"', src_cols, '"'), collapse = ", ")
      DBI::dbExecute(conn, paste0(
        'INSERT INTO archive."', tbl, '" (', qcols, ') ',
        'SELECT ', qcols, ' FROM "', tbl, '"'
      ))
      store_once_copied <- c(store_once_copied, tbl)
    }

    # Provenance row in _archive_meta
    .write_archive_meta(conn, pop, replicate)

    # ── Messages: archive summary ─────────────────────────────────────────
    if (length(store_and_reset) > 0L) {
      sar_detail <- paste(
        paste0(names(sar_counts), " (", sar_counts, " rows)"),
        collapse = ", "
      )
      message("Archived replicate ", replicate,
              " [stamped with replicate column]: ", sar_detail)
    }
    if (length(store_once_copied) > 0L)
      message("Copied to archive (first time, no replicate column): ",
              paste(store_once_copied, collapse = ", "))
    if (length(store_once_skipped) > 0L)
      message("Already in archive, skipped (store_once): ",
              paste(store_once_skipped, collapse = ", "))
  }

  # ── Phase 2: Reset working DB ─────────────────────────────────────────────
  # Runs only after all archive writes complete (on.exit handles DETACH)

  reset_tbls <- c(store_and_reset, reset_only)
  for (tbl in reset_tbls) {
    DBI::dbExecute(conn, paste0('DELETE FROM "', tbl, '"'))
  }
  if (length(reset_tbls) > 0L)
    message("Wiped from working DB: ", paste(reset_tbls, collapse = ", "))

  # ── Phase 3: Increment counter ────────────────────────────────────────────

  options(tidybreed.replicate = replicate + 1L)
  message("tidybreed.replicate auto-incremented to ",
          getOption("tidybreed.replicate"),
          " (set options(tidybreed.replicate = N) to override)")

  invisible(pop)
}


# ── Internal helpers ──────────────────────────────────────────────────────────

#' Resolve the archive DB file path from function arg + global options
#'
#' @keywords internal
.resolve_archive_path <- function(pop, archive_path) {

  # 1. Explicit argument wins
  if (!is.null(archive_path)) return(archive_path)

  opt_dir  <- getOption("tidybreed.archive_path",    NULL)
  opt_name <- getOption("tidybreed.db_name_archive", NULL)

  # 4. No db_name_archive → skip archiving
  if (is.null(opt_name)) return(NULL)

  # 2. Both archive_path option and db_name_archive set
  if (!is.null(opt_dir)) return(file.path(opt_dir, opt_name))

  # 3. Only db_name_archive set → same directory as working DB
  if (identical(pop$db_path, ":memory:"))
    stop(
      "Cannot archive from an in-memory database without an explicit ",
      "'archive_path' argument or 'tidybreed.archive_path' option.",
      call. = FALSE
    )

  file.path(dirname(pop$db_path), opt_name)
}


#' Filter a table list to those present in the working DB; warn about the rest
#'
#' @keywords internal
.filter_archive_tables <- function(tbls, db_tables, category) {
  missing_tbls <- setdiff(tbls, db_tables)
  if (length(missing_tbls) > 0L)
    message("archive_replicate(): skipping ", length(missing_tbls),
            " table(s) not found in the working DB (",
            category, "): ",
            paste(missing_tbls, collapse = ", "))
  intersect(tbls, db_tables)
}


#' List tables in the currently-attached 'archive' DB
#'
#' @keywords internal
.list_archive_tables <- function(conn) {
  tryCatch(
    DBI::dbGetQuery(
      conn,
      "SELECT table_name FROM duckdb_tables() WHERE database_name = 'archive'"
    )$table_name,
    error = function(e) character(0L)
  )
}


#' Ensure an archive table exists and its schema matches the source.
#'
#' Creates the archive table from the source schema on first call. On
#' subsequent calls checks that the schemas are identical (strict policy —
#' errors on any drift, including user-added columns that appear only in one
#' DB). The `replicate` column added by archive_replicate() is excluded from
#' the comparison.
#'
#' @return Character vector of source column names (in schema order), invisibly.
#' @keywords internal
.ensure_archive_table <- function(conn, tbl_name, add_replicate) {

  # Source schema
  src_info  <- DBI::dbGetQuery(conn, paste0('DESCRIBE "', tbl_name, '"'))
  src_cols  <- src_info$column_name
  src_types <- stats::setNames(src_info$column_type, src_cols)

  # Does the archive table already exist?
  arc_info <- tryCatch(
    DBI::dbGetQuery(conn, paste0('DESCRIBE archive."', tbl_name, '"')),
    error = function(e) NULL
  )

  if (!is.null(arc_info)) {
    # --- Schema compatibility check ---
    arc_cols  <- arc_info$column_name
    arc_types <- stats::setNames(arc_info$column_type, arc_cols)

    # Exclude the 'replicate' stamp column from comparison (we added it)
    arc_cols_cmp  <- setdiff(arc_cols, "replicate")
    arc_types_cmp <- arc_types[arc_cols_cmp]

    new_in_src <- setdiff(src_cols,      arc_cols_cmp)
    new_in_arc <- setdiff(arc_cols_cmp,  src_cols)
    msgs <- character(0L)
    if (length(new_in_src) > 0L)
      msgs <- c(msgs, paste0("column(s) in source but not in archive: ",
                             paste(new_in_src, collapse = ", ")))
    if (length(new_in_arc) > 0L)
      msgs <- c(msgs, paste0("column(s) in archive but not in source: ",
                             paste(new_in_arc, collapse = ", ")))
    if (length(msgs) > 0L)
      stop("Schema mismatch for archive table '", tbl_name, "': ",
           paste(msgs, collapse = "; "),
           ".\nResolve the column difference before calling archive_replicate().",
           call. = FALSE)

    # Type check for shared columns
    for (col in src_cols) {
      if (!identical(src_types[[col]], arc_types_cmp[[col]]))
        stop("Type mismatch for column '", col, "' in table '", tbl_name,
             "': source has '", src_types[[col]],
             "', archive has '", arc_types_cmp[[col]], "'.",
             call. = FALSE)
    }

  } else {
    # --- Create archive table from source schema ---
    col_defs <- paste(
      paste0('"', src_cols, '" ', src_types),
      collapse = ", "
    )
    if (add_replicate) col_defs <- paste0(col_defs, ', "replicate" INTEGER')
    DBI::dbExecute(conn,
      paste0('CREATE TABLE archive."', tbl_name, '" (', col_defs, ')'))
  }

  invisible(src_cols)
}


#' Append a provenance row to _archive_meta in the archive DB
#'
#' @keywords internal
.write_archive_meta <- function(conn, pop, replicate) {

  pkg_ver <- tryCatch(
    as.character(utils::packageVersion("tidybreed")),
    error = function(e) "unknown"
  )

  # Create _archive_meta on first call
  arc_existing <- .list_archive_tables(conn)
  if (!("_archive_meta" %in% arc_existing)) {
    DBI::dbExecute(conn, paste0(
      'CREATE TABLE archive."_archive_meta" (',
      '"replicate" INTEGER, ',
      '"pop_name" VARCHAR, ',
      '"working_db_path" VARCHAR, ',
      '"tidybreed_version" VARCHAR, ',
      '"archived_at" TIMESTAMP',
      ')'
    ))
  }

  # Escape single quotes in string values
  esc <- function(x) gsub("'", "''", x, fixed = TRUE)

  DBI::dbExecute(conn, paste0(
    'INSERT INTO archive."_archive_meta" ',
    '("replicate","pop_name","working_db_path","tidybreed_version","archived_at") ',
    'VALUES (',
    replicate, ',',
    "'", esc(pop$pop_name), "',",
    "'", esc(pop$db_path),  "',",
    "'", esc(pkg_ver),      "',",
    'CURRENT_TIMESTAMP',
    ')'
  ))

  invisible(NULL)
}
