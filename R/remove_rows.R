# ---- Internal helpers -------------------------------------------------------

# Delete rows from a single table using an exact composite-key match.
# key_df: data frame whose columns are the row keys (from TABLE_ROW_KEYS).
# Returns the number of rows deleted (invisibly).
delete_exact_rows <- function(conn, table_name, key_df,
                              dry_run = FALSE, verbose = TRUE) {

  if (nrow(key_df) == 0L) return(invisible(0L))

  key_cols  <- names(key_df)
  temp_name <- paste0("_dr_", table_name)

  DBI::dbWriteTable(conn, temp_name, key_df, overwrite = TRUE)
  on.exit(DBI::dbExecute(conn, paste0("DROP TABLE IF EXISTS ", temp_name)),
          add = TRUE)

  join_sql <- paste(
    paste0("t.", key_cols, " = f.", key_cols),
    collapse = " AND "
  )

  if (dry_run) {
    n <- DBI::dbGetQuery(
      conn,
      paste0("SELECT COUNT(*) AS n FROM ", table_name, " AS t ",
             "INNER JOIN ", temp_name, " AS f ON ", join_sql)
    )$n
    if (verbose) {
      message("[dry_run] Would delete ", n, " row", if (n != 1L) "s" else "",
              " from `", table_name, "`")
    }
    return(invisible(n))
  }

  n <- DBI::dbExecute(
    conn,
    paste0("DELETE FROM ", table_name, " AS t ",
           "USING ", temp_name, " AS f WHERE ", join_sql)
  )

  if (verbose) {
    message("Deleted ", n, " row", if (n != 1L) "s" else "",
            " from `", table_name, "`")
  }
  invisible(n)
}


# Delete rows from a table matching a set of id_ind values.
# Used by cross-table mode; always keys on id_ind.
delete_by_id_ind <- function(conn, table_name, id_ind_vals,
                             dry_run = FALSE, verbose = TRUE) {

  if (length(id_ind_vals) == 0L) return(invisible(0L))

  key_df    <- data.frame(id_ind = id_ind_vals, stringsAsFactors = FALSE)
  temp_name <- paste0("_dr_", table_name)
  n_ids     <- length(id_ind_vals)

  DBI::dbWriteTable(conn, temp_name, key_df, overwrite = TRUE)
  on.exit(DBI::dbExecute(conn, paste0("DROP TABLE IF EXISTS ", temp_name)),
          add = TRUE)

  join_sql <- "t.id_ind = f.id_ind"

  if (dry_run) {
    n <- DBI::dbGetQuery(
      conn,
      paste0("SELECT COUNT(*) AS n FROM ", table_name, " AS t ",
             "INNER JOIN ", temp_name, " AS f ON ", join_sql)
    )$n
    if (verbose) {
      message("[dry_run] Would delete ", n, " row", if (n != 1L) "s" else "",
              " from `", table_name, "` (",
              n_ids, " id_ind value", if (n_ids != 1L) "s" else "", ")")
    }
    return(invisible(n))
  }

  n <- DBI::dbExecute(
    conn,
    paste0("DELETE FROM ", table_name, " AS t ",
           "USING ", temp_name, " AS f WHERE ", join_sql)
  )

  if (verbose) {
    message("Deleted ", n, " row", if (n != 1L) "s" else "",
            " from `", table_name, "` (",
            n_ids, " id_ind value", if (n_ids != 1L) "s" else "", ")")
  }
  invisible(n)
}


# ---- Exported function ------------------------------------------------------

#' Remove rows from one or more population tables
#'
#' Removes rows from the database based on a filtered `tidybreed_table`.
#' Chain after [get_table()] and [filter()] to target specific rows. Returns
#' the population invisibly — consistent with all other action functions.
#'
#' A filter is **required** by default. Calling `remove_rows()` without a
#' prior [filter()] will stop with an error explaining how to opt in to
#' full-table deletion via `confirm_all = TRUE`.
#'
#' @param tbl A `tidybreed_table` from [get_table()], with a [filter()] applied.
#' @param tables `NULL` (default) to delete from the current table only;
#'   `"all"` to delete from every `ind_*` table plus `genome_haplotype` and
#'   `genome_genotype` that exist in the population; or a character vector of
#'   specific table names (each must have an `id_ind` column). When not `NULL`,
#'   the source table must also have an `id_ind` column.
#' @param confirm_all Logical. Set to `TRUE` to allow deletion of all rows when
#'   no filter has been applied. This is a deliberate safeguard — you must
#'   explicitly opt in to wiping an entire table. Default `FALSE`.
#' @param dry_run Logical. If `TRUE`, report what would be deleted without
#'   modifying the database. Default `FALSE`.
#' @param verbose Logical. If `TRUE`, print a message after each deletion
#'   reporting the row count removed. Default `TRUE`.
#'
#' @return The `tidybreed_pop` object (invisibly).
#'
#' @details
#' **Single-table mode** (`tables = NULL`): uses the composite row key from
#' the internal `TABLE_ROW_KEYS` registry to delete exactly the rows matched
#' by the filter — no more, no less. For example, filtering `ind_tbv` by
#' `trait_name == "ADG"` deletes only the ADG rows, not all TBV rows for those
#' animals.
#'
#' **Cross-table mode** (`tables != NULL`): extracts unique `id_ind` values
#' from the filtered table and issues a `DELETE ... WHERE id_ind IN (...)`
#' for each target table via a temp-table JOIN. `tables = "all"` targets every
#' `ind_*` table plus `genome_haplotype`, `genome_genotype`, and
#' `ind_true_index` that currently exist in the population (including
#' `ind_meta`).
#'
#' A temporary DuckDB table is used for both modes to avoid SQL injection risk
#' and to handle large ID sets efficiently.
#'
#' @examples
#' \dontrun{
#' # Delete specific phenotype records for culled animals
#' pop |>
#'   get_table("ind_phenotype") |>
#'   dplyr::filter(id_ind %in% culled_ids, trait_name == "litter_size") |>
#'   remove_rows()
#'
#' # Remove all data for culled animals across every individual table
#' pop |>
#'   get_table("ind_meta") |>
#'   dplyr::filter(id_ind %in% culled_ids) |>
#'   remove_rows(tables = "all")
#'
#' # Delete from an explicit subset of tables only
#' pop |>
#'   get_table("ind_meta") |>
#'   dplyr::filter(id_ind %in% culled_ids) |>
#'   remove_rows(tables = c("ind_phenotype", "ind_tbv"))
#'
#' # Preview before deleting
#' pop |>
#'   get_table("ind_meta") |>
#'   dplyr::filter(id_ind %in% culled_ids) |>
#'   remove_rows(tables = "all", dry_run = TRUE)
#'
#' # Wipe an entire table (requires confirm_all = TRUE)
#' pop |>
#'   get_table("ind_phenotype") |>
#'   remove_rows(confirm_all = TRUE)
#' }
#'
#' @export
remove_rows <- function(tbl, tables = NULL, confirm_all = FALSE,
                        dry_run = FALSE, verbose = TRUE) {

  # ---- A. Input validation --------------------------------------------------
  if (!inherits(tbl, "tidybreed_table")) {
    stop("'tbl' must be a tidybreed_table from get_table().", call. = FALSE)
  }

  pop        <- tbl$pop
  validate_tidybreed_pop(pop)
  conn       <- pop$db_conn
  table_name <- tbl$table_name

  # ---- B. Guard: require a filter unless confirm_all = TRUE -----------------
  if (length(tbl$pending_filter) == 0L && !isTRUE(confirm_all)) {
    n_total <- DBI::dbGetQuery(
      conn, paste0("SELECT COUNT(*) AS n FROM ", table_name)
    )$n
    stop(
      "remove_rows(): no filter applied — this would delete ALL ",
      n_total, " row", if (n_total != 1L) "s" else "",
      " from '", table_name, "'.\n",
      "  To delete an entire table, pass confirm_all = TRUE:\n",
      "    get_table(\"", table_name, "\") |> remove_rows(confirm_all = TRUE)\n",
      "  To preview what a filter would delete, use dry_run = TRUE:\n",
      "    get_table(\"", table_name, "\") |> filter(...) |> remove_rows(dry_run = TRUE)",
      call. = FALSE
    )
  }

  # ---- C. Collect filtered rows ---------------------------------------------
  collected_df <- dplyr::collect(tbl$tbl)

  if (nrow(collected_df) == 0L) {
    warning("filter() matched 0 rows in '", table_name,
            "'. Nothing to delete.", call. = FALSE)
    return(invisible(pop))
  }

  # ---- D. Single-table mode -------------------------------------------------
  if (is.null(tables)) {

    key_cols <- TABLE_ROW_KEYS[[table_name]]
    if (is.null(key_cols)) {
      stop("Cannot delete from '", table_name,
           "': table is not registered in TABLE_ROW_KEYS. ",
           "Use tables = c('", table_name, "') or tables = 'all' if ",
           "this table has an id_ind column.",
           call. = FALSE)
    }

    missing_cols <- setdiff(key_cols, names(collected_df))
    if (length(missing_cols) > 0L) {
      stop("Cannot delete from '", table_name,
           "': collected data is missing key column(s): ",
           paste(missing_cols, collapse = ", "),
           call. = FALSE)
    }

    key_df <- collected_df[, key_cols, drop = FALSE]
    delete_exact_rows(conn, table_name, key_df, dry_run, verbose)

    return(invisible(pop))
  }

  # ---- E. Cross-table mode --------------------------------------------------
  if (!"id_ind" %in% names(collected_df)) {
    stop("remove_rows() with tables != NULL requires the source table '",
         table_name, "' to have an 'id_ind' column.",
         call. = FALSE)
  }

  id_ind_vals <- unique(collected_df[["id_ind"]])

  if (identical(tables, "all")) {
    db_tables   <- DBI::dbListTables(conn)
    ind_tbls    <- db_tables[grepl("^ind_", db_tables)]
    candidates  <- c(ind_tbls, "genome_haplotype", "genome_genotype")
    target_tables <- intersect(IND_TABLE_ID_IND_COLS,
                               intersect(candidates, db_tables))
  } else {
    if (!is.character(tables) || length(tables) == 0L) {
      stop("'tables' must be NULL, \"all\", or a non-empty character vector.",
           call. = FALSE)
    }
    db_tables  <- DBI::dbListTables(conn)
    bad_exist  <- setdiff(tables, db_tables)
    if (length(bad_exist) > 0L) {
      stop("The following tables do not exist in this population: ",
           paste(bad_exist, collapse = ", "),
           call. = FALSE)
    }
    bad_known  <- setdiff(tables, IND_TABLE_ID_IND_COLS)
    if (length(bad_known) > 0L) {
      stop("The following tables do not have an 'id_ind' column: ",
           paste(bad_known, collapse = ", "),
           call. = FALSE)
    }
    target_tables <- tables
  }

  for (tgt in target_tables) {
    delete_by_id_ind(conn, tgt, id_ind_vals, dry_run, verbose)
  }

  invisible(pop)
}
