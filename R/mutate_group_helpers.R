#' Shared internal helpers for mutate_group_*() functions
#'
#' @keywords internal
#' @name mutate_group_helpers
NULL


#' Validate a tidybreed_table, col_name, and overwrite state before any RNG.
#'
#' Returns a named list with all needed state, or NULL when the filter matched
#' zero rows (callers must propagate the early return).
#'
#' @param tbl_obj A `tidybreed_table` from `get_table()`.
#' @param col_name Character scalar. Target column name.
#' @param output_type `"INTEGER"` or `"VARCHAR"`.
#' @param overwrite Logical. If `FALSE`, error when filtered rows already have
#'   non-NULL values in `col_name`.
#' @return Named list or NULL.
#' @keywords internal
.resolve_group_target <- function(tbl_obj, col_name, output_type, overwrite) {

  stopifnot(inherits(tbl_obj, "tidybreed_table"))

  table_name <- tbl_obj$table_name
  pop        <- tbl_obj$pop
  conn       <- pop$db_conn

  reserved <- TABLE_RESERVED_COLS[[table_name]]
  if (is.null(reserved)) reserved <- character(0)
  validate_sql_identifier(col_name, what = "col_name", reserved = reserved)

  pk_col <- TABLE_PRIMARY_KEYS[[table_name]]
  if (is.null(pk_col)) {
    stop(
      "Table '", table_name, "' has no registered primary key. ",
      "mutate_group_*() requires a table with a known primary key.",
      call. = FALSE
    )
  }

  # Collect filtered PKs — sort for stable, reproducible ordering
  pks     <- sort(dplyr::collect(dplyr::select(tbl_obj$tbl, dplyr::all_of(pk_col)))[[pk_col]])
  n_total <- length(pks)

  if (n_total == 0L) {
    warning(
      "filter() matched 0 rows in '", table_name, "'. ",
      "No group values assigned.",
      call. = FALSE
    )
    return(NULL)
  }

  existing_cols <- DBI::dbListFields(conn, table_name)
  col_exists    <- col_name %in% existing_cols

  if (col_exists) {
    existing_type <- toupper(DBI::dbGetQuery(
      conn,
      paste0(
        "SELECT data_type FROM information_schema.columns ",
        "WHERE table_name = '", table_name, "' ",
        "AND column_name = '", col_name, "'"
      )
    )$data_type)

    type_ok <- switch(output_type,
      "INTEGER" = existing_type %in% c("INTEGER", "BIGINT", "INT", "INT4", "INT8",
                                        "HUGEINT", "SMALLINT", "TINYINT"),
      "VARCHAR" = existing_type %in% c("VARCHAR"),
      FALSE
    )
    if (!type_ok) {
      stop(
        "Column '", col_name, "' already exists in '", table_name, "' ",
        "with type '", existing_type, "', which is incompatible with ",
        "the expected '", output_type, "' type.",
        call. = FALSE
      )
    }

    if (!overwrite) {
      n_existing <- count_nonnull_in_filter(conn, table_name, col_name, pk_col, pks)
      if (n_existing > 0L) {
        stop(
          n_existing, " row", if (n_existing != 1L) "s" else "",
          " in the filtered subset already have non-NULL values for '",
          col_name, "'. Set overwrite = TRUE to replace them.",
          call. = FALSE
        )
      }
    }
  }

  list(
    pop        = pop,
    conn       = conn,
    table_name = table_name,
    pk_col     = pk_col,
    pks        = pks,
    n_total    = n_total,
    col_exists = col_exists
  )
}


#' Create a column in a table if it does not already exist.
#'
#' @keywords internal
.ensure_col_exists <- function(conn, table_name, col_name, output_type) {
  if (!col_name %in% DBI::dbListFields(conn, table_name)) {
    DBI::dbExecute(conn, paste0(
      "ALTER TABLE ", table_name, " ADD COLUMN ", col_name, " ", output_type
    ))
  }
  invisible(NULL)
}


#' Write PK → value pairs to a table column via a temp-table UPDATE JOIN.
#'
#' Thin wrapper around `mutate_table_vector()` (defined in mutate_table.R).
#'
#' @keywords internal
.write_group_values <- function(conn, table_name, col_name, pk_col, pks, values) {
  mutate_table_vector(conn, table_name, col_name, pks, values, pk_col)
  invisible(NULL)
}


#' Emit group-level diagnostics after assignment.
#'
#' No-op when `quiet = TRUE`.
#'
#' @keywords internal
.summarize_group_assignment <- function(values, n_total, table_name, col_name, quiet) {
  if (quiet) return(invisible(NULL))

  n_assigned <- sum(!is.na(values))
  n_null     <- n_total - n_assigned

  msg <- paste0(
    "mutate_group: assigned ", n_assigned, " of ", n_total,
    " rows to '", col_name, "' in '", table_name, "'"
  )
  if (n_null > 0L) {
    msg <- paste0(msg, "; ", n_null, " row", if (n_null != 1L) "s" else "", " left NULL")
  }
  message(msg)

  if (n_assigned > 0L) {
    grp_counts <- table(values[!is.na(values)])
    message(
      "  Groups: n=", length(grp_counts),
      ", size min=", min(grp_counts),
      " mean=", round(mean(grp_counts), 1),
      " max=", max(grp_counts)
    )
  }

  invisible(NULL)
}


#' Largest-remainder integer rounding.
#'
#' Rounds `weights * n` to integers summing exactly to `n`. Tie-breaking is
#' stable left-to-right (earlier groups win).
#'
#' @param weights Numeric vector of proportional weights (must sum to 1).
#' @param n Integer. Total count to distribute.
#' @return Integer vector of the same length as `weights`.
#' @keywords internal
.largest_remainder <- function(weights, n) {
  exact     <- weights * n
  floors    <- floor(exact)
  remainder <- as.integer(round(n - sum(floors)))

  if (remainder > 0L) {
    fracs     <- exact - floors
    # order() is stable: equal fracs break left-to-right (lower index wins)
    order_idx <- order(-fracs, seq_along(fracs))
    bump      <- order_idx[seq_len(remainder)]
    floors[bump] <- floors[bump] + 1L
  }

  as.integer(floors)
}
