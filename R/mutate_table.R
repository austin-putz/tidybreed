#' S3 class wrapping a table reference for generic mutation
#'
#' @description
#' Returned by [get_table()]. Carries the parent population object, the
#' table name, a lazy dplyr tibble for reading, and a pending filter list.
#' Supports [dplyr::filter()], [collect()], and [mutate_table()].
#'
#' @param pop A `tidybreed_pop` object.
#' @param table_name Character scalar. Name of the table to reference.
#' @return A `tidybreed_table` S3 object.
#' @keywords internal
new_tidybreed_table <- function(pop, table_name) {
  structure(
    list(
      pop            = pop,
      table_name     = table_name,
      tbl            = dplyr::tbl(pop$db_conn, table_name),
      pending_filter = list()
    ),
    class = c("tidybreed_table", "list")
  )
}


#' Print method for tidybreed_table
#'
#' @description
#' Prints a summary header (rows × fields) plus a data preview of up to `n`
#' rows and `c` columns. Any active filter is shown below the header.
#'
#' @param x A `tidybreed_table` object.
#' @param n Maximum number of rows to preview (default 10).
#' @param c Maximum number of columns to preview (default 10).
#' @param ... Ignored.
#' @return `x` invisibly.
#' @export
print.tidybreed_table <- function(x, n = 10, c = 10, ...) {
  conn <- x$pop$db_conn

  n_rows <- DBI::dbGetQuery(
    conn, paste0("SELECT COUNT(*) AS n FROM ", x$table_name)
  )$n
  tbl_cols <- colnames(x$tbl)
  n_cols   <- length(tbl_cols)

  cat("<tidybreed_table: ", x$table_name, ">  [",
      n_rows, " row", if (n_rows != 1) "s" else "", " × ",
      n_cols, " field", if (n_cols != 1) "s" else "", "]\n", sep = "")

  if (length(x$pending_filter) > 0) {
    exprs <- vapply(x$pending_filter, rlang::as_label, character(1))
    cat("Filter: ", paste(exprs, collapse = " & "), "\n", sep = "")
  }

  preview_cols <- head(tbl_cols, c)
  preview_tbl  <- dplyr::collect(
    dplyr::select(x$tbl, dplyr::all_of(preview_cols))
  )
  print(utils::head(preview_tbl, n), n = n)

  if (n_cols > c) {
    remaining <- tbl_cols[(c + 1):n_cols]
    shown     <- min(5L, length(remaining))
    names_str <- paste(remaining[seq_len(shown)], collapse = ", ")
    if (length(remaining) > shown) names_str <- paste0(names_str, ", …")
    cat("# … with ", n_cols - c, " more field",
        if (n_cols - c != 1) "s" else "", ": ", names_str, "\n", sep = "")
  }

  invisible(x)
}


#' Filter method for tidybreed_table
#'
#' @description
#' Stashes filter predicates on the `tidybreed_table` object (for
#' [mutate_table()] to use in the SQL WHERE clause) AND applies them eagerly
#' to the underlying lazy dplyr tbl (so [collect()] continues to work
#' correctly for read-only queries).
#'
#' @param .data A `tidybreed_table` object returned by [get_table()].
#' @param ... Unquoted predicate expressions.
#' @param .preserve Not used; present for S3 signature compatibility.
#' @return The `tidybreed_table` object with the filter stashed and applied.
#' @export
filter.tidybreed_table <- function(.data, ..., .preserve = FALSE) {
  new_quosures <- rlang::enquos(...)
  .data$pending_filter <- c(.data$pending_filter, new_quosures)
  .data$tbl <- dplyr::filter(.data$tbl, !!!new_quosures)
  .data
}


#' Collect method for tidybreed_table
#'
#' @description
#' Delegates to the underlying lazy dplyr tbl, preserving backward
#' compatibility for all existing `get_table(...) |> collect()` patterns.
#'
#' @param x A `tidybreed_table` object.
#' @param ... Passed to [dplyr::collect()].
#' @return A tibble with the (filtered) table contents.
#' @importFrom dplyr collect
#' @export
collect.tidybreed_table <- function(x, ...) {
  dplyr::collect(x$tbl, ...)
}


#' Select method for tidybreed_table
#'
#' Applies column selection to the underlying lazy tbl and returns a modified
#' `tidybreed_table` so that further dplyr chains continue to work.
#'
#' @param .data A `tidybreed_table` object.
#' @param ... Column selection (tidy-select).
#' @return A `tidybreed_table` with the selection applied.
#' @export
select.tidybreed_table <- function(.data, ...) {
  .data$tbl <- dplyr::select(.data$tbl, ...)
  .data
}


#' Arrange method for tidybreed_table
#'
#' Applies row ordering to the underlying lazy tbl and returns a modified
#' `tidybreed_table` so that further dplyr chains continue to work.
#'
#' @param .data A `tidybreed_table` object.
#' @param ... Ordering expressions.
#' @return A `tidybreed_table` with the ordering applied.
#' @importFrom dplyr arrange
#' @export
arrange.tidybreed_table <- function(.data, ...) {
  .data$tbl <- dplyr::arrange(.data$tbl, ...)
  .data
}


#' Pull method for tidybreed_table
#'
#' Pulls a single column as a vector, delegating to the underlying lazy tbl.
#'
#' @param .data A `tidybreed_table` object.
#' @param var Column to pull (tidy-select).
#' @param name Optional column to use as names.
#' @param ... Passed to [dplyr::pull()].
#' @return A vector.
#' @importFrom dplyr pull
#' @export
pull.tidybreed_table <- function(.data, var = -1, name = NULL, ...) {
  dplyr::pull(.data$tbl, {{ var }}, ...)
}


#' Count method for tidybreed_table
#'
#' Counts rows (optionally by group), delegating to the underlying lazy tbl.
#'
#' @param x A `tidybreed_table` object.
#' @param ... Grouping columns (tidy-select).
#' @param wt Optional weighting column.
#' @param sort Logical; sort by descending count?
#' @param name Name for the count column.
#' @return A lazy tibble.
#' @importFrom dplyr count
#' @export
count.tidybreed_table <- function(x, ..., wt = NULL, sort = FALSE, name = "n") {
  dplyr::count(x$tbl, ..., sort = sort, name = name)
}


#' slice_max method for tidybreed_table
#'
#' Applies [dplyr::slice_max()] to the underlying lazy tbl and returns a
#' modified `tidybreed_table` so that further dplyr chains continue to work.
#'
#' @param .data A `tidybreed_table` object.
#' @param order_by Column or expression to order by (tidy-eval).
#' @param ... Passed to [dplyr::slice_max()].
#' @return A `tidybreed_table` with the slice applied.
#' @importFrom dplyr slice_max
#' @export
slice_max.tidybreed_table <- function(.data, order_by, ...) {
  .data$tbl <- dplyr::slice_max(.data$tbl, order_by = {{ order_by }}, ...)
  .data
}


#' slice_min method for tidybreed_table
#'
#' Applies [dplyr::slice_min()] to the underlying lazy tbl and returns a
#' modified `tidybreed_table` so that further dplyr chains continue to work.
#'
#' @param .data A `tidybreed_table` object.
#' @param order_by Column or expression to order by (tidy-eval).
#' @param ... Passed to [dplyr::slice_min()].
#' @return A `tidybreed_table` with the slice applied.
#' @importFrom dplyr slice_min
#' @export
slice_min.tidybreed_table <- function(.data, order_by, ...) {
  .data$tbl <- dplyr::slice_min(.data$tbl, order_by = {{ order_by }}, ...)
  .data
}


#' slice_head method for tidybreed_table
#'
#' Applies [dplyr::slice_head()] to the underlying lazy tbl and returns a
#' modified `tidybreed_table` so that further dplyr chains continue to work.
#'
#' @param .data A `tidybreed_table` object.
#' @param ... Passed to [dplyr::slice_head()] (e.g. `n = 10`).
#' @return A `tidybreed_table` with the slice applied.
#' @importFrom dplyr slice_head
#' @export
slice_head.tidybreed_table <- function(.data, ...) {
  .data$tbl <- dplyr::slice_head(.data$tbl, ...)
  .data
}


#' slice_tail method for tidybreed_table
#'
#' Applies [dplyr::slice_tail()] to the underlying lazy tbl and returns a
#' modified `tidybreed_table` so that further dplyr chains continue to work.
#'
#' @param .data A `tidybreed_table` object.
#' @param ... Passed to [dplyr::slice_tail()] (e.g. `n = 10`).
#' @return A `tidybreed_table` with the slice applied.
#' @importFrom dplyr slice_tail
#' @export
slice_tail.tidybreed_table <- function(.data, ...) {
  .data$tbl <- dplyr::slice_tail(.data$tbl, ...)
  .data
}


#' slice_sample method for tidybreed_table
#'
#' Applies [dplyr::slice_sample()] to the underlying lazy tbl and returns a
#' modified `tidybreed_table` so that further dplyr chains continue to work.
#'
#' @param .data A `tidybreed_table` object.
#' @param ... Passed to [dplyr::slice_sample()] (e.g. `n = 10`, `prop = 0.1`).
#' @return A `tidybreed_table` with the slice applied.
#' @importFrom dplyr slice_sample
#' @export
slice_sample.tidybreed_table <- function(.data, ...) {
  .data$tbl <- dplyr::slice_sample(.data$tbl, ...)
  .data
}


# ── Internal helpers ──────────────────────────────────────────────────────────

#' Count non-NULL values in a filtered set of rows
#'
#' Used by [mutate_table()] to report how many existing values will be
#' replaced vs. filled in before executing an UPDATE.
#'
#' @param conn DuckDB connection
#' @param table_name Target table name
#' @param field_name Column to inspect
#' @param pk_col Primary key column name
#' @param filter_ids Vector of PK values identifying the affected rows
#' @keywords internal
count_nonnull_in_filter <- function(conn, table_name, field_name, pk_col, filter_ids) {
  filter_df <- data.frame(pk = filter_ids, stringsAsFactors = FALSE)
  names(filter_df) <- pk_col
  DBI::dbWriteTable(conn, "_tb_cnt_nn", filter_df, overwrite = TRUE)
  n <- DBI::dbGetQuery(conn, paste0(
    "SELECT COUNT(*) AS n FROM ", table_name, " AS t ",
    "INNER JOIN _tb_cnt_nn AS f ON t.", pk_col, " = f.", pk_col,
    " WHERE t.", field_name, " IS NOT NULL"
  ))$n
  DBI::dbExecute(conn, "DROP TABLE _tb_cnt_nn")
  n
}


#' Execute a scalar UPDATE on a table column
#'
#' @param conn DuckDB connection
#' @param table_name Target table name
#' @param field_name Column to update
#' @param value Scalar R value
#' @param db_type DuckDB type string
#' @param filter_ids NULL = update all rows; vector of PK values = update only those rows
#' @param pk_col Primary key column name (required when filter_ids is not NULL)
#' @keywords internal
mutate_table_scalar <- function(conn, table_name, field_name, value, db_type,
                                filter_ids, pk_col) {
  sql_value <- format_sql_value(value, db_type)

  if (is.null(filter_ids)) {
    DBI::dbExecute(conn, paste0(
      "UPDATE ", table_name, " SET ", field_name, " = ", sql_value
    ))
    return(invisible(NULL))
  }

  filter_df  <- data.frame(pk = filter_ids, stringsAsFactors = FALSE)
  names(filter_df) <- pk_col
  temp_name  <- paste0("_tb_sc_", gsub("[^a-zA-Z0-9]", "_", field_name))
  DBI::dbWriteTable(conn, temp_name, filter_df, overwrite = TRUE)

  DBI::dbExecute(conn, paste0(
    "UPDATE ", table_name,
    " SET ", field_name, " = ", sql_value,
    " FROM ", temp_name, " AS f",
    " WHERE ", table_name, ".", pk_col, " = f.", pk_col
  ))
  DBI::dbExecute(conn, paste0("DROP TABLE ", temp_name))
  invisible(NULL)
}


#' Execute a vector UPDATE on a table column via a temp-table JOIN
#'
#' @param conn DuckDB connection
#' @param table_name Target table name
#' @param field_name Column to update
#' @param pks_ordered PK values in canonical row order (from `get_pks_in_order()`)
#' @param values Vector of new values, same length and order as `pks_ordered`
#' @param pk_col Primary key column name
#' @keywords internal
mutate_table_vector <- function(conn, table_name, field_name,
                                pks_ordered, values, pk_col) {
  update_df  <- data.frame(pk = pks_ordered, new_value = values,
                            stringsAsFactors = FALSE)
  names(update_df)[1] <- pk_col

  temp_name <- paste0("_tb_vec_", gsub("[^a-zA-Z0-9]", "_", field_name))
  DBI::dbWriteTable(conn, temp_name, update_df, overwrite = TRUE)

  DBI::dbExecute(conn, paste0(
    "UPDATE ", table_name,
    " SET ", field_name, " = t.new_value",
    " FROM ", temp_name, " AS t",
    " WHERE ", table_name, ".", pk_col, " = t.", pk_col
  ))
  DBI::dbExecute(conn, paste0("DROP TABLE ", temp_name))
  invisible(NULL)
}


# ── Exported function ─────────────────────────────────────────────────────────

#' Add or modify columns in any population database table
#'
#' @description
#' Generic replacement for the table-specific `mutate_ind_meta()` /
#' `mutate_genome_meta()` family. Chain after [get_table()] (and optionally
#' [dplyr::filter()]) to add new columns or update existing ones. `value` can
#' be a scalar broadcast to all (or filtered) rows, or a tibble keyed on the
#' table's primary key for per-row values — see `@examples` below for both,
#' plus the schema-pre-declaration pattern for empty tables.
#'
#' @param tbl_obj A `tidybreed_table` object returned by [get_table()].
#' @param ... Named arguments of the form `column_name = value`. `value` can
#'   be a scalar (applied to all affected rows) or a tibble/data frame
#'   containing the primary key column and a column named after the field.
#'   Plain vectors (length > 1) are not allowed; use a tibble instead.
#' @param .set_default Logical; if `TRUE`, creates new columns with a SQL
#'   DEFAULT constraint set to the provided value. The DEFAULT applies to future
#'   INSERT operations (e.g., from [add_founders()], [add_phenotype()]) when
#'   the column is not explicitly specified in `...`. Only valid for scalar
#'   values; an error is raised if used with tibbles. Has no effect
#'   on columns that already exist (DuckDB does not support modifying existing
#'   column defaults). Default: `FALSE`.
#'
#' @return The parent `tidybreed_pop` object (invisibly).
#'
#' @details
#' **Type inference** (`logical` → BOOLEAN, `integer` → INTEGER,
#' `numeric` → DOUBLE, `Date` → DATE, `POSIXct` → TIMESTAMP,
#' `character` → VARCHAR) is performed automatically via
#' `infer_duckdb_type()`.
#'
#' **Per-row values (tibble path)**: to set different values for different rows,
#' pass a tibble/data frame containing the primary key column plus a column named
#' after the field being updated:
#' ```
#' mutate_table(gen = tibble::tibble(id_ind = c("A_1", "A_2"),
#'                                    gen = c(1L, 2L)))
#' ```
#' The tibble is entirely self-contained; any upstream [dplyr::filter()] is ignored
#' for that field (but continues to apply to scalar fields in the same call).
#' All IDs in the tibble must exist in the table; an error is raised for unknown IDs.
#'
#' **Filtering**: when a [dplyr::filter()] is applied upstream, scalar values broadcast
#' to all matching rows. New columns created in this context will be `NULL` for
#' all non-matching rows.
#'
#' **Default constraints**: When `.set_default = TRUE`, new columns are created
#' with a SQL DEFAULT constraint. This means future INSERT operations from
#' [add_founders()], [add_phenotype()], or other `add_*()` functions will
#' automatically use the default value when the column is not explicitly
#' provided in `...`.
#'
#' For populated tables, existing rows are still updated via the standard UPDATE
#' mechanism; the DEFAULT only affects subsequent INSERT operations. For empty
#' tables (e.g., right after [open_pop()] / [define_genome()], before any
#' founders are added), the DEFAULT is set at the schema level with no data
#' operations.
#'
#' Note: DEFAULT constraints can only be added when creating new columns. If a
#' column already exists, `.set_default` is ignored (but the UPDATE proceeds
#' normally).
#'
#' @examples
#' \dontrun{
#' pop <- open_pop(pop_name = "sim", db_name = ":memory:") |>
#'   define_genome(n_loci = 100, n_chr = 5, chr_len_Mb = 50) |>
#'   define_founder_haplotypes(n_haplotypes = 100, line_name = "A")
#' pop <- pop |>
#'   get_table("founder_haplotypes") |>
#'   add_founders(n_males = 10, n_females = 100, line_name = "A")
#'
#' # Broadcast: assign gen = 1 to every individual
#' pop <- pop |>
#'   get_table("ind_meta") |>
#'   mutate_table(gen = 1L)
#'
#' # Filtered subset: bump gen only for males
#' pop <- pop |>
#'   get_table("ind_meta") |>
#'   dplyr::filter(sex == "M") |>
#'   mutate_table(gen = 2L)
#'
#' # Per-row values via a tibble keyed on the primary key column
#' pop <- pop |>
#'   get_table("ind_meta") |>
#'   mutate_table(
#'     gen = tibble::tibble(
#'       id_ind = c("A_1", "A_2"),
#'       gen    = c(5L, 6L)
#'     )
#'   )
#'
#' # Pre-declare a typed column schema before any data exists
#' pop2 <- open_pop(pop_name = "sim2", db_name = ":memory:") |>
#'   define_genome(n_loci = 100, n_chr = 5, chr_len_Mb = 50)
#' pop2 <- pop2 |>
#'   get_table("ind_meta") |>
#'   mutate_table(farm = NA_character_)
#'
#' # Pre-declare WITH a SQL DEFAULT so future add_founders() rows auto-fill
#' pop3 <- open_pop(pop_name = "sim3", db_name = ":memory:") |>
#'   define_genome(n_loci = 100, n_chr = 5, chr_len_Mb = 50) |>
#'   define_founder_haplotypes(n_haplotypes = 100, line_name = "A")
#' pop3 <- pop3 |>
#'   get_table("ind_meta") |>
#'   mutate_table(gen = 0L, active = TRUE, .set_default = TRUE)
#'
#' # Future founders automatically get gen = 0L and active = TRUE
#' pop3 <- pop3 |>
#'   get_table("founder_haplotypes") |>
#'   add_founders(n_males = 10, n_females = 100, line_name = "A")
#' pop3 |> get_table("ind_meta") |> dplyr::collect()  # all have defaults
#'
#' # Override the default by supplying an explicit value
#' pop3 <- pop3 |>
#'   get_table("founder_haplotypes") |>
#'   dplyr::filter(line_name == "A") |>
#'   add_founders(n_males = 5, n_females = 50, line_name = "B", gen = 1L)
#'
#' # Line B: gen = 1L (explicit), active = TRUE (default)
#' pop3 |>
#'   get_table("ind_meta") |>
#'   dplyr::filter(line_name == "B") |>
#'   dplyr::collect()
#' }
#'
#' @export
mutate_table <- function(tbl_obj, ..., .set_default = FALSE) {

  if (!inherits(tbl_obj, "tidybreed_table")) {
    stop(
      "mutate_table() must be called after get_table(). ",
      "Use: pop %>% get_table('table_name') %>% mutate_table(...)",
      call. = FALSE
    )
  }

  pop        <- tbl_obj$pop
  table_name <- tbl_obj$table_name

  validate_tidybreed_pop(pop)

  # Extract field values, excluding .set_default if it was passed in ...
  all_args <- list(...)
  if (".set_default" %in% names(all_args)) {
    .set_default <- all_args[[".set_default"]]
    field_values <- all_args[names(all_args) != ".set_default"]
  } else {
    field_values <- all_args
  }

  if (length(field_values) == 0) {
    warning("No fields specified. Returning population unchanged.", call. = FALSE)
    return(invisible(pop))
  }

  # Validate .set_default parameter
  if (!is.null(.set_default)) {
    stopifnot(
      is.logical(.set_default),
      length(.set_default) == 1,
      !is.na(.set_default)
    )
  }

  existing_cols <- DBI::dbListFields(pop$db_conn, table_name)

  reserved <- TABLE_RESERVED_COLS[[table_name]]
  if (is.null(reserved)) reserved <- character(0)

  pk_col <- TABLE_PRIMARY_KEYS[[table_name]]

  n_total <- DBI::dbGetQuery(
    pop$db_conn,
    paste0("SELECT COUNT(*) AS n FROM ", table_name)
  )$n

  if (n_total == 0) {
    # Table is empty — still create any new columns (schema-only ALTER TABLE).
    # This supports pre-declaring typed column schemas before data arrives, e.g.:
    #   pop |> get_table("ind_meta") |> mutate_table(gen = NA_integer_)
    for (field_name in names(field_values)) {
      value <- field_values[[field_name]]
      validate_sql_identifier(field_name, what = "field name", reserved = reserved)
      if (!field_name %in% existing_cols) {
        db_type <- infer_duckdb_type(value)

        # Build ALTER TABLE with optional DEFAULT
        alter_sql <- paste0(
          "ALTER TABLE ", table_name, " ADD COLUMN ", field_name, " ", db_type
        )
        if (!is.null(.set_default) && .set_default) {
          sql_value <- format_sql_value(value, db_type)
          alter_sql <- paste0(alter_sql, " DEFAULT ", sql_value)
        }

        DBI::dbExecute(pop$db_conn, alter_sql)

        # Message with DEFAULT indicator
        default_msg <- if (!is.null(.set_default) && .set_default) " with DEFAULT" else ""
        message("Added new column '", field_name, "' (", db_type, ") to `",
                table_name, "` (table is currently empty)", default_msg)
      } else {
        message("Column '", field_name, "' already exists in `", table_name,
                "` (table is currently empty; no rows updated)")
      }
    }
    return(invisible(pop))
  }

  # Resolve filter
  has_filter  <- length(tbl_obj$pending_filter) > 0
  filter_ids  <- NULL
  n_effective <- n_total

  if (has_filter) {
    if (is.null(pk_col)) {
      stop(
        "Cannot use filter() with mutate_table() on table '", table_name,
        "': primary key not registered for this table.",
        call. = FALSE
      )
    }
    filter_ids  <- dplyr::collect(dplyr::select(tbl_obj$tbl, dplyr::all_of(pk_col)))[[pk_col]]
    n_effective <- length(filter_ids)
    if (n_effective == 0) {
      warning(
        "filter() matched 0 rows in '", table_name, "'. No values updated.",
        call. = FALSE
      )
      return(invisible(pop))
    }
  }

  for (field_name in names(field_values)) {

    value <- field_values[[field_name]]

    validate_sql_identifier(field_name, what = "field name", reserved = reserved)

    # Dispatch on input type: tibble (data frame), scalar, or invalid vector
    if (is.data.frame(value)) {
      # --- Tibble / keyed path ---
      if (is.null(pk_col)) {
        stop(
          "Cannot use a tibble value in mutate_table() on table '", table_name,
          "': no primary key registered for this table.",
          call. = FALSE
        )
      }
      if (!pk_col %in% names(value)) {
        stop(
          "Tibble for field '", field_name, "' must contain a column named '",
          pk_col, "' (the primary key for '", table_name, "').",
          call. = FALSE
        )
      }
      if (!field_name %in% names(value)) {
        stop(
          "Tibble for field '", field_name, "' must contain a column named '",
          field_name, "'.",
          call. = FALSE
        )
      }

      # Validate all IDs exist in the table
      tbl_ids <- DBI::dbGetQuery(pop$db_conn,
                                  paste0("SELECT ", pk_col, " FROM ", table_name))[[pk_col]]
      unknown <- setdiff(value[[pk_col]], tbl_ids)
      if (length(unknown) > 0) {
        stop(
          "Tibble for field '", field_name, "' contains ", length(unknown),
          " ID(s) not found in '", table_name, "': ",
          paste(head(unknown, 5), collapse = ", "),
          if (length(unknown) > 5) ", …" else "",
          call. = FALSE
        )
      }

      is_vector      <- TRUE
      n_effective    <- nrow(value)
      db_type        <- infer_duckdb_type(value[[field_name]])
      pks_tibble     <- value[[pk_col]]
      vals_tibble    <- value[[field_name]]
      tibble_filter  <- pks_tibble

    } else if (length(value) == 1) {
      # --- Scalar path (unchanged) ---
      db_type <- infer_duckdb_type(value)
      is_vector <- FALSE
      tibble_filter <- NULL

    } else {
      # --- Plain vector (length > 1) — error ---
      stop(
        "Plain vectors are not allowed in mutate_table() — row order is not ",
        "guaranteed and can silently assign wrong values to wrong rows.\n\n",
        "Use a tibble with the primary key column instead:\n",
        "  mutate_table(", field_name, " = tibble::tibble(\n",
        "    ", pk_col %||% "pk", " = c(...),\n",
        "    ", field_name, " = c(...)\n",
        "  ))",
        call. = FALSE
      )
    }

    # Block .set_default with tibbles
    if (is_vector && !is.null(.set_default) && .set_default) {
      stop(
        "Cannot use .set_default = TRUE with tibble/vector values for field '",
        field_name, "'. DEFAULT constraints require scalar values.",
        call. = FALSE
      )
    }

    is_new_col <- !field_name %in% existing_cols
    n_replaced <- 0L
    n_filled   <- 0L

    if (is_new_col) {
      # Build ALTER TABLE with optional DEFAULT
      alter_sql <- paste0(
        "ALTER TABLE ", table_name, " ADD COLUMN ", field_name, " ", db_type
      )
      if (!is.null(.set_default) && .set_default) {
        sql_value <- format_sql_value(value, db_type)
        alter_sql <- paste0(alter_sql, " DEFAULT ", sql_value)
      }

      DBI::dbExecute(pop$db_conn, alter_sql)
      existing_cols <- c(existing_cols, field_name)
    } else {
      # Count non-NULL rows among the affected rows BEFORE updating so we can
      # report exactly what changed (replacements vs. NULL fills).
      # For tibble path, use the tibble IDs; for filter/scalar, use filter_ids
      affected_ids <- if (is_vector) tibble_filter else filter_ids
      if (!is.null(affected_ids)) {
        n_replaced <- count_nonnull_in_filter(
          pop$db_conn, table_name, field_name, pk_col, affected_ids
        )
      } else {
        n_replaced <- DBI::dbGetQuery(
          pop$db_conn,
          paste0("SELECT COUNT(*) AS n FROM ", table_name,
                 " WHERE ", field_name, " IS NOT NULL")
        )$n
      }
      n_filled <- n_effective - n_replaced
    }

    if (is_vector) {
      mutate_table_vector(pop$db_conn, table_name, field_name, pks_tibble, vals_tibble, pk_col)
    } else {
      mutate_table_scalar(pop$db_conn, table_name, field_name, value, db_type,
                          filter_ids, pk_col)
    }

    # Emit informational message
    if (is_new_col) {
      default_msg <- if (!is.null(.set_default) && .set_default) " with DEFAULT" else ""
      if (has_filter) {
        n_null <- n_total - n_effective
        message(
          "Added new column '", field_name, "' (", db_type, ") to `", table_name, "`",
          default_msg, "; ",
          n_effective, " row", if (n_effective != 1) "s" else "", " set, ",
          n_null, " row", if (n_null != 1) "s" else "", " NULL"
        )
      } else {
        message(
          "Added new column '", field_name, "' (", db_type, ") to `", table_name, "`",
          default_msg, "; ",
          n_total, " row", if (n_total != 1) "s" else "", " set"
        )
      }
    } else {
      scope_str <- if (has_filter) {
        paste0(" in `", table_name, "` [", n_effective, " of ", n_total, " rows]")
      } else {
        paste0(" in `", table_name, "` [", n_total,
               " row", if (n_total != 1) "s" else "", "]")
      }
      if (n_replaced > 0) {
        warning(
          "'", field_name, "': replaced ", n_replaced, " existing value",
          if (n_replaced != 1) "s" else "",
          scope_str,
          call. = FALSE
        )
      }
      if (n_filled > 0) {
        message(
          "'", field_name, "': filled ", n_filled, " NULL row",
          if (n_filled != 1) "s" else "",
          scope_str
        )
      }
    }
  }

  invisible(pop)
}
