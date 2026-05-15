#' Compute and write a derived column from a cross-table join
#'
#' @description
#' A pipe-friendly action function that:
#' 1. Collects the current (optionally filtered) primary table
#' 2. Left-joins an optional secondary table by `join_by`
#' 3. Applies a user-supplied `compute` function to produce a derived vector
#' 4. Writes the result to one or more destination columns across any tables
#'
#' This reduces ~25-line dplyr join + `mutate_table()` workflows for computed
#' fields (e.g. `puberty_date = birth_date + AP_value`) to a single pipe step.
#'
#' @param tbl_obj A `tidybreed_table` returned by [get_table()] with an
#'   optional [filter()] applied.
#' @param compute A function `function(df) -> vector`. `df` is the joined
#'   tibble (primary table left-joined to `join_table`). The returned vector
#'   must be a scalar (broadcast to all rows) or have the same length as
#'   `nrow(df)`. Primary table column names take precedence on name collisions;
#'   any colliding column from `join_table` is suffixed with `".<join_table>"`.
#' @param join_table Character scalar or `NULL`. Name of a secondary table to
#'   left-join with the primary table. Must exist in `pop$tables`. `NULL`
#'   (default) skips the join; `compute` receives only the primary table rows.
#' @param join_by Character scalar. Column name to join on; must exist in both
#'   tables when `join_table` is not `NULL`. Default `"id_ind"`.
#' @param write_to Named character vector. Names are destination table names;
#'   values are destination column names. Example:
#'   `c(ind_meta = "puberty_date", ind_phenotype = "pheno_date")`.
#'   Each destination table must be registered in the package primary-key map
#'   (`TABLE_PRIMARY_KEYS`). Reserved columns are blocked.
#'
#' @return The parent `tidybreed_pop` object (invisibly).
#'
#' @details
#' ## Writing back to the primary table
#'
#' Uses the table's own primary key for exact-row matching. A filtered
#' `ind_phenotype` (e.g. `filter(trait_name == "AP")`) only updates those
#' specific rows — other trait records for the same individual are untouched.
#'
#' ## Writing to a secondary table
#'
#' Requires at most one computed value per `join_by` ID. If the primary table
#' has multiple rows for the same individual that produce *different* computed
#' values, `compute_derived()` errors with a message asking you to pre-filter.
#' Identical values for the same `join_by` ID are silently deduplicated (first
#' occurrence wins).
#'
#' ## Column name collisions during the join
#'
#' If the primary and secondary tables share a column name other than `join_by`,
#' the secondary table's copy gains a `".<join_table>"` suffix inside `compute`.
#' Primary table columns are always unambiguous (no suffix).
#'
#' ## join_table row expansion
#'
#' If `join_table` itself contains multiple rows per `join_by` value, the
#' left-join silently expands the primary table rows. Pre-filter `join_table`
#' externally before passing it if this is a concern.
#'
#' ## New column type inference
#'
#' When a destination column does not yet exist, its DuckDB type is inferred
#' from `result_vec` via `infer_duckdb_type()`. R's `Date` arithmetic preserves
#' the `Date` class, so `birth_date + integer` correctly produces a `DATE`
#' column. Use typed `NA` vectors (e.g. `NA_real_`) if `result_vec` is all-NA.
#'
#' @examples
#' \dontrun{
#' # Age-at-puberty: puberty_date = birth_date (ind_meta) + AP value (days)
#' pop <- pop |>
#'   get_table("ind_phenotype") |>
#'   dplyr::filter(trait_name == "AP", pheno_number == 1L) |>
#'   compute_derived(
#'     compute    = \(df) df$birth_date + as.integer(df$value),
#'     join_table = "ind_meta",
#'     join_by    = "id_ind",
#'     write_to   = c(ind_meta = "puberty_date", ind_phenotype = "pheno_date")
#'   )
#'
#' # Compute entirely from the primary table (no join needed)
#' pop <- pop |>
#'   get_table("ind_phenotype") |>
#'   dplyr::filter(trait_name == "BW") |>
#'   compute_derived(
#'     compute  = \(df) df$value * 2.205,   # kg -> lb
#'     write_to = c(ind_phenotype = "value_lb")
#'   )
#' }
#'
#' @export
compute_derived <- function(tbl_obj, compute, join_table = NULL,
                             join_by = "id_ind", write_to) {

  # ── 1. Input validation ───────────────────────────────────────────────────

  if (!inherits(tbl_obj, "tidybreed_table")) {
    stop(
      "compute_derived() must be called after get_table(). ",
      "Use: pop |> get_table('table_name') |> compute_derived(...)",
      call. = FALSE
    )
  }

  pop        <- tbl_obj$pop
  table_name <- tbl_obj$table_name
  con        <- pop$db_conn
  validate_tidybreed_pop(pop)

  if (!is.function(compute)) {
    stop(
      "'compute' must be a function, e.g. \\(df) df$birth_date + as.integer(df$value)",
      call. = FALSE
    )
  }

  if (is.null(write_to) || !is.character(write_to) ||
      is.null(names(write_to)) || any(!nzchar(names(write_to)))) {
    stop(
      "'write_to' must be a named character vector. ",
      "Example: c(ind_meta = 'puberty_date', ind_phenotype = 'pheno_date')",
      call. = FALSE
    )
  }

  if (!is.null(join_table)) {
    if (!is.character(join_table) || length(join_table) != 1L || !nzchar(join_table)) {
      stop("'join_table' must be a single non-empty character string or NULL.", call. = FALSE)
    }
    if (!join_table %in% pop$tables) {
      stop("join_table '", join_table, "' does not exist in this population.", call. = FALSE)
    }
  }

  if (!is.character(join_by) || length(join_by) != 1L || !nzchar(join_by)) {
    stop("'join_by' must be a single non-empty character string.", call. = FALSE)
  }

  # Validate each destination: SQL safety, reserved cols, and PK registration
  for (i in seq_along(write_to)) {
    dest_table <- names(write_to)[i]
    dest_col   <- write_to[i]
    reserved   <- TABLE_RESERVED_COLS[[dest_table]]
    if (is.null(reserved)) reserved <- character(0)
    validate_sql_identifier(dest_col, what = "destination column", reserved = reserved)

    if (is.null(TABLE_PRIMARY_KEYS[[dest_table]])) {
      stop(
        "Destination table '", dest_table, "' is not registered in the primary-key ",
        "map (TABLE_PRIMARY_KEYS). Known tables: ",
        paste(names(TABLE_PRIMARY_KEYS), collapse = ", "),
        call. = FALSE
      )
    }
    if (!dest_table %in% pop$tables) {
      stop(
        "Destination table '", dest_table, "' does not exist in this population.",
        call. = FALSE
      )
    }
  }

  # ── 2. Collect filtered primary table ─────────────────────────────────────

  df_primary <- dplyr::collect(tbl_obj$tbl)

  if (nrow(df_primary) == 0L) {
    message("compute_derived(): primary table filter matched 0 rows; nothing written.")
    return(invisible(pop))
  }

  # join_by presence check (only needed when a join_table or secondary writes exist)
  needs_join_by <- !is.null(join_table) ||
    any(names(write_to) != table_name)

  if (needs_join_by && !join_by %in% names(df_primary)) {
    stop(
      "join_by column '", join_by, "' not found in primary table '", table_name, "'.",
      call. = FALSE
    )
  }

  # ── 3. Optional join with secondary table ─────────────────────────────────

  if (!is.null(join_table)) {
    join_table_cols <- DBI::dbListFields(con, join_table)
    if (!join_by %in% join_table_cols) {
      stop(
        "join_by column '", join_by, "' not found in join_table '", join_table, "'.",
        call. = FALSE
      )
    }
    join_ids <- unique(df_primary[[join_by]])
    df_join  <- dplyr::tbl(con, join_table) |>
      dplyr::filter(.data[[join_by]] %in% !!join_ids) |>
      dplyr::collect()

    df_joined <- dplyr::left_join(
      df_primary, df_join,
      by     = join_by,
      suffix = c("", paste0(".", join_table))
    )
  } else {
    df_joined <- df_primary
  }

  # ── 4. Apply compute function ─────────────────────────────────────────────

  result_vec <- tryCatch(
    compute(df_joined),
    error = function(e) {
      stop("Error in compute() function: ", conditionMessage(e), call. = FALSE)
    }
  )

  if (length(result_vec) == 1L) {
    result_vec <- rep(result_vec, nrow(df_joined))
  } else if (length(result_vec) != nrow(df_joined)) {
    stop(
      "compute() returned ", length(result_vec), " values but the joined data has ",
      nrow(df_joined), " rows. The function must return a scalar or a vector of the ",
      "same length as the (joined) primary table.",
      call. = FALSE
    )
  }

  # ── 5. Write to each destination ──────────────────────────────────────────

  for (i in seq_along(write_to)) {
    dest_table <- names(write_to)[i]
    dest_col   <- write_to[i]
    dest_pk    <- TABLE_PRIMARY_KEYS[[dest_table]]

    is_primary_table <- identical(dest_table, table_name)

    if (is_primary_table) {
      # Exact-row match using the primary table's own PK — correctly handles
      # multiple rows per id_ind (e.g. different traits in ind_phenotype).
      pk_col_primary <- TABLE_PRIMARY_KEYS[[table_name]]
      if (!pk_col_primary %in% names(df_primary)) {
        stop(
          "Primary key column '", pk_col_primary, "' not found in the collected ",
          "primary table '", table_name, "'. Ensure the PK column is not excluded ",
          "by a prior select().",
          call. = FALSE
        )
      }
      pk_vals         <- df_joined[[pk_col_primary]]
      values_to_write <- result_vec

    } else {
      # Secondary table: deduplicate to one computed value per join_by ID.
      if (!join_by %in% names(df_joined)) {
        stop(
          "Cannot write to '", dest_table, "': join_by column '", join_by,
          "' not found in joined data.",
          call. = FALSE
        )
      }

      join_ids_vec <- df_joined[[join_by]]
      uid          <- unique(join_ids_vec)

      # Build a data.frame to preserve type (e.g. Date) through dedup
      result_df         <- data.frame(join_id = join_ids_vec, stringsAsFactors = FALSE)
      result_df$.result <- result_vec   # preserves Date/POSIXct class

      # Conflict check: same join_by ID → different non-NA computed value
      if (nrow(result_df) > length(uid)) {
        conflict_ids <- character(0)
        for (jid in uid) {
          vals   <- result_df$.result[result_df$join_id == jid]
          non_na <- vals[!is.na(vals)]
          if (length(non_na) > 0L && length(unique(non_na)) > 1L) {
            conflict_ids <- c(conflict_ids, as.character(jid))
          }
        }
        if (length(conflict_ids) > 0L) {
          stop(
            "The same '", join_by, "' value maps to different computed results for: ",
            paste(head(conflict_ids, 5L), collapse = ", "),
            if (length(conflict_ids) > 5L) {
              paste0(" (and ", length(conflict_ids) - 5L, " more)")
            } else {
              ""
            },
            ". Filter the primary table to one row per ", join_by, " before ",
            "calling compute_derived(), or only write back to the primary table.",
            call. = FALSE
          )
        }
      }

      # Keep first occurrence per join_by ID
      dedup_df <- result_df[!duplicated(result_df$join_id), , drop = FALSE]

      # Fetch dest PK + join_by from the destination table
      dest_pk_df <- dplyr::tbl(con, dest_table) |>
        dplyr::select(dplyr::all_of(c(dest_pk, join_by))) |>
        dplyr::filter(.data[[join_by]] %in% !!uid) |>
        dplyr::collect()

      # Align: merge preserves Date/POSIXct class on .result column
      merged <- merge(
        dest_pk_df, dedup_df,
        by.x    = join_by,
        by.y    = "join_id",
        all.x   = TRUE,
        sort    = FALSE
      )
      values_to_write <- merged$.result
      pk_vals         <- merged[[dest_pk]]
    }

    # ── 5a. Add column if not present ─────────────────────────────────────────
    existing_cols <- DBI::dbListFields(con, dest_table)
    if (!dest_col %in% existing_cols) {
      db_type <- infer_duckdb_type(result_vec)
      DBI::dbExecute(con, paste0(
        "ALTER TABLE ", dest_table, " ADD COLUMN ", dest_col, " ", db_type
      ))
      message(
        "compute_derived(): added new column '", dest_col,
        "' (", db_type, ") to `", dest_table, "`"
      )
    }

    # ── 5b. Write via temp-table JOIN + UPDATE ────────────────────────────────
    mutate_table_vector(
      conn        = con,
      table_name  = dest_table,
      field_name  = dest_col,
      pks_ordered = pk_vals,
      values      = values_to_write,
      pk_col      = dest_pk
    )

    n_written <- sum(!is.na(values_to_write))
    message(
      "compute_derived(): wrote ", length(pk_vals), " rows to `",
      dest_table, "$", dest_col, "` (", n_written, " non-NA)"
    )
  }

  invisible(pop)
}
