#' Compute a selection index from a tidybreed table
#'
#' @description
#' Calculates a named selection index by multiplying each individual's values by
#' the index weights defined in [define_index()] and summing the products.
#' Results are appended to the `ind_index` table with an auto-incrementing
#' `index_number` so that successive runs are preserved.
#'
#' The first argument must be a `tidybreed_table` obtained from `get_table()`
#' (optionally filtered). The table must contain `id_ind`, a trait/phenotype
#' name column, and a numeric value column (`value_col`). For `ind_phenotype`
#' the name column is `phenotype_name`; for every other table (`ind_ebv`,
#' `ind_tbv`, or a user-defined table) it is `trait_name` — this is detected
#' automatically from the table's columns. Works with `ind_ebv` (EBVs),
#' `ind_phenotype` (phenotypes), `ind_tbv` (TBVs), or any user-defined table
#' with the same structure.
#'
#' ## Uniqueness requirement
#' There must be exactly one row per `(id_ind, trait_name)` after any filter is
#' applied. If duplicates remain, an error is thrown — filter the table down to
#' a single model, evaluation, or phenotype record before calling `add_index()`.
#'
#' ## Completeness requirement
#' Every individual must have a value for **every** trait in the index. If any
#' individual is missing a value for any required trait, an error is thrown
#' before any rows are written.
#'
#' @param tbl A `tidybreed_table` from `get_table()` (optionally filtered).
#'   Must contain `id_ind`, a name column (`trait_name`, or `phenotype_name`
#'   for `ind_phenotype`), and the column specified by `value_col`. Any table
#'   with this structure is accepted: `ind_ebv`, `ind_phenotype`, `ind_tbv`,
#'   or a custom table.
#' @param index_name Character scalar. Name of the index to compute; must
#'   already exist in `index_meta` (created via [define_index()]).
#' @param value_col Character scalar or `NULL`. The column in `tbl` that holds
#'   the numeric value to weight. When `NULL` (default), auto-detected from the
#'   table name: `ind_ebv` → `"ebv_value"`, `ind_phenotype` → `"pheno_value"`,
#'   `ind_tbv` → `"tbv_value"`. An error is thrown for unknown tables if
#'   `value_col` is not supplied.
#' @param overwrite_index Logical. If `TRUE`, all existing `ind_index` rows for
#'   this `index_name` are deleted before inserting; new rows receive
#'   `index_number = 1`. Default `FALSE`.
#' @param delete_all Logical. If `TRUE`, **all** rows in `ind_index` are deleted
#'   before inserting; new rows receive `index_number = 1`. Takes precedence
#'   over `overwrite_index`. Default `FALSE`.
#' @param ... Optional extra columns to add to `ind_index`. Scalar values only
#'   (broadcast to all inserted rows). Types are inferred via
#'   [infer_duckdb_type()].
#'
#' @return The `tidybreed_pop` (invisibly). Assign the result back.
#'
#' @seealso [define_index()], [get_table()], [add_ebv()], [add_phenotype()],
#'   [add_tbv()]
#'
#' @examples
#' \dontrun{
#' # Compute index from EBVs — filter to a specific model and eval_number
#' pop <- pop |>
#'   get_table("ind_ebv") |>
#'   dplyr::filter(model == "blup_v1", eval_number == 1L) |>
#'   add_index("terminal")
#'
#' # Compute phenotypic index — filter to first record per individual per trait
#' pop <- pop |>
#'   get_table("ind_phenotype") |>
#'   dplyr::filter(pheno_number == 1L) |>
#'   add_index("terminal")
#'
#' # Compute true-value index from TBVs (value_col auto-detected)
#' pop <- pop |>
#'   get_table("ind_tbv") |>
#'   add_index("terminal")
#'
#' # Explicit value_col for a user-defined table
#' pop <- pop |>
#'   get_table("my_scores") |>
#'   add_index("terminal", value_col = "composite_score")
#'
#' # Replace previous run with a fresh one
#' pop <- pop |>
#'   get_table("ind_ebv") |>
#'   dplyr::filter(model == "blup_v2", eval_number == 1L) |>
#'   add_index("terminal", overwrite_index = TRUE)
#' }
#' @export
add_index <- function(tbl,
                      index_name,
                      value_col       = NULL,
                      overwrite_index = FALSE,
                      delete_all      = FALSE,
                      ...) {

  # ---- Input validation ----
  stopifnot(inherits(tbl, "tidybreed_table"))

  pop <- tbl$pop
  validate_tidybreed_pop(pop)

  stopifnot(is.character(index_name), length(index_name) == 1)
  validate_sql_identifier(index_name, what = "index name")

  stopifnot(is.logical(overwrite_index), length(overwrite_index) == 1)
  stopifnot(is.logical(delete_all), length(delete_all) == 1)

  extra_cols <- list(...)
  for (nm in names(extra_cols)) {
    if (length(extra_cols[[nm]]) != 1L) {
      stop(
        "Extra field '", nm, "' passed to add_index() must be a scalar ",
        "(broadcast to all rows). Use mutate_table() with a tibble for per-row values.",
        call. = FALSE
      )
    }
  }

  # ---- Resolve value_col ----
  known_value_cols <- c(
    ind_ebv       = "ebv_value",
    ind_phenotype = "pheno_value",
    ind_tbv       = "tbv_value"
  )

  if (is.null(value_col)) {
    auto <- known_value_cols[tbl$table_name]
    if (is.na(auto)) {
      stop(
        "Cannot auto-detect value_col for table '", tbl$table_name, "'. ",
        "Supply value_col explicitly (e.g. value_col = \"my_value\").",
        call. = FALSE
      )
    }
    value_col <- unname(auto)
  } else {
    stopifnot(is.character(value_col), length(value_col) == 1L)
    validate_sql_identifier(value_col, what = "value_col")
  }

  # ---- Validate required columns exist ----
  # ind_phenotype uses phenotype_name; all other value tables use trait_name
  tbl_fields <- DBI::dbListFields(pop$db_conn, tbl$table_name)
  name_col   <- if (tbl$table_name == "ind_phenotype" &&
                    "phenotype_name" %in% tbl_fields) "phenotype_name"
                else "trait_name"

  required_cols <- c("id_ind", name_col, value_col)
  missing_cols  <- setdiff(required_cols, tbl_fields)
  if (length(missing_cols) > 0) {
    stop(
      "Table '", tbl$table_name, "' is missing required column(s): ",
      paste(missing_cols, collapse = ", "), ". ",
      "add_index() requires id_ind, ", name_col, ", and the value column ('",
      value_col, "').",
      call. = FALSE
    )
  }

  # ---- Look up index weights ----
  idx_weights <- DBI::dbGetQuery(
    pop$db_conn,
    paste0(
      "SELECT trait_name, index_weight FROM index_meta WHERE index_name = '",
      gsub("'", "''", index_name), "' ORDER BY trait_name"
    )
  )

  if (nrow(idx_weights) == 0) {
    stop(
      "Index '", index_name, "' not found in index_meta. ",
      "Call define_index() to register it first.",
      call. = FALSE
    )
  }

  index_traits <- idx_weights$trait_name
  index_wts    <- setNames(idx_weights$index_weight, idx_weights$trait_name)

  # ---- Print weights message ----
  wt_str <- paste(
    vapply(index_traits, function(t) {
      sprintf("%s (wt=%.4g)", t, index_wts[t])
    }, character(1)),
    collapse = ", "
  )
  message("Computing index '", index_name, "': ", wt_str)

  # ---- Collect value rows ----
  if (length(tbl$pending_filter) == 0) {
    warning(
      "No filter applied to '", tbl$table_name, "'. ",
      "If this table can have multiple rows per (id_ind, trait_name) ",
      "(e.g. ind_phenotype with pheno_number > 1, or ind_ebv with multiple ",
      "eval_number values), filter first to ensure exactly one row per pair.",
      call. = FALSE
    )
  }
  raw_df   <- dplyr::collect(tbl$tbl)[, c("id_ind", name_col, value_col),
                                       drop = FALSE]
  # Normalise to 'trait_name' internally (index_meta always uses trait_name)
  value_df <- raw_df
  if (name_col != "trait_name") {
    names(value_df)[names(value_df) == name_col] <- "trait_name"
  }

  # ---- Filter to only index traits ----
  value_df <- value_df[value_df$trait_name %in% index_traits, , drop = FALSE]

  # ---- Check for duplicate (id_ind, trait_name) rows ----
  dup_check <- table(paste(value_df$id_ind, value_df$trait_name, sep = "\x1f"))
  if (any(dup_check > 1)) {
    n_dups <- sum(dup_check > 1)
    stop(
      n_dups, " (id_ind, ", name_col, ") combination(s) have more than one row. ",
      "Filter '", tbl$table_name, "' to a single row per (id_ind, ", name_col, ") ",
      "before calling add_index(). ",
      "Example: dplyr::filter(model == \"blup_v1\", eval_number == 1L) for ind_ebv, ",
      "or dplyr::filter(pheno_number == 1L) for ind_phenotype.",
      call. = FALSE
    )
  }

  # ---- Check that all index traits have values for all individuals ----
  individuals <- unique(value_df$id_ind)
  n_ind       <- length(individuals)

  if (n_ind == 0) {
    stop(
      "No rows found for index traits in '", tbl$table_name, "' after filtering.",
      call. = FALSE
    )
  }

  missing_info <- character(0)
  for (tr in index_traits) {
    tr_inds   <- value_df$id_ind[value_df$trait_name == tr]
    n_missing <- n_ind - length(tr_inds)
    if (n_missing > 0) {
      missing_info <- c(missing_info,
        sprintf("  trait '%s': %d of %d individuals missing values",
                tr, n_missing, n_ind))
    }
  }

  if (length(missing_info) > 0) {
    stop(
      "Cannot compute index '", index_name, "' — some individuals are missing ",
      "values for required traits in '", tbl$table_name, "':\n",
      paste(missing_info, collapse = "\n"), "\n",
      "Ensure all index traits have values for all individuals before calling ",
      "add_index().",
      call. = FALSE
    )
  }

  # ---- Pivot wide and compute index ----
  individuals_sorted <- sort(individuals)

  # Build n_ind × n_traits value matrix via direct lookup
  value_mat <- vapply(index_traits, function(tr) {
    tr_rows <- value_df[value_df$trait_name == tr, , drop = FALSE]
    tr_rows[[value_col]][match(individuals_sorted, tr_rows$id_ind)]
  }, numeric(length(individuals_sorted)))

  index_values <- as.numeric(value_mat %*% index_wts[index_traits])

  # ---- Determine index_number ----
  conn <- pop$db_conn

  if (delete_all) {
    DBI::dbExecute(conn, "DELETE FROM ind_index")
    new_index_num <- rep(1L, n_ind)

  } else if (overwrite_index) {
    DBI::dbExecute(
      conn,
      paste0("DELETE FROM ind_index WHERE index_name = '",
             gsub("'", "''", index_name), "'")
    )
    new_index_num <- rep(1L, n_ind)

  } else {
    # Query max index_number per individual for this index_name
    max_df <- DBI::dbGetQuery(
      conn,
      paste0(
        "SELECT id_ind, MAX(index_number) AS max_num FROM ind_index ",
        "WHERE index_name = '", gsub("'", "''", index_name), "' ",
        "GROUP BY id_ind"
      )
    )
    max_map <- setNames(max_df$max_num, max_df$id_ind)
    new_index_num <- vapply(individuals_sorted, function(id) {
      prev <- max_map[id]
      if (is.na(prev)) 1L else as.integer(prev) + 1L
    }, integer(1))
  }

  # ---- Build insertion data frame ----
  start_id  <- next_int_id(conn, "ind_index", "id_index")
  result_df <- data.frame(
    id_index     = seq.int(start_id, start_id + n_ind - 1L),
    id_ind       = individuals_sorted,
    index_name   = index_name,
    index_number = new_index_num,
    index_value  = index_values,
    stringsAsFactors = FALSE
  )

  # ---- Attach extra columns ----
  if (length(extra_cols) > 0) {
    prepped <- prepare_extra_cols(extra_cols, nrow(result_df), "ind_index", conn)
    for (nm in names(prepped)) result_df[[nm]] <- prepped[[nm]]
  }

  # ---- Insert ----
  DBI::dbWriteTable(conn, "ind_index", result_df, append = TRUE)

  run_num <- new_index_num[1]
  message("Added index '", index_name, "' (run #", run_num, ") for ",
          n_ind, " individuals")

  invisible(pop)
}
