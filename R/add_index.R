#' Compute a selection index from EBVs
#'
#' @description
#' Calculates a named selection index by multiplying each individual's EBVs by
#' the index weights defined in [define_index()] and summing the products.
#' Results are appended to the `ind_index` table with an auto-incrementing
#' `index_number` so that successive runs are preserved.
#'
#' The first argument **must** be a `tidybreed_table` obtained from
#' `get_table("ind_ebv")` (optionally filtered). This is the only function in
#' the package that requires piping through `ind_ebv` rather than `ind_meta`.
#'
#' ## EBV resolution
#' If no filter is applied, `add_index()` automatically selects the row with the
#' highest `eval_number` per `(id_ind, trait_name)` and issues a warning. After
#' this auto-selection, if any `(id_ind, trait_name)` pair still has more than
#' one row (e.g. because multiple models are present), an error is thrown — you
#' must filter down to a single model first.
#'
#' ## Completeness requirement
#' Every individual must have an EBV for **every** trait in the index. If any
#' individual is missing an EBV for any required trait, an error is thrown
#' before any rows are written.
#'
#' @param tbl A `tidybreed_table` from `get_table("ind_ebv")` (optionally
#'   filtered). Must be the `ind_ebv` table.
#' @param index_name Character scalar. Name of the index to compute; must
#'   already exist in `index_meta` (created via [define_index()]).
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
#' @seealso [define_index()], [get_table()], [add_ebv()]
#'
#' @examples
#' \dontrun{
#' # Compute index using latest EBVs (warning issued — no filter applied)
#' pop <- pop |>
#'   get_table("ind_ebv") |>
#'   add_index("terminal")
#'
#' # Filter to a specific model and eval_number before computing (no warning)
#' pop <- pop |>
#'   get_table("ind_ebv") |>
#'   dplyr::filter(model == "blup_v1", eval_number == 1L) |>
#'   add_index("terminal")
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
                      overwrite_index = FALSE,
                      delete_all      = FALSE,
                      ...) {

  # ---- Input validation ----
  stopifnot(inherits(tbl, "tidybreed_table"))
  if (tbl$table_name != "ind_ebv") {
    stop(
      "add_index() requires a tidybreed_table from get_table(\"ind_ebv\"). ",
      "Got table '", tbl$table_name, "' instead. ",
      "Pipe through get_table(\"ind_ebv\") (and optionally filter()) first.",
      call. = FALSE
    )
  }

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
        "(broadcast to all rows). Use mutate_table() for per-row vectors.",
        call. = FALSE
      )
    }
  }

  # ---- Look up index weights ----
  idx_weights <- DBI::dbGetQuery(
    pop$db_conn,
    paste0(
      "SELECT trait_name, index_wt FROM index_meta WHERE index_name = '",
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
  index_wts    <- setNames(idx_weights$index_wt, idx_weights$trait_name)

  # ---- Print weights message ----
  wt_str <- paste(
    vapply(index_traits, function(t) {
      sprintf("%s (wt=%.4g)", t, index_wts[t])
    }, character(1)),
    collapse = ", "
  )
  message("Computing index '", index_name, "': ", wt_str)

  # ---- Resolve which EBV rows to use ----
  no_filter <- length(tbl$pending_filter) == 0

  if (no_filter) {
    warning(
      "No filter applied to ind_ebv. Assuming the latest eval_number per ",
      "individual per trait. Supply a filter (e.g. dplyr::filter(model == ..., ",
      "eval_number == ...)) if this is not intended.",
      call. = FALSE
    )
    # Auto-select latest eval_number per (id_ind, trait_name)
    ebv_df <- DBI::dbGetQuery(
      pop$db_conn,
      "SELECT e.id_ind, e.trait_name, e.ebv
       FROM ind_ebv e
       INNER JOIN (
         SELECT id_ind, trait_name, MAX(eval_number) AS max_eval
         FROM ind_ebv
         GROUP BY id_ind, trait_name
       ) m ON e.id_ind = m.id_ind
          AND e.trait_name = m.trait_name
          AND e.eval_number = m.max_eval"
    )
  } else {
    ebv_df <- dplyr::collect(tbl$tbl)[, c("id_ind", "trait_name", "ebv")]
  }

  # ---- Filter to only index traits ----
  ebv_df <- ebv_df[ebv_df$trait_name %in% index_traits, , drop = FALSE]

  # ---- Check for duplicate (id_ind, trait_name) rows ----
  dup_check <- table(paste(ebv_df$id_ind, ebv_df$trait_name, sep = "\x1f"))
  if (any(dup_check > 1)) {
    n_dups <- sum(dup_check > 1)
    stop(
      n_dups, " (id_ind, trait_name) combination(s) still have more than one ",
      "EBV row after resolution. Filter ind_ebv to a single model and/or ",
      "eval_number before calling add_index(). ",
      "Example: dplyr::filter(model == \"blup_v1\", eval_number == 1L)",
      call. = FALSE
    )
  }

  # ---- Check that all index traits have EBVs for all individuals ----
  individuals <- unique(ebv_df$id_ind)
  n_ind       <- length(individuals)

  if (n_ind == 0) {
    stop("No EBV rows found for index traits after filtering.", call. = FALSE)
  }

  missing_info <- character(0)
  for (tr in index_traits) {
    tr_inds    <- ebv_df$id_ind[ebv_df$trait_name == tr]
    n_missing  <- n_ind - length(tr_inds)
    if (n_missing > 0) {
      missing_info <- c(missing_info,
        sprintf("  trait '%s': %d of %d individuals missing EBVs",
                tr, n_missing, n_ind))
    }
  }

  if (length(missing_info) > 0) {
    stop(
      "Cannot compute index '", index_name, "' — some individuals are missing ",
      "EBVs for required traits:\n",
      paste(missing_info, collapse = "\n"), "\n",
      "Ensure all index traits have EBVs for all individuals before calling ",
      "add_index().",
      call. = FALSE
    )
  }

  # ---- Pivot wide and compute index ----
  individuals_sorted <- sort(individuals)

  # Build n_ind × n_traits EBV matrix via direct lookup — avoids stats::reshape()
  # which behaves differently on tibbles (from dplyr::collect) vs plain data.frames,
  # causing column name mangling that makes wide[[tr]] return NULL.
  ebv_mat <- vapply(index_traits, function(tr) {
    tr_rows <- ebv_df[ebv_df$trait_name == tr, , drop = FALSE]
    tr_rows$ebv[match(individuals_sorted, tr_rows$id_ind)]
  }, numeric(length(individuals_sorted)))

  index_values <- as.numeric(ebv_mat %*% index_wts[index_traits])

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
  result_df <- data.frame(
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
