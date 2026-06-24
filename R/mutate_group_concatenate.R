#' Create a composite group column by concatenating existing columns
#'
#' @description
#' Deterministically builds a `VARCHAR` group label for each filtered row by
#' concatenating the values of one or more existing columns in the same table,
#' separated by `sep`. Useful for contemporary groups (`sex_line_gen`), sib
#' groups (`sire_id_dam_id`), and any label derived from existing row data.
#'
#' Pipe `get_table("ind_meta")` — optionally through `dplyr::filter()` — into
#' this function. All `source_columns` must exist in the same table as the
#' target column. No randomness is involved; results are fully deterministic.
#'
#' @section NULL handling:
#' The `null_action` argument controls what happens when any source column
#' value is `NULL` (`NA`) for a given row:
#' - `"propagate"` (default): any `NULL` source produces a `NULL` result for
#'   that row. Safe for sib-group labels where a missing parent ID should
#'   produce no group, not a partial label.
#' - `"skip"`: `NULL` sources are silently omitted; remaining non-`NULL`
#'   values are concatenated. A row where all sources are `NULL` still yields
#'   a `NULL` result.
#' - `"literal"`: `NULL` sources are replaced with the string `"NA"` before
#'   concatenation.
#'
#' @param tbl A `tidybreed_table` from `get_table()` (optionally filtered).
#' @param col_name Character scalar. Name of the group column to create or
#'   update. Must be a valid SQL identifier and not a reserved column.
#' @param source_columns Character vector of one or more column names in the
#'   same table to concatenate. Order is preserved.
#' @param sep Character scalar. Separator string. Default `"_"`.
#' @param null_action One of `"propagate"` (default), `"skip"`, or
#'   `"literal"`. See section NULL handling above.
#' @param overwrite Logical. If `FALSE` (default), error when any filtered row
#'   already has a non-`NULL` value in `col_name`.
#' @param quiet Logical. If `TRUE`, suppress informational messages.
#'
#' @return The modified `tidybreed_pop` object (invisibly).
#'
#' @seealso [mutate_group_seq()], [mutate_group_named()], [mutate_table()]
#'
#' @examples
#' \dontrun{
#' # Contemporary group: sex × line_name
#' pop <- pop |>
#'   get_table("ind_meta") |>
#'   mutate_group_concatenate(
#'     col_name       = "cg",
#'     source_columns = c("sex", "line_name")
#'   )
#'
#' # Sib group from parent IDs — NULL propagates when a parent is unknown
#' pop <- pop |>
#'   get_table("ind_meta") |>
#'   mutate_group_concatenate(
#'     col_name       = "sib_group",
#'     source_columns = c("id_parent_1", "id_parent_2"),
#'     null_action    = "propagate"
#'   )
#'
#' # Generation × sex label, NULLs as literal "NA"
#' pop <- pop |>
#'   get_table("ind_meta") |>
#'   dplyr::filter(gen == 1L) |>
#'   mutate_group_concatenate(
#'     col_name       = "gen_sex",
#'     source_columns = c("gen", "sex"),
#'     sep            = ":",
#'     null_action    = "literal"
#'   )
#' }
#' @export
mutate_group_concatenate <- function(tbl,
                                      col_name,
                                      source_columns,
                                      sep         = "_",
                                      null_action = "propagate",
                                      overwrite   = FALSE,
                                      quiet       = FALSE) {

  # --- Validate structure ---
  info <- .resolve_group_target(tbl, col_name, "VARCHAR", overwrite)
  if (is.null(info)) return(invisible(tbl$pop))

  pop        <- info$pop
  conn       <- info$conn
  table_name <- info$table_name
  pk_col     <- info$pk_col
  pks        <- info$pks
  n_total    <- info$n_total

  null_action <- match.arg(null_action, c("propagate", "skip", "literal"))

  if (!is.character(sep) || length(sep) != 1L || is.na(sep)) {
    stop("'sep' must be a single non-NA character string.", call. = FALSE)
  }
  if (!is.character(source_columns) || length(source_columns) < 1L ||
      any(is.na(source_columns))) {
    stop("'source_columns' must be a non-empty character vector with no NA values.",
         call. = FALSE)
  }
  if (col_name %in% source_columns) {
    stop(
      "'col_name' ('", col_name, "') cannot be one of the 'source_columns'. ",
      "Source and target columns must be distinct.",
      call. = FALSE
    )
  }
  existing_cols <- DBI::dbListFields(conn, table_name)
  missing_cols  <- setdiff(source_columns, existing_cols)
  if (length(missing_cols) > 0L) {
    stop(
      "Source column", if (length(missing_cols) != 1L) "s" else "",
      " not found in '", table_name, "': ",
      paste(missing_cols, collapse = ", "),
      call. = FALSE
    )
  }

  # --- Collect source data for filtered rows (sorted by pk_col) ---
  collected <- dplyr::collect(
    dplyr::select(tbl$tbl, dplyr::all_of(c(pk_col, source_columns)))
  )
  collected <- collected[order(collected[[pk_col]]), ]

  src_df <- collected[, source_columns, drop = FALSE]

  # --- Build concatenated values in R ---
  values <- vapply(seq_len(nrow(src_df)), function(i) {
    raw <- unlist(src_df[i, ], use.names = FALSE)
    is_null <- is.na(raw)

    if (null_action == "propagate" && any(is_null)) return(NA_character_)

    row_chr <- as.character(raw)
    if (null_action == "skip") {
      row_chr <- row_chr[!is_null]
    } else if (null_action == "literal") {
      row_chr[is_null] <- "NA"
    }

    if (length(row_chr) == 0L) return(NA_character_)
    paste(row_chr, collapse = sep)
  }, character(1))

  # --- Write and return ---
  .ensure_col_exists(conn, table_name, col_name, "VARCHAR")
  .write_group_values(conn, table_name, col_name, pk_col, collected[[pk_col]], values)
  .summarize_group_assignment(values, n_total, table_name, col_name, quiet)

  invisible(pop)
}
