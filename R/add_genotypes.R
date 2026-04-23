#' Mark animals as genotyped on a SNP chip
#'
#' @description
#' Records which individuals have been genotyped on a named SNP chip by
#' writing or updating a BOOLEAN column `has_<chip_name>` in `ind_meta`.
#' Pipe a `tidybreed_table` (from [get_table()] and optionally [filter()]) as
#' the first argument to restrict which animals are marked.
#'
#' The operation is **additive**: animals already marked `TRUE` remain `TRUE`.
#' Only new animals are flipped. This mirrors real life — once an animal is
#' genotyped it stays genotyped.
#'
#' @param tbl A `tidybreed_table` object from [get_table()] (optionally piped
#'   through [filter()]). The table must contain an `id_ind` column when a
#'   filter is applied.
#' @param chip_name Character. Name of an existing SNP chip (must have an
#'   `is_<chip_name>` column in `genome_meta`, created by [define_chip()]).
#' @param col_name Character. Name of the BOOLEAN column to write in
#'   `ind_meta`. Default: `paste0("has_", chip_name)`.
#'
#' @return The modified `tidybreed_pop` (invisibly).
#'
#' @examples
#' \dontrun{
#' # Genotype all females in generation 1 on the 50k chip
#' pop <- pop |>
#'   get_table("ind_meta") |>
#'   dplyr::filter(sex == "F", gen == 1L) |>
#'   add_genotypes("50k")
#'
#' # Also genotype all generation 2 animals (additive — gen 1 females stay TRUE)
#' pop <- pop |>
#'   get_table("ind_meta") |>
#'   dplyr::filter(gen == 2L) |>
#'   add_genotypes("50k")
#' }
#' @export
add_genotypes <- function(tbl,
                         chip_name,
                         col_name = paste0("has_", chip_name)) {

  stopifnot(inherits(tbl, "tidybreed_table"))
  pop <- tbl$pop
  validate_tidybreed_pop(pop)
  stopifnot(is.character(chip_name), length(chip_name) == 1L, nzchar(chip_name))
  validate_sql_identifier(col_name, what = "col_name")

  chip_col <- paste0("is_", chip_name)
  if (!chip_col %in% DBI::dbListFields(pop$db_conn, "genome_meta")) {
    stop("Chip '", chip_name, "' not found in genome_meta. ",
         "Call define_chip() first.", call. = FALSE)
  }

  if (!"ind_meta" %in% pop$tables) {
    stop("ind_meta does not exist. Call add_founders() first.", call. = FALSE)
  }

  if (length(tbl$pending_filter) == 0) {
    subset_ids <- NULL
  } else {
    collected <- dplyr::collect(tbl)
    if (!"id_ind" %in% names(collected)) {
      stop("Filtered table '", tbl$table_name,
           "' must contain 'id_ind' to subset individuals for genotyping.",
           call. = FALSE)
    }
    subset_ids <- unique(collected[["id_ind"]])
  }

  # Ensure has_<chip_name> column exists (DEFAULT FALSE keeps existing TRUEs intact)
  if (!col_name %in% DBI::dbListFields(pop$db_conn, "ind_meta")) {
    DBI::dbExecute(
      pop$db_conn,
      paste0("ALTER TABLE ind_meta ADD COLUMN ", col_name, " BOOLEAN DEFAULT FALSE")
    )
  }

  # Targeted UPDATE — additive; never resets existing TRUEs
  if (is.null(subset_ids)) {
    DBI::dbExecute(
      pop$db_conn,
      paste0("UPDATE ind_meta SET ", col_name, " = TRUE")
    )
  } else {
    ids_sql <- paste0("'", subset_ids, "'", collapse = ", ")
    DBI::dbExecute(
      pop$db_conn,
      paste0("UPDATE ind_meta SET ", col_name, " = TRUE WHERE id_ind IN (", ids_sql, ")")
    )
  }

  n_genotyped <- DBI::dbGetQuery(
    pop$db_conn,
    paste0("SELECT COUNT(*) AS n FROM ind_meta WHERE ", col_name, " = TRUE")
  )$n

  message("Chip '", chip_name, "': ", n_genotyped, " animal(s) now marked as genotyped.")
  invisible(pop)
}
