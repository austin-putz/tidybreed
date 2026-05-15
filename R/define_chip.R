#' Define a SNP chip
#'
#' @description
#' Marks the loci in a filtered `genome_meta` table as members of a SNP chip,
#' creating a `BOOLEAN` column (default: `is_{chip_name}`) in `genome_meta`.
#'
#' Pipe `get_table("genome_meta")` — optionally piped through `dplyr::filter()`
#' — into this function. The unique `locus_id` values in the collected result
#' determine chip membership. All other loci receive `FALSE`.
#'
#' This replaces the old positional selection methods (`n`/`method`, `locus_tf`,
#' `locus_ids`, `locus_names`). Selection is now done entirely by the caller
#' using standard dplyr verbs on the `genome_meta` table, which is order-safe
#' and compatible with dynamic genomes.
#'
#' @param tbl A `tidybreed_table` from `get_table("genome_meta")` (optionally
#'   piped through `dplyr::filter()`). Must contain a `locus_id` column.
#' @param chip_name Character scalar. Name of the SNP chip; used in messages
#'   and to derive the default column name.
#' @param col_name Character scalar. Column name created in `genome_meta`.
#'   Default: `paste0("is_", chip_name)`. Must be a valid SQL identifier.
#'
#' @return The modified `tidybreed_pop` object (invisibly).
#'
#' @seealso [define_qtl()], [add_genotypes()], [extract_genotypes()]
#'
#' @examples
#' \dontrun{
#' # All loci on chromosomes 1-5
#' pop <- pop |>
#'   get_table("genome_meta") |>
#'   dplyr::filter(chr %in% 1:5) |>
#'   define_chip("chr1to5")
#'
#' # Random chip — sample locus names, then filter
#' selected <- pop |>
#'   get_table("genome_meta") |>
#'   dplyr::collect() |>
#'   dplyr::slice_sample(n = 500) |>
#'   dplyr::pull(locus_name)
#'
#' pop <- pop |>
#'   get_table("genome_meta") |>
#'   dplyr::filter(locus_name %in% selected) |>
#'   define_chip("50K")
#'
#' # Complement of an existing chip
#' pop <- pop |>
#'   get_table("genome_meta") |>
#'   dplyr::filter(is_50K == FALSE) |>
#'   define_chip("non50K")
#' }
#' @export
define_chip <- function(tbl, chip_name, col_name = paste0("is_", chip_name)) {

  stopifnot(inherits(tbl, "tidybreed_table"))
  pop <- tbl$pop
  validate_tidybreed_pop(pop)

  stopifnot(is.character(chip_name), length(chip_name) == 1L, nchar(chip_name) > 0L)
  validate_sql_identifier(col_name, what = "col_name")

  if (!"genome_meta" %in% pop$tables) {
    stop("genome_meta table does not exist. Call initialize_genome() first.",
         call. = FALSE)
  }

  # Collect filtered locus IDs
  collected <- dplyr::collect(tbl)
  if (!"locus_id" %in% names(collected)) {
    stop(
      "The filtered table must contain 'locus_id'. ",
      "Pipe get_table(\"genome_meta\") into define_chip().",
      call. = FALSE
    )
  }
  if (nrow(collected) == 0L) {
    stop(
      "No loci selected — your filter returned zero rows. ",
      "Check the filter before piping into define_chip().",
      call. = FALSE
    )
  }

  selected_ids <- collected$locus_id

  # Guard: error if column already exists
  existing_cols <- DBI::dbListFields(pop$db_conn, "genome_meta")
  if (col_name %in% existing_cols) {
    stop(
      "Column '", col_name, "' already exists in genome_meta. ",
      "Use a different chip_name or drop the column first.",
      call. = FALSE
    )
  }

  # Add column defaulting FALSE, then flip selected loci to TRUE
  DBI::dbExecute(
    pop$db_conn,
    paste0("ALTER TABLE genome_meta ADD COLUMN ", col_name, " BOOLEAN DEFAULT FALSE")
  )
  ids_sql <- paste(selected_ids, collapse = ", ")
  DBI::dbExecute(
    pop$db_conn,
    paste0("UPDATE genome_meta SET ", col_name, " = TRUE WHERE locus_id IN (", ids_sql, ")")
  )

  message(
    "Defined chip '", chip_name, "' with ", length(selected_ids),
    " SNPs in column '", col_name, "'"
  )

  invisible(pop)
}
