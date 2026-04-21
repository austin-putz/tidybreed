#' Extract genotype data for genotyped animals on a SNP chip
#'
#' @description
#' Returns a tibble of genotypes (0/1/2 encoding) for animals that have been
#' marked as genotyped via [add_genotypes()], restricted to loci on the named
#' chip. Follows the same `filter()` -> action pipe pattern as [add_phenotype()]:
#' pipe through `dplyr::filter()` to further restrict which animals are included.
#'
#' The returned set is the **intersection** of:
#' * Animals with `has_<chip_name> == TRUE` in `ind_meta`
#' * Animals matching any pending `filter()` predicates
#' * Loci with `is_<chip_name> == TRUE` in `genome_meta`
#'
#' @param pop A `tidybreed_pop` object.
#' @param chip_name Character. Name of a chip previously defined via
#'   [define_chip()] and applied to animals via [add_genotypes()].
#' @param col_name Character. Name of the BOOLEAN column in `ind_meta` that
#'   records genotyping status. Default: `paste0("has_", chip_name)`.
#'
#' @return A tibble with one row per genotyped animal and one column per chip
#'   locus (`id_ind` + `locus_N` columns). Locus columns use 0/1/2 integer
#'   encoding and are ordered by `locus_id`.
#'
#' @examples
#' \dontrun{
#' # All genotyped animals on the 50k chip
#' geno <- extract_genotypes(pop, "50k")
#'
#' # Only females genotyped on the HD chip
#' geno <- pop |>
#'   dplyr::filter(sex == "F") |>
#'   extract_genotypes("HD")
#' }
#' @export
extract_genotypes <- function(pop,
                             chip_name,
                             col_name = paste0("has_", chip_name)) {

  stopifnot(inherits(pop, "tidybreed_pop"))
  validate_tidybreed_pop(pop)
  stopifnot(is.character(chip_name), length(chip_name) == 1L, nzchar(chip_name))
  validate_sql_identifier(col_name, what = "col_name")

  chip_col <- paste0("is_", chip_name)
  if (!chip_col %in% DBI::dbListFields(pop$db_conn, "genome_meta")) {
    stop("Chip '", chip_name, "' not found in genome_meta. ",
         "Call define_chip() first.", call. = FALSE)
  }

  if (!col_name %in% DBI::dbListFields(pop$db_conn, "ind_meta")) {
    stop("Column '", col_name, "' not found in ind_meta. ",
         "Call add_genotypes('", chip_name, "') first.", call. = FALSE)
  }

  resolved   <- resolve_pending_filter(pop)
  pop        <- resolved$pop
  subset_ids <- resolved$ids

  # Resolve genotyped animal IDs, intersected with any pending filter
  if (is.null(subset_ids)) {
    has_ids <- DBI::dbGetQuery(
      pop$db_conn,
      paste0("SELECT id_ind FROM ind_meta WHERE ", col_name, " = TRUE")
    )$id_ind
  } else {
    ids_sql <- paste0("'", subset_ids, "'", collapse = ", ")
    has_ids <- DBI::dbGetQuery(
      pop$db_conn,
      paste0("SELECT id_ind FROM ind_meta WHERE ", col_name, " = TRUE ",
             "AND id_ind IN (", ids_sql, ")")
    )$id_ind
  }

  if (length(has_ids) == 0) {
    stop("No genotyped animals found for chip '", chip_name,
         "' in the filtered set. Call add_genotypes() first.", call. = FALSE)
  }

  # Resolve chip locus IDs -> genome_genotype column names
  chip_locus_ids <- DBI::dbGetQuery(
    pop$db_conn,
    paste0("SELECT locus_id FROM genome_meta WHERE ", chip_col,
           " = TRUE ORDER BY locus_id")
  )$locus_id

  if (length(chip_locus_ids) == 0) {
    stop("No loci found on chip '", chip_name, "'.", call. = FALSE)
  }

  chip_cols <- paste0("locus_", chip_locus_ids)

  # Dynamic SELECT: id_ind + chip locus columns, filtered to genotyped animals
  cols_sql <- paste(c("id_ind", chip_cols), collapse = ", ")
  ids_sql  <- paste0("'", has_ids, "'", collapse = ", ")

  df <- DBI::dbGetQuery(
    pop$db_conn,
    paste0("SELECT ", cols_sql,
           " FROM genome_genotype WHERE id_ind IN (", ids_sql, ")")
  )

  tibble::as_tibble(df)
}
