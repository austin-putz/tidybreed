#' Materialize genotype dosage values into `ind_genotype`
#'
#' @description
#' Computes 0/1/2 genotype dosages (sum of alleles across an individual's
#' haplotype strands) from the long `ind_haplotype` table and writes them to the
#' on-demand `ind_genotype` cache. `ind_genotype` is **never** auto-populated by
#' [add_founders()] or [add_offspring()] — call `add_dosage()` explicitly when
#' you need dosage values for marker-assisted selection, allele-frequency
#' queries, or other downstream analysis.
#'
#' Pipe a `tidybreed_table` (from [get_table()] and optionally [filter()]) as the
#' first argument to select individuals. As with [add_tbv()], the **distinct**
#' `id_ind` values in the collected table form the candidate set, so a table with
#' multiple rows per individual (e.g. `ind_phenotype`) does not multiply work.
#'
#' @section `add_dosage()` vs. [add_genotypes()]:
#' These are different operations despite similar names.
#' [add_genotypes()] marks animals as *physically genotyped* on a chip by
#' writing a BOOLEAN `has_<chip>` column to `ind_meta`; it touches no dosage
#' data. `add_dosage()` **materializes simulated dosage values** (ground truth
#' from `ind_haplotype`) into `ind_genotype`.
#'
#' @param tbl A `tidybreed_table` from [get_table()] (optionally filtered). Any
#'   table with an `id_ind` column is accepted.
#' @param chip_name Character or `NULL`. Name of a chip defined via
#'   [define_chip()] (which writes `is_<chip_name>` to `genome_meta`). Loci are
#'   restricted to `is_<chip_name> = TRUE`. Errors if the column is missing.
#' @param locus_names Character vector or `NULL`. Explicit loci to materialize.
#'   Takes precedence over `chip_name` when both are supplied.
#' @param overwrite_dosage Logical. When `TRUE`, existing `ind_genotype` rows for
#'   the candidate individuals are deleted before inserting (cache-scope reset).
#'   When `FALSE` (default), rows are upserted via `INSERT OR REPLACE` — dosage
#'   is fully determined by the haplotypes, so re-running is idempotent.
#'
#' @return The modified `tidybreed_pop` (invisibly).
#'
#' @seealso [add_genotypes()], [extract_genotypes()], [define_chip()]
#'
#' @examples
#' \dontrun{
#' # Dosage for all generation-5 candidates on the 50K chip
#' pop <- pop |>
#'   get_table("ind_meta") |>
#'   dplyr::filter(gen == 5L) |>
#'   add_dosage(chip_name = "50K")
#'
#' # Dosage at specific QTL only
#' pop <- pop |>
#'   get_table("ind_meta") |>
#'   add_dosage(locus_names = c("Locus_1", "Locus_42"))
#'
#' # Marker-assisted selection: homozygous-favorable at a QTL
#' ids <- pop |>
#'   get_table("ind_genotype") |>
#'   dplyr::filter(locus_name == "Locus_1", dosage_value == 2L) |>
#'   dplyr::pull(id_ind)
#' }
#' @export
add_dosage <- function(tbl, chip_name = NULL, locus_names = NULL,
                       overwrite_dosage = FALSE) {

  stopifnot(inherits(tbl, "tidybreed_table"))
  pop <- tbl$pop
  validate_tidybreed_pop(pop)
  assert_diploid_only(pop)

  # --- Resolve candidate individuals (distinct id_ind) ---
  if (length(tbl$pending_filter) == 0) {
    subset_ids <- NULL
  } else {
    collected <- dplyr::collect(tbl)
    if (!"id_ind" %in% names(collected)) {
      stop("Filtered table '", tbl$table_name,
           "' must contain 'id_ind' to select individuals for add_dosage().",
           call. = FALSE)
    }
    subset_ids <- unique(collected[["id_ind"]])
    if (length(subset_ids) == 0L) {
      warning("No individuals matched; no dosage computed.", call. = FALSE)
      return(invisible(pop))
    }
  }

  # --- Resolve locus filter (locus_names takes precedence over chip_name) ---
  locus_filter_sql <- ""
  if (!is.null(locus_names)) {
    stopifnot(is.character(locus_names), length(locus_names) >= 1L)
    name_sql <- sql_in_list(locus_names, what = "locus name")
    locus_filter_sql <- paste0(" AND h.locus_name IN (", name_sql, ")")
  } else if (!is.null(chip_name)) {
    stopifnot(is.character(chip_name), length(chip_name) == 1L, nzchar(chip_name))
    chip_col <- paste0("is_", chip_name)
    validate_sql_identifier(chip_col, what = "chip name")
    if (!chip_col %in% DBI::dbListFields(pop$db_conn, "genome_meta")) {
      stop("Chip '", chip_name, "' not found in genome_meta. ",
           "Call define_chip() first.", call. = FALSE)
    }
    locus_filter_sql <- paste0(
      " AND h.locus_name IN (SELECT locus_name FROM genome_meta WHERE ",
      chip_col, " = TRUE)")
  }

  id_filter_sql <- ""
  if (!is.null(subset_ids)) {
    ids_sql <- sql_in_list(subset_ids, what = "individual ID")
    id_filter_sql <- paste0(" AND h.id_ind IN (", ids_sql, ")")
  }

  # --- Optional cache-scope reset for the candidate individuals ---
  if (isTRUE(overwrite_dosage)) {
    if (is.null(subset_ids)) {
      DBI::dbExecute(pop$db_conn, "DELETE FROM ind_genotype")
    } else {
      ids_sql <- sql_in_list(subset_ids, what = "individual ID")
      DBI::dbExecute(pop$db_conn,
        paste0("DELETE FROM ind_genotype WHERE id_ind IN (", ids_sql, ")"))
    }
  }

  # --- Materialize dosage. INSERT OR REPLACE relies on the (id_ind, locus_id)
  #     primary key, so re-running is idempotent. ---
  n_before <- DBI::dbGetQuery(pop$db_conn,
    "SELECT COUNT(*) AS n FROM ind_genotype")$n
  DBI::dbExecute(pop$db_conn, paste0(
    "INSERT OR REPLACE INTO ind_genotype (id_ind, locus_id, locus_name, dosage_value) ",
    "SELECT h.id_ind, h.locus_id, h.locus_name, CAST(SUM(h.allele) AS UTINYINT) ",
    "FROM ind_haplotype h ",
    "WHERE 1 = 1", id_filter_sql, locus_filter_sql, " ",
    "GROUP BY h.id_ind, h.locus_id, h.locus_name"
  ))
  n_after <- DBI::dbGetQuery(pop$db_conn,
    "SELECT COUNT(*) AS n FROM ind_genotype")$n

  message("Materialized dosage into ind_genotype (", n_after,
          " total rows; ", n_after - n_before, " new).")

  invisible(pop)
}
