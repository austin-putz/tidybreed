#' Assert that a set of individuals all have ploidy = 2
#'
#' Ploidy-aware dosage and genotype extraction for real polyploidy (ploidy > 2)
#' is not yet supported. Sex-linked and organelle `chr_inheritance` rules
#' (reduced `from_parent_1`/`from_parent_2` copy counts) are fine here —
#' `SUM(allele)` over however many `ind_haplotype` rows exist per locus already
#' produces correct dosage regardless of row count, so this guard targets actual
#' organism ploidy, not chromosome copy-count shape.
#'
#' @param pop A `tidybreed_pop`.
#' @param ids Optional character vector of `id_ind` to check; `NULL` (default)
#'   checks every individual in `ind_meta`, and `character(0)` checks nobody.
#' @return Invisible `NULL` on success; errors otherwise.
#' @keywords internal
assert_ploidy_2 <- function(pop, ids = NULL) {
  count_sql <- "SELECT COUNT(*) AS n FROM ind_meta WHERE ploidy != 2"
  n_bad <- if (is.null(ids)) {
    DBI::dbGetQuery(pop$db_conn, count_sql)$n
  } else {
    .with_id_view(pop$db_conn, ids, function(view) {
      DBI::dbGetQuery(pop$db_conn, paste0(
        count_sql, " AND id_ind IN (SELECT id_ind FROM ", view, ")"))$n
    })
  }

  if (n_bad > 0) {
    stop(
      n_bad, " individual(s) have ploidy != 2. Ploidy-aware dosage and ",
      "genotype extraction for real polyploidy are not supported. ",
      "Sex-linked/organelle chr_inheritance rules (reduced from_parent_1/",
      "from_parent_2 copy counts) are fine here and do not trigger this error.",
      call. = FALSE
    )
  }

  invisible(NULL)
}
