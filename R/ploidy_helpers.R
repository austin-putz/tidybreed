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
#'   checks every individual in `ind_meta`.
#' @return Invisible `NULL` on success; errors otherwise.
#' @keywords internal
assert_ploidy_2 <- function(pop, ids = NULL) {
  if (is.null(ids)) {
    where_sql <- "WHERE ploidy != 2"
  } else {
    where_sql <- paste0(
      "WHERE id_ind IN (", sql_in_list(ids, what = "individual ID"), ") ",
      "AND ploidy != 2"
    )
  }

  n_bad <- DBI::dbGetQuery(
    pop$db_conn,
    paste0("SELECT COUNT(*) AS n FROM ind_meta ", where_sql)
  )$n

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
