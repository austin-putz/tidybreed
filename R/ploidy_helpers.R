#' Assert that `chr_meta` has no non-diploid chromosomes
#'
#' Ploidy-aware dosage and genotype extraction is Stage 4 (`define_chr()`)
#' scope. Until that ships, any `chr_meta` row with `copies_M`/`copies_F` other
#' than 2 would silently produce wrong dosage sums rather than erroring, since
#' the current dosage/export SQL assumes exactly two haplotype rows per
#' individual per locus. `define_chr()` does not exist yet, so this guard is
#' forward defense — it protects the day `chr_meta` gains non-default rows,
#' not a fix for a currently reachable bug.
#'
#' @param pop A `tidybreed_pop`.
#' @return Invisible `NULL` on success; errors otherwise.
#' @keywords internal
assert_diploid_only <- function(pop) {
  if (!"chr_meta" %in% pop$tables) {
    return(invisible(NULL))
  }

  n_non_diploid <- DBI::dbGetQuery(
    pop$db_conn,
    "SELECT COUNT(*) AS n FROM chr_meta WHERE copies_M != 2 OR copies_F != 2"
  )$n

  if (n_non_diploid > 0) {
    stop(
      "chr_meta contains ", n_non_diploid, " chromosome(s) with non-diploid ",
      "copy number (copies_M != 2 or copies_F != 2). Ploidy-aware dosage and ",
      "genotype extraction are not supported until Stage 4 (define_chr()) ",
      "ships; this function assumes diploid autosomes only.",
      call. = FALSE
    )
  }

  invisible(NULL)
}
