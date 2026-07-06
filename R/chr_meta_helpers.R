#' Resolve a chromosome's copy count for a given sex and ploidy
#'
#' Translates `chr_meta`'s relative `copy_mode` (`"full"`/`"half"`/`"none"`)
#' into an absolute copy count for one individual, given their own ploidy.
#' Shared by [add_founders()] and [add_offspring()] so both agree on expected
#' row counts per chromosome per sex — duplicating this mapping inline in both
#' files would risk them silently disagreeing.
#'
#' @param copy_mode Character scalar: `"full"`, `"half"`, or `"none"`.
#' @param ploidy Integer scalar. The individual's own ploidy (always `2L` in
#'   this version of tidybreed).
#' @return Integer copy count (`0L`, `1L`, or `2L` in this version).
#' @keywords internal
resolve_chr_copy_count <- function(copy_mode, ploidy = 2L) {
  switch(
    copy_mode,
    full = as.integer(ploidy),
    half = as.integer(ploidy %/% 2L),
    none = 0L,
    stop("Invalid copy_mode: '", copy_mode,
         "'. Must be 'full', 'half', or 'none'.", call. = FALSE)
  )
}


#' Load `chr_meta` as a named list keyed by `chr_name`
#'
#' Convenience reader used by [add_founders()] and [add_offspring()] to avoid
#' re-querying `chr_meta` per chromosome inside a loop.
#'
#' @param conn A DBI connection.
#' @return Named list (by `chr_name`), each element a 1-row data.frame with
#'   `copy_mode_M`, `copy_mode_F`, `hemi_parent`, `recombines_M`, `recombines_F`.
#' @keywords internal
get_chr_meta_map <- function(conn) {
  df <- DBI::dbGetQuery(conn, "SELECT * FROM chr_meta")
  stats::setNames(
    lapply(seq_len(nrow(df)), function(i) df[i, , drop = FALSE]),
    df$chr_name
  )
}


#' Error if any of the given loci sit on a non-"full" chromosome
#'
#' `rescale_effects_to_target()`'s Falconer variance formula
#' (`2 * p * (1-p) * a^2`) assumes every QTL is diploid/autosomal — it would
#' silently overstate the variance contribution of any QTL on a hemizygous
#' (`copy_mode != "full"`) locus, since a single-copy Bernoulli(p) draw has
#' genic variance `p*(1-p)*a^2`, not `2*p*(1-p)*a^2`. Rather than generalizing
#' the formula (an approximation that would need to average over the sexes
#' actually present in the population), `define_additive_effects()` fails
#' loudly instead when `scale_to_target = TRUE` and any selected locus is
#' sex-linked/organelle.
#'
#' @param conn A DBI connection.
#' @param locus_names Character vector of locus names to check.
#' @return Invisible `NULL` on success; errors otherwise.
#' @keywords internal
assert_qtl_autosomal <- function(conn, locus_names) {
  if (length(locus_names) == 0L) return(invisible(NULL))

  bad <- DBI::dbGetQuery(conn, paste0(
    "SELECT gm.locus_name, gm.chr_name, cm.copy_mode_M, cm.copy_mode_F ",
    "FROM genome_meta gm JOIN chr_meta cm ON gm.chr_name = cm.chr_name ",
    "WHERE gm.locus_name IN (", sql_in_list(locus_names, what = "locus name"), ") ",
    "AND (cm.copy_mode_M != 'full' OR cm.copy_mode_F != 'full') ",
    "ORDER BY gm.locus_name"
  ))

  if (nrow(bad) > 0) {
    example <- bad[1, ]
    stop(
      "Falconer variance scaling (scale_to_target = TRUE) assumes diploid/",
      "autosomal QTL; locus '", example$locus_name, "' is on chromosome '",
      example$chr_name, "' with copy_mode_M = '", example$copy_mode_M,
      "', copy_mode_F = '", example$copy_mode_F, "' (", nrow(bad),
      " affected locus/loci total). Pass scale_to_target = FALSE and supply ",
      "`effects` manually for sex-linked/organelle QTL, or exclude these loci ",
      "from the QTL set.",
      call. = FALSE
    )
  }

  invisible(NULL)
}


#' Is a chromosome "plain" (both sexes full copy, recombines in both sexes)?
#'
#' A chromosome is handled by the fast, unchanged autosome path in
#' [add_founders()]/[add_offspring()] iff both sexes get a full complement and
#' recombination occurs in both sexes — i.e. `chr_meta` has its default row for
#' it. Any other combination (sex-linked copy mode, or achiasmy in either sex)
#' routes through the special-chromosome path.
#'
#' @param chr_row One row of `chr_meta` (as returned by [get_chr_meta_map()]).
#' @return Logical scalar.
#' @keywords internal
is_plain_autosome <- function(chr_row) {
  identical(chr_row$copy_mode_M, "full") &&
    identical(chr_row$copy_mode_F, "full") &&
    isTRUE(chr_row$recombines_M) &&
    isTRUE(chr_row$recombines_F)
}
