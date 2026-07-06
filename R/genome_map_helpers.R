#' Resolve the genetic map to one position per locus
#'
#' @description
#' Internal helper. Returns exactly one genetic-map row per
#' `genome_meta.locus_id`, applying the sex/line fallback precedence for a gamete
#' produced by a parent of sex `S` and line `L`, within a chosen `map_name`.
#' This is the single source of genetic positions for founder LD, recombination,
#' and (later) the optimized offspring kernel.
#'
#' Per-locus precedence (first match wins), mirroring `genome_effects`'s
#' line-fallback:
#' `(sex=S, line=L)` -> `(sex=S, line NULL)` -> `(sex NULL, line=L)` ->
#' `(sex NULL, line NULL)`.
#'
#' @param conn A DuckDB connection.
#' @param sex Parent sex (`"M"`/`"F"`) or `NULL` to resolve only the
#'   sex-agnostic (`sex IS NULL`) rows.
#' @param line_name Parent line or `NULL` to resolve only the line-agnostic
#'   (`line_name IS NULL`) rows.
#' @param map_name Map identity; default `"default"`.
#'
#' @return A data frame with one row per locus —
#'   `locus_id, locus_name, chr, chr_name, pos_bp, pos_cM` — ordered by
#'   `locus_id`. Errors if any locus has no resolved map row, or if `pos_cM` is
#'   not monotonic non-decreasing within any chromosome after fallback.
#'
#' @keywords internal
#' @noRd
resolve_genome_map <- function(conn, sex = NULL, line_name = NULL,
                               map_name = "default") {
  sx <- DBI::dbQuoteLiteral(conn, if (is.null(sex))       NA_character_ else sex)
  ln <- DBI::dbQuoteLiteral(conn, if (is.null(line_name)) NA_character_ else line_name)
  mp <- DBI::dbQuoteLiteral(conn, map_name)

  # Priority: sex-specific (+2) beats sex-agnostic; line-specific (+1) beats
  # line-agnostic. A NULL S/L only matches the corresponding IS NULL rows,
  # because `col = NULL` is never true in SQL. Highest priority wins per locus.
  q <- sprintf("
    WITH cand AS (
      SELECT locus_id, pos_cM,
             (CASE WHEN sex IS NOT NULL THEN 2 ELSE 0 END) +
             (CASE WHEN line_name IS NOT NULL THEN 1 ELSE 0 END) AS prio
      FROM genome_map
      WHERE map_name = %s
        AND (sex IS NULL OR sex = %s)
        AND (line_name IS NULL OR line_name = %s)
    ),
    ranked AS (
      SELECT locus_id, pos_cM,
             ROW_NUMBER() OVER (PARTITION BY locus_id ORDER BY prio DESC) AS rn
      FROM cand
    )
    SELECT m.locus_id, m.locus_name, m.chr, m.chr_name, m.pos_bp, r.pos_cM
    FROM genome_meta m
    LEFT JOIN ranked r ON m.locus_id = r.locus_id AND r.rn = 1
    ORDER BY m.locus_id", mp, sx, ln)

  res <- DBI::dbGetQuery(conn, q)

  # Every locus must resolve to exactly one map position.
  if (any(is.na(res$pos_cM))) {
    miss <- res$locus_id[is.na(res$pos_cM)]
    stop(sprintf(
      "resolve_genome_map(): %d locus/loci have no genome_map row for (sex=%s, line_name=%s, map_name=%s). First missing locus_id: %s.",
      length(miss),
      if (is.null(sex)) "NULL" else sex,
      if (is.null(line_name)) "NULL" else line_name,
      map_name, miss[1L]),
      call. = FALSE)
  }

  # pos_cM must be monotonic non-decreasing within each chromosome, on the
  # RESOLVED map (not just the raw stored rows). Rows are locus_id-ordered, which
  # is physical/pos order by construction.
  non_mono <- tapply(res$pos_cM, res$chr, function(v) any(diff(v) < 0))
  if (any(non_mono)) {
    bad <- names(non_mono)[non_mono]
    stop(sprintf(
      "resolve_genome_map(): pos_cM is not monotonic non-decreasing within chromosome(s): %s.",
      paste(bad, collapse = ", ")),
      call. = FALSE)
  }

  res
}

#' Validate the structural integrity of the genome_map table
#'
#' @description
#' Internal helper. Run after every write to `genome_map`. Checks the logical
#' key, agreement with `genome_meta`, and column validity. This is where most
#' downstream map bugs are caught early.
#'
#' Checks:
#' * logical-key uniqueness `(locus_id, sex, line_name, map_name)` with explicit
#'   NULL normalization (SQL `NULL != NULL` would otherwise miss duplicate
#'   `(locus, NULL, NULL, 'default')` rows);
#' * every `locus_id`/`locus_name` agrees with `genome_meta` — no orphan loci,
#'   and one `locus_id` is never paired with two different `locus_name`s;
#' * every row has a valid `sex` (`NULL`/`'M'`/`'F'`) and a non-missing
#'   `map_name`.
#'
#' @param conn A DuckDB connection.
#' @return `invisible(TRUE)` on success; errors otherwise.
#'
#' @keywords internal
#' @noRd
validate_genome_map <- function(conn) {
  # 1. Duplicate logical keys (NULL-normalized).
  dups <- DBI::dbGetQuery(conn, "
    SELECT COUNT(*) AS n FROM (
      SELECT locus_id, COALESCE(sex, '<NA>') AS s,
             COALESCE(line_name, '<NA>') AS l, map_name
      FROM genome_map
      GROUP BY locus_id, COALESCE(sex, '<NA>'), COALESCE(line_name, '<NA>'), map_name
      HAVING COUNT(*) > 1
    )")$n
  if (dups > 0) {
    stop(sprintf("validate_genome_map(): %d duplicate logical key(s) (locus_id, sex, line_name, map_name) in genome_map.",
                 dups), call. = FALSE)
  }

  # 2a. Orphan loci: genome_map rows whose locus_id is absent from genome_meta.
  orphans <- DBI::dbGetQuery(conn, "
    SELECT COUNT(DISTINCT gm.locus_id) AS n
    FROM genome_map gm
    LEFT JOIN genome_meta m ON gm.locus_id = m.locus_id
    WHERE m.locus_id IS NULL")$n
  if (orphans > 0) {
    stop(sprintf("validate_genome_map(): %d genome_map locus_id(s) not present in genome_meta.",
                 orphans), call. = FALSE)
  }

  # 2b. locus_id <-> locus_name mismatch against genome_meta.
  mism <- DBI::dbGetQuery(conn, "
    SELECT COUNT(*) AS n
    FROM genome_map gm
    JOIN genome_meta m ON gm.locus_id = m.locus_id
    WHERE gm.locus_name IS DISTINCT FROM m.locus_name")$n
  if (mism > 0) {
    stop(sprintf("validate_genome_map(): %d genome_map row(s) have a locus_name that disagrees with genome_meta.",
                 mism), call. = FALSE)
  }

  # 3. Invalid sex or missing map_name.
  bad <- DBI::dbGetQuery(conn, "
    SELECT COUNT(*) AS n FROM genome_map
    WHERE (sex IS NOT NULL AND sex NOT IN ('M', 'F')) OR map_name IS NULL")$n
  if (bad > 0) {
    stop(sprintf("validate_genome_map(): %d genome_map row(s) have an invalid sex (not NULL/'M'/'F') or a NULL map_name.",
                 bad), call. = FALSE)
  }

  # 4. pos_cM must be present and finite for every stored row (genome_map rows are,
  #    by definition, loci with a genetic position).
  bad_pos <- DBI::dbGetQuery(conn, "
    SELECT COUNT(*) AS n FROM genome_map
    WHERE pos_cM IS NULL OR NOT isfinite(pos_cM)")$n
  if (bad_pos > 0) {
    stop(sprintf("validate_genome_map(): %d genome_map row(s) have a NULL or non-finite pos_cM.",
                 bad_pos), call. = FALSE)
  }

  # 5. Surrogate key id_genome_map must be unique (R-managed via next_int_id;
  #    no DB constraint, so enforce here).
  dup_id <- DBI::dbGetQuery(conn, "
    SELECT COUNT(*) AS n FROM (
      SELECT id_genome_map FROM genome_map
      GROUP BY id_genome_map HAVING COUNT(*) > 1
    )")$n
  if (dup_id > 0) {
    stop(sprintf("validate_genome_map(): %d duplicate id_genome_map value(s).",
                 dup_id), call. = FALSE)
  }

  invisible(TRUE)
}
