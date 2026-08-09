#' Resolve chromosome inheritance (copy counts) to one row per chromosome
#'
#' @description
#' Internal helper. Returns exactly one `chr_inheritance` row per
#' `genome_meta.chr_name`, applying the offspring-sex/line fallback precedence for
#' an offspring of sex `S` and line `L`. This is the single source of per-parent
#' copy counts for [add_founders()] (a founder is treated as an offspring of its
#' own line/sex) and [add_offspring()].
#'
#' Per-chromosome precedence (first match wins), mirroring `resolve_genome_map()`:
#' `(offspring_sex=S, line=L)` -> `(offspring_sex=S, NULL)` ->
#' `(NULL, line=L)` -> `(NULL, NULL)`.
#'
#' **Whose sex/line?** Call with the **offspring's** sex and line — the karyotype
#' an individual of that line/sex should carry.
#'
#' @param conn A DuckDB connection.
#' @param offspring_sex Offspring sex (`"M"`/`"F"`) or `NULL` to resolve only the
#'   sex-agnostic (`offspring_sex IS NULL`) rows.
#' @param line_name Offspring line or `NULL` to resolve only the line-agnostic
#'   (`line_name IS NULL`) rows.
#' @return A data frame `chr_name, from_parent_1, from_parent_2`, one row per
#'   chromosome, ordered by `chr_name`. Errors if any chromosome has no resolved
#'   row.
#' @keywords internal
#' @noRd
resolve_chr_inheritance <- function(conn, offspring_sex = NULL, line_name = NULL) {
  sx <- DBI::dbQuoteLiteral(conn, if (is.null(offspring_sex)) NA_character_ else offspring_sex)
  ln <- DBI::dbQuoteLiteral(conn, if (is.null(line_name))     NA_character_ else line_name)

  q <- sprintf("
    WITH cand AS (
      SELECT chr_name, from_parent_1, from_parent_2,
             (CASE WHEN offspring_sex IS NOT NULL THEN 2 ELSE 0 END) +
             (CASE WHEN line_name     IS NOT NULL THEN 1 ELSE 0 END) AS prio
      FROM chr_inheritance
      WHERE (offspring_sex IS NULL OR offspring_sex = %s)
        AND (line_name     IS NULL OR line_name     = %s)
    ),
    ranked AS (
      SELECT chr_name, from_parent_1, from_parent_2,
             ROW_NUMBER() OVER (PARTITION BY chr_name ORDER BY prio DESC) AS rn
      FROM cand
    )
    SELECT c.chr_name, r.from_parent_1, r.from_parent_2
    FROM (SELECT DISTINCT chr_name FROM genome_meta) c
    LEFT JOIN ranked r ON c.chr_name = r.chr_name AND r.rn = 1
    ORDER BY c.chr_name", sx, ln)

  res <- DBI::dbGetQuery(conn, q)

  unresolved <- is.na(res$from_parent_1) | is.na(res$from_parent_2)
  if (any(unresolved)) {
    miss <- res$chr_name[unresolved]
    stop(sprintf(
      "resolve_chr_inheritance(): %d chromosome(s) have no chr_inheritance row for (offspring_sex=%s, line_name=%s). First missing: '%s'.",
      length(miss),
      if (is.null(offspring_sex)) "NULL" else offspring_sex,
      if (is.null(line_name))     "NULL" else line_name,
      miss[1L]), call. = FALSE)
  }
  res$from_parent_1 <- as.integer(res$from_parent_1)
  res$from_parent_2 <- as.integer(res$from_parent_2)
  res
}


#' Resolve chromosome recombination to one row per chromosome
#'
#' @description
#' Internal helper. Returns exactly one `chr_recombination` row per
#' `genome_meta.chr_name`, applying the parent-sex/line fallback precedence for a
#' gamete produced by a parent of sex `S` and line `L`.
#'
#' **Whose sex/line?** Call with the **producing parent's** sex and line —
#' recombination is a property of the parent's germline, not the offspring's. In
#' an A x B cross the sire's gamete recombines under Line A's rule, the dam's
#' under Line B's. Passing the offspring line here would silently defeat pure-line
#' recombination defaults in crossbred offspring.
#'
#' @param conn A DuckDB connection.
#' @param parent_sex Producing-parent sex (`"M"`/`"F"`) or `NULL` for the
#'   sex-agnostic (`parent_sex IS NULL`) rows.
#' @param line_name Producing-parent line or `NULL` for the line-agnostic rows.
#' @return A data frame `chr_name, recombines`, one row per chromosome, ordered by
#'   `chr_name`. Errors if any chromosome has no resolved row.
#' @keywords internal
#' @noRd
resolve_chr_recombination <- function(conn, parent_sex = NULL, line_name = NULL) {
  sx <- DBI::dbQuoteLiteral(conn, if (is.null(parent_sex)) NA_character_ else parent_sex)
  ln <- DBI::dbQuoteLiteral(conn, if (is.null(line_name))  NA_character_ else line_name)

  q <- sprintf("
    WITH cand AS (
      SELECT chr_name, recombines,
             (CASE WHEN parent_sex IS NOT NULL THEN 2 ELSE 0 END) +
             (CASE WHEN line_name  IS NOT NULL THEN 1 ELSE 0 END) AS prio
      FROM chr_recombination
      WHERE (parent_sex IS NULL OR parent_sex = %s)
        AND (line_name  IS NULL OR line_name  = %s)
    ),
    ranked AS (
      SELECT chr_name, recombines,
             ROW_NUMBER() OVER (PARTITION BY chr_name ORDER BY prio DESC) AS rn
      FROM cand
    )
    SELECT c.chr_name, r.recombines
    FROM (SELECT DISTINCT chr_name FROM genome_meta) c
    LEFT JOIN ranked r ON c.chr_name = r.chr_name AND r.rn = 1
    ORDER BY c.chr_name", sx, ln)

  res <- DBI::dbGetQuery(conn, q)

  if (any(is.na(res$recombines))) {
    miss <- res$chr_name[is.na(res$recombines)]
    stop(sprintf(
      "resolve_chr_recombination(): %d chromosome(s) have no chr_recombination row for (parent_sex=%s, line_name=%s). First missing: '%s'.",
      length(miss),
      if (is.null(parent_sex)) "NULL" else parent_sex,
      if (is.null(line_name))  "NULL" else line_name,
      miss[1L]), call. = FALSE)
  }
  res$recombines <- as.logical(res$recombines)
  res
}


#' Within-table sex-vs-line shadowing check (deterministic)
#'
#' For each configured non-NULL line `L` and each supported sex `S in {M,F}`:
#' if a `(sex=S, line=NULL)` row and a `(sex=NULL, line=L)` row both exist and
#' their value columns **differ**, and no explicit `(sex=S, line=L)` row is
#' present, error — the user must materialize the explicit `(S, L)` row so the
#' intended value is unambiguous. Dormant until `line_name` is populated.
#'
#' @param rows Full table as a data frame with `chr_name`, `line_name`, the sex
#'   column, and the value columns.
#' @param sex_col Name of the sex column (`"offspring_sex"` / `"parent_sex"`).
#' @param value_cols Character vector of value columns to compare.
#' @param label Caller label for the error message.
#' @keywords internal
#' @noRd
.chr_shadow_check <- function(rows, sex_col, value_cols, label) {
  lines <- unique(rows$line_name[!is.na(rows$line_name)])
  if (!length(lines)) return(invisible(NULL))
  val_key <- function(r) paste(vapply(value_cols, function(cc) as.character(r[[cc]]), character(1)),
                               collapse = ",")
  for (cn in unique(rows$chr_name)) {
    crows <- rows[rows$chr_name == cn, , drop = FALSE]
    sx <- crows[[sex_col]]
    for (L in lines) {
      for (S in c("M", "F")) {
        r_s_null <- crows[!is.na(sx) & sx == S & is.na(crows$line_name), , drop = FALSE]
        r_null_l <- crows[is.na(sx) & !is.na(crows$line_name) & crows$line_name == L, , drop = FALSE]
        r_s_l    <- crows[!is.na(sx) & sx == S &
                            !is.na(crows$line_name) & crows$line_name == L, , drop = FALSE]
        if (nrow(r_s_null) == 1L && nrow(r_null_l) == 1L && nrow(r_s_l) == 0L &&
            val_key(r_s_null) != val_key(r_null_l)) {
          stop(sprintf(
            "%s: shadowing conflict on chromosome '%s' for sex '%s', line '%s'. The (sex=%s, line=NULL) rule and the (sex=NULL, line=%s) rule differ and no explicit (sex=%s, line=%s) row exists. Add an explicit (sex=%s, line=%s) row so the intended value is unambiguous.",
            label, cn, S, L, S, L, S, L, S, L), call. = FALSE)
        }
      }
    }
  }
  invisible(NULL)
}


#' Validate the structural integrity of the chr_inheritance table
#'
#' @description
#' Internal helper. Run after every write to `chr_inheritance` (inside the write
#' transaction). Checks: NULL-normalized logical-key uniqueness
#' `(chr_name, offspring_sex, line_name)`; no orphan `chr_name` (R-enforced FK to
#' `genome_meta`); valid `offspring_sex`; non-negative counts summing to `<= 2`
#' (the constant diploid-release invariant, not an `ind_meta.ploidy` lookup);
#' resolvability for both `M` and `F`; and the deterministic sex-vs-line shadowing
#' check.
#'
#' @param conn A DuckDB connection.
#' @return `invisible(TRUE)` on success; errors otherwise.
#' @keywords internal
#' @noRd
validate_chr_inheritance <- function(conn) {
  # 1. Duplicate logical keys (NULL-normalized).
  dups <- DBI::dbGetQuery(conn, "
    SELECT COUNT(*) AS n FROM (
      SELECT chr_name, COALESCE(offspring_sex, '<NA>') AS s,
             COALESCE(line_name, '<NA>') AS l
      FROM chr_inheritance
      GROUP BY chr_name, COALESCE(offspring_sex, '<NA>'), COALESCE(line_name, '<NA>')
      HAVING COUNT(*) > 1
    )")$n
  if (dups > 0) {
    stop(sprintf("validate_chr_inheritance(): %d duplicate logical key(s) (chr_name, offspring_sex, line_name).",
                 dups), call. = FALSE)
  }

  # 2. Orphan chr_name (R-enforced FK to genome_meta).
  orphans <- DBI::dbGetQuery(conn, "
    SELECT COUNT(DISTINCT ci.chr_name) AS n
    FROM chr_inheritance ci
    LEFT JOIN (SELECT DISTINCT chr_name FROM genome_meta) g ON ci.chr_name = g.chr_name
    WHERE g.chr_name IS NULL")$n
  if (orphans > 0) {
    stop(sprintf("validate_chr_inheritance(): %d chr_name(s) not present in genome_meta.",
                 orphans), call. = FALSE)
  }

  # 3. Invalid offspring_sex.
  bad_sex <- DBI::dbGetQuery(conn, "
    SELECT COUNT(*) AS n FROM chr_inheritance
    WHERE offspring_sex IS NOT NULL AND offspring_sex NOT IN ('M', 'F')")$n
  if (bad_sex > 0) {
    stop(sprintf("validate_chr_inheritance(): %d row(s) have an invalid offspring_sex (not NULL/'M'/'F').",
                 bad_sex), call. = FALSE)
  }

  # 4. Counts: non-negative, sum <= 2 (constant diploid-release invariant).
  bad_cnt <- DBI::dbGetQuery(conn, "
    SELECT COUNT(*) AS n FROM chr_inheritance
    WHERE from_parent_1 IS NULL OR from_parent_2 IS NULL
       OR from_parent_1 < 0 OR from_parent_2 < 0
       OR from_parent_1 + from_parent_2 > 2")$n
  if (bad_cnt > 0) {
    stop(sprintf("validate_chr_inheritance(): %d row(s) violate 0 <= from_parent_1 + from_parent_2 <= 2 (diploid release).",
                 bad_cnt), call. = FALSE)
  }

  # 5. Resolves for both offspring sexes (hard error via the resolver).
  resolve_chr_inheritance(conn, "M")
  resolve_chr_inheritance(conn, "F")

  # 6. Deterministic sex-vs-line shadowing check.
  rows <- DBI::dbGetQuery(conn,
    "SELECT chr_name, offspring_sex, line_name, from_parent_1, from_parent_2 FROM chr_inheritance")
  .chr_shadow_check(rows, "offspring_sex", c("from_parent_1", "from_parent_2"),
                    "validate_chr_inheritance()")

  invisible(TRUE)
}


#' Validate the structural integrity of the chr_recombination table
#'
#' @description
#' Internal helper. Run after every write to `chr_recombination` (inside the write
#' transaction). Same shape as `validate_chr_inheritance()`: NULL-normalized key
#' uniqueness `(chr_name, parent_sex, line_name)`; no orphan `chr_name`; valid
#' `parent_sex`; non-null `recombines`; resolvability for both `M` and `F`; and
#' the deterministic sex-vs-line shadowing check on `recombines`.
#'
#' @param conn A DuckDB connection.
#' @return `invisible(TRUE)` on success; errors otherwise.
#' @keywords internal
#' @noRd
validate_chr_recombination <- function(conn) {
  dups <- DBI::dbGetQuery(conn, "
    SELECT COUNT(*) AS n FROM (
      SELECT chr_name, COALESCE(parent_sex, '<NA>') AS s,
             COALESCE(line_name, '<NA>') AS l
      FROM chr_recombination
      GROUP BY chr_name, COALESCE(parent_sex, '<NA>'), COALESCE(line_name, '<NA>')
      HAVING COUNT(*) > 1
    )")$n
  if (dups > 0) {
    stop(sprintf("validate_chr_recombination(): %d duplicate logical key(s) (chr_name, parent_sex, line_name).",
                 dups), call. = FALSE)
  }

  orphans <- DBI::dbGetQuery(conn, "
    SELECT COUNT(DISTINCT cr.chr_name) AS n
    FROM chr_recombination cr
    LEFT JOIN (SELECT DISTINCT chr_name FROM genome_meta) g ON cr.chr_name = g.chr_name
    WHERE g.chr_name IS NULL")$n
  if (orphans > 0) {
    stop(sprintf("validate_chr_recombination(): %d chr_name(s) not present in genome_meta.",
                 orphans), call. = FALSE)
  }

  bad_sex <- DBI::dbGetQuery(conn, "
    SELECT COUNT(*) AS n FROM chr_recombination
    WHERE parent_sex IS NOT NULL AND parent_sex NOT IN ('M', 'F')")$n
  if (bad_sex > 0) {
    stop(sprintf("validate_chr_recombination(): %d row(s) have an invalid parent_sex (not NULL/'M'/'F').",
                 bad_sex), call. = FALSE)
  }

  bad_rec <- DBI::dbGetQuery(conn, "
    SELECT COUNT(*) AS n FROM chr_recombination WHERE recombines IS NULL")$n
  if (bad_rec > 0) {
    stop(sprintf("validate_chr_recombination(): %d row(s) have a NULL recombines.",
                 bad_rec), call. = FALSE)
  }

  resolve_chr_recombination(conn, "M")
  resolve_chr_recombination(conn, "F")

  rows <- DBI::dbGetQuery(conn,
    "SELECT chr_name, parent_sex, line_name, recombines FROM chr_recombination")
  .chr_shadow_check(rows, "parent_sex", "recombines", "validate_chr_recombination()")

  invisible(TRUE)
}


#' Build the resolved chromosome-rule lookup for one line context
#'
#' @description
#' Convenience reader used by [add_founders()] and [add_offspring()] to avoid
#' re-querying the chromosome tables per chromosome inside a loop. Resolves both
#' tables for `"M"` and `"F"` at the given `line_name` and returns a named list
#' keyed by `chr_name`.
#'
#' Each element is a list with:
#' * `from_parent_1` — named numeric `c(M = , F = )`, copies from parent_1 (sire)
#' * `from_parent_2` — named numeric `c(M = , F = )`, copies from parent_2 (dam)
#' * `recombines`    — named logical `c(M = , F = )`, does the parent-of-that-sex
#'   recombine
#'
#' For inheritance the `M`/`F` index is the **offspring** sex; for recombination
#' it is the **producing-parent** sex. Callers pass the correct entity's
#' `line_name` (offspring line for inheritance, parent line for recombination); in
#' the common no-line-specific case `line_name = NULL` resolves the defaults for
#' every individual.
#'
#' @param conn A DBI connection.
#' @param line_name Line context, or `NULL` for the line-agnostic default.
#' @return Named list (by `chr_name`) of resolved rules.
#' @keywords internal
#' @noRd
get_chr_rules_map <- function(conn, line_name = NULL) {
  inh_M <- resolve_chr_inheritance(conn, "M", line_name)
  inh_F <- resolve_chr_inheritance(conn, "F", line_name)
  rec_M <- resolve_chr_recombination(conn, "M", line_name)
  rec_F <- resolve_chr_recombination(conn, "F", line_name)

  chrs <- inh_M$chr_name   # all four are chr_name-ordered over the same set
  stats::setNames(lapply(seq_along(chrs), function(i) {
    list(
      from_parent_1 = c(M = inh_M$from_parent_1[i], F = inh_F$from_parent_1[i]),
      from_parent_2 = c(M = inh_M$from_parent_2[i], F = inh_F$from_parent_2[i]),
      recombines    = c(M = rec_M$recombines[i],    F = rec_F$recombines[i])
    )
  }), chrs)
}


#' Error if any of the given loci sit on a non-diploid-autosomal chromosome
#'
#' `rescale_effects_to_target()`'s Falconer variance formula
#' (`2 * p * (1-p) * a^2`) assumes every QTL is diploid/autosomal — it would
#' silently overstate the variance contribution of any QTL on a hemizygous locus
#' (a single-copy Bernoulli(p) draw has genic variance `p*(1-p)*a^2`, not
#' `2*p*(1-p)*a^2`). Rather than generalizing the formula,
#' `define_additive_effects()` fails loudly when `scale_to_target = TRUE` and any
#' selected locus is sex-linked/organelle.
#'
#' A chromosome is diploid-autosomal iff its **resolved** inheritance is
#' `from_parent_1 = 1, from_parent_2 = 1` for **both** offspring sexes (line-
#' agnostic default resolution — QTL variance scaling is population-level).
#'
#' @param conn A DBI connection.
#' @param locus_names Character vector of locus names to check.
#' @return Invisible `NULL` on success; errors otherwise.
#' @keywords internal
#' @noRd
assert_qtl_autosomal <- function(conn, locus_names) {
  if (length(locus_names) == 0L) return(invisible(NULL))

  inh_M <- resolve_chr_inheritance(conn, "M")
  inh_F <- resolve_chr_inheritance(conn, "F")   # same chr_name order as inh_M
  full_diploid <- inh_M$from_parent_1 == 1L & inh_M$from_parent_2 == 1L &
    inh_F$from_parent_1 == 1L & inh_F$from_parent_2 == 1L
  bad_chr <- inh_M$chr_name[!full_diploid]
  if (length(bad_chr) == 0L) return(invisible(NULL))

  bad <- DBI::dbGetQuery(conn, paste0(
    "SELECT locus_name, chr_name FROM genome_meta ",
    "WHERE locus_name IN (", sql_in_list(locus_names, what = "locus name"), ") ",
    "AND chr_name IN (", sql_in_list(bad_chr, what = "chromosome name"), ") ",
    "ORDER BY locus_name"
  ))

  if (nrow(bad) > 0) {
    example <- bad[1, ]
    stop(
      "Falconer variance scaling (scale_to_target = TRUE) assumes diploid/",
      "autosomal QTL; locus '", example$locus_name, "' is on chromosome '",
      example$chr_name, "', which is not diploid-autosomal (chr_inheritance is ",
      "not from_parent_1 = 1, from_parent_2 = 1 for both sexes; ", nrow(bad),
      " affected locus/loci total). Pass scale_to_target = FALSE and supply ",
      "`effects` manually for sex-linked/organelle QTL, or exclude these loci ",
      "from the QTL set.",
      call. = FALSE
    )
  }

  invisible(NULL)
}


#' Is a chromosome "plain" (both sexes full diploid copy, recombines in both)?
#'
#' A chromosome is handled by the fast, unchanged autosome path in
#' [add_founders()]/[add_offspring()] iff every offspring sex inherits one copy
#' from each parent (`from_parent_1 = from_parent_2 = 1`) and every producing-
#' parent sex recombines — i.e. the seeded default. Any other combination
#' (sex-linked copy count, or achiasmy in either sex) routes through the
#' special-chromosome path.
#'
#' @param chr_row One element of `get_chr_rules_map()` (resolved rule for one chr).
#' @return Logical scalar.
#' @keywords internal
#' @noRd
is_plain_autosome <- function(chr_row) {
  chr_row$from_parent_1[["M"]] == 1L && chr_row$from_parent_2[["M"]] == 1L &&
    chr_row$from_parent_1[["F"]] == 1L && chr_row$from_parent_2[["F"]] == 1L &&
    isTRUE(chr_row$recombines[["M"]]) && isTRUE(chr_row$recombines[["F"]])
}
