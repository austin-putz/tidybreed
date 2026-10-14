#' Define a chromosome's inheritance or recombination rule
#'
#' @description
#' Changes how one chromosome is passed on, for sex chromosomes (X/Y, Z/W,
#' X0, Z0) and organelles (mitochondria, plastids). [define_genome()] starts
#' every chromosome as a diploid autosome: one copy from each parent,
#' recombining in both sexes. Call `define_chromosome()` after
#' [define_genome()] and **before [add_founders()]**, once per rule that
#' differs from that default.
#'
#' Each call sets **one** of two things:
#'
#' * **Inheritance**: give `from_parent_1` and `from_parent_2`, the number of
#'   copies an offspring inherits from its sire (parent 1) and its dam
#'   (parent 2), optionally for one `offspring_sex`. Writes
#'   `chr_inheritance`.
#' * **Recombination**: give `recombines`, whether a parent recombines this
#'   chromosome when it makes gametes, optionally for one `parent_sex`.
#'   Writes `chr_recombination`.
#'
#' They are separate calls because "which sex" means different animals: the
#' **offspring** for inheritance, the **parent** for recombination.
#'
#' @details
#' **Copy counts** are absolute numbers of copies, at most 2 in total:
#'
#' | Chromosome | `offspring_sex` | `from_parent_1` (sire) | `from_parent_2` (dam) |
#' |---|---|---|---|
#' | Autosome (the default) | either | 1 | 1 |
#' | Mammalian X in males | `"M"` | 0 | 1 |
#' | Mammalian Y in males | `"M"` | 1 | 0 |
#' | Mammalian Y in females | `"F"` | 0 | 0 |
#' | Avian Z in females | `"F"` | 1 | 0 |
#' | Avian W in females | `"F"` | 0 | 1 |
#' | Avian W in males | `"M"` | 0 | 0 |
#' | Maternal mitochondria | either (`NULL`) | 0 | 1 |
#'
#' Only the sexes that differ from the default need a rule: a mammalian X in
#' females keeps the default `1, 1`. A chromosome a parent carries in one copy
#' is passed on unchanged; one it carries in two copies recombines as
#' `chr_recombination` says.
#'
#' **Which rule applies.** For an animal of sex S in line L, the most
#' specific stored rule wins: `(S, L)`, then `(S, any line)`, then
#' `(any sex, L)`, then the default `(any sex, any line)`. Inheritance rules
#' are looked up with the offspring's sex and line, recombination rules with
#' the parent's.
#'
#' **Timing.** Founders and offspring are written with the rules in force at
#' the time, so set every inheritance rule before [add_founders()]. Once
#' individuals exist, an inheritance call that would change the copies any of
#' them carry is refused (nothing is written): their stored haplotypes would
#' no longer match, and breeding from them would fail. A call that changes
#' nothing for existing individuals, such as a rule for a line with no
#' individuals yet, is still allowed. Recombination rules can be changed at
#' any time; they apply to gametes made afterwards.
#'
#' Every write is checked in the same transaction (one rule per key, both
#' sexes resolvable, at most 2 copies); a rejected call leaves the previous
#' rule in place. See the rules with `get_table(pop, "chr_inheritance")` and
#' `get_table(pop, "chr_recombination")`.
#'
#' @param pop A `tidybreed_pop` with a genome from [define_genome()].
#' @param chr_name Character scalar: a chromosome name from
#'   [define_genome()]'s `chr_names` (e.g. `"X"`).
#' @param offspring_sex Inheritance calls only: `"M"` or `"F"` for a rule that
#'   applies to offspring of that sex, or `NULL` (default) for both.
#' @param parent_sex Recombination calls only: `"M"` or `"F"` for a rule that
#'   applies when a parent of that sex makes gametes, or `NULL` (default) for
#'   both.
#' @param line_name Optional line name: the rule applies only to animals of
#'   that line (the offspring's line for inheritance, the parent's for
#'   recombination), for crossbreeding designs whose lines differ. `NULL`
#'   (default) applies to every line. A chromosome that any line's rule makes
#'   special is drawn on [add_offspring()]'s slower per-offspring path for
#'   every line.
#' @param from_parent_1,from_parent_2 Inheritance calls: whole numbers,
#'   always given together, of copies from the sire (parent 1) and the dam
#'   (parent 2). Their sum is 0, 1 or 2; see the table in *Details*.
#' @param recombines Recombination calls: `TRUE` if the chromosome recombines
#'   when a parent makes gametes, `FALSE` if it is passed on whole (a Y, a W).
#'   To stop one sex recombining on every chromosome, use `recombines_M` /
#'   `recombines_F` in [define_genome()] instead.
#' @param overwrite Logical. `TRUE` (default) replaces an existing rule with
#'   the same chromosome, sex and line. `FALSE` refuses instead, to guard
#'   against redefining a rule by accident. Every chromosome already has a
#'   default rule (both sexes, all lines) from [define_genome()], so
#'   `overwrite = FALSE` is useful only for sex- or line-specific rules.
#'
#' @return The `tidybreed_pop`, invisibly. Assign it back.
#'
#' @seealso [define_genome()], [add_founders()], [add_offspring()].
#'
#' @examples
#' # Mammalian X and Y: inheritance for the sexes that differ, and a
#' # non-recombining Y
#' pop <- open_pop(pop_name = "mammal", db_name = ":memory:") |>
#'   define_genome(n_loci = 300, n_chr = 3, chr_len_Mb = 100,
#'                 chr_names = c("1", "X", "Y"))
#' pop <- pop |>
#'   define_chromosome("X", offspring_sex = "M", from_parent_1 = 0, from_parent_2 = 1) |>
#'   define_chromosome("Y", offspring_sex = "M", from_parent_1 = 1, from_parent_2 = 0) |>
#'   define_chromosome("Y", offspring_sex = "F", from_parent_1 = 0, from_parent_2 = 0) |>
#'   define_chromosome("Y", recombines = FALSE)
#' get_table(pop, "chr_inheritance")
#'
#' # Founders are written with those rules: a male has one X copy and one Y
#' set.seed(1)
#' pop <- pop |> define_founder_haplotypes(n_haplotypes = 40)
#' pop <- pop |>
#'   get_table("founder_haplotypes") |>
#'   add_founders(n_males = 2, n_females = 2, line_name = "A")
#' close_pop(pop)
#'
#' # Avian Z and W: females (ZW) are the heterogametic sex
#' pop <- open_pop(pop_name = "chicken", db_name = ":memory:") |>
#'   define_genome(n_loci = 300, n_chr = 3, chr_len_Mb = 100,
#'                 chr_names = c("1", "Z", "W"))
#' pop <- pop |>
#'   define_chromosome("Z", offspring_sex = "F", from_parent_1 = 1, from_parent_2 = 0) |>
#'   define_chromosome("W", offspring_sex = "F", from_parent_1 = 0, from_parent_2 = 1) |>
#'   define_chromosome("W", offspring_sex = "M", from_parent_1 = 0, from_parent_2 = 0) |>
#'   define_chromosome("W", recombines = FALSE)
#' close_pop(pop)
#'
#' # Mitochondria: one maternal copy in both sexes, never recombining
#' pop <- open_pop(pop_name = "mito", db_name = ":memory:") |>
#'   define_genome(n_loci = 200, n_chr = 2, chr_len_Mb = c(100, 0.02),
#'                 chr_names = c("1", "MT"))
#' pop <- pop |>
#'   define_chromosome("MT", from_parent_1 = 0, from_parent_2 = 1) |>
#'   define_chromosome("MT", recombines = FALSE)
#' close_pop(pop)
#' @export
define_chromosome <- function(pop, chr_name,
                              offspring_sex = NULL,
                              parent_sex    = NULL,
                              line_name     = NULL,
                              from_parent_1 = NULL,
                              from_parent_2 = NULL,
                              recombines    = NULL,
                              overwrite     = TRUE) {

  validate_tidybreed_pop(pop)
  conn <- pop$db_conn

  stopifnot(is.character(chr_name), length(chr_name) == 1L, nchar(chr_name) > 0L)
  stopifnot(is.logical(overwrite), length(overwrite) == 1L, !is.na(overwrite))

  known_chrs <- DBI::dbGetQuery(conn,
    "SELECT DISTINCT chr_name FROM genome_meta")$chr_name
  if (!chr_name %in% known_chrs) {
    stop("chr_name '", chr_name, "' not found in genome_meta. ",
         "Define it via define_genome()'s chr_names argument first.",
         call. = FALSE)
  }

  chk_sex <- function(x, nm) {
    if (!is.null(x) && (!is.character(x) || length(x) != 1L || !x %in% c("M", "F"))) {
      stop(nm, " must be NULL, 'M', or 'F'.", call. = FALSE)
    }
  }
  if (!is.null(line_name) &&
      (!is.character(line_name) || length(line_name) != 1L || nchar(line_name) == 0L)) {
    stop("line_name must be NULL or a non-empty string.", call. = FALSE)
  }

  has_inheritance   <- !is.null(from_parent_1) || !is.null(from_parent_2)
  has_recombination <- !is.null(recombines)

  if (has_inheritance && has_recombination) {
    stop("define_chromosome() sets exactly one concern per call: supply either ",
         "from_parent_1/from_parent_2 (inheritance) OR recombines (recombination), ",
         "not both. Use two calls.", call. = FALSE)
  }
  if (!has_inheritance && !has_recombination) {
    stop("define_chromosome() needs a concern: supply from_parent_1 + ",
         "from_parent_2 (inheritance) or recombines (recombination).",
         call. = FALSE)
  }

  if (has_inheritance) {
    if (!is.null(parent_sex)) {
      stop("parent_sex is a recombination key; do not supply it on an ",
           "inheritance call (from_parent_* given). Set recombination in a ",
           "separate define_chromosome(..., recombines = ) call.", call. = FALSE)
    }
    if (is.null(from_parent_1) || is.null(from_parent_2)) {
      stop("from_parent_1 and from_parent_2 must be supplied together.",
           call. = FALSE)
    }
    ok_count <- function(v) is.numeric(v) && length(v) == 1L && !is.na(v) &&
      is.finite(v) && v >= 0 && v == as.integer(v)
    if (!ok_count(from_parent_1) || !ok_count(from_parent_2)) {
      stop("from_parent_1 and from_parent_2 must each be a single non-negative ",
           "integer.", call. = FALSE)
    }
    fp1 <- as.integer(from_parent_1)
    fp2 <- as.integer(from_parent_2)
    if (fp1 + fp2 > 2L) {
      stop("from_parent_1 + from_parent_2 must be <= 2 in this diploid release ",
           "(got ", fp1, " + ", fp2, ").", call. = FALSE)
    }
    chk_sex(offspring_sex, "offspring_sex")

    table       <- "chr_inheritance"
    sex_col     <- "offspring_sex"
    sex_val     <- offspring_sex
    insert_cols <- "chr_name, offspring_sex, line_name, from_parent_1, from_parent_2"
    df <- data.frame(
      chr_name      = chr_name,
      offspring_sex = if (is.null(offspring_sex)) NA_character_ else offspring_sex,
      line_name     = if (is.null(line_name))     NA_character_ else line_name,
      from_parent_1 = fp1,
      from_parent_2 = fp2,
      stringsAsFactors = FALSE
    )
    validator <- validate_chr_inheritance
    msg <- sprintf(
      "Defined chr_inheritance rule for '%s' (offspring_sex=%s, line_name=%s): from_parent_1=%d, from_parent_2=%d",
      chr_name, if (is.null(offspring_sex)) "NULL" else offspring_sex,
      if (is.null(line_name)) "NULL" else line_name, fp1, fp2)

  } else {
    if (!is.null(offspring_sex)) {
      stop("offspring_sex is an inheritance key; do not supply it on a ",
           "recombination call (recombines given). Set inheritance in a ",
           "separate define_chromosome(..., from_parent_1 = , from_parent_2 = ) ",
           "call.", call. = FALSE)
    }
    if (!is.logical(recombines) || length(recombines) != 1L || is.na(recombines)) {
      stop("recombines must be a single non-NA logical.", call. = FALSE)
    }
    chk_sex(parent_sex, "parent_sex")

    table       <- "chr_recombination"
    sex_col     <- "parent_sex"
    sex_val     <- parent_sex
    insert_cols <- "chr_name, parent_sex, line_name, recombines"
    df <- data.frame(
      chr_name   = chr_name,
      parent_sex = if (is.null(parent_sex)) NA_character_ else parent_sex,
      line_name  = if (is.null(line_name))  NA_character_ else line_name,
      recombines = recombines,
      stringsAsFactors = FALSE
    )
    validator <- validate_chr_recombination
    msg <- sprintf(
      "Defined chr_recombination rule for '%s' (parent_sex=%s, line_name=%s): recombines=%s",
      chr_name, if (is.null(parent_sex)) "NULL" else parent_sex,
      if (is.null(line_name)) "NULL" else line_name, recombines)
  }

  # NULL-safe logical-key predicate (IS NOT DISTINCT FROM matches NULL to NULL).
  chr_lit  <- DBI::dbQuoteLiteral(conn, chr_name)
  sex_lit  <- DBI::dbQuoteLiteral(conn, if (is.null(sex_val))   NA_character_ else sex_val)
  line_lit <- DBI::dbQuoteLiteral(conn, if (is.null(line_name)) NA_character_ else line_name)
  key_pred <- sprintf(
    "chr_name IS NOT DISTINCT FROM %s AND %s IS NOT DISTINCT FROM %s AND line_name IS NOT DISTINCT FROM %s",
    chr_lit, sex_col, sex_lit, line_lit)

  tmp_tbl <- "__define_chromosome_tmp"

  # Whole write is one transaction: delete-then-insert (RNG-neutral register +
  # INSERT) -> run the matching validator -> COMMIT. Any failure rolls back, so a
  # rejected replacement never leaves the config table corrupt.
  DBI::dbExecute(conn, "BEGIN TRANSACTION")
  committed <- FALSE
  on.exit({
    try(duckdb::duckdb_unregister(conn, tmp_tbl), silent = TRUE)
    if (!committed) try(DBI::dbExecute(conn, "ROLLBACK"), silent = TRUE)
  }, add = TRUE)

  exists_n <- DBI::dbGetQuery(conn,
    sprintf("SELECT COUNT(*) AS n FROM %s WHERE %s", table, key_pred))$n

  # Animals already written carry the copies the rules gave them, so a rule
  # that changes them would break every later mating. Compared after the
  # write, before COMMIT: a harmless call (an identical rule, a rule for a
  # line with no animals yet) still goes through.
  if (has_inheritance) karyo_before <- .dc_karyotypes(conn, chr_name)

  if (!overwrite && exists_n > 0L) {
    stop(sprintf(
      "define_chromosome(): a %s row already exists for (chr_name=%s, %s=%s, line_name=%s) and overwrite = FALSE.",
      table, chr_name, sex_col,
      if (is.null(sex_val)) "NULL" else sex_val,
      if (is.null(line_name)) "NULL" else line_name), call. = FALSE)
  }

  if (exists_n > 0L) {
    DBI::dbExecute(conn, sprintf("DELETE FROM %s WHERE %s", table, key_pred))
  }
  duckdb::duckdb_register(conn, tmp_tbl, df)
  DBI::dbExecute(conn, sprintf(
    "INSERT INTO %s (%s) SELECT %s FROM %s", table, insert_cols, insert_cols, tmp_tbl))

  validator(conn)
  if (has_inheritance) {
    .dc_refuse_changed_karyotypes(karyo_before, .dc_karyotypes(conn, chr_name),
                                  chr_name)
  }

  DBI::dbExecute(conn, "COMMIT")
  committed <- TRUE

  message(msg)
  invisible(pop)
}


#' Copy counts of one chromosome for every (sex, line) already in `ind_meta`
#'
#' Resolved exactly as [add_founders()] and [add_offspring()] resolve them, so
#' the comparison sees what those animals actually carry.
#'
#' @return Data frame `sex, line_name, from_parent_1, from_parent_2`, ordered.
#' @keywords internal
#' @noRd
.dc_karyotypes <- function(conn, chr_name) {
  groups <- DBI::dbGetQuery(conn, paste(
    "SELECT DISTINCT sex, line_name FROM ind_meta",
    "ORDER BY sex NULLS FIRST, line_name NULLS FIRST"))
  out <- groups
  out$from_parent_1 <- integer(nrow(groups))
  out$from_parent_2 <- integer(nrow(groups))
  for (i in seq_len(nrow(groups))) {
    sx <- if (is.na(groups$sex[i])) NULL else groups$sex[i]
    ln <- if (is.na(groups$line_name[i])) NULL else groups$line_name[i]
    r  <- resolve_chr_inheritance(conn, sx, ln)
    r  <- r[r$chr_name == chr_name, , drop = FALSE]
    out$from_parent_1[i] <- r$from_parent_1
    out$from_parent_2[i] <- r$from_parent_2
  }
  out
}

#' Refuse an inheritance rule that changes animals already written
#' @keywords internal
#' @noRd
.dc_refuse_changed_karyotypes <- function(before, after, chr_name) {
  changed <- before$from_parent_1 != after$from_parent_1 |
             before$from_parent_2 != after$from_parent_2
  if (!any(changed)) return(invisible(NULL))
  b <- before[changed, , drop = FALSE]
  a <- after[changed, , drop = FALSE]
  lab <- function(x) if (is.na(x)) "NULL" else x
  what <- vapply(seq_len(nrow(b)), function(i) sprintf(
    "sex %s, line %s: %d,%d -> %d,%d", lab(b$sex[i]), lab(b$line_name[i]),
    b$from_parent_1[i], b$from_parent_2[i],
    a$from_parent_1[i], a$from_parent_2[i]), character(1))
  stop("define_chromosome() would change the copies of chromosome '", chr_name,
       "' that existing individuals already carry (from_parent_1,",
       "from_parent_2): ", paste(utils::head(what, 5L), collapse = "; "),
       if (length(what) > 5L) "; ..." else "", ". Their haplotypes were ",
       "written under the old rule, so breeding from them would fail. Set ",
       "chromosome inheritance rules before add_founders(), in a new ",
       "population. Nothing was written.", call. = FALSE)
}
