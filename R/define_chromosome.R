#' Define a chromosome's inheritance or recombination rule
#'
#' @description
#' Sets a non-default rule for one chromosome — sex chromosomes (X/Y, Z/W, X0/Z0)
#' and organelles (mitochondria, plastids). [define_genome()] seeds every
#' chromosome with the default diploid-autosome rule (one copy from each parent,
#' recombining in both sexes); call `define_chromosome()` before [add_founders()]
#' to override that default.
#'
#' The rule is stored across two explicit long tables, and **each call sets
#' exactly one concern**:
#'
#' * **Inheritance** (copy count): supply `from_parent_1` and `from_parent_2` —
#'   the number of copies of this chromosome an offspring of `offspring_sex`
#'   inherits from parent_1 (sire) and parent_2 (dam). Writes `chr_inheritance`.
#' * **Recombination**: supply `recombines` — whether a parent of `parent_sex`
#'   recombines this chromosome when making gametes. Writes `chr_recombination`.
#'
#' Copy counts are **absolute** (correct at ploidy 2, which this version
#' enforces): an autosome is `from_parent_1 = 1, from_parent_2 = 1`; a male's X is
#' `0, 1`; a male's Y is `1, 0`; a female's Y is `0, 0`; maternal mitochondria are
#' `0, 1`. The two concerns are kept in separate calls because "which sex" means
#' the **offspring** for copy count but the **producing parent** for
#' recombination — mixing them in one call would make that ambiguous.
#'
#' @param pop A `tidybreed_pop` object.
#' @param chr_name Character scalar. Must already exist as a `genome_meta.chr_name`
#'   value (defined via [define_genome()]'s `chr_names` argument).
#' @param offspring_sex `NULL` (default; the sex-agnostic default row), `"M"`, or
#'   `"F"`. The **inheritance** row key — the offspring sex whose copy counts this
#'   call sets. Only valid on an inheritance call.
#' @param parent_sex `NULL` (default; both parent sexes), `"M"`, or `"F"`. The
#'   **recombination** row key — the producing-parent sex whose recombination this
#'   call sets. Only valid on a recombination call.
#' @param line_name `NULL` (default; all lines) or a line name. Reserved for
#'   line-specific rules (crossbreeding); resolves with the same precedence as the
#'   sex dimension.
#' @param from_parent_1,from_parent_2 Non-negative integers supplied **together**
#'   on an inheritance call. Copies inherited from parent_1 (sire) and parent_2
#'   (dam); their sum must be `<= 2` in this diploid release.
#' @param recombines Single logical supplied on a recombination call. `TRUE` if
#'   the chromosome recombines during gamete formation, `FALSE` for
#'   non-recombining chromosomes (Y, W, MT, most organelles). For a genome-wide
#'   achiasmatic sex, set `recombines_M`/`recombines_F` on [define_genome()]
#'   instead of one call per chromosome.
#' @param overwrite Logical. `TRUE` (default) upserts the row for this logical
#'   key. `FALSE` errors if the exact `(chr_name, sex, line_name)` key already
#'   exists (note: [define_genome()] pre-seeds the `(chr, NULL, NULL)` default
#'   row, so `overwrite = FALSE` on a default-keyed call always collides — it is
#'   meant for guarding *new* sex/line-specific rows).
#'
#' @return The `tidybreed_pop` (invisibly). Assign the result back.
#'
#' @seealso [define_genome()]
#'
#' @examples
#' \dontrun{
#' pop <- open_pop(pop_name = "test", db_name = ":memory:") |>
#'   define_genome(n_loci = 1000, n_chr = 3, chr_len_Mb = 100,
#'                 chr_names = c("1", "X", "Y"))
#'
#' # Mammalian sex chromosomes — inheritance (override only the deviating sexes)
#' pop <- pop |>
#'   define_chromosome("X", offspring_sex = "M", from_parent_1 = 0, from_parent_2 = 1) |>
#'   define_chromosome("Y", offspring_sex = "M", from_parent_1 = 1, from_parent_2 = 0) |>
#'   define_chromosome("Y", offspring_sex = "F", from_parent_1 = 0, from_parent_2 = 0) |>
#'   # Y is non-recombining (both parent sexes)
#'   define_chromosome("Y", recombines = FALSE)
#'
#' # Avian Z/W — females (ZW) are the heterogametic sex, males are ZZ
#' pop_zw <- open_pop(pop_name = "chicken", db_name = ":memory:") |>
#'   define_genome(n_loci = 1000, n_chr = 3, chr_len_Mb = 100,
#'                 chr_names = c("1", "Z", "W"))
#' pop_zw <- pop_zw |>
#'   define_chromosome("Z", offspring_sex = "F", from_parent_1 = 1, from_parent_2 = 0) |>
#'   define_chromosome("W", offspring_sex = "F", from_parent_1 = 0, from_parent_2 = 1) |>
#'   define_chromosome("W", offspring_sex = "M", from_parent_1 = 0, from_parent_2 = 0) |>
#'   # W is non-recombining (both parent sexes)
#'   define_chromosome("W", recombines = FALSE)
#'
#' # Mitochondria — maternal inheritance + non-recombining (separate calls)
#' pop2 <- open_pop(pop_name = "test2", db_name = ":memory:") |>
#'   define_genome(n_loci = 1000, n_chr = 2, chr_len_Mb = 100,
#'                 chr_names = c("1", "MT"))
#' pop2 <- pop2 |>
#'   define_chromosome("MT", from_parent_1 = 0, from_parent_2 = 1) |>
#'   define_chromosome("MT", recombines = FALSE)
#' }
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
      v >= 0 && v == as.integer(v)
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

  DBI::dbExecute(conn, "COMMIT")
  committed <- TRUE

  message(msg)
  invisible(pop)
}
