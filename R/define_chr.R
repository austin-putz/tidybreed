#' Define a chromosome's inheritance rule
#'
#' @description
#' Sets a non-default inheritance rule for one chromosome in `chr_meta` — sex
#' chromosomes (X/Y, Z/W, X0/Z0) and organelles (mitochondria, plastids).
#' `initialize_genome()`/`define_genome()` seeds every chromosome with the
#' default autosome rule (`copy_mode_M = copy_mode_F = "full"`,
#' `recombines = TRUE`); call `define_chr()` before [add_founders()] to
#' override that default for a specific chromosome.
#'
#' `copy_mode_M`/`copy_mode_F` are **relative to an individual's own ploidy**,
#' not an absolute copy count: `"full"` means the same as that individual's
#' ploidy for this chromosome, `"half"` means hemizygous (half of that
#' ploidy), and `"none"` means absent. In this version of tidybreed every
#' individual's ploidy is `2`, so `"full"` = 2 copies and `"half"` = 1 copy.
#'
#' @param pop A `tidybreed_pop` object.
#' @param chr_name Character scalar. Must already exist as a
#'   `genome_meta.chr_name` value (defined via [define_genome()]'s
#'   `chr_names` argument).
#' @param copy_mode_M Character scalar: `"full"` (default), `"half"`, or
#'   `"none"`. Copy count for males, relative to ploidy.
#' @param copy_mode_F Character scalar: `"full"` (default), `"half"`, or
#'   `"none"`. Copy count for females, relative to ploidy.
#' @param hemi_parent `NULL`, `"parent_1"`, or `"parent_2"`. Which parent
#'   supplies the reduced copy when either sex's copy_mode is not `"full"`.
#'   Must be `NULL` when both copy_modes are `"full"`; must be
#'   `"parent_1"`/`"parent_2"` otherwise.
#' @param recombines Logical. `TRUE` (default) if recombination occurs during
#'   gamete formation; `FALSE` for non-recombining chromosomes (Y, W, MT, most
#'   organelles).
#' @param overwrite Logical. If `TRUE` (default), re-calling `define_chr()`
#'   for the same `chr_name` updates the existing row in place — this is a
#'   normal edit workflow (e.g. fixing a typo'd `hemi_parent`), not append-only
#'   simulation output, so the default favors correction over silent no-op.
#'
#' @return The `tidybreed_pop` (invisibly). Assign the result back.
#'
#' @seealso [define_genome()]
#'
#' @examples
#' \dontrun{
#' # Mammalian sex chromosomes
#' pop <- pop |>
#'   define_chr("X", copy_mode_M = "half", copy_mode_F = "full",
#'              hemi_parent = "parent_2", recombines = TRUE) |>
#'   define_chr("Y", copy_mode_M = "half", copy_mode_F = "none",
#'              hemi_parent = "parent_1", recombines = FALSE)
#'
#' # Mitochondria (maternal, non-recombining)
#' pop <- pop |>
#'   define_chr("MT", copy_mode_M = "half", copy_mode_F = "half",
#'              hemi_parent = "parent_2", recombines = FALSE)
#' }
#' @export
define_chr <- function(pop, chr_name,
                       copy_mode_M = "full", copy_mode_F = "full",
                       hemi_parent = NULL, recombines = TRUE,
                       overwrite = TRUE) {

  validate_tidybreed_pop(pop)

  stopifnot(is.character(chr_name), length(chr_name) == 1L, nchar(chr_name) > 0L)

  valid_modes <- c("full", "half", "none")
  stopifnot(is.character(copy_mode_M), length(copy_mode_M) == 1L)
  stopifnot(is.character(copy_mode_F), length(copy_mode_F) == 1L)
  if (!copy_mode_M %in% valid_modes) {
    stop("copy_mode_M must be one of ", paste(valid_modes, collapse = ", "),
         ". Got '", copy_mode_M, "'.", call. = FALSE)
  }
  if (!copy_mode_F %in% valid_modes) {
    stop("copy_mode_F must be one of ", paste(valid_modes, collapse = ", "),
         ". Got '", copy_mode_F, "'.", call. = FALSE)
  }

  both_full <- copy_mode_M == "full" && copy_mode_F == "full"
  if (both_full && !is.null(hemi_parent)) {
    stop("hemi_parent must be NULL when copy_mode_M and copy_mode_F are both 'full'.",
         call. = FALSE)
  }
  if (!both_full && (is.null(hemi_parent) || length(hemi_parent) != 1L ||
                     !hemi_parent %in% c("parent_1", "parent_2"))) {
    stop("hemi_parent must be 'parent_1' or 'parent_2' when copy_mode_M or ",
         "copy_mode_F is not 'full'.", call. = FALSE)
  }

  stopifnot(is.logical(recombines), length(recombines) == 1L, !is.na(recombines))
  stopifnot(is.logical(overwrite), length(overwrite) == 1L, !is.na(overwrite))

  known_chrs <- DBI::dbGetQuery(pop$db_conn,
    "SELECT DISTINCT chr_name FROM genome_meta")$chr_name
  if (!chr_name %in% known_chrs) {
    stop("chr_name '", chr_name, "' not found in genome_meta. ",
         "Define it via define_genome()'s chr_names argument first.",
         call. = FALSE)
  }

  df <- data.frame(
    chr_name    = chr_name,
    copy_mode_M = copy_mode_M,
    copy_mode_F = copy_mode_F,
    hemi_parent = if (is.null(hemi_parent)) NA_character_ else hemi_parent,
    recombines  = recombines,
    stringsAsFactors = FALSE
  )

  # RNG-neutral write (duckdb_register + INSERT, not dbWriteTable/dbAppendTable
  # — define_chr() is typically called before add_founders(), so an RNG leak
  # here would shift the first stochastic draw of the simulation).
  tmp_tbl <- paste0("__define_chr_tmp_", as.integer(Sys.time()))
  duckdb::duckdb_register(pop$db_conn, tmp_tbl, df)
  on.exit(duckdb::duckdb_unregister(pop$db_conn, tmp_tbl), add = TRUE)

  conflict_action <- if (overwrite) {
    paste0(
      "DO UPDATE SET copy_mode_M = EXCLUDED.copy_mode_M, ",
      "copy_mode_F = EXCLUDED.copy_mode_F, ",
      "hemi_parent = EXCLUDED.hemi_parent, ",
      "recombines = EXCLUDED.recombines"
    )
  } else {
    "DO NOTHING"
  }

  DBI::dbExecute(pop$db_conn, paste0(
    "INSERT INTO chr_meta SELECT * FROM ", tmp_tbl,
    " ON CONFLICT (chr_name) ", conflict_action
  ))

  message("Defined chr_meta rule for '", chr_name, "': M=", copy_mode_M,
          ", F=", copy_mode_F, ", hemi_parent=",
          if (is.null(hemi_parent)) "NULL" else hemi_parent,
          ", recombines=", recombines)

  invisible(pop)
}
