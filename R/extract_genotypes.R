#' Extract genotype data for individuals, by chip and/or QTL loci
#'
#' @description
#' Returns a tibble of genotypes (0/1/2 encoding) for a set of individuals,
#' restricted to loci selected by a chip definition, a filtered
#' `genome_effects` table, or both. Pipe a `tidybreed_table` (from
#' [get_table()] and optionally [filter()]) as the first argument to restrict
#' which individuals are included.
#'
#' At least one of `chip_name` or `effects_tbl` must be provided. When both
#' are given the locus sets are **unioned** (deduplicated, ordered by
#' `locus_id`).
#'
#' **Chip path** — the returned individual set is the intersection of:
#' * Animals with `has_<chip_name> == TRUE` in `ind_meta`
#' * Animals matching any pending `filter()` predicates on `tbl`
#' * Loci with `is_<chip_name> == TRUE` in `genome_meta`
#'
#' **QTL path (`effects_tbl`)** — the returned individual set is:
#' * Animals matching any pending `filter()` predicates on `tbl`
#'   (or all individuals in `genome_genotype` when no filter is applied)
#' * Loci whose `locus_name` appears in the collected `effects_tbl`
#'
#' @param tbl A `tidybreed_table` object from [get_table()] (optionally piped
#'   through [filter()]). The table must contain an `id_ind` column when a
#'   filter is applied.
#' @param chip_name Character or `NULL`. Name of a chip previously defined via
#'   [define_chip()] and applied to animals via [add_genotypes()]. When `NULL`
#'   the chip path is skipped.
#' @param effects_tbl A `tidybreed_table` from
#'   `get_table(pop, "genome_effects")` (optionally filtered), or `NULL`.
#'   The collected table must contain a `locus_name` column. Use
#'   [dplyr::filter()] to restrict by `trait_name`, `genome_effect_type`,
#'   `genome_value`, `line_name`, etc. When `NULL` the QTL path is skipped.
#' @param col_name Character. Name of the BOOLEAN column in `ind_meta` that
#'   records chip genotyping status. Default: `paste0("has_", chip_name)`.
#'   Ignored when `chip_name` is `NULL`.
#'
#' @return A tibble with one row per individual and one column per selected
#'   locus (`id_ind` + `locus_N` columns). Locus columns use 0/1/2 integer
#'   encoding and are ordered by `locus_id`.
#'
#' @examples
#' \dontrun{
#' # All genotyped animals on the 50k chip (unchanged usage)
#' geno <- pop |> get_table("ind_meta") |> extract_genotypes("50k")
#'
#' # Only females genotyped on the HD chip
#' geno <- pop |>
#'   get_table("ind_meta") |>
#'   dplyr::filter(sex == "F") |>
#'   extract_genotypes("HD")
#'
#' # QTL loci for a trait (all individuals)
#' geno <- pop |>
#'   get_table("ind_meta") |>
#'   extract_genotypes(
#'     effects_tbl = get_table(pop, "genome_effects") |>
#'       dplyr::filter(trait_name == "ADG")
#'   )
#'
#' # Large-effect QTL only, females only
#' geno <- pop |>
#'   get_table("ind_meta") |>
#'   dplyr::filter(sex == "F") |>
#'   extract_genotypes(
#'     effects_tbl = get_table(pop, "genome_effects") |>
#'       dplyr::filter(trait_name == "ADG", abs(genome_value) > 0.5)
#'   )
#'
#' # Chip loci + QTL loci unioned
#' geno <- pop |>
#'   get_table("ind_meta") |>
#'   extract_genotypes(
#'     chip_name   = "50k",
#'     effects_tbl = get_table(pop, "genome_effects") |>
#'       dplyr::filter(trait_name == "ADG")
#'   )
#' }
#' @export
extract_genotypes <- function(tbl,
                              chip_name   = NULL,
                              effects_tbl = NULL,
                              col_name    = NULL) {

  stopifnot(inherits(tbl, "tidybreed_table"))
  pop <- tbl$pop
  validate_tidybreed_pop(pop)

  if (is.null(chip_name) && is.null(effects_tbl)) {
    stop("At least one of 'chip_name' or 'effects_tbl' must be provided.",
         call. = FALSE)
  }

  # --- chip_name validation ---
  if (!is.null(chip_name)) {
    stopifnot(is.character(chip_name), length(chip_name) == 1L, nzchar(chip_name))
    if (is.null(col_name)) col_name <- paste0("has_", chip_name)
    validate_sql_identifier(col_name, what = "col_name")

    chip_col <- paste0("is_", chip_name)
    if (!chip_col %in% DBI::dbListFields(pop$db_conn, "genome_meta")) {
      stop("Chip '", chip_name, "' not found in genome_meta. ",
           "Call define_chip() first.", call. = FALSE)
    }
    if (!col_name %in% DBI::dbListFields(pop$db_conn, "ind_meta")) {
      stop("Column '", col_name, "' not found in ind_meta. ",
           "Call add_genotypes('", chip_name, "') first.", call. = FALSE)
    }
  }

  # --- effects_tbl validation ---
  if (!is.null(effects_tbl)) {
    if (!inherits(effects_tbl, "tidybreed_table"))
      stop("'effects_tbl' must be a tidybreed_table from get_table('genome_effects').",
           call. = FALSE)
    if (!identical(effects_tbl$table_name, "genome_effects"))
      stop("'effects_tbl' must be from get_table('genome_effects'), got '",
           effects_tbl$table_name, "'.", call. = FALSE)
  }

  # --- Resolve pending individual filter from tbl ---
  if (length(tbl$pending_filter) == 0) {
    subset_ids <- NULL
  } else {
    collected <- dplyr::collect(tbl)
    if (!"id_ind" %in% names(collected)) {
      stop("Filtered table '", tbl$table_name,
           "' must contain 'id_ind' to subset individuals for genotype extraction.",
           call. = FALSE)
    }
    subset_ids <- unique(collected[["id_ind"]])
  }

  # --- Resolve final individual IDs ---
  if (!is.null(chip_name)) {
    # Chip path: intersect has_<chip> == TRUE with any pending filter
    if (is.null(subset_ids)) {
      has_ids <- DBI::dbGetQuery(
        pop$db_conn,
        paste0("SELECT id_ind FROM ind_meta WHERE ", col_name, " = TRUE")
      )$id_ind
    } else {
      ids_sql <- paste0("'", subset_ids, "'", collapse = ", ")
      has_ids <- DBI::dbGetQuery(
        pop$db_conn,
        paste0("SELECT id_ind FROM ind_meta WHERE ", col_name, " = TRUE ",
               "AND id_ind IN (", ids_sql, ")")
      )$id_ind
    }
    if (length(has_ids) == 0) {
      stop("No genotyped animals found for chip '", chip_name,
           "' in the filtered set. Call add_genotypes() first.", call. = FALSE)
    }
  } else {
    # QTL-only path: use pending filter or all individuals in genome_genotype
    if (is.null(subset_ids)) {
      has_ids <- DBI::dbGetQuery(
        pop$db_conn,
        "SELECT id_ind FROM genome_genotype"
      )$id_ind
    } else {
      has_ids <- subset_ids
    }
    if (length(has_ids) == 0) {
      stop("No individuals found in the filtered set.", call. = FALSE)
    }
  }

  # --- Resolve locus IDs (union of chip and/or QTL paths) ---
  locus_ids <- integer(0)

  if (!is.null(chip_name)) {
    chip_ids <- DBI::dbGetQuery(
      pop$db_conn,
      paste0("SELECT locus_id FROM genome_meta WHERE ", chip_col,
             " = TRUE ORDER BY locus_id")
    )$locus_id
    if (length(chip_ids) == 0)
      stop("No loci found on chip '", chip_name, "'.", call. = FALSE)
    locus_ids <- union(locus_ids, chip_ids)
  }

  if (!is.null(effects_tbl)) {
    effects_df <- dplyr::collect(effects_tbl)
    if (!"locus_name" %in% names(effects_df))
      stop("Collected 'effects_tbl' must contain a 'locus_name' column.",
           call. = FALSE)
    qtl_names <- unique(effects_df[["locus_name"]])
    if (length(qtl_names) == 0)
      stop("No QTL loci found in the filtered 'effects_tbl'. ",
           "Call add_additive_effects() first.", call. = FALSE)
    qtl_name_sql <- paste0("'", qtl_names, "'", collapse = ", ")
    qtl_ids <- DBI::dbGetQuery(
      pop$db_conn,
      paste0("SELECT locus_id FROM genome_meta WHERE locus_name IN (",
             qtl_name_sql, ") ORDER BY locus_id")
    )$locus_id
    locus_ids <- sort(union(locus_ids, qtl_ids))
  }

  locus_ids <- sort(locus_ids)

  # --- Dynamic SELECT: id_ind + locus columns, filtered to resolved individuals ---
  locus_cols <- paste0("locus_", locus_ids)
  cols_sql   <- paste(c("id_ind", locus_cols), collapse = ", ")
  ids_sql    <- paste0("'", has_ids, "'", collapse = ", ")

  df <- DBI::dbGetQuery(
    pop$db_conn,
    paste0("SELECT ", cols_sql,
           " FROM genome_genotype WHERE id_ind IN (", ids_sql, ")")
  )

  tibble::as_tibble(df)
}
