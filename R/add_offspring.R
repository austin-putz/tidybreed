#' Add offspring via recombination
#'
#' @description
#' Produces offspring by simulating chromosomal crossovers from explicit parent
#' pairs. Each row of `matings` defines exactly one offspring, giving full
#' control over the mating design.
#'
#' @param pop A `tidybreed_pop` object
#' @param matings A tibble or data.frame where each row = one offspring.
#'
#'   **Required columns** (canonical names match `ind_meta`):
#'   - `id_parent_1` — sire / parent 1 ID; must exist in `ind_meta`.
#'     Alias: `id_sire`
#'   - `id_parent_2` — dam / parent 2 ID; must exist in `ind_meta`.
#'     Alias: `id_dam`
#'   - `sex` — `"M"` or `"F"`
#'   - `line` — line name for offspring IDs (same format as `add_founders()`)
#'
#'   **Optional extra columns** (e.g. `gen = 2L`, `farm = "Iowa"`) are
#'   validated and written directly to `ind_meta`. If a column does not yet
#'   exist in `ind_meta` it is added automatically. Scalar values in a tibble
#'   are recycled to all rows before being passed to this function.
#'
#' @return The modified `tidybreed_pop` object (invisibly).
#'   Assign the result back: `pop <- add_offspring(pop, matings)`
#'
#' @details
#' **Mating design flexibility:**
#' - Multiple offspring per pair: repeat the pair row in `matings`
#' - Multiple sires per dam (pooled semen / polyspermy): include each sire
#'   as a separate row with the same `id_parent_2`
#' - Multiple dams per sire: include each dam as a separate row with the
#'   same `id_parent_1`
#' - Cross-line matings: use any `line` value in the `line` column
#'
#' **Column aliases:**
#' `id_sire` is accepted as an alias for `id_parent_1`, and `id_dam` for
#' `id_parent_2`. Both naming styles produce identical results.
#'
#' **Recombination model:**
#' Gametes are simulated using the Haldane map function. For chromosome i,
#' the number of crossovers ~ Poisson(chr_len_Mb\[i\] / 100), assuming
#' approximately 1 Morgan per 100 Mb. Crossover positions are uniform within
#' each chromosome, and the starting haplotype is chosen at random.
#'
#' **Offspring IDs:**
#' IDs follow the same `{line}-{n}` format as `add_founders()`. Numbering
#' continues from the current maximum for each line.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' pop <- initialize_genome(
#'   pop_name = "test", n_loci = 1000, n_chr = 5,
#'   chr_len_Mb = 100, n_haplotypes = 200
#' )
#' pop <- pop |>
#'   add_founders(n_males = 5, n_females = 5, line_name = "A") |>
#'   mutate_ind_meta(gen = 1L)
#'
#' # One offspring per mating, extra metadata column (gen)
#' matings <- tibble::tibble(
#'   id_parent_1 = rep("A-1", 5),
#'   id_parent_2 = paste0("A-", 6:10),
#'   sex         = c("M", "F", "M", "F", "M"),
#'   line        = "A",
#'   gen         = 2L
#' )
#' pop <- pop |> add_offspring(matings)
#'
#' # Animal-breeder-style aliases
#' matings2 <- tibble::tibble(
#'   id_sire = rep("A-1", 3),
#'   id_dam  = paste0("A-", 6:8),
#'   sex     = c("M", "F", "M"),
#'   line    = "A",
#'   gen     = 2L
#' )
#' pop <- pop |> add_offspring(matings2)
#' }
add_offspring <- function(pop, matings) {

  # ============================================================================
  # 1. Validate pop and matings types
  # ============================================================================

  stopifnot(inherits(pop, "tidybreed_pop"))
  validate_tidybreed_pop(pop)

  if (!is.data.frame(matings)) {
    stop("matings must be a data.frame or tibble", call. = FALSE)
  }

  if (nrow(matings) == 0) {
    stop("matings must have at least 1 row", call. = FALSE)
  }

  # ============================================================================
  # 2. Handle column aliases (id_sire -> id_parent_1, id_dam -> id_parent_2)
  # ============================================================================

  mat_names <- names(matings)

  if ("id_sire" %in% mat_names && !"id_parent_1" %in% mat_names) {
    names(matings)[names(matings) == "id_sire"] <- "id_parent_1"
  }
  if ("id_dam" %in% mat_names && !"id_parent_2" %in% mat_names) {
    names(matings)[names(matings) == "id_dam"] <- "id_parent_2"
  }

  # ============================================================================
  # 3. Check required columns are present
  # ============================================================================

  required_cols <- c("id_parent_1", "id_parent_2", "sex", "line")
  missing_req   <- setdiff(required_cols, names(matings))

  if (length(missing_req) > 0) {
    stop(
      "matings is missing required column(s): ",
      paste(missing_req, collapse = ", "), ". ",
      "Required: id_parent_1 (or id_sire), id_parent_2 (or id_dam), sex, line",
      call. = FALSE
    )
  }

  # ============================================================================
  # 4. Validate required column values
  # ============================================================================

  invalid_sex <- !matings$sex %in% c("M", "F")
  if (any(invalid_sex)) {
    bad <- paste(unique(matings$sex[invalid_sex]), collapse = ", ")
    stop(
      "Invalid sex value(s): ", bad, ". Must be 'M' or 'F'.",
      call. = FALSE
    )
  }

  invalid_line <- !grepl("^[a-zA-Z][a-zA-Z0-9_-]*$", matings$line)
  if (any(invalid_line)) {
    bad <- paste(unique(matings$line[invalid_line]), collapse = ", ")
    stop(
      "Invalid line value(s): ", bad, ". ",
      "Must start with a letter and contain only letters, numbers, underscores, or hyphens.",
      call. = FALSE
    )
  }

  # ============================================================================
  # 5. Identify and validate extra columns
  # ============================================================================

  extra_cols <- setdiff(names(matings), required_cols)

  for (col in extra_cols) {
    validate_field_name(col, existing_cols = character(0))
  }

  # ============================================================================
  # 6. Validate required tables exist
  # ============================================================================

  needed_tables <- c("ind_meta", "genome_meta", "genome_haplotype", "genome_genotype")
  missing_tables <- setdiff(needed_tables, pop$tables)

  if (length(missing_tables) > 0) {
    extra_hint <- if ("ind_meta" %in% missing_tables) " Call add_founders() first." else ""
    stop(
      "Required table(s) not found: ", paste(missing_tables, collapse = ", "), ".",
      extra_hint,
      call. = FALSE
    )
  }

  # ============================================================================
  # 7. Validate all parent IDs exist in ind_meta
  # ============================================================================

  unique_parents <- unique(c(matings$id_parent_1, matings$id_parent_2))
  parent_id_list <- paste(paste0("'", unique_parents, "'"), collapse = ", ")

  existing_ids <- DBI::dbGetQuery(
    pop$db_conn,
    paste0("SELECT id_ind FROM ind_meta WHERE id_ind IN (", parent_id_list, ")")
  )$id_ind

  missing_parents <- setdiff(unique_parents, existing_ids)
  if (length(missing_parents) > 0) {
    stop(
      "Parent ID(s) not found in ind_meta: ",
      paste(missing_parents, collapse = ", "),
      call. = FALSE
    )
  }

  # ============================================================================
  # 8. Load genome metadata (once)
  # ============================================================================

  n_loci    <- pop$metadata$n_loci
  locus_cols <- paste0("locus_", seq_len(n_loci))
  chr_len_Mb <- pop$metadata$chr_len_Mb

  genome_meta_df <- get_table(pop, "genome_meta") |>
    dplyr::select(locus_id, chr, pos_Mb) |>
    dplyr::arrange(locus_id) |>
    dplyr::collect()

  # ============================================================================
  # 9. Load parent haplotypes in a single batch query
  # ============================================================================

  parent_haps_raw <- DBI::dbGetQuery(
    pop$db_conn,
    paste0(
      "SELECT * FROM genome_haplotype WHERE id_ind IN (", parent_id_list,
      ") ORDER BY id_ind, parent_origin"
    )
  )

  parent_haps <- vector("list", length(unique_parents))
  names(parent_haps) <- unique_parents

  for (pid in unique_parents) {
    rows <- parent_haps_raw[parent_haps_raw$id_ind == pid, ]
    rows <- rows[order(rows$parent_origin), ]

    if (nrow(rows) != 2L) {
      stop(
        "Parent '", pid, "' does not have exactly 2 rows in genome_haplotype.",
        call. = FALSE
      )
    }

    hap_mat <- as.matrix(rows[, locus_cols])
    storage.mode(hap_mat) <- "integer"
    parent_haps[[pid]] <- hap_mat
  }

  # ============================================================================
  # 10. Determine offspring IDs (per-line sequential numbering)
  # ============================================================================

  unique_lines <- unique(matings$line)
  line_start   <- integer(length(unique_lines))
  names(line_start) <- unique_lines

  for (ln in unique_lines) {
    res <- DBI::dbGetQuery(
      pop$db_conn,
      paste0(
        "SELECT MAX(CAST(SUBSTRING(id_ind FROM POSITION('-' IN id_ind) + 1) AS INTEGER)) AS max_num ",
        "FROM ind_meta WHERE id_ind LIKE '", ln, "-%'"
      )
    )
    line_start[ln] <- if (is.na(res$max_num)) 1L else as.integer(res$max_num) + 1L
  }

  n_offspring  <- nrow(matings)
  offspring_ids <- character(n_offspring)
  line_counter  <- as.list(line_start)

  for (i in seq_len(n_offspring)) {
    ln             <- matings$line[i]
    offspring_ids[i] <- paste0(ln, "-", line_counter[[ln]])
    line_counter[[ln]] <- line_counter[[ln]] + 1L
  }

  # ============================================================================
  # 11. Generate gametes and build genomic matrices
  # ============================================================================

  hap_sire_mat <- matrix(0L, nrow = n_offspring, ncol = n_loci)
  hap_dam_mat  <- matrix(0L, nrow = n_offspring, ncol = n_loci)

  for (i in seq_len(n_offspring)) {
    hap_sire_mat[i, ] <- make_gamete(parent_haps[[matings$id_parent_1[i]]], genome_meta_df, chr_len_Mb)
    hap_dam_mat[i, ]  <- make_gamete(parent_haps[[matings$id_parent_2[i]]], genome_meta_df, chr_len_Mb)
  }

  # ============================================================================
  # 12. Build genome_haplotype data frame (2 rows per offspring)
  # ============================================================================

  sire_locus_df <- setNames(as.data.frame(hap_sire_mat), locus_cols)
  dam_locus_df  <- setNames(as.data.frame(hap_dam_mat),  locus_cols)

  genome_haplotype_df <- dplyr::bind_rows(
    dplyr::bind_cols(
      tibble::tibble(id_ind = offspring_ids, parent_origin = 1L),
      sire_locus_df
    ),
    dplyr::bind_cols(
      tibble::tibble(id_ind = offspring_ids, parent_origin = 2L),
      dam_locus_df
    )
  )

  # ============================================================================
  # 13. Build genome_genotype data frame (1 row per offspring)
  # ============================================================================

  geno_locus_df <- setNames(as.data.frame(hap_sire_mat + hap_dam_mat), locus_cols)

  genome_genotype_df <- dplyr::bind_cols(
    tibble::tibble(id_ind = offspring_ids),
    geno_locus_df
  )

  # ============================================================================
  # 14. Build ind_meta rows and handle extra columns
  # ============================================================================

  ind_meta_new <- tibble::tibble(
    id_ind      = offspring_ids,
    id_parent_1 = matings$id_parent_1,
    id_parent_2 = matings$id_parent_2,
    line        = matings$line,
    sex         = matings$sex
  )

  if (length(extra_cols) > 0) {
    existing_ind_meta_cols <- DBI::dbListFields(pop$db_conn, "ind_meta")

    for (col in extra_cols) {
      col_val <- matings[[col]]
      db_type <- infer_duckdb_type(col_val)

      if (!col %in% existing_ind_meta_cols) {
        DBI::dbExecute(pop$db_conn, paste0(
          "ALTER TABLE ind_meta ADD COLUMN ", col, " ", db_type
        ))
      }

      ind_meta_new[[col]] <- col_val
    }
  }

  # ============================================================================
  # 15. Write all tables to database
  # ============================================================================

  DBI::dbWriteTable(pop$db_conn, "ind_meta",         ind_meta_new,        append = TRUE)
  DBI::dbWriteTable(pop$db_conn, "genome_haplotype",  genome_haplotype_df, append = TRUE)
  DBI::dbWriteTable(pop$db_conn, "genome_genotype",   genome_genotype_df,  append = TRUE)

  # ============================================================================
  # 16. Update metadata and return
  # ============================================================================

  if (is.null(pop$metadata$n_individuals)) {
    pop$metadata$n_individuals <- n_offspring
  } else {
    pop$metadata$n_individuals <- pop$metadata$n_individuals + n_offspring
  }

  message("Added ", n_offspring, " offspring")

  invisible(pop)
}
