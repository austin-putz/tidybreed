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
#' IDs follow the same `{line}_{n}` format as `add_founders()`. Numbering
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
#'   add_founders(n_males = 5, n_females = 5, line_name = "A", gen = 1L)
#'
#' # One offspring per mating, extra metadata column (gen)
#' matings <- tibble::tibble(
#'   id_parent_1 = rep("A_1", 5),
#'   id_parent_2 = paste0("A_", 6:10),
#'   sex         = c("M", "F", "M", "F", "M"),
#'   line        = "A",
#'   gen         = 2L
#' )
#' pop <- pop |> add_offspring(matings)
#'
#' # Animal-breeder-style aliases
#' matings2 <- tibble::tibble(
#'   id_sire = rep("A_1", 3),
#'   id_dam  = paste0("A_", 6:8),
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

  required_cols <- c("id_parent_1", "id_parent_2", "sex", "line_name")
  missing_req   <- setdiff(required_cols, names(matings))

  if (length(missing_req) > 0) {
    stop(
      "matings is missing required column(s): ",
      paste(missing_req, collapse = ", "), ". ",
      "Required: id_parent_1 (or id_sire), id_parent_2 (or id_dam), sex, line_name",
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

  invalid_line <- !grepl("^[a-zA-Z][a-zA-Z0-9_]*$", matings$line_name)
  if (any(invalid_line)) {
    bad <- paste(unique(matings$line_name[invalid_line]), collapse = ", ")
    stop(
      "Invalid line value(s): ", bad, ". ",
      "Must start with a letter and contain only letters, numbers, or underscores.",
      call. = FALSE
    )
  }

  # ============================================================================
  # 5. Identify and validate extra columns
  # ============================================================================

  extra_cols <- setdiff(names(matings), required_cols)

  for (col in extra_cols) {
    validate_sql_identifier(col, what = "column name",
                           reserved = TABLE_RESERVED_COLS[["ind_meta"]])
  }

  # ============================================================================
  # 6. Validate required tables exist
  # ============================================================================

  needed_tables <- c("ind_meta", "genome_meta", "ind_haplotype")
  missing_tables <- setdiff(needed_tables, pop$tables)

  if (length(missing_tables) > 0) {
    stop(
      "Required table(s) not found: ", paste(missing_tables, collapse = ", "), ". ",
      "Ensure the population was created with initialize_genome().",
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

  genome_meta_df <- dplyr::tbl(pop$db_conn, "genome_meta") |>
    dplyr::select(locus_id, chr, pos_Mb) |>
    dplyr::arrange(locus_id) |>
    dplyr::collect()

  n_loci     <- nrow(genome_meta_df)
  locus_cols <- paste0("locus_", seq_len(n_loci))
  chr_info   <- build_chr_info(genome_meta_df)

  # ============================================================================
  # 9. Load parent haplotypes in a single batch query
  # ============================================================================

  # Load parent haplotypes from the long ind_haplotype table so each allele's
  # line_origin travels with it. Rows come back grouped by (id_ind,
  # parent_origin, locus_id); for a diploid parent that is exactly 2 x n_loci
  # rows, giving a 2 x n_loci allele matrix (row 1 = the parent's own
  # parent_origin 1 homolog, row 2 = its parent_origin 2 homolog) and a parallel
  # line_origin matrix. Reading via dbGetQuery does not advance R's RNG.
  parent_haps_raw <- DBI::dbGetQuery(
    pop$db_conn,
    paste0(
      "SELECT id_ind, parent_origin, locus_id, allele, line_origin ",
      "FROM ind_haplotype WHERE id_ind IN (", parent_id_list,
      ") ORDER BY id_ind, parent_origin, locus_id"
    )
  )

  parent_haps <- vector("list", length(unique_parents))
  parent_lo   <- vector("list", length(unique_parents))
  names(parent_haps) <- unique_parents
  names(parent_lo)   <- unique_parents

  for (pid in unique_parents) {
    rows <- parent_haps_raw[parent_haps_raw$id_ind == pid, ]
    if (nrow(rows) != 2L * n_loci) {
      stop(
        "Parent '", pid, "' does not have exactly 2 x n_loci rows in ind_haplotype.",
        call. = FALSE
      )
    }
    # byrow: the first n_loci rows (parent_origin 1, locus_id order) fill row 1,
    # the next n_loci (parent_origin 2) fill row 2.
    hap_mat <- matrix(as.integer(rows$allele), nrow = 2L, ncol = n_loci, byrow = TRUE)
    lo_mat  <- matrix(rows$line_origin,        nrow = 2L, ncol = n_loci, byrow = TRUE)
    parent_haps[[pid]] <- hap_mat
    parent_lo[[pid]]   <- lo_mat
  }

  # ============================================================================
  # 10. Determine offspring IDs (per-line sequential numbering)
  # ============================================================================

  unique_lines <- unique(matings$line_name)
  line_start   <- integer(length(unique_lines))
  names(line_start) <- unique_lines

  for (ln in unique_lines) {
    res <- DBI::dbGetQuery(
      pop$db_conn,
      paste0(
        "SELECT MAX(CAST(SUBSTRING(id_ind FROM POSITION('_' IN id_ind) + 1) AS INTEGER)) AS max_num ",
        "FROM ind_meta WHERE LEFT(id_ind, ", nchar(ln) + 1L, ") = '", ln, "_'"
      )
    )
    line_start[ln] <- if (is.na(res$max_num)) 1L else as.integer(res$max_num) + 1L
  }

  n_offspring  <- nrow(matings)
  offspring_ids <- character(n_offspring)
  line_counter  <- as.list(line_start)

  for (i in seq_len(n_offspring)) {
    ln             <- matings$line_name[i]
    offspring_ids[i] <- paste0(ln, "_", line_counter[[ln]])
    line_counter[[ln]] <- line_counter[[ln]] + 1L
  }

  # ============================================================================
  # 11. Generate gametes and build genomic matrices
  # ============================================================================

  hap_sire_mat    <- matrix(0L, nrow = n_offspring, ncol = n_loci)
  hap_dam_mat     <- matrix(0L, nrow = n_offspring, ncol = n_loci)
  hap_sire_lo_mat <- matrix(NA_character_, nrow = n_offspring, ncol = n_loci)
  hap_dam_lo_mat  <- matrix(NA_character_, nrow = n_offspring, ncol = n_loci)
  col_idx <- seq_len(n_loci)

  for (i in seq_len(n_offspring)) {
    sire <- matings$id_parent_1[i]
    dam  <- matings$id_parent_2[i]

    # make_gamete() returns allele + the per-locus homolog (1/2) it came from.
    # line_origin follows whichever parental homolog donated each locus, so it
    # is correct across F1/F2/backcross generations.
    gs <- make_gamete(parent_haps[[sire]], chr_info)
    hap_sire_mat[i, ]    <- gs$allele
    hap_sire_lo_mat[i, ] <- parent_lo[[sire]][cbind(gs$homolog, col_idx)]

    gd <- make_gamete(parent_haps[[dam]], chr_info)
    hap_dam_mat[i, ]    <- gd$allele
    hap_dam_lo_mat[i, ] <- parent_lo[[dam]][cbind(gd$homolog, col_idx)]
  }

  # ============================================================================
  # 12. Build the wide allele frame (transient; the vehicle for the UNPIVOT into
  #     the long ind_haplotype). Never written to a wide DB table.
  # ============================================================================

  hap_wide <- dplyr::bind_rows(
    dplyr::bind_cols(
      tibble::tibble(id_ind = offspring_ids, parent_origin = 1L),
      setNames(as.data.frame(hap_sire_mat), locus_cols)
    ),
    dplyr::bind_cols(
      tibble::tibble(id_ind = offspring_ids, parent_origin = 2L),
      setNames(as.data.frame(hap_dam_mat), locus_cols)
    )
  )

  # ============================================================================
  # 14. Build ind_meta rows and handle extra columns
  # ============================================================================

  ind_meta_new <- tibble::tibble(
    id_ind      = offspring_ids,
    id_parent_1 = matings$id_parent_1,
    id_parent_2 = matings$id_parent_2,
    line_name   = matings$line_name,
    sex         = matings$sex
  )

  if (length(extra_cols) > 0) {
    # Build a named list of per-offspring values from matings columns
    extra_list <- stats::setNames(
      lapply(extra_cols, function(col) matings[[col]]),
      extra_cols
    )
    prepped <- prepare_extra_cols(extra_list, n_offspring, "ind_meta", pop$db_conn)
    for (nm in names(prepped)) ind_meta_new[[nm]] <- prepped[[nm]]
  }

  # ============================================================================
  # 15. Write all tables to database
  # ============================================================================

  DBI::dbWriteTable(pop$db_conn, "ind_meta", ind_meta_new, append = TRUE)

  # Write the long ind_haplotype. A wide allele frame and a wide line_origin
  # frame are unpivoted and joined in DuckDB (per (id_ind, parent_origin,
  # locus)) so R never materializes the long frame. ind_genotype stays
  # on-demand (add_dosage()).
  lo_df <- dplyr::bind_rows(
    dplyr::bind_cols(
      tibble::tibble(id_ind = offspring_ids, parent_origin = 1L),
      setNames(as.data.frame(hap_sire_lo_mat, stringsAsFactors = FALSE), locus_cols)
    ),
    dplyr::bind_cols(
      tibble::tibble(id_ind = offspring_ids, parent_origin = 2L),
      setNames(as.data.frame(hap_dam_lo_mat, stringsAsFactors = FALSE), locus_cols)
    )
  )
  duckdb::duckdb_register(pop$db_conn, "__tmp_hap", as.data.frame(hap_wide))
  duckdb::duckdb_register(pop$db_conn, "__tmp_lo",  lo_df)
  DBI::dbExecute(pop$db_conn, paste0(
    "INSERT INTO ind_haplotype ",
    "(id_ind, parent_origin, strand, line_origin, locus_id, locus_name, allele) ",
    "SELECT a.id_ind, a.parent_origin, 1, l.line_origin, gm.locus_id, gm.locus_name, a.allele ",
    "FROM (UNPIVOT __tmp_hap ON COLUMNS(* EXCLUDE (id_ind, parent_origin)) ",
    "INTO NAME locus_col VALUE allele) a ",
    "JOIN (UNPIVOT __tmp_lo ON COLUMNS(* EXCLUDE (id_ind, parent_origin)) ",
    "INTO NAME locus_col VALUE line_origin) l ",
    "ON a.id_ind = l.id_ind AND a.parent_origin = l.parent_origin ",
    "AND a.locus_col = l.locus_col ",
    "JOIN genome_meta gm ON a.locus_col = 'locus_' || gm.locus_id"
  ))
  duckdb::duckdb_unregister(pop$db_conn, "__tmp_hap")
  duckdb::duckdb_unregister(pop$db_conn, "__tmp_lo")

  # ============================================================================
  # 16. Return
  # ============================================================================

  message("Added ", n_offspring, " offspring")

  invisible(pop)
}
