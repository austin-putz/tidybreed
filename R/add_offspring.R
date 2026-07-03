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
  # 7. Validate all parent IDs exist in ind_meta; fetch their sex and ploidy
  # ============================================================================

  unique_parents <- unique(c(matings$id_parent_1, matings$id_parent_2))
  parent_id_list <- paste(paste0("'", unique_parents, "'"), collapse = ", ")

  parent_meta <- DBI::dbGetQuery(
    pop$db_conn,
    paste0("SELECT id_ind, sex, ploidy FROM ind_meta WHERE id_ind IN (", parent_id_list, ")")
  )

  missing_parents <- setdiff(unique_parents, parent_meta$id_ind)
  if (length(missing_parents) > 0) {
    stop(
      "Parent ID(s) not found in ind_meta: ",
      paste(missing_parents, collapse = ", "),
      call. = FALSE
    )
  }

  parent_sex    <- stats::setNames(parent_meta$sex, parent_meta$id_ind)
  parent_ploidy <- stats::setNames(as.integer(parent_meta$ploidy), parent_meta$id_ind)

  # ============================================================================
  # 8. Load genome metadata (once); classify chromosomes as "plain autosome"
  #    (both sexes full copy, recombines — chr_meta's default row) vs.
  #    "special" (any sex-linked/organelle rule). Autosomes are handled by the
  #    unchanged fast path below; special chromosomes get their own branch,
  #    executed strictly AFTER the autosome path for each offspring so autosome
  #    RNG draws are never perturbed by chr_meta configuration.
  # ============================================================================

  genome_meta_df <- dplyr::tbl(pop$db_conn, "genome_meta") |>
    dplyr::select(locus_id, chr, chr_name, pos_Mb) |>
    dplyr::arrange(locus_id) |>
    dplyr::collect()

  n_loci     <- nrow(genome_meta_df)
  locus_cols <- paste0("locus_", seq_len(n_loci))
  chr_info   <- build_chr_info(genome_meta_df)

  chr_meta_map   <- get_chr_meta_map(pop$db_conn)
  chr_ids_all    <- sort(unique(genome_meta_df$chr))
  chr_id_keys    <- as.character(chr_ids_all)
  chr_name_by_id <- stats::setNames(
    vapply(chr_ids_all, function(x) genome_meta_df$chr_name[genome_meta_df$chr == x][1], character(1)),
    chr_id_keys
  )
  is_autosome_id <- vapply(chr_id_keys, function(k) is_plain_autosome(chr_meta_map[[chr_name_by_id[[k]]]]), logical(1))

  autosome_locus_idx <- which(genome_meta_df$chr %in% chr_ids_all[is_autosome_id])
  special_chr_names  <- unname(chr_name_by_id[chr_id_keys[!is_autosome_id]])

  # chr_info's locus_idx values reference positions in the FULL genome_meta_df
  # ordering. Parent haplotype matrices for the autosome path are built COMPACT
  # (n_autosome_loci wide, excluding special-chromosome loci entirely — see
  # step 9), so autosome_chr_info must be re-expressed with LOCAL indices into
  # that compact width. For the all-default case (no special chromosomes),
  # autosome_locus_idx is the identity 1:n_loci, so this remap is a no-op and
  # autosome_chr_info is byte-identical to the original chr_info — preserving
  # today's exact RNG sequence.
  full_to_local <- integer(nrow(genome_meta_df))
  full_to_local[autosome_locus_idx] <- seq_along(autosome_locus_idx)
  autosome_chr_info <- lapply(chr_info[chr_id_keys[is_autosome_id]], function(ci) {
    list(locus_idx = full_to_local[ci$locus_idx], pos_Mb = ci$pos_Mb, chr_len = ci$chr_len)
  })

  special_locus_ids <- stats::setNames(
    lapply(special_chr_names, function(cn) genome_meta_df$locus_id[genome_meta_df$chr_name == cn]),
    special_chr_names
  )
  # chr_info's locus_idx values reference positions in the FULL genome_meta_df
  # ordering; re-express each special chromosome's descriptor with LOCAL
  # (1..n_chr_loci) indices so it can be used against a compact per-chromosome
  # matrix (see parent_chr_haps below) without touching make_gamete() itself.
  special_chr_local_info <- stats::setNames(
    lapply(special_chr_names, function(cn) {
      orig <- chr_info[[chr_id_keys[chr_name_by_id == cn]]]
      list(list(locus_idx = seq_along(orig$locus_idx), pos_Mb = orig$pos_Mb, chr_len = orig$chr_len))
    }),
    special_chr_names
  )

  # ============================================================================
  # 9. Load parent haplotypes in a single batch query
  # ============================================================================

  # Load parent haplotypes from the long ind_haplotype table so each allele's
  # line_origin travels with it. Reading via dbGetQuery does not advance R's RNG.
  parent_haps_raw <- DBI::dbGetQuery(
    pop$db_conn,
    paste0(
      "SELECT id_ind, parent_origin, locus_id, allele, line_origin ",
      "FROM ind_haplotype WHERE id_ind IN (", parent_id_list,
      ") ORDER BY id_ind, parent_origin, locus_id"
    )
  )

  n_autosome_loci <- length(autosome_locus_idx)

  parent_haps     <- vector("list", length(unique_parents))
  parent_lo       <- vector("list", length(unique_parents))
  parent_chr_haps <- vector("list", length(unique_parents))
  names(parent_haps)     <- unique_parents
  names(parent_lo)       <- unique_parents
  names(parent_chr_haps) <- unique_parents

  for (pid in unique_parents) {
    rows <- parent_haps_raw[parent_haps_raw$id_ind == pid, ]

    # --- Autosome portion: every individual always has exactly 2 copies of
    # every autosome (autosomes are "full/full" by definition of the
    # classification above), regardless of sex. Built COMPACT
    # (2 x n_autosome_loci, special-chromosome loci excluded entirely) using
    # LOCAL column positions (full_to_local), matching autosome_chr_info's
    # local indexing above. ---
    auto_rows <- rows[rows$locus_id %in% genome_meta_df$locus_id[autosome_locus_idx], ]
    if (nrow(auto_rows) != 2L * n_autosome_loci) {
      stop(
        "Parent '", pid, "' has ", nrow(auto_rows), " autosomal ind_haplotype rows; ",
        "expected ", 2L * n_autosome_loci, ".", call. = FALSE
      )
    }
    hap_mat <- matrix(0L, nrow = 2L, ncol = n_autosome_loci)
    lo_mat  <- matrix(NA_character_, nrow = 2L, ncol = n_autosome_loci)
    for (po in c(1L, 2L)) {
      po_rows    <- auto_rows[auto_rows$parent_origin == po, ]
      local_cols <- full_to_local[po_rows$locus_id]
      hap_mat[po, local_cols] <- as.integer(po_rows$allele)
      lo_mat[po, local_cols]  <- po_rows$line_origin
    }
    parent_haps[[pid]] <- hap_mat
    parent_lo[[pid]]   <- lo_mat

    # --- Special-chromosome portions: compact k x n_chr_loci matrices, k
    # resolved from this parent's OWN sex (0, 1, or 2 in this version). ---
    chr_mats <- list()
    for (cn in special_chr_names) {
      chr_locus_ids <- special_locus_ids[[cn]]
      n_chr_loci    <- length(chr_locus_ids)
      cr <- rows[rows$locus_id %in% chr_locus_ids, ]
      cr <- cr[order(cr$parent_origin, cr$locus_id), ]

      cm_own <- if (parent_sex[[pid]] == "M") chr_meta_map[[cn]]$copy_mode_M else chr_meta_map[[cn]]$copy_mode_F
      k <- resolve_chr_copy_count(cm_own, ploidy = parent_ploidy[[pid]])

      if (nrow(cr) != k * n_chr_loci) {
        stop(
          "Parent '", pid, "' has ", nrow(cr), " ind_haplotype rows for chromosome '",
          cn, "'; expected ", k * n_chr_loci, " given chr_meta copy_mode rules for sex '",
          parent_sex[[pid]], "'.", call. = FALSE
        )
      }
      if (k == 0L) {
        chr_mats[[cn]] <- list(hap = matrix(integer(0), nrow = 0, ncol = n_chr_loci),
                               lo  = matrix(character(0), nrow = 0, ncol = n_chr_loci))
      } else {
        chr_mats[[cn]] <- list(
          hap = matrix(as.integer(cr$allele), nrow = k, ncol = n_chr_loci, byrow = TRUE),
          lo  = matrix(cr$line_origin,        nrow = k, ncol = n_chr_loci, byrow = TRUE)
        )
      }
    }
    parent_chr_haps[[pid]] <- chr_mats
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

  # Autosome path: BYTE-IDENTICAL to the pre-Stage-4 code — same single
  # make_gamete() call per parent, same chr_info-derived descriptor list. When
  # every chromosome is a plain autosome (chr_meta all-default),
  # autosome_chr_info == chr_info and special_chr_names is empty, so this is
  # the ONLY code executed per offspring, exactly as before.
  hap_sire_mat    <- matrix(0L, nrow = n_offspring, ncol = n_autosome_loci)
  hap_dam_mat     <- matrix(0L, nrow = n_offspring, ncol = n_autosome_loci)
  hap_sire_lo_mat <- matrix(NA_character_, nrow = n_offspring, ncol = n_autosome_loci)
  hap_dam_lo_mat  <- matrix(NA_character_, nrow = n_offspring, ncol = n_autosome_loci)
  col_idx <- seq_len(n_autosome_loci)

  # Long-format accumulator for special-chromosome (sex-linked/organelle) rows
  # — one entry per (offspring, parent_origin) written. Empty when there are
  # no special chromosomes.
  special_rows <- vector("list", 0L)

  for (i in seq_len(n_offspring)) {
    sire    <- matings$id_parent_1[i]
    dam     <- matings$id_parent_2[i]
    off_sex <- matings$sex[i]

    # make_gamete() returns allele + the per-locus homolog (1/2) it came from.
    # line_origin follows whichever parental homolog donated each locus, so it
    # is correct across F1/F2/backcross generations.
    gs <- make_gamete(parent_haps[[sire]], autosome_chr_info)
    hap_sire_mat[i, ]    <- gs$allele
    hap_sire_lo_mat[i, ] <- parent_lo[[sire]][cbind(gs$homolog, col_idx)]

    gd <- make_gamete(parent_haps[[dam]], autosome_chr_info)
    hap_dam_mat[i, ]    <- gd$allele
    hap_dam_lo_mat[i, ] <- parent_lo[[dam]][cbind(gd$homolog, col_idx)]

    # --- Special-chromosome path, executed strictly AFTER the autosome path
    # above for this offspring (never reordered ahead of it), so any RNG draws
    # here cannot perturb draws already made for the autosome gametes. ---
    for (cn in special_chr_names) {
      chr_row <- chr_meta_map[[cn]]
      cm_off  <- if (off_sex == "M") chr_row$copy_mode_M else chr_row$copy_mode_F
      if (cm_off == "none") next

      contributors <- if (cm_off == "half") {
        side <- chr_row$hemi_parent
        list(list(pid = if (side == "parent_1") sire else dam,
                  parent_origin = if (side == "parent_1") 1L else 2L))
      } else {
        list(list(pid = sire, parent_origin = 1L),
             list(pid = dam,  parent_origin = 2L))
      }

      for (contrib in contributors) {
        mats <- parent_chr_haps[[contrib$pid]][[cn]]
        k_p  <- nrow(mats$hap)
        if (k_p == 0L) {
          stop(
            "chr_meta misconfiguration: parent '", contrib$pid,
            "' has zero copies of chromosome '", cn,
            "' to pass on as its designated contribution.", call. = FALSE
          )
        }
        if (k_p == 2L && isTRUE(chr_row$recombines)) {
          g <- make_gamete(mats$hap, special_chr_local_info[[cn]])
          allele <- g$allele
          lo     <- mats$lo[cbind(g$homolog, seq_along(g$homolog))]
        } else {
          g <- pass_through_gamete(mats$hap, mats$lo)
          allele <- g$allele
          lo     <- g$line_origin
        }
        special_rows[[length(special_rows) + 1L]] <- data.frame(
          id_ind        = offspring_ids[i],
          parent_origin = contrib$parent_origin,
          locus_id      = special_locus_ids[[cn]],
          allele        = as.integer(allele),
          line_origin   = lo,
          stringsAsFactors = FALSE
        )
      }
    }
  }

  # ============================================================================
  # 12. Build the wide allele frame (transient; the vehicle for the UNPIVOT into
  #     the long ind_haplotype). Never written to a wide DB table. Restricted
  #     to autosome columns — special-chromosome loci are written separately
  #     from `special_rows` below.
  # ============================================================================

  autosome_locus_cols <- locus_cols[autosome_locus_idx]

  hap_wide <- dplyr::bind_rows(
    dplyr::bind_cols(
      tibble::tibble(id_ind = offspring_ids, parent_origin = 1L),
      setNames(as.data.frame(hap_sire_mat), autosome_locus_cols)
    ),
    dplyr::bind_cols(
      tibble::tibble(id_ind = offspring_ids, parent_origin = 2L),
      setNames(as.data.frame(hap_dam_mat), autosome_locus_cols)
    )
  )

  # ============================================================================
  # 14. Build ind_meta rows and handle extra columns
  # ============================================================================

  # Offspring ploidy is computed, not hardcoded, as the sum of each parent's
  # gamete contribution (own_ploidy / 2 per parent) — always 2L in this
  # version (both parents are always 2L), but written this way so this code
  # doesn't need revisiting when real polyploidy ships later.
  offspring_ploidy <- parent_ploidy[matings$id_parent_1] %/% 2L +
    parent_ploidy[matings$id_parent_2] %/% 2L
  stopifnot(all(offspring_ploidy == 2L))

  ind_meta_new <- tibble::tibble(
    id_ind      = offspring_ids,
    id_parent_1 = matings$id_parent_1,
    id_parent_2 = matings$id_parent_2,
    line_name   = matings$line_name,
    sex         = matings$sex,
    ploidy      = as.integer(unname(offspring_ploidy))
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

  # Write the long ind_haplotype (autosomes). A wide allele frame and a wide
  # line_origin frame are unpivoted and joined in DuckDB (per (id_ind,
  # parent_origin, locus)) so R never materializes the long frame. ind_genotype
  # stays on-demand (add_dosage()).
  lo_df <- dplyr::bind_rows(
    dplyr::bind_cols(
      tibble::tibble(id_ind = offspring_ids, parent_origin = 1L),
      setNames(as.data.frame(hap_sire_lo_mat, stringsAsFactors = FALSE), autosome_locus_cols)
    ),
    dplyr::bind_cols(
      tibble::tibble(id_ind = offspring_ids, parent_origin = 2L),
      setNames(as.data.frame(hap_dam_lo_mat, stringsAsFactors = FALSE), autosome_locus_cols)
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

  # Write special-chromosome (sex-linked/organelle) rows, already in long
  # format — no UNPIVOT needed. RNG-neutral (register + INSERT).
  if (length(special_rows) > 0L) {
    special_df <- do.call(rbind, special_rows)
    duckdb::duckdb_register(pop$db_conn, "__tmp_special_hap", special_df)
    DBI::dbExecute(pop$db_conn, paste0(
      "INSERT INTO ind_haplotype (id_ind, parent_origin, strand, line_origin, locus_id, locus_name, allele) ",
      "SELECT s.id_ind, s.parent_origin, 1, s.line_origin, gm.locus_id, gm.locus_name, s.allele ",
      "FROM __tmp_special_hap s JOIN genome_meta gm ON s.locus_id = gm.locus_id"
    ))
    duckdb::duckdb_unregister(pop$db_conn, "__tmp_special_hap")
  }

  # ============================================================================
  # 16. Return
  # ============================================================================

  message("Added ", n_offspring, " offspring")

  invisible(pop)
}
