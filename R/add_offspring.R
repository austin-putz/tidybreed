#' Add offspring via recombination
#'
#' @description
#' Produces offspring by simulating chromosomal crossovers from explicit parent
#' pairs. Each row of `matings` defines exactly one offspring, giving full
#' control over the mating design.
#'
#' @param pop A `tidybreed_pop` object
#' @param seed Optional integer base seed for the per-gamete recombination
#'   streams. Each gamete draws from its own deterministic `dqrng` sub-stream
#'   keyed on its global offspring index, so output is independent of
#'   `batch_size` and (later) thread count. `seed = NULL` (default) draws one
#'   base seed from base-R's RNG, so an upstream [set.seed()] still fully
#'   reproduces the call; pass an integer to fix the gamete streams independently
#'   of surrounding base-R state. The resolved base seed is reported in the
#'   message and attached to the returned object as `attr(pop, "base_seed")`; it
#'   is not written to any table.
#' @param store_crossovers Logical (default `FALSE`). When `TRUE`, every
#'   crossover drawn during meiosis is recorded as one row in `ind_crossover`
#'   (`id_ind`, `parent_origin`, `chr`, `chr_name`, `pos_cM`). Absence of a row
#'   for a `(id_ind, parent_origin, chr)` means that gamete's chromosome did not
#'   recombine. Includes crossovers on recombining special chromosomes. `FALSE`
#'   costs nothing (no buffer emitted, no `ind_crossover` write).
#' @param batch_size Optional integer. Offspring generated + written per batch.
#'   Bounds peak memory at roughly `batch_size x n_loci` long rows regardless of
#'   the total number of offspring (a single huge mating streams to disk batch by
#'   batch). Any value — including `batch_size = nrow(matings)` — produces
#'   byte-identical output for a fixed seed; it only trades peak memory against
#'   per-batch overhead. Overrides `max_batch_mem` when both are set.
#' @param max_batch_mem Optional per-batch memory budget as bytes (e.g. `512e6`)
#'   or a string (e.g. `"512MB"`, `"2GB"`). Used to derive `batch_size` when
#'   `batch_size` is `NULL`. When both are `NULL` (the default), the batch size
#'   is auto-picked from detected available system memory (falling back to a
#'   conservative fixed budget if detection is unavailable).
#' @param matings A tibble or data.frame where each row = one offspring.
#'
#'   **Required columns** (canonical names match `ind_meta`):
#'   - `id_parent_1` — sire / parent 1 ID; must exist in `ind_meta`.
#'     Alias: `id_sire`
#'   - `id_parent_2` — dam / parent 2 ID; must exist in `ind_meta`.
#'     Alias: `id_dam`
#'   - `sex` — `"M"` or `"F"`
#'   - `line_name` — line name for offspring IDs (same format as `add_founders()`)
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
#' - Cross-line matings: use any `line_name` value in the `line_name` column
#'
#' **Column aliases:**
#' `id_sire` is accepted as an alias for `id_parent_1`, and `id_dam` for
#' `id_parent_2`. Both naming styles produce identical results.
#'
#' **Recombination model:**
#' Gametes are simulated using the Haldane map function. Each gamete draws from
#' its own deterministic `dqrng` sub-stream keyed on `(base seed, offspring
#' index, parent role)`, so output is identical for any `batch_size` (and, in
#' later versions, any thread count) and is set by `seed` (see above). The
#' crossover count per chromosome ~ Poisson(chr_len_cM\[i\] / 100), i.e.
#' Poisson(genetic length in Morgans) since 1 Morgan = 100 cM, drawn by a
#' uniform-only inversion sampler. Genetic positions come from the resolved
#' genetic map (`genome_map`, per the gamete-producing parent's sex/line);
#' crossover positions are uniform in cM within each chromosome, and the starting
#' haplotype is chosen at random. This
#' applies to plain autosomes (`chr_meta`'s default `copy_mode_M`/`copy_mode_F
#' = "full"`, `recombines_M`/`recombines_F = TRUE`). Sex chromosomes and
#' organelles configured via [define_chr()] follow their own inheritance rule
#' instead: the contributing parent's copy recombines via the same Haldane model
#' only when that parent carries two copies of the chromosome AND the
#' contributor's own-sex `chr_meta.recombines_M`/`recombines_F` is `TRUE`;
#' otherwise the parent's single stored copy is passed straight through unchanged
#' (no recombination, no additional RNG draws).
#'
#' **Offspring IDs:**
#' IDs follow the same `"{line_name}_{n}"` format as `add_founders()`. Numbering
#' continues from the current maximum for each line.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' pop <- open_pop(pop_name = "test", db_name = ":memory:") |>
#'   define_genome(n_loci = 1000, n_chr = 5, chr_len_Mb = 100) |>
#'   define_founder_haplotypes(n_haplotypes = 200, line_name = "A")
#' pop <- pop |>
#'   get_table("founder_haplotypes") |>
#'   add_founders(n_males = 5, n_females = 5, line_name = "A", gen = 1L)
#'
#' # One offspring per mating, extra metadata column (gen)
#' matings <- tibble::tibble(
#'   id_parent_1 = rep("A_1", 5),
#'   id_parent_2 = paste0("A_", 6:10),
#'   sex         = c("M", "F", "M", "F", "M"),
#'   line_name   = "A",
#'   gen         = 2L
#' )
#' pop <- pop |> add_offspring(matings)
#'
#' # Animal-breeder-style aliases
#' matings2 <- tibble::tibble(
#'   id_sire   = rep("A_1", 3),
#'   id_dam    = paste0("A_", 6:8),
#'   sex       = c("M", "F", "M"),
#'   line_name = "A",
#'   gen       = 2L
#' )
#' pop <- pop |> add_offspring(matings2)
#' }
add_offspring <- function(pop, matings, seed = NULL, store_crossovers = FALSE,
                          batch_size = NULL, max_batch_mem = NULL) {

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
      "Ensure the population was created with open_pop() |> define_genome().",
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
    paste0("SELECT id_ind, sex, ploidy, line_name FROM ind_meta WHERE id_ind IN (", parent_id_list, ")")
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
  # Parent line — carried for the sex/line map-resolution seam (R9). v1 has only
  # the default map, so every gamete resolves to it regardless of parent line.
  parent_line   <- stats::setNames(parent_meta$line_name, parent_meta$id_ind)

  # ============================================================================
  # 8. Load genome metadata (once); classify chromosomes as "plain autosome"
  #    (both sexes full copy, recombines in both sexes — chr_meta's default row) vs.
  #    "special" (any sex-linked/organelle rule). Autosomes are handled by the
  #    unchanged fast path below; special chromosomes get their own branch,
  #    executed strictly AFTER the autosome path for each offspring so autosome
  #    RNG draws are never perturbed by chr_meta configuration.
  # ============================================================================

  # Recombination is driven by pos_cM (genetic distance), resolved PER GAMETE for
  # the gamete-producing parent's (sex, line) via resolve_genome_map()'s decision-#6
  # precedence. The chromosome CLASSIFICATION (which chromosomes are plain autosomes
  # vs special, and the compact autosome-index remap) is position-independent, so it
  # is computed once from a base resolution below; only the per-locus pos_cM/chr_len
  # differ per (sex, line). In v1 only the default map exists, so every (sex, line)
  # resolves to the same map and output is unchanged; adding sex/line-specific
  # genome_map rows changes recombination for the matching gametes with no further
  # code change.
  genome_meta_df <- resolve_genome_map(pop$db_conn)   # base/default map

  n_loci     <- nrow(genome_meta_df)

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
  # chr_info locus_idx values reference the FULL genome_meta_df ordering; the
  # autosome fast path uses COMPACT parent matrices (special loci excluded), so
  # autosome descriptors are re-expressed with LOCAL indices via full_to_local.
  full_to_local <- integer(n_loci)
  full_to_local[autosome_locus_idx] <- seq_along(autosome_locus_idx)
  special_locus_ids <- stats::setNames(
    lapply(special_chr_names, function(cn) genome_meta_df$locus_id[genome_meta_df$chr_name == cn]),
    special_chr_names
  )

  # Turn a resolved map into the gamete descriptors: autosome fast-path info +
  # per-special-chromosome LOCAL info. Locus structure is shared across maps;
  # positions (pos_cM/chr_len) come from the resolved map.
  build_gamete_info <- function(resolved_map) {
    ci <- build_chr_info(resolved_map)
    autosome <- lapply(ci[chr_id_keys[is_autosome_id]], function(x) {
      list(locus_idx = full_to_local[x$locus_idx], pos_cM = x$pos_cM, chr_len = x$chr_len)
    })
    special <- stats::setNames(
      lapply(special_chr_names, function(cn) {
        orig <- ci[[chr_id_keys[chr_name_by_id == cn]]]
        list(list(locus_idx = seq_along(orig$locus_idx), pos_cM = orig$pos_cM,
                  chr_len = orig$chr_len))
      }),
      special_chr_names
    )
    list(autosome = autosome, special = special)
  }

  # Resolve + cache one descriptor set per unique (parent sex, parent line) among the
  # parents, so a gamete picks its producer's map in O(1) inside the loop.
  gamete_info_key <- stats::setNames(
    vapply(unique_parents, function(pid) {
      paste0(parent_sex[[pid]], "\r",
             if (is.na(parent_line[[pid]])) "" else parent_line[[pid]])
    }, character(1)),
    unique_parents
  )
  gamete_info_cache <- list()
  for (pid in unique_parents) {
    key <- gamete_info_key[[pid]]
    if (is.null(gamete_info_cache[[key]])) {
      ln <- parent_line[[pid]]
      rm_ <- resolve_genome_map(pop$db_conn, sex = parent_sex[[pid]],
                                line_name = if (is.na(ln)) NULL else ln)
      gamete_info_cache[[key]] <- build_gamete_info(rm_)
    }
  }

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

  # One O(N) split of ROW INDICES replaces the old O(parents x rows) per-parent
  # scan (parent_haps_raw[parent_haps_raw$id_ind == pid, ] inside the loop).
  # Splitting indices (not the data.frame) keeps peak memory low — only one
  # parent's subset is materialized at a time, as before — while removing the
  # repeated full-table scan. split() preserves each group's original row order
  # (the query's ORDER BY id_ind, parent_origin, locus_id), so the per-parent
  # frames are identical — RNG- and output-neutral.
  parent_row_idx <- split(seq_len(nrow(parent_haps_raw)), parent_haps_raw$id_ind)

  for (pid in unique_parents) {
    rows <- parent_haps_raw[parent_row_idx[[pid]], ]

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
  # 11. Batch setup (parent-keyed descriptors + offspring metadata; once)
  # ============================================================================

  # ---- Packed, C++-ready kernel inputs (Stage 3a) ---------------------------
  # The gamete seam (make_gametes_batch_r / the C++ kernel) is a numeric function
  # over packed integer arrays — no string-keyed lists, no character line_origin.
  # Build them once here from the per-parent structures assembled in §9.
  n_par   <- length(unique_parents)
  par_row <- stats::setNames(seq_len(n_par), unique_parents)  # packed row per pid

  parent_allele <- matrix(0L, nrow = 2L * n_par, ncol = n_autosome_loci)
  parent_lo_str <- matrix(NA_character_, nrow = 2L * n_par, ncol = n_autosome_loci)
  for (pid in unique_parents) {
    p <- par_row[[pid]]
    parent_allele[c(2L * p - 1L, 2L * p), ] <- parent_haps[[pid]]
    parent_lo_str[c(2L * p - 1L, 2L * p), ] <- parent_lo[[pid]]
  }
  # line_origin dictionary (integer codes, first-seen order); decoded back to
  # VARCHAR at write time. Bijective, so the stored line_origin is unchanged.
  lo_levels      <- unique(as.vector(parent_lo_str))
  lo_levels      <- lo_levels[!is.na(lo_levels)]
  parent_lo_code <- matrix(match(parent_lo_str, lo_levels),
                           nrow = 2L * n_par, ncol = n_autosome_loci)

  # Packed-parent index per offspring (global; sliced per batch below).
  parent_1_idx_all <- unname(par_row[matings$id_parent_1])
  parent_2_idx_all <- unname(par_row[matings$id_parent_2])

  # Per-gamete resolved-map key (group-by-map seam). A gamete's map is its
  # PRODUCING parent's (sex, line) map, so a sire and a dam — or two lines — can
  # use different maps (sex-specific recombination, achiasmy, line maps). Each
  # distinct key's autosome descriptor is flattened to the kernel's flat arrays
  # once here; the batch loop groups gametes by key and calls the kernel per
  # group. In v1 with only the default map there is one key and one call.
  sire_mapkey_all <- unname(gamete_info_key[matings$id_parent_1])
  dam_mapkey_all  <- unname(gamete_info_key[matings$id_parent_2])

  chr_arrays_by_key <- stats::setNames(
    lapply(names(gamete_info_cache), function(key) {
      auto <- gamete_info_cache[[key]]$autosome
      pos  <- numeric(n_autosome_loci)
      for (ci in auto) pos[ci$locus_idx] <- ci$pos_cM
      list(
        chr_start  = vapply(auto, function(ci) min(ci$locus_idx), integer(1)),
        chr_end    = vapply(auto, function(ci) max(ci$locus_idx), integer(1)),
        chr_len_cM = vapply(auto, function(ci) ci$chr_len, numeric(1)),
        chr_pos_cM = pos
      )
    }),
    names(gamete_info_cache)
  )

  # Autosome locus_ids in LOCAL (compact) column order — column k of the kernel's
  # locus_idx (1..n_autosome_loci) maps to autosome_locus_ids[k].
  autosome_locus_ids <- genome_meta_df$locus_id[autosome_locus_idx]

  # chr number/name lookups for ind_crossover. The autosome kernel returns
  # `xover_chr` (position in the ascending autosome chr order); map it here.
  autosome_chr_nums  <- chr_ids_all[is_autosome_id]
  autosome_chr_names <- unname(chr_name_by_id[chr_id_keys[is_autosome_id]])
  # chr_name -> chr number, for special-chromosome crossover rows.
  chr_num_by_name    <- stats::setNames(chr_ids_all, unname(chr_name_by_id[chr_id_keys]))

  # Offspring ploidy is computed, not hardcoded, as the sum of each parent's
  # gamete contribution (own_ploidy / 2 per parent) — always 2L in this version
  # (both parents are always 2L), but written this way so this code doesn't need
  # revisiting when real polyploidy ships later. (RNG-neutral; no draws.)
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

  # ============================================================================
  # 12. Seed resolution + batched long write, all inside ONE transaction.
  #
  # RNG model (Stage 2): every gamete draws from its own dqrng sub-stream keyed
  # on its GLOBAL offspring index `o` (the original `matings` row), so output is
  # independent of `batch_size` and iteration order. `seed = NULL` draws one base
  # seed from base-R RNG, so an upstream set.seed() reproduces the whole run;
  # `seed = <int>` uses it directly. The resolved base seed is surfaced (message
  # + invisible attribute), never persisted (CLAUDE.md's no-metadata stance).
  # ============================================================================

  if (!is.null(seed)) {
    # Guard finiteness/range BEFORE as.integer(): a non-finite or out-of-int32
    # seed (e.g. Inf, 3e9) makes as.integer() return NA, and `NA != NA` would
    # yield `if (NA)` — a cryptic "missing value" error instead of this message.
    if (length(seed) != 1L || is.na(seed) || !is.finite(seed) ||
        seed != trunc(seed) || seed > .Machine$integer.max ||
        seed < -.Machine$integer.max) {
      stop("seed must be NULL or a single integer within +/- .Machine$integer.max.",
           call. = FALSE)
    }
    base_seed <- as.integer(seed)
  } else {
    base_seed <- sample.int(.Machine$integer.max, 1L)
  }
  if (length(store_crossovers) != 1L || is.na(store_crossovers) ||
      !is.logical(store_crossovers)) {
    stop("store_crossovers must be a single TRUE/FALSE.", call. = FALSE)
  }

  B    <- resolve_batch_size(n_offspring, n_loci, batch_size, max_batch_mem)
  conn <- pop$db_conn

  DBI::dbExecute(conn, "BEGIN TRANSACTION")
  committed <- FALSE
  on.exit(if (!committed) try(DBI::dbExecute(conn, "ROLLBACK"), silent = TRUE),
          add = TRUE)

  for (b_start in seq.int(1L, n_offspring, by = B)) {
    b_idx  <- b_start:min(b_start + B - 1L, n_offspring)
    b_n    <- length(b_idx)
    b_ids  <- offspring_ids[b_idx]
    b_sire <- matings$id_parent_1[b_idx]
    b_dam  <- matings$id_parent_2[b_idx]
    b_sex  <- matings$sex[b_idx]

    # --- Autosomes: gamete-flat kernel, grouped by resolved-map key (kind = 1
    # streams). Each gamete owns an independent stream keyed on its global
    # o = b_idx[i], so output is identical for any B and any map partitioning,
    # with or without special chromosomes. In v1 there is a single map group.
    gam_o      <- rep(b_idx, each = 2L)
    gam_origin <- rep(c(1L, 2L), times = b_n)
    gam_parent <- integer(2L * b_n)
    gam_parent[c(TRUE, FALSE)] <- parent_1_idx_all[b_idx]
    gam_parent[c(FALSE, TRUE)] <- parent_2_idx_all[b_idx]
    gam_key <- character(2L * b_n)
    gam_key[c(TRUE, FALSE)] <- sire_mapkey_all[b_idx]
    gam_key[c(FALSE, TRUE)] <- dam_mapkey_all[b_idx]

    auto_parts  <- list()
    xover_parts <- list()   # batch crossover accumulator (autosomes + special)
    if (n_autosome_loci > 0L) {
      for (key in unique(gam_key)) {
        sel <- which(gam_key == key)
        ca  <- chr_arrays_by_key[[key]]
        gg  <- make_gametes_batch_r(
          parent_allele, parent_lo_code,
          gam_parent[sel], gam_o[sel], gam_origin[sel],
          ca$chr_start, ca$chr_end, ca$chr_pos_cM, ca$chr_len_cM,
          base_seed, store_crossovers)

        # Rows are in this group's gamete order; each gamete is one n_loci block.
        row_o <- rep(gam_o[sel], each = n_autosome_loci)
        auto_parts[[length(auto_parts) + 1L]] <- data.frame(
          id_ind        = b_ids[match(row_o, b_idx)],
          parent_origin = gg$parent_origin,
          strand        = 1L,
          line_origin   = lo_levels[gg$line_origin_code],
          locus_id      = autosome_locus_ids[gg$locus_idx],
          allele        = gg$allele,
          stringsAsFactors = FALSE
        )
        if (store_crossovers && length(gg$xover_pos_cM)) {
          xover_parts[[length(xover_parts) + 1L]] <- data.frame(
            id_ind        = b_ids[match(gg$xover_gamete_o, b_idx)],
            parent_origin = gg$xover_parent_origin,
            chr           = autosome_chr_nums[gg$xover_chr],
            chr_name      = autosome_chr_names[gg$xover_chr],
            pos_cM        = gg$xover_pos_cM,
            stringsAsFactors = FALSE
          )
        }
      }
    }
    auto_long <- if (length(auto_parts) > 0L) do.call(rbind, auto_parts) else NULL

    # --- Special chromosomes (sex chromosomes / organelles): the R path, on the
    # independent kind = 2 stream. Seed ONCE per (offspring o, parent_origin r)
    # and draw all of that (gamete, r)'s special chromosomes in ascending-chr
    # order from the one stream — never re-seed per chromosome (that would
    # collide). Order-independent from autosomes: different stream kind. ---
    batch_special <- list()
    if (length(special_chr_names) > 0L) {
      sr <- vector("list", b_n * length(special_chr_names))
      sr_count <- 0L

      for (j in seq_len(b_n)) {
        o       <- b_idx[j]
        sire    <- b_sire[j]
        dam     <- b_dam[j]
        off_sex <- b_sex[j]

        for (r in c(1L, 2L)) {
          # Which special chromosomes does parent_origin r contribute to this
          # offspring (given its sex)? Ascending chr order (special_chr_names is
          # already chr-ascending).
          r_chrs <- character(0)
          for (cn in special_chr_names) {
            chr_row <- chr_meta_map[[cn]]
            cm_off  <- if (off_sex == "M") chr_row$copy_mode_M else chr_row$copy_mode_F
            if (cm_off == "none") next
            if (cm_off == "full") {
              r_chrs <- c(r_chrs, cn)                      # both r contribute
            } else {                                       # "half": hemi side only
              hemi_po <- if (chr_row$hemi_parent == "parent_1") 1L else 2L
              if (hemi_po == r) r_chrs <- c(r_chrs, cn)
            }
          }
          if (length(r_chrs) == 0L) next

          pid <- if (r == 1L) sire else dam
          .seed_gamete_stream(base_seed, .gamete_stream_id(o, r, 2L))

          for (cn in r_chrs) {
            chr_row <- chr_meta_map[[cn]]
            mats    <- parent_chr_haps[[pid]][[cn]]
            k_p     <- nrow(mats$hap)
            if (k_p == 0L) {
              stop(
                "chr_meta misconfiguration: parent '", pid,
                "' has zero copies of chromosome '", cn,
                "' to pass on as its designated contribution.", call. = FALSE
              )
            }
            contrib_recombines <- if (parent_sex[[pid]] == "M") {
              chr_row$recombines_M
            } else {
              chr_row$recombines_F
            }

            xpos <- numeric(0)
            if (k_p == 2L && isTRUE(contrib_recombines)) {
              si <- gamete_info_cache[[gamete_info_key[[pid]]]]$special[[cn]][[1L]]
              d  <- .draw_chr_recombination(si$pos_cM, si$chr_len, store_crossovers)
              col_idx <- seq_along(d$hap_idx_vec)
              allele  <- mats$hap[cbind(d$hap_idx_vec, col_idx)]
              lo      <- mats$lo[cbind(d$hap_idx_vec, col_idx)]
              xpos    <- d$xover_pos_cM
            } else {
              gg     <- pass_through_gamete(mats$hap, mats$lo)
              allele <- gg$allele
              lo     <- gg$line_origin
            }

            sr_count <- sr_count + 1L
            sr[[sr_count]] <- data.frame(
              id_ind        = b_ids[j],
              parent_origin = r,
              strand        = 1L,
              line_origin   = lo,
              locus_id      = special_locus_ids[[cn]],
              allele        = as.integer(allele),
              stringsAsFactors = FALSE
            )
            if (store_crossovers && length(xpos)) {
              xover_parts[[length(xover_parts) + 1L]] <- data.frame(
                id_ind        = b_ids[j],
                parent_origin = r,
                chr           = unname(chr_num_by_name[cn]),
                chr_name      = cn,
                pos_cM        = xpos,
                stringsAsFactors = FALSE
              )
            }
          }
        }
      }
      batch_special <- if (sr_count > 0L) sr[seq_len(sr_count)] else list()
    }

    # Assemble the batch's long rows: kernel autosome output (already long) plus
    # the special-chromosome rows. No wide matrices, no UNPIVOT.
    long_frame <- auto_long
    if (length(batch_special) > 0L) {
      special_long <- do.call(rbind, batch_special)
      long_frame <- if (is.null(long_frame)) special_long else
        rbind(long_frame, special_long)
    }

    .write_long_haplotypes(conn, long_frame)

    if (store_crossovers && length(xover_parts) > 0L) {
      .write_crossovers(conn, do.call(rbind, xover_parts))
    }

    rm(auto_parts, auto_long, long_frame, batch_special, xover_parts,
       gam_o, gam_origin, gam_parent, gam_key)
    gc(verbose = FALSE)
  }

  # ind_meta rows: extra user columns validated + type-inferred here (RNG-neutral
  # ALTER TABLE), then the single dbWriteTable append — issued after every gamete
  # draw so its RNG advance matches the pre-Stage-1 stream position exactly.
  if (length(extra_cols) > 0) {
    extra_list <- stats::setNames(
      lapply(extra_cols, function(col) matings[[col]]),
      extra_cols
    )
    prepped <- prepare_extra_cols(extra_list, n_offspring, "ind_meta", conn)
    for (nm in names(prepped)) ind_meta_new[[nm]] <- prepped[[nm]]
  }
  DBI::dbWriteTable(conn, "ind_meta", ind_meta_new, append = TRUE)

  DBI::dbExecute(conn, "COMMIT")
  committed <- TRUE

  # ============================================================================
  # 13. Return
  # ============================================================================

  message("Added ", n_offspring, " offspring (base_seed = ", base_seed, ")")

  # Surface the resolved base seed for note-taking / exact replay via
  # add_offspring(seed = <that value>). Not persisted to any table.
  attr(pop, "base_seed") <- base_seed
  invisible(pop)
}
