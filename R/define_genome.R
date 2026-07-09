#' Define the genome structure of a breeding population
#'
#' @description
#' Adds genome tables (`genome_meta`, `genome_map`, `ind_haplotype`,
#' `ind_genotype`, `chr_meta`) to a population opened with [open_pop()].
#' Pipe-friendly — accepts a `tidybreed_pop` and returns a `tidybreed_pop`.
#' Physical position (`pos_bp`, base pairs) lives in `genome_meta`; the genetic
#' map (`pos_cM`, centiMorgans) lives in the long `genome_map` table. Haplotypes
#' are stored in long format; `ind_genotype` is an on-demand dosage cache
#' populated by [add_dosage()].
#'
#' This is the second step in setting up a simulation:
#' ```r
#' pop <- open_pop(...) |>
#'   define_genome(n_loci = 50000, n_chr = 18, chr_len_Mb = 100)
#' ```
#'
#' To create the founder haplotype pool needed by [add_founders()], call
#' [define_founder_haplotypes()] after this function.
#'
#' @param pop A `tidybreed_pop` object from [open_pop()].
#' @param n_loci Integer scalar. Total number of loci to simulate.
#' @param n_chr Integer scalar. Number of chromosomes.
#' @param chr_len_Mb Numeric scalar or numeric vector of length `n_chr`.
#'   Chromosome length(s) in megabases. A single scalar applies the same
#'   length to all chromosomes; a vector specifies each chromosome separately.
#'   Determines the physical coordinate `genome_meta.pos_bp`.
#' @param cM_per_Mb Numeric scalar or numeric vector of length `n_chr`. Genetic
#'   map rate in centiMorgans per megabase, used to derive the genetic map
#'   position `genome_map.pos_cM` from `pos_bp` (`pos_cM = pos_bp/1e6 *
#'   cM_per_Mb`). Default `1.0`. Must be finite and strictly positive.
#' @param locus_names Character vector of length `n_loci` or `NULL`. Custom
#'   locus names. When `NULL` (default), names are auto-generated as
#'   `"Locus_1"`, `"Locus_2"`, etc.
#' @param chr_names Character vector of length `n_chr` or `NULL`. Custom
#'   chromosome names. When `NULL` (default), chromosomes are numbered
#'   `1, 2, ..., n_chr`.
#'
#' @return The input `tidybreed_pop` with genome tables added.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Simple setup
#' pop <- open_pop(pop_name = "A", db_name = ":memory:") |>
#'   define_genome(n_loci = 1000, n_chr = 10, chr_len_Mb = 100)
#'
#' # Different chromosome lengths (cattle)
#' pop <- open_pop(pop_name = "Cattle", db_name = ":memory:") |>
#'   define_genome(
#'     n_loci = 50000,
#'     n_chr  = 30,
#'     chr_len_Mb = c(158, 137, 121, 120, 121, 119, 112, 113, 105, 104,
#'                    107,  91,  84,  84,  85,  81,  75,  66,  64,  72,
#'                     71,  61,  52,  62,  42,  51,  45,  46,  51,  42)
#'   )
#'
#' # With founder haplotypes
#' pop <- open_pop(pop_name = "B", db_name = ":memory:") |>
#'   define_genome(n_loci = 1000, n_chr = 10, chr_len_Mb = 100) |>
#'   define_founder_haplotypes(n_haplotypes = 100,
#'                             min_allele_freq = 0.05,
#'                             max_allele_freq = 0.95)
#' }
define_genome <- function(pop,
                          n_loci,
                          n_chr,
                          chr_len_Mb,
                          cM_per_Mb   = 1.0,
                          locus_names = NULL,
                          chr_names   = NULL) {

  stopifnot(inherits(pop, "tidybreed_pop"))
  stopifnot(is.numeric(n_loci), length(n_loci) == 1, n_loci > 0)
  stopifnot(is.numeric(n_chr),  length(n_chr)  == 1, n_chr  > 0)
  if (!is.numeric(chr_len_Mb) || !length(chr_len_Mb) %in% c(1, n_chr) ||
      any(!is.finite(chr_len_Mb)) || any(chr_len_Mb <= 0)) {
    stop("chr_len_Mb must be a finite, strictly positive numeric scalar or ",
         "vector of length n_chr.", call. = FALSE)
  }
  if (!is.numeric(cM_per_Mb) || !length(cM_per_Mb) %in% c(1, n_chr) ||
      any(!is.finite(cM_per_Mb)) || any(cM_per_Mb <= 0)) {
    stop("cM_per_Mb must be a finite, strictly positive numeric scalar or ",
         "vector of length n_chr.", call. = FALSE)
  }

  db_conn <- pop$db_conn

  # Preflight: refuse to re-define a genome. define_genome() replaces genome_meta /
  # genome_map (CREATE OR REPLACE) but then plain-CREATEs ind_haplotype/ind_genotype
  # — a second call would partially clobber the genome before failing on the
  # existing haplotype tables. Fail cleanly up front instead.
  if (DBI::dbExistsTable(db_conn, "genome_meta") &&
      DBI::dbGetQuery(db_conn, "SELECT COUNT(*) AS n FROM genome_meta")$n > 0L) {
    stop("Genome already defined for this population (genome_meta is non-empty). ",
         "Call define_genome() once on a fresh population from open_pop().",
         call. = FALSE)
  }

  # Expand chr_len_Mb / cM_per_Mb if single value
  if (length(chr_len_Mb) == 1)
    chr_len_Mb <- rep(chr_len_Mb, n_chr)
  if (length(cM_per_Mb) == 1)
    cM_per_Mb <- rep(cM_per_Mb, n_chr)

  # Generate locus names if not provided
  if (is.null(locus_names)) {
    locus_names <- paste0("Locus_", seq_len(n_loci))
  } else {
    stopifnot(length(locus_names) == n_loci)
  }

  # locus_name must be unique — it is a denormalized key in the long
  # ind_haplotype/ind_genotype tables, and a duplicate would silently corrupt
  # joins to genome_effects. Validated in R rather than via a DB-level UNIQUE
  # constraint (kept simple; genome_meta carries no DB constraints).
  if (anyDuplicated(locus_names)) {
    dup <- unique(locus_names[duplicated(locus_names)])
    stop("locus_names must be unique. Duplicated: ",
         paste(utils::head(dup, 5L), collapse = ", "),
         if (length(dup) > 5L) ", ..." else "", call. = FALSE)
  }

  # Generate chromosome names if not provided
  if (is.null(chr_names)) {
    chr_names <- as.character(seq_len(n_chr))
  } else {
    stopifnot(length(chr_names) == n_chr)
  }

  # Assign loci to chromosomes (evenly distributed)
  loci_per_chr   <- diff(round(seq(0, n_loci, length.out = n_chr + 1)))
  chr_assignment <- rep(seq_len(n_chr), times = loci_per_chr)

  # Generate physical positions (pos_bp) within each chromosome (evenly spaced),
  # then derive genetic positions (pos_cM) from the sampled pos_bp. pos_bp is the
  # single physical source of truth; pos_cM = pos_bp/1e6 * cM_per_Mb.
  pos_bp <- numeric(n_loci)
  pos_cM <- numeric(n_loci)
  for (i in seq_len(n_chr)) {
    chr_loci   <- which(chr_assignment == i)
    n_chr_loci <- length(chr_loci)

    chr_len_bp_i <- round(chr_len_Mb[i] * 1e6)
    pos_bp_i     <- round(seq(0, chr_len_bp_i,
                              length.out = n_chr_loci + 2)[2:(n_chr_loci + 1)])

    # Validate BEFORE writing: physical positions must be positive, strictly
    # increasing, and duplicate-free. Very short chromosomes or very dense loci
    # make round(seq(...)) collide or hit 0. No silent repair.
    if (any(pos_bp_i <= 0) || (n_chr_loci > 1 && any(diff(pos_bp_i) <= 0))) {
      stop(sprintf(
        "chromosome %s: physical span %.0f bp is too short to place %d loci as positive, strictly-increasing integer positions; increase chr_len_Mb[%d] or reduce loci.",
        chr_names[i], chr_len_bp_i, n_chr_loci, i), call. = FALSE)
    }

    pos_bp[chr_loci] <- pos_bp_i
    pos_cM[chr_loci] <- (pos_bp_i / 1e6) * cM_per_Mb[i]
  }

  # ---- genome_meta (physical only) ---------------------------------------
  # Explicit typed CREATE so pos_bp is BIGINT (bulletproof for any genome size;
  # dbWriteTable would infer DOUBLE/INTEGER). Populate via duckdb_register +
  # INSERT (NOT dbWriteTable, which advances R's RNG via a random temp name and
  # would shift the seeded draw sequence).
  DBI::dbExecute(db_conn, paste0(
    "CREATE OR REPLACE TABLE genome_meta (",
    "locus_id INTEGER, locus_name VARCHAR, chr INTEGER, chr_name VARCHAR, ",
    "pos_bp BIGINT)"
  ))
  genome_meta_df <- tibble::tibble(
    locus_id   = seq_len(n_loci),
    locus_name = locus_names,
    chr        = chr_assignment,
    chr_name   = chr_names[chr_assignment],
    pos_bp     = pos_bp
  )
  duckdb::duckdb_register(db_conn, "__tmp_genome_meta",
                          as.data.frame(genome_meta_df))
  DBI::dbExecute(db_conn, paste0(
    "INSERT INTO genome_meta ",
    "SELECT locus_id, locus_name, chr, chr_name, pos_bp FROM __tmp_genome_meta"
  ))
  duckdb::duckdb_unregister(db_conn, "__tmp_genome_meta")

  # Create empty long ind_haplotype table (one row per individual x haplotype x locus).
  # No PRIMARY KEY: the 4-col ART index dominated bulk-insert cost (~53% of
  # add_founders runtime; benched ~2x faster to drop, and it makes remove_rows'
  # id_ind DELETE ~74x faster via zonemap pruning). Uniqueness of
  # (id_ind, parent_origin, strand, locus_id) is guaranteed R-side: id_ind is
  # MAX(id)+1 from ind_meta (which has its own PK) and each individual emits each
  # (parent_origin, strand, locus_id) tuple exactly once. Contrast ind_genotype,
  # which keeps its PK for INSERT OR REPLACE idempotency.
  DBI::dbExecute(db_conn, paste0(
    "CREATE TABLE ind_haplotype (",
    "id_ind VARCHAR, parent_origin UTINYINT, strand UTINYINT, ",
    "line_origin VARCHAR, locus_id INTEGER, locus_name VARCHAR, allele UTINYINT)"
  ))

  # Create empty long ind_genotype table (on-demand dosage cache; add_dosage())
  DBI::dbExecute(db_conn, paste0(
    "CREATE TABLE ind_genotype (",
    "id_ind VARCHAR, locus_id INTEGER, locus_name VARCHAR, dosage_value UTINYINT, ",
    "PRIMARY KEY (id_ind, locus_id))"
  ))

  # Create empty ind_crossover table (opt-in crossover-event storage). Created
  # here so restore_pop() sees it and it exists before the first insert; row
  # writes land with the Stage-2 kernel (add_offspring(store_crossovers = TRUE)).
  # id_crossover is BIGINT (not INTEGER like other id_ PKs): crossover events are
  # ~1-2 orders of magnitude more numerous than any other id_-keyed row (a
  # ~25-Morgan genome is ~50 events/offspring), so large multi-generation runs
  # would blow past INTEGER's 2.1B ceiling. Deliberate exception; empty now so
  # the wider PK costs nothing.
  DBI::dbExecute(db_conn, paste0(
    "CREATE TABLE ind_crossover (",
    "id_crossover BIGINT PRIMARY KEY, id_ind VARCHAR, parent_origin UTINYINT, ",
    "chr INTEGER, chr_name VARCHAR, pos_cM DOUBLE)"
  ))

  # ---- genome_map (genetic map, long) ------------------------------------
  # One row per (locus x sex x line x map) with a defined genetic position. v1
  # writes a single default map (sex NULL = both sexes, line_name NULL = all
  # lines, map_name = "default"). Explicit typed CREATE; populate RNG-safely.
  DBI::dbExecute(db_conn, paste0(
    "CREATE OR REPLACE TABLE genome_map (",
    "id_genome_map INTEGER, locus_id INTEGER, locus_name VARCHAR, ",
    "sex VARCHAR, line_name VARCHAR, map_name VARCHAR, pos_cM DOUBLE)"
  ))
  start_map_id  <- next_int_id(db_conn, "genome_map", "id_genome_map")
  genome_map_df <- tibble::tibble(
    id_genome_map = seq.int(start_map_id, length.out = n_loci),
    locus_id      = seq_len(n_loci),
    locus_name    = locus_names,
    sex           = NA_character_,
    line_name     = NA_character_,
    map_name      = "default",
    pos_cM        = pos_cM
  )
  duckdb::duckdb_register(db_conn, "__tmp_genome_map",
                          as.data.frame(genome_map_df))
  DBI::dbExecute(db_conn, paste0(
    "INSERT INTO genome_map ",
    "SELECT id_genome_map, locus_id, locus_name, sex, line_name, map_name, pos_cM ",
    "FROM __tmp_genome_map"
  ))
  duckdb::duckdb_unregister(db_conn, "__tmp_genome_map")
  validate_genome_map(db_conn)

  # Create chr_meta with default diploid-autosome rows. copy_mode_M/copy_mode_F
  # are relative to each individual's own ploidy ("full"/"half"/"none"), not an
  # absolute copy count — see define_chr() for non-default inheritance rules
  # (sex chromosomes, organelles). recombines_M/recombines_F are per-sex (matching
  # copy_mode_M/_F); FALSE gives whole-chromosome achiasmy for that sex.
  DBI::dbExecute(db_conn, paste0(
    "CREATE TABLE chr_meta (",
    "chr_name VARCHAR PRIMARY KEY, ",
    "copy_mode_M VARCHAR NOT NULL DEFAULT 'full', ",
    "copy_mode_F VARCHAR NOT NULL DEFAULT 'full', ",
    "hemi_parent VARCHAR, ",
    "recombines_M BOOLEAN NOT NULL DEFAULT TRUE, ",
    "recombines_F BOOLEAN NOT NULL DEFAULT TRUE)"
  ))
  chr_meta_df <- tibble::tibble(
    chr_name     = chr_names,
    copy_mode_M  = "full",
    copy_mode_F  = "full",
    hemi_parent  = NA_character_,
    recombines_M = TRUE,
    recombines_F = TRUE
  )
  # Use duckdb_register + INSERT (not dbAppendTable/dbWriteTable, which advance
  # R's RNG via an internal random temp name and would shift the simulation's
  # seeded draw sequence).
  duckdb::duckdb_register(db_conn, "__tmp_chr_meta", as.data.frame(chr_meta_df))
  DBI::dbExecute(db_conn, "INSERT INTO chr_meta SELECT * FROM __tmp_chr_meta")
  duckdb::duckdb_unregister(db_conn, "__tmp_chr_meta")

  # Register genome-table schema descriptions
  register_schema_meta(db_conn, .genome_table_descriptions())

  # Update pop$tables
  pop$tables <- c(pop$tables, "genome_meta", "genome_map",
                  "ind_haplotype", "ind_genotype", "ind_crossover", "chr_meta")

  chr_len_str <- if (length(unique(chr_len_Mb)) == 1) {
    paste0("all equal to ", chr_len_Mb[1], " Mb")
  } else {
    paste(paste0(chr_names, ":", chr_len_Mb), collapse = ", ")
  }
  message(
    "Defined genome: ", n_loci, " loci across ", n_chr, " chromosomes",
    " | chr lengths (Mb): ", chr_len_str,
    "\n  Tables written: genome_meta, genome_map, ind_haplotype, ind_genotype, ind_crossover, chr_meta"
  )

  pop
}
