#' Initialize a breeding population with genome
#'
#' @description
#' Creates a new `tidybreed_pop` object with a DuckDB backend and initializes
#' the genome tables. This is typically the first function called when setting
#' up a breeding program simulation.
#'
#' The function creates three core genome tables:
#' - `genome_meta`: Locus metadata (position, chromosome, etc.)
#' - `genome_haplotype`: Phased haplotypes for each individual (2 rows per individual)
#' - `genome_genotype`: Genotypes coded as 0/1/2 for each individual
#'
#' Optionally creates founder haplotypes for later sampling:
#' - `founder_haplotypes`: Pool of haplotypes to sample from when adding founders
#'
#' @param pop_name Character string naming the population
#' @param n_loci Integer. Total number of loci to simulate
#' @param n_chr Integer. Number of chromosomes
#' @param chr_len_Mb Numeric. Chromosome length in megabases. Can be:
#'   - Single value (all chromosomes same length)
#'   - Vector of length `n_chr` (chromosome-specific lengths)
#' @param db_path Character. Path to DuckDB database file. Default creates
#'   a file in the current directory named `{pop_name}_tidybreed.duckdb`.
#'   Use `:memory:` for in-memory database (not recommended for large sims).
#' @param locus_names Character vector of locus names. If NULL (default),
#'   names are generated as "Locus_1", "Locus_2", etc.
#' @param chr_names Character vector of chromosome names. If NULL (default),
#'   chromosomes are numbered 1, 2, ..., n_chr.
#' @param overwrite Logical. If TRUE and `db_path` exists, overwrite it.
#'   Default is FALSE.
#' @param n_haplotypes Integer. Number of founder haplotypes to generate.
#'   If NULL (default), no founder haplotypes are created.
#' @param allele_freq_dist Character. Distribution for sampling allele frequencies.
#'   Currently only "uniform" is supported. Ignored if `fixed_allele_freq` is provided.
#' @param min_allele_freq Numeric. Minimum allele frequency for uniform distribution.
#'   Default is 0.01 (1%).
#' @param max_allele_freq Numeric. Maximum allele frequency for uniform distribution.
#'   Default is 0.99 (99%).
#' @param fixed_allele_freq Numeric. If provided, use this fixed allele frequency
#'   for all loci (e.g., 0.5 for 50% frequency). Overrides `allele_freq_dist`,
#'   `min_allele_freq`, and `max_allele_freq`.
#'
#' @return A `tidybreed_pop` object containing:
#'   - DuckDB connection
#'   - Population metadata
#'   - Genome tables (genome_meta, genome_haplotype, genome_genotype)
#'   - Optionally: founder_haplotypes table (if `n_haplotypes` is provided)
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Simple initialization
#' pop <- initialize_genome(
#'   pop_name = "A",
#'   n_loci = 1000,
#'   n_chr = 10,
#'   chr_len_Mb = 100
#' )
#'
#' # Different chromosome lengths
#' pop <- initialize_genome(
#'   pop_name = "Cattle",
#'   n_loci = 50000,
#'   n_chr = 30,
#'   chr_len_Mb = c(158, 137, 121, 120, 121, 119, 112, 113, 105, 104,
#'                  107, 91, 84, 84, 85, 81, 75, 66, 64, 72,
#'                  71, 61, 52, 62, 42, 51, 45, 46, 51, 42)
#' )
#'
#' # In-memory (not recommended for large simulations)
#' pop <- initialize_genome(
#'   pop_name = "test",
#'   n_loci = 100,
#'   n_chr = 2,
#'   chr_len_Mb = 100,
#'   db_path = ":memory:"
#' )
#'
#' # With founder haplotypes (uniform allele frequencies)
#' pop <- initialize_genome(
#'   pop_name = "B",
#'   n_loci = 1000,
#'   n_chr = 10,
#'   chr_len_Mb = 100,
#'   n_haplotypes = 100,
#'   min_allele_freq = 0.05,
#'   max_allele_freq = 0.95
#' )
#'
#' # With founder haplotypes (fixed 50% allele frequency)
#' pop <- initialize_genome(
#'   pop_name = "C",
#'   n_loci = 1000,
#'   n_chr = 10,
#'   chr_len_Mb = 100,
#'   n_haplotypes = 50,
#'   fixed_allele_freq = 0.5
#' )
#' }
initialize_genome <- function(pop_name,
                              n_loci,
                              n_chr,
                              chr_len_Mb,
                              db_path = NULL,
                              locus_names = NULL,
                              chr_names = NULL,
                              overwrite = FALSE,
                              n_haplotypes = NULL,
                              allele_freq_dist = "uniform",
                              min_allele_freq = 0.01,
                              max_allele_freq = 0.99,
                              fixed_allele_freq = NULL) {

  # Input validation
  stopifnot(is.character(pop_name), length(pop_name) == 1)
  stopifnot(is.numeric(n_loci), length(n_loci) == 1, n_loci > 0)
  stopifnot(is.numeric(n_chr), length(n_chr) == 1, n_chr > 0)
  stopifnot(is.numeric(chr_len_Mb), length(chr_len_Mb) %in% c(1, n_chr))

  # Validate founder haplotype parameters
  if (!is.null(n_haplotypes)) {
    stopifnot(is.numeric(n_haplotypes), length(n_haplotypes) == 1, n_haplotypes > 0)
    stopifnot(allele_freq_dist %in% c("uniform"))

    if (!is.null(fixed_allele_freq)) {
      stopifnot(is.numeric(fixed_allele_freq), length(fixed_allele_freq) == 1)
      stopifnot(fixed_allele_freq >= 0, fixed_allele_freq <= 1)
    } else {
      stopifnot(is.numeric(min_allele_freq), length(min_allele_freq) == 1)
      stopifnot(is.numeric(max_allele_freq), length(max_allele_freq) == 1)
      stopifnot(min_allele_freq >= 0, min_allele_freq <= 1)
      stopifnot(max_allele_freq >= 0, max_allele_freq <= 1)
      stopifnot(min_allele_freq < max_allele_freq)
    }
  }

  # Expand chr_len_Mb if single value
  if (length(chr_len_Mb) == 1) {
    chr_len_Mb <- rep(chr_len_Mb, n_chr)
  }

  # Set default db_path
  if (is.null(db_path)) {
    db_path <- paste0(pop_name, "_tidybreed.duckdb")
  }

  # Check if file exists
  if (db_path != ":memory:" && file.exists(db_path)) {
    if (!overwrite) {
      stop(
        "Database file '", db_path, "' already exists. ",
        "Set overwrite = TRUE to replace it.",
        call. = FALSE
      )
    } else {
      message("Overwriting existing database: ", db_path)
      file.remove(db_path)
    }
  }

  # Create DuckDB connection
  drv <- duckdb::duckdb()
  db_conn <- DBI::dbConnect(drv, dbdir = db_path)

  # Generate locus names if not provided
  if (is.null(locus_names)) {
    locus_names <- paste0("Locus_", seq_len(n_loci))
  } else {
    stopifnot(length(locus_names) == n_loci)
  }

  # Generate chromosome names if not provided
  if (is.null(chr_names)) {
    chr_names <- as.character(seq_len(n_chr))
  } else {
    stopifnot(length(chr_names) == n_chr)
  }

  # Assign loci to chromosomes (evenly distributed)
  loci_per_chr <- diff(round(seq(0, n_loci, length.out = n_chr + 1)))
  chr_assignment <- rep(seq_len(n_chr), times = loci_per_chr)

  # Generate positions within each chromosome
  pos_Mb <- numeric(n_loci)
  for (i in seq_len(n_chr)) {
    chr_loci <- which(chr_assignment == i)
    n_chr_loci <- length(chr_loci)
    # Evenly space loci across chromosome
    pos_Mb[chr_loci] <- seq(0, chr_len_Mb[i], length.out = n_chr_loci + 2)[2:(n_chr_loci + 1)]
  }

  # Create genome_meta table
  genome_meta <- tibble::tibble(
    locus_id = seq_len(n_loci),
    locus_name = locus_names,
    chr = chr_assignment,
    chr_name = chr_names[chr_assignment],
    pos_Mb = pos_Mb
  )

  DBI::dbWriteTable(db_conn, "genome_meta", genome_meta, overwrite = TRUE)

  # Create empty genome_haplotype table
  # Structure: id_ind, parent_origin (1 or 2), then one column per locus
  hap_cols <- paste0("locus_", seq_len(n_loci))
  hap_schema <- paste(
    "id_ind VARCHAR,",
    "parent_origin INTEGER,",
    paste(hap_cols, "UTINYINT", collapse = ", ")
  )

  DBI::dbExecute(
    db_conn,
    paste0("CREATE TABLE genome_haplotype (", hap_schema, ")")
  )

  # Create empty genome_genotype table
  # Structure: id_ind, then one column per locus with genotype (0/1/2)
  geno_schema <- paste(
    "id_ind VARCHAR,",
    paste(hap_cols, "UTINYINT", collapse = ", ")
  )

  DBI::dbExecute(
    db_conn,
    paste0("CREATE TABLE genome_genotype (", geno_schema, ")")
  )

  # Eagerly create all core individual/trait tables (empty).
  # This lets users call get_table() and mutate_table() on any table
  # immediately after initialize_genome(), before any data is added.

  DBI::dbExecute(db_conn, "
    CREATE TABLE ind_meta (
      id_ind      VARCHAR PRIMARY KEY,
      id_parent_1 VARCHAR,
      id_parent_2 VARCHAR,
      line        VARCHAR,
      sex         VARCHAR
    )
  ")

  DBI::dbExecute(db_conn, "
    CREATE TABLE trait_effect_cov (
      effect_name VARCHAR,
      trait_1     VARCHAR,
      trait_2     VARCHAR,
      cov         DOUBLE,
      PRIMARY KEY (effect_name, trait_1, trait_2)
    )
  ")

  # Generate founder haplotypes if requested
  tables_created <- c("genome_meta", "genome_haplotype", "genome_genotype",
                      "ind_meta", "trait_effect_cov")
  founder_metadata <- list()

  if (!is.null(n_haplotypes)) {
    message("Generating ", n_haplotypes, " founder haplotypes...")

    # Determine allele frequencies for each locus
    if (!is.null(fixed_allele_freq)) {
      # Use fixed allele frequency for all loci
      allele_freqs <- rep(fixed_allele_freq, n_loci)
      message("  Using fixed allele frequency: ", fixed_allele_freq)
    } else {
      # Sample allele frequencies from uniform distribution
      if (allele_freq_dist == "uniform") {
        allele_freqs <- stats::runif(n_loci, min = min_allele_freq, max = max_allele_freq)
        message("  Allele frequencies sampled from uniform(", min_allele_freq, ", ", max_allele_freq, ")")
      }
    }

    # Add allele frequencies to genome_meta table
    genome_meta$founder_allele_freq <- allele_freqs
    DBI::dbWriteTable(db_conn, "genome_meta", genome_meta, overwrite = TRUE)

    # Generate haplotypes by sampling from Bernoulli(allele_freq) for each locus
    # Create matrix: rows = haplotypes, columns = loci
    haplotype_matrix <- matrix(0L, nrow = n_haplotypes, ncol = n_loci)

    for (j in seq_len(n_loci)) {
      # Sample alleles (0 or 1) from Bernoulli with probability = allele_freq
      haplotype_matrix[, j] <- stats::rbinom(n_haplotypes, size = 1, prob = allele_freqs[j])
    }

    # Create founder_haplotypes table
    founder_haplotypes <- tibble::tibble(
      hap_id = paste0("hap_", seq_len(n_haplotypes))
    )

    # Add locus columns
    for (j in seq_len(n_loci)) {
      founder_haplotypes[[paste0("locus_", j)]] <- haplotype_matrix[, j]
    }

    # Write to database
    DBI::dbWriteTable(db_conn, "founder_haplotypes", founder_haplotypes, overwrite = TRUE)

    # Update tables list and metadata
    tables_created <- c(tables_created, "founder_haplotypes")
    founder_metadata <- list(
      n_haplotypes = n_haplotypes,
      allele_freq_dist = if (!is.null(fixed_allele_freq)) "fixed" else allele_freq_dist,
      fixed_allele_freq = fixed_allele_freq,
      min_allele_freq = if (is.null(fixed_allele_freq)) min_allele_freq else NULL,
      max_allele_freq = if (is.null(fixed_allele_freq)) max_allele_freq else NULL
    )

    message("  Created founder_haplotypes table")
  }

  # Create tidybreed_pop object
  pop <- new_tidybreed_pop(
    db_conn = db_conn,
    pop_name = pop_name,
    db_path = db_path,
    tables = tables_created,
    metadata = c(
      list(
        n_loci = n_loci,
        n_chr = n_chr,
        chr_len_Mb = chr_len_Mb,
        chr_names = chr_names
      ),
      founder_metadata
    )
  )

  validate_tidybreed_pop(pop)

  # Create remaining core tables (trait_meta, trait_effects, trait_random_effects,
  # ind_phenotype, ind_tbv, ind_ebv) and register them in pop$tables.
  pop <- ensure_trait_tables(pop)

  message("Initialized population '", pop_name, "' with ", n_loci, " loci across ", n_chr, " chromosomes")
  if (!is.null(n_haplotypes)) {
    message("  Generated ", n_haplotypes, " founder haplotypes")
  }
  message("Database: ", db_path)

  return(pop)
}
