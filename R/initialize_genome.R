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
#' To create the founder haplotype pool needed by [add_founders()], call
#' [define_founder_haplotypes()] after this function.
#'
#' @param pop_name Character scalar. Name of the population; used as the prefix
#'   for the default `.duckdb` filename.
#' @param n_loci Integer scalar. Total number of loci to simulate.
#' @param n_chr Integer scalar. Number of chromosomes.
#' @param chr_len_Mb Numeric scalar or numeric vector of length `n_chr`.
#'   Chromosome length(s) in megabases. A single scalar applies the same
#'   length to all chromosomes; a vector specifies each chromosome separately.
#' @param db_path Character scalar. Path to the DuckDB database file. Default
#'   creates a file in the current directory named
#'   `{pop_name}_tidybreed.duckdb`. Use `":memory:"` for an in-memory
#'   database (not recommended for large simulations).
#' @param locus_names Character vector of length `n_loci` or `NULL`. Custom
#'   locus names. When `NULL` (default), names are auto-generated as
#'   `"Locus_1"`, `"Locus_2"`, etc.
#' @param chr_names Character vector of length `n_chr` or `NULL`. Custom
#'   chromosome names. When `NULL` (default), chromosomes are numbered
#'   `1, 2, ..., n_chr`.
#' @param overwrite Logical scalar. If `TRUE` and `db_path` already exists,
#'   the file is deleted and recreated. Default `FALSE`.
#'
#' @return A `tidybreed_pop` object containing:
#'   - DuckDB connection
#'   - Population metadata
#'   - Genome tables (genome_meta, genome_haplotype, genome_genotype)
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
#' # With founder haplotypes — call define_founder_haplotypes() after init
#' pop <- initialize_genome(
#'   pop_name = "B",
#'   n_loci = 1000,
#'   n_chr = 10,
#'   chr_len_Mb = 100
#' ) |>
#'   define_founder_haplotypes(n_haplotypes = 100,
#'                             min_allele_freq = 0.05,
#'                             max_allele_freq = 0.95)
#' }
initialize_genome <- function(pop_name,
                              n_loci,
                              n_chr,
                              chr_len_Mb,
                              db_path = NULL,
                              locus_names = NULL,
                              chr_names = NULL,
                              overwrite = FALSE) {

  # Input validation
  stopifnot(is.character(pop_name), length(pop_name) == 1)
  stopifnot(is.numeric(n_loci), length(n_loci) == 1, n_loci > 0)
  stopifnot(is.numeric(n_chr), length(n_chr) == 1, n_chr > 0)
  stopifnot(is.numeric(chr_len_Mb), length(chr_len_Mb) %in% c(1, n_chr))

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

  # Create _schema_meta system table (must exist before any register_schema_meta calls)
  DBI::dbExecute(db_conn, "
    CREATE TABLE _schema_meta (
      id_schema_meta INTEGER PRIMARY KEY,
      object_type    VARCHAR NOT NULL,
      table_name     VARCHAR NOT NULL,
      column_name    VARCHAR,
      description    VARCHAR NOT NULL,
      notes          VARCHAR
    )
  ")

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
      line_name   VARCHAR,
      sex         VARCHAR
    )
  ")

  DBI::dbExecute(db_conn, "
    CREATE TABLE trait_effect_cov (
      id_trait_effect_cov INTEGER PRIMARY KEY,
      effect_name         VARCHAR,
      trait_name_1        VARCHAR,
      trait_name_2        VARCHAR,
      cov_value           DOUBLE
    )
  ")

  DBI::dbExecute(db_conn, "
    CREATE TABLE genome_effects (
      id_genome_effect   INTEGER PRIMARY KEY,
      locus_name         VARCHAR NOT NULL,
      line_name          VARCHAR,
      trait_name         VARCHAR NOT NULL,
      genome_effect_type VARCHAR NOT NULL,
      genome_value       DOUBLE  NOT NULL,
      base_allele_freq   DOUBLE
    )
  ")

  # Create new phenotype-layer tables (eagerly, empty)
  DBI::dbExecute(db_conn, "
    CREATE TABLE phenotype_meta (
      id_phenotype_meta        INTEGER PRIMARY KEY,
      phenotype_name           VARCHAR UNIQUE NOT NULL,
      type                     VARCHAR,
      mean                     DOUBLE DEFAULT 0,
      expressed_sex            VARCHAR DEFAULT 'both',
      repeatable               BOOLEAN DEFAULT FALSE,
      min_value                DOUBLE,
      max_value                DOUBLE,
      prevalence               DOUBLE,
      thresholds               VARCHAR,
      cat_values               VARCHAR,
      cat_names                VARCHAR,
      store_liability          BOOLEAN DEFAULT FALSE,
      missing_component_action VARCHAR DEFAULT 'skip',
      formula_tbv              VARCHAR,
      formula                  VARCHAR
    )
  ")

  DBI::dbExecute(db_conn, "
    CREATE TABLE phenotype_components (
      id_phenotype_comp   INTEGER PRIMARY KEY,
      phenotype_name      VARCHAR NOT NULL,
      source_trait_name   VARCHAR NOT NULL,
      contributor_type    VARCHAR NOT NULL,
      group_column        VARCHAR,
      group_table         VARCHAR DEFAULT 'ind_meta',
      aggregation         VARCHAR DEFAULT 'sum',
      weight              DOUBLE  DEFAULT 1.0,
      weight_type         VARCHAR DEFAULT 'fixed',
      covariate_name      VARCHAR,
      covariate_table     VARCHAR,
      poly_order          INTEGER,
      poly_scale_min      DOUBLE,
      poly_scale_max      DOUBLE,
      genome_effect_types VARCHAR DEFAULT 'additive',
      missing_action      VARCHAR DEFAULT 'skip',
      contributor_filter  VARCHAR
    )
  ")

  DBI::dbExecute(db_conn, "
    CREATE TABLE phenotype_residual_cov (
      id_residual_cov  INTEGER PRIMARY KEY,
      phenotype_name_1 VARCHAR NOT NULL,
      phenotype_name_2 VARCHAR NOT NULL,
      cov_value        DOUBLE  NOT NULL,
      condition_column VARCHAR,
      condition_table  VARCHAR DEFAULT 'ind_meta',
      condition_level  VARCHAR,
      weight_type      VARCHAR DEFAULT 'fixed',
      poly_order       INTEGER
    )
  ")

  # Register genome-layer schema descriptions
  register_schema_meta(db_conn, .genome_layer_descriptions())

  tables_created <- c("_schema_meta",
                      "genome_meta", "genome_haplotype", "genome_genotype",
                      "ind_meta", "trait_effect_cov", "genome_effects",
                      "phenotype_meta", "phenotype_components",
                      "phenotype_residual_cov")

  # Create tidybreed_pop object
  pop <- new_tidybreed_pop(
    db_conn  = db_conn,
    pop_name = pop_name,
    db_path  = db_path,
    tables   = tables_created
  )

  validate_tidybreed_pop(pop)

  # Create remaining core tables (trait_meta, trait_effects, trait_random_effects,
  # ind_phenotype, ind_tbv, ind_ebv) and register them in pop$tables.
  pop <- ensure_trait_tables(pop)

  message("Initialized population '", pop_name, "' with ", n_loci, " loci across ", n_chr, " chromosomes")
  message("Database: ", db_path)

  return(pop)
}
