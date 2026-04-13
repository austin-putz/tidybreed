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
#'   names are generated as "SNP_1", "SNP_2", etc.
#' @param chr_names Character vector of chromosome names. If NULL (default),
#'   chromosomes are numbered 1, 2, ..., n_chr.
#' @param overwrite Logical. If TRUE and `db_path` exists, overwrite it.
#'   Default is FALSE.
#'
#' @return A `tidybreed_pop` object containing:
#'   - DuckDB connection
#'   - Population metadata
#'   - Three genome tables (genome_meta, genome_haplotype, genome_genotype)
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

  # Generate locus names if not provided
  if (is.null(locus_names)) {
    locus_names <- paste0("SNP_", seq_len(n_loci))
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
  # Structure: ind_id, parent_origin (1 or 2), then one column per locus
  hap_cols <- paste0("locus_", seq_len(n_loci))
  hap_schema <- paste(
    "ind_id VARCHAR,",
    "parent_origin INTEGER,",
    paste(hap_cols, "INTEGER", collapse = ", ")
  )

  DBI::dbExecute(
    db_conn,
    paste0("CREATE TABLE genome_haplotype (", hap_schema, ")")
  )

  # Create empty genome_genotype table
  # Structure: ind_id, then one column per locus with genotype (0/1/2)
  geno_schema <- paste(
    "ind_id VARCHAR,",
    paste(hap_cols, "INTEGER", collapse = ", ")
  )

  DBI::dbExecute(
    db_conn,
    paste0("CREATE TABLE genome_genotype (", geno_schema, ")")
  )

  # Create tidybreed_pop object
  pop <- new_tidybreed_pop(
    db_conn = db_conn,
    pop_name = pop_name,
    db_path = db_path,
    tables = c("genome_meta", "genome_haplotype", "genome_genotype"),
    metadata = list(
      n_loci = n_loci,
      n_chr = n_chr,
      chr_len_Mb = chr_len_Mb,
      chr_names = chr_names
    )
  )

  validate_tidybreed_pop(pop)

  message("Initialized population '", pop_name, "' with ", n_loci, " loci across ", n_chr, " chromosomes")
  message("Database: ", db_path)

  return(pop)
}
