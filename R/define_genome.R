#' Define the genome structure of a breeding population
#'
#' @description
#' Adds genome tables (`genome_meta`, `genome_haplotype`, `genome_genotype`) to
#' a population opened with [open_pop()]. Pipe-friendly — accepts a
#' `tidybreed_pop` and returns a `tidybreed_pop`.
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
                          locus_names = NULL,
                          chr_names   = NULL) {

  stopifnot(inherits(pop, "tidybreed_pop"))
  stopifnot(is.numeric(n_loci), length(n_loci) == 1, n_loci > 0)
  stopifnot(is.numeric(n_chr),  length(n_chr)  == 1, n_chr  > 0)
  stopifnot(is.numeric(chr_len_Mb), length(chr_len_Mb) %in% c(1, n_chr))

  db_conn <- pop$db_conn

  # Expand chr_len_Mb if single value
  if (length(chr_len_Mb) == 1)
    chr_len_Mb <- rep(chr_len_Mb, n_chr)

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
  loci_per_chr   <- diff(round(seq(0, n_loci, length.out = n_chr + 1)))
  chr_assignment <- rep(seq_len(n_chr), times = loci_per_chr)

  # Generate positions within each chromosome (evenly spaced)
  pos_Mb <- numeric(n_loci)
  for (i in seq_len(n_chr)) {
    chr_loci   <- which(chr_assignment == i)
    n_chr_loci <- length(chr_loci)
    pos_Mb[chr_loci] <- seq(0, chr_len_Mb[i], length.out = n_chr_loci + 2)[2:(n_chr_loci + 1)]
  }

  # Write genome_meta
  genome_meta_df <- tibble::tibble(
    locus_id   = seq_len(n_loci),
    locus_name = locus_names,
    chr        = chr_assignment,
    chr_name   = chr_names[chr_assignment],
    pos_Mb     = pos_Mb
  )
  DBI::dbWriteTable(db_conn, "genome_meta", genome_meta_df, overwrite = TRUE)

  # Create empty genome_haplotype table (wide: one column per locus)
  hap_cols   <- paste0("locus_", seq_len(n_loci))
  hap_schema <- paste(
    "id_ind VARCHAR,",
    "parent_origin INTEGER,",
    paste(hap_cols, "UTINYINT", collapse = ", ")
  )
  DBI::dbExecute(db_conn, paste0("CREATE TABLE genome_haplotype (", hap_schema, ")"))

  # Create empty genome_genotype table (wide: one column per locus, 0/1/2 encoding)
  geno_schema <- paste(
    "id_ind VARCHAR,",
    paste(hap_cols, "UTINYINT", collapse = ", ")
  )
  DBI::dbExecute(db_conn, paste0("CREATE TABLE genome_genotype (", geno_schema, ")"))

  # Register genome-table schema descriptions
  register_schema_meta(db_conn, .genome_table_descriptions())

  # Update pop$tables
  pop$tables <- c(pop$tables, "genome_meta", "genome_haplotype", "genome_genotype")

  chr_len_str <- if (length(unique(chr_len_Mb)) == 1) {
    paste0("all equal to ", chr_len_Mb[1], " Mb")
  } else {
    paste(paste0(chr_names, ":", chr_len_Mb), collapse = ", ")
  }
  message(
    "Defined genome: ", n_loci, " loci across ", n_chr, " chromosomes",
    " | chr lengths (Mb): ", chr_len_str,
    "\n  Tables written: genome_meta, genome_haplotype, genome_genotype"
  )

  pop
}
