#' Define founder haplotypes for a tidybreed population
#'
#' @description
#' Generates a pool of founder haplotypes from per-locus Bernoulli distributions
#' and stores them in the `founder_haplotypes` table. This table is required by
#' [add_founders()] to assign phased alleles to new individuals.
#'
#' Call this function after [initialize_genome()] and before [add_founders()].
#' Also writes a `founder_allele_freq` column to `genome_meta` recording the
#' per-locus allele frequency used during sampling.
#'
#' @param pop A `tidybreed_pop` object returned by [initialize_genome()].
#' @param n_haplotypes Integer scalar. Number of founder haplotypes to generate.
#'   Must be positive.
#' @param fixed_allele_freq Numeric scalar in [0, 1] or `NULL`. When supplied,
#'   every locus is assigned this allele frequency exactly. Overrides
#'   `allele_freq_dist`, `min_allele_freq`, and `max_allele_freq`.
#' @param allele_freq_dist Character scalar. Distribution for sampling per-locus
#'   allele frequencies. Currently only `"uniform"` is supported. Ignored when
#'   `fixed_allele_freq` is supplied.
#' @param min_allele_freq Numeric scalar in (0, 1). Minimum allele frequency for
#'   the uniform distribution. Default `0.01`.
#' @param max_allele_freq Numeric scalar in (0, 1). Maximum allele frequency for
#'   the uniform distribution. Default `0.99`.
#'
#' @return The `tidybreed_pop` object (invisibly), with `founder_haplotypes`
#'   registered in `pop$tables`.
#'
#' @seealso [initialize_genome()], [add_founders()]
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Uniform allele frequencies
#' pop <- initialize_genome(
#'   pop_name = "A",
#'   n_loci = 1000,
#'   n_chr = 10,
#'   chr_len_Mb = 100
#' ) |>
#'   define_founder_haplotypes(
#'     n_haplotypes  = 100,
#'     min_allele_freq = 0.05,
#'     max_allele_freq = 0.95
#'   )
#'
#' # Fixed allele frequency
#' pop <- initialize_genome(
#'   pop_name = "B",
#'   n_loci = 1000,
#'   n_chr = 10,
#'   chr_len_Mb = 100
#' ) |>
#'   define_founder_haplotypes(
#'     n_haplotypes      = 50,
#'     fixed_allele_freq = 0.5
#'   )
#' }
define_founder_haplotypes <- function(pop,
                                      n_haplotypes,
                                      fixed_allele_freq = NULL,
                                      allele_freq_dist  = "uniform",
                                      min_allele_freq   = 0.01,
                                      max_allele_freq   = 0.99) {

  stopifnot(inherits(pop, "tidybreed_pop"))
  validate_tidybreed_pop(pop)

  # Validate n_haplotypes
  stopifnot(is.numeric(n_haplotypes), length(n_haplotypes) == 1, n_haplotypes > 0)
  n_haplotypes <- as.integer(n_haplotypes)

  # Validate allele frequency parameters
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

  # Guard: genome_meta must exist
  if (!"genome_meta" %in% pop$tables) {
    stop(
      "genome_meta table does not exist. ",
      "Call initialize_genome() before define_founder_haplotypes().",
      call. = FALSE
    )
  }

  # Guard: no silent overwrite of existing table
  if ("founder_haplotypes" %in% pop$tables ||
      DBI::dbExistsTable(pop$db_conn, "founder_haplotypes")) {
    stop(
      "founder_haplotypes table already exists in this population. ",
      "Drop it before calling define_founder_haplotypes() again.",
      call. = FALSE
    )
  }

  # Derive n_loci from genome_meta
  n_loci <- DBI::dbGetQuery(pop$db_conn, "SELECT COUNT(*) AS n FROM genome_meta")$n
  n_loci <- as.integer(n_loci)
  if (n_loci == 0L) {
    stop("genome_meta is empty. Cannot generate founder haplotypes.", call. = FALSE)
  }

  # Read genome_meta to update founder_allele_freq column
  genome_meta <- DBI::dbGetQuery(pop$db_conn, "SELECT * FROM genome_meta")

  message("Generating ", n_haplotypes, " founder haplotypes...")

  # Determine allele frequencies for each locus
  if (!is.null(fixed_allele_freq)) {
    allele_freqs <- rep(fixed_allele_freq, n_loci)
    message("  Using fixed allele frequency: ", fixed_allele_freq)
  } else {
    if (allele_freq_dist == "uniform") {
      allele_freqs <- stats::runif(n_loci, min = min_allele_freq, max = max_allele_freq)
      message("  Allele frequencies sampled from uniform(", min_allele_freq, ", ", max_allele_freq, ")")
    }
  }

  # Write founder_allele_freq to genome_meta
  genome_meta$founder_allele_freq <- allele_freqs
  DBI::dbWriteTable(pop$db_conn, "genome_meta", genome_meta, overwrite = TRUE)

  # Generate haplotypes: rows = haplotypes, columns = loci
  haplotype_matrix <- matrix(0L, nrow = n_haplotypes, ncol = n_loci)
  for (j in seq_len(n_loci)) {
    haplotype_matrix[, j] <- stats::rbinom(n_haplotypes, size = 1, prob = allele_freqs[j])
  }

  # Build founder_haplotypes tibble
  founder_haplotypes <- tibble::tibble(
    hap_id = paste0("hap_", seq_len(n_haplotypes))
  )
  for (j in seq_len(n_loci)) {
    founder_haplotypes[[paste0("locus_", j)]] <- haplotype_matrix[, j]
  }

  # Write to database
  DBI::dbWriteTable(pop$db_conn, "founder_haplotypes", founder_haplotypes, overwrite = FALSE)

  # Register in pop$tables
  pop$tables <- unique(c(pop$tables, "founder_haplotypes"))

  # Register schema metadata
  register_schema_meta(pop$db_conn, rbind(
    .sm_tbl("founder_haplotypes",
            "Pool of founder haplotypes sampled from per-locus Bernoulli distributions. Used by add_founders() to assign phased alleles."),
    .sm_col("founder_haplotypes", "hap_id",
            "Haplotype identifier (e.g. 'hap_1', 'hap_2')")
  ))

  message("  Created founder_haplotypes table")

  invisible(pop)
}
