# Internal helpers for define_founder_haplotypes(). None are exported.

# ---------------------------------------------------------------------------
# Allele frequency generators (independent methods)
# ---------------------------------------------------------------------------

.gen_allele_freqs_fixed <- function(n_loci, allele_freq) {
  rep(allele_freq, n_loci)
}

.gen_allele_freqs_uniform <- function(n_loci, min_allele_freq, max_allele_freq) {
  stats::runif(n_loci, min = min_allele_freq, max = max_allele_freq)
}

.gen_allele_freqs_beta <- function(n_loci, shape1, shape2) {
  stats::rbeta(n_loci, shape1 = shape1, shape2 = shape2)
}

.gen_allele_freqs_balding_nichols <- function(n_loci, fst, mean_freq) {
  shape1 <- mean_freq * (1 - fst) / fst
  shape2 <- (1 - mean_freq) * (1 - fst) / fst
  if (shape1 <= 0 || shape2 <= 0) {
    stop(
      "Balding-Nichols parameterization produced invalid Beta shape parameters. ",
      "Ensure fst is in (0, 1) and mean_freq is in (0, 1).",
      call. = FALSE
    )
  }
  stats::rbeta(n_loci, shape1 = shape1, shape2 = shape2)
}

# ---------------------------------------------------------------------------
# Haplotype matrix generators
# ---------------------------------------------------------------------------

.gen_haplotypes_from_freqs <- function(n_haplotypes, allele_freqs) {
  n_loci <- length(allele_freqs)
  haplotype_matrix <- matrix(0L, nrow = n_haplotypes, ncol = n_loci)
  for (j in seq_len(n_loci)) {
    haplotype_matrix[, j] <- stats::rbinom(n_haplotypes, size = 1L,
                                            prob = allele_freqs[j])
  }
  haplotype_matrix
}

.gen_haplotypes_mosaic <- function(n_haplotypes, n_loci, genome_meta,
                                    n_templates, switch_rate) {
  chrs <- sort(unique(genome_meta$chr))

  # Build traversal order (loci in chr/pos order) and d_Mb (distance from
  # previous locus within the same chromosome; 0 at each chromosome's start)
  traversal_order <- integer(n_loci)
  d_Mb            <- numeric(n_loci)
  k <- 0L
  for (chr_id in chrs) {
    chr_mask        <- genome_meta$chr == chr_id
    chr_loci        <- which(chr_mask)
    chr_pos         <- genome_meta$pos_Mb[chr_mask]
    ord             <- order(chr_pos)
    chr_loci_sorted <- chr_loci[ord]
    chr_pos_sorted  <- chr_pos[ord]
    n_chr_loci      <- length(chr_loci_sorted)
    idx             <- seq(k + 1L, k + n_chr_loci)
    traversal_order[idx]              <- chr_loci_sorted
    d_Mb[chr_loci_sorted[1L]]         <- 0
    if (n_chr_loci > 1L) {
      d_Mb[chr_loci_sorted[-1L]] <- diff(chr_pos_sorted)
    }
    k <- k + n_chr_loci
  }

  # Per-locus switch probability (Haldane-like)
  p_switch <- 1 - exp(-switch_rate * d_Mb)

  # Identify the first locus of each chromosome (template resets here)
  chr_first_loci <- vapply(chrs, function(chr_id) {
    chr_mask <- genome_meta$chr == chr_id
    chr_loci <- which(chr_mask)
    chr_pos  <- genome_meta$pos_Mb[chr_mask]
    chr_loci[which.min(chr_pos)]
  }, integer(1L))
  is_chr_first <- logical(n_loci)
  is_chr_first[chr_first_loci] <- TRUE

  # Generate n_templates template haplotypes from uniform allele frequencies
  template_freqs  <- stats::runif(n_loci, 0.01, 0.99)
  template_matrix <- matrix(0L, nrow = n_templates, ncol = n_loci)
  for (j in seq_len(n_loci)) {
    template_matrix[, j] <- stats::rbinom(n_templates, size = 1L,
                                           prob = template_freqs[j])
  }

  # Pre-generate all switch-decision uniform draws (vectorised over haplotypes)
  U <- matrix(stats::runif(n_haplotypes * n_loci),
               nrow = n_haplotypes, ncol = n_loci)

  # Build each haplotype as a mosaic of templates
  haplotype_matrix <- matrix(0L, nrow = n_haplotypes, ncol = n_loci)
  for (i in seq_len(n_haplotypes)) {
    current_template <- sample.int(n_templates, 1L)
    for (step in seq_len(n_loci)) {
      j <- traversal_order[step]
      if (is_chr_first[j]) {
        current_template <- sample.int(n_templates, 1L)
      } else if (U[i, step] < p_switch[j]) {
        current_template <- sample.int(n_templates, 1L)
      }
      haplotype_matrix[i, j] <- template_matrix[current_template, j]
    }
  }

  haplotype_matrix
}

.gen_haplotypes_gaussian_copula <- function(n_haplotypes, n_loci, genome_meta,
                                             decay_rate) {
  chrs <- sort(unique(genome_meta$chr))

  # Sample per-locus allele frequencies and convert to latent thresholds
  allele_freqs <- stats::runif(n_loci, 0.01, 0.99)
  thresholds   <- stats::qnorm(allele_freqs)

  # Build traversal order, d_Mb, and rho per locus
  traversal_order <- integer(n_loci)
  d_Mb            <- numeric(n_loci)
  k <- 0L
  for (chr_id in chrs) {
    chr_mask        <- genome_meta$chr == chr_id
    chr_loci        <- which(chr_mask)
    chr_pos         <- genome_meta$pos_Mb[chr_mask]
    ord             <- order(chr_pos)
    chr_loci_sorted <- chr_loci[ord]
    chr_pos_sorted  <- chr_pos[ord]
    n_chr_loci      <- length(chr_loci_sorted)
    idx             <- seq(k + 1L, k + n_chr_loci)
    traversal_order[idx]              <- chr_loci_sorted
    d_Mb[chr_loci_sorted[1L]]         <- 0
    if (n_chr_loci > 1L) {
      d_Mb[chr_loci_sorted[-1L]] <- diff(chr_pos_sorted)
    }
    k <- k + n_chr_loci
  }

  rho <- exp(-decay_rate * d_Mb)

  # Force rho = 0 at each chromosome's first locus so the AR(1) restarts
  # independently (z_curr = eps ~ N(0,1) when rho = 0)
  chr_first_loci <- vapply(chrs, function(chr_id) {
    chr_mask <- genome_meta$chr == chr_id
    chr_loci <- which(chr_mask)
    chr_pos  <- genome_meta$pos_Mb[chr_mask]
    chr_loci[which.min(chr_pos)]
  }, integer(1L))
  rho[chr_first_loci] <- 0

  # Pre-generate noise matrix: n_haplotypes × n_loci
  eps <- matrix(stats::rnorm(n_haplotypes * n_loci),
                nrow = n_haplotypes, ncol = n_loci)

  # Walk loci in traversal order; vectorise over haplotypes at each step
  haplotype_matrix <- matrix(0L, nrow = n_haplotypes, ncol = n_loci)
  z_prev <- stats::rnorm(n_haplotypes)   # initial latent vector

  for (step in seq_len(n_loci)) {
    j      <- traversal_order[step]
    rho_j  <- rho[j]
    scale  <- sqrt(max(0, 1 - rho_j^2))
    z_curr <- rho_j * z_prev + scale * eps[, step]
    haplotype_matrix[, j] <- as.integer(z_curr < thresholds[j])
    z_prev <- z_curr
  }

  list(haplotype_matrix = haplotype_matrix,
       allele_freqs     = allele_freqs)
}

# ---------------------------------------------------------------------------
# Shared DB write (called by all methods)
# ---------------------------------------------------------------------------

.write_founder_haplotypes <- function(pop, haplotype_matrix, allele_freqs) {
  n_haplotypes <- nrow(haplotype_matrix)
  n_loci       <- ncol(haplotype_matrix)

  # Write founder_allele_freq back to genome_meta
  genome_meta <- DBI::dbGetQuery(pop$db_conn, "SELECT * FROM genome_meta")
  genome_meta$founder_allele_freq <- allele_freqs
  DBI::dbWriteTable(pop$db_conn, "genome_meta", genome_meta, overwrite = TRUE)

  # Build founder_haplotypes tibble: hap_id + locus_1..locus_n
  fh <- tibble::tibble(hap_id = paste0("hap_", seq_len(n_haplotypes)))
  for (j in seq_len(n_loci)) {
    fh[[paste0("locus_", j)]] <- haplotype_matrix[, j]
  }

  DBI::dbWriteTable(pop$db_conn, "founder_haplotypes", fh, overwrite = FALSE)

  pop$tables <- unique(c(pop$tables, "founder_haplotypes"))

  register_schema_meta(pop$db_conn, rbind(
    .sm_tbl("founder_haplotypes",
            "Pool of founder haplotypes sampled by define_founder_haplotypes(). Rows are sampled by add_founders() to assign phased alleles."),
    .sm_col("founder_haplotypes", "hap_id",
            "Haplotype identifier (e.g. 'hap_1', 'hap_2')")
  ))

  pop
}
