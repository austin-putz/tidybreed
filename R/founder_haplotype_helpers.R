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

.gen_haplotypes_mosaic <- function(n_haplotypes, n_loci, locus_map,
                                    n_templates, switch_rate) {
  chrs <- sort(unique(locus_map$chr))

  # Build traversal order (loci in chr/pos order) and d_cM (genetic distance from
  # previous locus within the same chromosome; 0 at each chromosome's start).
  # switch_rate is per cM; locus_map is the resolved default genetic map.
  traversal_order <- integer(n_loci)
  d_cM            <- numeric(n_loci)
  k <- 0L
  for (chr_id in chrs) {
    chr_mask        <- locus_map$chr == chr_id
    chr_loci        <- which(chr_mask)
    chr_pos         <- locus_map$pos_cM[chr_mask]
    ord             <- order(chr_pos)
    chr_loci_sorted <- chr_loci[ord]
    chr_pos_sorted  <- chr_pos[ord]
    n_chr_loci      <- length(chr_loci_sorted)
    idx             <- seq(k + 1L, k + n_chr_loci)
    traversal_order[idx]              <- chr_loci_sorted
    d_cM[chr_loci_sorted[1L]]         <- 0
    if (n_chr_loci > 1L) {
      d_cM[chr_loci_sorted[-1L]] <- diff(chr_pos_sorted)
    }
    k <- k + n_chr_loci
  }

  # Per-locus switch probability (Haldane-like), per cM
  p_switch <- 1 - exp(-switch_rate * d_cM)

  # Identify the first locus of each chromosome (template resets here)
  chr_first_loci <- vapply(chrs, function(chr_id) {
    chr_mask <- locus_map$chr == chr_id
    chr_loci <- which(chr_mask)
    chr_pos  <- locus_map$pos_cM[chr_mask]
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

.gen_haplotypes_gaussian_copula <- function(n_haplotypes, n_loci, locus_map,
                                             decay_rate) {
  chrs <- sort(unique(locus_map$chr))

  # Sample per-locus allele frequencies and convert to latent thresholds
  allele_freqs <- stats::runif(n_loci, 0.01, 0.99)
  thresholds   <- stats::qnorm(allele_freqs)

  # Build traversal order, d_cM (genetic distance), and rho per locus.
  # decay_rate is per cM; locus_map is the resolved default genetic map.
  traversal_order <- integer(n_loci)
  d_cM            <- numeric(n_loci)
  k <- 0L
  for (chr_id in chrs) {
    chr_mask        <- locus_map$chr == chr_id
    chr_loci        <- which(chr_mask)
    chr_pos         <- locus_map$pos_cM[chr_mask]
    ord             <- order(chr_pos)
    chr_loci_sorted <- chr_loci[ord]
    chr_pos_sorted  <- chr_pos[ord]
    n_chr_loci      <- length(chr_loci_sorted)
    idx             <- seq(k + 1L, k + n_chr_loci)
    traversal_order[idx]              <- chr_loci_sorted
    d_cM[chr_loci_sorted[1L]]         <- 0
    if (n_chr_loci > 1L) {
      d_cM[chr_loci_sorted[-1L]] <- diff(chr_pos_sorted)
    }
    k <- k + n_chr_loci
  }

  rho <- exp(-decay_rate * d_cM)

  # Force rho = 0 at each chromosome's first locus so the AR(1) restarts
  # independently (z_curr = eps ~ N(0,1) when rho = 0)
  chr_first_loci <- vapply(chrs, function(chr_id) {
    chr_mask <- locus_map$chr == chr_id
    chr_loci <- which(chr_mask)
    chr_pos  <- locus_map$pos_cM[chr_mask]
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

.write_founder_haplotypes <- function(pop, haplotype_matrix, allele_freqs,
                                       line_name = NULL) {
  n_haplotypes <- nrow(haplotype_matrix)
  n_loci       <- ncol(haplotype_matrix)

  # Write founder_allele_freq into genome_meta via ALTER + UPDATE (last call wins
  # for multi-line). NEVER rewrite the table with dbWriteTable(overwrite): that
  # re-infers column types from R and would demote pos_bp from BIGINT. The UPDATE
  # is keyed on locus_id via a registered temp relation (RNG-safe; no dbWriteTable).
  if (!"founder_allele_freq" %in% DBI::dbListFields(pop$db_conn, "genome_meta")) {
    DBI::dbExecute(pop$db_conn,
      "ALTER TABLE genome_meta ADD COLUMN founder_allele_freq DOUBLE")
  }
  # allele_freqs is in locus_id order (matching haplotype_matrix columns); align it
  # positionally to the sorted locus_id list.
  lid <- DBI::dbGetQuery(pop$db_conn,
    "SELECT locus_id FROM genome_meta ORDER BY locus_id")$locus_id
  freq_df <- data.frame(locus_id = lid, founder_allele_freq = allele_freqs)
  duckdb::duckdb_register(pop$db_conn, "__tmp_founder_freq", freq_df)
  DBI::dbExecute(pop$db_conn,
    "UPDATE genome_meta AS g SET founder_allele_freq = t.founder_allele_freq
       FROM __tmp_founder_freq AS t WHERE g.locus_id = t.locus_id")
  duckdb::duckdb_unregister(pop$db_conn, "__tmp_founder_freq")

  # Locus names in locus_id order (haplotype_matrix columns are in locus_id order)
  gm_order <- DBI::dbGetQuery(
    pop$db_conn, "SELECT locus_name FROM genome_meta ORDER BY locus_id")$locus_name

  # Write LONG founder_haplotypes (one row per haplotype x locus), batched over
  # haplotypes so the full n_hap x n_loci frame is never materialized at once —
  # peak memory is bounded by ~B x n_loci. The pool generator's draws already
  # happened above; only the WRITE is batched, so seeded pools are byte-identical.
  #
  # RNG discipline: dbWriteTable advances R's RNG by a fixed, size-independent
  # amount (temp-name generation); register+INSERT is RNG-neutral. To reproduce
  # the historical single-dbWriteTable advance exactly (preserving downstream
  # seeded simulations), the FIRST batch uses dbWriteTable — which also creates
  # the table when new — and any remaining batches use register+INSERT. One
  # transaction makes the whole pool write atomic.
  is_new_table <- !DBI::dbExistsTable(pop$db_conn, "founder_haplotypes")
  lname <- if (!is.null(line_name)) line_name else NA_character_
  B     <- resolve_batch_size(n_haplotypes, n_loci)
  conn  <- pop$db_conn

  DBI::dbExecute(conn, "BEGIN TRANSACTION")
  committed <- FALSE
  on.exit(if (!committed) try(DBI::dbExecute(conn, "ROLLBACK"), silent = TRUE),
          add = TRUE)

  first_batch <- TRUE
  for (h_start in seq.int(1L, n_haplotypes, by = B)) {
    h_idx <- h_start:min(h_start + B - 1L, n_haplotypes)
    sub   <- haplotype_matrix[h_idx, , drop = FALSE]      # length(h_idx) x n_loci
    fh <- data.frame(
      line_name    = lname,
      haplotype_id = rep(h_idx, times = n_loci),
      locus_name   = rep(gm_order, each = length(h_idx)),
      allele       = as.integer(as.vector(sub)),
      stringsAsFactors = FALSE
    )
    if (first_batch) {
      DBI::dbWriteTable(conn, "founder_haplotypes", fh, append = TRUE)
      first_batch <- FALSE
    } else {
      duckdb::duckdb_register(conn, "__tmp_fh", fh)
      DBI::dbExecute(conn, paste0(
        "INSERT INTO founder_haplotypes ",
        "(line_name, haplotype_id, locus_name, allele) ",
        "SELECT line_name, haplotype_id, locus_name, allele FROM __tmp_fh"))
      duckdb::duckdb_unregister(conn, "__tmp_fh")
    }
    rm(sub, fh)
    gc(verbose = FALSE)
  }

  DBI::dbExecute(conn, "COMMIT")
  committed <- TRUE

  pop$tables <- unique(c(pop$tables, "founder_haplotypes"))

  if (is_new_table) {
    register_schema_meta(pop$db_conn, rbind(
      .sm_tbl("founder_haplotypes",
              "Pool of founder haplotypes in long format (one row per haplotype x locus) sampled by add_founders() to assign phased alleles."),
      .sm_col("founder_haplotypes", "line_name",
              "Founder line label matching add_founders() line_name. NULL = shared pool for all lines."),
      .sm_col("founder_haplotypes", "haplotype_id",
              "Sequential haplotype identifier within the pool (unique per line_name)"),
      .sm_col("founder_haplotypes", "locus_name",
              "Locus name; FK to genome_meta.locus_name"),
      .sm_col("founder_haplotypes", "allele",
              "Allele on this haplotype at this locus: 0 or 1")
    ))
  }

  pop
}
