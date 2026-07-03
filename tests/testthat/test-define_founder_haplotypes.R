# Helper: minimal pop (no founders needed)
make_fh_pop <- function(pop_name = "fh_test", n_loci = 100, n_chr = 2) {
  open_pop(pop_name = pop_name, db_name = ":memory:") |>
    define_genome(n_loci = n_loci, n_chr = n_chr, chr_len_Mb = 50)
}

# founder_haplotypes is stored long. Reconstruct a haplotype x locus matrix
# (rows = one per (line_name, haplotype_id), cols in locus_id order).
fh_wide <- function(pop) {
  long <- get_table(pop, "founder_haplotypes") |> dplyr::collect()
  gm <- DBI::dbGetQuery(pop$db_conn,
    "SELECT locus_id, locus_name FROM genome_meta ORDER BY locus_id")
  long$locus_id <- gm$locus_id[match(long$locus_name, gm$locus_name)]
  key_df <- unique(long[, c("line_name", "haplotype_id")])
  key_df <- key_df[order(key_df$line_name, key_df$haplotype_id), , drop = FALSE]
  keys <- paste(key_df$line_name, key_df$haplotype_id, sep = "\r")
  rk   <- paste(long$line_name, long$haplotype_id, sep = "\r")
  m <- matrix(0L, nrow = length(keys), ncol = nrow(gm))
  m[cbind(match(rk, keys), long$locus_id)] <- as.integer(long$allele)
  m
}

# Number of distinct haplotypes in the pool.
fh_n_haps <- function(pop) {
  DBI::dbGetQuery(pop$db_conn,
    "SELECT COUNT(*) AS n FROM (SELECT DISTINCT line_name, haplotype_id FROM founder_haplotypes)")$n
}

# ============================================================
# Existing methods (updated to new API)
# ============================================================

test_that("method = 'uniform' creates haplotypes with per-locus frequencies in range", {
  pop <- make_fh_pop("fh_uniform") |>
    define_founder_haplotypes(
      n_haplotypes    = 50,
      method          = "uniform",
      min_allele_freq = 0.1,
      max_allele_freq = 0.9
    )

  expect_true("founder_haplotypes" %in% pop$tables)

  fh <- get_table(pop, "founder_haplotypes") |> dplyr::collect()
  expect_equal(fh_n_haps(pop), 50L)
  expect_setequal(unique(fh$haplotype_id), 1:50)
  expect_true("line_name" %in% colnames(fh))
  expect_true(all(is.na(fh$line_name)))

  alleles <- fh_wide(pop)
  expect_equal(nrow(alleles), 50L)
  expect_equal(ncol(alleles), 100L)
  expect_true(all(alleles %in% c(0, 1)))

  gm <- get_table(pop, "genome_meta") |> dplyr::collect()
  expect_true("founder_allele_freq" %in% colnames(gm))
  expect_equal(length(gm$founder_allele_freq), 100L)
  expect_true(all(gm$founder_allele_freq >= 0.1))
  expect_true(all(gm$founder_allele_freq <= 0.9))

  close_pop(pop)
})


test_that("method = 'fixed' creates haplotypes with uniform allele frequency", {
  pop <- make_fh_pop("fh_fixed", n_loci = 50) |>
    define_founder_haplotypes(
      n_haplotypes = 30,
      method       = "fixed",
      allele_freq  = 0.5
    )

  expect_true("founder_haplotypes" %in% pop$tables)

  gm <- get_table(pop, "genome_meta") |> dplyr::collect()
  expect_true(all(gm$founder_allele_freq == 0.5))

  expect_equal(fh_n_haps(pop), 30L)

  close_pop(pop)
})


test_that("method = 'fixed' defaults allele_freq to 0.5", {
  pop <- make_fh_pop("fh_fixed_default", n_loci = 20) |>
    define_founder_haplotypes(n_haplotypes = 10, method = "fixed")

  gm <- get_table(pop, "genome_meta") |> dplyr::collect()
  expect_true(all(gm$founder_allele_freq == 0.5))

  close_pop(pop)
})


# ============================================================
# method = "beta"
# ============================================================

test_that("method = 'beta' generates valid haplotypes with correct dimensions", {
  set.seed(1)
  pop <- make_fh_pop("fh_beta") |>
    define_founder_haplotypes(n_haplotypes = 40, method = "beta")

  expect_equal(fh_n_haps(pop), 40L)
  alleles <- fh_wide(pop)
  expect_equal(nrow(alleles), 40L)
  expect_equal(ncol(alleles), 100L)
  expect_true(all(alleles %in% c(0, 1)))

  gm <- get_table(pop, "genome_meta") |> dplyr::collect()
  expect_true("founder_allele_freq" %in% colnames(gm))
  expect_true(all(gm$founder_allele_freq >= 0))
  expect_true(all(gm$founder_allele_freq <= 1))

  close_pop(pop)
})


test_that("method = 'beta' with shape1 = shape2 = 2 produces frequencies near 0.5", {
  set.seed(42)
  pop <- make_fh_pop("fh_beta_sym", n_loci = 200) |>
    define_founder_haplotypes(n_haplotypes = 50, method = "beta",
                              beta_shape1 = 2, beta_shape2 = 2)

  gm <- get_table(pop, "genome_meta") |> dplyr::collect()
  # Beta(2,2) is symmetric around 0.5; most freqs between 0.1 and 0.9
  expect_true(mean(gm$founder_allele_freq > 0.1 & gm$founder_allele_freq < 0.9) > 0.8)

  close_pop(pop)
})


# ============================================================
# method = "balding_nichols"
# ============================================================

test_that("method = 'balding_nichols' generates valid haplotypes", {
  set.seed(7)
  pop <- make_fh_pop("fh_bn") |>
    define_founder_haplotypes(n_haplotypes = 40, method = "balding_nichols",
                              fst = 0.1, mean_freq = 0.5)

  expect_equal(fh_n_haps(pop), 40L)
  alleles <- fh_wide(pop)
  expect_equal(nrow(alleles), 40L)
  expect_equal(ncol(alleles), 100L)
  expect_true(all(alleles %in% c(0, 1)))

  gm <- get_table(pop, "genome_meta") |> dplyr::collect()
  expect_true(all(gm$founder_allele_freq >= 0))
  expect_true(all(gm$founder_allele_freq <= 1))

  close_pop(pop)
})


test_that("method = 'balding_nichols' validates fst boundaries", {
  pop0 <- make_fh_pop("fh_bn_val")
  on.exit(close_pop(pop0))

  expect_error(
    define_founder_haplotypes(pop0, n_haplotypes = 10,
                              method = "balding_nichols", fst = 0),
    "fst"
  )
  expect_error(
    define_founder_haplotypes(pop0, n_haplotypes = 10,
                              method = "balding_nichols", fst = 1),
    "fst"
  )
})


test_that("method = 'balding_nichols' validates mean_freq boundaries", {
  pop0 <- make_fh_pop("fh_bn_mf")
  on.exit(close_pop(pop0))

  expect_error(
    define_founder_haplotypes(pop0, n_haplotypes = 10,
                              method = "balding_nichols", mean_freq = 0),
    "mean_freq"
  )
  expect_error(
    define_founder_haplotypes(pop0, n_haplotypes = 10,
                              method = "balding_nichols", mean_freq = 1),
    "mean_freq"
  )
})


# ============================================================
# method = "mosaic"
# ============================================================

test_that("method = 'mosaic' generates valid haplotypes", {
  set.seed(3)
  pop <- make_fh_pop("fh_mosaic") |>
    define_founder_haplotypes(n_haplotypes = 40, method = "mosaic")

  expect_equal(fh_n_haps(pop), 40L)
  alleles <- fh_wide(pop)
  expect_equal(nrow(alleles), 40L)
  expect_equal(ncol(alleles), 100L)
  expect_true(all(alleles %in% c(0, 1)))

  close_pop(pop)
})


test_that("method = 'mosaic' stores empirical colMeans as founder_allele_freq", {
  set.seed(5)
  pop <- make_fh_pop("fh_mosaic_freq") |>
    define_founder_haplotypes(n_haplotypes = 50, method = "mosaic")

  empirical <- colMeans(fh_wide(pop))

  gm <- get_table(pop, "genome_meta") |> dplyr::collect()
  expect_equal(gm$founder_allele_freq, unname(empirical), tolerance = 1e-9)

  close_pop(pop)
})


test_that("method = 'mosaic' with switch_rate = 0 and n_chr = 1 produces at most n_templates unique haplotypes", {
  set.seed(99)
  n_templates_used <- 2L
  pop <- open_pop(pop_name = "fh_mosaic_sr0", db_name = ":memory:") |>
    define_genome(n_loci = 20, n_chr = 1, chr_len_Mb = 50) |>
    define_founder_haplotypes(
      n_haplotypes = 30,
      method       = "mosaic",
      n_templates  = n_templates_used,
      switch_rate  = 0
    )

  n_unique <- nrow(unique(fh_wide(pop)))

  # With no switching on a single chromosome, each haplotype is an exact
  # copy of its starting template, so unique count <= n_templates
  expect_lte(n_unique, n_templates_used)

  close_pop(pop)
})


# ============================================================
# method = "gaussian_copula"
# ============================================================

test_that("method = 'gaussian_copula' generates valid haplotypes", {
  set.seed(11)
  pop <- make_fh_pop("fh_gc") |>
    define_founder_haplotypes(n_haplotypes = 40, method = "gaussian_copula")

  expect_equal(fh_n_haps(pop), 40L)
  alleles <- fh_wide(pop)
  expect_equal(nrow(alleles), 40L)
  expect_equal(ncol(alleles), 100L)
  expect_true(all(alleles %in% c(0, 1)))

  close_pop(pop)
})


test_that("method = 'gaussian_copula' stores empirical colMeans as founder_allele_freq", {
  set.seed(13)
  pop <- make_fh_pop("fh_gc_freq") |>
    define_founder_haplotypes(n_haplotypes = 50, method = "gaussian_copula")

  empirical <- colMeans(fh_wide(pop))

  gm <- get_table(pop, "genome_meta") |> dplyr::collect()
  expect_equal(gm$founder_allele_freq, unname(empirical), tolerance = 1e-9)

  close_pop(pop)
})


test_that("method = 'gaussian_copula' with very high decay_rate gives near-independent loci", {
  set.seed(17)
  pop <- open_pop(pop_name = "fh_gc_ind", db_name = ":memory:") |>
    define_genome(n_loci = 50, n_chr = 1, chr_len_Mb = 100) |>
    define_founder_haplotypes(n_haplotypes = 2000, method = "gaussian_copula",
                              decay_rate = 1000)

  mat <- fh_wide(pop)

  # With rho ≈ 0, correlations between adjacent loci should be near zero
  r_adj <- cor(mat[, 1:(ncol(mat) - 1)], mat[, 2:ncol(mat)])
  # Take diagonal (adjacent-locus correlations)
  adj_cors <- diag(r_adj)
  expect_true(mean(abs(adj_cors)) < 0.15)

  close_pop(pop)
})


# ============================================================
# Wrong-method argument validation
# ============================================================

test_that("passing allele_freq with method = 'uniform' errors with informative message", {
  pop0 <- make_fh_pop("fh_wrong1")
  on.exit(close_pop(pop0))

  expect_error(
    define_founder_haplotypes(pop0, n_haplotypes = 10,
                              method = "uniform", allele_freq = 0.5),
    regexp = "allele_freq.*belongs to method.*fixed"
  )
})


test_that("passing fst with method = 'beta' errors with informative message", {
  pop0 <- make_fh_pop("fh_wrong2")
  on.exit(close_pop(pop0))

  expect_error(
    define_founder_haplotypes(pop0, n_haplotypes = 10,
                              method = "beta", fst = 0.1),
    regexp = "fst.*belongs to method.*balding_nichols"
  )
})


test_that("passing decay_rate with method = 'fixed' errors with informative message", {
  pop0 <- make_fh_pop("fh_wrong3")
  on.exit(close_pop(pop0))

  expect_error(
    define_founder_haplotypes(pop0, n_haplotypes = 10,
                              method = "fixed", decay_rate = 1.0),
    regexp = "decay_rate.*belongs to method.*gaussian_copula"
  )
})


test_that("passing mosaic arg with method = 'gaussian_copula' errors", {
  pop0 <- make_fh_pop("fh_wrong4")
  on.exit(close_pop(pop0))

  expect_error(
    define_founder_haplotypes(pop0, n_haplotypes = 10,
                              method = "gaussian_copula", switch_rate = 2.0),
    regexp = "switch_rate.*belongs to method.*mosaic"
  )
})


test_that("invalid method string raises match.arg error", {
  pop0 <- make_fh_pop("fh_bad_method")
  on.exit(close_pop(pop0))

  expect_error(
    define_founder_haplotypes(pop0, n_haplotypes = 10, method = "invalid_method"),
    regexp = "should be one of"
  )
})


# ============================================================
# Existing guard tests (unchanged)
# ============================================================

test_that("open_pop() + define_genome() alone does not create founder_haplotypes", {
  pop <- make_fh_pop("fh_no_founders")
  on.exit(close_pop(pop))

  expect_false("founder_haplotypes" %in% pop$tables)
  gm <- get_table(pop, "genome_meta") |> dplyr::collect()
  expect_false("founder_allele_freq" %in% colnames(gm))
})


test_that("define_founder_haplotypes() errors if called twice with same (NULL) line_name", {
  pop <- make_fh_pop("fh_twice", n_loci = 10) |>
    define_founder_haplotypes(n_haplotypes = 10)

  expect_error(
    define_founder_haplotypes(pop, n_haplotypes = 10),
    "line_name = NULL"
  )

  close_pop(pop)
})


test_that("define_founder_haplotypes() errors if called twice with same named line_name", {
  pop <- make_fh_pop("fh_twice_named", n_loci = 10) |>
    define_founder_haplotypes(n_haplotypes = 10, line_name = "A")

  expect_error(
    define_founder_haplotypes(pop, n_haplotypes = 10, line_name = "A"),
    "line 'A'"
  )

  close_pop(pop)
})


test_that("define_founder_haplotypes() succeeds for two distinct line_names", {
  pop <- make_fh_pop("fh_two_lines", n_loci = 20) |>
    define_founder_haplotypes(n_haplotypes = 30, line_name = "A") |>
    define_founder_haplotypes(n_haplotypes = 40, line_name = "B")

  expect_equal(fh_n_haps(pop), 70L)
  haps <- DBI::dbGetQuery(pop$db_conn,
    "SELECT DISTINCT line_name, haplotype_id FROM founder_haplotypes")
  expect_equal(sort(unique(haps$line_name)), c("A", "B"))
  expect_equal(sum(haps$line_name == "A"), 30L)
  expect_equal(sum(haps$line_name == "B"), 40L)
  # haplotype_id restarts per line (1..n).
  expect_setequal(haps$haplotype_id[haps$line_name == "A"], 1:30)
  expect_setequal(haps$haplotype_id[haps$line_name == "B"], 1:40)

  close_pop(pop)
})


test_that("define_founder_haplotypes() invalid line_name format errors", {
  pop0 <- make_fh_pop("fh_bad_line")
  on.exit(close_pop(pop0))

  expect_error(
    define_founder_haplotypes(pop0, n_haplotypes = 10, line_name = "1bad"),
    "must start with a letter"
  )
  expect_error(
    define_founder_haplotypes(pop0, n_haplotypes = 10, line_name = "bad-name"),
    "must start with a letter"
  )
})


test_that("define_founder_haplotypes() is pipe-friendly and returns tidybreed_pop", {
  pop <- make_fh_pop("fh_pipe", n_loci = 10) |>
    define_founder_haplotypes(n_haplotypes = 5)

  expect_true(inherits(pop, "tidybreed_pop"))
  expect_true("founder_haplotypes" %in% pop$tables)

  close_pop(pop)
})
