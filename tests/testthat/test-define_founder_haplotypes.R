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


test_that("method = 'fixed' realizes the target frequency exactly at every locus", {
  set.seed(11)
  pop <- make_fh_pop("fh_fixed_exact", n_loci = 200) |>
    define_founder_haplotypes(
      n_haplotypes = 100,
      method       = "fixed",
      allele_freq  = 0.3
    )

  m <- fh_wide(pop)
  expect_equal(nrow(m), 100L)
  # Exact allocation: every column carries exactly 30 one-alleles.
  expect_true(all(colSums(m) == 30L))
  expect_true(all(colMeans(m) == 0.3))

  # Independent per-locus permutations, so loci are not carbon copies of
  # each other (that would be perfect LD).
  expect_gt(length(unique(apply(m, 2, paste0, collapse = ""))), 1L)

  close_pop(pop)
})


test_that("method = 'fixed' handles the 0 and 1 boundaries exactly", {
  pop0 <- make_fh_pop("fh_fixed_0", n_loci = 20) |>
    define_founder_haplotypes(n_haplotypes = 10, method = "fixed",
                              allele_freq = 0)
  expect_true(all(fh_wide(pop0) == 0L))
  close_pop(pop0)

  pop1 <- make_fh_pop("fh_fixed_1", n_loci = 20) |>
    define_founder_haplotypes(n_haplotypes = 10, method = "fixed",
                              allele_freq = 1)
  expect_true(all(fh_wide(pop1) == 1L))
  close_pop(pop1)
})


test_that("method = 'fixed' warns when allele_freq is off the 1/n_haplotypes grid", {
  pop <- make_fh_pop("fh_fixed_grid", n_loci = 20)
  expect_warning(
    pop <- define_founder_haplotypes(pop, n_haplotypes = 7, method = "fixed",
                                     allele_freq = 0.5),
    regexp = "not exactly representable"
  )

  m <- fh_wide(pop)
  # round(0.5 * 7) = 4 (R rounds half to even at .5 -> 4 here)
  expect_true(all(colSums(m) == 4L))
  gm <- get_table(pop, "genome_meta") |> dplyr::collect()
  expect_true(all(gm$founder_allele_freq == 4 / 7))

  close_pop(pop)
})


test_that("method = 'fixed' defaults allele_freq to 0.5", {
  pop <- make_fh_pop("fh_fixed_default", n_loci = 20) |>
    define_founder_haplotypes(n_haplotypes = 10, method = "fixed")

  gm <- get_table(pop, "genome_meta") |> dplyr::collect()
  expect_true(all(gm$founder_allele_freq == 0.5))
  expect_true(all(colSums(fh_wide(pop)) == 5L))

  close_pop(pop)
})


test_that("every method reproduces byte-identically from the same seed", {
  methods <- list(
    fixed           = list(method = "fixed"),
    uniform         = list(method = "uniform"),
    beta            = list(method = "beta"),
    balding_nichols = list(method = "balding_nichols"),
    mosaic          = list(method = "mosaic"),
    gaussian_copula = list(method = "gaussian_copula")
  )

  run_one <- function(nm, args) {
    set.seed(4242)
    # suppressWarnings: mosaic legitimately warns about monomorphic loci at a
    # small n_templates (see its own test below). Irrelevant to reproducibility.
    pop <- suppressWarnings(do.call(
      define_founder_haplotypes,
      c(list(make_fh_pop(paste0("fh_seed_", nm), n_loci = 60, n_chr = 3),
             n_haplotypes = 40), args)
    ))
    out <- list(
      hap  = fh_wide(pop),
      freq = get_table(pop, "genome_meta") |> dplyr::collect() |>
               dplyr::arrange(locus_id) |> dplyr::pull(founder_allele_freq),
      # The RNG stream position after the call must also be seed-determined,
      # or downstream add_founders()/add_offspring() would not reproduce.
      tail = stats::runif(3)
    )
    close_pop(pop)
    out
  }

  for (nm in names(methods)) {
    a <- run_one(nm, methods[[nm]])
    b <- run_one(nm, methods[[nm]])
    expect_identical(a, b, info = nm)
  }
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
                              fst = 0.1, mean_allele_freq = 0.5)

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


test_that("method = 'balding_nichols' validates mean_allele_freq boundaries", {
  pop0 <- make_fh_pop("fh_bn_mf")
  on.exit(close_pop(pop0))

  expect_error(
    define_founder_haplotypes(pop0, n_haplotypes = 10,
                              method = "balding_nichols", mean_allele_freq = 0),
    "mean_allele_freq"
  )
  expect_error(
    define_founder_haplotypes(pop0, n_haplotypes = 10,
                              method = "balding_nichols", mean_allele_freq = 1),
    "mean_allele_freq"
  )
})


# ============================================================
# method = "mosaic"
# ============================================================

test_that("method = 'mosaic' generates valid haplotypes", {
  set.seed(3)
  # suppressWarnings: the monomorphic-loci warning is inherent to mosaic at a
  # small n_templates and has its own test below.
  pop <- suppressWarnings(
    make_fh_pop("fh_mosaic") |>
      define_founder_haplotypes(n_haplotypes = 40, method = "mosaic"))

  expect_equal(fh_n_haps(pop), 40L)
  alleles <- fh_wide(pop)
  expect_equal(nrow(alleles), 40L)
  expect_equal(ncol(alleles), 100L)
  expect_true(all(alleles %in% c(0, 1)))

  close_pop(pop)
})


test_that("method = 'mosaic' stores empirical colMeans as founder_allele_freq", {
  set.seed(5)
  pop <- suppressWarnings(
    make_fh_pop("fh_mosaic_freq") |>
      define_founder_haplotypes(n_haplotypes = 50, method = "mosaic"))

  empirical <- colMeans(fh_wide(pop))

  gm <- get_table(pop, "genome_meta") |> dplyr::collect()
  expect_equal(gm$founder_allele_freq, unname(empirical), tolerance = 1e-9)

  close_pop(pop)
})


test_that("method = 'mosaic' with template_switch_rate = 0 and n_chr = 1 produces at most n_templates unique haplotypes", {
  set.seed(99)
  n_templates_used <- 2L
  pop <- suppressWarnings(
    open_pop(pop_name = "fh_mosaic_sr0", db_name = ":memory:") |>
      define_genome(n_loci = 20, n_chr = 1, chr_len_Mb = 50) |>
      define_founder_haplotypes(
        n_haplotypes = 30,
        method       = "mosaic",
        n_templates  = n_templates_used,
        template_switch_rate = 0
      ))

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


test_that("method = 'gaussian_copula' with very high ld_decay_rate gives near-independent loci", {
  set.seed(17)
  pop <- open_pop(pop_name = "fh_gc_ind", db_name = ":memory:") |>
    define_genome(n_loci = 50, n_chr = 1, chr_len_Mb = 100) |>
    define_founder_haplotypes(n_haplotypes = 2000, method = "gaussian_copula",
                              ld_decay_rate = 1000)

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


test_that("passing ld_decay_rate with method = 'fixed' errors with informative message", {
  pop0 <- make_fh_pop("fh_wrong3")
  on.exit(close_pop(pop0))

  expect_error(
    define_founder_haplotypes(pop0, n_haplotypes = 10,
                              method = "fixed", ld_decay_rate = 1.0),
    regexp = "ld_decay_rate.*belongs to method.*gaussian_copula"
  )
})


test_that("passing mosaic arg with method = 'gaussian_copula' errors", {
  pop0 <- make_fh_pop("fh_wrong4")
  on.exit(close_pop(pop0))

  expect_error(
    define_founder_haplotypes(pop0, n_haplotypes = 10,
                              method = "gaussian_copula", template_switch_rate = 2.0),
    regexp = "template_switch_rate.*belongs to method.*mosaic"
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


# ============================================================
# Strict argument validation
# ============================================================

test_that("n_haplotypes rejects fractional, sub-1, over-range, and non-numeric values", {
  pop0 <- make_fh_pop("fh_nhap")
  on.exit(close_pop(pop0))

  # 0 < n < 1 used to coerce to 0L and fail deep inside the write path with
  # "invalid '(to - from)/by'".
  expect_error(define_founder_haplotypes(pop0, n_haplotypes = 0.5),
               "must be at least 1")
  expect_error(define_founder_haplotypes(pop0, n_haplotypes = 2.7),
               "must be a whole number")
  # Used to overflow silently to NA and die at matrix(nrow = NA).
  expect_error(define_founder_haplotypes(pop0, n_haplotypes = 1e10),
               "too large to represent as an integer")
  expect_error(define_founder_haplotypes(pop0, n_haplotypes = "10"),
               "must be numeric")
  expect_error(define_founder_haplotypes(pop0, n_haplotypes = NA_real_),
               "must not be NA")
})


test_that("n_templates validates type before coercion and is bounded by n_haplotypes", {
  pop0 <- make_fh_pop("fh_ntmpl")
  on.exit(close_pop(pop0))

  # as.integer() used to run BEFORE is.numeric(), making the type check
  # tautological so a string sailed through.
  expect_error(
    define_founder_haplotypes(pop0, n_haplotypes = 10, method = "mosaic",
                              n_templates = "5"),
    "must be numeric"
  )
  expect_error(
    define_founder_haplotypes(pop0, n_haplotypes = 10, method = "mosaic",
                              n_templates = 2.7),
    "must be a whole number"
  )
  expect_error(
    define_founder_haplotypes(pop0, n_haplotypes = 10, method = "mosaic",
                              n_templates = 50),
    "must not exceed"
  )
})


test_that("non-finite shape and rate parameters are rejected", {
  pop0 <- make_fh_pop("fh_inf")
  on.exit(close_pop(pop0))

  # rbeta(n, Inf, Inf) silently returns exactly 0.5 at every locus.
  expect_error(
    define_founder_haplotypes(pop0, n_haplotypes = 10, method = "beta",
                              beta_shape1 = Inf),
    "must be finite"
  )
  expect_error(
    define_founder_haplotypes(pop0, n_haplotypes = 10, method = "mosaic",
                              template_switch_rate = Inf),
    "must be finite"
  )
  # Inf here used to write silent NA alleles when two loci shared a pos_cM.
  expect_error(
    define_founder_haplotypes(pop0, n_haplotypes = 10, method = "gaussian_copula",
                              ld_decay_rate = Inf),
    "must be finite"
  )
})


test_that("ld_decay_rate = 0 is allowed and gives complete within-chromosome LD", {
  set.seed(21)
  # Mirrors the template_switch_rate = 0 test: both are the same maximum-LD
  # limit, but ld_decay_rate = 0 used to be rejected while the mosaic
  # equivalent, template_switch_rate = 0, was accepted.
  pop <- suppressWarnings(
    open_pop(pop_name = "fh_decay0", db_name = ":memory:") |>
      define_genome(n_loci = 20, n_chr = 1, chr_len_Mb = 50) |>
      define_founder_haplotypes(n_haplotypes = 30, method = "gaussian_copula",
                                ld_decay_rate = 0))

  m <- fh_wide(pop)

  # rho = 1 throughout, so a haplotype's latent z is constant along the whole
  # chromosome and its allele at locus j is just (z < threshold_j). Every row is
  # therefore determined by one scalar, which makes the alleles perfectly
  # *nested* rather than perfectly correlated: sorting loci by threshold, a 1 at
  # one locus implies a 1 at every looser locus. Two consequences to assert:
  #   (a) at most n_loci + 1 distinct rows are reachable (one per z interval),
  #       versus 2^n_loci under independence;
  #   (b) no pair of loci shows all four 0/1 combinations (no "recombinants").
  expect_lte(nrow(unique(m)), ncol(m) + 1L)

  poly <- which(apply(m, 2, function(x) length(unique(x)) > 1L))
  if (length(poly) >= 2L) {
    combos <- nrow(unique(m[, poly[1:2], drop = FALSE]))
    expect_lte(combos, 3L)
  }

  close_pop(pop)
})


test_that("mosaic warns when a small n_templates leaves many loci monomorphic", {
  set.seed(31)
  expect_warning(
    pop <- open_pop(pop_name = "fh_mono", db_name = ":memory:") |>
      define_genome(n_loci = 200, n_chr = 2, chr_len_Mb = 50) |>
      define_founder_haplotypes(n_haplotypes = 50, method = "mosaic",
                                n_templates = 2),
    "monomorphic"
  )

  # Templates are the only source of variation, so ~2/(K+1) of loci are
  # monomorphic no matter how many haplotypes are drawn.
  f <- colMeans(fh_wide(pop))
  expect_gt(mean(f == 0 | f == 1), 0.4)

  close_pop(pop)
})


# ============================================================
# exact_freq
# ============================================================

test_that("exact_freq = TRUE realizes drawn frequencies exactly for distribution methods", {
  for (m in c("uniform", "beta", "balding_nichols")) {
    set.seed(77)
    pop <- make_fh_pop(paste0("fh_exact_", m), n_loci = 150) |>
      define_founder_haplotypes(n_haplotypes = 50, method = m,
                                exact_freq = TRUE)

    realised <- colMeans(fh_wide(pop))
    stored <- get_table(pop, "genome_meta") |> dplyr::collect() |>
      dplyr::arrange(locus_id) |> dplyr::pull(founder_allele_freq)

    # Every locus sits exactly on the 1/n_haplotypes grid...
    expect_equal(realised * 50, round(realised * 50), tolerance = 1e-12,
                 info = m)
    # ...and founder_allele_freq is the realised value, not the pre-rounding target.
    expect_equal(stored, unname(realised), tolerance = 1e-12, info = m)

    close_pop(pop)
  }
})


test_that("exact_freq defaults to FALSE for distribution methods, TRUE for fixed", {
  set.seed(5)
  pop_b <- make_fh_pop("fh_exact_default_b", n_loci = 300) |>
    define_founder_haplotypes(n_haplotypes = 50, method = "uniform")
  realised_b <- colMeans(fh_wide(pop_b))
  stored_b <- get_table(pop_b, "genome_meta") |> dplyr::collect() |>
    dplyr::arrange(locus_id) |> dplyr::pull(founder_allele_freq)
  # Binomial sampling: realised scatters away from the stored target.
  expect_gt(max(abs(realised_b - stored_b)), 0)
  close_pop(pop_b)

  set.seed(5)
  pop_f <- make_fh_pop("fh_exact_default_f", n_loci = 100) |>
    define_founder_haplotypes(n_haplotypes = 50, method = "fixed",
                              allele_freq = 0.4)
  expect_true(all(colMeans(fh_wide(pop_f)) == 0.4))
  close_pop(pop_f)
})


test_that("exact_freq is rejected for the LD methods", {
  pop0 <- make_fh_pop("fh_exact_ld")
  on.exit(close_pop(pop0))

  for (m in c("mosaic", "gaussian_copula")) {
    expect_error(
      define_founder_haplotypes(pop0, n_haplotypes = 10, method = m,
                                exact_freq = TRUE),
      "belongs to method"
    )
  }
})


test_that("LD methods resolve the genetic map for their own line_name", {
  # A line-specific map with distances 100x the default changes which loci fall
  # in the same LD block, so the two pools must differ. Before this fix both
  # branches called resolve_genome_map() with no line_name and silently used the
  # default map even when the user named a line.
  build <- function(nm) {
    pop <- open_pop(pop_name = nm, db_name = ":memory:") |>
      define_genome(n_loci = 60, n_chr = 1, chr_len_Mb = 50)
    dflt <- DBI::dbGetQuery(pop$db_conn,
      "SELECT locus_id, locus_name, pos_cM FROM genome_map ORDER BY locus_id")
    start_id <- next_int_id(pop$db_conn, "genome_map", "id_genome_map")
    line_map <- data.frame(
      id_genome_map = seq.int(start_id, length.out = nrow(dflt)),
      locus_id      = dflt$locus_id,
      locus_name    = dflt$locus_name,
      sex           = NA_character_,
      line_name     = "A",
      map_name      = "default",
      pos_cM        = dflt$pos_cM * 100
    )
    DBI::dbWriteTable(pop$db_conn, "genome_map", line_map, append = TRUE)
    pop
  }

  set.seed(8)
  pop_a <- suppressWarnings(
    build("fh_map_a") |>
      define_founder_haplotypes(n_haplotypes = 30, method = "mosaic",
                                n_templates = 4, line_name = "A"))
  hap_a <- fh_wide(pop_a)
  close_pop(pop_a)

  set.seed(8)
  pop_b <- suppressWarnings(
    build("fh_map_b") |>
      define_founder_haplotypes(n_haplotypes = 30, method = "mosaic",
                                n_templates = 4, line_name = "B"))
  hap_b <- fh_wide(pop_b)
  close_pop(pop_b)

  # Same seed, same parameters -- the only difference is which map was resolved.
  expect_false(identical(hap_a, hap_b))
})
