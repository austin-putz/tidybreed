make_effects_pop <- function(pop_name = "eff", n_ind = 500, n_loci = 500) {
  pop <- initialize_genome(
    pop_name        = pop_name,
    n_loci          = n_loci,
    n_chr           = 5,
    chr_len_Mb      = 100,
    n_haplotypes    = 200,
    db_path         = ":memory:",
    fixed_allele_freq = 0.5
  )
  pop <- add_founders(pop, n_males = n_ind / 2, n_females = n_ind / 2,
                      line_name = "A")
  pop
}


test_that("set_qtl_effects() rescales to target_add_var within tolerance", {
  set.seed(42)
  pop <- make_effects_pop("eff_scale", n_ind = 600, n_loci = 600)

  pop <- add_trait(pop, "ADG", target_add_var = 0.5, residual_var = 0.5)
  pop <- define_qtl(pop, "ADG", n = 100, method = "random")
  pop <- set_qtl_effects(pop, "ADG", distribution = "normal", seed = 1)

  genome <- dplyr::collect(get_table(pop, "genome_meta"))
  a <- genome$add_ADG
  a[is.na(a)] <- 0

  geno <- DBI::dbGetQuery(pop$db_conn, "SELECT * FROM genome_genotype")
  locus_cols <- grep("^locus_", names(geno), value = TRUE)
  locus_cols <- locus_cols[order(as.integer(sub("^locus_", "", locus_cols)))]
  X <- as.matrix(geno[, locus_cols])
  realised <- var(as.numeric(X %*% a))

  expect_equal(realised, 0.5, tolerance = 0.05)
  close_pop(pop)
})


test_that("TBV mean is approximately 0 for founder population", {
  set.seed(7)
  pop <- make_effects_pop("eff_mean", n_ind = 500, n_loci = 500)

  pop <- add_trait(pop, "ADG", target_add_var = 100, residual_var = 100,
                   target_add_mean = 0)
  pop <- define_qtl(pop, "ADG", n = 200, method = "random")
  pop <- set_qtl_effects(pop, "ADG", distribution = "normal",
                         base = "founder_haplotypes", seed = 3)

  pop <- pop |> get_table("ind_meta") |> add_tbv("ADG")

  tbv_df <- dplyr::collect(get_table(pop, "ind_tbv"))
  expect_equal(nrow(tbv_df), 500)
  # mean TBV should be close to 0 (within ~2 SE = 2*sqrt(100/500) ≈ 0.9)
  expect_equal(mean(tbv_df$tbv), 0, tolerance = 2.0)
  # var TBV should be close to target_add_var = 100
  expect_equal(var(tbv_df$tbv), 100, tolerance = 15)

  close_pop(pop)
})


test_that("base_allele_freq_{trait} column is written to genome_meta", {
  pop <- make_effects_pop("eff_base_col")

  pop <- add_trait(pop, "ADG", target_add_var = 1, residual_var = 1)
  pop <- define_qtl(pop, "ADG", n = 50, method = "random")
  pop <- set_qtl_effects(pop, "ADG", distribution = "normal")

  genome <- dplyr::collect(get_table(pop, "genome_meta"))
  expect_true("base_allele_freq_ADG" %in% names(genome))
  # all values should be between 0 and 1
  freqs <- genome$base_allele_freq_ADG
  expect_true(all(freqs >= 0 & freqs <= 1))

  close_pop(pop)
})


test_that("base = 'current_pop' via tidybreed_table pipe works", {
  set.seed(17)
  pop <- make_effects_pop("eff_currpop", n_ind = 200, n_loci = 300)

  pop <- get_table(pop, "ind_meta") |> mutate_table(gen = 0L)
  pop <- add_trait(pop, "ADG", target_add_var = 50, residual_var = 50)
  pop <- define_qtl(pop, "ADG", n = 100, method = "random")

  # Pass tidybreed_table as first arg with base = "current_pop"
  pop <- pop |>
    get_table("ind_meta") |>
    dplyr::filter(gen == 0L) |>
    set_qtl_effects("ADG", base = "current_pop", distribution = "normal", seed = 5)

  genome <- dplyr::collect(get_table(pop, "genome_meta"))
  expect_true("add_ADG" %in% names(genome))
  expect_true("base_allele_freq_ADG" %in% names(genome))

  # TBV mean should be ≈ 0
  pop <- pop |> get_table("ind_meta") |> add_tbv("ADG")
  tbv_df <- dplyr::collect(get_table(pop, "ind_tbv"))
  expect_equal(mean(tbv_df$tbv), 0, tolerance = 3.0)

  close_pop(pop)
})


test_that("set_qtl_effects() accepts manual effects", {
  pop <- make_effects_pop("eff_manual")
  pop <- add_trait(pop, "ADG", target_add_var = 1, residual_var = 1)
  pop <- define_qtl(pop, "ADG", n = 10, method = "even")

  # 10 QTL, all effects = 2.0
  pop <- set_qtl_effects(pop, "ADG", effects = rep(2.0, 10))
  genome <- dplyr::collect(get_table(pop, "genome_meta"))
  expect_equal(sum(!is.na(genome$add_ADG)), 10)
  expect_true(all(genome$add_ADG[!is.na(genome$add_ADG)] == 2.0))

  close_pop(pop)
})


test_that("set_qtl_effects_multi() hits target variances per trait", {
  set.seed(123)
  pop <- make_effects_pop("eff_multi", n_ind = 800, n_loci = 600)

  pop <- add_trait(pop, "ADG", target_add_var = 0.25, residual_var = 0.75)
  pop <- add_trait(pop, "BW",  target_add_var = 0.50, residual_var = 0.50)

  # Same QTL for both traits (full pleiotropy)
  pop <- define_qtl(pop, "ADG", n = 150, method = "random")
  qtl_tf <- dplyr::collect(get_table(pop, "genome_meta"))$is_QTL_ADG
  pop <- define_qtl(pop, "BW", locus_tf = qtl_tf)

  G <- matrix(c(0.25, 0.10, 0.10, 0.50), 2, 2)
  pop <- set_qtl_effects_multi(pop, c("ADG", "BW"), G = G, seed = 7)

  genome <- dplyr::collect(get_table(pop, "genome_meta"))
  aA <- genome$add_ADG; aA[is.na(aA)] <- 0
  aB <- genome$add_BW;  aB[is.na(aB)] <- 0

  geno <- DBI::dbGetQuery(pop$db_conn, "SELECT * FROM genome_genotype")
  locus_cols <- grep("^locus_", names(geno), value = TRUE)
  locus_cols <- locus_cols[order(as.integer(sub("^locus_", "", locus_cols)))]
  X <- as.matrix(geno[, locus_cols])

  bv_A <- as.numeric(X %*% aA)
  bv_B <- as.numeric(X %*% aB)

  expect_equal(var(bv_A), 0.25, tolerance = 0.06)
  expect_equal(var(bv_B), 0.50, tolerance = 0.10)

  # Realised correlation should be positive (target ~0.28)
  expect_gt(cor(bv_A, bv_B), 0.05)

  close_pop(pop)
})
