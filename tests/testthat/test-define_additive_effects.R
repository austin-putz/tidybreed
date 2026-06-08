make_effects_pop <- function(pop_name = "eff", n_ind = 500, n_loci = 500) {
  pop <- initialize_genome(
    pop_name   = pop_name,
    n_loci     = n_loci,
    n_chr      = 5,
    chr_len_Mb = 100,
    db_path    = ":memory:"
  ) |>
    define_founder_haplotypes(n_haplotypes = 200, fixed_allele_freq = 0.5)
  pop <- add_founders(pop, n_males = n_ind / 2, n_females = n_ind / 2,
                      line_name = "A")
  pop
}


test_that("define_additive_effects() rescales to target_add_var within tolerance", {
  set.seed(42)
  pop <- make_effects_pop("eff_scale", n_ind = 600, n_loci = 600)

  pop <- define_trait(pop, "ADG", target_add_var = 0.5)
  sel <- pop |> get_table("genome_meta") |> dplyr::collect() |>
    dplyr::slice_sample(n = 100) |> dplyr::pull(locus_name)
  pop <- pop |>
    get_table("genome_meta") |>
    dplyr::filter(locus_name %in% sel) |>
    define_additive_effects("ADG", distribution = "normal", seed = 1)

  # Effects are now in genome_effects, not genome_meta columns
  eff <- DBI::dbGetQuery(pop$db_conn,
    "SELECT locus_name, genome_value FROM genome_effects WHERE trait_name = 'ADG'")
  locus_order <- DBI::dbGetQuery(pop$db_conn,
    "SELECT locus_id, locus_name FROM genome_meta ORDER BY locus_id")
  a <- rep(0, nrow(locus_order))
  idx <- match(eff$locus_name, locus_order$locus_name)
  a[idx] <- eff$genome_value

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

  pop <- define_trait(pop, "ADG", target_add_var = 100, target_add_mean = 0)
  sel <- pop |> get_table("genome_meta") |> dplyr::collect() |>
    dplyr::slice_sample(n = 200) |> dplyr::pull(locus_name)
  pop <- pop |>
    get_table("genome_meta") |>
    dplyr::filter(locus_name %in% sel) |>
    define_additive_effects("ADG", distribution = "normal",
                          base = "founder_haplotypes", seed = 3)

  pop <- pop |> get_table("ind_meta") |> add_tbv("ADG")

  tbv_df <- dplyr::collect(get_table(pop, "ind_tbv"))
  expect_equal(nrow(tbv_df), 500)
  # mean TBV should be close to 0 (within ~2 SE = 2*sqrt(100/500) ≈ 0.9)
  expect_equal(mean(tbv_df$tbv_value), 0, tolerance = 2.0)
  # var TBV should be close to target_add_var = 100
  expect_equal(var(tbv_df$tbv_value), 100, tolerance = 15)

  close_pop(pop)
})


test_that("base_allele_freq written to genome_effects, not genome_meta", {
  pop <- make_effects_pop("eff_base_col")

  pop <- define_trait(pop, "ADG", target_add_var = 1)
  sel <- pop |> get_table("genome_meta") |> dplyr::collect() |>
    dplyr::slice_sample(n = 50) |> dplyr::pull(locus_name)
  pop <- pop |>
    get_table("genome_meta") |>
    dplyr::filter(locus_name %in% sel) |>
    define_additive_effects("ADG", distribution = "normal")

  # base_allele_freq is in genome_effects, not genome_meta
  genome_cols <- DBI::dbListFields(pop$db_conn, "genome_meta")
  expect_false("base_allele_freq_ADG" %in% genome_cols)
  expect_false("add_ADG"              %in% genome_cols)
  expect_false("is_QTL_ADG"          %in% genome_cols)

  eff <- DBI::dbGetQuery(pop$db_conn,
    "SELECT base_allele_freq FROM genome_effects WHERE trait_name = 'ADG'")
  expect_equal(nrow(eff), 50)
  expect_true(all(eff$base_allele_freq >= 0 & eff$base_allele_freq <= 1))

  close_pop(pop)
})


test_that("base = 'current_pop' via base_tbl argument works", {
  set.seed(17)
  pop <- make_effects_pop("eff_currpop", n_ind = 200, n_loci = 300)

  pop <- get_table(pop, "ind_meta") |> mutate_table(gen = 0L)
  pop <- define_trait(pop, "ADG", target_add_var = 50)
  sel <- pop |> get_table("genome_meta") |> dplyr::collect() |>
    dplyr::slice_sample(n = 100) |> dplyr::pull(locus_name)

  gen0_tbl <- get_table(pop, "ind_meta") |> dplyr::filter(gen == 0L)
  pop <- pop |>
    get_table("genome_meta") |>
    dplyr::filter(locus_name %in% sel) |>
    define_additive_effects("ADG", base = "current_pop", base_tbl = gen0_tbl,
                          distribution = "normal", seed = 5)

  eff <- DBI::dbGetQuery(pop$db_conn,
    "SELECT locus_name, genome_value, base_allele_freq FROM genome_effects WHERE trait_name = 'ADG'")
  expect_equal(nrow(eff), 100)

  # TBV mean should be ≈ 0
  pop <- pop |> get_table("ind_meta") |> add_tbv("ADG")
  tbv_df <- dplyr::collect(get_table(pop, "ind_tbv"))
  expect_equal(mean(tbv_df$tbv_value), 0, tolerance = 3.0)

  close_pop(pop)
})


test_that("define_additive_effects() accepts manual effects", {
  pop <- make_effects_pop("eff_manual")
  pop <- define_trait(pop, "ADG", target_add_var = 1)
  sel <- pop |> get_table("genome_meta") |> dplyr::collect() |>
    dplyr::slice_sample(n = 10) |> dplyr::pull(locus_name)
  pop <- pop |>
    get_table("genome_meta") |>
    dplyr::filter(locus_name %in% sel) |>
    define_additive_effects("ADG", effects = rep(2.0, 10))

  eff <- DBI::dbGetQuery(pop$db_conn,
    "SELECT genome_value FROM genome_effects WHERE trait_name = 'ADG'")
  expect_equal(nrow(eff), 10)
  expect_true(all(eff$genome_value == 2.0))

  close_pop(pop)
})


test_that("re-calling define_additive_effects() replaces existing rows", {
  pop <- make_effects_pop("eff_replace")
  pop <- define_trait(pop, "ADG", target_add_var = 1)
  sel <- pop |> get_table("genome_meta") |> dplyr::collect() |>
    dplyr::slice_sample(n = 20) |> dplyr::pull(locus_name)

  pop <- pop |>
    get_table("genome_meta") |>
    dplyr::filter(locus_name %in% sel) |>
    define_additive_effects("ADG", effects = rep(1.0, 20))

  n_before <- DBI::dbGetQuery(pop$db_conn,
    "SELECT COUNT(*) AS n FROM genome_effects WHERE trait_name = 'ADG'")$n
  expect_equal(n_before, 20L)

  # Call again with different loci set
  sel2 <- pop |> get_table("genome_meta") |> dplyr::collect() |>
    dplyr::slice_sample(n = 30) |> dplyr::pull(locus_name)
  pop <- pop |>
    get_table("genome_meta") |>
    dplyr::filter(locus_name %in% sel2) |>
    define_additive_effects("ADG", effects = rep(3.0, 30))

  n_after <- DBI::dbGetQuery(pop$db_conn,
    "SELECT COUNT(*) AS n FROM genome_effects WHERE trait_name = 'ADG'")$n
  expect_equal(n_after, 30L)
  eff_vals <- DBI::dbGetQuery(pop$db_conn,
    "SELECT genome_value FROM genome_effects WHERE trait_name = 'ADG'")$genome_value
  expect_true(all(eff_vals == 3.0))

  close_pop(pop)
})


test_that("define_additive_effects() hits target variances per trait (multi-trait)", {
  set.seed(123)
  pop <- make_effects_pop("eff_multi", n_ind = 800, n_loci = 600)

  pop <- define_trait(pop, "ADG", target_add_var = 0.25)
  pop <- define_trait(pop, "BW",  target_add_var = 0.50)

  # Same QTL for both traits (full pleiotropy via method = "shared")
  sel <- pop |> get_table("genome_meta") |> dplyr::collect() |>
    dplyr::slice_sample(n = 150) |> dplyr::pull(locus_name)

  G <- matrix(c(0.25, 0.10, 0.10, 0.50), 2, 2)
  pop <- pop |>
    get_table("genome_meta") |>
    dplyr::filter(locus_name %in% sel) |>
    define_additive_effects(trait_name = c("ADG", "BW"), G = G,
                             method = "shared", seed = 7)

  locus_order <- DBI::dbGetQuery(pop$db_conn,
    "SELECT locus_id, locus_name FROM genome_meta ORDER BY locus_id")
  n_loci <- nrow(locus_order)

  load_eff <- function(t) {
    e <- DBI::dbGetQuery(pop$db_conn, paste0(
      "SELECT locus_name, genome_value FROM genome_effects WHERE trait_name = '", t, "'"))
    a <- rep(0, n_loci)
    idx <- match(e$locus_name, locus_order$locus_name)
    a[idx] <- e$genome_value
    a
  }
  aA <- load_eff("ADG")
  aB <- load_eff("BW")

  geno <- DBI::dbGetQuery(pop$db_conn, "SELECT * FROM genome_genotype")
  locus_cols <- grep("^locus_", names(geno), value = TRUE)
  locus_cols <- locus_cols[order(as.integer(sub("^locus_", "", locus_cols)))]
  X <- as.matrix(geno[, locus_cols])

  bv_A <- as.numeric(X %*% aA)
  bv_B <- as.numeric(X %*% aB)

  expect_equal(var(bv_A), 0.25, tolerance = 0.06)
  expect_equal(var(bv_B), 0.50, tolerance = 0.10)
  expect_gt(cor(bv_A, bv_B), 0.05)

  close_pop(pop)
})


test_that("define_additive_effects() errors on bare tidybreed_pop", {
  pop <- make_effects_pop("eff_err_pop")
  pop <- define_trait(pop, "ADG", target_add_var = 1)

  expect_error(
    define_additive_effects(pop, "ADG"),
    "tidybreed_table"
  )
  close_pop(pop)
})


test_that("define_additive_effects() errors when filter returns zero rows", {
  pop <- make_effects_pop("eff_err_empty")
  pop <- define_trait(pop, "ADG", target_add_var = 1)

  expect_error(
    pop |> get_table("genome_meta") |> dplyr::filter(locus_name == "NONEXISTENT") |>
      define_additive_effects("ADG"),
    "zero rows"
  )
  close_pop(pop)
})
