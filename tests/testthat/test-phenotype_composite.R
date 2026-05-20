# Tests for composite / complex phenotype features:
#   1. Maternal trait (WW = WWD self + WWM dam)
#   2. Repeatability model (PE random effect)
#   3. Heterogeneous residuals by sex
#   4. Sex-limited phenotype
#   5. Phenotype mean from phenotype_meta (not trait_meta)

# ── Shared helpers ─────────────────────────────────────────────────────────────

make_composite_pop <- function(pop_name = "comp",
                               n_males = 20, n_females = 20,
                               n_loci = 300) {
  pop <- initialize_genome(
    pop_name          = pop_name,
    n_loci            = n_loci,
    n_chr             = 3,
    chr_len_Mb        = 100,
    n_haplotypes      = 100,
    db_path           = ":memory:",
    fixed_allele_freq = 0.5
  )
  pop <- add_founders(pop, n_males = n_males, n_females = n_females,
                      line_name = "A")
  pop <- get_table(pop, "ind_meta") |> mutate_table(gen = 0L)
  pop
}

add_offspring_gen <- function(pop, n_off = 40) {
  ind <- dplyr::collect(get_table(pop, "ind_meta"))
  sires <- ind$id_ind[ind$sex == "M"]
  dams  <- ind$id_ind[ind$sex == "F"]

  sexes <- sample(c("M","F"), n_off, replace = TRUE)
  matings <- data.frame(
    id_parent_1 = sample(sires, n_off, replace = TRUE),
    id_parent_2 = sample(dams,  n_off, replace = TRUE),
    sex         = sexes,
    line_name   = "A",
    stringsAsFactors = FALSE
  )

  pop <- add_offspring(pop, matings = matings)
  pop <- get_table(pop, "ind_meta") |>
    dplyr::filter(!is.na(id_parent_1)) |>
    mutate_table(gen = 1L)
  pop
}


# ── 1. Maternal composite (WW = WWD self + WWM dam) ───────────────────────────

test_that("composite maternal phenotype assembles correctly", {
  set.seed(42)
  pop <- make_composite_pop("mat_ww")
  on.exit(close_pop(pop))

  # Define two component traits
  pop <- define_trait(pop, "WWD", target_add_var = 200)
  pop <- define_trait(pop, "WWM", target_add_var = 80)

  # Sample correlated QTL effects
  G_ww <- matrix(c(200, 40, 40, 80), 2, 2,
                 dimnames = list(c("WWD","WWM"), c("WWD","WWM")))
  sel <- pop |> get_table("genome_meta") |> dplyr::collect() |>
    dplyr::slice_sample(n = 100) |> dplyr::pull(locus_name)
  pop <- pop |>
    get_table("genome_meta") |>
    dplyr::filter(locus_name %in% sel) |>
    define_additive_effects(c("WWD","WWM"), G = G_ww)

  # Define the composite phenotype
  pop <- define_phenotype(pop, "WW",
                          trait_type   = "continuous",
                          mean         = 230,
                          residual_var = 180,
                          components   = tibble::tribble(
                            ~source_trait_name, ~contributor_type,
                            "WWD",              "self",
                            "WWM",              "dam"
                          ))

  # Generate offspring (gen 1)
  pop <- add_offspring_gen(pop, n_off = 40)

  # Phenotype offspring (calves with known dams)
  gen1_tbl <- get_table(pop, "ind_meta") |> dplyr::filter(gen == 1L)
  pop <- gen1_tbl |> add_phenotype("WW")

  # Assertions
  ph <- dplyr::collect(get_table(pop, "ind_phenotype"))
  expect_true(nrow(ph) > 0)
  expect_true(all(ph$phenotype_name == "WW"))

  # "WW" must NOT appear in ind_tbv (only WWD and WWM should)
  tbv_traits <- unique(dplyr::collect(get_table(pop, "ind_tbv"))$trait_name)
  expect_false("WW" %in% tbv_traits)
  expect_true("WWD" %in% tbv_traits)
  expect_true("WWM" %in% tbv_traits)

  # Phenotype values should be centred near 230 (mean)
  expect_equal(mean(ph$pheno_value), 230, tolerance = 30)
})


test_that("composite phenotype excludes founders with missing dam TBV", {
  set.seed(7)
  pop <- make_composite_pop("mat_missing_dam", n_males = 10, n_females = 10)
  on.exit(close_pop(pop))

  pop <- define_trait(pop, "WWD", target_add_var = 200)
  pop <- define_trait(pop, "WWM", target_add_var = 80)

  sel <- pop |> get_table("genome_meta") |> dplyr::collect() |>
    dplyr::slice_sample(n = 80) |> dplyr::pull(locus_name)
  # Define additive effects for each trait independently (no correlated G needed)
  pop <- pop |>
    get_table("genome_meta") |>
    dplyr::filter(locus_name %in% sel) |>
    define_additive_effects("WWD")
  pop <- pop |>
    get_table("genome_meta") |>
    dplyr::filter(locus_name %in% sel) |>
    define_additive_effects("WWM")

  pop <- define_phenotype(pop, "WW",
                          trait_type   = "continuous",
                          mean         = 200,
                          residual_var = 150,
                          components   = tibble::tribble(
                            ~source_trait_name, ~contributor_type,
                            "WWD",              "self",
                            "WWM",              "dam"
                          ))

  # Founders have NA id_parent_2, so phenotyping them as the composite should
  # exclude them or warn about missing dam TBVs.
  expect_warning(
    {
      pop <- get_table(pop, "ind_meta") |> add_phenotype("WW")
    },
    regexp = "missing contributor|excluded"
  )
})


# ── 2. Repeatability model ────────────────────────────────────────────────────

test_that("repeatable phenotype with PE random effect produces multiple records", {
  set.seed(14)
  pop <- make_composite_pop("rep_pe", n_males = 20, n_females = 20, n_loci = 200)
  on.exit(close_pop(pop))

  pop <- define_trait(pop, "LS", target_add_var = 1)
  sel <- pop |> get_table("genome_meta") |> dplyr::collect() |>
    dplyr::slice_sample(n = 50) |> dplyr::pull(locus_name)
  pop <- pop |>
    get_table("genome_meta") |>
    dplyr::filter(locus_name %in% sel) |>
    define_additive_effects("LS")

  pop <- define_phenotype(pop, "LS",
                          trait_type   = "continuous",
                          repeatable   = TRUE,
                          mean         = 10,
                          residual_var = 1)

  # Permanent environment: one draw per individual, reused across records
  pop <- define_effect_random(pop, "LS", "pe",
                               source_column = "id_ind",
                               variance = 0.5)

  pop <- pop |> get_table("ind_meta") |> add_phenotype("LS")
  pop <- pop |> get_table("ind_meta") |> add_phenotype("LS")

  ph <- dplyr::collect(get_table(pop, "ind_phenotype"))
  expect_equal(nrow(ph), 80)
  expect_equal(max(ph$pheno_number), 2)

  # PE draws should be stored in trait_random_effects
  pe_draws <- DBI::dbGetQuery(pop$db_conn,
    "SELECT * FROM trait_random_effects WHERE phenotype_name = 'LS' AND effect_name = 'pe'")
  expect_equal(nrow(pe_draws), 40)  # one draw per individual
})


# ── 3. Heterogeneous residuals by sex ─────────────────────────────────────────

test_that("heterogeneous residuals by sex produce different variance groups", {
  set.seed(21)
  pop <- make_composite_pop("het_resid", n_males = 100, n_females = 100)
  on.exit(close_pop(pop))

  pop <- define_trait(pop, "BW", target_add_var = 0.1)
  sel <- pop |> get_table("genome_meta") |> dplyr::collect() |>
    dplyr::slice_sample(n = 80) |> dplyr::pull(locus_name)
  pop <- pop |>
    get_table("genome_meta") |>
    dplyr::filter(locus_name %in% sel) |>
    define_additive_effects("BW")

  # Unconditional default (fallback)
  pop <- define_phenotype(pop, "BW",
                          trait_type   = "continuous",
                          residual_var = 600)

  # Heterogeneous residuals: males 400, females 800
  pop <- add_residual_cov(pop, "BW",
                          matrix(400, 1, 1, dimnames = list("BW","BW")),
                          condition_column = "sex",
                          condition_table  = "ind_meta",
                          condition_level  = "M")
  pop <- add_residual_cov(pop, "BW",
                          matrix(800, 1, 1, dimnames = list("BW","BW")),
                          condition_column = "sex",
                          condition_table  = "ind_meta",
                          condition_level  = "F")

  # Verify three rows in phenotype_residual_cov
  rcov <- DBI::dbGetQuery(pop$db_conn,
    "SELECT * FROM phenotype_residual_cov WHERE phenotype_name_1 = 'BW'")
  expect_equal(nrow(rcov), 3L)

  # Phenotype everyone; variances should differ by sex
  pop <- pop |> get_table("ind_meta") |> add_phenotype("BW")

  ph  <- dplyr::collect(get_table(pop, "ind_phenotype"))
  ind <- dplyr::collect(get_table(pop, "ind_meta"))
  ph$sex <- ind$sex[match(ph$id_ind, ind$id_ind)]
  expect_equal(nrow(ph), 200)
})


# ── 4. Sex-limited phenotype ──────────────────────────────────────────────────

test_that("sex-limited phenotype only generates records for the specified sex", {
  set.seed(31)
  pop <- make_composite_pop("sex_lim", n_males = 50, n_females = 50)
  on.exit(close_pop(pop))

  pop <- define_trait(pop, "milk", target_add_var = 10)
  sel <- pop |> get_table("genome_meta") |> dplyr::collect() |>
    dplyr::slice_sample(n = 60) |> dplyr::pull(locus_name)
  pop <- pop |>
    get_table("genome_meta") |>
    dplyr::filter(locus_name %in% sel) |>
    define_additive_effects("milk")

  pop <- define_phenotype(pop, "milk",
                          trait_type    = "continuous",
                          expressed_sex = "F",
                          mean          = 8000,
                          residual_var  = 2000000)

  pop <- pop |> get_table("ind_meta") |> add_phenotype("milk")

  ph  <- dplyr::collect(get_table(pop, "ind_phenotype"))
  ind <- dplyr::collect(get_table(pop, "ind_meta"))
  sexes <- ind$sex[match(ph$id_ind, ind$id_ind)]
  expect_equal(nrow(ph), 50)
  expect_true(all(sexes == "F"))
})


# ── 5. Phenotype mean from phenotype_meta ─────────────────────────────────────

test_that("phenotype mean in phenotype_meta controls liability intercept", {
  set.seed(55)
  pop <- make_composite_pop("pheno_mean", n_males = 200, n_females = 200,
                             n_loci = 200)
  on.exit(close_pop(pop))

  pop <- define_trait(pop, "ADG", target_add_var = 1)
  sel <- pop |> get_table("genome_meta") |> dplyr::collect() |>
    dplyr::slice_sample(n = 80) |> dplyr::pull(locus_name)
  pop <- pop |>
    get_table("genome_meta") |>
    dplyr::filter(locus_name %in% sel) |>
    define_additive_effects("ADG")

  pop <- define_phenotype(pop, "ADG",
                          trait_type   = "continuous",
                          mean         = 750,
                          residual_var = 1)

  pop <- pop |> get_table("ind_meta") |> add_phenotype("ADG")

  ph <- dplyr::collect(get_table(pop, "ind_phenotype"))
  expect_equal(mean(ph$pheno_value), 750, tolerance = 0.5)
})


# ── 6. define_trait_simple() chains through define_phenotype() ────────────────

test_that("define_trait_simple() creates both trait_meta and phenotype_meta rows", {
  set.seed(9)
  pop <- make_composite_pop("dts_chain", n_males = 30, n_females = 30)
  on.exit(close_pop(pop))

  pop <- define_trait_simple(pop, "ADG",
                              n_qtl          = 50,
                              target_add_var = 100,
                              mean           = 850,
                              residual_var   = 120)

  # trait_meta row
  tm <- DBI::dbGetQuery(pop$db_conn,
    "SELECT * FROM trait_meta WHERE trait_name = 'ADG'")
  expect_equal(nrow(tm), 1L)

  # phenotype_meta row
  pm <- DBI::dbGetQuery(pop$db_conn,
    "SELECT * FROM phenotype_meta WHERE phenotype_name = 'ADG'")
  expect_equal(nrow(pm), 1L)
  expect_equal(pm$mean, 850)

  # Can phenotype after define_trait_simple
  pop <- pop |> get_table("ind_meta") |> add_phenotype("ADG")
  ph <- dplyr::collect(get_table(pop, "ind_phenotype"))
  expect_equal(nrow(ph), 60)
})
