make_pheno_pop <- function(pop_name = "ph", n_ind = 200, n_loci = 400) {
  pop <- initialize_genome(
    pop_name   = pop_name,
    n_loci     = n_loci,
    n_chr      = 4,
    chr_len_Mb = 100,
    db_path    = ":memory:"
  ) |>
    define_founder_haplotypes(n_haplotypes = 100, fixed_allele_freq = 0.5)
  pop <- add_founders(pop, n_males = n_ind / 2, n_females = n_ind / 2,
                      line_name = "A")
  pop <- get_table(pop, "ind_meta") |> mutate_table(gen = 0L)
  pop
}

# Helper: sample n_qtl loci and assign additive effects for each trait
apply_random_qtl <- function(pop, trait_names, n_qtl, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  sel <- pop |> get_table("genome_meta") |> dplyr::collect() |>
    dplyr::slice_sample(n = n_qtl) |> dplyr::pull(locus_name)
  for (t in trait_names) {
    pop <- pop |>
      get_table("genome_meta") |>
      dplyr::filter(locus_name %in% sel) |>
      define_additive_effects(t)
  }
  pop
}

# Helper: define trait + phenotype + QTL for a single continuous trait
setup_simple_trait <- function(pop, trait_name, target_add_var = 0.25,
                                residual_var = 0.75, mean = 0,
                                n_qtl = 50, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  pop <- define_trait(pop, trait_name, target_add_var = target_add_var)
  pop <- apply_random_qtl(pop, trait_name, n_qtl = n_qtl)
  pop <- define_phenotype(pop, trait_name,
                          type         = "continuous",
                          mean         = mean,
                          residual_var = residual_var)
  pop
}


test_that("add_phenotype() writes records and TBVs for continuous trait", {
  set.seed(99)
  pop <- make_pheno_pop("ph_cont")
  pop <- setup_simple_trait(pop, "ADG",
                             target_add_var = 0.25, residual_var = 0.75,
                             mean = 10, n_qtl = 80, seed = 1)
  pop <- pop |> get_table("ind_meta") |> add_phenotype("ADG")

  ph <- dplyr::collect(get_table(pop, "ind_phenotype"))
  expect_equal(nrow(ph), 200)
  expect_true(all(ph$phenotype_name == "ADG"))
  expect_equal(mean(ph$pheno_value), 10, tolerance = 0.3)

  tbv <- dplyr::collect(get_table(pop, "ind_tbv"))
  expect_equal(nrow(tbv), 200)

  close_pop(pop)
})


test_that("get_table() |> filter() restricts phenotyped subset", {
  pop <- make_pheno_pop("ph_filter")
  pop <- setup_simple_trait(pop, "ADG", n_qtl = 50)

  pop <- pop |>
    get_table("ind_meta") |>
    dplyr::filter(sex == "F") |>
    add_phenotype("ADG")

  ph  <- dplyr::collect(get_table(pop, "ind_phenotype"))
  expect_equal(nrow(ph), 100)
  ind <- dplyr::collect(get_table(pop, "ind_meta"))
  expect_true(all(ind$sex[match(ph$id_ind, ind$id_ind)] == "F"))

  close_pop(pop)
})


test_that("categorical trait with prevalence respects target rate approximately", {
  set.seed(7)
  pop <- initialize_genome(
    pop_name   = "ph_categorical",
    n_loci     = 300,
    n_chr      = 3,
    chr_len_Mb = 100,
    db_path    = ":memory:"
  ) |>
    define_founder_haplotypes(n_haplotypes = 200, fixed_allele_freq = 0.5)
  pop <- add_founders(pop, n_males = 500, n_females = 500, line_name = "A")
  pop <- define_trait(pop, "mort", target_add_var = 1)
  pop <- apply_random_qtl(pop, "mort", n_qtl = 50)
  pop <- define_phenotype(pop, "mort",
                          type            = "categorical",
                          prevalence      = 0.1,
                          cat_values      = c(0, 1),
                          cat_names       = c("Alive", "Dead"),
                          store_liability = TRUE,
                          residual_var    = 1)
  pop <- pop |> get_table("ind_meta") |> add_phenotype("mort")

  ph <- dplyr::collect(get_table(pop, "ind_phenotype"))
  rate <- mean(ph$pheno_value)
  expect_true(abs(rate - 0.1) <= 0.04)
  expect_true(all(ph$pheno_value %in% c(0, 1)))
  expect_true("liability_value" %in% names(ph))
  expect_true("cat_name"        %in% names(ph))
  expect_true(all(ph$cat_name %in% c("Alive", "Dead")))

  close_pop(pop)
})


test_that("categorical trait with explicit thresholds produces correct categories", {
  set.seed(42)
  pop <- initialize_genome(
    pop_name   = "ph_cat_thresh",
    n_loci     = 200,
    n_chr      = 2,
    chr_len_Mb = 100,
    db_path    = ":memory:"
  ) |>
    define_founder_haplotypes(n_haplotypes = 100, fixed_allele_freq = 0.5)
  pop <- add_founders(pop, n_males = 250, n_females = 250, line_name = "A")
  pop <- define_trait(pop, "body_score", target_add_var = 1)
  pop <- apply_random_qtl(pop, "body_score", n_qtl = 40)
  pop <- define_phenotype(pop, "body_score",
                          type            = "categorical",
                          thresholds      = c(-1, 0, 1),
                          cat_values      = c(1, 2, 3, 4),
                          cat_names       = c("thin", "fair", "good", "excellent"),
                          store_liability = TRUE,
                          residual_var    = 1)
  pop <- pop |> get_table("ind_meta") |> add_phenotype("body_score")

  ph <- dplyr::collect(get_table(pop, "ind_phenotype"))
  expect_true(all(ph$pheno_value %in% c(1, 2, 3, 4)))
  expect_true("liability_value" %in% names(ph))
  expect_true("cat_name"        %in% names(ph))
  expect_true(all(ph$cat_name %in% c("thin", "fair", "good", "excellent")))

  close_pop(pop)
})


test_that("count trait clips to min/max", {
  set.seed(3)
  pop <- make_pheno_pop("ph_count", n_ind = 200)
  pop <- define_trait(pop, "litter", target_add_var = 4, target_add_mean = 10)
  pop <- apply_random_qtl(pop, "litter", n_qtl = 60)
  pop <- define_phenotype(pop, "litter",
                          type         = "count",
                          mean         = 10,
                          min_value    = 0,
                          max_value    = 20,
                          residual_var = 8)
  pop <- pop |> get_table("ind_meta") |> add_phenotype("litter")

  ph <- dplyr::collect(get_table(pop, "ind_phenotype"))
  expect_true(all(ph$pheno_value >= 0))
  expect_true(all(ph$pheno_value <= 20))
  expect_true(all(ph$pheno_value == round(ph$pheno_value)))

  close_pop(pop)
})


test_that("user_values override bypasses the model", {
  pop <- make_pheno_pop("ph_user", n_ind = 20)
  pop <- setup_simple_trait(pop, "ADG", n_qtl = 10)

  ids    <- dplyr::pull(dplyr::collect(get_table(pop, "ind_meta")), id_ind)
  custom <- stats::setNames(seq_along(ids) * 1.0, ids)

  pop <- pop |>
    get_table("ind_meta") |>
    dplyr::filter(.data$id_ind %in% !!ids) |>
    add_phenotype("ADG", user_values = list(ADG = custom))

  ph <- dplyr::collect(get_table(pop, "ind_phenotype"))
  expect_equal(length(ph$pheno_value), length(ids))
  expect_equal(
    sort(ph$pheno_value[match(ids, ph$id_ind)]),
    sort(unname(custom))
  )

  close_pop(pop)
})


test_that("fixed-effect class shifts phenotype by level", {
  set.seed(11)
  pop <- make_pheno_pop("ph_covariate", n_ind = 400)
  pop <- define_trait(pop, "ADG", target_add_var = 0.01)
  pop <- apply_random_qtl(pop, "ADG", n_qtl = 20)
  pop <- define_phenotype(pop, "ADG",
                          type         = "continuous",
                          residual_var = 0.01)
  pop <- pop |>
    define_effect_fixed_class("ADG", "sex",
                              source_column = "sex",
                              levels = c(M = 10, F = -10))
  pop <- pop |> get_table("ind_meta") |> add_phenotype("ADG")

  ph  <- dplyr::collect(get_table(pop, "ind_phenotype"))
  ind <- dplyr::collect(get_table(pop, "ind_meta"))
  ph$sex <- ind$sex[match(ph$id_ind, ind$id_ind)]

  mean_m <- mean(ph$pheno_value[ph$sex == "M"])
  mean_f <- mean(ph$pheno_value[ph$sex == "F"])
  expect_gt(mean_m - mean_f, 15)

  close_pop(pop)
})


# ─── ... forwarding tests ─────────────────────────────────────────────────────

make_pheno_pop_with_trait <- function() {
  set.seed(42)
  pop <- make_pheno_pop("ph_extra", n_ind = 20, n_loci = 100)
  pop <- setup_simple_trait(pop, "ADG",
                             target_add_var = 0.25, residual_var = 0.75,
                             n_qtl = 20, seed = 1)
  pop
}

test_that("add_phenotype() accepts scalar ... fields and writes to ind_phenotype", {
  pop <- make_pheno_pop_with_trait()

  pop <- pop |>
    get_table("ind_meta") |>
    add_phenotype("ADG", test_env = "barn_A")

  result <- get_table(pop, "ind_phenotype") |> dplyr::collect()
  expect_true("test_env" %in% names(result))
  expect_true(all(result$test_env == "barn_A"))

  close_pop(pop)
})

test_that("add_phenotype() ... errors for non-scalar value", {
  pop <- make_pheno_pop_with_trait()

  expect_error(
    pop |> get_table("ind_meta") |> add_phenotype("ADG", test_env = c("a", "b")),
    "scalar"
  )

  close_pop(pop)
})

test_that("add_phenotype() defaults to all phenotypes in phenotype_meta", {
  set.seed(5)
  pop <- make_pheno_pop("ph_default_pheno", n_ind = 40, n_loci = 200)
  pop <- setup_simple_trait(pop, "ADG", n_qtl = 30)
  pop <- setup_simple_trait(pop, "BW",  target_add_var = 2, residual_var = 2, n_qtl = 30)

  pop <- pop |> get_table("ind_meta") |> add_phenotype()  # no trait_name

  ph <- dplyr::collect(get_table(pop, "ind_phenotype"))
  expect_true("ADG" %in% ph$phenotype_name)
  expect_true("BW"  %in% ph$phenotype_name)

  close_pop(pop)
})


test_that("sex-limited phenotype only records correct sex", {
  set.seed(21)
  pop <- make_pheno_pop("ph_sex_lim", n_ind = 100)
  pop <- define_trait(pop, "milk", target_add_var = 10)
  pop <- apply_random_qtl(pop, "milk", n_qtl = 40)
  pop <- define_phenotype(pop, "milk",
                          type          = "continuous",
                          expressed_sex = "F",
                          mean          = 30,
                          residual_var  = 5)
  pop <- pop |> get_table("ind_meta") |> add_phenotype("milk")

  ph  <- dplyr::collect(get_table(pop, "ind_phenotype"))
  ind <- dplyr::collect(get_table(pop, "ind_meta"))
  sexes <- ind$sex[match(ph$id_ind, ind$id_ind)]
  expect_true(all(sexes == "F"))
  expect_equal(nrow(ph), 50)  # 50 females out of 100

  close_pop(pop)
})


test_that("non-repeatable phenotype refuses second record for same individual", {
  set.seed(13)
  pop <- make_pheno_pop("ph_nonrep", n_ind = 40)
  pop <- setup_simple_trait(pop, "ADG", n_qtl = 30)

  pop <- pop |> get_table("ind_meta") |> add_phenotype("ADG")
  expect_warning(
    pop <- pop |> get_table("ind_meta") |> add_phenotype("ADG"),
    "not repeatable"
  )

  ph <- dplyr::collect(get_table(pop, "ind_phenotype"))
  expect_equal(nrow(ph), 40)  # second call produced 0 new records

  close_pop(pop)
})


test_that("repeatable phenotype allows multiple records per individual", {
  set.seed(14)
  pop <- make_pheno_pop("ph_rep", n_ind = 20)
  pop <- define_trait(pop, "LS", target_add_var = 1)
  pop <- apply_random_qtl(pop, "LS", n_qtl = 20)
  pop <- define_phenotype(pop, "LS",
                          type         = "continuous",
                          repeatable   = TRUE,
                          residual_var = 1)
  pop <- define_effect_random(pop, "LS", "pe",
                               source_column = "id_ind",
                               variance = 0.5)

  pop <- pop |> get_table("ind_meta") |> add_phenotype("LS")
  pop <- pop |> get_table("ind_meta") |> add_phenotype("LS")

  ph <- dplyr::collect(get_table(pop, "ind_phenotype"))
  expect_equal(nrow(ph), 40)
  expect_true(max(ph$pheno_number) == 2)

  close_pop(pop)
})
