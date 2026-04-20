make_pheno_pop <- function(pop_name = "ph", n_ind = 200, n_loci = 400) {
  pop <- initialize_genome(
    pop_name        = pop_name,
    n_loci          = n_loci,
    n_chr           = 4,
    chr_len_Mb      = 100,
    n_haplotypes    = 100,
    db_path         = ":memory:",
    fixed_allele_freq = 0.5
  )
  pop <- add_founders(pop, n_males = n_ind / 2, n_females = n_ind / 2,
                      line_name = "A")
  pop <- mutate_ind_meta(pop, gen = 0L)
  pop
}


test_that("add_phenotype() writes records and TBVs for continuous trait", {
  set.seed(99)
  pop <- make_pheno_pop("ph_cont")

  pop <- pop |>
    add_trait("ADG", target_add_var = 0.25, residual_var = 0.75, mean = 10) |>
    define_qtl("ADG", n = 80, method = "random") |>
    set_qtl_effects("ADG", seed = 1) |>
    add_phenotype("ADG")

  ph <- dplyr::collect(get_table(pop, "ind_phenotype"))
  expect_equal(nrow(ph), 200)
  expect_true(all(ph$trait_name == "ADG"))
  expect_equal(mean(ph$value), 10, tolerance = 0.3)

  tbv <- dplyr::collect(get_table(pop, "ind_tbv"))
  expect_equal(nrow(tbv), 200)

  close_pop(pop)
})


test_that("filter.tidybreed_pop() restricts phenotyped subset", {
  pop <- make_pheno_pop("ph_filter")

  pop <- pop |>
    add_trait("ADG", target_add_var = 1, residual_var = 1) |>
    define_qtl("ADG", n = 50) |>
    set_qtl_effects("ADG")

  pop <- pop |>
    dplyr::filter(sex == "F") |>
    add_phenotype("ADG")

  ph <- dplyr::collect(get_table(pop, "ind_phenotype"))
  expect_equal(nrow(ph), 100)
  # All phenotyped individuals should be female
  ind <- dplyr::collect(get_table(pop, "ind_meta"))
  expect_true(all(ind$sex[match(ph$id_ind, ind$id_ind)] == "F"))

  # pending_filter should be empty after consumption
  expect_equal(length(pop$pending_filter), 0)

  close_pop(pop)
})


test_that("binary trait respects prevalence approximately", {
  set.seed(7)
  pop <- initialize_genome(
    pop_name        = "ph_binary",
    n_loci          = 300,
    n_chr           = 3,
    chr_len_Mb      = 100,
    n_haplotypes    = 200,
    db_path         = ":memory:",
    fixed_allele_freq = 0.5
  )
  pop <- add_founders(pop, n_males = 500, n_females = 500, line_name = "A")

  pop <- pop |>
    add_trait("mort", trait_type = "binary", prevalence = 0.1,
              target_add_var = 1, residual_var = 1) |>
    define_qtl("mort", n = 50) |>
    set_qtl_effects("mort") |>
    add_phenotype("mort")

  ph <- dplyr::collect(get_table(pop, "ind_phenotype"))
  rate <- mean(ph$value)
  expect_equal(rate, 0.1, tolerance = 0.03)
  expect_true(all(ph$value %in% c(0, 1)))

  close_pop(pop)
})


test_that("count trait clips to min/max", {
  set.seed(3)
  pop <- make_pheno_pop("ph_count", n_ind = 200)

  pop <- pop |>
    add_trait("litter", trait_type = "count",
              min_value = 0, max_value = 20,
              target_add_var = 4, residual_var = 8, mean = 10) |>
    define_qtl("litter", n = 60) |>
    set_qtl_effects("litter") |>
    add_phenotype("litter")

  ph <- dplyr::collect(get_table(pop, "ind_phenotype"))
  expect_true(all(ph$value >= 0))
  expect_true(all(ph$value <= 20))
  expect_true(all(ph$value == round(ph$value)))

  close_pop(pop)
})


test_that("user_values override bypasses the model", {
  pop <- make_pheno_pop("ph_user", n_ind = 20)

  pop <- pop |>
    add_trait("ADG", target_add_var = 1, residual_var = 1) |>
    define_qtl("ADG", n = 10) |>
    set_qtl_effects("ADG")

  ids <- dplyr::pull(dplyr::collect(get_table(pop, "ind_meta")), id_ind)
  custom <- stats::setNames(seq_along(ids) * 1.0, ids)

  pop <- pop |>
    dplyr::filter(id_ind %in% !!ids) |>
    add_phenotype("ADG", user_values = list(ADG = custom))

  ph <- dplyr::collect(get_table(pop, "ind_phenotype"))
  expect_equal(length(ph$value), length(ids))
  # Exact match by id_ind
  expect_equal(
    sort(ph$value[match(ids, ph$id_ind)]),
    sort(unname(custom))
  )

  close_pop(pop)
})


test_that("fixed-effect covariate shifts phenotype by level", {
  set.seed(11)
  pop <- make_pheno_pop("ph_covariate", n_ind = 400)

  pop <- pop |>
    add_trait("ADG", target_add_var = 0.01, residual_var = 0.01, mean = 0) |>
    define_qtl("ADG", n = 20) |>
    set_qtl_effects("ADG") |>
    add_trait_covariate("ADG", "sex",
                        effect_class = "fixed", source_column = "sex",
                        levels = c(M = 10, F = -10)) |>
    add_phenotype("ADG")

  ph <- dplyr::collect(get_table(pop, "ind_phenotype"))
  ind <- dplyr::collect(get_table(pop, "ind_meta"))
  ph$sex <- ind$sex[match(ph$id_ind, ind$id_ind)]

  mean_m <- mean(ph$value[ph$sex == "M"])
  mean_f <- mean(ph$value[ph$sex == "F"])
  expect_gt(mean_m - mean_f, 15)      # difference is ~20 with small noise

  close_pop(pop)
})
