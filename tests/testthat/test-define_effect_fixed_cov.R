# Tests for define_effect_fixed_cov()
#
# Contribution model applied by add_phenotype():
#   slope * (source_column - center)^poly_order

make_cov_pop <- function(pop_name = "fc", n_ind = 10, ages = NULL) {
  pop <- open_pop(pop_name = pop_name, db_name = ":memory:") |>
    define_genome(n_loci = 200, n_chr = 2, chr_len_Mb = 100) |>
    define_founder_haplotypes(n_haplotypes = 40, method = "beta")
  pop <- pop |>
    get_table("founder_haplotypes") |>
    add_founders(n_males = n_ind / 2, n_females = n_ind / 2, line_name = "A")

  ids <- get_table(pop, "ind_meta") |> dplyr::pull(id_ind)
  if (is.null(ages)) ages <- seq(10, by = 10, length.out = length(ids))
  pop <- get_table(pop, "ind_meta") |>
    mutate_table(age_days = tibble::tibble(id_ind = ids, age_days = ages))
  pop
}

# Trait + QTL + phenotype with a deterministic (zero-residual) phenotype model
setup_cov_trait <- function(pop, mean = 100, n_qtl = 30) {
  pop <- define_trait(pop, "ADG", target_add_var = 0.25)
  sel <- get_table(pop, "genome_meta") |>
    dplyr::collect() |>
    dplyr::slice_sample(n = n_qtl) |>
    dplyr::pull(locus_name)
  pop <- pop |>
    get_table("genome_meta") |>
    dplyr::filter(locus_name %in% sel) |>
    define_additive_effects("ADG")
  define_phenotype(pop, "ADG", type = "continuous",
                   mean = mean, residual_var = 0)
}

# Join phenotype records to age and TBV for exact-arithmetic assertions
cov_frame <- function(pop) {
  ph <- dplyr::collect(get_table(pop, "ind_phenotype"))
  im <- dplyr::collect(get_table(pop, "ind_meta"))
  tb <- dplyr::collect(get_table(pop, "ind_tbv"))
  ph$age <- im$age_days[match(ph$id_ind, im$id_ind)]
  ph$tbv <- tb$tbv_value[match(ph$id_ind, tb$id_ind)]
  ph
}


test_that("define_effect_fixed_cov() writes a correct phenotype_effects row", {
  set.seed(1)
  pop <- make_cov_pop("fc_row")
  pop <- setup_cov_trait(pop)
  pop <- define_effect_fixed_cov(pop, "ADG", "age",
                                 source_column = "age_days",
                                 slope = 2, center = 0)

  row <- DBI::dbGetQuery(pop$db_conn,
    "SELECT * FROM phenotype_effects WHERE effect_name = 'age'")

  expect_equal(nrow(row), 1)
  expect_equal(row$phenotype_name, "ADG")
  expect_equal(row$effect_class,   "fixed_cov")
  expect_equal(row$source_column,  "age_days")
  expect_equal(row$source_table,   "ind_meta")
  expect_equal(row$slope,          2)
  expect_equal(row$center,         0)
  expect_equal(row$poly_order,     1L)
  expect_true(is.na(row$levels_json))
  expect_true(is.na(row$distribution))

  close_pop(pop)
})


test_that("linear covariate contributes slope * (x - center) exactly", {
  set.seed(2)
  pop <- make_cov_pop("fc_linear")
  pop <- setup_cov_trait(pop, mean = 100)
  pop <- define_effect_fixed_cov(pop, "ADG", "age",
                                 source_column = "age_days",
                                 slope = 2, center = 0)
  pop <- pop |> get_table("ind_meta") |> add_phenotype("ADG")

  ph <- cov_frame(pop)
  expect_equal(ph$pheno_value, 100 + ph$tbv + 2 * ph$age)

  close_pop(pop)
})


test_that("center is subtracted before scaling by slope", {
  set.seed(3)
  pop <- make_cov_pop("fc_center")
  pop <- setup_cov_trait(pop, mean = 0)
  pop <- define_effect_fixed_cov(pop, "ADG", "age",
                                 source_column = "age_days",
                                 slope = 3, center = 55)
  pop <- pop |> get_table("ind_meta") |> add_phenotype("ADG")

  ph <- cov_frame(pop)
  expect_equal(ph$pheno_value, ph$tbv + 3 * (ph$age - 55))

  close_pop(pop)
})


test_that("poly_order = 2 gives a quadratic contribution", {
  set.seed(4)
  pop <- make_cov_pop("fc_quad")
  pop <- setup_cov_trait(pop, mean = 0)
  pop <- define_effect_fixed_cov(pop, "ADG", "age",
                                 source_column = "age_days",
                                 slope = 0.5, center = 50, poly_order = 2L)
  pop <- pop |> get_table("ind_meta") |> add_phenotype("ADG")

  ph <- cov_frame(pop)
  expect_equal(ph$pheno_value, ph$tbv + 0.5 * (ph$age - 50)^2)

  close_pop(pop)
})


test_that("NULL center is auto-computed as the column mean", {
  set.seed(5)
  # ages 10..100 -> mean 55
  pop <- make_cov_pop("fc_autocenter")
  pop <- setup_cov_trait(pop, mean = 0)

  expect_message(
    pop <- define_effect_fixed_cov(pop, "ADG", "age",
                                   source_column = "age_days",
                                   slope = 1, center = NULL),
    "Auto-computed center"
  )

  row <- DBI::dbGetQuery(pop$db_conn,
    "SELECT center FROM phenotype_effects WHERE effect_name = 'age'")
  expect_equal(row$center, 55)

  # And the stored centre is the one actually applied
  pop <- pop |> get_table("ind_meta") |> add_phenotype("ADG")
  ph  <- cov_frame(pop)
  expect_equal(ph$pheno_value, ph$tbv + (ph$age - 55))

  close_pop(pop)
})


test_that("auto-center errors when the source column is all NULL", {
  set.seed(6)
  pop <- make_cov_pop("fc_allnull")
  pop <- setup_cov_trait(pop)
  pop <- get_table(pop, "ind_meta") |> mutate_table(empty_col = NA_real_)

  expect_error(
    define_effect_fixed_cov(pop, "ADG", "age",
                            source_column = "empty_col", slope = 1),
    "Could not auto-compute center"
  )

  close_pop(pop)
})


test_that("multiple covariate effects are additive", {
  set.seed(7)
  pop <- make_cov_pop("fc_multi")
  pop <- setup_cov_trait(pop, mean = 0)

  ids <- get_table(pop, "ind_meta") |> dplyr::pull(id_ind)
  pop <- get_table(pop, "ind_meta") |>
    mutate_table(weight = tibble::tibble(
      id_ind = ids,
      weight = seq(5, by = 5, length.out = length(ids))
    ))

  pop <- define_effect_fixed_cov(pop, "ADG", "age",
                                 source_column = "age_days",
                                 slope = 2, center = 0)
  pop <- define_effect_fixed_cov(pop, "ADG", "wt",
                                 source_column = "weight",
                                 slope = -1, center = 0)
  pop <- pop |> get_table("ind_meta") |> add_phenotype("ADG")

  ph <- cov_frame(pop)
  im <- dplyr::collect(get_table(pop, "ind_meta"))
  ph$wt <- im$weight[match(ph$id_ind, im$id_ind)]

  expect_equal(ph$pheno_value, ph$tbv + 2 * ph$age - ph$wt)

  close_pop(pop)
})


test_that("duplicate effect_name requires overwrite = TRUE", {
  set.seed(8)
  pop <- make_cov_pop("fc_dup")
  pop <- setup_cov_trait(pop)
  pop <- define_effect_fixed_cov(pop, "ADG", "age",
                                 source_column = "age_days",
                                 slope = 2, center = 0)

  expect_error(
    define_effect_fixed_cov(pop, "ADG", "age",
                            source_column = "age_days",
                            slope = 9, center = 0),
    "already exists"
  )

  pop <- define_effect_fixed_cov(pop, "ADG", "age",
                                 source_column = "age_days",
                                 slope = 9, center = 0, overwrite = TRUE)

  row <- DBI::dbGetQuery(pop$db_conn,
    "SELECT * FROM phenotype_effects WHERE effect_name = 'age'")
  expect_equal(nrow(row), 1)   # replaced, not appended
  expect_equal(row$slope, 9)

  close_pop(pop)
})


test_that("define_effect_fixed_cov() rejects invalid input", {
  set.seed(9)
  pop <- make_cov_pop("fc_invalid")
  pop <- setup_cov_trait(pop)

  # Unknown phenotype
  expect_error(
    define_effect_fixed_cov(pop, "NOPE", "age",
                            source_column = "age_days", slope = 1),
    "not found in phenotype_meta"
  )

  # Unknown source table
  expect_error(
    define_effect_fixed_cov(pop, "ADG", "age",
                            source_column = "age_days", slope = 1,
                            source_table = "no_such_table"),
    "does not exist"
  )

  # poly_order must be >= 1
  expect_error(
    define_effect_fixed_cov(pop, "ADG", "age",
                            source_column = "age_days", slope = 1,
                            poly_order = 0L),
    "poly_order"
  )

  # SQL-unsafe identifiers
  expect_error(
    define_effect_fixed_cov(pop, "ADG", "3bad",
                            source_column = "age_days", slope = 1),
    "Invalid effect name"
  )

  # slope must be a numeric scalar
  expect_error(
    define_effect_fixed_cov(pop, "ADG", "age",
                            source_column = "age_days", slope = c(1, 2))
  )

  close_pop(pop)
})
