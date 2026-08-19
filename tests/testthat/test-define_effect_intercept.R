# Tests for define_effect_intercept()
#
# Updates phenotype_meta.mean for an existing phenotype. That value is the
# population mean added to every individual's liability by add_phenotype().

make_int_pop <- function(pop_name = "ic", n_ind = 10) {
  pop <- open_pop(pop_name = pop_name, db_name = ":memory:") |>
    define_genome(n_loci = 200, n_chr = 2, chr_len_Mb = 100) |>
    define_founder_haplotypes(n_haplotypes = 40, method = "beta")
  pop |>
    get_table("founder_haplotypes") |>
    add_founders(n_males = n_ind / 2, n_females = n_ind / 2, line_name = "A")
}

setup_int_trait <- function(pop, mean = 0, n_qtl = 30) {
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

get_pheno_mean <- function(pop, phenotype_name = "ADG") {
  DBI::dbGetQuery(pop$db_conn, paste0(
    "SELECT mean FROM phenotype_meta WHERE phenotype_name = '",
    phenotype_name, "'"))$mean
}


test_that("define_effect_intercept() updates mean in phenotype_meta", {
  set.seed(1)
  pop <- make_int_pop("ic_update")
  pop <- setup_int_trait(pop, mean = 0)
  expect_equal(get_pheno_mean(pop), 0)

  expect_message(
    pop <- define_effect_intercept(pop, "ADG", mean = 850),
    "Set mean for phenotype"
  )
  expect_equal(get_pheno_mean(pop), 850)

  close_pop(pop)
})


test_that("the updated intercept is what add_phenotype() applies", {
  set.seed(2)
  pop <- make_int_pop("ic_applied")
  pop <- setup_int_trait(pop, mean = 0)
  pop <- define_effect_intercept(pop, "ADG", mean = 500)
  pop <- pop |> get_table("ind_meta") |> add_phenotype("ADG")

  ph <- dplyr::collect(get_table(pop, "ind_phenotype"))
  tb <- dplyr::collect(get_table(pop, "ind_tbv"))
  ph$tbv <- tb$tbv_value[match(ph$id_ind, tb$id_ind)]

  # residual_var = 0, no other effects -> pheno == mean + tbv exactly
  expect_equal(ph$pheno_value, 500 + ph$tbv)

  close_pop(pop)
})


test_that("define_effect_intercept() overrides the define_phenotype() mean", {
  set.seed(3)
  pop <- make_int_pop("ic_override")
  pop <- setup_int_trait(pop, mean = 100)
  expect_equal(get_pheno_mean(pop), 100)

  pop <- define_effect_intercept(pop, "ADG", mean = -25)
  expect_equal(get_pheno_mean(pop), -25)

  pop <- pop |> get_table("ind_meta") |> add_phenotype("ADG")
  ph <- dplyr::collect(get_table(pop, "ind_phenotype"))
  tb <- dplyr::collect(get_table(pop, "ind_tbv"))
  ph$tbv <- tb$tbv_value[match(ph$id_ind, tb$id_ind)]
  expect_equal(ph$pheno_value, -25 + ph$tbv)

  close_pop(pop)
})


test_that("define_effect_intercept() only touches the named phenotype", {
  set.seed(4)
  pop <- make_int_pop("ic_scoped")
  pop <- setup_int_trait(pop, mean = 10)

  pop <- define_trait(pop, "BW", target_add_var = 0.1)
  sel <- get_table(pop, "genome_meta") |>
    dplyr::collect() |> dplyr::slice_sample(n = 20) |> dplyr::pull(locus_name)
  pop <- pop |>
    get_table("genome_meta") |>
    dplyr::filter(locus_name %in% sel) |>
    define_additive_effects("BW")
  pop <- define_phenotype(pop, "BW", type = "continuous",
                          mean = 20, residual_var = 0)

  pop <- define_effect_intercept(pop, "ADG", mean = 999)

  expect_equal(get_pheno_mean(pop, "ADG"), 999)
  expect_equal(get_pheno_mean(pop, "BW"),   20)   # untouched

  close_pop(pop)
})


test_that("define_effect_intercept() accepts negative and zero means", {
  set.seed(5)
  pop <- make_int_pop("ic_signs")
  pop <- setup_int_trait(pop, mean = 100)

  pop <- define_effect_intercept(pop, "ADG", mean = 0)
  expect_equal(get_pheno_mean(pop), 0)

  pop <- define_effect_intercept(pop, "ADG", mean = -3.5)
  expect_equal(get_pheno_mean(pop), -3.5)

  close_pop(pop)
})


test_that("define_effect_intercept() rejects invalid input", {
  set.seed(6)
  pop <- make_int_pop("ic_invalid")
  pop <- setup_int_trait(pop)

  expect_error(define_effect_intercept(pop, "NOPE", mean = 1),
               "not found in phenotype_meta")

  expect_error(define_effect_intercept(pop, "3bad", mean = 1),
               "Invalid phenotype name")

  expect_error(define_effect_intercept(pop, "ADG", mean = c(1, 2)))
  expect_error(define_effect_intercept(pop, "ADG", mean = NA_real_))
  expect_error(define_effect_intercept(pop, "ADG", mean = "high"))

  close_pop(pop)
})


test_that("define_effect_intercept() is not fooled by SQL injection in the name", {
  set.seed(7)
  pop <- make_int_pop("ic_inject")
  pop <- setup_int_trait(pop, mean = 10)

  # Rejected by the identifier validator, so phenotype_meta is untouched
  expect_error(
    define_effect_intercept(pop, "ADG'; DROP TABLE phenotype_meta; --", mean = 1)
  )
  expect_true("phenotype_meta" %in% DBI::dbListTables(pop$db_conn))
  expect_equal(get_pheno_mean(pop), 10)

  close_pop(pop)
})
