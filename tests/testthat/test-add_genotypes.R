make_chip_pop <- function(pop_name = "geno_test", n_loci = 100, chip_name = "50k",
                          n_males = 5, n_females = 5) {
  pop <- open_pop(pop_name = pop_name, db_name = ":memory:") |>
    define_genome(n_loci = n_loci, n_chr = 2, chr_len_Mb = 100) |>
    define_founder_haplotypes(n_haplotypes = 20) |>
    get_table("founder_haplotypes") |>
    add_founders(n_males = n_males, n_females = n_females, line_name = "A")

  # Filter to first chromosome as chip loci (stable, order-independent)
  pop <- pop |>
    get_table("genome_meta") |>
    dplyr::filter(chr == 1L) |>
    define_chip(chip_name)

  pop
}


test_that("add_genotypes() adds has_<chip> column to ind_meta", {
  pop <- make_chip_pop("test_ag_col")

  pop <- pop |> get_table("ind_meta") |> add_genotypes("50k")

  expect_true("has_50k" %in% DBI::dbListFields(pop$db_conn, "ind_meta"))

  ind <- dplyr::collect(get_table(pop, "ind_meta"))
  expect_true(all(ind$has_50k))  # all animals genotyped (no filter applied)

  close_pop(pop)
})

test_that("add_genotypes() with filter only marks filtered animals", {
  pop <- make_chip_pop("test_ag_filter")

  pop <- pop |>
    get_table("ind_meta") |>
    dplyr::filter(sex == "F") |>
    add_genotypes("50k")

  ind <- dplyr::collect(get_table(pop, "ind_meta"))
  expect_true(all(ind$has_50k[ind$sex == "F"]))
  expect_true(all(!ind$has_50k[ind$sex == "M"]))

  close_pop(pop)
})

test_that("add_genotypes() is additive (union semantics)", {
  pop <- make_chip_pop("test_ag_additive")

  pop <- pop |>
    get_table("ind_meta") |>
    dplyr::filter(sex == "F") |>
    add_genotypes("50k")

  pop <- pop |>
    get_table("ind_meta") |>
    dplyr::filter(sex == "M") |>
    add_genotypes("50k")

  ind <- dplyr::collect(get_table(pop, "ind_meta"))
  expect_true(all(ind$has_50k))

  close_pop(pop)
})

test_that("add_genotypes() with filter returns tidybreed_pop", {
  pop <- make_chip_pop("test_ag_class")

  pop <- pop |>
    get_table("ind_meta") |>
    dplyr::filter(sex == "F") |>
    add_genotypes("50k")

  expect_s3_class(pop, "tidybreed_pop")
  close_pop(pop)
})

test_that("add_genotypes() errors if chip not defined", {
  pop <- open_pop(pop_name = "test_ag_noship", db_name = ":memory:") |>
    define_genome(n_loci = 100, n_chr = 2, chr_len_Mb = 100) |>
    define_founder_haplotypes(n_haplotypes = 20) |>
    get_table("founder_haplotypes") |>
    add_founders(n_males = 5, n_females = 5, line_name = "A")

  expect_error(
    pop |> get_table("ind_meta") |> add_genotypes("nonexistent"),
    "not found in genome_meta"
  )
  close_pop(pop)
})

test_that("add_genotypes() supports col_name override", {
  pop <- make_chip_pop("test_ag_colname", n_males = 3, n_females = 3)

  pop <- pop |> get_table("ind_meta") |> add_genotypes("50k", col_name = "genotyped_50k")

  expect_true("genotyped_50k" %in% DBI::dbListFields(pop$db_conn, "ind_meta"))
  close_pop(pop)
})
