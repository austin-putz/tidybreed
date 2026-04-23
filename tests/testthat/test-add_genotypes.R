test_that("add_genotypes() adds has_<chip> column to ind_meta", {
  pop <- initialize_genome(
    pop_name = "test", n_loci = 100, n_chr = 2,
    chr_len_Mb = 100, n_haplotypes = 20, db_path = ":memory:"
  ) |>
    add_founders(n_males = 5, n_females = 5, line_name = "A") |>
    define_chip("50k", n = 80, method = "random")

  pop <- pop |> get_table("ind_meta") |> add_genotypes("50k")

  expect_true("has_50k" %in% DBI::dbListFields(pop$db_conn, "ind_meta"))

  ind <- dplyr::collect(get_table(pop, "ind_meta"))
  expect_true(all(ind$has_50k))  # all animals genotyped (no filter applied)

  close_pop(pop)
})

test_that("add_genotypes() with filter only marks filtered animals", {
  pop <- initialize_genome(
    pop_name = "test", n_loci = 100, n_chr = 2,
    chr_len_Mb = 100, n_haplotypes = 20, db_path = ":memory:"
  ) |>
    add_founders(n_males = 5, n_females = 5, line_name = "A") |>
    define_chip("50k", n = 80, method = "random")

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
  pop <- initialize_genome(
    pop_name = "test", n_loci = 100, n_chr = 2,
    chr_len_Mb = 100, n_haplotypes = 20, db_path = ":memory:"
  ) |>
    add_founders(n_males = 5, n_females = 5, line_name = "A") |>
    define_chip("50k", n = 80, method = "random")

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
  pop <- initialize_genome(
    pop_name = "test", n_loci = 100, n_chr = 2,
    chr_len_Mb = 100, n_haplotypes = 20, db_path = ":memory:"
  ) |>
    add_founders(n_males = 5, n_females = 5, line_name = "A") |>
    define_chip("50k", n = 80, method = "random")

  pop <- pop |>
    get_table("ind_meta") |>
    dplyr::filter(sex == "F") |>
    add_genotypes("50k")

  expect_s3_class(pop, "tidybreed_pop")
  close_pop(pop)
})

test_that("add_genotypes() errors if chip not defined", {
  pop <- initialize_genome(
    pop_name = "test", n_loci = 100, n_chr = 2,
    chr_len_Mb = 100, n_haplotypes = 20, db_path = ":memory:"
  ) |>
    add_founders(n_males = 5, n_females = 5, line_name = "A")

  expect_error(
    pop |> get_table("ind_meta") |> add_genotypes("nonexistent"),
    "not found in genome_meta"
  )
  close_pop(pop)
})

test_that("add_genotypes() supports col_name override", {
  pop <- initialize_genome(
    pop_name = "test", n_loci = 100, n_chr = 2,
    chr_len_Mb = 100, n_haplotypes = 20, db_path = ":memory:"
  ) |>
    add_founders(n_males = 3, n_females = 3, line_name = "A") |>
    define_chip("50k", n = 60, method = "random")

  pop <- pop |> get_table("ind_meta") |> add_genotypes("50k", col_name = "genotyped_50k")

  expect_true("genotyped_50k" %in% DBI::dbListFields(pop$db_conn, "ind_meta"))
  close_pop(pop)
})
