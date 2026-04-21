test_that("extract_genotypes() returns a tibble with correct dimensions", {
  n_chip <- 60
  pop <- initialize_genome(
    pop_name = "test", n_loci = 100, n_chr = 2,
    chr_len_Mb = 100, n_haplotypes = 20, db_path = ":memory:"
  ) |>
    add_founders(n_males = 5, n_females = 5, line_name = "A") |>
    define_chip("50k", n = n_chip, method = "random") |>
    add_genotypes("50k")

  geno <- extract_genotypes(pop, "50k")

  expect_s3_class(geno, "tbl_df")
  expect_equal(nrow(geno), 10)
  expect_equal(ncol(geno), n_chip + 1L)  # id_ind + chip loci
  expect_true("id_ind" %in% names(geno))

  close_pop(pop)
})

test_that("extract_genotypes() respects filter()", {
  pop <- initialize_genome(
    pop_name = "test", n_loci = 100, n_chr = 2,
    chr_len_Mb = 100, n_haplotypes = 20, db_path = ":memory:"
  ) |>
    add_founders(n_males = 5, n_females = 5, line_name = "A") |>
    define_chip("50k", n = 60, method = "random") |>
    add_genotypes("50k")

  geno <- pop |>
    dplyr::filter(sex == "F") |>
    extract_genotypes("50k")

  expect_equal(nrow(geno), 5)
  close_pop(pop)
})

test_that("extract_genotypes() intersects filter with has_<chip> == TRUE", {
  pop <- initialize_genome(
    pop_name = "test", n_loci = 100, n_chr = 2,
    chr_len_Mb = 100, n_haplotypes = 20, db_path = ":memory:"
  ) |>
    add_founders(n_males = 5, n_females = 5, line_name = "A") |>
    define_chip("50k", n = 60, method = "random")

  # Only genotype females
  pop <- pop |>
    dplyr::filter(sex == "F") |>
    add_genotypes("50k")

  # Extract all — should only return females (males not genotyped)
  geno <- extract_genotypes(pop, "50k")
  expect_equal(nrow(geno), 5)

  close_pop(pop)
})

test_that("extract_genotypes() returns only chip loci columns", {
  n_chip <- 40
  pop <- initialize_genome(
    pop_name = "test", n_loci = 100, n_chr = 2,
    chr_len_Mb = 100, n_haplotypes = 20, db_path = ":memory:"
  ) |>
    add_founders(n_males = 3, n_females = 3, line_name = "A") |>
    define_chip("HD", n = n_chip, method = "even") |>
    add_genotypes("HD")

  geno <- extract_genotypes(pop, "HD")

  locus_cols <- grep("^locus_", names(geno), value = TRUE)
  expect_equal(length(locus_cols), n_chip)

  # All genotype values must be 0, 1, or 2
  geno_vals <- unlist(geno[, locus_cols])
  expect_true(all(geno_vals %in% c(0L, 1L, 2L)))

  close_pop(pop)
})

test_that("extract_genotypes() errors if add_genotypes() not called", {
  pop <- initialize_genome(
    pop_name = "test", n_loci = 100, n_chr = 2,
    chr_len_Mb = 100, n_haplotypes = 20, db_path = ":memory:"
  ) |>
    add_founders(n_males = 5, n_females = 5, line_name = "A") |>
    define_chip("50k", n = 60, method = "random")

  expect_error(extract_genotypes(pop, "50k"), "not found in ind_meta")
  close_pop(pop)
})

test_that("extract_genotypes() errors if chip not defined", {
  pop <- initialize_genome(
    pop_name = "test", n_loci = 100, n_chr = 2,
    chr_len_Mb = 100, n_haplotypes = 20, db_path = ":memory:"
  ) |>
    add_founders(n_males = 5, n_females = 5, line_name = "A")

  expect_error(extract_genotypes(pop, "nonexistent"), "not found in genome_meta")
  close_pop(pop)
})
