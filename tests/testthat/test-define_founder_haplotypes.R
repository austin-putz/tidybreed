test_that("define_founder_haplotypes() creates haplotypes with uniform frequencies", {

  temp_db <- tempfile(fileext = ".duckdb")

  pop <- initialize_genome(
    pop_name = "test_founders",
    n_loci = 100,
    n_chr = 2,
    chr_len_Mb = 50,
    db_path = temp_db
  ) |>
    define_founder_haplotypes(
      n_haplotypes    = 50,
      min_allele_freq = 0.1,
      max_allele_freq = 0.9
    )

  expect_true("founder_haplotypes" %in% pop$tables)

  founder_haps <- get_table(pop, "founder_haplotypes") |> dplyr::collect()
  expect_equal(nrow(founder_haps), 50)
  expect_equal(ncol(founder_haps), 101)  # hap_id + 100 loci

  expect_equal(founder_haps$hap_id, paste0("hap_", 1:50))

  locus_cols <- paste0("locus_", 1:100)
  alleles <- as.matrix(founder_haps[, locus_cols])
  expect_true(all(alleles %in% c(0, 1)))

  genome_meta <- get_table(pop, "genome_meta") |> dplyr::collect()
  expect_true("founder_allele_freq" %in% colnames(genome_meta))
  expect_equal(length(genome_meta$founder_allele_freq), 100)
  expect_true(all(genome_meta$founder_allele_freq >= 0.1))
  expect_true(all(genome_meta$founder_allele_freq <= 0.9))

  close_pop(pop)
  unlink(temp_db)
})


test_that("define_founder_haplotypes() creates haplotypes with fixed frequency", {

  temp_db <- tempfile(fileext = ".duckdb")

  pop <- initialize_genome(
    pop_name = "test_fixed",
    n_loci = 50,
    n_chr = 2,
    chr_len_Mb = 50,
    db_path = temp_db
  ) |>
    define_founder_haplotypes(
      n_haplotypes      = 30,
      fixed_allele_freq = 0.5
    )

  expect_true("founder_haplotypes" %in% pop$tables)

  genome_meta <- get_table(pop, "genome_meta") |> dplyr::collect()
  expect_true(all(genome_meta$founder_allele_freq == 0.5))

  expect_equal(nrow(dplyr::collect(get_table(pop, "founder_haplotypes"))), 30L)

  close_pop(pop)
  unlink(temp_db)
})


test_that("initialize_genome() alone does not create founder_haplotypes", {

  pop <- initialize_genome(
    pop_name = "test_no_founders",
    n_loci = 100,
    n_chr = 2,
    chr_len_Mb = 50,
    db_path = ":memory:"
  )

  expect_false("founder_haplotypes" %in% pop$tables)

  genome_meta <- get_table(pop, "genome_meta") |> dplyr::collect()
  expect_false("founder_allele_freq" %in% colnames(genome_meta))

  close_pop(pop)
})


test_that("define_founder_haplotypes() errors if called twice", {

  pop <- initialize_genome(
    pop_name = "t",
    n_loci = 10,
    n_chr = 1,
    chr_len_Mb = 10,
    db_path = ":memory:"
  ) |>
    define_founder_haplotypes(n_haplotypes = 10)

  expect_error(
    define_founder_haplotypes(pop, n_haplotypes = 10),
    "already exists"
  )

  close_pop(pop)
})


test_that("define_founder_haplotypes() is pipe-friendly and returns tidybreed_pop", {

  pop <- initialize_genome(
    pop_name = "t",
    n_loci = 10,
    n_chr = 1,
    chr_len_Mb = 10,
    db_path = ":memory:"
  ) |>
    define_founder_haplotypes(n_haplotypes = 5)

  expect_true(inherits(pop, "tidybreed_pop"))
  expect_true("founder_haplotypes" %in% pop$tables)

  close_pop(pop)
})
