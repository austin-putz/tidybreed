test_that("print.tidybreed_pop() on a bare pop shows header + Database + hint only", {
  pop <- open_pop(pop_name = "bare", db_name = ":memory:")
  on.exit(close_pop(pop))

  out <- paste(capture.output(print(pop)), collapse = "\n")

  # Header rule with the population name, connection status, and hint
  expect_true(grepl("tidybreed population: bare", out))
  expect_true(grepl("in-memory", out))
  expect_true(grepl("\\[connected\\]", out))
  expect_true(grepl("schema(pop)", out, fixed = TRUE))

  # No data sections yet — every count is zero / table empty
  expect_false(grepl("Genome", out))
  expect_false(grepl("Model", out))
  expect_false(grepl("Individuals", out))
  expect_false(grepl("Records", out))
})


test_that("print.tidybreed_pop() shows a Genome section after define_genome()", {
  pop <- open_pop(pop_name = "geno", db_name = ":memory:") |>
    define_genome(n_loci = 200, n_chr = 4, chr_len_Mb = 100)
  on.exit(close_pop(pop))

  out <- paste(capture.output(print(pop)), collapse = "\n")

  expect_true(grepl("Genome", out))
  expect_true(grepl("200 loci", out))
  expect_true(grepl("4 chr", out))
  expect_true(grepl("Mb", out))
  # Still no individuals or founder pool
  expect_false(grepl("Individuals", out))
  expect_false(grepl("founder pool", out))
})


test_that("print.tidybreed_pop() shows Model / Individuals / Records for a full pop", {
  pop <- open_pop(pop_name = "full", db_name = ":memory:") |>
    define_genome(n_loci = 200, n_chr = 4, chr_len_Mb = 100) |>
    define_founder_haplotypes(n_haplotypes = 50, line_name = "Duroc") |>
    define_founder_haplotypes(n_haplotypes = 50, line_name = "Landrace")
  on.exit(close_pop(pop))

  pop <- pop |>
    get_table("founder_haplotypes") |> dplyr::filter(line_name == "Duroc") |>
    add_founders(n_males = 5, n_females = 10, line_name = "Duroc")
  pop <- pop |>
    get_table("founder_haplotypes") |> dplyr::filter(line_name == "Landrace") |>
    add_founders(n_males = 5, n_females = 10, line_name = "Landrace")

  pop <- pop |> define_trait("ADG", target_add_var = 0.25)
  pop <- pop |>
    get_table("genome_meta") |> dplyr::filter(chr %in% 1:2) |>
    define_additive_effects("ADG")
  pop <- pop |> define_phenotype("ADG", residual_var = 1)
  pop <- pop |> get_table("ind_meta") |> add_phenotype("ADG")

  out <- paste(capture.output(print(pop)), collapse = "\n")

  # Genome + founder pool (two lines of 50 → 100 haplotypes)
  expect_true(grepl("founder pool: 100 haplotypes", out))
  # Model section: trait, phenotype, QTL
  expect_true(grepl("Model", out))
  expect_true(grepl("1 trait", out))
  expect_true(grepl("1 phenotype", out))
  expect_true(grepl("QTL", out))
  # Individuals with sex and line breakdowns (two lines → by line shown)
  expect_true(grepl("Individuals", out))
  expect_true(grepl("by sex", out))
  expect_true(grepl("by line", out))
  expect_true(grepl("Duroc", out))
  expect_true(grepl("Landrace", out))
  # Records section
  expect_true(grepl("Records", out))
  expect_true(grepl("phenotypes", out))
})


test_that("print.tidybreed_pop() omits 'by line' for a single-line pop", {
  pop <- open_pop(pop_name = "oneline", db_name = ":memory:") |>
    define_genome(n_loci = 100, n_chr = 2, chr_len_Mb = 50) |>
    define_founder_haplotypes(n_haplotypes = 50, line_name = "A")
  pop <- pop |>
    get_table("founder_haplotypes") |>
    add_founders(n_males = 5, n_females = 5, line_name = "A")
  on.exit(close_pop(pop))

  out <- paste(capture.output(print(pop)), collapse = "\n")

  expect_true(grepl("Individuals", out))
  expect_true(grepl("by sex", out))
  expect_false(grepl("by line", out))
})


test_that("print.tidybreed_pop() reports a closed connection as [disconnected]", {
  pop <- open_pop(pop_name = "closed", db_name = ":memory:") |>
    define_genome(n_loci = 100, n_chr = 2, chr_len_Mb = 50)
  close_pop(pop)

  out <- paste(capture.output(res <- print(pop)), collapse = "\n")

  expect_true(grepl("\\[disconnected\\]", out))
  # No queries run against the dead connection; no data sections
  expect_false(grepl("Genome", out))
  # Returns the object invisibly
  expect_identical(res, pop)
})
