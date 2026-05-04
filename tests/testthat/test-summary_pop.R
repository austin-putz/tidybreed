make_summary_pop <- function() {
  pop <- initialize_genome(
    pop_name = "SumTest",
    n_loci = 100,
    n_chr = 2,
    chr_len_Mb = 50,
    n_haplotypes = 30,
    db_path = ":memory:",
    fixed_allele_freq = 0.5
  )
  add_founders(pop, n_males = 10, n_females = 10, line_name = "A", gen = 0L)
}


test_that("summary.tidybreed_pop returns tidybreed_summary object", {
  pop  <- make_summary_pop()
  summ <- summary(pop)
  expect_s3_class(summ, "tidybreed_summary")
  expect_true(is.list(summ$tables))
  expect_true("ind_meta" %in% names(summ$tables))
  close_pop(pop)
})


test_that("summary reports correct row and column counts for ind_meta", {
  pop  <- make_summary_pop()
  summ <- summary(pop)
  meta <- summ$tables[["ind_meta"]]
  expect_equal(meta$n_rows, 20L)
  expect_gte(meta$n_cols, 6L)  # 5 core + gen
  close_pop(pop)
})


test_that("genome_haplotype is flagged as wide and locus cols not summarized", {
  pop  <- make_summary_pop()
  summ <- summary(pop)
  hap  <- summ$tables[["genome_haplotype"]]
  expect_true(hap$is_wide_genome)
  expect_gte(hap$n_locus_cols, 100L)
  # Only id_ind and parent_origin should be summarized
  expect_true(all(hap$col_summaries$col_name %in% c("id_ind", "parent_origin")))
  close_pop(pop)
})


test_that("genome_genotype is also flagged as wide", {
  pop  <- make_summary_pop()
  summ <- summary(pop)
  geno <- summ$tables[["genome_genotype"]]
  expect_true(geno$is_wide_genome)
  expect_gte(geno$n_locus_cols, 100L)
  close_pop(pop)
})


test_that("low-cardinality VARCHAR column gets frequency table display", {
  pop  <- make_summary_pop()
  summ <- summary(pop)
  cs   <- summ$tables[["ind_meta"]]$col_summaries
  sex_display <- cs$display[cs$col_name == "sex"]
  expect_match(sex_display, "F:")
  expect_match(sex_display, "M:")
  close_pop(pop)
})


test_that("DOUBLE column gets 5-number summary display", {
  pop <- make_summary_pop()
  pop <- pop |>
    add_trait("ADG", target_add_var = 1, residual_var = 1) |>
    define_qtl("ADG", n = 20) |>
    set_qtl_effects("ADG", seed = 1)
  pop <- pop |> get_table("ind_meta") |> add_phenotype("ADG")
  summ    <- summary(pop, tables = "ind_phenotype")
  cs      <- summ$tables[["ind_phenotype"]]$col_summaries
  val_row <- cs[cs$col_name == "value", ]
  expect_match(val_row$display, "Min:")
  expect_match(val_row$display, "Q1:")
  expect_match(val_row$display, "Med:")
  close_pop(pop)
})


test_that("tables parameter restricts output to named subset", {
  pop  <- make_summary_pop()
  summ <- summary(pop, tables = c("ind_meta", "genome_meta"))
  expect_equal(sort(names(summ$tables)), sort(c("ind_meta", "genome_meta")))
  close_pop(pop)
})


test_that("unknown table name in tables raises a warning", {
  pop <- make_summary_pop()
  expect_warning(
    summary(pop, tables = c("ind_meta", "does_not_exist")),
    "not found"
  )
  close_pop(pop)
})


test_that("empty table shows 0 rows with no per-column breakdown", {
  pop  <- initialize_genome(
    pop_name = "EmptyTest",
    n_loci = 50,
    n_chr = 1,
    chr_len_Mb = 50,
    db_path = ":memory:"
  )
  summ <- summary(pop, tables = "ind_meta")
  meta <- summ$tables[["ind_meta"]]
  expect_equal(meta$n_rows, 0L)
  expect_equal(nrow(meta$col_summaries), 0L)
  close_pop(pop)
})


test_that("max_values = 1L forces 5-number summary for 2-unique INTEGER column", {
  # Build a pop with gen = 0 and gen = 1 so gen has 2 unique INTEGER values
  pop <- initialize_genome(
    pop_name = "SumTest2",
    n_loci = 50,
    n_chr = 1,
    chr_len_Mb = 50,
    n_haplotypes = 20,
    db_path = ":memory:",
    fixed_allele_freq = 0.5
  )
  pop <- add_founders(pop, n_males = 5, n_females = 5, line_name = "A", gen = 0L)
  pop <- add_founders(pop, n_males = 5, n_females = 5, line_name = "A", gen = 1L)

  # With max_values = 2, gen (2 unique) gets a frequency table
  summ_freq <- summary(pop, tables = "ind_meta", max_values = 2L)
  cs_freq   <- summ_freq$tables[["ind_meta"]]$col_summaries
  gen_freq  <- cs_freq[cs_freq$col_name == "gen", ]
  expect_match(gen_freq$display, "0:")

  # With max_values = 1, gen (2 unique > 1) gets a 5-number summary
  summ_5num <- summary(pop, tables = "ind_meta", max_values = 1L)
  cs_5num   <- summ_5num$tables[["ind_meta"]]$col_summaries
  gen_5num  <- cs_5num[cs_5num$col_name == "gen", ]
  expect_match(gen_5num$display, "Min:")

  close_pop(pop)
})


test_that("column with NAs shows NA count in display", {
  pop  <- make_summary_pop()
  summ <- summary(pop, tables = "ind_meta")
  cs   <- summ$tables[["ind_meta"]]$col_summaries
  # id_parent_1 is NA for all founders
  parent_row <- cs[cs$col_name == "id_parent_1", ]
  expect_match(parent_row$display, "NAs")
  close_pop(pop)
})


test_that("print.tidybreed_summary produces output without error", {
  pop  <- make_summary_pop()
  summ <- summary(pop)
  expect_output(print(summ), "tidybreed Population")
  expect_output(print(summ), "ind_meta")
  close_pop(pop)
})
