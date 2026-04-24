make_tiny_pop <- function(pop_name = "t") {
  pop <- initialize_genome(
    pop_name        = pop_name,
    n_loci          = 500,
    n_chr           = 5,
    chr_len_Mb      = 100,
    n_haplotypes    = 100,
    db_path         = ":memory:",
    fixed_allele_freq = 0.5
  )
  pop <- add_founders(pop, n_males = 25, n_females = 25, line_name = "A")
  pop <- get_table(pop, "ind_meta") |> mutate_table(gen = 0L)
  pop
}


test_that("add_trait() creates the trait tables and inserts the row", {
  pop <- make_tiny_pop("trait_basic")

  pop <- add_trait(pop,
                   trait_name      = "ADG",
                   units           = "g/day",
                   trait_type      = "continuous",
                   target_add_var  = 0.25,
                   residual_var    = 0.75,
                   target_add_mean = 850)

  tables <- DBI::dbListTables(pop$db_conn)
  for (tbl in c("trait_meta", "trait_effects", "trait_effect_cov",
                "ind_phenotype", "ind_tbv", "ind_ebv")) {
    expect_true(tbl %in% tables)
  }

  row <- DBI::dbGetQuery(pop$db_conn,
    "SELECT * FROM trait_meta WHERE trait_name = 'ADG'")
  expect_equal(nrow(row), 1)
  expect_equal(row$trait_type, "continuous")
  expect_equal(row$target_add_mean, 850)

  # Variances are stored in trait_effect_cov
  expect_equal(get_effect_var(pop, "gen_add", "ADG"), 0.25)
  expect_equal(get_effect_var(pop, "residual", "ADG"), 0.75)

  close_pop(pop)
})


test_that("add_trait() refuses duplicate names without overwrite", {
  pop <- make_tiny_pop("trait_dup")
  pop <- add_trait(pop, "ADG", target_add_var = 0.25, residual_var = 0.75)

  expect_error(add_trait(pop, "ADG"), "already exists")

  # overwrite = TRUE replaces the row; variances update in trait_effect_cov
  pop <- add_trait(pop, "ADG", target_add_var = 0.1, residual_var = 0.9,
                   overwrite = TRUE)
  expect_equal(get_effect_var(pop, "gen_add", "ADG"), 0.1)
  expect_equal(get_effect_var(pop, "residual", "ADG"), 0.9)

  close_pop(pop)
})


test_that("add_trait() validates trait_type-specific args", {
  pop <- make_tiny_pop("trait_validate")

  # binary requires prevalence
  expect_error(
    add_trait(pop, "mort", trait_type = "binary"),
    "prevalence"
  )
  # categorical requires thresholds
  expect_error(
    add_trait(pop, "score", trait_type = "categorical"),
    "thresholds"
  )

  close_pop(pop)
})


test_that("add_trait() rejects SQL-unsafe names", {
  pop <- make_tiny_pop("trait_names")
  expect_error(add_trait(pop, "3bad"), "Invalid trait name")
  expect_error(add_trait(pop, "SELECT"), "reserved keyword")
  close_pop(pop)
})
