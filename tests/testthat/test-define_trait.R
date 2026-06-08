make_tiny_pop <- function(pop_name = "t") {
  pop <- initialize_genome(
    pop_name   = pop_name,
    n_loci     = 500,
    n_chr      = 5,
    chr_len_Mb = 100,
    db_path    = ":memory:"
  ) |>
    define_founder_haplotypes(n_haplotypes = 100, fixed_allele_freq = 0.5)
  pop <- add_founders(pop, n_males = 25, n_females = 25, line_name = "A")
  pop <- get_table(pop, "ind_meta") |> mutate_table(gen = 0L)
  pop
}


test_that("define_trait() creates the trait tables and inserts the row", {
  pop <- make_tiny_pop("trait_basic")

  pop <- define_trait(pop,
                      trait_name     = "ADG",
                      units          = "g/day",
                      target_add_var = 0.25,
                      target_add_mean = 850)

  tables <- DBI::dbListTables(pop$db_conn)
  for (tbl in c("trait_meta", "trait_effects", "trait_effect_cov",
                "ind_phenotype", "ind_tbv", "ind_ebv",
                "phenotype_meta", "phenotype_components",
                "phenotype_residual_cov")) {
    expect_true(tbl %in% tables, info = paste("missing table:", tbl))
  }

  row <- DBI::dbGetQuery(pop$db_conn,
    "SELECT * FROM trait_meta WHERE trait_name = 'ADG'")
  expect_equal(nrow(row), 1)
  expect_equal(row$units, "g/day")
  expect_equal(row$target_add_mean, 850)

  # Additive genetic variance stored in trait_effect_cov
  expect_equal(get_effect_var(pop, "gen_add", "ADG"), 0.25)

  close_pop(pop)
})


test_that("define_trait() refuses duplicate names without overwrite", {
  pop <- make_tiny_pop("trait_dup")
  pop <- define_trait(pop, "ADG", target_add_var = 0.25)

  expect_error(define_trait(pop, "ADG"), "already exists")

  pop <- define_trait(pop, "ADG", target_add_var = 0.1, overwrite = TRUE)
  expect_equal(get_effect_var(pop, "gen_add", "ADG"), 0.1)

  close_pop(pop)
})


test_that("define_trait() rejects SQL-unsafe names", {
  pop <- make_tiny_pop("trait_names")
  expect_error(define_trait(pop, "3bad"), "Invalid trait name")
  expect_error(define_trait(pop, "SELECT"), "reserved keyword")
  close_pop(pop)
})


test_that("define_trait() stores target_add_mean", {
  pop <- make_tiny_pop("trait_mean")
  pop <- define_trait(pop, "ADG", target_add_var = 1.0, target_add_mean = 500)

  row <- DBI::dbGetQuery(pop$db_conn,
    "SELECT target_add_mean FROM trait_meta WHERE trait_name = 'ADG'")
  expect_equal(row$target_add_mean, 500)

  close_pop(pop)
})


test_that("define_trait() stores expressed_parent", {
  pop <- make_tiny_pop("trait_ep")
  pop <- define_trait(pop, "mat", target_add_var = 1.0,
                      expressed_parent = "parent_2")

  row <- DBI::dbGetQuery(pop$db_conn,
    "SELECT expressed_parent FROM trait_meta WHERE trait_name = 'mat'")
  expect_equal(row$expressed_parent, "parent_2")

  close_pop(pop)
})


test_that("define_trait() inserts a global index_meta row", {
  pop <- make_tiny_pop("trait_idx_row")
  on.exit(close_pop(pop))

  pop <- define_trait(pop, "ADG", target_add_var = 0.25)

  rows <- DBI::dbGetQuery(pop$db_conn,
    "SELECT * FROM index_meta WHERE index_name IS NULL AND trait_name = 'ADG'")
  expect_equal(nrow(rows), 1L)
})


test_that("define_trait() overwrite = TRUE replaces trait_meta row", {
  pop <- make_tiny_pop("trait_overwrite")
  on.exit(close_pop(pop))

  pop <- define_trait(pop, "ADG", target_add_var = 0.25, units = "g/day")
  pop <- define_trait(pop, "ADG", target_add_var = 0.10, units = "kg/day",
                      overwrite = TRUE)

  row <- DBI::dbGetQuery(pop$db_conn,
    "SELECT * FROM trait_meta WHERE trait_name = 'ADG'")
  expect_equal(nrow(row), 1L)
  expect_equal(row$units, "kg/day")
  expect_equal(get_effect_var(pop, "gen_add", "ADG"), 0.10)
})


test_that("define_trait() can be called without target_add_var", {
  pop <- make_tiny_pop("trait_no_var")
  on.exit(close_pop(pop))

  pop <- define_trait(pop, "latent", description = "no QTL yet")
  row <- DBI::dbGetQuery(pop$db_conn,
    "SELECT * FROM trait_meta WHERE trait_name = 'latent'")
  expect_equal(nrow(row), 1L)
})
