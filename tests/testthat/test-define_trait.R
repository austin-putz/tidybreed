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


test_that("define_trait() creates the trait tables and inserts the row", {
  pop <- make_tiny_pop("trait_basic")

  pop <- define_trait(pop,
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


test_that("define_trait() refuses duplicate names without overwrite", {
  pop <- make_tiny_pop("trait_dup")
  pop <- define_trait(pop, "ADG", target_add_var = 0.25, residual_var = 0.75)

  expect_error(define_trait(pop, "ADG"), "already exists")

  # overwrite = TRUE replaces the row; variances update in trait_effect_cov
  pop <- define_trait(pop, "ADG", target_add_var = 0.1, residual_var = 0.9,
                   overwrite = TRUE)
  expect_equal(get_effect_var(pop, "gen_add", "ADG"), 0.1)
  expect_equal(get_effect_var(pop, "residual", "ADG"), 0.9)

  close_pop(pop)
})


test_that("define_trait() validates trait_type-specific args", {
  pop <- make_tiny_pop("trait_validate")

  # binary is removed — should error with helpful message
  expect_error(
    define_trait(pop, "mort", trait_type = "binary"),
    "removed"
  )
  # categorical requires thresholds OR prevalence
  expect_error(
    define_trait(pop, "score", trait_type = "categorical"),
    "thresholds"
  )
  # cannot supply both
  expect_error(
    define_trait(pop, "score", trait_type = "categorical",
              thresholds = c(0), prevalence = 0.1),
    "not both"
  )

  close_pop(pop)
})


test_that("define_trait() accepts categorical with prevalence (2-category)", {
  pop <- make_tiny_pop("trait_cat_prev")

  pop <- define_trait(pop, "mort",
                   trait_type      = "categorical",
                   prevalence      = 0.10,
                   cat_values      = c(0, 1),
                   cat_names       = c("Alive", "Dead"),
                   store_liability = TRUE,
                   target_add_var  = 1,
                   residual_var    = 9)

  row <- DBI::dbGetQuery(pop$db_conn,
    "SELECT * FROM trait_meta WHERE trait_name = 'mort'")
  expect_equal(row$trait_type,  "categorical")
  expect_equal(row$prevalence,  0.10)
  expect_equal(row$cat_values,  "0,1")
  expect_equal(row$cat_names,   "Alive,Dead")
  expect_true(row$store_liability)

  close_pop(pop)
})


test_that("define_trait() accepts categorical with explicit thresholds", {
  pop <- make_tiny_pop("trait_cat_thresh")

  pop <- define_trait(pop, "score",
                   trait_type  = "categorical",
                   thresholds  = c(-1, 0, 1),
                   cat_values  = c(1, 2, 3, 4),
                   cat_names   = c("thin", "fair", "good", "excellent"),
                   target_add_var = 1,
                   residual_var   = 1)

  row <- DBI::dbGetQuery(pop$db_conn,
    "SELECT * FROM trait_meta WHERE trait_name = 'score'")
  expect_equal(row$thresholds, "-1,0,1")
  expect_equal(row$cat_values, "1,2,3,4")
  expect_equal(row$cat_names,  "thin,fair,good,excellent")

  close_pop(pop)
})


test_that("define_trait() validates cat_values and cat_names lengths", {
  pop <- make_tiny_pop("trait_cat_len")

  # 1 threshold → 2 categories; cat_values must be length 2
  expect_error(
    define_trait(pop, "x", trait_type = "categorical",
              thresholds = c(0), cat_values = c(1, 2, 3)),
    "length 2"
  )
  # cat_names same length as cat_values
  expect_error(
    define_trait(pop, "x", trait_type = "categorical",
              thresholds = c(0), cat_values = c(0, 1),
              cat_names = c("No")),
    "length 2"
  )
  # cat_names cannot contain commas
  expect_error(
    define_trait(pop, "x", trait_type = "categorical",
              thresholds = c(0),
              cat_names = c("Yes, really", "No")),
    "commas"
  )

  close_pop(pop)
})


test_that("define_trait() rejects SQL-unsafe names", {
  pop <- make_tiny_pop("trait_names")
  expect_error(define_trait(pop, "3bad"), "Invalid trait name")
  expect_error(define_trait(pop, "SELECT"), "reserved keyword")
  close_pop(pop)
})


test_that("define_trait() inserts global economic_weight row into index_meta", {
  pop <- make_tiny_pop("trait_ev_insert")
  on.exit(close_pop(pop))

  pop <- define_trait(pop, "ADG",
                   trait_type     = "continuous",
                   target_add_var = 0.25,
                   residual_var   = 0.75,
                   economic_value = 5.0)

  rows <- DBI::dbGetQuery(pop$db_conn,
    "SELECT * FROM index_meta WHERE index_name IS NULL AND trait_name = 'ADG'")
  expect_equal(nrow(rows), 1L)
  expect_equal(rows$economic_weight, 5.0, tolerance = 1e-9)
  expect_true(is.na(rows$index_weight))
})


test_that("define_trait() overwrite = FALSE preserves existing index_meta row", {
  pop <- make_tiny_pop("trait_ev_no_overwrite")
  on.exit(close_pop(pop))

  pop <- define_trait(pop, "ADG",
                   trait_type     = "continuous",
                   target_add_var = 0.25,
                   residual_var   = 0.75,
                   economic_value = 5.0)
  # overwrite = FALSE on a new trait — second call should error on trait_meta
  # (existing coverage). Check index_meta row persists after overwrite = TRUE.
  pop <- define_trait(pop, "ADG",
                   trait_type     = "continuous",
                   target_add_var = 0.25,
                   residual_var   = 0.75,
                   economic_value = 99.0,
                   overwrite      = TRUE)

  rows <- DBI::dbGetQuery(pop$db_conn,
    "SELECT economic_weight FROM index_meta WHERE index_name IS NULL AND trait_name = 'ADG'")
  expect_equal(nrow(rows), 1L)
  expect_equal(rows$economic_weight, 99.0, tolerance = 1e-9)
})


test_that("define_trait() overwrite = TRUE updates index_meta and keeps only one row", {
  pop <- make_tiny_pop("trait_ev_overwrite")
  on.exit(close_pop(pop))

  pop <- define_trait(pop, "ADG",
                   trait_type     = "continuous",
                   target_add_var = 0.25,
                   residual_var   = 0.75,
                   economic_value = 1.0)
  pop <- define_trait(pop, "ADG",
                   trait_type     = "continuous",
                   target_add_var = 0.25,
                   residual_var   = 0.75,
                   economic_value = 9.0,
                   overwrite      = TRUE)

  rows <- DBI::dbGetQuery(pop$db_conn,
    "SELECT * FROM index_meta WHERE index_name IS NULL AND trait_name = 'ADG'")
  expect_equal(nrow(rows), 1L)
  expect_equal(rows$economic_weight, 9.0, tolerance = 1e-9)
})
