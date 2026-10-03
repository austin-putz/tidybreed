make_tiny_pop <- function(pop_name = "t") {
  pop <- open_pop(pop_name = pop_name, db_name = ":memory:") |>
    define_genome(n_loci = 500, n_chr = 5, chr_len_Mb = 100) |>
    define_founder_haplotypes(n_haplotypes = 100, method = "fixed")
  pop <- pop |>
    get_table("founder_haplotypes") |>
    add_founders( n_males = 25, n_females = 25, line_name = "A")
  pop <- get_table(pop, "ind_meta") |> mutate_table(gen = 0L)
  pop
}


test_that("define_trait() creates the trait tables and inserts the row", {
  pop <- make_tiny_pop("trait_basic")

  pop <- define_trait(pop,
                      trait_name     = "ADG",
                      units          = "g/day")

  tables <- DBI::dbListTables(pop$db_conn)
  for (tbl in c("trait_meta", "phenotype_effects", "trait_var_comp",
                "ind_phenotype", "ind_tbv", "ind_ebv",
                "phenotype_meta", "phenotype_components",
                "phenotype_var_comp")) {
    expect_true(tbl %in% tables, info = paste("missing table:", tbl))
  }

  row <- DBI::dbGetQuery(pop$db_conn,
    "SELECT * FROM trait_meta WHERE trait_name = 'ADG'")
  expect_equal(nrow(row), 1)
  expect_equal(row$units, "g/day")

  # define_trait() writes no target: targets enter only through
  # trait_var_comp's writers (define_effect_cov_matrix(), generator G =).
  expect_true(is.na(get_trait_var(pop, "additive", "ADG")))

  close_pop(pop)
})


test_that("define_trait() refuses duplicate names without overwrite", {
  pop <- make_tiny_pop("trait_dup")
  pop <- with_additive_target(pop, "ADG", 0.25)

  expect_error(define_trait(pop, "ADG"), "already exists")

  # overwrite replaces the trait_meta row and leaves the stored target alone.
  pop <- define_trait(pop, "ADG", overwrite = TRUE)
  expect_equal(get_trait_var(pop, "additive", "ADG"), 0.25)

  close_pop(pop)
})


test_that("define_trait() rejects SQL-unsafe names", {
  pop <- make_tiny_pop("trait_names")
  expect_error(define_trait(pop, "3bad"), "Invalid trait name")
  expect_error(define_trait(pop, "SELECT"), "reserved keyword")
  close_pop(pop)
})


test_that("define_trait() has no target arguments and trait_meta no target column", {
  # 0.73.0: targets enter only through trait_var_comp (one entry path, §6C);
  # target_add_mean was written and never read.
  pop <- make_tiny_pop("trait_no_target")
  on.exit(close_pop(pop))
  expect_false(any(c("target_add_var", "target_add_mean") %in%
                     names(formals(define_trait))))
  expect_error(define_trait(pop, "ADG", target_add_var = 1), "unused argument")
  pop <- define_trait(pop, "ADG")
  cols <- DBI::dbGetQuery(pop$db_conn, "SELECT * FROM trait_meta LIMIT 0")
  expect_false("target_add_mean" %in% names(cols))
  expect_false(exists("define_trait_simple", envir = asNamespace("tidybreed")))
})


test_that("trait_meta carries no expressed_parent column", {
  # Imprinting is a property of an effect's scope, not of a trait: it is one
  # origin row on the member (define_additive_effects(parent_origin = )), which
  # can differ per locus, per line and per effect owner. The trait-wide flag
  # could express none of that and is gone, argument and column both.
  pop <- make_tiny_pop("trait_ep")

  cols <- DBI::dbGetQuery(pop$db_conn, "SELECT * FROM trait_meta LIMIT 0")
  expect_false("expressed_parent" %in% names(cols))
  expect_false("expressed_parent" %in% names(formals(define_trait)))
  expect_error(define_trait(pop, "mat", expressed_parent = "parent_2"),
               "unused argument")

  close_pop(pop)
})


test_that("define_trait() inserts a global index_meta row", {
  pop <- make_tiny_pop("trait_idx_row")
  on.exit(close_pop(pop))

  pop <- with_additive_target(pop, "ADG", 0.25)

  rows <- DBI::dbGetQuery(pop$db_conn,
    "SELECT * FROM index_meta WHERE index_name IS NULL AND trait_name = 'ADG'")
  expect_equal(nrow(rows), 1L)
})


test_that("define_trait() overwrite = TRUE replaces trait_meta row", {
  pop <- make_tiny_pop("trait_overwrite")
  on.exit(close_pop(pop))

  pop <- with_additive_target(pop, "ADG", 0.25, units = "g/day")
  pop <- define_trait(pop, "ADG", units = "kg/day", overwrite = TRUE)

  row <- DBI::dbGetQuery(pop$db_conn,
    "SELECT * FROM trait_meta WHERE trait_name = 'ADG'")
  expect_equal(nrow(row), 1L)
  expect_equal(row$units, "kg/day")

  # Changing a stored target is an explicit remove-then-write, never an
  # overwrite.
  expect_error(suppressMessages(
    define_effect_cov_matrix(pop, "additive", 0.10, trait_name = "ADG")),
    "already stored")
  pop <- get_table(pop, "trait_var_comp") |>
    dplyr::filter(effect_name == "additive", trait_name_1 == "ADG") |>
    remove_rows()
  pop <- suppressMessages(
    define_effect_cov_matrix(pop, "additive", 0.10, trait_name = "ADG"))
  expect_equal(get_trait_var(pop, "additive", "ADG"), 0.10)
})


test_that("define_trait() can be called with only a name and description", {
  pop <- make_tiny_pop("trait_no_var")
  on.exit(close_pop(pop))

  pop <- define_trait(pop, "latent", description = "no QTL yet")
  row <- DBI::dbGetQuery(pop$db_conn,
    "SELECT * FROM trait_meta WHERE trait_name = 'latent'")
  expect_equal(nrow(row), 1L)
})
