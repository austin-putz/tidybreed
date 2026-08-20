# Guards against drift between the schema registries in R/sql_utils.R
# (TABLE_RESERVED_COLS, TABLE_ROW_KEYS, TABLE_PRIMARY_KEYS) and the columns the
# tables actually have. Drift here fails at a distance: mutate_table() silently
# stops protecting a renamed column, and remove_rows() builds SQL against a
# column that no longer exists.


# Columns deliberately listed before they exist. These are added later by
# ALTER TABLE, and are reserved up front so a user cannot create a conflicting
# column of the same name in the meantime.
#   replicate                 — stamped by archive_replicate()
#   liability_value, cat_name — added by add_phenotype() per phenotype type
DEFERRED_COLS <- list(
  ind_meta       = "replicate",
  ind_phenotype  = c("liability_value", "cat_name", "replicate"),
  ind_tbv        = "replicate",
  ind_ebv        = "replicate",
  ind_index      = "replicate",
  ind_true_index = "replicate"
)

# Population exercising every table the registries reference.
make_pop_all_tables <- function() {
  pop <- make_test_pop(n_males = 2, n_females = 2, n_loci = 50, n_chr = 2)

  pop <- define_trait(pop, "ADG", target_add_var = 100)
  pop <- pop |> get_table("genome_meta") |> define_additive_effects("ADG")
  pop <- define_phenotype(pop, "ADG", mean = 500, residual_var = 50)

  pop <- define_effect_fixed_class(pop, "ADG", effect_name = "sex",
                                   source_column = "sex",
                                   levels = c(M = 10, F = 0))
  pop <- define_effect_random(pop, "ADG", effect_name = "pen",
                              source_column = "line_name", variance = 5)
  pop <- define_index(pop, "IDX", trait_names = "ADG", index_wts = 1)

  pop |> get_table("ind_meta") |> add_phenotype("ADG")
}


test_that("every column in TABLE_RESERVED_COLS exists on its table", {
  pop <- make_pop_all_tables()
  on.exit(close_pop(pop), add = TRUE)
  live <- DBI::dbListTables(pop$db_conn)

  for (tbl in intersect(names(TABLE_RESERVED_COLS), live)) {
    actual   <- DBI::dbListFields(pop$db_conn, tbl)
    expected <- setdiff(TABLE_RESERVED_COLS[[tbl]], DEFERRED_COLS[[tbl]])
    expect_setequal(intersect(expected, actual), expected)
  }
})

test_that("every column in TABLE_ROW_KEYS and TABLE_PRIMARY_KEYS exists on its table", {
  pop <- make_pop_all_tables()
  on.exit(close_pop(pop), add = TRUE)
  live <- DBI::dbListTables(pop$db_conn)

  for (tbl in intersect(names(TABLE_ROW_KEYS), live)) {
    actual <- DBI::dbListFields(pop$db_conn, tbl)
    expect_setequal(intersect(TABLE_ROW_KEYS[[tbl]], actual), TABLE_ROW_KEYS[[tbl]])
  }
  for (tbl in intersect(names(TABLE_PRIMARY_KEYS), live)) {
    actual <- DBI::dbListFields(pop$db_conn, tbl)
    expect_true(TABLE_PRIMARY_KEYS[[tbl]] %in% actual)
  }
})

test_that("trait_effects is registered against phenotype_name, not trait_name", {
  # Regression: both registries listed trait_name after the column was renamed
  # to phenotype_name, so remove_rows() on trait_effects could never work.
  expect_true("phenotype_name" %in% TABLE_RESERVED_COLS$trait_effects)
  expect_false("trait_name"    %in% TABLE_RESERVED_COLS$trait_effects)
  expect_identical(TABLE_ROW_KEYS$trait_effects,
                   c("phenotype_name", "effect_name"))
})

test_that("TABLE_RESERVED_COLS$trait_effects covers the whole table", {
  pop <- make_pop_all_tables()
  on.exit(close_pop(pop), add = TRUE)
  expect_setequal(DBI::dbListFields(pop$db_conn, "trait_effects"),
                  TABLE_RESERVED_COLS$trait_effects)
})

test_that("remove_rows() deletes an exact row from trait_effects", {
  # This is the user-visible failure the stale row key caused.
  pop <- make_pop_all_tables()
  on.exit(close_pop(pop), add = TRUE)

  before <- get_table(pop, "trait_effects") |> dplyr::collect()
  expect_setequal(before$effect_name, c("sex", "pen"))

  pop |>
    get_table("trait_effects") |>
    dplyr::filter(effect_name == "sex") |>
    remove_rows(verbose = FALSE)

  after <- get_table(pop, "trait_effects") |> dplyr::collect()
  expect_identical(after$effect_name, "pen")
  expect_identical(after$phenotype_name, "ADG")
})

test_that("mutate_table() blocks the renamed trait_effects reserved columns", {
  pop <- make_pop_all_tables()
  on.exit(close_pop(pop), add = TRUE)

  expect_error(
    pop |> get_table("trait_effects") |> mutate_table(phenotype_name = "X"),
    regexp = "reserved"
  )
  expect_error(
    pop |> get_table("trait_effects") |> mutate_table(null_class_action = "error"),
    regexp = "reserved"
  )
})
