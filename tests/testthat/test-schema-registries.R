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

test_that("phenotype_effects is registered against phenotype_name, not trait_name", {
  # Regression: both registries listed trait_name after the column was renamed
  # to phenotype_name, so remove_rows() on phenotype_effects could never work.
  expect_true("phenotype_name" %in% TABLE_RESERVED_COLS$phenotype_effects)
  expect_false("trait_name"    %in% TABLE_RESERVED_COLS$phenotype_effects)
  expect_identical(TABLE_ROW_KEYS$phenotype_effects,
                   c("phenotype_name", "effect_name"))
})

test_that("TABLE_RESERVED_COLS$phenotype_effects covers the whole table", {
  pop <- make_pop_all_tables()
  on.exit(close_pop(pop), add = TRUE)
  expect_setequal(DBI::dbListFields(pop$db_conn, "phenotype_effects"),
                  TABLE_RESERVED_COLS$phenotype_effects)
})

test_that("remove_rows() deletes an exact row from phenotype_effects", {
  # This is the user-visible failure the stale row key caused.
  pop <- make_pop_all_tables()
  on.exit(close_pop(pop), add = TRUE)

  before <- get_table(pop, "phenotype_effects") |> dplyr::collect()
  expect_setequal(before$effect_name, c("sex", "pen"))

  pop |>
    get_table("phenotype_effects") |>
    dplyr::filter(effect_name == "sex") |>
    remove_rows(verbose = FALSE)

  after <- get_table(pop, "phenotype_effects") |> dplyr::collect()
  expect_identical(after$effect_name, "pen")
  expect_identical(after$phenotype_name, "ADG")
})

test_that("mutate_table() blocks the renamed phenotype_effects reserved columns", {
  pop <- make_pop_all_tables()
  on.exit(close_pop(pop), add = TRUE)

  expect_error(
    pop |> get_table("phenotype_effects") |> mutate_table(phenotype_name = "X"),
    regexp = "reserved"
  )
  expect_error(
    pop |> get_table("phenotype_effects") |> mutate_table(null_class_action = "error"),
    regexp = "reserved"
  )
})


test_that("package-managed tables reserve every column they have", {
  # Regression: founder_haplotypes, phenotype_components, and
  # phenotype_random_effects had no TABLE_RESERVED_COLS entry at all, so
  # mutate_table() blocked nothing on them and a user could overwrite e.g.
  # phenotype_random_effects.draw_value or phenotype_components.contributor_type
  # in place. Tables written exclusively by a define_*/add_* function reserve
  # all of their columns; only entity-shaped tables leave room for user columns.
  pop <- make_pop_all_tables()
  on.exit(close_pop(pop), add = TRUE)

  managed <- c("founder_haplotypes", "phenotype_components",
               "phenotype_random_effects", "phenotype_effects", "genome_effects",
               "phenotype_meta", "ind_true_index")

  for (tbl in managed) {
    expect_true(tbl %in% names(TABLE_RESERVED_COLS), info = tbl)
    expected <- setdiff(TABLE_RESERVED_COLS[[tbl]], DEFERRED_COLS[[tbl]])
    expect_setequal(DBI::dbListFields(pop$db_conn, tbl), expected)
  }
})

test_that("mutate_table() blocks writes to the newly registered tables", {
  pop <- make_pop_all_tables()
  on.exit(close_pop(pop), add = TRUE)

  expect_error(
    pop |> get_table("founder_haplotypes") |> mutate_table(allele = 1L),
    regexp = "reserved"
  )
  expect_error(
    pop |> get_table("phenotype_random_effects") |> mutate_table(draw_value = 1),
    regexp = "reserved"
  )
  expect_error(
    pop |> get_table("phenotype_components") |> mutate_table(contributor_type = "dam"),
    regexp = "reserved"
  )
})


test_that("every system table either has a row key or is explicitly excluded", {
  # Completeness, not just correctness: a missing TABLE_ROW_KEYS entry made
  # remove_rows() hard-error on phenotype_meta, phenotype_components,
  # phenotype_random_effects and founder_haplotypes. Adding a new table must now
  # force a decision — register a key, or say why deletion is refused.
  expect_setequal(
    c(names(TABLE_ROW_KEYS), names(TABLE_NO_ROW_DELETE)),
    SYSTEM_TABLES
  )
  # A table cannot be both deletable and refused.
  expect_length(intersect(names(TABLE_ROW_KEYS), names(TABLE_NO_ROW_DELETE)), 0L)
  # Every refusal carries a reason for the user.
  expect_true(all(nzchar(TABLE_NO_ROW_DELETE)))
})

test_that("remove_rows() refuses _schema_meta with a reason, not a generic error", {
  pop <- make_pop_all_tables()
  on.exit(close_pop(pop), add = TRUE)

  expect_error(
    pop |>
      get_table("_schema_meta") |>
      dplyr::filter(table_name == "ind_meta") |>
      remove_rows(verbose = FALSE),
    regexp = "define_schema_description"
  )
})

test_that("remove_rows() deletes from the newly registered observation tables", {
  pop <- make_pop_all_tables()
  on.exit(close_pop(pop), add = TRUE)

  pop |>
    get_table("phenotype_random_effects") |>
    dplyr::filter(effect_name == "pen") |>
    remove_rows(verbose = FALSE)
  expect_equal(
    nrow(dplyr::collect(get_table(pop, "phenotype_random_effects"))), 0L
  )

  pop |>
    get_table("phenotype_meta") |>
    dplyr::filter(phenotype_name == "ADG") |>
    remove_rows(verbose = FALSE)
  expect_equal(nrow(dplyr::collect(get_table(pop, "phenotype_meta"))), 0L)
})

test_that("remove_rows() deletes rows whose key columns are NULL", {
  # Regression: delete_exact_rows() joined with `=`, and `NULL = NULL` is NULL,
  # so the default chr_inheritance / chr_recombination rows seeded by
  # define_genome() (offspring_sex and line_name both NULL) matched nothing.
  # The delete reported "Deleted 0 rows" as a success while changing nothing.
  pop <- make_pop_all_tables()
  on.exit(close_pop(pop), add = TRUE)

  before <- get_table(pop, "chr_inheritance") |> dplyr::collect()
  expect_true(all(is.na(before$offspring_sex)))
  expect_gt(nrow(before), 1L)
  target <- before$chr_name[[1L]]

  pop |>
    get_table("chr_inheritance") |>
    dplyr::filter(chr_name == target) |>
    remove_rows(verbose = FALSE)

  after <- get_table(pop, "chr_inheritance") |> dplyr::collect()
  expect_equal(nrow(after), nrow(before) - 1L)
  expect_false(target %in% after$chr_name)
})

test_that("remove_rows() deletes a whole line's founder_haplotypes pool", {
  # The shared pool has line_name IS NULL, so this also exercises the NULL-safe
  # join on a table whose key is composite and partly nullable.
  pop <- make_pop_all_tables()
  on.exit(close_pop(pop), add = TRUE)

  before <- get_table(pop, "founder_haplotypes") |> dplyr::collect()
  expect_gt(nrow(before), 0L)
  expect_true(all(is.na(before$line_name)))

  pop |>
    get_table("founder_haplotypes") |>
    dplyr::filter(haplotype_id == 1L) |>
    remove_rows(verbose = FALSE)

  after <- get_table(pop, "founder_haplotypes") |> dplyr::collect()
  expect_false(1L %in% after$haplotype_id)
  expect_equal(nrow(after), nrow(before) - sum(before$haplotype_id == 1L))
})
