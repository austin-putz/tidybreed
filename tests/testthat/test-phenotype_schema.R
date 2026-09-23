# Phase 1 of plans/sample_correlated_effects.md: the residual-realization
# columns live in the base ind_phenotype DDL, liability_value / cat_name are no
# longer added by on-demand ALTER TABLE, and phenotype_meta carries
# condition_change_action.

IND_PHENOTYPE_BASE_COLS <- c(
  "id_phenotype", "id_ind", "phenotype_name", "pheno_value", "pheno_number",
  "liability_value", "cat_name", "residual_value", "residual_condition_level"
)

test_that("ind_phenotype and phenotype_meta carry the new columns from open_pop()", {
  pop <- open_pop(pop_name = "ph_schema", db_name = ":memory:")
  on.exit(close_pop(pop), add = TRUE)

  expect_identical(DBI::dbListFields(pop$db_conn, "ind_phenotype"),
                   IND_PHENOTYPE_BASE_COLS)
  expect_true("condition_change_action" %in%
                DBI::dbListFields(pop$db_conn, "phenotype_meta"))

  # Every base column is reserved and described.
  expect_true(all(IND_PHENOTYPE_BASE_COLS %in% TABLE_RESERVED_COLS$ind_phenotype))
  expect_true("condition_change_action" %in% TABLE_RESERVED_COLS$phenotype_meta)
  d <- describe_table(pop, "ind_phenotype")
  expect_true(all(IND_PHENOTYPE_BASE_COLS %in% d$column_name))
  expect_true(all(nchar(d$description[d$column_name %in% IND_PHENOTYPE_BASE_COLS]) > 0L))
  dm <- describe_table(pop, "phenotype_meta")
  expect_true(nchar(dm$description[dm$column_name == "condition_change_action"]) > 0L)
})

test_that("store_liability and cat_names populate base columns without ALTER TABLE", {
  set.seed(11)
  pop <- make_test_pop(n_males = 50, n_females = 50, n_loci = 100, n_chr = 2)
  on.exit(close_pop(pop), add = TRUE)
  pop <- define_trait(pop, "mort", target_add_var = 1)
  pop <- pop |> get_table("genome_meta") |> define_additive_effects("mort")
  pop <- define_phenotype(pop, "mort",
                          type            = "categorical",
                          prevalence      = 0.2,
                          cat_values      = c(0, 1),
                          cat_names       = c("Alive", "Dead"),
                          store_liability = TRUE,
                          residual_var    = 1)

  cols_before <- DBI::dbListFields(pop$db_conn, "ind_phenotype")
  pop <- pop |> get_table("ind_meta") |> add_phenotype("mort")
  cols_after  <- DBI::dbListFields(pop$db_conn, "ind_phenotype")

  # No ALTER TABLE: the column set is unchanged by the write.
  expect_identical(cols_after, cols_before)

  ph <- dplyr::collect(get_table(pop, "ind_phenotype"))
  expect_equal(nrow(ph), 100L)
  expect_false(anyNA(ph$liability_value))
  expect_true(all(ph$cat_name %in% c("Alive", "Dead")))
  # Liability and category agree: category 1 ("Dead") is the upper tail.
  expect_true(all(ph$liability_value[ph$pheno_value == 1] >
                    max(ph$liability_value[ph$pheno_value == 0])))
  # The liability-scale residual is stored; the unconditional R has no level.
  expect_false(anyNA(ph$residual_value))
  expect_true(all(is.na(ph$residual_condition_level)))
})

test_that("a continuous phenotype leaves liability_value and cat_name NULL", {
  set.seed(12)
  pop <- make_test_pop(n_males = 10, n_females = 10, n_loci = 100, n_chr = 2)
  on.exit(close_pop(pop), add = TRUE)
  pop <- define_trait(pop, "ADG", target_add_var = 100)
  pop <- pop |> get_table("genome_meta") |> define_additive_effects("ADG")
  pop <- define_phenotype(pop, "ADG", mean = 500, residual_var = 50)
  pop <- pop |> get_table("ind_meta") |> add_phenotype("ADG")

  ph <- dplyr::collect(get_table(pop, "ind_phenotype"))
  expect_equal(nrow(ph), 20L)
  expect_true(all(is.na(ph$liability_value)))
  expect_true(all(is.na(ph$cat_name)))
  # The base columns are reserved: a user cannot shadow them.
  expect_error(pop |> get_table("ind_phenotype") |> mutate_table(residual_value = 1),
               "reserved")
})

test_that("define_phenotype() stores condition_change_action and validates it", {
  pop <- open_pop(pop_name = "ph_cca", db_name = ":memory:")
  on.exit(close_pop(pop), add = TRUE)

  pop <- define_phenotype(pop, "A", residual_var = 1)
  pop <- define_phenotype(pop, "B", residual_var = 1,
                          condition_change_action = "independent")
  expect_error(define_phenotype(pop, "C", condition_change_action = "warn"),
               "should be one of")

  pm <- DBI::dbGetQuery(pop$db_conn,
    "SELECT phenotype_name, condition_change_action FROM phenotype_meta ORDER BY phenotype_name")
  expect_identical(pm$condition_change_action, c("error", "independent"))

  # overwrite = TRUE replaces the stored value like any other phenotype_meta field
  pop <- define_phenotype(pop, "A", condition_change_action = "independent",
                          overwrite = TRUE)
  v <- DBI::dbGetQuery(pop$db_conn,
    "SELECT condition_change_action FROM phenotype_meta WHERE phenotype_name = 'A'")
  expect_identical(v$condition_change_action, "independent")
})
