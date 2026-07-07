test_that("_schema_meta table exists after open_pop() |> define_genome()", {
  pop <- open_pop(pop_name = "schema_test", db_name = ":memory:") |>
    define_genome(n_loci = 100, n_chr = 2, chr_len_Mb = 50)
  on.exit(close_pop(pop))

  expect_true("_schema_meta" %in% pop$tables)
  expect_true(DBI::dbExistsTable(pop$db_conn, "_schema_meta"))
})


test_that("_schema_meta is populated with table and column entries after open_pop() |> define_genome()", {
  pop <- open_pop(pop_name = "schema_test", db_name = ":memory:") |>
    define_genome(n_loci = 100, n_chr = 2, chr_len_Mb = 50)
  on.exit(close_pop(pop))

  n_tbl <- DBI::dbGetQuery(
    pop$db_conn,
    "SELECT COUNT(*) AS n FROM _schema_meta WHERE object_type = 'table'"
  )$n
  expect_gt(n_tbl, 0L)

  n_col <- DBI::dbGetQuery(
    pop$db_conn,
    "SELECT COUNT(*) AS n FROM _schema_meta WHERE object_type = 'column'"
  )$n
  expect_gt(n_col, 0L)
})


test_that("ensure_trait_tables() registers trait-layer descriptions in _schema_meta", {
  pop <- open_pop(pop_name = "schema_test", db_name = ":memory:") |>
    define_genome(n_loci = 100, n_chr = 2, chr_len_Mb = 50)
  on.exit(close_pop(pop))

  n_trait_tbl <- DBI::dbGetQuery(
    pop$db_conn,
    "SELECT COUNT(*) AS n FROM _schema_meta WHERE object_type = 'table' AND table_name = 'ind_tbv'"
  )$n
  expect_equal(n_trait_tbl, 1L)

  n_trait_col <- DBI::dbGetQuery(
    pop$db_conn,
    "SELECT COUNT(*) AS n FROM _schema_meta WHERE object_type = 'column' AND table_name = 'ind_tbv'"
  )$n
  expect_gt(n_trait_col, 0L)
})


test_that("schema() returns correct class and columns", {
  pop <- open_pop(pop_name = "schema_test", db_name = ":memory:") |>
    define_genome(n_loci = 100, n_chr = 2, chr_len_Mb = 50)
  on.exit(close_pop(pop))

  s <- schema(pop)
  expect_s3_class(s, "tidybreed_schema")
  expect_true(all(c("table_name", "n_rows", "n_cols", "description") %in% names(s)))
})


test_that("schema() covers all user-visible tables and excludes _schema_meta", {
  pop <- open_pop(pop_name = "schema_test", db_name = ":memory:") |>
    define_genome(n_loci = 100, n_chr = 2, chr_len_Mb = 50)
  on.exit(close_pop(pop))

  s <- schema(pop)
  expected_tables <- setdiff(pop$tables, "_schema_meta")
  expect_setequal(s$table_name, expected_tables)
  expect_false("_schema_meta" %in% s$table_name)
})


test_that("schema() includes descriptions for core tables", {
  pop <- open_pop(pop_name = "schema_test", db_name = ":memory:") |>
    define_genome(n_loci = 100, n_chr = 2, chr_len_Mb = 50)
  on.exit(close_pop(pop))

  s <- schema(pop)
  genome_meta_row <- s[s$table_name == "genome_meta", ]
  expect_equal(nrow(genome_meta_row), 1L)
  expect_true(nchar(genome_meta_row$description) > 0L)

  ind_tbv_row <- s[s$table_name == "ind_tbv", ]
  expect_equal(nrow(ind_tbv_row), 1L)
  expect_true(nchar(ind_tbv_row$description) > 0L)
})


test_that("describe_table() returns correct class and columns", {
  pop <- open_pop(pop_name = "schema_test", db_name = ":memory:") |>
    define_genome(n_loci = 100, n_chr = 2, chr_len_Mb = 50)
  on.exit(close_pop(pop))

  d <- describe_table(pop, "ind_meta")
  expect_s3_class(d, "tidybreed_table_desc")
  expect_true(all(c("column_name", "column_type", "description", "notes") %in% names(d)))
})


test_that("describe_table() includes core columns of ind_meta", {
  pop <- open_pop(pop_name = "schema_test", db_name = ":memory:") |>
    define_genome(n_loci = 100, n_chr = 2, chr_len_Mb = 50)
  on.exit(close_pop(pop))

  d <- describe_table(pop, "ind_meta")
  expect_true("id_ind" %in% d$column_name)
  expect_true("sex"    %in% d$column_name)

  id_ind_row <- d[d$column_name == "id_ind", ]
  expect_true(nchar(id_ind_row$description) > 0L)
})


test_that("describe_table() errors on unknown table", {
  pop <- open_pop(pop_name = "schema_test", db_name = ":memory:") |>
    define_genome(n_loci = 100, n_chr = 2, chr_len_Mb = 50)
  on.exit(close_pop(pop))

  expect_error(describe_table(pop, "does_not_exist"),
               regexp = "does not exist")
})


test_that("define_schema_description() upserts a table-level description", {
  pop <- open_pop(pop_name = "schema_test", db_name = ":memory:") |>
    define_genome(n_loci = 100, n_chr = 2, chr_len_Mb = 50)
  on.exit(close_pop(pop))

  pop |>
    get_table("ind_meta") |>
    define_schema_description(description = "Custom description for testing")
  row <- DBI::dbGetQuery(pop$db_conn,
    "SELECT description FROM _schema_meta
     WHERE object_type = 'table' AND table_name = 'ind_meta' AND column_name IS NULL")
  expect_equal(row$description, "Custom description for testing")
})


test_that("define_schema_description() upserts a column-level description", {
  pop <- open_pop(pop_name = "schema_test", db_name = ":memory:") |>
    define_genome(n_loci = 100, n_chr = 2, chr_len_Mb = 50)
  on.exit(close_pop(pop))

  pop |>
    get_table("ind_meta") |>
    define_schema_description("sex", "Sex of the individual (M/F)")

  d <- describe_table(pop, "ind_meta")
  sex_row <- d[d$column_name == "sex", ]
  expect_equal(sex_row$description, "Sex of the individual (M/F)")
})


test_that("define_schema_description() upsert is idempotent (replaces prior value)", {
  pop <- open_pop(pop_name = "schema_test", db_name = ":memory:") |>
    define_genome(n_loci = 100, n_chr = 2, chr_len_Mb = 50)
  on.exit(close_pop(pop))

  pop |> get_table("ind_meta") |> define_schema_description("sex", "First description")
  pop |> get_table("ind_meta") |> define_schema_description("sex", "Second description")

  d <- describe_table(pop, "ind_meta")
  sex_row <- d[d$column_name == "sex", ]
  expect_equal(nrow(sex_row), 1L)
  expect_equal(sex_row$description, "Second description")
})


test_that("define_schema_description() errors when column does not exist", {
  pop <- open_pop(pop_name = "schema_test", db_name = ":memory:") |>
    define_genome(n_loci = 100, n_chr = 2, chr_len_Mb = 50)
  on.exit(close_pop(pop))

  expect_error(
    pop |> get_table("ind_meta") |>
      define_schema_description("nonexistent_col", "Bad desc"),
    regexp = "not found"
  )
})


test_that("define_schema_description() can be chained for multiple columns", {
  pop <- open_pop(pop_name = "schema_test", db_name = ":memory:") |>
    define_genome(n_loci = 100, n_chr = 2, chr_len_Mb = 50)
  on.exit(close_pop(pop))

  pop |>
    get_table("ind_meta") |>
    define_schema_description("sex",       "Sex of individual")    |>
    define_schema_description("id_ind",    "Unique identifier")    |>
    define_schema_description("line_name", "Genetic line")

  d <- describe_table(pop, "ind_meta")
  expect_equal(d$description[d$column_name == "sex"],       "Sex of individual")
  expect_equal(d$description[d$column_name == "id_ind"],    "Unique identifier")
  expect_equal(d$description[d$column_name == "line_name"], "Genetic line")
})


test_that("print.tidybreed_pop() hides _schema_meta from table list", {
  pop <- open_pop(pop_name = "schema_test", db_name = ":memory:") |>
    define_genome(n_loci = 100, n_chr = 2, chr_len_Mb = 50)
  on.exit(close_pop(pop))

  output <- capture.output(print(pop))
  table_line <- output[grep("^Tables:", output)]
  expect_false(grepl("_schema_meta", table_line))
})


test_that("print.tidybreed_pop() shows schema() hint", {
  pop <- open_pop(pop_name = "schema_test", db_name = ":memory:") |>
    define_genome(n_loci = 100, n_chr = 2, chr_len_Mb = 50)
  on.exit(close_pop(pop))

  output <- paste(capture.output(print(pop)), collapse = "\n")
  expect_true(grepl("schema(pop)", output, fixed = TRUE))
})


test_that("define_chip() auto-registers column description in _schema_meta", {
  pop <- open_pop(pop_name = "schema_test", db_name = ":memory:") |>
    define_genome(n_loci = 100, n_chr = 2, chr_len_Mb = 50)
  on.exit(close_pop(pop))

  pop <- pop |>
    get_table("genome_meta") |>
    dplyr::filter(chr == 1L) |>
    define_chip("testchip")

  row <- DBI::dbGetQuery(pop$db_conn,
    "SELECT description FROM _schema_meta
     WHERE object_type = 'column' AND table_name = 'genome_meta'
     AND column_name = 'is_testchip'")
  expect_equal(nrow(row), 1L)
  expect_true(grepl("testchip", row$description))
})


test_that("summary.tidybreed_pop() pulls descriptions from _schema_meta (not hard-coded)", {
  pop <- open_pop(pop_name = "schema_test", db_name = ":memory:") |>
    define_genome(n_loci = 100, n_chr = 2, chr_len_Mb = 50)
  on.exit(close_pop(pop))

  # Override the genome_meta table description so we can verify it was pulled from DB
  pop |>
    get_table("genome_meta") |>
    define_schema_description(description = "CUSTOM_TEST_DESCRIPTION_XYZ")

  summ <- summary(pop)
  genome_desc <- summ$tables[["genome_meta"]]$description
  expect_equal(genome_desc, "CUSTOM_TEST_DESCRIPTION_XYZ")
})


test_that("summary.tidybreed_pop() excludes _schema_meta from results", {
  pop <- open_pop(pop_name = "schema_test", db_name = ":memory:") |>
    define_genome(n_loci = 100, n_chr = 2, chr_len_Mb = 50)
  on.exit(close_pop(pop))

  summ <- summary(pop)
  expect_false("_schema_meta" %in% names(summ$tables))
})
