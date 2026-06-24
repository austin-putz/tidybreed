test_that("archive_replicate() errors on non-tidybreed_pop input", {
  expect_error(archive_replicate("not a pop"), "'pop' must be a tidybreed_pop")
})

test_that("archive_replicate() errors on invalid replicate values", {
  pop <- make_test_pop()
  arc <- tempfile(fileext = ".duckdb")
  on.exit(unlink(arc), add = TRUE)

  expect_error(archive_replicate(pop, replicate = 0L,  archive_path = arc), "positive integer")
  expect_error(archive_replicate(pop, replicate = NA,  archive_path = arc), "positive integer")
  expect_error(archive_replicate(pop, replicate = -1L, archive_path = arc), "positive integer")
})

test_that("archive_replicate() errors when a table appears in two categories", {
  pop <- make_test_pop()
  arc <- tempfile(fileext = ".duckdb")
  on.exit(unlink(arc), add = TRUE)

  expect_error(
    archive_replicate(pop,
                      replicate = 1L,
                      archive_path = arc,
                      store_and_reset = c("ind_meta", "ind_phenotype"),
                      store_once      = c("ind_meta", "genome_meta")),
    "more than one category"
  )
})

test_that("archive_replicate() creates store_and_reset tables with replicate column", {
  pop <- make_test_pop()
  arc <- tempfile(fileext = ".duckdb")
  on.exit({
    DBI::dbDisconnect(pop$db_conn, shutdown = TRUE)
    unlink(arc)
  }, add = TRUE)

  pop <- archive_replicate(pop,
                           replicate       = 1L,
                           archive_path    = arc,
                           store_and_reset = "ind_meta",
                           store_once      = character(0L),
                           reset_only      = character(0L))

  # Archive DB should have ind_meta with replicate column
  arc_conn <- DBI::dbConnect(duckdb::duckdb(), dbdir = arc, read_only = TRUE)
  on.exit(DBI::dbDisconnect(arc_conn, shutdown = TRUE), add = TRUE)

  expect_true("ind_meta" %in% DBI::dbListTables(arc_conn))
  arc_cols <- DBI::dbListFields(arc_conn, "ind_meta")
  expect_true("replicate" %in% arc_cols)
  expect_true("id_ind"    %in% arc_cols)

  rows <- DBI::dbGetQuery(arc_conn, "SELECT DISTINCT replicate FROM ind_meta")
  expect_equal(rows$replicate, 1L)
})

test_that("archive_replicate() creates store_once tables WITHOUT replicate column", {
  pop <- make_test_pop()
  pop <- define_trait(pop, "ADG", target_add_var = 100)

  arc <- tempfile(fileext = ".duckdb")
  on.exit({
    DBI::dbDisconnect(pop$db_conn, shutdown = TRUE)
    unlink(arc)
  }, add = TRUE)

  pop <- archive_replicate(pop,
                           replicate       = 1L,
                           archive_path    = arc,
                           store_and_reset = character(0L),
                           store_once      = "trait_meta",
                           reset_only      = character(0L))

  arc_conn <- DBI::dbConnect(duckdb::duckdb(), dbdir = arc, read_only = TRUE)
  on.exit(DBI::dbDisconnect(arc_conn, shutdown = TRUE), add = TRUE)

  expect_true("trait_meta" %in% DBI::dbListTables(arc_conn))
  arc_cols <- DBI::dbListFields(arc_conn, "trait_meta")
  expect_false("replicate" %in% arc_cols)
})

test_that("second call appends new replicate rows without overwriting prior ones", {
  pop <- make_test_pop()
  arc <- tempfile(fileext = ".duckdb")
  on.exit({
    DBI::dbDisconnect(pop$db_conn, shutdown = TRUE)
    unlink(arc)
  }, add = TRUE)

  # First replicate
  pop <- archive_replicate(pop,
                           replicate       = 1L,
                           archive_path    = arc,
                           store_and_reset = "ind_meta",
                           store_once      = character(0L),
                           reset_only      = character(0L))

  # Add new founders and archive a second replicate
  pop <- get_table(pop, "founder_haplotypes") |>
    add_founders(n_males = 3, n_females = 3, line_name = "B")

  pop <- archive_replicate(pop,
                           replicate       = 2L,
                           archive_path    = arc,
                           store_and_reset = "ind_meta",
                           store_once      = character(0L),
                           reset_only      = character(0L))

  arc_conn <- DBI::dbConnect(duckdb::duckdb(), dbdir = arc, read_only = TRUE)
  on.exit(DBI::dbDisconnect(arc_conn, shutdown = TRUE), add = TRUE)

  reps <- DBI::dbGetQuery(arc_conn, "SELECT DISTINCT replicate FROM ind_meta ORDER BY replicate")$replicate
  expect_equal(reps, c(1L, 2L))
})

test_that("store_once tables are NOT duplicated on second call", {
  pop <- make_test_pop()
  pop <- define_trait(pop, "ADG", target_add_var = 100)

  arc <- tempfile(fileext = ".duckdb")
  on.exit({
    DBI::dbDisconnect(pop$db_conn, shutdown = TRUE)
    unlink(arc)
  }, add = TRUE)

  pop <- archive_replicate(pop,
                           replicate       = 1L,
                           archive_path    = arc,
                           store_and_reset = "ind_meta",
                           store_once      = "trait_meta",
                           reset_only      = character(0L))

  pop <- get_table(pop, "founder_haplotypes") |>
    add_founders(n_males = 3, n_females = 3, line_name = "B")

  pop <- archive_replicate(pop,
                           replicate       = 2L,
                           archive_path    = arc,
                           store_and_reset = "ind_meta",
                           store_once      = "trait_meta",
                           reset_only      = character(0L))

  arc_conn <- DBI::dbConnect(duckdb::duckdb(), dbdir = arc, read_only = TRUE)
  on.exit(DBI::dbDisconnect(arc_conn, shutdown = TRUE), add = TRUE)

  # trait_meta should have only the original rows — no duplication from rep 2
  n_traits <- DBI::dbGetQuery(arc_conn, "SELECT COUNT(*) AS n FROM trait_meta")$n
  expect_equal(n_traits, 1L)  # only "ADG"
})

test_that("store_and_reset tables are empty in working DB after success", {
  pop <- make_test_pop()
  arc <- tempfile(fileext = ".duckdb")
  on.exit({
    DBI::dbDisconnect(pop$db_conn, shutdown = TRUE)
    unlink(arc)
  }, add = TRUE)

  n_before <- DBI::dbGetQuery(pop$db_conn, "SELECT COUNT(*) AS n FROM ind_meta")$n
  expect_gt(n_before, 0L)

  pop <- archive_replicate(pop,
                           replicate       = 1L,
                           archive_path    = arc,
                           store_and_reset = "ind_meta",
                           store_once      = character(0L),
                           reset_only      = character(0L))

  n_after <- DBI::dbGetQuery(pop$db_conn, "SELECT COUNT(*) AS n FROM ind_meta")$n
  expect_equal(n_after, 0L)
})

test_that("reset_only tables are empty in working DB after success", {
  pop <- make_test_pop()
  arc <- tempfile(fileext = ".duckdb")
  on.exit({
    DBI::dbDisconnect(pop$db_conn, shutdown = TRUE)
    unlink(arc)
  }, add = TRUE)

  n_hap <- DBI::dbGetQuery(pop$db_conn, "SELECT COUNT(*) AS n FROM genome_haplotype")$n
  expect_gt(n_hap, 0L)

  pop <- archive_replicate(pop,
                           replicate       = 1L,
                           archive_path    = arc,
                           store_and_reset = character(0L),
                           store_once      = character(0L),
                           reset_only      = "genome_haplotype")

  n_after <- DBI::dbGetQuery(pop$db_conn, "SELECT COUNT(*) AS n FROM genome_haplotype")$n
  expect_equal(n_after, 0L)
})

test_that("tidybreed.replicate increments after success", {
  pop <- make_test_pop()
  arc <- tempfile(fileext = ".duckdb")
  on.exit({
    DBI::dbDisconnect(pop$db_conn, shutdown = TRUE)
    unlink(arc)
  }, add = TRUE)

  old_opt <- getOption("tidybreed.replicate")
  on.exit(options(tidybreed.replicate = old_opt), add = TRUE)

  options(tidybreed.replicate = 5L)
  archive_replicate(pop,
                    replicate       = 5L,
                    archive_path    = arc,
                    store_and_reset = character(0L),
                    store_once      = character(0L),
                    reset_only      = character(0L))

  expect_equal(getOption("tidybreed.replicate"), 6L)
})

test_that("missing tables skipped with message, not error", {
  pop <- make_test_pop()
  arc <- tempfile(fileext = ".duckdb")
  on.exit({
    DBI::dbDisconnect(pop$db_conn, shutdown = TRUE)
    unlink(arc)
  }, add = TRUE)

  # A table that will never exist in any pop; should skip with message, not error
  expect_message(
    archive_replicate(pop,
                      replicate       = 1L,
                      archive_path    = arc,
                      store_and_reset = c("ind_meta", "this_table_does_not_exist"),
                      store_once      = character(0L),
                      reset_only      = character(0L)),
    "skipping"
  )
})

test_that("db_name_archive = NULL skips archiving but still resets", {
  old_opts <- options(tidybreed.db_name_archive = NULL)
  on.exit(options(old_opts), add = TRUE)

  pop <- make_test_pop()
  on.exit(DBI::dbDisconnect(pop$db_conn, shutdown = TRUE), add = TRUE)

  n_before <- DBI::dbGetQuery(pop$db_conn, "SELECT COUNT(*) AS n FROM ind_meta")$n
  expect_gt(n_before, 0L)

  pop <- archive_replicate(pop,
                           replicate       = 1L,
                           archive_path    = NULL,
                           store_and_reset = "ind_meta",
                           store_once      = character(0L),
                           reset_only      = character(0L))

  n_after <- DBI::dbGetQuery(pop$db_conn, "SELECT COUNT(*) AS n FROM ind_meta")$n
  expect_equal(n_after, 0L)
})

test_that("replicate column collision guard errors before any write", {
  pop <- make_test_pop()
  # Add a 'replicate' column to ind_meta to trigger the guard
  DBI::dbExecute(pop$db_conn,
    "ALTER TABLE ind_meta ADD COLUMN replicate INTEGER DEFAULT 1")

  arc <- tempfile(fileext = ".duckdb")
  on.exit({
    DBI::dbDisconnect(pop$db_conn, shutdown = TRUE)
    unlink(arc)
  }, add = TRUE)

  expect_error(
    archive_replicate(pop,
                      replicate       = 1L,
                      archive_path    = arc,
                      store_and_reset = "ind_meta",
                      store_once      = character(0L),
                      reset_only      = character(0L)),
    "already contains a column named 'replicate'"
  )

  # Working DB rows should be untouched since we errored before writing
  n_rows <- DBI::dbGetQuery(pop$db_conn, "SELECT COUNT(*) AS n FROM ind_meta")$n
  expect_gt(n_rows, 0L)
})

test_that("_archive_meta provenance table is written and accumulates rows", {
  pop <- make_test_pop()
  arc <- tempfile(fileext = ".duckdb")
  on.exit({
    DBI::dbDisconnect(pop$db_conn, shutdown = TRUE)
    unlink(arc)
  }, add = TRUE)

  pop <- archive_replicate(pop,
                           replicate       = 1L,
                           archive_path    = arc,
                           store_and_reset = "ind_meta",
                           store_once      = character(0L),
                           reset_only      = character(0L))

  pop <- get_table(pop, "founder_haplotypes") |>
    add_founders(n_males = 2, n_females = 2, line_name = "C")

  pop <- archive_replicate(pop,
                           replicate       = 2L,
                           archive_path    = arc,
                           store_and_reset = "ind_meta",
                           store_once      = character(0L),
                           reset_only      = character(0L))

  arc_conn <- DBI::dbConnect(duckdb::duckdb(), dbdir = arc, read_only = TRUE)
  on.exit(DBI::dbDisconnect(arc_conn, shutdown = TRUE), add = TRUE)

  expect_true("_archive_meta" %in% DBI::dbListTables(arc_conn))
  meta <- DBI::dbGetQuery(arc_conn, 'SELECT * FROM "_archive_meta" ORDER BY replicate')
  expect_equal(nrow(meta), 2L)
  expect_equal(meta$replicate, c(1L, 2L))
  expect_true(all(!is.na(meta$tidybreed_version)))
})

test_that("archive_replicate() returns pop invisibly", {
  pop <- make_test_pop()
  arc <- tempfile(fileext = ".duckdb")
  on.exit({
    DBI::dbDisconnect(pop$db_conn, shutdown = TRUE)
    unlink(arc)
  }, add = TRUE)

  result <- archive_replicate(pop,
                              replicate       = 1L,
                              archive_path    = arc,
                              store_and_reset = character(0L),
                              store_once      = character(0L),
                              reset_only      = character(0L))

  expect_true(inherits(result, "tidybreed_pop"))
})
