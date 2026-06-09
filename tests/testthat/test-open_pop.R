test_that("open_pop() returns tidybreed_pop with correct class and pop_name", {
  pop <- open_pop(pop_name = "my_sim", db_name = ":memory:")
  on.exit(close_pop(pop))

  expect_s3_class(pop, "tidybreed_pop")
  expect_equal(pop$pop_name, "my_sim")
  expect_equal(pop$db_path, ":memory:")
})

test_that("open_pop() creates all core tables before define_genome()", {
  pop <- open_pop(db_name = ":memory:")
  on.exit(close_pop(pop))

  core_tables <- c("_schema_meta", "ind_meta", "trait_var_comp",
                   "genome_effects", "phenotype_meta", "phenotype_components",
                   "phenotype_var_comp")
  expect_true(all(core_tables %in% pop$tables))

  # Genome tables are NOT present yet
  expect_false("genome_meta" %in% pop$tables)
  expect_false("genome_haplotype" %in% pop$tables)
})

test_that("open_pop() with :memory: sets run_dirs to character(0)", {
  pop <- open_pop(db_name = ":memory:")
  on.exit(close_pop(pop))

  expect_equal(pop$run_dirs, character(0))
  expect_equal(length(pop$run_dirs), 0)
})

test_that("open_pop() on disk creates layer 2-3 folders and run_dirs base", {
  tmp <- tempfile()
  dir.create(tmp)
  on.exit(unlink(tmp, recursive = TRUE))

  pop <- open_pop(pop_name   = "t",
                  base_dir   = tmp,
                  output_dir = "out",
                  scenario_dir = "sc1",
                  tools      = NULL,
                  db_name    = "sim.duckdb")
  on.exit(close_pop(pop), add = TRUE)

  expect_true(dir.exists(file.path(tmp, "out", "sc1")))
  expect_true(file.exists(file.path(tmp, "out", "sc1", "sim.duckdb")))
  expect_equal(pop$run_dirs[["base"]], file.path(tmp, "out", "sc1"))
})

test_that("open_pop() creates tool subdirs when tools is provided", {
  tmp <- tempfile()
  dir.create(tmp)
  on.exit(unlink(tmp, recursive = TRUE))

  pop <- open_pop(pop_name     = "t",
                  base_dir     = tmp,
                  output_dir   = "out",
                  scenario_dir = "sc1",
                  tools        = c("blupf90", "plink"),
                  db_name      = "sim.duckdb")
  on.exit(close_pop(pop), add = TRUE)

  expect_true(dir.exists(file.path(tmp, "out", "sc1", "blupf90")))
  expect_true(dir.exists(file.path(tmp, "out", "sc1", "plink")))
  expect_equal(names(pop$run_dirs), c("base", "blupf90", "plink"))
  expect_equal(pop$run_dirs[["blupf90"]],
               file.path(tmp, "out", "sc1", "blupf90"))
})

test_that("open_pop() auto-generates scenario_dir when scenario_dir is NULL", {
  tmp <- tempfile()
  dir.create(tmp)
  on.exit(unlink(tmp, recursive = TRUE))

  pop <- open_pop(pop_name     = "t",
                  base_dir     = tmp,
                  output_dir   = "out",
                  scenario_dir = NULL,
                  tools        = NULL,
                  db_name      = "sim.duckdb")
  on.exit(close_pop(pop), add = TRUE)

  # Scenario dir should be a timestamp folder matching YYYYMMDD_HHMMSS
  scenario_dirs <- list.dirs(file.path(tmp, "out"),
                             recursive = FALSE, full.names = FALSE)
  expect_length(scenario_dirs, 1)
  expect_match(scenario_dirs[1], "^[0-9]{8}_[0-9]{6}$")
})

test_that("open_pop() clean=TRUE deletes existing db before recreating", {
  tmp <- tempfile()
  dir.create(tmp)
  on.exit(unlink(tmp, recursive = TRUE))

  # First run
  pop1 <- open_pop(pop_name     = "t",
                   base_dir     = tmp,
                   output_dir   = "out",
                   scenario_dir = "sc1",
                   tools        = NULL,
                   db_name      = "sim.duckdb",
                   clean        = FALSE)

  # Write a dummy table to first run's DB so we can check it was replaced
  DBI::dbExecute(pop1$db_conn, "CREATE TABLE _sentinel (x INTEGER)")
  close_pop(pop1)

  # Second run with clean = TRUE — should delete and recreate DB
  pop2 <- open_pop(pop_name     = "t",
                   base_dir     = tmp,
                   output_dir   = "out",
                   scenario_dir = "sc1",
                   tools        = NULL,
                   db_name      = "sim.duckdb",
                   clean        = TRUE)
  on.exit(close_pop(pop2), add = TRUE)

  existing <- DBI::dbListTables(pop2$db_conn)
  expect_false("_sentinel" %in% existing)
  expect_true("ind_meta" %in% existing)
})

test_that(".create_run_dir() creates a timestamped subdirectory", {
  tmp <- tempfile()
  dir.create(tmp)
  on.exit(unlink(tmp, recursive = TRUE))

  pop <- open_pop(pop_name     = "t",
                  base_dir     = tmp,
                  output_dir   = "out",
                  scenario_dir = "sc1",
                  tools        = "blupf90",
                  db_name      = "sim.duckdb")
  on.exit(close_pop(pop), add = TRUE)

  run_dir <- tidybreed:::.create_run_dir(pop, "blupf90")
  expect_true(dir.exists(run_dir))
  expect_true(grepl("^[0-9]{8}_[0-9]{6}_[a-z0-9]{6}$", basename(run_dir)))
})

test_that(".create_run_dir() errors for undeclared tool", {
  tmp <- tempfile()
  dir.create(tmp)
  on.exit(unlink(tmp, recursive = TRUE))

  pop <- open_pop(pop_name     = "t",
                  base_dir     = tmp,
                  output_dir   = "out",
                  scenario_dir = "sc1",
                  tools        = "blupf90",
                  db_name      = "sim.duckdb")
  on.exit(close_pop(pop), add = TRUE)

  expect_error(
    tidybreed:::.create_run_dir(pop, "plink"),
    regexp = "not registered"
  )
})

test_that(".create_run_dir() errors when run_dirs is empty (in-memory pop)", {
  pop <- open_pop(db_name = ":memory:")
  on.exit(close_pop(pop))

  expect_error(
    tidybreed:::.create_run_dir(pop, "blupf90"),
    regexp = "not registered"
  )
})

test_that("ensure_trait_tables() tables are present after open_pop()", {
  pop <- open_pop(db_name = ":memory:")
  on.exit(close_pop(pop))

  trait_tables <- c("trait_meta", "trait_effects", "ind_phenotype",
                    "ind_tbv", "ind_ebv", "index_meta", "ind_index")
  expect_true(all(trait_tables %in% pop$tables))
})
