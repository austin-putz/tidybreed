# Helper: minimal pop (no founders needed for table tests)
make_dt_pop <- function() {
  initialize_genome(
    pop_name   = "dt_test",
    n_loci     = 50,
    n_chr      = 1,
    chr_len_Mb = 100,
    db_path    = ":memory:"
  )
}


# ============================================================
# Basic creation
# ============================================================

test_that("define_table() creates table with correct schema", {
  pop <- make_dt_pop()
  on.exit(close_pop(pop))

  pop <- pop |> define_table(
    "sim_timing",
    run_id       = NA_integer_,
    duration_sec = NA_real_,
    label        = NA_character_
  )

  # Table appears in pop$tables
  expect_true("sim_timing" %in% pop$tables)

  # Table exists in database
  expect_true(DBI::dbExistsTable(pop$db_conn, "sim_timing"))

  # Schema matches
  cols <- DBI::dbListFields(pop$db_conn, "sim_timing")
  expect_equal(cols, c("run_id", "duration_sec", "label"))
})


test_that("define_table() supports BOOLEAN and DATE columns", {
  pop <- make_dt_pop()
  on.exit(close_pop(pop))

  pop <- pop |> define_table(
    "misc_types",
    flag   = FALSE,
    dob    = as.Date(NA),
    score  = NA_real_
  )

  info <- DBI::dbGetQuery(pop$db_conn,
    "PRAGMA table_info('misc_types')")
  types <- setNames(info$type, info$name)

  expect_equal(types[["flag"]],  "BOOLEAN")
  expect_equal(types[["dob"]],   "DATE")
  expect_equal(types[["score"]], "DOUBLE")
})


# ============================================================
# primary_key
# ============================================================

test_that("define_table() declares PRIMARY KEY correctly", {
  pop <- make_dt_pop()
  on.exit(close_pop(pop))

  pop <- pop |> define_table(
    "runs",
    run_id = NA_integer_,
    notes  = NA_character_,
    primary_key = "run_id"
  )

  info <- DBI::dbGetQuery(pop$db_conn, "PRAGMA table_info('runs')")
  pk_col <- info$name[info$pk > 0]
  expect_equal(pk_col, "run_id")
})


test_that("define_table() with no primary_key creates no PK constraint", {
  pop <- make_dt_pop()
  on.exit(close_pop(pop))

  pop <- pop |> define_table("nokey", x = NA_integer_, y = NA_real_)

  info <- DBI::dbGetQuery(pop$db_conn, "PRAGMA table_info('nokey')")
  expect_true(all(info$pk == 0))
})


# ============================================================
# overwrite
# ============================================================

test_that("define_table() errors when table exists and overwrite = FALSE", {
  pop <- make_dt_pop()
  on.exit(close_pop(pop))

  pop <- pop |> define_table("dup_tbl", x = NA_integer_)

  expect_error(
    pop |> define_table("dup_tbl", x = NA_integer_),
    "already exists"
  )
})


test_that("define_table() overwrites when overwrite = TRUE", {
  pop <- make_dt_pop()
  on.exit(close_pop(pop))

  pop <- pop |> define_table("dup_tbl", x = NA_integer_)

  # Insert a row so we can verify it's gone after overwrite
  DBI::dbAppendTable(pop$db_conn, "dup_tbl",
                     data.frame(x = 99L))

  pop <- pop |> define_table("dup_tbl", x = NA_integer_, y = NA_real_,
                              overwrite = TRUE)

  # Table recreated with new schema
  cols <- DBI::dbListFields(pop$db_conn, "dup_tbl")
  expect_equal(cols, c("x", "y"))

  # Old data gone
  n <- DBI::dbGetQuery(pop$db_conn, "SELECT COUNT(*) AS n FROM dup_tbl")$n
  expect_equal(n, 0L)

  # Still registered
  expect_true("dup_tbl" %in% pop$tables)
})


# ============================================================
# Input validation errors
# ============================================================

test_that("define_table() errors on empty ... (no columns)", {
  pop <- make_dt_pop()
  on.exit(close_pop(pop))

  expect_error(pop |> define_table("empty_cols"), "at least one column")
})


test_that("define_table() errors on reserved system table name", {
  pop <- make_dt_pop()
  on.exit(close_pop(pop))

  expect_error(
    pop |> define_table("ind_meta", x = NA_integer_),
    "system-managed"
  )
  expect_error(
    pop |> define_table("genome_meta", x = NA_integer_),
    "system-managed"
  )
})


test_that("define_table() errors on invalid SQL identifier for table name", {
  pop <- make_dt_pop()
  on.exit(close_pop(pop))

  expect_error(pop |> define_table("123bad",  x = NA_integer_), "Invalid")
  expect_error(pop |> define_table("bad-name", x = NA_integer_), "Invalid")
})


test_that("define_table() errors on invalid SQL identifier for column name", {
  pop <- make_dt_pop()
  on.exit(close_pop(pop))

  expect_error(
    pop |> define_table("my_tbl", `bad col` = NA_integer_),
    "Invalid"
  )
})


test_that("define_table() errors when primary_key not in defined columns", {
  pop <- make_dt_pop()
  on.exit(close_pop(pop))

  expect_error(
    pop |> define_table("runs", run_id = NA_integer_, primary_key = "nonexistent"),
    "not among the defined columns"
  )
})


# ============================================================
# Integration: get_table() + DBI round-trip
# ============================================================

test_that("get_table() works immediately after define_table()", {
  pop <- make_dt_pop()
  on.exit(close_pop(pop))

  pop <- pop |> define_table(
    "sim_timing",
    run_id       = NA_integer_,
    duration_sec = NA_real_,
    label        = NA_character_,
    primary_key  = "run_id"
  )

  # Append a row
  DBI::dbAppendTable(
    pop$db_conn, "sim_timing",
    tibble::tibble(run_id = 1L, duration_sec = 3.14, label = "test")
  )

  # Read back via get_table()
  result <- get_table(pop, "sim_timing") |> dplyr::collect()

  expect_equal(nrow(result), 1L)
  expect_equal(result$run_id,       1L)
  expect_equal(result$duration_sec, 3.14, tolerance = 1e-9)
  expect_equal(result$label,        "test")
})


test_that("mutate_table() adds columns to a user-defined table", {
  pop <- make_dt_pop()
  on.exit(close_pop(pop))

  pop <- pop |> define_table("custom", x = NA_integer_)
  pop <- pop |> get_table("custom") |> mutate_table(y = NA_character_)

  cols <- DBI::dbListFields(pop$db_conn, "custom")
  expect_true("y" %in% cols)
})
