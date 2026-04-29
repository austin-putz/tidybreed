# Tests for .set_default parameter in mutate_table()

# Helper to create test population
create_test_pop <- function(n_haplotypes = 50) {
  initialize_genome(
    pop_name = "test_defaults",
    n_loci = 100,
    n_chr = 5,
    chr_len_Mb = 50,
    n_haplotypes = n_haplotypes,
    db_path = ":memory:"
  )
}

# ── Basic functionality ──────────────────────────────────────────────────────

test_that(".set_default creates DEFAULT on empty table", {
  pop <- create_test_pop()

  # Add column with DEFAULT
  pop <- pop |>
    get_table("ind_meta") |>
    mutate_table(gen = 0L, .set_default = TRUE)

  # Verify DEFAULT in schema
  schema_info <- DBI::dbGetQuery(
    pop$db_conn,
    "SELECT column_name, column_default FROM information_schema.columns
     WHERE table_name = 'ind_meta' AND column_name = 'gen'"
  )

  expect_equal(nrow(schema_info), 1)
  expect_false(is.na(schema_info$column_default))

  close_pop(pop)
})

test_that("future INSERTs use DEFAULT value", {
  pop <- create_test_pop()

  # Set defaults
  pop <- pop |>
    get_table("ind_meta") |>
    mutate_table(gen = 0L, active = TRUE, .set_default = TRUE)

  # Add founders WITHOUT specifying gen or active
  pop <- pop |>
    add_founders(n_males = 5, n_females = 5, line_name = "A")

  # Verify defaults applied
  result <- pop |> get_table("ind_meta") |> collect()
  expect_true(all(result$gen == 0L))
  expect_true(all(result$active == TRUE))

  close_pop(pop)
})

test_that("explicit values override DEFAULT", {
  pop <- create_test_pop()

  # Set default
  pop <- pop |>
    get_table("ind_meta") |>
    mutate_table(gen = 0L, .set_default = TRUE)

  # Explicit override
  pop <- pop |>
    add_founders(n_males = 5, n_females = 5, line_name = "A", gen = 1L)

  result <- pop |> get_table("ind_meta") |> collect()
  expect_true(all(result$gen == 1L))  # Override worked

  close_pop(pop)
})

test_that(".set_default works on populated table", {
  pop <- create_test_pop()

  # Add founders first
  pop <- pop |>
    add_founders(n_males = 5, n_females = 5, line_name = "A")

  # Add column with DEFAULT to populated table
  pop <- pop |>
    get_table("ind_meta") |>
    mutate_table(farm = "Iowa", .set_default = TRUE)

  # Existing rows get UPDATE
  result1 <- pop |> get_table("ind_meta") |> collect()
  expect_true(all(result1$farm == "Iowa"))

  # New rows get DEFAULT
  pop <- pop |>
    add_founders(n_males = 2, n_females = 2, line_name = "B")

  result2 <- pop |>
    get_table("ind_meta") |>
    dplyr::filter(line == "B") |>
    collect()
  expect_true(all(result2$farm == "Iowa"))

  close_pop(pop)
})

# ── Type coverage ────────────────────────────────────────────────────────────

test_that("BOOLEAN type works with DEFAULT", {
  pop <- create_test_pop()

  pop <- pop |>
    get_table("ind_meta") |>
    mutate_table(active = TRUE, .set_default = TRUE)

  pop <- pop |> add_founders(n_males = 5, n_females = 5, line_name = "A")

  result <- pop |> get_table("ind_meta") |> collect()
  expect_true(all(result$active == TRUE))

  close_pop(pop)
})

test_that("INTEGER type works with DEFAULT", {
  pop <- create_test_pop()

  pop <- pop |>
    get_table("ind_meta") |>
    mutate_table(gen = 42L, .set_default = TRUE)

  pop <- pop |> add_founders(n_males = 5, n_females = 5, line_name = "A")

  result <- pop |> get_table("ind_meta") |> collect()
  expect_true(all(result$gen == 42L))

  close_pop(pop)
})

test_that("DOUBLE type works with DEFAULT", {
  pop <- create_test_pop()

  pop <- pop |>
    get_table("ind_meta") |>
    mutate_table(weight = 3.14, .set_default = TRUE)

  pop <- pop |> add_founders(n_males = 5, n_females = 5, line_name = "A")

  result <- pop |> get_table("ind_meta") |> collect()
  expect_true(all(abs(result$weight - 3.14) < 1e-10))

  close_pop(pop)
})

test_that("VARCHAR type works with DEFAULT", {
  pop <- create_test_pop()

  pop <- pop |>
    get_table("ind_meta") |>
    mutate_table(farm = "TestFarm", .set_default = TRUE)

  pop <- pop |> add_founders(n_males = 5, n_females = 5, line_name = "A")

  result <- pop |> get_table("ind_meta") |> collect()
  expect_true(all(result$farm == "TestFarm"))

  close_pop(pop)
})

test_that("DATE type works with DEFAULT", {
  pop <- create_test_pop()

  test_date <- as.Date("2024-01-01")
  pop <- pop |>
    get_table("ind_meta") |>
    mutate_table(birth_date = test_date, .set_default = TRUE)

  pop <- pop |> add_founders(n_males = 5, n_females = 5, line_name = "A")

  result <- pop |> get_table("ind_meta") |> collect()
  expect_true(all(result$birth_date == test_date))

  close_pop(pop)
})

test_that("NA values produce DEFAULT NULL", {
  pop <- create_test_pop()

  pop <- pop |>
    get_table("ind_meta") |>
    mutate_table(gen = NA_integer_, .set_default = TRUE)

  # Schema should have DEFAULT NULL
  schema_info <- DBI::dbGetQuery(
    pop$db_conn,
    "SELECT column_default FROM information_schema.columns
     WHERE table_name = 'ind_meta' AND column_name = 'gen'"
  )

  # DuckDB may return NA or "NULL" string
  expect_true(is.na(schema_info$column_default) ||
              toupper(trimws(schema_info$column_default)) == "NULL")

  close_pop(pop)
})

# ── Error cases ──────────────────────────────────────────────────────────────

test_that(".set_default with vector errors", {
  pop <- create_test_pop()

  # Add some rows
  pop <- pop |> add_founders(n_males = 5, n_females = 5, line_name = "A")

  expect_error(
    pop |>
      get_table("ind_meta") |>
      mutate_table(gen = 1:10, .set_default = TRUE),
    "Cannot use .set_default = TRUE with vector values"
  )

  close_pop(pop)
})

test_that("invalid .set_default value errors", {
  pop <- create_test_pop()

  expect_error(
    pop |>
      get_table("ind_meta") |>
      mutate_table(gen = 0L, .set_default = NA),
    "is.logical\\(.set_default\\) is not TRUE|!is.na\\(.set_default\\) is not TRUE"
  )

  expect_error(
    pop |>
      get_table("ind_meta") |>
      mutate_table(gen = 0L, .set_default = "yes"),
    "is.logical\\(.set_default\\) is not TRUE"
  )

  close_pop(pop)
})

# ── Edge cases ───────────────────────────────────────────────────────────────

test_that("multiple columns with .set_default", {
  pop <- create_test_pop()

  pop <- pop |>
    get_table("ind_meta") |>
    mutate_table(
      gen = 0L,
      active = TRUE,
      farm = "Iowa",
      .set_default = TRUE
    )

  # Verify all have DEFAULT
  schema_info <- DBI::dbGetQuery(
    pop$db_conn,
    "SELECT column_name, column_default FROM information_schema.columns
     WHERE table_name = 'ind_meta' AND column_name IN ('gen', 'active', 'farm')"
  )

  expect_equal(nrow(schema_info), 3)
  expect_true(all(!is.na(schema_info$column_default)))

  # Verify INSERTs use defaults
  pop <- pop |> add_founders(n_males = 5, n_females = 5, line_name = "A")
  result <- pop |> get_table("ind_meta") |> collect()
  expect_true(all(result$gen == 0L))
  expect_true(all(result$active == TRUE))
  expect_true(all(result$farm == "Iowa"))

  close_pop(pop)
})

test_that(".set_default on existing column is ignored", {
  pop <- create_test_pop()

  # Create column without DEFAULT
  pop <- pop |>
    get_table("ind_meta") |>
    mutate_table(gen = 0L)

  # Add founders
  pop <- pop |> add_founders(n_males = 5, n_females = 5, line_name = "A")

  # Try to add DEFAULT to existing column (should be ignored, UPDATE still works)
  pop <- pop |>
    get_table("ind_meta") |>
    mutate_table(gen = 1L, .set_default = TRUE)

  # UPDATE happened
  result <- pop |> get_table("ind_meta") |> collect()
  expect_true(all(result$gen == 1L))

  close_pop(pop)
})

test_that("partial override of defaults", {
  pop <- create_test_pop()

  # Set multiple defaults
  pop <- pop |>
    get_table("ind_meta") |>
    mutate_table(
      gen = 0L,
      active = TRUE,
      farm = "FarmA",
      .set_default = TRUE
    )

  # Override only gen
  pop <- pop |>
    add_founders(n_males = 5, n_females = 5, line_name = "A", gen = 1L)

  result <- pop |> get_table("ind_meta") |> collect()
  expect_true(all(result$gen == 1L))       # Overridden
  expect_true(all(result$active == TRUE))   # Default
  expect_true(all(result$farm == "FarmA"))  # Default

  close_pop(pop)
})

# ── Integration tests ────────────────────────────────────────────────────────

test_that("DEFAULT works with add_phenotype()", {
  pop <- create_test_pop()

  # Add trait and founders
  pop <- pop |>
    add_trait_simple("ADG", n_qtl = 10, target_add_var = 100, residual_var = 50)

  pop <- pop |>
    add_founders(n_males = 5, n_females = 5, line_name = "A")

  # Set default for phenotype location
  pop <- pop |>
    get_table("ind_phenotype") |>
    mutate_table(location = "Lab1", .set_default = TRUE)

  # Add phenotypes without location
  pop <- pop |>
    get_table("ind_meta") |>
    add_phenotype("ADG")

  # Verify default applied
  phenos <- pop |> get_table("ind_phenotype") |> collect()
  expect_true(all(phenos$location == "Lab1"))

  close_pop(pop)
})

test_that("DEFAULT persists across multiple add_founders() calls", {
  pop <- create_test_pop()

  # Set defaults
  pop <- pop |>
    get_table("ind_meta") |>
    mutate_table(gen = 0L, farm = "MainFarm", .set_default = TRUE)

  # First batch
  pop <- pop |>
    add_founders(n_males = 5, n_females = 5, line_name = "A")

  # Second batch (different line, defaults still apply)
  pop <- pop |>
    add_founders(n_males = 3, n_females = 3, line_name = "B")

  result <- pop |> get_table("ind_meta") |> collect()
  expect_equal(nrow(result), 16)
  expect_true(all(result$gen == 0L))
  expect_true(all(result$farm == "MainFarm"))

  close_pop(pop)
})

test_that(".set_default = FALSE is default behavior", {
  pop <- create_test_pop()

  # Don't specify .set_default (should default to FALSE)
  pop <- pop |>
    get_table("ind_meta") |>
    mutate_table(gen = 0L)

  # Check schema - should have no DEFAULT
  schema_info <- DBI::dbGetQuery(
    pop$db_conn,
    "SELECT column_default FROM information_schema.columns
     WHERE table_name = 'ind_meta' AND column_name = 'gen'"
  )

  # Should be NULL (no default)
  expect_true(is.na(schema_info$column_default))

  close_pop(pop)
})
