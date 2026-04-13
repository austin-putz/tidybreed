# Note: These tests will work once add_founders() is implemented
# For now, we'll create mock ind_meta tables to test mutate_ind_meta()

# Helper function to create a mock population with ind_meta table
create_mock_pop_with_ind_meta <- function(n_males = 10, n_females = 100) {
  temp_db <- tempfile(fileext = ".duckdb")

  pop <- initialize_genome(
    pop_name = "test_pop",
    n_loci = 100,
    n_chr = 2,
    chr_len_Mb = 50,
    db_path = temp_db
  )

  # Create ind_meta table manually (simulating what add_founders will do)
  n_ind <- n_males + n_females
  ind_meta <- tibble::tibble(
    ind_id = paste0("ind_", seq_len(n_ind)),
    parent_1 = "0",
    parent_2 = "0",
    population = "A",
    sex = c(rep("M", n_males), rep("F", n_females))
  )

  DBI::dbWriteTable(pop$db_conn, "ind_meta", ind_meta, overwrite = TRUE)

  # Update pop object tables list
  pop$tables <- c(pop$tables, "ind_meta")

  return(pop)
}


# 1. Basic Functionality Tests ----

test_that("mutate_ind_meta adds single scalar column", {
  pop <- create_mock_pop_with_ind_meta()

  pop <- pop %>% mutate_ind_meta(gen = 0)

  ind_meta <- get_table(pop, "ind_meta") %>% dplyr::collect()

  expect_true("gen" %in% colnames(ind_meta))
  expect_true(all(ind_meta$gen == 0))

  close_pop(pop)
  unlink(pop$db_path)
})


test_that("mutate_ind_meta adds multiple scalar columns", {
  pop <- create_mock_pop_with_ind_meta()

  pop <- pop %>% mutate_ind_meta(
    gen = 0,
    farm = "A",
    is_selected = FALSE
  )

  ind_meta <- get_table(pop, "ind_meta") %>% dplyr::collect()

  expect_true(all(c("gen", "farm", "is_selected") %in% colnames(ind_meta)))
  expect_true(all(ind_meta$gen == 0))
  expect_true(all(ind_meta$farm == "A"))
  expect_true(all(ind_meta$is_selected == FALSE))

  close_pop(pop)
  unlink(pop$db_path)
})


test_that("mutate_ind_meta adds vector column", {
  pop <- create_mock_pop_with_ind_meta(n_males = 10, n_females = 100)
  n_ind <- 110

  weight_values <- rnorm(n_ind, mean = 80, sd = 10)
  pop <- pop %>% mutate_ind_meta(weight_kg = weight_values)

  ind_meta <- get_table(pop, "ind_meta") %>% dplyr::collect()

  expect_true("weight_kg" %in% colnames(ind_meta))
  expect_equal(nrow(ind_meta), n_ind)
  expect_equal(ind_meta$weight_kg, weight_values)

  close_pop(pop)
  unlink(pop$db_path)
})


test_that("mutate_ind_meta updates existing column", {
  pop <- create_mock_pop_with_ind_meta()

  # Add column with value 0
  pop <- pop %>% mutate_ind_meta(gen = 0)

  # Update to value 1
  pop <- pop %>% mutate_ind_meta(gen = 1)

  ind_meta <- get_table(pop, "ind_meta") %>% dplyr::collect()

  expect_true(all(ind_meta$gen == 1))

  close_pop(pop)
  unlink(pop$db_path)
})


# 2. Type Handling Tests ----

test_that("mutate_ind_meta handles logical/BOOLEAN type", {
  pop <- create_mock_pop_with_ind_meta()

  pop <- pop %>% mutate_ind_meta(is_selected = TRUE)

  ind_meta <- get_table(pop, "ind_meta") %>% dplyr::collect()

  expect_true("is_selected" %in% colnames(ind_meta))
  expect_type(ind_meta$is_selected, "logical")
  expect_true(all(ind_meta$is_selected == TRUE))

  close_pop(pop)
  unlink(pop$db_path)
})


test_that("mutate_ind_meta handles integer/INTEGER type", {
  pop <- create_mock_pop_with_ind_meta()

  pop <- pop %>% mutate_ind_meta(gen = 0L)

  ind_meta <- get_table(pop, "ind_meta") %>% dplyr::collect()

  expect_true("gen" %in% colnames(ind_meta))
  expect_type(ind_meta$gen, "integer")

  close_pop(pop)
  unlink(pop$db_path)
})


test_that("mutate_ind_meta handles numeric/DOUBLE type", {
  pop <- create_mock_pop_with_ind_meta()

  pop <- pop %>% mutate_ind_meta(weight = 75.5)

  ind_meta <- get_table(pop, "ind_meta") %>% dplyr::collect()

  expect_true("weight" %in% colnames(ind_meta))
  expect_type(ind_meta$weight, "double")
  expect_equal(ind_meta$weight[1], 75.5)

  close_pop(pop)
  unlink(pop$db_path)
})


test_that("mutate_ind_meta handles Date/DATE type", {
  pop <- create_mock_pop_with_ind_meta()

  test_date <- as.Date("2024-01-15")
  pop <- pop %>% mutate_ind_meta(date_birth = test_date)

  ind_meta <- get_table(pop, "ind_meta") %>% dplyr::collect()

  expect_true("date_birth" %in% colnames(ind_meta))
  expect_s3_class(ind_meta$date_birth, "Date")
  expect_equal(ind_meta$date_birth[1], test_date)

  close_pop(pop)
  unlink(pop$db_path)
})


test_that("mutate_ind_meta handles character/VARCHAR type", {
  pop <- create_mock_pop_with_ind_meta()

  pop <- pop %>% mutate_ind_meta(farm = "FarmA")

  ind_meta <- get_table(pop, "ind_meta") %>% dplyr::collect()

  expect_true("farm" %in% colnames(ind_meta))
  expect_type(ind_meta$farm, "character")
  expect_equal(ind_meta$farm[1], "FarmA")

  close_pop(pop)
  unlink(pop$db_path)
})


# 3. Error Handling Tests ----

test_that("mutate_ind_meta errors if ind_meta doesn't exist", {
  temp_db <- tempfile(fileext = ".duckdb")

  pop <- initialize_genome(
    pop_name = "test_pop",
    n_loci = 100,
    n_chr = 2,
    chr_len_Mb = 50,
    db_path = temp_db
  )

  expect_error(
    pop %>% mutate_ind_meta(gen = 0),
    "ind_meta table does not exist"
  )

  close_pop(pop)
  unlink(temp_db)
})


test_that("mutate_ind_meta errors on vector length mismatch", {
  pop <- create_mock_pop_with_ind_meta(n_males = 10, n_females = 100)

  # Provide vector of wrong length
  expect_error(
    pop %>% mutate_ind_meta(weight = rnorm(50)),
    "Length of 'weight' \\(50\\) does not match number of individuals \\(110\\)"
  )

  close_pop(pop)
  unlink(pop$db_path)
})


test_that("mutate_ind_meta errors on reserved column names", {
  pop <- create_mock_pop_with_ind_meta()

  expect_error(
    pop %>% mutate_ind_meta(ind_id = "new_id"),
    "Cannot modify reserved column 'ind_id'"
  )

  expect_error(
    pop %>% mutate_ind_meta(parent_1 = "parent"),
    "Cannot modify reserved column 'parent_1'"
  )

  expect_error(
    pop %>% mutate_ind_meta(parent_2 = "parent"),
    "Cannot modify reserved column 'parent_2'"
  )

  expect_error(
    pop %>% mutate_ind_meta(population = "B"),
    "Cannot modify reserved column 'population'"
  )

  expect_error(
    pop %>% mutate_ind_meta(sex = "F"),
    "Cannot modify reserved column 'sex'"
  )

  close_pop(pop)
  unlink(pop$db_path)
})


test_that("mutate_ind_meta errors on invalid field names", {
  pop <- create_mock_pop_with_ind_meta()

  # Field name starting with number
  expect_error(
    pop %>% mutate_ind_meta(`1gen` = 0),
    "Invalid field name"
  )

  # Field name with special characters
  expect_error(
    pop %>% mutate_ind_meta(`my-field` = 0),
    "Invalid field name"
  )

  # Field name with spaces
  expect_error(
    pop %>% mutate_ind_meta(`my field` = 0),
    "Invalid field name"
  )

  close_pop(pop)
  unlink(pop$db_path)
})


test_that("mutate_ind_meta errors on SQL reserved words", {
  pop <- create_mock_pop_with_ind_meta()

  expect_error(
    pop %>% mutate_ind_meta(SELECT = 0),
    "SQL reserved keyword"
  )

  expect_error(
    pop %>% mutate_ind_meta(WHERE = 0),
    "SQL reserved keyword"
  )

  close_pop(pop)
  unlink(pop$db_path)
})


test_that("mutate_ind_meta errors on unsupported types", {
  pop <- create_mock_pop_with_ind_meta()

  # List type
  expect_error(
    pop %>% mutate_ind_meta(mylist = list(1, 2, 3)),
    "Unsupported type"
  )

  close_pop(pop)
  unlink(pop$db_path)
})


# 4. Edge Cases Tests ----

test_that("mutate_ind_meta handles NA scalar values", {
  pop <- create_mock_pop_with_ind_meta()

  expect_warning(
    pop <- pop %>% mutate_ind_meta(field1 = NA),
    "Cannot infer type from NA"
  )

  ind_meta <- get_table(pop, "ind_meta") %>% dplyr::collect()

  expect_true("field1" %in% colnames(ind_meta))
  expect_true(all(is.na(ind_meta$field1)))

  close_pop(pop)
  unlink(pop$db_path)
})


test_that("mutate_ind_meta handles empty ind_meta table", {
  temp_db <- tempfile(fileext = ".duckdb")

  pop <- initialize_genome(
    pop_name = "test_pop",
    n_loci = 100,
    n_chr = 2,
    chr_len_Mb = 50,
    db_path = temp_db
  )

  # Create empty ind_meta table
  DBI::dbExecute(
    pop$db_conn,
    "CREATE TABLE ind_meta (ind_id VARCHAR, parent_1 VARCHAR, parent_2 VARCHAR, population VARCHAR, sex VARCHAR)"
  )
  pop$tables <- c(pop$tables, "ind_meta")

  expect_warning(
    pop %>% mutate_ind_meta(gen = 0),
    "ind_meta table is empty"
  )

  close_pop(pop)
  unlink(temp_db)
})


test_that("mutate_ind_meta handles no fields specified", {
  pop <- create_mock_pop_with_ind_meta()

  expect_warning(
    result <- pop %>% mutate_ind_meta(),
    "No fields specified"
  )

  close_pop(pop)
  unlink(pop$db_path)
})


test_that("mutate_ind_meta handles special characters in strings", {
  pop <- create_mock_pop_with_ind_meta()

  # String with single quotes
  pop <- pop %>% mutate_ind_meta(farm = "O'Brien's Farm")

  ind_meta <- get_table(pop, "ind_meta") %>% dplyr::collect()

  expect_equal(ind_meta$farm[1], "O'Brien's Farm")

  close_pop(pop)
  unlink(pop$db_path)
})


test_that("mutate_ind_meta handles mixed scalar and vector in single call", {
  pop <- create_mock_pop_with_ind_meta(n_males = 10, n_females = 100)
  n_ind <- 110

  weight_values <- rnorm(n_ind, mean = 80, sd = 10)

  pop <- pop %>% mutate_ind_meta(
    gen = 0,                    # scalar
    weight_kg = weight_values   # vector
  )

  ind_meta <- get_table(pop, "ind_meta") %>% dplyr::collect()

  expect_true(all(ind_meta$gen == 0))
  expect_equal(ind_meta$weight_kg, weight_values)

  close_pop(pop)
  unlink(pop$db_path)
})


# 5. Integration Tests ----

test_that("mutate_ind_meta works in pipe workflow", {
  temp_db <- tempfile(fileext = ".duckdb")

  # Simulate complete workflow (without actual add_founders)
  pop <- initialize_genome(
    pop_name = "test_pop",
    n_loci = 100,
    n_chr = 2,
    chr_len_Mb = 50,
    db_path = temp_db
  )

  # Manually create ind_meta
  ind_meta <- tibble::tibble(
    ind_id = paste0("ind_", 1:110),
    parent_1 = "0",
    parent_2 = "0",
    population = "A",
    sex = c(rep("M", 10), rep("F", 100))
  )
  DBI::dbWriteTable(pop$db_conn, "ind_meta", ind_meta, overwrite = TRUE)
  pop$tables <- c(pop$tables, "ind_meta")

  # Test piping
  pop <- pop %>%
    mutate_ind_meta(
      gen = 0,
      farm = "A"
    ) %>%
    mutate_ind_meta(
      date_birth = Sys.Date()
    )

  ind_meta_result <- get_table(pop, "ind_meta") %>% dplyr::collect()

  expect_true(all(c("gen", "farm", "date_birth") %in% colnames(ind_meta_result)))

  close_pop(pop)
  unlink(temp_db)
})


test_that("mutate_ind_meta preserves existing data", {
  pop <- create_mock_pop_with_ind_meta()

  # Add first set of fields
  pop <- pop %>% mutate_ind_meta(
    gen = 0,
    farm = "A",
    weight = 75.5
  )

  # Add more fields
  pop <- pop %>% mutate_ind_meta(
    date_birth = Sys.Date(),
    is_selected = TRUE
  )

  ind_meta <- get_table(pop, "ind_meta") %>% dplyr::collect()

  # Check all fields present
  expect_true(all(c("gen", "farm", "weight", "date_birth", "is_selected") %in% colnames(ind_meta)))

  # Check original values preserved
  expect_true(all(ind_meta$gen == 0))
  expect_true(all(ind_meta$farm == "A"))
  expect_true(all(ind_meta$weight == 75.5))

  close_pop(pop)
  unlink(pop$db_path)
})


test_that("mutate_ind_meta can be called multiple times", {
  pop <- create_mock_pop_with_ind_meta()

  pop <- pop %>% mutate_ind_meta(gen = 0)
  pop <- pop %>% mutate_ind_meta(farm = "A")
  pop <- pop %>% mutate_ind_meta(weight = 75.5)

  ind_meta <- get_table(pop, "ind_meta") %>% dplyr::collect()

  expect_true(all(c("gen", "farm", "weight") %in% colnames(ind_meta)))
  expect_equal(nrow(ind_meta), 110)

  close_pop(pop)
  unlink(pop$db_path)
})
