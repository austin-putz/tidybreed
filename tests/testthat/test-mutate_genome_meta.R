test_that("mutate_genome_meta adds scalar field", {
  # Initialize genome
  pop <- initialize_genome(
    pop_name = "test_mgm_scalar",
    n_loci = 100,
    n_chr = 5,
    chr_len_Mb = 100,
    db_path = ":memory:"
  )

  # Add scalar field
  pop <- mutate_genome_meta(pop, is_QTL = FALSE)

  # Check column exists
  cols <- DBI::dbListFields(pop$db_conn, "genome_meta")
  expect_true("is_QTL" %in% cols)

  # Check all values are FALSE
  result <- get_table(pop, "genome_meta") %>%
    dplyr::select(is_QTL) %>%
    dplyr::collect()

  expect_equal(nrow(result), 100)
  expect_true(all(result$is_QTL == FALSE))

  close_pop(pop)
})


test_that("mutate_genome_meta adds vector field", {
  pop <- initialize_genome(
    pop_name = "test_mgm_vector",
    n_loci = 100,
    n_chr = 5,
    chr_len_Mb = 100,
    db_path = ":memory:"
  )

  # Create vector
  maf_values <- runif(100, 0.01, 0.5)

  # Add vector field
  pop <- mutate_genome_meta(pop, maf = maf_values)

  # Check column exists
  cols <- DBI::dbListFields(pop$db_conn, "genome_meta")
  expect_true("maf" %in% cols)

  # Check values match
  result <- get_table(pop, "genome_meta") %>%
    dplyr::select(locus_id, maf) %>%
    dplyr::arrange(locus_id) %>%
    dplyr::collect()

  expect_equal(nrow(result), 100)
  expect_equal(result$maf, maf_values, tolerance = 1e-10)

  close_pop(pop)
})


test_that("mutate_genome_meta adds multiple fields", {
  pop <- initialize_genome(
    pop_name = "test_mgm_multi",
    n_loci = 100,
    n_chr = 5,
    chr_len_Mb = 100,
    db_path = ":memory:"
  )

  # Add multiple fields
  pop <- mutate_genome_meta(
    pop,
    is_50k = FALSE,
    is_HD = TRUE,
    maf = 0.25
  )

  # Check all columns exist
  cols <- DBI::dbListFields(pop$db_conn, "genome_meta")
  expect_true(all(c("is_50k", "is_HD", "maf") %in% cols))

  # Check values
  result <- get_table(pop, "genome_meta") %>%
    dplyr::select(is_50k, is_HD, maf) %>%
    dplyr::collect()

  expect_true(all(result$is_50k == FALSE))
  expect_true(all(result$is_HD == TRUE))
  expect_true(all(result$maf == 0.25))

  close_pop(pop)
})


test_that("mutate_genome_meta updates existing field", {
  pop <- initialize_genome(
    pop_name = "test_mgm_update",
    n_loci = 100,
    n_chr = 5,
    chr_len_Mb = 100,
    db_path = ":memory:"
  )

  # Add field
  pop <- mutate_genome_meta(pop, is_QTL = FALSE)

  # Update field
  pop <- mutate_genome_meta(pop, is_QTL = TRUE)

  # Check updated value
  result <- get_table(pop, "genome_meta") %>%
    dplyr::select(is_QTL) %>%
    dplyr::collect()

  expect_true(all(result$is_QTL == TRUE))

  close_pop(pop)
})


test_that("mutate_genome_meta handles different types", {
  pop <- initialize_genome(
    pop_name = "test_mgm_types",
    n_loci = 100,
    n_chr = 5,
    chr_len_Mb = 100,
    db_path = ":memory:"
  )

  # Add fields of different types
  pop <- mutate_genome_meta(
    pop,
    is_functional = TRUE,                      # logical
    effect_size = 0.5,                         # numeric
    gene_count = 3L,                           # integer
    gene_name = "GENE1"                        # character
  )

  # Check all columns exist and have correct types
  result <- get_table(pop, "genome_meta") %>%
    dplyr::select(is_functional, effect_size, gene_count, gene_name) %>%
    dplyr::collect()

  expect_type(result$is_functional, "logical")
  expect_type(result$effect_size, "double")
  expect_type(result$gene_count, "integer")
  expect_type(result$gene_name, "character")

  close_pop(pop)
})


test_that("mutate_genome_meta errors on reserved columns", {
  pop <- initialize_genome(
    pop_name = "test_mgm_reserved",
    n_loci = 100,
    n_chr = 5,
    chr_len_Mb = 100,
    db_path = ":memory:"
  )

  # Try to modify reserved columns
  expect_error(
    mutate_genome_meta(pop, locus_id = 999),
    "Cannot modify reserved column 'locus_id'"
  )

  expect_error(
    mutate_genome_meta(pop, locus_name = "BAD"),
    "Cannot modify reserved column 'locus_name'"
  )

  expect_error(
    mutate_genome_meta(pop, chr = 99),
    "Cannot modify reserved column 'chr'"
  )

  expect_error(
    mutate_genome_meta(pop, chr_name = "BAD"),
    "Cannot modify reserved column 'chr_name'"
  )

  expect_error(
    mutate_genome_meta(pop, pos_Mb = 999),
    "Cannot modify reserved column 'pos_Mb'"
  )

  close_pop(pop)
})


test_that("mutate_genome_meta errors on vector length mismatch", {
  pop <- initialize_genome(
    pop_name = "test_mgm_length",
    n_loci = 100,
    n_chr = 5,
    chr_len_Mb = 100,
    db_path = ":memory:"
  )

  # Try to add vector with wrong length
  expect_error(
    mutate_genome_meta(pop, maf = runif(50)),  # Only 50 values for 100 loci
    "does not match number of loci"
  )

  expect_error(
    mutate_genome_meta(pop, maf = runif(200)),  # Too many values
    "does not match number of loci"
  )

  close_pop(pop)
})


test_that("mutate_genome_meta errors on invalid field names", {
  pop <- initialize_genome(
    pop_name = "test_mgm_invalid",
    n_loci = 100,
    n_chr = 5,
    chr_len_Mb = 100,
    db_path = ":memory:"
  )

  # Invalid characters
  expect_error(
    mutate_genome_meta(pop, `my field` = 1),
    "Invalid field name"
  )

  # Starts with number
  expect_error(
    mutate_genome_meta(pop, `1field` = 1),
    "Invalid field name"
  )

  # SQL keyword
  expect_error(
    mutate_genome_meta(pop, SELECT = 1),
    "SQL reserved keyword"
  )

  close_pop(pop)
})


test_that("mutate_genome_meta warns on no fields", {
  pop <- initialize_genome(
    pop_name = "test_mgm_nofields",
    n_loci = 100,
    n_chr = 5,
    chr_len_Mb = 100,
    db_path = ":memory:"
  )

  # No fields specified
  expect_warning(
    mutate_genome_meta(pop),
    "No fields specified"
  )

  close_pop(pop)
})


test_that("mutate_genome_meta errors if genome_meta doesn't exist", {
  # Create pop without genome_meta (not possible with normal workflow, but test edge case)
  pop <- list(
    db_conn = NULL,
    tables = c(),
    metadata = list()
  )
  class(pop) <- "tidybreed_pop"

  expect_error(
    mutate_genome_meta(pop, test = 1),
    "genome_meta table does not exist"
  )
})
