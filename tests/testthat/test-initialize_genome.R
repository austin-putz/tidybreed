test_that("initialize_genome() emits deprecation warning", {
  expect_warning(
    pop <- initialize_genome("t", 100, 2, 100, db_path = ":memory:"),
    regexp = "deprecated"
  )
  close_pop(pop)
})


test_that("initialize_genome() deprecated wrapper produces correct tables", {
  pop <- suppressWarnings(
    initialize_genome("t", 100, 2, 100, db_path = ":memory:")
  )
  on.exit(close_pop(pop))

  expect_s3_class(pop, "tidybreed_pop")
  expect_true("genome_meta" %in% pop$tables)
  expect_true("ind_meta"    %in% pop$tables)

  gm <- get_table(pop, "genome_meta") |> dplyr::collect()
  expect_equal(nrow(gm), 100)
})


test_that("initialize_genome() deprecated wrapper produces same tables as open_pop() |> define_genome()", {
  pop_new <- open_pop(pop_name = "t", db_name = ":memory:") |>
    define_genome(n_loci = 100, n_chr = 2, chr_len_Mb = 100)
  on.exit(close_pop(pop_new))

  pop_old <- suppressWarnings(
    initialize_genome("t", 100, 2, 100, db_path = ":memory:")
  )
  on.exit(close_pop(pop_old), add = TRUE)

  # Both should expose the same set of tables
  expect_setequal(pop_new$tables, pop_old$tables)

  # Both should have 100 loci in genome_meta
  gm_new <- get_table(pop_new, "genome_meta") |> dplyr::collect()
  gm_old <- get_table(pop_old, "genome_meta") |> dplyr::collect()
  expect_equal(nrow(gm_new), nrow(gm_old))
  expect_equal(ncol(gm_new), ncol(gm_old))
})


test_that("initialize_genome() deprecated wrapper works with overwrite = TRUE", {
  tmp <- tempfile(fileext = ".duckdb")
  on.exit(unlink(tmp))

  pop1 <- suppressWarnings(
    initialize_genome("t", 50, 1, 100, db_path = tmp)
  )
  close_pop(pop1)

  # overwrite = TRUE -> clean = TRUE in open_pop
  expect_warning(
    pop2 <- initialize_genome("t", 80, 1, 100, db_path = tmp, overwrite = TRUE),
    regexp = "deprecated"
  )
  on.exit(close_pop(pop2), add = TRUE)

  gm <- get_table(pop2, "genome_meta") |> dplyr::collect()
  expect_equal(nrow(gm), 80)
})
