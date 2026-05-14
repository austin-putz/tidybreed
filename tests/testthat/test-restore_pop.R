test_that("restore_pop() reconstructs metadata from a freshly initialized DB", {
  tmp <- tempfile(fileext = ".duckdb")
  pop <- initialize_genome(
    pop_name     = "test_pop",
    n_loci       = 100,
    n_chr        = 2,
    chr_len_Mb   = 50,
    n_haplotypes = 30,
    db_path      = tmp
  )
  close_pop(pop)

  pop2 <- restore_pop(tmp, pop_name = "test_pop")

  expect_s3_class(pop2, "tidybreed_pop")
  expect_equal(pop2$pop_name,            "test_pop")
  expect_equal(pop2$db_path,             tmp)
  expect_equal(pop2$metadata$n_loci,     100L)
  expect_equal(pop2$metadata$n_chr,      2L)
  expect_equal(length(pop2$metadata$chr_len_Mb), 2L)
  expect_equal(pop2$metadata$chr_names,  c("1", "2"))
  expect_equal(pop2$metadata$n_individuals, 0L)

  close_pop(pop2)
  unlink(tmp)
})


test_that("restore_pop() derives n_individuals from ind_meta", {
  tmp <- tempfile(fileext = ".duckdb")
  pop <- initialize_genome(
    pop_name = "count_test", n_loci = 50, n_chr = 2,
    chr_len_Mb = 50, n_haplotypes = 20, db_path = tmp
  )
  pop <- add_founders(pop, n_males = 10, n_females = 20, line_name = "A")
  close_pop(pop)

  pop2 <- restore_pop(tmp)

  expect_equal(pop2$metadata$n_individuals, 30L)

  close_pop(pop2)
  unlink(tmp)
})


test_that("restore_pop() preserves heterogeneous chr_len_Mb", {
  tmp <- tempfile(fileext = ".duckdb")
  pop <- initialize_genome(
    pop_name     = "cattle",
    n_loci       = 90,
    n_chr        = 3,
    chr_len_Mb   = c(158, 137, 121),
    n_haplotypes = 20,
    db_path      = tmp
  )
  close_pop(pop)

  pop2 <- restore_pop(tmp)

  expect_equal(pop2$metadata$n_chr, 3L)
  expect_equal(pop2$metadata$chr_names, c("1", "2", "3"))
  # chr_len_Mb from max(pos_Mb) will be <= original; just check ordering
  expect_true(pop2$metadata$chr_len_Mb[1] >= pop2$metadata$chr_len_Mb[2])
  expect_true(pop2$metadata$chr_len_Mb[2] >= pop2$metadata$chr_len_Mb[3])

  close_pop(pop2)
  unlink(tmp)
})


test_that("add_founders() works after restore_pop() on empty DB", {
  tmp <- tempfile(fileext = ".duckdb")
  pop <- initialize_genome(
    pop_name = "fresh", n_loci = 50, n_chr = 2,
    chr_len_Mb = 50, n_haplotypes = 20, db_path = tmp
  )
  close_pop(pop)

  pop2 <- restore_pop(tmp)

  expect_no_error({
    pop2 <- add_founders(pop2, n_males = 5, n_females = 5, line_name = "A")
  })

  ind_meta <- get_table(pop2, "ind_meta") |> dplyr::collect()
  expect_equal(nrow(ind_meta), 10L)

  close_pop(pop2)
  unlink(tmp)
})


test_that("add_offspring() works after restore; IDs continue without collision", {
  tmp <- tempfile(fileext = ".duckdb")
  pop <- initialize_genome(
    pop_name = "offspring_test", n_loci = 50, n_chr = 2,
    chr_len_Mb = 50, n_haplotypes = 20, db_path = tmp
  )
  pop <- add_founders(pop, n_males = 5, n_females = 5, line_name = "A")
  close_pop(pop)

  pop2 <- restore_pop(tmp)

  matings <- tibble::tibble(
    id_parent_1 = "A_1",
    id_parent_2 = "A_6",
    sex         = "M",
    line        = "A"
  )

  expect_no_error({
    pop2 <- add_offspring(pop2, matings)
  })

  ind_meta <- get_table(pop2, "ind_meta") |> dplyr::collect()
  expect_equal(nrow(ind_meta), 11L)
  expect_equal(length(unique(ind_meta$id_ind)), 11L)
  expect_true("A_11" %in% ind_meta$id_ind)

  close_pop(pop2)
  unlink(tmp)
})


test_that("pop_name inferred from _tidybreed.duckdb suffix", {
  tmp_dir  <- tempdir()
  db_path  <- file.path(tmp_dir, "mysim_tidybreed.duckdb")

  pop <- initialize_genome(
    pop_name = "mysim", n_loci = 50, n_chr = 2,
    chr_len_Mb = 50, n_haplotypes = 10, db_path = db_path
  )
  close_pop(pop)

  pop2 <- restore_pop(db_path)
  expect_equal(pop2$pop_name, "mysim")

  close_pop(pop2)
  unlink(db_path)
})


test_that("pop_name override parameter works", {
  tmp <- tempfile(fileext = ".duckdb")
  pop <- initialize_genome(
    pop_name = "original", n_loci = 50, n_chr = 2,
    chr_len_Mb = 50, n_haplotypes = 10, db_path = tmp
  )
  close_pop(pop)

  pop2 <- restore_pop(tmp, pop_name = "replicate_01")
  expect_equal(pop2$pop_name, "replicate_01")

  close_pop(pop2)
  unlink(tmp)
})


test_that("restore_pop() errors when file does not exist", {
  expect_error(
    restore_pop("/tmp/definitely_does_not_exist_xyz123.duckdb"),
    "does not exist"
  )
})


test_that("restore_pop() errors for ':memory:' path", {
  expect_error(
    restore_pop(":memory:"),
    "cannot restore in-memory databases"
  )
})


test_that("restore_pop() errors when genome_meta table is missing", {
  tmp  <- tempfile(fileext = ".duckdb")
  drv  <- duckdb::duckdb()
  conn <- DBI::dbConnect(drv, dbdir = tmp)
  DBI::dbExecute(conn, "CREATE TABLE dummy (x INTEGER)")
  DBI::dbDisconnect(conn, shutdown = TRUE)

  expect_error(
    restore_pop(tmp),
    "genome_meta"
  )

  unlink(tmp)
})
