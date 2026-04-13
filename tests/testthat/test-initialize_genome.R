test_that("initialize_genome creates population object", {

  # Create temporary database
  temp_db <- tempfile(fileext = ".duckdb")

  # Initialize genome
  pop <- initialize_genome(
    pop_name = "test_pop",
    n_loci = 100,
    n_chr = 2,
    chr_len_Mb = 50,
    db_path = temp_db
  )

  # Check object class
  expect_s3_class(pop, "tidybreed_pop")

  # Check components
  expect_equal(pop$pop_name, "test_pop")
  expect_equal(pop$db_path, temp_db)
  expect_true(file.exists(temp_db))

  # Check tables exist
  expect_true(all(c("genome_meta", "genome_haplotype", "genome_genotype") %in% pop$tables))

  # Check genome_meta has correct number of rows
  genome_meta <- get_table(pop, "genome_meta") %>% dplyr::collect()
  expect_equal(nrow(genome_meta), 100)
  expect_equal(ncol(genome_meta), 5)  # locus_id, locus_name, chr, chr_name, pos_Mb

  # Check chromosome assignment
  expect_equal(unique(genome_meta$chr), c(1, 2))

  # Clean up
  close_pop(pop)
  unlink(temp_db)
})


test_that("initialize_genome handles different chromosome lengths", {

  temp_db <- tempfile(fileext = ".duckdb")

  pop <- initialize_genome(
    pop_name = "test_pop",
    n_loci = 100,
    n_chr = 3,
    chr_len_Mb = c(100, 50, 75),
    db_path = temp_db
  )

  genome_meta <- get_table(pop, "genome_meta") %>% dplyr::collect()

  # Check max position per chromosome
  chr_max <- genome_meta %>%
    dplyr::group_by(chr) %>%
    dplyr::summarise(max_pos = max(pos_Mb)) %>%
    dplyr::arrange(chr)

  expect_true(all(chr_max$max_pos <= c(100, 50, 75)))
  expect_true(all(chr_max$max_pos > c(90, 45, 70)))  # Close to max

  close_pop(pop)
  unlink(temp_db)
})


test_that("initialize_genome prevents overwriting without permission", {

  temp_db <- tempfile(fileext = ".duckdb")

  # Create first population
  pop1 <- initialize_genome(
    pop_name = "test1",
    n_loci = 10,
    n_chr = 1,
    chr_len_Mb = 100,
    db_path = temp_db
  )

  close_pop(pop1)

  # Try to create second without overwrite
  expect_error(
    initialize_genome(
      pop_name = "test2",
      n_loci = 20,
      n_chr = 1,
      chr_len_Mb = 100,
      db_path = temp_db,
      overwrite = FALSE
    ),
    "already exists"
  )

  # Should work with overwrite = TRUE
  expect_message(
    pop2 <- initialize_genome(
      pop_name = "test2",
      n_loci = 20,
      n_chr = 1,
      chr_len_Mb = 100,
      db_path = temp_db,
      overwrite = TRUE
    ),
    "Overwriting"
  )

  # Check new population has correct loci
  genome_meta <- get_table(pop2, "genome_meta") %>% dplyr::collect()
  expect_equal(nrow(genome_meta), 20)

  close_pop(pop2)
  unlink(temp_db)
})


test_that("initialize_genome works with custom locus names", {

  temp_db <- tempfile(fileext = ".duckdb")

  custom_names <- paste0("rs", 1:50)

  pop <- initialize_genome(
    pop_name = "test",
    n_loci = 50,
    n_chr = 1,
    chr_len_Mb = 100,
    db_path = temp_db,
    locus_names = custom_names
  )

  genome_meta <- get_table(pop, "genome_meta") %>% dplyr::collect()

  expect_equal(genome_meta$locus_name, custom_names)

  close_pop(pop)
  unlink(temp_db)
})


test_that("get_table returns lazy tibble", {

  temp_db <- tempfile(fileext = ".duckdb")

  pop <- initialize_genome(
    pop_name = "test",
    n_loci = 100,
    n_chr = 2,
    chr_len_Mb = 100,
    db_path = temp_db
  )

  genome <- get_table(pop, "genome_meta")

  expect_s3_class(genome, "tbl_duckdb_connection")

  # Test that dplyr verbs work
  chr1_loci <- genome %>%
    dplyr::filter(chr == 1) %>%
    dplyr::collect()

  expect_true(all(chr1_loci$chr == 1))
  expect_true(nrow(chr1_loci) > 0)

  close_pop(pop)
  unlink(temp_db)
})


test_that("print method works", {

  temp_db <- tempfile(fileext = ".duckdb")

  pop <- initialize_genome(
    pop_name = "TestPop",
    n_loci = 100,
    n_chr = 2,
    chr_len_Mb = 100,
    db_path = temp_db
  )

  expect_output(print(pop), "tidybreed population")
  expect_output(print(pop), "TestPop")
  expect_output(print(pop), "Loci: 100")

  close_pop(pop)
  unlink(temp_db)
})
