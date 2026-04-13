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


test_that("initialize_genome creates founder haplotypes with uniform frequencies", {

  temp_db <- tempfile(fileext = ".duckdb")

  pop <- initialize_genome(
    pop_name = "test_founders",
    n_loci = 100,
    n_chr = 2,
    chr_len_Mb = 50,
    db_path = temp_db,
    n_haplotypes = 50,
    min_allele_freq = 0.1,
    max_allele_freq = 0.9
  )

  # Check founder_haplotypes table exists
  expect_true("founder_haplotypes" %in% pop$tables)

  # Check table dimensions
  founder_haps <- get_table(pop, "founder_haplotypes") %>% dplyr::collect()
  expect_equal(nrow(founder_haps), 50)
  expect_equal(ncol(founder_haps), 101)  # hap_id + 100 loci

  # Check haplotype IDs
  expect_equal(founder_haps$hap_id, paste0("hap_", 1:50))

  # Check alleles are 0 or 1
  locus_cols <- paste0("locus_", 1:100)
  alleles <- as.matrix(founder_haps[, locus_cols])
  expect_true(all(alleles %in% c(0, 1)))

  # Check genome_meta has founder_allele_freq column
  genome_meta <- get_table(pop, "genome_meta") %>% dplyr::collect()
  expect_true("founder_allele_freq" %in% colnames(genome_meta))
  expect_equal(length(genome_meta$founder_allele_freq), 100)

  # Check allele frequencies are within bounds
  expect_true(all(genome_meta$founder_allele_freq >= 0.1))
  expect_true(all(genome_meta$founder_allele_freq <= 0.9))

  # Check metadata
  expect_equal(pop$metadata$n_haplotypes, 50)
  expect_equal(pop$metadata$allele_freq_dist, "uniform")

  close_pop(pop)
  unlink(temp_db)
})


test_that("initialize_genome creates founder haplotypes with fixed frequency", {

  temp_db <- tempfile(fileext = ".duckdb")

  pop <- initialize_genome(
    pop_name = "test_fixed",
    n_loci = 50,
    n_chr = 2,
    chr_len_Mb = 50,
    db_path = temp_db,
    n_haplotypes = 30,
    fixed_allele_freq = 0.5
  )

  # Check founder_haplotypes table exists
  expect_true("founder_haplotypes" %in% pop$tables)

  # Check genome_meta has founder_allele_freq column with all 0.5
  genome_meta <- get_table(pop, "genome_meta") %>% dplyr::collect()
  expect_true(all(genome_meta$founder_allele_freq == 0.5))

  # Check metadata
  expect_equal(pop$metadata$n_haplotypes, 30)
  expect_equal(pop$metadata$allele_freq_dist, "fixed")
  expect_equal(pop$metadata$fixed_allele_freq, 0.5)

  close_pop(pop)
  unlink(temp_db)
})


test_that("initialize_genome works without founder haplotypes", {

  temp_db <- tempfile(fileext = ".duckdb")

  pop <- initialize_genome(
    pop_name = "test_no_founders",
    n_loci = 100,
    n_chr = 2,
    chr_len_Mb = 50,
    db_path = temp_db
  )

  # Check founder_haplotypes table does NOT exist
  expect_false("founder_haplotypes" %in% pop$tables)

  # Check genome_meta does NOT have founder_allele_freq column
  genome_meta <- get_table(pop, "genome_meta") %>% dplyr::collect()
  expect_false("founder_allele_freq" %in% colnames(genome_meta))

  # Check metadata does not have founder info
  expect_null(pop$metadata$n_haplotypes)

  close_pop(pop)
  unlink(temp_db)
})
