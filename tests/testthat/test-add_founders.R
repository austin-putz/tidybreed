test_that("add_founders creates ind_meta table with correct structure", {
  # Create test population
  pop <- initialize_genome(
    pop_name = "test",
    n_loci = 100,
    n_chr = 5,
    chr_len_Mb = 100,
    n_haplotypes = 50,
    db_path = ":memory:"
  )

  # Add founders
  pop <- add_founders(pop, n_males = 10, n_females = 20, line_name = "A")

  # Check ind_meta table exists
  expect_true("ind_meta" %in% pop$tables)
  expect_true("ind_meta" %in% DBI::dbListTables(pop$db_conn))

  # Read ind_meta
  ind_meta <- get_table(pop, "ind_meta") %>% dplyr::collect()

  # Check structure
  expect_equal(nrow(ind_meta), 30)
  expect_equal(ncol(ind_meta), 5)
  expect_true(all(c("id_ind", "id_parent_1", "id_parent_2", "line", "sex") %in% colnames(ind_meta)))

  # Check data types
  expect_type(ind_meta$id_ind, "character")
  expect_type(ind_meta$id_parent_1, "character")
  expect_type(ind_meta$id_parent_2, "character")
  expect_type(ind_meta$line, "character")
  expect_type(ind_meta$sex, "character")

  close_pop(pop)
})


test_that("add_founders correctly assigns IDs and sex", {
  pop <- initialize_genome(
    pop_name = "test",
    n_loci = 50,
    n_chr = 2,
    chr_len_Mb = 100,
    n_haplotypes = 20,
    db_path = ":memory:"
  )

  pop <- add_founders(pop, n_males = 5, n_females = 10, line_name = "A")

  ind_meta <- get_table(pop, "ind_meta") %>% dplyr::collect()

  # Check IDs
  expected_ids <- paste0("A-", 1:15)
  expect_equal(ind_meta$id_ind, expected_ids)

  # Check line
  expect_true(all(ind_meta$line == "A"))

  # Check sex assignment
  expect_equal(sum(ind_meta$sex == "M"), 5)
  expect_equal(sum(ind_meta$sex == "F"), 10)
  expect_equal(ind_meta$sex[1:5], rep("M", 5))
  expect_equal(ind_meta$sex[6:15], rep("F", 10))

  # Check parents are NULL
  expect_true(all(is.na(ind_meta$id_parent_1)))
  expect_true(all(is.na(ind_meta$id_parent_2)))

  close_pop(pop)
})


test_that("add_founders creates correct genome_haplotype table", {
  pop <- initialize_genome(
    pop_name = "test",
    n_loci = 50,
    n_chr = 2,
    chr_len_Mb = 100,
    n_haplotypes = 20,
    db_path = ":memory:"
  )

  pop <- add_founders(pop, n_males = 5, n_females = 5, line_name = "A")

  haps <- get_table(pop, "genome_haplotype") %>% dplyr::collect()

  # Check 2 rows per individual
  expect_equal(nrow(haps), 20)  # 10 individuals × 2

  # Check structure: id_ind, parent_origin, locus_1, ..., locus_50
  expect_equal(ncol(haps), 52)  # id_ind + parent_origin + 50 loci
  expect_true("id_ind" %in% colnames(haps))
  expect_true("parent_origin" %in% colnames(haps))

  # Check locus columns exist
  locus_cols <- paste0("locus_", 1:50)
  expect_true(all(locus_cols %in% colnames(haps)))

  # Check parent_origin values
  expect_true(all(haps$parent_origin %in% c(1L, 2L)))

  # Check each individual has exactly 2 rows
  for (i in 1:10) {
    ind_id <- paste0("A-", i)
    ind_haps <- haps %>% dplyr::filter(id_ind == !!ind_id)
    expect_equal(nrow(ind_haps), 2)
    expect_equal(sort(ind_haps$parent_origin), c(1L, 2L))
  }

  # Check haplotype values are 0 or 1
  hap_matrix <- as.matrix(haps[, locus_cols])
  expect_true(all(hap_matrix %in% c(0, 1)))

  close_pop(pop)
})


test_that("add_founders creates correct genome_genotype table", {
  pop <- initialize_genome(
    pop_name = "test",
    n_loci = 50,
    n_chr = 2,
    chr_len_Mb = 100,
    n_haplotypes = 20,
    db_path = ":memory:"
  )

  pop <- add_founders(pop, n_males = 5, n_females = 5, line_name = "A")

  genos <- get_table(pop, "genome_genotype") %>% dplyr::collect()

  # Check 1 row per individual
  expect_equal(nrow(genos), 10)

  # Check structure: id_ind, locus_1, ..., locus_50
  expect_equal(ncol(genos), 51)  # id_ind + 50 loci
  expect_true("id_ind" %in% colnames(genos))

  # Check locus columns exist
  locus_cols <- paste0("locus_", 1:50)
  expect_true(all(locus_cols %in% colnames(genos)))

  # Check genotype values are 0, 1, or 2
  geno_matrix <- as.matrix(genos[, locus_cols])
  expect_true(all(geno_matrix %in% c(0, 1, 2)))

  close_pop(pop)
})


test_that("genotypes equal sum of haplotypes", {
  pop <- initialize_genome(
    pop_name = "test",
    n_loci = 50,
    n_chr = 2,
    chr_len_Mb = 100,
    n_haplotypes = 20,
    db_path = ":memory:"
  )

  pop <- add_founders(pop, n_males = 5, n_females = 5, line_name = "A")

  haps <- get_table(pop, "genome_haplotype") %>% dplyr::collect()
  genos <- get_table(pop, "genome_genotype") %>% dplyr::collect()

  # Check for each individual
  for (i in 1:10) {
    ind_id <- paste0("A-", i)

    hap_rows <- haps %>% dplyr::filter(id_ind == !!ind_id)
    geno_row <- genos %>% dplyr::filter(id_ind == !!ind_id)

    # Check each locus
    for (j in 1:50) {
      locus_name <- paste0("locus_", j)

      hap1 <- hap_rows[[locus_name]][1]
      hap2 <- hap_rows[[locus_name]][2]
      geno <- geno_row[[locus_name]][1]

      expect_equal(geno, hap1 + hap2)
    }
  }

  close_pop(pop)
})


test_that("add_founders works with multiple lines", {
  pop <- initialize_genome(
    pop_name = "test",
    n_loci = 50,
    n_chr = 2,
    chr_len_Mb = 100,
    n_haplotypes = 20,
    db_path = ":memory:"
  )

  # Add line A
  pop <- add_founders(pop, n_males = 5, n_females = 5, line_name = "A")

  # Add line B
  pop <- add_founders(pop, n_males = 3, n_females = 7, line_name = "B")

  ind_meta <- get_table(pop, "ind_meta") %>% dplyr::collect()

  # Check total individuals
  expect_equal(nrow(ind_meta), 20)

  # Check line A IDs
  line_a <- ind_meta %>% dplyr::filter(line == "A")
  expect_equal(nrow(line_a), 10)
  expect_equal(line_a$id_ind, paste0("A-", 1:10))

  # Check line B IDs
  line_b <- ind_meta %>% dplyr::filter(line == "B")
  expect_equal(nrow(line_b), 10)
  expect_equal(line_b$id_ind, paste0("B-", 1:10))

  # Check sex assignment
  expect_equal(sum(line_a$sex == "M"), 5)
  expect_equal(sum(line_a$sex == "F"), 5)
  expect_equal(sum(line_b$sex == "M"), 3)
  expect_equal(sum(line_b$sex == "F"), 7)

  close_pop(pop)
})


test_that("sequential additions to same line continue numbering", {
  pop <- initialize_genome(
    pop_name = "test",
    n_loci = 50,
    n_chr = 2,
    chr_len_Mb = 100,
    n_haplotypes = 20,
    db_path = ":memory:"
  )

  # First batch
  pop <- add_founders(pop, n_males = 5, n_females = 5, line_name = "A")

  # Second batch
  pop <- add_founders(pop, n_males = 2, n_females = 3, line_name = "A")

  ind_meta <- get_table(pop, "ind_meta") %>% dplyr::collect()

  # Check total individuals
  expect_equal(nrow(ind_meta), 15)

  # Check IDs are sequential
  expected_ids <- paste0("A-", 1:15)
  expect_equal(ind_meta$id_ind, expected_ids)

  # Check all belong to line A
  expect_true(all(ind_meta$line == "A"))

  # Check metadata counter
  expect_equal(pop$metadata$n_individuals, 15)

  close_pop(pop)
})


test_that("add_founders works with only males", {
  pop <- initialize_genome(
    pop_name = "test",
    n_loci = 50,
    n_chr = 2,
    chr_len_Mb = 100,
    n_haplotypes = 20,
    db_path = ":memory:"
  )

  pop <- add_founders(pop, n_males = 10, n_females = 0, line_name = "A")

  ind_meta <- get_table(pop, "ind_meta") %>% dplyr::collect()

  expect_equal(nrow(ind_meta), 10)
  expect_true(all(ind_meta$sex == "M"))

  close_pop(pop)
})


test_that("add_founders works with only females", {
  pop <- initialize_genome(
    pop_name = "test",
    n_loci = 50,
    n_chr = 2,
    chr_len_Mb = 100,
    n_haplotypes = 20,
    db_path = ":memory:"
  )

  pop <- add_founders(pop, n_males = 0, n_females = 10, line_name = "A")

  ind_meta <- get_table(pop, "ind_meta") %>% dplyr::collect()

  expect_equal(nrow(ind_meta), 10)
  expect_true(all(ind_meta$sex == "F"))

  close_pop(pop)
})


test_that("add_founders errors if founder_haplotypes doesn't exist", {
  pop <- initialize_genome(
    pop_name = "test",
    n_loci = 50,
    n_chr = 2,
    chr_len_Mb = 100,
    db_path = ":memory:"
    # Note: no n_haplotypes parameter
  )

  expect_error(
    add_founders(pop, n_males = 5, n_females = 5, line_name = "A"),
    "founder_haplotypes table does not exist"
  )

  close_pop(pop)
})


test_that("add_founders errors if pop not valid tidybreed_pop", {
  expect_error(
    add_founders(list(), n_males = 5, n_females = 5, line_name = "A")
  )
})


test_that("add_founders errors if n_males + n_females = 0", {
  pop <- initialize_genome(
    pop_name = "test",
    n_loci = 50,
    n_chr = 2,
    chr_len_Mb = 100,
    n_haplotypes = 20,
    db_path = ":memory:"
  )

  expect_error(
    add_founders(pop, n_males = 0, n_females = 0, line_name = "A"),
    "At least one founder must be specified"
  )

  close_pop(pop)
})


test_that("add_founders errors if line_name has invalid format", {
  pop <- initialize_genome(
    pop_name = "test",
    n_loci = 50,
    n_chr = 2,
    chr_len_Mb = 100,
    n_haplotypes = 20,
    db_path = ":memory:"
  )

  # Starts with number
  expect_error(
    add_founders(pop, n_males = 5, n_females = 5, line_name = "1A"),
    "line_name must start with letter"
  )

  # Contains invalid characters
  expect_error(
    add_founders(pop, n_males = 5, n_females = 5, line_name = "A B"),
    "line_name must start with letter"
  )

  close_pop(pop)
})


test_that("add_founders errors if n_males or n_females invalid", {
  pop <- initialize_genome(
    pop_name = "test",
    n_loci = 50,
    n_chr = 2,
    chr_len_Mb = 100,
    n_haplotypes = 20,
    db_path = ":memory:"
  )

  # Negative values
  expect_error(
    add_founders(pop, n_males = -1, n_females = 5, line_name = "A")
  )

  expect_error(
    add_founders(pop, n_males = 5, n_females = -1, line_name = "A")
  )

  # Non-numeric
  expect_error(
    add_founders(pop, n_males = "five", n_females = 5, line_name = "A")
  )

  close_pop(pop)
})


test_that("integration: initialize_genome -> add_founders -> mutate_table", {
  # Full pipeline
  pop <- initialize_genome(
    pop_name = "test",
    n_loci = 100,
    n_chr = 5,
    chr_len_Mb = 100,
    n_haplotypes = 50,
    db_path = ":memory:"
  ) %>%
    add_founders(n_males = 10, n_females = 100, line_name = "A") %>%
    get_table("ind_meta") %>%
    mutate_table(
      gen = 0,
      farm = "FarmA",
      date_birth = as.Date("2024-01-01")
    )

  # Check ind_meta has all expected columns
  ind_meta <- get_table(pop, "ind_meta") %>% dplyr::collect()

  expect_equal(nrow(ind_meta), 110)
  expect_true(all(c("id_ind", "id_parent_1", "id_parent_2", "line", "sex", "gen", "farm", "date_birth") %in% colnames(ind_meta)))

  # Check custom metadata
  expect_true(all(ind_meta$gen == 0))
  expect_true(all(ind_meta$farm == "FarmA"))
  expect_true(all(ind_meta$date_birth == as.Date("2024-01-01")))

  close_pop(pop)
})


test_that("haplotypes are sampled with replacement", {
  # Create genome with small number of haplotypes
  pop <- initialize_genome(
    pop_name = "test",
    n_loci = 10,
    n_chr = 1,
    chr_len_Mb = 100,
    n_haplotypes = 5,  # Small pool
    db_path = ":memory:"
  )

  # Add many founders (will require sampling with replacement)
  pop <- add_founders(pop, n_males = 20, n_females = 20, line_name = "A")

  # Should not error (would error if sampling without replacement)
  # Just check it worked
  ind_meta <- get_table(pop, "ind_meta") %>% dplyr::collect()
  expect_equal(nrow(ind_meta), 40)

  close_pop(pop)
})


test_that("add_founders handles large lines efficiently", {
  skip_on_cran()

  pop <- initialize_genome(
    pop_name = "test",
    n_loci = 1000,
    n_chr = 10,
    chr_len_Mb = 100,
    n_haplotypes = 100,
    db_path = ":memory:"
  )

  # Add 1000 founders
  expect_no_error({
    pop <- add_founders(pop, n_males = 500, n_females = 500, line_name = "A")
  })

  ind_meta <- get_table(pop, "ind_meta") %>% dplyr::collect()
  expect_equal(nrow(ind_meta), 1000)

  haps <- get_table(pop, "genome_haplotype") %>% dplyr::collect()
  expect_equal(nrow(haps), 2000)  # 2 per individual

  genos <- get_table(pop, "genome_genotype") %>% dplyr::collect()
  expect_equal(nrow(genos), 1000)

  close_pop(pop)
})


# ─── ... forwarding tests ────────────────────────────────────────────────────

make_pop_for_extra <- function() {
  initialize_genome(pop_name = "t", n_loci = 20, n_chr = 1,
                    chr_len_Mb = 10, n_haplotypes = 20, db_path = ":memory:")
}

test_that("add_founders() accepts scalar ... fields and writes them to ind_meta", {
  pop <- make_pop_for_extra()
  pop <- add_founders(pop, n_males = 5, n_females = 5, line_name = "A",
                      gen = 0L, farm = "Iowa")

  result <- get_table(pop, "ind_meta") |> dplyr::collect()
  expect_equal(nrow(result), 10L)
  expect_true("gen"  %in% names(result))
  expect_true("farm" %in% names(result))
  expect_true(all(result$gen  == 0L))
  expect_true(all(result$farm == "Iowa"))

  close_pop(pop)
})

test_that("add_founders() ... uses correct DuckDB type from R type", {
  pop <- make_pop_for_extra()
  pop <- add_founders(pop, n_males = 3, n_females = 3, line_name = "A",
                      gen_int = 1L, gen_dbl = 1.0, flag = TRUE)

  col_types <- DBI::dbGetQuery(pop$db_conn,
    "SELECT column_name, data_type FROM information_schema.columns
     WHERE table_name = 'ind_meta'
     AND column_name IN ('gen_int', 'gen_dbl', 'flag')")
  types <- setNames(col_types$data_type, col_types$column_name)

  expect_equal(types[["gen_int"]], "INTEGER")
  expect_equal(types[["gen_dbl"]], "DOUBLE")
  expect_equal(types[["flag"]],    "BOOLEAN")

  close_pop(pop)
})

test_that("add_founders() ... accepts a vector of correct length", {
  pop <- make_pop_for_extra()
  scores <- c(1L, 2L, 3L, 4L, 5L)
  pop <- add_founders(pop, n_males = 5, n_females = 0, line_name = "A",
                      rank = scores)

  result <- get_table(pop, "ind_meta") |> dplyr::collect()
  expect_equal(result$rank, scores)

  close_pop(pop)
})

test_that("add_founders() ... errors for vector of wrong length", {
  pop <- make_pop_for_extra()
  expect_error(
    add_founders(pop, n_males = 5, n_females = 5, line_name = "A",
                 rank = c(1L, 2L, 3L)),
    "must equal the number of rows"
  )
  close_pop(pop)
})

test_that("add_founders() ... blocks reserved column names", {
  pop <- make_pop_for_extra()
  expect_error(
    add_founders(pop, n_males = 3, n_females = 3, line_name = "A", sex = "X"),
    "reserved"
  )
  close_pop(pop)
})

test_that("add_founders() ... preserves values across two cohort calls", {
  pop <- make_pop_for_extra()
  pop <- add_founders(pop, n_males = 3, n_females = 3, line_name = "A", gen = 0L)
  pop <- add_founders(pop, n_males = 2, n_females = 2, line_name = "B", gen = 1L)

  result <- get_table(pop, "ind_meta") |> dplyr::collect()
  expect_equal(sum(result$gen == 0L), 6L)
  expect_equal(sum(result$gen == 1L), 4L)

  close_pop(pop)
})
