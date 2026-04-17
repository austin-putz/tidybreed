# Helper: minimal population with founders
make_test_pop <- function(n_loci = 100, n_chr = 2, n_males = 5, n_females = 5) {
  pop <- initialize_genome(
    pop_name    = "test",
    n_loci      = n_loci,
    n_chr       = n_chr,
    chr_len_Mb  = 100,
    n_haplotypes = 50,
    db_path     = ":memory:"
  )
  add_founders(pop, n_males = n_males, n_females = n_females, line_name = "A")
}


test_that("add_offspring creates rows in all three tables", {
  pop <- make_test_pop()

  matings <- tibble::tibble(
    id_parent_1 = "A-1",
    id_parent_2 = "A-6",
    sex         = "M",
    line        = "A"
  )

  pop <- add_offspring(pop, matings)

  ind_meta <- get_table(pop, "ind_meta") |> dplyr::collect()
  haps     <- get_table(pop, "genome_haplotype") |> dplyr::collect()
  genos    <- get_table(pop, "genome_genotype")  |> dplyr::collect()

  # 10 founders + 1 offspring
  expect_equal(nrow(ind_meta), 11)
  # 10 founders × 2 + 1 offspring × 2
  expect_equal(nrow(haps), 22)
  expect_equal(nrow(genos), 11)

  close_pop(pop)
})


test_that("add_offspring writes correct parent IDs and sex to ind_meta", {
  pop <- make_test_pop()

  matings <- tibble::tibble(
    id_parent_1 = "A-1",
    id_parent_2 = "A-6",
    sex         = "F",
    line        = "A"
  )

  pop <- add_offspring(pop, matings)

  offspring <- get_table(pop, "ind_meta") |>
    dplyr::filter(id_ind == "A-11") |>
    dplyr::collect()

  expect_equal(nrow(offspring), 1)
  expect_equal(offspring$id_parent_1, "A-1")
  expect_equal(offspring$id_parent_2, "A-6")
  expect_equal(offspring$sex,         "F")
  expect_equal(offspring$line,        "A")

  close_pop(pop)
})


test_that("add_offspring assigns sequential IDs continuing from founders", {
  pop <- make_test_pop(n_males = 3, n_females = 3)

  matings <- tibble::tibble(
    id_parent_1 = c("A-1", "A-1", "A-2"),
    id_parent_2 = c("A-4", "A-5", "A-6"),
    sex         = c("M",   "F",   "M"),
    line        = "A"
  )

  pop <- add_offspring(pop, matings)

  offspring_ids <- get_table(pop, "ind_meta") |>
    dplyr::filter(!is.na(id_parent_1)) |>
    dplyr::pull(id_ind) |>
    sort()

  expect_equal(offspring_ids, c("A-7", "A-8", "A-9"))

  close_pop(pop)
})


test_that("add_offspring: multiple offspring per dam (same dam repeated)", {
  pop <- make_test_pop()

  matings <- tibble::tibble(
    id_parent_1 = rep("A-1", 5),
    id_parent_2 = rep("A-6", 5),
    sex         = rep("F",   5),
    line        = "A"
  )

  pop <- add_offspring(pop, matings)

  ind_meta <- get_table(pop, "ind_meta") |> dplyr::collect()
  expect_equal(nrow(ind_meta), 15)  # 10 founders + 5 offspring

  n_from_dam6 <- get_table(pop, "ind_meta") |>
    dplyr::filter(id_parent_2 == "A-6") |>
    dplyr::count() |>
    dplyr::pull(n)
  expect_equal(n_from_dam6, 5)

  close_pop(pop)
})


test_that("add_offspring: multiple sires per dam (pooled semen pattern)", {
  pop <- make_test_pop()

  matings <- tibble::tibble(
    id_parent_1 = c("A-1", "A-2", "A-3"),
    id_parent_2 = rep("A-6", 3),
    sex         = c("M", "F", "M"),
    line        = "A"
  )

  pop <- add_offspring(pop, matings)

  offspring <- get_table(pop, "ind_meta") |>
    dplyr::filter(!is.na(id_parent_1)) |>
    dplyr::collect()

  expect_equal(nrow(offspring), 3)
  expect_equal(sort(offspring$id_parent_1), c("A-1", "A-2", "A-3"))
  expect_true(all(offspring$id_parent_2 == "A-6"))

  close_pop(pop)
})


test_that("add_offspring: cross-line offspring with independent ID sequences", {
  pop <- initialize_genome(
    pop_name    = "test",
    n_loci      = 100,
    n_chr       = 2,
    chr_len_Mb  = 100,
    n_haplotypes = 50,
    db_path     = ":memory:"
  )
  pop <- add_founders(pop, n_males = 2, n_females = 2, line_name = "A")
  pop <- add_founders(pop, n_males = 2, n_females = 2, line_name = "B")

  # Cross A sire × B dam, assign to line "C"
  matings <- tibble::tibble(
    id_parent_1 = c("A-1", "A-2"),
    id_parent_2 = c("B-3", "B-4"),
    sex         = c("M",   "F"),
    line        = "C"
  )

  pop <- add_offspring(pop, matings)

  offspring_ids <- get_table(pop, "ind_meta") |>
    dplyr::filter(line == "C") |>
    dplyr::pull(id_ind) |>
    sort()

  expect_equal(offspring_ids, c("C-1", "C-2"))

  close_pop(pop)
})


test_that("add_offspring writes extra column (gen) to ind_meta", {
  pop <- make_test_pop()

  matings <- tibble::tibble(
    id_parent_1 = "A-1",
    id_parent_2 = "A-6",
    sex         = "M",
    line        = "A",
    gen         = 2L
  )

  pop <- add_offspring(pop, matings)

  ind_meta <- get_table(pop, "ind_meta") |> dplyr::collect()

  expect_true("gen" %in% colnames(ind_meta))

  offspring_gen <- dplyr::filter(ind_meta, id_ind == "A-11")$gen
  expect_equal(offspring_gen, 2L)

  close_pop(pop)
})


test_that("add_offspring adds a new extra column not previously in ind_meta", {
  pop <- make_test_pop()

  # ind_meta currently has no 'farm' column
  existing_cols <- DBI::dbListFields(pop$db_conn, "ind_meta")
  expect_false("farm" %in% existing_cols)

  matings <- tibble::tibble(
    id_parent_1 = c("A-1", "A-2"),
    id_parent_2 = c("A-6", "A-7"),
    sex         = c("M",   "F"),
    line        = "A",
    farm        = c("Iowa", "Iowa")
  )

  pop <- add_offspring(pop, matings)

  updated_cols <- DBI::dbListFields(pop$db_conn, "ind_meta")
  expect_true("farm" %in% updated_cols)

  farm_vals <- get_table(pop, "ind_meta") |>
    dplyr::filter(!is.na(id_parent_1)) |>
    dplyr::pull(farm) |>
    sort()
  expect_equal(farm_vals, c("Iowa", "Iowa"))

  close_pop(pop)
})


test_that("add_offspring accepts id_sire / id_dam aliases", {
  pop <- make_test_pop()

  matings_alias <- tibble::tibble(
    id_sire = "A-1",
    id_dam  = "A-6",
    sex     = "M",
    line    = "A"
  )

  pop <- add_offspring(pop, matings_alias)

  offspring <- get_table(pop, "ind_meta") |>
    dplyr::filter(!is.na(id_parent_1)) |>
    dplyr::collect()

  expect_equal(nrow(offspring), 1)
  expect_equal(offspring$id_parent_1, "A-1")
  expect_equal(offspring$id_parent_2, "A-6")

  close_pop(pop)
})


test_that("add_offspring: genome_haplotype has 2 rows per offspring with correct structure", {
  pop <- make_test_pop(n_loci = 50)

  matings <- tibble::tibble(
    id_parent_1 = c("A-1", "A-2"),
    id_parent_2 = c("A-6", "A-7"),
    sex         = c("M",   "F"),
    line        = "A"
  )

  pop <- add_offspring(pop, matings)

  offspring_ids <- c("A-11", "A-12")

  haps <- get_table(pop, "genome_haplotype") |>
    dplyr::filter(id_ind %in% offspring_ids) |>
    dplyr::collect()

  expect_equal(nrow(haps), 4)  # 2 offspring × 2 haplotypes
  expect_true(all(haps$parent_origin %in% c(1L, 2L)))

  for (oid in offspring_ids) {
    ind_haps <- dplyr::filter(haps, id_ind == oid)
    expect_equal(nrow(ind_haps), 2)
    expect_equal(sort(ind_haps$parent_origin), c(1L, 2L))
  }

  locus_cols <- paste0("locus_", 1:50)
  hap_vals   <- unlist(haps[, locus_cols])
  expect_true(all(hap_vals %in% c(0L, 1L)))

  close_pop(pop)
})


test_that("add_offspring: genotype equals sum of the two offspring haplotypes", {
  pop <- make_test_pop(n_loci = 50)

  matings <- tibble::tibble(
    id_parent_1 = "A-1",
    id_parent_2 = "A-6",
    sex         = "M",
    line        = "A"
  )

  pop <- add_offspring(pop, matings)

  offspring_id <- "A-11"
  locus_cols   <- paste0("locus_", 1:50)

  haps <- get_table(pop, "genome_haplotype") |>
    dplyr::filter(id_ind == offspring_id) |>
    dplyr::arrange(parent_origin) |>
    dplyr::collect()

  geno <- get_table(pop, "genome_genotype") |>
    dplyr::filter(id_ind == offspring_id) |>
    dplyr::collect()

  hap1_vals  <- as.integer(haps[1, locus_cols])
  hap2_vals  <- as.integer(haps[2, locus_cols])
  geno_vals  <- as.integer(geno[1, locus_cols])
  expected   <- hap1_vals + hap2_vals

  expect_equal(geno_vals, expected)

  close_pop(pop)
})


test_that("add_offspring recombination produces variation across offspring", {
  set.seed(42)
  pop <- initialize_genome(
    pop_name    = "test",
    n_loci      = 500,
    n_chr       = 5,
    chr_len_Mb  = 100,
    n_haplotypes = 100,
    db_path     = ":memory:"
  )
  pop <- add_founders(pop, n_males = 1, n_females = 1, line_name = "A")

  # Produce 50 offspring from the same pair
  matings <- tibble::tibble(
    id_parent_1 = rep("A-1", 50),
    id_parent_2 = rep("A-2", 50),
    sex         = rep("M",   50),
    line        = "A"
  )

  pop <- add_offspring(pop, matings)

  locus_cols <- paste0("locus_", 1:500)

  offspring_ids <- get_table(pop, "ind_meta") |>
    dplyr::filter(!is.na(id_parent_1)) |>
    dplyr::pull(id_ind)

  haps <- get_table(pop, "genome_haplotype") |>
    dplyr::filter(id_ind %in% offspring_ids, parent_origin == 1L) |>
    dplyr::collect()

  hap_mat <- as.matrix(haps[, locus_cols])

  # Not all 50 gametes should be identical
  unique_gametes <- nrow(unique(hap_mat))
  expect_gt(unique_gametes, 1)

  close_pop(pop)
})


# ============================================================================
# Error handling tests
# ============================================================================

test_that("add_offspring errors on non-data.frame matings", {
  pop <- make_test_pop()
  expect_error(add_offspring(pop, list(id_parent_1 = "A-1")), "data.frame")
  close_pop(pop)
})


test_that("add_offspring errors on empty matings", {
  pop  <- make_test_pop()
  empty <- tibble::tibble(id_parent_1 = character(), id_parent_2 = character(),
                          sex = character(), line = character())
  expect_error(add_offspring(pop, empty), "at least 1 row")
  close_pop(pop)
})


test_that("add_offspring errors on missing required columns", {
  pop <- make_test_pop()

  no_sex <- tibble::tibble(id_parent_1 = "A-1", id_parent_2 = "A-6", line = "A")
  expect_error(add_offspring(pop, no_sex), "sex")

  no_line <- tibble::tibble(id_parent_1 = "A-1", id_parent_2 = "A-6", sex = "M")
  expect_error(add_offspring(pop, no_line), "line")

  close_pop(pop)
})


test_that("add_offspring errors on invalid sex values", {
  pop <- make_test_pop()

  bad <- tibble::tibble(id_parent_1 = "A-1", id_parent_2 = "A-6",
                        sex = "X", line = "A")
  expect_error(add_offspring(pop, bad), "sex")

  close_pop(pop)
})


test_that("add_offspring errors on invalid line names", {
  pop <- make_test_pop()

  bad <- tibble::tibble(id_parent_1 = "A-1", id_parent_2 = "A-6",
                        sex = "M", line = "1badline")
  expect_error(add_offspring(pop, bad), "line")

  close_pop(pop)
})


test_that("add_offspring errors on non-existent parent IDs", {
  pop <- make_test_pop()

  bad <- tibble::tibble(id_parent_1 = "Z-999", id_parent_2 = "A-6",
                        sex = "M", line = "A")
  expect_error(add_offspring(pop, bad), "Z-999")

  close_pop(pop)
})


test_that("add_offspring errors on reserved extra column name (id_ind)", {
  pop <- make_test_pop()

  bad <- tibble::tibble(id_parent_1 = "A-1", id_parent_2 = "A-6",
                        sex = "M", line = "A", id_ind = "fake")
  expect_error(add_offspring(pop, bad), "id_ind")

  close_pop(pop)
})
