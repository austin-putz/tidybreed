test_that("define_chip works with random selection", {
  pop <- initialize_genome(
    pop_chip_name = "test_chip_random",
    n_loci = 1000,
    n_chr = 10,
    chr_len_Mb = 100,
    db_path = ":memory:"
  )

  # Define chip with random selection
  pop <- define_chip(pop, chip_name = "50k", n = 500, method = "random")

  # Check column exists
  cols <- DBI::dbListFields(pop$db_conn, "genome_meta")
  expect_true("is_50k" %in% cols)

  # Check number of selected SNPs
  result <- get_table(pop, "genome_meta") %>%
    dplyr::select(is_50k) %>%
    dplyr::collect()

  expect_equal(sum(result$is_50k), 500)

  close_pop(pop)
})


test_that("define_chip works with even spacing", {
  pop <- initialize_genome(
    pop_chip_name = "test_chip_even",
    n_loci = 1000,
    n_chr = 10,
    chr_len_Mb = 100,
    db_path = ":memory:"
  )

  # Define chip with even spacing
  pop <- define_chip(pop, chip_name = "HD", n = 100, method = "even")

  # Check column exists
  cols <- DBI::dbListFields(pop$db_conn, "genome_meta")
  expect_true("is_HD" %in% cols)

  # Check number of selected SNPs
  result <- get_table(pop, "genome_meta") %>%
    dplyr::filter(is_HD == TRUE) %>%
    dplyr::collect()

  # Should be close to 100 (might be slightly less due to rounding/duplicates)
  expect_true(nrow(result) >= 95 && nrow(result) <= 105)

  close_pop(pop)
})


test_that("define_chip works with chromosome_even spacing", {
  pop <- initialize_genome(
    pop_chip_name = "test_chip_chreve",
    n_loci = 1000,
    n_chr = 10,
    chr_len_Mb = 100,
    db_path = ":memory:"
  )

  # Define chip with chromosome-even spacing
  pop <- define_chip(pop, chip_name = "10k", n = 500, method = "chromosome_even")

  # Check column exists
  cols <- DBI::dbListFields(pop$db_conn, "genome_meta")
  expect_true("is_10k" %in% cols)

  # Check total number of selected SNPs
  result <- get_table(pop, "genome_meta") %>%
    dplyr::select(chr, is_10k) %>%
    dplyr::collect()

  total_selected <- sum(result$is_10k)
  expect_true(total_selected >= 490 && total_selected <= 510)  # Allow some rounding variation

  # Check that SNPs are distributed across chromosomes
  snps_per_chr <- result %>%
    dplyr::filter(is_10k == TRUE) %>%
    dplyr::group_by(chr) %>%
    dplyr::summarise(n = dplyr::n(), .groups = "drop")

  # Should have SNPs on multiple chromosomes
  expect_true(nrow(snps_per_chr) > 5)

  close_pop(pop)
})


test_that("define_chip works with locus IDs", {
  pop <- initialize_genome(
    pop_chip_name = "test_chip_ids",
    n_loci = 1000,
    n_chr = 10,
    chr_len_Mb = 100,
    db_path = ":memory:"
  )

  # Select specific locus IDs
  selected_ids <- c(1, 10, 50, 100, 200, 500, 1000)

  # Define chip
  pop <- define_chip(pop, chip_name = "custom", locus_ids = selected_ids)

  # Check selected loci
  result <- get_table(pop, "genome_meta") %>%
    dplyr::filter(is_custom == TRUE) %>%
    dplyr::select(locus_id) %>%
    dplyr::collect()

  expect_equal(nrow(result), length(selected_ids))
  expect_true(all(result$locus_id %in% selected_ids))

  close_pop(pop)
})


test_that("define_chip works with locus names", {
  pop <- initialize_genome(
    pop_chip_name = "test_chip_names",
    n_loci = 1000,
    n_chr = 10,
    chr_len_Mb = 100,
    db_path = ":memory:"
  )

  # Select specific locus names
  selected_names <- c("Locus_1", "Locus_10", "Locus_50", "Locus_100")

  # Define chip
  pop <- define_chip(pop, chip_name = "manifest", locus_names = selected_names)

  # Check selected loci
  result <- get_table(pop, "genome_meta") %>%
    dplyr::filter(is_manifest == TRUE) %>%
    dplyr::select(locus_name) %>%
    dplyr::collect()

  expect_equal(nrow(result), length(selected_names))
  expect_true(all(result$locus_name %in% selected_names))

  close_pop(pop)
})


test_that("define_chip works with locus_tf logical vector", {
  pop <- initialize_genome(
    pop_chip_name = "test_chip_tf",
    n_loci = 1000,
    n_chr = 10,
    chr_len_Mb = 100,
    db_path = ":memory:"
  )

  # Define a 50k chip first
  pop <- define_chip(pop, chip_name = "50k", n = 500, method = "random")

  # Pull the logical vector for the 50k chip
  chip_tf <- get_table(pop, "genome_meta") %>%
    dplyr::pull(is_50k)

  # Define complement chip using !chip_tf
  pop <- define_chip(pop, chip_name = "non50k", locus_tf = !chip_tf)

  result <- get_table(pop, "genome_meta") %>%
    dplyr::select(is_50k, is_non50k) %>%
    dplyr::collect()

  # Complement: every locus is on exactly one chip
  expect_true(all(result$is_50k != result$is_non50k))
  expect_equal(sum(result$is_non50k), 500)

  close_pop(pop)
})


test_that("define_chip errors with non-logical locus_tf", {
  pop <- initialize_genome(
    pop_chip_name = "test_chip_tf_bad_type",
    n_loci = 100,
    n_chr = 5,
    chr_len_Mb = 100,
    db_path = ":memory:"
  )

  expect_error(
    define_chip(pop, chip_name = "bad", locus_tf = 1:100),
    "locus_tf must be a logical vector"
  )

  close_pop(pop)
})


test_that("define_chip errors with wrong-length locus_tf", {
  pop <- initialize_genome(
    pop_chip_name = "test_chip_tf_bad_len",
    n_loci = 100,
    n_chr = 5,
    chr_len_Mb = 100,
    db_path = ":memory:"
  )

  expect_error(
    define_chip(pop, chip_name = "bad", locus_tf = rep(TRUE, 50)),
    "locus_tf length"
  )

  close_pop(pop)
})


test_that("define_chip respects custom column name", {
  pop <- initialize_genome(
    pop_chip_name = "test_chip_colname",
    n_loci = 1000,
    n_chr = 10,
    chr_len_Mb = 100,
    db_path = ":memory:"
  )

  # Define chip with custom column name
  pop <- define_chip(
    pop,
    chip_name = "bovine_50k",
    n = 500,
    col_chip_name = "SNP_50k"
  )

  # Check custom column name exists
  cols <- DBI::dbListFields(pop$db_conn, "genome_meta")
  expect_true("SNP_50k" %in% cols)
  expect_false("is_bovine_50k" %in% cols)  # Default name should NOT exist

  close_pop(pop)
})


test_that("define_chip can define multiple chips", {
  pop <- initialize_genome(
    pop_chip_name = "test_chip_multi",
    n_loci = 1000,
    n_chr = 10,
    chr_len_Mb = 100,
    db_path = ":memory:"
  )

  # Define multiple chips
  pop <- pop %>%
    define_chip(chip_name = "50k", n = 500, method = "random") %>%
    define_chip(chip_name = "HD", n = 900, method = "even")

  # Check both columns exist
  cols <- DBI::dbListFields(pop$db_conn, "genome_meta")
  expect_true("is_50k" %in% cols)
  expect_true("is_HD" %in% cols)

  # Check number of SNPs
  result <- get_table(pop, "genome_meta") %>%
    dplyr::select(is_50k, is_HD) %>%
    dplyr::collect()

  expect_equal(sum(result$is_50k), 500)
  expect_true(sum(result$is_HD) >= 890 && sum(result$is_HD) <= 910)

  close_pop(pop)
})


test_that("define_chip errors with no selection method", {
  pop <- initialize_genome(
    pop_chip_name = "test_chip_nomethod",
    n_loci = 1000,
    n_chr = 10,
    chr_len_Mb = 100,
    db_path = ":memory:"
  )

  # No selection method specified
  expect_error(
    define_chip(pop, chip_name = "50k"),
    "Must specify one selection method: n, locus_tf, locus_ids, or locus_names"
  )

  close_pop(pop)
})


test_that("define_chip errors with multiple selection methods", {
  pop <- initialize_genome(
    pop_chip_name = "test_chip_multimethod",
    n_loci = 1000,
    n_chr = 10,
    chr_len_Mb = 100,
    db_path = ":memory:"
  )

  # Multiple selection methods
  expect_error(
    define_chip(pop, chip_name = "50k", n = 500, locus_ids = c(1, 2, 3)),
    "Cannot specify multiple selection methods"
  )

  close_pop(pop)
})


test_that("define_chip errors with invalid method", {
  pop <- initialize_genome(
    pop_chip_name = "test_chip_badmethod",
    n_loci = 1000,
    n_chr = 10,
    chr_len_Mb = 100,
    db_path = ":memory:"
  )

  # Invalid method
  expect_error(
    define_chip(pop, chip_name = "50k", n = 500, method = "bad_method"),
    "Invalid method"
  )

  close_pop(pop)
})


test_that("define_chip errors when n exceeds total loci", {
  pop <- initialize_genome(
    pop_chip_name = "test_chip_toomany",
    n_loci = 100,
    n_chr = 5,
    chr_len_Mb = 100,
    db_path = ":memory:"
  )

  # n exceeds total loci
  expect_error(
    define_chip(pop, chip_name = "50k", n = 500),
    "exceeds total loci"
  )

  close_pop(pop)
})


test_that("define_chip errors with invalid locus IDs", {
  pop <- initialize_genome(
    pop_chip_name = "test_chip_badids",
    n_loci = 1000,
    n_chr = 10,
    chr_len_Mb = 100,
    db_path = ":memory:"
  )

  # Invalid locus IDs
  expect_error(
    define_chip(pop, chip_name = "custom", locus_ids = c(1, 2, 9999)),
    "Invalid locus_ids"
  )

  close_pop(pop)
})


test_that("define_chip errors with invalid locus names", {
  pop <- initialize_genome(
    pop_chip_name = "test_chip_badnames",
    n_loci = 1000,
    n_chr = 10,
    chr_len_Mb = 100,
    db_path = ":memory:"
  )

  # Invalid locus names
  expect_error(
    define_chip(pop, chip_name = "custom", locus_names = c("Locus_1", "BadName")),
    "Invalid locus_names"
  )

  close_pop(pop)
})


test_that("define_chip integration: chip + effects", {
  pop <- initialize_genome(
    pop_chip_name = "test_chip_integration",
    n_loci = 1000,
    n_chr = 10,
    chr_len_Mb = 100,
    db_path = ":memory:"
  )

  # Define chip
  pop <- define_chip(pop, chip_name = "50k", n = 500, method = "random")

  # Get genome metadata
  genome <- get_table(pop, "genome_meta") %>% dplyr::collect()
  n_loci <- nrow(genome)

  # Create effect vector (effects only for chip SNPs)
  set.seed(123)
  effect_vec <- ifelse(
    genome$is_50k,
    rnorm(n_loci, 0, 1),
    0
  )

  # Add effects
  pop <- mutate_genome_meta(pop, effect_50k = effect_vec)

  # Verify
  result <- get_table(pop, "genome_meta") %>%
    dplyr::select(is_50k, effect_50k) %>%
    dplyr::collect()

  # All non-chip SNPs should have 0 effect
  non_chip_effects <- result$effect_50k[!result$is_50k]
  expect_true(all(non_chip_effects == 0))

  # Chip SNPs should have non-zero effects (mostly)
  chip_effects <- result$effect_50k[result$is_50k]
  expect_true(sum(chip_effects != 0) > 450)  # Most should be non-zero

  close_pop(pop)
})
