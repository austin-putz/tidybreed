make_extract_pop <- function(pop_name = "extract_test", chip_name = "50k",
                             n_males = 5, n_females = 5) {
  pop <- open_pop(pop_name = pop_name, db_name = ":memory:") |>
    define_genome(n_loci = 100, n_chr = 2, chr_len_Mb = 100) |>
    define_founder_haplotypes(n_haplotypes = 20) |>
    get_table("founder_haplotypes") |>
    add_founders(n_males = n_males, n_females = n_females, line_name = "A")

  # Use chr == 1 loci as chip; with n_loci = 100, n_chr = 2 this is 50 loci
  pop <- pop |>
    get_table("genome_meta") |>
    dplyr::filter(chr == 1L) |>
    define_chip(chip_name)

  pop
}

make_qtl_pop <- function(pop_name = "qtl_test", n_males = 5, n_females = 5) {
  pop <- open_pop(pop_name = pop_name, db_name = ":memory:") |>
    define_genome(n_loci = 100, n_chr = 2, chr_len_Mb = 100) |>
    define_founder_haplotypes(n_haplotypes = 20) |>
    get_table("founder_haplotypes") |>
    add_founders(n_males = n_males, n_females = n_females, line_name = "A") |>
    define_trait("ADG", target_add_mean = 100, target_add_var = 400)

  # chr 1 loci are QTL (50 loci)
  pop <- pop |>
    get_table("genome_meta") |>
    dplyr::filter(chr == 1L) |>
    define_additive_effects("ADG")

  pop
}


# ---------------------------------------------------------------------------
# Existing chip tests (unchanged)
# ---------------------------------------------------------------------------

test_that("extract_genotypes() returns a tibble with correct dimensions", {
  pop <- make_extract_pop("test_eg_dim")
  n_chip <- sum(dplyr::collect(get_table(pop, "genome_meta"))$is_50k)

  pop <- pop |> get_table("ind_meta") |> add_genotypes("50k")
  geno <- pop |> get_table("ind_meta") |> extract_genotypes("50k")

  expect_s3_class(geno, "tbl_df")
  expect_equal(nrow(geno), 10)
  expect_equal(ncol(geno), n_chip + 1L)  # id_ind + chip loci
  expect_true("id_ind" %in% names(geno))

  close_pop(pop)
})

test_that("extract_genotypes() respects filter()", {
  pop <- make_extract_pop("test_eg_filter")
  pop <- pop |> get_table("ind_meta") |> add_genotypes("50k")

  geno <- pop |>
    get_table("ind_meta") |>
    dplyr::filter(sex == "F") |>
    extract_genotypes("50k")

  expect_equal(nrow(geno), 5)
  close_pop(pop)
})

test_that("extract_genotypes() intersects filter with has_<chip> == TRUE", {
  pop <- make_extract_pop("test_eg_intersect")

  # Only genotype females
  pop <- pop |>
    get_table("ind_meta") |>
    dplyr::filter(sex == "F") |>
    add_genotypes("50k")

  # Extract all — should only return females (males not genotyped)
  geno <- pop |> get_table("ind_meta") |> extract_genotypes("50k")
  expect_equal(nrow(geno), 5)

  close_pop(pop)
})

test_that("extract_genotypes() returns only chip loci columns with correct encoding", {
  pop <- open_pop(pop_name = "test_eg_cols", db_name = ":memory:") |>
    define_genome(n_loci = 100, n_chr = 2, chr_len_Mb = 100) |>
    define_founder_haplotypes(n_haplotypes = 20) |>
    get_table("founder_haplotypes") |>
    add_founders(n_males = 3, n_females = 3, line_name = "A")

  pop <- pop |>
    get_table("genome_meta") |>
    dplyr::filter(chr == 2L) |>
    define_chip("HD")

  n_chip <- sum(dplyr::collect(get_table(pop, "genome_meta"))$is_HD)

  pop <- pop |> get_table("ind_meta") |> add_genotypes("HD")
  geno <- pop |> get_table("ind_meta") |> extract_genotypes("HD")

  locus_cols <- grep("^locus_", names(geno), value = TRUE)
  expect_equal(length(locus_cols), n_chip)

  # All genotype values must be 0, 1, or 2
  geno_vals <- unlist(geno[, locus_cols])
  expect_true(all(geno_vals %in% c(0L, 1L, 2L)))

  close_pop(pop)
})

test_that("extract_genotypes() errors if add_genotypes() not called", {
  pop <- make_extract_pop("test_eg_nogenotype")

  expect_error(
    pop |> get_table("ind_meta") |> extract_genotypes("50k"),
    "not found in ind_meta"
  )
  close_pop(pop)
})

test_that("extract_genotypes() errors if chip not defined", {
  pop <- open_pop(pop_name = "test_eg_nochip", db_name = ":memory:") |>
    define_genome(n_loci = 100, n_chr = 2, chr_len_Mb = 100) |>
    define_founder_haplotypes(n_haplotypes = 20) |>
    get_table("founder_haplotypes") |>
    add_founders(n_males = 5, n_females = 5, line_name = "A")

  expect_error(
    pop |> get_table("ind_meta") |> extract_genotypes("nonexistent"),
    "not found in genome_meta"
  )
  close_pop(pop)
})


# ---------------------------------------------------------------------------
# New effects_tbl tests
# ---------------------------------------------------------------------------

test_that("extract_genotypes() with effects_tbl returns QTL loci for all individuals", {
  pop <- make_qtl_pop("test_eg_qtl_basic")

  geno <- pop |>
    get_table("ind_meta") |>
    extract_genotypes(
      effects_tbl = get_table(pop, "genome_effects") |>
        dplyr::filter(trait_name == "ADG")
    )

  expect_s3_class(geno, "tbl_df")
  expect_equal(nrow(geno), 10)
  expect_true("id_ind" %in% names(geno))

  locus_cols <- grep("^locus_", names(geno), value = TRUE)
  expect_equal(length(locus_cols), 50L)  # chr 1 has 50 QTL loci

  close_pop(pop)
})

test_that("extract_genotypes() effects_tbl respects pending tbl filter", {
  pop <- make_qtl_pop("test_eg_qtl_filter")

  geno <- pop |>
    get_table("ind_meta") |>
    dplyr::filter(sex == "F") |>
    extract_genotypes(
      effects_tbl = get_table(pop, "genome_effects") |>
        dplyr::filter(trait_name == "ADG")
    )

  expect_equal(nrow(geno), 5)

  close_pop(pop)
})

test_that("extract_genotypes() chip + effects_tbl unions loci (no overlap)", {
  pop <- make_qtl_pop("test_eg_union_nolap")

  # chip on chr 2 (50 loci), QTL on chr 1 (50 loci) → 100 total
  pop <- pop |>
    get_table("genome_meta") |>
    dplyr::filter(chr == 2L) |>
    define_chip("HD")
  pop <- pop |> get_table("ind_meta") |> add_genotypes("HD")

  geno <- pop |>
    get_table("ind_meta") |>
    extract_genotypes(
      chip_name   = "HD",
      effects_tbl = get_table(pop, "genome_effects") |>
        dplyr::filter(trait_name == "ADG")
    )

  locus_cols <- grep("^locus_", names(geno), value = TRUE)
  expect_equal(length(locus_cols), 100L)

  close_pop(pop)
})

test_that("extract_genotypes() chip + effects_tbl deduplicates overlapping loci", {
  pop <- make_qtl_pop("test_eg_union_overlap")

  # chip on chr 1 (50 loci) — same loci as the QTL → union should still be 50
  pop <- pop |>
    get_table("genome_meta") |>
    dplyr::filter(chr == 1L) |>
    define_chip("50k")
  pop <- pop |> get_table("ind_meta") |> add_genotypes("50k")

  geno <- pop |>
    get_table("ind_meta") |>
    extract_genotypes(
      chip_name   = "50k",
      effects_tbl = get_table(pop, "genome_effects") |>
        dplyr::filter(trait_name == "ADG")
    )

  locus_cols <- grep("^locus_", names(geno), value = TRUE)
  expect_equal(length(locus_cols), 50L)

  close_pop(pop)
})

test_that("extract_genotypes() errors when effects_tbl filter yields no loci", {
  pop <- make_qtl_pop("test_eg_qtl_empty")

  expect_error(
    pop |>
      get_table("ind_meta") |>
      extract_genotypes(
        effects_tbl = get_table(pop, "genome_effects") |>
          dplyr::filter(trait_name == "NONEXISTENT")
      ),
    "No QTL loci found"
  )

  close_pop(pop)
})

test_that("extract_genotypes() errors when effects_tbl is wrong table", {
  pop <- make_qtl_pop("test_eg_qtl_wrongtbl")

  expect_error(
    pop |>
      get_table("ind_meta") |>
      extract_genotypes(
        effects_tbl = get_table(pop, "ind_meta")
      ),
    "genome_effects"
  )

  close_pop(pop)
})

test_that("extract_genotypes() errors when both chip_name and effects_tbl are NULL", {
  pop <- make_qtl_pop("test_eg_both_null")

  expect_error(
    pop |> get_table("ind_meta") |> extract_genotypes(),
    "At least one of"
  )

  close_pop(pop)
})
