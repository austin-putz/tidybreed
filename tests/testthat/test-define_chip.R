make_genome_pop <- function(pop_name = "chip_test", n_loci = 1000, n_chr = 10) {
  initialize_genome(
    pop_name     = pop_name,
    n_loci       = n_loci,
    n_chr        = n_chr,
    chr_len_Mb   = 100,
    n_haplotypes = 20,
    db_path      = ":memory:"
  )
}


test_that("define_chip() creates is_{chip_name} column with correct TRUE count", {
  pop <- make_genome_pop("chip_basic")
  set.seed(1)
  selected <- pop |> get_table("genome_meta") |> dplyr::collect() |>
    dplyr::slice_sample(n = 500) |> dplyr::pull(locus_name)

  pop <- pop |>
    get_table("genome_meta") |>
    dplyr::filter(locus_name %in% selected) |>
    define_chip("50K")

  cols <- DBI::dbListFields(pop$db_conn, "genome_meta")
  expect_true("is_50K" %in% cols)

  result <- pop |> get_table("genome_meta") |> dplyr::collect()
  expect_equal(sum(result$is_50K), 500)
  expect_equal(sum(!result$is_50K), 500)

  close_pop(pop)
})


test_that("define_chip() correctly marks non-selected loci FALSE", {
  pop <- make_genome_pop("chip_false_check", n_loci = 100, n_chr = 5)

  pop <- pop |>
    get_table("genome_meta") |>
    dplyr::filter(chr == 1L) |>
    define_chip("chr1")

  result <- pop |> get_table("genome_meta") |> dplyr::collect()
  chr1_loci    <- result$is_chr1[result$chr == 1L]
  nonchr1_loci <- result$is_chr1[result$chr != 1L]

  expect_true(all(chr1_loci))
  expect_true(all(!nonchr1_loci))

  close_pop(pop)
})


test_that("define_chip() supports custom col_name", {
  pop <- make_genome_pop("chip_colname")
  set.seed(42)
  selected <- pop |> get_table("genome_meta") |> dplyr::collect() |>
    dplyr::slice_sample(n = 200) |> dplyr::pull(locus_name)

  pop <- pop |>
    get_table("genome_meta") |>
    dplyr::filter(locus_name %in% selected) |>
    define_chip("bovine_50k", col_name = "SNP_50k")

  cols <- DBI::dbListFields(pop$db_conn, "genome_meta")
  expect_true("SNP_50k" %in% cols)
  expect_false("is_bovine_50k" %in% cols)

  close_pop(pop)
})


test_that("define_chip() can define multiple chips on the same pop", {
  pop <- make_genome_pop("chip_multi")
  set.seed(7)
  sel_a <- pop |> get_table("genome_meta") |> dplyr::collect() |>
    dplyr::slice_sample(n = 500) |> dplyr::pull(locus_name)
  set.seed(8)
  sel_b <- pop |> get_table("genome_meta") |> dplyr::collect() |>
    dplyr::slice_sample(n = 900) |> dplyr::pull(locus_name)

  pop <- pop |>
    get_table("genome_meta") |>
    dplyr::filter(locus_name %in% sel_a) |>
    define_chip("50K")

  pop <- pop |>
    get_table("genome_meta") |>
    dplyr::filter(locus_name %in% sel_b) |>
    define_chip("HD")

  cols <- DBI::dbListFields(pop$db_conn, "genome_meta")
  expect_true("is_50K" %in% cols)
  expect_true("is_HD"  %in% cols)

  res <- pop |> get_table("genome_meta") |> dplyr::collect()
  expect_equal(sum(res$is_50K), 500)
  expect_equal(sum(res$is_HD),  900)

  close_pop(pop)
})


test_that("define_chip() complement via filter on existing chip column", {
  pop <- make_genome_pop("chip_complement")
  set.seed(2)
  sel <- pop |> get_table("genome_meta") |> dplyr::collect() |>
    dplyr::slice_sample(n = 500) |> dplyr::pull(locus_name)

  pop <- pop |>
    get_table("genome_meta") |>
    dplyr::filter(locus_name %in% sel) |>
    define_chip("50K")

  pop <- pop |>
    get_table("genome_meta") |>
    dplyr::filter(is_50K == FALSE) |>
    define_chip("non50K")

  res <- pop |> get_table("genome_meta") |> dplyr::collect()
  expect_true(all(res$is_50K != res$is_non50K))
  expect_equal(sum(res$is_non50K), 500)

  close_pop(pop)
})


test_that("define_chip() errors on empty filter result", {
  pop <- make_genome_pop("chip_empty_filter", n_loci = 100, n_chr = 5)

  expect_error(
    pop |>
      get_table("genome_meta") |>
      dplyr::filter(chr == 999L) |>
      define_chip("bad"),
    "No loci selected"
  )

  close_pop(pop)
})


test_that("define_chip() errors when column already exists", {
  pop <- make_genome_pop("chip_dup", n_loci = 100, n_chr = 5)
  set.seed(1)
  sel <- pop |> get_table("genome_meta") |> dplyr::collect() |>
    dplyr::slice_sample(n = 50) |> dplyr::pull(locus_name)

  pop <- pop |>
    get_table("genome_meta") |>
    dplyr::filter(locus_name %in% sel) |>
    define_chip("chip_a")

  expect_error(
    pop |>
      get_table("genome_meta") |>
      dplyr::filter(locus_name %in% sel) |>
      define_chip("chip_a"),
    "already exists"
  )

  close_pop(pop)
})


test_that("define_chip() errors when tbl has no locus_id column", {
  pop <- make_genome_pop("chip_wrong_tbl", n_loci = 100, n_chr = 5)
  pop <- add_founders(pop, n_males = 5, n_females = 5, line_name = "A")

  expect_error(
    pop |>
      get_table("ind_meta") |>
      define_chip("chip_bad"),
    "locus_id"
  )

  close_pop(pop)
})


test_that("define_chip() errors when first arg is not tidybreed_table", {
  pop <- make_genome_pop("chip_wrong_arg", n_loci = 100, n_chr = 5)

  expect_error(
    define_chip(pop, "bad_chip"),
    "inherits"
  )

  close_pop(pop)
})


test_that("define_chip() integration: chip column used downstream by add_genotypes", {
  pop <- initialize_genome(
    pop_name         = "chip_integration",
    n_loci           = 200,
    n_chr            = 4,
    chr_len_Mb       = 100,
    n_haplotypes     = 50,
    db_path          = ":memory:",
    fixed_allele_freq = 0.5
  )
  pop <- add_founders(pop, n_males = 10, n_females = 10, line_name = "A")

  set.seed(99)
  sel <- pop |> get_table("genome_meta") |> dplyr::collect() |>
    dplyr::slice_sample(n = 100) |> dplyr::pull(locus_name)

  pop <- pop |>
    get_table("genome_meta") |>
    dplyr::filter(locus_name %in% sel) |>
    define_chip("test_chip")

  pop <- pop |> get_table("ind_meta") |> add_genotypes("test_chip")

  geno <- pop |> get_table("ind_meta") |> extract_genotypes("test_chip")
  expect_equal(nrow(geno), 20)
  expect_equal(ncol(geno) - 1L, 100)  # id_ind + 100 locus columns

  close_pop(pop)
})
