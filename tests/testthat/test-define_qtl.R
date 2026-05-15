make_qtl_pop <- function(pop_name = "qtl", n_loci = 500, n_chr = 5) {
  pop <- initialize_genome(
    pop_name          = pop_name,
    n_loci            = n_loci,
    n_chr             = n_chr,
    chr_len_Mb        = 100,
    n_haplotypes      = 100,
    db_path           = ":memory:",
    fixed_allele_freq = 0.5
  )
  pop <- add_founders(pop, n_males = 20, n_females = 20, line_name = "A")
  pop <- add_trait(pop, "ADG", target_add_var = 0.25, residual_var = 0.75)
  pop
}


test_that("define_qtl() writes is_QTL_{trait} with correct TRUE count", {
  pop <- make_qtl_pop("qtl_basic")
  set.seed(1)
  sel <- pop |> get_table("genome_meta") |> dplyr::collect() |>
    dplyr::slice_sample(n = 50) |> dplyr::pull(locus_name)

  pop <- pop |>
    get_table("genome_meta") |>
    dplyr::filter(locus_name %in% sel) |>
    define_qtl("ADG")

  cols <- DBI::dbListFields(pop$db_conn, "genome_meta")
  expect_true("is_QTL_ADG" %in% cols)

  result <- pop |> get_table("genome_meta") |> dplyr::collect()
  expect_equal(sum(result$is_QTL_ADG), 50)
  expect_equal(sum(!result$is_QTL_ADG), 450)

  close_pop(pop)
})


test_that("define_qtl() with multiple trait_names writes one column per trait", {
  pop <- make_qtl_pop("qtl_multi_trait")
  pop <- add_trait(pop, "BW", target_add_var = 0.5, residual_var = 0.5)

  set.seed(5)
  sel <- pop |> get_table("genome_meta") |> dplyr::collect() |>
    dplyr::slice_sample(n = 60) |> dplyr::pull(locus_name)

  pop <- pop |>
    get_table("genome_meta") |>
    dplyr::filter(locus_name %in% sel) |>
    define_qtl(c("ADG", "BW"))

  cols <- DBI::dbListFields(pop$db_conn, "genome_meta")
  expect_true("is_QTL_ADG" %in% cols)
  expect_true("is_QTL_BW"  %in% cols)

  result <- pop |> get_table("genome_meta") |> dplyr::collect()
  expect_equal(sum(result$is_QTL_ADG), 60)
  expect_equal(sum(result$is_QTL_BW),  60)
  # Same loci selected for both traits
  expect_true(all(result$is_QTL_ADG == result$is_QTL_BW))

  close_pop(pop)
})


test_that("define_qtl() defaults trait_name to all traits in trait_meta", {
  pop <- make_qtl_pop("qtl_default_trait")
  pop <- add_trait(pop, "BW", target_add_var = 0.5, residual_var = 0.5)

  set.seed(3)
  sel <- pop |> get_table("genome_meta") |> dplyr::collect() |>
    dplyr::slice_sample(n = 40) |> dplyr::pull(locus_name)

  pop <- pop |>
    get_table("genome_meta") |>
    dplyr::filter(locus_name %in% sel) |>
    define_qtl()  # no trait_name

  cols <- DBI::dbListFields(pop$db_conn, "genome_meta")
  expect_true("is_QTL_ADG" %in% cols)
  expect_true("is_QTL_BW"  %in% cols)

  result <- pop |> get_table("genome_meta") |> dplyr::collect()
  expect_equal(sum(result$is_QTL_ADG), 40)
  expect_equal(sum(result$is_QTL_BW),  40)

  close_pop(pop)
})


test_that("define_qtl() shared QTL via filter on existing is_QTL column", {
  pop <- make_qtl_pop("qtl_shared")
  pop <- add_trait(pop, "BW", target_add_var = 0.5, residual_var = 0.5)

  set.seed(10)
  sel <- pop |> get_table("genome_meta") |> dplyr::collect() |>
    dplyr::slice_sample(n = 50) |> dplyr::pull(locus_name)

  pop <- pop |>
    get_table("genome_meta") |>
    dplyr::filter(locus_name %in% sel) |>
    define_qtl("ADG")

  # BW gets exactly the same QTL as ADG
  pop <- pop |>
    get_table("genome_meta") |>
    dplyr::filter(is_QTL_ADG == TRUE) |>
    define_qtl("BW")

  result <- pop |> get_table("genome_meta") |> dplyr::collect()
  expect_equal(sum(result$is_QTL_BW), 50)
  expect_true(all(result$is_QTL_ADG == result$is_QTL_BW))

  close_pop(pop)
})


test_that("define_qtl() pipeline into add_additive_effects() works", {
  pop <- make_qtl_pop("qtl_pipeline")

  set.seed(42)
  sel <- pop |> get_table("genome_meta") |> dplyr::collect() |>
    dplyr::slice_sample(n = 100) |> dplyr::pull(locus_name)

  pop <- pop |>
    get_table("genome_meta") |>
    dplyr::filter(locus_name %in% sel) |>
    define_qtl("ADG") |>
    add_additive_effects("ADG", distribution = "normal", seed = 7L)

  cols <- DBI::dbListFields(pop$db_conn, "genome_meta")
  expect_true("add_ADG" %in% cols)

  gm <- pop |> get_table("genome_meta") |> dplyr::collect()
  expect_equal(sum(!is.na(gm$add_ADG)), 100)

  close_pop(pop)
})


test_that("define_qtl() errors without a prior add_trait()", {
  pop <- initialize_genome("qtl_no_trait", n_loci = 100, n_chr = 2,
                           chr_len_Mb = 100, db_path = ":memory:")

  expect_error(
    pop |>
      get_table("genome_meta") |>
      define_qtl("ADG"),
    "No traits defined|not found"
  )

  close_pop(pop)
})


test_that("define_qtl() errors when trait_name is not in trait_meta", {
  pop <- make_qtl_pop("qtl_bad_trait")

  expect_error(
    pop |>
      get_table("genome_meta") |>
      dplyr::filter(chr == 1L) |>
      define_qtl("NonExistentTrait"),
    "not found in trait_meta"
  )

  close_pop(pop)
})


test_that("define_qtl() errors on empty filter result", {
  pop <- make_qtl_pop("qtl_empty_filter")

  expect_error(
    pop |>
      get_table("genome_meta") |>
      dplyr::filter(chr == 999L) |>
      define_qtl("ADG"),
    "No loci selected"
  )

  close_pop(pop)
})


test_that("define_qtl() errors when QTL column already exists", {
  pop <- make_qtl_pop("qtl_dup")
  set.seed(2)
  sel <- pop |> get_table("genome_meta") |> dplyr::collect() |>
    dplyr::slice_sample(n = 30) |> dplyr::pull(locus_name)

  pop <- pop |>
    get_table("genome_meta") |>
    dplyr::filter(locus_name %in% sel) |>
    define_qtl("ADG")

  expect_error(
    pop |>
      get_table("genome_meta") |>
      dplyr::filter(locus_name %in% sel) |>
      define_qtl("ADG"),
    "already exists"
  )

  close_pop(pop)
})
