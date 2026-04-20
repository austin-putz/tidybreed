make_trait_pop <- function(pop_name = "qtl", n_loci = 500, n_chr = 5) {
  pop <- initialize_genome(
    pop_name        = pop_name,
    n_loci          = n_loci,
    n_chr           = n_chr,
    chr_len_Mb      = 100,
    n_haplotypes    = 100,
    db_path         = ":memory:",
    fixed_allele_freq = 0.5
  )
  pop <- add_founders(pop, n_males = 20, n_females = 20, line_name = "A")
  pop <- add_trait(pop, "ADG", target_add_var = 0.25, residual_var = 0.75)
  pop
}


test_that("define_qtl() writes is_QTL_{trait} with requested count", {
  pop <- make_trait_pop("qtl_random")

  pop <- define_qtl(pop, "ADG", n = 50, method = "random")

  cols <- DBI::dbListFields(pop$db_conn, "genome_meta")
  expect_true("is_QTL_ADG" %in% cols)

  result <- dplyr::collect(get_table(pop, "genome_meta"))
  expect_equal(sum(result$is_QTL_ADG), 50)

  close_pop(pop)
})


test_that("define_qtl() errors without a prior add_trait()", {
  pop <- initialize_genome("qtl_no_trait", n_loci = 100, n_chr = 2,
                           chr_len_Mb = 100, db_path = ":memory:")
  expect_error(define_qtl(pop, "ADG", n = 10),
               "No traits defined|not found")
  close_pop(pop)
})


test_that("define_qtl() supports all selection methods", {
  pop <- make_trait_pop("qtl_methods", n_loci = 400, n_chr = 4)

  # even
  pop <- define_qtl(pop, "ADG", n = 40, method = "even")
  res <- dplyr::collect(get_table(pop, "genome_meta"))
  expect_true(sum(res$is_QTL_ADG) >= 35 && sum(res$is_QTL_ADG) <= 45)

  # by locus ids
  pop <- add_trait(pop, "BW", target_add_var = 0.5, residual_var = 0.5)
  pop <- define_qtl(pop, "BW", locus_ids = c(1L, 100L, 200L, 300L, 400L))
  res <- dplyr::collect(get_table(pop, "genome_meta"))
  expect_equal(sum(res$is_QTL_BW), 5)

  # by logical vector
  pop <- add_trait(pop, "FE", target_add_var = 0.5, residual_var = 0.5)
  tf <- rep(FALSE, 400); tf[seq(1, 400, 10)] <- TRUE
  pop <- define_qtl(pop, "FE", locus_tf = tf)
  res <- dplyr::collect(get_table(pop, "genome_meta"))
  expect_equal(sum(res$is_QTL_FE), 40)

  close_pop(pop)
})
