# assert_diploid_only() — Stage 3 forward-defense guard. chr_meta today only
# ever has default diploid rows (define_chr() is Stage 4 and doesn't exist
# yet), so tests mutate chr_meta directly via DBI to exercise the non-default
# path.

test_that("assert_diploid_only() passes silently for default all-diploid chr_meta", {
  pop <- make_pop_base(n_loci = 10, n_chr = 2, chr_len_Mb = 100)
  on.exit(close_pop(pop), add = TRUE)

  expect_null(assert_diploid_only(pop))
})

test_that("assert_diploid_only() errors when a chromosome has non-diploid copies_M", {
  pop <- make_pop_base(n_loci = 10, n_chr = 2, chr_len_Mb = 100)
  on.exit(close_pop(pop), add = TRUE)

  DBI::dbExecute(pop$db_conn, "UPDATE chr_meta SET copies_M = 1 WHERE chr_name = '1'")

  expect_error(assert_diploid_only(pop), "Stage 4")
})

test_that("assert_diploid_only() errors when a chromosome has non-diploid copies_F", {
  pop <- make_pop_base(n_loci = 10, n_chr = 2, chr_len_Mb = 100)
  on.exit(close_pop(pop), add = TRUE)

  DBI::dbExecute(pop$db_conn, "UPDATE chr_meta SET copies_F = 4 WHERE chr_name = '2'")

  expect_error(assert_diploid_only(pop), "Stage 4")
})
