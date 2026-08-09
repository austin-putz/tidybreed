# assert_ploidy_2() guards against real (>2) ploidy for dosage/export. Stage 4
# added ind_meta.ploidy and chr_inheritance (sex chromosomes, organelles) —
# sex-linked chr_inheritance rows must NOT trip this guard; only ploidy != 2
# should.

test_that("assert_ploidy_2() passes silently when all individuals have ploidy 2", {
  pop <- make_pop_base(n_loci = 10, n_chr = 2, chr_len_Mb = 100)
  on.exit(close_pop(pop), add = TRUE)

  expect_null(assert_ploidy_2(pop))
})

test_that("assert_ploidy_2() errors when an individual has ploidy != 2", {
  pop <- make_test_pop(n_loci = 10, n_chr = 2, chr_len_Mb = 100,
                       n_males = 2, n_females = 2)
  on.exit(close_pop(pop), add = TRUE)

  ids <- DBI::dbGetQuery(pop$db_conn, "SELECT id_ind FROM ind_meta")$id_ind
  DBI::dbExecute(pop$db_conn, paste0(
    "UPDATE ind_meta SET ploidy = 4 WHERE id_ind = '", ids[1], "'"))

  expect_error(assert_ploidy_2(pop), "ploidy")
})

test_that("assert_ploidy_2() scoped to `ids` ignores bad individuals outside that set", {
  pop <- make_test_pop(n_loci = 10, n_chr = 2, chr_len_Mb = 100,
                       n_males = 2, n_females = 2)
  on.exit(close_pop(pop), add = TRUE)

  ids <- DBI::dbGetQuery(pop$db_conn, "SELECT id_ind FROM ind_meta ORDER BY id_ind")$id_ind
  bad_id  <- ids[1]
  good_id <- ids[2]
  DBI::dbExecute(pop$db_conn, paste0(
    "UPDATE ind_meta SET ploidy = 4 WHERE id_ind = '", bad_id, "'"))

  expect_null(assert_ploidy_2(pop, ids = good_id))
  expect_error(assert_ploidy_2(pop, ids = bad_id), "ploidy")
  expect_error(assert_ploidy_2(pop, ids = c(bad_id, good_id)), "ploidy")
})

test_that("assert_ploidy_2() does not trip on a sex-linked chr_inheritance row", {
  pop <- make_pop_base(n_loci = 10, n_chr = 2, chr_len_Mb = 100)
  on.exit(close_pop(pop), add = TRUE)

  suppressMessages(
    pop <- define_chromosome(pop, "1", offspring_sex = "M",
                             from_parent_1 = 0, from_parent_2 = 1))

  expect_null(assert_ploidy_2(pop))
})
