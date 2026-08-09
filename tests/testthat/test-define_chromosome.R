# define_chromosome(): sets a non-default chromosome rule across the two long
# tables (chr_inheritance / chr_recombination). One concern per call. Uses
# make_pop_base() (genome only, no founders needed).

.inh_row <- function(pop, chr, sex = NULL) {
  pred <- if (is.null(sex)) "offspring_sex IS NULL" else sprintf("offspring_sex = '%s'", sex)
  DBI::dbGetQuery(pop$db_conn, sprintf(
    "SELECT * FROM chr_inheritance WHERE chr_name = '%s' AND %s", chr, pred))
}
.rec_row <- function(pop, chr, sex = NULL) {
  pred <- if (is.null(sex)) "parent_sex IS NULL" else sprintf("parent_sex = '%s'", sex)
  DBI::dbGetQuery(pop$db_conn, sprintf(
    "SELECT * FROM chr_recombination WHERE chr_name = '%s' AND %s", chr, pred))
}

test_that("inheritance call writes chr_inheritance keyed by offspring_sex", {
  pop <- make_pop_base(n_loci = 10, n_chr = 2, chr_len_Mb = 100)
  on.exit(close_pop(pop), add = TRUE)

  suppressMessages(
    pop <- define_chromosome(pop, "1", offspring_sex = "M",
                             from_parent_1 = 0, from_parent_2 = 1))
  row <- .inh_row(pop, "1", "M")
  expect_equal(row$from_parent_1, 0L)
  expect_equal(row$from_parent_2, 1L)
  expect_true(is.na(row$line_name))

  # Nothing written to chr_recombination, and the seeded default rows survive.
  expect_equal(nrow(.rec_row(pop, "1", NULL)), 1L)
  expect_equal(.inh_row(pop, "1", NULL)$from_parent_1, 1L)   # seed untouched
  expect_equal(nrow(.inh_row(pop, "2", "M")), 0L)            # sibling chr untouched
})

test_that("recombination call writes chr_recombination keyed by parent_sex", {
  pop <- make_pop_base(n_loci = 10, n_chr = 2, chr_len_Mb = 100)
  on.exit(close_pop(pop), add = TRUE)

  suppressMessages(pop <- define_chromosome(pop, "1", recombines = FALSE))
  row <- .rec_row(pop, "1", NULL)
  expect_false(row$recombines)
  # No inheritance row added; the seeded NULL default inheritance is intact.
  expect_equal(nrow(.inh_row(pop, "1", "M")), 0L)
})

test_that("define_chromosome() errors when both concerns are supplied", {
  pop <- make_pop_base(n_loci = 10, n_chr = 2, chr_len_Mb = 100)
  on.exit(close_pop(pop), add = TRUE)
  expect_error(
    define_chromosome(pop, "1", from_parent_1 = 0, from_parent_2 = 1,
                      recombines = FALSE),
    "exactly one concern")
})

test_that("define_chromosome() errors when no concern is supplied", {
  pop <- make_pop_base(n_loci = 10, n_chr = 2, chr_len_Mb = 100)
  on.exit(close_pop(pop), add = TRUE)
  expect_error(define_chromosome(pop, "1"), "needs a concern")
})

test_that("define_chromosome() rejects the other concern's sex key", {
  pop <- make_pop_base(n_loci = 10, n_chr = 2, chr_len_Mb = 100)
  on.exit(close_pop(pop), add = TRUE)
  # parent_sex on an inheritance call
  expect_error(
    define_chromosome(pop, "1", parent_sex = "M",
                      from_parent_1 = 0, from_parent_2 = 1),
    "parent_sex is a recombination key")
  # offspring_sex on a recombination call
  expect_error(
    define_chromosome(pop, "1", offspring_sex = "M", recombines = FALSE),
    "offspring_sex is an inheritance key")
})

test_that("define_chromosome() enforces from_parent_1 + from_parent_2 <= 2", {
  pop <- make_pop_base(n_loci = 10, n_chr = 2, chr_len_Mb = 100)
  on.exit(close_pop(pop), add = TRUE)
  expect_error(
    define_chromosome(pop, "1", from_parent_1 = 2, from_parent_2 = 2),
    "<= 2")
  expect_error(
    define_chromosome(pop, "1", from_parent_1 = 1),   # not supplied together
    "supplied together")
})

test_that("define_chromosome() rejects an unknown chr_name", {
  pop <- make_pop_base(n_loci = 10, n_chr = 2, chr_len_Mb = 100)
  on.exit(close_pop(pop), add = TRUE)
  expect_error(
    define_chromosome(pop, "Z9", from_parent_1 = 0, from_parent_2 = 1),
    "not found in genome_meta")
})

test_that("overwrite = TRUE upserts the same logical key in place", {
  pop <- make_pop_base(n_loci = 10, n_chr = 2, chr_len_Mb = 100)
  on.exit(close_pop(pop), add = TRUE)

  suppressMessages({
    pop <- define_chromosome(pop, "1", offspring_sex = "M",
                             from_parent_1 = 0, from_parent_2 = 1)
    pop <- define_chromosome(pop, "1", offspring_sex = "M",
                             from_parent_1 = 1, from_parent_2 = 0)
  })
  row <- .inh_row(pop, "1", "M")
  expect_equal(nrow(row), 1L)
  expect_equal(row$from_parent_1, 1L)
  expect_equal(row$from_parent_2, 0L)
})

test_that("overwrite = FALSE errors on an existing NULL-safe key across all four combos", {
  # One chromosome per combo so the shadowing validator can't couple them.
  pop <- make_pop_base(n_loci = 16, n_chr = 4, chr_len_Mb = 100)
  on.exit(close_pop(pop), add = TRUE)

  # (sex=NULL, line=NULL) on chr "1": collides with the define_genome() seed.
  expect_error(
    define_chromosome(pop, "1", from_parent_1 = 1, from_parent_2 = 1,
                      overwrite = FALSE),
    "already exists")

  # (sex=M, line=NULL) on chr "2": first insert OK, second collides.
  suppressMessages(
    pop <- define_chromosome(pop, "2", offspring_sex = "M",
                             from_parent_1 = 0, from_parent_2 = 1, overwrite = FALSE))
  expect_error(
    define_chromosome(pop, "2", offspring_sex = "M",
                      from_parent_1 = 1, from_parent_2 = 0, overwrite = FALSE),
    "already exists")

  # (sex=NULL, line=L) on chr "3": first insert OK, second collides.
  suppressMessages(
    pop <- define_chromosome(pop, "3", line_name = "LA",
                             from_parent_1 = 0, from_parent_2 = 1, overwrite = FALSE))
  expect_error(
    define_chromosome(pop, "3", line_name = "LA",
                      from_parent_1 = 1, from_parent_2 = 0, overwrite = FALSE),
    "already exists")

  # (sex=M, line=L) on chr "4": distinct key, inserts without collision.
  expect_error(
    suppressMessages(define_chromosome(pop, "4", offspring_sex = "M", line_name = "LA",
                                       from_parent_1 = 0, from_parent_2 = 1,
                                       overwrite = FALSE)),
    NA)
})

test_that("a failed validation rolls back, preserving the prior valid row", {
  pop <- make_pop_base(n_loci = 10, n_chr = 2, chr_len_Mb = 100)
  on.exit(close_pop(pop), add = TRUE)

  # Establish a valid sex-specific row.
  suppressMessages(
    pop <- define_chromosome(pop, "1", offspring_sex = "M",
                             from_parent_1 = 0, from_parent_2 = 1))
  before <- .inh_row(pop, "1", "M")

  # Force a shadowing conflict: a (sex=M, line=NULL) rule and a (sex=NULL,
  # line=L) rule that differ, with no explicit (M, L) row. The (NULL, L) write
  # must be rejected AND rolled back.
  expect_error(
    define_chromosome(pop, "1", line_name = "LB",
                      from_parent_1 = 1, from_parent_2 = 1),
    "shadowing conflict")

  # The (sex=NULL, line=LB) row must NOT have been left behind.
  left <- DBI::dbGetQuery(pop$db_conn,
    "SELECT COUNT(*) AS n FROM chr_inheritance WHERE chr_name = '1' AND line_name = 'LB'")$n
  expect_equal(left, 0L)
  # The prior (M) row is intact.
  after <- .inh_row(pop, "1", "M")
  expect_equal(after$from_parent_1, before$from_parent_1)
  expect_equal(after$from_parent_2, before$from_parent_2)
})

test_that("define_chromosome() does not perturb the RNG stream", {
  pop <- make_pop_base(n_loci = 10, n_chr = 2, chr_len_Mb = 100)
  on.exit(close_pop(pop), add = TRUE)

  set.seed(42)
  before <- runif(1)

  set.seed(42)
  suppressMessages(
    pop <- define_chromosome(pop, "1", offspring_sex = "M",
                             from_parent_1 = 0, from_parent_2 = 1))
  after <- runif(1)

  expect_equal(before, after)
})
