# define_chr(): sets a non-default chr_meta inheritance rule (sex chromosomes,
# organelles). Uses make_pop_base() (genome only, no founders needed).

test_that("define_chr() writes copy_mode/hemi_parent/recombines correctly", {
  pop <- make_pop_base(n_loci = 10, n_chr = 2, chr_len_Mb = 100)
  on.exit(close_pop(pop), add = TRUE)

  pop <- pop |> define_chr("1",
    copy_mode_M = "half", copy_mode_F = "full",
    hemi_parent = "parent_2", recombines = TRUE)

  row <- DBI::dbGetQuery(pop$db_conn, "SELECT * FROM chr_meta WHERE chr_name = '1'")
  expect_equal(row$copy_mode_M, "half")
  expect_equal(row$copy_mode_F, "full")
  expect_equal(row$hemi_parent, "parent_2")
  expect_true(row$recombines)

  # Untouched chromosome keeps the default row
  row2 <- DBI::dbGetQuery(pop$db_conn, "SELECT * FROM chr_meta WHERE chr_name = '2'")
  expect_equal(row2$copy_mode_M, "full")
  expect_equal(row2$copy_mode_F, "full")
  expect_true(is.na(row2$hemi_parent))
})

test_that("define_chr() rejects an unknown chr_name", {
  pop <- make_pop_base(n_loci = 10, n_chr = 2, chr_len_Mb = 100)
  on.exit(close_pop(pop), add = TRUE)

  expect_error(
    pop |> define_chr("Z9", copy_mode_M = "half", copy_mode_F = "none",
                      hemi_parent = "parent_1", recombines = FALSE),
    "not found in genome_meta"
  )
})

test_that("define_chr() rejects invalid copy_mode values", {
  pop <- make_pop_base(n_loci = 10, n_chr = 2, chr_len_Mb = 100)
  on.exit(close_pop(pop), add = TRUE)

  expect_error(
    pop |> define_chr("1", copy_mode_M = "triploid", copy_mode_F = "full",
                      hemi_parent = "parent_1"),
    "copy_mode_M"
  )
  expect_error(
    pop |> define_chr("1", copy_mode_M = "full", copy_mode_F = "bogus",
                      hemi_parent = "parent_1"),
    "copy_mode_F"
  )
})

test_that("define_chr() rejects hemi_parent supplied when both copy_modes are 'full'", {
  pop <- make_pop_base(n_loci = 10, n_chr = 2, chr_len_Mb = 100)
  on.exit(close_pop(pop), add = TRUE)

  expect_error(
    pop |> define_chr("1", copy_mode_M = "full", copy_mode_F = "full",
                      hemi_parent = "parent_1"),
    "hemi_parent must be NULL"
  )
})

test_that("define_chr() rejects missing/invalid hemi_parent when a copy_mode is non-full", {
  pop <- make_pop_base(n_loci = 10, n_chr = 2, chr_len_Mb = 100)
  on.exit(close_pop(pop), add = TRUE)

  expect_error(
    pop |> define_chr("1", copy_mode_M = "half", copy_mode_F = "full"),
    "hemi_parent must be"
  )
  expect_error(
    pop |> define_chr("1", copy_mode_M = "half", copy_mode_F = "full",
                      hemi_parent = "dad"),
    "hemi_parent must be"
  )
})

test_that("define_chr(overwrite = TRUE) updates an existing row in place", {
  pop <- make_pop_base(n_loci = 10, n_chr = 2, chr_len_Mb = 100)
  on.exit(close_pop(pop), add = TRUE)

  pop <- pop |> define_chr("1", copy_mode_M = "half", copy_mode_F = "full",
                           hemi_parent = "parent_2", recombines = TRUE)
  pop <- pop |> define_chr("1", copy_mode_M = "none", copy_mode_F = "half",
                           hemi_parent = "parent_1", recombines = FALSE,
                           overwrite = TRUE)

  row <- DBI::dbGetQuery(pop$db_conn, "SELECT * FROM chr_meta WHERE chr_name = '1'")
  expect_equal(row$copy_mode_M, "none")
  expect_equal(row$copy_mode_F, "half")
  expect_equal(row$hemi_parent, "parent_1")
  expect_false(row$recombines)
})

test_that("define_chr() does not perturb the RNG stream", {
  pop <- make_pop_base(n_loci = 10, n_chr = 2, chr_len_Mb = 100)
  on.exit(close_pop(pop), add = TRUE)

  set.seed(42)
  before <- runif(1)

  set.seed(42)
  pop <- pop |> define_chr("1", copy_mode_M = "half", copy_mode_F = "full",
                           hemi_parent = "parent_2")
  after <- runif(1)

  expect_equal(before, after)
})
