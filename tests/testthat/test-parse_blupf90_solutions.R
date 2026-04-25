test_that("parse_blupf90_solutions reads solutions.orig with header and alignment", {
  tmp <- tempfile()
  dir.create(tmp)
  on.exit(unlink(tmp, recursive = TRUE))

  # Reproduce the aligned, header-bearing format blupf90+ writes when OPTION origID is active.
  # Format: trait  effect  level  original_id  solution  (5 columns, multi-space aligned)
  writeLines(c(
    "trait effect level original_id solution",
    "   1   2         1 A-1           1.23456789",
    "   1   2         2 A-2          -0.98765432",
    "   1   2         3 A-3           0.00000000",
    "   1   1         1 10            5.00000000"   # mu (effect 1) — should be excluded
  ), file.path(tmp, "solutions.orig"))

  result <- tidybreed:::parse_blupf90_solutions(
    eval_dir          = tmp,
    trait             = c("ADG"),
    animal_effect_num = 2L,
    all_ped_ids       = c("A-1", "A-2", "A-3"),
    model             = "test_model",
    date_calc         = as.Date("2026-01-01")
  )

  expect_s3_class(result, "tbl_df")
  expect_equal(nrow(result), 3L)
  expect_equal(result$id_ind,     c("A-1", "A-2", "A-3"))
  expect_equal(result$trait_name, rep("ADG", 3))
  expect_equal(result$model,      rep("test_model", 3))
  expect_equal(result$ebv,        c(1.23456789, -0.98765432, 0.00000000),
               tolerance = 1e-6)
  expect_equal(result$date_calc,  rep(as.Date("2026-01-01"), 3))
})

test_that("parse_blupf90_solutions handles numeric-looking original IDs", {
  tmp <- tempfile()
  dir.create(tmp)
  on.exit(unlink(tmp, recursive = TRUE))

  writeLines(c(
    "trait effect level original_id solution",
    "   1   2         1 10           1.00000000",
    "   1   2         2 20          -2.00000000"
  ), file.path(tmp, "solutions.orig"))

  result <- tidybreed:::parse_blupf90_solutions(
    eval_dir          = tmp,
    trait             = c("BW"),
    animal_effect_num = 2L,
    all_ped_ids       = c("10", "20"),
    model             = "m1",
    date_calc         = as.Date("2026-01-01")
  )

  expect_equal(nrow(result), 2L)
  expect_equal(result$id_ind, c("10", "20"))
  expect_equal(result$ebv, c(1.0, -2.0), tolerance = 1e-6)
})

test_that("parse_blupf90_solutions errors when solutions.orig is missing", {
  tmp <- tempfile()
  dir.create(tmp)
  on.exit(unlink(tmp, recursive = TRUE))

  expect_error(
    tidybreed:::parse_blupf90_solutions(tmp, "ADG", 2L, character(), "m", Sys.Date()),
    "solutions\\.orig"
  )
})

test_that("parse_blupf90_solutions returns empty tibble when no IDs match", {
  tmp <- tempfile()
  dir.create(tmp)
  on.exit(unlink(tmp, recursive = TRUE))

  writeLines(c(
    "trait effect level original_id solution",
    "   1   2         1 A-1           1.00000000"
  ), file.path(tmp, "solutions.orig"))

  expect_warning(
    result <- tidybreed:::parse_blupf90_solutions(tmp, "ADG", 2L, c("B-999"), "m", Sys.Date()),
    "No animal EBVs matched"
  )
  expect_equal(nrow(result), 0L)
})
