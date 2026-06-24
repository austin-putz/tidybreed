# Helper: small pop with 20 individuals (10M, 10F) and a gen column
make_group_pop <- function(n_males = 10, n_females = 10) {
  pop <- make_test_pop(n_males = n_males, n_females = n_females)
  pop <- pop |>
    get_table("ind_meta") |>
    mutate_table(gen = 1L)
  pop
}


# ============================================================
# mutate_group_seq()
# ============================================================

test_that("mutate_group_seq() creates INTEGER column with correct group counts (scalar n_per_group)", {
  pop <- make_group_pop()
  pop <- pop |>
    get_table("ind_meta") |>
    mutate_group_seq(col_name = "pen", n_per_group = 5L, quiet = TRUE)

  result <- pop |> get_table("ind_meta") |> dplyr::collect()
  expect_true("pen" %in% names(result))
  expect_type(result$pen, "integer")

  # 20 animals / 5 per group = 4 full groups, 0 leftover
  expect_true(all(!is.na(result$pen)))
  expect_equal(sort(unique(result$pen)), 1L:4L)
  expect_true(all(table(result$pen) == 5L))

  close_pop(pop)
})


test_that("mutate_group_seq() leaves remainder NULL when include_leftover = FALSE", {
  pop <- make_group_pop(n_males = 7, n_females = 6)  # 13 animals
  expect_warning(
    pop <- pop |>
      get_table("ind_meta") |>
      mutate_group_seq(col_name = "pen", n_per_group = 5L, quiet = TRUE),
    "left unassigned"
  )

  result <- pop |> get_table("ind_meta") |> dplyr::collect()
  n_null    <- sum(is.na(result$pen))
  n_full    <- sum(!is.na(result$pen))
  expect_equal(n_full, 10L)   # 2 full groups of 5
  expect_equal(n_null, 3L)    # 3 leftover
  expect_equal(sort(unique(result$pen[!is.na(result$pen)])), 1L:2L)

  close_pop(pop)
})


test_that("mutate_group_seq() includes remainder in last group when include_leftover = TRUE", {
  pop <- make_group_pop(n_males = 7, n_females = 6)  # 13 animals
  pop <- pop |>
    get_table("ind_meta") |>
    mutate_group_seq(col_name = "pen", n_per_group = 5L, include_leftover = TRUE, quiet = TRUE)

  result <- pop |> get_table("ind_meta") |> dplyr::collect()
  expect_true(all(!is.na(result$pen)))
  expect_equal(sort(unique(result$pen)), 1L:3L)
  cnt <- as.integer(table(result$pen))
  expect_true(all(cnt[1:2] == 5L))
  expect_equal(cnt[3], 3L)

  close_pop(pop)
})


test_that("mutate_group_seq() same seed gives identical assignment", {
  pop1 <- make_group_pop()
  pop2 <- make_group_pop()

  pop1 <- pop1 |> get_table("ind_meta") |>
    mutate_group_seq(col_name = "pen", n_per_group = 5L, seed = 99L, quiet = TRUE)
  pop2 <- pop2 |> get_table("ind_meta") |>
    mutate_group_seq(col_name = "pen", n_per_group = 5L, seed = 99L, quiet = TRUE)

  r1 <- pop1 |> get_table("ind_meta") |> dplyr::arrange(id_ind) |> dplyr::collect()
  r2 <- pop2 |> get_table("ind_meta") |> dplyr::arrange(id_ind) |> dplyr::collect()
  expect_equal(r1$pen, r2$pen)

  close_pop(pop1); close_pop(pop2)
})


test_that("mutate_group_seq() auto-continues from MAX + 1", {
  pop <- make_group_pop()

  # Gen 1: groups 1-2 (males only, 10 / 5 = 2 groups)
  pop <- pop |>
    get_table("ind_meta") |>
    dplyr::filter(sex == "M") |>
    mutate_group_seq(col_name = "pen", n_per_group = 5L, quiet = TRUE)

  # Gen 2 style: all rows still NULL for females → auto-continues from 3
  pop <- pop |>
    get_table("ind_meta") |>
    dplyr::filter(sex == "F") |>
    mutate_group_seq(col_name = "pen", n_per_group = 5L, quiet = TRUE)

  result <- pop |> get_table("ind_meta") |> dplyr::collect()
  expect_equal(sort(unique(result$pen)), 1L:4L)

  close_pop(pop)
})


test_that("mutate_group_seq() start = 1L reuses labels in different filtered rows without error", {
  pop <- make_group_pop()

  pop <- pop |>
    get_table("ind_meta") |>
    dplyr::filter(sex == "M") |>
    mutate_group_seq(col_name = "pen", n_per_group = 5L, start = 1L, quiet = TRUE)

  # Females start at 1 too — this is label reuse, not overwrite
  pop <- pop |>
    get_table("ind_meta") |>
    dplyr::filter(sex == "F") |>
    mutate_group_seq(col_name = "pen", n_per_group = 5L, start = 1L, quiet = TRUE)

  result <- pop |> get_table("ind_meta") |> dplyr::collect()
  expect_true(all(!is.na(result$pen)))
  expect_equal(sort(unique(result$pen)), 1L:2L)

  close_pop(pop)
})


test_that("mutate_group_seq() n_groups mode: balanced split with include_leftover", {
  # 11 animals, 5 groups: floor(11/5) = 2, leftover = 1
  pop <- make_group_pop(n_males = 6, n_females = 5)  # 11 total
  pop <- pop |>
    get_table("ind_meta") |>
    mutate_group_seq(col_name = "grp", n_groups = 5L, include_leftover = TRUE, quiet = TRUE)

  result <- pop |> get_table("ind_meta") |> dplyr::collect()
  expect_true(all(!is.na(result$grp)))
  cnt <- as.integer(sort(table(result$grp)))
  expect_equal(sum(cnt), 11L)
  # First group gets the extra animal: one group of 3, four groups of 2
  expect_equal(sort(cnt), c(2L, 2L, 2L, 2L, 3L))

  close_pop(pop)
})


test_that("mutate_group_seq() n_groups mode: without include_leftover leaves remainder NULL", {
  pop <- make_group_pop(n_males = 6, n_females = 5)  # 11 total
  expect_warning(
    pop <- pop |>
      get_table("ind_meta") |>
      mutate_group_seq(col_name = "grp", n_groups = 5L, quiet = TRUE),
    "left unassigned"
  )

  result <- pop |> get_table("ind_meta") |> dplyr::collect()
  expect_equal(sum(!is.na(result$grp)), 10L)  # 5 groups × 2
  expect_equal(sum(is.na(result$grp)), 1L)

  close_pop(pop)
})


test_that("mutate_group_seq() variable n_per_group range assigns correctly", {
  pop <- make_group_pop()
  pop <- pop |>
    get_table("ind_meta") |>
    mutate_group_seq(
      col_name = "pen", n_per_group = c(3L, 7L),
      include_leftover = TRUE, seed = 1L, quiet = TRUE
    )

  result <- pop |> get_table("ind_meta") |> dplyr::collect()
  expect_true(all(!is.na(result$pen)))
  expect_true(all(result$pen >= 1L))

  close_pop(pop)
})


test_that("mutate_group_seq() overwrite = FALSE errors when filtered rows are already filled", {
  pop <- make_group_pop()
  pop <- pop |>
    get_table("ind_meta") |>
    mutate_group_seq(col_name = "pen", n_per_group = 5L, quiet = TRUE)

  expect_error(
    pop |> get_table("ind_meta") |>
      mutate_group_seq(col_name = "pen", n_per_group = 5L, overwrite = FALSE),
    "already have non-NULL values"
  )

  close_pop(pop)
})


test_that("mutate_group_seq() overwrite = TRUE replaces existing values", {
  pop <- make_group_pop()
  pop <- pop |>
    get_table("ind_meta") |>
    mutate_group_seq(col_name = "pen", n_per_group = 5L, quiet = TRUE)

  expect_no_error(
    pop <- pop |>
      get_table("ind_meta") |>
      mutate_group_seq(col_name = "pen", n_per_group = 5L,
                       overwrite = TRUE, seed = 7L, quiet = TRUE)
  )

  result <- pop |> get_table("ind_meta") |> dplyr::collect()
  expect_true(all(!is.na(result$pen)))

  close_pop(pop)
})


test_that("mutate_group_seq() empty filter warns and returns unchanged", {
  pop <- make_group_pop()
  expect_warning(
    pop <- pop |>
      get_table("ind_meta") |>
      dplyr::filter(sex == "X") |>
      mutate_group_seq(col_name = "pen", n_per_group = 5L),
    "0 rows"
  )
  expect_false("pen" %in% DBI::dbListFields(pop$db_conn, "ind_meta"))

  close_pop(pop)
})


test_that("mutate_group_seq() errors on incompatible existing column type", {
  pop <- make_group_pop()
  pop <- pop |>
    get_table("ind_meta") |>
    mutate_table(pen = "string_pen")  # VARCHAR

  expect_error(
    pop |> get_table("ind_meta") |>
      mutate_group_seq(col_name = "pen", n_per_group = 5L),
    "incompatible"
  )

  close_pop(pop)
})


test_that("mutate_group_seq() errors when n_groups > n_rows", {
  pop <- make_group_pop(n_males = 3, n_females = 2)  # 5 total
  expect_error(
    pop |> get_table("ind_meta") |>
      mutate_group_seq(col_name = "pen", n_groups = 10L),
    "exceeds"
  )
  close_pop(pop)
})


# ============================================================
# mutate_group_named()
# ============================================================

test_that("mutate_group_named() creates VARCHAR column with correct counts", {
  pop <- make_group_pop()
  pop <- pop |>
    get_table("ind_meta") |>
    mutate_group_named(
      col_name    = "farm",
      group_names = c("farm_A", "farm_B"),
      counts      = c(12L, 8L),
      quiet       = TRUE
    )

  result <- pop |> get_table("ind_meta") |> dplyr::collect()
  expect_true("farm" %in% names(result))
  expect_type(result$farm, "character")
  expect_equal(sum(result$farm == "farm_A", na.rm = TRUE), 12L)
  expect_equal(sum(result$farm == "farm_B", na.rm = TRUE), 8L)

  close_pop(pop)
})


test_that("mutate_group_named() proportions exact method gives correct split", {
  pop <- make_group_pop()  # 20 animals
  pop <- pop |>
    get_table("ind_meta") |>
    mutate_group_named(
      col_name    = "farm",
      group_names = c("A", "B"),
      proportions = c(0.5, 0.5),
      quiet       = TRUE
    )

  result <- pop |> get_table("ind_meta") |> dplyr::collect()
  expect_equal(sum(result$farm == "A", na.rm = TRUE), 10L)
  expect_equal(sum(result$farm == "B", na.rm = TRUE), 10L)

  close_pop(pop)
})


test_that("mutate_group_named() accepts proportions within tolerance", {
  pop <- make_group_pop()
  expect_no_error(
    pop |> get_table("ind_meta") |>
      mutate_group_named(
        col_name    = "farm",
        group_names = c("A", "B", "C"),
        proportions = c(1/3, 1/3, 1/3),
        quiet       = TRUE
      )
  )
  close_pop(pop)
})


test_that("mutate_group_named() largest-remainder is stable left-to-right", {
  # 10 animals, 1/3 each => floors = c(3,3,3), remainder 1 goes to first group
  pop <- make_group_pop(n_males = 5, n_females = 5)  # 10
  pop <- pop |>
    get_table("ind_meta") |>
    mutate_group_named(
      col_name    = "grp",
      group_names = c("X", "Y", "Z"),
      proportions = c(1/3, 1/3, 1/3),
      seed        = 1L,
      quiet       = TRUE
    )

  result <- pop |> get_table("ind_meta") |> dplyr::collect()
  cnt <- table(result$grp)
  expect_equal(as.integer(cnt["X"]), 4L)   # wins the tie
  expect_equal(as.integer(cnt["Y"]), 3L)
  expect_equal(as.integer(cnt["Z"]), 3L)

  close_pop(pop)
})


test_that("mutate_group_named() prop_method = 'random' assigns all rows", {
  pop <- make_group_pop()
  pop <- pop |>
    get_table("ind_meta") |>
    mutate_group_named(
      col_name    = "farm",
      group_names = c("A", "B"),
      proportions = c(0.6, 0.4),
      prop_method = "random",
      seed        = 5L,
      quiet       = TRUE
    )

  result <- pop |> get_table("ind_meta") |> dplyr::collect()
  expect_true(all(!is.na(result$farm)))
  expect_true(all(result$farm %in% c("A", "B")))

  close_pop(pop)
})


test_that("mutate_group_named() same seed gives identical assignments", {
  pop1 <- make_group_pop()
  pop2 <- make_group_pop()

  pop1 <- pop1 |> get_table("ind_meta") |>
    mutate_group_named(col_name = "farm", group_names = c("A", "B"),
                       proportions = c(0.6, 0.4), seed = 42L, quiet = TRUE)
  pop2 <- pop2 |> get_table("ind_meta") |>
    mutate_group_named(col_name = "farm", group_names = c("A", "B"),
                       proportions = c(0.6, 0.4), seed = 42L, quiet = TRUE)

  r1 <- pop1 |> get_table("ind_meta") |> dplyr::arrange(id_ind) |> dplyr::collect()
  r2 <- pop2 |> get_table("ind_meta") |> dplyr::arrange(id_ind) |> dplyr::collect()
  expect_equal(r1$farm, r2$farm)

  close_pop(pop1); close_pop(pop2)
})


test_that("mutate_group_named() counts leftover warns and leaves NULL", {
  pop <- make_group_pop()  # 20 animals
  expect_warning(
    pop <- pop |>
      get_table("ind_meta") |>
      mutate_group_named(
        col_name    = "grp",
        group_names = c("A", "B"),
        counts      = c(8L, 8L),
        quiet       = TRUE
      ),
    "unassigned"
  )

  result <- pop |> get_table("ind_meta") |> dplyr::collect()
  expect_equal(sum(is.na(result$grp)), 4L)

  close_pop(pop)
})


test_that("mutate_group_named() counts leftover added to last group when include_leftover = TRUE", {
  pop <- make_group_pop()  # 20 animals
  pop <- pop |>
    get_table("ind_meta") |>
    mutate_group_named(
      col_name         = "grp",
      group_names      = c("A", "B"),
      counts           = c(8L, 8L),
      include_leftover = TRUE,
      quiet            = TRUE
    )

  result <- pop |> get_table("ind_meta") |> dplyr::collect()
  expect_true(all(!is.na(result$grp)))
  expect_equal(sum(result$grp == "A"), 8L)
  expect_equal(sum(result$grp == "B"), 12L)

  close_pop(pop)
})


test_that("mutate_group_named() underfill_action = 'error' stops before sampling", {
  pop <- make_group_pop(n_males = 3, n_females = 2)  # 5 animals
  expect_error(
    pop |> get_table("ind_meta") |>
      mutate_group_named(
        col_name    = "farm",
        group_names = c("A", "B"),
        counts      = c(10L, 10L)
      ),
    "fewer than the total"
  )
  close_pop(pop)
})


test_that("mutate_group_named() underfill_action = 'proportional' rescales counts", {
  pop <- make_group_pop(n_males = 3, n_females = 2)  # 5 animals
  pop <- pop |>
    get_table("ind_meta") |>
    mutate_group_named(
      col_name         = "farm",
      group_names      = c("A", "B"),
      counts           = c(4L, 4L),
      underfill_action = "proportional",
      quiet            = TRUE
    )

  result <- pop |> get_table("ind_meta") |> dplyr::collect()
  expect_true(all(!is.na(result$farm)))
  expect_equal(nrow(result), 5L)

  close_pop(pop)
})


test_that("mutate_group_named() underfill_action = 'sequential' fills in order", {
  pop <- make_group_pop(n_males = 3, n_females = 2)  # 5 animals
  pop <- pop |>
    get_table("ind_meta") |>
    mutate_group_named(
      col_name         = "farm",
      group_names      = c("A", "B", "C"),
      counts           = c(3L, 3L, 3L),
      underfill_action = "sequential",
      quiet            = TRUE
    )

  result <- pop |> get_table("ind_meta") |> dplyr::collect()
  # C gets 0 (runs out), B gets 2, A gets 3 — or A gets 3, B gets 2, C gets 0
  expect_true(all(result$farm[!is.na(result$farm)] %in% c("A", "B", "C")))
  expect_equal(sum(!is.na(result$farm)), 5L)
  expect_equal(sum(result$farm == "A", na.rm = TRUE), 3L)
  expect_equal(sum(result$farm == "B", na.rm = TRUE), 2L)
  expect_equal(sum(result$farm == "C", na.rm = TRUE), 0L)

  close_pop(pop)
})


test_that("mutate_group_named() overwrite = FALSE errors when rows filled", {
  pop <- make_group_pop()
  pop <- pop |>
    get_table("ind_meta") |>
    mutate_group_named(
      col_name    = "farm",
      group_names = c("A", "B"),
      counts      = c(10L, 10L),
      quiet       = TRUE
    )

  expect_error(
    pop |> get_table("ind_meta") |>
      mutate_group_named(
        col_name    = "farm",
        group_names = c("A", "B"),
        counts      = c(10L, 10L)
      ),
    "already have non-NULL values"
  )

  close_pop(pop)
})


test_that("mutate_group_named() empty filter warns and returns unchanged", {
  pop <- make_group_pop()
  expect_warning(
    pop |>
      get_table("ind_meta") |>
      dplyr::filter(sex == "X") |>
      mutate_group_named(
        col_name    = "farm",
        group_names = c("A"),
        counts      = c(1L)
      ),
    "0 rows"
  )
  expect_false("farm" %in% DBI::dbListFields(pop$db_conn, "ind_meta"))

  close_pop(pop)
})


test_that("mutate_group_named() errors when proportions do not sum near 1", {
  pop <- make_group_pop()
  expect_error(
    pop |> get_table("ind_meta") |>
      mutate_group_named(
        col_name    = "farm",
        group_names = c("A", "B"),
        proportions = c(0.3, 0.3)
      ),
    "sum to approximately 1"
  )
  close_pop(pop)
})


# ============================================================
# mutate_group_concatenate()
# ============================================================

test_that("mutate_group_concatenate() creates VARCHAR column from two columns", {
  pop <- make_group_pop()
  pop <- pop |>
    get_table("ind_meta") |>
    mutate_group_concatenate(
      col_name       = "cg",
      source_columns = c("sex", "line_name"),
      quiet          = TRUE
    )

  result <- pop |> get_table("ind_meta") |> dplyr::collect()
  expect_true("cg" %in% names(result))
  expect_type(result$cg, "character")
  # All animals have sex and line_name → no NAs
  expect_true(all(!is.na(result$cg)))
  expect_true(all(grepl("_", result$cg)))

  close_pop(pop)
})


test_that("mutate_group_concatenate() uses custom separator", {
  pop <- make_group_pop()
  pop <- pop |>
    get_table("ind_meta") |>
    mutate_group_concatenate(
      col_name       = "cg",
      source_columns = c("sex", "line_name"),
      sep            = ":",
      quiet          = TRUE
    )

  result <- pop |> get_table("ind_meta") |> dplyr::collect()
  expect_true(all(grepl(":", result$cg, fixed = TRUE)))
  expect_false(any(grepl("_", result$cg, fixed = TRUE)))

  close_pop(pop)
})


test_that("mutate_group_concatenate() null_action = 'propagate' produces NA for rows with NULL source", {
  pop <- make_group_pop()
  # founders have id_parent_1 = NA; propagate should yield NA
  pop <- pop |>
    get_table("ind_meta") |>
    mutate_group_concatenate(
      col_name       = "sib_grp",
      source_columns = c("id_parent_1", "id_parent_2"),
      null_action    = "propagate",
      quiet          = TRUE
    )

  result <- pop |> get_table("ind_meta") |> dplyr::collect()
  # Founders have NA parents → result should be NA
  expect_true(all(is.na(result$sib_grp)))

  close_pop(pop)
})


test_that("mutate_group_concatenate() null_action = 'skip' omits NULL sources", {
  pop <- make_group_pop()
  # Founders have id_parent_1 = NA — skipping should produce just id_parent_2 value
  # But since both parents are NA for founders, all-NA row → NA result
  pop <- pop |>
    get_table("ind_meta") |>
    mutate_group_concatenate(
      col_name       = "sib_grp",
      source_columns = c("sex", "id_parent_1"),
      null_action    = "skip",
      quiet          = TRUE
    )

  result <- pop |> get_table("ind_meta") |> dplyr::collect()
  # sex is always non-NULL; id_parent_1 is NULL → skipped; result = just sex value
  expect_true(all(!is.na(result$sib_grp)))
  expect_true(all(result$sib_grp %in% c("M", "F")))

  close_pop(pop)
})


test_that("mutate_group_concatenate() null_action = 'literal' replaces NULL with 'NA'", {
  pop <- make_group_pop()
  pop <- pop |>
    get_table("ind_meta") |>
    mutate_group_concatenate(
      col_name       = "lbl",
      source_columns = c("sex", "id_parent_1"),
      null_action    = "literal",
      quiet          = TRUE
    )

  result <- pop |> get_table("ind_meta") |> dplyr::collect()
  expect_true(all(!is.na(result$lbl)))
  # Should contain "NA" literal where parent is missing
  expect_true(any(grepl("NA", result$lbl)))

  close_pop(pop)
})


test_that("mutate_group_concatenate() respects filter", {
  pop <- make_group_pop()
  pop <- pop |>
    get_table("ind_meta") |>
    dplyr::filter(sex == "M") |>
    mutate_group_concatenate(
      col_name       = "lbl",
      source_columns = c("sex", "line_name"),
      quiet          = TRUE
    )

  result <- pop |> get_table("ind_meta") |> dplyr::collect()
  expect_equal(sum(!is.na(result$lbl)), 10L)  # only males assigned
  expect_true(all(is.na(result$lbl[result$sex == "F"])))

  close_pop(pop)
})


test_that("mutate_group_concatenate() errors on missing source column", {
  pop <- make_group_pop()
  expect_error(
    pop |> get_table("ind_meta") |>
      mutate_group_concatenate(
        col_name       = "cg",
        source_columns = c("sex", "nonexistent_col")
      ),
    "not found"
  )
  close_pop(pop)
})


test_that("mutate_group_concatenate() errors when col_name is in source_columns", {
  pop <- make_group_pop()
  # First create the column
  pop <- pop |>
    get_table("ind_meta") |>
    mutate_group_concatenate(
      col_name       = "cg",
      source_columns = c("sex", "line_name"),
      quiet          = TRUE
    )
  expect_error(
    pop |> get_table("ind_meta") |>
      mutate_group_concatenate(
        col_name       = "cg",
        source_columns = c("sex", "cg"),
        overwrite      = TRUE
      ),
    "cannot be one of the"
  )
  close_pop(pop)
})


test_that("mutate_group_concatenate() overwrite = FALSE errors when rows filled", {
  pop <- make_group_pop()
  pop <- pop |>
    get_table("ind_meta") |>
    mutate_group_concatenate(
      col_name       = "cg",
      source_columns = c("sex", "line_name"),
      quiet          = TRUE
    )

  expect_error(
    pop |> get_table("ind_meta") |>
      mutate_group_concatenate(
        col_name       = "cg",
        source_columns = c("sex", "line_name")
      ),
    "already have non-NULL values"
  )

  close_pop(pop)
})


test_that("mutate_group_concatenate() empty filter warns and returns unchanged", {
  pop <- make_group_pop()
  expect_warning(
    pop |>
      get_table("ind_meta") |>
      dplyr::filter(sex == "X") |>
      mutate_group_concatenate(
        col_name       = "cg",
        source_columns = c("sex", "line_name")
      ),
    "0 rows"
  )
  expect_false("cg" %in% DBI::dbListFields(pop$db_conn, "ind_meta"))

  close_pop(pop)
})


# ============================================================
# Integration: group column used with define_effect_random()
# ============================================================

test_that("pen column from mutate_group_seq() is compatible with define_effect_random()", {
  pop <- make_group_pop()
  pop <- pop |>
    get_table("ind_meta") |>
    mutate_group_seq(col_name = "pen", n_per_group = 5L, quiet = TRUE)

  # Verify column type is INTEGER — required for source_column in define_effect_random
  result <- pop |> get_table("ind_meta") |> dplyr::collect()
  expect_type(result$pen, "integer")
  expect_true("pen" %in% DBI::dbListFields(pop$db_conn, "ind_meta"))

  close_pop(pop)
})
