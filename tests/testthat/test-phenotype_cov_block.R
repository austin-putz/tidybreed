# Covariance blocks in phenotype_var_comp — D1 (whole-block declaration), the
# residual strata rules, D3 (realization lock), D5 (diagonal writers), D6
# (condition_change_action agreement) and §5.6 (named-effect compatibility).
# See plans/sample_correlated_effects.md and phase_2.md.

sym <- function(names, values) {
  n <- length(names)
  matrix(values, n, n, byrow = TRUE, dimnames = list(names, names))
}

# A population with traits A, B, C (QTL placed) but no phenotypes defined.
make_block_pop <- function(pop_name = "blk", traits = c("A", "B", "C")) {
  pop <- suppressMessages(make_test_pop(pop_name, n_loci = 60, n_chr = 1,
                                        n_males = 6, n_females = 6))
  for (t in traits) {
    pop <- suppressMessages(define_trait(pop, t, target_add_var = 1))
    pop <- suppressMessages(
      pop |> get_table("genome_meta") |> define_additive_effects(t))
  }
  pop
}

resid_rows <- function(pop, where = "") {
  DBI::dbGetQuery(pop$db_conn, paste0(
    "SELECT * FROM phenotype_var_comp WHERE effect_name = 'residual' ", where,
    " ORDER BY condition_level NULLS FIRST, phenotype_name_1, phenotype_name_2"))
}

effect_rows <- function(pop, effect) {
  DBI::dbGetQuery(pop$db_conn, paste0(
    "SELECT * FROM phenotype_var_comp WHERE effect_name = '", effect,
    "' ORDER BY phenotype_name_1, phenotype_name_2"))
}


# ── D1: a block is declared in one call, as a complete matrix ───────────────

test_that("D1: a fragment of a block errors, names the missing pair, and changes nothing", {
  pop <- make_block_pop("d1_fragment")
  on.exit(close_pop(pop))

  pop <- define_residual_cov(pop, c("A", "B"), sym(c("A", "B"), c(1, .3, .3, 2)))
  before <- resid_rows(pop)
  expect_equal(nrow(before), 4L)

  expect_error(
    define_residual_cov(pop, c("B", "C"), sym(c("B", "C"), c(2, .1, .1, 1))),
    "omits: A \\(paired with B\\)")
  expect_error(
    define_residual_cov(pop, c("B", "C"), sym(c("B", "C"), c(2, .1, .1, 1))),
    'define_residual_cov\\(pop, c\\("A", "B", "C"\\), R\\)')
  expect_identical(resid_rows(pop), before)
})

test_that("D1: the complete block succeeds and an explicit zero is a stored pair row", {
  pop <- make_block_pop("d1_complete")
  on.exit(close_pop(pop))

  R <- sym(c("A", "B", "C"), c(1, .5, 0,
                               .5, 2, .1,
                               0, .1, 3))
  set.seed(11); seed_before <- .Random.seed
  pop <- define_residual_cov(pop, c("A", "B", "C"), R)
  expect_identical(.Random.seed, seed_before)   # writer is RNG-neutral

  rows <- resid_rows(pop)
  expect_equal(nrow(rows), 9L)
  ac <- rows[rows$phenotype_name_1 == "A" & rows$phenotype_name_2 == "C", ]
  expect_equal(nrow(ac), 1L)
  expect_equal(ac$cov_value, 0)
  expect_true(all(is.na(rows$condition_column)))
  expect_true(all(is.na(rows$condition_table)))

  # Values land in the right cells whatever the matrix's dimname order
  R_perm <- R[c("C", "A", "B"), c("C", "A", "B")]
  pop <- define_residual_cov(pop, c("A", "B", "C"), R_perm)
  rows <- resid_rows(pop)
  expect_equal(nrow(rows), 9L)
  expect_equal(rows$cov_value[rows$phenotype_name_1 == "B" &
                              rows$phenotype_name_2 == "C"], .1)
})

test_that("D1: non-PSD, asymmetric, non-finite and misnamed matrices are rejected at define time", {
  pop <- make_block_pop("d1_matrix")
  on.exit(close_pop(pop))

  expect_error(define_residual_cov(pop, c("A", "B"), sym(c("A", "B"), c(1, 2, 2, 1))),
               "positive semi-definite")
  expect_error(define_residual_cov(pop, c("A", "B"), sym(c("A", "B"), c(1, .2, .3, 1))),
               "symmetric")
  expect_error(define_residual_cov(pop, c("A", "B"), sym(c("A", "B"), c(1, NA, NA, 1))),
               "finite")
  expect_error(define_residual_cov(pop, "A", matrix(-1, 1, 1, dimnames = list("A", "A"))),
               "non-negative")
  expect_error(define_residual_cov(pop, c("A", "B"), sym(c("A", "X"), c(1, 0, 0, 1))),
               "dimnames")
  expect_equal(nrow(resid_rows(pop)), 0L)

  # PSD with a zero eigenvalue is fine (perfect correlation)
  pop <- define_residual_cov(pop, c("A", "B"), sym(c("A", "B"), c(1, 1, 1, 1)))
  expect_equal(nrow(resid_rows(pop)), 4L)
})

test_that("D1: a strict subset of an existing block errors and names the omitted member", {
  pop <- make_block_pop("d1_subset")
  on.exit(close_pop(pop))

  pop <- define_residual_cov(pop, c("A", "B", "C"),
                             sym(c("A", "B", "C"), c(1, 0, 0, 0, 1, 0, 0, 0, 1)))
  expect_error(define_residual_cov(pop, c("A", "B"), sym(c("A", "B"), c(1, 0, 0, 1))),
               "omits: C")
  expect_error(define_residual_cov(pop, "A", matrix(5, 1, 1, dimnames = list("A", "A"))),
               "block is \\{A, B, C\\}")
  expect_equal(nrow(resid_rows(pop)), 9L)
})

test_that("D1: redeclaring the same block replaces its rows exactly once", {
  pop <- make_block_pop("d1_redeclare")
  on.exit(close_pop(pop))

  pop <- define_residual_cov(pop, c("A", "B"), sym(c("A", "B"), c(1, .3, .3, 2)))
  pop <- define_residual_cov(pop, c("A", "B"), sym(c("A", "B"), c(4, .6, .6, 5)))
  rows <- resid_rows(pop)
  expect_equal(nrow(rows), 4L)
  expect_equal(rows$cov_value[rows$phenotype_name_1 == "A" &
                              rows$phenotype_name_2 == "A"], 4)
})

test_that("D1: two singleton blocks merge into one block when declared together", {
  # The README / vignette pattern: residual_var per phenotype, then a full R.
  pop <- make_block_pop("d1_merge")
  on.exit(close_pop(pop))

  pop <- suppressMessages(define_phenotype(pop, "A", residual_var = 1))
  pop <- suppressMessages(define_phenotype(pop, "B", residual_var = 2))
  expect_equal(nrow(resid_rows(pop)), 2L)

  pop <- suppressMessages(
    define_effect_cov_matrix(pop, "residual", sym(c("A", "B"), c(1, .4, .4, 2))))
  rows <- resid_rows(pop)
  expect_equal(nrow(rows), 4L)
  expect_equal(sort(unique(c(rows$phenotype_name_1, rows$phenotype_name_2))),
               c("A", "B"))
})

test_that("define_effect_cov_matrix() honours trait_names for the residual route", {
  pop <- make_block_pop("d1_trait_names")
  on.exit(close_pop(pop))

  M <- matrix(c(1, .2, .2, 1), 2, 2)
  pop <- suppressMessages(
    define_effect_cov_matrix(pop, "residual", M, trait_names = c("A", "B")))
  expect_equal(nrow(resid_rows(pop)), 4L)
  expect_error(define_effect_cov_matrix(pop, "residual", M), "row names")
})


# ── Strata ──────────────────────────────────────────────────────────────────

test_that("strata: every stratum names the same phenotypes and one condition column", {
  pop <- make_block_pop("strata")
  on.exit(close_pop(pop))

  R_ab <- sym(c("A", "B"), c(1, .3, .3, 2))
  pop <- define_residual_cov(pop, c("A", "B"), R_ab)

  # {A}-only conditional stratum of an {A, B} block: D1 (N != U)
  expect_error(
    define_residual_cov(pop, "A", matrix(1, 1, 1, dimnames = list("A", "A")),
                        condition_column = "sex", condition_level = "M"),
    "omits: B")

  # {A, B} per level succeeds
  pop <- define_residual_cov(pop, c("A", "B"), 2 * R_ab,
                             condition_column = "sex", condition_level = "M")
  pop <- define_residual_cov(pop, c("A", "B"), 3 * R_ab,
                             condition_column = "sex", condition_level = "F")
  rows <- resid_rows(pop)
  expect_equal(nrow(rows), 12L)
  expect_equal(sum(is.na(rows$condition_level)), 4L)
  expect_equal(unique(rows$condition_table[!is.na(rows$condition_column)]), "ind_meta")

  # A second condition column on the same block
  expect_error(
    define_residual_cov(pop, c("A", "B"), R_ab,
                        condition_column = "line_name", condition_level = "A"),
    "at most one condition column")

  # Column without level, level without column
  expect_error(define_residual_cov(pop, c("A", "B"), R_ab, condition_column = "sex"),
               "both `condition_column` and `condition_level`")
  expect_error(define_residual_cov(pop, c("A", "B"), R_ab, condition_level = "M"),
               "both `condition_column` and `condition_level`")

  # Redeclaring one stratum leaves the others alone
  pop <- define_residual_cov(pop, c("A", "B"), 5 * R_ab,
                             condition_column = "sex", condition_level = "M")
  rows <- resid_rows(pop)
  expect_equal(nrow(rows), 12L)
  expect_equal(rows$cov_value[rows$phenotype_name_1 == "A" & rows$phenotype_name_2 == "A" &
                              !is.na(rows$condition_level) & rows$condition_level == "M"], 5)
  expect_equal(rows$cov_value[rows$phenotype_name_1 == "A" & rows$phenotype_name_2 == "A" &
                              !is.na(rows$condition_level) & rows$condition_level == "F"], 3)
  expect_equal(rows$cov_value[rows$phenotype_name_1 == "A" & rows$phenotype_name_2 == "A" &
                              is.na(rows$condition_level)], 1)
})

test_that("strata: conditional strata alone (no unconditional stratum) are a valid block", {
  pop <- make_block_pop("strata_cond_only")
  on.exit(close_pop(pop))

  R1 <- matrix(400, 1, 1, dimnames = list("A", "A"))
  pop <- define_residual_cov(pop, "A", R1, condition_column = "sex", condition_level = "M")
  pop <- define_residual_cov(pop, "A", 2 * R1, condition_column = "sex", condition_level = "F")
  rows <- resid_rows(pop)
  expect_equal(nrow(rows), 2L)
  expect_true(all(!is.na(rows$condition_level)))
})

test_that("strata: growing a block that already has several strata requires clearing it", {
  pop <- make_block_pop("strata_grow")
  on.exit(close_pop(pop))

  R_ab <- sym(c("A", "B"), c(1, .3, .3, 2))
  pop <- define_residual_cov(pop, c("A", "B"), R_ab)
  pop <- define_residual_cov(pop, c("A", "B"), R_ab,
                             condition_column = "sex", condition_level = "M")

  R_abc <- sym(c("A", "B", "C"), c(1, .3, 0, .3, 2, 0, 0, 0, 1))
  expect_error(define_residual_cov(pop, c("A", "B", "C"), R_abc),
               "every stratum of a block names the same phenotypes")
  expect_error(define_residual_cov(pop, c("A", "B", "C"), R_abc),
               'get_table\\("phenotype_var_comp"\\)')
  expect_equal(nrow(resid_rows(pop)), 8L)

  # The recipe from the message
  pop <- suppressMessages(
    pop |> get_table("phenotype_var_comp") |>
      dplyr::filter(effect_name == "residual", phenotype_name_1 %in% c("A", "B")) |>
      remove_rows())
  pop <- define_residual_cov(pop, c("A", "B", "C"), R_abc)
  pop <- define_residual_cov(pop, c("A", "B", "C"), R_abc,
                             condition_column = "sex", condition_level = "M")
  expect_equal(nrow(resid_rows(pop)), 18L)
})


# ── Writer rollback ─────────────────────────────────────────────────────────

test_that("a failure after the DELETE rolls the writer back", {
  pop <- make_block_pop("rollback")
  on.exit(close_pop(pop))

  pop <- define_residual_cov(pop, c("A", "B"), sym(c("A", "B"), c(1, .3, .3, 2)))
  before <- resid_rows(pop)

  testthat::local_mocked_bindings(
    next_int_id = function(...) stop("simulated failure after DELETE"),
    .package = "tidybreed")
  expect_error(
    define_residual_cov(pop, c("A", "B"), sym(c("A", "B"), c(9, 0, 0, 9))),
    "simulated failure")
  expect_identical(resid_rows(pop), before)

  # The connection is usable again (no dangling transaction)
  expect_equal(DBI::dbGetQuery(pop$db_conn, "SELECT 1 AS x")$x, 1L)
})


# ── D3: realization lock ────────────────────────────────────────────────────

test_that("D3 (named effect): a block with stored draws cannot be redefined until they are removed", {
  pop <- make_block_pop("d3_named", traits = c("A", "B"))
  on.exit(close_pop(pop))

  pop <- suppressMessages(define_phenotype(pop, "A", residual_var = 1))
  pop <- suppressMessages(define_phenotype(pop, "B", residual_var = 1))
  R_pen <- sym(c("A", "B"), c(1, .2, .2, 1))
  pop <- suppressMessages(define_effect_cov_matrix(pop, "pen", R_pen))
  pop <- suppressMessages(define_effect_random(pop, "A", "pen", source_column = "sex"))
  pop <- suppressMessages(define_effect_random(pop, "B", "pen", source_column = "sex"))
  pop <- suppressMessages(pop |> get_table("ind_meta") |> add_phenotype("A"))

  n_draws <- DBI::dbGetQuery(pop$db_conn,
    "SELECT COUNT(*) AS n FROM phenotype_random_effects WHERE effect_name = 'pen'")$n
  expect_gt(n_draws, 0L)

  err <- tryCatch(define_effect_cov_matrix(pop, "pen", 2 * R_pen),
                  error = function(e) conditionMessage(e))
  expect_match(err, "realized draw")
  expect_match(err, 'get_table\\("phenotype_random_effects"\\)')
  expect_match(err, "remove_rows\\(\\)")
  expect_match(err, "not removed by that call")
  expect_equal(effect_rows(pop, "pen")$cov_value[1L], 1)

  pop <- suppressMessages(
    pop |> get_table("phenotype_random_effects") |>
      dplyr::filter(effect_name == "pen", phenotype_name %in% c("A", "B")) |>
      remove_rows())
  pop <- suppressMessages(define_effect_cov_matrix(pop, "pen", 2 * R_pen))
  expect_equal(effect_rows(pop, "pen")$cov_value[1L], 2)
})

test_that("D3 (residual): the records add_phenotype() writes lock the block; remove_rows() on residual_value clears it", {
  pop <- make_block_pop("d3_resid", traits = c("A", "B"))
  on.exit(close_pop(pop))

  R <- sym(c("A", "B"), c(1, .3, .3, 2))
  pop <- define_residual_cov(pop, c("A", "B"), R)
  pop <- suppressMessages(define_phenotype(pop, "A"))
  pop <- suppressMessages(define_phenotype(pop, "B"))
  pop <- suppressMessages(pop |> get_table("ind_meta") |> add_phenotype("A"))
  expect_equal(DBI::dbGetQuery(pop$db_conn,
    "SELECT COUNT(*) AS n FROM ind_phenotype WHERE residual_value IS NOT NULL")$n, 12L)

  err <- tryCatch(define_residual_cov(pop, c("A", "B"), 2 * R),
                  error = function(e) conditionMessage(e))
  expect_match(err, "has 12 realized draws in ind_phenotype")
  expect_match(err, "!is.na\\(residual_value\\)")
  expect_match(err, "remove_rows\\(\\)")
  expect_equal(resid_rows(pop)$cov_value[1L], 1)

  # A user_values record carries no residual and does not lock
  pop <- suppressMessages(
    pop |> get_table("ind_phenotype") |>
      dplyr::filter(phenotype_name %in% c("A", "B"), !is.na(residual_value)) |>
      remove_rows())
  pop <- suppressMessages(pop |> get_table("ind_meta") |>
    add_phenotype("B", user_values = rep(1, 12)))
  pop <- define_residual_cov(pop, c("A", "B"), 3 * R)
  expect_equal(resid_rows(pop)$cov_value[1L], 3)
})


# ── D5: the diagonal writers ────────────────────────────────────────────────

test_that("D5: define_phenotype(residual_var = ) writes a singleton, overwrites an unrealized one, and errors inside a block", {
  pop <- make_block_pop("d5_pheno")
  on.exit(close_pop(pop))

  # 1. no block -> singleton
  pop <- suppressMessages(define_phenotype(pop, "C", residual_var = 7))
  rows <- resid_rows(pop, "AND phenotype_name_1 = 'C'")
  expect_equal(nrow(rows), 1L)
  expect_equal(rows$cov_value, 7)

  # unrealized singleton -> overwritten
  pop <- suppressMessages(define_phenotype(pop, "C", residual_var = 8, overwrite = TRUE))
  expect_equal(resid_rows(pop, "AND phenotype_name_1 = 'C'")$cov_value, 8)

  # 2. inside {A, B} -> error naming the block and define_residual_cov()
  pop <- define_residual_cov(pop, c("A", "B"), sym(c("A", "B"), c(1, .3, .3, 2)))
  err <- tryCatch(define_phenotype(pop, "A", mean = 50, residual_var = 5),
                  error = function(e) conditionMessage(e))
  expect_match(err, "define_phenotype\\(residual_var = \\)")
  expect_match(err, "block is \\{A, B\\}")
  expect_match(err, "define_residual_cov\\(pop, c\\(\"A\", \"B\"\\), R\\)")
  # ... and phenotype_meta was not written
  expect_equal(DBI::dbGetQuery(pop$db_conn,
    "SELECT COUNT(*) AS n FROM phenotype_meta WHERE phenotype_name = 'A'")$n, 0L)

  # Without residual_var the phenotype is defined and the block is untouched
  pop <- suppressMessages(define_phenotype(pop, "A", mean = 50))
  expect_equal(nrow(resid_rows(pop, "AND phenotype_name_1 IN ('A','B')")), 4L)

  # An overwrite that errors on residual_var keeps the old phenotype_meta row
  expect_error(define_phenotype(pop, "A", mean = 60, residual_var = 5, overwrite = TRUE))
  expect_equal(DBI::dbGetQuery(pop$db_conn,
    "SELECT mean FROM phenotype_meta WHERE phenotype_name = 'A'")$mean, 50)
})

test_that("D5: define_phenotype(overwrite = TRUE) without residual_var leaves phenotype_var_comp untouched", {
  pop <- make_block_pop("d5_overwrite")
  on.exit(close_pop(pop))

  pop <- define_residual_cov(pop, c("A", "B"), sym(c("A", "B"), c(1, .3, .3, 2)))
  pop <- suppressMessages(define_phenotype(pop, "A", mean = 1))
  pop <- suppressMessages(define_phenotype(pop, "B", mean = 1))
  before <- resid_rows(pop)
  pop <- suppressMessages(define_phenotype(pop, "A", mean = 2, overwrite = TRUE))
  expect_identical(resid_rows(pop), before)
  expect_equal(DBI::dbGetQuery(pop$db_conn,
    "SELECT mean FROM phenotype_meta WHERE phenotype_name = 'A'")$mean, 2)
})

test_that("D5: define_phenotype(residual_var = ) on a realized singleton errors with the remove_rows() recipe", {
  pop <- make_block_pop("d5_locked", traits = "A")
  on.exit(close_pop(pop))

  pop <- suppressMessages(define_phenotype(pop, "A", residual_var = 1))
  pop <- suppressMessages(pop |> get_table("ind_meta") |> add_phenotype("A"))

  err <- tryCatch(define_phenotype(pop, "A", residual_var = 2, overwrite = TRUE),
                  error = function(e) conditionMessage(e))
  expect_match(err, "realized draws in ind_phenotype")
  expect_match(err, "remove_rows\\(\\)")
  expect_equal(resid_rows(pop)$cov_value, 1)
})

test_that("D5: define_effect_random(variance = ) mirrors the residual rules", {
  pop <- make_block_pop("d5_random", traits = c("A", "B"))
  on.exit(close_pop(pop))

  pop <- suppressMessages(define_phenotype(pop, "A", residual_var = 1))
  pop <- suppressMessages(define_phenotype(pop, "B", residual_var = 1))

  # No block -> singleton written
  pop <- suppressMessages(
    define_effect_random(pop, "A", "pen", source_column = "sex", variance = 3))
  expect_equal(effect_rows(pop, "pen")$cov_value, 3)

  # Overwrite with a new variance on an unrealized singleton
  pop <- suppressMessages(
    define_effect_random(pop, "A", "pen", source_column = "sex", variance = 4,
                         overwrite = TRUE))
  expect_equal(effect_rows(pop, "pen")$cov_value, 4)

  # No variance and none stored -> error
  expect_error(define_effect_random(pop, "B", "herd", source_column = "sex"),
               "No variance found")

  # Multi-member block -> variance here is an error naming define_effect_cov_matrix()
  pop <- suppressMessages(
    define_effect_cov_matrix(pop, "pen", sym(c("A", "B"), c(4, .5, .5, 2))))
  err <- tryCatch(
    define_effect_random(pop, "B", "pen", source_column = "sex", variance = 9),
    error = function(e) conditionMessage(e))
  expect_match(err, "define_effect_random\\(variance = \\)")
  expect_match(err, 'define_effect_cov_matrix\\(pop, "pen", R\\)')
  expect_equal(nrow(effect_rows(pop, "pen")), 4L)
  # ... and no phenotype_effects row was written (transaction rolled back)
  expect_equal(DBI::dbGetQuery(pop$db_conn,
    "SELECT COUNT(*) AS n FROM phenotype_effects WHERE phenotype_name = 'B'")$n, 0L)

  # Without variance it joins the block
  pop <- suppressMessages(define_effect_random(pop, "B", "pen", source_column = "sex"))
  expect_equal(DBI::dbGetQuery(pop$db_conn,
    "SELECT COUNT(*) AS n FROM phenotype_effects WHERE effect_name = 'pen'")$n, 2L)
})

test_that("define_effect_random(overwrite = TRUE) discards that phenotype's draws and is one transaction", {
  pop <- make_block_pop("random_overwrite", traits = "A")
  on.exit(close_pop(pop))

  pop <- suppressMessages(define_phenotype(pop, "A", residual_var = 1))
  pop <- suppressMessages(
    define_effect_random(pop, "A", "pen", source_column = "sex", variance = 3))
  pop <- suppressMessages(pop |> get_table("ind_meta") |> add_phenotype("A"))
  expect_gt(DBI::dbGetQuery(pop$db_conn,
    "SELECT COUNT(*) AS n FROM phenotype_random_effects")$n, 0L)

  # Without overwrite: the existing-effect error, nothing changed
  expect_error(define_effect_random(pop, "A", "pen", source_column = "sex", variance = 5),
               "overwrite = TRUE")
  expect_equal(effect_rows(pop, "pen")$cov_value, 3)

  # With overwrite: draws gone, variance replaced
  pop <- suppressMessages(
    define_effect_random(pop, "A", "pen", source_column = "sex", variance = 5,
                         overwrite = TRUE))
  expect_equal(DBI::dbGetQuery(pop$db_conn,
    "SELECT COUNT(*) AS n FROM phenotype_random_effects")$n, 0L)
  expect_equal(effect_rows(pop, "pen")$cov_value, 5)

  # A failure inside the call restores the row that overwrite deleted
  testthat::local_mocked_bindings(
    next_int_id = function(...) stop("simulated failure"), .package = "tidybreed")
  expect_error(define_effect_random(pop, "A", "pen", source_column = "sex",
                                    variance = 6, overwrite = TRUE),
               "simulated failure")
  expect_equal(DBI::dbGetQuery(pop$db_conn,
    "SELECT COUNT(*) AS n FROM phenotype_effects WHERE effect_name = 'pen'")$n, 1L)
  expect_equal(effect_rows(pop, "pen")$cov_value, 5)
})


# ── §5.6: named-effect block compatibility ──────────────────────────────────

test_that("§5.6: a block of two or more requires normal, random, source-compatible members", {
  pop <- make_block_pop("s56", traits = c("A", "B"))
  on.exit(close_pop(pop))

  pop <- suppressMessages(define_phenotype(pop, "A", residual_var = 1))
  pop <- suppressMessages(define_phenotype(pop, "B", residual_var = 1))
  R_pen <- sym(c("A", "B"), c(1, .2, .2, 1))
  pop <- suppressMessages(define_effect_cov_matrix(pop, "pen", R_pen))
  pop <- suppressMessages(define_effect_random(pop, "A", "pen", source_column = "sex"))

  # define_effect_random() side
  expect_error(define_effect_random(pop, "B", "pen", source_column = "line_name"),
               "different grouping columns")
  expect_error(define_effect_random(pop, "B", "pen", source_column = "sex",
                                    source_table = "ind_phenotype"),
               "different grouping columns")
  expect_error(define_effect_random(pop, "B", "pen", source_column = "sex",
                                    distribution = "gamma"),
               'distribution = "normal"')
  expect_equal(DBI::dbGetQuery(pop$db_conn,
    "SELECT COUNT(*) AS n FROM phenotype_effects WHERE effect_name = 'pen'")$n, 1L)
  pop <- suppressMessages(define_effect_random(pop, "B", "pen", source_column = "sex"))

  # define_effect_cov_matrix() side: rows already disagree -> the block is refused
  pop <- suppressMessages(
    define_effect_random(pop, "A", "herd", source_column = "sex", variance = 1))
  pop <- suppressMessages(
    define_effect_random(pop, "B", "herd", source_column = "line_name", variance = 1))
  expect_error(define_effect_cov_matrix(pop, "herd", R_pen), "different grouping columns")
  expect_equal(nrow(effect_rows(pop, "herd")), 2L)   # two singletons, unchanged

  # A fixed effect sharing the name cannot be joined
  pop <- suppressMessages(define_effect_fixed_class(pop, "A", "batch", source_column = "sex",
                                                    levels = c(M = 1, F = 0)))
  pop <- suppressMessages(
    define_effect_random(pop, "B", "batch", source_column = "sex", variance = 1))
  expect_error(define_effect_cov_matrix(pop, "batch", R_pen), "can only join random effects")
})

test_that("§5.6: a 1 x 1 gamma effect is legal until something joins it to a second phenotype", {
  pop <- make_block_pop("s56_gamma", traits = c("A", "B"))
  on.exit(close_pop(pop))

  pop <- suppressMessages(define_phenotype(pop, "A", residual_var = 1))
  pop <- suppressMessages(define_phenotype(pop, "B", residual_var = 1))
  pop <- suppressMessages(
    define_effect_random(pop, "A", "litter", source_column = "sex", variance = 1,
                         distribution = "gamma"))
  pop <- suppressMessages(pop |> get_table("ind_meta") |> add_phenotype("A"))
  draws <- DBI::dbGetQuery(pop$db_conn,
    "SELECT draw_value FROM phenotype_random_effects WHERE effect_name = 'litter'")
  expect_equal(nrow(draws), 2L)
  expect_true(all(draws$draw_value > 0))   # gamma, not normal

  expect_error(
    define_effect_cov_matrix(pop, "litter", sym(c("A", "B"), c(1, 0, 0, 1))),
    'A uses "gamma"')
})


# ── D6: condition_change_action agreement ───────────────────────────────────

test_that("D6: condition_change_action must agree across a residual block at definition time", {
  pop <- make_block_pop("d6", traits = c("A", "B"))
  on.exit(close_pop(pop))

  # Block declared before the phenotypes
  pop <- define_residual_cov(pop, c("A", "B"), sym(c("A", "B"), c(1, .3, .3, 2)))
  pop <- suppressMessages(define_phenotype(pop, "A"))
  err <- tryCatch(define_phenotype(pop, "B", condition_change_action = "independent"),
                  error = function(e) conditionMessage(e))
  expect_match(err, "A = 'error', B = 'independent'")
  expect_equal(DBI::dbGetQuery(pop$db_conn,
    "SELECT COUNT(*) AS n FROM phenotype_meta WHERE phenotype_name = 'B'")$n, 0L)

  # Fix A first, then B goes in
  pop <- suppressMessages(
    define_phenotype(pop, "A", condition_change_action = "independent", overwrite = TRUE))
  pop <- suppressMessages(
    define_phenotype(pop, "B", condition_change_action = "independent"))

  # Flipping one member back is refused
  expect_error(define_phenotype(pop, "A", condition_change_action = "error",
                                overwrite = TRUE),
               "must agree")

  # define_residual_cov() checks the block members that exist
  pop2 <- make_block_pop("d6_writer", traits = c("A", "B"))
  on.exit(close_pop(pop2), add = TRUE)
  pop2 <- suppressMessages(define_phenotype(pop2, "A"))
  pop2 <- suppressMessages(
    define_phenotype(pop2, "B", condition_change_action = "independent"))
  expect_error(define_residual_cov(pop2, c("A", "B"), sym(c("A", "B"), c(1, 0, 0, 1))),
               "must agree")
  expect_equal(nrow(resid_rows(pop2)), 0L)
  # Singletons never need agreement
  pop2 <- define_residual_cov(pop2, "A", matrix(1, 1, 1, dimnames = list("A", "A")))
  pop2 <- define_residual_cov(pop2, "B", matrix(1, 1, 1, dimnames = list("B", "B")))
  expect_equal(nrow(resid_rows(pop2)), 2L)
})
