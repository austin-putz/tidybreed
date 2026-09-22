# D7 — what a failed add_phenotype() call leaves behind
# (plans/sample_correlated_effects.md D7, §7 "Reproducibility and integrity"
# items 3 and 5; Phase 7).
#
# The contract: the database is atomic, the RNG is not. An error anywhere in
# the call — a Stage-1 rejection, a Stage-2 error after some draws, or a
# failed Stage-3 write — leaves `ind_phenotype` and `phenotype_random_effects`
# exactly as they were, while `.Random.seed` stays advanced by exactly the
# draws made before the error. Nothing in the package restores the seed, so a
# retry after an error draws different values (re-seed to reproduce).
#
# Exact assertions replay the resolver's contract (`n * m` standard normals
# per block, entities in sorted order, `sqrt(v) * z` for a 1 x 1 block) and
# the Stage-2 order: named effects (byte-sorted effect, blocks by first
# member, levels sorted) before residuals (blocks by first member).

# Twelve founders A_1..A_12; one continuous repeatable phenotype per trait,
# residual variance declared by the test through `...`. TBVs are written up
# front so the `add_tbv()` upsert inside `.ap_plan()` rewrites the same values
# and a whole database snapshot can be compared before and after a failed call.
make_d7_pop <- function(pop_name, traits = "A", seed = 11, ...) {
  set.seed(seed)
  pop <- suppressMessages(make_test_pop(pop_name, n_loci = 60, n_chr = 1,
                                        n_males = 6, n_females = 6))
  for (t in traits) {
    pop <- suppressMessages(define_trait(pop, t, target_add_var = 1))
    pop <- suppressMessages(
      pop |> get_table("genome_meta") |> define_additive_effects(t))
    pop <- suppressMessages(define_phenotype(
      pop, t, type = "continuous", mean = 10, repeatable = TRUE, ...))
  }
  suppressMessages(pop |> get_table("ind_meta") |> add_tbv(trait_name = traits))
}

set_col <- function(pop, col, value, ids = NULL) {
  tbl <- get_table(pop, "ind_meta")
  if (!is.null(ids)) tbl <- dplyr::filter(tbl, .data$id_ind %in% !!ids)
  args <- stats::setNames(list(value), col)
  suppressWarnings(suppressMessages(do.call(mutate_table, c(list(tbl), args))))
}

add_ph <- function(pop, phenos, ...) {
  suppressMessages(pop |> get_table("ind_meta") |> add_phenotype(phenos, ...))
}

# Every base table, each sorted by all its columns, so two snapshots are
# identical iff the database holds the same rows.
db_snapshot <- function(pop) {
  conn <- pop$db_conn
  tbls <- sort(DBI::dbGetQuery(conn, paste0(
    "SELECT table_name FROM information_schema.tables ",
    "WHERE table_type = 'BASE TABLE' AND table_schema = 'main'"))$table_name)
  stats::setNames(lapply(tbls, function(t)
    DBI::dbGetQuery(conn, paste0("SELECT * FROM \"", t, "\" ORDER BY ALL"))),
    tbls)
}

# `ind_tbv` is the one table a failed call still touches: `.ap_plan()` runs
# the add_tbv() upsert, which rewrites the same values. Those values are
# bit-identical, so this compares every table the same way -- row for row,
# `expect_identical()`. (Phase 7 had to give `ind_tbv` a tolerance, because
# the evaluator's parallel `SUM()` re-ordered the summation; Phase 8's exact
# accumulator removed the reason. See test-genome-effects-determinism.R.)
expect_db_unchanged <- function(pop, before) {
  expect_identical(db_snapshot(pop), before)
}

resid_of <- function(pop, t, pheno_number = 1L) {
  DBI::dbGetQuery(pop$db_conn, paste0(
    "SELECT residual_value FROM ind_phenotype WHERE phenotype_name = '", t,
    "' AND pheno_number = ", pheno_number, " ORDER BY id_ind"))$residual_value
}

draws_of <- function(pop, t, eff) {
  r <- DBI::dbGetQuery(pop$db_conn, paste0(
    "SELECT level, draw_value FROM phenotype_random_effects ",
    "WHERE phenotype_name = '", t, "' AND effect_name = '", eff, "' ORDER BY level"))
  stats::setNames(r$draw_value, r$level)
}

# The RNG state after `expr`, starting from `seed`.
state_after <- function(seed, expr) {
  set.seed(seed)
  force(expr)
  .Random.seed
}


# ── Stage 3: a failed write ─────────────────────────────────────────────────

test_that("D7 — a failed Stage-3 write leaves the database unchanged and the RNG advanced by exactly the Stage-2 draws; a retry draws new values", {
  pop <- make_d7_pop("d7_stage3", residual_var = 1)
  on.exit(close_pop(pop))
  n <- 12L
  pop <- set_col(pop, "pen", "P1")
  pop <- suppressMessages(define_effect_random(
    pop, "A", "pen", source_column = "pen", variance = 4))
  before <- db_snapshot(pop)

  # A reserved extra column is rejected inside the commit, after the
  # phenotype_random_effects INSERT has already run in the same transaction.
  s_fail <- state_after(41, expect_error(add_ph(pop, "A", pheno_value = 1),
                                         "reserved"))

  # Database: unchanged, every table. RNG: one pen level plus n residuals.
  expect_db_unchanged(pop, before)
  expect_identical(s_fail, state_after(41, stats::rnorm(1L + n)))

  # The schema is part of "unchanged": an extra column that was added by
  # ALTER TABLE earlier in the same transaction (prepare_extra_cols() runs
  # per field, so 'farm_x' is added before 'pheno_value' is refused) rolls
  # back with it.
  cols <- DBI::dbListFields(pop$db_conn, "ind_phenotype")
  expect_error(suppressMessages(add_ph(pop, "A", farm_x = "F1", pheno_value = 1)),
               "reserved")
  expect_identical(DBI::dbListFields(pop$db_conn, "ind_phenotype"), cols)
  expect_db_unchanged(pop, before)

  # A retry without re-seeding continues the stream — it does not reproduce
  # the failed call's draws.
  set.seed(41)
  z <- stats::rnorm(2L * (1L + n))
  set.seed(41)
  expect_error(add_ph(pop, "A", pheno_value = 1), "reserved")
  add_ph(pop, "A")
  expect_equal(unname(draws_of(pop, "A", "pen")), 2 * z[1L + n + 1L])
  expect_equal(resid_of(pop, "A"), z[(1L + n + 2L):(2L * (1L + n))])
  expect_false(isTRUE(all.equal(unname(draws_of(pop, "A", "pen")), 2 * z[1L])))

  # Re-seeding is how a retry reproduces the failed call
  pop2 <- make_d7_pop("d7_stage3_b", residual_var = 1)
  on.exit(close_pop(pop2), add = TRUE)
  pop2 <- set_col(pop2, "pen", "P1")
  pop2 <- suppressMessages(define_effect_random(
    pop2, "A", "pen", source_column = "pen", variance = 4))
  expect_error(add_ph(pop2, "A", seed = 41, pheno_value = 1), "reserved")
  add_ph(pop2, "A", seed = 41)
  expect_equal(unname(draws_of(pop2, "A", "pen")), 2 * z[1L])
  expect_equal(resid_of(pop2, "A"), z[2L:(1L + n)])
})


# ── Stage 2: an error after some draws ──────────────────────────────────────

test_that("D7 — a Stage-2 error in the named-effect adapter keeps the draws already made and writes nothing", {
  pop <- make_d7_pop("d7_stage2_named", residual_var = 1)
  on.exit(close_pop(pop))
  pop <- set_col(pop, "herd", "H1")
  pop <- set_col(pop, "herd", "H2", ids = c("A_2", "A_4", "A_6"))
  pop <- set_col(pop, "pen", "P1")
  pop <- suppressMessages(define_effect_random(
    pop, "A", "herd", source_column = "herd", variance = 1))
  pop <- suppressMessages(define_effect_random(
    pop, "A", "pen", source_column = "pen", variance = 4))
  # 'herd' < 'pen' in byte order: herd's two levels are drawn, then pen's
  # missing variance (rows removed by hand) is found.
  DBI::dbExecute(pop$db_conn,
    "DELETE FROM phenotype_var_comp WHERE effect_name = 'pen'")
  before <- db_snapshot(pop)

  s_fail <- state_after(42, expect_error(add_ph(pop, "A"),
    "No variance stored for random effect 'pen' / phenotype 'A'"))
  expect_db_unchanged(pop, before)
  expect_identical(s_fail, state_after(42, stats::rnorm(2L)))
})

test_that("D7 — a Stage-2 error in the residual adapter (D2) discards every earlier block's draws from the write but not from the stream", {
  R_BC <- matrix(c(1, .5, .5, 1), 2, 2, dimnames = list(c("B", "C"), c("B", "C")))
  pop <- make_d7_pop("d7_stage2_resid", traits = c("A", "B", "C"))
  on.exit(close_pop(pop))
  n <- 12L
  pop <- suppressMessages(define_residual_cov(pop, "A", matrix(1, dimnames = list("A", "A"))))
  pop <- set_col(pop, "farm", "F1")
  pop <- suppressMessages(define_residual_cov(
    pop, c("B", "C"), R_BC, condition_column = "farm", condition_level = "F1"))
  pop <- suppressMessages(define_residual_cov(
    pop, c("B", "C"), 4 * R_BC, condition_column = "farm", condition_level = "F2"))
  pop <- set_col(pop, "pen", "P1")
  pop <- suppressMessages(define_effect_random(
    pop, "A", "pen", source_column = "pen", variance = 4))

  # B stored under F1 for everyone; then two animals move farm, so C's draw
  # for them would condition on a B residual from another stratum (D2).
  add_ph(pop, "B", seed = 43)
  pop <- set_col(pop, "farm", "F2", ids = c("A_3", "A_5"))
  before <- db_snapshot(pop)

  # Stage-2 order: pen (1 new level), residual block {A} (n), then {B, C}
  # errors before its first draw. A's fully resolved records are never written.
  s_fail <- state_after(44, expect_error(add_ph(pop, c("A", "C")),
    "2 record\\(s\\) in block \\{B, C\\}.*stored under 'F1', now 'F2'"))
  expect_db_unchanged(pop, before)
  expect_identical(s_fail, state_after(44, stats::rnorm(1L + n)))
  expect_equal(nrow(before$ind_phenotype), n)
  expect_equal(nrow(before$phenotype_random_effects), 0L)
})


test_that("D7 — the Stage-1 add_tbv() upsert is the one write a failed call leaves behind", {
  # The contract is "no phenotype record and no draw", not "no write at all":
  # Stage 1 materializes the TBVs the plan reads, and they stay. They cost no
  # RNG and the retry rewrites them, so nothing stochastic survives.
  pop <- make_d7_pop("d7_tbv", residual_var = 1)
  on.exit(close_pop(pop))
  DBI::dbExecute(pop$db_conn, "DELETE FROM ind_tbv")
  expect_equal(nrow(dplyr::collect(get_table(pop, "ind_tbv"))), 0L)

  expect_error(add_ph(pop, "A", pheno_value = 1), "reserved")
  expect_equal(nrow(dplyr::collect(get_table(pop, "ind_tbv"))), 12L)
  expect_equal(nrow(dplyr::collect(get_table(pop, "ind_phenotype"))), 0L)
  expect_equal(nrow(dplyr::collect(get_table(pop, "phenotype_random_effects"))), 0L)
})


# ── Rejections before any draw ──────────────────────────────────────────────

test_that("D7 — a call rejected before its first draw leaves the RNG untouched", {
  pop <- make_d7_pop("d7_reject", traits = c("A", "B"))
  on.exit(close_pop(pop))
  # A singular residual block: A has zero variance
  pop <- suppressMessages(define_residual_cov(
    pop, c("A", "B"), matrix(c(0, 0, 0, 1), 2, 2,
                             dimnames = list(c("A", "B"), c("A", "B")))))
  before <- db_snapshot(pop)

  # Stage 1: an unknown phenotype
  s <- state_after(45, expect_error(add_ph(pop, "Z"), "not found in phenotype_meta"))
  expect_identical(s, state_after(45, NULL))

  # Stage 2, before the adapters: user_residual is validated against the plan
  s <- state_after(46, expect_error(add_ph(pop, "A", user_residual = c(1, 2)),
                                    "user_residual length"))
  expect_identical(s, state_after(46, NULL))

  # Stage 2, inside the resolver: a fixed value outside the support of the
  # singular block is refused before B's first normal
  s <- state_after(47, expect_error(
    add_ph(pop, c("A", "B"), user_residual = list(A = rep(1, 12))),
    "outside the support"))
  expect_identical(s, state_after(47, NULL))

  expect_db_unchanged(pop, before)
})
