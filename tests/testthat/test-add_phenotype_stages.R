# The three stages of add_phenotype(): PLAN (no RNG, no writes), RESOLVE
# (RNG, no writes), COMMIT (writes, no RNG). See ?add_phenotype_stages and
# plans/sample_correlated_effects.md §5.5 / §7 "Record planning (Stage 1)".

sym <- function(names, values) {
  n <- length(names)
  matrix(values, n, n, byrow = TRUE, dimnames = list(names, names))
}

# Twelve founders A_1..A_12 (lexicographic order A_1, A_10, A_11, A_12, A_2, ...),
# traits with QTL, one continuous repeatable phenotype per trait. Seeded so two
# pops built with the same seed have identical genomes and TBVs.
make_stage_pop <- function(pop_name, traits = "A", seed = 11) {
  set.seed(seed)
  pop <- suppressMessages(make_test_pop(pop_name, n_loci = 60, n_chr = 1,
                                        n_males = 6, n_females = 6))
  for (t in traits) {
    pop <- suppressMessages(define_trait(pop, t, target_add_var = 1))
    pop <- suppressMessages(
      pop |> get_table("genome_meta") |> define_additive_effects(t))
    pop <- suppressMessages(define_phenotype(
      pop, t, type = "continuous", mean = 10, residual_var = 1,
      repeatable = TRUE))
  }
  pop
}

set_col <- function(pop, col, value, ids = NULL) {
  tbl <- get_table(pop, "ind_meta")
  if (!is.null(ids)) tbl <- dplyr::filter(tbl, .data$id_ind %in% !!ids)
  args <- stats::setNames(list(value), col)
  suppressWarnings(suppressMessages(do.call(mutate_table, c(list(tbl), args))))
}

phen <- function(pop, ...) {
  DBI::dbGetQuery(pop$db_conn, paste0(
    "SELECT * FROM ind_phenotype ", ..., " ORDER BY id_phenotype"))
}

re_rows <- function(pop) {
  DBI::dbGetQuery(pop$db_conn,
    "SELECT * FROM phenotype_random_effects ORDER BY phenotype_name, effect_name, level")
}

all_ids <- function(pop) sort(DBI::dbGetQuery(pop$db_conn,
                                              "SELECT id_ind FROM ind_meta")$id_ind)

# The RNG state after `expr`, starting from `seed`.
state_after <- function(seed, expr) {
  set.seed(seed)
  force(expr)
  .Random.seed
}


# ── Stage 1: exclusions consume no RNG and leave no stochastic state ────────

test_that("null_class_action = 'skip': seeded output is identical with or without the excluded animal, and its random-effect level is never drawn", {
  build <- function(name) {
    pop <- make_stage_pop(name)
    pop <- set_col(pop, "farm", "F1")
    pop <- set_col(pop, "farm", NA_character_, ids = "A_3")
    pop <- set_col(pop, "pen", "P1")
    pop <- set_col(pop, "pen", "P9", ids = "A_3")   # a level only A_3 touches
    pop <- suppressMessages(define_effect_fixed_class(
      pop, "A", "farm", source_column = "farm", levels = c(F1 = 2)))
    pop <- suppressMessages(define_effect_random(
      pop, "A", "pen", source_column = "pen", variance = 2))
    pop
  }
  pop_with    <- build("st_skip_with")
  pop_without <- build("st_skip_without")
  on.exit({ close_pop(pop_with); close_pop(pop_without) })

  set.seed(1)
  expect_warning(
    suppressMessages(pop_with |> get_table("ind_meta") |> add_phenotype("A")),
    "null_class_action = 'skip'")
  set.seed(1)
  suppressMessages(pop_without |> get_table("ind_meta") |>
    dplyr::filter(.data$id_ind != "A_3") |> add_phenotype("A"))

  a <- phen(pop_with); b <- phen(pop_without)
  expect_false("A_3" %in% a$id_ind)
  expect_identical(a$id_ind, b$id_ind)
  expect_equal(a$pheno_value, b$pheno_value)

  # The skipped animal's pen level was never drawn: no stochastic state
  expect_identical(re_rows(pop_with)$level, "P1")
  expect_identical(re_rows(pop_with)$draw_value, re_rows(pop_without)$draw_value)
})

test_that("composite and formula_tbv exclusions consume no RNG either", {
  build <- function(name, how) {
    pop <- make_stage_pop(name, traits = c("A_direct", "A_social"))
    pop <- set_col(pop, "pen", "P1")
    pop <- set_col(pop, "pen", "P2", ids = c("A_2", "A_4", "A_6"))
    pop <- set_col(pop, "pen", NA_character_, ids = "A_5")   # no group → excluded
    if (how == "components") {
      pop <- suppressMessages(define_phenotype(
        pop, "A_obs", type = "continuous", mean = 10, residual_var = 1,
        components = tibble::tribble(
          ~source_trait_name, ~contributor_type, ~group_column,
          "A_direct",         "self",            NA_character_,
          "A_social",         "group",           "pen")))
    } else {
      pop <- suppressMessages(define_phenotype(
        pop, "A_obs", type = "continuous", mean = 10, residual_var = 1,
        formula_tbv = "A_direct + group_sum(A_social, pen)"))
    }
    pop
  }
  for (how in c("components", "formula_tbv")) {
    pop_with    <- build(paste0("st_excl_with_", how), how)
    pop_without <- build(paste0("st_excl_without_", how), how)

    set.seed(3)
    expect_warning(
      suppressMessages(pop_with |> get_table("ind_meta") |> add_phenotype("A_obs")),
      "missing components")
    set.seed(3)
    suppressMessages(pop_without |> get_table("ind_meta") |>
      dplyr::filter(.data$id_ind != "A_5") |> add_phenotype("A_obs"))

    a <- phen(pop_with); b <- phen(pop_without)
    expect_false("A_5" %in% a$id_ind)
    expect_identical(a$id_ind, b$id_ind)
    expect_equal(a$pheno_value, b$pheno_value)
    close_pop(pop_with); close_pop(pop_without)
  }
})

test_that("the repeatable guard excludes before any draw", {
  pop1 <- make_stage_pop("st_rep_1")
  pop2 <- make_stage_pop("st_rep_2")
  on.exit({ close_pop(pop1); close_pop(pop2) })
  pop1 <- suppressMessages(define_phenotype(pop1, "A", mean = 10, residual_var = 1,
                                            repeatable = FALSE, overwrite = TRUE))
  pop2 <- suppressMessages(define_phenotype(pop2, "A", mean = 10, residual_var = 1,
                                            repeatable = FALSE, overwrite = TRUE))

  set.seed(5)
  suppressMessages(pop1 |> get_table("ind_meta") |>
    dplyr::filter(.data$id_ind %in% c("A_1", "A_2")) |> add_phenotype("A"))
  set.seed(6)
  expect_warning(
    suppressMessages(pop1 |> get_table("ind_meta") |> add_phenotype("A")),
    "not repeatable")
  set.seed(6)
  suppressMessages(pop2 |> get_table("ind_meta") |>
    dplyr::filter(!.data$id_ind %in% c("A_1", "A_2")) |> add_phenotype("A"))

  a <- phen(pop1, "WHERE id_ind NOT IN ('A_1', 'A_2')")
  b <- phen(pop2)
  expect_identical(a$id_ind, b$id_ind)
  expect_equal(a$pheno_value, b$pheno_value)
})


# ── Stage 1: pheno_number and record order ──────────────────────────────────

test_that("pheno_number assigned in Stage 1 is what Stage 3 writes", {
  pop <- make_stage_pop("st_pn")
  on.exit(close_pop(pop))
  suppressMessages(pop |> get_table("ind_meta") |>
    dplyr::filter(.data$id_ind %in% c("A_1", "A_2", "A_3")) |> add_phenotype("A"))
  suppressMessages(pop |> get_table("ind_meta") |>
    dplyr::filter(.data$id_ind %in% c("A_1")) |> add_phenotype("A"))

  tbl  <- pop |> get_table("ind_meta") |>
    dplyr::filter(.data$id_ind %in% c("A_1", "A_2", "A_5"))
  plan <- suppressMessages(.ap_plan(tbl, "A"))
  e    <- plan$entries[["A"]]
  expect_identical(e$path, "model")
  expect_identical(e$id_ind, c("A_1", "A_2", "A_5"))
  expect_identical(e$pheno_number, c(3L, 2L, 1L))

  suppressMessages(add_phenotype(tbl, "A"))
  written <- phen(pop, "WHERE id_phenotype > 4")
  expect_identical(written$id_ind, e$id_ind)
  expect_identical(written$pheno_number, e$pheno_number)
})

test_that("records are planned and written in id_ind order, not physical row order", {
  pop <- make_stage_pop("st_order")
  on.exit(close_pop(pop))
  suppressMessages(pop |> get_table("ind_meta") |> add_phenotype("A"))
  ids <- all_ids(pop)
  expect_identical(phen(pop)$id_ind, ids)   # A_1, A_10, A_11, A_12, A_2, ...
  expect_identical(phen(pop)$id_phenotype, seq_along(ids))
})

# §7 "Reproducibility and integrity" item 4: shuffle ind_meta's physical order
# between two seeded runs and compare. (The test above asserts the *output*
# order; this one asserts the draws themselves do not depend on the input
# order, which is what makes the id_ind sort load-bearing rather than
# cosmetic. Both adapters run: a named effect with three levels and a
# two-phenotype residual block.)
test_that("seeded output is identical when ind_meta's physical row order is reversed", {
  build <- function(name, reverse) {
    pop <- make_stage_pop(name, traits = c("A", "B"))
    pop <- suppressMessages(define_residual_cov(
      pop, c("A", "B"), sym(c("A", "B"), c(2, .8, .8, 1))))
    pop <- set_col(pop, "pen", "P1")
    pop <- set_col(pop, "pen", "P2", ids = c("A_2", "A_4", "A_7"))
    pop <- set_col(pop, "pen", "P3", ids = c("A_11"))
    for (t in c("A", "B")) {
      pop <- suppressMessages(define_effect_random(
        pop, t, "pen", source_column = "pen", variance = 4))
    }
    if (reverse) {
      # CREATE TABLE AS keeps the column types; nothing declares a SQL
      # foreign key to ind_meta, so the table can be rebuilt in place.
      DBI::dbExecute(pop$db_conn,
        "CREATE TABLE __shuf AS SELECT * FROM ind_meta ORDER BY id_ind DESC")
      DBI::dbExecute(pop$db_conn, "DROP TABLE ind_meta")
      DBI::dbExecute(pop$db_conn, "ALTER TABLE __shuf RENAME TO ind_meta")
    }
    suppressMessages(pop |> get_table("ind_meta") |>
                       add_phenotype(c("A", "B"), seed = 77))
    pop
  }
  ordered  <- build("st_phys_a", FALSE)
  on.exit(close_pop(ordered))
  shuffled <- build("st_phys_b", TRUE)
  on.exit(close_pop(shuffled), add = TRUE)

  # The physical order really did change, and the records still come out sorted
  expect_false(identical(
    DBI::dbGetQuery(ordered$db_conn,  "SELECT id_ind FROM ind_meta")$id_ind,
    DBI::dbGetQuery(shuffled$db_conn, "SELECT id_ind FROM ind_meta")$id_ind))

  cols <- c("id_phenotype", "id_ind", "phenotype_name", "pheno_value",
            "pheno_number", "residual_value")
  expect_equal(phen(ordered)[cols], phen(shuffled)[cols])
  expect_equal(re_rows(ordered), re_rows(shuffled))
})

test_that("user_residual is positional over the planned (id_ind-ordered) records", {
  pop <- make_stage_pop("st_ures")
  on.exit(close_pop(pop))
  ids <- all_ids(pop)
  resid <- seq_along(ids) / 10
  suppressMessages(pop |> get_table("ind_meta") |>
    add_phenotype("A", user_residual = resid))
  ph  <- phen(pop)
  tbv <- DBI::dbGetQuery(pop$db_conn,
    "SELECT id_ind, tbv_value FROM ind_tbv WHERE trait_name = 'A'")
  tbv <- stats::setNames(tbv$tbv_value, tbv$id_ind)[ph$id_ind]
  expect_equal(ph$pheno_value - 10 - unname(tbv), resid)
  expect_error(
    pop |> get_table("ind_meta") |> add_phenotype("A", user_residual = 1:3),
    "must equal 12")
})


# ── Stage 2 / 3: RNG accounting and atomicity ───────────────────────────────

test_that("a call consumes exactly its draws: n residuals plus one normal per new random-effect level, nothing else", {
  pop <- make_stage_pop("st_rng")
  on.exit(close_pop(pop))
  n <- 12L

  # No random effects: exactly n normals (add_tbv, planning, commit are silent)
  s_call <- state_after(21, suppressMessages(
    pop |> get_table("ind_meta") |> add_phenotype("A")))
  expect_identical(s_call, state_after(21, stats::rnorm(n)))

  # A random effect with three new levels: 3 + n normals, levels drawn first
  pop <- set_col(pop, "pen", "P1")
  pop <- set_col(pop, "pen", "P2", ids = c("A_2", "A_4"))
  pop <- set_col(pop, "pen", "P3", ids = c("A_6"))
  pop <- suppressMessages(define_effect_random(
    pop, "A", "pen", source_column = "pen", variance = 4))
  s_call <- state_after(22, suppressMessages(
    pop |> get_table("ind_meta") |> add_phenotype("A")))
  expect_identical(s_call, state_after(22, stats::rnorm(3 + n)))
  set.seed(22)
  expect_equal(re_rows(pop)$draw_value, stats::rnorm(3, sd = 2))  # sorted levels

  # Levels already drawn are reused: only the n residuals are consumed
  s_call <- state_after(23, suppressMessages(
    pop |> get_table("ind_meta") |> add_phenotype("A")))
  expect_identical(s_call, state_after(23, stats::rnorm(n)))
  expect_equal(nrow(re_rows(pop)), 3L)

  # user_values skips the model: the call is RNG-neutral
  s_call <- state_after(24, suppressMessages(
    pop |> get_table("ind_meta") |> add_phenotype("A", user_values = rep(1, n))))
  expect_identical(s_call, state_after(24, NULL))
})

# (The full D7 contract — every table unchanged, seed advanced by exactly the
# draws made, for failures in every stage — is test-add_phenotype_failure_contract.R.)
test_that("Stage 3 is atomic: a failed record write leaves no random-effect draws behind", {
  pop <- make_stage_pop("st_atomic")
  on.exit(close_pop(pop))
  pop <- set_col(pop, "pen", "P1")
  pop <- suppressMessages(define_effect_random(
    pop, "A", "pen", source_column = "pen", variance = 4))

  # A reserved extra column is rejected inside the commit, after the
  # phenotype_random_effects INSERT has already run in the same transaction.
  expect_error(
    suppressMessages(pop |> get_table("ind_meta") |> add_phenotype("A", pheno_value = 1)),
    "reserved")
  expect_equal(nrow(phen(pop)), 0L)
  expect_equal(nrow(re_rows(pop)), 0L)

  # The connection is usable and a clean call succeeds afterwards
  suppressMessages(pop |> get_table("ind_meta") |> add_phenotype("A"))
  expect_equal(nrow(phen(pop)), 12L)
  expect_equal(nrow(re_rows(pop)), 1L)
})


# ── Derived formulas over records of the same call ──────────────────────────

test_that("a derived phenotype consumes a feeder planned in the same call", {
  pop <- make_stage_pop("st_derived")
  on.exit(close_pop(pop))
  pop <- suppressMessages(define_phenotype(pop, "D", type = "derived_formula",
                                           formula = "2 * A"))
  suppressMessages(pop |> get_table("ind_meta") |> add_phenotype(c("D", "A")))
  ph <- phen(pop)
  a  <- ph[ph$phenotype_name == "A", ]
  d  <- ph[ph$phenotype_name == "D", ]
  expect_equal(nrow(d), 12L)
  expect_true(all(a$id_phenotype < min(d$id_phenotype)))   # topological order
  expect_equal(d$pheno_value, 2 * a$pheno_value[match(d$id_ind, a$id_ind)])
})


# ── Stage 1: residual stratum lookup contract ───────────────────────────────

test_that("the condition table must have exactly one row per planned individual", {
  pop <- make_stage_pop("st_stratum", traits = c("A", "B"))
  on.exit(close_pop(pop))
  conn <- pop$db_conn
  R <- sym(c("A", "B"), c(1, .3, .3, 1))
  pop <- suppressMessages(define_residual_cov(pop, c("A", "B"), R))
  pop <- suppressMessages(define_residual_cov(pop, c("A", "B"), 2 * R,
    condition_column = "farm", condition_table = "ind_env", condition_level = "F1"))
  pop <- suppressMessages(define_residual_cov(pop, c("A", "B"), 3 * R,
    condition_column = "farm", condition_table = "ind_env", condition_level = "F2"))

  ids <- all_ids(pop)
  env <- data.frame(id_ind = c(ids, "A_1", "A_2"),
                    farm = c(rep(c("F1", "F2"), 6), "F2", "F1"),
                    stringsAsFactors = FALSE)
  DBI::dbExecute(conn, "CREATE TABLE ind_env (id_ind VARCHAR, farm VARCHAR)")
  duckdb::duckdb_register(conn, "__env", env)
  DBI::dbExecute(conn, "INSERT INTO ind_env SELECT * FROM __env")
  duckdb::duckdb_unregister(conn, "__env")

  expect_error(
    suppressMessages(pop |> get_table("ind_meta") |> add_phenotype(c("A", "B"))),
    "'ind_env'.*'farm'.*2 have several rows.*A_1, A_2")
  expect_equal(nrow(phen(pop)), 0L)

  DBI::dbExecute(conn, "DELETE FROM ind_env WHERE rowid >= 12")
  DBI::dbExecute(conn, "DELETE FROM ind_env WHERE id_ind = 'A_7'")
  expect_error(
    suppressMessages(pop |> get_table("ind_meta") |> add_phenotype(c("A", "B"))),
    "1 planned individual\\(s\\) have no row.*A_7")

  DBI::dbExecute(conn, "INSERT INTO ind_env VALUES ('A_7', 'F3')")
  expect_warning(
    suppressMessages(pop |> get_table("ind_meta") |> add_phenotype(c("A", "B"))),
    "1 planned record\\(s\\) of block \\{A, B\\} have a farm value matching no residual stratum \\('F3'\\).*A_7")
  expect_equal(nrow(phen(pop)), 24L)

  plan <- suppressMessages(.ap_plan(pop |> get_table("ind_meta"), c("A", "B")))
  e <- plan$entries[["A"]]
  expect_identical(e$condition_table, "ind_env")
  expect_identical(e$condition_column, "farm")
  expect_identical(e$id_ind, ids)
  expect_identical(e$condition_value,
                   ifelse(ids == "A_7", "F3", rep(c("F1", "F2"), 6)))
})


# ── user_values: named vectors are checked against the planned set ──────────

test_that("named user_values must name planned individuals, each once, and are written in id_ind order", {
  pop <- make_stage_pop("st_uv")
  on.exit(close_pop(pop))

  expect_error(
    pop |> get_table("ind_meta") |>
      add_phenotype("A", user_values = c(A_2 = 1, ghost = 2, other = 3)),
    "2 unknown \\(e.g. ghost, other\\)")
  expect_error(
    pop |> get_table("ind_meta") |>
      add_phenotype("A", user_values = c(A_2 = 1, A_2 = 2)),
    "Duplicate names")
  expect_equal(nrow(phen(pop)), 0L)

  suppressMessages(pop |> get_table("ind_meta") |>
    add_phenotype("A", user_values = c(A_3 = 30, A_1 = 10, A_2 = 20)))
  ph <- phen(pop)
  expect_identical(ph$id_ind, c("A_1", "A_2", "A_3"))
  expect_identical(ph$pheno_value, c(10, 20, 30))
})

test_that("a derived phenotype has no model terms: effects declared on it are ignored and draw nothing", {
  pop <- make_stage_pop("st_derived_fx")
  on.exit(close_pop(pop))
  pop <- set_col(pop, "pen", "P1")
  pop <- suppressMessages(define_phenotype(pop, "D", type = "derived_formula",
                                           formula = "2 * A"))
  pop <- suppressMessages(define_effect_random(
    pop, "D", "pen", source_column = "pen", variance = 4))

  s_call <- state_after(31, suppressMessages(
    pop |> get_table("ind_meta") |> add_phenotype(c("A", "D"))))
  expect_identical(s_call, state_after(31, stats::rnorm(12)))   # A's residuals only
  expect_equal(nrow(re_rows(pop)), 0L)
  ph <- phen(pop)
  a <- ph[ph$phenotype_name == "A", ]; d <- ph[ph$phenotype_name == "D", ]
  expect_equal(d$pheno_value, 2 * a$pheno_value[match(d$id_ind, a$id_ind)])
})
