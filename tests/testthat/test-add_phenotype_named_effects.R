# The named-effect adapter of add_phenotype() (Phase 6 of
# plans/sample_correlated_effects.md): one persistent draw per (effect,
# level), conditional on the coordinates the level has already realized for
# the block's other phenotypes, in this call or an earlier one. See §5.6,
# §5.8 and §7 "Named-effect blocks", "Defect 3".
#
# Exact tests replay the resolver's contract: `n * m` standard normals in
# entity (byte-sorted level) order and coordinate order within an entity,
# applied through chol(C) (draws = z %*% U). Named-effect draws precede the
# residual draws of a call.

sym <- function(names, values) {
  n <- length(names)
  matrix(values, n, n, byrow = TRUE, dimnames = list(names, names))
}

# Founders A_1..A_n with traits (QTL placed) and one continuous repeatable
# phenotype per trait with residual variance 1, plus a `pen` column.
make_ne_pop <- function(pop_name, traits = c("A", "B"), n_males = 6,
                        n_females = 6, seed = 11) {
  set.seed(seed)
  pop <- suppressMessages(make_test_pop(pop_name, n_loci = 60, n_chr = 1,
                                        n_males = n_males, n_females = n_females))
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

# Assign pens round-robin over sorted ids: A_1 -> pens[1], A_2 -> pens[2], ...
assign_pens <- function(pop, pens, ids = all_ids(pop)) {
  pen <- rep(pens, length.out = length(ids))
  for (p in unique(pen)) pop <- set_col(pop, "pen", p, ids = ids[pen == p])
  pop
}

all_ids <- function(pop) sort(DBI::dbGetQuery(pop$db_conn,
                                              "SELECT id_ind FROM ind_meta")$id_ind)

sex_of <- function(pop) {
  m <- DBI::dbGetQuery(pop$db_conn, "SELECT id_ind, sex FROM ind_meta")
  stats::setNames(m$sex, m$id_ind)
}


pen_of <- function(pop) {
  m <- DBI::dbGetQuery(pop$db_conn, "SELECT id_ind, pen FROM ind_meta ORDER BY id_ind")
  stats::setNames(as.character(m$pen), m$id_ind)
}

# draw_value of one (phenotype, effect), named by level, in byte-sorted order
draws_of <- function(pop, t, eff) {
  r <- DBI::dbGetQuery(pop$db_conn, paste0(
    "SELECT level, draw_value FROM phenotype_random_effects ",
    "WHERE phenotype_name = '", t, "' AND effect_name = '", eff, "'"))
  r <- r[order(r$level, method = "radix"), , drop = FALSE]
  stats::setNames(r$draw_value, r$level)
}

# pheno_value - mean - tbv - residual: the summed random-effect contribution
# of every record of one phenotype at one pheno_number, named by id_ind
contrib_of <- function(pop, t, pheno_number = 1L) {
  r <- DBI::dbGetQuery(pop$db_conn, paste0(
    "SELECT p.id_ind, p.pheno_value - 10 - v.tbv_value - p.residual_value AS c ",
    "FROM ind_phenotype AS p JOIN ind_tbv AS v ",
    "ON v.id_ind = p.id_ind AND v.trait_name = p.phenotype_name ",
    "WHERE p.phenotype_name = '", t, "' AND p.pheno_number = ", pheno_number,
    " ORDER BY p.id_ind"))
  stats::setNames(r$c, r$id_ind)
}

resid_of <- function(pop, t, pheno_number = 1L) {
  r <- DBI::dbGetQuery(pop$db_conn, paste0(
    "SELECT id_ind, residual_value FROM ind_phenotype WHERE phenotype_name = '",
    t, "' AND pheno_number = ", pheno_number, " ORDER BY id_ind"))
  stats::setNames(r$residual_value, r$id_ind)
}

state_after <- function(seed, expr) {
  set.seed(seed)
  force(expr)
  .Random.seed
}

add_ph <- function(pop, phenos, tbl = get_table(pop, "ind_meta"), ...) {
  suppressMessages(add_phenotype(tbl, phenos, ...))
}

set_pen_effect <- function(pop, phenos, R, source_column = "pen") {
  pop <- suppressMessages(define_effect_cov_matrix(pop, "pen", R))
  for (t in phenos) {
    pop <- suppressMessages(define_effect_random(
      pop, t, "pen", source_column = source_column))
  }
  pop
}

re_rows_all <- function(pop) DBI::dbGetQuery(pop$db_conn,
  "SELECT * FROM phenotype_random_effects ORDER BY phenotype_name, effect_name, level")

# Standard normals in resolver order for n entities x m coordinates.
z_mat <- function(n, m) matrix(stats::rnorm(n * m), n, m, byrow = TRUE)

# §5.8: pen on ADG (A) and BF (B), negatively correlated
R_pen <- sym(c("A", "B"), c(150, -12, -12, 4))
U_pen <- chol(R_pen)
PENS  <- c("P01", "P02", "P03", "P04")


# ── §5.8: the worked example, across calls ──────────────────────────────────

test_that("day 0 / day 100 / day 200: a level's draw is realized once, the second phenotype conditions on it, and a reuse draws nothing", {
  pop <- make_ne_pop("ne_seq", n_males = 6, n_females = 6)
  on.exit(close_pop(pop))
  pop <- assign_pens(pop, PENS)
  pop <- suppressMessages(define_effect_cov_matrix(pop, "pen", R_pen))
  pop <- suppressMessages(define_effect_random(pop, "A", "pen", source_column = "pen"))
  pop <- suppressMessages(define_effect_random(pop, "B", "pen", source_column = "pen"))
  n <- 12L

  # Day 0: A only. Four pen normals (sorted levels) then n residual normals.
  s <- state_after(1, pop <- add_ph(pop, "A", seed = 1))
  expect_identical(s, state_after(1, stats::rnorm(4 + n)))
  set.seed(1)
  z <- stats::rnorm(4)
  a <- draws_of(pop, "A", "pen")
  expect_identical(names(a), PENS)
  expect_equal(unname(a), z * sqrt(150))                    # marginal
  expect_equal(unname(resid_of(pop, "A")), stats::rnorm(n))  # residuals after
  # B stays latent: nothing drawn, nothing stored
  expect_length(draws_of(pop, "B", "pen"), 0L)
  # Every record carries its pen's draw
  expect_equal(unname(contrib_of(pop, "A")), unname(a[pen_of(pop)]))

  # Day 100: B only. Each pen conditions on its stored A.
  s <- state_after(2, pop <- add_ph(pop, "B", seed = 2))
  expect_identical(s, state_after(2, stats::rnorm(4 + n)))
  set.seed(2)
  z <- stats::rnorm(4)
  b <- draws_of(pop, "B", "pen")
  expect_identical(names(b), PENS)
  expect_equal(unname(b), (-12 / 150) * unname(a) + sqrt(4 - 144 / 150) * z)
  expect_equal(draws_of(pop, "A", "pen"), a)                 # A untouched
  expect_equal(unname(contrib_of(pop, "B")), unname(b[pen_of(pop)]))

  # Day 200: A again for the same pens. No pen draw; only the n residuals.
  s <- state_after(3, pop <- add_ph(pop, "A", seed = 3))
  expect_identical(s, state_after(3, stats::rnorm(n)))
  expect_equal(draws_of(pop, "A", "pen"), a)
  expect_equal(unname(contrib_of(pop, "A", 2L)), unname(a[pen_of(pop)]))
})

test_that("mixed patterns: stored and new levels in one call are one entity list, conditional where stored, marginal where not", {
  pop <- make_ne_pop("ne_mixed", n_males = 8, n_females = 8)
  on.exit(close_pop(pop))
  ids <- all_ids(pop)
  first <- ids[1:8]; second <- ids[9:16]
  pop <- assign_pens(pop, PENS, ids = first)
  pop <- set_pen_effect(pop, c("A", "B"), R_pen)

  # A on the first batch: P01..P04 get an A draw
  pop <- add_ph(pop, "A", tbl = get_table(pop, "ind_meta") |>
                  dplyr::filter(.data$id_ind %in% !!first), seed = 4)
  a <- draws_of(pop, "A", "pen")
  expect_identical(names(a), PENS)

  # B on everyone, with the second batch in new pens P05..P08
  new_pens <- c("P05", "P06", "P07", "P08")
  pop <- assign_pens(pop, new_pens, ids = second)
  s <- state_after(5, pop <- add_ph(pop, "B", seed = 5))
  expect_identical(s, state_after(5, stats::rnorm(8 + 16)))
  set.seed(5)
  z <- stats::rnorm(8)                                       # entity order
  b <- draws_of(pop, "B", "pen")
  expect_identical(names(b), c(PENS, new_pens))
  expect_equal(unname(b[PENS]), (-12 / 150) * unname(a) + sqrt(4 - 144 / 150) * z[1:4])
  expect_equal(unname(b[new_pens]), 2 * z[5:8])              # marginal sd 2
  expect_length(draws_of(pop, "A", "pen"), 4L)               # P05..P08 A latent
})

test_that("Defect 3 closed: both phenotypes in one call reuse a stored A exactly and draw B conditional on it; new levels draw jointly", {
  pop <- make_ne_pop("ne_joint", n_males = 6, n_females = 6)
  on.exit(close_pop(pop))
  ids <- all_ids(pop)
  pop <- assign_pens(pop, c("P01", "P02", "P03"), ids = ids[1:9])
  pop <- assign_pens(pop, "P04", ids = ids[10:12])
  pop <- set_pen_effect(pop, c("A", "B"), R_pen)

  pop <- add_ph(pop, "A", tbl = get_table(pop, "ind_meta") |>
                  dplyr::filter(.data$pen != "P04"), seed = 6)
  a <- draws_of(pop, "A", "pen")
  expect_identical(names(a), c("P01", "P02", "P03"))

  # Sample-set groups in byte order: {A, B} (P04) before {B} (P01..P03).
  s <- state_after(7, pop <- add_ph(pop, c("A", "B"), seed = 7))
  expect_identical(s, state_after(7, stats::rnorm(2 + 3 + 24)))
  set.seed(7)
  z_joint <- z_mat(1, 2)
  z_cond  <- stats::rnorm(3)
  joint   <- z_joint %*% U_pen
  a2 <- draws_of(pop, "A", "pen")
  b2 <- draws_of(pop, "B", "pen")
  expect_equal(a2[c("P01", "P02", "P03")], a)                # not redrawn
  expect_equal(unname(a2["P04"]), unname(joint[1, 1]))
  expect_equal(unname(b2["P04"]), unname(joint[1, 2]))
  expect_equal(unname(b2[c("P01", "P02", "P03")]),
               (-12 / 150) * unname(a) + sqrt(4 - 144 / 150) * z_cond)
  # Distribution of the joint draw: A = sqrt(150) z1, B = -12/sqrt(150) z1 + ... z2
  expect_equal(unname(U_pen[1, 1]), sqrt(150))
  expect_equal(unname(U_pen[1, 2]), -12 / sqrt(150))
})

test_that("a three-phenotype block conditions each level on whatever subset it has stored", {
  pop <- make_ne_pop("ne_abc", traits = c("A", "B", "C"), n_males = 6, n_females = 6)
  on.exit(close_pop(pop))
  ids <- all_ids(pop)
  R <- sym(c("A", "B", "C"), c(1, .5, .2, .5, 1, .4, .2, .4, 1))
  pop <- assign_pens(pop, c("P01", "P02", "P03"), ids = ids)
  pop <- set_pen_effect(pop, c("A", "B", "C"), R)

  # P01, P02, P03 all get A; P01 and P02 also get B (C's animals filtered)
  pop <- add_ph(pop, "A", seed = 8)
  pop <- add_ph(pop, "B", tbl = get_table(pop, "ind_meta") |>
                  dplyr::filter(.data$pen != "P03"), seed = 9)
  a <- draws_of(pop, "A", "pen"); b <- draws_of(pop, "B", "pen")
  expect_identical(names(b), c("P01", "P02"))

  # C for everyone: P01/P02 condition on (A, B); P03 on A alone. One
  # resolver call (one sample set {C}); patterns split inside it, so the
  # stream is one normal per level in level order.
  s <- state_after(10, pop <- add_ph(pop, "C", seed = 10))
  expect_identical(s, state_after(10, stats::rnorm(3 + 12)))
  set.seed(10)
  z <- stats::rnorm(3)
  cc <- draws_of(pop, "C", "pen")
  # E[C | A, B] and V[C | A, B]
  S_ab <- R[c("A", "B"), c("A", "B")]; s_c <- R["C", c("A", "B")]
  for (k in 1:2) {
    p <- c("P01", "P02")[k]
    mu <- sum(solve(S_ab, s_c) * c(a[p], b[p]))
    v  <- 1 - sum(s_c * solve(S_ab, s_c))
    expect_equal(unname(cc[p]), mu + sqrt(v) * z[k])
  }
  expect_equal(unname(cc["P03"]), .2 * unname(a["P03"]) + sqrt(1 - .04) * z[3])
})


# ── Ordering, persistence, singletons ───────────────────────────────────────

test_that("effects draw in byte-sorted effect_name order, then blocks, then residuals; a NULL level draws nothing and contributes 0", {
  pop <- make_ne_pop("ne_order", traits = "A", n_males = 6, n_females = 6)
  on.exit(close_pop(pop))
  ids <- all_ids(pop)
  pop <- assign_pens(pop, c("P01", "P02"), ids = ids[1:10])   # A_11, A_12: NULL pen
  pop <- set_col(pop, "herd", "H1")
  pop <- set_col(pop, "herd", "H2", ids = ids[7:12])
  pop <- suppressMessages(define_effect_random(pop, "A", "pen", source_column = "pen", variance = 4))
  pop <- suppressMessages(define_effect_random(pop, "A", "herd", source_column = "herd", variance = 9))

  s <- state_after(11, pop <- add_ph(pop, "A", seed = 11))
  expect_identical(s, state_after(11, stats::rnorm(2 + 2 + 12)))
  set.seed(11)
  z_herd <- stats::rnorm(2); z_pen <- stats::rnorm(2)
  expect_equal(unname(draws_of(pop, "A", "herd")), 3 * z_herd)   # 'herd' < 'pen'
  expect_equal(unname(draws_of(pop, "A", "pen")),  2 * z_pen)
  expect_equal(unname(resid_of(pop, "A")), stats::rnorm(12))

  herd <- DBI::dbGetQuery(pop$db_conn, "SELECT id_ind, herd FROM ind_meta ORDER BY id_ind")
  expected <- 3 * z_herd[match(herd$herd, c("H1", "H2"))] +
    ifelse(is.na(pen_of(pop)), 0, 2 * z_pen[match(pen_of(pop), c("P01", "P02"))])
  expect_equal(unname(contrib_of(pop, "A")), unname(expected))
  expect_equal(nrow(re_rows_all(pop)), 4L)
})

test_that("a phenotype whose planned levels are all NULL draws nothing, alone or in a block, and contributes 0", {
  pop <- make_ne_pop("ne_allnull", n_males = 4, n_females = 4)
  on.exit(close_pop(pop))
  ids <- all_ids(pop)
  pop <- set_col(pop, "pen", NA_character_)                 # declare the column
  pop <- suppressMessages(define_effect_random(pop, "A", "pen", source_column = "pen", variance = 4))
  s <- state_after(30, pop <- add_ph(pop, "A", seed = 30))   # 1 x 1 block
  expect_identical(s, state_after(30, stats::rnorm(8)))      # residuals only
  expect_equal(nrow(re_rows_all(pop)), 0L)
  expect_equal(unname(contrib_of(pop, "A")), rep(0, 8))

  # In a block with B: A_1..A_4 (the males) have no pen, A_5..A_8 are in
  # P01/P02. A expressed in males only, so A plans no entity at all while
  # B plans both pens: A's coordinate is planned by nobody and stays
  # latent; B draws marginally -- 2 normals, then 4 + 8 residuals.
  pop <- suppressMessages(define_phenotype(
    pop, "A", type = "continuous", mean = 10, expressed_sex = "M",
    repeatable = TRUE, overwrite = TRUE))
  pop <- assign_pens(pop, c("P01", "P02"), ids = ids[5:8])
  pop <- suppressMessages(define_effect_cov_matrix(pop, "pen", R_pen))  # A's 1 x 1 joins B
  pop <- suppressMessages(define_effect_random(pop, "B", "pen", source_column = "pen"))
  expect_identical(unname(sex_of(pop)[ids[1:4]]), rep("M", 4))
  s <- state_after(31, pop <- add_ph(pop, c("A", "B"), seed = 31))
  expect_identical(s, state_after(31, stats::rnorm(2 + 4 + 8)))
  r <- re_rows_all(pop)
  expect_identical(r$phenotype_name, c("B", "B"))
  expect_identical(r$level, c("P01", "P02"))
  set.seed(31)
  expect_equal(r$draw_value, 2 * stats::rnorm(2))
  expect_equal(unname(contrib_of(pop, "A", 2L)), rep(0, 4))
})

test_that("a permanent-environment effect (source_column = id_ind) is one draw per animal reused across repeated records", {
  pop <- make_ne_pop("ne_pe", traits = "A", n_males = 4, n_females = 4)
  on.exit(close_pop(pop))
  pop <- suppressMessages(define_effect_random(pop, "A", "pe", source_column = "id_ind", variance = 2.25))

  pop <- add_ph(pop, "A", seed = 12)
  pe <- draws_of(pop, "A", "pe")
  expect_identical(names(pe), all_ids(pop))
  set.seed(12)
  expect_equal(unname(pe), 1.5 * stats::rnorm(8))
  expect_equal(unname(contrib_of(pop, "A", 1L)), unname(pe))

  s <- state_after(13, pop <- add_ph(pop, "A", seed = 13))
  expect_identical(s, state_after(13, stats::rnorm(8)))        # residuals only
  expect_equal(unname(contrib_of(pop, "A", 2L)), unname(pe))
  expect_equal(draws_of(pop, "A", "pe"), pe)
})

test_that("a 1 x 1 gamma or uniform block keeps its marginal sampler and still draws only new levels", {
  pop <- make_ne_pop("ne_gamma", traits = c("A", "B"), n_males = 6, n_females = 6)
  on.exit(close_pop(pop))
  ids <- all_ids(pop)
  pop <- assign_pens(pop, c("P01", "P02", "P03"))
  pop <- suppressMessages(define_effect_random(
    pop, "A", "pen", source_column = "pen", variance = 4, distribution = "gamma"))
  pop <- suppressMessages(define_effect_random(
    pop, "B", "pen", source_column = "pen", variance = 3, distribution = "uniform"))

  pop <- add_ph(pop, "A", tbl = get_table(pop, "ind_meta") |>
                  dplyr::filter(.data$pen != "P03"), seed = 14)
  set.seed(14)
  expect_equal(unname(draws_of(pop, "A", "pen")),
               stats::rgamma(2, shape = 1, rate = 1 / 2))
  # P03 is new; P01/P02 reused
  pop <- add_ph(pop, "A", seed = 15)
  set.seed(15)
  g <- stats::rgamma(1, shape = 1, rate = 1 / 2)
  expect_equal(unname(draws_of(pop, "A", "pen")["P03"]), g)
  expect_length(draws_of(pop, "A", "pen"), 3L)

  pop <- add_ph(pop, "B", seed = 16)
  set.seed(16)
  expect_equal(unname(draws_of(pop, "B", "pen")),
               stats::runif(3, min = -3, max = 3))
})

test_that("a block member without a random term for the effect is a latent coordinate; its stored draws still condition", {
  pop <- make_ne_pop("ne_latent", n_males = 6, n_females = 6)
  on.exit(close_pop(pop))
  pop <- assign_pens(pop, PENS)
  pop <- set_pen_effect(pop, c("A", "B"), R_pen)
  pop <- add_ph(pop, "A", seed = 17)
  a <- draws_of(pop, "A", "pen")

  # Drop A's pen term: A no longer draws pen, but its stored draws remain
  # and condition B's.
  pop <- suppressMessages(define_effect_random(
    pop, "A", "pen_dummy", source_column = "pen", variance = 0))
  DBI::dbExecute(pop$db_conn,
    "DELETE FROM phenotype_effects WHERE phenotype_name = 'A' AND effect_name = 'pen'")
  expect_equal(draws_of(pop, "A", "pen"), a)

  s <- state_after(18, pop <- add_ph(pop, c("A", "B"), seed = 18))
  # 'pen' (4 B draws) before 'pen_dummy' (4 levels; variance 0 still
  # consumes 4 normals), then 24 residuals
  expect_identical(s, state_after(18, stats::rnorm(4 + 4 + 24)))
  set.seed(18)
  z <- stats::rnorm(4)
  b <- draws_of(pop, "B", "pen")
  expect_equal(unname(b), (-12 / 150) * unname(a) + sqrt(4 - 144 / 150) * z)
  expect_equal(unname(contrib_of(pop, "A", 2L)), rep(0, 12))
})


# ── Backstops and integrity ─────────────────────────────────────────────────

test_that("§5.6 backstop: a row edited after definition is caught at add_phenotype()", {
  pop <- make_ne_pop("ne_backstop", n_males = 4, n_females = 4)
  on.exit(close_pop(pop))
  pop <- assign_pens(pop, c("P01", "P02"))
  pop <- set_col(pop, "herd", "H1")
  pop <- set_pen_effect(pop, c("A", "B"), R_pen)

  DBI::dbExecute(pop$db_conn,
    "UPDATE phenotype_effects SET distribution = 'gamma' WHERE phenotype_name = 'B' AND effect_name = 'pen'")
  expect_error(add_ph(pop, "A"),
               "add_phenotype\\(\\): 'pen' covariance block \\{A, B\\} requires distribution = \"normal\".*B uses \"gamma\"")
  DBI::dbExecute(pop$db_conn,
    "UPDATE phenotype_effects SET distribution = 'normal', source_column = 'herd' WHERE phenotype_name = 'B' AND effect_name = 'pen'")
  expect_error(add_ph(pop, "A"), "different grouping columns.*B reads ind_meta.herd")
  DBI::dbExecute(pop$db_conn,
    "UPDATE phenotype_effects SET source_column = 'pen', effect_class = 'fixed_class' WHERE phenotype_name = 'B' AND effect_name = 'pen'")
  expect_error(add_ph(pop, "A"), "can only join random effects")

  # Nothing was drawn or written by the rejected calls
  expect_equal(nrow(re_rows_all(pop)), 0L)
  expect_equal(DBI::dbGetQuery(pop$db_conn, "SELECT COUNT(*) AS n FROM ind_phenotype")$n, 0L)
})

test_that("a random term whose variance rows are gone, or a named block with strata, is an error before any draw", {
  pop <- make_ne_pop("ne_novar", n_males = 4, n_females = 4)
  on.exit(close_pop(pop))
  pop <- assign_pens(pop, c("P01", "P02"))
  pop <- set_pen_effect(pop, c("A", "B"), R_pen)

  DBI::dbExecute(pop$db_conn,
    "UPDATE phenotype_var_comp SET condition_table = 'ind_meta', condition_column = 'sex', condition_level = 'M' WHERE effect_name = 'pen'")
  s <- state_after(19, expect_error(add_ph(pop, "A"),
    "'pen' covariance block \\{A, B\\} has conditional strata; condition_column is residual-only"))
  expect_identical(s, state_after(19, NULL))

  DBI::dbExecute(pop$db_conn, "DELETE FROM phenotype_var_comp WHERE effect_name = 'pen'")
  s <- state_after(20, expect_error(add_ph(pop, "A"),
    "No variance stored for random effect 'pen' / phenotype 'A'"))
  expect_identical(s, state_after(20, NULL))
  expect_equal(DBI::dbGetQuery(pop$db_conn, "SELECT COUNT(*) AS n FROM ind_phenotype")$n, 0L)
})

test_that("D3 is live for named effects: the first call locks the block, remove_rows() on the draws clears it", {
  pop <- make_ne_pop("ne_lock", n_males = 4, n_females = 4)
  on.exit(close_pop(pop))
  pop <- assign_pens(pop, c("P01", "P02"))
  pop <- set_pen_effect(pop, c("A", "B"), R_pen)
  pop <- add_ph(pop, "A", seed = 21)
  expect_error(define_effect_cov_matrix(pop, "pen", 2 * R_pen), "realized draw")
  pop <- suppressMessages(
    pop |> get_table("phenotype_random_effects") |>
      dplyr::filter(.data$effect_name == "pen") |> remove_rows())
  pop <- suppressMessages(define_effect_cov_matrix(pop, "pen", 2 * R_pen))
  pop <- add_ph(pop, "B", seed = 22)
  set.seed(22)
  expect_equal(unname(draws_of(pop, "B", "pen")), sqrt(8) * stats::rnorm(2))  # marginal, no A stored
})
