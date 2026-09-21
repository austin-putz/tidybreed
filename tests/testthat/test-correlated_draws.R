# The pure resolver resolve_correlated_draws() and the block loader
# find_covariance_blocks(). See plans/sample_correlated_effects.md §5.2–5.4,
# §7 "Numerical", and phase_3.md.
#
# Every resolver test runs without a database. The RNG contract under test:
# a successful call consumes exactly n * m standard normals from rnorm(), in
# entity order; a rejected call and a call with nothing to draw consume none.

sym <- function(names, values) {
  n <- length(names)
  matrix(values, n, n, byrow = TRUE, dimnames = list(names, names))
}

R3 <- sym(c("A", "B", "C"), c(1.0, 0.5, 0.2,
                              0.5, 2.0, 0.3,
                              0.2, 0.3, 1.5))

# The normals a call will consume, in the order the resolver consumes them.
replay_z <- function(seed, n, m) {
  set.seed(seed)
  matrix(stats::rnorm(n * m), nrow = n, ncol = m, byrow = TRUE)
}

# Reference conditional moments by the textbook formula with a pseudoinverse.
cond_moments <- function(R, s, o, e_o) {
  if (length(o) == 0L) return(list(mean = rep(0, length(s)), cov = R[s, s, drop = FALSE]))
  W <- R[s, o, drop = FALSE] %*% MASS::ginv(R[o, o, drop = FALSE])
  list(mean = as.numeric(W %*% e_o),
       cov  = R[s, s, drop = FALSE] - W %*% R[o, s, drop = FALSE])
}


# ── One-dimensional cases ───────────────────────────────────────────────────

test_that("1 x 1 block, nothing observed: draws are sqrt(v) * z", {
  R <- sym("A", 4)
  set.seed(11)
  x <- resolve_correlated_draws(R, "A", 1:5)
  z <- replay_z(11, 5, 1)
  expect_equal(dim(x), c(5L, 1L))
  expect_identical(colnames(x), "A")
  expect_equal(as.numeric(x), 2 * as.numeric(z))
})

test_that("one sample coordinate conditional on one observed coordinate", {
  R <- R3[c("A", "B"), c("A", "B")]
  e <- c(0.7, -1.2, 0, 2.5)
  set.seed(12)
  x <- resolve_correlated_draws(R, "B", 1:4, observed = cbind(A = e))
  z <- replay_z(12, 4, 1)
  mu    <- R["B", "A"] / R["A", "A"] * e
  c_var <- R["B", "B"] - R["B", "A"]^2 / R["A", "A"]
  expect_equal(as.numeric(x), mu + sqrt(c_var) * as.numeric(z))
})

test_that("one entity, several sample coordinates, several observed", {
  set.seed(13)
  x <- resolve_correlated_draws(R3, c("B", "C"), "ind_1",
                                observed = cbind(A = 0.9))
  z <- replay_z(13, 1, 2)
  ref <- cond_moments(R3, c("B", "C"), "A", 0.9)
  expect_equal(as.numeric(x), ref$mean + as.numeric(z %*% chol(ref$cov)))
})


# ── Mixed patterns in one call ──────────────────────────────────────────────

test_that("entities with different observed patterns are resolved per pattern, in entity order", {
  obs <- rbind(c(NA,   NA),    # nothing observed
               c(0.5,  NA),    # A observed
               c(NA,  -0.3),   # B observed
               c(1.0,  0.4),   # A and B observed
               c(-0.8, NA))    # A observed again
  colnames(obs) <- c("A", "B")
  set.seed(14)
  x <- resolve_correlated_draws(R3, "C", 1:5, observed = obs)
  z <- replay_z(14, 5, 1)

  for (i in 1:5) {
    o   <- colnames(obs)[!is.na(obs[i, ])]
    ref <- cond_moments(R3, "C", o, obs[i, o])
    expect_equal(unname(x[i, "C"]), ref$mean + sqrt(as.numeric(ref$cov)) * z[i, 1],
                 info = paste("entity", i))
  }
})

test_that("a column that is NA for every entity is the same as not supplying it", {
  set.seed(15)
  x1 <- resolve_correlated_draws(R3, c("B", "C"), 1:3,
                                 observed = cbind(A = c(NA_real_, NA, NA)))
  set.seed(15)
  x2 <- resolve_correlated_draws(R3, c("B", "C"), 1:3)
  expect_equal(x1, x2)
  # a logical all-NA matrix (what cbind(A = c(NA, NA, NA)) gives) is accepted too
  set.seed(15)
  x3 <- resolve_correlated_draws(R3, c("B", "C"), 1:3, observed = cbind(A = c(NA, NA, NA)))
  expect_equal(x3, x2)
})

test_that("a latent coordinate is never realized and does not change the draw", {
  # Block {A, B, C}: sample B given A; C is latent. Must equal the 2-block draw.
  set.seed(16)
  x3 <- resolve_correlated_draws(R3, "B", 1:4, observed = cbind(A = c(1, 2, 3, 4)))
  set.seed(16)
  x2 <- resolve_correlated_draws(R3[c("A", "B"), c("A", "B")], "B", 1:4,
                                 observed = cbind(A = c(1, 2, 3, 4)))
  expect_equal(x3, x2)
  expect_identical(colnames(x3), "B")
})

test_that("a data frame is accepted for `observed` and for `entity_keys`", {
  set.seed(17)
  x1 <- resolve_correlated_draws(R3, "C", 1:2, observed = data.frame(A = c(1, NA)))
  set.seed(17)
  x2 <- resolve_correlated_draws(R3, "C", 1:2, observed = cbind(A = c(1, NA)))
  expect_equal(x1, x2)
  # a data frame of keys counts rows, not columns
  set.seed(17)
  x3 <- resolve_correlated_draws(R3, "C",
                                 data.frame(id_ind = c("A_1", "A_2"), pheno_number = c(1L, 1L),
                                            extra = 0, more = 0),
                                 observed = cbind(A = c(1, NA)))
  expect_equal(x3, x2)
})


# ── RNG contract ────────────────────────────────────────────────────────────

test_that("a call consumes exactly n * m normals, whatever the patterns", {
  obs <- cbind(A = c(NA, 0.5, NA, 1.0), B = c(NA, NA, -0.3, 0.4))
  set.seed(21)
  invisible(resolve_correlated_draws(R3, "C", 1:4, observed = obs))
  s_call <- get(".Random.seed", envir = globalenv())
  set.seed(21)
  invisible(stats::rnorm(4 * 1))
  expect_identical(s_call, get(".Random.seed", envir = globalenv()))

  set.seed(22)
  invisible(resolve_correlated_draws(R3, c("B", "C"), 1:7))
  s_call <- get(".Random.seed", envir = globalenv())
  set.seed(22)
  invisible(stats::rnorm(7 * 2))
  expect_identical(s_call, get(".Random.seed", envir = globalenv()))
})

test_that("zero entities and zero sample coordinates are RNG-neutral and shaped", {
  set.seed(23)
  before <- get(".Random.seed", envir = globalenv())

  x0 <- resolve_correlated_draws(R3, c("A", "B"), character(0))
  expect_equal(dim(x0), c(0L, 2L))
  expect_identical(colnames(x0), c("A", "B"))
  expect_identical(get(".Random.seed", envir = globalenv()), before)

  x0 <- resolve_correlated_draws(R3, c("A", "B"), NULL)
  expect_equal(dim(x0), c(0L, 2L))
  expect_identical(get(".Random.seed", envir = globalenv()), before)

  x0 <- resolve_correlated_draws(R3, c("A", "B"),
                                 data.frame(id_ind = character(0), pheno_number = integer(0)))
  expect_equal(dim(x0), c(0L, 2L))
  expect_identical(get(".Random.seed", envir = globalenv()), before)

  x1 <- resolve_correlated_draws(R3, character(0), 1:5,
                                 observed = cbind(A = 1:5, B = 1:5))
  expect_equal(dim(x1), c(5L, 0L))
  expect_identical(get(".Random.seed", envir = globalenv()), before)

  x2 <- resolve_correlated_draws(R3, character(0), character(0))
  expect_equal(dim(x2), c(0L, 0L))
  expect_identical(get(".Random.seed", envir = globalenv()), before)
})

test_that("a rejected call leaves the RNG untouched", {
  set.seed(24)
  before <- get(".Random.seed", envir = globalenv())
  v  <- c(A = 1, B = 1, C = 0.2)
  P3 <- outer(v, v); P3["C", "C"] <- 1.5
  expect_error(resolve_correlated_draws(P3, "C", 1:3,
                                        observed = cbind(A = c(1, 1, 2), B = c(1, 1, 2 + 1e-3))),
               "outside the support")
  expect_error(resolve_correlated_draws(R3, "Z", 1:3), "not in the block")
  expect_identical(get(".Random.seed", envir = globalenv()), before)
})

test_that("zero conditional variance still consumes its normals and returns the mean exactly", {
  P <- sym(c("A", "B"), c(1, 1, 1, 1))
  set.seed(25)
  x <- resolve_correlated_draws(P, "B", 1:3, observed = cbind(A = c(1, -2, 0.25)))
  expect_equal(as.numeric(x), c(1, -2, 0.25))
  s_call <- get(".Random.seed", envir = globalenv())
  set.seed(25)
  invisible(stats::rnorm(3))
  expect_identical(s_call, get(".Random.seed", envir = globalenv()))
})


# ── Singular and near-singular blocks ───────────────────────────────────────

test_that("perfect correlation: the sample coordinate is determined by the observed one", {
  P <- sym(c("A", "B"), c(4, 2, 2, 1))          # corr = 1, sd 2 and 1
  set.seed(31)
  x <- resolve_correlated_draws(P, "B", 1:3, observed = cbind(A = c(2, -4, 1)))
  expect_equal(as.numeric(x), c(1, -2, 0.5))
})

test_that("near-perfect correlation gives a tiny but positive conditional variance", {
  rho <- 1 - 1e-6
  P <- sym(c("A", "B"), c(1, rho, rho, 1))
  set.seed(32)
  x <- resolve_correlated_draws(P, "B", 1:2, observed = cbind(A = c(1, 1)))
  z <- replay_z(32, 2, 1)
  expect_equal(as.numeric(x), rho * 1 + sqrt(1 - rho^2) * as.numeric(z))
})

test_that("a zero-variance coordinate is a valid block member", {
  Z <- sym(c("A", "B"), c(1, 0, 0, 0))
  # sampled: always exactly 0
  set.seed(33)
  x <- resolve_correlated_draws(Z, c("A", "B"), 1:3)
  z <- replay_z(33, 3, 2)
  expect_equal(x[, "B"], c(0, 0, 0))
  expect_equal(x[, "A"], z[, 1])
  # observed at 0: consistent, conditions nothing
  set.seed(34)
  x1 <- resolve_correlated_draws(Z, "A", 1:3, observed = cbind(B = c(0, 0, 0)))
  set.seed(34)
  x2 <- resolve_correlated_draws(Z, "A", 1:3)
  expect_equal(x1, x2)
})

test_that("support consistency: off-support observations error, in-support proceed", {
  Z <- sym(c("A", "B"), c(1, 0, 0, 0))
  expect_error(resolve_correlated_draws(Z, "A", 1:2, observed = cbind(B = c(0, 0.1))),
               "outside the support.*1 of 2 entities.*entity 2: B = 0.1")
  expect_error(resolve_correlated_draws(Z, "A", c("A_7", "A_9"), observed = cbind(B = c(0, 0.1))),
               "entity A_9: B = 0.1")
  expect_error(resolve_correlated_draws(Z, "A", data.frame(id_ind = c("A_7", "A_9"), pheno_number = 1:2),
                                        observed = cbind(B = c(0, 0.1))),
               "entity id_ind = A_9, pheno_number = 2: B = 0.1")
  # numerically-zero noise on a zero-variance coordinate is inside the support
  expect_silent(resolve_correlated_draws(Z, "A", 1:2, observed = cbind(B = c(0, 1e-12))))

  # A and B perfectly correlated (B = sqrt(2) * A), C partly correlated with both.
  v  <- c(A = 1, B = sqrt(2), C = 0.2)
  P3 <- outer(v, v); P3["C", "C"] <- 1.5
  expect_silent(resolve_correlated_draws(P3, "C", 1:2,
                                         observed = cbind(A = c(1, -1), B = sqrt(2) * c(1, -1))))
  expect_error(resolve_correlated_draws(P3, "C", 1:2,
                                        observed = cbind(A = c(1, 1), B = c(sqrt(2), 2))),
               "outside the support.*1 of 2 entities.*entity 2")
  # the same values are fine when only one of the two is observed
  expect_silent(resolve_correlated_draws(P3, "C", 1:2, observed = cbind(B = c(sqrt(2), 2))))
})

test_that("the tolerance is scale-aware: the same relative perturbation at variance 1e-6 and 1e6", {
  for (v in c(1e-6, 1, 1e6)) {
    sd <- sqrt(v)
    P  <- sym(c("A", "B", "C"), c(v, v, 0,
                                  v, v, 0,
                                  0, 0, v))
    ok  <- cbind(A = sd, B = sd * (1 + 1e-10))
    bad <- cbind(A = sd, B = sd * (1 + 1e-6))
    expect_silent(resolve_correlated_draws(P, "C", 1, observed = ok))
    expect_error(resolve_correlated_draws(P, "C", 1, observed = bad),
                 "outside the support", info = paste("variance", v))
  }
})

test_that("tiny negative eigenvalues are absorbed; materially negative are rejected", {
  # Rank-1 PSD matrix perturbed by rounding-scale noise off the diagonal.
  u <- c(1, 2, 3); names(u) <- c("A", "B", "C")
  P <- outer(u, u)
  P_tiny <- P; P_tiny["A", "C"] <- P_tiny["C", "A"] <- P["A", "C"] + 1e-13
  ev <- eigen(P_tiny, symmetric = TRUE, only.values = TRUE)$values
  expect_true(min(ev) < 0)                          # it really is slightly indefinite
  set.seed(41)
  expect_silent(x <- resolve_correlated_draws(P_tiny, c("B", "C"), 1:3,
                                              observed = cbind(A = c(1, 2, 3))))
  expect_equal(x[, "B"], 2 * c(1, 2, 3), tolerance = 1e-6)
  expect_equal(x[, "C"], 3 * c(1, 2, 3), tolerance = 1e-6)

  P_bad <- P; P_bad["A", "C"] <- P_bad["C", "A"] <- P["A", "C"] + 0.5
  expect_error(resolve_correlated_draws(P_bad, "B", 1, observed = cbind(A = 1)),
               "not positive semi-definite")
})

test_that("extreme but valid variance scales are stable", {
  for (v in c(1e-12, 1e12)) {
    R <- R3 * v
    set.seed(42)
    x <- resolve_correlated_draws(R, c("B", "C"), 1:3, observed = cbind(A = sqrt(v) * c(1, -1, 0.5)))
    z <- replay_z(42, 3, 2)
    for (i in 1:3) {
      ref <- cond_moments(R, c("B", "C"), "A", sqrt(v) * c(1, -1, 0.5)[i])
      expect_equal(x[i, ], ref$mean + as.numeric(z[i, ] %*% chol(ref$cov)),
                   tolerance = 1e-10, ignore_attr = TRUE, info = paste(v, i))
    }
  }
})


# ── Distributional check ────────────────────────────────────────────────────

test_that("conditional draws reproduce the conditional covariance (joint over many entities)", {
  n <- 20000L
  set.seed(51)
  A <- stats::rnorm(n, sd = 1)
  x <- resolve_correlated_draws(R3, c("B", "C"), seq_len(n), observed = cbind(A = A))
  ref <- cond_moments(R3, c("B", "C"), "A", 0)
  resid <- x - cbind(R3["B", "A"] * A, R3["C", "A"] * A)
  S <- stats::cov(resid)
  # sampling sd of a covariance estimate ~ sqrt((s_ii s_jj + s_ij^2) / n)
  se <- sqrt((outer(diag(ref$cov), diag(ref$cov)) + ref$cov^2) / n)
  expect_true(all(abs(S - ref$cov) < 4 * se))
  expect_true(all(abs(colMeans(resid)) < 4 * sqrt(diag(ref$cov) / n)))
})


# ── Validation ──────────────────────────────────────────────────────────────

test_that("resolver input validation", {
  expect_error(resolve_correlated_draws(matrix(1, 2, 2), "A", 1),
               "identical, unique, non-empty row and column names")
  M <- R3; colnames(M) <- c("A", "B", "X")
  expect_error(resolve_correlated_draws(M, "A", 1), "identical, unique")
  M <- R3; M["A", "B"] <- NA
  expect_error(resolve_correlated_draws(M, "A", 1), "finite")
  M <- R3; M["A", "B"] <- 0.9
  expect_error(resolve_correlated_draws(M, "A", 1), "symmetric within tolerance")
  expect_error(resolve_correlated_draws(R3, c("A", "A"), 1), "duplicates")
  expect_error(resolve_correlated_draws(R3, "Z", 1), "not in the block: \\{Z\\}")
  expect_error(resolve_correlated_draws(R3, 1, 1), "character vector")
  expect_error(resolve_correlated_draws(R3, "A", 1:3, observed = cbind(B = 1:2)),
               "2 rows but `entity_keys` has 3")
  expect_error(resolve_correlated_draws(R3, "A", 1:2, observed = cbind(Z = 1:2)),
               "`observed` columns not in the block: \\{Z\\}")
  expect_error(resolve_correlated_draws(R3, "A", 1:2, observed = cbind(A = 1:2)),
               "both observed and sampled: \\{A\\}")
  expect_error(resolve_correlated_draws(R3, "A", 1:2, observed = cbind(B = c(1, Inf))),
               "finite or NA")
  expect_error(resolve_correlated_draws(R3, "A", 1:2, observed = matrix(1:2, 2, 1)),
               "unique column names")
  expect_error(resolve_correlated_draws(R3, "A", 1, tolerance = -1),
               "non-negative")
  expect_error(resolve_correlated_draws(R3, "A", 1:2, observed = cbind(B = c("a", "b"))),
               "numeric matrix")
})


# ── find_covariance_blocks() ────────────────────────────────────────────────

# The loader is the one function in this file that needs a database. A bare
# open_pop() is enough: phenotype_var_comp exists from the start and the
# writers do not require phenotype_meta rows.
make_cov_pop <- function(name) {
  suppressMessages(open_pop(pop_name = name, db_name = ":memory:"))
}

test_that("find_covariance_blocks(): components, strata, singletons, and absent phenotypes", {
  pop <- make_cov_pop("fcb_basic")
  on.exit(close_pop(pop))

  # Declared as {B, A}: the block is returned in sorted order.
  pop <- define_residual_cov(pop, c("B", "A"), sym(c("B", "A"), c(2, .3, .3, 1)))
  pop <- define_residual_cov(pop, c("A", "B"), sym(c("A", "B"), c(1.5, .4, .4, 2.5)),
                             condition_column = "sex", condition_level = "M")
  pop <- define_residual_cov(pop, c("A", "B"), sym(c("A", "B"), c(0.5, .1, .1, 1.0)),
                             condition_column = "sex", condition_level = "F")
  pop <- define_residual_cov(pop, "C", sym("C", 4))
  pop <- suppressMessages(
    define_effect_cov_matrix(pop, "pen", sym(c("A", "C"), c(1, .1, .1, 1))))

  blocks <- find_covariance_blocks(pop$db_conn, "residual", c("C", "A", "D"))
  expect_length(blocks, 2L)
  expect_null(names(blocks))
  expect_identical(blocks[[1]]$phenotypes, c("A", "B"))
  expect_identical(blocks[[2]]$phenotypes, "C")

  ab <- blocks[[1]]
  expect_identical(ab$effect_name, "residual")
  expect_identical(ab$condition_table, "ind_meta")
  expect_identical(ab$condition_column, "sex")
  expect_equal(ab$unconditional, sym(c("A", "B"), c(1, .3, .3, 2)))
  expect_identical(names(ab$conditional), c("F", "M"))
  expect_equal(ab$conditional$M, sym(c("A", "B"), c(1.5, .4, .4, 2.5)))
  expect_equal(ab$conditional$F, sym(c("A", "B"), c(0.5, .1, .1, 1.0)))

  cc <- blocks[[2]]
  expect_null(cc$condition_column)
  expect_null(cc$condition_table)
  expect_equal(cc$unconditional, sym("C", 4))
  expect_identical(cc$conditional, list())

  # Asking for one member returns the whole component.
  expect_identical(find_covariance_blocks(pop$db_conn, "residual", "B")[[1]]$phenotypes,
                   c("A", "B"))
  # A phenotype with no rows is in no block; a different effect is a different graph.
  expect_identical(find_covariance_blocks(pop$db_conn, "residual", "D"), list())
  expect_identical(find_covariance_blocks(pop$db_conn, "residual", character(0)), list())
  pen <- find_covariance_blocks(pop$db_conn, "pen", "C")
  expect_length(pen, 1L)
  expect_identical(pen[[1]]$phenotypes, c("A", "C"))
  expect_identical(find_covariance_blocks(pop$db_conn, "pen", "B"), list())
  # Callers find the absent phenotypes with setdiff.
  expect_identical(setdiff(c("C", "A", "D"),
                           unlist(lapply(blocks, `[[`, "phenotypes"))), "D")
})

test_that("find_covariance_blocks(): a block with only conditional strata", {
  pop <- make_cov_pop("fcb_cond_only")
  on.exit(close_pop(pop))
  pop <- define_residual_cov(pop, "A", sym("A", 9), condition_column = "farm",
                             condition_table = "ind_meta", condition_level = "x")
  b <- find_covariance_blocks(pop$db_conn, "residual", "A")[[1]]
  expect_null(b$unconditional)
  expect_identical(b$condition_column, "farm")
  expect_equal(b$conditional$x, sym("A", 9))
})

test_that("find_covariance_blocks(): the stored matrix is exactly symmetric even when the input was only symmetric within tolerance", {
  pop <- make_cov_pop("fcb_sym")
  on.exit(close_pop(pop))
  R <- sym(c("A", "B"), c(1, .3, .3, 2))
  R["A", "B"] <- 0.3 + 1e-12          # accepted by the writer's tolerance
  pop <- define_residual_cov(pop, c("A", "B"), R)
  M <- find_covariance_blocks(pop$db_conn, "residual", "A")[[1]]$unconditional
  expect_identical(M, t(M))
  expect_equal(M["A", "B"], 0.3 + 5e-13)
})

test_that("find_covariance_blocks(): hand-removed rows are reported, not guessed around", {
  pop <- make_cov_pop("fcb_broken")
  on.exit(close_pop(pop))
  pop <- define_residual_cov(pop, c("A", "B"), sym(c("A", "B"), c(1, .3, .3, 2)))
  pop <- define_residual_cov(pop, c("A", "B"), sym(c("A", "B"), c(1, .3, .3, 2)),
                             condition_column = "sex", condition_level = "M")

  # one pair row gone from the unconditional stratum
  suppressMessages(
    pop |> get_table("phenotype_var_comp") |>
      dplyr::filter(effect_name == "residual", phenotype_name_1 == "A",
                    phenotype_name_2 == "B", is.na(condition_column)) |>
      remove_rows())
  expect_error(find_covariance_blocks(pop$db_conn, "residual", "B"),
               "block \\{A, B\\} is incomplete in the unconditional stratum: no row for \\(A, B\\)")
  expect_error(find_covariance_blocks(pop$db_conn, "residual", "B"),
               'define_residual_cov\\(pop, c\\("A", "B"\\), R\\)')

  # the whole unconditional stratum gone: the conditional one still loads
  suppressMessages(
    pop |> get_table("phenotype_var_comp") |>
      dplyr::filter(effect_name == "residual", is.na(condition_column)) |>
      remove_rows())
  b <- find_covariance_blocks(pop$db_conn, "residual", "A")[[1]]
  expect_null(b$unconditional)
  expect_identical(names(b$conditional), "M")

  # B's rows gone from the conditional stratum too: {A, B} is still one block
  # by A's remaining (A, B) row, but the stratum is incomplete
  suppressMessages(
    pop |> get_table("phenotype_var_comp") |>
      dplyr::filter(effect_name == "residual", phenotype_name_1 == "B") |>
      remove_rows())
  expect_error(find_covariance_blocks(pop$db_conn, "residual", "A"),
               "incomplete in the ind_meta.sex = 'M' stratum: no row for \\(B, A\\), \\(B, B\\)")
})

test_that("find_covariance_blocks(): two condition columns on one block is an error", {
  pop <- make_cov_pop("fcb_two_cols")
  on.exit(close_pop(pop))
  pop <- define_residual_cov(pop, "A", sym("A", 1), condition_column = "sex",
                             condition_level = "M")
  # The writer refuses a second column; forge the state directly.
  DBI::dbExecute(pop$db_conn,
    "INSERT INTO phenotype_var_comp
       (id_phenotype_var_comp, effect_name, phenotype_name_1, phenotype_name_2,
        cov_value, condition_column, condition_table, condition_level)
     VALUES (999, 'residual', 'A', 'A', 2, 'farm', 'ind_meta', 'x')")
  expect_error(find_covariance_blocks(pop$db_conn, "residual", "A"),
               "conditioned on more than one column \\(ind_meta.farm, ind_meta.sex\\)")

  # a conditional row with no level is reported, not silently dropped
  DBI::dbExecute(pop$db_conn,
    "UPDATE phenotype_var_comp SET condition_column = 'sex', condition_level = NULL
     WHERE id_phenotype_var_comp = 999")
  expect_error(find_covariance_blocks(pop$db_conn, "residual", "A"),
               "NULL condition_table or condition_level")
})

test_that("find_covariance_blocks() does not touch the RNG", {
  pop <- make_cov_pop("fcb_rng")
  on.exit(close_pop(pop))
  pop <- define_residual_cov(pop, c("A", "B"), sym(c("A", "B"), c(1, .3, .3, 2)))
  set.seed(61)
  before <- get(".Random.seed", envir = globalenv())
  invisible(find_covariance_blocks(pop$db_conn, "residual", "A"))
  expect_identical(get(".Random.seed", envir = globalenv()), before)
})

test_that("loader output feeds the resolver directly", {
  pop <- make_cov_pop("fcb_feed")
  on.exit(close_pop(pop))
  pop <- define_residual_cov(pop, c("A", "B", "C"), R3)
  b <- find_covariance_blocks(pop$db_conn, "residual", "B")[[1]]
  set.seed(62)
  x1 <- resolve_correlated_draws(b$unconditional, c("B", "C"), 1:3,
                                 observed = cbind(A = c(1, NA, 0.5)))
  set.seed(62)
  x2 <- resolve_correlated_draws(R3, c("B", "C"), 1:3,
                                 observed = cbind(A = c(1, NA, 0.5)))
  expect_equal(x1, x2)
})

test_that(".blupf90_residual_cov(): block-diagonal over the unconditional strata; strata-only and absent traits error", {
  pop <- make_cov_pop("fcb_blup")
  on.exit(close_pop(pop))
  pop <- define_residual_cov(pop, c("A", "B"), sym(c("A", "B"), c(2, .3, .3, 1)))
  pop <- define_residual_cov(pop, "C", sym("C", 4))
  pop <- define_residual_cov(pop, "D", sym("D", 9),
                             condition_column = "sex", condition_level = "M")

  R <- .blupf90_residual_cov(pop, c("C", "B", "A"))
  expect_equal(R, matrix(c(4, 0, 0, 0, 1, .3, 0, .3, 2), 3, 3,
                         dimnames = list(c("C", "B", "A"), c("C", "B", "A"))))
  expect_error(.blupf90_residual_cov(pop, c("A", "E")),
               "Residual covariance matrix not found for traits: E")
  expect_error(.blupf90_residual_cov(pop, c("A", "D")),
               "block \\{D\\} has only conditional strata")
})
