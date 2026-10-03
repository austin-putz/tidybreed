# The ported QTL-effect method (R/qtl_congruence.R), tested on the internals.
# Tests "paper-N" are the source project's tests/test_qtl_effects_paper.R,
# numbered and named as there, pointed at the tidybreed port (never at the
# frozen paper file). Gates A1, A2, A4 and A5 of
# plans/import_qtl_effect_methods.md §11 are on the same internals here; the
# generator-level gates are in test-define_additive_effects-anchor.R.

test_that("paper-1: exactness under 'realised'", {
  G <- qtl_paper_G(); X <- qtl_paper_X()
  fit <- .qtl_sim_effects(X, G, anchor = "realised")
  expect_lt(max(abs(fit$G_anchor - G)), 1e-9)
  expect_lt(max(abs(stats::cov(qtl_breeding_values(X, fit)) - G)), 1e-9)
})

test_that("paper-2: all-heterozygote panel keeps every locus (Codex F1)", {
  X_het <- matrix(1, 300, 5)
  fit <- .qtl_sim_effects(X_het, matrix(4, 1, 1), anchor = "genic",
                          warn_bounds = NULL)   # no founder variance: warns by design
  expect_equal(fit$n_causal, 5L)
})

test_that("paper-3: rank-1 target is feasible on one locus (Codex F2)", {
  set.seed(1002)
  X_one <- matrix(stats::rbinom(600, 2, 0.3), 600, 1)
  G_rank1 <- matrix(c(4, 6, 6, 9), 2, 2)
  fit <- .qtl_sim_effects(X_one, G_rank1, anchor = "genic")
  expect_lt(max(abs(fit$G_anchor - G_rank1)), 1e-9)
})

test_that("paper-4: matrix-aware diagnostics see through a cancelling trace (F4)", {
  set.seed(33)
  n_c <- 12000
  blk <- function(n) {
    a <- stats::rbinom(n, 1, 0.5)
    b <- ifelse(stats::rbinom(n, 1, 0.1) == 1, 1 - a, a)
    cbind(a, b)
  }
  Xc4 <- cbind(blk(n_c), blk(n_c)) + cbind(blk(n_c), blk(n_c))
  B_c <- cbind(c(sqrt(5), sqrt(5), 0, 0), c(0, 0, sqrt(5), -sqrt(5)))
  G_c <- stats::cov(scale(Xc4, TRUE, FALSE) %*% B_c)
  expect_warning(
    fit <- .qtl_sim_effects(Xc4, G_c, architecture = B_c, anchor = "realised"),
    "Realised variance departs")
  expect_lt(abs(fit$trace_ratio_equilibrium - 1), 0.08)
  expect_lt(min(fit$equilibrium_relative_eigen), 0.7)
  expect_gt(max(fit$equilibrium_relative_eigen), 3)
})

test_that("paper-5: dual anchor is exact at founders and equilibrium", {
  skip_on_cran()
  set.seed(37)
  X_d <- qtl_mkhap(1500) + qtl_mkhap(1500)
  G <- qtl_paper_G()
  fit <- .qtl_sim_effects(X_d, G, anchor = "dual")
  expect_lt(max(abs(fit$G_founder - G)), 1e-8)
  expect_lt(max(abs(fit$G_equilibrium - G)), 1e-8)
})

test_that("paper-6: custom reference anchor, weights and genotypes", {
  G <- qtl_paper_G(); X <- qtl_paper_X()
  set.seed(4242)
  p_r <- colMeans(X) / 2
  w_self <- 2 * p_r * (1 - p_r) * 1.5
  fit_w <- .qtl_sim_effects(X, G, anchor = "reference", reference = w_self,
                            warn_bounds = NULL)
  expect_lt(max(abs(crossprod(fit_w$B * sqrt(w_self)) - G)), 1e-9)

  X_ref <- matrix(stats::rbinom(900 * 30, 2, 0.45), 900, 30)
  fit_ref <- .qtl_sim_effects(X, G, anchor = "reference", reference = X_ref,
                              warn_bounds = NULL)
  expect_lt(max(abs(stats::cov(sweep(X_ref, 2, colMeans(X_ref), "-") %*%
                                 fit_ref$B) - G)), 1e-9)
})

test_that("paper-7: singular target under the dual anchor (q, not k)", {
  skip_on_cran()
  set.seed(37)
  X_d <- qtl_mkhap(1500) + qtl_mkhap(1500)
  G_sing <- matrix(c(4, 2, 2, 1), 2, 2)
  fit <- .qtl_sim_effects(X_d, G_sing, anchor = "dual")
  expect_equal(fit$rank_G, 1L)
  expect_lt(max(abs(fit$G_founder - G_sing)), 1e-8)
  expect_lt(max(abs(fit$G_equilibrium - G_sing)), 1e-8)
})

test_that("paper-8: invalid input is rejected for the right reason", {
  G <- qtl_paper_G(); X <- qtl_paper_X()
  set.seed(37)
  X_d <- qtl_mkhap(200) + qtl_mkhap(200)
  expect_error(.qtl_sim_effects(X, G, ploidy = 2.5), "`ploidy` must be an integer")
  expect_error(.qtl_sim_effects(X, G, ploidy = c(2, 2)), "`ploidy` must be a single")
  expect_error(.qtl_sim_effects(X, G, rel_tol = -1), "`rel_tol` must lie in")
  expect_error(.qtl_sim_effects(X, G, warn_bounds = c(1.2, 0.8)),
               "0 < lower <= upper")
  expect_error(.qtl_sim_effects(X, G, anchor = "reference", reference = rep(1, 5)),
               "one weight per locus")
  expect_error(.qtl_sim_effects(X, G, anchor = "reference",
                                reference = c(rep(1, 29), -1)),
               "non-negative")
  expect_error(.qtl_sim_effects(X, G, anchor = "reference",
                                reference = matrix(1, 100, 29)),
               "one column per locus")
  expect_error(.qtl_sim_effects(X, G, architecture = matrix(0, 30, 2)),
               "must be exactly 30 x 3")
  expect_error(.qtl_sim_effects(X_d, G, anchor = "dual",
                                architecture = matrix(1, ncol(X_d), 3)),
               "cannot be used with `anchor = \"dual\"`")
})

test_that("paper-9: inbred panel with repulsion LD is endpoint-feasible (dual)", {
  skip_on_cran()
  n_i <- 4000
  h_in <- rbind(c(0, 1), c(1, 0))
  H_in <- h_in[rep(1:2, each = n_i / 2), ]
  X_in2 <- H_in + H_in
  fit <- .qtl_sim_effects(X_in2, matrix(4, 1, 1), anchor = "dual",
                          warn_bounds = NULL)
  expect_lt(abs(fit$G_founder[1, 1] - 4), 1e-6)
  expect_lt(abs(fit$G_equilibrium[1, 1] - 4), 1e-6)

  # paper-10: the same solution misses generation 1 (Theorem 4). For a fully
  # inbred panel S1 - Delta = (S0 - diag(S0)) / 2.
  p_in <- colMeans(X_in2) / 2
  S0_in <- stats::cov(X_in2)
  S1_in <- diag(2 * p_in * (1 - p_in)) + (S0_in - diag(diag(S0_in))) / 2
  v1 <- drop(crossprod(fit$B, S1_in %*% fit$B))
  expect_gt(abs(v1 - 4), 1)
})

test_that("paper-11: dual path, non-HWE founders miss later generations, HWE do not", {
  skip_on_cran()
  # The source's own free-recombination transmission, verbatim: this is a
  # property of the dual anchor, which is internal only (Q6).
  advance <- function(A1, A2, n_off) {
    mk <- function() {
      par <- sample.int(nrow(A1), n_off, TRUE)
      cf  <- matrix(stats::rbinom(n_off * ncol(A1), 1, 0.5), n_off, ncol(A1))
      ifelse(cf == 1, A1[par, , drop = FALSE], A2[par, , drop = FALSE])
    }
    list(mk(), mk())
  }
  path <- function(B, A1, A2, p0, n_gen) {
    out <- numeric(n_gen + 1)
    for (g in 0:n_gen) {
      out[g + 1] <- stats::var(drop(sweep(A1 + A2, 2, 2 * p0, "-") %*% B))
      nx <- advance(A1, A2, nrow(A1)); A1 <- nx[[1]]; A2 <- nx[[2]]
    }
    out
  }
  n_i <- 4000
  H_in <- rbind(c(0, 1), c(1, 0))[rep(1:2, each = n_i / 2), ]
  X_in2 <- H_in + H_in
  fit_in <- .qtl_sim_effects(X_in2, matrix(4, 1, 1), anchor = "dual",
                             warn_bounds = NULL)
  set.seed(52)
  pi_nonhwe <- path(fit_in$B, H_in, H_in, colMeans(X_in2) / 2, 3)
  expect_lt(abs(pi_nonhwe[1] - 4), 0.1)
  expect_gt(abs(pi_nonhwe[2] - 4), 1)

  set.seed(53)
  Hh1 <- qtl_mkhap(3000); Hh2 <- qtl_mkhap(3000); X_h <- Hh1 + Hh2
  fit_h <- .qtl_sim_effects(X_h, matrix(4, 1, 1), anchor = "dual",
                            warn_bounds = NULL)
  set.seed(54)
  pi_hwe <- path(fit_h$B, Hh1, Hh2, colMeans(X_h) / 2, 4)
  expect_lt(max(abs(pi_hwe - 4)), 0.6)
})

# ── Gates on the internals ───────────────────────────────────────────────────

test_that("A1: k = 1 genic congruence is the scalar rescale", {
  for (s in c(11, 12, 13)) {
    set.seed(s)
    p  <- stats::runif(50, 0.05, 0.95)
    w  <- 2 * p * (1 - p)
    b0 <- matrix(stats::rnorm(50), 50, 1)
    for (g in c(1e-6, 4, 1e6)) {
      cg <- .qtl_congruence(b0, .qtl_target_eigen(matrix(g)),
                            .qtl_anchor_diag(w))
      expect_equal(cg$B, b0 * sqrt(g / sum(w * b0^2)), tolerance = 1e-12)
    }
  }
})

test_that("A2: k = 2 genic is exact and matches a dense C^-1/2 G^1/2 oracle", {
  G_full  <- matrix(c(4, 1.2, 1.2, 9), 2, 2)
  G_rank1 <- matrix(c(4, 6, 6, 9), 2, 2)
  for (s in c(21, 22, 23)) for (m in c(2L, 200L)) for (G in list(G_full, G_rank1)) {
    set.seed(s)
    p  <- stats::runif(m, 0.05, 0.95)
    w  <- 2 * p * (1 - p)
    B0 <- matrix(stats::rnorm(m * 2), m, 2)
    cg <- .qtl_congruence(B0, .qtl_target_eigen(G), .qtl_anchor_diag(w))
    expect_lt(max(abs(crossprod(cg$B * sqrt(w)) - G)), 1e-10)
    C  <- crossprod(B0 * sqrt(w))
    oracle <- B0 %*% qtl_sym_pow(C, -1 / 2) %*% qtl_sym_pow(G, 1 / 2)
    expect_lt(max(abs(cg$B - oracle)), 1e-10 * max(1, max(abs(oracle))))
  }
})

test_that("A4: rank-deficient architecture, rank-1 target in another direction", {
  set.seed(31)
  m  <- 40
  p  <- stats::runif(m, 0.1, 0.9)
  w  <- 2 * p * (1 - p)
  b  <- stats::rnorm(m)
  B0 <- cbind(b, 2 * b)                      # proportional columns: rank 1
  G  <- matrix(c(1, -2, -2, 4), 2, 2)        # rank 1, direction (1, -2)
  cg <- .qtl_congruence(B0, .qtl_target_eigen(G), .qtl_anchor_diag(w))
  expect_lt(max(abs(crossprod(cg$B * sqrt(w)) - G)), 1e-10)
})

test_that("A5: the two rank errors are distinct and name the ranks", {
  G <- matrix(c(4, 1, 1, 9), 2, 2)
  # (a) the anchor: one segregating locus cannot carry a rank-2 target.
  B0 <- matrix(stats::rnorm(6), 3, 2)
  expect_error(
    .qtl_congruence(B0, .qtl_target_eigen(G), .qtl_anchor_diag(c(0.5, 0, 0))),
    "anchor cannot carry the target.*rank\\(G\\) = 2.*rank 1")
  # (b) the architecture: proportional columns on a rank-20 anchor.
  set.seed(32)
  b  <- stats::rnorm(20)
  w  <- stats::runif(20, 0.1, 0.5)
  expect_error(
    .qtl_congruence(cbind(b, 3 * b), .qtl_target_eigen(G), .qtl_anchor_diag(w)),
    "drawn architecture cannot carry the target.*= 1 but rank\\(G\\) = 2.*rank 20")
})

test_that("design anchor rank and covariance never form an m x m matrix", {
  set.seed(41)
  X  <- matrix(stats::rbinom(100 * 30, 2, 0.4), 100, 30)
  Xc <- sweep(X, 2, colMeans(X), "-")
  a  <- .qtl_anchor_design(Xc, 99)
  B  <- matrix(stats::rnorm(60), 30, 2)
  expect_equal(a$cov(B), crossprod(B, stats::cov(X) %*% B), tolerance = 1e-12)
  expect_equal(a$rank(1e-10), qr(Xc)$rank)
})
