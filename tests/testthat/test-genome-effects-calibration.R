# Step-5b internals (plans/import_qtl_effect_methods_phase_5_plan.md, 5b.9):
# the source project's non-additive suite ported to `.na_calibrate()` and its
# helpers, plus gates G4 (the dominance-mean solver) and G6 (per-trait units).
# No database. Architectures are drawn here and supplied, the source suite's
# test-18 path; the source generator (helper-nonadd-generator-oracle.R) is a
# cross-check on what must agree, never on coefficients.

relF <- function(A, B) norm(as.matrix(A - B), "F") / max(norm(as.matrix(B), "F"), 1e-300)

# Founders with LD and some inbred rows (the source suite's mk_panel()).
na_panel <- function(n, m, n_anc = 60, F_extra = 0) {
  anc <- t(replicate(n_anc, stats::rbinom(m, 1, stats::runif(m, 0.1, 0.9))))
  hap <- function() {
    h <- anc[sample(n_anc, 1), ]
    for (b in sort(sample(2:m, 4))) { s <- anc[sample(n_anc, 1), ]; h[b:m] <- s[b:m] }
    h
  }
  H1 <- t(replicate(n, hap())); H2 <- t(replicate(n, hap()))
  if (F_extra > 0) { sel <- stats::runif(n) < F_extra; H2[sel, ] <- H1[sel, ] }
  X <- H1 + H2
  X[, colMeans(X) > 0 & colMeans(X) < 2, drop = FALSE]
}

set.seed(1)
NA_X <- na_panel(300, 120, F_extra = 0.3)
NA_M <- ncol(NA_X)
NA_PAIRS <- matrix(sample(NA_M, 60), ncol = 2)

# Supplied architectures, drawn as the source draws them.
na_draw <- function(m, k, r) {
  list(B_a = matrix(stats::rnorm(m * k), m, k),
       z   = matrix(stats::rnorm(m * k), m, k),
       B_aa = if (r > 0) matrix(stats::rnorm(r * k), r, k))
}

na_anchors_for <- function(X, anchor, pairs = NULL) {
  a <- .na_anchors(anchor, p = colMeans(X) / 2, X = X)
  if (!is.null(pairs)) a <- .na_aa_anchor(a, pairs)
  a
}

named <- function(G) {
  G <- as.matrix(G)
  nm <- paste0("T", seq_len(nrow(G)))
  dimnames(G) <- list(nm, nm)
  G
}

run_na <- function(X, G_A, G_D = NULL, G_AA = NULL, pairs = NULL,
                   anchor = "genic", dr = NULL, ...) {
  k <- nrow(as.matrix(G_A))
  if (is.null(dr)) dr <- na_draw(ncol(X), k, if (is.null(G_AA)) 0 else nrow(pairs))
  .na_calibrate(na_anchors_for(X, anchor, if (!is.null(G_AA)) pairs),
                named(G_A), if (!is.null(G_D)) named(G_D),
                if (!is.null(G_AA)) named(G_AA), B_a = dr$B_a, z = dr$z,
                B_aa = if (!is.null(G_AA)) dr$B_aa, pairs = pairs, ...)
}

# The functional model's decomposition under the oracle's own covariates.
oracle_blocks <- function(X, cal, pairs, anchor) {
  dec <- nonadd_oracle$nonadd_decompose(X, cal$B_a, cal$B_d, cal$B_aa,
                                        if (!is.null(cal$B_aa)) pairs)
  if (anchor == "realised") {
    list(A = dec$real_A, D = dec$real_D, AA = dec$real_AA, id = dec$identity_error)
  } else {
    list(A = dec$genic_A, D = dec$genic_D, AA = dec$genic_AA, id = dec$identity_error)
  }
}

rand_pd <- function(k, s = 1) s * (crossprod(matrix(stats::rnorm(k * k), k)) + diag(k))


test_that("1. three blocks exact under both anchors, k = 1, 2, 3, and the oracle agrees", {
  set.seed(101)
  for (k in 1:3) for (anc in c("genic", "realised")) {
    GA <- rand_pd(k); GD <- rand_pd(k, 0.3); GAA <- rand_pd(k, 0.2)
    dr <- na_draw(NA_M, k, nrow(NA_PAIRS))
    cal <- run_na(NA_X, GA, GD, GAA, NA_PAIRS, anc, dr)
    expect_lt(relF(cal$delivered$A, named(GA)), 1e-10)
    expect_lt(relF(cal$delivered$D, named(GD)), 1e-10)
    expect_lt(relF(cal$delivered$AA, named(GAA)), 1e-10)
    # The source's own measuring instrument on our functional (a, d, e).
    ob <- oracle_blocks(NA_X, cal, NA_PAIRS, anc)
    expect_lt(relF(ob$A, GA), 1e-10)
    expect_lt(relF(ob$D, GD), 1e-10)
    expect_lt(relF(ob$AA, GAA), 1e-10)
    expect_lt(ob$id, 1e-10)
    # The source on the same supplied architectures delivers the same targets.
    src <- suppressWarnings(nonadd_oracle$sim_qtl_effects_nonadd(
      NA_X, GA, GD, GAA, NA_PAIRS, anchor = anc, B_a = dr$B_a,
      B_d = (0.19 + 0.097 * dr$z) * abs(dr$B_a), B_aa = dr$B_aa))
    expect_lt(relF(src$delivered$G_A, cal$delivered$A), 1e-10)
    expect_lt(relF(src$delivered$G_D, cal$delivered$D), 1e-10)
  }
})

test_that("2 / C6. k = 1 genic reproduces Zeng et al. (2013) Appendix A", {
  set.seed(102)
  p  <- stats::runif(300, 0.1, 0.9)
  Xh <- matrix(stats::rbinom(500 * 300, 2, rep(p, each = 500)), 500, 300)
  a0 <- stats::rnorm(300); z <- stats::rnorm(300)
  d0 <- (0.19 + 0.097 * z) * abs(a0)
  zz <- nonadd_oracle$zeng_appendix_A(colMeans(Xh) / 2, a0, d0, V_A = 4, V_D = 1)
  cal <- .na_calibrate(.na_anchors("genic", p = colMeans(Xh) / 2), named(4),
                       named(1), B_a = matrix(a0), z = matrix(z))
  expect_lt(relF(cal$B_a, matrix(zz$a)), 1e-12)
  expect_lt(relF(cal$B_d, matrix(zz$d)), 1e-12)
})

test_that("3 / C5. G_A below the floor is refused, naming the floor and the anchor", {
  set.seed(103)
  for (anc in c("genic", "realised")) {
    err <- tryCatch(run_na(NA_X, 1e-6, 5, anchor = anc), error = identity)
    expect_s3_class(err, "error")
    expect_match(conditionMessage(err), "below the additive floor for this sampled architecture")
    expect_match(conditionMessage(err), paste0("\"", anc, "\" anchor"))
    expect_match(conditionMessage(err), "at least T1 = [0-9.e+-]+")
    expect_match(conditionMessage(err), "Raise `G_A` or lower `G_D` / `G_AA`")
  }
})

test_that("4 / C4 (a). No non-zero D or A x A is Part A's congruence, nothing else", {
  set.seed(104)
  GA <- named(matrix(c(4, 2.4, 2.4, 9), 2))
  for (anc in c("genic", "realised")) {
    dr <- na_draw(NA_M, 2, 0)
    an <- na_anchors_for(NA_X, anc)
    cal <- .na_calibrate(an, GA, B_a = dr$B_a)
    ref <- .qtl_calibrate(dr$B_a, .qtl_target_std(unname(GA)), an$A)$B
    expect_identical(cal$route, "additive")
    expect_identical(cal$B_alpha, ref)
    # The plain congruence on the same normalised architecture agrees to 1e-12.
    v0 <- diag(an$A$cov(dr$B_a))
    cg <- .qtl_congruence(sweep(dr$B_a, 2, sqrt(v0), "/"),
                          .qtl_target_std(unname(GA))$eigen, an$A)$B %*%
      diag(sqrt(diag(GA)))
    expect_lt(relF(cal$B_alpha, cg), 1e-12)
    # Explicit zero D and A x A blocks take the same route and produce zeros.
    z0 <- matrix(0, 2, 2)
    cal0 <- .na_calibrate(.na_aa_anchor(an, NA_PAIRS), GA, named(z0), named(z0),
                          B_a = dr$B_a, z = dr$z, B_aa = matrix(1, nrow(NA_PAIRS), 2),
                          pairs = NA_PAIRS)
    expect_identical(cal0$B_alpha, ref)
    expect_true(all(cal0$B_d == 0) && all(cal0$B_aa == 0))
  }
})

test_that("5 / C7 / D4. Inbreeding depression: exact for one trait, reported for several", {
  set.seed(105)
  p  <- stats::runif(300, 0.1, 0.9)
  Xh <- matrix(stats::rbinom(500 * 300, 2, rep(p, each = 500)), 500, 300)
  for (anc in c("genic", "realised")) for (id in c(3, -2)) {
    cal <- run_na(Xh, 4, 1, anchor = anc,
                  inbreeding_depression = c(T1 = id))
    expect_lt(abs(cal$delivered$A[1, 1] - 4), 1e-10)
    expect_lt(abs(cal$delivered$D[1, 1] - 1), 1e-10)
    expect_lt(abs(cal$inbreeding$delivered[["T1"]] - id), 1e-10)
  }
  # k = 2: the blocks stay exact, the depression is approximate (D4 (a)).
  set.seed(106)
  cal <- run_na(Xh, diag(c(4, 6)), diag(c(1, 0.5)),
                inbreeding_depression = c(T1 = 3, T2 = 1))
  expect_lt(relF(cal$delivered$A, named(diag(c(4, 6)))), 1e-12)
  expect_lt(relF(cal$delivered$D, named(diag(c(1, 0.5)))), 1e-12)
  expect_true(all(is.finite(cal$inbreeding$delivered)))
  expect_gt(max(abs(cal$inbreeding$delivered - c(3, 1))), 1e-8)
  expect_error(run_na(Xh, 4, 1, inbreeding_depression = c(T1 = 1e6)),
               "inbreeding_depression` for trait 'T1'")
})

test_that("5 / G5. Only the named traits' degree means are solved", {
  set.seed(107)
  dr <- na_draw(NA_M, 2, 0)
  an <- na_anchors_for(NA_X, "genic")
  GA <- named(diag(c(4, 6))); GD <- named(diag(c(1, 0.5)))
  one <- .na_calibrate(an, GA, GD, B_a = dr$B_a, z = dr$z,
                       inbreeding_depression = c(T2 = 1))
  none <- .na_calibrate(an, GA, GD, B_a = dr$B_a, z = dr$z)
  expect_true(is.na(one$inbreeding$mean[["T1"]]))
  expect_false(is.na(one$inbreeding$mean[["T2"]]))
  # Before the joint calibration, trait 1's column is the drawn degrees.
  u <- abs(dr$B_a)
  B_d1 <- (0.19 + 0.097 * dr$z[, 1]) * u[, 1]
  expect_identical(one$inbreeding$mean[["T1"]], NA_real_)
  expect_false(identical(one$B_d, none$B_d))
  expect_true(is.finite(sum(B_d1)))
})

test_that("6 / D3 (a). Fixed loci keep their effects and add no variance under the anchor", {
  set.seed(108)
  Xf <- cbind(NA_X, 0, 2)
  m <- ncol(Xf)
  pr <- rbind(NA_PAIRS, c(1, m))          # a pair with a fixed member
  GA <- matrix(c(4, 2.4, 2.4, 9), 2)
  for (anc in c("genic", "realised")) {
    dr <- na_draw(m, 2, nrow(pr))
    cal <- run_na(Xf, GA, 0.3 * GA, 0.2 * GA, pr, anc, dr)
    expect_true(all(cal$B_a[(m - 1):m, ] != 0))
    expect_true(all(cal$B_d[(m - 1):m, ] != 0))
    expect_true(all(cal$B_aa[nrow(pr), ] != 0))
    expect_lt(relF(cal$delivered$A, named(GA)), 1e-10)
    # Per locus and per pair, never the whole model: the fixed loci's own
    # contrast columns and the pair containing one carry no variance.
    an <- na_anchors_for(Xf, anc, pr)
    w_or_var <- function(a) a$locus_variance
    expect_true(all(w_or_var(an$A)[(m - 1):m] == 0))
    expect_true(all(w_or_var(an$D)[(m - 1):m] == 0))
    expect_true(w_or_var(an$AA)[nrow(pr)] == 0)
  }
})

test_that("9. Rank-one D and A x A targets through the whole calibration", {
  set.seed(109)
  for (anc in c("genic", "realised")) {
    cal <- run_na(NA_X, diag(2) + 1, 0.2 * tcrossprod(c(1, 1)),
                  0.1 * tcrossprod(c(1, -2)), NA_PAIRS, anc)
    expect_lt(relF(cal$delivered$A, named(diag(2) + 1)), 1e-10)
    expect_lt(relF(cal$delivered$D, named(0.2 * tcrossprod(c(1, 1)))), 1e-10)
    expect_lt(relF(cal$delivered$AA, named(0.1 * tcrossprod(c(1, -2)))), 1e-10)
  }
})

test_that("10. A rank-deficient additive architecture is refused by the stage", {
  set.seed(110)
  b <- stats::rnorm(NA_M)
  an <- na_anchors_for(NA_X, "genic")
  C <- matrix(stats::rnorm(NA_M * 2), NA_M, 2)
  expect_error(.na_additive_stage(cbind(b, 2 * b), C,
                                  .qtl_target_std(diag(2)), an$A, "genic",
                                  c("T1", "T2")),
               "drawn additive architecture has rank 1 under the \"genic\" anchor but there are 2 traits")
})

test_that("11. Target scale invariance over 24 orders of magnitude", {
  set.seed(111)
  GA <- matrix(c(4, 2.4, 2.4, 9), 2); GD <- matrix(c(1, .3, .3, 2), 2)
  GAA <- matrix(c(.5, .1, .1, .8), 2)
  dr <- na_draw(NA_M, 2, nrow(NA_PAIRS))
  worst <- 0
  for (s in c(1e-12, 1e-6, 1, 1e6, 1e12)) for (anc in c("genic", "realised")) {
    cal <- run_na(NA_X, s * GA, s * GD, s * GAA, NA_PAIRS, anc, dr)
    worst <- max(worst, relF(cal$delivered$A, named(s * GA)),
                 relF(cal$delivered$D, named(s * GD)),
                 relF(cal$delivered$AA, named(s * GAA)))
  }
  expect_lt(worst, 1e-10)
})

test_that("12. Randomised sweep over k, anchor, blocks and ranks", {
  set.seed(112)
  worst <- 0
  for (i in 1:16) {
    k <- sample(1:3, 1); anc <- sample(c("genic", "realised"), 1)
    GA <- rand_pd(k)
    rd <- sample(0:k, 1)
    GD <- if (stats::runif(1) < 0.7) 0.3 * tcrossprod(matrix(stats::rnorm(k * rd), k))
    ra <- sample(1:k, 1)
    GAA <- if (stats::runif(1) < 0.7) 0.2 * tcrossprod(matrix(stats::rnorm(k * ra), k))
    cal <- run_na(NA_X, GA, GD, GAA, if (!is.null(GAA)) NA_PAIRS, anc)
    e <- relF(cal$delivered$A, named(GA))
    if (!is.null(GD)) e <- max(e, if (rd == 0) max(abs(cal$delivered$D)) else
      relF(cal$delivered$D, named(GD)))
    if (!is.null(GAA)) e <- max(e, relF(cal$delivered$AA, named(GAA)))
    worst <- max(worst, e)
  }
  expect_lt(worst, 1e-10)
})

test_that("13 / D3 (a). All-heterozygous loci: segregating under genic, no variance under realised", {
  set.seed(113)
  Xf <- NA_X[, 1:60]; Xf[, 1:20] <- 1
  GA <- matrix(c(4, 2.4, 2.4, 9), 2); GD <- matrix(c(1, .3, .3, 2), 2)
  rg <- run_na(Xf, GA, GD, anchor = "genic")
  expect_true(all(rowSums(rg$B_a[1:20, ]^2) > 0))
  expect_lt(relF(rg$delivered$A, named(GA)), 1e-10)
  rr <- run_na(Xf, GA, GD, anchor = "realised")
  an <- na_anchors_for(Xf, "realised")
  expect_true(all(an$A$locus_variance[1:20] == 0))
  expect_true(all(an$D$locus_variance[1:20] == 0))
  expect_true(all(rowSums(rr$B_a[1:20, ]^2) > 0))
  expect_lt(relF(rr$delivered$A, named(GA)), 1e-10)
  expect_lt(relF(rr$delivered$D, named(GD)), 1e-10)
})

test_that("14. Zero-rank targets", {
  set.seed(114)
  GA <- matrix(c(4, 2.4, 2.4, 9), 2)
  z0 <- matrix(0, 2, 2)
  r0 <- run_na(NA_X, GA, z0, z0, NA_PAIRS)
  expect_true(max(abs(r0$B_d)) == 0 && max(abs(r0$B_aa)) == 0)
  expect_lt(relF(r0$delivered$A, named(GA)), 1e-10)
  rz <- run_na(NA_X, z0)
  expect_true(max(abs(rz$B_a)) == 0)
  r1 <- run_na(NA_X, GA, 0.2 * tcrossprod(c(1, 1)))
  expect_lt(relF(r1$delivered$D, named(0.2 * tcrossprod(c(1, 1)))), 1e-10)
  expect_lt(relF(r1$delivered$A, named(GA)), 1e-10)
})

test_that("15. Trait names on every delivered block", {
  set.seed(115)
  cal <- run_na(NA_X, diag(2), 0.3 * diag(2), 0.2 * diag(2), NA_PAIRS)
  for (b in c("A", "D", "AA")) {
    expect_identical(dimnames(cal$delivered[[b]]), list(c("T1", "T2"), c("T1", "T2")))
  }
  expect_identical(names(cal$inbreeding$delivered), c("T1", "T2"))
})

test_that("17. Edge sizes: one locus with A + D, one pair", {
  set.seed(117)
  r1 <- run_na(NA_X[, 1, drop = FALSE], 2, 0.5)
  expect_lt(abs(r1$delivered$A[1, 1] - 2), 1e-10)
  expect_lt(abs(r1$delivered$D[1, 1] - 0.5), 1e-10)
  r3 <- run_na(NA_X[, 1:2], 1, NULL, 0.1, matrix(1:2, 1))
  expect_lt(abs(r3$delivered$AA[1, 1] - 0.1), 1e-10)
  expect_lt(abs(r3$delivered$A[1, 1] - 1), 1e-10)
})

test_that("18. The functional additive effects stay in the supplied architecture's span", {
  set.seed(118)
  dr <- na_draw(NA_M, 2, 0)
  cal <- run_na(NA_X, matrix(c(4, 2.4, 2.4, 9), 2), 0.3 * diag(2), dr = dr)
  T_hat <- qr.solve(dr$B_a, cal$B_a)
  expect_lt(max(abs(cal$B_a - dr$B_a %*% T_hat)), 1e-10 * max(abs(cal$B_a)))
})

test_that("Anchor objects: cross() is the dense M, and the design rank is cached", {
  set.seed(119)
  B1 <- matrix(stats::rnorm(NA_M * 2), NA_M); B2 <- matrix(stats::rnorm(NA_M * 3), NA_M)
  w <- stats::runif(NA_M)
  expect_lt(relF(.qtl_anchor_diag(w)$cross(B1, B2), crossprod(B1, diag(w) %*% B2)), 1e-13)
  Xc <- sweep(NA_X, 2, colMeans(NA_X))
  an <- .qtl_anchor_design(Xc, nrow(Xc) - 1)
  expect_lt(relF(an$cross(B1, B2), crossprod(B1, crossprod(Xc) %*% B2) / (nrow(Xc) - 1)), 1e-12)
  expect_identical(an$rank(1e-10), an$rank(1e-10))
  expect_identical(an$rank(1e-10), qr(Xc)$rank)
})

test_that("The coupling accumulates hubs exactly like the source's loop", {
  set.seed(120)
  m <- 10; k <- 2
  pr <- rbind(c(1, 2), c(1, 3), c(1, 4), c(5, 6))
  B_d <- matrix(stats::rnorm(m * k), m); B_aa <- matrix(stats::rnorm(nrow(pr) * k), nrow(pr))
  b <- stats::runif(m, -1, 1); cc <- stats::runif(m, -1, 1)
  ref <- b * B_d
  for (r in seq_len(nrow(pr))) {
    ref[pr[r, 1], ] <- ref[pr[r, 1], ] + cc[pr[r, 2]] * B_aa[r, ]
    ref[pr[r, 2], ] <- ref[pr[r, 2], ] + cc[pr[r, 1]] * B_aa[r, ]
  }
  expect_lt(max(abs(.na_coupling(m, k, b, cc, B_d, B_aa, pr) - ref)), 1e-14)
})


# -- G4: the dominance-degree mean solver -------------------------------------

dd_check <- function(r, u, v, w, M, rho, sd) {
  d  <- r$mu * u + sd * v
  id <- sum(w * d); vd <- M$cov(matrix(d))[1, 1]
  expect_true(is.finite(r$mu))
  expect_gt(vd, 0)
  expect_lt(abs(id / sqrt(vd) - rho), 1e-8 * max(1, abs(rho)))
}

test_that("G4. Degenerate coefficients: one genic locus, rho = +/-1", {
  p <- 0.3; w <- 2 * p * (1 - p); M <- .qtl_anchor_diag(w^2)
  for (rho in c(1, -1)) {
    r <- .na_solve_dd_mean(1, 0.5, w, M, rho, 0.097, 0.19)
    expect_identical(r$branch, "degenerate")
    dd_check(r, 1, 0.5, w, M, rho, 0.097)
  }
})

test_that("G4. A = 0 with B != 0 takes the linear root (mu = 0.1)", {
  # A = rho^2 u'Mu - (w'u)^2 vanishes when rho is the ratio of d = u itself
  # (x -> infinity). With w_1 = w_2 the ratio is symmetric in loci 1 and 2,
  # so d = (2, 1, 1) has the ratio of u = (1, 2, 1) without being
  # proportional to it: at mu = 0.1, sd = 0.1 that is v = (1, -1, 0).
  p <- c(0.2, 0.8, 0.5); w <- 2 * p * (1 - p); M <- .qtl_anchor_diag(w^2)
  u <- c(1, 2, 1); v <- c(1, -1, 0)
  rho <- sum(w * u) / sqrt(sum(w^2 * u^2))
  r <- .na_solve_dd_mean(u, v, w, M, rho, 0.1, 0.19)
  expect_identical(r$branch, "linear")
  expect_equal(r$mu, 0.1, tolerance = 1e-12)
  dd_check(r, u, v, w, M, rho, 0.1)
})

test_that("G4. B = 0 two-locus genic case: rho = +/-1 gives mu = +/-0.1", {
  p <- c(0.5, 0.5); w <- 2 * p * (1 - p); M <- .qtl_anchor_diag(w^2)
  r <- .na_solve_dd_mean(c(1, 1), c(1, -1), w, M, 1, 0.1, 0.19)
  expect_equal(r$mu, 0.1, tolerance = 1e-12)
  r <- .na_solve_dd_mean(c(1, 1), c(1, -1), w, M, -1, 0.1, 0.19)
  expect_equal(r$mu, -0.1, tolerance = 1e-12)
})

test_that("G4. The preferred mean already on ID = 0 moves one sd past the root", {
  p <- 0.3; w <- 2 * p * (1 - p); M <- .qtl_anchor_diag(w^2)
  r <- .na_solve_dd_mean(1, -1.9, w, M, 1, 0.1, 0.19)
  expect_equal(r$mu, 0.29, tolerance = 1e-12)
  dd_check(r, 1, -1.9, w, M, 1, 0.1)
  r <- .na_solve_dd_mean(1, -1.9, w, M, -1, 0.1, 0.19)
  expect_equal(r$mu, 0.09, tolerance = 1e-12)
})

dd_fixture <- function() {
  set.seed(121)
  u <- abs(stats::rnorm(20)); v <- stats::rnorm(20) * u
  p <- stats::runif(20, 0.1, 0.9); w <- 2 * p * (1 - p)
  M <- .qtl_anchor_diag(w^2)
  Q <- M$cov(cbind(u, v)); g <- c(sum(w * u), sum(w * v))
  cf <- solve(Q, g)
  # The line d = x u + v reaches the Cauchy-Schwarz bound on one side only:
  # at direction (cf1, cf2), which has positive weight on v iff cf2 > 0.
  list(u = u, v = v, w = w, M = M, Q = Q, g = g, cf = cf,
       bound = sign(cf[2]) * sqrt(sum(g * cf)))
}

test_that("G4. A repeated root at x = 0 with positive V_D there", {
  f <- dd_fixture()
  v0 <- f$v + (f$cf[1] / f$cf[2]) * f$u    # moves the extremum to x = 0
  r <- .na_solve_dd_mean(f$u, v0, f$w, f$M, f$bound, 0.097, 0.19)
  expect_identical(r$branch, "quadratic")
  # A double root's location is ill-conditioned (the ratio is flat there), so
  # only its neighbourhood is asserted; the ratio itself is checked exactly.
  expect_lt(abs(r$x), 1e-3)
  dd_check(r, f$u, v0, f$w, f$M, f$bound, 0.097)
})

test_that("G4. rho = 0 whose only root has zero dominance variance is refused", {
  p <- c(0.3, 0.6); w <- 2 * p * (1 - p); M <- .qtl_anchor_diag(w^2)
  # d = sd (x u + v) with v = -u vanishes at x = 1, the root of ID.
  expect_error(.na_solve_dd_mean(c(1, 2), -c(1, 2), w, M, 0, 0.1, 0.19),
               "positive dominance variance")
})

test_that("G4. rho = 0, negative rho, and rho just inside / outside the bound", {
  f <- dd_fixture()
  for (rho in c(0, -1, 1, f$bound * (1 - 1e-6))) {
    r <- .na_solve_dd_mean(f$u, f$v, f$w, f$M, rho, 0.097, 0.19)
    dd_check(r, f$u, f$v, f$w, f$M, rho, 0.097)
  }
  expect_error(.na_solve_dd_mean(f$u, f$v, f$w, f$M, f$bound * (1 + 1e-6), 0.097, 0.19),
               "Cauchy-Schwarz bound")
})

test_that("G4. An unattainable one-locus request is refused", {
  p <- 0.3; w <- 2 * p * (1 - p); M <- .qtl_anchor_diag(w^2)
  expect_error(.na_solve_dd_mean(1, 0.5, w, M, 2, 0.097, 0.19),
               "no dominance-degree mean")
})

test_that("G4. The solution in x = mu / sd does not depend on sd", {
  f <- dd_fixture()
  out <- lapply(10^seq(-8, 8, 4), function(sd) {
    .na_solve_dd_mean(f$u, f$v, f$w, f$M, 1, sd, 0.19 * sd)
  })
  expect_true(all(vapply(out, function(r) r$branch, "") == "quadratic"))
  x <- vapply(out, function(r) r$x, 0)
  expect_lt(max(abs(x - x[3])), 1e-10 * abs(x[3]))
})


# -- G6: per-trait units and the floor boundary -------------------------------

test_that("G6. Per-trait units: S = diag(1e-6, 1e6) transforms consistently", {
  set.seed(130)
  GA <- matrix(c(4, 2.4, 2.4, 9), 2); GD <- matrix(c(1, .3, .3, 2), 2)
  GAA <- matrix(c(.5, .1, .1, .8), 2)
  S <- diag(c(1e-6, 1e6))
  dr <- na_draw(NA_M, 2, nrow(NA_PAIRS))
  ds <- list(B_a = dr$B_a %*% S, z = dr$z, B_aa = dr$B_aa %*% S)
  for (anc in c("genic", "realised")) {
    c0 <- run_na(NA_X, GA, GD, GAA, NA_PAIRS, anc, dr)
    c1 <- run_na(NA_X, S %*% GA %*% S, S %*% GD %*% S, S %*% GAA %*% S,
                 NA_PAIRS, anc, ds)
    std <- .qtl_target_std(S %*% GA %*% S)
    expect_lt(.qtl_target_error(c1$delivered$A, std), 1e-10)
    # The small-unit trait is preserved, and the effects transform with S.
    expect_lt(relF(c1$B_a[, 1], c0$B_a[, 1] * 1e-6), 1e-8)
    expect_lt(relF(c1$B_a[, 2], c0$B_a[, 2] * 1e6), 1e-8)
  }
  # Same feasibility decision: a target below the floor is refused in both units.
  expect_error(run_na(NA_X, 1e-6 * GA, GD, anchor = "genic", dr = dr), "floor")
  expect_error(run_na(NA_X, S %*% (1e-6 * GA) %*% S, S %*% GD %*% S,
                      anchor = "genic", dr = ds), "floor")
})

test_that("G6. Targets on, just above, and below the independently computed floor", {
  set.seed(131)
  dr <- na_draw(NA_M, 2, nrow(NA_PAIRS))
  GD <- matrix(c(4, 1, 1, 3), 2); GAA <- matrix(c(2, .5, .5, 2), 2)
  big <- run_na(NA_X, diag(c(100, 100)), GD, GAA, NA_PAIRS, "genic", dr)
  # The floor from the residual coefficients, computed here.
  p <- colMeans(NA_X) / 2; w <- 2 * p * (1 - p)
  C <- (1 - 2 * p) * big$B_d
  for (r in seq_len(nrow(NA_PAIRS))) {
    k <- NA_PAIRS[r, 1]; l <- NA_PAIRS[r, 2]
    C[k, ] <- C[k, ] + (2 * p[l] - 1) * big$B_aa[r, ]
    C[l, ] <- C[l, ] + (2 * p[k] - 1) * big$B_aa[r, ]
  }
  M <- diag(w); B <- dr$B_a
  Cres <- C - B %*% solve(t(B) %*% M %*% B, t(B) %*% M %*% C)
  floor <- t(Cres) %*% M %*% Cres
  floor <- (floor + t(floor)) / 2
  on <- run_na(NA_X, floor, GD, GAA, NA_PAIRS, "genic", dr)
  expect_lt(relF(on$delivered$A, named(floor)), 1e-8)
  above <- run_na(NA_X, floor * (1 + 1e-6), GD, GAA, NA_PAIRS, "genic", dr)
  expect_lt(relF(above$delivered$A, named(floor * (1 + 1e-6))), 1e-10)
  expect_error(run_na(NA_X, floor * 0.9, GD, GAA, NA_PAIRS, "genic", dr),
               "below the additive floor")
  # Per-trait units: trait 1 in tiny units below its own floor, trait 2 in
  # huge units far above it. A floor test on the raw scale sees only trait 2.
  low <- diag(c(0.5 * floor[1, 1], 100 * floor[2, 2]))
  S <- diag(c(1e-6, 1e6))
  ds <- list(B_a = dr$B_a %*% S, z = dr$z, B_aa = dr$B_aa %*% S)
  expect_error(run_na(NA_X, low, GD, GAA, NA_PAIRS, "genic", dr),
               "below the additive floor")
  expect_error(run_na(NA_X, S %*% low %*% S, S %*% GD %*% S, S %*% GAA %*% S,
                      NA_PAIRS, "genic", ds),
               "below the additive floor")
})

test_that("G6. One trait, G_A = 4, original-unit floor 1: the standardised residual is 0.75", {
  M <- .qtl_anchor_diag(c(1, 1))
  st <- .na_additive_stage(matrix(c(1, 0)), matrix(c(0, 1)),
                           .qtl_target_std(matrix(4)), M, "genic", "T1")
  expect_equal(st$floor[1, 1], 1, tolerance = 1e-14)
  expect_equal(st$floor_s[1, 1], 0.25, tolerance = 1e-14)
  expect_equal(st$G_tilde_min, 0.75, tolerance = 1e-14)
  expect_equal(M$cov(st$B_alpha)[1, 1], 4, tolerance = 1e-14)
})


test_that("The source oracle is isolated and runs on its own", {
  set.seed(140)
  fns <- ls(nonadd_oracle, all.names = TRUE)
  expect_true(all(vapply(fns, function(f) {
    identical(environment(get(f, envir = nonadd_oracle)), nonadd_oracle)
  }, logical(1))))
  expect_identical(parent.env(nonadd_oracle), baseenv())
  src <- suppressWarnings(nonadd_oracle$sim_qtl_effects_nonadd(
    NA_X[, 1:40], diag(2), 0.3 * diag(2), 0.2 * diag(2), matrix(1:20, ncol = 2)))
  expect_lt(relF(src$delivered$G_D, 0.3 * diag(2)), 1e-10)
})
