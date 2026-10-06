# Test oracle for gate B2 of plans/import_qtl_effect_methods_phase_4_plan.md
# (decision D4): the source project's measuring instrument, copied verbatim
# so the test runs on every machine. Source:
#   /Users/austinputz/Claude/simulate_qtl_effects
#   non-additive/R/qtl_effects_nonadd.R, nonadd_covariates() and
#   nonadd_decompose(), unchanged from commit 8f8a97c (checked 2026-10-05).
# Do not edit these two functions: they are the independent reference that
# extract_genetic_variance() is checked against. They share no code with the
# package's conversion or projection helpers.

## Covariates and anchors. Returns the NOIA design blocks and the per-component anchor matrices.
nonadd_covariates <- function(X, pairs = NULL, anchor = c("genic", "realised")) {
  anchor <- match.arg(anchor)
  n <- nrow(X); m <- ncol(X)
  p <- colMeans(X) / 2; q <- 1 - p
  W <- (X == 1) * 1
  Z_A <- sweep(X, 2, 2 * p)
  vx  <- colSums(Z_A^2) / n                      # observed dosage variance (population form)
  cwx <- colSums(sweep(W, 2, colMeans(W)) * Z_A) / n
  if (anchor == "realised") {
    b <- ifelse(vx > 0, cwx / vx, 0)
  } else {
    b <- q - p                                    # HWE value of Cov(w,x)/Var(x)
  }
  Z_D <- sweep(W, 2, colMeans(W)) - sweep(Z_A, 2, b, "*")
  cc  <- 2 * p - 1                                # E[x - 1] = p - q
  Z_AA <- NULL
  if (!is.null(pairs)) {
    Z_AA <- Z_A[, pairs[, 1], drop = FALSE] * Z_A[, pairs[, 2], drop = FALSE]
    Z_AA <- sweep(Z_AA, 2, colMeans(Z_AA))
  }
  if (anchor == "genic") {
    M_A  <- diag(2 * p * q, m)
    M_D  <- diag((2 * p * q)^2, m)
    M_AA <- if (!is.null(pairs)) diag(4 * p[pairs[, 1]] * q[pairs[, 1]] * p[pairs[, 2]] * q[pairs[, 2]], nrow(pairs)) else NULL
  } else {
    M_A  <- crossprod(Z_A) / (n - 1)
    M_D  <- crossprod(Z_D) / (n - 1)
    M_AA <- if (!is.null(pairs)) crossprod(Z_AA) / (n - 1) else NULL
  }
  list(p = p, b = b, c = cc, w_id = 2 * p * q, Z_A = Z_A, Z_D = Z_D, Z_AA = Z_AA,
       M_A = M_A, M_D = M_D, M_AA = M_AA)
}

## Exact NOIA decomposition of the functional model on a genotype matrix: components,
## realised (by-individual) and genic (by-locus, HWE+LE) variance matrices, alpha, and the
## inbreeding depression. This is the measuring instrument every comparison uses.
nonadd_decompose <- function(X, B_a, B_d = NULL, B_aa = NULL, pairs = NULL) {
  X <- as.matrix(X); n <- nrow(X); m <- ncol(X); k <- ncol(B_a)
  cv <- nonadd_covariates(X, pairs, anchor = "realised")
  C <- matrix(0, m, k)
  if (!is.null(B_d))  C <- C + cv$b * B_d
  if (!is.null(B_aa)) for (r in seq_len(nrow(pairs))) {
    C[pairs[r, 1], ] <- C[pairs[r, 1], ] + cv$c[pairs[r, 2]] * B_aa[r, ]
    C[pairs[r, 2], ] <- C[pairs[r, 2], ] + cv$c[pairs[r, 1]] * B_aa[r, ]
  }
  B_alpha <- B_a + C
  BV <- cv$Z_A %*% B_alpha
  DD <- if (!is.null(B_d))  cv$Z_D  %*% B_d  else matrix(0, n, k)
  AA <- if (!is.null(B_aa)) cv$Z_AA %*% B_aa else matrix(0, n, k)
  ## functional genotypic value, computed directly, to prove the identity g = const + BV + DD + AA
  g <- (X - 1) %*% B_a
  if (!is.null(B_d))  g <- g + (X == 1) %*% B_d
  if (!is.null(B_aa)) g <- g + ((X[, pairs[, 1], drop = FALSE] - 1) * (X[, pairs[, 2], drop = FALSE] - 1)) %*% B_aa
  gc <- sweep(g, 2, colMeans(g))
  p <- cv$p; q <- 1 - p
  ## genic (HWE + LE) statistical variances at this population's allele frequencies
  bH <- q - p
  Ch <- matrix(0, m, k)
  if (!is.null(B_d))  Ch <- Ch + bH * B_d
  if (!is.null(B_aa)) for (r in seq_len(nrow(pairs))) {
    Ch[pairs[r, 1], ] <- Ch[pairs[r, 1], ] + cv$c[pairs[r, 2]] * B_aa[r, ]
    Ch[pairs[r, 2], ] <- Ch[pairs[r, 2], ] + cv$c[pairs[r, 1]] * B_aa[r, ]
  }
  alphaH <- B_a + Ch
  gen_A  <- crossprod(alphaH, (2 * p * q) * alphaH)
  gen_D  <- if (!is.null(B_d))  crossprod(B_d, (2 * p * q)^2 * B_d) else NULL
  gen_AA <- if (!is.null(B_aa)) crossprod(B_aa, (4 * p[pairs[, 1]] * q[pairs[, 1]] * p[pairs[, 2]] * q[pairs[, 2]]) * B_aa) else NULL
  list(BV = BV, DD = DD, AA = AA, g = g, alpha = B_alpha,
       identity_error = max(abs(gc - (BV + DD + AA))),
       real_A = stats::var(BV), real_D = stats::var(DD), real_AA = stats::var(AA),
       real_G = stats::var(g), cov_A_D = stats::cov(BV, DD), cov_A_AA = stats::cov(BV, AA),
       genic_A = gen_A, genic_D = gen_D, genic_AA = gen_AA,
       inbr_depr = if (!is.null(B_d)) as.vector(crossprod(2 * p * q, B_d)) else NULL,
       p = p)
}
