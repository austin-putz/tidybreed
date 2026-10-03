# Fixtures for the ported QTL-effect method suite (test-qtl-congruence.R).
# Ported from the source project's tests/test_qtl_effects_paper.R, where they
# were built inline; their seeds and shapes are unchanged.

# The paper suite's shared target and panel: 3 traits, 1200 x 30 dosages.
qtl_paper_G <- function() matrix(c(4, 2, 1, 2, 9, 3, 1, 3, 16), 3, 3)

qtl_paper_X <- function() {
  set.seed(1001)
  matrix(stats::rbinom(1200 * 30, 2, 0.3), 1200, 30)
}

qtl_breeding_values <- function(X, fit) sweep(X, 2, fit$centre, "-") %*% fit$B

# Haplotypes in 6 two-locus LD blocks (switch rate 0.08, every second block in
# repulsion). Draws from the current RNG stream: seed before calling.
qtl_mkhap <- function(nh, nb = 6) {
  do.call(cbind, lapply(seq_len(nb), function(b) {
    a  <- stats::rbinom(nh, 1, 0.5)
    bb <- ifelse(stats::rbinom(nh, 1, 0.08) == 1, 1 - a, a)
    if (b %% 2 == 0) bb <- 1 - bb
    cbind(a, bb)
  }))
}

# Matrix symmetric square root and inverse square root (full rank), for the
# dense oracle of gate A2.
qtl_sym_pow <- function(S, power) {
  e <- eigen((S + t(S)) / 2, symmetric = TRUE)
  e$vectors %*% diag(pmax(e$values, 0)^power, nrow = length(e$values)) %*%
    t(e$vectors)
}
