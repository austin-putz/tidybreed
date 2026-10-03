# Exact multivariate QTL-effect calibration: the congruence B = B0 A with
# B' M B = G for a named reference-population covariance M (the anchor).
#
# Ported from the source project `simulate_qtl_effects/R/qtl_effects.R` (the
# current file, never the frozen `qtl_effects_paper.R`); see
# plans/import_qtl_effect_methods.md §1.1, §1A and §7. The algebra is the
# source's, character for character, with the `.congruence` closure lifted to
# `.qtl_congruence()` so the generators can call it. Two changes:
#
# * A second rank error. The source raises one, `rank(B0' M B0) < rank(G)`.
#   That blames the drawn architecture even when the anchor itself has too few
#   independent directions to carry the target for *any* effects (Theorem 1).
#   `.qtl_congruence()` checks `rank(G) <= rank(M)` first, so the two causes
#   get distinct messages (gate A5).
# * The anchor is an object (`.qtl_anchor_diag()`, `.qtl_anchor_design()`)
#   rather than closure state, so a generator can build it from base allele
#   frequencies without any genotype matrix.
#
# `.qtl_sim_effects()` mirrors the source's `sim_qtl_effects()` without the
# architecture samplers and `marginal =` (deferred, Q5). It is internal and
# exists so the source's test suite ports with its fixtures unchanged; the
# public surface is `define_additive_effects(anchor = )`. The `dual`,
# `marginal_observed` and `reference` anchors are internal only (Q6).

#' @noRd
.qtl_max_abs <- function(x) if (length(x)) max(abs(x)) else 0

#' @noRd
.qtl_validate_numeric_matrix <- function(x, name, min_rows = 1L, min_cols = 1L) {
  if (!is.matrix(x) || !is.numeric(x))
    stop("`", name, "` must be a numeric matrix.", call. = FALSE)
  if (nrow(x) < min_rows || ncol(x) < min_cols)
    stop("`", name, "` must have at least ", min_rows, " row(s) and ",
         min_cols, " column(s).", call. = FALSE)
  if (anyNA(x) || any(!is.finite(x)))
    stop("`", name, "` must contain only finite, non-missing values.",
         call. = FALSE)
  invisible(TRUE)
}

#' @noRd
.qtl_validate_dosage <- function(x, ploidy, name, tol) {
  if (min(x) < -tol || max(x) > ploidy + tol)
    stop("`", name, "` must contain allele dosages in [0, ploidy].", call. = FALSE)
  invisible(TRUE)
}

#' @noRd
.qtl_validate_scalar <- function(x, name, integral = FALSE, min = -Inf, max = Inf) {
  if (!is.numeric(x) || length(x) != 1L || !is.finite(x))
    stop("`", name, "` must be a single finite numeric value.", call. = FALSE)
  if (integral && abs(x - round(x)) > .Machine$double.eps^0.5)
    stop("`", name, "` must be an integer.", call. = FALSE)
  if (x < min || x > max)
    stop("`", name, "` must lie in [", min, ", ", max, "].", call. = FALSE)
  invisible(TRUE)
}

#' @noRd
.qtl_validate_exact_dim <- function(x, name, nr, nc) {
  if (nrow(x) != nr || ncol(x) != nc)
    stop("`", name, "` must be exactly ", nr, " x ", nc, "; got ",
         nrow(x), " x ", ncol(x), ".", call. = FALSE)
  invisible(TRUE)
}

#' @noRd
.qtl_validate_warn_bounds <- function(warn_bounds) {
  if (is.null(warn_bounds)) return(invisible(TRUE))
  if (!is.numeric(warn_bounds) || length(warn_bounds) != 2L ||
      any(!is.finite(warn_bounds)))
    stop("`warn_bounds` must be NULL or two finite numeric values.", call. = FALSE)
  if (warn_bounds[1] <= 0 || warn_bounds[1] > warn_bounds[2])
    stop("`warn_bounds` must satisfy 0 < lower <= upper.", call. = FALSE)
  invisible(TRUE)
}

#' Eigen-decompose a matrix that must be positive semidefinite
#'
#' Negative eigenvalues within `rel_tol` of the largest absolute eigenvalue are
#' rounding and are clamped to zero; anything more negative is an error naming
#' the smallest eigenvalue. `rank` counts eigenvalues above `rel_tol` times the
#' largest.
#' @noRd
.qtl_psd_eigen <- function(x, name, rel_tol) {
  x <- (x + t(x)) / 2
  e <- eigen(x, symmetric = TRUE)
  scale <- if (length(e$values)) max(abs(e$values)) else 0
  if (scale > 0 && min(e$values) < -rel_tol * scale)
    stop("`", name, "` must be positive semi-definite; smallest eigenvalue is ",
         format(min(e$values), digits = 5), ".", call. = FALSE)
  values <- pmax(e$values, 0)
  rank <- if (scale == 0) 0L else sum(values > max(values) * rel_tol)
  list(values = values, vectors = e$vectors, rank = rank, scale = scale)
}

#' Validate a target covariance and decompose it
#'
#' Symmetry is checked relative to the matrix's own scale, as in the source.
#' @return `.qtl_psd_eigen()`'s list, plus `G` (the symmetrised input).
#' @noRd
.qtl_target_eigen <- function(G, name = "G", rel_tol = 1e-10) {
  .qtl_validate_numeric_matrix(G, name)
  if (nrow(G) != ncol(G)) stop("`", name, "` must be square.", call. = FALSE)
  g_scale <- .qtl_max_abs(G)
  if (.qtl_max_abs(G - t(G)) > rel_tol * max(g_scale, .Machine$double.eps))
    stop("`", name, "` must be symmetric within the relative tolerance.",
         call. = FALSE)
  G <- (G + t(G)) / 2
  c(.qtl_psd_eigen(G, name, rel_tol), list(G = G))
}

#' Per-contrast spectrum of a candidate covariance relative to the target
#'
#' The eigenvalues of `G^{-1/2} candidate G^{-1/2}` on the target's range. All
#' ones means the candidate equals the target; a trace ratio can hide a pair of
#' contrasts that miss in opposite directions, these cannot.
#' @noRd
.qtl_relative_spectrum <- function(target_eigen, candidate) {
  q <- target_eigen$rank
  if (q == 0L) return(numeric())
  V <- target_eigen$vectors[, seq_len(q), drop = FALSE]
  gamma <- target_eigen$values[seq_len(q)]
  W <- V %*% diag(1 / sqrt(gamma), nrow = q)
  R <- crossprod(W, candidate %*% W)
  sort(pmax(eigen((R + t(R)) / 2, symmetric = TRUE, only.values = TRUE)$values, 0))
}

#' @noRd
.qtl_spectrum_text <- function(x) {
  if (!length(x)) return("not defined (zero target)")
  sprintf("[%.3g, %.3g]", min(x), max(x))
}


# ---------------------------------------------------------------------------
# Anchors
# ---------------------------------------------------------------------------

#' A diagonal anchor: M = diag(w)
#'
#' The genic anchor is `w = n_eligible * p * q`. Never forms an m x m matrix.
#' @noRd
.qtl_anchor_diag <- function(w) {
  w <- as.numeric(w)
  list(kind = "diagonal", weights = w,
       cov = function(B) crossprod(B * sqrt(w)),
       locus_variance = w,
       rank = function(rel_tol) {
         if (!length(w) || .qtl_max_abs(w) == 0) 0L
         else sum(w > .qtl_max_abs(w) * rel_tol)
       })
}

#' A design anchor: M = Xc' Xc / denominator
#'
#' `Xc` is the column-centred n x m design. `cov()` is `O(nmk)`: it forms
#' `Xc %*% B`, never `Xc' Xc`. `rank()` uses the singular values of `Xc`
#' (squared, so the cut-off is on the eigenvalue scale `.qtl_psd_eigen()` uses).
#' @noRd
.qtl_anchor_design <- function(Xc, denominator) {
  list(kind = "design", design = Xc, denominator = denominator,
       cov = function(B) crossprod(Xc %*% B) / denominator,
       locus_variance = colSums(Xc^2) / denominator,
       rank = function(rel_tol) {
         if (!length(Xc)) return(0L)
         d2 <- svd(Xc, nu = 0L, nv = 0L)$d^2
         if (!length(d2) || max(d2) == 0) 0L else sum(d2 > max(d2) * rel_tol)
       })
}


# ---------------------------------------------------------------------------
# The congruence
# ---------------------------------------------------------------------------

#' Right-transform an architecture so that B' M B = G exactly
#'
#' `B = B0 A`, `A = Uc Dc^{-1/2} Q Gamma^{1/2} V'`, with `C = B0' M B0 =
#' Uc Dc Uc'`, `G = V Gamma V'` (rank-q parts) and `Q` the polar factor of a
#' thin SVD of `Uc' V` (Proposition 2). `Q` is what makes rank-deficient
#' architectures correct (Codex F2). For k = 1 this is the scalar rescale
#' `b0 * sqrt(G / sum(w b0^2))`: the polar factor cancels both eigenvector
#' signs.
#'
#' Two feasibility errors, checked in this order (gate A5):
#' 1. `rank(G) > rank(M)`: the anchor cannot carry the target for any effects
#'    (Theorem 1).
#' 2. `rank(B0' M B0) < rank(G)`: the anchor could, but this architecture
#'    cannot (e.g. proportional columns against a full-rank target).
#'
#' @param B0 m x k architecture.
#' @param target `.qtl_target_eigen()` result.
#' @param anchor `.qtl_anchor_diag()` / `.qtl_anchor_design()` result.
#' @return list(B, A, C, rank_C, rank_M).
#' @noRd
.qtl_congruence <- function(B0, target, anchor, rel_tol = 1e-10) {
  q <- target$rank
  k <- ncol(B0)
  if (q == 0L) {
    return(list(B = B0 * 0, A = matrix(0, k, k), C = matrix(0, k, k),
                rank_C = 0L, rank_M = NA_integer_))
  }
  rank_M <- anchor$rank(rel_tol)
  if (rank_M < q) {
    stop("The anchor cannot carry the target for any effects: rank(G) = ", q,
         " but the reference covariance M has rank ", rank_M, " at the ",
         "selected loci (too few independent segregating directions). Select ",
         "more loci, or loci that segregate in the base population.",
         call. = FALSE)
  }
  C <- anchor$cov(B0)
  C <- (C + t(C)) / 2
  eC <- .qtl_psd_eigen(C, "B0' M B0", rel_tol)
  rk <- eC$rank
  if (rk < q) {
    stop("The drawn architecture cannot carry the target: rank(B0' M B0) = ",
         rk, " but rank(G) = ", q, ", although the anchor has rank ", rank_M,
         ". The architecture's effect columns are (nearly) linearly ",
         "dependent on the segregating loci.", call. = FALSE)
  }
  Uc <- eC$vectors[, seq_len(rk), drop = FALSE]; dc <- eC$values[seq_len(rk)]
  Vg <- target$vectors[, seq_len(q), drop = FALSE]; gamma <- target$values[seq_len(q)]
  Qal <- svd(crossprod(Uc, Vg), nu = q, nv = q)
  A <- Uc %*% diag(1 / sqrt(dc), nrow = rk) %*% (Qal$u %*% t(Qal$v)) %*%
       diag(sqrt(gamma), nrow = q) %*% t(Vg)
  list(B = B0 %*% A, A = A, C = C, rank_C = rk, rank_M = rank_M)
}


# ---------------------------------------------------------------------------
# Internal mirror of the source's sim_qtl_effects() (tests and diagnostics)
# ---------------------------------------------------------------------------

#' @noRd
.qtl_call_architecture <- function(architecture, m, k, info) {
  if (is.null(architecture)) {
    B <- matrix(0, m, k)
    pool <- which(info$anchor_variable)
    if (length(pool)) B[pool, ] <- stats::rnorm(length(pool) * k)
    return(B)
  }
  if (is.matrix(architecture)) return(architecture)
  if (!is.function(architecture))
    stop("`architecture` must be NULL, an m x k numeric matrix, or a function.",
         call. = FALSE)
  fml <- formals(architecture)
  accepts_info <- !is.null(fml) && ("..." %in% names(fml) || "info" %in% names(fml))
  if (accepts_info) architecture(m, k, info) else architecture(m, k)
}

#' Isotropic basis for the dual anchor (Theorem 3)
#'
#' `q` is rank(G), not the number of trait columns. Internal only (Q6).
#' @noRd
.qtl_dual_basis <- function(Xc, weights, q, rel_tol) {
  n <- nrow(Xc)
  seg <- weights > .qtl_max_abs(weights) * rel_tol
  m_seg <- sum(seg)
  if (m_seg < q)
    stop("anchor = \"dual\" needs at least rank(G) = ", q, " segregating loci.",
         call. = FALSE)

  Y <- sweep(Xc[, seg, drop = FALSE], 2, sqrt(weights[seg]), "/")

  if (m_seg <= n) {
    R <- crossprod(Y) / (n - 1)
    eR <- eigen((R + t(R)) / 2, symmetric = TRUE)
    U <- eR$vectors
    lambda <- eR$values
    r <- m_seg
  } else {
    K <- tcrossprod(Y) / (n - 1)
    eK <- eigen((K + t(K)) / 2, symmetric = TRUE)
    r <- sum(eK$values > .qtl_max_abs(eK$values) * rel_tol)
    if (r == 0L) stop("anchor = \"dual\": founder panel has no variation.", call. = FALSE)
    lambda <- eK$values[seq_len(r)]
    U <- crossprod(Y, eK$vectors[, seq_len(r), drop = FALSE]) %*%
      diag(1 / sqrt((n - 1) * lambda), nrow = r)
  }

  n_pos <- sum(lambda > 1 + rel_tol)
  n_zero <- sum(abs(lambda - 1) <= rel_tol)
  n_neg <- sum(lambda < 1 - rel_tol) + (m_seg - r)
  max_dim <- n_zero + min(n_pos, n_neg)
  if (q > max_dim)
    stop("anchor = \"dual\" is infeasible for a target of rank ", q,
         ". Witt index of W = S0 - Delta is ", max_dim, ".", call. = FALSE)

  idx_zero <- which(abs(lambda - 1) <= rel_tol)
  idx_pos <- which(lambda > 1 + rel_tol)
  idx_neg <- which(lambda < 1 - rel_tol)
  n_use_zero <- min(q, length(idx_zero))
  Uiso <- if (n_use_zero) U[, idx_zero[seq_len(n_use_zero)], drop = FALSE] else matrix(0, m_seg, 0)
  need <- q - n_use_zero

  if (need > 0L) {
    n_free <- m_seg - r
    n_from_free <- min(need, n_free)
    Zfree <- matrix(0, m_seg, 0)
    if (n_from_free > 0L) {
      Z <- matrix(stats::rnorm(m_seg * n_from_free), m_seg, n_from_free)
      Z <- Z - U %*% crossprod(U, Z)
      Zfree <- qr.Q(qr(Z))[, seq_len(n_from_free), drop = FALSE]
    }
    for (a in seq_len(need)) {
      i <- idx_pos[a]
      pos_vec <- U[, i] / sqrt(lambda[i] - 1)
      neg_vec <- if (a <= n_from_free) Zfree[, a] else {
        j <- idx_neg[a - n_from_free]
        U[, j] / sqrt(1 - lambda[j])
      }
      Uiso <- cbind(Uiso, pos_vec + neg_vec)
    }
  }
  V <- matrix(0, ncol(Xc), q)
  V[seg, ] <- sweep(Uiso, 1, sqrt(weights[seg]), "/")
  list(V = V, inertia = c(pos = n_pos, neg = n_neg, zero = n_zero), max_dim = max_dim)
}

#' Internal mirror of the source's `sim_qtl_effects()`
#'
#' Same arguments and return list as the source, minus `marginal`,
#' `marginal_iters` and `warn_on_mismatch` (`warn_bounds = NULL` turns the
#' warning off, §7.4). Used by the ported test suite.
#' @noRd
.qtl_sim_effects <- function(X, G, architecture = NULL,
                             anchor = c("genic", "marginal_observed", "realised",
                                        "dual", "reference"),
                             reference = NULL, ploidy = 2, rel_tol = 1e-10,
                             warn_bounds = c(0.8, 1.25)) {
  .qtl_validate_numeric_matrix(X, "X", min_rows = 2L)
  .qtl_validate_numeric_matrix(G, "G")
  if (nrow(G) != ncol(G)) stop("`G` must be square.", call. = FALSE)
  .qtl_validate_scalar(ploidy, "ploidy", integral = TRUE, min = 1,
                       max = .Machine$integer.max)
  ploidy <- as.integer(round(ploidy))
  .qtl_validate_scalar(rel_tol, "rel_tol", min = 0, max = 1)
  .qtl_validate_warn_bounds(warn_bounds)
  dosage_tol <- sqrt(.Machine$double.eps) * max(1, ploidy)
  .qtl_validate_dosage(X, ploidy, "X", dosage_tol)

  anchor <- match.arg(anchor)
  if (anchor != "reference" && !is.null(reference))
    stop("`reference` is used only for `anchor = \"reference\"`.", call. = FALSE)
  if (anchor == "reference" && is.null(reference))
    stop("Supply `reference` when `anchor = \"reference\"`.", call. = FALSE)
  if (anchor == "dual" && !is.null(architecture))
    stop("`architecture` cannot be used with `anchor = \"dual\"`: the isotropic ",
         "basis is determined by W = S0 - Delta, not sampled.", call. = FALSE)

  eG <- .qtl_target_eigen(G, "G", rel_tol)
  G_requested <- eG$G
  q <- eG$rank
  k <- ncol(G_requested)
  G_target <- if (q == 0L) matrix(0, k, k) else {
    Vg <- eG$vectors[, seq_len(q), drop = FALSE]
    Vg %*% diag(eG$values[seq_len(q)], nrow = q) %*% t(Vg)
  }
  dimnames(G_target) <- dimnames(G_requested)

  n <- nrow(X); m <- ncol(X); centre <- colMeans(X); Xc <- sweep(X, 2, centre, "-")
  p <- pmin(pmax(centre / ploidy, 0), 1)
  v_hwe <- ploidy * p * (1 - p)
  v_observed <- colSums(Xc^2) / (n - 1)

  anc <- if (anchor == "genic" || anchor == "dual") {
    .qtl_anchor_diag(v_hwe)
  } else if (anchor == "marginal_observed") {
    .qtl_anchor_diag(v_observed)
  } else if (anchor == "realised") {
    .qtl_anchor_design(Xc, n - 1)
  } else if (is.numeric(reference) && is.null(dim(reference))) {
    if (length(reference) != m)
      stop("A numeric `reference` must give one weight per locus: expected ", m,
           " values, got ", length(reference), ".", call. = FALSE)
    if (anyNA(reference) || any(!is.finite(reference)))
      stop("`reference` weights must be finite and non-missing.", call. = FALSE)
    if (any(reference < 0))
      stop("`reference` weights must be non-negative (they are locus variances).",
           call. = FALSE)
    .qtl_anchor_diag(reference)
  } else {
    .qtl_validate_numeric_matrix(reference, "reference", min_rows = 2L, min_cols = 1L)
    if (ncol(reference) != m)
      stop("`reference` must have one column per locus (", m, "); got ",
           ncol(reference), ".", call. = FALSE)
    .qtl_validate_dosage(reference, ploidy, "reference", dosage_tol)
    .qtl_anchor_design(sweep(reference, 2, colMeans(reference), "-"),
                       nrow(reference) - 1)
  }

  info <- list(p = p, v_hwe = v_hwe, v_observed = v_observed,
               segregating = v_hwe > 0,
               anchor_variable = anc$locus_variance > 0)
  trait_names <- colnames(G_requested)
  if (is.null(trait_names)) trait_names <- rownames(G_requested)
  if (is.null(trait_names)) trait_names <- paste0("Trait", seq_len(k))
  dual_inertia <- NULL

  if (q == 0L) {
    if (!is.null(architecture))
      .qtl_validate_exact_dim(.qtl_call_architecture(architecture, m, k, info),
                              "architecture", m, k)
    B0 <- matrix(0, m, k); B <- B0; A <- matrix(0, k, k); rank_C <- 0L
  } else if (anchor == "dual") {
    db <- .qtl_dual_basis(Xc, v_hwe, q, rel_tol)
    dual_inertia <- db$inertia
    B0 <- db$V
    C <- crossprod(B0 * sqrt(v_hwe))
    eC <- .qtl_psd_eigen(C, "dual basis Gram matrix", rel_tol)
    rank_C <- eC$rank
    if (rank_C < q)
      stop("Dual basis is rank deficient under the genic anchor.", call. = FALSE)
    Uc <- eC$vectors[, seq_len(rank_C), drop = FALSE]; dc <- eC$values[seq_len(rank_C)]
    Vg <- eG$vectors[, seq_len(q), drop = FALSE]; gamma <- eG$values[seq_len(q)]
    C_inv_half <- Uc %*% diag(1 / sqrt(dc), nrow = rank_C) %*% t(Uc)
    A <- C_inv_half %*% diag(sqrt(gamma), nrow = q) %*% t(Vg)
    B <- B0 %*% A
  } else {
    B0 <- .qtl_call_architecture(architecture, m, k, info)
    .qtl_validate_numeric_matrix(B0, "architecture output")
    .qtl_validate_exact_dim(B0, "architecture output", m, k)
    cg <- .qtl_congruence(B0, eG, anc, rel_tol)
    B <- cg$B; A <- cg$A; rank_C <- cg$rank_C
  }

  dimnames(B) <- list(colnames(X), trait_names)
  dimnames(B0) <- list(colnames(X),
                       if (ncol(B0) == k) trait_names else paste0("iso", seq_len(ncol(B0))))

  G_anchor <- anc$cov(B)
  G_founder <- crossprod(Xc %*% B) / (n - 1)
  G_equilibrium <- crossprod(B * sqrt(v_hwe))
  founder_relative_eigen <- .qtl_relative_spectrum(eG, G_founder)
  equilibrium_relative_eigen <- .qtl_relative_spectrum(eG, G_equilibrium)

  if (q > 0L && !is.null(warn_bounds)) {
    if (min(founder_relative_eigen) < warn_bounds[1] ||
        max(founder_relative_eigen) > warn_bounds[2] ||
        min(equilibrium_relative_eigen) < warn_bounds[1] ||
        max(equilibrium_relative_eigen) > warn_bounds[2]) {
      warning("Realised variance departs from `G`. Founder spectrum: ",
              .qtl_spectrum_text(founder_relative_eigen),
              ". Equilibrium spectrum: ",
              .qtl_spectrum_text(equilibrium_relative_eigen), ".", call. = FALSE)
    }
  }

  row_norm2 <- rowSums(B^2)
  n_causal <- sum(row_norm2 > max(row_norm2, 0) * rel_tol^2)
  trace_ratio <- if (sum(diag(G_target)) > 0)
    sum(diag(G_equilibrium)) / sum(diag(G_target)) else NA_real_

  list(B = B, B0 = B0, transform = A, centre = centre, anchor = anchor,
       G_requested = G_requested, G_target = G_target, G_anchor = G_anchor,
       G_founder = G_founder, G_equilibrium = G_equilibrium,
       rank_G = q, rank_architecture = rank_C, n_causal = n_causal,
       trace_ratio_equilibrium = trace_ratio,
       founder_relative_eigen = founder_relative_eigen,
       equilibrium_relative_eigen = equilibrium_relative_eigen,
       dual_inertia = dual_inertia)
}
