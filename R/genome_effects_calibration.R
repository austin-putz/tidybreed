# Calibration of additive, dominance and additive-by-additive effects to three
# targets at once: the engine of define_genome_effects() (plan §1.3, §9.4;
# phase-5 plan 5b.4).
#
# Ported from the source project `simulate_qtl_effects`,
# non-additive/R/qtl_effects_nonadd.R at commit 318e54f: the body of
# sim_qtl_effects_nonadd(), solve_additive_stage() and solve_dd_mean(). What
# changed, and why:
#
# * Pure. Every architecture is supplied (B_a, the standard-normal dominance
#   degree deviations z, B_aa); nothing is drawn and nothing is written.
# * One congruence. The A x A and dominance stages call Part A's
#   `.qtl_calibrate()` (correlation-scale rank and PSD judgement, verified
#   against the target as stored), not the source's `congruence()`; the two
#   have the same algebra.
# * Anchors are objects (`.qtl_anchor_diag()`, `.qtl_anchor_design()`), never
#   dense m x m or r x r matrices.
# * The additive stage works on the correlation scale of `G_A`, and its floor
#   is computed from the residual coupling `C_res`, which is PSD by
#   construction, with a rounding budget relative to the operands. The
#   source's floor test refused a target exactly on the floor.
# * The dominance-degree mean solver keeps the source's equation and root
#   policy and replaces its arithmetic, which failed on feasible inputs
#   (`.na_solve_dd_mean()`).
# * No eligibility zeroing (decision D3 (a)): a locus or pair that does not
#   segregate under the anchor keeps its effect; it contributes no variance
#   there, and other populations get whatever their frequencies give.

#' The three anchors, and the coupling's `b` and `c`
#'
#' * `"genic"`: from base allele frequencies `p` alone.
#'   `M_A = diag(2pq)`, `M_D = diag((2pq)^2)`, `M_AA = diag(2p_k q_k 2p_l q_l)`,
#'   `b = q - p`, `c = 2p - 1`.
#' * `"realised"`: from the base individuals' dosages `X` (n x m). `Z_A = X -
#'   2p`, `Z_D` the heterozygote indicator made orthogonal to dosage within
#'   each locus by the observed regression `b`, `Z_AA` the centred products of
#'   each pair's `Z_A` columns; every `M_c = Z_c' Z_c / (n - 1)`.
#'
#' The A x A anchor needs the pairs, which may be drawn later: add it with
#' `.na_aa_anchor()`.
#'
#' @return list(anchor, p, w = 2pq, b, cc, A, D, AA = NULL, Z_A = NULL).
#' @keywords internal
#' @noRd
.na_anchors <- function(anchor, p = NULL, X = NULL) {
  if (anchor == "genic") {
    w <- 2 * p * (1 - p)
    return(list(anchor = anchor, p = p, w = w, b = (1 - p) - p, cc = 2 * p - 1,
                A = .qtl_anchor_diag(w), D = .qtl_anchor_diag(w^2), AA = NULL,
                Z_A = NULL))
  }
  n   <- nrow(X)
  p   <- colMeans(X) / 2
  Z_A <- sweep(X, 2L, 2 * p, "-")
  W   <- (X == 1) * 1
  Wc  <- sweep(W, 2L, colMeans(W), "-")
  vx  <- colSums(Z_A^2) / n
  cwx <- colSums(Wc * Z_A) / n
  b   <- ifelse(vx > 0, cwx / vx, 0)
  Z_D <- Wc - sweep(Z_A, 2L, b, "*")
  list(anchor = anchor, p = p, w = 2 * p * (1 - p), b = b, cc = 2 * p - 1,
       A = .qtl_anchor_design(Z_A, n - 1L), D = .qtl_anchor_design(Z_D, n - 1L),
       AA = NULL, Z_A = Z_A)
}

#' Add the A x A anchor for a set of pairs
#'
#' @param pairs Two-column integer matrix of positions in the anchors' loci.
#' @keywords internal
#' @noRd
.na_aa_anchor <- function(anchors, pairs) {
  if (anchors$anchor == "genic") {
    anchors$AA <- .qtl_anchor_diag(anchors$w[pairs[, 1]] * anchors$w[pairs[, 2]])
  } else {
    Z <- anchors$Z_A[, pairs[, 1], drop = FALSE] *
      anchors$Z_A[, pairs[, 2], drop = FALSE]
    anchors$AA <- .qtl_anchor_design(sweep(Z, 2L, colMeans(Z), "-"),
                                     nrow(Z) - 1L)
  }
  anchors
}

#' The coupling `C = diag(b) B_d + E_c`
#'
#' `E_c` puts `c_l e_r` on locus `k` and `c_k e_r` on locus `l` for each pair
#' `r = (k, l)`: the additive effect a pair induces at a locus whose partner's
#' mean dosage is not 1. Accumulated with one `rowsum()`, never a loop over
#' pairs; a locus in several pairs (a hub) collects all of them.
#'
#' @return m x k matrix.
#' @keywords internal
#' @noRd
.na_coupling <- function(m, k, b, cc, B_d = NULL, B_aa = NULL, pairs = NULL) {
  C <- matrix(0, m, k)
  if (!is.null(B_d)) C <- b * B_d
  if (!is.null(B_aa) && nrow(pairs) > 0L) {
    s <- rowsum(rbind(cc[pairs[, 2]] * B_aa, cc[pairs[, 1]] * B_aa),
                c(pairs[, 1], pairs[, 2]))
    at <- as.integer(rownames(s))
    C[at, ] <- C[at, ] + s
  }
  C
}

#' Calibrate one block, naming it and the anchor in any error
#' @keywords internal
#' @noRd
.na_block_calibrate <- function(B0, std, anchor, what, anchor_name = NULL) {
  tryCatch(.qtl_calibrate(B0, std, anchor)$B, error = function(e) {
    stop("The '", what, "' block",
         if (!is.null(anchor_name)) paste0(" (\"", anchor_name, "\" anchor)"),
         ": ", conditionMessage(e), call. = FALSE)
  })
}

#' Calibrate A (+ D, + A x A) to their targets, highest order first
#'
#' Stage 1 calibrates `B_aa` to `G_AA`. Stage 2 builds the dominance
#' architecture `B_d = (mean + sd z) |B_a|` (each trait named in
#' `inbreeding_depression` gets its own solved mean, the others keep the drawn
#' degrees, decision D4) and calibrates it to `G_D`. Stage 3 fixes the
#' coupling `C` and solves `(B_a T + C)' M_A (B_a T + C) = G_A`
#' (`.na_additive_stage()`). With no non-zero dominance or A x A target it is
#' Part A's congruence alone (`.qtl_calibrate()`), never the quadratic with
#' `C = 0` (gate C4).
#'
#' Every present block's delivered covariance is verified against its target
#' as stored, at `QTL_CALIBRATION_TOL` on the correlation scale; a miss is an
#' error.
#'
#' @param anchors `.na_anchors()` result (with `.na_aa_anchor()` when `G_AA`).
#' @param G_A,G_D,G_AA Named k x k targets; `NULL` = absent block.
#' @param B_a m x k additive architecture.
#' @param z m x k standard-normal dominance degree deviations (with `G_D`).
#' @param B_aa r x k A x A architecture (with `G_AA`).
#' @param pairs r x 2 positions into the loci (with `G_AA`).
#' @param inbreeding_depression `NULL` or a numeric vector named by trait.
#' @return list(route, B_a = functional additive effects a, B_d, B_aa,
#'   B_alpha = statistical additive effects under the anchor, C, floor,
#'   floor_s, G_tilde_min, delivered (list A, D, AA), inbreeding (list
#'   requested, delivered, mean)).
#' @keywords internal
#' @noRd
.na_calibrate <- function(anchors, G_A, G_D = NULL, G_AA = NULL, B_a, z = NULL,
                          B_aa = NULL, pairs = NULL,
                          dominance_degree_mean = 0.19,
                          dominance_degree_sd = 0.097,
                          inbreeding_depression = NULL, rel_tol = 1e-10) {
  m <- nrow(B_a); k <- ncol(B_a)
  traits <- rownames(G_A)
  if (is.null(traits)) traits <- paste0("Trait", seq_len(k))
  std_A  <- .qtl_target_std(unname(G_A), "G_A")
  std_D  <- if (!is.null(G_D))  .qtl_target_std(unname(G_D), "G_D")
  std_AA <- if (!is.null(G_AA)) .qtl_target_std(unname(G_AA), "G_AA")
  live <- function(s) !is.null(s) && s$rank > 0L
  delivered <- list(A = NULL, D = NULL, AA = NULL)
  ib <- list(requested = inbreeding_depression, delivered = NULL,
             mean = stats::setNames(rep(NA_real_, k), traits))

  if (!live(std_D) && !live(std_AA)) {
    B <- .na_block_calibrate(B_a, std_A, anchors$A, "additive", anchors$anchor)
    delivered$A <- anchors$A$cov(B)
    if (!is.null(G_D)) {
      B_d <- matrix(0, m, k)
      delivered$D <- matrix(0, k, k)
    }
    if (!is.null(G_AA)) {
      B_aa <- matrix(0, nrow(pairs), k)
      delivered$AA <- matrix(0, k, k)
    }
    return(list(route = "additive", B_a = B,
                B_d = if (!is.null(G_D)) B_d, B_aa = if (!is.null(G_AA)) B_aa,
                B_alpha = B, C = matrix(0, m, k), floor = NULL, floor_s = NULL,
                G_tilde_min = NA_real_, delivered = .na_names(delivered, traits),
                inbreeding = ib))
  }

  # Stage 1: A x A.
  if (!is.null(G_AA)) {
    B_aa <- .na_block_calibrate(B_aa, std_AA, anchors$AA, "additive_by_additive",
                                anchors$anchor)
    delivered$AA <- anchors$AA$cov(B_aa)
  }

  # Stage 2: dominance, from the degrees of the architecture |B_a|.
  B_d <- NULL
  if (!is.null(G_D)) {
    u <- abs(B_a)
    B_d <- (dominance_degree_mean + dominance_degree_sd * z) * u
    for (t in names(inbreeding_depression)) {
      j   <- match(t, traits)
      rho <- inbreeding_depression[[t]] / sqrt(std_D$d[j])
      v   <- z[, j] * u[, j]
      r <- tryCatch(
        .na_solve_dd_mean(u[, j], v, anchors$w, anchors$D, rho,
                          dominance_degree_sd, dominance_degree_mean, rel_tol),
        error = function(e) stop("`inbreeding_depression` for trait '", t,
                                 "': ", conditionMessage(e), call. = FALSE))
      B_d[, j] <- r$mu * u[, j] + dominance_degree_sd * v
      ib$mean[t] <- r$mu
    }
    B_d <- .na_block_calibrate(B_d, std_D, anchors$D, "dominance", anchors$anchor)
    delivered$D <- anchors$D$cov(B_d)
    ib$delivered <- stats::setNames(as.vector(crossprod(anchors$w, B_d)), traits)
  }

  # Stage 3: additive, through the coupling alpha = a + b d + sum e c.
  C  <- .na_coupling(m, k, anchors$b, anchors$cc, B_d,
                     if (!is.null(G_AA)) B_aa, pairs)
  st <- .na_additive_stage(B_a, C, std_A, anchors$A, anchors$anchor, traits,
                           rel_tol)
  delivered$A <- anchors$A$cov(st$B_alpha)
  err <- .qtl_target_error((delivered$A + t(delivered$A)) / 2, std_A)
  if (err > QTL_CALIBRATION_TOL) {
    stop("The calibration did not reach the 'additive' target: the delivered ",
         "covariance differs from it by ", format(err, digits = 3), " on the ",
         "correlation scale (tolerance ", QTL_CALIBRATION_TOL, "). The drawn ",
         "architecture is numerically ill-conditioned against this anchor; ",
         "select more segregating loci, or try another seed. Nothing was ",
         "written.", call. = FALSE)
  }
  list(route = "nonadditive", B_a = st$B_alpha - C, B_d = B_d,
       B_aa = if (!is.null(G_AA)) B_aa, B_alpha = st$B_alpha, C = C,
       floor = st$floor, floor_s = st$floor_s, G_tilde_min = st$G_tilde_min,
       delivered = .na_names(delivered, traits), inbreeding = ib)
}

#' @keywords internal
#' @noRd
.na_names <- function(delivered, traits) {
  lapply(delivered, function(x) {
    if (is.null(x)) return(NULL)
    x <- (x + t(x)) / 2
    dimnames(x) <- list(traits, traits)
    x
  })
}

#' The additive stage: solve `(B_a T + C)' M (B_a T + C) = G_A` for `T`
#'
#' Completing the square on the correlation scale of `G_A`
#' (`S = diag(sqrt(diag(G_A)))`, positive definite on this route by decision
#' D5), with every quantity in standardised coordinates:
#'
#' ```
#' R_A     = S^-1 G_A S^-1
#' B_s     = B_a S^-1 ;  C_s = C S^-1
#' P_s     = B_s' M B_s ;  Q_s = B_s' M C_s
#' C_res_s = C_s - B_s P_s^-1 Q_s        # the coupling the span cannot cancel
#' floor_s = C_res_s' M C_res_s          # PSD by construction
#' G~_s    = R_A - floor_s
#' ```
#'
#' Feasible iff `G~_s` is PSD. A negative eigenvalue within
#' `rel_tol * max(1, ||floor_s||)` is rounding (a target exactly on the floor)
#' and is clamped to 0. The solution is `B_s Y + C_res_s` with
#' `Y' P_s Y = G~_s`; among the orthogonal freedom, the source's choice that
#' maximises `tr(T)` (for k = 1, the larger root). Rescaled to original units
#' once, at the end; `floor = S floor_s S` is for messages only.
#'
#' The floor is conditional on the sampled architecture: the least additive
#' covariance reachable within the drawn additive span, given the calibrated
#' dominance and A x A effects.
#'
#' @param std_A `.qtl_target_std()` of `G_A`.
#' @param M The additive anchor object.
#' @param anchor `"genic"` or `"realised"`, for messages.
#' @return list(B_alpha, floor, floor_s, G_tilde_min).
#' @keywords internal
#' @noRd
.na_additive_stage <- function(B_a, C, std_A, M, anchor, traits,
                               rel_tol = 1e-10) {
  k <- ncol(B_a)
  s <- sqrt(std_A$d)
  R_A <- std_A$G / outer(s, s)
  B_s <- sweep(B_a, 2L, s, "/")
  C_s <- sweep(C, 2L, s, "/")

  P_s <- M$cov(B_s)
  P_s <- (P_s + t(P_s)) / 2
  ep  <- eigen(P_s, symmetric = TRUE)
  top <- max(ep$values, 0)
  rk  <- if (top == 0) 0L else sum(ep$values > rel_tol * top)
  if (rk < k) {
    stop("The drawn additive architecture has rank ", rk, " under the \"",
         anchor, "\" anchor but there are ", k, " traits: the additive stage ",
         "solves for a k x k factor of the architecture and needs it at full ",
         "rank. Select more loci that segregate in the base.", call. = FALSE)
  }
  V <- ep$vectors; lam <- ep$values
  Pinv_Q  <- V %*% (crossprod(V, M$cross(B_s, C_s)) / lam)
  C_res_s <- C_s - B_s %*% Pinv_Q
  floor_s <- M$cov(C_res_s)
  floor_s <- (floor_s + t(floor_s)) / 2

  Gt <- R_A - floor_s
  Gt <- (Gt + t(Gt)) / 2
  eg <- eigen(Gt, symmetric = TRUE)
  budget <- rel_tol * max(1, max(abs(eigen(floor_s, symmetric = TRUE,
                                           only.values = TRUE)$values)))
  floor <- floor_s * outer(s, s)
  dimnames(floor) <- list(traits, traits)
  if (min(eg$values) < -budget) {
    stop("`G_A` is below the additive floor for this sampled architecture ",
         "under the \"", anchor, "\" anchor. The dominance and ",
         "additive-by-additive effects already induce additive variance that ",
         "the drawn additive architecture cannot cancel: at least ",
         paste0(traits, " = ", signif(diag(floor), 6), collapse = ", "),
         " (the floor's diagonal; G_A - floor has smallest eigenvalue ",
         signif(min(eg$values), 4), " on the correlation scale). Raise `G_A` ",
         "or lower `G_D` / `G_AA`. The floor belongs to this draw: another ",
         "seed gives another floor. Nothing was written.", call. = FALSE)
  }
  Gth   <- eg$vectors %*% (sqrt(pmax(eg$values, 0)) * t(eg$vectors))
  Pinvh <- V %*% (t(V) / sqrt(lam))
  sv    <- svd(Gth %*% Pinvh)
  Y     <- Pinvh %*% (sv$v %*% t(sv$u)) %*% Gth
  B_alpha <- sweep(B_s %*% Y + C_res_s, 2L, s, "*")
  list(B_alpha = B_alpha, floor = floor, floor_s = floor_s,
       G_tilde_min = min(eg$values))
}

#' The dominance-degree mean for a requested inbreeding depression
#'
#' The source's `solve_dd_mean()` equation with robust arithmetic (phase-5
#' plan 5b.4). Degrees are `d = sd (x u + v)`, `x = mu / sd`, so
#' `ID = sd (x w'u + w'v)` and `V_D = sd^2 (x^2 u'Mu + 2x u'Mv + v'Mv)`. The
#' target ratio `rho = ID / sqrt(V_D)` gives
#' `rho^2 V_D - ID^2 = sd^2 (A x^2 + B x + C)`, which does not depend on `sd`.
#'
#' Branches: all three coefficients negligible (the ratio holds wherever the
#' sign fits and `V_D > 0`: candidates in a fixed order); `A` negligible (the
#' linear root); `A`, `B` negligible and `C` not (a limit no finite mean
#' attains); otherwise the stable quadratic with `sign+(0) = 1`, and the
#' vertex when the discriminant is within rounding of 0. Every
#' candidate must give finite `x`, `ID` with the sign of `rho`, and a
#' **positive** `V_D`; the largest accepted quadratic root wins, as in the
#' source. The result is verified.
#'
#' @param u `|B_a[, t]|`; `v` `z_t * u`; `w` `2pq`; `M` the dominance anchor.
#' @return list(mu, x, branch, ratio).
#' @keywords internal
#' @noRd
.na_solve_dd_mean <- function(u, v, w, M, rho, sd, mean, tol = 1e-10) {
  Q   <- M$cov(cbind(u, v))
  uMu <- Q[1, 1]; uMv <- (Q[1, 2] + Q[2, 1]) / 2; vMv <- Q[2, 2]
  wu  <- sum(w * u); wv <- sum(w * v)
  r2  <- rho^2
  A <- r2 * uMu - wu^2
  B <- 2 * (r2 * uMv - wu * wv)
  C <- r2 * vMv - wv^2
  s <- r2 * (uMu + 2 * abs(uMv) + vMv) + (abs(wu) + abs(wv))^2
  small <- function(x) abs(x) <= tol * s
  id <- function(x) x * wu + wv
  vd <- function(x) x^2 * uMu + 2 * x * uMv + vMv
  failed <- character(0)
  accept <- function(x) {
    if (!is.finite(x)) return(FALSE)
    i  <- id(x)
    ib <- tol * (abs(x * wu) + abs(wv))
    ok_sign <- if (rho == 0) abs(i) <= ib else
      (abs(i) > ib && sign(i) == sign(rho))
    v  <- vd(x)
    ok_var <- is.finite(v) && v > tol * (x^2 * uMu + 2 * abs(x * uMv) + vMv)
    if (!ok_sign) failed <<- c(failed, "sign")
    if (!ok_var)  failed <<- c(failed, "variance")
    ok_sign && ok_var
  }
  bound <- function() {
    g <- c(wu, wv)
    Qi <- tryCatch(solve(Q), error = function(e) NULL)
    if (is.null(Qi)) return("")
    paste0(" (attainable |ratio| at most ",
           signif(sqrt(max(0, sum(g * (Qi %*% g)))), 6), ")")
  }

  if (small(A) && small(B) && small(C)) {
    branch <- "degenerate"
    xp <- mean / sd
    cand <- if (abs(wu) <= tol * (abs(wu) + abs(wv))) {
      c(xp, xp + 1, xp - 1)
    } else {
      x0 <- -wv / wu
      c(xp, 2 * x0 - xp, x0 + sign(rho) * sign(wu))
    }
    x <- NA_real_
    for (cx in cand) if (accept(cx)) { x <- cx; break }
  } else if (small(A) && !small(B)) {
    branch <- "linear"
    cand <- -C / B
    x <- if (accept(cand)) cand else NA_real_
  } else if (small(A)) {
    stop("the requested ratio ID / sqrt(V_D) = ", signif(rho, 6), " is a ",
         "limit no finite dominance-degree mean attains for this ",
         "architecture", bound(), ".", call. = FALSE)
  } else {
    branch <- "quadratic"
    disc <- B^2 - 4 * A * C
    if (disc < -tol * (B^2 + 4 * abs(A) * abs(C))) {
      stop("the ratio ID / sqrt(V_D) = ", signif(rho, 6), " exceeds the ",
           "Cauchy-Schwarz bound for this architecture", bound(), ".",
           call. = FALSE)
    }
    # A discriminant within rounding of 0 is a repeated root, taken at the
    # vertex: sqrt() of the residue would move both roots by ~sqrt(eps), and
    # C / q of two residues is noise. This includes rho = 0, where the
    # quadratic is -(ID)^2 and its one root is the root of ID.
    cand <- if (abs(disc) <= tol * (B^2 + 4 * abs(A) * abs(C))) {
      -B / (2 * A)
    } else {
      qq <- -(B + (if (B >= 0) 1 else -1) * sqrt(disc)) / 2
      c(qq / A, C / qq)
    }
    ok <- vapply(cand, accept, logical(1))
    x <- if (any(ok)) max(cand[ok]) else NA_real_
  }
  if (is.na(x)) {
    stop("no dominance-degree mean gives ",
         if ("sign" %in% failed && !"variance" %in% failed)
           "the requested sign of inbreeding depression"
         else if ("variance" %in% failed && !"sign" %in% failed)
           "a positive dominance variance at the requested ratio"
         else "the requested sign of inbreeding depression with a positive dominance variance",
         " for this architecture.", call. = FALSE)
  }
  ratio <- id(x) / sqrt(vd(x))
  if (abs(ratio - rho) > QTL_CALIBRATION_TOL * max(1, abs(rho))) {
    stop("the solved dominance-degree mean delivers the ratio ",
         signif(ratio, 8), " instead of ", signif(rho, 8), ".", call. = FALSE)
  }
  list(mu = sd * x, x = x, branch = branch, ratio = ratio)
}
