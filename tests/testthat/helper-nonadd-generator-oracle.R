# Test oracle for the step-5b calibration gates
# (plans/import_qtl_effect_methods_phase_5_plan.md, 5b.9 "Oracle isolation").
#
# The source project's non-additive generator, copied verbatim:
#   /Users/austinputz/Claude/simulate_qtl_effects, commit 318e54f,
#   non-additive/R/qtl_effects_nonadd.R
# Entry points used by the tests: sim_qtl_effects_nonadd(), nonadd_covariates(),
# congruence(), solve_additive_stage(), solve_dd_mean(), zeng_appendix_A().
# Their transitive dependencies, checked with codetools::findGlobals() when the
# copy was made (2026-10-09): nonadd_decompose() (called unconditionally for
# the founder diagnostics) and the source's own .na_* helpers .na_max_abs(),
# .na_validate_numeric_matrix(), .na_validate_scalar(),
# .na_validate_exact_dim(), .na_psd_eigen(), .na_target(),
# .na_relative_spectrum(), .na_spectrum_text(), .na_trace_ratio().
#
# They are evaluated in their own environment whose parent is baseenv(): the
# package's internals (which also use the .na_ prefix) cannot shadow them, and
# they cannot silently reach a package function or .GlobalEnv. `stats::` calls
# are explicit in the source. Do not edit the functions: they are the
# independent reference. Production intentionally diverges (correlation-scale
# congruence, no eligibility zeroing, a robust degree-mean solver), so tests
# compare what must agree -- delivered covariances, feasibility, identities --
# not coefficients.

nonadd_oracle <- local({
  env <- new.env(parent = baseenv())
  eval(quote({
  .na_max_abs <- function(x) if (length(x)) max(abs(x)) else 0

  .na_validate_numeric_matrix <- function(x, name, min_rows = 1L, min_cols = 1L) {
    if (!is.matrix(x) || !is.numeric(x))
      stop("`", name, "` must be a numeric matrix.", call. = FALSE)
    if (nrow(x) < min_rows || ncol(x) < min_cols)
      stop("`", name, "` must have at least ", min_rows, " row(s) and ", min_cols,
           " column(s).", call. = FALSE)
    if (anyNA(x) || any(!is.finite(x)))
      stop("`", name, "` must contain only finite, non-missing values.", call. = FALSE)
    invisible(TRUE)
  }

  .na_validate_scalar <- function(x, name, integral = FALSE, min = -Inf, max = Inf) {
    if (!is.numeric(x) || length(x) != 1L || !is.finite(x))
      stop("`", name, "` must be a single finite numeric value.", call. = FALSE)
    if (integral && abs(x - round(x)) > .Machine$double.eps^0.5)
      stop("`", name, "` must be an integer.", call. = FALSE)
    if (x < min || x > max)
      stop("`", name, "` must lie in [", min, ", ", max, "].", call. = FALSE)
    invisible(TRUE)
  }

  .na_validate_exact_dim <- function(x, name, nr, nc) {
    if (nrow(x) != nr || ncol(x) != nc)
      stop("`", name, "` must be exactly ", nr, " x ", nc, "; got ", nrow(x), " x ", ncol(x),
           ".", call. = FALSE)
    invisible(TRUE)
  }

  .na_psd_eigen <- function(x, name, rel_tol) {
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

  ## A k x k target: square, symmetric within the relative tolerance, PSD. Returns the
  ## symmetrised matrix, its eigen-structure and rank.
  .na_target <- function(G, name, k, rel_tol) {
    G <- as.matrix(G)
    .na_validate_numeric_matrix(G, name)
    if (nrow(G) != ncol(G)) stop("`", name, "` must be square.", call. = FALSE)
    if (!is.null(k) && nrow(G) != k)
      stop("`", name, "` must be ", k, " x ", k, " to match `G_A`; got ", nrow(G), " x ",
           ncol(G), ".", call. = FALSE)
    if (.na_max_abs(G - t(G)) > rel_tol * max(.na_max_abs(G), .Machine$double.eps))
      stop("`", name, "` must be symmetric within the relative tolerance.", call. = FALSE)
    Gs <- (G + t(G)) / 2
    e <- .na_psd_eigen(Gs, name, rel_tol)
    list(G = Gs, eigen = e, rank = e$rank)
  }

  ## Relative-variance spectrum (Codex F4): eigenvalues of Gam^-1/2 V' G_cand V Gam^-1/2 over the
  ## support of the target. Trace cancellation cannot hide a per-contrast miss from this.
  .na_relative_spectrum <- function(target_eigen, candidate) {
    q <- target_eigen$rank
    if (q == 0L) return(numeric())
    V <- target_eigen$vectors[, seq_len(q), drop = FALSE]
    gamma <- target_eigen$values[seq_len(q)]
    W <- V %*% diag(1 / sqrt(gamma), nrow = q)
    R <- crossprod(W, candidate %*% W)
    sort(pmax(eigen((R + t(R)) / 2, symmetric = TRUE, only.values = TRUE)$values, 0))
  }

  .na_spectrum_text <- function(x) {
    if (!length(x)) return("not defined (zero target)")
    sprintf("[%.3g, %.3g]", min(x), max(x))
  }

  .na_trace_ratio <- function(candidate, target) {
    tt <- sum(diag(target)); if (tt > 0) sum(diag(candidate)) / tt else NA_real_
  }

  sim_qtl_effects_nonadd <- function(X, G_A, G_D = NULL, G_AA = NULL, pairs = NULL,
                                     anchor = c("genic", "realised"),
                                     B_a = NULL, B_d = NULL, B_aa = NULL,
                                     dd_mean = 0.19, dd_sd = 0.097,
                                     inbr_depr = NULL, rel_tol = 1e-8,
                                     warn_bounds = c(0.8, 1.25), warn_on_mismatch = TRUE) {
    anchor <- match.arg(anchor)

    ## ---- validation ------------------------------------------------------------------------
    .na_validate_numeric_matrix(X, "X", min_rows = 2L)
    storage.mode(X) <- "double"
    n <- nrow(X); m <- ncol(X)
    dosage_tol <- sqrt(.Machine$double.eps) * 2
    if (min(X) < -dosage_tol || max(X) > 2 + dosage_tol)
      stop("`X` must contain diploid allele dosages in [0, 2]; the non-additive construction ",
           "(heterozygote indicator, (2pq)^2 dominance weights) is defined for diploids only.",
           call. = FALSE)
    .na_validate_scalar(rel_tol, "rel_tol", min = 0, max = 1)
    if (!is.null(warn_bounds)) {
      if (!is.numeric(warn_bounds) || length(warn_bounds) != 2L || any(!is.finite(warn_bounds)))
        stop("`warn_bounds` must be NULL or two finite numeric values.", call. = FALSE)
      if (warn_bounds[1] <= 0 || warn_bounds[1] > warn_bounds[2])
        stop("`warn_bounds` must satisfy 0 < lower <= upper.", call. = FALSE)
    }
    if (!is.logical(warn_on_mismatch) || length(warn_on_mismatch) != 1L || is.na(warn_on_mismatch))
      stop("`warn_on_mismatch` must be a single TRUE or FALSE.", call. = FALSE)

    tA <- .na_target(G_A, "G_A", NULL, rel_tol); G_A <- tA$G; k <- nrow(G_A)
    trait_names <- colnames(G_A); if (is.null(trait_names)) trait_names <- rownames(G_A)
    if (is.null(trait_names)) trait_names <- paste0("Trait", seq_len(k))
    dimnames(G_A) <- list(trait_names, trait_names)
    has_d <- !is.null(G_D); has_aa <- !is.null(G_AA)
    tD <- tAA <- NULL
    if (has_d)  { tD  <- .na_target(G_D,  "G_D",  k, rel_tol); G_D  <- tD$G;  dimnames(G_D)  <- dimnames(G_A) }
    if (has_aa) { tAA <- .na_target(G_AA, "G_AA", k, rel_tol); G_AA <- tAA$G; dimnames(G_AA) <- dimnames(G_A) }

    if (!has_aa && !is.null(pairs)) stop("`pairs` is used only with `G_AA`.", call. = FALSE)
    if (has_aa) {
      if (is.null(pairs)) stop("Supply `pairs` (an npair x 2 matrix of locus indices) with `G_AA`.", call. = FALSE)
      pairs <- as.matrix(pairs)
      .na_validate_numeric_matrix(pairs, "pairs")
      if (ncol(pairs) != 2L) stop("`pairs` must have exactly 2 columns; got ", ncol(pairs), ".", call. = FALSE)
      if (any(abs(pairs - round(pairs)) > .Machine$double.eps^0.5))
        stop("`pairs` must hold integer locus indices.", call. = FALSE)
      if (min(pairs) < 1 || max(pairs) > m)
        stop("`pairs` must index loci in 1..", m, " (the columns of `X`).", call. = FALSE)
      if (any(pairs[, 1] == pairs[, 2]))
        stop("`pairs` must not pair a locus with itself (rows ",
             paste(which(pairs[, 1] == pairs[, 2]), collapse = ", "), ").", call. = FALSE)
      key <- paste(pmin(pairs[, 1], pairs[, 2]), pmax(pairs[, 1], pairs[, 2]))
      if (anyDuplicated(key))
        stop("`pairs` must not repeat a pair (in either order); duplicated rows: ",
             paste(which(duplicated(key)), collapse = ", "), ".", call. = FALSE)
      storage.mode(pairs) <- "integer"
    }
    npair <- if (has_aa) nrow(pairs) else 0L

    if (!has_d  && !is.null(B_d))  stop("`B_d` is used only with `G_D`.",   call. = FALSE)
    if (!has_aa && !is.null(B_aa)) stop("`B_aa` is used only with `G_AA`.", call. = FALSE)
    if (!is.null(B_a))  { B_a  <- as.matrix(B_a);  .na_validate_numeric_matrix(B_a,  "B_a");  .na_validate_exact_dim(B_a,  "B_a",  m, k) }
    if (!is.null(B_d))  { B_d  <- as.matrix(B_d);  .na_validate_numeric_matrix(B_d,  "B_d");  .na_validate_exact_dim(B_d,  "B_d",  m, k) }
    if (!is.null(B_aa)) { B_aa <- as.matrix(B_aa); .na_validate_numeric_matrix(B_aa, "B_aa"); .na_validate_exact_dim(B_aa, "B_aa", npair, k) }
    .na_validate_scalar(dd_mean, "dd_mean")
    .na_validate_scalar(dd_sd, "dd_sd", min = 0)
    if (!is.null(inbr_depr)) {
      if (!has_d) stop("`inbr_depr` needs `G_D`: inbreeding depression is a property of the dominance effects.", call. = FALSE)
      if (!is.numeric(inbr_depr) || length(inbr_depr) != k || any(!is.finite(inbr_depr)))
        stop("`inbr_depr` must be a finite numeric vector of length ", k, " (one target per trait).", call. = FALSE)
      if (!is.null(B_d))
        stop("`inbr_depr` re-derives `B_d` from `B_a` and the dominance degrees; do not also supply `B_d`.", call. = FALSE)
      if (dd_sd <= 0)
        stop("`inbr_depr` needs `dd_sd` > 0: with a constant dominance degree the ratio ID / sqrt(V_D) is fixed by `B_a` and cannot be targeted.", call. = FALSE)
    }

    ## ---- anchors and eligibility -------------------------------------------------------------
    cov0 <- nonadd_covariates(X, pairs = if (has_aa) pairs else NULL, anchor = anchor)
    M_A <- cov0$M_A; M_D <- cov0$M_D; M_AA <- cov0$M_AA; b <- cov0$b; cc <- cov0$c
    ## Eligibility is the anchor's own locus variance per component (the additive F1 rule), not a
    ## frequency test: under "realised" an all-heterozygous locus has zero dosage variance and
    ## zero dominance-contrast variance, so it gets nothing; under "genic" it is a segregating
    ## locus and keeps its effect. Fixed loci (2pq = 0) are ineligible under either.
    elig_A  <- diag(M_A) > rel_tol * max(diag(M_A), 1e-300)
    elig_D  <- if (has_d)  diag(M_D)  > rel_tol * max(diag(M_D),  1e-300) else NULL
    elig_AA <- if (has_aa) diag(M_AA) > rel_tol * max(diag(M_AA), 1e-300) else NULL

    ## ---- functional draws (RNG order unchanged from the original) ----------------------------
    supplied <- list(B_a = !is.null(B_a), B_d = !is.null(B_d), B_aa = !is.null(B_aa))
    if (is.null(B_a)) { B_a <- matrix(stats::rnorm(m * k), m, k); B_a[!elig_A, ] <- 0 }
    h <- NULL
    if (has_d) {
      if (is.null(B_d)) {
        h   <- matrix(stats::rnorm(m * k, dd_mean, dd_sd), m, k)   # dominance degrees d = h|a|
        B_d <- h * abs(B_a); B_d[!elig_D, ] <- 0
      }
    }
    if (has_aa) {
      if (is.null(B_aa)) { B_aa <- matrix(stats::rnorm(npair * k), npair, k); B_aa[!elig_AA, ] <- 0 }
    }

    ## ---- stage 1: A x A ----------------------------------------------------------------------
    if (has_aa) B_aa <- congruence(B_aa, M_AA, G_AA, rel_tol, what = "B_aa", target = "G_AA")

    ## ---- stage 2: dominance (optionally with an inbreeding-depression target) ----------------
    id_info <- NULL
    if (has_d) {
      w <- cov0$w_id                                     # 2 p q: weight of d_j in the inbreeding depression
      if (!is.null(inbr_depr)) {
        ## solve per trait for the dominance-degree mean that gives ID / sqrt(V_D) its target ratio;
        ## exact for k = 1 and for diagonal G_D (the k x k congruence below is then diagonal)
        id_info <- vector("list", k)
        for (t in seq_len(k)) {
          if (G_D[t, t] <= 0)
            stop("`inbr_depr` for trait ", t, " needs a positive dominance variance `G_D[", t, ", ", t, "]`.", call. = FALSE)
          u <- abs(B_a[, t]); z <- if (is.null(h)) stats::rnorm(m) else (h[, t] - dd_mean) / dd_sd
          v <- z * u; rho <- inbr_depr[t] / sqrt(G_D[t, t])
          r <- solve_dd_mean(u, v, w, M_D, rho, dd_sd)
          B_d[, t] <- (r$mu * u + dd_sd * v); id_info[[t]] <- r
        }
        B_d[!elig_D, ] <- 0
      }
      B_d <- congruence(B_d, M_D, G_D, rel_tol, what = "B_d", target = "G_D")
    }

    ## ---- stage 3: additive, through the coupling alpha = a + b d + sum e c -------------------
    C <- matrix(0, m, k)
    if (has_d)  C <- C + b * B_d
    if (has_aa) for (r in seq_len(npair)) {
      C[pairs[r, 1], ] <- C[pairs[r, 1], ] + cc[pairs[r, 2]] * B_aa[r, ]
      C[pairs[r, 2], ] <- C[pairs[r, 2], ] + cc[pairs[r, 1]] * B_aa[r, ]
    }
    sol <- solve_additive_stage(B_a, C, M_A, G_A, rel_tol, anchor = anchor)
    B_a <- sol$B_a; B_alpha <- B_a + C

    ## ---- names -------------------------------------------------------------------------------
    locus_names <- colnames(X)
    dimnames(B_a) <- list(locus_names, trait_names); dimnames(B_alpha) <- dimnames(B_a)
    if (has_d)  dimnames(B_d)  <- dimnames(B_a)
    if (has_aa) dimnames(B_aa) <- list(rownames(pairs), trait_names)

    ## ---- diagnostics: what was delivered, under the anchor, in the founders, at equilibrium ---
    delivered <- list(
      G_A  = crossprod(B_alpha, M_A %*% B_alpha),
      G_D  = if (has_d)  crossprod(B_d,  M_D  %*% B_d)  else NULL,
      G_AA = if (has_aa) crossprod(B_aa, M_AA %*% B_aa) else NULL)
    founders <- nonadd_decompose(X, B_a, if (has_d) B_d else NULL, if (has_aa) B_aa else NULL,
                                 if (has_aa) pairs else NULL)
    comp <- list(A = list(t = tA, founder = founders$real_A, equil = founders$genic_A))
    if (has_d)  comp$D  <- list(t = tD,  founder = founders$real_D,  equil = founders$genic_D)
    if (has_aa) comp$AA <- list(t = tAA, founder = founders$real_AA, equil = founders$genic_AA)
    diagnostics <- lapply(comp, function(z) list(
      rank_target = z$t$rank,
      founder_relative_eigen     = .na_relative_spectrum(z$t$eigen, z$founder),
      equilibrium_relative_eigen = .na_relative_spectrum(z$t$eigen, z$equil),
      trace_ratio_founder        = .na_trace_ratio(z$founder, z$t$G),
      trace_ratio_equilibrium    = .na_trace_ratio(z$equil,   z$t$G)))
    for (nm in names(diagnostics)) {
      dimnames(delivered[[paste0("G_", nm)]]) <- list(trait_names, trait_names)
    }
    if (isTRUE(warn_on_mismatch) && !is.null(warn_bounds)) {
      off <- vapply(names(diagnostics), function(nm) {
        d <- diagnostics[[nm]]
        sp <- c(d$founder_relative_eigen, d$equilibrium_relative_eigen)
        length(sp) > 0 && (min(sp) < warn_bounds[1] || max(sp) > warn_bounds[2])
      }, logical(1))
      if (any(off)) {
        msg <- vapply(names(diagnostics)[off], function(nm) {
          d <- diagnostics[[nm]]
          sprintf("G_%s: founders %s, equilibrium %s", nm,
                  .na_spectrum_text(d$founder_relative_eigen), .na_spectrum_text(d$equilibrium_relative_eigen))
        }, character(1))
        warning("Realised variance departs from the target under anchor \"", anchor, "\". ",
                "Relative-variance spectra -- ", paste(msg, collapse = "; "), ".", call. = FALSE)
      }
    }
    n_causal <- list(A = sum(rowSums(B_a^2) > 0))
    if (has_d)  n_causal$D  <- sum(rowSums(B_d^2)  > 0)
    if (has_aa) n_causal$AA <- sum(rowSums(B_aa^2) > 0)

    out <- list(B_a = B_a, B_d = if (has_d) B_d else NULL, B_aa = if (has_aa) B_aa else NULL,
                pairs = if (has_aa) pairs else NULL, B_alpha = B_alpha,
                anchor = anchor, T = sol$T, G_tilde = sol$G_tilde, min_G_A = sol$min_G_A,
                inbr_depr_info = id_info, b = b, c = cc,
                G_target = list(G_A = G_A, G_D = if (has_d) G_D else NULL, G_AA = if (has_aa) G_AA else NULL),
                delivered = delivered, founders = founders, diagnostics = diagnostics,
                eligible = list(A = elig_A, D = elig_D, AA = elig_AA), n_causal = n_causal,
                supplied = supplied,
                inbr_depr = if (has_d) as.vector(crossprod(cov0$w_id, B_d)) else NULL)
    class(out) <- "qtl_effects_nonadd"
    out
  }

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

  ## Congruence correction: returns B0 A with (B0 A)' M (B0 A) = G. Same algebra as the additive
  ## implementation (R/qtl_effects.R, manuscript Proposition 2): C = B0' M B0 = Uc Dc Uc',
  ## G = Vg Gam Vg', A = Uc Dc^-1/2 Q Gam^1/2 Vg' with Q the polar factor of Uc' Vg from a thin SVD.
  ## Q is what makes the reduced-rank case exact. The earlier symmetric-root form
  ## B0 C^-1/2 G^1/2 delivers G^1/2 P_range(C) G^1/2, which equals G only when range(G) lies
  ## inside range(C): for a rank-deficient architecture and a target pointing elsewhere it passed
  ## the rank check and returned the wrong covariance silently (found 2026-09-19; test 9).
  congruence <- function(B0, M, G, rel_tol = 1e-8, what = "B0", target = "G") {
    k <- ncol(B0)
    Cm <- crossprod(B0, M %*% B0)
    ec <- eigen((Cm + t(Cm)) / 2, symmetric = TRUE)
    eg <- eigen((G + t(G)) / 2, symmetric = TRUE)
    rk <- sum(ec$values > rel_tol * max(ec$values, 1e-300))
    q  <- sum(eg$values > rel_tol * max(eg$values, 1e-300))
    if (rk < q)
      stop("`", what, "` has rank ", rk, " under this anchor but `", target, "` has rank ", q,
           ": the target is not attainable from this architecture.", call. = FALSE)
    if (q == 0L) return(B0 * 0)
    Uc <- ec$vectors[, seq_len(rk), drop = FALSE]; dc <- ec$values[seq_len(rk)]
    Vg <- eg$vectors[, seq_len(q),  drop = FALSE]; gam <- eg$values[seq_len(q)]
    Qal <- svd(crossprod(Uc, Vg), nu = q, nv = q)
    A <- Uc %*% diag(1 / sqrt(dc), nrow = rk) %*% (Qal$u %*% t(Qal$v)) %*%
         diag(sqrt(gam), nrow = q) %*% t(Vg)
    B0 %*% A
  }

  ## Additive stage. Solve (B0 T + C)' M (B0 T + C) = G for T (k x k).
  ## Completing the square: (T + P^-1 Q)' P (T + P^-1 Q) = G - R + Q' P^-1 Q =: G_tilde,
  ## P = B0' M B0, Q = B0' M C, R = C' M C. Feasible iff G_tilde is PSD; the minimum additive
  ## covariance reachable by any T is R - Q' P^-1 Q (the part of the dominance/epistatic
  ## contribution to alpha that the additive architecture cannot cancel).
  solve_additive_stage <- function(B0, C, M, G, rel_tol = 1e-8, anchor = "this") {
    k <- ncol(B0)
    P <- crossprod(B0, M %*% B0); Q <- crossprod(B0, M %*% C); R <- crossprod(C, M %*% C)
    P <- (P + t(P)) / 2
    ep <- eigen(P, symmetric = TRUE)
    rk <- sum(ep$values > rel_tol * max(ep$values, 1e-300))
    if (rk < k)
      stop("`B_a` has rank ", rk, " under the \"", anchor, "\" anchor but there are ", k,
           " traits: the additive stage solves a matrix quadratic in the k x k right factor and ",
           "needs a full-rank additive architecture (a rank-deficient `G_A` is fine; `B_a` is not).",
           call. = FALSE)
    Pinv  <- ep$vectors %*% diag(1 / ep$values, k) %*% t(ep$vectors)
    Pinvh <- ep$vectors %*% diag(1 / sqrt(ep$values), k) %*% t(ep$vectors)
    min_G_A <- R - t(Q) %*% Pinv %*% Q
    Gt <- G - min_G_A; Gt <- (Gt + t(Gt)) / 2
    eg <- eigen(Gt, symmetric = TRUE)
    if (min(eg$values) < -rel_tol * max(abs(eg$values), 1e-300))
      stop(sprintf(paste0("`G_A` is infeasible for this architecture: G_A - (R - Q'P^-1 Q) is not PSD ",
                          "(smallest eigenvalue %.3g). The dominance/epistatic effects already imply at least ",
                          "this much additive (co)variance under the \"%s\" anchor (see `$min_G_A`); raise `G_A` ",
                          "or lower `G_D` / `G_AA`."), min(eg$values), anchor), call. = FALSE)
    Gth <- eg$vectors %*% diag(sqrt(pmax(eg$values, 0)), k) %*% t(eg$vectors)
    ## orthogonal freedom U: pick the one that maximises tr(T) (for k = 1: the larger root)
    Mx <- Gth %*% Pinvh
    sv <- svd(Mx); U <- sv$v %*% t(sv$u)
    Tm <- -Pinv %*% Q + Pinvh %*% U %*% Gth
    list(B_a = B0 %*% Tm, T = Tm, G_tilde = Gt, min_G_A = min_G_A)
  }

  ## Dominance-degree mean mu such that ID(mu) / sqrt(V_D(mu)) = rho, with d(mu) = mu u + sd v.
  ## ID = w'd (linear), V_D = d' M d (quadratic): rho^2 V_D = ID^2 is a quadratic in mu.
  solve_dd_mean <- function(u, v, w, M, rho, sd) {
    Mu <- M %*% u; Mv <- M %*% v
    uMu <- sum(u * Mu); uMv <- sum(u * Mv); vMv <- sum(v * Mv)
    wu <- sum(w * u); wv <- sum(w * v)
    A <- rho^2 * uMu - wu^2
    B <- 2 * sd * (rho^2 * uMv - wu * wv)
    Cc <- sd^2 * (rho^2 * vMv - wv^2)
    disc <- B^2 - 4 * A * Cc
    if (disc < 0) stop(sprintf("`inbr_depr` is infeasible: the ratio ID / sqrt(V_D) = %.3f exceeds the Cauchy-Schwarz bound for this architecture.", rho), call. = FALSE)
    roots <- (-B + c(-1, 1) * sqrt(disc)) / (2 * A)
    id <- roots * wu + sd * wv
    ok <- sign(id) == sign(rho) | rho == 0
    if (!any(ok)) stop("`inbr_depr`: no dominance-degree mean gives the requested sign of inbreeding depression.", call. = FALSE)
    mu <- max(roots[ok])
    list(mu = mu, roots = roots, disc = disc)
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

  ## Zeng et al. (2013) Appendix A: k = 1, HWE + LE ("by locus"): d <- s d, then a <- t a with t the
  ## larger root of  t^2 S1 + t S2 + S3 = V_A.  Returned alongside our stage-3 T for the same inputs.
  zeng_appendix_A <- function(p, a0, d0, V_A, V_D) {
    q <- 1 - p
    s <- sqrt(V_D / sum((2 * p * q * d0)^2)); d <- s * d0
    S1 <- sum(2 * p * q * a0^2); S2 <- sum(2 * p * q * 2 * (q - p) * a0 * d); S3 <- sum(2 * p * q * (q - p)^2 * d^2)
    disc <- S2^2 - 4 * S1 * (S3 - V_A)
    if (disc < 0) stop("Zeng Appendix A: no real root")
    t <- (-S2 + sqrt(disc)) / (2 * S1)
    list(a = t * a0, d = d, t = t, s = s, disc = disc)
  }
  }), env)
  env
})
