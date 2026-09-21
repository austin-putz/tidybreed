#' Correlated draws: block loading and the conditional resolver
#'
#' @description
#' The two database-independent halves of sequential correlated sampling
#' (`plans/sample_correlated_effects.md` §5.2–5.4):
#'
#' * [find_covariance_blocks()] reads `phenotype_var_comp` once and returns the
#'   complete covariance block(s) touching a set of phenotypes, one validated
#'   matrix per stratum. It is the only place a stored block becomes a matrix.
#' * [resolve_correlated_draws()] takes one such matrix and, for a set of
#'   entities with possibly different observed coordinates, returns draws of the
#'   requested coordinates from the Gaussian conditional on what each entity has
#'   already realized. It has no notion of individuals, strata or tables — the
#'   adapters in `add_phenotype()` own those (the residual adapter in
#'   `R/add_phenotype_stages.R`; the named-effect adapter lands in Phase 6).
#'
#' Neither function writes. The resolver is the only one that consumes RNG, and
#' it does so only after every check has passed and only when there is
#' something to draw.
#'
#' @name correlated_draws
#' @keywords internal
NULL


# ── Block loading ───────────────────────────────────────────────────────────

#' Load the covariance block(s) touching a set of phenotypes
#'
#' A *block* is a connected component of the graph whose edges are stored pair
#' rows in `phenotype_var_comp` for `effect_name` (a `cov_value` of `0` is an
#' edge). The requested phenotypes may fall into several components; each is
#' returned with every stratum stored for it. A phenotype with no row at all
#' belongs to no block and is simply absent from the result — callers decide
#' what that means (for the residual effect it is "no residual variance").
#'
#' The writers in `R/phenotype_cov_block.R` guarantee that every stratum of a
#' block is a complete, symmetric matrix over the same phenotypes and that a
#' block has one `(condition_table, condition_column)`. Rows can still be
#' removed by hand with `remove_rows()`, so the loader re-checks those
#' invariants and errors, naming the block, when they no longer hold.
#'
#' @param conn A DBI connection.
#' @param effect_name `"residual"` or a named random effect.
#' @param phenotype_names Character vector of the phenotypes whose blocks are
#'   wanted.
#' @return A list of blocks, ordered by their first member; possibly empty.
#'   Each block is a list with:
#'   \describe{
#'     \item{`effect_name`}{as supplied.}
#'     \item{`phenotypes`}{sorted character vector of members.}
#'     \item{`condition_table`, `condition_column`}{the block's condition
#'       column, or `NULL` when it has only the unconditional stratum.}
#'     \item{`unconditional`}{the unconditional matrix, or `NULL` if none is
#'       stored.}
#'     \item{`conditional`}{a list of matrices named by `condition_level`
#'       (empty when there are no conditional strata).}
#'   }
#'   Every matrix is `phenotypes` x `phenotypes` with dimnames, exactly
#'   symmetric, in the sorted member order.
#' @keywords internal
find_covariance_blocks <- function(conn, effect_name, phenotype_names) {
  stopifnot(is.character(effect_name), length(effect_name) == 1L)
  targets <- unique(as.character(phenotype_names))
  if (length(targets) == 0L) return(list())

  members <- .pvc_block_members(conn, effect_name, targets)
  rows <- .pvc_block_rows(conn, effect_name, members)
  if (nrow(rows) == 0L) return(list())

  # Connected components over pair-row existence (all strata together).
  comp <- stats::setNames(seq_along(members), members)
  find_root <- function(i) {
    while (comp[[i]] != i) i <- comp[[i]]
    i
  }
  for (k in seq_len(nrow(rows))) {
    a <- find_root(match(rows$phenotype_name_1[k], members))
    b <- find_root(match(rows$phenotype_name_2[k], members))
    if (a != b) comp[[max(a, b)]] <- min(a, b)
  }
  roots <- vapply(seq_along(members), find_root, integer(1))
  present <- unique(c(rows$phenotype_name_1, rows$phenotype_name_2))
  groups <- split(members, roots)
  groups <- Filter(function(g) any(g %in% present), groups)
  groups <- groups[order(vapply(groups, `[[`, character(1), 1L))]

  unname(lapply(groups, function(g) .cd_build_block(effect_name, sort(g), rows)))
}

# Assemble one block from its rows, checking the D1 invariants.
.cd_build_block <- function(effect_name, members, rows) {
  n <- length(members)
  rows <- rows[rows$phenotype_name_1 %in% members &
               rows$phenotype_name_2 %in% members, , drop = FALSE]
  label <- paste0("The '", effect_name, "' covariance block ", .pvc_set(members))
  redeclare <- paste0(
    "Redeclare the block:\n\n  ", .pvc_writer_call(effect_name, members), "\n")

  is_uncond <- is.na(rows$condition_column)
  cond <- rows[!is_uncond, , drop = FALSE]
  cond_key <- character(0)
  if (nrow(cond) > 0L) {
    if (anyNA(cond$condition_table) || anyNA(cond$condition_level)) {
      stop(label, " has a conditional row with a NULL condition_table or ",
           "condition_level. ", redeclare, call. = FALSE)
    }
    cond_key <- sort(unique(paste(cond$condition_table, cond$condition_column,
                                  sep = ".")))
    if (length(cond_key) != 1L) {
      stop(label, " is conditioned on more than one column (",
           paste(cond_key, collapse = ", "), "). A block has one condition ",
           "column. ", redeclare, call. = FALSE)
    }
  }
  build <- function(sub, what) {
    if (nrow(sub) != n * n) {
      pairs <- paste(sub$phenotype_name_1, sub$phenotype_name_2)
      want  <- paste(rep(members, each = n), rep(members, times = n))
      if (anyDuplicated(pairs)) {
        stop(label, " has duplicate rows in the ", what, " stratum. ",
             redeclare, call. = FALSE)
      }
      missing <- setdiff(want, pairs)
      stop(label, " is incomplete in the ", what, " stratum: no row for ",
           paste0("(", sub("^(\\S+) (\\S+)$", "\\1, \\2", missing), ")",
                  collapse = ", "),
           ". ", redeclare, call. = FALSE)
    }
    M <- matrix(NA_real_, n, n, dimnames = list(members, members))
    M[cbind(match(sub$phenotype_name_1, members),
            match(sub$phenotype_name_2, members))] <- sub$cov_value
    if (anyNA(M)) {
      stop(label, " has duplicate rows in the ", what, " stratum. ",
           redeclare, call. = FALSE)
    }
    if (!all(is.finite(M))) {
      stop(label, " has a non-finite cov_value in the ", what, " stratum. ",
           redeclare, call. = FALSE)
    }
    if (!identical(M, t(M))) {
      stop(label, " is not symmetric in the ", what, " stratum. ",
           redeclare, call. = FALSE)
    }
    M
  }

  unconditional <- NULL
  if (any(is_uncond)) {
    unconditional <- build(rows[is_uncond, , drop = FALSE], "unconditional")
  }
  conditional <- list()
  if (nrow(cond) > 0L) {
    levels <- sort(unique(cond$condition_level))
    conditional <- lapply(levels, function(lv) {
      build(cond[cond$condition_level == lv, , drop = FALSE],
            paste0(cond_key, " = '", lv, "'"))
    })
    names(conditional) <- levels
  }

  list(
    effect_name      = effect_name,
    phenotypes       = members,
    condition_table  = if (nrow(cond) > 0L) cond$condition_table[[1L]]  else NULL,
    condition_column = if (nrow(cond) > 0L) cond$condition_column[[1L]] else NULL,
    unconditional    = unconditional,
    conditional      = conditional
  )
}


# ── The resolver ────────────────────────────────────────────────────────────

#' Draw coordinates of a Gaussian block conditional on observed coordinates
#'
#' @description
#' For `n` entities sharing one covariance matrix `covariance` over a block of
#' coordinates, draws the `sample_coordinates` of every entity from the
#' multivariate normal conditional on that entity's observed coordinates. An
#' entity's observed set is whatever is non-`NA` in its row of `observed`;
#' entities with different observed sets are grouped by pattern and the
#' conditional mean coefficients and covariance are computed once per pattern.
#' Coordinates that are neither observed nor sampled are latent and are never
#' realized.
#'
#' The function is pure apart from RNG: it never reads or writes a database and
#' knows nothing about individuals, strata or tables.
#'
#' @section Numerical contract:
#' * `tolerance` is *relative*. The absolute tolerance is
#'   `tolerance * lambda_max(covariance)`, so the same relative perturbation is
#'   treated the same way at variance `1e-6` and `1e6`. The default is
#'   `nrow(covariance) * sqrt(.Machine$double.eps)`.
#' * `covariance` must be square, finite, symmetric within tolerance, with
#'   identical unique row and column names, and positive semi-definite within
#'   tolerance. Perfect correlation and zero-variance coordinates are valid.
#' * The observed block `R_oo` is inverted by Cholesky when it is positive
#'   definite and by an eigen pseudoinverse when it is only semi-definite. In
#'   the singular case every observed vector must lie in the support of its
#'   Gaussian: the component of `e_o` outside the range of `R_oo` must be
#'   within `tolerance * max(||e_o||, sqrt(lambda_max))` of zero. A vector
#'   outside the support (a non-zero value for a zero-variance coordinate, or
#'   two perfectly correlated coordinates that disagree) has probability zero
#'   and defines no conditional distribution, so it is an error rather than a
#'   pseudoinverse guess.
#' * The conditional covariance is symmetrized; eigenvalues in
#'   `[-tolerance * lambda_max, 0)` are projected to zero and anything more
#'   negative is an error. It is factored by Cholesky when positive definite
#'   (platform-deterministic) and by `V sqrt(D)` otherwise. A coordinate with
#'   zero conditional variance is returned at its conditional mean exactly.
#' * **RNG.** Every check runs before the first random number is drawn, so a
#'   rejected call leaves `.Random.seed` untouched. A successful call consumes
#'   exactly `n * length(sample_coordinates)` standard normals from
#'   [stats::rnorm()], in entity order and `sample_coordinates` order within
#'   an entity, whatever the observed patterns and even when some conditional
#'   variances are zero. Zero entities or zero sample coordinates consume
#'   nothing.
#'
#' @param covariance Numeric matrix with dimnames: one stratum of a block, as
#'   returned by [find_covariance_blocks()].
#' @param sample_coordinates Character vector of the coordinates to draw, a
#'   subset of the block's names with no duplicates. May be empty.
#' @param entity_keys Opaque entity identifiers: a vector or list (one entity
#'   per element) or a data frame (one entity per row). Only the count is used,
#'   and it fixes `n`. Rows of the result follow its order.
#' @param observed `NULL`, or a numeric matrix / data frame with `n` rows whose
#'   column names are block coordinates not in `sample_coordinates`. `NA`
#'   means "not observed for this entity"; every non-`NA` value must be
#'   finite. Stored and caller-fixed values both go here.
#' @param tolerance Relative tolerance; see the numerical contract. `NULL`
#'   uses the default.
#' @return A numeric `n` x `length(sample_coordinates)` matrix with the sample
#'   coordinates as column names, rows in `entity_keys` order.
#' @keywords internal
resolve_correlated_draws <- function(covariance, sample_coordinates,
                                     entity_keys, observed = NULL,
                                     tolerance = NULL) {

  # ── validation ───────────────────────────────────────────────────────────
  if (!is.matrix(covariance) || !is.numeric(covariance) ||
      nrow(covariance) != ncol(covariance) || nrow(covariance) < 1L) {
    stop("`covariance` must be a square numeric matrix.", call. = FALSE)
  }
  coords <- rownames(covariance)
  if (is.null(coords) || is.null(colnames(covariance)) ||
      !identical(coords, colnames(covariance)) || anyDuplicated(coords) ||
      anyNA(coords) || !all(nzchar(coords))) {
    stop("`covariance` must have identical, unique, non-empty row and column ",
         "names.", call. = FALSE)
  }
  n_block <- length(coords)
  if (!all(is.finite(covariance))) {
    stop("`covariance` must be finite (no NA, NaN or Inf).", call. = FALSE)
  }
  if (is.null(tolerance)) {
    tolerance <- n_block * sqrt(.Machine$double.eps)
  }
  if (!is.numeric(tolerance) || length(tolerance) != 1L ||
      !is.finite(tolerance) || tolerance < 0) {
    stop("`tolerance` must be a single non-negative number.", call. = FALSE)
  }
  lambda_max <- max(abs(covariance))
  tol_abs <- tolerance * lambda_max
  if (max(abs(covariance - t(covariance))) > tol_abs) {
    stop("`covariance` must be symmetric within tolerance (max discrepancy ",
         format(max(abs(covariance - t(covariance)))), ", tolerance ",
         format(tol_abs), ").", call. = FALSE)
  }
  covariance <- (covariance + t(covariance)) / 2
  if (lambda_max > 0) {
    ev <- eigen(covariance, symmetric = TRUE, only.values = TRUE)$values
    lambda_max <- max(ev)
    tol_abs <- tolerance * lambda_max
    if (min(ev) < -tol_abs) {
      stop("`covariance` is not positive semi-definite (smallest eigenvalue ",
           format(min(ev)), ", tolerance ", format(-tol_abs), ").",
           call. = FALSE)
    }
  }

  if (!is.character(sample_coordinates)) {
    stop("`sample_coordinates` must be a character vector.", call. = FALSE)
  }
  if (anyDuplicated(sample_coordinates)) {
    stop("`sample_coordinates` must not contain duplicates.", call. = FALSE)
  }
  bad <- setdiff(sample_coordinates, coords)
  if (length(bad) > 0L) {
    stop("`sample_coordinates` not in the block: ", .pvc_set(bad),
         ". The block is ", .pvc_set(coords), ".", call. = FALSE)
  }
  n <- if (is.data.frame(entity_keys)) nrow(entity_keys) else length(entity_keys)
  m <- length(sample_coordinates)

  if (is.null(observed)) {
    observed <- matrix(numeric(0), nrow = n, ncol = 0L)
  } else {
    if (is.data.frame(observed)) observed <- as.matrix(observed)
    if (is.matrix(observed) && is.logical(observed) && all(is.na(observed))) {
      storage.mode(observed) <- "double"     # cbind(A = c(NA, NA)) is logical
    }
    if (!is.matrix(observed) ||
        !(is.numeric(observed) || (ncol(observed) == 0L))) {
      stop("`observed` must be NULL or a numeric matrix / data frame.",
           call. = FALSE)
    }
    if (nrow(observed) != n) {
      stop("`observed` has ", nrow(observed), " rows but `entity_keys` has ",
           n, " entities.", call. = FALSE)
    }
    if (ncol(observed) > 0L) {
      ocn <- colnames(observed)
      if (is.null(ocn) || anyDuplicated(ocn)) {
        stop("`observed` must have unique column names.", call. = FALSE)
      }
      bad <- setdiff(ocn, coords)
      if (length(bad) > 0L) {
        stop("`observed` columns not in the block: ", .pvc_set(bad), ".",
             call. = FALSE)
      }
      both <- intersect(ocn, sample_coordinates)
      if (length(both) > 0L) {
        stop("Coordinates cannot be both observed and sampled: ",
             .pvc_set(both), ".", call. = FALSE)
      }
      if (any(is.infinite(observed) | is.nan(observed))) {
        stop("`observed` values must be finite or NA.", call. = FALSE)
      }
      storage.mode(observed) <- "double"
    }
  }

  out <- matrix(NA_real_, nrow = n, ncol = m,
                dimnames = list(NULL, sample_coordinates))
  if (n == 0L || m == 0L) return(out)

  # ── per-pattern conditional moments, all before any RNG ──────────────────
  # Pattern = the set of observed (non-NA) coordinates of an entity.
  if (ncol(observed) == 0L) {
    pattern_key <- rep("", n)
  } else {
    seen <- !is.na(observed)
    pattern_key <- apply(seen, 1L, function(r) paste(colnames(observed)[r],
                                                     collapse = "\r"))
  }
  groups <- split(seq_len(n), pattern_key)
  groups <- groups[order(names(groups))]

  R_ss <- covariance[sample_coordinates, sample_coordinates, drop = FALSE]
  plans <- lapply(seq_along(groups), function(g) {
    key <- names(groups)[[g]]          # "" = nothing observed
    idx <- groups[[g]]
    o   <- if (nzchar(key)) strsplit(key, "\r", fixed = TRUE)[[1L]] else character(0)
    if (length(o) == 0L) {
      return(list(idx = idx, mean = matrix(0, length(idx), m),
                  factor = .cd_factor(R_ss, tol_abs, sample_coordinates)))
    }
    E_o  <- observed[idx, o, drop = FALSE]
    R_oo <- covariance[o, o, drop = FALSE]
    R_so <- covariance[sample_coordinates, o, drop = FALSE]

    eig  <- eigen(R_oo, symmetric = TRUE)
    keep <- eig$values > tol_abs
    if (all(keep)) {
      # positive definite: Cholesky solve
      U <- chol(R_oo)
      W <- t(backsolve(U, backsolve(U, t(R_so), transpose = TRUE)))
    } else {
      # semi-definite (eigenvalues in (0, tol_abs] count as zero): the
      # observed vector must lie in the range of R_oo
      V_null <- eig$vectors[, !keep, drop = FALSE]
      off    <- sqrt(rowSums((E_o %*% V_null)^2))
      scale  <- pmax(sqrt(rowSums(E_o^2)), sqrt(lambda_max))
      bad    <- which(off > tolerance * scale)
      if (length(bad) > 0L) {
        stop("Observed values for ", .pvc_set(o), " lie outside the support ",
             "of their distribution for ", length(bad), " of ", length(idx),
             " entities in this pattern (e.g. entity ",
             .cd_entity_label(entity_keys, idx[bad[[1L]]]), ": ",
             paste0(o, " = ", format(E_o[bad[[1L]], ]), collapse = ", "),
             "). The covariance of ", .pvc_set(o), " is singular, so a ",
             "value off its support has probability zero and no conditional ",
             "distribution.", call. = FALSE)
      }
      V_k <- eig$vectors[, keep, drop = FALSE]
      pinv <- V_k %*% (t(V_k) / eig$values[keep])
      W <- R_so %*% pinv
    }
    C <- R_ss - W %*% t(R_so)
    list(idx = idx, mean = E_o %*% t(W),
         factor = .cd_factor(C, tol_abs, sample_coordinates))
  })

  # ── draw ─────────────────────────────────────────────────────────────────
  z <- matrix(stats::rnorm(n * m), nrow = n, ncol = m, byrow = TRUE)
  for (p in plans) {
    out[p$idx, ] <- p$mean + z[p$idx, , drop = FALSE] %*% p$factor
  }
  out
}

# A printable name for entity i, for error messages only.
.cd_entity_label <- function(entity_keys, i) {
  if (is.data.frame(entity_keys)) {
    paste0(names(entity_keys), " = ", vapply(entity_keys[i, , drop = FALSE],
                                             function(v) format(v[[1L]]),
                                             character(1)), collapse = ", ")
  } else if (is.atomic(entity_keys)) {
    format(entity_keys[[i]])
  } else {
    as.character(i)
  }
}

# Upper-triangular-style factor F with t(F) %*% F = C, so that z %*% F has
# covariance C. Cholesky when PD; otherwise eigen with tiny negatives zeroed.
.cd_factor <- function(C, tol_abs, coords) {
  C  <- (C + t(C)) / 2
  ev <- eigen(C, symmetric = TRUE)
  if (min(ev$values) < -tol_abs) {
    stop("The conditional covariance of ", .pvc_set(coords), " is materially ",
         "indefinite (smallest eigenvalue ", format(min(ev$values)),
         ", tolerance ", format(-tol_abs), ").", call. = FALSE)
  }
  if (min(ev$values) > tol_abs) {
    return(chol(C))
  }
  # Fix each eigenvector's sign (largest-magnitude component positive) so the
  # factor, and hence a seeded draw, does not depend on the LAPACK build.
  V <- ev$vectors
  s <- apply(V, 2L, function(v) sign(v[[which.max(abs(v))]]))
  V <- sweep(V, 2L, s, "*")
  d <- pmax(ev$values, 0)
  t(V %*% diag(sqrt(d), nrow = length(d)))
}
