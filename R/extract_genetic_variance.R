#' Measure the genetic covariance blocks of a cohort
#'
#' @description
#' Reports the genetic (co)variances a selected group of individuals actually
#' carries, in the shape of `trait_var_comp`: one row per
#' `(effect_name, trait_name_1, trait_name_2)`, every trait pair in both
#' orientations. So "did the model deliver its target?" is one
#' [dplyr::inner_join()] against the stored targets. Read-only: nothing is
#' written, and `ind_tgv` is not touched.
#'
#' `tbl` selects the individuals, exactly as for [add_tgv()]: any filtered
#' table with `id_ind` (`ind_meta` by generation or date, `ind_phenotype`,
#' `ind_ebv`, ...).
#'
#' @section Two anchors, two definitions:
#' * `anchor = "realised"` (default) measures the selected individuals: sample
#'   covariances (divisor `n - 1`) of each block's values. Under linkage
#'   disequilibrium (LD) or departures from Hardy-Weinberg equilibrium (HWE)
#'   the blocks correlate; those cross-block covariances are reported in
#'   `between_components`, so the block rows plus `between_components` sum
#'   exactly to `total`.
#' * `anchor = "genic"` is the HWE + linkage-equilibrium expectation at
#'   `base_tbl`'s allele frequencies. The blocks are orthogonal by
#'   construction, there is no `between_components` row, and `total` is the
#'   sum of the blocks. Available for fully decomposable models only (case 1
#'   below).
#'
#' `base_tbl` is like the generators' `base_tbl` (see [extract_allele_freq()])
#' with two differences: its `NULL` default is **the selected individuals**
#' (their whole-genotype frequencies, whatever table `tbl` is), not the founder
#' pool, so it follows the cohort as frequencies drift and does not reproduce a
#' generation target; and it is genic only (a non-`NULL` `base_tbl` with
#' `"realised"` is an error). An explicit `base_tbl` keeps
#' [extract_allele_freq()]'s semantics: an `ind_haplotype` filter selects allele
#' copies. To compare with a generation target, pass the generation base.
#'
#' @section What the realised blocks are, under LD:
#' The realised blocks are the sample covariances of the NOIA **contrast
#' components** at the cohort's frequencies (the source project's
#' `nonadd_decompose()`): each locus' heterozygosity is regressed on its own
#' dosage only, and the additive block is `(X - 2p) alpha` with
#' `alpha_j = a_j + b_j d_j + sum_l e_jl (2 p_l - 1)`. Under LE these are the
#' orthogonal least-squares components. Under LD they are not: dominance at
#' another locus, or a pair product, can still regress on the dosages, so the
#' `additive` row is **not** the variance of the cohort's joint least-squares
#' additive projection (the breeding values a regression of the genetic value
#' on all dosages would give), and `between_components = 0` does not certify
#' that it is. For example, one pair `e (g_1 - 1)(g_2 - 1)` on a cohort with
#' both loci at `p = 0.5` and HWE margins but in LD reports `additive = 0`,
#' while `lm(g ~ g_1 + g_2)` explains part of its variance. `"full"` means
#' every stored term has a supported shape, not that the cohort's additive
#' breeding-value variance has been recovered.
#'
#' @section What the function can decompose:
#' Each trait's stored model falls in one of three cases, reported in the
#' `decomposition` column.
#' 1. `"full"`: every term is common-scope (no line or parent-of-origin scope)
#'    on diploid-autosomal loci, and is an order-one `additive`, order-one
#'    `dominance`, order-one diploid `indicator` (any of the three genotype
#'    states), or two-member `additive x additive` term, of any owner and at
#'    any centre. These are converted to functional effects and re-projected
#'    on the NOIA (natural and orthogonal interactions) model at the cohort's
#'    frequencies, so both codings of one model give the same report.
#' 2. `"additive_only"`: every term is an order-one `additive` term and some
#'    have a line or parent-of-origin scope (crossbreeding, imprinting). The
#'    `additive` row is the **evaluated additive variance**, the covariance of
#'    the evaluated additive component, not a re-projection.
#' 3. `"partial"`: anything else. The covered terms are decomposed as in case
#'    1; everything else goes to `unpartitioned`, the variance of those terms'
#'    value alone (not `residual`, which in this package is environmental
#'    noise). In this case the `additive` row is the additive contrast
#'    component of the covered terms, not the breeding value of the whole
#'    model.
#'
#' Scope variants of one term compete (the most specific matching scope
#' wins), so a term family with any scoped variant is uncovered as a whole.
#' An off-diagonal row of two traits in different cases carries the less
#' complete label.
#'
#' @section Which rows appear:
#' Rows follow the statistical decomposition, not the stored term kinds.
#' `additive` is reported for every trait with a decomposed term (a one-locus
#' genotype surface equal to the dosage is all additive; dominance and
#' interaction terms induce additive effects too). `dominance` appears when
#' the converted model has a non-zero dominance coefficient,
#' `additive_by_additive` when it has a non-zero pair. A coefficient whose
#' stored contributions cancel to within their floating-point rounding (a
#' genotype surface that is exactly linear, written in any row order) counts
#' as zero. A block whose variance
#' is 0 in this cohort is reported as 0; a missing row means the model has no
#' such block, so a stored target for it falls out of the `inner_join()` and
#' into the [dplyr::anti_join()]. An off-diagonal block row appears only when
#' both traits have the block. `total` is always reported.
#'
#' A locus fixed in the cohort contributes nothing to the dominance or
#' interaction blocks, but an interaction with a fixed partner is still a real
#' additive effect (`e (g_1 - 1)(g_2 - 1) = e (g_1 - 1)` when `g_2 = 2`), and is
#' kept in `additive`.
#'
#' @section Comparing with a target:
#' The effect names are `trait_var_comp`'s (the interaction block is
#' `additive_by_additive`, not the `ind_tgv` component name `interaction`).
#' A mismatch with the stored target is information, never an error. The
#' comparison is like-for-like only when the trait's terms are all
#' `"generated"` **and** the anchor, the reference population, the measured
#' model and the calibrated scope match the generation call: a generated line
#' variant competes with common fallback terms, so a mixed cohort measures the
#' active combination; a restored or later-edited founder pool is not the
#' generation-time base. Stored diagonals are generation targets under their
#' calibration contract (genic or realised), not universally a genic total.
#'
#' `define_phenotype(prevalence = )` sums the stored diagonals. On a cohort
#' with LD, departures from HWE or drifted frequencies the realised `total`
#' differs from that sum by `between_components` and the drift in each block,
#' which this function shows. Even a matching variance does not guarantee the
#' prevalence: the threshold also assumes a near-normal liability, which a
#' skewed finite-locus model need not give.
#'
#' @param tbl A `tidybreed_table` from [get_table()], optionally filtered,
#'   selecting the individuals.
#' @param trait_name Character vector of traits. `NULL` (default) is every
#'   trait with stored genome-effect terms, in `trait_meta` order. A named
#'   trait without terms is an error.
#' @param base_tbl Genic anchor only: a filtered `tidybreed_table` whose allele
#'   frequencies define the expectation (`founder_haplotypes`, `ind_haplotype`
#'   copies, or individuals). `NULL` (default) uses the selected individuals.
#' @param anchor `"realised"` (default) or `"genic"`.
#'
#' @return A tibble with columns `effect_name`, `trait_name_1`,
#'   `trait_name_2`, `cov_value`, `n_ind`, `decomposition`, `anchor`, sorted by
#'   effect (`additive`, `dominance`, `additive_by_additive`, `unpartitioned`,
#'   `between_components`, `total`) and trait order. A `message()` names the
#'   population the estimates describe; the `anchor` column travels with the
#'   rows, but which cohort or base produced them is only in that message.
#'
#' @seealso [add_tgv()], [define_effect_cov_matrix()], [extract_allele_freq()].
#'
#' @examples
#' \dontrun{
#' # Realised blocks of generation 5
#' meas <- get_table(pop, "ind_meta") |>
#'   dplyr::filter(gen == 5L) |>
#'   extract_genetic_variance()
#'
#' # Target vs measured, population-wide targets (a line: line_name == "A")
#' targets <- get_table(pop, "trait_var_comp") |>
#'   dplyr::filter(is.na(line_name)) |>
#'   dplyr::collect() |>
#'   dplyr::select(effect_name, trait_name_1, trait_name_2, target = cov_value)
#' dplyr::inner_join(targets, meas,
#'                   by = c("effect_name", "trait_name_1", "trait_name_2"))
#' dplyr::anti_join(targets, meas,      # targets with no measured block
#'                  by = c("effect_name", "trait_name_1", "trait_name_2"))
#'
#' # Genic: a common-scope model generated with anchor = "genic" on the founder
#' # pool, measured at the founder pool's frequencies
#' get_table(pop, "ind_meta") |>
#'   extract_genetic_variance(anchor = "genic",
#'                            base_tbl = get_table(pop, "founder_haplotypes"))
#' }
#' @export
extract_genetic_variance <- function(tbl, trait_name = NULL, base_tbl = NULL,
                                     anchor = c("realised", "genic")) {
  if (!inherits(tbl, "tidybreed_table")) {
    stop("'tbl' must be a tidybreed_table from get_table() selecting ",
         "individuals.", call. = FALSE)
  }
  anchor <- match.arg(anchor)
  pop  <- tbl$pop
  validate_tidybreed_pop(pop)
  conn <- pop$db_conn
  if (!is.null(base_tbl)) {
    if (anchor == "realised") {
      stop("'base_tbl' applies to anchor = \"genic\" only. Under \"realised\" ",
           "the projection uses the frequencies and regressions of the ",
           "individuals 'tbl' selects; drop 'base_tbl', or use ",
           "anchor = \"genic\".", call. = FALSE)
    }
    .validate_base_tbl(base_tbl, pop)
  }

  ids <- resolve_subset_ids(tbl, "extract_genetic_variance()",
                            all_if_null = TRUE)
  n <- length(ids)
  if (n < 2L) {
    stop("extract_genetic_variance() needs at least 2 individuals (a sample ",
         "covariance needs two); 'tbl' selects ", n, ".", call. = FALSE)
  }

  traits <- .egv_traits(conn, trait_name)
  model  <- .gev_read_model(conn, traits)
  cls    <- .egv_classify(conn, model, traits)

  if (anchor == "genic" && any(cls$case != "full")) {
    bad <- traits[cls$case != "full"]
    stop("anchor = \"genic\" needs a fully decomposable model (common-scope ",
         "additive / dominance / one-locus indicator / additive x additive ",
         "terms on diploid-autosomal loci). Trait(s) ",
         paste0("'", bad, "'", collapse = ", "), " are '",
         paste(unique(cls$case[cls$case != "full"]), collapse = "', '"),
         "'. Use anchor = \"realised\", which measures any model.",
         call. = FALSE)
  }

  # Knowable sizes, before any evaluation.
  proj  <- cls$case != "additive_only"
  cov_terms <- model$terms$id_genome_effect[
    cls$covered & model$terms$trait_name %in% traits[proj]]
  cov_loci <- sort(unique(model$members$locus_id[
    model$members$id_genome_effect %in% cov_terms]))
  if (anchor == "realised" && length(cov_loci) > 0L) {
    .dosage_guard(n, length(cov_loci), "extract_genetic_variance()",
                  .EGV_SIZE_FIX)
  }

  res <- if (anchor == "realised") {
    .egv_realised(conn, ids, traits, model, cls, cov_loci)
  } else {
    .egv_genic(conn, ids, traits, model, cls, cov_loci, base_tbl)
  }

  message(.egv_message(anchor, n, base_tbl, traits, cls))
  .egv_rows(res, traits, cls, n, anchor)
}

.EGV_SIZE_FIX <- paste0(
  "Select fewer individuals with 'tbl', or use anchor = \"genic\" ",
  "(fully decomposable models).")

.EGV_EFFECT_ORDER <- c("additive", "dominance", "additive_by_additive",
                       "unpartitioned", "between_components", "total")

.EGV_CASE_RANK <- c(full = 1L, additive_only = 2L, partial = 3L)


# -- Traits -------------------------------------------------------------------

#' The traits to report: named ones (each must have terms), or every trait
#' with terms in `trait_meta` order
#'
#' Deliberately not `.gev_resolve_traits(conn, NULL)`, which returns every
#' trait (with or without terms) and which `add_tgv()` relies on.
#'
#' @keywords internal
#' @noRd
.egv_traits <- function(conn, trait_name) {
  with_terms <- DBI::dbGetQuery(conn, paste0(
    "SELECT t.trait_name FROM trait_meta t WHERE EXISTS (SELECT 1 FROM ",
    "genome_effects e WHERE e.trait_name = t.trait_name) ORDER BY t.id_trait"
  ))$trait_name
  if (is.null(trait_name)) {
    if (length(with_terms) == 0L) {
      stop("No trait has genome-effect terms. Call define_additive_effects() ",
           "or define_genome_effect_terms() first.", call. = FALSE)
    }
    return(with_terms)
  }
  trait_name <- .gev_resolve_traits(conn, trait_name)
  if (anyDuplicated(trait_name)) {
    stop("'trait_name' names a trait twice.", call. = FALSE)
  }
  none <- setdiff(trait_name, with_terms)
  if (length(none) > 0L) {
    stop("No genome-effect terms for trait(s) ",
         paste0("'", none, "'", collapse = ", "), ", so there is no genetic ",
         "variance to measure.", call. = FALSE)
  }
  trait_name
}


# -- Classification -----------------------------------------------------------

#' Locus ids on chromosomes that are `1, 1` for both sexes in every line
#'
#' @keywords internal
#' @noRd
.egv_diploid_loci <- function(conn) {
  lines <- DBI::dbGetQuery(conn, paste0(
    "SELECT DISTINCT line_name FROM chr_inheritance ",
    "WHERE line_name IS NOT NULL"))$line_name
  ok <- NULL
  for (ln in c(list(NULL), as.list(lines))) {
    for (sx in c("M", "F")) {
      inh <- resolve_chr_inheritance(conn, sx, ln)
      this <- stats::setNames(inh$from_parent_1 == 1L & inh$from_parent_2 == 1L,
                              inh$chr_name)
      ok <- if (is.null(ok)) this else ok & this[names(ok)]
    }
  }
  chrs <- names(ok)[ok]
  if (length(chrs) == 0L) return(integer(0))
  DBI::dbGetQuery(conn, paste0(
    "SELECT locus_id FROM genome_meta WHERE chr_name IN (",
    sql_in_list(chrs, what = "chromosome name"), ")"))$locus_id
}

#' Covered terms and the case of each trait
#'
#' A term is **covered** when it has one of the decomposable shapes, every
#' member locus is diploid-autosomal, and it has no origin rows. A family
#' (`family_key`: trait, owner, member loci and states) is covered only as a
#' whole: its scope variants compete, so a common variant must never be
#' projected apart from scoped siblings.
#'
#' @return `list(covered, case)`: `covered` logical per row of `model$terms`,
#'   `case` character per trait (`"full"`, `"additive_only"`, `"partial"`).
#' @keywords internal
#' @noRd
.egv_classify <- function(conn, model, traits) {
  terms <- model$terms
  mem   <- model$members
  dip   <- .egv_diploid_loci(conn)
  mid   <- as.character(mem$id_genome_effect)
  n_mem <- tabulate(match(mem$id_genome_effect, terms$id_genome_effect),
                    nbins = nrow(terms))
  all_dip <- tapply(mem$locus_id %in% dip, mid, all)
  n_add   <- tapply(mem$contrast_name == "additive", mid, sum)
  key <- as.character(terms$id_genome_effect)
  m1  <- mem[!duplicated(mem$id_genome_effect), , drop = FALSE]
  m1  <- m1[match(terms$id_genome_effect, m1$id_genome_effect), , drop = FALSE]

  shape <- (n_mem == 1L & m1$contrast_name %in% c("additive", "dominance")) |
    (n_mem == 1L & m1$contrast_name == "indicator" &
       !is.na(m1$copy_count_value) & m1$copy_count_value == 2L &
       !is.na(m1$dosage_value) & m1$dosage_value %in% 0:2) |
    (n_mem == 2L & n_add[key] == 2L)
  unscoped <- !terms$id_genome_effect %in% model$origins$id_genome_effect
  term_ok  <- shape & as.logical(all_dip[key]) & unscoped
  fam_ok   <- tapply(term_ok, terms$family_key, all)
  covered  <- unname(as.logical(fam_ok[terms$family_key]))

  additive1 <- n_mem == 1L & m1$contrast_name == "additive"
  case <- vapply(traits, function(t) {
    w <- terms$trait_name == t
    if (all(covered[w])) "full"
    else if (all(additive1[w])) "additive_only"
    else "partial"
  }, character(1), USE.NAMES = FALSE)
  list(covered = covered, case = case)
}


# -- Canonical coefficients ---------------------------------------------------

#' One trait set's covered terms as functional coefficient matrices
#'
#' `.stored_to_functional()` per trait, owners summed, laid out on the shared
#' covered-locus index: `A`, `D` are `m x k`, `pairs` the union of canonical
#' pairs (column indices into `loci`) with `E` its `(r, k)` coefficients. A
#' trait not projected (case 2) has zero columns.
#'
#' @keywords internal
#' @noRd
.egv_coefficients <- function(model, cls, traits, loci) {
  k <- length(traits); m <- length(loci)
  A <- D <- matrix(0, m, k, dimnames = list(NULL, traits))
  pk <- list()
  proj <- cls$case != "additive_only"
  has_cov <- stats::setNames(logical(k), traits)
  for (j in seq_len(k)) {
    if (!proj[j]) next
    w <- cls$covered & model$terms$trait_name == traits[j]
    if (!any(w)) next
    has_cov[j] <- TRUE
    tt <- model$terms[w, , drop = FALSE]
    mm <- model$members[model$members$id_genome_effect %in%
                          tt$id_genome_effect, , drop = FALSE]
    f <- .stored_to_functional(tt, mm)
    idx <- match(as.integer(names(f$a)), loci)
    A[idx, j] <- f$a
    D[idx, j] <- f$d
    if (nrow(f$pairs) > 0L) {
      f$pairs$trait <- j
      pk[[length(pk) + 1L]] <- f$pairs
    }
  }
  pr <- if (length(pk) > 0L) do.call(rbind, pk) else
    data.frame(locus_id_1 = integer(0), locus_id_2 = integer(0),
               e = numeric(0), trait = integer(0))
  key <- paste(pr$locus_id_1, pr$locus_id_2)
  ukey <- unique(key[order(pr$locus_id_1, pr$locus_id_2)])
  E <- matrix(0, length(ukey), k, dimnames = list(NULL, traits))
  if (nrow(pr) > 0L) E[cbind(match(key, ukey), pr$trait)] <- pr$e
  ks <- strsplit(ukey, " ", fixed = TRUE)
  pairs <- cbind(match(as.integer(vapply(ks, `[`, "", 1)), loci),
                 match(as.integer(vapply(ks, `[`, "", 2)), loci))
  if (length(ukey) == 0L) pairs <- matrix(integer(0), 0, 2)
  list(A = A, D = D, E = E, pairs = pairs, has_cov = has_cov)
}

#' The induced additive coefficients: alpha = a + b d + sum_l e c_l
#'
#' @keywords internal
#' @noRd
.egv_alpha <- function(co, b, cc) {
  alpha <- co$A + b * co$D
  if (nrow(co$pairs) > 0L) {
    for (r in seq_len(nrow(co$pairs))) {
      k <- co$pairs[r, 1]; l <- co$pairs[r, 2]
      alpha[k, ] <- alpha[k, ] + cc[l] * co$E[r, ]
      alpha[l, ] <- alpha[l, ] + cc[k] * co$E[r, ]
    }
  }
  alpha
}

#' A x A values in deterministic pair chunks
#'
#' `sum_r e_r * centre(z_k z_l)` accumulated chunk by chunk in canonical pair
#' order, so no `n x r` matrix exists. The source's `nonadd_covariates()`
#' builds that matrix whole; with all pairs of 500 loci it is ~2 GB.
#'
#' @param Z_A Centred dosages, `n x m`.
#' @param pairs Two-column matrix of column indices into `Z_A`, one row per pair.
#' @param E `(r, k)` coefficients.
#' @param chunk Pairs per chunk.
#' @return `n x k` matrix.
#' @keywords internal
#' @noRd
.egv_aa_values <- function(Z_A, pairs, E, chunk) {
  out <- matrix(0, nrow(Z_A), ncol(E), dimnames = list(NULL, colnames(E)))
  r <- nrow(pairs)
  if (r == 0L) return(out)
  chunk <- max(1L, as.integer(chunk))
  for (s in seq(1L, r, by = chunk)) {
    idx <- s:min(r, s + chunk - 1L)
    Z <- Z_A[, pairs[idx, 1], drop = FALSE] * Z_A[, pairs[idx, 2], drop = FALSE]
    Z <- sweep(Z, 2L, colMeans(Z), "-")
    out <- out + Z %*% E[idx, , drop = FALSE]
  }
  out
}

#' Pairs per chunk: one chunk's `n x chunk` matrix stays under the cell limit
#'
#' @keywords internal
#' @noRd
.egv_pair_chunk <- function(n) {
  max(1L, as.integer(floor(QTL_REALISED_MAX_CELLS / max(1L, n))))
}


# -- Realised -----------------------------------------------------------------

#' Per-individual block values and totals for the selected cohort
#'
#' @return `list(blocks, avail, G)`: `blocks` a named list of centred `n x k`
#'   value matrices, `avail` a logical block-by-trait matrix, `G` the centred
#'   evaluated totals.
#' @keywords internal
#' @noRd
.egv_realised <- function(conn, ids, traits, model, cls, cov_loci) {
  n <- length(ids); k <- length(traits)
  res <- .gev_evaluate(conn, ids, traits, model)
  G   <- .egv_by_trait(res, ids, traits)
  miss <- vapply(traits, function(t) sum(!ids %in% res$id_ind[res$trait_name == t]),
                 integer(1))
  if (any(miss > 0L)) {
    stop("Some selected individuals have no genetic value for a trait (no ",
         "term matches any of their allele copies): ",
         paste0(traits[miss > 0L], " (", miss[miss > 0L], ")", collapse = ", "),
         ". Every row is computed on the same individuals; select a narrower ",
         "'tbl', or name fewer traits.", call. = FALSE)
  }

  zero  <- matrix(0, n, k, dimnames = list(NULL, traits))
  blocks <- list(additive = zero, dominance = zero, additive_by_additive = zero,
                 unpartitioned = zero)
  avail <- matrix(FALSE, length(blocks), k,
                  dimnames = list(names(blocks), traits))

  if (length(cov_loci) > 0L) {
    tmp <- paste0("_egv_ids_", as.character(round(as.numeric(Sys.time()) * 1000)))
    duckdb::duckdb_register(conn, tmp, data.frame(id_ind = ids,
                                                  stringsAsFactors = FALSE))
    on.exit(try(duckdb::duckdb_unregister(conn, tmp), silent = TRUE),
            add = TRUE)
    dx <- .collect_dosages(
      conn, paste0("SELECT id_ind FROM ", tmp), cov_loci,
      who = "extract_genetic_variance()", where = "'tbl'",
      size_fix = .EGV_SIZE_FIX,
      incomplete = paste0(
        "some selected individuals lack both allele copies at a decomposed ",
        "locus (individuals without haplotypes?)."))
    X <- dx$X[match(ids, dx$id_ind), , drop = FALSE]
    co <- .egv_coefficients(model, cls, traits, dx$locus_id)

    # The source's nonadd_covariates(), realised anchor, without its dense
    # anchor matrices.
    p   <- colMeans(X) / 2
    Z_A <- sweep(X, 2L, 2 * p, "-")
    Wc  <- sweep((X == 1) * 1, 2L, colMeans((X == 1) * 1), "-")
    vx  <- colSums(Z_A^2) / n
    cwx <- colSums(Wc * Z_A) / n
    b   <- ifelse(vx > 0, cwx / vx, 0)
    cc  <- 2 * p - 1
    Z_D <- Wc - sweep(Z_A, 2L, b, "*")

    alpha <- .egv_alpha(co, b, cc)
    blocks$additive  <- Z_A %*% alpha
    blocks$dominance <- Z_D %*% co$D
    blocks$additive_by_additive <- .egv_aa_values(Z_A, co$pairs, co$E,
                                                  .egv_pair_chunk(n))
    avail["additive", ]  <- co$has_cov
    avail["dominance", ] <- co$has_cov & colSums(co$D != 0) > 0
    avail["additive_by_additive", ] <- co$has_cov & colSums(co$E != 0) > 0
  }

  for (j in seq_len(k)) {
    t <- traits[j]
    if (cls$case[j] == "additive_only") {
      a <- res[res$trait_name == t & res$component_name == "additive", ,
               drop = FALSE]
      v <- numeric(n)
      v[match(a$id_ind, ids)] <- a$tgv_value
      blocks$additive[, j] <- v - mean(v)
      avail["additive", j] <- TRUE
    } else if (cls$case[j] == "partial") {
      w <- !cls$covered & model$terms$trait_name == t
      sub <- .egv_submodel(model, w)
      u <- .gev_evaluate(conn, ids, t, sub)
      U <- .egv_by_trait(u, ids, t)   # absent individuals contribute 0
      blocks$unpartitioned[, j] <- U[, 1]
      avail["unpartitioned", j] <- TRUE
    }
  }

  # The identities say the blocks add up to the centred evaluated total. If
  # they do not, the canonicalisation is wrong: refuse to report.
  S <- Reduce(`+`, blocks)
  for (j in seq_len(k)) {
    tol <- 1e-10 + 1e-8 * max(abs(G[, j]))
    if (max(abs(S[, j] - G[, j])) > tol) {
      stop("Internal error: the decomposed blocks of trait '", traits[j],
           "' do not add up to its evaluated genetic value (max difference ",
           signif(max(abs(S[, j] - G[, j])), 3), "). Please report this.",
           call. = FALSE)
    }
  }
  list(blocks = blocks, avail = avail, G = G)
}

#' Evaluator rows -> centred `n x k` totals, aligned to `ids` by name
#'
#' Components are added in a fixed order; an individual with no row is 0.
#'
#' @keywords internal
#' @noRd
.egv_by_trait <- function(res, ids, traits) {
  out <- matrix(0, length(ids), length(traits), dimnames = list(NULL, traits))
  res <- res[order(res$trait_name, res$id_ind, res$component_name), ,
             drop = FALSE]
  for (j in seq_along(traits)) {
    r <- res[res$trait_name == traits[j], , drop = FALSE]
    if (nrow(r) == 0L) next
    tot <- rowsum(r$tgv_value, r$id_ind, reorder = FALSE)
    out[match(rownames(tot), ids), j] <- tot[, 1]
    out[, j] <- out[, j] - mean(out[, j])
  }
  out
}

#' Restrict a `.gev_read_model()` result to some terms
#'
#' @keywords internal
#' @noRd
.egv_submodel <- function(model, keep) {
  ids <- model$terms$id_genome_effect[keep]
  list(terms   = model$terms[keep, , drop = FALSE],
       members = model$members[model$members$id_genome_effect %in% ids, ,
                               drop = FALSE],
       origins = model$origins[model$origins$id_genome_effect %in% ids, ,
                               drop = FALSE])
}


# -- Genic --------------------------------------------------------------------

#' Closed-form genic blocks at the base frequencies (case 1 only)
#'
#' @return `list(cov, avail)`: `cov` a named list of `k x k` matrices.
#' @keywords internal
#' @noRd
.egv_genic <- function(conn, ids, traits, model, cls, cov_loci, base_tbl) {
  p  <- .egv_base_freq(conn, ids, cov_loci, base_tbl)
  co <- .egv_coefficients(model, cls, traits, cov_loci)
  q  <- 1 - p
  w  <- 2 * p * q
  alpha <- .egv_alpha(co, q - p, 2 * p - 1)
  wAA <- if (nrow(co$pairs) > 0L) w[co$pairs[, 1]] * w[co$pairs[, 2]] else
    numeric(0)
  cv <- list(
    additive  = crossprod(alpha, w * alpha),
    dominance = crossprod(co$D, w^2 * co$D),
    additive_by_additive = crossprod(co$E, wAA * co$E))
  cv$total <- cv$additive + cv$dominance + cv$additive_by_additive
  avail <- rbind(additive  = co$has_cov,
                 dominance = co$has_cov & colSums(co$D != 0) > 0,
                 additive_by_additive = co$has_cov & colSums(co$E != 0) > 0)
  list(cov = cv, avail = avail)
}

#' Base allele frequencies at the covered loci
#'
#' `base_tbl = NULL`: the selected individuals' whole-genotype frequencies,
#' whatever table selected them (an `ind_haplotype` filter picks individuals
#' here, never copies). Otherwise [extract_allele_freq()] with its own
#' semantics. A locus without copies in the base is an error.
#'
#' @keywords internal
#' @noRd
.egv_base_freq <- function(conn, ids, loci, base_tbl) {
  if (is.null(base_tbl)) {
    tmp <- paste0("_egv_base_", as.character(round(as.numeric(Sys.time()) * 1000)))
    duckdb::duckdb_register(conn, tmp, data.frame(id_ind = ids,
                                                  stringsAsFactors = FALSE))
    on.exit(try(duckdb::duckdb_unregister(conn, tmp), silent = TRUE),
            add = TRUE)
    f <- DBI::dbGetQuery(conn, paste0(
      "SELECT h.locus_id, AVG(CAST(h.allele AS DOUBLE)) AS allele_freq ",
      "FROM ind_haplotype h JOIN ", tmp, " i USING (id_ind) ",
      "WHERE h.locus_id IN (", paste(as.integer(loci), collapse = ", "), ") ",
      "GROUP BY h.locus_id"))
  } else {
    f <- as.data.frame(extract_allele_freq(base_tbl))
  }
  p <- f$allele_freq[match(loci, f$locus_id)]
  bad <- is.na(p)
  if (any(bad)) {
    nm <- DBI::dbGetQuery(conn, paste0(
      "SELECT locus_name FROM genome_meta WHERE locus_id IN (",
      paste(as.integer(loci[bad]), collapse = ", "), ") ORDER BY locus_id"))
    stop("The base has no allele copies at decomposed locus/loci ",
         paste0("'", utils::head(nm$locus_name, 5), "'", collapse = ", "),
         if (sum(bad) > 5L) paste0(" (", sum(bad), " total)") else "",
         ". A missing frequency is not 0 or 1; pass a base_tbl that covers ",
         "every locus with an effect.", call. = FALSE)
  }
  p
}


# -- Output -------------------------------------------------------------------

#' Assemble the long tibble
#'
#' @keywords internal
#' @noRd
.egv_rows <- function(res, traits, cls, n, anchor) {
  k <- length(traits)
  rank <- .EGV_CASE_RANK[cls$case]
  dec <- function(i, j) names(.EGV_CASE_RANK)[max(rank[i], rank[j])]
  rows <- list()
  add <- function(effect, M, keep = NULL) {
    for (i in seq_len(k)) for (j in seq_len(k)) {
      if (!is.null(keep) && !(keep[i] && keep[j])) next
      rows[[length(rows) + 1L]] <<- data.frame(
        effect_name = effect, trait_name_1 = traits[i],
        trait_name_2 = traits[j], cov_value = unname(M[i, j]),
        n_ind = as.integer(n), decomposition = dec(i, j), anchor = anchor,
        stringsAsFactors = FALSE)
    }
  }
  if (anchor == "realised") {
    bl <- res$blocks
    present <- names(bl)[rowSums(res$avail) > 0]
    for (b in present) {
      add(b, stats::cov(bl[[b]], bl[[b]]), res$avail[b, ])
    }
    # Every ordered pair of different blocks, both trait orientations,
    # computed directly from the value matrices (never total minus blocks).
    btw <- matrix(0, k, k)
    for (b in present) for (b2 in present) {
      if (b != b2) btw <- btw + stats::cov(bl[[b]], bl[[b2]])
    }
    add("between_components", btw)
    add("total", stats::cov(res$G, res$G))
  } else {
    for (b in rownames(res$avail)) {
      if (any(res$avail[b, ])) add(b, res$cov[[b]], res$avail[b, ])
    }
    add("total", res$cov$total)
  }
  out <- do.call(rbind, rows)
  out <- out[order(match(out$effect_name, .EGV_EFFECT_ORDER),
                   match(out$trait_name_1, traits),
                   match(out$trait_name_2, traits)), , drop = FALSE]
  rownames(out) <- NULL
  tibble::as_tibble(out)
}

#' The population the estimates describe (decision D1)
#'
#' @keywords internal
#' @noRd
.egv_message <- function(anchor, n, base_tbl, traits, cls) {
  head <- if (anchor == "realised") {
    paste0("Realised genetic covariances of the ", n, " selected individuals ",
           "(frequencies and regressions are the cohort's own).")
  } else if (is.null(base_tbl)) {
    paste0("Genic (HWE + linkage equilibrium) expectation at the ",
           "whole-genotype allele frequencies of the ", n, " selected ",
           "individuals (base_tbl = NULL follows the cohort; pass the ",
           "generation base to compare with a target).")
  } else {
    paste0("Genic (HWE + linkage equilibrium) expectation at the allele ",
           "frequencies of '", base_tbl$table_name, "'",
           if (length(base_tbl$pending_filter) > 0L) " (filtered)" else "",
           ", an explicit base (", if (base_tbl$table_name %in%
             c("founder_haplotypes", "ind_haplotype")) "allele copies" else
             "individuals", ").")
  }
  ao <- traits[cls$case == "additive_only"]
  pa <- traits[cls$case == "partial"]
  paste0(head,
         if (length(ao)) paste0(" Evaluated additive variance (scoped ",
                                "additive model, not re-projected): ",
                                paste(ao, collapse = ", "), ".") else "",
         if (length(pa)) paste0(" Partly decomposed; the rest is in ",
                                "'unpartitioned': ", paste(pa, collapse = ", "),
                                ".") else "")
}
