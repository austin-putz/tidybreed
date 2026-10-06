#' Build `terms` rows for a functional or Cockerham (a, d) pair
#'
#' @description
#' Expands one locus's additive and dominance coefficients into the two member
#' rows [define_genome_effect_terms()] takes. There is no third table and no separate
#' storage mode — functional coding is `additive` with `center_value = 0.5` plus
#' an `indicator` on the heterozygous state, Cockerham coding is `additive` plus
#' `dominance`, both centred at `p`.
#'
#' * `coding = "functional"` — the genotypic value is `a * (g - 1) + d * 1[g = 1]`,
#'   where `g` is the dosage of allele 1. Values are **absolute**, so the
#'   genetic component has a non-zero mean by construction.
#' * `coding = "cockerham"` — `a` is the average effect `alpha` and `d` the
#'   dominance deviation `delta` on the HWE-orthogonal contrast
#'   (`-2p^2`, `2pq`, `-2q^2`). Values are deviations, with mean 0.
#'
#' @section The reported mean is not written anywhere:
#' For functional coding the implied mean of the **genetic component** is
#' `mu = a(p - q) + 2pq*d`. It is reported and stored nowhere — in particular
#' not in `phenotype_meta.mean`, which would double-count it once non-additive
#' genetic values reach the phenotype layer, because the raw genetic values
#' already have expectation `mu`. The running total across loci is reported only
#' because every term here is a single-locus main effect; no such total exists
#' for an epistatic term, whose expectation depends on joint genotype
#' frequencies and LD.
#'
#' @param locus_name Character vector of locus names.
#' @param a Numeric additive coefficient(s), recycled to `locus_name`.
#' @param d Numeric dominance coefficient(s), recycled. Default `0`.
#' @param p Numeric allele-1 frequency per locus, recycled. Required for
#'   `"cockerham"`; used for the reported mean under `"functional"`.
#' @param coding `"functional"` (default) or `"cockerham"`.
#' @param effect_name Optional label carried onto both members of each locus.
#' @param report Logical; print the implied genetic mean. Default `TRUE`.
#'
#' @return A `terms` data frame with the fixed builder column set (see
#'   [aa_terms()]): one row per non-zero coefficient. `term_id` encodes the
#'   builder, the locus and `a` / `d`, so it never collides with another
#'   builder's ids when outputs are `rbind()`-ed.
#'
#' @seealso [define_genome_effect_terms()], [aa_terms()], [genotype_terms()].
#'
#' @examples
#' \dontrun{
#' tt <- ad_terms("Locus_10", a = 0.4, d = 0.2, p = 0.3)
#' pop <- pop |> define_genome_effect_terms("ADG", tt, effect_owner = "functional")
#' }
#' @export
ad_terms <- function(locus_name, a, d = 0, p,
                     coding      = c("functional", "cockerham"),
                     effect_name = NULL,
                     report      = TRUE) {
  coding <- match.arg(coding)
  if (!is.character(locus_name) || length(locus_name) == 0L) {
    stop("'locus_name' must be a non-empty character vector.", call. = FALSE)
  }
  n <- length(locus_name)
  if (anyDuplicated(locus_name)) {
    stop("'locus_name' must not repeat a locus: ",
         paste0("'", unique(locus_name[duplicated(locus_name)]), "'",
                collapse = ", "), ".", call. = FALSE)
  }
  a <- .ad_recycle(a, n, "a")
  d <- .ad_recycle(d, n, "d")
  if (missing(p)) {
    stop("'p' (the allele-1 frequency at each locus) is required: it is the ",
         "centre of the additive contrast under Cockerham coding and the ",
         "basis of the reported mean under functional coding.", call. = FALSE)
  }
  p <- .ad_recycle(p, n, "p")
  if (any(p < 0 | p > 1)) stop("'p' must be between 0 and 1.", call. = FALSE)
  lab <- if (is.null(effect_name)) NA_character_ else
    .ad_recycle_chr(effect_name, n, "effect_name")

  centre <- if (coding == "functional") rep(0.5, n) else p
  add <- .terms_frame(
    term_id       = .term_id("ad", as.list(locus_name), "a"),
    locus_name    = locus_name,
    contrast_name = "additive",
    center_value  = centre,
    genome_value  = a,
    effect_name   = lab)
  dom <- if (coding == "functional") {
    .terms_frame(term_id = .term_id("ad", as.list(locus_name), "d"),
                 locus_name = locus_name, contrast_name = "indicator",
                 copy_count_value = 2L, dosage_value = 1L,
                 genome_value = d, effect_name = lab)
  } else {
    .terms_frame(term_id = .term_id("ad", as.list(locus_name), "d"),
                 locus_name = locus_name, contrast_name = "dominance",
                 center_value = p, genome_value = d, effect_name = lab)
  }
  out <- rbind(add[a != 0, , drop = FALSE], dom[d != 0, , drop = FALSE])
  if (nrow(out) == 0L) {
    stop("Every a and d is zero, so there is no term to write.", call. = FALSE)
  }

  if (isTRUE(report)) {
    mu <- if (coding == "functional") a * (2 * p - 1) + 2 * p * (1 - p) * d
          else rep(0, n)
    message("Implied genetic mean (", coding, " coding), reported only — ",
            "written to no table: ",
            paste0(locus_name, " mu = ", format(round(mu, 6)),
                   collapse = "; "),
            ". Running total over these single-locus main effects: ",
            format(round(sum(mu), 6)), ".")
  }
  rownames(out) <- NULL
  out
}

#' @keywords internal
#' @noRd
.ad_recycle <- function(x, n, what, of = "locus_name") {
  if (!is.numeric(x) || anyNA(x) || any(!is.finite(x))) {
    stop("'", what, "' must be numeric, finite and free of NA.", call. = FALSE)
  }
  if (length(x) == 1L) return(rep(as.numeric(x), n))
  if (length(x) != n) {
    stop("'", what, "' must have length 1 or length(", of, ") (", n, ").",
         call. = FALSE)
  }
  as.numeric(x)
}

#' @keywords internal
#' @noRd
.ad_recycle_chr <- function(x, n, what, of = "locus_name") {
  if (!is.character(x)) stop("'", what, "' must be character.", call. = FALSE)
  if (length(x) == 1L) return(rep(x, n))
  if (length(x) != n) {
    stop("'", what, "' must have length 1 or length(", of, ") (", n, ").",
         call. = FALSE)
  }
  x
}

#' Build `terms` rows from a table of genotype values
#'
#' @description
#' Turns a genotype-by-value table into `indicator` terms: one term per row of
#' `genotypes`, one member per locus column. This is indicator completeness made
#' concrete — an arbitrary function of a finite set of genotype states *is* a
#' sum of indicator terms, so a hand-entered surface needs no second
#' representation, only rows.
#'
#' A cell you leave out is a term you did not write, contributing zero. Ask for
#' the opposite with `define_genome_effect_terms(require_complete = TRUE)`, which
#' then demands every reachable `(copy_count, dosage)` state.
#'
#' @param genotypes A data frame whose columns are named by `locus_name` and
#'   hold **dosages** of allele 1, one row per genotype combination. A
#'   `copy_count` list may accompany it for variable-copy loci.
#' @param value Numeric vector of genotypic values, one per row of `genotypes`.
#' @param copy_count Optional named list or vector giving `copy_count_value` per
#'   locus column. Omitted, it is inferred by the writer at diploid-autosomal
#'   loci; a variable-copy locus needs it, or a `copy_count` column-shaped data
#'   frame matching `genotypes`.
#' @param drop_zero Logical; drop rows whose `value` is exactly 0, since they
#'   contribute nothing. Default `TRUE`.
#' @param effect_name Optional label carried onto every term.
#'
#' @return A `terms` data frame with the fixed builder column set (see
#'   [aa_terms()]) and `nrow(genotypes) * ncol(genotypes)` rows (before
#'   `drop_zero`). Each row of `genotypes` is one term; its `term_id` encodes
#'   the surface's loci and the row number, so surfaces over different loci
#'   never share a `term_id` when bound together.
#'
#' @seealso [define_genome_effect_terms()], [ad_terms()], [aa_terms()].
#'
#' @examples
#' \dontrun{
#' cells <- expand.grid(Locus_10 = 0:2, Locus_44 = 0:2)
#' vals  <- c(0, 0, 0,  0, 1.4, 2.1,  0, 2.1, 3.6)
#' pop <- pop |> define_genome_effect_terms(
#'   "ADG", genotype_terms(cells, vals), effect_owner = "epistasis_AxA")
#' }
#' @export
genotype_terms <- function(genotypes, value, copy_count = NULL,
                           drop_zero = TRUE, effect_name = NULL) {
  if (!is.data.frame(genotypes) || ncol(genotypes) == 0L ||
      nrow(genotypes) == 0L) {
    stop("'genotypes' must be a data frame with at least one row and one ",
         "locus column.", call. = FALSE)
  }
  if (!is.numeric(value) || length(value) != nrow(genotypes) || anyNA(value)) {
    stop("'value' must be a numeric vector with one non-missing entry per row ",
         "of 'genotypes' (", nrow(genotypes), ").", call. = FALSE)
  }
  loci <- names(genotypes)
  if (anyDuplicated(loci)) {
    stop("'genotypes' names the same locus column twice.", call. = FALSE)
  }
  for (lc in loci) {
    g <- genotypes[[lc]]
    if (!is.numeric(g) || anyNA(g) || any(g < 0) || any(g != trunc(g))) {
      stop("Column '", lc, "' of 'genotypes' must hold non-negative whole ",
           "dosages of allele 1.", call. = FALSE)
    }
  }
  cc <- stats::setNames(rep(list(rep(NA_integer_, nrow(genotypes))),
                             length(loci)), loci)
  if (!is.null(copy_count)) {
    if (is.null(names(copy_count))) {
      stop("'copy_count' must be named by locus.", call. = FALSE)
    }
    extra <- setdiff(names(copy_count), loci)
    if (length(extra) > 0L) {
      stop("'copy_count' names ", paste0("'", extra, "'", collapse = ", "),
           ", which is not a column of 'genotypes'.", call. = FALSE)
    }
    for (nm in names(copy_count)) {
      cc[[nm]] <- .ad_recycle(copy_count[[nm]], nrow(genotypes),
                              paste0("copy_count$", nm))
    }
  }

  keep <- if (isTRUE(drop_zero)) value != 0 else rep(TRUE, length(value))
  if (!any(keep)) {
    stop("Every value is zero, so there is no term to write.", call. = FALSE)
  }
  # One id per row from the row's own loci, so two surfaces bound together
  # cannot share a term_id unless they are the same row of the same loci.
  ids <- .term_id("geno", rep(list(.sort_c(loci)), nrow(genotypes)),
                  seq_len(nrow(genotypes)))
  lab <- if (is.null(effect_name)) NA_character_ else effect_name
  out <- do.call(rbind, lapply(loci, function(lc) {
    .terms_frame(
      term_id          = ids[keep],
      locus_name       = lc,
      contrast_name    = "indicator",
      copy_count_value = as.integer(cc[[lc]])[keep],
      dosage_value     = as.integer(genotypes[[lc]])[keep],
      genome_value     = value[keep],
      effect_name      = lab)
  }))
  rownames(out) <- NULL
  out
}


#' Build `terms` rows for additive-by-additive pairs
#'
#' @description
#' Expands one coefficient per locus pair into the two-member `additive x
#' additive` term [define_genome_effect_terms()] takes: two rows per pair, both
#' `contrast_name = "additive"`, carrying the same coefficient.
#'
#' * `coding = "functional"` — the pair contributes `e * (g_1 - 1) * (g_2 - 1)`
#'   (centres 0.5), where `g` is the dosage of allele 1. Values are absolute,
#'   so the term has a non-zero mean.
#' * `coding = "cockerham"` — the pair contributes
#'   `e * (g_1 - 2 p_1) * (g_2 - 2 p_2)` (centres `p_1`, `p_2`), which has mean
#'   0 under HWE and linkage equilibrium.
#'
#' Within a pair the two loci are put in C-locale order (each `p` moves with its
#' locus), so `(L10, L44)` and `(L44, L10)` give the same rows. `e` is never
#' rescaled. A pair with `e = 0` is dropped, and every `e` being zero is an
#' error, as [ad_terms()] does for `a` and `d`. Every input is validated
#' before zero pairs are dropped.
#'
#' @section Every builder returns the same columns:
#' [ad_terms()], `aa_terms()` and [genotype_terms()] return exactly
#' `term_id`, `locus_name`, `contrast_name`, `center_value`,
#' `copy_count_value`, `dosage_value`, `genome_value`, `effect_name`, in that
#' order and with those types (character, character, character, double,
#' integer, integer, double, character), with a typed `NA` where a column does
#' not apply. So their outputs `rbind()` in any combination. A `term_id`
#' encodes the builder and the term's loci unambiguously, whatever characters
#' the locus names contain, so bound outputs never share a `term_id` unless
#' they describe the same term on the same loci. (The writer may still refuse
#' overlapping definitions on the same loci; that is its family rule, not an
#' id collision.)
#'
#' @section Functional `a` is not Cockerham `alpha` once pairs exist:
#' `ad_terms(coding = "cockerham")` takes the average effect `alpha`, not the
#' functional `a`. Already without pairs `alpha = a + (q - p) d`; with pairs it
#' also needs `sum_l e_jl (2 p_l - 1)`, which `ad_terms()` never sees. Feeding
#' functional `a` to the Cockerham branch therefore changes the **genotypic
#' values**, not just how the variance is split between components. For a
#' hand-written model with pairs, use functional coding for both builders.
#'
#' @section The reported mean is not written anywhere:
#' Under functional coding a `message()` reports each pair's share of the
#' implied genetic mean, `e (2 p_1 - 1)(2 p_2 - 1)` (its expectation under HWE
#' and linkage equilibrium), and their running total. It is written to no
#' table, for the reason given in [ad_terms()].
#'
#' @param locus_name_1,locus_name_2 Character vectors of the same length, one
#'   pair per position. A pair's two loci must differ.
#' @param e Numeric coefficient per pair, recycled.
#' @param p_1,p_2 Allele-1 frequency of each pair's first and second locus,
#'   recycled. Required under both codings: the centres under Cockerham, the
#'   basis of the reported mean under functional.
#' @param coding `"functional"` (default) or `"cockerham"`.
#' @param effect_name Optional label carried onto every term.
#' @param report Logical; print the implied genetic mean (functional coding).
#'   Default `TRUE`.
#'
#' @return A `terms` data frame with the fixed builder column set: two rows per
#'   pair with a non-zero `e`.
#'
#' @seealso [define_genome_effect_terms()], [ad_terms()], [genotype_terms()].
#'
#' @examples
#' \dontrun{
#' pop |>
#'   define_genome_effect_terms(
#'     trait_name = "ADG",
#'     terms = rbind(
#'       ad_terms(locus_name = c("Locus_10", "Locus_44"),
#'                a = c(0.30, -0.12), d = c(0.10, 0.05), p = c(0.35, 0.60)),
#'       aa_terms(locus_name_1 = "Locus_10", locus_name_2 = "Locus_44",
#'                e = 0.08, p_1 = 0.35, p_2 = 0.60)))
#' }
#' @export
aa_terms <- function(locus_name_1, locus_name_2, e, p_1, p_2,
                     coding      = c("functional", "cockerham"),
                     effect_name = NULL,
                     report      = TRUE) {
  coding <- match.arg(coding)
  for (nm in c("locus_name_1", "locus_name_2")) {
    x <- get(nm)
    if (!is.character(x) || length(x) == 0L || anyNA(x) || any(!nzchar(x))) {
      stop("'", nm, "' must be a non-empty character vector of locus names.",
           call. = FALSE)
    }
  }
  n <- length(locus_name_1)
  if (length(locus_name_2) != n) {
    stop("'locus_name_1' and 'locus_name_2' must have the same length (got ",
         n, " and ", length(locus_name_2), ").", call. = FALSE)
  }
  if (missing(p_1) || missing(p_2)) {
    stop("'p_1' and 'p_2' (the allele-1 frequency at each pair's loci) are ",
         "required: they are the centres under Cockerham coding and the basis ",
         "of the reported mean under functional coding.", call. = FALSE)
  }
  e   <- .ad_recycle(e,   n, "e",   of = "locus_name_1")
  p_1 <- .ad_recycle(p_1, n, "p_1", of = "locus_name_1")
  p_2 <- .ad_recycle(p_2, n, "p_2", of = "locus_name_1")
  if (any(c(p_1, p_2) < 0 | c(p_1, p_2) > 1)) {
    stop("'p_1' and 'p_2' must be between 0 and 1.", call. = FALSE)
  }
  lab <- if (is.null(effect_name)) rep(NA_character_, n) else
    .ad_recycle_chr(effect_name, n, "effect_name", of = "locus_name_1")

  same <- locus_name_1 == locus_name_2
  if (any(same)) {
    stop("A pair names the same locus twice: ",
         paste0("'", unique(locus_name_1[same]), "'", collapse = ", "),
         ". An additive-by-additive term needs two different loci.",
         call. = FALSE)
  }
  # Canonical order within each pair; each p moves with its locus.
  swap <- vapply(seq_len(n), function(i) {
    .sort_c(c(locus_name_1[i], locus_name_2[i]))[1] != locus_name_1[i]
  }, logical(1))
  l1 <- ifelse(swap, locus_name_2, locus_name_1)
  l2 <- ifelse(swap, locus_name_1, locus_name_2)
  q1 <- ifelse(swap, p_2, p_1)
  q2 <- ifelse(swap, p_1, p_2)
  ids <- .term_id("aa", lapply(seq_len(n), function(i) c(l1[i], l2[i])))
  if (anyDuplicated(ids)) {
    dup <- which(duplicated(ids))
    stop("Pair(s) repeated: ",
         paste0("(", l1[dup], ", ", l2[dup], ")", collapse = ", "),
         ". Give each pair once, with its total coefficient.", call. = FALSE)
  }

  keep <- e != 0
  if (!any(keep)) {
    stop("Every e is zero, so there is no term to write.", call. = FALSE)
  }
  c1 <- if (coding == "functional") rep(0.5, n) else q1
  c2 <- if (coding == "functional") rep(0.5, n) else q2
  k  <- which(keep)
  out <- .terms_frame(
    term_id       = rep(ids[k], each = 2L),
    locus_name    = as.vector(rbind(l1[k], l2[k])),
    contrast_name = "additive",
    center_value  = as.vector(rbind(c1[k], c2[k])),
    genome_value  = rep(e[k], each = 2L),
    effect_name   = rep(lab[k], each = 2L))

  if (isTRUE(report) && coding == "functional") {
    mu <- e[k] * (2 * q1[k] - 1) * (2 * q2[k] - 1)
    message("Implied genetic mean (functional coding), reported only — ",
            "written to no table: ",
            paste0("(", l1[k], ", ", l2[k], ") mu = ", format(round(mu, 6)),
                   collapse = "; "),
            ". Running total over these pairs: ",
            format(round(sum(mu), 6)), ".")
  }
  rownames(out) <- NULL
  out
}


# -- Shared builder plumbing --------------------------------------------------

#' Sort strings in C-locale byte order, the same on every platform
#'
#' @keywords internal
#' @noRd
.sort_c <- function(x) x[order(x, method = "radix")]

#' Deterministic, collision-free `term_id`s for the term builders
#'
#' `<builder>:<len>:<name>|<len>:<name>...#<suffix>`. Length-prefixing each
#' locus name makes the locus list decode uniquely whatever characters the
#' names contain (`c("A", "B")` is `1:A|1:B`, the single locus `"AxB"` is
#' `3:AxB`), and the builder prefix keeps builders apart. Pure: no counters,
#' no randomness, so a builder's output is a function of its input alone.
#'
#' @param builder `"ad"`, `"aa"` or `"geno"`.
#' @param loci List of character vectors, one per term, in the order to encode.
#' @param suffix Optional vector recycled to `length(loci)`.
#' @keywords internal
#' @noRd
.term_id <- function(builder, loci, suffix = NULL) {
  enc <- vapply(loci, function(l) {
    paste0(nchar(l, type = "chars"), ":", l, collapse = "|")
  }, character(1))
  out <- paste0(builder, ":", enc)
  if (!is.null(suffix)) out <- paste0(out, "#", suffix)
  out
}

#' The one `terms` frame every builder returns
#'
#' Fixed columns, order and types, with a typed `NA` where a column does not
#' apply, so builder outputs `rbind()` in any combination (gate C16).
#'
#' @keywords internal
#' @noRd
.terms_frame <- function(term_id, locus_name, contrast_name,
                         center_value = NA_real_, copy_count_value = NA_integer_,
                         dosage_value = NA_integer_, genome_value,
                         effect_name = NA_character_) {
  n <- length(locus_name)
  rec <- function(x) if (length(x) == 1L) rep(x, n) else x
  data.frame(
    term_id          = as.character(rec(term_id)),
    locus_name       = as.character(locus_name),
    contrast_name    = as.character(rec(contrast_name)),
    center_value     = as.numeric(rec(center_value)),
    copy_count_value = as.integer(rec(copy_count_value)),
    dosage_value     = as.integer(rec(dosage_value)),
    genome_value     = as.numeric(rec(genome_value)),
    effect_name      = as.character(rec(effect_name)),
    stringsAsFactors = FALSE)
}


# -- The NOIA conversion pair (plan Q13) ---------------------------------------
#
# Functional coding: a (g - 1) + d 1[g = 1] + sum e_kl (g_k - 1)(g_l - 1).
# Statistical (NOIA) coding at frequencies p: alpha (g - 2p) + d x_D(p) +
# sum e_kl (g_k - 2p_k)(g_l - 2p_l), with x_D the Cockerham dominance contrast
# (-2p^2, 2pq, -2q^2) the evaluator uses. The two describe the same genotypic
# values up to a constant. Both functions are pure: no database access.

#' Stored common-scope terms -> functional (a, d, e)
#'
#' Converts one trait's covered terms, as `.gev_read_model()` returns them,
#' to functional coefficients keyed by `locus_id`. Sign convention: **stored
#' term = functional term + kappa**, summed over terms.
#'
#' | Stored term | Functional | kappa |
#' |---|---|---|
#' | `additive` v, centre c | a += v | v (1 - 2c) |
#' | `dominance` v, centre c | d += v; a -= v (1 - 2c) | -v (c^2 + (1 - c)^2) |
#' | `indicator (2, 1)` v | d += v | 0 |
#' | `indicator (2, 2)` v | a += v/2; d -= v/2 | v/2 |
#' | `indicator (2, 0)` v | a -= v/2; d -= v/2 | v/2 |
#' | `additive x additive` v, c_k, c_l | e += v; a_k += v (1 - 2c_l); a_l += v (1 - 2c_k) | v (1 - 2c_k)(1 - 2c_l) |
#'
#' @param terms,members The `terms` and `members` of a `.gev_read_model()`
#'   result, restricted to covered terms (every term one of the shapes above).
#' @return `list(a, d, pairs, kappa)`: `a`, `d` named numeric vectors (names
#'   are `locus_id`), every locus a term touches present in both; `pairs` a
#'   data frame `locus_id_1 < locus_id_2, e`, one row per canonical pair;
#'   `kappa` one number.
#' @keywords internal
#' @noRd
.stored_to_functional <- function(terms, members) {
  loci <- as.character(sort(unique(members$locus_id)))
  a <- d <- stats::setNames(numeric(length(loci)), loci)
  kappa <- 0
  pk <- character(0); pe <- numeric(0)
  mem_by <- split(members, members$id_genome_effect)
  for (i in seq_len(nrow(terms))) {
    v  <- terms$genome_value[i]
    mm <- mem_by[[as.character(terms$id_genome_effect[i])]]
    mm <- mm[order(mm$locus_id), , drop = FALSE]
    if (nrow(mm) == 1L) {
      l <- as.character(mm$locus_id)
      c <- mm$center_value
      switch(mm$contrast_name,
        additive = {
          a[l] <- a[l] + v
          kappa <- kappa + v * (1 - 2 * c)
        },
        dominance = {
          d[l] <- d[l] + v
          a[l] <- a[l] - v * (1 - 2 * c)
          kappa <- kappa - v * (c^2 + (1 - c)^2)
        },
        indicator = {
          if (!identical(as.integer(mm$copy_count_value), 2L)) {
            stop("Internal error: .stored_to_functional() got an indicator ",
                 "that is not a diploid state.", call. = FALSE)
          }
          switch(as.character(mm$dosage_value),
            "1" = { d[l] <- d[l] + v },
            "2" = { a[l] <- a[l] + v / 2; d[l] <- d[l] - v / 2
                    kappa <- kappa + v / 2 },
            "0" = { a[l] <- a[l] - v / 2; d[l] <- d[l] - v / 2
                    kappa <- kappa + v / 2 })
        })
    } else if (nrow(mm) == 2L && all(mm$contrast_name == "additive")) {
      lk <- as.character(mm$locus_id[1]); ll <- as.character(mm$locus_id[2])
      ck <- mm$center_value[1];           cl <- mm$center_value[2]
      a[lk] <- a[lk] + v * (1 - 2 * cl)
      a[ll] <- a[ll] + v * (1 - 2 * ck)
      kappa <- kappa + v * (1 - 2 * ck) * (1 - 2 * cl)
      pk <- c(pk, paste(lk, ll)); pe <- c(pe, v)
    } else {
      stop("Internal error: .stored_to_functional() got a term outside the ",
           "covered shapes.", call. = FALSE)
    }
  }
  pairs <- data.frame(locus_id_1 = integer(0), locus_id_2 = integer(0),
                      e = numeric(0))
  if (length(pk) > 0L) {
    # Owners (and duplicate variants) sum into one coefficient per pair.
    e  <- tapply(pe, pk, sum)
    ks <- strsplit(names(e), " ", fixed = TRUE)
    pairs <- data.frame(
      locus_id_1 = as.integer(vapply(ks, `[`, "", 1)),
      locus_id_2 = as.integer(vapply(ks, `[`, "", 2)),
      e = as.numeric(e))
    pairs <- pairs[order(pairs$locus_id_1, pairs$locus_id_2), , drop = FALSE]
    rownames(pairs) <- NULL
  }
  list(a = a, d = d, pairs = pairs, kappa = kappa)
}

#' Functional (a, d, e) -> statistical (NOIA) coefficients at frequencies p
#'
#' The forward map of plan §3:
#' `alpha_j = a_j + (q_j - p_j) d_j + sum_l e_jl (2 p_l - 1)`, with `d` and
#' `e` unchanged and every centre at `p`. Sign convention: **functional model
#' = statistical model + mu**, where `mu` is the functional model's mean under
#' HWE and linkage equilibrium at `p`. A round trip through
#' `.stored_to_functional()` therefore returns `kappa = -mu`.
#'
#' Returns coefficient data, not a writer frame; `.noia_terms()` builds the
#' frame.
#'
#' @param a,d Named numeric vectors (same names, any locus key).
#' @param pairs Data frame `locus_1`, `locus_2`, `e` using the same keys.
#' @param p Named numeric vector of allele-1 frequencies covering every key.
#' @return `list(alpha, d, p, pairs, mu)`; `pairs` gains `p_1`, `p_2`.
#' @keywords internal
#' @noRd
.noia_to_stored <- function(a, d, pairs, p) {
  keys <- names(a)
  if (is.null(keys) || !identical(sort(keys), sort(names(d)))) {
    stop("Internal error: 'a' and 'd' must be named by the same loci.",
         call. = FALSE)
  }
  need <- unique(c(keys, as.character(pairs$locus_1),
                   as.character(pairs$locus_2)))
  if (!all(need %in% names(p)) || anyNA(p[need])) {
    stop("Internal error: 'p' must cover every locus.", call. = FALSE)
  }
  d <- d[keys]
  pk <- p[keys]
  alpha <- a + (1 - 2 * pk) * d
  mu <- sum(a * (2 * pk - 1) + 2 * pk * (1 - pk) * d)
  if (nrow(pairs) > 0L) {
    p1 <- unname(p[as.character(pairs$locus_1)])
    p2 <- unname(p[as.character(pairs$locus_2)])
    for (r in seq_len(nrow(pairs))) {
      k <- as.character(pairs$locus_1[r]); l <- as.character(pairs$locus_2[r])
      alpha[k] <- alpha[k] + pairs$e[r] * (2 * p2[r] - 1)
      alpha[l] <- alpha[l] + pairs$e[r] * (2 * p1[r] - 1)
    }
    mu <- mu + sum(pairs$e * (2 * p1 - 1) * (2 * p2 - 1))
    pairs$p_1 <- p1
    pairs$p_2 <- p2
  }
  list(alpha = alpha, d = d, p = pk, pairs = pairs, mu = mu)
}

#' A `.noia_to_stored()` result as one writer-ready `terms` frame
#'
#' Cockerham [ad_terms()] for `alpha` / `d` plus Cockerham [aa_terms()] for the
#' pairs, keyed by locus name. Zero coefficients are dropped by the builders.
#'
#' @param stat A `.noia_to_stored()` result whose keys are locus names.
#' @keywords internal
#' @noRd
.noia_terms <- function(stat) {
  out <- list()
  keep <- stat$alpha != 0 | stat$d != 0
  if (any(keep)) {
    out[[1]] <- ad_terms(names(stat$alpha)[keep], a = unname(stat$alpha[keep]),
                         d = unname(stat$d[keep]), p = unname(stat$p[keep]),
                         coding = "cockerham", report = FALSE)
  }
  pr <- stat$pairs
  if (nrow(pr) > 0L && any(pr$e != 0)) {
    out[[length(out) + 1L]] <- aa_terms(
      as.character(pr$locus_1), as.character(pr$locus_2), e = pr$e,
      p_1 = pr$p_1, p_2 = pr$p_2, coding = "cockerham", report = FALSE)
  }
  if (length(out) == 0L) {
    stop("Internal error: every coefficient is zero.", call. = FALSE)
  }
  do.call(rbind, out)
}
