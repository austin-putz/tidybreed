#' Build `terms` rows for a functional or Cockerham (a, d) pair
#'
#' @description
#' Expands one locus's additive and dominance coefficients into the two member
#' rows [define_genome_effects()] takes. There is no third table and no separate
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
#' @return A `terms` data frame: two rows per locus with a non-zero
#'   coefficient, `term_id` `"<locus>_a"` / `"<locus>_d"`.
#'
#' @seealso [define_genome_effects()], [genotype_terms()].
#'
#' @examples
#' \dontrun{
#' tt <- ad_terms("Locus_10", a = 0.4, d = 0.2, p = 0.3)
#' pop <- pop |> define_genome_effects("ADG", tt, effect_owner = "functional")
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
  add <- data.frame(
    term_id       = paste0(locus_name, "_a"),
    locus_name    = locus_name,
    contrast_name = "additive",
    center_value  = centre,
    genome_value  = a,
    effect_name   = lab,
    stringsAsFactors = FALSE
  )
  dom <- if (coding == "functional") {
    data.frame(term_id = paste0(locus_name, "_d"), locus_name = locus_name,
               contrast_name = "indicator", copy_count_value = 2L,
               dosage_value = 1L, genome_value = d, effect_name = lab,
               stringsAsFactors = FALSE)
  } else {
    data.frame(term_id = paste0(locus_name, "_d"), locus_name = locus_name,
               contrast_name = "dominance", center_value = p,
               genome_value = d, effect_name = lab, stringsAsFactors = FALSE)
  }
  out <- .ad_rbind_fill(add[a != 0, , drop = FALSE], dom[d != 0, , drop = FALSE])
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
.ad_recycle <- function(x, n, what) {
  if (!is.numeric(x) || anyNA(x)) {
    stop("'", what, "' must be numeric and free of NA.", call. = FALSE)
  }
  if (length(x) == 1L) return(rep(as.numeric(x), n))
  if (length(x) != n) {
    stop("'", what, "' must have length 1 or length(locus_name) (", n, ").",
         call. = FALSE)
  }
  as.numeric(x)
}

#' @keywords internal
#' @noRd
.ad_recycle_chr <- function(x, n, what) {
  if (!is.character(x)) stop("'", what, "' must be character.", call. = FALSE)
  if (length(x) == 1L) return(rep(x, n))
  if (length(x) != n) {
    stop("'", what, "' must have length 1 or length(locus_name) (", n, ").",
         call. = FALSE)
  }
  x
}

#' rbind two frames with different column sets, filling the gaps with NA
#'
#' @keywords internal
#' @noRd
.ad_rbind_fill <- function(x, y) {
  cols <- union(names(x), names(y))
  fill <- function(d) {
    for (cl in setdiff(cols, names(d))) d[[cl]] <- rep(NA, nrow(d))
    d[, cols, drop = FALSE]
  }
  rbind(fill(x), fill(y))
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
#' the opposite with `define_genome_effects(require_complete = TRUE)`, which
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
#' @return A `terms` data frame with `nrow(genotypes) * ncol(genotypes)` rows
#'   (before `drop_zero`).
#'
#' @seealso [define_genome_effects()], [ad_terms()].
#'
#' @examples
#' \dontrun{
#' cells <- expand.grid(Locus_10 = 0:2, Locus_44 = 0:2)
#' vals  <- c(0, 0, 0,  0, 1.4, 2.1,  0, 2.1, 3.6)
#' pop <- pop |> define_genome_effects(
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
  ids <- seq_len(nrow(genotypes))
  out <- do.call(rbind, lapply(loci, function(lc) {
    data.frame(
      term_id          = ids[keep],
      locus_name       = lc,
      contrast_name    = "indicator",
      copy_count_value = as.integer(cc[[lc]])[keep],
      dosage_value     = as.integer(genotypes[[lc]])[keep],
      genome_value     = value[keep],
      effect_name      = if (is.null(effect_name)) NA_character_ else effect_name,
      stringsAsFactors = FALSE
    )
  }))
  rownames(out) <- NULL
  out
}
