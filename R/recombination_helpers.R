#' Pre-compute per-chromosome locus metadata for gamete simulation
#'
#' Called once per `add_offspring()` invocation. The returned list is passed to
#' every `make_gamete()` call, avoiding O(n_loci × n_chr) masking work inside
#' the per-offspring loop.
#'
#' @param locus_map Data frame with columns `chr` and `pos_cM` — the resolved
#'   genetic map (from `resolve_genome_map()`), one row per locus, ordered by
#'   `locus_id`.
#' @return Named list, one element per chromosome. Each element is a list with
#'   `locus_idx` (integer indices into the locus dimension), `pos_cM`
#'   (genetic positions in centiMorgans), and `chr_len` (derived as
#'   `MAX(pos_cM)` for the chromosome; crossovers beyond the last locus are
#'   invisible to inheritance).
#' @keywords internal
build_chr_info <- function(locus_map) {
  chrs <- sort(unique(locus_map$chr))
  stats::setNames(
    lapply(chrs, function(chr_id) {
      mask <- locus_map$chr == chr_id
      list(
        locus_idx = which(mask),
        pos_cM    = locus_map$pos_cM[mask],
        chr_len   = max(locus_map$pos_cM[mask])
      )
    }),
    as.character(chrs)
  )
}

#' Produce a gamete from a parent's two haplotypes via chromosomal crossovers
#'
#' Simulates chromosomal recombination using the Haldane map function.
#' Crossover count per chromosome ~ Poisson(chr_len_cM / 100), i.e.
#' Poisson(genetic length in Morgans) since 1 Morgan = 100 cM. Crossover
#' positions are uniform in genetic (cM) distance within each chromosome.
#'
#' @param hap_matrix 2 x n_loci integer matrix. Row 1 = haplotype from
#'   parent_origin 1, row 2 = haplotype from parent_origin 2.
#' @param chr_info List returned by `build_chr_info()`. Contains pre-computed
#'   locus indices, positions, and lengths per chromosome.
#' @return A list with two length-`n_loci` integer vectors:
#'   \describe{
#'     \item{`allele`}{the gamete alleles (0 or 1).}
#'     \item{`homolog`}{which of the parent's two homologs (1 or 2) donated the
#'       allele at each locus. Used by [add_offspring()] to inherit `line_origin`
#'       from the contributing parental segment.}
#'   }
#'   The random-draw sequence (`rpois`, `sample`, `runif`) is identical to the
#'   allele-only version — `homolog` is derived from values already computed, so
#'   seeded output is unchanged.
#' @keywords internal
make_gamete <- function(hap_matrix, chr_info) {

  n_loci  <- ncol(hap_matrix)
  gamete  <- integer(n_loci)
  homolog <- integer(n_loci)

  for (ci in chr_info) {
    chr_locus_idx <- ci$locus_idx
    chr_pos       <- ci$pos_cM
    chr_len       <- ci$chr_len

    n_cross     <- rpois(1L, lambda = chr_len / 100)
    current_hap <- sample(1L:2L, 1L)

    if (n_cross == 0L) {
      # Whole chromosome comes from one homolog; expand the scalar to a vector.
      hap_idx_vec <- rep(current_hap, length(chr_locus_idx))
    } else {
      cross_pos   <- sort(runif(n_cross, min = 0, max = chr_len))
      # findInterval counts crossovers strictly before each locus position;
      # parity of that count determines how many times the haplotype toggled.
      n_toggles   <- findInterval(chr_pos, cross_pos)
      hap_idx_vec <- (current_hap - 1L + n_toggles %% 2L) %% 2L + 1L
    }

    gamete[chr_locus_idx]  <- hap_matrix[cbind(hap_idx_vec, chr_locus_idx)]
    homolog[chr_locus_idx] <- hap_idx_vec
  }

  list(allele = gamete, homolog = homolog)
}

#' Generate autosome gametes for a batch of offspring (the swappable seam)
#'
#' Stage-0 extraction of the per-offspring autosome gamete generation out of
#' [add_offspring()]. This is the **single swappable seam** that the later
#' dqrng R reference (Stage 2) and C++ kernel (Stage 3) drop into; its signature
#' is deliberately dependency-free — no `tidybreed_pop`, no DBI connection, no
#' `data.frame` — so it is drop-in swappable and unit-testable in isolation.
#'
#' In Stage 0 it keeps the base-R internals ([make_gamete()], base RNG) and the
#' wide-matrix output shape, so seeded `ind_haplotype` output is unchanged. The
#' long-output / `base_seed` / `store_crossovers` signature arrives in Stage 2.
#'
#' **RNG order:** offspring are processed in the given order, **parent_1 (sire)
#' gamete then parent_2 (dam) gamete** — the identical `make_gamete()` call
#' sequence [add_offspring()] used inline, so the global RNG stream is untouched.
#' When special chromosomes are configured, [add_offspring()] calls this per
#' offspring (batch of one) so each offspring's special-chromosome draws remain
#' interleaved directly after its autosome draws (base-R global-stream order).
#'
#' @param sire_ids,dam_ids Character vectors of parent IDs, length `n` (one per
#'   offspring). `sire_ids[i]`/`dam_ids[i]` name the parent_1/parent_2 of
#'   offspring `i`; both index into `parent_haps`/`parent_lo`.
#' @param parent_haps Named list; `parent_haps[[pid]]` is the parent's
#'   `2 x n_autosome_loci` integer allele matrix (homolog x locus).
#' @param parent_lo Named list; `parent_lo[[pid]]` is the parallel
#'   `2 x n_autosome_loci` character `line_origin` matrix.
#' @param autosome_info_by_parent Named list; `autosome_info_by_parent[[pid]]` is
#'   the parent's resolved autosome `chr_info` (from `build_chr_info()`, LOCAL
#'   indices) for its `(sex, line)` map.
#' @param n_autosome_loci Integer. Number of autosome loci (matrix column count).
#' @return A list of four matrices, each `n x n_autosome_loci`:
#'   `sire_allele`, `dam_allele` (integer alleles) and `sire_lo`, `dam_lo`
#'   (character `line_origin`, gathered from the donating homolog).
#' @keywords internal
make_gametes_batch <- function(sire_ids, dam_ids, parent_haps, parent_lo,
                               autosome_info_by_parent, n_autosome_loci) {
  n       <- length(sire_ids)
  col_idx <- seq_len(n_autosome_loci)

  sire_allele <- matrix(0L, nrow = n, ncol = n_autosome_loci)
  dam_allele  <- matrix(0L, nrow = n, ncol = n_autosome_loci)
  sire_lo     <- matrix(NA_character_, nrow = n, ncol = n_autosome_loci)
  dam_lo      <- matrix(NA_character_, nrow = n, ncol = n_autosome_loci)

  for (i in seq_len(n)) {
    sire <- sire_ids[i]
    dam  <- dam_ids[i]

    # make_gamete() returns allele + the per-locus homolog (1/2) it came from;
    # line_origin follows the donating homolog so it stays correct across
    # F1/F2/backcross generations.
    gs <- make_gamete(parent_haps[[sire]], autosome_info_by_parent[[sire]])
    sire_allele[i, ] <- gs$allele
    sire_lo[i, ]     <- parent_lo[[sire]][cbind(gs$homolog, col_idx)]

    gd <- make_gamete(parent_haps[[dam]], autosome_info_by_parent[[dam]])
    dam_allele[i, ] <- gd$allele
    dam_lo[i, ]     <- parent_lo[[dam]][cbind(gd$homolog, col_idx)]
  }

  list(sire_allele = sire_allele, dam_allele = dam_allele,
       sire_lo = sire_lo, dam_lo = dam_lo)
}

#' Pass a chromosome copy through to a gamete without recombination
#'
#' Used for non-recombining chromosomes (the contributing parent's own-sex
#' `chr_meta.recombines_M`/`recombines_F` is `FALSE` — Y, W, MT, most organelles)
#' and for any chromosome where the contributing parent's own copy count for it is
#' 1 (nothing to recombine against regardless of the `recombines_*` flags).
#' Deliberately not routed through [make_gamete()] —
#' that function's whole contract is "simulate crossovers," and a chromosome
#' that never recombines has no crossover model at all.
#'
#' @param hap_matrix `k x n_chr_loci` integer matrix, `k` in `{1, 2}` — the
#'   contributing parent's own stored copies for this chromosome's loci only.
#' @param lo_matrix `k x n_chr_loci` character matrix of `line_origin`,
#'   parallel to `hap_matrix`.
#' @return A list with `allele` and `line_origin`, each a length-`n_chr_loci`
#'   vector. When `k == 2`, one whole row is chosen uniformly at random (a
#'   single `sample.int(2, 1)` draw) and passed through unchanged — no
#'   crossover, no `homolog` vector needed. When `k == 1`, the only row is
#'   passed through with **no** random draw at all — there is nothing to
#'   choose, so strictly hemizygous, non-recombining inheritance (e.g.
#'   patrilineal Y, matrilineal MT) never consumes RNG state.
#' @keywords internal
pass_through_gamete <- function(hap_matrix, lo_matrix) {
  k <- nrow(hap_matrix)
  chosen <- if (k == 1L) 1L else sample.int(k, size = 1L)
  list(
    allele      = hap_matrix[chosen, ],
    line_origin = lo_matrix[chosen, ]
  )
}
