#' Pre-compute per-chromosome locus metadata for gamete simulation
#'
#' Called once per `add_offspring()` invocation. The returned list is passed to
#' every `make_gamete()` call, avoiding O(n_loci × n_chr) masking work inside
#' the per-offspring loop.
#'
#' @param genome_meta_df Data frame with columns `chr` and `pos_Mb`.
#' @return Named list, one element per chromosome. Each element is a list with
#'   `locus_idx` (integer indices into the locus dimension), `pos_Mb`
#'   (positions), and `chr_len` (derived as `MAX(pos_Mb)` for the chromosome;
#'   crossovers beyond the last locus are invisible to inheritance).
#' @keywords internal
build_chr_info <- function(genome_meta_df) {
  chrs <- sort(unique(genome_meta_df$chr))
  stats::setNames(
    lapply(chrs, function(chr_id) {
      mask <- genome_meta_df$chr == chr_id
      list(
        locus_idx = which(mask),
        pos_Mb    = genome_meta_df$pos_Mb[mask],
        chr_len   = max(genome_meta_df$pos_Mb[mask])
      )
    }),
    as.character(chrs)
  )
}

#' Produce a gamete from a parent's two haplotypes via chromosomal crossovers
#'
#' Simulates chromosomal recombination using the Haldane map function.
#' Crossover count per chromosome ~ Poisson(chr_len_Mb / 100), assuming
#' approximately 1 Morgan per 100 Mb. Crossover positions are uniform within
#' each chromosome.
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
    chr_pos       <- ci$pos_Mb
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

#' Pass a chromosome copy through to a gamete without recombination
#'
#' Used for `chr_meta.recombines = FALSE` chromosomes (Y, W, MT, most
#' organelles) and for any chromosome where the contributing parent's own
#' copy count for it is 1 (nothing to recombine against regardless of the
#' `recombines` flag). Deliberately not routed through [make_gamete()] —
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
