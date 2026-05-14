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
#' @return Integer vector of length n_loci representing the gamete (0 or 1)
#' @keywords internal
make_gamete <- function(hap_matrix, chr_info) {

  n_loci <- ncol(hap_matrix)
  gamete  <- integer(n_loci)

  for (ci in chr_info) {
    chr_locus_idx <- ci$locus_idx
    chr_pos       <- ci$pos_Mb
    chr_len       <- ci$chr_len

    n_cross     <- rpois(1L, lambda = chr_len / 100)
    current_hap <- sample(1L:2L, 1L)

    if (n_cross == 0L) {
      gamete[chr_locus_idx] <- hap_matrix[current_hap, chr_locus_idx]
    } else {
      cross_pos   <- sort(runif(n_cross, min = 0, max = chr_len))
      # findInterval counts crossovers strictly before each locus position;
      # parity of that count determines how many times the haplotype toggled.
      n_toggles   <- findInterval(chr_pos, cross_pos)
      hap_indices <- (current_hap - 1L + n_toggles %% 2L) %% 2L + 1L
      gamete[chr_locus_idx] <- hap_matrix[cbind(hap_indices, chr_locus_idx)]
    }
  }

  gamete
}
