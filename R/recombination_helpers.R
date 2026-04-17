#' Produce a gamete from a parent's two haplotypes via chromosomal crossovers
#'
#' Simulates chromosomal recombination using the Haldane map function.
#' Crossover count per chromosome ~ Poisson(chr_len_Mb / 100), assuming
#' approximately 1 Morgan per 100 Mb. Crossover positions are uniform within
#' each chromosome.
#'
#' @param hap_matrix 2 x n_loci integer matrix. Row 1 = haplotype from
#'   parent_origin 1, row 2 = haplotype from parent_origin 2.
#' @param genome_meta_df Data frame with columns `chr` (integer) and `pos_Mb`
#'   (double), one row per locus. Row j must correspond to locus_j (sorted by
#'   locus_id ascending).
#' @param chr_len_Mb Numeric vector indexed by chromosome number (1-based).
#'   Element i = length of chromosome i in megabases.
#' @return Integer vector of length n_loci representing the gamete (0 or 1)
#' @keywords internal
make_gamete <- function(hap_matrix, genome_meta_df, chr_len_Mb) {

  n_loci <- ncol(hap_matrix)
  gamete  <- integer(n_loci)
  chrs    <- sort(unique(genome_meta_df$chr))

  for (chr_id in chrs) {
    chr_mask      <- genome_meta_df$chr == chr_id
    chr_locus_idx <- which(chr_mask)               # column indices in gamete
    chr_pos       <- genome_meta_df$pos_Mb[chr_mask]
    chr_len       <- chr_len_Mb[chr_id]

    n_cross     <- rpois(1L, lambda = chr_len / 100)
    current_hap <- sample(1L:2L, 1L)               # random starting haplotype

    if (n_cross == 0L) {
      gamete[chr_locus_idx] <- hap_matrix[current_hap, chr_locus_idx]
    } else {
      cross_pos <- sort(runif(n_cross, min = 0, max = chr_len))
      cross_ptr <- 1L

      for (k in seq_along(chr_locus_idx)) {
        locus_pos <- chr_pos[k]

        # Toggle haplotype for each crossover passed before this locus
        while (cross_ptr <= n_cross && cross_pos[cross_ptr] < locus_pos) {
          current_hap <- 3L - current_hap            # 1 <-> 2
          cross_ptr   <- cross_ptr + 1L
        }

        gamete[chr_locus_idx[k]] <- hap_matrix[current_hap, chr_locus_idx[k]]
      }
    }
  }

  gamete
}
