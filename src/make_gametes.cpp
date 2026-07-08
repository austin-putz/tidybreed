// Stage 3b — C++17 recombination kernel (autosome gametes).
//
// A bit-for-bit port of the R reference make_gametes_batch_r()
// (R/recombination_helpers.R) and its per-gamete core (make_gamete /
// .draw_chr_recombination / .draw_poisson_dqrng). Randomness is driven by the
// SAME dqrng engine the R side uses (Scheme A: seed the generator from
// convert_seed<uint64_t>(c(base_seed, sid)), draw uniform01()), so a given base
// seed produces identical output in R and C++. The pre-port spike
// (dev/spikes/dqrng_parity) proved the three primitives used here match the R
// dqrng API bit-for-bit: uniform01() == dqrunif(n), L*uniform01() ==
// dqrunif(n,0,L), and the log-accumulation Poisson count.
//
// Gamete-flat / one-map-per-call: the caller groups gametes by resolved-map key
// and calls this once per group (see add_offspring()). Autosome-only; special
// chromosomes stay on the R path.

#include <Rcpp.h>
#include <dqrng_generator.h>   // random_64bit_wrapper, xoroshiro128plusplus
#include <convert_seed.h>      // dqrng::convert_seed<uint64_t>(IntegerVector)
#include <algorithm>
#include <cmath>
using namespace Rcpp;

// The default dqrng generator (matches R's dqRNGkind("Xoroshiro128++")).
typedef dqrng::random_64bit_wrapper<dqrng::xoroshiro128plusplus> rng_t;

// Fold (o, parent_origin, kind = 1 autosome) to a stream id, identical to R's
// .gamete_stream_id(o, r, 1L) = ((o*2 + (r-1))*2 + 0).
static inline int gamete_stream_id(int o, int r) {
  return (o * 2 + (r - 1)) * 2;
}

// Log-accumulation inversion Poisson count on the current stream — the EXACT
// algorithm of .draw_poisson_dqrng(). lambda <= 0 draws nothing and returns 0;
// lambda > 30 errors (documented ceiling), matching the R guard.
static inline int draw_poisson(rng_t& rng, double lambda) {
  if (lambda <= 0.0) return 0;
  if (lambda > 30.0) {
    Rcpp::stop("Recombination rate lambda = %s exceeds the supported ceiling "
               "(30 Morgans / 3000 cM per chromosome).", lambda);
  }
  double target = -lambda, acc = 0.0;
  int k = 0;
  for (;;) {
    acc += std::log(rng.uniform01());
    if (acc <= target) break;
    ++k;
  }
  return k;
}

//' Autosome gamete kernel (C++), one resolved-map group
//'
//' Bit-for-bit port of [make_gametes_batch_r()]. Same arguments, same list
//' return; see the R reference for the contract. Not called directly — dispatched
//' to by the internal selector `make_gametes_batch()`.
//'
//' @param parent_allele Integer matrix (2*n_parents) x n_autosome_loci.
//' @param parent_lo_code Integer matrix, same shape: line_origin codes.
//' @param gamete_parent_idx Integer length G: 1-based packed parent row.
//' @param gamete_o Integer length G: global offspring index per gamete.
//' @param gamete_origin Integer length G: parent_origin (1/2) per gamete.
//' @param chr_start,chr_end Integer per-chromosome 1-based inclusive ranges.
//' @param chr_pos_cM Numeric per-locus positions (cM), length n_autosome_loci.
//' @param chr_len_cM Numeric per-chromosome length (cM).
//' @param base_seed Integer base seed.
//' @param store_crossovers Logical.
//' @return Named list matching make_gametes_batch_r().
//' @keywords internal
// [[Rcpp::export]]
List make_gametes_batch_cpp(IntegerMatrix parent_allele,
                            IntegerMatrix parent_lo_code,
                            IntegerVector gamete_parent_idx,
                            IntegerVector gamete_o,
                            IntegerVector gamete_origin,
                            IntegerVector chr_start,
                            IntegerVector chr_end,
                            NumericVector chr_pos_cM,
                            NumericVector chr_len_cM,
                            int base_seed,
                            bool store_crossovers) {
  const int G      = gamete_parent_idx.size();
  const int n_loci = parent_allele.ncol();
  const int n_chr  = chr_start.size();

  // Contiguity contract (guards future out-of-order locus_id additions).
  for (int c = 1; c < n_chr; ++c) {
    if (chr_start[c] != chr_end[c - 1] + 1) {
      Rcpp::stop("make_gametes_batch_cpp(): chromosomes must be contiguous "
                 "locus-index blocks.");
    }
  }

  const R_xlen_t total = (R_xlen_t)G * n_loci;
  IntegerVector parent_origin(total);
  IntegerVector locus_idx(total);
  IntegerVector allele(total);
  IntegerVector lo_code(total);

  // Ragged crossover buffers (grown only when store_crossovers).
  std::vector<int>    xo_o, xo_r, xo_chr;
  std::vector<double> xo_pos;

  std::vector<double> cross_pos;   // reused per chromosome
  R_xlen_t pos = 0;

  for (int i = 0; i < G; ++i) {
    const int o    = gamete_o[i];
    const int r    = gamete_origin[i];
    const int pidx = gamete_parent_idx[i];        // 1-based packed parent row
    const int row0 = 2 * (pidx - 1);              // 0-based homolog-1 row
    const int row1 = row0 + 1;                    // homolog-2 row

    rng_t rng;
    rng.seed(dqrng::convert_seed<uint64_t>(
      IntegerVector::create(base_seed, gamete_stream_id(o, r))));

    for (int c = 0; c < n_chr; ++c) {
      const int lo = chr_start[c] - 1;            // 0-based first locus
      const int hi = chr_end[c] - 1;              // 0-based last locus
      const double len = chr_len_cM[c];

      // (1) crossover count, (2) start homolog, (3) sorted positions — the fixed
      // canonical draw order; identical to .draw_chr_recombination().
      const int n_cross = draw_poisson(rng, len / 100.0);
      const int start   = (rng.uniform01() < 0.5) ? 1 : 2;

      cross_pos.clear();
      if (n_cross > 0) {
        cross_pos.resize(n_cross);
        for (int k = 0; k < n_cross; ++k) cross_pos[k] = len * rng.uniform01();
        std::sort(cross_pos.begin(), cross_pos.end());
      }

      for (int L = lo; L <= hi; ++L) {
        int hap_idx;
        if (n_cross == 0) {
          hap_idx = start;
        } else {
          // findInterval: count of crossovers <= this locus's position.
          const int n_tog = (int)(std::upper_bound(cross_pos.begin(),
                                                    cross_pos.end(),
                                                    chr_pos_cM[L]) - cross_pos.begin());
          hap_idx = ((start - 1 + (n_tog % 2)) % 2) + 1;   // {1,2}
        }
        const int src_row = (hap_idx == 1) ? row0 : row1;
        parent_origin[pos] = r;
        locus_idx[pos]     = L + 1;               // 1-based
        allele[pos]        = parent_allele(src_row, L);
        lo_code[pos]       = parent_lo_code(src_row, L);
        ++pos;
      }

      if (store_crossovers && n_cross > 0) {
        for (int k = 0; k < n_cross; ++k) {
          xo_o.push_back(o);
          xo_r.push_back(r);
          xo_chr.push_back(c + 1);                // 1-based chromosome position
          xo_pos.push_back(cross_pos[k]);
        }
      }
    }
  }

  return List::create(
    _["parent_origin"]       = parent_origin,
    _["locus_idx"]           = locus_idx,
    _["allele"]              = allele,
    _["line_origin_code"]    = lo_code,
    _["xover_gamete_o"]      = wrap(xo_o),
    _["xover_parent_origin"] = wrap(xo_r),
    _["xover_chr"]           = wrap(xo_chr),
    _["xover_pos_cM"]        = wrap(xo_pos));
}
