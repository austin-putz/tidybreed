// dqrng R <-> C++ uniform-parity spike (THROWAWAY, dev-only).
//
// Purpose: prove that an independent C++ dqrng generator, seeded from a
// per-gamete stream key, reproduces the SAME uniform sequence that the dqrng
// R API (dqset.seed + dqrunif) produces from the same key. This is the
// foundational assumption behind the per-gamete-stream design in
// plans/optimize_add_offspring.md (decisions #2/#3). If it fails, that design
// needs rework before the C++ kernel (Stage 3) is built.
//
// NOT part of the package build. Compiled ad hoc via Rcpp::sourceCpp() from
// dev/spikes/dqrng_parity/run.R. dqrng/BH/sitmo are header-only LinkingTo deps
// here; they are NOT added to DESCRIPTION in this increment.

// [[Rcpp::depends(dqrng, BH, sitmo)]]
#include <Rcpp.h>
#include <dqrng_generator.h>   // random_64bit_wrapper, xoroshiro128plusplus (default)
#include <convert_seed.h>      // dqrng::convert_seed<uint64_t>(IntegerVector)
using namespace Rcpp;

// The default dqrng generator (matches R's dqRNGkind("Xoroshiro128++")).
typedef dqrng::random_64bit_wrapper<dqrng::xoroshiro128plusplus> rng_t;

// Scheme A — seed the generator directly from the whole key vector, letting
// dqrng's own convert_seed() fold it to a uint64. Mirrors, on the R side,
//   dqrng::dqset.seed(key)         # key = c(base_seed, o, r, kind)
//   dqrng::dqrunif(n)
// Both paths route the identical integer vector through the identical
// convert_seed(), so the generator seed is bit-identical by construction.
// [[Rcpp::export]]
NumericVector spike_unif_seedvec(IntegerVector key, int n) {
  uint64_t s = dqrng::convert_seed<uint64_t>(key);
  rng_t rng;
  rng.seed(s);
  NumericVector out(n);
  for (int i = 0; i < n; ++i) out[i] = rng.uniform01();
  return out;
}

// Scheme B — seed + stream. Mirrors, on the R side,
//   dqrng::dqset.seed(seed, stream)   # stream = c(o, r, kind)
//   dqrng::dqrunif(n)
// For xoroshiro128plusplus, seed(s, st) does gen.seed(s) then gen.jump(st).
// [[Rcpp::export]]
NumericVector spike_unif_seed_stream(IntegerVector seed, IntegerVector stream, int n) {
  uint64_t s  = dqrng::convert_seed<uint64_t>(seed);
  uint64_t st = dqrng::convert_seed<uint64_t>(stream);
  rng_t rng;
  rng.seed(s, st);
  NumericVector out(n);
  for (int i = 0; i < n; ++i) out[i] = rng.uniform01();
  return out;
}

// Raw uint64 -> [0,1) transform check: expose uniform01() alongside the raw
// draw so we can confirm dqrunif uses (x >> 11) * 2^-53 (the generator's
// uniform01), not the boost uniform_real path, for the [0,1) case.
// [[Rcpp::export]]
NumericVector spike_unif_from_seedvec_first(IntegerVector key, int n) {
  return spike_unif_seedvec(key, n);
}

// ---------------------------------------------------------------------------
// Stage-3 pre-port additions: prove the TWO transforms the recombination
// kernel layers on top of uniform01() also match R bit-for-bit, since the
// original spike only covered the 1-arg [0,1) uniform.
// ---------------------------------------------------------------------------

// Ranged uniform via L * uniform01() (the crossover-position transform). If
// this matches dqrng::dqrunif(n, min = 0, max = L), the C++ kernel can draw
// positions as chr_len * uniform01() and stay bit-identical to the R reference.
// [[Rcpp::export]]
NumericVector spike_ranged(IntegerVector key, int n, double L) {
  uint64_t s = dqrng::convert_seed<uint64_t>(key);
  rng_t rng;
  rng.seed(s);
  NumericVector out(n);
  for (int i = 0; i < n; ++i) out[i] = L * rng.uniform01();
  return out;
}

// Log-accumulation inversion Poisson count (the crossover-count sampler), the
// EXACT algorithm .draw_poisson_dqrng() uses in R. Must return an identical
// count for the same stream.
// [[Rcpp::export]]
int spike_poisson(IntegerVector key, double lambda) {
  uint64_t s = dqrng::convert_seed<uint64_t>(key);
  rng_t rng;
  rng.seed(s);
  if (lambda <= 0) return 0;
  double target = -lambda, acc = 0.0;
  int k = 0;
  for (;;) {
    acc += std::log(rng.uniform01());
    if (acc <= target) break;
    ++k;
  }
  return k;
}
