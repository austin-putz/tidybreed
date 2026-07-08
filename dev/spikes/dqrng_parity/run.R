#!/usr/bin/env Rscript
# dqrng R <-> C++ uniform-parity spike driver (THROWAWAY, dev-only).
#
#   Rscript dev/spikes/dqrng_parity/run.R
#
# Compiles spike.cpp and checks, for several per-gamete stream keys, whether
# the dqrng R API (dqset.seed + dqrunif) and an independent C++ dqrng generator
# produce BIT-IDENTICAL uniform sequences. See spike.cpp for the design.
#
# Verdict (which derivation is bit-parity-safe) is recorded in
# plans/optimize_add_offspring.md's stream-derivation section.

stopifnot(requireNamespace("Rcpp", quietly = TRUE),
          requireNamespace("dqrng", quietly = TRUE))

Rcpp::sourceCpp(file.path("dev", "spikes", "dqrng_parity", "spike.cpp"))

# Match the C++ default generator (dqrng::xoroshiro128plusplus).
dqrng::dqRNGkind("Xoroshiro128++")

# dqset.seed() accepts a scalar or 2-int seed, and a SCALAR stream (probed:
# a length-3+ vector errors "out-of-range seed"). So a per-gamete stream key
# (base_seed, o, r, kind) must fold (o, r, kind) into a single stream id.
# kind = 1L autosome / 2L special; r = parent_origin in {1,2}.
sid <- function(o, r, kind) ((as.integer(o) * 2L + (as.integer(r) - 1L)) * 2L +
                             (as.integer(kind) - 1L))   # unique, fits 32-bit for o < 5e8

keys <- list(
  list(base = 42L,         o = 1L,    r = 1L, kind = 1L),
  list(base = 42L,         o = 1L,    r = 2L, kind = 1L),
  list(base = 42L,         o = 2L,    r = 1L, kind = 1L),
  list(base = 42L,         o = 1L,    r = 1L, kind = 2L),   # special stream, same (o,r)
  list(base = 7L,          o = 9999L, r = 2L, kind = 1L),
  list(base = 2147483647L, o = 500L,  r = 1L, kind = 1L)    # near INT_MAX base seed
)
ns <- c(1L, 5L, 50L)

cat("=== dqrng version:", as.character(packageVersion("dqrng")), "===\n\n")

report <- function(scheme, ok_vec, supported) {
  if (!supported) { cat(sprintf("  %-28s UNSUPPORTED\n", scheme)); return(invisible()) }
  cat(sprintf("  %-28s %s\n", scheme,
              if (all(ok_vec)) "ALL IDENTICAL" else
                sprintf("MISMATCH in %d/%d cases", sum(!ok_vec), length(ok_vec))))
}

# ---- Scheme A: 2-int seed c(base, sid) through convert_seed ------------------
a_ok <- logical(0); a_supported <- TRUE
for (k in keys) {
  seedvec <- c(as.integer(k$base), sid(k$o, k$r, k$kind))
  for (n in ns) {
    r_out <- tryCatch({ dqrng::dqset.seed(seedvec); dqrng::dqrunif(n) },
                      error = function(e) { a_supported <<- FALSE; e })
    if (!a_supported) break
    a_ok <- c(a_ok, identical(r_out, spike_unif_seedvec(seedvec, n)))
  }
  if (!a_supported) break
}
report("Scheme A (2-int seedvec):", a_ok, a_supported)

# ---- Scheme B: seed = base, scalar stream = sid(o,r,kind) --------------------
b_ok <- logical(0); b_supported <- TRUE
for (k in keys) {
  base <- as.integer(k$base); stream <- sid(k$o, k$r, k$kind)
  for (n in ns) {
    r_out <- tryCatch({ dqrng::dqset.seed(base, stream); dqrng::dqrunif(n) },
                      error = function(e) { b_supported <<- FALSE; e })
    if (!b_supported) break
    b_ok <- c(b_ok, identical(r_out, spike_unif_seed_stream(base, stream, n)))
  }
  if (!b_supported) break
}
report("Scheme B (seed + stream):", b_ok, b_supported)

# ---- Stage-3 pre-port: ranged-uniform parity dqrunif(n,0,L) vs L*uniform01() -
# The kernel draws crossover positions as uniforms scaled to (0, chr_len). If
# dqrng's ranged form is min + (max-min)*uniform01() (one uniform01 per value),
# C++ `L*uniform01()` is bit-identical and Stage-2 output is preserved.
ranged_ok <- logical(0)
for (k in keys) {
  seedvec <- c(as.integer(k$base), sid(k$o, k$r, k$kind))
  for (L in c(1, 50, 137.035)) for (n in ns) {
    dqrng::dqset.seed(seedvec)
    r_out <- dqrng::dqrunif(n, min = 0, max = L)
    ranged_ok <- c(ranged_ok, identical(r_out, spike_ranged(seedvec, n, L)))
  }
}
report("ranged uniform (0,L):", ranged_ok, TRUE)

# ---- Stage-3 pre-port: Poisson-count parity (log-accumulation inversion) -----
# The EXACT R sampler .draw_poisson_dqrng() vs the C++ mirror, same stream.
r_poisson <- function(seedvec, lambda) {
  dqrng::dqset.seed(seedvec)
  if (lambda <= 0) return(0L)
  s <- 0; k <- 0L
  repeat { s <- s + log(dqrng::dqrunif(1L)); if (s <= -lambda) break; k <- k + 1L }
  k
}
pois_ok <- logical(0)
for (k in keys) {
  seedvec <- c(as.integer(k$base), sid(k$o, k$r, k$kind))
  for (lam in c(0.3, 1, 5, 25)) {
    pois_ok <- c(pois_ok, identical(as.integer(r_poisson(seedvec, lam)),
                                    as.integer(spike_poisson(seedvec, lam))))
  }
}
report("poisson count (log-accum):", pois_ok, TRUE)

# ---- Independence sanity: distinct keys must give distinct streams ----------
mk <- function(base, o, r, kind) spike_unif_seedvec(c(as.integer(base), sid(o, r, kind)), 5L)
streams <- list(mk(42, 1, 1, 1), mk(42, 1, 2, 1), mk(42, 2, 1, 1), mk(42, 1, 1, 2))
cat(sprintf("  %-28s %s\n", "independence (4 keys):",
            if (!any(duplicated(streams))) "ALL DISTINCT" else "COLLISION (bad)"))

cat("\n=== VERDICT ===\n")
if (a_supported && length(a_ok) && all(a_ok)) {
  cat("Scheme A (2-int seed c(base, sid) via convert_seed) is BIT-PARITY-SAFE.\n",
      "  sid(o,r,kind) = ((o*2 + (r-1))*2 + (kind-1))\n",
      "  Stage 2 R ref: dqrng::dqset.seed(c(base_seed, sid)); dqrng::dqrunif(...)\n",
      "  Stage 3 C++:   dqrng::convert_seed<uint64_t>(c(base_seed, sid)) -> xoroshiro128plusplus.uniform01()\n",
      sep = "")
} else if (b_supported && length(b_ok) && all(b_ok)) {
  cat("Scheme B (seed=base_seed, stream=sid; jump) is BIT-PARITY-SAFE (A unavailable).\n")
} else {
  cat("NEITHER scheme is bit-parity-safe as written — revisit derivation before Stage 3.\n")
}
