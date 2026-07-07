# Optimize `add_offspring()` — Increment 1 summary (Spike + Phase 0 + Stage 0)

**Date:** 2026-07-07
**Branch:** `feat/position-coordinates-genome-map`
**Scope executed:** the first increment of
[`plans/optimize_add_offspring.md`](optimize_add_offspring.md) — the dqrng R↔C++
parity **Spike**, **Phase 0** profiling + committed baseline, and **Stage 0**
(output-neutral seam extraction + cheap structural wins). Stages 1–4 were **not**
started; they are re-planned from the Phase-0 numbers below.

**Status: complete and green.** No `DESCRIPTION`/`NEWS` bump and no commit yet
(reserved for when requested). `dqrng`/`Rcpp` were **not** added to `DESCRIPTION`
— the spike is dev-only throwaway.

---

## Headline finding (reshapes the roadmap)

**Gamete generation is only ~2–6% of `add_offspring()` wall time.** The
dqrng/C++/parallel gamete kernel (Stages 2–4) — the design's headline effort —
optimizes a **<5% slice**. The real costs are:

1. the **O(P×N) parent-matrix scan** (now fixed in Stage 0), and
2. the **wide→`UNPIVOT` DB write**, which now dominates everywhere
   (34–70% of `add_offspring()`, ~69% of `add_founders()`, ~96% of the founder
   pool write) and drove ~43 GB of allocation churn for just 200 offspring.

**Recommendation for the regroup:** reorder so **Stage 1 (direct long write +
batching + one-transaction-per-call, shared across all three write paths)** is the
next priority — it targets the now-dominant cost. Stages 2–4 (the C++ kernel)
become a *later, smaller* optimization. Because the parity spike is already done,
Stages 2–3 stay unblocked whenever we want them.

---

## Workstream A — dqrng R↔C++ uniform-parity spike ✅ green

**Files:** `dev/spikes/dqrng_parity/spike.cpp`, `dev/spikes/dqrng_parity/run.R`
(dev-only, throwaway; not part of the package build).

**Goal:** prove an independent C++ dqrng generator reproduces the *bit-identical*
uniform sequence the dqrng R API (`dqset.seed` + `dqrunif`) produces from the same
per-gamete stream key — the assumption Stages 2–3 rest on.

**Result:** the risk is retired, and the derivation is **simpler than the design's
original hash64 proposal**. `dqrng::dqset.seed()` accepts only a scalar/2-int seed
and a scalar stream (a length-≥3 vector errors `"out-of-range seed"`), so the
per-gamete key `(base_seed, o, r, kind)` folds to a single stream id
`sid = ((o*2 + (r-1))*2 + (kind-1))` (`kind`: 1 = autosome, 2 = special;
`r` = `parent_origin` ∈ {1,2}). **No hand-rolled `splitmix64`/`hash64` is needed —
dqrng's own `convert_seed()` does the 64-bit folding**, so feeding the identical
integer vector through it in R and C++ yields a bit-identical generator seed by
construction.

**Adopted (Scheme A):**
- **R reference (Stage 2):**
  `dqrng::dqRNGkind("Xoroshiro128++"); dqrng::dqset.seed(c(base_seed, sid)); dqrng::dqrunif(n)`
- **C++ kernel (Stage 3):** `dqrng::convert_seed<uint64_t>(c(base_seed, sid))` →
  `dqrng::random_64bit_wrapper<dqrng::xoroshiro128plusplus>` → `uniform01()`
  (`= (x >> 11) * 2^-53`, which the spike confirmed `dqrunif` uses for `[0,1)`).
- **Bound:** `sid` fits signed-32-bit for `o < ~5e8` (>> realistic `n_offspring`);
  if ever exceeded, split `sid` into a 2-int stream — `convert_seed` handles it.
- Scheme B (`dqset.seed(base_seed, sid)` scalar stream → C++ `seed(base)+jump(sid)`)
  was **also** proven bit-identical and is the fallback; Scheme A is preferred.

**Spike output (dqrng 0.4.1):** Scheme A and Scheme B both `ALL IDENTICAL` across
6 tuples × {1, 5, 50} draws; distinct keys → distinct streams (`ALL DISTINCT`);
`dqrunif` matches C++ `uniform01()`. Verdict recorded in the design doc's
stream-derivation section and its implementation-order table.

---

## Workstream B — Phase 0 profiler + committed baseline ✅

**File:** `dev/benchmarks/profile_haplotype_paths.R` (new; distinct from the
broad-throughput `benchmark_haplotype_scale.R`). Base-R `Rprof(line.profiling=)`
attributes self-time to source lines, bucketed into named steps; peak memory =
max Vcells since `gc(reset=TRUE)`. `profvis`/`bench` are optional, `requireNamespace`-
guarded, and **not** package `Suggests`. Scenarios: founder creation, litter/hatch,
1-sire × N, and a `define_chr()` sex-chromosome run. The committed baseline table +
conclusion live at the bottom of that script.

**Baseline (small: 5000 loci, 10 chr, 600-hap pool, 250M+250F/line, 500 offspring;
read the % split, not the machine-specific seconds):**

| Path / call | total | dominant steps |
|---|---|---|
| `define_founder_haplotypes()` A | 0.80s | pool long-materialize + write **95.8%** |
| `add_founders()` A | 14.2s | DB-write (UNPIVOT) **68.8%**; read-pool + rebuild-matrix **27.9%** |
| `add_offspring()` litter (250 parents × 500) | 40.1s | **parent-matrix O(P×N) scan 58.6%**; DB-write 34.3%; gamete-gen **1.7%** |
| `add_offspring()` 1-sire × 500 | 18.7s | DB-write 60.5%; O(P×N) scan 29.8%; gamete-gen 3.2% |
| `add_offspring()` + sex chr (100 off) | 1.3s | DB-write 50.4%; O(P×N) scan 35.1%; gamete-gen 1.9% |

`bench`: `add_offspring(200)` `mem_alloc` = **43.2 GB** (dense wide-matrix churn).

`large` (50k loci, 30 chr, 10k offspring) is env-gated (`TIDYBREED_BENCH_LARGE=1`)
and intentionally not run on the unmodified tree — the dense wide matrices (~2 GB
each × 4) are expected to thrash/OOM, itself the Stage-1 motivation. Re-run large
after Stage 1 to confirm the ceiling.

---

## Workstream C — Stage 0 refactor ✅ output-neutral

Pure structural refactor; seeded `ind_haplotype` output unchanged.

- **C1 — seam extraction.** New `make_gametes_batch()` in
  `R/recombination_helpers.R`: the dependency-free (no `pop`/DBI/`data.frame`)
  autosome gamete generator that Stages 2–3 will swap into. Stage 0 keeps base-R
  `make_gamete()` and the wide-matrix output shape. `R/add_offspring.R` now calls it
  as a **true batch** in the autosome-only fast path, and **per-offspring**
  (interleaved with the special-chromosome path) when a sex chromosome is
  configured — preserving base-R global-stream draw order in both cases.
- **C2 — parent scan.** Replaced the per-parent full `data.frame` scan
  (`parent_haps_raw[parent_haps_raw$id_ind == pid, ]`, O(parents × rows)) with one
  `split()` of **row indices** + a single subset per parent. Splitting *indices*
  (not the frame) keeps peak memory low (one parent's subset at a time) while
  removing the repeated scan.
- **C3 — `special_rows`.** Preallocated to its known upper bound and filled via a
  counter, replacing the quadratic `list[[length + 1]]` growth.

### Verification

- **Autosome output-neutrality:** `test-parity.R` golden files (2-line → F1 → F2)
  green — byte-identical haplotypes/dosage/TBV/export.
- **Special-path output-neutrality:** proved `make_gametes_batch()` reproduces the
  inline `make_gamete()` draw sequence exactly for the same seed (single-offspring
  sire-then-dam **and** batch-of-N == N sequential calls). The special path routes
  only through the seam plus a non-RNG accumulator change ⇒ byte-identical.
- **Full suite:** `testthat::test_local()` → **1229 pass / 0 fail / 0 error**
  (9 pre-existing warnings, unrelated).

### Measured Stage-0 wins (litter, 250 parents × 500 offspring)

| metric | before | after |
|---|---|---|
| O(P×N) parent scan | 58.6% (19.1s) | ~10% (1.5s) — **~12×** on that step |
| total time | 40.0s | 12.5s — **~3.2×** |
| `mem_alloc` (200 offspring) | 43.2 GB | 1.55 GB — **~28×** less churn |
| peak Vcells (litter) | 402 MB | 418 MB (index-split; a data.frame-split variant briefly hit 659 MB and was reverted) |

`gamete-generation` stays ~2–6% (unchanged), confirming Stage 0 is a pure refactor
plus the two structural wins — and reconfirming the headline finding that the
gamete kernel is a small slice.

---

## Files touched

**Package code (Stage 0, output-neutral):**
- `R/recombination_helpers.R` — added `make_gametes_batch()` seam.
- `R/add_offspring.R` — call the seam; index-`split()` parent scan; preallocated
  `special_rows`.

**Dev-only (not part of the package build; no `DESCRIPTION` change):**
- `dev/spikes/dqrng_parity/spike.cpp`, `dev/spikes/dqrng_parity/run.R` — the spike.
- `dev/benchmarks/profile_haplotype_paths.R` — Phase 0 profiler + committed baseline.

**Design doc updated:**
- `plans/optimize_add_offspring.md` — stream-derivation section, open-question #1,
  and the implementation-order Spike row marked done with the adopted scheme.

---

## Open items for the regroup

1. **Reprioritize Stage 1 ahead of the C++ kernel** (per the headline finding):
   direct long write + batching + one-transaction-per-call, shared across
   `add_offspring()`, `add_founders()`, and `define_founder_haplotypes()`'s pool
   write. This is where the now-dominant DB-write / churn cost lives.
   - **Memory decisions folded into the plan (2026-07-07):** batch size is
     **exposed** (`batch_size`/`max_batch_mem`) **and defaults to a RAM-aware
     auto-pick** (detect free memory on Windows / macOS Intel+ARM / Ubuntu,
     conservative fallback). Batching bounds peak independent of total offspring, so a
     single high-fecundity mega-mating (aquaculture/poultry) is handled with no
     separate mechanism. A count-based `matings` spec is **deferred but non-foreclosing**
     (keep `matings` expansion at the top of `add_offspring()`). See Stage 1 (1b) +
     the "High-fecundity / low-memory contract" in `optimize_add_offspring.md`.
2. **Run the gated large-scale profile** to confirm the wide-matrix ceiling (and,
   after Stage 1, the improvement).
3. **Commit Stage 0** with the `DESCRIPTION` version + `NEWS.md` entry when ready
   (deferred until requested).
