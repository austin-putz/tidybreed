# Optimize the haplotype write path — Stage 3 summary (C++17 kernel + R↔C++ parity)

**Date:** 2026-07-07
**Branch:** `feat/position-coordinates-genome-map`
**Scope executed:** Stage 3 of
[`plans/optimize_add_offspring.md`](optimize_add_offspring.md) — port the `add_offspring()`
autosome gamete kernel to **compiled C++17** (`Rcpp` + the `dqrng` C++ headers), driven by
the same `dqrng` engine as the Stage-2 R reference, so a given seed produces
**byte-identical** output in R and C++. Seeded output is **unchanged from 0.52.0**
(within-version R↔C++ parity + forward determinism; no cross-version guarantee, per
CLAUDE.md). Stage 4 (across-individual `RcppParallel`) was **not** started.

**Split into two commits** (user decision): **3a** an output-neutral R-seam refactor to a
C++-ready layout, then **3b** the C++ port. The pure-R kernel is **kept permanently** as the
parity oracle and a no-compiler fallback.

**Status: complete and green.** Full suite **1531 pass / 0 fail / 0 error** (C++ is the
default path). `R CMD check`: **0 errors**, `checking compiled code ... OK`,
`checking tests ... OK` (remaining 4 warnings / 4 notes are pre-existing, unrelated).
Version bumped to **0.53.0**; `NEWS.md` consolidates the Stage 2 (dqrng) and Stage 3 (C++)
work. Commits `cf6e415` (3a) and `49f1f33` (3b).

---

## What changed

### Pre-port micro-spike — `dev/spikes/dqrng_parity/`
The original spike proved only the 1-arg uniform (`uniform01() == dqrunif(n)`). Extended it
to confirm the **two** transforms the kernel adds: the **ranged uniform**
(`dqrunif(n,0,L) == L*uniform01()`) and the **log-accumulation Poisson count** — both
**bit-identical** (dqrng 0.4.1). This de-risked the port with **no R-reference fallback
needed** and let the kernel keep the **cM unit** (not the plan's "Morgans" — a cosmetic ÷100
that would only risk ULP drift).

### Stage 3a — gamete-flat, packed/coded/long-native R seam (output-neutral)
`R/recombination_helpers.R`, `R/add_offspring.R`.
- **`make_gametes_batch_r()`** reshaped to a dependency-free numeric function over **packed
  integer arrays** (`parent_allele` / `parent_lo_code`, homolog-major), `line_origin`
  **integer codes**, and **long parallel output vectors** — no string-keyed lists, no
  character `line_origin`, no wide matrices. The per-gamete core (`make_gamete()`,
  `.draw_chr_recombination()`, `.draw_poisson_dqrng()`) is **unchanged**, so streams and
  alleles are identical.
- **Gamete-flat / one-map-per-call (group-by-map).** A gamete's map is its *producing
  parent's* resolved `(sex, line)` map, so a sire and a dam of the same offspring can use
  **different** maps (sex-specific / achiasmy / line maps). `add_offspring()` groups gametes
  by map key and calls the kernel once per group. This corrected an initial single-map
  collapse that regressed the existing sex-specific-map test.
- **`add_offspring()`** builds packed parents + the `line_origin` dictionary + per-key chr
  arrays once; the batch loop consumes the kernel's long output directly (no wide→long
  reshape). Special-chromosome path and `ind_crossover` writes unchanged.
- **Output-neutral:** byte-identical `ind_haplotype`; parity golden green with **no
  regeneration**; suite back to **1332**.

### Stage 3b — C++17 kernel — `src/make_gametes.cpp`
- **`make_gametes_batch_cpp()`**, a bit-for-bit port of `make_gametes_batch_r()`. Same
  **log-accumulation inversion Poisson sampler** (not `std::poisson_distribution`),
  `chr_len * uniform01()` positions, `std::upper_bound` == R `findInterval` parity, ragged
  crossover buffer. Randomness from the same `dqrng` engine via Scheme A:
  `convert_seed<uint64_t>(c(base_seed, sid))` → `xoroshiro128plusplus.uniform01()`.
  **Autosome-only** — special chromosomes stay on the R path (wrapper-level parity).
- **Runtime selector `make_gametes_batch()`** — prefers the compiled kernel, falls back to
  the R oracle; `getOption("tidybreed.kernel")` / `TIDYBREED_KERNEL` (`"auto"` default,
  `"r"` forces R). `add_offspring()` calls the dispatcher.

### Packaging (first compiled code in the package)
- `DESCRIPTION`: `Imports: Rcpp`, `LinkingTo: Rcpp, dqrng, BH, sitmo`,
  `SystemRequirements: C++17`, `withr` added to `Suggests`. Version **0.53.0**.
- `R/tidybreed-package.R`: `@useDynLib tidybreed, .registration = TRUE` +
  `@importFrom Rcpp sourceCpp`; `NAMESPACE` regenerated. `Rcpp::compileAttributes()` →
  `src/RcppExports.cpp` + `R/RcppExports.R`. No `Makevars` needed (code is C++11-compatible;
  `SystemRequirements: C++17` declared).
- `.gitignore` excludes `src/*.o|so|dll`; `.Rbuildignore` excludes the standalone manual dev
  scripts (`tests/test_*.R`) and stray output dirs.

---

## Tests

- **`tests/testthat/test-make_gametes_parity.R`** (new): R↔C++ kernel parity via
  `expect_identical` across several seeds, chr configs (single / many-short), and both
  `store_crossovers` states; `lambda`-ceiling error parity; **wrapper-level** `add_offspring()`
  parity (C++ vs forced-R) **including a configured sex chromosome** — both `ind_haplotype`
  and `ind_crossover` identical; and the runtime selector.
- **Parity golden** (`parity_golden/*.rds`) unchanged — 3a is output-neutral, so **no
  regeneration**.
- Existing nets stay green: `test-genome_map.R` (sex-specific map), `test-recombination_dqrng.R`
  (R core unchanged), `test-sex_chromosomes.R`, `test-add_offspring.R`, `test-long-writers.R`.

**Full suite:** `testthat::test_local()` → **1531 pass / 0 fail / 0 error** (41 files).
**`R CMD check`:** compiled code OK, tests OK, **0 errors** (4 pre-existing warnings / 4 notes).

---

## Files touched

**Package code:**
- `R/recombination_helpers.R` — gamete-flat packed seam, `.chr_info_from_arrays()`, selector.
- `R/add_offspring.R` — packed inputs, group-by-map dispatch, long-output consumption.
- `R/tidybreed-package.R` — `@useDynLib` / `@importFrom Rcpp sourceCpp`.
- `src/make_gametes.cpp` — **new** C++ kernel; `src/RcppExports.cpp`, `R/RcppExports.R` — generated.
- `DESCRIPTION`, `NAMESPACE`, `NEWS.md`, `.gitignore`, `.Rbuildignore`, `man/*` — packaging/docs.

**Tests / dev:**
- `tests/testthat/test-make_gametes_parity.R` — **new**.
- `dev/spikes/dqrng_parity/{spike.cpp,run.R}` — ranged-uniform + Poisson parity checks.

**New dependencies:** `Rcpp` (Imports); `Rcpp, dqrng, BH, sitmo` (LinkingTo).

### Post-stage correctness audit (2026-07-08)

A follow-up read-through found **no silent correctness bug** — the per-gamete `(sex, line)`
map resolution, the compact-locus contiguity, the `xover_chr` mapping, and the
`line_origin` code round-trip all verified correct. It closed test-coverage gaps and did
light cleanup (all output-neutral):
- **Map coverage** (`test-genome_map.R`): new **line-specific-map** and **cross-line
  distinct-map** tests through `add_offspring()`, plus a dam-side assertion on the
  sex-specific test — the cross-line test proves each gamete resolves its own parent's map
  (line-A sire changes, line-B dam byte-identical). Non-vacuity confirmed by a broken-key
  experiment.
- **Crossover batch-invariance** (`test-recombination_dqrng.R`): `ind_crossover` **events**
  identical across `batch_size` (surrogate `id_crossover` excluded — it is write-order).
- **Cleanup:** corrected two stale comments in `add_offspring.R`; added an int32
  stream-id overflow guard (loud error instead of silent R/C++ divergence at absurd
  offspring counts); C++ error format `%s`→`%g`.

Full suite **1537 pass / 0 fail** on both the C++ and forced-R paths; `R CMD check` 0 errors.

---

## Deferred to Stage 4 (explicitly not in this stage)

Across-individual parallelism (`RcppParallel::parallelFor` over the offspring range, with
`std::thread` as a zero-dependency fallback) and its thread-count-invariance test. The
per-gamete `dqrng` streams already make this reproducible; the single-threaded C++ kernel and
its parity oracle are the foundation it builds on.

## Open items

1. **Stage 4** — parallelize across individuals; thread-count invariance test.
2. **Pre-existing check hygiene** (out of Stage 3 scope): non-ASCII characters in several
   older `R/*.R` files, a `%||%`-style Rd `\name` warning, and undeclared deps in the dev
   vignette/scripts (`tidyr`, `ggplot2`, `patchwork`). None block install or the compiled build.
