# Stage 2 — Line-Origin TBV (COMPLETE)

**Status**: Done (tidybreed 0.48.0, 2026-07-02)
**Parent plan**: [refactor_haplotype.md](refactor_haplotype.md) (Stage 2 of 5)
**Also closes out**: the Stage 1 benchmark exit criterion deferred in
[refactor_haplotype_stage_1.md](refactor_haplotype_stage_1.md) (long-format scale validation,
since wide-vs-long comparison is moot post-greenfield-rewrite).

---

## Outcome

- **Full test suite: 0 failures, 0 errors, 1024 passing** (up from ~1005 at Stage 1), including
  15 new tests in `tests/testthat/test-add_tbv.R`.
- **Parity preserved**: `tests/testthat/test-parity.R` still matches its Stage 1 golden files
  within `tolerance = 1e-8` — no re-baseline needed, since the parity simulation only defines
  population-wide effects and the new SQL is algebraically identical to the old matrix math in
  that case.
- **Benchmark script** added and smoke-tested at small scale (500 ind × 5,000 loci); large-scale
  (5,000 × 50,000) run is opt-in (`TIDYBREED_BENCH_LARGE=1`) and takes considerably longer — not
  run automatically.

---

## What shipped

### `add_tbv()` rewrite (`R/add_tbv.R`)

Replaced the per-trait matrix-multiply computation (which only read population-wide effects,
`genome_effects.line_name IS NULL`) with a centered SQL query that joins `ind_haplotype` to
`genome_effects` on `(locus_name, line_origin)`, with a `NOT EXISTS`-based per-locus fallback to
the population-wide row when no line-specific row exists:

```sql
SELECT h.id_ind,
       SUM((h.allele - COALESCE(e.base_allele_freq, 0)) * e.genome_value) AS tbv_value
FROM ind_haplotype h
JOIN genome_effects e
  ON  h.locus_name          = e.locus_name
  AND e.trait_name          = '{t}'
  AND e.genome_effect_type  = 'additive'
  AND ( e.line_name = h.line_origin
        OR (e.line_name IS NULL AND NOT EXISTS (
              SELECT 1 FROM genome_effects e2
              WHERE e2.locus_name = h.locus_name AND e2.trait_name = '{t}'
                AND e2.genome_effect_type = 'additive' AND e2.line_name = h.line_origin
            )) )
WHERE h.id_ind IN (...) [AND h.parent_origin = {1|2}]   -- imprinting only
GROUP BY h.id_ind
```

- `COALESCE(e.base_allele_freq, 0)` was added on top of the plan's draft SQL — `base_allele_freq`
  is nullable, and without the `COALESCE`, SQL NULL propagation would zero out the *entire* summed
  term (both centering and the raw effect) for that row, not just the centering constant.
- Imprinting (`expressed_parent != "both"`) restricts to one `parent_origin`, composed with the
  same line-specific/fallback logic.
- The "no additive effects" existence check now counts both population-wide and line-specific
  rows (previously only checked `line_name IS NULL`), since line-specific-only traits are valid.
- `index_names`/`type`/`overwrite_index` (true-index computation from already-computed `ind_tbv`
  rows) is untouched.

### Dead code removed

`get_genotype_matrix()` (`R/define_additive_effects.R`) and `get_haplotype_matrix()`
(`R/phenotype_helpers.R`) were each used only by the old `add_tbv()` body. Deleted along with
their generated `man/*.Rd` files.

### Tests (`tests/testthat/test-add_tbv.R`, new file)

15 tests covering every Stage 2 exit-criterion case from the parent plan: F1 crossbred TBV
(hand-verified), F2 recombination (line_origin follows the mosaic genome through two
generations), backcross, `line_origin IS NULL` fallback, partial line-specific coverage (some
loci only), a line-specific row coexisting with a population-wide row for the same locus
(line-specific wins, no double-counting), `base_allele_freq` differing by line, and imprinted +
line-specific interaction. Verification uses an independent R-side re-implementation
(`independent_tbv()` helper in the test file) of the fallback/centering logic, cross-checked
against `add_tbv()`'s SQL output — not a copy of the SQL under test.

### Benchmark (deferred Stage 1 exit criterion)

`dev/benchmarks/benchmark_haplotype_scale.R` — standalone `Rscript`, excluded from
`R CMD check`/`testthat` (`^dev$` added to `.Rbuildignore`). Times `add_founders()`/
`add_offspring()` insert throughput, `add_tbv()` for both population-wide and line-specific
effects, and `extract_genotypes()` PIVOT export, using only existing package functions (no
reinvented setup). Wide-vs-long comparison from the original plan was dropped — the wide tables
no longer exist post-Stage-1 (greenfield rewrite), so that decision is moot; this instead
validates the long format's own throughput at scale.

Small-scale smoke test (500 ind × 5,000 loci, in-memory DB, single run on a dev laptop — not a
tuned or repeated measurement):

| Stage | Rows/individuals | Elapsed | Rate |
|-------|------------------:|--------:|------|
| `add_founders()` | 10,000,000 rows | 23.1s | ~433K rows/sec |
| `add_offspring()` (250 offspring) | 2,500,000 rows | 33.6s | ~74K rows/sec |
| `add_tbv()` population-wide (5,000 QTL) | 1,250 ind | 0.74s | ~1,680 ind/sec |
| `add_tbv()` line-specific w/ fallback (500 QTL) | 1,250 ind | 0.35s | ~3,531 ind/sec |
| `extract_genotypes()` PIVOT (500-locus chip) | 1,250 ind × 500 loci | 0.20s | — |

The large/realistic scale (5,000 × 50,000, `TIDYBREED_BENCH_LARGE=1`) was not run in this
session — extrapolating from the small-scale founder throughput, it involves on the order of 1
billion rows for `add_founders()` alone and would take tens of minutes. Run it manually when
real numbers at that scale are needed.

---

## Deviations from the Stage 2 plan (intentional)

- None functionally — the shipped SQL matches the plan's draft with the one addition
  (`COALESCE` on `base_allele_freq`) noted above and in `R/add_tbv.R`'s inline comment.

---

## What's next (not done in Stage 2)

- **Stage 3** — dosage cache / export hardening beyond current `add_dosage()`.
- **Stage 4** — `define_chr()`, sex chromosomes, organelles, polyploidy.
- **Stage 5** — sparse mutation storage (`add_mutation()`).
- **Large-scale benchmark run** (5,000 × 50,000) — opt-in, not yet executed; see
  `dev/benchmarks/benchmark_haplotype_scale.R`.
