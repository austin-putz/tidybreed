# Import QTL-effect methods — Step 5a implementation review

Reviewed 2026-10-09: tidybreed **0.75.3**, commit **`9c08715`**, against its
0.75.2 parent **`0276271`**. This reviews the implemented 5a changes described in
[import_qtl_effect_methods_phase_5.md](import_qtl_effect_methods_phase_5.md),
the 5a requirements and downstream dependencies in
[import_qtl_effect_methods_phase_5_plan.md](import_qtl_effect_methods_phase_5_plan.md),
and the production code and tests. It does not repeat the earlier review of the
unimplemented 5b design.

**Assessment: the efficiency change is effective, and I found no blocking regression
for the planned 5b workloads. Proceed with 5b. There is one low-priority precision
defect in the new grouping helper; its small guard is worth adding before building
further on it.** Both writer performance gates pass in an independent rerun, and
all 12 targeted test files pass unchanged.

## Finding 1 — Low: guard the numeric pair key against double-precision collisions

**Location:** [`R/define_genome_effect_terms.R:440`](../R/define_genome_effect_terms.R#L440),
`.ge_pair_key()`. Used by `.ge_build()` and `.ge_frame_term_messages()` for
duplicate loci and slots.

The helper computes `x * (max(y) + 1) + y`. It falls back to string keys for
missing, negative or fractional `y`, but never checks that the resulting integers
remain exactly representable. Distinct pairs can therefore receive the same key
above `2^53`, despite both input columns being valid integers.

Reproduced with the current implementation:

```r
k <- tidybreed:::.ge_pair_key(
  c(4194304L, 4194304L, 1L),
  c(3L,       4L,       2147483647L)
)
format(k, digits = 22)
# "9007199254740996" "9007199254740996" "4294967295"
k[1] == k[2]
# TRUE: (4194304, 3) and (4194304, 4) are distinct pairs.
```

The third row sets the radix to `2147483648`; the first two keys round to the
same double. At a sufficiently large model with sparse high locus IDs, this can
turn two distinct loci into a false repeated-locus violation and reject a valid
write. The comment's `max(x) * max(y)` condition is also insufficient: the actual
expression includes the extra radix increment and final `y`.

**Recommended change:** retain the fast path only when a conservative upper bound
on the entire encoded key stays below `2^53`; otherwise use the existing string
fallback. Add a hand-authored uniqueness test at this boundary and an ordinary
small-key case. These test the helper's mathematical contract, not old output.

**Severity and limit:** this is a directly reproduced helper defect, not an
end-to-end failure reproduced with millions of terms. Its callers use compact
group positions, so reaching this example needs over four million groups as well
as a very high locus ID. The planned 124,750-pair / 500-locus workload is safely
below the boundary, including with integer-max locus IDs. This does **not** block
the planned 5b work.

## What the code review verified

| Area | Assessment |
|---|---|
| Writer scalar checks | `.ge_build()` chooses the first failing input term, then its first rule. It builds members in term order and ascending locus-ID order with integer slots. The new test pins cross-rule precedence and user labels. |
| Member and origin messages | The vectorised rules interleave messages in the previous row/rule order. Unknown origin match types still report only their own rule. |
| Structural validation | Grouped slot, locus and origin checks retain ordered messages on the tested frames. Family predicates use pre-split rows, with family and variant ordering preserved for unique stored term IDs. |
| Atomic writes | `.ge_commit()` and its transaction are unchanged. Candidate validation and **whole stored table** validation still run; target writes still use the hook inside the transaction. The conflict/rollback test passes. |
| Reverse conversion | `.stored_to_functional()` sorts contributions by term and within-term position before accumulation. Dominance cancellation still uses absolute contributions and counts; pair summation and cancellation remain intact. Mixed-shape numerical comparisons were bit-identical. |
| Scope and target helpers | The additional `.dae_*` and `.gev_*` refactors use position-based grouping while retaining shape-based target classification and scope matching. Additive, removal and prevalence tests pass. |
| Evaluator preparation | `.gev_variant_map()` and `.gev_preflight()` reach common-family data by position. Origin signatures remain aligned to term positions. Evaluation SQL and the exact deterministic accumulator are unchanged. |
| Scope of the change | No API, schema, `NAMESPACE` or manual-page changes. The only test-file change is the added writer-order contract file; no existing test file was edited. |

The extra evaluator and additive-helper changes are justified extensions of 5a:
they remove the same repeated-lookup costs from paths that the upcoming generator
will use. The roxygen comment repair is unrelated to effect math.

The contract test was committed with the refactor rather than in a preceding
commit. The results document explains that it was run against 0.75.2 first.
That historical execution is an implementation report, not something this review
independently certifies; the current contract tests pass.

## Independent verification

### Targeted tests

Executed from the repository root:

```r
pkgload::load_all(".", quiet = TRUE)
testthat::test_local(
  ".",
  filter = "genome-effect|extract_genetic_variance|define_additive_effects|remove_generated_effects|prevalence-threshold",
  reporter = "summary"
)
```

**Result: exit 0; all 12 selected files passed, with no reported failures,
warnings or skips.** Coverage includes the writer and its error-order contract,
builders, schema, hand-derived fixtures, evaluator, thread determinism, extractor,
both additive-generator files, prevalence and generated-effect removal.
Log: `/private/tmp/tidybreed_phase5a_targeted.log`.

### Development comparisons

The scratch script sources the five pre-change R files into a separate
environment and compares old/current results using `identical()`. It remains
outside the repository, respecting the prohibition on committed golden-from-old
tests. Its cases cover shuffled/nonconsecutive IDs, additive/dominance/indicator/
A×A shapes, coefficients across several scales, high integer locus IDs,
malformed slots and scalar values, scoped origins, unusual line names, target
classification, replacement modes and evaluator maps/preflight.

**Result: 7,020 comparisons identical, zero mismatches** (300 conversions,
900 structural checks, 900 target/scope-label checks, 3,600 replacement checks,
1,200 map/preflight checks and 120 writer builds/errors). The separate precision
boundary probe reproduces finding 1.

Script: `/private/tmp/tidybreed_phase5a_review_probes.R`.
Log: `/private/tmp/tidybreed_phase5a_review_probes.log`.

### Writer benchmark

Ran the committed `dev/benchmarks/benchmark_genome_effect_writer.R` at every
planned size. R 4.5.3; Darwin 24.6.0; x86_64; DuckDB 1.5.5, 16 threads; 500 loci
on five chromosomes. The sandbox did not expose a core count (`detectCores()`
returned `NA`). Targeted tests overlapped this run, so these are development
measurements rather than isolated timing estimates.

| Pairs | Setup (`aa_terms`) | Write | Time per pair | Peak R heap |
|---:|---:|---:|---:|---:|
| 1,000 | 0.024 s | 0.244 s | 244 µs | 165 MB |
| 4,000 | 0.152 s | 0.470 s | 118 µs | 179 MB |
| 16,000 | 0.430 s | 0.747 s | 47 µs | 194 MB |
| 64,000 | 1.653 s | 3.601 s | 56 µs | 243 MB |
| 124,750 | 3.260 s | 7.193 s | 58 µs | 390 MB |

- Time-per-pair growth from 4,000 to 64,000: **0.48×**, below the 2× gate.
- All-pairs write: **7.2 s**, below the 60 s gate.
- Replacement of 16,000 pairs into 20,200 stored terms, three traits/owners and
  100 line-scoped variants: **1.2 s**, peak 280 MB.
- Reverse conversion of those 16,000 terms: **0.17 s**.

Log: `/private/tmp/tidybreed_phase5a_writer_benchmark.log`.

## Readiness for 5b and limits

The shared writer, conversion and evaluator setup are ready for the planned
generator. Preserve the existing `.ge_build()` / `.dae_stack()` / `.ge_commit()`
path, including full-table validation and atomic target writes. The planned
additive-only route must still share Part A's calibration **and storage**, since
5a does not change its retained zero-row behavior.

The pre-5b grep gate also passes: no `define_genome_effects(` references in
`R/`, `man/`, `tests/` or `vignettes/`. The generator, new solver and rank-note
helper are still 5b work; this review does not certify their unimplemented gates.

Two performance qualifications should remain explicit:

- The speed evidence is for growing numbers of common A×A families. Family
  validation still compares pairs of variants within a family, and exceptional
  scoped paths retain some per-term work. It is not a proof that every possible
  high-order or heavily scoped model has linear runtime.
- Large realised evaluation remains expensive. I did not rerun the extractor's
  2,000-individual all-pairs benchmark; its timings and in-memory spill finding
  are reported by the implementation document. Keep ordinary 5b/5c correctness
  fixtures small, and profile forward `.noia_to_stored()` in 5b as planned.

I did not rerun the full package suite, regenerate documentation or run pkgdown.
The implementation document's 4,295-expectation full-suite result is separate
from the independent targeted checks above. This review changes only this
Markdown file; it makes no production or test-code changes.
