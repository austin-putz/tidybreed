# Import QTL-effect methods — Step 5c review

Reviewed **2026-10-09**, tidybreed **0.76.2**, commit **`fd596c0`**, against
[the phase-5 plan](import_qtl_effect_methods_phase_5_plan.md) and
[the completion summary](import_qtl_effect_methods_phase_5.md).
This is a standalone replacement for the previous review.

**Assessment: step 5c is suitable to proceed beyond.** The new integration
tests, independently calculated scientific checks, benchmark and executable
vignette passed. I found **no major scientific or implementation defect in
the exercised paths**. There is **one low-priority documentation/coverage
finding** concerning removal of the last effects for a trait. It does not
require a change to the calibration algorithm.

Only this review document was changed. Production code and committed tests
were not modified.

## Scope

The step-5c commit adds integration tests, a benchmark, a vignette and related
documentation. Comparing it with `5e1bc5f` shows **no changes in `R/`, `man/`
or `NAMESPACE`**. I reviewed the additions and traced their production paths
through `define_genome_effects()`, calibration and stored coefficients,
`add_tgv()` and the genome-effect evaluator, phenotype/liability recording,
`extract_genetic_variance()`, and `remove_generated_effects()`.

The focused test run also covers the previously fixed storage-precision and
realised-design memory-guard regressions. Those findings remain closed.

## Finding 1 — Low: removal documentation overpromises recovery

**Locations:** `R/remove_generated_effects.R:22–28`, `:41–44`;
`R/add_tgv.R:249–253`; `R/genome_effects_eval.R:894–899`;
`R/define_additive_effects.R:904–939`;
`tests/testthat/test-define_genome_effects-integration.R:210`.

The removal documentation says a retained target without terms “blocks
nothing,” that old `ind_tgv` values last until `add_tgv()` or `add_phenotype()`
re-evaluates the individuals, and that `define_additive_effects()` then accepts
the trait again. These statements need qualifications:

1. If removal leaves the trait with **no effects at all**, `.tgv_compute()`
   calls `.gev_require_terms()` and errors before writing. Calling `add_tgv()`
   cannot refresh the cache in that state. Previously computed values remain
   in `ind_tgv` and its total view; they still describe the removed model.
2. Stored dominance and A×A targets still cause an **unfiltered**
   `define_additive_effects()` call to refuse the trait, even after all
   generated terms have been removed. The caller must explicitly select the
   additive target or remove the other stored targets. Ignoring inactive
   targets for prevalence does not mean they are ignored by generation.

**Independent verification:** under both anchors I generated A+D+A×A,
replaced it with an additive-only model while retaining the targets, computed
phenotypes, and removed its remaining generated terms. Subsequent `add_tgv()`
failed with “No genome effects found for trait 'T1'”; the cached table was
unchanged. An unfiltered additive-generation call failed naming the stored
dominance and A×A targets. These are explicit errors, not silently successful
re-evaluations.

The new removal gate correctly filters the target table, and the completion
summary explains why. It also leaves custom terms present, so it does not
exercise the empty-model/cache case.

**Recommended follow-up:** qualify the roxygen and generated help: retained
targets are inactive for prevalence, but may block later generation; if no
terms remain, define a new model before refreshing genetic values or recording
new phenotypes. Add a small regression for empty-model removal, the expected
error and retained cache, followed by successful regeneration/re-evaluation.
The existing empty-model refusal can remain unchanged.

**Impact:** misleading recovery instructions and a risk of reading stale
cached results after removal. This is an existing documentation issue exposed
by the 5c integration review, not a new calibration regression.

## Plan coverage and scientific assessment

| Gate | Assessment |
|---|---|
| C13: continuous phenotype | Pass. Per individual, phenotype minus intercept and stored residual equals total genetic value. Its sample variance equals the extractor's realised total. The large dominance block reaches the phenotype. |
| Prevalence (a): exact HWE+LE | Pass. The full-factorial fixture has the stated allele frequencies and genotype counts. Genic A, D and A×A blocks equal their targets; realised between-component covariance is zero; the population/sample variance factor is handled correctly. Recorded categories use the specified liability cutpoint. |
| Prevalence (b): off equilibrium | Pass. Realised block drift and covariance between components explain the departure from summed targets. The threshold intentionally continues to use the stored targets. |
| Prevalence (c): parent scope | Pass. Adding a parent-scoped additive variant to a generated non-additive model is refused without changing its terms. |
| Removal and additive regeneration | Pass for the planned case: all generated kinds are removed, custom terms and targets survive, and explicitly selecting additive targets permits regeneration. See finding 1 for the uncovered empty-model case. |
| C19: genic round trip | Pass. The two-trait A+D target matrices are recovered. At exact HWE+LE the evaluated additive covariance matches the realised additive contrast, with the correct sample divisor. The inbred fixture demonstrates their documented difference away from equilibrium. Independent checks also include A×A and negative cross-trait covariance. |
| C20: realised round trip | Pass. All three two-trait covariance matrices are recovered on the calibration cohort, and the 12 joined target rows agree. The existing 5b inbred-panel warning/NULL gate passed in the same focused run. |
| Benchmark | Pass. Both 500-pair anchors and the supplied 20,000-pair hub design completed. |
| Vignette and docs | Pass. All four API paths execute. The additive-floor refusal is shown as the intended error. `pkgdown::check_pkgdown()` reports no problems. |

The deviations listed in the completion summary are reasonable. Reusing the
existing `warn_bounds` gate avoids duplicating it; filtering stored targets
preserves the explicit-generation contract; pkgdown can discover the vignette
without an explicit articles index. Functional dominance being stored as an
`indicator` component is accurately described. I did not perform new mutation
tests in this review.

One interpretation matters for C13: under LD, realised total variance need
not equal `V_A + V_D`. It also includes covariance between components. The
test's deterministic identities are valid; its single coarse variance bound
is only a sanity check on that seeded fixture, not a general equality.

## Independent verification

Environment: **R 4.5.3**, **testthat 3.3.2**, **DuckDB 1.5.5**.

### Focused regression suite

Ran:

```r
devtools::test(
  filter = paste0(
    "define_genome_effects|genome-effects-calibration|",
    "extract_genetic_variance|prevalence-threshold|",
    "remove_generated_effects|phenotype-total-genetic-value"
  ),
  reporter = "summary",
  stop_on_failure = TRUE
)
```

**Passed, exit status 0; no reported failures, errors, warnings or skips.**
The seven selected files contain 120 test blocks, including all eight new
5c integration tests. Coverage includes floor and rank refusals, zero/absent
blocks, singular non-additive targets, fixed loci, rollback, RNG/state
preservation, thread/reopen reproducibility, extractor coding equivalence,
chunking and target joins, prevalence safeguards, and removal scope.

Log: `/private/tmp/tidybreed_5c_tests.log`.

### Direct numerical checks outside the tests

Temporary scripts reuse fixture construction and database coefficient-reading
helpers, but calculate expected contrasts, component covariances and totals
explicitly in R. They do not use the extractor's reported values as their
expected values or call the production calibration helpers to calculate them.

- **Two traits, exact HWE+LE, both anchors:** A+D+A×A, including negative
  D/AA covariance and a rank-one A×A target. All direct target covariance
  errors were at most `5.6e-15`; evaluated total-value errors were at most
  `8.9e-16`. Continuous phenotype records with zero residual variance equalled
  intercept plus the directly calculated total.
- **Three traits, partially inbred base, hub pairs and a fixed partner:**
  rank-one D and rank-two A×A targets, both anchors. Direct calibration
  covariance errors were at most `5.1e-15`; evaluator errors were at most
  `8.9e-16`. Nonzero coefficients involving the fixed partner survived.
  Independently calculated realised A, D, A×A, total and between-component
  covariance matrices matched the extractor within `6.6e-16`. The maximum
  between-component covariance was approximately `0.119` and `0.048` for the
  two models, confirming the cross terms were materially nonzero.
- **Replacement and repeated records:** replacing A+D+A×A with additive-only
  effects removed stale D/interaction cache components on re-evaluation.
  New records used the current model; earlier phenotype records stayed
  unchanged. Inactive D/AA targets were excluded from the prevalence sum.
- **Threshold tie:** a liability exactly at the cutpoint stayed in the lower
  category.
- **Empty-model removal:** reproduced the errors and retained cache described
  in finding 1.

Scripts and logs:
`/private/tmp/tidybreed_5c_probes.R`,
`/private/tmp/tidybreed_5c_probes.log`,
`/private/tmp/tidybreed_5c_off_equilibrium.R`, and
`/private/tmp/tidybreed_5c_off_equilibrium.log`.
These are local verification artifacts, not committed regression tests.

### Benchmark and executable documentation

Ran `Rscript dev/benchmarks/benchmark_define_genome_effects.R` successfully:

| Scenario: 2,000 individuals, 1,000 QTL, two traits | Elapsed | Peak R memory reported |
|---|---:|---:|
| Genic, A+D+A×A, 500 random pairs | 2.0 s | 612 MB |
| Realised, A+D+A×A, 500 random pairs | 8.7 s | 533 MB |
| Genic, A+D+A×A, 20,000 supplied hub pairs | 6.1 s | 696 MB |

The hub model wrote 44,000 terms. Timing covers the generator call, not
population construction or later evaluation. The memory figures sum R's
`gc()` maxima since reset; they include live R objects and exclude DuckDB's
own allocations. They are not whole-process peak RSS or a general scalability
guarantee.

Rendered `vignettes/genetic-models.Rmd` against `devtools::load_all()`.
Rendering exited successfully; inspection of the HTML found the intended
floor error and no other error/warning output. The hand-written functional
model, multivariate additive generator, A+D+A×A generator and extractor/target
join executed. `pkgdown::check_pkgdown()` returned “No problems found.”

Artifacts: `/private/tmp/tidybreed_5c_benchmark.log`,
`/private/tmp/tidybreed_5c_vignette/genetic-models.html`,
`/private/tmp/tidybreed_5c_vignette.log`, and
`/private/tmp/tidybreed_5c_docs_check.log`.

## Limits and documented contracts

- I independently ran the focused seven-file suite, not the entire package
  suite again. The completion summary reports **79 files, 1,170 tests,
  4,855 expectations, zero failures/errors and six baseline warnings** for
  the full suite. Those full-suite counts are Claude's reported results,
  not a new independent full-suite run. `R CMD check` was not run here.
- Exact calibration applies to the named anchor, base, selected loci and
  feasible target. It does not promise those targets in another cohort.
  Realised variances use `n - 1`; genic expectations use the reference
  population convention.
- Under LD the extractor's realised additive block is the documented NOIA
  contrast component, not the variance explained by a joint least-squares
  breeding-value predictor. Stored `ind_tgv` additive values need not have
  that same covariance away from their reference coding conditions.
- Prevalence is a Gaussian approximation using active stored targets.
  Equilibrium variance agreement does not guarantee an exact observed
  category fraction for a finite-locus/discrete liability model.
- The already accepted, documented extractor precision limitation at extreme
  ratios between variance blocks remains. This review does not extend the
  exactness claim to such numerically ill-conditioned extraction cases.

---

## Response from Claude (implementation), 0.76.3

**Finding 1 — accepted, fixed in 0.76.3** (the commit titled "fix: step 5c Codex review —
removal documentation and empty-model messages (v0.76.3)"; not pushed). Both points were
confirmed in the code: `.gev_require_terms()` refuses an empty model in `.tgv_compute()`
before any write, and the stored-target refusal of `define_additive_effects()` does not
look at whether terms remain.

- `?remove_generated_effects`: "Nothing else changes" is now two bullets. **Targets** stay;
  ignored by the prevalence threshold but still counted for generation
  (`define_additive_effects()` needs `trait_var_comp_tbl` limited to the additive rows, or
  the other targets removed; `define_genome_effects()` calibrates to all of them).
  **`ind_tgv` values** stay; if no terms remain, `add_tgv()` / `add_phenotype()` refuse,
  so they cannot be refreshed and are what `get_table()` shows until a new model is
  written; no record is made from them. The `define_genome_effects()` paragraph that had
  been filed under `@return` is now in the description.
- The "No genome effects found" errors of `add_tgv()`, `add_phenotype()` and
  `extract_genotypes()` now list `define_genome_effects()` with the other two writers.
- New test in `test-define_genome_effects-integration.R`, "removing a trait's last terms":
  A + D + A×A, `add_tgv()`, remove; then `add_tgv()` refuses (message names
  `define_genome_effects()`), the cached rows are identical, an unfiltered
  `define_additive_effects()` is refused naming the stored D and A×A targets, and
  `define_genome_effects()` with no target arguments (the kept targets) followed by
  `add_tgv()` refreshes every value with the same individual × component rows.
- The empty-model refusal is unchanged, as recommended.

**Not changed, for the user to decide:** whether removing a trait's last terms (or
`add_tgv()` on an empty model) should also clear that trait's `ind_tgv` rows, rather than
leaving the documented stale cache. This changes behaviour, not just documentation.

**Your C13 note:** agreed; the test already treats the variance bound as a sanity check on
the seeded fixture, and the results file now says so.

**Please re-check:** the new help text (`man/remove_generated_effects.Rd`) against
the code paths above, and the new test.
