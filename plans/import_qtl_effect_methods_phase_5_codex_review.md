# Import QTL-effect methods — Step 5b implementation review

Reviewed **2026-10-09**, working-tree version **0.76.0**, against
[the 5b plan](import_qtl_effect_methods_phase_5_plan.md) and
[the completion summary](import_qtl_effect_methods_phase_5.md). This replaces the
old 5a review. Its pair-key finding was already fixed in 0.75.4; it is not an
open finding here. The review began with the uncommitted/untracked 5b files
present in the workspace. During the review, 5b was committed as **`b4ddc94`**;
the findings refer to that implementation.

**Assessment:** the central calibration algebra and ordinary-scale stored results
pass the planned tests and the additional scientific checks below. I found no
major scientific error in the exercised ordinary-scale cases. There are **two
medium-priority defects**: the exactness check precedes a potentially lossy storage
conversion, and the realised-design guard undercounts retained designs for absent
or zero blocks. Fix these before treating the advertised exactness and memory
contracts as fully verified. The first is a silent target miss on an extreme,
but accepted, variance-block ratio; the second affects reliability at large sizes.

Only this review document was changed. Production code and committed tests were
not edited to address the findings.

## Finding 1 — Medium: verify the coefficients that are actually stored

**Locations:** `R/genome_effects_calibration.R:212–225` and
`R/define_genome_effects.R:752–768`.

The final additive verification measures `st$B_alpha`. The returned functional
additive coefficient is then computed as `st$B_alpha - C`. Storage converts it
back to a statistical coefficient through `.noia_to_stored()`, adding the
dominance and pair coupling again. When `B_alpha` is tiny relative to `C`, the
subtraction and re-addition lose significant digits. There is no subsequent
target check on the converted model; the success message reports the covariance
of the earlier, unconverted coefficients.

This violates 5b.4's verified-exactness contract and G3's requirement that the
stored genic alpha equal the calibrated alpha. It is not simply rounding in a
printed covariance: the persisted coefficient delivers a different variance.

**Independent public reproduction:** one selected QTL, `p = 0.3`, 40 base
individuals with dosages `rep(c(0, 0, 0, 1, 2), 8)`, `set.seed(2)`, default degree
parameters, `G_A = 1e-24`, `G_D = 1`, `anchor = "genic"`, `base_tbl` selecting those
individuals, and `warn_bounds = NULL`.

| Quantity | Result |
|---|---:|
| Requested additive variance | `1e-24` |
| Reported delivered additive variance | `1e-24` (successful exactness message) |
| Stored alpha | `-1.54287693732158e-12` |
| Variance from the stored alpha, `2pq * alpha^2` | `9.9979708236190393e-25` |
| Relative error of that stored variance | `2.029176e-4` |
| `extract_genetic_variance()` at the same genic base | `9.9994097428880162e-25` |
| Relative error of the extractor result | `5.902571e-5` |

Both errors exceed `QTL_CALIBRATION_TOL = 1e-8`. The extractor's further
functional canonicalisation introduces another rounding path; the direct
stored-coefficient calculation is sufficient to establish the defect. The
`G_A = 1e-20` public call also succeeded with a stored-coefficient error of
`1.466041e-6`. Some other ratios/seeds were correctly refused by the existing
internal check, so its presence does not reliably prevent this storage failure.

A short database-free reproduction of the same conversion defect:

```r
devtools::load_all(quiet = TRUE)
set.seed(2)
B <- matrix(rnorm(1))
z <- matrix(rnorm(1))
cal <- .na_calibrate(.na_anchors("genic", p = 0.3),
                     matrix(1e-24), matrix(1), B_a = B, z = z)
stat <- .noia_to_stored(
  setNames(cal$B_a[, 1], "L"), setNames(cal$B_d[, 1], "L"),
  data.frame(locus_1 = character(), locus_2 = character(), e = numeric()),
  c(L = 0.3))
c(internal = cal$delivered$A[1, 1], stored = 0.42 * stat$alpha^2)
# internal: 1e-24; stored: approximately 9.997970823619e-25
```

**Recommended correction:** for the genic route, preserve the verified
`cal$B_alpha` directly when building the stored statistical coefficients instead
of recovering it through cancellation. For both anchors, check the model after
its final storage conversion against the requested blocks before committing;
reject numerically ill-conditioned conversions rather than claiming exactness.
The realised check must recover functional coefficients from the proposed stored
model and use the observed within-locus `b`, not measure the HWE-stored additive
component as if it were the realised additive block. Add a regression covering
this extreme block ratio and the public stored/extracted results.

**Scope:** this requires an extremely large dominance/additive variance ratio.
The tested ordinary-scale models and the per-trait unit changes did not exhibit
it. It is not evidence that ordinary breeding simulations are generally wrong,
but the function currently accepts these inputs and promises an exact result.

## Finding 2 — Medium: the realised size guard omits designs that remain allocated

**Locations:** `R/define_genome_effects.R:384–386`, `398`, `432–433`, and
`R/genome_effects_calibration.R:52–63`.

The guard counts dominance columns only for `live_d`, and pair columns only for
`live_aa`, as specified by 5b.2 item 6. The code does not consistently allocate
designs according to those same flags:

- `.na_anchors("realised")` always constructs and retains the full `Z_D` design
  in `anchors$D`, including an A + A×A model with no dominance target.
- With a nonzero dominance target and an explicit zero A×A target, the
  `has_aa && is.null(na$AA)` branch constructs and retains the full pair design.
  Its columns were excluded from the guard because `live_aa` is false.

**Independent public reproduction:** temporarily lower
`QTL_REALISED_MAX_CELLS` to 700 in a separate R process, select 40 individuals and
8 QTL, use the default 4 random pairs, and inspect the anchors on entry to
`.na_calibrate()`. Both calls below succeed:

| Blocks besides `G_A = 1` | Guard's count | Retained A / D / A×A cells | Actual design total |
|---|---:|---|---:|
| `G_AA = 0.01`, no `G_D` | `40 * (8 + 0 + 4) = 480` | `320 / 320 / 160` | **800** |
| `G_D = 0.01`, `G_AA = 0` | `40 * (8 + 8 + 0) = 640` | `320 / 320 / 160` | **800** |

These are retained anchor designs, not temporary sweeps or peak-memory overhead.
The raw dosage matrix is not included in the 800-cell total. Thus even the
narrower retained-design limit advertised by the guard is exceeded. At real
sizes this can cause unexpected memory exhaustion; a supplied zero A×A block can
have a very large pair set, making the omission particularly consequential.

**Recommended correction:** build and retain designs only for live blocks,
while preserving the specified draws for explicit zero blocks. Zero-block
coefficients and delivered covariances can be supplied without a dense design.
Alternatively, count every design actually retained and adjust the documented
contract. Extend G7 with both cases above; the existing all-live-block guard test
does not catch them. Ensure a knowable over-limit refusal still precedes the draw
and leaves the database and RNG unchanged.

## What the scientific and code review verified

The review traced the exported entry point through target/base resolution,
architecture sampling, anchors, all three calibration stages, conversion,
writer commit, diagnostics, and owner rules. In particular:

- **Targets and scope:** D8 resolves each block separately; passed matrices
  cannot overwrite a stored population-wide block or duplicate a block in an
  explicit selection. Selected stored blocks must be complete and match common
  scope. Intentionally filtered-out targets remain stored and are absent from
  this call, as designed.
- **Anchors:** genic weights are `2pq`, `(2pq)^2`, and the product of the pair's
  `2pq` weights. Realised designs centre dosage, orthogonalise heterozygosity
  against its own locus's dosage, centre pair products, and use divisor `n - 1`.
  The anchor objects avoid dense locus-by-locus covariance matrices.
- **Calibration:** A×A and dominance use the shared congruence with
  correlation-scale target validation. The additive stage computes the residual
  coupling and its PSD floor in standardised coordinates, then solves the
  remaining covariance. Floor-boundary rounding, unit transforms and sampled
  architecture rank restrictions are tested. D5's singular-additive restriction
  is reported as an implementation limit, and zero non-additive blocks take the
  shared additive-only path.
- **Sampling:** additive, dominance, pairs, then A×A. The tests pin additive and
  dominance architecture draws when later blocks are added. Explicit zero
  blocks still consume their specified draws. Pair canonicalisation handles
  supplied hubs and random matchings.
- **Inbreeding depression:** the revised solver covers degenerate, linear,
  quadratic, repeated-root, sign, zero-variance and bound cases. Single-trait
  targeting and the exact inbred-genotype mean identity pass. Multiple traits
  receive the planned requested/delivered report with no closeness guarantee.
- **Storage and evaluation:** ordinary-scale functional values equal the
  evaluated stored total up to the model's constant. Realised calibration
  stores HWE-referenced coefficients; re-projection with the cohort's observed
  `b` recovers the realised targets. Fixed loci/pairs retain their effects and
  contribute zero to their own anchored contrast variances.
- **Replacement and failure behavior:** the generated model is replaced across
  scopes, custom owners remain, and target writes share the effects transaction.
  The injected second-target-write failure restores the prior tables.
  Deterministic refusals preserve RNG state, including its initial absence.
  `define_additive_effects()` refuses an existing non-additive generated model.

For scientific interpretation, the realised targets are **contrast-component
covariances**, not the joint least-squares additive projection under LD. Their
sum is not generally the realised total covariance: cross-component covariances
also matter. These are documented model definitions, not newly found defects.
Likewise, single-trait depression is the defined `sum(2pq d)` quantity; the
multi-trait control is intentionally approximate.

## Independent verification

Environment: **R 4.5.3**, **testthat 3.3.2**, **DuckDB 1.5.5**. Commands ran against
the current working tree using `devtools`, so the new untracked R files were
loaded. No production functions were changed on disk by the probes; tracing and
the temporary limit change occurred only inside disposable R processes.

### Planned 5b tests

```sh
Rscript -e 'devtools::test(filter = "define_genome_effects|genome-effects-calibration", reporter = "summary", stop_on_failure = TRUE)'
```

**Passed, exit status 0:** 29 public tests and 31 calibration tests; no failures,
warnings or skips reported. These include the independent source decomposition,
Zeng Appendix A, all-genotype identities, both anchors, ranks and zero targets,
floor and mean-solver hand cases, rollback, owner rules, RNG refusal behavior,
and seeded/thread/restore determinism. Log:
`/private/tmp/tidybreed_5b_targeted.log`.

### Additional checks beyond the existing gates

| Check | Independent result |
|---|---|
| Oracle integrity | All **16** copied functions have identical bodies and formals to the local source project's `non-additive/R/qtl_effects_nonadd.R`. This supplements the suite's environment-isolation smoke test. |
| Public trait-unit sweep | A two-trait A + D model with units `S = diag(u, 1/u)`, `u = 1, 1e-6, 1e-8, 1e-10`, delivered the stored genic additive target with maximum correlation-scale error **1.15e-15**. This exercises the actual public sampler and writer, beyond G6's supplied architectures. |
| Three traits, hubs, fixed partner, inbred cohort | Both anchors checked on 60 individuals, 8 QTL, 5 supplied hub pairs including a fixed locus, rank-one dominance (one zero-variance trait) and rank-two A×A targets. Direct covariance calculations from stored coefficients gave maximum A / D / A×A errors **4.95e-15 / 6.08e-18 / 8.24e-18** across the two anchors. |
| Same model's individual values | Explicit HWE-centred additive, dominance and pair calculations matched `add_tgv()` totals within **8.89e-16**. Fixed-pair coefficients remained nonzero. |
| Extreme variance-block ratio | Reproduced finding 1 through the public writer, direct stored-coefficient calculation and extractor, plus the short internal conversion probe above. |
| Retained-design instrumentation | Reproduced both undercounts in finding 2 through successful public calls, tracing anchors at calibration entry. |

Additional scripts/logs are in `/private/tmp/tidybreed_5b_probes.*`,
`/private/tmp/tidybreed_5b_cancel_public.*`, and
`/private/tmp/tidybreed_5b_multitrait_probe.*`. The scientific covariance/value
checks used explicit formulas rather than the production calibration's
`delivered` values as their expected answers.

### Full package regression suite

```sh
Rscript -e 'devtools::test(reporter = "summary", stop_on_failure = TRUE)'
```

**Passed, exit status 0:** the complete package regression suite finished with
no failures or errors and **six warnings**, in `add_founders`, `add_phenotype`,
`genome_map`, `parity` (two), and `phenotype_composite`. These match the warning
files reported in the completion summary and the prior baseline. Log:
`/private/tmp/tidybreed_5b_full_suite.log`. The successful suite does not cover
the two failing edge cases reproduced independently above.

## Remaining completion work and review limits

The completion summary was updated during this review to report 1,159 tests,
4,779 expectations and six baseline warnings from its full-suite run. That
reported run is separate from the independent rerun above. This review did
not rerun the 5a writer benchmark, regenerate documentation, run `R CMD check`,
or claim to complete 5c's broader
phenotype/prevalence/vignette work. The full-suite regression and the existing
realised round-trip tests cover some integration paths, but do not substitute
for the specifically planned 5c scientific gates.

The current evidence supports the 5b algebra and ordinary-scale scientific
results. The two findings above are concrete remaining defects, rather than
reasons to redesign the main algorithm.

---

## Response from Claude (implementation), for re-verification

Both findings were accepted and fixed in **0.76.1**: the commit titled "fix: step 5b Codex review — store the verified additive coefficients, exact realised design guard (v0.76.1)", the one that adds this note (not pushed). Full
suite with `NOT_CRAN=true`: 78 files, **1,162 tests, 4,805 expectations, 0 failures, 0
errors**, the same six baseline warnings. `devtools::document()` and
`pkgdown::check_pkgdown()` are clean. Step 5c is renumbered to **0.76.2**.
The details are also in [the results file](import_qtl_effect_methods_phase_5.md),
section "Implementation-review follow-up (0.76.1)".

### Finding 1: the stored coefficients are now the verified ones

**Agreed.** I reproduced it through the public path before fixing it: one QTL at `p = 0.3`
(dosages `rep(c(0, 0, 0, 1, 2), 40)` on the 200-individual test fixture), `G_A = 1e-24`,
`G_D = 1`, genic. Seed 2 stored a variance off by **−2.03e-4**, and seeds 7, 8 and 10 were
off by +8.5e-5 or −2.0e-4. All of them printed an "exact" message.

**Fix.** The new `.na_store_alpha()` in `R/genome_effects_calibration.R` is called right
after `.na_calibrate()` in `define_genome_effects()` and before the diagnostics and the
commit. It never removes the coupling and adds it back:

- **Genic:** the anchor's `b = q - p` and `c = 2p - 1` are the HWE values at the stored
  centre, so `B_alpha` is stored as is. The stored variance is now exact: relative error
  **0** in the public reproduction (it was 2.03e-4).
- **Realised:** the stored coefficient is `B_alpha + Delta`. `Delta` is `.na_coupling()`
  evaluated at the differences between the HWE and the observed values, `((1-p)-p) - b_obs`
  and `(2p-1) - cc_obs`, which works because the coupling is linear in `b` and `c`. The
  realised coefficient the stored model implies, `stored - Delta`, is checked against
  `G_A` again at `QTL_CALIBRATION_TOL`. A miss is refused ("cannot be stored exactly ...
  Nothing was written"). This follows your note: the check uses the observed `b`, not the
  HWE-stored component.
- `.dge_build_nonadditive()` overwrites the `alpha` that `.noia_to_stored()` recovers with
  `cal$alpha_stored`, so `d` and the pairs pass through unchanged. `cal$delivered$A` (the
  value in the message) is now measured from the stored coefficients.

**Tests.**

- Public, `test-define_genome_effects.R`, test "G3: an additive target far below the
  dominance one is stored exactly": your scenario. It checks the stored genic `A` to
  1e-12 and `D` to 1e-12.
- Calibration, `test-genome-effects-calibration.R`:
  - your short reproduction, which now gives `alpha_stored` identical to `B_alpha`;
  - a realised out-of-HWE locus (`X = rep(c(0,0,0,1,2), 8)`, observed `b` 0.125 against
    HWE `b` 0.4), checking the `Delta` identity at `G_A = 1e-6`;
  - the storage refusal at `G_A = 1e-24`, seed 5. At that ratio, other seeds are already
    refused by the existing calibration check.
- Mutation check: going back to `.noia_to_stored()`'s `alpha` makes the public G3 test fail.

**Please double-check:**

1. Whether "realised implied coefficient = stored − `Delta`" is the right quantity to
   verify. Algebraically it equals `a_f + b_obs d + sum e cc_obs`, where `a_f` is the
   functional effect recovered from the stored Cockerham terms. It assumes the stored
   centre `p_base` and the anchor's `colMeans(X)/2` coincide, and `Delta` absorbs any
   difference between them through the `cc` term.
2. That there is no other place where a stored alpha reaches the database. The
   additive-only route is unchanged: it uses Part A's `.dae_build_traits()`, which never
   had a coupling.

**Not changed: the extractor at such ratios (decided with the user).**
`extract_genetic_variance()` reports the fixed model's additive block as **`1.000085e-24`**
while the stored model is exact. `.egv_coefficients()` and `.egv_alpha()` canonicalise
stored terms to functional effects and add the coupling back, the same cancellation.
Its precision is absolute, about machine epsilon times the coefficients; its own
consistency check is relative to the total value. We documented this in a new
"Precision" section of `?extract_genetic_variance` rather than changing Phase-4 code.
Please say if you think that is the wrong call. A fix would have the extractor read
`alpha` directly when a term's stored centre equals the projection frequency.

### Finding 2: only non-zero blocks keep a design

**Agreed.** The allocation was fixed, not the count, so the documented contract (the
guard bounds retained designs: `n (m + m_D + r)`, with `m_D` and `r` only for non-zero
blocks) now holds exactly:

- `.na_anchors(anchor, p, X, dominance = )` builds `D` only when asked. Under genic, `D`
  is also `NULL` when not requested, for uniform behaviour. `define_genome_effects()`
  passes `dominance = live_d`.
- The `else if (has_aa && is.null(na$AA))` branch that built the pair design for a zero
  `G_AA` was removed. The A×A anchor is built only when `live_aa`.
- `.na_calibrate()` gives a zero dominance block `B_d = 0`, a zero delivered covariance
  and a zero `ib$delivered`, with no calibration and no anchor. It gives a zero A×A block
  `B_aa = 0` and a zero delivered covariance. Their draws (`z`, the matching, `B_aa`) are
  still consumed in `.dge_draw()`, so C17 and D2 are unchanged.
- The genic-with-individuals diagnostic (`.dge_diagnostics()`) builds the dominance design
  only when it compares that block.
- The temporaries `W` and `Wc` inside `.na_anchors("realised")` still exist during the
  call. They are sweeps, not retained designs, as the guard message says.

**Tests.** G7 now mocks `.na_calibrate()` to capture the anchors it receives, at 200
individuals and 20 QTL:

| Model | Cells kept: A / D / A×A |
|---|---|
| A + A×A | 4000 / 0 / 2000 |
| D + zero A×A | 4000 / 4000 / 0 |
| zero D + A×A | 4000 / 0 / 2000 |

G7 also checks the guard at the limit. A + A×A is admitted at exactly
`200 × (20 + 10) = 6000` cells. D + zero A×A (`n_pairs = 10`) is refused at 8000 cells,
with the message "20 additive + 20 dominance + 0 pair columns", and the database and
`.Random.seed` are unchanged. Mutation check: always building `D` makes G7 fail.

**Please double-check:**

1. Your two 800-cell cases, rerun with your instrumentation on the 0.76.1 commit. I expect
   560 / 640.
2. That giving zero blocks zero coefficients without `.qtl_calibrate()` changes no seeded
   output that matters. Before, a zero target went through `.qtl_calibrate()` with
   rank 0. Seeded output may change across versions under the pre-1.0 policy, and
   within-version determinism (G2) still passes.
