# Review after QTL-effect methods Step 2 (0.73.1)

**Reviewed:** 2026-10-03. **Purpose:** decisions and fixes to consider before Step 3 and before the non-additive generator. This is a review, not an implementation record.

## Scope and assessment

I traced `plans/import_qtl_effect_methods.md` and the Step 2 results against the additive generator, congruence, target storage and readers, base-frequency selection, founder sampling, genome-effect writer, phenotype path and the new tests. I also checked the planned Step 3–5 contracts. The core congruence is a useful port: for well-conditioned feasible targets, the tests cover full and rank-one matrices, realised and genic anchors, atomic target/term writes, and repeatability. The base selection API makes the intended population explicit. The issues below matter because the interface presently promises more exactness and provenance than the implementation can establish in some cases.

**Phase count:** the main plan has **three** unfinished steps after Step 2: Step 3 (consolidation and phenotypes), Step 4 (variance extraction), and Step 5 (additive/dominance/A×A generator). It is not one remaining phase. Step 3 must land before the general generator, or phenotypes will omit non-additive genetic values.

### What “target this G in that population” currently means

| User intent | Current route | Exactness claim to make |
|---|---|---|
| HWE/LE limit at a chosen pool or cohort's frequencies | Filter `base_tbl`, use `anchor = "genic"`, `method = "shared"` | Exact under the **genic** covariance, subject to feasibility and finding 1; finite animals need not hit it. |
| The observed covariance in named diploid individuals | Filter an individual table as `base_tbl`, use `anchor = "realised"` | Exact in those individuals, including their LD, if their QTL genotypes are complete and the size limit permits it. |
| A line's own allele-copy effect variant | Set `line_name` and choose/accept its line base | Calibrates that variant. The full active model in a line or crossbred can also contain fallback terms, so its observed covariance must be measured. |
| Different per-trait QTL sets or a mixed parent-origin model | `method = "union"` or separate scoped calls | No general exact matrix guarantee; the covariance constraints depend on overlap and which inherited copies express each trait. |

This table should appear, in shorter form, in the Step 5 vignette. It distinguishes a target *parameter* from a realised population statistic without asking users to understand the internal term tables first.

## Findings, in priority order

### 1. A valid, poorly conditioned target is silently changed (scientific accuracy; fix before advertising “exact G”)

`R/qtl_congruence.R` `.qtl_psd_eigen()` treats an eigenvalue below `1e-10 × largest` as zero (lines 84–104), and `.qtl_congruence()` constructs `B` from only those retained eigenvectors (lines 202–238). `R/define_additive_effects.R` then labels the result “exact” (lines 437–455, 499–503), while `.tvc_write_block()` stores the **original** matrix (lines 335–369). For example, with `G = diag(c(1, 1e-11))`, identity `B0`, and `M = I`, the code reports rank 1 and delivers `diag(c(1, 0))` while storing `diag(c(1, 1e-11))`. I reproduced this with a direct call to the internal functions. The small eigenvalue is far above double precision at scale 1, so this is a real lost target, not rounding at machine precision. The same tolerance also accepts a slightly indefinite matrix, such as `matrix(c(1, 1, 1, 1 - 2e-11), 2)`, as a rank-one target while storing its negative eigenvalue. This tolerance was inherited from the source method; the integration issue is the unchecked “exact” claim against the original stored matrix.

**Recommendation:** separate PSD tolerance from numerical rank. If a positive eigenvalue is too small for stable calibration, error with its value and the effective rank, before writing the target; otherwise retain it and verify `B' M B` against the exact stored `G` under a stated mixed absolute/relative tolerance. Add an ill-conditioned positive-definite gate at the exported generator level. Apply the same policy to Step 5's three blocks and its additive floor.

### 2. Founder-pool “expectation” is biased by its divisor (scientific diagnostic)

`R/define_additive_effects.R` `.dae_diagnostics()` uses `2 * crossprod(Hc %*% B) / (nrow(H) - 1)` (line 1066). `R/add_founders.R` draws each of the two founder haplotypes **with replacement** (lines 207–215). Conditional on the stored pool, the expected variance of a new founder dosage is `2 * crossprod(Hc %*% B) / nrow(H)`. The present diagnostic is too large by `n_h/(n_h - 1)`. With the two-haplotype pool `{0, 1}` and effect 1, it reports 1; the actual draw distribution is `Binomial(2, 0.5)` with variance 0.5. This is most visible for exactly the small pools Q22 discusses.

**Recommendation:** use the population divisor `n_h` for the stated with-replacement pool expectation and add a tiny exact-enumeration test. Retain `n-1` for *observed* covariance of selected individuals and for a realised anchor. Update §7.4 and Q22's interpretation and example rates if those rates used the old divisor.

### 3. `G =` bypasses the stored non-additive-target check (target selection; fix in Step 3)

The default additive call is documented to error when the call's traits already have a dominance or A×A target, unless `trait_var_comp_tbl` explicitly filters those blocks out (`R/define_additive_effects.R` lines 63–79; main plan §6C). But `.dae_resolve_target()` returns immediately on a passed `G` (lines 819–837). Its `other <- ...` check is reached only when `G = NULL` (lines 865–874). I reproduced a call with a stored `dominance = 0.2` target and passed `G = 1`: it succeeded, wrote additive effects, and left both target rows. The dominance block was silently ignored. This matters as soon as users stage Step 5 targets in advance.

**Recommendation:** resolve the non-additive presence before the passed-`G` return. A fresh `G` should require an explicit choice to omit those blocks, with an API that does not collide with the current `G` + `trait_var_comp_tbl` refusal, or it should refuse and point to the general generator once it exists. Decide that API in Step 3, then add tests for passed `G` with stored `dominance` and `additive_by_additive` targets.

### 4. “Generated owner proves calibration” is not yet a safe premise for prevalence (Step 3 design blocker)

Main plan §6A/Q21 says any `generated` terms prove calibration to the stored target. The public `define_genome_effect_terms(..., effect_owner = "generated", allow_reserved_owner = TRUE)` can write arbitrary coefficients under that owner (`R/define_genome_effect_terms.R` lines 161–220). `remove_rows()` can remove a subset of generated terms or target rows after calibration. Neither path changes ownership. A future prevalence check based only on owner and target presence can therefore compute a threshold from a target unrelated to the active model. Step 3's proposed `define_effect_cov_matrix()` refusal only blocks one way to *write* a target; it does not protect these edits or a target already present. Also, `method = "union"` explicitly leaves nonzero off-diagonal targets approximate, so “calibrated” must not be taken to mean that every covariance in the stored block is attained.

**Recommendation:** revise §6A/Q21 before implementing PH7. Either make the reserved owner unforgeable through exported writers and guard `remove_rows()` on target/generated tables, or verify the active model against the target when deriving prevalence. The latter must use the stated anchor and handle scoped and non-additive models; if it cannot, require `thresholds =`. A practical first release can require explicit `thresholds =` for any model whose provenance cannot be proven. Add a regression where manual coefficients are written under `generated` with the override flag and a regression where one generated term is removed. Do not describe owner alone as proof while those public paths remain.

### 5. `method = "union"` can store a target for a trait receiving no effects (target contract)

When an existing per-trait QTL set is empty, `R/define_additive_effects.R` warns but continues (lines 382–395). It calibrates only nonempty masks (lines 446–453), skips the empty trait in the build (lines 474–482), yet writes the complete passed target block (lines 484–488). A positive target variance for that trait is then unattained, and `exact` can even be set to `TRUE` when all off-diagonals are zero (line 454). The completion message presents a calibration to a target that the active model does not contain.

**Recommendation:** with `scale_to_target = TRUE`, refuse an empty trait QTL set whenever its target diagonal is positive. Decide explicitly whether a zero target with no terms is acceptable; if so, message that it was not generated. Add a two-trait union gate with one absent per-trait set.

### 6. `parent_origin` silently truncates fractional inputs (wrong biological scope)

`.dae_parent_origin()` converts to `as.integer()` before checking membership in `{1, 2}` (`R/define_additive_effects.R` lines 546–577). I reproduced `.dae_parent_origin(1.9, "T")` returning parent 1. A typo can silently turn a biologically invalid value into a paternal effect.

**Recommendation:** validate finite, scalar/vector, exact integer values before coercion, including named entries. Add tests for `1.9`, `Inf`, duplicate names, and missing values.

### 7. An allowed projected pool can fail after the target and terms commit (atomicity and usability)

`extract_allele_freq()`'s shared `.validate_base_tbl()` accepts a `founder_haplotypes` selection exposing only `locus_name` and `allele` (`R/extract_allele_freq.R` lines 145–169), and the generator's `base_tbl` documentation says this is the accepted shape. But the default diagnostic later queries `line_name` and `haplotype_id` from that **projected subquery** (`R/define_additive_effects.R` lines 1046–1060). I reproduced this with `select(locus_name, allele)`: calibration and `.ge_commit()` finished first (line 488), then the diagnostic raised a DuckDB missing-column error. `trait_var_comp` contained one target row after the failed call. The caller sees an error even though the target and effects were committed. This is especially confusing when retrying, because the target now exists and the retry is refused.

**Recommendation:** collect or validate all columns the diagnostic needs before commit, or run a diagnostic directly on the physical pool keyed through the filtered selection. Keep the advertised base shape true. Test a projected founder selection, including the database state after a diagnostic failure.

### 8. A scientifically feasible cross-scope covariance is refused (scope limitation; document or extend)

`R/define_additive_effects.R` lines 332–347 refuse **any** difference in `parent_origin` across traits. The explanation is true for paternal-only versus maternal-only traits under random mating: their covariance is zero. It is not true for common-scope versus paternal-only. At one locus, `Cov(X_p + X_m, X_p) = pq`, while their variances are `2pq` and `pq`; a nonzero covariance is feasible, with a constrained maximum correlation. The current single `n_eligible` and common congruence cannot represent this mixed-scope anchor, so refusal is understandable, but the message and the plan's “any target G” framing are overbroad.

**Recommendation:** either specify that one call requires one identical parent scope and give the actual reason, or introduce a block anchor that represents shared paternal/maternal copies and checks feasibility for mixed scopes. This is a later extension; do not imply that all mixed-scope off-diagonals are biologically zero.

### 9. The seed guarantee is broader than the code can deliver (documentation and reproducibility)

The Step 2 note says “seed is applied after every check” and “No refusal touches `.Random.seed`.” In fact `set.seed(seed)` precedes `.qtl_congruence()` (`R/define_additive_effects.R` lines 418–443), where both anchor-rank and drawn-architecture-rank errors occur. A request for rank-two `G` on one segregating QTL can therefore fail **after** resetting and consuming the RNG. The latter error necessarily depends on the draw, but the anchor-rank error is knowable before it.

**Recommendation:** move the anchor-rank feasibility check before `set.seed()` and narrow the guarantee to “input and anchor feasibility errors do not touch RNG”; document that stochastic architecture failure can consume RNG. Add a seed-state gate for the former and a clear message for the latter.

## Changes to the remaining plan

1. **Before Step 3:** resolve findings 1–3 and 5–7 as Step 2 corrections, with gates. They change current results and messages, so record the version and tests in a small correction phase note. The Q22 fix should be included in the pool correction.
2. **Step 3:** replace Q21's owner-only proof with an enforceable contract (finding 4). `PH7` must test public mutation paths, not only ordinary generator calls. For a binary trait, a wrong liability variance directly changes prevalence, making this a scientific behavior gate. Keep the planned move to `ind_tgv_total`, the `formula_tgv` rename, exact component sums, and composite-phenotype checks; they are prerequisites for non-additive effects to reach observations.
3. **Step 4:** report the reference population and interpretation of each variance estimate prominently. §8's `anchor = "genic"` on a selected cohort is a projection at `base_tbl` frequencies, while `anchor = "realised"` measures the selected cohort. The default `base_tbl = NULL` intentionally drifts with the cohort and is not a check of the original generation target. Add a worked target-vs-measured example that filters `trait_var_comp` by `line_name` and calls out this distinction. For case 2, label the result “evaluated additive variance,” not a new-population breeding-value variance; §8 already says this but the output should make it hard to miss.
4. **Step 5:** apply the corrected target-rank rule to `G_A`, `G_D`, `G_AA`, and test non-additive effects in *phenotypes*, not only `ind_tgv` or covariance calculations. A passed block should not bypass a stored optional block without an explicit user decision. The planned `define_genome_effects()` has no `seed` argument while `define_additive_effects()` does; decide whether the first-release API should be consistent. The current plan explicitly chooses no seed (Q12), so this is a usability decision, not a correctness defect.
5. **Scope promise:** state precisely what is guaranteed. A common additive generator hits a feasible, numerically resolvable `G` under its named anchor and chosen QTL set. A line-specific call calibrates its own *variant*; common fallback terms, other origins, LD, finite sampling and selection can change the variance of the active model in a line or crossbred cohort. `method = "union"` generally does not hit a nonzero off-diagonal. The extractor should be the documented way to measure the target population the user actually means.

## Verification and limitations of this review

- Direct R probes reproduced the lost small eigenvalue and fractional `parent_origin` acceptance. The pool divisor follows exactly from the `sample(..., replace = TRUE)` call and the covariance of independent draws from the empirical pool.
- I ran the focused `qtl-congruence` and `define_additive_effects-anchor` tests independently during this review: both passed with process exit 0. These tests do not exercise the cases above.
- I ran the full package suite at 0.73.1 with `devtools::test(reporter = "summary")`: process exit 0, no failed tests or errors, 11 warnings shown by the reporter. The new findings remain uncovered by that suite.
- A package-level reproduction of the projected founder selection raised a missing-`line_name` binder error after one target row had been committed.
- A package-level reproduction of passed `G` with an existing dominance target succeeded and left both target rows, despite the default-call contract.
- The Step 2 note reports a full suite pass at 0.73.0; a passing suite does not cover the cases above. I did not modify package code or the main plan in this review.

**Focused and full-suite results:** passed. The full-suite warning list is in `/tmp/tidybreed_phase2_review_tests.log` for this run.

---

## Response (Claude, 2026-10-03, 0.73.2)

Every finding was checked against the code. Findings 1–3 and 5–9 hold and are fixed in
0.73.2 with gates. Finding 4 is partly right and becomes open question **Q23** in the main
plan. Plan changes are in the main plan's Step 2 "Corrections, 0.73.2" block and its Step 3–5
lists.

| # | Verdict | What changed |
|---|---|---|
| 1 | **Agree, and broader.** The cut-off depends on the traits' units (`diag(1, 1e-11)` is a full-rank target recorded in small units, not a rounding artefact). Fixing only the target exposed the same truncation on the architecture: `MVN(0, G)` gives the small trait a tiny `B0` column, so `rank(B0' M B0)` dropped it too | `.qtl_target_std()` decides rank and PSD on the correlation scale. `.qtl_calibrate()` normalises `B0` columns to unit anchor variance, runs the congruence, rescales, and **verifies** `B' M B` against the stored `G` at `1e-8` (correlation scale), erroring before any write. The same verified check sets "exact" under `"union"`. **Not changed:** `matrix(c(1, 1, 1, 1 - 2e-11), 2)` is still accepted as rank one. On the correlation scale its off-diagonal is `1 + 1e-11`, which is numerical precision. The verified delivery differs from the stored value by `1e-11 < 1e-8`, so "exact" holds to the stated tolerance. Gates R1 (unit level, unit invariance, generator) and R1b (a missed calibration writes nothing) |
| 2 | Agree | Divisor `n_h`. Gate R2: a two-haplotype pool `{0, 1}` gives spectrum `[1, 1]`, matching exact enumeration of `Binomial(2, ½)`. §7.4 updated. Q22's rates are annotated, not re-run (5% shift at 20 haplotypes; the decision does not depend on them) |
| 3 | Agree | Refused now, with the store-then-`trait_var_comp_tbl` route in the error; no new argument. Whether a passed `G` should get an explicit opt-out is listed for step 3/5. Gate R3 (both `dominance` and `additive_by_additive`, and the route works) |
| 4 | **Partly.** `remove_rows()` cannot remove generated terms: it refuses all three `genome_effect*` tables (pinned). `"union"` delivers the diagonal exactly, which is all a prevalence threshold reads. `allow_reserved_owner = TRUE` is a real forge path | Q23 decided (a), built in 0.73.2: the exported writer has no override. PH7 keeps the regression |
| 5 | Agree. Also, "exact" was wrong for a *zero* target covariance with overlapping QTL sets | Positive-variance trait without QTL: error. Zero-variance trait without QTL: message. Exactness computed. Existing test updated (it encoded the bug); gates R5, R5b |
| 6 | Agree | Validated before coercion, names checked. Tests in `test-genome-effects-writer.R` |
| 7 | Agree | Diagnostics computed before the commit, on the un-projected filtered pool; reported after. Gate R7 (projected base works; a diagnostic failure leaves target and terms unwritten) |
| 8 | Agree that the message was overbroad | The refusal stays, which is your recommendation's first option. The message and roxygen now give the real reason (one anchor per call) and say that common vs one-parent covariance is non-zero. A mixed-scope anchor is a later extension. Gate 44 extended |
| 9 | Agree | Anchor-rank check before `set.seed()`. `seed` roxygen names the draw-dependent failures. Gate R9 |

**Phase count:** agreed, three steps remain (3, 4, 5).
