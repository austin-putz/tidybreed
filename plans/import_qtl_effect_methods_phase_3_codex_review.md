# Review of import QTL-effect methods — Step 3

**Date:** 2026-10-05.  
**Reviewed version:** 0.74.4, commit `76aa0df`, including the 3a, 3b, 3b review, 3c and 3c follow-up commits since `be43687` (0.73.2).  
**Scope:** [Step 3 implementation plan](import_qtl_effect_methods_phase_3_plan.md), [implementation summary](import_qtl_effect_methods_phase_3.md), and the relevant contracts in [the main plan](import_qtl_effect_methods.md). Scientific correctness was the priority. Package implementation, tests and existing plans were left unchanged; this review is the only repository file added.

## Assessment

The central consolidation is sound: `ind_tgv` is the component table, consumers can read its total, and dominance, indicators and interactions now reach phenotypes. The independent mixed-model oracle is useful, and the change to preserve custom fields during re-evaluation is sensible. The component and formula validation also close real problems.

**The scientific contract around stored targets is still incomplete.** I reproduced public API paths where all terms have owner `generated`, but the variance used for a categorical threshold does not describe their genetic values. There is also a pre-existing omission of environmental variance from those thresholds. These are substantive errors in the variance supplied to the threshold; they are separate from finite sampling, selection, LD, or an accepted Gaussian approximation.

I would correct findings 1–4 before relying on Step 3's claim that generated ownership establishes a trustworthy liability variance. Findings 5 and the scientific qualifications below concern interpretation and the scope of that claim. Findings 6–9 are additional implementation and boundary-contract issues found during the review.

| Finding | Priority | Observed consequence | Relationship to Step 3 |
|---|---|---|---|
| 1. Explicit target block can belong to another scope | High | Common effects deliver genic variance 100; threshold reads target 1 | Existing target-selection branch bypasses the new protection |
| 2. Zero-target `union` trait retains old effects | High | Stored target becomes 0 while six old effects still deliver variance 1 | Existing replacement edge case violates the strengthened generated-owner contract |
| 3. Reciprocal parent scopes are counted as one variance | High | Each parent's terms deliver variance 1; total variance is 2; threshold uses 1 | New prevalence reader does not distinguish whole-model variance from scope calibration |
| 4. Random environmental variance is omitted | High | Requested prevalence 0.1 gives 0.3442 in a Gaussian example | Pre-existing threshold issue, still present in the changed path |
| 5. Structural `additive` can silently be used as a breeding value | Medium, interpretation | A genetically additive indicator surface produces an all-zero default true index | Chosen structural API, with a scientific safeguard removed |
| 6. Non-finite formula results are stored | Medium | `formula_tgv = "T / 0"` writes infinite phenotypes | Pre-existing numerical issue survives the DSL rewrite |
| 7. Constant genetic-value formula fails at evaluation | Low | `formula_tgv = "2"` is accepted, then fails on a multi-animal subset | Accepted grammar and evaluator disagree |
| 8. Failed phenotype overwrite can erase the old definition | Medium | `mean = NULL, overwrite = TRUE` errors after deletion | Pre-existing atomicity issue only partly addressed by 3a |
| 9. Threshold ties disagree with the strict-exceedance oracle | Medium, contract | Every liability exactly on a cutpoint enters the upper category | Pre-existing convention is untested by the new continuous-liability fixtures |

The severities concern effects on results, not whether the relevant branch was first introduced in this release. I checked the baseline diff to distinguish old defects from new behavior.

## 1. Explicit target selection can break the target-to-scope association — high

**Locations:** [target selection](../R/define_additive_effects.R#L892), [single-block check](../R/define_additive_effects.R#L940), [generated-term scope inference](../R/define_effect_cov_matrix.R#L419), and [prevalence variance lookup](../R/add_phenotype_stages.R#L1196).

`.dae_resolve_target()` validates that `trait_var_comp_tbl` comes from the same database, contains additive rows, and represents one complete covariance matrix. It checks that the rows have only one `line_name`, but does **not** check whether that line is the target scope associated with the effects being generated.

This permits the following entirely public sequence:

```r
pop <- define_trait(pop, "T")
pop <- define_effect_cov_matrix(pop, "additive", 1, trait_name = "T")
pop <- define_effect_cov_matrix(pop, "additive", 100,
                                trait_name = "T", line_name = "A")

line_target <- get_table(pop, "trait_var_comp") |>
  dplyr::filter(effect_name == "additive", line_name == "A")

pop <- get_table(pop, "genome_meta") |>
  define_additive_effects("T", trait_var_comp_tbl = line_target,
                          seed = 11, warn_bounds = NULL)
```

The call has `line_name = NULL`, so it writes **common** terms, with no origin rows. It calibrates them to 100, while the population-wide stored target remains 1.

**Measured in the review:**

- `sum(2 * p_j * (1 - p_j) * a_j^2)` from the written common coefficients: **100**.
- Population-wide additive target: **1**.
- `.ap_prevalence_genetic_var(pop, "T")`: **1**.
- Origin rows: **0**.
- A subsequent categorical phenotype with `prevalence = 0.1, residual_var = 1` succeeds; the particular 60-animal probe produced an affected fraction of **0.4**. The deterministic variance mismatch is the evidence; the observed fraction is illustrative.

There is a second consequence: after deleting line A's target, `define_effect_cov_matrix(..., line_name = "A")` accepts a replacement. The new protection thinks the common terms depend on the population-wide block. Their actual calibration source was line A's block, so the association inferred by `.tvc_generated_terms()` is wrong.

The same problem can occur when explicitly choosing population-wide rows despite an existing target for the effect's line, or choosing another line's rows. Ownership establishes that a calibration happened; it does not identify **which** stored target was used.

**Recommended correction:** validate the selected block's scope against the generator's target-resolution contract before any RNG use or write. If arbitrary source scopes are intentional, persist the actual calibration association and use it in every later reader and replacement safeguard. Inferring it from the origin predicates is insufficient once explicit selection can override it.

**Required gates:** common call with line rows; line call with another line's rows; line call explicitly choosing a population-wide block when its own block exists; valid same-scope selection; valid line-to-population fallback. Refusals should leave RNG, targets, terms and evaluated values unchanged.

## 2. A zero-target trait in `union` can keep old effects outside the candidate loci — high

**Locations:** [empty per-trait QTL selection](../R/define_additive_effects.R#L385) and [replacement loop](../R/define_additive_effects.R#L483).

For `method = "union"`, an empty per-trait QTL set is allowed when the target diagonal is zero. That is reasonable. However, the replacement loop executes `if (!any(rows)) next` **before resolving deletions for the trait's scope**.

Consequently, “no effects to insert” also becomes “do not delete the old effects.” The calibration and delivered-covariance calculation only examine the candidate mask, so the surviving effects are invisible to the call's verification.

**Public reproduction:**

1. Generate trait T on loci 1–6 with target variance 1.
2. Generate trait U on loci 7–12 with target variance 1.
3. Remove the old additive target rows.
4. Select loci 1–6 and call:

```r
G <- diag(c(1, 0))
dimnames(G) <- list(c("T", "U"), c("T", "U"))

pop <- get_table(pop, "genome_meta") |>
  dplyr::filter(locus_id <= 6L) |>
  define_additive_effects(c("T", "U"), G = G,
                          method = "union", seed = 13,
                          warn_bounds = NULL)
```

**Measured:** U has six surviving generated terms. Its genic variance is **1 before and after** the call, but its stored diagonal and prevalence variance are both **0**. The call succeeds.

This differs from the already acknowledged `union` limitation: an approximate off-diagonal is allowed, but the per-trait diagonals are supposed to be calibrated. Here the diagonal itself is wrong for the retained model.

**Recommended correction:** determine the scope deletions for every requested trait, including traits with no new terms, and commit those deletions with the target. If preserving old terms outside the candidate set is intended instead, the call must refuse a zero target that they contradict and must measure them in its delivery check. The existing replacement semantics favor clearing the scope.

A resulting trait with no terms also needs an explicit policy: either retain a validated zero-valued architecture or let the existing “no genome effects” error apply. It must not keep the old nonzero model merely to avoid that decision.

**Required gate:** start with nonzero effects outside the new candidate loci, then rerun with a zero diagonal. Assert on the stored terms and an independently calculated whole-scope variance, not only on the call's masked `delivered` matrix.

## 3. Two independently calibrated parent scopes do not make one calibrated trait — high

**Locations:** [eligible-copy variance convention](../R/define_additive_effects.R#L611), [scope construction](../R/define_additive_effects.R#L622), [parent fallback warning](../R/define_additive_effects.R#L712), and [one-diagonal-per-kind prevalence calculation](../R/add_phenotype_stages.R#L1196).

The generator correctly calibrates a one-parent additive model using one eligible copy per locus. But the same trait can have both paternal-only and maternal-only terms, each generated against the same stored target. The prevalence reader counts the additive kind once, regardless of how many contributing scopes were calibrated separately.

**Public reproduction:** store population-wide additive target 1, then generate the same trait twice:

```r
pop <- get_table(pop, "genome_meta") |>
  define_additive_effects("T", parent_origin = 1,
                          seed = 31, warn_bounds = NULL)
pop <- get_table(pop, "genome_meta") |>
  define_additive_effects("T", parent_origin = 2,
                          seed = 32, warn_bounds = NULL)
```

Both sets remain, by the documented scope replacement rules. Their parents are disjoint, so this does not trigger the nested-parent fallback warning.

At HWE and LE, with independent paternal and maternal gametes:

```text
V_paternal = sum_j p_j q_j a_paternal,j^2 = 1
V_maternal = sum_j p_j q_j a_maternal,j^2 = 1
V_total    = V_paternal + V_maternal      = 2
```

**Measured from the written terms:** the two sums are **1 and 1**, while `.ap_prevalence_genetic_var()` returns **1**. Neither generator call emitted a warning in this probe. Every term has the reserved owner, every individual calibration is correct, and the mismatch exists in the reference equilibrium population.

For a normal approximation with residual variance 1, the current threshold is based on variance 2 instead of the actual liability variance 3. This error precedes any questions about selection or LD.

**Recommended correction:** the prevalence eligibility check needs a whole-model scope contract. A conservative first fix is to refuse automatic prevalence for multi-scope/imprinted combinations that a single stored diagonal cannot describe and name `thresholds =`. A broader solution can represent and calibrate the combined model, including the covariance of its contributing copies. Counting or summing target rows mechanically is insufficient for overlapping fallback variants.

**Required gates:** reciprocal parent-only effects at one common reference; a common-plus-parent fallback; same target but different locus coverage; target replacement followed by completion of the documented scope regeneration sequence.

The 0.74.2 decision to **warn** when a target change strands another parent scope was explicitly accepted in the plan. This review is not treating that chosen warning policy as an accidental regression. The example above is stronger: it needs no target change and remains wrong even after both scopes have been calibrated to the current target. The documentation should distinguish calibration of a variant from calibration of the active trait.

## 4. The categorical threshold excludes random environmental variance — high, pre-existing

**Locations:** [liability assembly](../R/add_phenotype_stages.R#L627), [threshold calculation](../R/add_phenotype_stages.R#L1248), and [named random-effect model](../R/define_effect_random.R#L1).

The liability includes `random_contrib[[t]]`, but the threshold uses only genetic variance and the unconditional residual diagonal. Named random effects contribute neither their variance nor their distribution to that threshold.

A clean counterexample eliminates all uncertainty about genetic architecture or non-normality:

- Generated genetic variance: **0**.
- Residual variance: **1**.
- Normal permanent-environment effect, one independent level per `id_ind`: **9**.
- Intercept: **0**.
- Requested prevalence: **0.1**.

```r
pop <- define_phenotype(pop, "T", type = "categorical",
                        prevalence = 0.1, residual_var = 1,
                        store_liability = TRUE)
pop <- define_effect_random(pop, "T", "pe",
                            source_column = "id_ind", variance = 9)
pop <- get_table(pop, "ind_meta") |> add_phenotype("T", seed = 15)
```

The actual liability is exactly Gaussian with variance 10. The implemented cutoff is `qnorm(0.9)`, whereas the appropriate marginal cutoff is `qnorm(0.9) * sqrt(10)`.

```r
pnorm(qnorm(0.9) / sqrt(10), lower.tail = FALSE)
# 0.3426424
```

**Measured in 10,000 individuals:** affected fraction **0.3442**, liability sample variance **10.12696**. This agrees with the wrong cutoff's exact probability. It is not an unexplained sampling deviation from 0.1.

The normal threshold depends on the mean and variance of the complete liability distribution. This is the model underlying the normal-quantile cutoff in [Wray et al. (2010), Methods: disease traits](https://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1000864). The numerical counterexample above is an independent calculation using the package's declared random effects.

**Recommended correction:** include the variances of independent, centered normal named effects when automatic prevalence represents a marginal reference population. For a conditional interpretation, define the conditioning explicitly and use the appropriate conditional mean and variance. Require explicit thresholds for model combinations whose reference distribution cannot be derived, including non-normal random effects, rather than silently omitting them.

**Related residual issue:** a stored “unconditional” residual row is also the sampler's **fallback stratum**. If records actually select conditional strata with different variances, that fallback diagonal is not automatically the population's marginal residual variance. Calling `var_unconditional` “marginal” in the threshold comments overstates what the stored row establishes. A mixture requires reference stratum weights/distribution or explicit thresholds. This follows from the residual selection code; I did not run a separate mixture simulation.

**Required gate:** the zero-genetic-variance Gaussian example above, checking the deterministic cutoff as well as the resulting distribution. Add eligibility tests for normal, gamma/uniform and conditional-stratum models.

## 5. The structural component API can silently supply an incorrect breeding-value interpretation — medium

**Locations:** [breeding-value description and index default](../R/add_tgv.R#L27), [missing-component zero rule](../R/add_tgv.R#L297), and [default true-index computation](../R/add_tgv.R#L349). Step 3 deletes `.gev_warn_tbv_stale()`.

The implementation deliberately stores declared model structure. Thus `component_name = "additive"` is not a statistical additive projection of arbitrary stored genotypic values. This distinction is valid and documented in parts of the API. The risk is the default true index and wording that equate the structural component with a breeding value without a model check.

**Public probe:** define an indicator surface with values `0, 1, 2` for dosages `0, 1, 2` at one autosomal locus. Define a one-trait index with weight 1, then call `add_tgv("T", index_names = "I")`.

- The genotypic values vary over **0, 1, 2**.
- Every default true-index value is **0**.

The surface is genetically additive: its average effect is 1 and its centered breeding value is `dosage - 2p`. It was merely entered through `genotype_terms()`, so there is no row declared `additive`. The zero-fill rule loses all selection information if that result is interpreted as a breeding-value index.

This is **not a failure of the documented structural matrix multiplication**. It is an interpretation hazard made more important by removing the old warning and moving the breeding-value API onto a function that evaluates arbitrary models. Equivalent functional representations can now produce different default indices.

**Recommended safeguard:** retain the structural component selector, but warn or refuse when a default breeding-value index is requested for a model whose structural additive rows are not known to be its statistical breeding values. Provide a clear route to a projection when available. `component_name = "total"` is an explicitly different selection objective and should not be presented as a general substitute for a breeding-value projection.

A single common statistical base with the required equilibrium assumptions makes the familiar partition meaningful; arbitrary functional surfaces, mismatched centers, interactions and LD require care. [Vitezica et al. (2017)](https://pmc.ncbi.nlm.nih.gov/articles/PMC5500131/) explains the population dependence and orthogonality requirements. The recommendation here concerns safe interpretation of the existing API, not adding an unplanned projection implementation during Step 3.

**Suggested gates:** an indicator surface equivalent to a purely additive function; functional A+D with `p != 0.5`; a generated additive model supplemented by a functional surface. Verify either a correct projection or an explicit diagnostic, rather than treating the structural output as the scientific oracle.

## 6. The formula DSL accepts and stores infinite phenotypes — medium, pre-existing

**Locations:** [formula evaluator](../R/formula_helpers.R#L420) and [formula-path exclusions](../R/add_phenotype_stages.R#L455).

`formula_tgv = "T / 0"` passes the closed grammar, which is expected for an arithmetic expression. However, evaluation does not check `is.finite()`, and exclusions only check `is.na()`.

**Measured:** all **60 of 60** continuous phenotype rows in the probe were written with non-finite `pheno_value`. This used a correctly generated additive model and residual variance 0. The derived-phenotype evaluator already has an explicit `Inf`/`NaN` policy; the genetic-value DSL lacks an equivalent one.

Whitelisting functions prevents arbitrary execution but does not make their numeric results valid. Overflow in `exp(T)` or division by a zero component can reach the same path. On categorical traits, an infinite liability can instead silently force an extreme category.

**Recommended correction:** validate evaluated genetic values before residual/random draws and record writes. Apply a documented error/skip policy consistently to `NA`, `NaN` and infinities, with examples of affected individuals. Keep the D7 record-table contract.

**Required gates:** division by zero, overflow, and a domain error such as a logarithm of a negative value; test both continuous and categorical outputs.

## 7. A constant `formula_tgv` passes definition but fails on multiple individuals — low

**Location:** [formula result naming](../R/formula_helpers.R#L424).

The grammar accepts numeric constants, and `define_phenotype(..., formula_tgv = "2")` succeeds. The evaluator receives a scalar, then assigns all individual names with `setNames()`, which errors instead of broadcasting it.

**Measured:** a 60-animal call fails with:

```text
'names' attribute [60] must be the same length as the vector [1]
```

**Recommended correction:** either broadcast a scalar to the selected individuals, or require at least one trait reference at definition time and document that restriction. The accepted grammar should correspond to what can be evaluated. The existing grammar and constant-derived-formula support favor broadcasting.

**Required gate:** a constant genetic-value expression over more than one individual, plus a formula containing both a trait and a numeric offset.

## 8. The phenotype-definition atomicity fix is incomplete — medium, pre-existing

**Locations:** [overwrite deletion](../R/define_phenotype.R#L519) and [late mean validation](../R/define_phenotype.R#L542).

Moving component validation before writes fixes the reported bad-component cases. Other validation still occurs after deletion of the existing definition, with no transaction around the definition replacement.

**Public probe:** create phenotype T, then call:

```r
define_phenotype(pop, "T", mean = NULL, overwrite = TRUE)
```

The call errors, but the old `phenotype_meta` row has already been deleted. **Measured:** zero phenotype-definition rows remain. Any old components are deleted too, while residual configuration can remain behind.

**Recommended correction:** validate all inputs before mutation and wrap deletion, metadata insertion, component insertion and residual changes in one transaction where necessary. Coordinate transaction ownership with the residual writer rather than nesting exported transactional calls. The summary should describe the actual validation guarantee until this is corrected; fixing component validation does not establish atomicity for every refusal.

**Required gate:** snapshot an existing definition, its components and residual rows; try invalid mean inputs and a controlled insertion failure; verify exact preservation.

## 9. Exact threshold ties use a different convention from the scientific oracle — medium

**Locations:** [categorical conversion](../R/phenotype_helpers.R#L43), the `prevalence` description in [define_phenotype()](../R/define_phenotype.R#L39), and the PH7 end-to-end test in `test-phenotype-total-genetic-value.R`.

The prevalence description refers to the fraction **above** a threshold, and PH7's oracle uses `liability_value > thr`. `liability_to_categorical()` instead uses the default `findInterval()` convention, which places a liability **equal** to a cutpoint in the upper category.

With continuous, nondegenerate residuals, exact ties have probability zero and the test cannot expose this difference. Discrete genotype surfaces with residual variance 0 are supported, so ties can be a substantial part of the population.

**Public probe:** the same dosage surface `0, 1, 2` as finding 5, explicit threshold 1, and residual variance 0. Of 60 individuals, **33** had liability exactly 1; all 33 entered category 2. The upper-category fraction was **0.8333333**, while the fraction strictly above the cutoff was **0.2833333**.

Either boundary convention can define an ordered categorical model. The defect is the discrepancy between the documented/tested strict-exceedance rule and the implemented inclusive rule. This is independent of whether a Gaussian approximation is suitable.

**Recommended correction:** choose and document the boundary convention, align the classifier and scientific oracle, and explicitly test exact ties. If retaining the current inclusive convention, qualify the “above” wording and update the oracle to reflect it. Also define the behavior of a prevalence request on a liability with zero total variance; no cutoff can create an arbitrary fraction in a point-mass distribution.

**Related input-validation gap, found by inspection:** `define_phenotype()` checks that cutpoints are present, but does not validate a finite, non-missing numeric vector. For example, `thresholds = c(1, NA_real_)` declares three categories; the later `sort(thresholds)` drops the missing value, leaving only two intervals. Category values/names can then describe a different model from the classifier. Reject invalid cutpoints before writing the definition and add a partial-`NA` regression gate. This case was not separately simulated during the review.

## Scientific qualifications that should be explicit

### Normality is an assumption even at the correct reference population

A calibrated covariance does not establish a normal liability. The current formula

```text
threshold = mean + qnorm(1 - prevalence) * sqrt(V_genetic + V_residual)
```

is a Gaussian approximation unless the full liability distribution is normal. The documentation emphasizes finite samples, selection and line-specific populations, but a few large QTL can defeat the requested prevalence in an infinite HWE/LE reference population too. The normality premise is explicit in [Wray et al. (2010)](https://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1000864).

An exact enumeration illustrates the boundary:

- One additive locus at `p = 0.5`.
- Additive coefficient `sqrt(2)`, giving exactly `V_A = 1`.
- HWE genotype probabilities `0.25, 0.5, 0.25`.
- Normal residual variance `0.01`.
- Requested prevalence `0.1`.

Summing the three genotype-specific normal tail probabilities at the implemented cutoff gives **0.224163**, not 0.1. This is an independent mathematical calculation, not a package Monte Carlo estimate.

Document `prevalence` as a reference Gaussian approximation, including this condition. Exact cutpoints can be derived from a fixed named reference distribution and passed with `thresholds =`; the plan's rejection of a new empirical quantile for every batch is reasonable and need not change. This qualification does not excuse findings 1–4, where the variance input itself is wrong.

### Sum of component variances needs orthogonality

For total `A + D + AA`, the general variance is

```text
V_total = V_A + V_D + V_AA
          + 2 Cov(A, D) + 2 Cov(A, AA) + 2 Cov(D, AA).
```

Cockerham contrasts at a consistent HWE/LE base support an orthogonal reference partition. LD and departures from the coding reference can introduce cross-component covariance; [Vitezica et al. (2017)](https://pmc.ncbi.nlm.nih.gov/articles/PMC5500131/) discusses these limits. Step 3's sum of stored diagonals is therefore a reference-model rule, not a measurement of total variance in arbitrary evaluated individuals.

The main plan already recognizes cross-component covariance for Step 4. Keep that qualification consistent in Step 3's roxygen and summary. Before Step 5 supplies non-additive generators, verify that a generated A+D+AA model's target reference really satisfies the assumptions needed by the prevalence rule. An arbitrary realised anchor should not inherit an orthogonal interpretation solely from the owner label.

### Deterministic accumulation and numerical accuracy are different guarantees

The ordered floating sum for up to four components is a reasonable way to preserve the one-component value and remove dependence on thread scheduling. It is deterministic, but does not make arbitrary floating-point cancellation mathematically exact. The existing decimal term/group accumulator also has its documented finite domain and rounding floor.

I found no new arithmetic defect in the ordinary-scale component sums. The switch should be described as deterministic summation, and the scientific tolerances should remain explicit. The implementation summary already distinguishes parity agreement to tolerance from literal equality with the old evaluator.

### Cached true indices can become stale

`overwrite_index = FALSE` intentionally preserves existing index rows after genetic values have been recomputed. The summary identifies this pre-existing behavior, and the API documents it. It is not a hidden new bug. Workflows that regenerate effects or change index weights must explicitly refresh the true indices; an index row is a cached calculation, not a live view of current `ind_tgv`.

### The parent-fallback warning suggests recovery operations that are refused

`.dae_warn_parent_only()` recommends `mode = "replace_owner"` through `define_genome_effect_terms()`, or removal of the unwanted variant. For the generated variants that trigger this warning, the exported writer refuses owner `generated`, and `remove_rows()` refuses the three genome-effect storage tables. Those advertised operations therefore do not remove the generated variant. This is a code-inspection finding, not an additional runtime probe. Update the warning to give a supported recovery route, or supply a controlled generator-scope reset operation. The reserved-owner safeguard should stay intact.

## What the implementation does well

- **One evaluator and one component table.** The consolidation removes the old disagreement about which owners contributed to a genetic value. Re-evaluation deletes obsolete component rows for the selected individual/trait combinations and upserts surviving ones.
- **Phenotype routing.** Simple, composite and DSL paths read totals by default, with explicit component selection where requested. The PH1/PH2/PH3 gates cover meaningful cases, including parents and group contributors.
- **Independent arithmetic check.** `cons_oracle()` in `test-tgv-consolidation.R` computes common-scope additive, dominance, interaction and indicator terms directly from allele copies. This is substantially stronger than checking one reader against another reader of the same result.
- **Mean semantics.** Keeping `mean` as an intercept is correct for the raw-sum storage contract. A functional surface's offset is not silently erased. PH4 verifies the record assembly identity.
- **Reproducibility.** Ordered component sums and the preserved deterministic group accumulator support the thread-count contract. The existing larger-group tests are appropriate for exposing parallel reduction differences.
- **Definition-time DSL checks.** The named-only optional arguments, per-reference placeholders, trait/table/column checks and closed call grammar are coherent. The rewrite removes the arbitrary execution path described in the implementation summary.
- **Target protection on the ordinary path.** The new refusal catches direct target rewrites under generated terms and the previously missed fallback-line case. Removing manual/unscaled generator arguments closes the original uncalibrated-owner path.
- **Schema and compatibility.** The component vocabulary, true-index component dimension and explicit old-database refusals match the planned breaking release. Custom-column preservation improves repeated phenotype generation.

## Test coverage gaps

The new tests give good confidence in storage and reader routing. They do not establish the strongest scientific claims in the summary:

1. **PH7's non-additive fixtures deliberately bypass calibration.** They plant coefficients under the reserved owner through an internal test writer. The A+D+AA end-to-end test verifies that classification uses the intended sum of target rows; it does not establish that those target rows equal the variance of its constructed liability. This is a legitimate routing test, but cannot serve as a calibration/prevalence oracle.
2. **No explicit target-scope mismatch gate.** Testing “one selected block” is weaker than checking that its scope is the one downstream readers associate with the generated terms.
3. **Zero-target coverage does not start with effects outside the candidates.** The deletion bug requires that initial state; a fresh empty trait cannot expose it.
4. **No whole-model reciprocal-parent prevalence gate.** Individual scoped calibrations can all pass while their combination has another variance.
5. **No full-liability environmental-variance gate.** Normal random effects make a deterministic analytic cutoff available, without requiring a stochastic tolerance to decide correctness.
6. **The grammar checks do not cover scalar/non-finite result contracts.** Numeric validity and vector shape require evaluation tests in addition to AST tests.

For the correction phase, use independent enumeration or variance calculations alongside the existing equality checks. Test public API state transitions. Keep internal reserved-owner fixtures for routing tests, but label them as such.

## Proposed correction plan

This is a recommendation for a follow-up phase, not an edit to the implementation or the accepted main plan.

1. **Repair target association and replacement before expanding non-additive generation.** Resolve finding 1, then finding 2. Add gates that compare the stored target with every retained effect in the relevant scope, including empty new QTL sets.
2. **Define automatic-prevalence eligibility for the complete genetic model.** Resolve finding 3 by a conservative refusal or a represented combined-scope calibration. Preserve the user-approved parent-scope warning policy while ensuring the phenotype reader cannot interpret an unvalidated combination as one calibrated diagonal. Give its warning a supported recovery route.
3. **Include the full environmental model or require explicit cutpoints.** Resolve finding 4 and define the residual-stratum interpretation. Validate these checks in PLAN before any random draws or record writes.
4. **Restore scientific diagnostics for default true indices.** Address finding 5 without changing the meaning of explicit structural component selection. Qualify breeding-value and prevalence claims by coding/reference/model assumptions.
5. **Finish numerical, boundary and atomicity contracts.** Correct findings 6–9. Ensure failed calls preserve the documented tables, accepted scalar formulas evaluate consistently, and threshold ties follow the declared convention.
6. **Strengthen later scientific gates.** Step 4 should retain all cross-component covariance terms and identify its population. Step 5 needs an actually calibrated A+D+AA model followed through to phenotypes, plus equilibrium and non-equilibrium examples. Reuse the independent sum-plus-cross-covariance identity from the main plan.
7. **Update the phase summary after corrections.** Replace the unrestricted implication “generated proves the target describes the active model” with the guarantee that can actually be checked. Distinguish scope calibration, whole-model variance, covariance delivery and the Gaussian prevalence approximation.

## Verification and limitations

The small probes used seed 3401, 12 autosomal loci on one chromosome, 100 uniform founder haplotypes and 30 male plus 30 female founders in line A. The environmental-variance probe used 5,000 founders of each sex. Each probe had a fresh in-memory population. This public setup reproduces the common fixture without the package's test helpers:

```r
set.seed(3401)
pop <- open_pop(pop_name = "review", db_name = ":memory:") |>
  define_genome(n_loci = 12, n_chr = 1, chr_len_Mb = 100) |>
  define_founder_haplotypes(n_haplotypes = 100) |>
  get_table("founder_haplotypes") |>
  add_founders(n_males = 30, n_females = 30, line_name = "A")
```

- Reviewed the five implementation/review commits from 0.73.2 through 0.74.4, the implementation plans and summary, relevant consumers, schema/restore changes and new scientific gates.
- Loaded 0.74.4 with `pkgload::load_all(..., recompile = FALSE)` under R 4.5.3, testthat 3.3.2 and DuckDB 1.5.5. Probe databases were in memory; no package source or test files were edited.
- Ran public API probes for findings 1–9. Only the diagnostic inspection used internal readers/DBI; the states producing the findings were created through public model APIs. The exact one-locus Gaussian-mixture calculation was evaluated separately.
- Probe scripts and logs are temporary review artifacts: `/tmp/tidybreed_phase3_review_probes.R`, `/tmp/tidybreed_phase3_review_probes.log`, `/tmp/tidybreed_phase3_review_more.R`, `/tmp/tidybreed_phase3_review_more.log`, `/tmp/tidybreed_phase3_review_ties.R` and `/tmp/tidybreed_phase3_review_ties.log`.
- **Full suite passed:** `NOT_CRAN=true` with `devtools::test(reporter = "summary", stop_on_failure = TRUE, recompile = FALSE)` completed with process exit 0, no failed tests or errors, and no reported skips. The six warnings match those documented for 0.74.4: named-pool `hap_id`, repeatability skips in two tests, a monomorphic pool, and the two pooled-base notices in parity. The log is `/tmp/tidybreed_phase3_review_tests.log`. The reproduced findings are not covered by that passing suite.
- This review does not certify the unimplemented Step 4 extractor or Step 5 generator. The findings concerning their integration identify requirements for later gates.
- The existing modification to `plans/consolidate_genetic_values.md` was present before the review and was left untouched.
