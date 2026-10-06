# Review of import QTL-effect methods — Phase 4 implementation

**Date:** 2026-10-05.

**Reviewed version:** tidybreed 0.75.1, commit `3cde9d1`.

**Implementation commits:** `aa76666` (4a) and `3cde9d1` (4b), compared with `7c6e2b7`.

**Scope:** all changes belonging to [the Phase 4 implementation](import_qtl_effect_methods_phase_4.md), including the builders, conversion helpers, extractor, shared dosage collector and its additive-generator caller, tests/oracle, benchmark, exports, generated documentation, and related plan/API documentation. The existing evaluator, writer, inheritance resolver, allele-frequency extractor and cohort resolver were inspected as integration dependencies.

**This is an implementation review for Claude. No package code was changed.** The earlier pre-implementation plan review is preserved below under “Historical plan review”; its findings describe the old baseline and should not be mistaken for outstanding implementation defects.

## Assessment

The implemented algebra largely follows the revised plan correctly. The genic formulas, inverse contrast mappings, constant sign conventions, monomorphic-partner main effects, family-level scope handling, and ordered cross-component covariance accounting are sound under the stated component definitions. The focused existing checks passed: **446 expectations, 0 failures, 0 warnings, 0 skips**.

I found three implementation/validation issues and one scientific interpretation issue. The latter is inherited from the source/plan rather than a failure to port their algebra: under LD, the reported `additive` block is not necessarily the cohort's least-squares additive projection. That distinction needs to be explicit before these measurements are interpreted as breeding-value variance or used to justify the later breeding-value export.

| # | Priority | Finding | Provenance |
|---|---|---|---|
| 1 | High scientific concern | Under LD, “additive projection” overstates what the marginal-contrast decomposition measures | Source/plan interpretation, carried into new documentation |
| 2 | Medium | Roundoff changes block availability for equivalent one-locus surfaces | New canonicalisation/support logic |
| 3 | Medium | Pair loci absent from `a`/`d` produce `NA` coefficients; an empty-main-effect model then silently loses induced main effects | New forward conversion and writer-frame helper |
| 4 | Medium | Fractional copy counts are silently truncated into different genotype states | Retained defect in the modified `genotype_terms()` builder; already present in `7c6e2b7` |

## 1. Qualify the realised additive block's meaning under LD

**Locations:** `R/extract_genetic_variance.R:44–55`, `.egv_alpha()` and `.egv_realised()` at `R/extract_genetic_variance.R:479–483`; related breeding-value claims in the main plan.

**Classification:** scientific interpretation/documentation concern, not a disagreement with the implemented accounting identity or the requested source port.

The implementation computes

```text
alpha_j = a_j + b_j d_j + sum_l e_jl (2 p_l - 1)
A = centred_dosages %*% alpha
```

Here `b_j` only orthogonalises heterozygosity against dosage **at the same locus**. It does not project the complete genetic value onto the joint dosage space. Under LD, dominance at another locus or the centred pair products can still have an additive regression on those dosages. Adding their cross-block covariances back into `total` proves exact reconstruction, but does not turn `A` into that joint additive projection.

**Reproduced counterexample with HWE at both loci.** Use one functional pair with `e = 1`, no main effects, and these 32 observed two-locus genotypes:

| Dosage at L1 / L2 | 0 | 1 | 2 |
|---|---:|---:|---:|
| 0 | 3 | 2 | 3 |
| 1 | 3 | 9 | 4 |
| 2 | 2 | 5 | 1 |

Each locus has counts `(8, 16, 8)`, hence exact HWE margins and `p = 0.5`. The joint distribution is not LE. The independent regression is reproducible without a database:

```r
cells <- expand.grid(L1 = 0:2, L2 = 0:2)
X <- cells[rep(1:9, c(3, 3, 2, 2, 9, 5, 3, 4, 1)), ]
g <- (X$L1 - 1) * (X$L2 - 1)
var(fitted(lm(g ~ L1 + L2, data = X)))  # 0.02099937
```

After assigning those dosages to the in-memory population and writing the functional pair through `aa_terms()`, the actual extractor reports:

```text
additive               0
additive_by_additive    0.28931452
between_components     0
total                  0.28931452
decomposition          full
```

However, `lm(g ~ L1 + L2)` gives dosage slopes approximately `(-0.07450980, -0.19215686)` and sample variance of fitted values **0.02099937**. Thus, even zero `between_components` does not certify that the additive block equals the cohort's additive regression variance. This is an interpretive counterexample, not a request to replace the existing source oracle with `lm()` for every block.

The literature distinguishes the marginal NOIA construction under LE from a multilocus regression under LD. The NOIA extension described by [Álvarez-Castro and Yang (2011)](https://pmc.ncbi.nlm.nih.gov/articles/PMC3247674/) assumes LE for its multilocus orthogonality. [Álvarez-Castro and Crujeiras (2019), “Mean and Additive Component” and “Regression Procedures”](https://www.frontiersin.org/journals/genetics/articles/10.3389/fgene.2019.00054/full) explains why marginal regressions do not suffice for an orthogonal decomposition under LD. The numerical comparison above is my independent inference/check of this implementation against the joint additive regression definition.

**Recommendation for Claude:** preserve the planned marginal-contrast accounting if that is the intended estimator, but explicitly describe the realised rows as covariances of those contrast components. State that `full` means all stored term shapes are supported, not that the cohort's additive breeding-value projection has been recovered. Qualify “additive projection” for both full and partial models under LD, and carry that qualification into the later breeding-value plan. If the intended scientific quantity is instead the cohort's joint regression breeding value, that requires a separately specified projection; agreement with `nonadd_decompose()` alone cannot validate that interpretation.

**Regression gate:** retain the above HWE-margin/LD panel as an independent interpretation fixture. Test the documented distinction explicitly, rather than requiring the current estimator to match a different estimator without changing its contract.

## 2. Exact-zero support checks make equivalent surfaces report different blocks

**Locations:** `R/genome_effect_terms_builders.R:491–535` (`.stored_to_functional()` accumulation), `R/extract_genetic_variance.R:485–486` and `578–579` (availability).

The converter incrementally accumulates `d`, then both anchors decide dominance support using the exact test `co$D != 0`. Cancellation of decimal coefficients can leave a floating-point residue. This affects **which rows exist**, not just their last numerical digits.

Using the existing extractor test fixtures, I wrote:

```r
genotype_terms(data.frame(L1 = 0:2), c(0.3, 0.2, 0.1))
genotype_terms(data.frame(L1 = c(0, 2, 1)), c(0.3, 0.1, 0.2))
ad_terms("L1", a = -0.1, d = 0, p = 0.5, report = FALSE)
```

These specify the same linear one-locus surface, up to the last builder's irrelevant constant. The first conversion yields `d = 1.387779e-17` and reports a dominance row (observed variance approximately `4.24e-35`). The reordered surface yields `d = 0` and omits that row, as does `ad_terms()`. I reproduced this through the real writer and extractor; reversing all surface rows also changes the small residual.

This violates the new report's structural coding-invariance promise. A dominance target can join a spurious “measured” row for one surface representation while falling into the `anti_join()` for another. Existing B14/B16 use values whose cancellations happen to be exact, so they miss this case.

**Recommendation:** make canonical support detection robust to accumulation error. Use a local error bound based on the contributions being cancelled, or another justified canonicalisation rule. Avoid a fixed absolute cutoff that would erase a legitimately small uncancelled effect merely because of the trait's units. Make support decisions consistently for `d` and pair coefficients and for both anchors.

**Regression gate:** the three representations above, including multiple surface-row permutations, should have identical block availability under both anchors. Include a small but genuinely nonzero coefficient to ensure the chosen rule preserves it.

## 3. The forward converter can silently drop pair-induced additive effects

**Locations:** `R/genome_effect_terms_builders.R:574–595` (`.noia_to_stored()`) and `606–612` (`.noia_terms()`).

The frequency validator explicitly considers the union of main-effect and pair loci, but `alpha`, `d`, and returned `p` are initialised only on `names(a)`. If a pair names a locus absent from `a`/`d`, `alpha[k] <- alpha[k] + ...` starts with `NA`. The function returns malformed coefficient data instead of expanding the zero main effects or rejecting the input.

The pair-only case is worse than an error:

```r
z <- setNames(numeric(), character())
s <- .noia_to_stored(
  a = z, d = z,
  pairs = data.frame(locus_1 = "L1", locus_2 = "L2", e = 1),
  p = c(L1 = 0.3, L2 = 0.6)
)
# s$alpha: L1 = NA, L2 = NA
# s$d and s$p: named numeric(0)
tt <- .noia_terms(s)
# Returns only the Cockerham pair; no additive main effects.
```

In `.noia_terms()`, `stat$alpha != 0 | stat$d != 0` becomes `logical(0)` because `d` is empty. The additive branch is skipped, and the pair is written successfully. The correct induced additive coefficients are **L1 = 0.2, L2 = -0.4**, with `mu = -0.08`. In the real evaluator probe, the functional model minus the returned statistical model ranged from **-0.32 to 0.68**, rather than being the constant `-0.08`.

**Impact:** this internal helper is a planned foundation for step 5. Existing N2 always provides explicit main-effect entries for every pair locus, so its round trip does not exercise this input. The current public extractor uses the inverse helper and is not directly affected by this forward-conversion defect.

**Recommendation:** initialise aligned main-effect vectors on the full locus union, treating omitted main coefficients as zero if sparse input is supported. Alternatively, enforce and document that every pair locus must be explicitly present in both vectors, and reject missing keys before performing arithmetic. In either policy, `.noia_terms()` should reject missing coefficients or inconsistent vector dimensions rather than silently emitting an incomplete model.

**Regression gate:** a pair-only functional model at unequal, non-0.5 frequencies; a pair with only one endpoint in the main vectors; and permuted named vectors. Assert that converted and functional evaluator values differ by `mu` alone, or that unsupported input is rejected at the converter boundary.

## 4. `genotype_terms()` truncates fractional copy counts before the writer can reject them

**Locations:** `R/genome_effect_terms_builders.R:212–215` and `232`.

**Provenance:** confirmed present in the baseline builder as well; this is a retained bug in code refactored by phase 4, not a newly introduced regression.

```r
tt <- genotype_terms(
  data.frame(L1 = 1), value = 1,
  copy_count = c(L1 = 2.9)
)
# tt$copy_count_value is 2L, without an error or warning.
```

`.ad_recycle()` checks finite numeric input, but not integrality or non-negativity. The builder calls `as.integer()` before handing the rows to `define_genome_effect_terms()`. Consequently, the writer's whole-number validation sees `2`, accepts the model, and the extractor classifies it as a fully supported diploid heterozygote term. A supplied copy count of `-0.5` similarly becomes `0L`; with dosage 0 it can become a valid copy-absence indicator.

**Recommendation:** validate raw copy-count values as non-negative whole numbers before integer conversion, before dropping zero-valued rows. Retain missing/inferred counts only where the builder's existing contract permits them. The writer cannot recover the original invalid value after coercion.

**Regression gate:** refuse `2.9`, `-0.5`, and malformed copy counts on a row whose coefficient is zero; continue to accept genuine `0`, `1`, and `2` states where the inheritance model supports them.

## Checks that held up

- The stored additive, dominance, all three diploid indicator, and A×A inverse identities have the correct coefficients and constant signs. `kappa = -mu` for a matching forward/inverse conversion is correct.
- The genic weights `2pq`, `(2pq)^2`, and `4p_kq_kp_lq_l`, with the induced additive coefficients, are correct for the stated HWE+LE reference. They measure that expectation, not the observed cohort variance.
- Fixed pair partners retain their induced main effect; a fixed heterozygote and fixed homozygotes are handled differently as required.
- The realised decomposition uses sample covariance with divisor `n - 1`. The divisor used in the within-locus regression cancels appropriately.
- Cross-component accounting includes D–A×A, both trait orientations, and uncovered values. Family classification precedes owner pooling, preserving scoped fallback competition.
- Default genic frequencies use the resolved individuals' whole genotypes, while an explicit copy-filtered base keeps allele-copy semantics. Genic extraction intentionally does not require an evaluated genetic value for each selected individual; that follows the revised plan and is not a missing-value defect.
- The shared dosage collector preserves the additive generator's completeness and size guards. Its SQL-order row labels correctly avoid the previous R/SQL collation mismatch.
- Pair chunking avoids unbounded pair-value/anchor matrices. The reported benchmark honestly uses 19,900 pairs rather than the planned 124,750. The existing writer bottleneck remains a step-5 scalability limitation, not an additional extractor correctness finding.

## Verification and limits

Ran the existing focused suites with `NOT_CRAN=true` and failure stopping enabled:

```r
testthat::test_local(
  ".",
  filter = "^(extract_genetic_variance|genome-effect-terms-builders|define_additive_effects-anchor|genome-effects-determinism|extract_allele_freq)$",
  stop_on_failure = TRUE
)
```

Result: **446 passed; 0 failed, warned, or skipped; exit status 0**. Also ran independent in-memory probes through the actual writer/evaluator/extractor for the findings above, and inspected the source project's measurement functions directly. The LD counterexample uses every two-locus dosage combination and exact single-locus HWE margins, so it does not depend on a missing genotype state or HWE departure.

I did not rerun the full package suite, the large benchmark, documentation generation, or mutation tests. The implementation results' reported full-suite/mutation/benchmark outcomes remain author-reported evidence; the focused results and counterexamples above were independently run for this review.

Temporary reproducibility artifacts for this session:

- `/tmp/tidybreed_phase4_implementation_probes.R`
- `/tmp/tidybreed_phase4_implementation_probes.log`
- `/tmp/tidybreed_phase4_implementation_tests.log`

**Claude follow-up:** resolve the scientific estimator wording/contract first, then add the cancellation and sparse-conversion gates before step 5 builds on these helpers. Fix the retained copy-count validation bug in the builder. No implementation changes have been made by this review.

---

## Historical plan review

# Review of import QTL-effect methods — Phase 4 plan

**Date:** 2026-10-05.  
**Reviewed baseline:** tidybreed 0.74.5, commit `7c6e2b7`.  
**Scope:** [Phase 4 plan](import_qtl_effect_methods_phase_4_plan.md), the mapping, measurement contract, gates and relevant questions in [the main plan](import_qtl_effect_methods.md), the current implementation and tests, and the source project's `nonadd_covariates()` / `nonadd_decompose()`. This reviews a proposed implementation, not completed Phase 4 code. Package code and existing plans were left unchanged; this review is the only repository file added.

## Assessment

**I agree with the overall direction and the 4a/4b split, but would revise several contracts before implementing them.** The extractor should use the existing evaluator for totals, canonicalise supported common-scope models, report every ordered cross-component covariance, and retain scoped families together. Those choices fit the codebase and the scientific accounting.

The main corrections are the rule for deciding which block rows exist, collision-free builder IDs, and the monomorphic-locus gate. The resource guard also needs to cover the arrays used by the projection, rather than only its input dosage matrix. These are substantive changes to the proposed behavior or verification; they are not reasons to redesign the whole phase.

| Item | Assessment | When to resolve |
|---|---|---|
| Block rows based on stored term kinds | Can omit induced additive variance and violate coding invariance | Before 4b |
| D5's locus names joined with `x` | Can silently merge different surfaces into one interaction | Before 4a |
| Monomorphic-locus equivalence in B12 | Incorrect if removing a pair also removes its induced main effect | Before 4b |
| Projection resource guard | Does not bound pair arrays or dense source covariance matrices | Before 4b |
| Default base selection | Correct stated intent; requires an explicit whole-individual route | Implementation contract for 4b |
| `trait_name = NULL` | The named shared resolver does not implement the promised default | Before 4b |
| Target/prevalence interpretation | Generated ownership alone does not make a comparison like-for-like | Documentation correction |
| Mutation checks | The proposed subtraction mutation need not fail | Verification correction |

## 1. Block availability must follow the statistical decomposition

**Plan locations:** 4b.4's rule that a block appears only when the model has terms of that kind; B8, B10 and B14.  
**Code:** [structural component classification](../R/genome_effects_eval.R#L90), [target-kind classification](../R/genome_effects_eval.R#L114), and [genotype builder](../R/genome_effect_terms_builders.R#L190).

The stored component names describe term structure. The proposed extractor deliberately reports a different, statistical decomposition. Presence rules cannot simply reuse `component_name`, `.gev_target_kind()`, or the presence of an order-one additive term.

A one-locus surface with values `c(0, 1, 2)` is exactly dosage. `genotype_terms()` stores only indicators, but its statistical additive variance is the entire variance of dosage. In a ten-individual probe with dosages `c(0,0,0,0,0,1,1,1,2,2)`, the current evaluator returned the expected values under component `indicator`. There were **no stored additive members**, while the source decomposition gave additive variance and total variance both **0.6777778**. D2/B14 correctly require this to be case 1 and equivalent to a functional additive model. Suppressing the additive row because no additive term was written would contradict that requirement.

The same issue arises for heterozygote-only models and pair-only models. Dominance can induce an additive coefficient through `b*d`; A×A can induce one through the partners' means. A two-member term containing additive members does not, by itself, tell us which projected main effects are present. Projection at a new reference can also create additive coefficients that were zero at the stored reference.

**Recommendation:** define statistical block availability after canonicalisation and projection. Always retain an induced additive contribution. It is reasonable to report an additive row for every covered model, including a zero value; alternatively, specify a canonical support rule that never suppresses an induced contribution. Retain the intended B8 rule that an absent dominance family does not manufacture a dominance measurement, using canonical `d` support rather than the literal stored contrast name. Distinguish a missing block from a supported block whose variance happens to be zero in this cohort.

Specify availability per trait and per trait pair. A full square is straightforward when both traits support a block. Mixed models need an explicit rule for missing blocks and zero-filled value columns, so one trait's dominance terms do not make an absent dominance target in another trait look measured. Cross-trait covariances between different blocks still belong in `between_components`.

**Additional gates:** dosage-only indicator surface; heterozygote-only model with unequal genotype frequencies; functional pair-only model; Cockerham pair-only model measured after frequency drift; and two traits with different block support. Test report contents as well as total accounting and coding invariance.

## 2. D5's proposed IDs are not collision-free

**Plan locations:** D5, 4a.1 and C16.  
**Code:** [custom locus-name validation](../R/define_genome.R#L177), [term grouping](../R/define_genome_effect_terms.R#L281), and [genotype builder](../R/genome_effect_terms_builders.R#L190).

The current numeric IDs do collide when surfaces are bound. Fixing that is necessary. However, joining arbitrary locus names with `x` does not uniquely encode their list. `define_genome()` requires nonempty, unique locus names; it does not reserve `x`, underscores, or other delimiters.

For example, the surfaces on loci `c("A", "B")` and on the single locus `"AxB"` both receive proposed ID `"AxB_1"` for their first row. Their loci are disjoint. When both coefficients are 1, the writer does **not** reject the combined frame: it groups the rows into one valid three-member indicator term. I reproduced this against the current writer by assigning exactly D5's proposed IDs. The write succeeded with **one term and three members**, changing the intended sum into a product.

The proposed `aa_terms()` IDs have the same problem: pairs `(A, BxC)` and `(AxB, C)` both become `AxBxC`. Cross-builder collisions also need consideration because locus names may contain the suffixes used by other builders.

**Recommendation:** use a shared deterministic encoding with a builder/type prefix and unambiguous encoding of each locus name, such as length-prefixed names or properly escaped structured strings. Canonicalise locus order where it matters. Do not solve this with random IDs or session counters: the builders should remain pure and reproducible.

D5's deliberate refusal when incompatible surfaces over the same loci share IDs is a separate policy. It does not establish uniqueness for different surfaces. State the composability promise precisely; overlapping definitions may still be refused by the writer's family rules.

**Additional gates:** the `A/B` versus `AxB` example; two ambiguous A×A pair names; names containing the chosen delimiter; and collisions across builder types. Assert the written term count, member sets and evaluated sum, not just unique-looking IDs in the builder output.

## 3. A monomorphic pair member can still change additive variance

**Plan locations:** 4b.2 rule 5 and B12 in the main plan.  
**Source:** `nonadd_covariates()` centres each dosage; `nonadd_decompose()` adds `c_l * e_kl` to its partner's additive coefficient.

Setting the within-locus regression to zero when dosage variance is zero is correct. The centred A×A value is also zero when either member is fixed. But the *functional interaction* need not disappear from the total model.

If locus 2 is fixed at dosage 2, then

```text
e * (g1 - 1) * (g2 - 1) = e * (g1 - 1).
```

It is an additive effect at locus 1. If locus 2 is fixed at dosage 0, the induced effect is `-e`. The projection captures this through `c2 = 2*p2 - 1`. Dropping the pair before computing the induced main effect loses real variance.

I evaluated a functional pair-only model with `e = 1`, the ten dosages above at locus 1, and dosage 2 at locus 2. The current evaluator's total variance was **0.6777778**. The source returned additive variance **0.6777778**, A×A variance **0**, and total variance **0.6777778**. Removing the pair without replacing its induced additive effect would instead give zero.

**Recommendation:** interpret “pairs contribute nothing” as referring to the centred A×A block only. Canonicalise first, retain the `e*c` contributions, and then handle zero-variance columns. Rewrite B12's equivalence to compare against a reduced model with the appropriate induced main effect and dropped constant. Cover fixed dosages 0, 1 and 2; a fixed heterozygote has `c = 0` and behaves differently from either fixed homozygote.

This correction also strengthens finding 1: a pair-only model may need an additive row even when its A×A variance is zero.

## 4. Bound the projection's allocations and avoid dense source anchors

**Plan locations:** 4b.2 rule 4, 4b.4 and the proposed benchmark.  
**Code:** [current dosage guard](../R/define_additive_effects.R#L1109).  
**Source:** `nonadd_covariates()` constructs `M_A`, `M_D` and `M_AA` as dense covariance or diagonal matrices.

The current limit is **20 million dosage cells**. It bounds `n*m`, not all the arrays needed by a literal port. In particular, it does not bound `n*r` for `r` interaction pairs or the source's `m*m` and `r*r` anchor matrices.

For example, 2,000 individuals and 500 loci use one million dosage cells. All 124,750 distinct pairs would produce **249.5 million** pair-value cells, about **2 GB** for one double matrix before temporary copies. The source's pair covariance matrix would require about **124 GB**. Even without many pairs, `n = 2` and `m = 100,000` pass the dosage guard while a dense locus covariance matrix requires about 80 GB.

The extractor does not need those dense anchors. It needs component values per individual and small trait covariance matrices. Likewise, the genic formulas should use row weights directly, as `nonadd_decompose()` does for its genic summaries, rather than materialising `diag(weights)`.

**Recommendation:** port the algebra, not the source's unnecessary anchor allocations. Compute A×A values in deterministic pair chunks, accumulating an `n*k` value matrix, or explicitly bound pair allocations before creating them. Keep a stated budget for the other dense work arrays. Perform knowable size checks after resolving the model and cohort and before expensive evaluation/allocation. Preserve the existing generator's size and completeness contract when extracting its dosage helper.

The planned 2,000-by-500 benchmark should state its **pair count and trait count**. Add a separate many-pairs resource gate; the current benchmark and B9 cannot expose this risk if they use only a few pairs.

## 5. Default genic frequencies must use the resolved individuals

**Plan locations:** 4b.1 and 4b.2.  
**Code:** [individual selection](../R/sql_utils.R#L491) and [frequency selection semantics](../R/extract_allele_freq.R#L69).

The stated contract is right: `tbl` selects individuals, and `base_tbl = NULL` uses those individuals' frequencies. But calling `extract_allele_freq(tbl)` directly does not always implement it.

When `tbl` is `ind_haplotype`, `resolve_subset_ids()` selects the distinct individuals represented by its filtered rows. `extract_allele_freq()` deliberately treats that same table as a selection of **allele copies**. A filter on `parent_origin`, `line_origin`, or `locus_id` therefore changes the frequency population if the original table is reused.

In the probe, `ind_haplotype |> filter(parent_origin == 1L)` selected all ten individuals. Its direct allele-copy frequency at locus A was **0.5**; the selected individuals' whole-dosage frequency was **0.35**. Both helpers behaved correctly under their existing contracts.

**Recommendation:** resolve the cohort IDs once and use a whole-individual selection for the default base, independent of which table supplied those IDs. An explicitly passed `base_tbl = ind_haplotype |> filter(...)` should keep the existing copy-selection semantics. Do not change `extract_allele_freq()` to resolve this distinction.

Also require finite, nonmissing frequencies at every retained covered locus. The helper intentionally returns `NA` where the base has no copies; that is different from observed frequency 0 or 1.

**Additional gates:** the same cohort selected through `ind_meta`, paternal haplotype rows, and repeated phenotype records gives identical default-base output; an explicit copy-filtered base gives the appropriately different projection; a base missing one required locus errors with its name.

## 6. The shared trait resolver does not implement the promised default

**Plan location:** 4b.1 says `trait_name = NULL` means every trait with terms, using `.gev_resolve_traits()`.  
**Code:** [trait resolver](../R/genome_effects_eval.R#L829).

`.gev_resolve_traits(conn, NULL)` currently returns **every trait in `trait_meta`**, in `id_trait` order. It does not filter to traits with terms. The probe had traits `T` and `Empty`, with terms only for T; the resolver returned both.

**Recommendation:** implement the extractor's advertised default explicitly, preserving trait order and selecting only traits with active stored terms. An explicitly requested trait with no terms should give a clear error. An entirely empty model should also have a defined error. Avoid changing a shared resolver's defaults incidentally, since other callers rely on it.

## 7. Tighten the target-comparison and prevalence explanation

**Plan location:** 4b.5.  
**Code:** [anchor/base restrictions](../R/define_additive_effects.R#L788), [target resolution](../R/define_additive_effects.R#L820), and [current prevalence variance readers](../R/add_phenotype_stages.R).

Filtering targets by `line_name` and using both `inner_join()` and `anti_join()` is a useful example. However, “all terms are generated” is only one condition for a like-for-like comparison. The anchor, reference selection, measured model, and calibrated scope must also match. A generated line variant can compete with common fallback terms; measuring a mixed cohort measures the active combination, rather than that variant alone. A restored or subsequently edited founder pool need not be the generation-time base.

The paragraph saying the sum of stored diagonals “is the genic total under HWE and LE at the base” is too broad. The additive generator already supports **realised** calibration, and the planned general generator does too. A realised target is not automatically a genic covariance. Scoped targets are also calibrations of variants, rather than a universally valid whole-cohort decomposition.

**Recommendation:** show a fully covered, common-scope model calibrated with `anchor = "genic"` for the genic equality example, and a realised model measured on its actual calibration individuals for the realised equality example. Describe stored diagonals as generation targets under their calibration contract. Explain that prevalence still assumes an appropriate liability distribution; matching variance alone does not guarantee prevalence for a skewed finite-locus genetic model.

For the scoped case, retain the existing “evaluated additive variance” wording and explain why a mismatch with a chosen line target can be meaningful. Do not imply that custom ownership is the only cause of target disagreement.

## 8. Correct the mutation-check expectation

**Plan location:** Verification item 4, requiring B3 to fail if `between_components` is calculated as total minus blocks.

If the component identities are correct, then `total - sum(blocks)` and the ordered sum of cross-component covariances are mathematically equal. An independent test of that ordered sum can therefore pass after this mutation. This remains true even when the total is independently evaluated. A numerical rounding difference is not a reliable mutation detector.

**Recommendation:** keep the direct calculation if its independent accounting is the implementation requirement, but do not promise that this algebraically equivalent mutation must fail. Mutate omission of D–A×A covariance, omission of `unpartitioned` cross terms, or use of only one off-diagonal orientation instead. Those mutations change the specified result and should fail a deliberately asymmetric two-trait fixture.

B2 should independently check block values against the copied oracle; B3 should check cross-block accounting against independent value matrices. Together these are stronger evidence than a claim about detecting an equivalent formula.

## Decisions D1–D7

| Decision | My recommendation | Reason / qualification |
|---|---|---|
| **D1 — `anchor` column plus message** | **Agree, with a provenance limitation stated.** | The anchor must travel with the rows. It does not identify which filtered cohort or base produced them. A saved tibble is self-describing about the estimator, not the full reference population. The message should explicitly distinguish cohort-derived frequencies from an explicit copy/pool base. If saved reference provenance is a requirement, add a reference column; a transient message cannot supply it. No persistent database metadata is needed. |
| **D2 — all three diploid indicators** | **Agree.** | Any one-locus diploid surface has an exact constant + additive + heterozygote representation. The proposed coefficients are correct. Add block-support gates from finding 1. |
| **D3 — project covered terms; retain the rest in `unpartitioned`** | **Agree.** | The family rule is essential because variants compete within an owner; owners then sum. Label the additive block as the additive projection of the covered part, not the whole model's breeding value. |
| **D4 — copy the source measurement oracle into tests** | **Agree.** | I verified that the named functions exist at commit `8f8a97c`. A self-contained oracle is better than a machine-dependent skip. Preserve provenance and keep it independent of the production conversion/projection helpers. Add hand-derived fixtures as well, since agreement with the same source algebra is not independent proof of its scientific interpretation. |
| **D5 — change genotype term IDs** | **Agree with the need; disagree with the proposed encoding.** | Use a deterministic, collision-free, type-prefixed encoding across all builders. The proposed `x` encoding silently merges valid different surfaces in a reproduced case. |
| **D6 — two commits, 0.75.0 and 0.75.1** | **Agree.** | Builders and conversion have their own exhaustive evaluator gates and are independently reviewable. The extractor then builds on tested algebra. Version allocation is bookkeeping in this pre-1.0 package. |
| **D7 — drop zero pairs; error if all are zero** | **Agree.** | It matches `ad_terms()` and avoids adding causality from inert coefficients. Validate lengths, names, finite coefficients/frequencies and duplicate pairs before dropping rows, so zero values do not hide malformed input. |

D3's family rule should use the evaluator's actual `family_key`, including owner and member states. Centres are intentionally excluded from that key. Sum owners during canonicalisation **after** preserving each owner's competition semantics; pooling coefficients before deciding covered/uncovered families could change which terms apply.

For a partial model, align uncovered values to the full cohort and fill an absent uncovered contribution with zero. A selected individual may match a covered family and no uncovered family. The missing-value error should apply to the **whole trait model**, not to the uncovered submodel alone. Skip dosage collection entirely when there are no covered loci.

## Recommendations in the main plan's questions

The questions directly relevant to this phase are Q9, Q13, Q16 and Q17. The remaining unresolved separation in Q20 also benefits from the extractor, but should stay a separate change.

- **Q9 — extract functional effects, do not persist them:** agree. The pure inverse built here is the appropriate foundation. A later export must specify the same coverage restrictions rather than pretending every scoped or higher-order model has one `(a,d,e)` representation.
- **Q13 — one internal conversion pair:** agree, including placing it alongside the builders. I also agree with Phase 4's decision not to route `ad_terms()` through a whole-model converter: its Cockerham argument already is alpha, and it has no pairs. Share the algebra and test its no-pair mean equivalence without inventing a conversion the function does not perform.
- **Q16 — one `between_components` row:** agree. Include every ordered pair and `unpartitioned`. Qualify the rationale: random mating can restore single-locus HWE while LD remains, so cross-block covariance is guaranteed zero by the **HWE + LE genic model**, not by random mating alone.
- **Q17 — functional manual models first; whole-model Cockerham builder later:** agree. A caution in §9.3 needs correction when preparing docs: feeding functional `a` directly as Cockerham alpha generally changes genotype differences, not merely the stored component split. The total differs only by a constant **after the required alpha conversion**. For example, with `a = 0`, `d = 1`, `p = 0.3`, the heterozygote indicator cannot be replaced by the centred dominance contrast alone: the missing additive coefficient is 0.4.
- **Q20 — separate evaluation parameters from generation targets and measured truth:** agree with the separation, the named evaluation model and a transparent default. The current BLUPF90 writer still reads the generation covariance, and `update_covars` still describes an unimplemented write-back into targets. The extractor should not become a writer into either generation table. The exact evaluation-table key, effect-to-component mapping and REML-history policy deserve a separate plan; I would not settle that schema or add overloaded writer behavior during Phase 4.

The other settled/deferred questions do not need reopening to build this phase. In particular, a per-individual breeding-value export, scoped genic decomposition, additional non-additive blocks and a separate cross-block covariance export can remain out of scope.

## Conversion details and gates worth retaining

The proposed per-term inverse is algebraically correct. Complete the constant entries and make the two constant conventions explicit:

```text
stored term = functional term + kappa

additive(v, c):       kappa = v * (1 - 2*c)
dominance(v, c):      kappa = -v * (c^2 + (1-c)^2)
indicator(2, 1, v):   kappa = 0
indicator(2, 0/2, v): kappa = v/2
AA(v, c1, c2):       kappa = v * (1-2*c1) * (1-2*c2)

functional model = converted statistical model + mu
```

For a matching forward/inverse round trip, `kappa = -mu`. Avoid using the same undocumented sign convention for both returned constants. Define whether `.noia_to_stored()` returns a writer-ready frame or coefficient data, and how it becomes the model accepted by `.stored_to_functional()`; the proposed N2 call otherwise leaves that interface implicit. Align coefficients, pairs and frequencies by locus identity, especially at loci that appear only in pairs or indicators.

N1 against the real evaluator is particularly valuable. Use all single-locus states and all nine pair states, arbitrary non-0.5 centres, unequal centres, and multiple owners. N2/N3, coding invariance under both anchors, read-only table comparisons, and restore determinism should stay. Strengthen B9 with the existing thread-count style of determinism test. Keep the observed sample covariance divisor `n-1`; the within-locus regression can use population-form moments because the common divisor cancels.

The writer already accepts the relevant `NA` fields: `.ge_check_member_fields()` requires indicator centres to be missing and non-indicator states to be missing, and `.ge_infer_copy_counts()` fills omitted indicator copy counts where unambiguous. The probe successfully wrote additive + functional dominance with those missing fields. Introduce the fixed typed frame and test it before considering any relaxation of writer validation.

Use a per-trait, mixed absolute/relative tolerance for the centred-value identity. A large raw genetic offset should not by itself enlarge the allowed decomposition error. The planned O(1) fixtures and explicit tolerance notes are appropriate.

## Verification and limits

- Read the builders, writer validation and inference, family keys, scope predicates and evaluator, inheritance checks, target/base resolution, dosage collection, frequency extraction, total view, existing test fixtures, and relevant phenotype/evaluation readers.
- Read the source measurement functions directly, including their form at commit `8f8a97c`; no external source checkout was altered.
- Ran package-level probes on in-memory populations using the current evaluator/writer. They reproduced the induced additive block, monomorphic-partner main effect, successful silent merge under D5's proposed ID, trait-default discrepancy and whole-cohort versus copy-base frequency distinction.
- Loaded the package with `compile = FALSE`; no package sources or generated documentation were rebuilt.
- Ran `test-genome-effects-writer.R`, `test-extract_allele_freq.R` and `test-define_additive_effects-anchor.R` with `NOT_CRAN=true` and `stop_on_failure = TRUE`: all passed, process exit 0, with no test warnings reported. Phase 4's new functions do not exist yet, so these checks establish current integration behavior rather than validate the proposed implementation. I did not run the full suite or rebuild documentation for this review.

**Review artifacts:** the temporary probe script and output are `/tmp/tidybreed_phase4_review_probes.R` and `/tmp/tidybreed_phase4_review_probes.log`; focused-test output is `/tmp/tidybreed_phase4_review_tests.log`.
