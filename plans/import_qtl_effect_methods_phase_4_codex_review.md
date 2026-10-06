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
