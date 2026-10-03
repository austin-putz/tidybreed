# Review of `import_qtl_effect_methods.md` (current draft, 2026-10-02)

**Scope.** Full reread of the current 2,003-line plan and prior review, relevant tidybreed code, and the source project's `nonadd_decompose()`, with three independent reviews of plan logic, code alignment, and edge cases. This reviews the proposal; none of the proposed releases is implemented.

## Assessment

The plan is substantially stronger than the September draft. The old review's main questions are answered: the coding constant and phenotype intercept are explicit; the extractor has supported model classes and cross-block accounting; manual effects and pairs have named inputs; target writes are atomic; additive exactness has limits; and phenotypes use the total genetic value. The five-release order and acceptance gates are detailed.

**Two material design issues remain, plus several smaller implementation traps.** The main issue is that `trait_var_comp` is described both as the single source of truth for the *active* model and as a catalogue whose blocks a call may filter out. The proposed binary-prevalence threshold makes this a user-visible bug. Settle these before Step 2. I do not see evidence that the core congruence or statistical-coding approach needs redesign.

## Findings to resolve in the plan

### 1. Stored targets do not identify the active model — high

Section 6C says the table describes the model a generator builds ([plan, lines 808–811](import_qtl_effect_methods.md)), but then permits filtering out stored A×A or dominance targets while leaving those rows in the table (lines 863–886; C17, line 1625). The generator replaces `generated` terms, so the excluded block is no longer active. A later default call reads it again. The target-versus-extracted `inner_join()` (lines 1225–1231) can compare a target for a block absent from the model. Q15's “always honest” join claim (lines 1815–1837) is too strong.

Two related cases need the same rule:

- **Partial trait sets.** A stored full matrix for `ADG` and `BF` looks like a valid 1×1 block if a later call selects only `ADG` and filters rows to the call's traits. Replacing `ADG` alone may invalidate the retained `ADG`–`BF` target. The current table has no block identifier ([`R/open_pop.R`, lines 277–284](../R/open_pop.R)). Check the whole stored block, or explicitly support partial-model semantics.
- **Manual and custom terms.** Section 7.1 intentionally ignores stored targets for manual `effects` and `scale_to_target = FALSE` (lines 960–967); §5 lets custom-owner terms sum with generated terms (lines 488–519). The extractor sums owners (lines 1195–1212). Its whole-model covariance need not equal a target calibrated for the generated owner alone. This is valid behavior, but the target comparison must name the owner/model it describes.

**Recommendation:** decide whether `trait_var_comp` holds *available targets* or the *active model specification*. For available targets, narrow the single-source claim and supply an active-block/owner-aware comparison path. For an active specification, filtering out a block must deactivate or remove it atomically. Add gates for a two-trait block followed by a one-trait rerun, an excluded stored block, and a generated-plus-custom model.

### 2. Binary prevalence can use variance from inactive effects — high

P2 calculates the prevalence threshold from the sum of **stored** additive, dominance, and A×A diagonals (lines 680–690). Section 6C permits the active model to omit stored blocks, and PH6 permits a hand-written model without stored targets. For example, filter an A+D+A×A target to A+D, rerun the generator, then call `add_phenotype(prevalence = ...)`: the threshold still includes A×A. This changes prevalence even if generation and evaluation work correctly. PH7 (line 1650) tests only all-stored A+D and additive-only cases.

**Recommendation:** derive the threshold from the active model and declared reference, or require an explicit liability variance/threshold when targets cannot describe it (manual or custom terms, composite traits). State the population/model behind any approximation. Add PH7 cases for a filtered-out target and a hand-written model with no target. Decide this with finding 1, before Step 3.

### 3. Group-contributor sums still need deterministic accumulation — medium-high

P2 fixes the `ind_tgv_total` view's floating `SUM()` so `self` totals are bit-identical across DuckDB thread counts (lines 667–676). The existing group-mate query separately does `SUM(t.tbv_value)` across mates ([`R/contributor_tbv.R`, lines 131–143](../R/contributor_tbv.R)). Replacing its input with `ind_tgv_total` leaves that second floating reduction. `group_sum()` and `group_mean()` can still vary with thread count. Apply the same bounded exact accumulation at the group aggregation point, including component-specific lookup. Extend PH5 with several mates and both thread settings.

### 4. The source oracle omits a realised cross covariance — medium

Section 8 correctly defines `between_components` from **every** ordered cross-block covariance (lines 1183–1191). The source project's `nonadd_decompose()` returns `cov_A_D` and `cov_A_AA`, but not `cov_D_AA` (`simulate_qtl_effects/non-additive/R/qtl_effects_nonadd.R:480`). D–A×A covariance can be nonzero under LD or non-HWE. B2/B3 (lines 1596–1599) should compute it directly from the source's `DD` and `AA` matrices, including both orientations for off-diagonal trait pairs. Keep the independent block-plus-cross-terms-equals-total test.

### 5. Define extractor cohort and missing-contribution behavior — medium

The realised API accepts any individual selection (§8, lines 1150–1176), but sample covariance is undefined for fewer than two contributing individuals. The current evaluator omits individuals with no contribution ([`R/genome_effects_eval.R`, lines 628–639](../R/genome_effects_eval.R)). If traits have different contributing subsets, covariance on the union, per-trait subsets, and complete pairs differ, as does `n_ind`. Specify one aligned-cohort rule for every trait pair, whether missing contributions error or exclude individuals, and behavior for `n_ind < 2`. Also specify how a monomorphic locus or zero observed variance is handled in realised NOIA projection. Add a small gate for missing contributions and one selected individual.

### 6. Distinguish anchor rank from drawn-architecture rank — medium

Section 1.1 says a target is feasible iff `rank(G) <= rank(M)` (lines 224–240), an existence claim for *some* effects. Congruence on a fixed drawn `B0` can deliver rank r only when `rank(B0' M B0) >= r`. A rank-deficient `B0` can fail even with full-rank `M`. A4 explicitly tests a rank-one `B0` and rank-one target; A5 should distinguish an insufficient anchor rank from an insufficient architecture rank. This improves failure messages for small QTL sets.

### 7. Qualify the realised-anchor variance wording — low

Section 6B rule 3 says generated `ind_tgv` additive variance **is** the stored target and a join says something true (lines 750–757). Section 4.3.1 correctly says stored HWE-coded additive values under a *realised* anchor can differ from that target (lines 459–466; example 4.134 versus 4). Qualify rule 3 as a genic covariance at the generation base for a genic-anchored generator; use the extractor's realised projection for a realised-anchor target check. Variance is a property of a population, not of “each individual's” value.

## Earlier review: current disposition

| September concern | Current plan |
|---|---|
| Functional/statistical total constant and phenotype intercept | Resolved in §§3 and 6A, C2/PH4. |
| Cross-block variance accounting | Resolved in §8 and B3, subject to finding 4's oracle. |
| Arbitrary stored models cannot all be NOIA-decomposed | Resolved by §8's three cases; finding 5 fills an input edge case. |
| Named manual effects, pair validation, exact writer path | Resolved in §§7.1 and 9.1–9.3, A12/C15/C16. |
| Absent, stored and zero targets | Resolved at call level by §6C/C17; finding 1 concerns persistence after the call. |
| Atomic target and term writes | Resolved in §§6C/7.1/9.2, A10/C11. |
| Additive exactness under shared, union, manual and unscaled paths | Resolved in §7.1–7.3, A7/A19. |
| Diagnostics, size guard, owner naming, phenotype paths and DSL rename | Resolved in §§5–7, 9–11. |

The old review's API names are obsolete: `define_genome_effect_terms()` is now the writer, and `define_genome_effects()` the new generator. Its four-step edit order is obsolete; §10 specifies five releases.

## Suggested changes before implementation

1. Settle the active-target contract, including multi-trait blocks, custom owners, and the meaning of a target-versus-realised comparison.
2. Use that decision to specify binary prevalence for filtered and hand-written models.
3. Add deterministic group aggregation, extractor cohort rules, the complete cross-block oracle, and rank-specific errors to the relevant sections and gates.
4. Correct §6B rule 3's realised-anchor wording.

After these changes, the remaining work looks like implementation and verification risk rather than a missing mathematical design. The source methods and storage mapping are coherent, and most formerly open questions now have explicit decisions.

## Disposition (2026-10-02)

All seven findings were checked against the plan, the code and `$SRC` and accepted. Each is now in `import_qtl_effect_methods.md`:

| Finding | Where |
|---|---|
| 1. Active vs stored targets | §6C: the table holds *available* targets; a stored block whose traits extend beyond the call errors under the default `NULL`. §8: block rows only for components the model has, and the generated-only caveat. Q15 wording narrowed. Gates A16, B8. |
| 2. Prevalence from inactive effects | §6A: active-block rule; errors naming `thresholds =` for untargeted terms and composite phenotypes. Gate PH7. Today's silent-zero bug is §10 step 0b, B-2. |
| 3. Group sums | §6A bullet, PH5. Live bug: §10 step 0b, B-1. |
| 4. Missing `Cov(D, AA)` | §8, gates B2/B3. |
| 5. Extractor cohort | §8 "Cohort" bullet, gate B12 (missing values error; `n_ind < 2` errors; monomorphic loci give `b = 0`). |
| 6. Anchor vs architecture rank | §1.1, §7.1, gate A5 (a)/(b). |
| 7. Realised-anchor wording | §6B rule 3. |
