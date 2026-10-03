# Review: import QTL-effect methods, phases 0b and 1

**Reviewed:** 2026-10-03, at `a07e1a9` (tidybreed 0.72.2).  
**Scope:** The implemented phase 0b commit `4e22d7f`, phase 1 commit `bd54462`, their follow-up commits `fd8c053` and `a07e1a9`, the two phase summaries, and the relevant requirements in `import_qtl_effect_methods.md`. Phases 2–5 have **not** been implemented and are not treated as missing code here. The unrelated working-tree change to `plans/consolidate_genetic_values.md` was left untouched.

## Findings

### 1. High — phase 0b still permits a silently wrong prevalence threshold when a target exists but the active effects were not calibrated to it

`.ap_check_prevalence()` checks only that an `additive` diagonal exists (`R/add_phenotype_stages.R:1152–1169`). `.ap_liability_records()` then uses that stored number in `sqrt(Va + Ve)` (`R/add_phenotype_stages.R:1183–1207`). However, `define_additive_effects(effects = ...)` skips scaling even when `trait_var_comp` already holds a target (`R/define_additive_effects.R:225–283`). `scale_to_target = FALSE` has the same issue. The check therefore establishes *presence of a row*, not whether the row describes the liability being generated.

**Reproduction on the current code:** With `set.seed(71)`, define 20 loci and 200 founders from fixed founder haplotypes. Define a simple trait with `target_add_var = 1`, assign the first ten QTL `effects = rep(10, 10)`, define a categorical phenotype with `prevalence = 0.1` and `residual_var = 1`, then phenotype the founders. The call succeeds. The stored additive target remains 1, the observed TBV variance is approximately 552.76, and 45.5% of records fall in the upper category. The exact rate is seed and sample dependent; the large mismatch is the point. This is a silent modeling error, not a failure of the new no-target test.

The main plan schedules an **active-model/target rule in phase 3** (`import_qtl_effect_methods.md`, §6A, prevalence bullet). As written, that rule tests which component kinds have active terms and whether each has a stored target. It still cannot distinguish calibrated effects from manual or unscaled effects when an old target row is present. In fact, §7.1 says manual effects ignore stored targets while keeping those rows available. The phase 3 design therefore needs a way to establish that a target describes the *active* terms, such as a comparison with delivered variance at the reference base or explicit calibration provenance. Add a regression covering *target present + unscaled/manual effects*. Until then, phase 0b fixes the missing-row and composite cases only; the phase summary and NEWS could say so explicitly. An explicit `thresholds =` is the current workaround for this model.

### 2. Medium — phase 1 does not reject pre-0.72 databases despite declaring their stored strings unreadable

Phase 1 changes stored `effect_name = 'gen_add'` to `'additive'` and reserved `effect_owner` to `'generated'`; `NEWS.md` says earlier databases are not readable and no migration exists. `restore_pop()` (`R/restore_pop.R:77–142`) checks old **column shapes** but does not check either old string. Thus it successfully restores a database with the previous strings. Current `get_trait_var(pop, 'additive', trait)` then returns `NA` (`R/define_effect_cov_matrix.R:179–187`), while the old target row is still present under `gen_add`. Writing a new additive target deletes only rows with the new name (`R/define_effect_cov_matrix.R:127–151`), so old and new vocabularies can coexist in one database.

**Reproduction:** Create a file-backed population, define trait `T` with target 1, change that row's `effect_name` to `gen_add` to model the prior release, close it, and call `restore_pop()` under 0.72.2. Restore succeeds; `get_trait_var(..., 'additive', 'T')` returns `NA`, while `get_trait_var(..., 'gen_add', 'T')` returns 1. No migration is needed to fix this: detect the old target or owner strings at restore time and stop with a rebuild message. Alternatively, document explicitly that restore can succeed but later operations can misinterpret the old model. The former matches `restore_pop()`'s existing stale-schema policy.

### 3. Low — the phase 0b regression test proves the formula route, but does not exercise the `components` route end to end

The B-1 implementation is in the shared `.group_mate_tbv()` helper (`R/contributor_tbv.R:119–150`), and code inspection confirms both routes call it. `test-group-contributor-determinism.R` covers the helper's sum/mean and an end-to-end `formula_tbv` phenotype, but not a `phenotype_components` group contributor. This is a coverage gap, not an observed defect. A small components-route assertion using the existing large-pen fixture would protect the second caller when phase 3 retargets the helper to `ind_tgv_total`.

## Checks that passed

- **B-1 implementation:** The mate query casts each `tbv_value` to `GEV_ACC_TYPE` before `SUM`, then casts the exact sum to `DOUBLE`. `COUNT` and `total / n_mates` preserve the documented mean and no-mate behavior. The group helper's exception handler uses the evaluator's accumulator error formatter. The existing large-pen test compares one- and eight-thread results with `expect_identical()`; its end-to-end formula test was checked against the old code in the phase summary.
- **B-2 missing-target and composite cases:** `define_phenotype()` refuses composite prevalence before writing metadata (`R/define_phenotype.R:248–266`). The PLAN-stage backstop catches a manually inserted composite row, and the missing-target check runs before TBV writes or RNG draws (`R/add_phenotype_stages.R:148–171`). It correctly exempts `user_values`. Tests cover both refusal paths, the no-write/no-RNG result, and explicit thresholds.
- **Phase 1 rename:** The exported writer, `NAMESPACE`, Rd, `_pkgdown.yml`, call sites, reserved owner, genetic covariance readers/writers, BLUPF90 reader, tests and current documentation use the new names consistently. A repository scan of `R`, `tests`, `man`, `vignettes`, `dev`, `tools`, `NAMESPACE`, `README.md`, `CLAUDE.md`, `package_summary.md`, `_pkgdown.yml` and `.claude/skills` found no remaining old writer/owner/variance names or quoted `epistasis`. The remaining `trait_names =` hits belong to `define_index()`, which is outside this phase. The phase 1 test diff changes names and one fixture label, with no assertion or tolerance change.
- **Phase boundaries:** `ind_tbv`, `formula_tbv`, and additive-only phenotype assembly still exist by design. The total-genetic-value conversion, active-target contract, and new `define_genome_effects()` generator belong to phases 2–5.

## Verification and limits

- Read the plans, phase summaries, current implementation and commit diffs. Verified the old-name census with `rg` and the export/reference entries directly.
- Ran two isolated R reproductions with `devtools::load_all()`: the prevalence mismatch in an in-memory population and acceptance of a file-backed old-string database by `restore_pop()`. Both completed successfully and produced the results above.
- The phase summaries report a passing full suite after phases 0b and 1 (964 tests, 3521 expectations at 0.72.0). A new full-suite run under 0.72.2 was stopped after reaching the `add_phenotype_failure_contract` file; it did not produce a full-suite result. The focused `group-contributor-determinism` test file passed under 0.72.2.

## Recommended order

1. Reject old stored strings in `restore_pop()` before a partially compatible database can be used.
2. Make the phase 0b documentation describe its actual missing-target/composite guarantee; use explicit thresholds for currently unscaled models.
3. Resolve the calibration-provenance gap in the main plan before phase 3, then implement its active-model prevalence check and target-present/unscaled regression.
4. Add components-route coverage when retargeting group lookup in phase 3.

---

## Response (Claude, 2026-10-03, tidybreed 0.72.3)

All three findings are accepted. Two are fixed in code, and one is a design question for
the user.

### 1. Prevalence with a stored but uncalibrated target — **accepted; design gap, now open as Q21**

The reproduction is correct. This is not a step-0b regression: before 0.71.2 the same
call used the same stale target. But the phase summary's wording invited reading B-2 as
covering it, and §6A's active-block rule as written would not have caught it either.
Done in 0.72.3:
- the `define_phenotype(prevalence = )` roxygen says the target is taken as given, and
  names `thresholds =` for uncalibrated effects;
- NEWS and an addendum in `_phase_0b.md` state B-2's actual guarantee: missing target
  and composite phenotype only.

No code check was added now. Each candidate rule commits the API in a different
direction, so the rule waits for the user's decision. Q21 sets out the options:
- (a) the owner is the provenance; generators write only calibrated terms (recommended);
- (b) a check against the genic variance of the active terms (no API change, but it
  fails correctly calibrated `anchor = "realised"` models);
- (c) stored provenance.

§6A, step 3 and gate PH7 now point to Q21. PH7 gains the target-present/unscaled
regression, using this reproduction.

### 2. `restore_pop()` accepts pre-0.72 strings — **fixed**

`restore_pop()` refuses a file whose `trait_var_comp.effect_name` holds `'gen_add'` or
`'epistasis'`, or whose `genome_effects.effect_owner` holds `'generated_additive_tbv'`.
The message names the strings found. It goes through the existing `stop_stale()`, so the
connection is closed and the message says to rebuild. Two tests are in
`test-restore_pop.R`, one for the old target and one for the old owner.

### 3. Components-route coverage — **fixed now rather than in step 3**

`test-group-contributor-determinism.R` has a `phenotype_components` test with the same
contributors as the formula test (self `D`, group-sum `S`, group-mean `D`), compared
with `expect_identical()` at 1 and 8 threads. It was checked against reverted code:
with `.group_mate_tbv()` back on a plain `SUM(t.tbv_value)`, both the formula and the
components tests fail, and both pass on the fix. A first version with only one group
contributor passed on the broken code, so it was replaced.

**Decision on finding 1 (user, 2026-10-03): Q21 option (a).** In step 3,
`define_additive_effects()` loses `effects` and `scale_to_target` and always calibrates.
Exact coefficients go through `define_genome_effect_terms()`. The prevalence threshold
then trusts a stored target only when every active term is `generated`, and
`define_effect_cov_matrix()` refuses to rewrite a target under existing generated terms.
