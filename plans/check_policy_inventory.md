# Check policy — inventory appendix

Raw inventories behind `plans/check_policy.md` §4, one per area, compiled by four
read-only agents against 6136b89 (v0.76.12). Line numbers are from that commit.
The plan is the source of decisions; this file is evidence. Claims marked
"verified" in the plan were re-checked by hand; the rest are as reported.

Areas: A. genome effects, B. phenotype and evaluation, C. genome and breeding,
D. infrastructure.


---

## Check inventory: genome effects, generators, evaluator, add_tgv

Scope: R/define_additive_effects.R, R/define_genome_effect_terms.R,
R/define_genome_effects.R, R/genome_effect_terms_builders.R, R/qtl_congruence.R,
R/genome_effects_calibration.R, R/genome_effects_eval.R,
R/genome_effects_helpers.R, R/remove_generated_effects.R,
R/extract_genetic_variance.R, R/add_tgv.R, R/contributor_tgv.R, plus
`assert_qtl_autosomal()` (R/chr_meta_helpers.R:388), which the generators call.
Line numbers are from the tree at commit 6136b89 (v0.76.12).

Labels used below:
- **1-Q21**: refusing is the only way to keep the reserved owner `"generated"`
  meaning "calibrated to the stored target". Bypassing would write generated
  terms that do not deliver the target. Most of these are mathematically
  infeasible requests (rank, floor, bound), not assumptions. Never bypassable.
- **R (resource)**: a size or memory guard. It is neither an assumption nor an
  integrity rule (see notes).
- **pass**: a re-throw (`stop(e)`) or a wrapper that adds context to an error
  raised elsewhere. Not a check of its own.

### 1. Category 0 counts per file

Each count covers `stop()`/`warning()` call sites that are only argument-shape
checks. Everything else is listed in the table.

| file | stop/warning sites | category 0 |
|---|---|---|
| define_additive_effects.R | 42 | 24 |
| define_genome_effect_terms.R | 35 | 30 (includes the aggregated row-local messages at 585 and 637) |
| define_genome_effects.R | 35 | 25 |
| genome_effect_terms_builders.R | 31 | 24 (the other 7 are internal assertions) |
| qtl_congruence.R | 35 | 27 (5 more are in `.qtl_sim_effects()` paths reached only by tests, listed below) |
| genome_effects_calibration.R | 10 | 0 |
| genome_effects_eval.R | 12 | 2 |
| genome_effects_helpers.R | 4 | 1 |
| remove_generated_effects.R | 4 | 4 |
| extract_genetic_variance.R | 10 | 6 |
| add_tgv.R | 14 (+1 `stopifnot` at 150) | 9 (+ the `stopifnot`) |
| contributor_tgv.R | 4 | 3 |
| chr_meta_helpers.R (`assert_qtl_autosomal` only) | 1 | 0 |

Rough totals outside category 0: about 30 category 1 (including 16 1-Q21 and 8
internal assertions), 9 category 2, 7 category 3, 4 R, 4 pass, and 11 warnings.

### 2. Category 1, 2 and 3 checks, and every warning

Entry points: DAE = `define_additive_effects()`, DGE = `define_genome_effects()`,
DGET = `define_genome_effect_terms()`, RGE = `remove_generated_effects()`,
EGV = `extract_genetic_variance()`, AP = `add_phenotype()` (Stage 1 calls the
`add_tgv()` engine through `.tgv_compute_ids()`).

### Category 3: modelling assumptions (where the `allow` design matters)

| file:line | function (entry) | what it refuses | cat | proposed `allow` | default | downstream behavior needed if bypassed | trusted by |
|---|---|---|---|---|---|---|---|
| define_genome_effect_terms.R:860 | `.ge_assert_complete` (DGET, only when `require_complete = TRUE`) | An indicator genotype surface that leaves out a reachable (copy_count, dosage) state. | 3 (today it is opt-in) | `incomplete_genotype_surface` | **refuse** (flip the default: the policy names "a missing genotype state" as refuse-level, but today the default `require_complete = FALSE` lets it pass silently) | Already defined: the evaluator scores an unmatched state as 0 (`HAVING` every slot matched). EGV treats each indicator term separately, so nothing there assumes completeness. | none |
| genome_effects_eval.R:917 | `.gev_require_contribution` (`add_tgv`, AP Stage 1) | An individual for whom no term of a trait matches any allele copy, which would give a silent zero (e.g. Y-linked QTL in females, a paternal-only X effect in males). | 3 | `no_matching_copy` (shared with EGV:470) | refuse | `.gev_write_tgv` must write an explicit 0 row (component `additive`, or one per expected component) so `.tgv_read`/`ind_tgv_total` return 0, not NA. Otherwise phenotype paths would treat it as a missing component (`missing_component_action`), which is different from 0. | `.tgv_true_index` (add_tgv.R:453 refuses NA), composite/formula phenotypes |
| extract_genetic_variance.R:470 | `.egv_realised` (EGV) | Same condition at EGV time: a selected individual has no genetic value for a trait. | 3 | `no_matching_copy` | refuse | Treat the missing value as 0 in `G` and in the blocks. The internal sum check at :545 must agree. | none |
| genome_effects_eval.R:899 | `.gev_require_terms` (`add_tgv`, AP) | A trait with no genome-effect terms at all (genetic value identically 0). | 3 | `trait_without_effects` | refuse | The same 0-row write as `no_matching_copy`. EGV must then report 0 variance rather than refuse at :274 (category 0 today). | AP Stage 1 |
| define_additive_effects.R:952, :965, :1023 | `.dae_resolve_target` (DAE) | A stored `dominance`/`additive_by_additive` target exists, and this additive-only call would ignore it. | 3 | `ignore_stored_target` | refuse (the stored target describes a model that will not be generated) | Already defined: the prevalence threshold ignores a target with no terms of its kind (remove_generated_effects.R roxygen; `define_phenotype`). This is the same state the explicit `trait_var_comp_tbl` filter route reaches today, so the allow would be sugar for that filter. | `define_phenotype(prevalence = )`, which re-checks it ("target with no terms ignored") |
| contributor_tgv.R:49 (missing-row half) | `.read_one_per_id` (AP: covariate weights add_phenotype.R:355, residual conditions add_phenotype_stages.R:545, group lookups via `.register_focal_groups`) | A planned individual has no row in the covariate, condition or group table. | 3 | `missing_contributor_row` | refuse | Return NA, which already means "missing component" (a `NULL` group value gives NA today), then `missing_component_action`. For a residual condition, NA must map to a defined condition level, or stay refused. | none |
| define_additive_effects.R:1507 (also called from DGE :430) | `.dae_require_base_at` | The base population has no allele copies at a selected QTL. | 3 (weak; arguably 0) | `qtl_outside_base` | refuse | Drop those loci from the QTL set before the draw (a defined behavior that changes the QTL count). Never treat p as 0 or 1. | genic weights `n_eligible*p*q` |

### Category 2: not implemented

| file:line | function (entry) | what it refuses | cat | if implemented, possible `allow` / note | trusted by |
|---|---|---|---|---|---|
| chr_meta_helpers.R:388 | `assert_qtl_autosomal` (DAE :451, DGE :432) | A QTL on a chromosome that is not 1,1 for both offspring sexes (X, Y, organelle). Genic weights assume 2pq. | 2 | Bypassing would write miscalibrated `"generated"` terms (a Q21 breach), so it is not category 3 until sex-linked calibration exists. Future name: `non_autosomal_qtl`. The message names DAE even when DGE calls it. | DAE/DGE calibration, prevalence |
| genome_effects_helpers.R:348 (message assembled at :376, raised via :325 / define_genome_effect_terms.R:1109) | `.ge_validate_dominance_ploidy` (DGET; every write and RGE via `validate_genome_effects`) | A `dominance` member at a locus whose chromosome is not 1,1 for both sexes, unless that member's own origin multiset demands 2 copies. | 2 | Could become category 3 `non_diploid_dominance` if the evaluator defines "dominance contrast = 0 when copy count ≠ 2" (a hemizygous male has no dominance deviation). Both this check and eval :770 would have to honour it. With no stored flag, the stored term itself must carry the behavior (see notes). | eval :770 |
| genome_effects_eval.R:770 | `.gev_evaluate` (`add_tgv`, AP, EGV) | At evaluation, a dominance contrast is selected for an individual carrying ≠ 2 copies. | 2 (also 1: haplotype data disagree with `chr_inheritance`) | The re-check of the row above. If `non_diploid_dominance` is added, this is where the 0 is produced. | — |
| define_genome_effect_terms.R:508 | `.ge_infer_copy_counts` (DGET) | `copy_count_value` cannot be inferred for an indicator on a chromosome with more than one reachable copy count. | 2 (no defined default; the fix is explicit input) | Not bypassable. Uses line-agnostic `.ge_reachable_copy_counts` (see notes). | — |
| define_additive_effects.R:377 | DAE | `parent_origin` differs across traits in one call (a mixed-scope anchor). | 2 | — | — |
| define_additive_effects.R:914 | `.dae_check_realised` (DAE; DGE :358) | `anchor = "realised"` with `line_name` or `parent_origin`. | 2 | — | — |
| define_genome_effects.R:392 | DGE | A singular `G_A` together with live dominance or epistasis ("this release"). | 2 | — | — |
| define_genome_effects.R:455 | DGE | `n_pairs` > floor(m/2): random matching with a locus in several pairs is not implemented (pass `pairs`). | 2 | — | — |
| define_genome_effects.R:625 | `.dge_check_inbreeding` (DGE) | `inbreeding_depression` with `dominance_degree_sd = 0` (the ratio is fixed, so it cannot be targeted). | 2 | — | — |
| genome_effects_helpers.R:683 | `.collect_dosages` (DAE/DGE realised anchor; EGV realised :489; swallowed into a message in diagnostics at define_additive_effects.R:1339 and define_genome_effects.R:933) | A base or selected individual lacks both copies at every locus (crossbred copy subsets, non-diploid, no haplotypes). | 2 | Cov(X) of partial genotypes is not defined. | — |
| extract_genetic_variance.R:206 | EGV | `anchor = "genic"` on a model that is not fully decomposable (scoped, higher-order, or off-diploid terms). | 2 | — | — |
| extract_genetic_variance.R:644 | `.egv_base_freq` (EGV genic) | The base has no copies at a decomposed locus. | 2 | No defined p. The locus has effects, so it cannot be dropped. | — |

### Category 1: integrity and safety (never bypassable)

| file:line | function (entry) | what it refuses | cat |
|---|---|---|---|
| define_genome_effect_terms.R:245 | `.ge_write_terms` (DGET) | A user writing under the reserved owner `"generated"`. | 1 (Q21) |
| define_genome_effect_terms.R:939 | `.ge_resolve_deletes` (DGET `mode = "replace_trait"`) | `replace_trait` deleting `"generated"` terms. | 1 (Q21) |
| define_genome_effect_terms.R:1096, genome_effects_helpers.R:325 | `.ge_commit` / `validate_genome_effects` (DGET, DAE, DGE, RGE) | Cross-row violations: orphan members or origins, duplicate family+scope, overlapping-but-incomparable scopes (no variant could be selected), bad member_slot or origin_slot, multiset sum ≠ state, plus the dominance-ploidy messages (category 2 above). | 1 |
| genome_effects_eval.R:608 | `.gev_select_variant` (`add_tgv`, AP, EGV) | Ambiguous specificity at evaluation (re-checks the write-time rule). | 1 |
| genome_effects_eval.R:292 | `.gev_additive_alphabet` | A line name containing `,` or `@`, which would break the label encoding. | 1 |
| genome_effects_eval.R:695 | `.gev_accumulator_error` | Overflow of the exact DECIMAL accumulator (`GEV_ACC_TYPE`); bypassing would break bit-identical sums. | 1 |
| genome_effects_eval.R:94; add_tgv.R:322; extract_genetic_variance.R:545; genome_effect_terms_builders.R:578, 581, 710, 715, 720, 761, 779 | various | Internal assertions. | 1 (internal) |
| define_additive_effects.R:777 | `.dae_refuse_nonadditive_model` (DAE) | Replacing part of a jointly calibrated DGE model. | 1 (Q21) |
| define_additive_effects.R:984; define_genome_effects.R:694 | `.dae_resolve_target` / `.dge_resolve_targets` | A new `G` when a block is already stored (stored targets are never overwritten). | 1 (Q21) |
| define_additive_effects.R:996 | `.dae_resolve_target` (DAE) | A new population-wide target that line-fallback generated terms (calibrated to the old one) would silently keep missing. | 1 (Q21) |
| define_additive_effects.R:1093 | `.tvc_block_from_rows` (DAE, DGE) | Calibrating a subset of a stored block that links other traits (would break the stored covariance). | 1 (Q21) |
| define_additive_effects.R:1119 | `.tvc_block_from_rows` | `trait_var_comp_tbl` selecting a block other than the one this scope reads. | 1 (Q21) |
| define_additive_effects.R:427 | DAE `method = "union"` | A positive target for a trait with no existing QTL. | 1-Q21 (infeasible) |
| define_genome_effects.R:401 | DGE | Both dominance-degree mean and sd are 0 with a non-zero `G_D`. | 1-Q21 (infeasible) |
| define_genome_effects.R:463 | DGE | `G_AA` rank > number of pairs. | 1-Q21 |
| define_genome_effects.R:473 | DGE realised | Target rank > n_ind − 1. | 1-Q21 |
| define_genome_effects.R:527 | DGE | Drawn random pair design is rank-deficient. | 1-Q21 |
| qtl_congruence.R:167 | `.qtl_anchor_rank_check` (DAE, DGE) | The anchor rank is below the target rank. | 1-Q21 |
| qtl_congruence.R:274 | `.qtl_calibrate` (DAE, DGE via calibration.R:117) | The calibration misses the target beyond tolerance. | 1-Q21 |
| qtl_congruence.R:416 | `.qtl_congruence` | The drawn architecture rank is below the target rank. | 1-Q21 |
| genome_effects_calibration.R:233, 287, 360, 379, 468, 475, 493, 503 | `.na_calibrate` / `.na_store_alpha` / `.na_additive_stage` / `.na_solve_dd_mean` (DGE) | Additive target missed; Cockerham storage inexact; additive architecture rank deficient; `G_A` below the additive floor; the inbreeding-depression ratio cannot be attained (limit, Cauchy-Schwarz, sign/variance) or does not verify. | 1-Q21 |
| genome_effects_calibration.R:117, 215 | wrappers | Add block or trait context to the errors above. | pass |
| add_tgv.R:489, 557; genome_effects_eval.R:694 | rollback re-throws | Rollback then re-throw (failure contract). | pass |
| contributor_tgv.R:49 (several-rows half) | `.read_one_per_id` (AP) | Several rows per individual in a covariate, condition or group table. Picking one would depend on row order (breaks "no implicit ordering"). | 1 |

### Resource guards

| file:line | function | what it refuses | note |
|---|---|---|---|
| define_genome_effects.R:798 | `.dge_dosage_guard` (DGE realised) | Design arrays above `QTL_REALISED_MAX_CELLS` (2e7). | R. The out-of-core realised path is not implemented, so it behaves like category 2. |
| genome_effects_helpers.R:631 | `.dosage_guard` (DAE/DGE realised, EGV) | Dosage matrix above 2e7 cells. | R |
| genome_effects_eval.R:661 | `.gev_preflight` (`add_tgv`, AP, EGV) | Label-vector enumeration above `options(tidybreed.label_vector_max)` (1e6). | R. **Already bypassable through `options()`**, which is a second relaxation mechanism (see notes). |

### Warnings

| file:line | function (entry) | what it warns | right level? |
|---|---|---|---|
| define_additive_effects.R:360 | DAE | `distribution` is ignored for multi-trait calls (always MVN). | **No.** It is an ignored argument, so it should be a category 0 refusal; it fires only when the user explicitly passed a non-normal distribution. |
| define_additive_effects.R:497 | DAE `method = "union"` | Covariances are not calibrated; the delivered covariance is reported. | Yes, warn (the user asked for union). But the terms are written as `"generated"` while missing the stored off-diagonals (see notes on Q21). |
| define_additive_effects.R:533 | DAE | Other-parent-scope generated terms at the target's line are now stale against the new target. | Yes. Refusing would deadlock (per the code comment), and `define_phenotype(prevalence=)` re-checks it by refusing two parent scopes at one line. |
| define_additive_effects.R:855 | `.dae_warn_parent_only` (DAE) | Two additive variants that differ only in parent (a legal fallback pair). | Yes, warn. |
| define_additive_effects.R:1377; define_genome_effects.R:972 | diagnostics (DAE, DGE) | The realised spectrum seen by another population falls outside `warn_bounds`. | Yes, warn. Its own off switch (`warn_bounds = NULL`) is a second relaxation mechanism. |
| define_additive_effects.R:1464 | `.dae_default_base` (DAE, DGE) | The default founder base pools more than one line pool (Wahlund). | Yes, warn. Arguably refuse, because calibration is to pooled heterozygosity, but the user can pass `base_tbl` on purpose. |
| genome_effects_eval.R:669 | `.gev_preflight` | Label-vector enumeration above `tidybreed.label_vector_warn`. | Resource warning; fine. |
| add_tgv.R:197 | `add_tgv` | No individuals matched; nothing computed. | Fine (a no-op). |
| add_tgv.R:369 | `.tgv_warn_structural_additive` (`add_tgv(index_names=, component_name="additive")`) | A true index on `additive` for traits with hand-written indicator or interaction terms, which can give an all-zero or partial index. | **Borderline.** It is literally what was asked, but the result would be silently wrong as a breeding value. Candidate for refuse with `allow = "structural_additive_index"`. |
| qtl_congruence.R:642 | `.qtl_sim_effects` (test-only) | The realised variance departs from `G`. | Dead path in R/. |

Also test-only (`.qtl_sim_effects()`, no caller in R/): qtl_congruence.R:460, 475, 486 (the dual anchor is infeasible) and 613 (the dual basis is rank deficient).

### 3. Notes

**Line-specific chromosome rules are ignored by every ploidy check here (gap).**
`assert_qtl_autosomal()`, `.ge_validate_dominance_ploidy()` and
`.ge_reachable_copy_counts()` (copy-count inference and `require_complete`) all
call `resolve_chr_inheritance(conn, sex)` with `line_name = NULL`. That resolves
only rows with `line_name IS NULL`. After v0.76.12 `add_offspring()` honours
line-specific rules, so a line with a haploid chromosome passes the generators'
autosomal check (miscalibration) and the dominance writer check. It is caught
only later by eval :770. `.egv_diploid_loci()` already intersects over all
lines, and these checks should do the same.

**`define_chromosome()` does not re-run `validate_genome_effects()`.** A rule
added after a dominance term exists is caught only at evaluation (eval :770).
That is the re-check pattern the policy asks for, but the error then blames
haplotypes rather than the new rule.

**A shared `allow` needs a shared name in two places:**
- `no_matching_copy`: eval :917 (`add_tgv`/AP) and EGV :470.
- `non_diploid_dominance`: helpers :348 (write) and eval :770 (evaluate).

Without stored flags, an allow given at write time cannot reach `add_tgv()`
(often called from `add_phenotype()`). So the write must store a term whose
evaluation is defined. One option is an explicit origin, or a contrast whose
value at copy count ≠ 2 is defined as 0. Otherwise the evaluator re-check would
refuse again.

**`ignore_stored_target` duplicates an existing route.** Filtering
`trait_var_comp_tbl` to the additive rows is already the named opt-out, so the
allow would only be sugar for that filter.

**Trust dependencies:**
- `define_phenotype(prevalence = )` trusts every 1-Q21 / Q21 row above
  (245, 939, 777, 984/694, 996, 1093, 1119 and the calibration failures). It
  re-checks that every term is `"generated"` and refuses two parent scopes at
  one line. That re-check is what makes the :533 and :855 warnings safe.
- `.gev_select_variant` (:608) re-checks the write-time incomparable-scope rule.
- EGV re-derives decomposability itself (`.egv_classify` / `.egv_diploid_loci`)
  rather than trusting the writers.
- `.stored_to_functional` internal errors trust `.egv_classify`.
- The generators' calibration trusts `assert_qtl_autosomal()`, which has the
  line gap above.

**Possibly miscategorised today:**
- `require_complete = FALSE` default: should be refuse-by-default with an allow
  (the policy example is "a missing genotype state").
- `method = "union"` (:497): it writes `"generated"` terms whose covariances miss
  the stored block. The diagonals are exact, so the prevalence threshold, which
  is per trait, is still correct. But any consumer that later trusts a stored
  covariance for generated terms would be wrong. Either refuse union under a
  non-zero stored covariance, or document the exception to Q21.
- The :360 warning should be a refusal.
- add_tgv.R:369 is borderline (see the warnings table).

**Two existing relaxation mechanisms sit outside `allow`:**
`options(tidybreed.label_vector_max / _warn)` and `warn_bounds = NULL`. Decide
whether these become `allow` names or stay as resource and verbosity knobs.

**Hard to classify:**
- Infeasible-target refusals (1-Q21): they are not assumptions; they are math
  with no solution.
- Resource guards.
- remove_generated_effects.R:95 ("no generated effects at this scope; nothing
  removed") is a no-op refusal. I counted it as category 0; it could be a
  warning.
- "Every coefficient is zero, nothing to write" in the builders (:124, :267,
  :435) is also a no-op refusal (category 0).
- An explicit indicator `copy_count_value` that is unreachable on its
  chromosome is **not** checked: the term silently never fires. That would be a
  category 3 warn candidate (`unreachable_genotype_state`).

---

## Check inventory — phenotype, effects, index, EBV files

Repo state: main @ 6136b89 (v0.76.12). Line numbers are the line of the `stop(` / `warning(` call.
Counts are `stop(` calls only (205 in these 19 files); `stopifnot()` / `match.arg()` are all category 0 and not counted.
Warnings: 19.

Totals: **cat 0 = 134, cat 1 = 23, cat 2 = 26, cat 3 = 22** (of which 5 already have a user-choice argument).

### 1. Per-file category 0 counts

| file | cat 0 | cat 1 | cat 2 | cat 3 | warnings |
|---|---|---|---|---|---|
| R/add_phenotype_stages.R | 14 | 2 | 5 | 9 | 6 |
| R/add_phenotype.R | 5 | 0 | 1 | 1 | 1 |
| R/define_phenotype.R | 18 | 1 | 3 | 1 | 0 |
| R/phenotype_cov_block.R | 8 | 2 | 6 | 3 | 0 |
| R/phenotype_helpers.R | 2 | 1 | 1 | 2 | 0 |
| R/correlated_draws.R | 14 | 7 | 2 | 1 | 0 |
| R/formula_helpers.R | 10 | 7 | 1 | 1 | 2 |
| R/define_effect_cov_matrix.R | 12 | 2 | 2 | 2 | 0 |
| R/mutate_derived.R | 13 | 0 | 2 | 0 | 0 |
| R/define_effect_random.R | 2 | 1 | 0 | 0 | 0 |
| R/define_effect_fixed_cov.R | 3 | 0 | 0 | 0 | 0 |
| R/define_effect_fixed_class.R | 2 | 0 | 0 | 0 | 0 |
| R/effect_helpers.R | 1 | 0 | 0 | 1 | 0 |
| R/define_residual_cov.R | 2 | 0 | 0 | 0 | 0 |
| R/define_condition_change_action.R | 2 | 0 | 0 | 0 | 0 |
| R/define_index.R | 2 | 0 | 0 | 0 | 0 |
| R/add_index.R | 5 | 0 | 1 | 1 | 1 |
| R/add_ebv.R | 10 | 0 | 1 | 0 | 8 |
| R/blupf90_helpers.R | 9 | 0 | 1 | 0 | 1 |

Entry-point shorthand: AP = `add_phenotype()`, DP = `define_phenotype()`, DRC = `define_residual_cov()`,
DECM = `define_effect_cov_matrix()`, DER = `define_effect_random()`, DEFC = `define_effect_fixed_class()`,
DEFV = `define_effect_fixed_cov()`, AEBV = `add_ebv()`, AI = `add_index()`, MD = `mutate_derived()`.

### 2. Category 1, 2, 3 checks and all warnings

### Category 3 — modelling assumptions (candidates for `allow`)

| file:line | function (entry) | refuses | cat | proposed `allow` | default | downstream behavior needed if bypassed | trusted by |
|---|---|---|---|---|---|---|---|
| add_phenotype_stages.R:154 | `.ap_plan` (AP) | simple phenotype whose trait has no genome-effect terms | 3 | `phenotype_without_genetic_effects` | refuse (usually a forgotten generator call) | TGV = 0 for every record: `add_tgv()` writes no `ind_tgv` rows, so `.tgv_read()`/`.tgv_by_id()` must read "no rows" as 0, not NA (today NA would drop/NA the record). Prevalence path re-checks via :1174. | `.ap_materialize_tgvs`, `.ap_prevalence_genetic_var` (:1174 is unreachable today because this runs first) |
| add_phenotype_stages.R:462 | `.ap_plan_phenotype` (AP) | formula_tgv: individuals with a missing contributor (only when `missing_component_action = "error"`) | 3 | existing arg `missing_component_action` | skip+warn (existing) | already defined: excluded, no record | — |
| add_phenotype.R:380 | `.assemble_composite_tgv` (AP) | components: individuals with a missing contributor (only when `"error"`) | 3 | existing `missing_component_action` | skip+warn | already defined | — |
| add_phenotype_stages.R:1059 | `.ap_residual_block` (AP) | stored residual of a block member drawn under a different condition stratum than the new record resolves to | 3 | existing arg `condition_change_action` (`"independent"`) | refuse (existing default `"error"`) | already defined: mismatched stored coordinates dropped from conditioning, warning :1051 | — |
| add_phenotype_stages.R:1174 | `.ap_prevalence_genetic_var` (AP) | prevalence on a trait with no genome effects | 3 (dead code) | (fold into `phenotype_without_genetic_effects`) | refuse | Vg = 0 is well defined; threshold from Vr+Ve | unreachable: :154 fires first for the same state |
| add_phenotype_stages.R:1180 | `.ap_prevalence_genetic_var` (AP) | prevalence when any term is not owner `"generated"` (user terms from `define_genome_effect_terms()`) | 3 | `prevalence_from_stored_target` | refuse (silently wrong prevalence) | threshold = mean + z·sqrt(stored target + Vr + Ve) even though target does not describe the terms; realized prevalence ≠ requested. Defined math, no later function reads the threshold. | relies on Q21 (DECM :469, generators always calibrate) |
| add_phenotype_stages.R:1193 | same | prevalence when model has terms outside additive / dominance / A×A (indicator surfaces, other interactions) | 3 | `prevalence_from_stored_target` | refuse | same: ignore unmodelled variance | — |
| add_phenotype_stages.R:1219 | same | prevalence when one kind has generated variants for >1 parent-of-origin scope at one line | 3 | `prevalence_from_stored_target` | refuse | same: use one target although real variance is a sum | — |
| add_phenotype_stages.R:1264 | `.ap_prevalence_env_var` (AP) | prevalence when the residual has conditional strata | 3 | `prevalence_from_stored_target` (or `prevalence_unconditional_residual`) | refuse | use the unconditional stratum as Ve (exists? if not, :1342 still refuses — cat 2) | — |
| add_phenotype_stages.R:1282 | same | prevalence with a gamma random effect (non-centred, mean sqrt(var)) | 3 | `prevalence_from_stored_target` | refuse | add its variance, ignore its mean shift | — |
| phenotype_cov_block.R:144 | `.pvc_check_matrix` (DP residual_var, DRC, DECM residual/named, DER variance) | covariance block not positive semi-definite | 3 | `non_psd_covariance` | refuse | **none today**: `resolve_correlated_draws` re-checks at correlated_draws.R:278 and `.cd_factor` :424, so a bypass would just move the failure to AP. Needs a defined projection (eigen-clamp to PSD, sign-fixed for determinism) in `resolve_correlated_draws`, and BLUPF90 R would need the same projection. Until then effectively cat 2. | correlated_draws.R:278/:424 (re-check), `.blupf90_residual_cov` (no re-check) |
| correlated_draws.R:278 | `resolve_correlated_draws` (AP) | the stored block is not PSD at draw time | 3 (re-check) | same `non_psd_covariance` | refuse | same projection | — (this is the downstream re-check of :144) |
| phenotype_cov_block.R:239 | `.pvc_check_realization_lock` (DP residual_var, DRC, DECM "residual") | redefine a residual block that already has realized `residual_value`s | 3 | `redefine_realized_block` | refuse (silently mixes two models) | AP already handles it mathematically: new records condition on stored residuals under the new R (`resolve_correlated_draws`). Needs: (a) if the new R is singular, stored values may be off its support → correlated_draws.R:379 refuses (fine, cat 2); (b) old records keep values drawn under old R, nothing records which R — acceptable only because no stored provenance is promised. | `.ap_residual_block` trusts stored residuals were drawn under the current R |
| phenotype_cov_block.R:257 | same (DECM named, DER variance) | redefine a named random-effect block with realized levels in `phenotype_random_effects` | 3 | `redefine_realized_block` | refuse | `.ap_named_effect_block` reuses stored level values and conditions new coordinates on them under the new matrix; same caveats. Note `overwrite = TRUE` in DER already bypasses this by deleting the draws (see notes). | `.ap_named_effect_block` |
| phenotype_helpers.R:139 | `.ap_covariate_terms` (AP) | NULL class in a fixed_class source column (only `null_class_action = "error"`) | 3 | existing arg `null_class_action` (default `"skip"`) | skip+warn | defined (skip / zero) | — |
| phenotype_helpers.R:161 | `.ap_covariate_terms` (AP) | a class level present in data but absent from `levels` | 3 | `unlisted_class_level` | refuse (typo in a level = silently wrong) | shift 0 (reference level) — the very next line (`shifts[is.na(shifts)] <- 0`) already implements it | — |
| formula_helpers.R:437 | `.eval_formula_tgv` (AP) | formula_tgv producing Inf/-Inf/NaN (div by 0, domain) | 3 | `non_finite_formula_value` | refuse | treat as missing → route through `missing_component_action` (exclude). Must not write Inf into `ind_tgv`/liability. | — |
| define_effect_cov_matrix.R:364 | `.tvc_write_block` (DECM genetic; generators' `G =`) | writing a genetic target where one is already stored (even identical) | 3 | `replace_stored_target` | refuse | replace rows in place. Safe only where :469 also passes (no generated terms at that scope). Generators read target at next run; BLUPF90 reads it as G prior; prevalence reads it but :469 guards. | `define_additive_effects`/`define_genome_effects` (calibration), `.ap_prevalence_genetic_var`, `write_renum_par` |
| define_effect_cov_matrix.R:469 | `.tvc_refuse_under_generated` (DECM) | a genetic target under existing generated terms of that kind/scope (Q21) | 3 → **keep unbypassable for now** | (`target_under_generated_terms`) | refuse | Bypass leaves "generated" terms not calibrated to the stored target; `.ap_prevalence_genetic_var` checks only the owner (:1180) and would trust it → silently wrong. Per the "re-check, no flags" rule, prevalence would have to verify calibration itself (e.g. recompute the model's variance from terms + base allele freqs and compare to target). Until then this is effectively cat 1. | `.ap_prevalence_genetic_var` (trusts owner = calibrated), `write_renum_par` G prior |
| define_phenotype.R:457 | `define_phenotype` | redefine an existing phenotype | 3 | existing arg `overwrite` | refuse | defined (replace meta + components). Gap: no check for existing `ind_phenotype` records (type/thresholds/mean change after realization is not locked, unlike cov blocks). | — |
| effect_helpers.R:38 | `.handle_effect_overwrite` (DER, DEFC, DEFV) | redefine an existing phenotype effect | 3 | existing arg `overwrite` | refuse | defined; but it deletes this phenotype's `phenotype_random_effects` draws (see notes) | — |
| add_index.R:248 | `add_index` | some individuals lack a value for an index trait | 3 | `index_missing_values` | refuse (partial index is silently wrong) | write `index_value = NA` for those individuals (or omit them); downstream selection on `ind_index` must accept NA | selection filters on `ind_index` |

### Category 1 — integrity / safety (never bypassable)

| file:line | function (entry) | refuses | cat | notes |
|---|---|---|---|---|
| add_phenotype_stages.R:719 | `.ap_named_effect_block` (AP) | a named-effect block carrying conditional strata (writers never store one; hand-edited) | 1 | defensive re-check of stored invariant |
| add_phenotype_stages.R:1026 | `.ap_residual_block` (AP) | duplicate `(id_ind, phenotype_name, pheno_number)` rows with residuals | 1 | key violation; conditioning undefined |
| define_phenotype.R:374 | `define_phenotype` | commas in `cat_names` | 1 | stored comma-joined; would corrupt decoding |
| phenotype_cov_block.R:166 | `.pvc_check_complete` (DP, DRC, DECM, DER) | partial redeclaration of a block (omits members already joined) | 1 | CLAUDE.md hard rule "declared whole". Could be argued cat 3 (fill omitted cross-covariances with 0) but that can make the block non-PSD; keep 1 |
| phenotype_cov_block.R:301 | `.check_condition_change_agreement` (DP, DRC, AP) | block members disagree on `condition_change_action` | 1 | block-scoped value; disagreement makes the action undefined. Re-checked at AP (:1034) — good example of re-check |
| phenotype_helpers.R:12 | `encode_levels_json` (DEFC) | level names containing `"` `,` `:` | 1 | hand-rolled JSON encoding would corrupt |
| correlated_draws.R:107, :113, :123, :127, :136, :140, :144 | `.cd_build_block` via `find_covariance_blocks` (AP, AEBV/BLUPF90) | stored block rows: NULL condition fields, >1 condition column, duplicate rows, incomplete, non-finite, not bit-symmetric | 1 | read-time re-checks of writer invariants (7 checks) |
| formula_helpers.R:150, :159, :168 | `.check_derived_formula` (DP, AP re-check before eval) | derived `formula` symbol that is an operator/whitelist name, disallowed call, non-numeric constant | 1 | closed grammar (eval security) |
| formula_helpers.R:286 | `dsl_call`→`bad()` in `.walk_formula_tgv_ast` (DP, AP) | malformed contributor-function arguments in formula_tgv | 1 | grammar of the DSL; one `stop` serves many messages |
| formula_helpers.R:341, :351, :362 | `.walk_formula_tgv_ast` (DP, AP) | function name used as trait, disallowed call, non-numeric constant | 1 | closed grammar |
| define_effect_cov_matrix.R:215 | `.check_effect_name_input` (DECM, DER, DEFC, DEFV, `write_phenotype_cov_block`) | defining a derived output name (`total`, …) | 1 | reserved vocabulary |
| define_effect_cov_matrix.R:219 | same | genetic name used for a phenotype-level effect | 1 | reserved vocabulary |
| define_effect_random.R:97 | `define_effect_random` | `effect_name = "residual"` | 1 | reserved; residual has its own writer |

### Category 2 — not implemented

| file:line | function (entry) | refuses | cat | notes |
|---|---|---|---|---|
| add_phenotype_stages.R:995 | `.ap_residual_block` (AP) | records resolving to no residual stratum when the block has no unconditional fallback | 2 | no variance to draw from |
| add_phenotype_stages.R:1230 | `.ap_prevalence_genetic_var` (AP) | prevalence when a term kind has no population-wide target row | 2 | no number; could be cat 3 if "missing target = 0" were defined, but that is silently wrong |
| add_phenotype_stages.R:1314 | `.ap_check_prevalence` (AP) | prevalence on composite (components/formula_tgv) phenotype | 2 | duplicate of define_phenotype.R:356 for rows defined before that check existed / hand edits |
| add_phenotype_stages.R:1342 | `.ap_liability_records` (AP) | prevalence with no unconditional residual (reachable when `user_residual` supplies residuals for a phenotype with no stored block) | 2 | |
| add_phenotype_stages.R:1354 | `.ap_liability_records` (AP) | prevalence with total liability variance 0 | 2 | mathematically undefined |
| add_phenotype.R:358 | `.assemble_composite_tgv` (AP) | stored `weight_type` other than fixed/covariate | 2 | re-check of define_phenotype.R:518 |
| define_phenotype.R:356 | `define_phenotype` | `prevalence` on a composite phenotype | 2 | |
| define_phenotype.R:401 | `define_phenotype` | `residual_var` on a `derived_formula` phenotype | 2 | a derived phenotype + measurement error is a real model; no path today |
| define_phenotype.R:518 | `define_phenotype` | `weight_type` other than fixed/covariate | 2 | |
| phenotype_cov_block.R:191 | `.pvc_check_strata` (DRC) | a second condition column on one residual block | 2 | |
| phenotype_cov_block.R:214 | `.pvc_check_strata` (DRC) | a stratum declared over a different phenotype set than the block | 2 | |
| phenotype_cov_block.R:350 | `validate_named_effect_block` (DECM, DER, AP backstop) | a named block joining non-random effects | 2 | |
| phenotype_cov_block.R:362 | same | non-normal distribution in a multi-member named block | 2 | |
| phenotype_cov_block.R:373 | same | block members reading different grouping columns | 2 | correlated levels on different factors undefined |
| phenotype_cov_block.R:412 | `validate_phenotype_cov_block` | condition strata on a non-residual effect | 2 | |
| phenotype_helpers.R:124 | `.ap_covariate_terms` (AP) | source table with duplicate `id_ind` rows for the subset | 2 | no rule picks a row (could also be read as cat 0 bad data) |
| correlated_draws.R:379 | `resolve_correlated_draws` (AP) | observed/fixed values off the support of a singular covariance | 2 | reachable with `user_residual` on a singular block, or after a relaxed realization lock |
| correlated_draws.R:424 | `.cd_factor` (AP) | conditional covariance materially indefinite | 2 | numerical; should be unreachable when the block is PSD |
| formula_helpers.R:634 | `.topo_sort_phenotypes` (AP) | cyclic derived_formula dependencies | 2 | |
| define_effect_cov_matrix.R:142 | `define_effect_cov_matrix` | `line_name` on a non-genetic effect | 2 | line-specific phenotype effects not implemented (schema bias says reserve the dimension) |
| define_effect_cov_matrix.R:210 | `.check_effect_name_input` | `GENETIC_EFFECT_NAMES_FUTURE` names | 2 | |
| mutate_derived.R:147 | `mutate_derived` | destination table not in `TABLE_PRIMARY_KEYS` | 2 | no key to write back by |
| mutate_derived.R:276 | `mutate_derived` | same key maps to different computed values | 2 | ambiguous write-back |
| add_index.R:215 | `add_index` | >1 row per (id_ind, trait) | 2 | no rule picks the row |
| add_ebv.R:197 | `add_ebv` | software other than blupf90 | 2 | |
| blupf90_helpers.R:279 | `.blupf90_residual_cov` (AEBV) | residual block with only conditional strata | 2 | BLUPF90 takes one R |

### Warnings

| file:line | function (entry) | warns about | right level? |
|---|---|---|---|
| add_phenotype_stages.R:112 | `.ap_plan` (AP) | filter matched no individuals; nothing written | yes |
| add_phenotype_stages.R:202 | `.ap_plan` (AP) | non-repeatable phenotype: already-phenotyped individuals skipped | yes (literally what `repeatable = FALSE` means; could be a message) |
| add_phenotype_stages.R:434 | `.ap_plan_phenotype` (AP) | individuals skipped by `null_class_action = "skip"` | yes (user-chosen action, mirrors missing_component_action) |
| add_phenotype_stages.R:462 | `.ap_plan_phenotype` (AP) | formula_tgv missing component → excluded (`skip`) | yes (precedent) |
| add_phenotype.R:381 | `.assemble_composite_tgv` (AP) | components missing → excluded (`skip`) | yes |
| add_phenotype_stages.R:1006 | `.ap_residual_block` (AP) | condition level with no stratum → drawn from unconditional R | yes (documented fallback; literally asked) |
| add_phenotype_stages.R:1051 | `.ap_residual_block` (AP) | `condition_change_action = "independent"` dropped coordinates | yes |
| formula_helpers.R:217 | `.validate_derived_formula` (DP) | derived formula names phenotypes not yet defined | ok at define time, **but the promised AP-time hard error is incomplete**: formula_helpers.R:537 fires only if *no* symbol has records; an undefined symbol next to a defined one becomes an all-NA column (`.pivot_pheno_wide`) and silently writes NA records. AP should refuse unknown symbols (cat 0). |
| formula_helpers.R:558 | `.eval_derived_formula` (AP) | derived formula Inf/NaN → NA, records written with NA | **inconsistent** with formula_helpers.R:437 (formula_tgv refuses). Make both the same check (`non_finite_formula_value`), default refuse. |
| add_index.R:192 | `add_index` | no filter applied to the table | redundant: :215 hard-refuses the only harmful case. Drop it (or message). |
| add_ebv.R:184 | `add_ebv` | `phenotype` filter ignored in parent_avg mode | too weak: an argument silently ignored → cat 0 refuse |
| add_ebv.R:191 | `add_ebv` | `phenotype` filter matched no records | redundant: blupf90_helpers.R:132 then refuses; make it the refusal |
| add_ebv.R:229 | `add_ebv` | no candidates matched | yes |
| add_ebv.R:330 | `ebv_parent_avg` (AEBV) | parents have several evaluations, latest used | modelling assumption that is the normal case after repeated evaluations → message, or an explicit eval selector |
| add_ebv.R:338 | `ebv_parent_avg` (AEBV) | no parents; all EBVs NA | yes |
| add_ebv.R:369 | `ebv_parent_avg` (AEBV) | one warning **per individual × trait** with a missing parent EBV | wrong shape: aggregate into one counted warning (R truncates >50) |
| add_ebv.R:425, :428 | `ebv_blupf90` (AEBV) | `eval_id` / `run_dir` ignored under managed folders | arguably cat 0 refuse (ignored argument); warn acceptable |
| blupf90_helpers.R:567 | `parse_blupf90_solutions` (AEBV) | no solutions matched pedigree IDs; returns empty | too weak: AEBV then writes nothing silently — should refuse (external failure) |

### 3. Notes

**Hard to classify**
- Missing variance for a random effect / residual (add_phenotype_stages.R:675, :861; define_effect_random.R:143; :1289) — classed cat 0 (missing config). It is *not* "treat as 0": the user can declare 0 explicitly, so no `allow` is needed.
- phenotype_cov_block.R:144 PSD is textbook cat 3 but has no downstream behavior today (two re-checks at draw time). It is only bypassable once a deterministic PSD projection exists in `resolve_correlated_draws` **and** BLUPF90's R uses the same projection. Tolerances differ between writer (`tol` 1e-8/1e-9 relative to max eigenvalue) and reader (`n·sqrt(eps)·max|cov|`); a near-singular block accepted by one could be refused by the other.
- define_effect_cov_matrix.R:469 (Q21) — cat 3 in spirit but must stay unbypassable until `.ap_prevalence_genetic_var` can verify calibration without trusting the owner. This is the clearest case of "later code trusts an earlier check without re-checking".
- add_phenotype_stages.R:1174 is dead (:154 fires first for the same state).
- correlated_draws.R:379 is cat 2, but it is the landing point for two cat-3 relaxations (`redefine_realized_block`, user_residual on a singular block); leave it as the backstop.
- formula_helpers.R:537 (derived formula, no feeder records) — classed cat 0 (ordering mistake), but it is inconsistent with the partial case which silently writes NA (see warning row :217).

**Shared `allow` names**
- `prevalence_from_stored_target`: add_phenotype_stages.R:1180, :1193, :1219, :1264, :1282 — all mean "place the threshold from stored variances even though they do not describe this liability". Keep :1230, :1342, :1354, :1314 refused (no number exists). Note define_phenotype() itself does none of these checks; they all fire at AP time, so the `allow` would go on `add_phenotype()` (or be re-checked there if given on DP).
- `redefine_realized_block`: phenotype_cov_block.R:239 and :257 (residual and named).
- `non_psd_covariance`: phenotype_cov_block.R:144 and correlated_draws.R:278 (the re-check must honour the same name, passed through AP — no stored flag).
- `non_finite_formula_value`: formula_helpers.R:437 (refuse today) and :558 (warn today) — unify.
- `phenotype_without_genetic_effects`: add_phenotype_stages.R:154 (+ dead :1174).

**Trust dependencies (X assumes Y passed)**
- `.ap_prevalence_genetic_var` trusts that every `"generated"` term was calibrated to the stored target (generators + DECM :469 + `.tvc_write_block` :364). It re-checks owner, kinds and parent scopes but cannot re-check calibration.
- `write_renum_par` (blupf90_helpers.R:326) uses the population-wide `additive` target as BLUPF90's G prior, with **no** owner/calibration re-check — same trust as prevalence, unguarded. It also pulls residual blocks keyed by *phenotype* names using the *trait* names (assumes simple phenotype = trait name).
- `.ap_residual_block` / `.ap_named_effect_block` trust stored realizations were drawn under the current block (realization lock :239/:257). `.handle_effect_overwrite` (effect_helpers.R:44-51) deletes `phenotype_random_effects` for **only this phenotype** when `overwrite = TRUE`, then DER's `.pvc_write_block` lock check passes — an existing bypass of the lock. In a multi-member block, other members' draws stay while this member's are gone.
- AP re-runs `validate_named_effect_block` (:726) and `.check_condition_change_agreement` (:1034), and `.check_derived_formula` / `.walk_formula_tgv_ast` before every eval — good re-check precedent.
- `.cd_build_block` (correlated_draws.R:107-144) re-checks every writer invariant of `phenotype_var_comp` on read (cat 1 backstops).

**Existing user-choice arguments that already act like `allow`**
- `missing_component_action` (skip|error), `null_class_action` (skip|error|zero), `condition_change_action` (error|independent): these choose *between defined behaviors*, stored per phenotype/effect, and are model configuration rather than check relaxation. Recommend they **stay** as named arguments (they pick a behavior, not just silence a refusal). `null_class_action = "zero"` and `condition_change_action = "independent"` have no `allow` equivalent.
- `overwrite` (define_phenotype.R:457, effect_helpers.R:38): a replace-vs-refuse switch; keep as is, but DP overwrite has no realization lock (type/thresholds/mean can change after records exist) while cov blocks do — inconsistent. Either lock DP overwrite when `ind_phenotype` has records (with `allow = "redefine_realized_phenotype"`) or document why not.
- `repeatable = FALSE` (warning :202) also behaves as a skip-and-warn action.

**Possibly miscategorised today**
- formula_helpers.R:558 should refuse (consistent with :437), not warn-and-write-NA.
- blupf90_helpers.R:567, add_ebv.R:184 and :191 should be refusals.
- add_ebv.R:369 should be one aggregated warning; add_ebv.R:330 a message.
- add_index.R:192 is noise given :215.
- define_phenotype.R:401 (`residual_var` on derived_formula) is a genuine missing model (cat 2), not a malformed input.
- add_phenotype_stages.R:154's message points users to `type = "derived_formula"` for a phenotype "with no genetic architecture", but derived formulas only combine other phenotypes; a purely environmental phenotype has no path today — this is what `phenotype_without_genetic_effects` would provide.

---

## Check inventory: genome / chromosome / breeding files

Scope: R/add_offspring.R, add_founders.R, define_chromosome.R, chr_meta_helpers.R,
define_genome.R, define_founder_haplotypes.R, founder_haplotype_helpers.R,
genome_map_helpers.R, recombination_helpers.R, haplotype_write_helpers.R,
ploidy_helpers.R, extract_genotypes.R, extract_allele_freq.R, define_chip.R.
Line numbers are from HEAD 6136b89. Counted: every `stop()` and `stopifnot()`
written in the file (calls to shared validators such as `validate_tidybreed_pop()`,
`.check_scalar()` and `validate_sql_identifier()` are not counted again).

### 1. Category 0 counts per file

| File | cat 0 | cat 1/2/3 + warnings listed below |
|---|---|---|
| add_offspring.R | 10 | 6 |
| add_founders.R | 9 | 2 |
| define_chromosome.R | 13 | 2 |
| chr_meta_helpers.R | 0 | 12 |
| define_genome.R | 9 | 2 |
| define_founder_haplotypes.R | 7 | 4 stops + 2 warnings |
| founder_haplotype_helpers.R | 0 | 1 |
| genome_map_helpers.R | 0 | 8 |
| recombination_helpers.R | 0 | 2 |
| haplotype_write_helpers.R | 3 | 0 |
| ploidy_helpers.R | 0 | 1 |
| extract_genotypes.R | 17 | 0 (calls `assert_ploidy_2`, listed under ploidy_helpers.R) |
| extract_allele_freq.R | 7 | 1 |
| define_chip.R | 5 | 1 |
| **Total** | **80** | **42 stops + 2 warnings** |

By category (listed rows): cat 1 = 28, cat 2 = 11 (incl. assert_qtl_autosomal and define_chip:90), cat 3 = 3 (chr_meta_helpers:159, define_founder_haplotypes:385, recombination_helpers:88), warnings = 2. Five more cat 3 checks are proposed in the notes.

### 2. Category 1, 2, 3 checks and warnings

Columns: location | function (exported entry point) | refuses | cat | allow name (cat 3) | default | behaviour needed if bypassed | trusted by.

| file:line | function (entry) | what it refuses | cat | allow | default | downstream behaviour if bypassed | trusted by |
|---|---|---|---|---|---|---|---|
| add_offspring.R:401 | add_offspring | parent's autosomal `ind_haplotype` row count != 2 x autosomal loci | 1 | – | – | – (stored haplotypes disagree with rules) | kernel packing at 406-415 |
| add_offspring.R:435 | add_offspring | parent's row count on a special chromosome != copies its (sex, line) rule implies | 1 | – | – | – | special path 736-784 |
| add_offspring.R:483 | add_offspring | > ~5.4e8 offspring in one call (int32 gamete stream id; R would give NA, C++ wraps) | 2 | – | – | bypass would break R/C++ parity | `.gamete_stream_id`, C++ kernel |
| add_offspring.R:567 | add_offspring (`stopifnot`) | offspring ploidy != 2 | 2 | – | – | – | everything downstream (strand always 1) |
| add_offspring.R:724 | add_offspring | offspring rule takes > 1 copy from one parent (`2,0`, uniparental disomy) | 2 | – | – | needs a transmission rule (which two copies, recombined?) | – |
| add_offspring.R:740 | add_offspring | the parent that must pass a chromosome carries 0 copies (e.g. a female in `id_parent_1` for a male offspring that takes Y from parent 1; or a cross whose lines have incompatible karyotype rules) | 1 | – | – | writing the offspring without that copy would break its karyotype vs the rules | – (root cause is a missing parent-role check, see notes) |
| add_founders.R:114 | add_founders | `ploidy != 2` | 2 | – | – | – | `assert_ploidy_2` (makes it unreachable), add_offspring:567 |
| add_founders.R:312 | add_founders | a founder rule needing > 1 copy from one parent | 2 | – | – | no founder mechanism (1 pool draw per slot) | – |
| define_chromosome.R:201 | define_chromosome | `from_parent_1 + from_parent_2 > 2` | 2 | – | – | – (also a SQL CHECK, define_genome.R:420) | validate_chr_inheritance:227 |
| define_chromosome.R:352 | define_chromosome | inheritance rule that changes the copies existing individuals carry | 1 | – | – | – | add_offspring:435/740 (re-checks row counts) |
| chr_meta_helpers.R:55 | resolve_chr_inheritance (define_chromosome, define_genome, add_founders, add_offspring, define_additive_effects, define_genome_effects, define_genome_effect_terms, extract_genetic_variance) | a chromosome has no inheritance rule for the requested (sex, line) | 1 | – | – | – | every caller |
| chr_meta_helpers.R:116 | resolve_chr_recombination (same writers + add_founders/add_offspring) | a chromosome has no recombination rule | 1 | – | – | – | every caller |
| chr_meta_helpers.R:159 | .chr_shadow_check (define_chromosome via validators) | a (sex, NULL-line) rule and a (NULL-sex, line) rule differ with no explicit (sex, line) row | 3 | `ambiguous_chr_rule` | refuse | resolvers already define it: sex-specific (prio 2) beats line-specific (prio 1). Bypass is safe only if `.chr_shadow_check` stays out of every later path (it does: only validators call it) | none at run time; resolvers apply precedence |
| chr_meta_helpers.R:196 | validate_chr_inheritance | duplicate logical key | 1 | – | – | – | resolvers (ROW_NUMBER picks one arbitrarily otherwise -> non-determinism) |
| chr_meta_helpers.R:207 | validate_chr_inheritance | rule names a chromosome not in genome_meta | 1 | – | – | – | – |
| chr_meta_helpers.R:216 | validate_chr_inheritance | offspring_sex not NULL/M/F (backstop; SQL CHECK too) | 1 | – | – | – | – |
| chr_meta_helpers.R:227 | validate_chr_inheritance | negative counts or sum > 2 | 2 | – | – | – | add_founders/add_offspring copy logic |
| chr_meta_helpers.R:268 | validate_chr_recombination | duplicate logical key | 1 | – | – | – | resolvers |
| chr_meta_helpers.R:278 | validate_chr_recombination | orphan chr_name | 1 | – | – | – | – |
| chr_meta_helpers.R:286 | validate_chr_recombination | invalid parent_sex (backstop) | 1 | – | – | – | – |
| chr_meta_helpers.R:293 | validate_chr_recombination | NULL recombines (backstop; NOT NULL in DDL) | 1 | – | – | – | – |
| chr_meta_helpers.R:388 | assert_qtl_autosomal (define_additive_effects:451, define_genome_effects:432) | generated QTL on a chromosome that is not 1,1 for both sexes (line-agnostic resolution only) | 2 (see notes; other agent covers) | – | – | would need a hemizygous calibration formula (`p(1-p)a^2`) to keep Q21's "generated = calibrated" | prevalence trust in define_phenotype (Q21) |
| define_genome.R:163 | define_genome | `n_loci < n_chr` (chromosome with zero loci) | 2 | – | – | empty chromosome has no genome_meta rows, so chr tables would orphan | validators |
| define_genome.R:192 | define_genome | genome already defined | 1 | – | – | – | – |
| define_founder_haplotypes.R:281 | define_founder_haplotypes | founder_haplotypes has no `line_name` column ("old format") | 1 (legacy) | – | – | – | – (should be deleted, see notes) |
| define_founder_haplotypes.R:293 | define_founder_haplotypes | a NULL-line pool already exists | 1 | – | – | a second pool reuses haplotype_id 1..n under the same key; add_founders' key matrix would merge two haplotypes per key | add_founders:181-190 |
| define_founder_haplotypes.R:309 | define_founder_haplotypes | a pool for this line already exists | 1 | – | – | same | add_founders |
| define_founder_haplotypes.R:385 | define_founder_haplotypes | mosaic `n_templates > n_haplotypes` | 3 | `excess_templates` | warn | none needed: the mosaic generator is defined for any K >= 2 (unused templates only cost RNG/time) | nothing downstream |
| define_founder_haplotypes.R:423 | warning() | fixed `allele_freq` not representable on the 1/n grid | warn | – | warn (right level) | – | – |
| define_founder_haplotypes.R:455 | warning() | > 10% of loci monomorphic in an LD-method pool | warn | – | warn (right level; result is what was asked) | – | – |
| founder_haplotype_helpers.R:228 | .write_founder_haplotypes | frequency vector length != loci (internal) | 1 | – | – | – | – |
| genome_map_helpers.R:63 | resolve_genome_map (add_offspring, define_founder_haplotypes LD methods) | a locus has no map position for (sex, line) | 1 | – | – | – | kernel |
| genome_map_helpers.R:78 | resolve_genome_map | resolved pos_cM decreases within a chromosome | 1 | – | – | crossover placement assumes sorted positions | kernel, `.draw_chr_recombination` |
| genome_map_helpers.R:119/130/141/150/160/172 | validate_genome_map (define_genome only) | duplicate key / orphan locus / locus_name mismatch / bad sex or NULL map_name / non-finite pos_cM / duplicate id | 1 (x6) | – | – | – | resolve_genome_map |
| recombination_helpers.R:88 (+ src/make_gametes.cpp:38) | .draw_poisson_dqrng (add_offspring) | a chromosome longer than 3000 cM (lambda > 30) | 3 | `long_genetic_map` | refuse (almost always a units error) | the inversion sampler is exact for any lambda (O(lambda) cost); the allow must reach BOTH the R sampler and the C++ kernel or R/C++ parity breaks | – |
| recombination_helpers.R:301 (+ cpp:89) | make_gametes_batch_r | chromosomes not contiguous locus blocks (internal) | 1 | – | – | – | – |
| ploidy_helpers.R:31 | assert_ploidy_2 (extract_genotypes:189, add_dosage:77) | individuals with ploidy != 2 | 2 | – | – | – | – (unreachable via API; see notes) |
| extract_allele_freq.R:111 | extract_allele_freq (and define_additive_effects / define_genome_effect_terms `base_tbl`) | computed frequency outside [0,1] (allele not 0/1) | 1 | – | – | – | effect generators |
| define_chip.R:90 | define_chip | chip column already exists | 2 (+1) | – | – | no redefine path; redefining in place would change what existing `has_<chip>` stamps mean | add_genotypes / extract_genotypes |

### 3. Notes

### Confirmed: define_chromosome() does not re-run genome-effect checks
`define_chromosome()` (R/define_chromosome.R:298-302) runs only the chr validator
and `.dc_refuse_changed_karyotypes()`. It calls none of `assert_qtl_autosomal()`,
`.ge_validate_dominance_ploidy()` or `validate_genome_effects()`. The karyotype
guard only looks at `ind_meta`, so before founders exist a user can:
1. generate additive effects on chromosome X (passes `assert_qtl_autosomal`, X is 1,1),
2. then `define_chromosome("X", offspring_sex = "M", 0, 1)`. Allowed: no individuals yet.
Result: `"generated"` terms calibrated with `2pq` on a now-hemizygous
chromosome (Q21's "generated means calibrated" is now false, and
`define_phenotype(prevalence=)` trusts it); `dominance` members on a non-diploid
chromosome that `.ge_validate_dominance_ploidy` would have refused; indicator
`copy_count_value` inferred by `.ge_reachable_copy_counts` now stale. Fix:
`define_chromosome()` inheritance writes should call the same genome-effect
validators before COMMIT (re-check, not a flag).

### Other configuration-order gaps
- **Line-specific rules escape `assert_qtl_autosomal`.** It resolves line-agnostic only
  (chr_meta_helpers.R:372-373). A rule `define_chromosome("5", line_name = "A", ...)`
  making chr 5 hemizygous in line A is never seen by generators, even for a
  line-A-scoped effect. Same for `.ge_validate_dominance_ploidy` and
  `.ge_reachable_copy_counts` ("line-agnostic: population-level").
- **`remove_rows()` can delete `chr_inheritance`, `chr_recombination`, `genome_map`,
  `genome_meta`, `founder_haplotypes` and single `ind_haplotype` rows**
  (R/sql_utils.R:174-193 registers them; only genome_effect* are in
  `TABLE_NO_ROW_DELETE`). It runs no validator, so it bypasses
  `.dc_refuse_changed_karyotypes`, the resolvability checks and the genome-effect
  checks. Deleting a line rule changes existing animals' karyotypes; deleting a default
  rule makes resolvers error later. Should be refused (cat 1) or routed through the validators.
- **Map length is checked at mating time, not at definition time.** define_genome
  accepts `chr_len_Mb * cM_per_Mb > 3000`. The failure appears only at the first
  `add_offspring()`, after founders were built. The lambda check should also run in define_genome.
- **Rule changes after founder pools / effects but before founders** silently
  reclassify chromosomes as plain/special (`.plain_autosome_chrs`). This is documented and
  what was asked, so no check is needed, but the docs should say that adding any line rule
  changes seeded offspring of every line.
- `define_founder_haplotypes()` LD methods use the map resolved at pool time. If maps
  ever become writable after pools exist, pools would not match. This cannot happen yet
  because nothing but define_genome writes genome_map.

### Missing checks worth adding (by proposed category)
- **add_offspring parent roles (cat 3, `parent_sex_role`, refuse when a sex
  chromosome rule depends on it, otherwise warn).** Nothing checks that
  `id_parent_1` is male / `id_parent_2` female. With autosomes only this is
  harmless (hermaphrodites, plants). With XY/ZW rules it surfaces as the cryptic
  add_offspring.R:740 "misconfiguration" error, or worse, a female sire passes an X
  where a Y was meant if copy counts happen to line up. Defined behaviour if
  allowed: transmit by `parent_origin` slot as now.
- **add_founders partial pool (cat 1).** `hap_data_matrix` is initialised to 0
  (add_founders.R:188). A `founder_haplotypes` filter on `locus_name` gives founders allele
  0 at every unselected locus, silently. A pool `locus_name` missing from genome_meta
  gives an NA subscript error. Refuse unless every locus is covered for every selected haplotype.
- **add_founders mixed / foreign pools (cat 3, `founder_pool_other_line`, warn).**
  The docs (define_founder_haplotypes.R:9-12) say add_founders "will sample only from the
  pool matching its own line_name". The code samples whatever the filter returns, so an
  unfiltered multi-line table mixes pools and stamps `line_origin = line_name`. That is
  literally what was asked, so warn. Fix the docs.
- **extract_genotypes on non-diploid chromosomes (cat 3, `non_diploid_dosage`, warn).**
  Dosage is `SUM(allele)` with 0 fill (extract_genotypes.R:256-259). An absent copy
  (a female on Y, a missing MT copy) reads as dosage 0, the same value as homozygous
  reference, and hemizygous loci are on a 0/1 scale. Defined behaviour if relaxed:
  return as now. Better default: NA for absent.
- **Line-name ID collision (cat 1, a bug, verified).** The line regex allows `_`, but
  the next-id query (add_founders.R:197-200, add_offspring.R:463-469) casts everything
  after the first `_`. A second `add_founders(line_name = "A_B")` fails with a DuckDB
  conversion error, as does line `A` once `A_B_*` exists (`LEFT(id,2)='A_'` matches).
  Either forbid `_` in line names or parse the suffix after the last `_` with the exact prefix.

### Probably miscategorised / should change
- define_founder_haplotypes.R:281 ("old format" founder_haplotypes without `line_name`)
  is a legacy compatibility path, which CLAUDE.md forbids before 1.0.0. Delete it.
- `assert_ploidy_2` (ploidy_helpers.R:31) is unreachable: add_founders.R:114 and
  add_offspring.R:567 guarantee ploidy 2 and mutate_table reserves `ploidy`. Keep one
  cat 2 "polyploidy not implemented" check at the writers and drop the reader-side one,
  or keep it as the documented re-check. Either way it is the same allow-less cat 2.
- add_offspring.R:567 is a `stopifnot` with no message for a cat 2 condition; it
  should be a named cat 2 refusal like add_founders.R:114.
- define_founder_haplotypes.R:385 refuses a well-defined case. It is cat 3 and should be
  a warning.
- `assert_qtl_autosomal` reads as a cat 3 modelling check, but under Q21 a generator
  must calibrate, so bypassing means a calibration formula that does not exist yet.
  That makes it cat 2 until a hemizygous genic-variance formula is implemented, at which
  point it disappears rather than becoming an `allow`.
- define_chip.R:90: also not transactional (ALTER at 98, UPDATE at 103). A failed UPDATE
  leaves the column, and every retry then hits this refusal.

### Shared `allow` names
- `ambiguous_chr_rule`: chr_meta_helpers.R:159 for both inheritance and recombination
  (one name, both validators).
- `long_genetic_map`: R sampler and C++ kernel together, plus the proposed define_genome
  pre-check.
- `non_diploid_dosage`: extract_genotypes and add_dosage (same `SUM(allele)` semantics).

### Trust dependencies
- add_offspring trusts define_chromosome:352 for karyotype consistency but re-checks
  row counts itself (401/435/740). That is good and is the pattern the policy wants.
- The generators and define_phenotype(prevalence) trust assert_qtl_autosomal at write
  time; nothing re-checks after rule changes (gap above).
- add_founders trusts define_founder_haplotypes:293/309 for unique (line, haplotype_id)
  keys and does not re-check them.

### Hard-rule observations (not checks)
- add_offspring.R:234-238 and 365-369 paste parent `id_ind`s unescaped into SQL text.
  This violates "Individual ids never appear in SQL text", and an id containing `'`
  breaks the query. extract_genotypes.R:163/244 and assert_ploidy_2 use `sql_in_list`
  (escaped, but still ids in SQL text).
- genome_map supports sex- and line-specific maps in the resolver, but no exported
  writer exists, so validate_genome_map runs only inside define_genome.

---

## Check inventory — infrastructure files

Scope: `stop()` and `warning()` calls in R/sql_utils.R, remove_rows.R, arg_checks.R,
archive_replicate.R, mutate_table.R, mutate_group_{seq,named,concatenate,helpers}.R,
tidybreed_pop.R, define_table.R, restore_pop.R, schema.R, open_pop.R, summary_pop.R.
Line numbers are from the current tree (HEAD 6136b89). `stopifnot()` calls are not
counted (all are category 0).

### 1. Category 0 counts per file

| file | cat 0 `stop()` | lines |
|---|---|---|
| sql_utils.R | 6 | 47, 390, 436, 457, 460, 499 |
| remove_rows.R | 5 | 191, 245, 259, 274, 280 |
| arg_checks.R | 9 | 22, 25, 29, 32, 37, 44, 50, 53, 64 (all of the file) |
| archive_replicate.R | 5 | 124, 127, 133, 139, 170 |
| mutate_table.R | 3 | 482, 603, 610 |
| mutate_group_seq.R | 6 | 102, 108, 122, 128, 131, 139 |
| mutate_group_named.R | 5 | 120, 125, 134, 185, 193 |
| mutate_group_concatenate.R | 4 | 93, 97, 101, 110 |
| mutate_group_helpers.R | 0 | — |
| tidybreed_pop.R | 5 | 51, 56, 61, 266, 318 |
| define_table.R | 4 | 98, 107, 121, 124 |
| restore_pop.R | 3 | 63, 69, 83 |
| schema.R | 3 | 1033, 1282, 1497 |
| open_pop.R | 2 | 219, 223 (both in internal `.create_run_dir()`) |
| summary_pop.R | 1 | 38 |
| **total** | **62** | |

### 2. Category 1, 2, 3 checks and every `warning()`

Exported entry point is the user-facing function that reaches the check.
"allow" / default / downstream columns are filled for category 3 only.

| file:line | function | what it refuses / warns | cat | proposed `allow` | default | downstream behavior if bypassed |
|---|---|---|---|---|---|---|
| sql_utils.R:50 | `mutate_table()`, `mutate_group_*()`, every `add_*(...)` extra-column path (via `prepare_extra_cols()`) | writing a reserved (package-managed) column | 1 | — | refuse | — |
| sql_utils.R:55 | same + `define_table()` | identifier not `^[A-Za-z][A-Za-z0-9_]*$` (pasted raw into SQL) | 1 (SQL injection) | — | refuse | — |
| sql_utils.R:61 | same | identifier is a SQL keyword | 1 | — | refuse | — |
| sql_utils.R:321 | `mutate_table()`, `define_table()`, `add_*(...)` extras | **warning**: logical `NA` cannot be typed, column created as VARCHAR | 3 | `untyped_na_column` | refuse (recommend) | column is VARCHAR forever; later numeric/logical writes are stored as text with no further warning — silently wrong, hence refuse |
| sql_utils.R:339 | same | **warning**: all-NA non-logical vector of unknown class -> VARCHAR | 3 | `untyped_na_column` | refuse | as above |
| sql_utils.R:348 | same | value class has no DuckDB mapping (factor, list, ...) | 2 | — | refuse | — |
| remove_rows.R:204 | `remove_rows()` | unfiltered delete of a whole table unless `confirm_all = TRUE` | 3 (destructive-intent guard) | `delete_all_rows` (replaces `confirm_all`) | refuse | empty table; every reader already handles an empty table. But see note N3: emptying `ind_meta`/`trait_meta`/`genome_meta` in single-table mode orphans dependants |
| remove_rows.R:220 | `remove_rows()` | **warning**: filter matched 0 rows, nothing deleted | 3 | `empty_selection` | warn | no-op is literally what was asked |
| remove_rows.R:229 | `remove_rows()` | delete from `TABLE_NO_ROW_DELETE` tables (3 genome_effect* tables, views, `_schema_meta`) | 1 | — | refuse | — |
| remove_rows.R:236 | `remove_rows()` | table not in `TABLE_ROW_KEYS` (includes every `define_table()` user table, even with a primary key) | 2 | — | refuse | — |
| remove_rows.R:286 | `remove_rows(tables=)` | target not in the static `IND_TABLE_ID_IND_COLS` list (message wrongly says "do not have an id_ind column" for user tables that do) | 2 | — | refuse | — |
| archive_replicate.R:159 | `archive_replicate()` | a store_and_reset table already has a `replicate` column | 1 (stamp column collision) | — | refuse | — |
| archive_replicate.R:280 | `archive_replicate()` | in-memory pop with only `db_name_archive` set: no directory to derive | 2 | — | refuse | — |
| archive_replicate.R:360 | `archive_replicate()` | working-table columns differ from the existing archive table | 2 (no schema-union path; could become a NULL-fill feature) | — | refuse | — |
| archive_replicate.R:368 | `archive_replicate()` | column type differs from archive | 1 (would cast / corrupt archived values) | — | refuse | — |
| mutate_table.R:504 | `mutate_table()` | **warning**: no fields given | 0 (empty call) | — | — | recommend plain `stop()` or silent no-op |
| mutate_table.R:569 | `mutate_table()` | `filter()` on a table with no registered PK (static `TABLE_PRIMARY_KEYS`; `define_table(primary_key=)` tables are not registered) | 2 | — | refuse | — |
| mutate_table.R:578 | `mutate_table()` | **warning**: filter matched 0 rows | 3 | `empty_selection` | warn | no-op |
| mutate_table.R:596 | `mutate_table()` | tibble value on a table with no registered PK | 2 | — | refuse | — |
| mutate_table.R:622 | `mutate_table()` | tibble contains ids not in the table | 3 | `unknown_ids` | refuse (typo'd ids are silently dropped) | unmatched rows are skipped by the UPDATE JOIN; nothing else changes |
| mutate_table.R:646 | `mutate_table()` | plain vector of length > 1 (positional assignment) | 1 (Design principle 7, no implicit ordering) | — | refuse | — |
| mutate_table.R:660 | `mutate_table()` | `.set_default = TRUE` with a tibble value | 2 (no single DEFAULT exists) | — | refuse | — |
| mutate_table.R:735 | `mutate_table()` | **warning**: replaced N existing non-NULL values | 3 | `replace_values` | warn (literally asked) | values overwritten; same outcome as `mutate_group_*(overwrite = TRUE)` |
| mutate_group_seq.R:112 | `mutate_group_seq()` | `n_groups` > filtered rows | 3 | `empty_groups` | refuse (asked for N groups, gets fewer; without `include_leftover` every row is NULL) | groups of size 0 never appear; size-1 groups when `include_leftover = TRUE`. Readers of the column see only the labels that exist |
| mutate_group_seq.R:227 | `mutate_group_seq()` | **warning**: rows left unassigned (NULL) because `include_leftover = FALSE` | 3 | `unassigned_rows` | warn | NULL group column; `add_phenotype()` effects already route NULL classes via `null_class_action`, components via `missing_component_action` |
| mutate_group_named.R:145 | `mutate_group_named()` | fewer rows than requested `counts` (`underfill_action = "error"`) | 3 (already has an action enum) | keep `underfill_action` (it selects a behavior, not a bypass) | refuse | "proportional" / "sequential" already defined |
| mutate_group_named.R:168 | `mutate_group_named()` | **warning**: extra rows left NULL | 3 | `unassigned_rows` | warn | as seq.R:227 |
| mutate_group_helpers.R:37 | all `mutate_group_*()` | table has no registered PK | 2 | — | refuse | — |
| mutate_group_helpers.R:49 | all `mutate_group_*()` | **warning**: filter matched 0 rows | 3 | `empty_selection` | warn | no-op |
| mutate_group_helpers.R:77 | all `mutate_group_*()` | existing column has an incompatible type | 1 (type safety, principle 5) | — | refuse | — |
| mutate_group_helpers.R:88 | all `mutate_group_*()` | filtered rows already non-NULL and `overwrite = FALSE` | 3 | `replace_values` (replaces `overwrite`) | refuse | values overwritten |
| tidybreed_pop.R:68 | every exported fn via `validate_tidybreed_pop()` | **warning**: `pop$tables` lists tables missing from the DB | 1 (registry out of sync; later calls fail with raw SQL errors) | — | recommend refuse | — |
| tidybreed_pop.R:296 | `close_pop(results_dir=)` | **warning**: in-memory DB cannot be moved; issued *after* disconnect | 0 | — | — | recommend `stop()` before disconnecting |
| tidybreed_pop.R:301 | `close_pop(results_dir=)` | **warning**: DB file vanished, nothing moved | 0 (environment) | — | — | keep warning |
| tidybreed_pop.R:312 | `close_pop(results_dir=)` | destination file already exists | 3 (destructive) | `replace_file` | refuse | destination replaced. Note the connection is already closed when this fires |
| define_table.R:86 | `define_table()` | redefining a system table | 1 | — | refuse | — |
| define_table.R:138 | `define_table()` | table exists and `overwrite = FALSE` | 3 | `replace_table` (or keep `overwrite`, see §3) | refuse | table dropped and recreated; user tables have no dependants in package code |
| restore_pop.R:55 | `restore_pop()` | `":memory:"` path | 2 (nothing on disk) | — | refuse | — |
| restore_pop.R:95 | `restore_pop()` | **warning**: `genome_meta` empty | — (informational) | — | — | recommend `message()` |
| restore_pop.R:100 | `restore_pop()` | **warning**: no `genome_meta` (normal before `define_genome()`) | — (informational) | — | — | recommend `message()` |
| restore_pop.R:115 (`stop_stale`, called at 127, 139, 170, 182, 187, 197, 203, 221, 230) | `restore_pop()` | file has a pre-v0.65 / v0.68 / v0.72 / v0.73 / v0.74 / v0.74.3 schema shape or stored names | 2 (no migration) | — | refuse | — |
| restore_pop.R:264 | `restore_pop()` | **warning**: no `_schema_meta` (pre-v0.36 file) | — | — | — | likely unreachable: such a file also lacks `trait_var_comp.line_name` and stops at 182 first. Candidate for deletion (legacy path) |
| schema.R:1045 | `schema(order = "size")` | in-memory pop has no blocks to size | 2 | — | refuse | — |
| schema.R:1049 | `schema(sizes = TRUE)` | **warning**: sizes ignored for in-memory pop | 2 (degraded) | — | warn | keep, or `message()` |
| summary_pop.R:48 | `summary()` | **warning**: unknown table names, silently dropped | 0 | — | — | recommend `stop()` — `get_table()` and `describe_table()` refuse the same input |

Totals: cat 1 = 10, cat 2 = 13 rows (the 9 `stop_stale` sites counted as one row), cat 3 = 15 (8 refusals + 7 warnings), plus 7 informational/cat-0 warnings.

### 3. Existing opt-out-like arguments (whole package)

| file:line | argument | what it bypasses / does | recommendation |
|---|---|---|---|
| remove_rows.R:186 | `confirm_all = FALSE` | refusal of an unfiltered whole-table delete (L204) | fold into `allow = "delete_all_rows"` |
| remove_rows.R:187 | `dry_run = FALSE` | not a bypass; preview mode | keep |
| define_table.R:77 | `overwrite = FALSE` | refusal when the table exists (L138); `TRUE` drops it | fold into `allow = "replace_table"` or keep as a verb-like mode — decide once for all `overwrite` |
| mutate_group_seq.R:87, mutate_group_named.R:105, mutate_group_concatenate.R:77 (checked in mutate_group_helpers.R:88) | `overwrite = FALSE` | refusal to replace non-NULL values | fold into `allow = "replace_values"`; **inconsistent** with `mutate_table()`, which replaces and only warns (L735) |
| mutate_group_named.R:103 | `underfill_action = "error"` | underfill refusal (L145); alternatives "proportional"/"sequential" | keep (behavior selector, `_action` family) |
| mutate_group_seq.R:85 / named.R:102 | `include_leftover = FALSE` | not a bypass, chooses where leftovers go; drives warnings L227 / L168 | keep |
| mutate_group_concatenate.R:76 | `null_action = "propagate"` | NULL handling selector | keep, but rename to `null_class_action`-style? it means something different; leave |
| mutate_table.R:479 | `.set_default = FALSE` | not a bypass | keep |
| open_pop.R:90 | `clean = TRUE` | **silently deletes an existing DB file and tool run dirs** (no check at all; default is destructive) | rename/invert: default should refuse an existing file; `allow = "replace_file"` to overwrite (same name as close_pop:312) |
| define_trait.R:64 | `overwrite = FALSE` | refusal on an existing trait (L81); TRUE replaces and clears `phenotype_effects` | fold into `allow = "replace_trait"`, or keep `overwrite` uniformly |
| define_phenotype.R:307 | `overwrite = FALSE` | refusal on existing phenotype (L456) | as define_trait |
| define_effect_fixed_cov.R:55, define_effect_fixed_class.R:50, define_effect_random.R:87 (via effect_helpers.R:28 `.handle_effect_overwrite`) | `overwrite = FALSE` | refusal on an existing effect name | as above |
| define_index.R:73 | `overwrite = FALSE` | **not a refusal**: `FALSE` is a silent no-op (ON CONFLICT ... DO NOTHING, L138-153) | inconsistent with every other `overwrite`; make `FALSE` refuse |
| define_chromosome.R:142 | `overwrite = TRUE` | **default TRUE**; `FALSE` refuses an existing rule (L283) | inconsistent default; align with the rest |
| define_phenotype.R:305 | `missing_component_action = c("skip","error")` | skip+warn vs refuse for missing components (stored in `phenotype_meta`) | keep (stored behavior, per memory note) |
| define_phenotype.R:306, define_condition_change_action.R:45 | `condition_change_action = c("error","independent")` | refusal when residual strata differ | keep (stored behavior) |
| define_effect_fixed_class.R:49 | `null_class_action = c("skip","error","zero")` | NULL class handling | keep |
| define_genome_effect_terms.R:202 | `require_complete = FALSE` | opt-**in** strictness (indicator tables must cover all states) | keep; or invert into the policy (refuse by default, `allow = "incomplete_indicator"`) |
| define_genome_effect_terms.R:196 | `mode = "append"` / `replace_*` | the only sanctioned replacement path for genome-effect rows | keep (cat 1 boundary) |
| define_genome_effect_terms.R:219 | `allow_reserved_owner = FALSE` (internal `.ge_write_terms()` only) | writing under the reserved owner `"generated"` | keep internal; **name clash risk** with the new `allow =` — rename to e.g. `.generator_owner = TRUE` |
| define_additive_effects.R:315, define_genome_effects.R:310, qtl_congruence.R:530 | `warn_bounds = c(0.8, 1.25)` (`NULL` turns off) | silences the calibration-diagnostic warning | keep as a tuning parameter, but `NULL`-to-silence overlaps `allow`; consider `allow = "calibration_drift"` for the off switch |
| correlated_draws.R:240 | `tolerance = NULL` (internal) | numeric PSD/symmetry tolerance | keep internal |
| define_phenotype.R:693 | `strict` in `.split_component_names()` (internal) | only changes the message suffix | not an opt-out |
| add_phenotype.R:247, define_additive_effects.R:316, define_genome_effects.R:311, mutate_group_seq.R:86, mutate_group_named.R:104 | `seed = NULL` | calls `set.seed()` — writes `.Random.seed` from R/ | not an opt-out; flag vs. CLAUDE.md "nothing in R/ touches .Random.seed" |

### 4. Notes

### N1. Shared validator — yes, put `.check_allow()` in R/arg_checks.R
`arg_checks.R` holds two generic, semantics-free helpers: `.check_scalar(x, name, min, max,
exclusive, whole, finite)` (14 call sites) and `.check_flag(x, name)` (4). Other validators are
scattered and domain-specific (`validate_sql_identifier()` in sql_utils.R,
`validate_tidybreed_pop()` in tidybreed_pop.R, `.check_effect_name_input()`,
`.check_phenotype_exists()`, `.pvc_check_*`, `.qtl_validate_*`, ...).
`.check_allow(allow, known)` fits the file's charter exactly (generic, uniform wording) and should
live there, with a companion that builds the message suffix, e.g.
`.allow_hint(name)` -> `"  To proceed anyway: allow = \"name\""`, and a predicate
`.allowed(allow, name)`. `.check_allow()` should refuse unknown names (cat 0) and list the
function's `known` set, so a typo cannot silently fail to relax a check. Keep a single package-wide
registry of allow names (a named character vector with one-line meanings, like
`TABLE_NO_ROW_DELETE`) so the same name means the same thing in every function and docs can be
generated from it.

### N2. Miscategorised / inconsistent today
- `mutate_table()` replaces values with a warning (L735); `mutate_group_*()` refuses the same
  thing unless `overwrite = TRUE` (helpers L88). One name (`replace_values`), one default.
- `define_index(overwrite = FALSE)` is a silent no-op; every other `overwrite = FALSE` refuses.
  `define_chromosome()` defaults `overwrite = TRUE`.
- `open_pop(clean = TRUE)` deletes an existing database silently by default — the most
  destructive default in the package and has no check.
- `summary_pop.R:48` warns on unknown table names while `get_table()` / `describe_table()` stop.
- `close_pop()` L296 warns about an impossible argument after the connection is already closed.
- `restore_pop.R:95/100` are normal states, not warnings; `restore_pop.R:264` is likely dead
  legacy code (CLAUDE.md: no legacy paths).
- `tidybreed_pop.R:68` (registry out of sync) is an integrity signal issued as a warning.
- `mutate_group_named.R:193` refuses proportions not summing to 1 within ~1.5e-6 and then
  normalises (L198) — the normalisation is effectively dead; either drop the check (cat 0 ->
  accept relative weights) or drop the normalisation.
- `remove_rows.R:286` message is wrong for user tables that do have `id_ind`.
- `remove_rows.R:245` fires when the user `select()`ed away key columns; it should re-read the
  physical table the way `resolve_subset_ids()` does (sql_utils.R:497-519) rather than refuse.

### N3. Missing category-1 checks (no refusal today)
- `remove_rows()` single-table mode deletes from `ind_meta`, `genome_meta`, `trait_meta`,
  `phenotype_meta`, `index_meta` with no dependant check: orphans in `ind_haplotype`/`ind_tgv`,
  `genome_effect_members` (locus_id), `genome_effects` (trait_name), `phenotype_components`, etc.
  Only the genome_effect* tables are protected (sql_utils.R:226).
- `archive_replicate()`: no guard against archiving the same `replicate` number twice (rows
  duplicate silently); archive writes (L188-216) are not in a transaction, so a failure at the
  2nd table (e.g. L360) leaves a partial archive, and a retry duplicates the first table's rows;
  `store_once` tables changed between replicates are skipped with only a `message()`;
  `ATTACH '<path>'` (L181) does not escape quotes in the path.
- `.create_run_dir()` (open_pop.R:233, reached from `add_ebv()` BLUPF90 path) calls `sample()`
  for a directory suffix, consuming the R RNG stream — a run-dir name changes every later
  seeded draw. Determinism issue (cat 1).
- `slice_sample.tidybreed_table()` (mutate_table.R:251) uses DuckDB `random()`, not the R seed,
  and `slice_*` is not recorded in `pending_filter`, so `resolve_subset_ids()` /
  `mutate_table()` filtered paths ignore it. Either refuse slices on tidybreed_table or record them.
- `remove_rows()`/`mutate_table()` write temp tables with `dbWriteTable()`; CLAUDE.md notes that
  advances the RNG. Harmless today only because they sit outside the `add_phenotype()` path.
- `mutate_table()` does not check that a scalar's inferred type matches an existing column's
  type (helpers L77 does this for `mutate_group_*()`); DuckDB casts or errors raw.
