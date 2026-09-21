# Correlated effects — Phase 5 results

**Spec:** `plans/sample_correlated_effects.md` (v3.5), §5.3, §5.4 (adapter
responsibilities), §5.5 Stage 2, §5.7, D2, D6, §7 "Residual blocks",
"Defect 4", "Fixed coordinates", "Repeated records", and §8 Phase 5.
**Branch:** `feat/genome-effects-v49`.
**Follows:** `plans/sample_correlated_effects_phase_4.md`.
**Status:** complete. Full suite green.
**Date:** 2026-09-21.

Phase 5 is the first phase that changes what `add_phenotype()` *computes*.
Residuals are now drawn through the Phase 3 resolver, one covariance block
at a time, each record conditional on the residuals its individual has
already realized for the block's other phenotypes — in this call, in an
earlier call, or supplied by the caller. Defects 1 (residuals never
persisted), 2 (joint draw gated on same-call identical sets) and 4
(heterogeneous variance ignored on single-phenotype calls) are closed.
The stage boundaries from Phase 4 did not move.

## What shipped

| Piece | Where | Notes |
|---|---|---|
| `.ap_resolve_residuals()` | `R/add_phenotype_stages.R` | Entry point of the residual adapter. Validates `user_residual` against the plan, loads the blocks touching the model-path phenotypes (`find_covariance_blocks()`), errors "No residual variance found" for a phenotype in no block, and runs `.ap_residual_block()` per block in loader order. Returns per phenotype `value`, `level` (the selected stratum per record) and `var_unconditional`. |
| `.ap_residual_block()` | same | One block: planned `(entity, phenotype)` rows → sorted entities `(id_ind, pheno_number)` → stratum per entity (D2 selection, fallback, error) → stored coordinates of every block member at the same `pheno_number` via a registered view → D6 (every call) then D2 on the stored set → `observed` = stored + fixed → one `resolve_correlated_draws()` call per `(stratum, sample set)` group in sorted group order → values back in planned-record order. |
| `.ap_fixed_residuals()` | same | The `user_residual` contract (§5.3): plain vector iff exactly one model-path phenotype; otherwise a named list over any subset of the model-path phenotypes; length = planned records; finite; no per-`id_ind` names. |
| `.ap_resolve()` | same | Reordered: every named-effect draw (pre-draw, then marginal new levels in plan order) **before** the residual adapter; record assembly afterwards in plan order with no RNG. `records` gain `residual_value` / `residual_condition_level` on the model path. |
| `.ap_liability_records()` | same | Takes the adapter's element instead of `resid_info`; writes the two residual columns; the `prevalence` threshold errors when there is no unconditional variance. |
| `add_phenotype()` | `R/add_phenotype.R` | `user_values` + `user_residual` together is an error. Roxygen: the residual model paragraph and the `user_residual` argument rewritten for the new contract. |
| Deleted | `R/add_phenotype_stages.R`, `R/phenotype_helpers.R` | `.ap_joint_residuals()` (ex-§8.5 and its equal-planned-sets restriction), the independent `rnorm` residual branch, `sample_residuals()` (the last `MASS::mvrnorm()` on the residual path), `get_residual_cov()` (and with it the `condition_column == ""` tolerance noted in Phase 4). Three `man/` pages go with them. |
| Tests | `tests/testthat/test-add_phenotype_residuals.R` (new) | 87 expectations; see below. Three existing tests that pinned the pre-Phase-5 `NULL` residuals updated (`test-phenotype_schema.R`, two in `test-phenotype_cov_block.R`), the Defect 4 assertions added to the composite test, one warning text in `test-add_phenotype_stages.R`. |
| Docs | `man/dot-ap_resolve_residuals.Rd`, `dot-ap_residual_block.Rd`, `dot-ap_fixed_residuals.Rd` new; `add_phenotype.Rd`, `add_phenotype_stages.Rd`, `dot-ap_resolve.Rd`, `dot-ap_liability_records.Rd`, `define_residual_cov.Rd` ("How the block is sampled"), `correlated_draws.Rd` updated; `get_residual_cov.Rd`, `sample_residuals.Rd`, `dot-ap_joint_residuals.Rd` removed | `NAMESPACE` unchanged. CLAUDE.md: the `add_phenotype()` Stage 2 description and the two `ind_phenotype` column notes. The Phase 8 roxygen pass still owns the rest of the user-facing narrative. |

## What a user sees differently

- **Residuals are correlated across calls.** `add_phenotype("A")` then,
  a hundred simulated days and a culling later, `add_phenotype("B")`: the
  survivors' `B` residual is `E[B | A] + conditional noise`, the culled get
  no `B` row, nothing else is touched. Same for `B` first, for a mixed
  subset where only some animals have an `A`, for three-phenotype blocks
  with different observation patterns per animal, and for two phenotypes
  with different planned sets in one call (which used to fall back to
  independent draws).
- **Heterogeneous residual variance works on every call.** A phenotype
  with `condition_column = "sex"` strata drew everyone from the
  unconditional variance on a single-phenotype call (Defect 4). Each record
  now draws from its stratum; a phenotype with *only* conditional strata
  can be sampled at all (it used to error "No residual variance found").
- **`residual_value` and `residual_condition_level` are populated** for
  every model-generated record, liability scale, with the *selected*
  stratum (`NULL` for the unconditional `R`). `user_values` and
  `derived_formula` records leave both `NULL`.
- **The realization lock is live.** `define_residual_cov()` /
  `define_phenotype(residual_var = )` on a block that has any model-path
  record now error with the `remove_rows()` recipe — the D3 rule from
  Phase 2, which until now could only be triggered by an `UPDATE`.
- **Stratum fallback is explicit.** A `NULL` condition value falls back to
  the unconditional `R` silently; a non-`NULL` value matching no stratum
  falls back with a warning (count, levels, ≤ 5 ids); with no unconditional
  stratum either case is an error naming the levels and the count. The
  pre-Phase-5 "residual set to 0" fallback is gone.
- **D2 / D6 at sampling.** A stored residual under a different stratum
  than the current record resolves to: error by default with up to five
  `id (coordinate: stored under X, now Y)` examples, or a warning naming the
  dropped coordinates under `condition_change_action = "independent"`;
  disagreeing actions across the block error before either.
- **`user_residual` is checked and richer.** A plain vector is accepted
  only when exactly one phenotype in the call is model-generated (it used to
  be applied to every phenotype); a list may name a subset and the rest are
  drawn conditional on it; naming a derived phenotype, per-`id_ind` names,
  non-finite values and combining with `user_values` are errors; a value off
  the support of a singular covariance is an error from the resolver.
  Supplied values are stored and condition later calls.
- **Categorical `prevalence` with only conditional strata errors** (the
  cut-point needs the marginal variance; use `thresholds =` or add an
  unconditional stratum). This path could not run before Phase 5.
- **RNG order.** Named-effect draws now all precede the residual draws; on
  a call with new random-effect levels the seeded values differ from Phase
  4 (allowed and expected).

## Why the shape is what it is

**One adapter, one block loop, one resolver call per group.** The adapter
follows §5.4's nine responsibilities in order and nothing else. Entities are
`(id_ind, pheno_number)`; the coordinates are the block's phenotypes,
including members not in the call (they can be stored, so they condition;
they are never sampled). Grouping is by `(stratum, sample set)` only — the
observed-pattern split is the resolver's (v3.3) — so a mixed subset is one
call, and the stream is entity order regardless of who has what.

**Sorted everything.** Entities by `(id_ind, pheno_number)` in byte order;
groups by a key that puts the unconditional stratum first, then levels in
byte order, then the sample set; blocks in loader order (first member).
Every exact test in the new file replays this order from `rnorm()`.

**Named effects before residuals.** The pre-Phase-4 code interleaved the
marginal named-effect draws with the residual draw per phenotype; the plan's
Stage 2 puts every named-effect block before the residual blocks. Doing the
reorder now, rather than in Phase 6, means Phase 6 replaces bodies without
touching the order again — one seeded-output change instead of two.

**Fixed coordinates condition in the same call.** `user_residual` values go
into `observed` alongside stored values, so a supplied `A` conditions a
generated `B` in one call (§5.3's "fixed" state), and both are written with
the stratum the record resolved to.

**Delete, don't rewrite, `get_residual_cov()`.** The plan said "rewrite
around strata"; `find_covariance_blocks()` already is that function with
the D1 invariants re-checked. The old reader's `condition_column == ""`
tolerance and first-column-wins stratum selection go with it. The only
other consumer was the prevalence threshold, which now reads the block's
unconditional diagonal.

**A warning for an unmatched non-`NULL` level.** D2 says both `NULL` and
"matches nothing" fall back silently. `NULL` is the documented "no group"
state and should be quiet; a farm code with no stratum is more often a
typo, and the pre-Phase-5 code warned on it. Keeping a warning for the
non-`NULL` case only was the smallest honest reading, and the plan header
records it as a refinement.

**Prevalence needs a marginal variance.** With strata there is no single
`V_E`; the old code fell to `0` (via `NA → 0`), which puts the threshold in
the wrong place silently. An error with the two fixes is better than a
guess, and this path was unreachable before Phase 5 anyway.

## Tests (`test-add_phenotype_residuals.R`, 87 expectations)

Exact tests replay the resolver's contract (`n × m` normals in entity
order, applied through `chol(C)`), so for `R = [1 .8; .8 1]` a conditioned
`B` is `.8 A + .6 z`. Distributional checks use `residual_value` directly.

**Residual blocks (§7 items 1–9).** Sequential `A → B` on 400 animals
(exact, `cor ≈ .8`, `var ≈ 1`); `A → B` on a culled subset (exact on the
survivors, the culled have no `B`, the `A` rows are byte-identical); `A` and
`B` in one call with `B` sex-limited (the `A`-only group draws first, then
the joint group; `n_A + n_B` normals consumed); `B` for a mixed group
(conditioned where `A` exists, marginal elsewhere, exact through one
resolver call); `B` first then `A`; an `{A, B, C}` block where males have
`(A, B)` stored and females `A` only, `C` exact against the hand-computed
conditional moments of each pattern; two disconnected blocks in a call in
non-block order (block `{A, B}` draws first, then `{C}`; `3n` normals);
an explicit zero covariance stays one block and conditioning on `A` leaves
`B`'s stream unchanged; a phenotype in no block still errors.

**Defect 4 (items 2–4).** Sex strata with no unconditional stratum: the
call succeeds, `residual_condition_level` equals the sex, `F` entities draw
before `M` with the right scale. A `farm` column with an unmatched value
and a `NULL`: without an unconditional stratum the error names
`farm = 'F9', NULL`, the count and the ids; with one, the unmatched value
warns, both fall back with `NULL` stored, and the fallback group draws
first. Item 1 is the strengthened composite test.

**Fixed coordinates (items 1–5).** Supplied `A` + generated `B` in one call
(exact, both stored, only `n` normals consumed); both supplied (RNG-neutral,
both stored, a later `C` exact against the two-coordinate conditional);
validation (plain vector with two model phenotypes, unknown name, derived
name, wrong length with the planned count, non-finite, unnamed list,
per-`id_ind` names, `user_values` together — all errors, nothing written;
a plain vector with one model phenotype and a derived one is fine); a
non-zero value for a zero-variance coordinate errors in the resolver before
any draw, `0` is accepted.

**Repeated records (items 1–3).** `A` twice then `B` once: `A(2)` is a
marginal draw, `B(1)` conditions on `A(1)`; unequal counts (`A` for all
then for males; `B` for all then for males): `B(1)` pairs with `A(1)` for
everyone and `B(2)` with `A(2)` for the males.

**D2 / D6 at sampling.** A `farm` block with `F1`/`F2` strata: after `A`
under `F1`, moving two animals to `F2` makes `B` error naming
`A_3 (A: stored under 'F1', now 'F2')`; with `'independent'` on both
phenotypes the same call warns, names the dropped coordinate `{A}`, and
draws the moved animals marginally under `4R` after the conditioned `F1`
group (exact). A hand-flipped `condition_change_action` on one member
errors with the D6 message. Sex strata across two calls never trigger D2.

**Lock.** After `add_phenotype("A")` on an `{A, B}` block,
`define_residual_cov()` errors with "12 realized draws"; the `remove_rows()`
recipe unlocks it.

## Verification

- `test-add_phenotype_residuals.R`: 87/87.
- `test-add_phenotype_stages.R` 57, `test-add_phenotype.R` 45,
  `test-phenotype_composite.R` 39 (4 added), `test-phenotype_schema.R` 21,
  `test-phenotype_cov_block.R` 109 (D3 residual test rewritten around the
  live lock), `test-correlated_draws.R` 115, `test-formula_phenotype.R` 41,
  `test-define_phenotype.R` 38, `test-mutate_derived.R` 37 — all green.
- Full suite: 3304 passed, 0 failed, 1 skipped (a mid-work snapshot ran
  3303/0/1; the one added expectation is the named-vector guard).
- `roxygen2::roxygenise(roclets = "rd")` clean; `NAMESPACE` unchanged.

## Review pass

- The first draft built the per-entity sample set with an R-level loop
  over planned rows (`sample_of[[r]] <- c(...)`); replaced by a logical
  `entity × in-call phenotype` matrix and a vectorized key, so the adapter
  does no per-row R work beyond the resolver.
- Entity keys handed to the resolver are the `(id_ind, pheno_number)` data
  frame, not the internal `"id\rpn"` string, so a resolver error names
  `id_ind = A_3, pheno_number = 1`.
- A named plain `user_residual` vector was silently treated positionally
  (the roxygen said names were unsupported); it now errors.
- Checked against the plan's adapter list (§5.4, nine responsibilities)
  and the §5.5 Stage 2 pseudocode line by line: the only deliberate
  deviations are the two recorded in the v3.5 header (warning on unmatched
  non-`NULL` levels; the prevalence-threshold error). The D6 check runs
  once per block on every call, as the pseudocode has it (a first draft
  skipped it when nothing was stored).
- Found and recorded, not fixed: the definition-time D6 check freezes
  `condition_change_action` once every member of a block is defined
  (plan header and §8 Phase 8).

## Plan bookkeeping

- Header: v3.5, Phases 0–5 shipped; "What changed from v3.4 to v3.5"
  (adapter, RNG order, fallback reporting, `user_residual` contract,
  prevalence error, live lock, the D6 open item).
- §5.3: `user_residual` as shipped; the case list is now the test list.
- §5.5: heading and v3.4 note updated; equal-set restriction marked
  deleted; stratum-lookup paragraph points at the adapter.
- D2, D6: v3.5 notes on where the rules live.
- §7: "Residual blocks 1–9", "Defect 4", "Fixed coordinates", "Repeated
  records" marked covered.
- §8: Phase 5 ✅; Phase 6 note (the residual adapter is the template, and
  the draw order is already right); Phase 8 gains the D6 mutability item.
- `NEWS.md` under `0.71.0 (in development)` → Changed; CLAUDE.md updated.
- No version bump (deferred to Phase 8).

## Next

Phase 6 — the named-effect adapter: replace `.ap_predraw_named_effects()`
and the normal branch of `.ap_resolve_random_terms()` with the same block
loop over `phenotype_random_effects` (entity `(effect_name, level)`, no
strata, no fixed coordinates), keep the gamma/uniform marginal sampler for
1 × 1 blocks, and add the `validate_named_effect_block()` backstop in
Stage 2.
