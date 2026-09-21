# Correlated effects — Phase 2 results

**Spec:** `plans/sample_correlated_effects.md` (v3.2), §5.2, §5.6, §5.9, D1, D3,
D5, D6 and §8 Phase 2.
**Branch:** `feat/genome-effects-v49`.
**Follows:** `plans/sample_correlated_effects_phase_1.md`.
**Status:** complete. Full suite green.
**Date:** 2026-09-20.

Phase 2 centralizes the *definition* side of the feature. Every write to
`phenotype_var_comp` now goes through one validator and one writer, so by the
time the sampler (Phases 3–6) reads a covariance block it can assume the block
is complete, symmetric, PSD, consistent across strata, and unchanged since its
draws were taken. **This is the first phase with user-visible behaviour
change**: declarations that used to be silently accepted — fragments of a block,
a diagonal rewrite inside a joint block, a non-PSD matrix, a redefinition after
draws exist — are errors now, each with the call that fixes it in the message.
The sampler itself is untouched; `add_phenotype()` behaves exactly as before.

---

## What shipped

| File | Change |
|---|---|
| `R/phenotype_cov_block.R` | **new** (~500 lines). `validate_phenotype_cov_block()` — D1, strata, D6, §5.6, D3 in that order against the current table state. `write_phenotype_cov_block()` (own transaction) and `.pvc_write_block()` (caller's transaction): validate → `DELETE` the stratum's rows for the block → `duckdb_register()` + `INSERT … SELECT`. Helpers `.pvc_block_members()` (fixpoint closure over pair-row existence, all strata), `.check_condition_change_agreement()` (D6, with an optional pending `phenotype_meta` row), `validate_named_effect_block()` (§5.6, with an optional pending `phenotype_effects` row) |
| `R/define_residual_cov.R` | body replaced by argument checks + one `write_phenotype_cov_block("residual", …)` call. `condition_column` and `condition_level` must come together; `condition_table` is `NULL` on unconditional rows. No `dbWriteTable()` — RNG-neutral. Roxygen documents D1, strata, D3, D6 |
| `R/define_effect_cov_matrix.R` | one argument-validation path for all routes (`trait_names` now also applies to `"residual"`); genetic route unchanged; `"residual"` → `define_residual_cov()`; any other name → `write_phenotype_cov_block()`. `write_phenotype_var_diag()` **deleted**. `get_phenotype_var()` drops the legacy `condition_column = ''` predicate |
| `R/define_effect_random.R` | one transaction around: overwrite-delete of the old row and its draws → `.pvc_write_block()` for a supplied `variance` (or a stored-variance lookup for `NULL`) → `validate_named_effect_block()` with the pending row → register + `INSERT` of the `phenotype_effects` row. Rejects `effect_name = "residual"`. A supplied `variance` is always written (it used to be ignored when a stored value existed) |
| `R/define_phenotype.R` | `residual_var` is validated through `validate_phenotype_cov_block()` and D6 through `.check_condition_change_agreement()` (with the pending row) **before** the `phenotype_meta` write; the overwrite path **no longer deletes** the phenotype's unconditional residual rows; the residual write still delegates to `define_residual_cov()` |
| `tests/testthat/test-phenotype_cov_block.R` | **new** — 110 expectations, see below |
| `man/` | `define_residual_cov.Rd`, `define_effect_cov_matrix.Rd`, `define_effect_random.Rd`, `define_phenotype.Rd` regenerated; `write_phenotype_var_diag.Rd` removed; internal pages for the new helpers |
| `CLAUDE.md` | `phenotype_var_comp` section gains the block rules; `define_phenotype()` / `define_residual_cov()` / `define_effect_random()` sections updated |
| `NEWS.md` | Phase 2 entry under **0.71.0 (in development)** — Breaking |
| `R/set_residual_cov.R`, `R/set_random_effect_cov.R` | **deleted** — two-line tombstones for functions removed in v0.7.0 (found in the review pass) |

## What a user sees differently

| Call | Before | Now |
|---|---|---|
| `define_residual_cov(c("A","B"), R)` then `define_residual_cov(c("B","C"), R2)` | accepted; `Cov(A, C)` undefined; the sampler's `load_*` returned `NULL` for `{A, B, C}` and fell back to independence | **error**: names the block `{A, B, C}`, the omitted `A (paired with B)`, and the complete-block call to make |
| `define_residual_cov("A", 1×1)` or `define_phenotype("A", residual_var = )` when `A ∈ {A, B}` | rewrote `Var(A)` in place, silently changing the correlation and possibly breaking PSD | **error** naming the block and `define_residual_cov(pop, c("A", "B"), R)` |
| `define_effect_random("B", "pen", variance = 9)` when `B ∈ {A, B}` for `pen` | `variance` silently ignored | **error** naming `define_effect_cov_matrix(pop, "pen", R)` |
| `define_effect_cov_matrix("pen", R)` after `add_phenotype()` has drawn pen levels | accepted; later draws conditioned under a different matrix than earlier ones | **error** with the `remove_rows()` recipe on `phenotype_random_effects` |
| non-PSD residual or named-effect matrix | accepted (diagonal ≥ 0 was the only check); `MASS::mvrnorm()` failed later, in `add_phenotype()` | **error** at define time |
| `define_phenotype("A", overwrite = TRUE)` with `A ∈ {A, B}` | **deleted** `(A, A)` and `(A, B)` from the unconditional block, leaving `{A, B}` half-defined | leaves `phenotype_var_comp` untouched |
| `define_phenotype("B", condition_change_action = "independent")` with `A ∈ {A, B}` set to `"error"` | accepted | **error** naming both values (D6) |
| `define_effect_random("B", "pen", source_column = "herd")` when `A`'s `pen` reads `sex`, `{A, B}` a `pen` block | accepted; a pen id and a herd id were about to be treated as one random vector | **error** (§5.6) |
| `define_effect_random(..., distribution = "gamma")` into a block of ≥ 2 | accepted; drawn normal by the correlated path anyway | **error**; a 1 × 1 gamma effect stays legal |
| `define_residual_cov(c("A","B"), R, condition_column = "sex")` (no level) | wrote rows with a column and a `NULL` level | **error**: both or neither |
| `define_effect_cov_matrix("residual", M, trait_names = )` | `trait_names` ignored for the residual route | honoured |

## Why the shape is what it is

**One writer, validated before the `DELETE`.** The plan's §5.9 asked for the
validator to run "inside the transaction, before `COMMIT`", by analogy with
`validate_genome_effects()`, which checks the table *after* the write. D1's
algorithm is stated on the pre-write state (find touched blocks, `N == U`, lock,
*then* delete and insert), and its messages need the omitted coordinates —
easiest to name before anything moves. So `.pvc_write_block()` validates first,
then deletes, then inserts, all in the caller's transaction. The guarantee is the
same (a rejected call changes nothing) and the transaction still matters: the
test suite proves it by mocking `next_int_id()` to fail after the `DELETE` and
asserting the prior rows are back.

**Block discovery is a fixpoint over pair-row existence across all strata.**
`.pvc_block_members()` repeats `phenotype_name_1 IN (S) OR phenotype_name_2 IN (S)`
until `S` stops growing. It deliberately ignores `cov_value` (an explicit `0` is
an edge) and `condition_*` (a block spans its strata), which is what makes the
strata rule checkable in one place: after `N == U`, every other stratum's name
set must equal `N`.

**Growing a multi-stratum block requires clearing it.** Found while writing the
strata check: with `{A, B}` unconditional and `{A, B}` for `M` on disk, no
sequence of calls reaches `{A, B, C}` in both strata, because whichever stratum
is written first is checked against the one that still says `{A, B}`. The
alternatives were to allow a transient inconsistency (a stratum that is a strict
subset of the block, which §5.2 forbids the sampler from having to reason about)
or to accept the deadlock and give the way out in the message. The message gives
the `remove_rows()` call on `phenotype_var_comp`, and the test runs that recipe
and redeclares. Blocks come from REML output and are declared once; this is not
a workflow worth a special case. Recorded as a v3.2 note on D1.

**`define_effect_random()` is one transaction and always honours `variance`.**
The old body ran the variance write *before* `.handle_effect_overwrite()`, so
with the D3 lock in place an `overwrite = TRUE` call on a realized singleton
would have failed on the draws it was about to delete. Reordering (overwrite →
variance → §5.6 → row) inside one transaction fixes that and makes a §5.6
rejection restore the row the overwrite deleted. Making a supplied `variance`
always write is the honest reading of D5 ("the D1 algorithm with
`N = {phenotype}`"): silently discarding a user's number was the old behaviour,
not a feature. `overwrite = TRUE` discarding the phenotype's stored draws for
that effect is unchanged and now documented on the argument; it is the explicit,
per-effect act D3 wants, not a covariance `force`.

**`define_phenotype()` validates before it writes `phenotype_meta`.** The
function is not itself transactional (four tables, `ensure_trait_tables()`,
`dbWriteTable()` for components) and making it so was out of scope. Running the
residual-block validator and the D6 check first — cheap, read-only — means a
rejected `residual_var` or a disagreeing `condition_change_action` leaves the
old `phenotype_meta` row in place. `define_residual_cov()` then validates again
inside its own transaction; two validations of a 1 × 1 block cost nothing.

**The old overwrite `DELETE` had to go.** `define_phenotype(overwrite = TRUE)`
deleted `phenotype_var_comp` rows with `phenotype_name_1 = <this>` and
`condition_column IS NULL` — for a phenotype in `{A, B}` that is `(A, A)` and
`(A, B)`, leaving `(B, A)` and `(B, B)`: a block with one direction of one pair.
D5's table says an overwrite without `residual_var` leaves the table untouched;
the plan called that "as today", which it was not. Corrected in the plan.

**`condition_table` is `NULL` on unconditional rows.** The two writers disagreed
(`'ind_meta'` vs `NULL`); a table name without a column is meaningless; the
reader only ever consulted it on conditional rows. Unified on `NULL`.

**Unchanged on purpose.** The genetic route of `define_effect_cov_matrix()`
(`trait_var_comp`) keeps its bare `DELETE` + `INSERT` — it is not part of this
plan and `MASS::mvrnorm()` already rejects a non-PSD `G`. `ensure_trait_tables()`
callers in the three `define_*` functions stay (noted in Phase 0 as follow-up).
`add_phenotype()` is not touched; its `all_equal` gate, `MASS::mvrnorm()` call
and the named-effect §7.5 block go in Phases 5–6.

## Tests (`test-phenotype_cov_block.R`, 110 expectations)

Plan §7 items covered: Residual blocks 10–15; Named effects 5–6 at both writer
sites and the writer half of 7; Diagonal writers 1–5; Reproducibility and
integrity 6–8.

- **D1** — fragment errors name the missing pair and the complete call, table
  unchanged; complete 3 × 3 with an explicit `0` stores 9 rows and the zero is a
  row; values land by name regardless of dimname order; non-PSD, asymmetric,
  `NA`, negative variance and wrong dimnames rejected; a zero eigenvalue (perfect
  correlation) accepted; strict subset and 1 × 1 inside a block name the block;
  redeclaring replaces exactly; two singletons merge when declared together (the
  README pattern); `trait_names` on the residual route.
- **Strata** — `{A}` conditional under `{A, B}` errors; `{A, B}` per level
  succeeds and each stratum can be redeclared alone; second condition column
  errors; column-without-level and level-without-column error; conditional-only
  blocks (no unconditional stratum) are valid; growing a multi-stratum block
  errors with the recipe, and the recipe followed by redeclaration works.
- **Rollback** — `next_int_id()` mocked to fail after the `DELETE`: prior rows
  intact, connection usable.
- **D3** — named: draws from `add_phenotype()` lock the block, message carries
  the `phenotype_random_effects` recipe and the "records not removed" sentence,
  redefinition succeeds after the recipe. Residual: `NULL` `residual_value` rows
  (everything today) do **not** lock; one row set by SQL locks with the
  `!is.na(residual_value)` recipe; the recipe removes exactly the locking rows.
- **D5** — `define_phenotype(residual_var = )`: singleton written, unrealized
  singleton overwritten, inside a block errors naming `define_residual_cov()`
  and `phenotype_meta` is not written (fresh and overwrite cases); overwrite
  without `residual_var` leaves a 2-block intact; realized singleton errors with
  the recipe. `define_effect_random(variance = )`: singleton, overwrite,
  no-variance-no-stored error, inside a block errors naming
  `define_effect_cov_matrix()` with no `phenotype_effects` row written, joins
  without `variance`; `overwrite = TRUE` discards draws; a failure mid-call
  restores the deleted row.
- **§5.6** — source column, source table and distribution mismatches refused at
  `define_effect_random()`; pre-existing disagreement refused at
  `define_effect_cov_matrix()`; a fixed effect sharing the name cannot be
  joined; a 1 × 1 gamma effect draws gamma (all draws positive), persists, and
  joining it errors naming `A uses "gamma"`.
- **D6** — block before phenotypes: second phenotype with a different value
  errors and is not written; fixing the first then defining the second works;
  flipping one back is refused; `define_residual_cov()` over disagreeing
  phenotypes errors; singletons never need agreement.
- **RNG** — `.Random.seed` identical across a `define_residual_cov()` call.

## Verification

- `test-phenotype_cov_block.R`: 110 passed.
- Affected files re-run individually: `test-define_phenotype.R` (38),
  `test-add_phenotype.R` (45), `test-phenotype_composite.R` (35),
  `test-archive_replicate.R` (35), `test-phenotype_schema.R` (21),
  `test-schema-registries.R` (110), `test-define_effect_fixed_cov.R` (26),
  `test-formula_phenotype.R` (41) — all green, no edits needed. No existing test
  relied on the behaviours that changed.
- Full suite: **3042 passed, 0 failed, 0 errors, 1 skipped** (2932 before Phase 2 + 110 new).
- `roxygen2::roxygenise(roclets = "rd")` clean; `NAMESPACE` unchanged.

## Review pass (after the suite)

A second read of the module plus adversarial probes not in the test file, all
behaving as intended: two different effects on overlapping phenotypes are
independent blocks; a three-way chain `{B, C, D}` after `{A, B}` names `A`;
four residual singletons merge into one 4-block and cannot be split back;
`define_effect_random(overwrite = TRUE)` on one member of a realized named
block discards only that member's draws and the block stays locked by the
other member's; a numeric `condition_level` round-trips and redeclares its own
stratum; errors inside `define_effect_random()` leave no open transaction;
blocks survive `close_pop()` / `restore_pop()`; a 6 × 6 block with permuted
dimnames round-trips **bit-exactly** (the old `format(x, scientific = FALSE)`
`INSERT` text rounded to seven significant digits — register + `INSERT` fixes
that quietly).

Two things found and handled outside the phase's scope:

- `R/set_residual_cov.R` and `R/set_random_effect_cov.R` were two-line
  "removed in v0.7.0" tombstones. Deleted (CLAUDE.md: leftover compatibility
  files are debt).
- **Pre-existing sampler crash**: a multi-phenotype `add_phenotype()` call
  whose individuals are all skipped reaches `MASS::mvrnorm(n = 0)` in
  `sample_residuals()` and errors with "non-conformable arguments". Untouched by
  Phase 2 (the file is not in the diff); recorded under Defect 2 in the plan,
  covered by Phase 3's empty-entity-set contract and gone with Phase 5.

## Plan bookkeeping

- Header: status → v3.2, Phases 0–2 shipped; new "What changed from v3.1 to
  v3.2" list (eight implementation notes, no decision changes).
- §0 table, §2 third-writer observation, §5.2 and §5.7 line references
  refreshed.
- §5.6: writer sites marked shipped; `add_phenotype()` backstop assigned to
  Phase 6 and named as the same function.
- §5.9 retitled ✅ shipped and rewritten around the shipped code.
- D1: v3.2 "growing" bullet. D5: "as today" corrected; implementation paragraph
  rewritten as shipped. D6: definition-time sites marked shipped.
- §7: coverage note pointing at the new test file and at which halves wait for
  Phases 5–6. §8 Phase 2 marked ✅ shipped.

## Next

**Phase 3 — pure resolver.** `find_covariance_block()` and
`resolve_correlated_draws()` in a new `R/correlated_draws.R`, tested with no
database access across every coordinate pattern and numerical edge case,
including RNG-neutrality of the nothing-to-draw cases. `.pvc_block_members()`
already does block discovery at the SQL level; Phase 3's `find_covariance_block()`
is the in-memory counterpart over a loaded set of rows (or a thin wrapper), and
its precondition — a complete, validated, PSD matrix — is what Phase 2 now
guarantees.
