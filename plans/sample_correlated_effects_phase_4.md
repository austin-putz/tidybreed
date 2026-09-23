# Correlated effects — Phase 4 results

**Spec:** `plans/sample_correlated_effects.md` (v3.4), §5.3, §5.5, §7 "Record
planning (Stage 1)" and §8 Phase 4.
**Branch:** `feat/genome-effects-v49`.
**Follows:** `plans/sample_correlated_effects_phase_3.md`.
**Status:** complete. Full suite green.
**Date:** 2026-09-21.

Phase 4 restructures `add_phenotype()` into the three stages the plan
describes — PLAN (no RNG, no writes), RESOLVE (RNG, no writes), COMMIT
(writes, no RNG) — without changing what the function computes. It is the
"largest single step" the plan warned about, and it lands as one commit. No
user-visible argument, column or message changed; what changed is *when*
things happen and what the RNG stream depends on.

## What shipped

| Piece | Where | Notes |
|---|---|---|
| `.ap_plan()` — Stage 1 | `R/add_phenotype_stages.R` (new) | Subset (sorted), metadata, path classification, topological sort, sex expression, repeatable guard, `add_tbv()` prerequisite, then per phenotype: fixed-effect terms + `null_class_action` skip, TBV (simple / composite / formula_tbv) + exclusions, `pheno_number`, random-effect level per record, residual condition value per record. Returns an in-memory plan. |
| `.ap_plan_phenotype()`, `.ap_materialize_tbvs()`, `.ap_condition_values()`, `.ap_read_by_id()`, `.ap_phenotyped_ids()` | same | Stage-1 helpers. `.ap_read_by_id()` is the one registered-view join for "rows of table X for these ids". |
| `.ap_resolve()` — Stage 2 | same | Joint named-effect pre-draw (ex-§7.5), joint residual draw (ex-§8.5), then per phenotype the marginal random-effect draws for new levels and the independent residual draw, liability, type conversion. All in memory. Draws for new random-effect levels accumulate in a pending frame that every later step consults. |
| `.ap_predraw_named_effects()`, `.ap_resolve_random_terms()`, `.ap_joint_residuals()`, `.ap_liability_records()`, `.ap_existing_draws()` | same | Stage-2 helpers. Bodies are the pre-Phase-4 draw paths, re-homed; Phases 5–6 replace them. |
| `.ap_commit()` — Stage 3 | same | `BEGIN` → register + `INSERT` into `phenotype_random_effects` → `ALTER TABLE` for extra columns → register + `INSERT` into `ind_phenotype` → `COMMIT`; `ROLLBACK` on any failure. "Wrote …" messages after the commit. |
| `add_phenotype()` | `R/add_phenotype.R` | Now ~60 lines: validation, seed, `.ap_plan()` → `.ap_resolve()` → `.ap_commit()`. Roxygen gains a "How a call runs" paragraph. `write_user_phenotype_values()` deleted (folded into the `user_values` path of Stage 1 / 2). |
| `.ap_covariate_terms()` | `R/phenotype_helpers.R` | Replaces `compute_covariate_contribution()`. Fixed effects evaluated, skip mask applied; **random effects are not drawn** — the level per record is returned for Stage 2. `phenotype_effects` rows read in `effect_name` order. |
| `next_pheno_numbers()` | same | Registered view + `JOIN` instead of `dbWriteTable()` of a temp table — RNG-neutral. |
| `.eval_derived_formula()` | `R/formula_helpers.R` | Gains `pending`: records planned earlier in the same call, treated as if written. Ids registered, not pasted. |
| Tests | `tests/testthat/test-add_phenotype_stages.R` (new) | 57 expectations; see below. |
| Docs | `man/add_phenotype_stages.Rd`, `man/dot-ap_*.Rd`, updated `add_phenotype.Rd`, `next_pheno_numbers.Rd`, `dot-eval_derived_formula.Rd`; `compute_covariate_contribution.Rd` and `write_user_phenotype_values.Rd` removed | All internal. `NAMESPACE` unchanged. |

## What a user sees differently

- **Seeded output no longer depends on physical row order.** The subset is
  sorted by `id_ind` before planning; random-effect levels and effect rows
  are sorted before any draw. Records are written in that order too
  (`id_phenotype` increases with sorted `id_ind` within a phenotype).
- **`user_values` / `user_residual` positional matching** is over that
  sorted order. The roxygen already said "`id_ind` order"; now it is true.
- **An excluded individual costs nothing.** With the same seed, the values
  written for the animals that *do* get records are identical whether or
  not an animal skipped by `null_class_action`, a missing component, or the
  repeatable guard was in the input — and a random-effect level touched
  only by an excluded animal is never drawn or stored.
- **A failed call writes nothing.** Both tables roll back. Before, a
  failure between phenotypes left the earlier phenotypes' records and any
  named-effect draws behind.
- **Condition-table lookup is strict.** A `condition_table` with zero or
  several rows for a planned individual is an error naming the table, the
  column and up to five ids; before, the first match was taken silently.
- **Named `user_values` are checked.** Every name must be a planned
  individual for that phenotype, each once; the old code wrote whatever
  names it was given, including ids that do not exist.
- **Effects on a `derived_formula` phenotype are ignored**, as the docs
  always said; the old loop evaluated them (and drew random effects for
  them) by accident of ordering.
- Seeded numeric values differ from pre-Phase-4 output (the `dbWriteTable()`
  RNG advances are gone, draws are in a fixed order, records are sorted).
  Allowed and expected — the reproducibility contract is within current
  code only.

Everything else — arguments, messages, warnings, error texts, columns —
is unchanged.

## Why the shape is what it is

**Stage 3 arrived early.** The plan scheduled the single transaction across
Phases 5–7. Once every draw is in memory the commit is ~40 lines, and
building it now meant the pre-Phase-4 draw paths could be moved into Stage 2
*as is* rather than being rewritten around their own writes and then
rewritten again. Phase 7 shrinks to the D7 test.

**Random-effect levels come from planned records, not the input subset.**
The old marginal path collected levels from the sex-filtered subset before
the skip, so a pen containing only a skipped animal still got a draw. The §7
"no stochastic state" test requires the opposite, so Stage 1 records the
level per *planned* record and Stage 2 draws only those. This is also why the
retained §7.5 pre-draw now reads its levels from the plan.

**The joint residual path compares planned sets.** The old `all_equal`
check compared pre-exclusion subsets and then indexed draws by id, so an
excluded animal still consumed a joint draw. Comparing post-exclusion sets
keeps the RNG property at the cost of a transient narrowing: two phenotypes
whose exclusions differ draw independently for that call. Phase 5 deletes
the restriction entirely, so this is a temporary cost on a path being
removed, recorded in the plan header.

**Derived formulas see pending records.** Writing everything in Stage 3
means a derived phenotype can no longer read its feeder from disk when both
are in the same call. `.eval_derived_formula()` takes the in-memory records
of the phenotypes planned before it and unions them with the disk rows;
since `pheno_number` is assigned in Stage 1, the union is exactly what the
disk would have held.

**`.ap_read_by_id()` is the one lookup idiom.** The plan's rule that planned
ids never enter SQL text is now the mechanism for the repeatable guard, the
simple-TBV read, effect source tables, the stratum lookup, `pheno_number`
and derived-formula feeders. The composite/formula TBV assembly helpers
still paste contributor ids; they are TBV code, not this plan's, and are
noted in §5.5 as out of scope.

**Stage 1 is not quite "no writes".** `add_tbv()` runs inside `.ap_plan()`
because the plan reads `ind_tbv`. It is idempotent, RNG-neutral and would be
run by the call anyway; the roxygen says so. Everything after it in Stage 1
is pure reads.

**`find_covariance_blocks()` does the stratum discovery.** Stage 1 asks the
Phase 3 loader which phenotypes sit in a conditional block and on which
`(table, column)`, then reads the values with the strict contract. The
retained §8.5 path consumes those values instead of re-querying. Phase 5
will consume the same field for the real per-entity stratum selection.

## Tests (`test-add_phenotype_stages.R`, 57 expectations)

Twelve founders, one QTL trait per phenotype, seeded construction so two
populations built the same way have identical TBVs.

**Exclusions consume no RNG (§7 items 1–2, plus the repeatable guard).**
Three paired-population tests: `null_class_action = "skip"` (with a random
effect whose level `P9` only the skipped animal touches — no `P9` row is
written, and the other level's draw is identical in both pops), composite
and formula_tbv group exclusions, and a non-repeatable second call. In each,
the records of the retained animals are `expect_equal` between the
population with the excluded animal and the one without it.

**`pheno_number` (§7 item 3).** After two prior calls, `.ap_plan()` reports
`c(3L, 2L, 1L)` for `A_1, A_2, A_5`; `add_phenotype()` writes exactly that.

**Order.** Records come back `A_1, A_10, A_11, A_12, A_2, …` with
`id_phenotype = 1:12`. `user_residual = seq_len(12)/10` recovers
`pheno_value − mean − tbv` in that order; wrong length errors with the
planned count.

**RNG accounting.** `.Random.seed` after a plain call equals the state
after `rnorm(n)`; with three new random-effect levels, after `rnorm(3 + n)`
with the level draws first and in sorted-level order; on the next call the
levels are reused and only `rnorm(n)` is consumed; a `user_values` call
leaves the seed untouched.

**Atomicity.** A reserved extra column (`pheno_value = 1`) is rejected
inside the transaction after the `phenotype_random_effects` insert has run:
both tables are empty afterwards and the next clean call succeeds.

**Derived from the same call.** `add_phenotype(c("D", "A"))` with
`D = 2 * A` writes `A` first and `D` equal to twice this call's `A` records.

**Stratum lookup contract.** A user table `ind_env` with two duplicated ids
errors naming `'ind_env'`, `'farm'`, "2 have several rows" and the ids; a
missing row errors with "1 planned individual(s) have no row … A_7"; an
unmatched level falls through to the retained unconditional-fallback
warning; `.ap_plan()` reports `condition_table`, `condition_column` and the
per-record values.

## Verification

- `test-add_phenotype_stages.R`: 57/57.
- Every phenotype-related file re-run after the last code change:
  `test-add_phenotype.R` 45, `test-formula_phenotype.R` 41,
  `test-phenotype_composite.R` 35, `test-define_effect_fixed_cov.R` 26,
  `test-phenotype_cov_block.R` 110, `test-mutate_derived.R` 37,
  `test-define_phenotype.R` 38, `test-phenotype_schema.R` 21,
  `test-correlated_draws.R` 115 — all green, no changes to any existing
  test.
- Full suite before the review pass: 3206 passed, 0 failed, 1 skipped.
  After it (final code): 3214 passed, 0 failed, 1 skipped.
- `roxygen2::roxygenise(roclets = "rd")` clean; `NAMESPACE` unchanged.

## Review pass

Second pass (full read-through for bugs and pre-1.0 cruft), on top of the
first:

- **Named `user_values` were unchecked.** The old path wrote a record for
  every name given — ids outside the subset, ids that do not exist,
  duplicates (which then shared a `pheno_number`). Names must now be
  planned individuals for that phenotype, each once; unknown or duplicate
  names error, and the named records are written in the planned
  (`id_ind`-sorted) order. Two tests added.
- **Derived phenotypes evaluated model terms.** An effect declared on a
  `derived_formula` phenotype was evaluated by the old per-phenotype loop
  as a side effect of ordering: fixed-class `NULL`s could skip individuals
  and a random effect was drawn and stored for a phenotype the docs say has
  "no mean/fixed/random contribution". The derived path now plans only
  `id_ind` / `pheno_number` and Stage 2 evaluates the formula alone. Test
  added (`.Random.seed` shows only the feeder's residuals were drawn).
- **Locale-dependent order.** `order()` / `sort()` on ids and levels used
  the session collation; `method = "radix"` (byte order) everywhere a sort
  precedes a draw, so seeded output is the same across locales.
- **Legacy guards removed** in the code this phase touched: `col %in%
  names(row)` checks for `formula_tbv`, `formula` and
  `missing_component_action` (all base columns since v0.31), `!is.null()`
  on data-frame-row fields (`type`, `mean`, `thresholds`, `cat_values`,
  `cat_names`, `null_class_action`, `poly_order`), the second
  `get_phenotype_var()` lookup in the prevalence-threshold branch (the same
  unconditional diagonal `get_residual_cov()` already returned), and the
  `requireNamespace("MASS")` check (`MASS` is in `Imports`).
- **`sample_residuals()`** lost its dead independent-draw branch and the
  `residual_var` argument: it is now `sample_residuals(n, R)`, called only
  by `.ap_joint_residuals()`, and goes entirely in Phase 5.
- `upsert_ind_tbv()` moved from `add_phenotype.R` to `add_tbv.R`, its only
  caller.
- Not touched, noted for a later cleanup: the `"NA"`-string tolerance in
  the TBV/formula helpers (`x != "NA"`), which predates `NULL` parents and
  is shared with `mutate_group_concatenate()`'s documented literal; and
  `get_residual_cov()`'s `condition_column == ""` tolerance, which Phase 5
  deletes with the function.

First pass:

- `records[[t]] <- NULL` would have *removed* the list element and shifted
  plan-order positions; the zero-record case now leaves the pre-allocated
  `NULL` in place, and the derived-formula `pending` frame is built from the
  whole records list rather than by position.
- `user_values` calls ran the named-effect pre-draw and the joint residual
  draw for nothing (the old code returned before both). Stage 2 now skips
  both when no entry is on the model path; the RNG test asserts the call is
  neutral.
- `.ap_read_by_id()` with `source_column = "id_ind"` on a non-`ind_meta`
  table would have selected `id_ind` twice; the column list drops it.
- `phenotype_name` is `unique()`d after validation — a duplicated name used
  to write the phenotype twice.
- Float-order note: TBVs are SQL aggregates, so the same animal's TBV can
  differ in the last bit between two populations built with different
  subsets; the paired tests compare with `expect_equal`, not `identical`.

## Plan bookkeeping

- Header → v3.4; "What changed from v3.3 to v3.4" (seven bullets).
- §5.5 marked shipped for Stage 1 and the stage boundaries, with a note on
  what Stage 2 currently holds; the `dbWriteTable`, id-in-SQL and stratum
  contract paragraphs annotated as done; §5.6 and §5.8 helper names
  refreshed; D2's "set to 0" note updated.
- §7 "Record planning" marked covered with the additional properties
  tested.
- §8 Phase 4 marked shipped; Phase 7 note that the transaction exists.
- §0–§2 file:line references declared historical.

## Next

**Phase 5 — residual integration.** Replace `.ap_joint_residuals()` and the
independent-draw branch of `.ap_resolve()` with the residual adapter over
`find_covariance_blocks()` + `resolve_correlated_draws()`: stratum per entity
from the plan's `condition_value`, stored `residual_value` lookup by
registered join, `user_residual` as fixed coordinates, D2/D6 checks, one
resolver call per `(block, stratum, sample set)`; write `residual_value` and
`residual_condition_level`; delete `sample_residuals()` (and its
`mvrnorm(n = 0)` crash), the equal-sets restriction and the unconditional
fallback warning; rewrite `get_residual_cov()` around strata or retire it.
Stage 1 and Stage 3 should not need to change.
