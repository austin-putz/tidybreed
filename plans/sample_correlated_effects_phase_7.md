# Correlated effects — Phase 7 results

**Spec:** `plans/sample_correlated_effects.md` (v3.7), D7 (RNG state on
failure), §5.5 (Stage-3 transaction), §7 "Reproducibility and integrity"
items 3–5, §8 Phase 7.
**Branch:** `feat/genome-effects-v49`.
**Follows:** `plans/sample_correlated_effects_phase_6.md`.
**Status:** complete. Full suite green.
**Date:** 2026-09-21; review pass 2026-09-22.

Phase 7 is the smallest phase: D7 was decided as *option 1* ("database
atomic; RNG advances on failure"), which is the behaviour every RNG-consuming
function in the package already has. The Stage-3 transaction has existed
since Phase 4. What was missing was the *statement* of the contract in the
docs and the *two-part test* the plan asks for — database unchanged, seed
advanced by exactly the draws made — including the Stage-2-failure half that
no test covered. **No code under `R/` changed behaviour.**

## What shipped

- **Verified, not built.** Nothing in `R/` reads or writes `.Random.seed`,
  and there is no `RNGkind()`, `withr::with_seed()` or `with_preserve_seed()`
  anywhere; the only `set.seed()` calls are the documented `seed =`
  arguments of `add_phenotype()`, `define_additive_effects()`,
  `define_trait_simple()` and the two `mutate_group_*()` helpers, and
  `recombination_helpers.R` seeds `dqrng` *from* the base stream. No
  `dbWriteTable()` is reachable from `add_phenotype()` (the §5.5 rule; the
  Stage-1 `add_tbv()` upsert is a column-named `INSERT … ON CONFLICT`, and
  a probe confirms it leaves `.Random.seed` identical). `.ap_commit()` is
  the only writer after Stage 1, and it is one `BEGIN … COMMIT` with
  `ROLLBACK` in `on.exit`. So D7 option 1 holds by construction; Phase 7
  pins it.
- **Docs.** `?add_phenotype` gains an "If a call fails" paragraph;
  `?add_phenotype_stages` gains an "On failure" paragraph naming the three
  failure points (Stage-1 rejection, Stage-2 error after some draws, failed
  Stage-3 write), the atomic write, the non-rewound stream, the retry
  consequence and the one write that remains (the Stage-1 `add_tbv()`
  upsert). CLAUDE.md's `add_phenotype()` entry gets a "Failure contract
  (D7)" paragraph with the rule *never add seed restoration to one
  function*.
- **Tests.** New `tests/testthat/test-add_phenotype_failure_contract.R`
  (38 expectations). Every integrity assertion is a **whole-database
  snapshot**: `db_snapshot()` reads every base table in `main` ordered by
  all columns; `expect_db_unchanged()` requires identity for every table
  except `ind_tbv`, which is compared with tolerance (see "Observation"
  below). Every RNG assertion replays the draws made before the error with
  `stats::rnorm()` and compares `.Random.seed`.

## What a user sees differently

Nothing at runtime. The documented promise is now explicit:

- A failed `add_phenotype()` never writes a phenotype record or a
  random-effect draw, whatever stage it fails in — including a block that
  was fully resolved before a later block failed.
- Re-running after a failure draws **different** values unless the call is
  re-seeded (`seed =` or `set.seed()`); the docs say so and say why (no
  function in the package rewinds the stream).
- The TBVs the failed call materialized stay in `ind_tbv` — the contract
  is "no phenotype record and no draw", not "no write at all". They cost no
  RNG and the retry rewrites them (up to summation order — see below).

## Tests (`test-add_phenotype_failure_contract.R`, 38 expectations)

| Test | What it pins | Draws before the error |
|---|---|---|
| Stage-3 write failure (reserved extra column, rejected after the `phenotype_random_effects` INSERT ran in the same transaction) | every table identical; seed = `rnorm(1 + 12)` (pen `P1`, then 12 residuals); a retry **without** re-seeding gets pen = `2·z[14]` and residuals `z[15:26]`, not `2·z[1]`; a retry **with** `seed = 41` reproduces `2·z[1]`, `z[2:13]` on a fresh pop | 1 + 12 |
| Stage-2 error in the **named-effect adapter** (`herd` variance 1 with two new levels; `pen` variance rows deleted by hand; `'herd' < 'pen'`) | error "No variance stored for random effect 'pen'"; every table identical; seed = `rnorm(2)` — the two herd levels were drawn and stay consumed | 2 |
| Stage-2 error in the **residual adapter** (D2: `B` stored under farm `F1`, two animals moved to `F2`, then `add_phenotype(c("A", "C"))` with a new pen level on `A`) | error names the 2 records in block `{B, C}`; every table identical — **the pen level and all 12 of block `{A}`'s residuals were drawn, `A`'s records were fully assembled, and none of it was written**; seed = `rnorm(1 + 12)` | 1 + 12 |
| Schema rollback (same test as row 1: `...` names a new column `farm_x` **and** the reserved `pheno_value`; `prepare_extra_cols()` runs per field, so the `ALTER TABLE ADD COLUMN` is issued before the refusal) | `ind_phenotype` has the same columns afterwards — DuckDB's DDL is transactional, so "unchanged" covers the schema | 1 + 12 |
| The one write a failed call leaves | with `ind_tbv` emptied first, a failed call leaves 12 TBV rows and still no records and no draws — the Stage-1 upsert is outside Stage 3 by design, costs no RNG, and the retry rewrites it | 1 + 12 |
| Rejections before any draw | unknown phenotype (Stage 1); `user_residual` of the wrong length (Stage 2, before the adapters); a fixed value off the support of a singular `{A, B}` block (inside the resolver, before `B`'s first normal) — each leaves the seed *identical* and the database unchanged | 0 |

The third test is the one the plan's v3.4 note flagged as missing ("the
Stage-2-failure half of D7"): it shows the atomicity is per *call*, not per
block — the residual adapter's block loop resolves `{A}` completely, then
`{B, C}` errors, and `{A}`'s records never reach the database.

**Pre-change check.** Stronger than a run against the previous commit: the
whole `R/` diff of this phase is roxygen comments
(`git diff R/ | grep -v "^[-+]#'"` is empty), so the tests provably assert
behaviour that already held at `40bc6d9` — that is the point of the phase. The Stage-3 test in
`test-add_phenotype_stages.R` ("Stage 3 is atomic") is left as it was, with
a pointer comment to the new file, which is the D7 home.

**Also added, in `test-add_phenotype_stages.R`:** "seeded output is
identical when ind_meta's physical row order is reversed" — §7 item 4,
which the review found was *not* covered (see below). Two pops built from
the same seed, one with `ind_meta` rebuilt in reverse `id_ind` order
(`CREATE TABLE AS … ORDER BY id_ind DESC`, then drop and rename — nothing
declares a SQL foreign key to `ind_meta`), the same seeded call on both,
records and draws compared. Both adapters run.

## Observation outside the plan's layer: `add_tbv()` is not bit-reproducible

Writing `db_snapshot()` surfaced this: `ind_tbv` is the one table a failed
call touches (the Stage-1 `add_tbv()` upsert), and its `tbv_value`s differ
across identical runs at the last bit — three consecutive `add_tbv("A")`
calls on the same 12 animals / 60 loci gave three different vectors,
max |Δ| ≈ 5e-16, and become identical after `SET threads = 1`. The cause is
DuckDB's parallel `SUM()` in the genome-effects evaluator: floating-point
summation order changes with the thread schedule.

This is not Phase 7's concern and was not changed. It is recorded in the
main plan under §8 Phase 8 as an open item, because it bears on CLAUDE.md's
"same seed reproduces within the current code" contract: `pheno_value`
inherits the wobble through the TBV term (existing tests use
`expect_equal`, so nothing fails). The options are an ordered reduction in
the evaluator (a single-threaded aggregate, or an R-side ordered sum — an
`ORDER BY` inside the subquery is not enough, DuckDB does not promise
aggregation order) or stating the contract as "identical within tolerance".
`expect_db_unchanged()` compares `ind_tbv` with tolerance for this reason.

## Adversarial probes (review pass)

Four probes beyond the tests, each run against the shipped code:

| Probe | Result |
|---|---|
| A new `...` column is ALTERed in, then a second `...` field is refused inside the same transaction | Column gone afterwards, later calls work → **tested** (DDL rolls back) |
| A failed call on a pop with **no** TBVs yet | `ind_tbv` gains 12 rows, `ind_phenotype` and `phenotype_random_effects` stay empty → **tested and documented** (the contract is "no record, no draw", not "no write") |
| Temporary registered views (`__ap_records`, `__ap_levels`, …) after a failure | None leaked — every `duckdb_register()` has its `on.exit` unregister |
| `add_tbv()` RNG-neutrality | `.Random.seed` identical across the call — the reason the Stage-1 write does not break the seed accounting |

No probe found a defect; two turned into tests because they pin facts the
docs now assert.

## Verification

- `test-add_phenotype_failure_contract.R`: 38 pass.
- `test-add_phenotype_stages.R`: 60 pass (57 + the physical-order shuffle).
- Full suite (`suite5.R`, `load_all` + `test_dir`): **3434 passed, 0
  failed, 1 skipped** on the final code (3422 after the first pass, +9 from
  the two review tests and +3 from the physical-order shuffle).
- `roxygen2::roxygenise(roclets = "rd")`: `man/add_phenotype.Rd` and
  `man/add_phenotype_stages.Rd` regenerated; `NAMESPACE` unchanged.

## Review against the plan

Checked item by item:

- **D7 option 1** — "database atomic, RNG advances on failure; documented
  and tested". Documented in two Rd pages and CLAUDE.md; tested for every
  failure point. No seed restoration anywhere. ✅
- **§7 item 3** — Stage-3 failure rolls back residuals, named draws and
  records together **and** seed advanced by exactly the Stage-2 draws, "both
  halves asserted". ✅ (whole-database identity + `rnorm(1 + n)` replay).
- **§7 item 5** — RNG accounting. Already asserted in Phase 4
  (`test-add_phenotype_stages.R`); Phase 7 adds the accounting *through* a
  failure. ✅
- **§7 item 4** — row-order independence. **Was not covered.** My first
  pass marked it ✅ against Phase 4's "records are planned and written in
  id_ind order" test, which asserts the *output* order and never shuffles
  the input — the item explicitly asks to "shuffle `ind_meta` physical
  order between two seeded runs and compare". Corrected: the shuffle test
  is written (above) and the plan now cites it and says the Phase 4 test
  is not the same thing. ✅
- **§8 Phase 7 (v3.4 note)** — "what remains is the Stage-2-failure half of
  D7 and the seed-advanced assertion". Both shipped. ✅
- **§8 Phase 8 roxygen list** — "D7 RNG contract" is done here; marked in
  the plan so Phase 8 does not redo it.
- **§5.5 idiom** — `.ap_commit()` follows the `define_chromosome.R`
  transaction idiom; unchanged. Re-verified alongside it: no
  `dbWriteTable()` is reachable from `add_phenotype()`, and no planned id
  reaches SQL text — the two other §5.5 rules the seed accounting depends
  on. ✅
- **§5.5 "what the transaction does not cover"** — annotated in the plan
  with what Phase 7 found it *does* cover (the schema) and the one write
  outside it (Stage 1's `add_tbv()`). ✅

One item (§7.4) was found uncovered during the review and is now covered;
nothing else in Phase 7's scope is outstanding. Two things were deliberately
*not* done: adding a hint to the error messages ("the RNG has advanced") —
the contract is a property of every failed R call and belongs in the docs,
not in every `stop()`; and making the Stage-1 `add_tbv()` write part of the
Stage-3 transaction — it is a prerequisite the plan places in Stage 1 on
purpose (§5.5), it is RNG-neutral and idempotent, and moving it would put a
long-running write inside the transaction that holds the phenotype tables.

## Plan bookkeeping

`plans/sample_correlated_effects.md` → v3.7: status header (Phases 0–7
shipped); "What changed from v3.6 to v3.7"; D7 heading ✅ with the shipped
note; §7 items 3–5 ✅ with test names and item 4 corrected (see above);
§8 Phase 7 ✅; §8 Phase 8 roxygen
list marks "D7 RNG contract" and the persistence note as done; new Phase 8
open item for `add_tbv()` bit-reproducibility. `NEWS.md` Phase 7 entry
under 0.71.0 → Changed. CLAUDE.md `add_phenotype()` failure-contract
paragraph. Version stays `0.71.0` (the in-development heading NEWS
accumulates under).

## Next

Phase 8 — documentation, housekeeping, performance: the remaining roxygen
items (culling example, `user_residual` subset-list contract — check what
Phases 5–6 already cover), the D6 `condition_change_action` mutability
decision, the `dev/benchmarks/` script, the `add_tbv()` reproducibility
decision above, and the final 0.71.0 bookkeeping.
