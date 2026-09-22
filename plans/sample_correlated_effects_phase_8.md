# Correlated effects — Phase 8 results

**Spec:** `plans/sample_correlated_effects.md` (v3.8), §8 Phase 8, D6
(mutability), D8 (new).
**Branch:** `feat/genome-effects-v49`.
**Follows:** `plans/sample_correlated_effects_phase_7.md`.
**Status:** complete. Full suite green. **This closes the plan** — Phases 0–8
are shipped and D1–D8 are decided, implemented and tested.
**Date:** 2026-09-22.

Phase 8 is the closing phase: documentation, housekeeping, performance. Most
of its roxygen list had already landed with the phase that introduced each
rule, so the real work was the two open **decisions** the plan deferred here,
and the benchmark that had to exist before anyone could claim the feature
scales. Both decisions went against the framing the plan offered them in, for
reasons recorded in D6 and D8.

## What shipped

**D8 — seeded output is bit-identical (new decision, code change).**

Phase 7 recorded an observation it could not act on: `add_tbv()`/`add_tgv()`
values differed by ~1e-15 between two runs of an identical population, so the
Phase 7 integrity test had to give `ind_tbv` a tolerance. The plan left it as
"decide whether the contract means bit-identical or within tolerance."

The cause was localized rather than guessed. Every step of the evaluator's one
statement is bit-stable on its own — `add_red` sums at most two allele copies
per group (floating-point addition *is* commutative; only associativity fails,
and two summands raise no associativity question), `state` sums integers,
`string_agg` carries an explicit `ORDER BY`, `product()` runs over a fixed
member count. The final `SUM(genome_value * prod)` over a trait's terms is the
only reduction with enough summands for order to matter, and DuckDB's parallel
hash aggregate combines partial sums in thread-completion order.

Decided **bit-identical**, implemented as exact `DECIMAL` accumulation
(`GEV_ACC_TYPE` in [genome_effects_eval.R](../R/genome_effects_eval.R)):
integer arithmetic, therefore associative, therefore a function of the stored
model alone. Plus `.gev_accumulator_error()`, which rethrows DuckDB's bare
conversion error as a tidybreed message when a model exceeds the accumulator's
1e20 range.

**D6 — `condition_change_action` mutability (open item, new function).**

New exported `define_condition_change_action(pop, phenotype_name, action)` in
[define_condition_change_action.R](../R/define_condition_change_action.R):
block-scoped, one transaction, one column. The D6 error message now names it
instead of the recipe it used to name, which could not work.

**Docs.** The culling example on `?add_phenotype`. CLAUDE.md: the
bit-identical paragraph in the reproducibility contract, the new function in
the `define_residual_cov()` section, the block-scoped note on the
`phenotype_meta.condition_change_action` row. `?define_phenotype` now says the
argument only sets the value while the phenotype is a block of one.

**Benchmark.** [benchmark_phenotype_scale.R](../dev/benchmarks/benchmark_phenotype_scale.R).

**Bookkeeping.** `NEWS.md` — three Phase 8 entries and the `0.71.0` heading
loses "(in development)", since the feature has now shipped. `DESCRIPTION`
stays at `0.71.0`: the version was already bumped for this feature and the
whole plan lands under it. `_pkgdown.yml` gains the new export — and, while
there, the five topics that were already missing from the reference index
(`ad_terms`, `add_tgv`, `define_genome_effects`, `extract_allele_freq`,
`genotype_terms`, all from the genome-effects work), because
`pkgdown::check_pkgdown()` fails on the whole set or none of it, so the new
entry could not be verified in isolation while they were outstanding.

## What a user sees differently

| Before | After |
|---|---|
| Two runs of the same seed on the same population gave `tbv_value`, `tgv_value` and `pheno_value` differing in the last bits | Bit-identical. `expect_identical()` holds |
| A model with astronomically scaled `genome_value` produced a huge or `Inf` genetic value | Errors above 1e20, naming `genome_value` as the usual cause |
| `condition_change_action` was frozen once its residual block had two defined members | `define_condition_change_action()` changes it on the whole block |
| The D6 error suggested `define_phenotype(..., overwrite = TRUE)` on every member — a sequence that cannot succeed | It names the writer that works |

## Why not the options the plan offered

**D8.** "Identical within tolerance" (option 1) was tempting only while the
wobble looked cosmetic. It is not: a breeding simulation is chaotic downstream
of selection, so once `select_parents()` exists a last-bit difference in a TBV
will occasionally flip a truncation ranking between two nearly tied animals,
and the two runs then produce different pedigrees rather than different last
bits — reported as a bug in selection, not traced to a `SUM()`. Of the three
deterministic options, `SET threads = 1` serializes the expensive join along
with the cheap aggregate, and an ordered `list_reduce()` changes the aggregate
state from one value per group to one per row, regressing the memory profile
of the package's hottest query and conflicting with design principle 2
("larger than RAM"). Exact accumulation has neither cost.

**D6.** The plan offered a "set on all members" form of `define_phenotype()`
or a relaxation letting the *last* member's change through. The relaxation is
order-dependent and leaves the database in the disagreeing state between calls,
where any `add_phenotype()` fails — the mid-sequence invalid state D1's
whole-block rule exists to prevent. Overloading `define_phenotype()` looks
smaller than it is: `overwrite = TRUE` replaces the whole `phenotype_meta` row,
so flipping one flag means restating `type`, `mean`, `expressed_sex` and the
rest, and forgetting one silently resets it. A one-column change must not be
spelled as a re-registration. The writer matches the scope of the property
instead.

## Measurements

Determinism, `add_tbv()` on an unchanged population, 6 runs each:

| loci | individuals | before | after |
|---|---|---|---|
| 60 | 12 | diverged, max abs 8.9e-16 | identical |
| 500 | 100 | diverged, max abs 1.3e-15 | identical |
| 2000 | 200 | diverged, max abs 2.2e-15 | identical |

Under `SET threads = 1` the old path was already identical, which is what
pointed at the parallel aggregate.

Cost of the exact accumulator against the plain `SUM()` (min of 3, the same
registered inputs):

| loci | individuals | `SUM` | exact | ratio |
|---|---|---|---|---|
| 500 | 200 | 0.046s | 0.050s | 1.09× |
| 2000 | 500 | 0.196s | 0.202s | 1.03× |
| 5000 | 1000 | 1.046s | 1.006s | 0.96× |

The join dominates, so the difference is inside run-to-run noise. An ordered
`list_reduce()` measured 1.26–1.32× on the same inputs *and* holds one value
per row instead of one per group — the reason it was not chosen.

`add_phenotype()` scale (`TIDYBREED_BENCH_LARGE=1`, 1000 loci, 200 QTL,
16 threads). Total cost per individual, call 1 / call 2:

| shape | n = 1,000 | n = 4,000 | n = 16,000 |
|---|---|---|---|
| independent | 0.47 / 0.90 | 0.133 / 0.108 | 0.095 / 0.099 |
| correlated | 0.31 / 0.35 | 0.112 / 0.134 | 0.085 / 0.090 |
| conditional | 0.34 / 0.37 | 0.139 / 0.131 | 0.101 / 0.107 |
| named_effect | 0.34 / 0.38 | 0.117 / 0.160 | 0.102 / 0.098 |

(ms/ind; the n = 1,000 row carries the process's warm-up, which is why the
first cell of the first shape is the largest number in the table.)

**No optimization was needed.** Per-individual cost falls with population size
in every shape — roughly 3–4× from 1,000 to 16,000 — and the stage split says
why nothing here is worth tuning:

- `plan` dominates the total everywhere (1.2–1.5 s of a ~1.6 s call at
  n = 16,000). It is the stage with no RNG and no writes, so it is also the
  stage an optimization could touch most safely; it is simply not the
  bottleneck anyone feared.
- `commit` is 17–28 ms at n = 16,000 against 8–13 ms at n = 1,000 — a 16×
  row increase for under a 3× cost. The single transaction is batching.
- `resolve` — which contains the observation-pattern query the plan singled
  out — is 0.055 s (independent) to 0.26 s (conditional) at n = 16,000, and
  the call-2-minus-call-1 gap stays in the tens of milliseconds rather than
  growing per entity. Reading back already-realized residuals is a constant
  number of statements, as designed.

## Tests

| File | What it pins |
|---|---|
| `test-genome-effects-determinism.R` (new, 18 expectations) | `add_tbv()` and `add_tgv()` bit-identical over 8 repeated runs of one population; a seeded `add_phenotype()` bit-identical between two identically built populations, asserting `residual_value` (never the problem) and `pheno_value` (was) separately; the accumulator-range error |
| `test-phenotype_cov_block.R` (109 \u2192 125 expectations, four new tests) | The D6 deadlock itself — both single-member flips refused — then the block-scoped write, both directions; that only `condition_change_action` changes and every other `phenotype_meta` column survives; idempotence; that realized draws lock the matrix but not the action; unknown phenotype, bad action, non-scalar name; a block of one |
| `test-add_phenotype_failure_contract.R` (modified) | `expect_db_unchanged()` drops the `ind_tbv` tolerance Phase 7 needed and compares every table with `expect_identical()` |

**Pre-change check** (the plan's "tests before the fix" rule, §8): with
`R/genome_effects_eval.R` stashed, `test-genome-effects-determinism.R` reports
**16 failed / 2 passed**; with it restored, 18/18. The two that pass without
the change are the `expect_length()` guard and the `residual_value` assertion —
both true before and after, which is the point of asserting them separately.

**Verification.** Full suite on the final code: **3463 passed, 0 failed,
1 skipped** (3434 at Phase 7; +18 determinism, +16 cov-block, −5 from the
failure-contract helper collapsing two expectations into one).
`pkgdown::check_pkgdown()` clean.

## Adversarial probes

| Probe | Result |
|---|---|
| Does a plain `SUM()` over many rows per group vary run-to-run on synthetic data? | **No** — a single registered chunk partitions deterministically. The synthetic reproduction failed, which is why the real query was instrumented stage by stage instead |
| Which CTE is order-dependent? | Only the final `SUM`. `add_red`, `state`, `product()` and the R-side ordered sum over the same per-term rows are all bit-stable |
| Is single-threaded `SUM` the same value as an R-side ordered sum? | Equal, **not** identical — both are *an* order. Bit-identity needs a canonical reduction, not merely a serial one |
| What does the exact accumulator do at its limits? | 1e19 fine; 1e20, an overflowing accumulation, `Inf` and `NaN` all raise DuckDB conversion errors — hence `.gev_accumulator_error()`. `NULL` passes through unchanged, so the dominance-mismatch `n_bad` path is untouched |
| Would a conservative pre-check be better than a rethrow? | **No.** The only cheap bound is `\|genome_value\| × 2^n_members` summed over terms, which rejects a legal 70-member interaction. Not pre-checked |
| Does `--check` mode actually catch a semantic change? | It did, immediately — on itself. The first version left `build_pop()` unseeded, so the founder haplotypes differed per run and `pheno_value` moved by whole units. Seeded; three runs now diff clean |

## Review against the plan

Every §8 Phase 8 bullet, item by item:

- **Roxygen list (11 items).** Ten had already landed with the phase that
  introduced the rule — verified in the source, not assumed: the
  `user_residual` subset-list contract at `add_phenotype.R:114-129`
  (Phase 5), the persistent-level note at `define_effect_random.R:10`
  (Phase 6), D7 (Phase 7), D1/D2/D3/D5/D6 in `define_phenotype.R` and
  `define_residual_cov.R` (Phases 1–2). The **culling example** was the one
  genuine gap — culling was described in prose at `add_phenotype.R:48` but had
  no runnable example. ✅
- **CLAUDE.md (4 items).** All four sections read as current; added the three
  Phase 8 items. ✅
- **`R/schema.R` descriptions / `test-schema-print.R`.** Unchanged since
  Phase 1; green. ✅
- **`NEWS.md` / `DESCRIPTION`.** ✅ (Version held at `0.71.0` deliberately —
  see "What shipped".)
- **D6 mutability.** ✅ Decided against both options the plan listed, with the
  reasoning written into D6 rather than only here.
- **Benchmark.** ✅ And the "without changing RNG semantics" clause is now
  enforceable rather than aspirational: `--check` is the diff target.
- **`add_tbv()` bit-reproducibility.** ✅ Resolved as D8 rather than
  documented as a limitation.

## Plan bookkeeping

`plans/sample_correlated_effects.md` → **v3.8**: status header (complete,
Phases 0–8, D1–D8); a "What changed from v3.7 to v3.8" block; D6 marked ✅ with
the new **Mutability** subsection; new **D8** section with cause, options,
decision and consequences; D7's closing note updated to say the `add_tbv()`
upsert is bit-identical too; §8 Phase 8 marked ✅ with its shipped list and the
original scope kept for the record.

## Review pass over Phases 0-8

A full re-read after Phase 8 closed. Two questions drove it: is anything dead,
and does the suite actually prove residual covariances are reproduced?

### The test gap that mattered

**Every multi-phenotype residual test used a unit-diagonal matrix.** `R_AB` is
`c(1, .8, .8, 1)`; the three-phenotype matrix has 1s on the diagonal; the
stratified tests use `R_AB` and `4 * R_AB`, which are scaled correlation
matrices. In all of them the conditional slope

```
beta = sigma_AB / var_A
```

is numerically equal to the correlation `sigma_AB / sqrt(var_A var_B)`. A
resolver that conditioned on the correlation instead of the covariance would
have passed **every one of those tests**. The 1x1 strata (variances 4 and 9)
covered `sqrt(v)` scaling for a single coordinate, but nothing covered the
two-coordinate case where the two quantities differ.

The implementation turned out to be **correct** — checked directly with
`R = [[4, 1.8], [1.8, 9]]`, where `beta = 0.45` and the correlation is `0.30`:
the stored residuals match `0.45 * A + sqrt(8.19) * z` exactly and do not match
the correlation form. But it was correct untested. Four tests now close it:

| Test | What it distinguishes |
|---|---|
| unequal variances, sequential | exact `beta` and conditional sd, plus an explicit `expect_false` that the correlation form was used |
| unequal variances, joint call | the same block through one `chol()` in a single call |
| three phenotypes, three variances | the whole matrix recovered by the sequential path, 4 SE per entry |
| strata with different unequal matrices | each stratum conditions on its own covariance, not a shared or scaled one |

### Edge cases probed (no defects)

| Probe | Result |
|---|---|
| Rank-1 block `[[4,4],[4,4]]` | `B == A` exactly (max abs difference 0), and the call still consumes its `n` normals — stream position is independent of the matrix's rank |
| Zero-variance member | residual is exactly `0`; the other coordinate keeps its variance |
| Non-PSD matrix | refused at definition time, naming the smallest eigenvalue |
| `dbWriteTable()` reachable from `add_phenotype()`? | No. Stage 1's `add_tbv()` upserts, Stage 3 registers + `INSERT`s |
| ids pasted into SQL in the phenotype layer? | No. Every `paste` touching `id_ind` builds an error message |
| `set.seed()` placement | Before Stage 1, deliberately: any RNG use added to planning would shift the draws and fail the "consumes exactly its draws" test loudly. Moving it after Stage 1 would hide that, so it stays |
| Duplicate DDL (Phase 0's claim) | Still one `CREATE TABLE` per table |
| Stale references to functions Phases 0-6 deleted | Only in `NEWS.md`, where they belong |

### Dead code removed

Found by `codetools::checkUsagePackage()` and a namespace walk for functions
and formals never referenced from any body:

| Removed | Why |
|---|---|
| `phenotype_components.missing_action` | Write-only. Duplicates the decided single `phenotype_meta.missing_component_action`; `.ap_missing_action()` reads the latter |
| `phenotype_components.contributor_filter` | Write-only, reserved for a spatial lookup with no counterpart anywhere |
| `update_covars_from_blupf90()` + its man page | A stub that ignored all three arguments and printed a message. `add_ebv()` prints it directly now |
| `keep` parameter of `.create_run_dir()` | Never read; documented as "reserved for a future release" |
| `effects_df` parameter of `write_meta_file()` | Never read |
| `user_pos` in `.schema_group_of()`, `pks` in `mutate_group_concatenate()` | Dead locals (the sibling `mutate_group_*` writers do use `pks`) |

`phenotype_components.component_names` was **kept** — it is the one reserved
column here whose counterpart, `ind_tgv.component_name`, is already written by
`add_tgv()`.

### Errors moved to the call that causes them

`define_phenotype(components = )` validated neither `group_column` nor
`group_table`, though both reach SQL through `.group_members_sql()` and
`define_effect_random()` / `define_effect_fixed_cov()` validate their
equivalents. The read-time existence checks in `.read_one_per_id()` close the
injection risk, so this was an asymmetry rather than a hole — but the error
arrived at `add_phenotype()` time instead of at the definition. Both are now
validated, `source_table` is validated like `source_column` in both effect
writers, and a `weight_type` other than `"fixed"`/`"covariate"` is rejected at
definition rather than accepted into the table and refused later.

`poly_scale_min` / `poly_scale_max` are read nowhere, since they exist for the
rejected `"legendre"` / `"raw_poly"` weight types. They were left in place:
`poly_order` is a coherent reserved dimension across three tables for the
random-regression work in this plan's §10, and removing half of it would be
worse than leaving it whole. Flagged rather than cut.

## Next

Nothing in this plan. The correlated-effects feature is complete.

Two things it leaves pointing forward, neither blocking:

1. **Selection must not break ties on raw float comparison.** D8 removes the
   wobble that would have made this bite, but `select_parents()` should still
   order on an explicit, total key rather than on `>` between two doubles.
2. **`plans/sample_correlated_effects.md` §9 and §10** stand as written: the v1
   limitations to document (deletion/regeneration coherence, informative
   missingness, count traits as transformed Gaussian, persistent level
   identity, concurrency) and the Layer 2/3 models deliberately out of scope.
