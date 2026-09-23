# Correlated effects — Phase 6 results

**Spec:** `plans/sample_correlated_effects.md` (v3.6), §5.4 (adapter
responsibilities), §5.5 Stage 2, §5.6, §5.8, §7 "Named effects",
"Defect 3", and §8 Phase 6.
**Branch:** `feat/genome-effects-v49`.
**Follows:** `plans/sample_correlated_effects_phase_5.md`.
**Status:** complete. Full suite green.
**Date:** 2026-09-21.

Phase 6 puts named random effects through the machinery Phase 5 built for
residuals, with a different entity key. A level of a random effect — a pen,
a herd, an `id_ind` for a permanent-environment effect — is a persistent
entity: its draw is realized once, the first time a planned record touches
it, and reused by every later record with that level. In a covariance block
a level's draw for one phenotype is now conditional on the draws it already
has stored for the block's other phenotypes, whether the phenotypes arrive
in one call or a season apart. Defect 3 (a level with a stored `ADG` draw
had its `(ADG, BF)` pair redrawn jointly and half of it discarded) is
closed. Nothing of the pre-Phase-4 draw code remains in `add_phenotype()`.

## What shipped

| Piece | Where | Notes |
|---|---|---|
| `.ap_named_effect_targets()` | `R/add_phenotype_stages.R` | Stage 1: the model-path phenotypes whose planned records carry a random term for each effect, keyed by byte-sorted `effect_name`. |
| `.ap_plan()` | same | Loads the named-effect blocks once — one `find_covariance_blocks()` call per effect, targeting only those phenotypes — into `plan$named_blocks`, beside `residual_blocks`. |
| `.ap_resolve_named_effects()` | same | Stage 2 entry point of the named-effect adapter. Effects in sorted order; a phenotype with a random term but no block is "No variance stored for random effect"; per block `.ap_named_effect_block()`; returns the summed contribution per planned record and the new `phenotype_random_effects` rows. |
| `.ap_named_effect_block()` | same | One block: integrity check (a named block has exactly one, unconditional stratum) → `validate_named_effect_block(caller = "add_phenotype()")` (the §5.6 backstop) → planned levels of each in-call phenotype → sorted entities (levels) → stored draws of every block member at those levels via a registered view → sample set per level (planned and not stored) → one `resolve_correlated_draws()` call per sample-set group in sorted order, or the marginal `rgamma`/`runif` sampler for a 1 × 1 non-normal block → new rows and per-record contribution (`0` for a `NULL` level). |
| `.ap_resolve()` | same | The pre-draw and the per-phenotype marginal loop are replaced by one `.ap_resolve_named_effects()` call; still before the residual adapter. |
| `.blupf90_residual_cov()` | `R/blupf90_helpers.R` | `write_renum_par()`'s residual matrix, from `find_covariance_blocks()`: block-diagonal across independent blocks, error for a trait in no block or with only conditional strata. |
| Deleted | `R/add_phenotype_stages.R`, `R/define_effect_cov_matrix.R` | `.ap_predraw_named_effects()` (§7.5, the last `MASS::mvrnorm()` on the phenotype path), `.ap_resolve_random_terms()`, `.ap_existing_draws()`, `.ap_append_random_effects()`, `load_phenotype_cov()`. Four `man/` pages go with them. |
| Tests | `tests/testthat/test-add_phenotype_named_effects.R` (new) | 74 expectations (65 at the first commit, +9 in the review pass); see below. One test added to `test-correlated_draws.R` for `.blupf90_residual_cov()`. No existing test changed. |
| Docs | `man/dot-ap_resolve_named_effects.Rd`, `dot-ap_named_effect_block.Rd`, `dot-ap_named_effect_targets.Rd`, `dot-blupf90_residual_cov.Rd` new; `add_phenotype.Rd` (random-shift bullet), `add_phenotype_stages.Rd` (Stage 2 = two adapters), `dot-ap_resolve.Rd`, `correlated_draws.Rd`, `define_effect_random.Rd` (the §5.8 persistence note, the conditional-sampling paragraph, the `distribution` semantics), `define_effect_cov_matrix.Rd`, `write_renum_par.Rd` updated | `NAMESPACE` unchanged. CLAUDE.md: `phenotype_random_effects` (persistent entity, written only by `add_phenotype()`), the Stage 2 description, the `define_effect_random()` bullet. |

## What a user sees differently

- **Named effects are correlated across calls.** `add_phenotype("ADG")`
  today draws pen `P1`'s `ADG` marginally and stores it; `BF` stays latent.
  `add_phenotype("BF")` next season draws `P1`'s `BF` as
  `E[BF | ADG] + conditional noise` from the declared covariance. Before,
  the two single-phenotype calls drew independently — the covariance was
  used only when both phenotypes were in one call.
- **A stored draw is never redrawn.** A call naming both phenotypes on a
  level that already has one of them keeps the stored value exactly and
  conditions the other on it (Defect 3). Levels that have neither draw
  jointly; the two groups are separate resolver calls.
- **Mixed and partial patterns.** Levels with a stored `A`, levels with a
  stored `(A, B)`, and brand-new levels can all appear in one call for `C`;
  each conditions on what it has.
- **The §5.6 checks run at sampling time too.** A `phenotype_effects` row
  changed after the block was declared (distribution, `effect_class`,
  `source_column`) is refused by `add_phenotype()` with the same message
  as the writers, before any draw.
- **Persistence is documented where the level is declared.**
  `define_effect_random()` now says that a level is forever and that a
  per-batch pen is `pen_batch`, next to its example.
- **BLUPF90 residual prior.** `write_renum_par()` reads the unconditional
  residual blocks; two traits in different blocks get an explicit `0`
  covariance (they are independent by declaration), and a trait with only
  conditional strata is an error rather than whichever row the old
  unfiltered query returned last.
- Nothing else. Fixed effects, `user_values`, derived formulas, the
  residual adapter and Stage 3 are untouched; the Phase 4 RNG-accounting
  test (`3 + n` normals, sorted levels, `rnorm(3, sd = 2)`) passes
  unchanged because a 1 × 1 normal block through the resolver *is*
  `sqrt(v) · z`.

## Why the shape is what it is

**The residual adapter, minus what does not apply.** Entity = level;
coordinates = the block's phenotypes (including members not in the call
and members with no random term for the effect — they can be stored, so
they condition; they are never sampled); stored = `phenotype_random_effects`
joined on a registered view of the planned levels; no strata (a named block
is always unconditional, and finding one with strata is an integrity
error); no fixed coordinates; no D2/D6. Grouping is by sample set only.
The function is about half the length of `.ap_residual_block()` for that
reason, and reads the same way.

**Targets are the phenotypes with a random term, not the phenotypes of the
call.** `find_covariance_blocks()` is asked about the phenotypes whose
planned records carry the effect; a call phenotype that is in the block
but has no `phenotype_effects` row for the effect is then a latent
coordinate like any absent member. That is what makes "drop the term, keep
the draws" (the latent-coordinate test) behave: the stored draws still
condition, nothing new is drawn for it.

**Blocks loaded once, in Stage 1.** The same reason as the residual blocks:
Stage 2 should not be reading tables to decide what to draw. One loader
call per effect, effects byte-sorted, so the RNG order is a function of
the model and the plan.

**`normal` 1 × 1 goes through the resolver.** The plan allowed the marginal
sampler for any 1 × 1 block; keeping it only for `gamma`/`uniform` means
there is exactly one Gaussian path. For one coordinate and nothing
observed the resolver draws `chol(v) · z = sqrt(v) · z`, which is what
`rnorm(sd = sqrt(v))` is, so the change is invisible to seeded output.

**No pending-draw bookkeeping.** The old code merged on-disk and
this-call draws per `(phenotype, effect)` because the pre-draw and the
marginal path could both touch one phenotype. Each `(effect, block)` is
now visited once and a phenotype is in one block per effect, so a level is
drawn at most once per call by construction and the merge is gone.

**Retire `load_phenotype_cov()` now.** The plan tied its removal to Phase 6.
Its last caller was the BLUPF90 parameter writer, and the reader was
wrong for the current schema (no `condition_column IS NULL` filter, so a
stratified residual would feed a conditional row's value into the prior).
The block loader is the one reader of stored matrices; the writer now
uses it.

## Tests (`test-add_phenotype_named_effects.R`, 74 expectations)

Exact tests replay the resolver's contract: `n × m` standard normals in
entity (byte-sorted level) order and coordinate order within an entity,
`draws = z %*% chol(C)`; named-effect draws precede the residual draws of a
call, so the residual stream is checked as `rnorm(n)` *after* the level
normals. Contributions are checked as
`pheno_value − mean − tbv − residual_value` per record against the stored
level draw.

| Test | Covers |
|---|---|
| day 0 / day 100 / day 200 (§5.8 exactly) | marginal draw for `A` with `B` latent; `B` conditional on stored `A` (`−12/150 · a + √3.04 · z`); a reuse draws nothing; contributions equal the level draw; RNG accounting for all three calls |
| mixed patterns | stored and new levels in one `B` call: conditional for `P01–P04`, marginal `2z` for `P05–P08`, one entity list, `A` still latent for the new pens |
| Defect 3 closed | both phenotypes in one call: `{A, B}` group (new level, joint through `chol(R)`) before `{B}` group (stored `A` reused exactly, `B` conditional); group order by byte-sorted sample key |
| three-phenotype block | `C` conditional on `(A, B)` for two levels and on `A` alone for a third, in one resolver call; hand-computed `E[C | ·]` and `V[C | ·]` |
| effect order, `NULL` level | `'herd'` before `'pen'` before residuals; a `NULL` pen draws nothing and contributes `0`; four rows stored |
| all-`NULL` levels | a phenotype none of whose planned records has a level: nothing drawn for it as a 1 × 1 block, and in a block with a partner that has levels the partner draws marginally while its coordinate stays latent (review-pass regression) |
| permanent environment | `source_column = "id_ind"`: one draw per animal, reused on the second record, second call draws residuals only |
| 1 × 1 gamma / uniform | `rgamma(shape = 1, rate = 1/√v)` and `runif(±√3v)` exactly; a new level on the next call is one new gamma draw, stored levels reused |
| latent member | the term dropped from `A` after its draws exist: `A` no longer draws, contributes `0`, and its stored draws still condition `B` |
| §5.6 backstop | `distribution`, `source_column` and `effect_class` each edited by SQL after definition → the writer's message from `add_phenotype()`; nothing drawn or written |
| integrity | strata on a named block, and variance rows deleted → errors before any draw, RNG untouched |
| D3 live | the first call locks the block; `remove_rows()` on the draws clears it; the redeclared block draws marginally |

`test-correlated_draws.R` gains `.blupf90_residual_cov()`: block-diagonal
assembly in the caller's trait order, the two error cases.

**Pre-change check** (§7 "tests before the fix"): the new file run against
`ebabe89` (the Phase 5 end state) fails the day-100 conditional draw, the
Defect 3 test, the mixed-pattern and three-phenotype tests, the backstop
test and the strata integrity test, and passes the ordering, permanent-
environment, gamma/uniform and D3 tests — which is the expected split (the
old code already persisted per level and drew marginally in sorted order;
it never conditioned).

## Verification

- `test-add_phenotype_named_effects.R`: 74/74.
- `test-add_phenotype_stages.R`, `test-add_phenotype.R`,
  `test-add_phenotype_residuals.R`, `test-phenotype_composite.R`,
  `test-phenotype_cov_block.R`, `test-correlated_draws.R` (1 added),
  `test-archive_replicate.R`, `test-schema-registries.R` — all green.
- Full suite at the first commit (`70b9909`): 3384 passed, 0 failed,
  1 skipped (3316 at the Phase 5 end state + 65 + 3). After the second
  review pass: 3393 passed, 0 failed, 1 skipped.
- `roxygen2::roxygenise(roclets = "rd")` clean; `NAMESPACE` unchanged.

## Review pass

Checked against §5.4's adapter responsibilities (1 discover, 2 stratum —
n/a, 3 stable order, 4 stored, 5 fixed — n/a, 6 group, 7 resolve, 8 merge,
9 persist in the outer transaction) and the CLAUDE.md rules (no ids in
SQL, no `dbWriteTable()`, sort before every draw, no old code). Fixes
applied before the commit:

- `.ap_named_effect_targets()` was computed in Stage 1 and again in
  Stage 2; it is now stored on the plan (`plan$named_targets`) beside the
  blocks it was used to load, so Stage 2 reads the plan and nothing else.
- The 1 × 1 dispatch read the first in-call term's `distribution` before
  knowing the block was 1 × 1; it now reads it only then, and the
  `normal`/`NA` case falls through to the resolver in one place.
- `get_phenotype_var()` (still `define_effect_random()`'s reader of the
  unconditional diagonal) escaped its literals with `gsub()`; it now uses
  `dbQuoteLiteral()` like every other reader on this path.
- `R/schema.R`'s `phenotype_random_effects` descriptions said "populated on
  first use"; they now state the persistence rule and the conditional
  draw.
- Confirmed nothing references the deleted functions outside historical
  `NEWS.md` entries and the plan's "before Phase 6" narrative.

Second pass (after the commit `70b9909`), against seven adversarial
probes — integer levels read from a non-`ind_meta` source table, every
planned level `NULL`, a sex-expressed phenotype whose pens are a subset of
the block partner's, a categorical phenotype in a block, a gamma 1 × 1 on
a reuse call, two effects with one partially stored, and conditioning
after `close_pop()` / `restore_pop()`:

- **Bug: a phenotype whose planned levels are all `NULL` crashed.**
  `.ap_named_effect_block()` built its planned-level frame with
  `data.frame(level = <length 0>, phenotype_name = t)`, which R refuses
  ("differing number of rows: 0, 1"). The pre-Phase-6 marginal path
  handled this case silently (nothing to draw). Fixed by recycling the
  name to the level count; the case is now a test ("all NULL … alone or
  in a block": no draw, RNG untouched by the effect, contribution `0`,
  and in a block the partner still draws for its own pens while the
  all-`NULL` phenotype's coordinate stays latent).
- `.ap_plan()`'s roxygen now lists `named_targets` / `named_blocks` in
  the returned plan.
- The other six probes passed unchanged; the categorical probe was
  wrong about the package (a `prevalence` phenotype writes category
  indices `1`/`2`, not `0`/`1`), not the adapter.

## Plan bookkeeping

- Header: v3.6, Phases 0–6 shipped; "What changed from v3.5 to v3.6"
  (adapter, Stage-1 block loading, final RNG order, backstop, 1 × 1
  dispatch, `load_phenotype_cov()` retired, docs).
- Defect 3: ✅ closed. §5.2 / D1: the "until Phase 6" parentheticals on
  `load_phenotype_cov()` removed. §5.3: named-effect test note.
- §5.5: heading "shipped in full"; the v3.4 note updated.
- §5.6: dispatch and backstop marked as shipped.
- §5.8: heading ✅; "What happens today" → "What happened before Phase 6";
  the persistence-note obligation marked done.
- §7: "Named effects" ✅ with the item → test mapping.
- §8: Phase 6 ✅ with the shipped summary; original scope retained.
- `NEWS.md` under `0.71.0 (in development)` → Changed; CLAUDE.md updated.
- `DESCRIPTION` bumped to `0.71.0` (the version `NEWS.md` accumulates
  under); Phase 8 keeps the final release bookkeeping.

## Next

Phase 7 — the D7 transaction/RNG boundary test: a forced Stage-2 failure
leaves the database unchanged and `.Random.seed` advanced by exactly the
Stage-2 draws made before the failure; a forced Stage-3 failure rolls back
named draws and records together. Then Phase 8 (docs pass, the D6
mutability item, `DESCRIPTION` bump, benchmark).
