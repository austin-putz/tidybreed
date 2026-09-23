# The three stages of `add_phenotype()`

[`add_phenotype()`](https://austin-putz.github.io/tidybreed/reference/add_phenotype.md)
runs in three stages with a strict boundary between them:

- **Stage 1 — PLAN**
  ([`.ap_plan()`](https://austin-putz.github.io/tidybreed/reference/dot-ap_plan.md)):
  decide the final record list for every phenotype. Sex expression, the
  repeatable guard, fixed-effect contributions and
  `null_class_action = "skip"`, formula/composite TBV evaluation and
  `missing_component_action`, path classification, `pheno_number`
  assignment, the residual condition value of every record, and the
  random-effect level every record touches. **No random number is drawn
  and nothing is written** (the one prerequisite write is
  [`add_tbv()`](https://austin-putz.github.io/tidybreed/reference/add_tbv.md),
  which materializes the TBVs the plan reads and is RNG-neutral).

- **Stage 2 — RESOLVE**
  ([`.ap_resolve()`](https://austin-putz.github.io/tidybreed/reference/dot-ap_resolve.md)):
  every random draw of the call, in a fixed order, and the liability /
  type conversion — all in memory. Nothing is written. Two adapters over
  [`find_covariance_blocks()`](https://austin-putz.github.io/tidybreed/reference/find_covariance_blocks.md)
  and
  [`resolve_correlated_draws()`](https://austin-putz.github.io/tidybreed/reference/resolve_correlated_draws.md)
  draw everything: the named-effect adapter
  ([`.ap_resolve_named_effects()`](https://austin-putz.github.io/tidybreed/reference/dot-ap_resolve_named_effects.md),
  entity `(effect_name, level)`, stored in `phenotype_random_effects`,
  one draw per level reused forever), then the residual adapter
  ([`.ap_resolve_residuals()`](https://austin-putz.github.io/tidybreed/reference/dot-ap_resolve_residuals.md),
  entity `(id_ind, pheno_number)`, stored in `ind_phenotype`, one draw
  per record). In both, an entity draws its planned coordinates from the
  block's Gaussian conditional on the coordinates it has already
  realized — stored on disk from an earlier call, or (residuals only)
  fixed by `user_residual`.

- **Stage 3 — COMMIT**
  ([`.ap_commit()`](https://austin-putz.github.io/tidybreed/reference/dot-ap_commit.md)):
  one transaction that inserts the new `phenotype_random_effects` rows
  and the `ind_phenotype` records via `duckdb_register()` + `INSERT`. No
  random number is drawn, so the stream a call consumes is a function of
  the model and the plan only.

The consequence is that a record that will not exist — an individual
skipped by `null_class_action`, excluded by a missing component, or
refused by the repeatable guard — never consumes RNG and never leaves
stochastic state behind. Planned records are ordered by `id_ind` within
a phenotype; that order is what `user_values` and `user_residual` match
positionally, and it does not depend on physical row order in the
database.

**On failure** (`plans/sample_correlated_effects.md` D7). The database
is atomic and the RNG is not. An error anywhere in the call — a Stage-1
rejection, a Stage-2 error after some draws (a missing variance, a
residual stratum change under `condition_change_action = "error"`), or a
failed Stage-3 write — leaves `ind_phenotype` and
`phenotype_random_effects` exactly as they were: Stage 3 is the only
writer and it rolls back as a whole, so a block resolved before the
failing one is never written on its own, and a column
[`prepare_extra_cols()`](https://austin-putz.github.io/tidybreed/reference/prepare_extra_cols.md)
added by `ALTER TABLE` earlier in the same transaction goes with it.
`.Random.seed` is left advanced by exactly the draws made before the
error; nothing in the package restores it, so a retry draws different
values unless the caller re-seeds. (The
[`add_tbv()`](https://austin-putz.github.io/tidybreed/reference/add_tbv.md)
upsert of Stage 1 is the one write that remains; it does not depend on
the RNG and the retry rewrites it.)

See `plans/sample_correlated_effects.md` §5.5 and D7.
