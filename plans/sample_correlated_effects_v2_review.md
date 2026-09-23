# Review of `sample_correlated_effects.md` v2

**Status:** Strong direction, but not yet implementation-ready  
**Reviewed:** `plans/sample_correlated_effects.md` and the current phenotype/effect code

## Executive assessment

The core decision is correct:

- persist realized liability-scale residuals;
- keep named random-effect draws in their existing per-level table;
- lazily sample missing coordinates from a conditional multivariate normal; and
- use one database-independent Gaussian resolver for residual and named-effect
  adapters.

That model handles the motivating cases: phenotypes recorded at different times,
culling between records, partially overlapping subsets, and pen/effect vectors
completed in later calls. The observed/requested/latent distinction and the
decision not to pre-draw lifetime trajectories are especially good.

I would retain the overall architecture. Before implementation, however, the
plan needs several corrections. The largest issue is not the conditional-MVN
formula; it is defining exactly which entities and covariance definition are in
force before any random number is consumed.

## Implementation blockers

### 1. Final record eligibility must be resolved before stochastic sampling

The proposed control flow says to resolve all stochastic state before the write
loop, but the current `add_phenotype()` code determines the final population
later:

- `compute_covariate_contribution()` can remove individuals through
  `null_class_action = "skip"`;
- formula-TBV evaluation can remove individuals with missing components;
- composite-TBV assembly can also remove individuals; and
- derived and direct-value paths bypass residual generation entirely.

If residuals and named effects are resolved from the preliminary
`subset_by_pheno`, the implementation can consume RNG for records that will not
exist. A residual cannot even be persisted under Option A unless its
`ind_phenotype` row exists. This also makes seeded output depend on later
filtering details.

Replace the proposed flow with three explicit stages:

1. **Plan records without RNG or writes.** Resolve final eligible IDs, TBVs,
   covariate inputs, phenotype paths, condition levels, named-effect levels, and
   the next `pheno_number` for every phenotype.
2. **Resolve stochastic state in memory.** Discover blocks, load stored
   coordinates, incorporate supplied residuals, and draw only coordinates that
   will be used by planned records.
3. **Commit once.** Insert new named-effect coordinates and phenotype records,
   including residual metadata, in one transaction.

This is a real restructuring of `add_phenotype()`, not just replacing sections
7.5 and 8.5. The rough estimate of approximately 250 net lines is probably too
optimistic unless record planning is first extracted into a clean internal
representation.

### 2. Covariance-block identity and replacement semantics need an exact rule

The plan can recover block membership from the existing row-oriented schema,
but only if the graph is defined by **the existence of a stored pair row**, not
by `cov_value != 0`. This distinction is essential: an explicitly stored zero
off-diagonal still declares two coordinates to be in the same block.

The definition functions also need exact replacement behavior. For a new matrix
over coordinate set `N`:

1. Find every existing block for the same effect and covariance stratum that
   intersects `N`.
2. Let `U` be `N` plus every coordinate in those touched blocks.
3. Require `N` to equal `U`; otherwise error and name the omitted coordinates
   and pairs.
4. Validate the complete new matrix.
5. Delete every row for the touched blocks and insert the new `N x N` rows in
   one transaction.

Without that rule, replacing `{A, B}` inside an existing `{A, B, C}` block can
leave a syntactically complete but globally unvalidated mixture of old and new
entries.

This invariant must apply to **all covariance writers**, not only
`define_residual_cov()` and `define_effect_cov_matrix()`.
`define_phenotype(residual_var = ...)` currently writes a scalar residual
diagonal, and `define_effect_random(variance = ...)` writes a scalar named-effect
diagonal. Either scalar writer can mutate one cell of an existing multivariate
block and invalidate its PSD property. Route all four entry points through one
validated block-update operation, or reject scalar updates when the coordinate
already belongs to a multivariate block.

For conditional residual covariance, define a covariance stratum as at least:

```text
(effect_name, condition_table, condition_column, condition_level)
```

The plan should decide whether every condition level must have exactly the same
coordinate set as the unconditional block. Requiring identical coordinate sets
in v1 is the clearest and safest contract. Cell-by-cell fallback between a
partial conditional matrix and an unconditional matrix should not be allowed.

### 3. `force = TRUE` makes future conditional draws incoherent

D3 correctly rejects covariance redefinition after realization, but its
`force = TRUE` escape hatch defeats the same guarantee.

Suppose `A` was drawn under `R1`, the block is forcibly replaced by `R2`, and
`B` is later requested. The proposed sampler will condition the old `A` value
using `R2`. This is not merely an old phenotype that is “distributionally
inconsistent”; the newly generated pair has no declared joint distribution.
Named effects have no covariance-version field, and residuals store a condition
level rather than a covariance-definition version, so the sampler cannot detect
or isolate the two eras.

Recommended v1 rule: **remove `force = TRUE` once any block coordinate has been
realized.** Model-building remains editable until the first draw. After that,
users must start a new population or explicitly remove/regenerate all dependent
records and draws through a future purpose-built operation.

If force is considered indispensable, add a covariance-definition ID to every
realization and specify that cross-version conditioning is forbidden. That is
substantially more machinery and still cannot repair phenotype values already
computed from the old draws.

### 4. Singular Gaussian conditioning needs a support-consistency check

The numerical section correctly allows positive-semidefinite covariance and a
pseudoinverse. One condition is missing. For singular `R_oo`, an observed vector
must lie in the support of its Gaussian distribution:

```text
(I - R_oo R_oo^+) e_o = 0
```

within scale-aware tolerance. This matters for perfect correlation,
zero-variance coordinates, user-supplied residuals, and any externally modified
database. For example, a nonzero supplied residual for a zero-variance
coordinate has probability zero and does not define a valid conditional
distribution.

The resolver should error when the observed vector lies materially outside the
support. Applying a pseudoinverse without this check returns a number, but not a
valid conditional draw.

Use a scale-aware eigenvalue tolerance such as a function of matrix dimension,
largest eigenvalue, and machine precision. A single absolute tolerance passed
unchanged across variances of very different magnitude will be unreliable.

### 5. Supplied residuals need an explicit “fixed current coordinate” state

The three-state model covers stored, newly sampled, and latent coordinates, but
`user_residual` introduces a fourth case: a current coordinate whose value is
fixed by the caller rather than sampled.

If `A` and `B` are requested together, the caller supplies `A`'s residual, and
`B` is model-generated, then `B` must be drawn conditionally on the supplied
`A` value in the same call. If both are supplied, neither is drawn, but both are
stored and can condition a later `C` coordinate.

Define the state model as:

- **stored:** prior realized coordinate;
- **fixed:** current caller-supplied coordinate;
- **sample:** current missing coordinate to draw;
- **latent:** neither requested nor realized.

The conditional observed set is `stored + fixed`; only `sample` consumes RNG.
The plan should also state whether a named list may supply residuals for only a
subset of current phenotypes. Supporting that naturally follows from the model
above and is more useful than the current all-or-nothing branch.

### 6. Heterogeneous-residual policy is underspecified for a block

`condition_change_action` is proposed as a per-phenotype setting, while the
decision is made for a multivariate block. The plan does not say what happens if
coordinates in the same block have different settings.

Choose one rule:

- require all phenotypes in a residual block to use the same action; or
- make the strictest action win and document that behavior.

Requiring agreement is easier to reason about. Because covariance can currently
be declared before phenotype metadata, enforce agreement when sampling as well
as when all metadata are available at definition time.

For blocks larger than two coordinates, `independent` also needs precise
semantics. The clean rule is: discard only stored coordinates whose actual
covariance stratum differs from the current record, but continue conditioning on
compatible stored coordinates. Then draw the requested coordinates from the
current stratum's conditional distribution. Warn by entity and identify the
dropped coordinates.

The selected stratum must be recorded, not merely the raw value in the condition
column. If an unmatched or missing condition uses the unconditional fallback,
store `NULL` as the actual selected stratum. If neither a matching conditional
matrix nor a complete unconditional fallback exists, error. The current behavior
that can set residuals to zero is not an acceptable fallback.

## Important corrections to the document

### Use the actual table name

The current table is `phenotype_random_effects`, not
`trait_random_effects`. The latter appears throughout the plan in the schema,
pseudocode, test discussion, and implementation phases. All of those references
should be corrected before the plan is used as an implementation checklist.

### Phase 0 is stale

The current `ensure_trait_tables()` is already pure idempotent DDL and its
documentation already says that old databases are not migrated. The migration
blocks listed in Phase 0 are no longer present at the cited lines. Remove Phase
0 from this plan and renumber the remaining phases.

The cited fallback in `add_phenotype()` does still exist, so its removal belongs
in residual integration rather than in a broad migration-cleanup phase.

### “No block” versus a one-dimensional block

A phenotype still needs its marginal variance. In the proposed representation,
a scalar variance is naturally a complete one-dimensional block. The phrase “a
phenotype in no block keeps the existing marginal `rnorm()` path” is ambiguous:
no covariance rows at all should still be an error, while a valid 1x1 block uses
the marginal Gaussian path. State this explicitly.

### Pen identity is persistent unless occasion is part of the level

The named-effect example correctly says a stored pen draw is reused forever.
That is a strong biological assumption that should be highlighted near the
user-facing example.

If `pen_id = "P1"` is reused over several batches or seasons, the model says the
same realized pen effect persists across all of them. A temporary pen or
pen-by-batch effect requires a level such as `interaction(pen_id, batch_id)` or a
dedicated `pen_occasion_id`. This is separate from correlating the ADG and BF
coordinates of the same pen-level vector.

This distinction is especially important because “effects at different times”
can mean either completing one persistent multivariate pen vector or modeling a
new time-indexed pen realization. This plan supports the first, not the second.

## Revised internal contract

The shared resolver can remain small if database concerns stay in adapters:

```r
resolve_correlated_draws <- function(covariance,
                                     sample_coordinates,
                                     entity_keys,
                                     observed,
                                     tolerance) {
  # observed contains both prior stored and current fixed coordinates
  # return values for sample_coordinates only, in normalized stable order
}
```

Adapter responsibilities should be explicit:

1. Discover the complete block by stored pair-row existence.
2. Select the exact covariance stratum for every planned entity.
3. Produce stable entity and coordinate order.
4. Load prior values and validate their covariance-definition compatibility.
5. Add fixed current values such as `user_residual`.
6. Group entities by `(stratum, observed coordinates, sample coordinates)`.
7. Call the resolver once per pattern.
8. Merge fixed and sampled values into planned phenotype rows.
9. Persist only inside the outer transaction.

The resolver responsibilities should be limited to:

- validating dimensions and finite values;
- checking support consistency for singular observed covariance;
- computing conditional means and covariances;
- projecting only numerically tiny negative eigenvalues to zero;
- sampling in stable order; and
- consuming no RNG for empty or fully resolved requests.

## Transaction and reproducibility boundary

The proposed database transaction is necessary, but it does not make the whole
operation atomic: an R error after drawing advances `.Random.seed` even if the
database rolls back.

Choose and document one of these contracts:

- database state is atomic, but RNG state advances on failure; or
- capture the pre-resolution RNG state and restore it on any error before a
  successful commit.

The second gives cleaner retry behavior, especially when `seed` is supplied.
Whichever rule is chosen, test it directly. Do not describe a rolled-back write
as fully reproducible unless RNG rollback is included.

## Additional tests to add

The existing test matrix is strong. Add these cases:

1. A scalar variance writer attempts to alter one diagonal of an existing
   multivariate block.
2. A replacement matrix names only a strict subset of a touched block.
3. An all-zero off-diagonal remains part of the declared block because the pair
   row exists.
4. Conditional covariance levels have different coordinate sets and are
   rejected.
5. A skipped individual consumes no residual RNG and gets no orphan stochastic
   state.
6. A current supplied `A` residual conditions a generated `B` residual in the
   same call.
7. A supplied value violates the support of a singular covariance matrix and
   errors.
8. Different `condition_change_action` values within one block are rejected.
9. In a three-coordinate block, mismatched stored coordinates are dropped under
   `independent` while compatible coordinates still condition the draw.
10. Missing condition level with no unconditional fallback errors rather than
    producing a zero residual.
11. Redefinition after realization has no override in v1.
12. A forced write failure verifies both database rollback and the documented
    RNG-state behavior.
13. Reusing `P1` reuses the same persistent pen draw, while `P1:batch2` is a new
    level.

## Recommended implementation sequence

1. **Centralize covariance definitions.** Implement block discovery by pair-row
   existence, strict replacement semantics, PSD validation, and the realization
   lock for every matrix and scalar writer.
2. **Stabilize schema and names.** Add residual columns and the condition action,
   update schema registries and archive behavior, and use
   `phenotype_random_effects` consistently.
3. **Build the pure resolver.** Include singular-support validation and tests
   independent of DuckDB.
4. **Extract record planning.** Determine final eligible entities,
   `pheno_number`, conditions, and effect levels without RNG or writes.
5. **Integrate residuals.** Treat supplied residuals as fixed coordinates and
   model-generated residuals as sampled coordinates.
6. **Integrate named effects.** Use the same resolver with persistent per-level
   entity identity and compatible source validation.
7. **Add the transaction/RNG boundary.** Persist planned rows atomically and
   implement the chosen RNG-on-error contract.
8. **Document scope.** Emphasize ordinal residual pairing, liability scale,
   persistent versus occasion-specific pen identity, and the exclusion of
   longitudinal covariance models.

## Bottom line

The plan has the right statistical center: lazy conditional Gaussian completion
is the correct way to preserve requested covariance when effects are realized at
different pipeline stages. I would approve that direction.

I would not start integration from the current phase list yet. First revise the
definition/update invariants, remove post-realization covariance forcing, plan
final records before drawing, and make supplied/singular cases explicit. Those
changes turn a good mathematical proposal into a database design that cannot
quietly create incompatible stochastic histories.
