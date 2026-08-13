# Codex review: sampling correlated stochastic effects across simulated time

**Status:** Design review and implementation recommendations  
**Plans reviewed:**

- `plans/sample_correlated_effects.md`
- `plans/sample_correlated_effects_gemini.md`

## Executive summary

Both plans identify the correct underlying defects and converge on the correct
general solution: persist realized stochastic deviations and lazily sample
missing coordinates from their conditional multivariate normal distribution.

The recommended direction is **Option A**:

- Store liability-scale residuals on `ind_phenotype` as `residual_value`.
- Continue storing named random-effect draws in `trait_random_effects`.
- Share one internal conditional-Gaussian sampling engine across both storage
  paths.

Before implementation, the design needs to be strengthened in several areas:

1. A one-phenotype call must discover correlated phenotypes outside the current
   request.
2. The implementation must distinguish stored, currently requested, and still
   latent coordinates.
3. Conditional covariance calculations must support positive-semidefinite and
   numerically near-singular matrices.
4. Changing a heterogeneous covariance condition between records must be an
   error unless cross-condition covariance is explicitly defined.
5. Repeated-record pairing semantics must be explicit.
6. User-supplied residuals should be persisted, while direct phenotype values
   and derived phenotypes should not invent latent residuals.
7. Stochastic draws and phenotype records should be written atomically.
8. Missing or malformed declared covariance matrices should fail loudly rather
   than silently reverting to independent draws.

This proposal should be described as support for **fixed multivariate Gaussian
stochastic blocks sampled across pipeline stages**. It is not, by itself, a
general longitudinal, random-regression, survival, or multivariate non-Gaussian
phenotype framework. The storage and resolver boundaries recommended below are
chosen so those models can be added without redefining what a realized residual
means.

The longer original plan is materially stronger than the Gemini proposal. It
addresses repeated records, heterogeneous covariance, incomplete matrices,
diagnostics, and non-Gaussian effects. The Gemini plan is a useful summary, but
it omits several decisions required for a safe implementation.

## Points of agreement

### Residual covariance is currently lost across calls

`sample_residuals()` returns an in-memory matrix, and the realized values are
discarded after phenotype construction. A later `add_phenotype()` call therefore
has no value on which to condition.

This is a storage problem, not merely a missing sampling branch.

### The identical-ID restriction is incorrect

The current joint residual path requires correlated phenotypes to be requested
in the same call for identical sets of individuals. This excludes ordinary
simulation workflows involving culling, selection, sex limitation, or records
taken at different stages.

Different subsets should be treated as different per-individual observation
patterns, not as a reason to discard the covariance.

### Partial named-effect redraws do not preserve covariance

For a named effect such as pen, herd, or litter, one phenotype's level may already
have a stored draw while another phenotype's coordinate is still missing. Drawing
a new unconditional joint vector and discarding the already-stored coordinate
does not correlate the newly stored coordinate with the earlier realization.

The missing coordinate must instead be drawn conditionally on the stored value.

### Conditional MVN sampling is the correct mechanism

For new coordinates `n` and already observed coordinates `o`:

```text
e_n | e_o ~ N(
  R_no R_oo^-1 e_o,
  R_nn - R_no R_oo^-1 R_on
)
```

This supports joint calls, sequential calls, partial overlap, culling, and partial
named-effect vectors through the same mathematical operation.

### Lazy sampling is preferable to lifetime pre-drawing

Eager pre-drawing requires knowing all future individuals, records, and relevant
phenotypes in advance. It also creates deviations for records that might never
exist.

Lazy conditioning draws only the values that are needed and naturally handles
new individuals and branching simulation pipelines.

### Residuals and named effects should not share a physical table

A residual is a per-record quantity. A named random effect is a per-level
quantity. A unified table would obscure that distinction, weaken `id_ind` typing,
and require a record-number field that is meaningless for most named effects.

The mathematical resolver should be shared; the schemas should remain separate.

### Correlated non-Gaussian effects require explicit rejection

The marginal named-effect path supports normal, gamma, and uniform distributions,
while the correlated path currently always uses an MVN. Silently changing a
gamma or uniform effect into a normal one is incorrect.

For the proposed implementation, every member of a correlated named-effect block
must use the normal distribution. Other distributions should produce a clear
error until an explicit copula or other multivariate distribution model exists.

### Sequential and joint draws need not be numerically identical

Different call sequences consume random numbers differently. Tests should require:

- exact reproducibility for the same seed and the same call sequence; and
- distributional agreement between joint and sequential sampling strategies.

They should not require identical realized values between different call
sequences.

## Required design improvements

### 1. Discover the complete covariance block

This is the most important missing implementation detail.

When a user calls:

```r
add_phenotype("off_test_wt")
```

the requested phenotype vector contains only `off_test_wt`. The sampler must
still discover that `on_test_wt` belongs to the same residual covariance block.
Restricting covariance lookup to phenotypes named in the current call would
preserve the existing sequential-call bug.

Introduce an internal operation such as:

```r
find_covariance_block(pop, effect_name, target_phenotype)
```

It should return the complete connected component containing the target
phenotype. For example, if A is correlated with B and B is correlated with C,
the block for A is A/B/C even if A/C has zero covariance.

The same discovery rule should be used for residual and named-effect blocks.

### 2. Model coordinate state explicitly

For each individual or named-effect level, distinguish:

- **observed:** coordinates already stored from earlier calls;
- **requested:** missing coordinates being generated in the current call; and
- **latent:** other block coordinates that have not been requested and should
  remain undrawn.

Only requested coordinates should be sampled and stored. Observed coordinates
condition the distribution. Latent coordinates should not be eagerly realized.

This state model handles:

- A and B requested together for the same individuals;
- A requested for more individuals than B;
- B requested later for a culled subset;
- B requested for a mixture of individuals with and without stored A residuals;
- named-effect levels for which only some phenotype coordinates exist; and
- covariance blocks containing more than two phenotypes.

The implementation should group entities by identical `(observed, requested)`
patterns so that each conditional covariance is computed once per group.

### 3. Make the Gaussian conditioning numerically robust

The draft helper's direct use of `solve(Roo)` is insufficient. Valid covariance
matrices can be positive-semidefinite rather than strictly positive-definite,
including blocks with perfect correlation or zero-variance coordinates.

The covariance definition path should validate:

- matching and unique row/column names;
- symmetry within a documented tolerance;
- finite entries;
- nonnegative diagonal values; and
- positive semidefiniteness within a documented tolerance.

The conditional sampler should:

- use a Cholesky solve when `Roo` is positive-definite;
- use a tolerance-aware eigen or pseudoinverse path when it is only
  positive-semidefinite;
- symmetrize the calculated conditional covariance;
- clamp only tiny negative eigenvalues caused by floating-point error; and
- reject materially invalid conditional covariance matrices.

It must also explicitly handle:

- zero entities;
- one entity;
- one requested coordinate;
- zero observed coordinates; and
- zero conditional variance.

The plan's helper sketch also references an undefined `n` in its unconditional
branch. The production helper should take entity count explicitly or derive it
unambiguously from its inputs.

### 4. Reject condition-level changes without a cross-condition model

The original plan proposes using the new record's condition-level covariance
matrix and warning if an earlier residual was sampled under a different level.
This is not mathematically coherent.

If phenotype A was drawn under `R_a` and phenotype B is drawn under `R_b`, neither
matrix specifies `Cov(A under a, B under b)`. Applying `R_b` retroactively does
not create a valid joint model.

Recommended v1 behavior:

- Store the covariance condition level used for each residual realization,
  potentially as `residual_condition_level`.
- Condition only when prior and current coordinates resolve to the same level.
- Error clearly when levels differ.
- Treat cross-condition covariance as a future feature requiring explicit model
  metadata.

Immutable conditions such as sex will naturally pass this check.

### 5. Define repeated-record pairing semantics

`pheno_number` provides storage identity, but it does not automatically provide a
biological time or occasion identity. Pairing A record 2 with B record 2 is only
meaningful when the ordinal records represent corresponding occasions.

Recommended narrow v1 contract:

- Residual covariance applies across distinct phenotype names.
- Prior residual lookup uses the same `pheno_number` across those phenotypes.
- Repeated records of the same phenotype have independent residuals.
- Persistent within-animal covariance should be represented with a permanent-
  environment named effect using `source_column = "id_ind"`.

Documentation should state that `pheno_number` supplies ordinal pairing, not
simulated-time matching. Irregular longitudinal data will eventually require an
explicit occasion or time model.

### 6. Persist `user_residual`, but not unexplained phenotype values

A value supplied through `user_residual` is an explicit liability-scale residual.
It should be stored in `residual_value` so later correlated coordinates can
condition on it.

By contrast:

- `user_values` provides a final phenotype value without identifying its residual
  decomposition; and
- `derived_formula` computes a phenotype without sampling a residual.

Those paths should leave `residual_value` as `NULL` and should not constrain later
stochastic draws.

The documentation should make this distinction explicit.

### 7. Make stochastic resolution and record writes atomic

The new workflow can involve:

- reading prior deviations;
- drawing new residuals;
- drawing new named-effect coordinates;
- adding schema columns; and
- writing phenotype records.

A failure after some draws are stored but before their phenotype records are
written can leave simulation state inconsistent. Resolve all required draws in
memory first, then persist the new draws and phenotype records in a database
transaction.

At minimum, the insertion of correlated draws and the phenotype records that use
them must succeed or fail together.

The transaction also reduces the chance that two callers independently resolve
the same missing coordinate before either writes it.

### 8. Distinguish absent, incomplete, and explicitly independent covariance

`load_phenotype_cov()` currently returns `NULL` when no rows exist and when any
requested matrix cell is missing. Those cases have different meanings:

- no covariance block exists: use the ordinary marginal sampler;
- a covariance block was declared but is incomplete: error;
- a complete block contains an explicit zero covariance: preserve the declared
  independence within that block.

Covariance completeness and validity should primarily be checked by
`define_residual_cov()` and `define_effect_cov_matrix()`, rather than being
discovered deep inside phenotype generation.

Sampling code should never silently interpret a malformed declared block as an
instruction to draw independently.

## Recommended schema

### Residual storage

Add the following column to the base `ind_phenotype` schema:

```sql
residual_value DOUBLE
```

Adding it to the base schema is preferable to altering the table on first use:

- ordinary generated phenotypes always have residuals;
- queries and tests see a stable schema;
- migrations are explicit;
- concurrent connections avoid conditional DDL; and
- users can reliably inspect the realized simulation components.

If heterogeneous residual covariance is retained, also add:

```sql
residual_condition_level VARCHAR
```

Populate `residual_value` for:

- model-generated residuals;
- unconditional residuals; and
- `user_residual` values.

Leave it `NULL` for:

- direct `user_values`; and
- derived-formula phenotype records.

### Named random-effect storage

Keep the existing logical identity:

```text
(phenotype_name, effect_name, level)
```

and continue storing realized values in `trait_random_effects.draw_value`.

Before treating coordinates as members of the same named-effect vector, validate
that the participating phenotype definitions use compatible source semantics.
At minimum, source table and source column must describe the same kind of level.
A pen identifier must not be conditionally paired with a herd identifier merely
because both effects share the same `effect_name` string.

## Recommended internal architecture

Use separate storage adapters around one mathematical resolver.

```text
Residual adapter
  entity identity: id_ind + pheno_number
  coordinate:      phenotype_name
  storage:         ind_phenotype.residual_value

Named-effect adapter
  entity identity: effect_name + level
  coordinate:      phenotype_name
  storage:         trait_random_effects.draw_value
```

The shared resolver can have a conceptual interface such as:

```r
resolve_correlated_draws <- function(
  covariance,
  requested_coordinates,
  entity_ids,
  stored_draws,
  tolerance
)
```

It should:

1. Receive a complete, validated covariance block.
2. Determine observed, requested, and latent coordinates for every entity.
3. Group entities by identical observation/request patterns.
4. Calculate conditional coefficients and covariance once per pattern.
5. Draw only missing requested coordinates.
6. Return a normalized entity-by-coordinate result without writing to the
   database.
7. Allow the caller to persist all changes atomically.

This shares the difficult logic while preserving the genuine schema distinction
between per-record residuals and per-level named effects.

## Same-call subset behavior

The implementation must also correct differing subsets within one
`add_phenotype()` call.

For example, if A is requested for individuals 1–100 and B for individuals
51–100:

- individuals 1–50 receive a marginal A draw;
- individuals 51–100 receive a joint A/B draw; and
- no B draw is created for individuals 1–50.

If A for individuals 1–50 is written before the overlapping group is processed,
B must be conditionally drawn from those exact A realizations. Prefer resolving
the full current-call entity state before writing, so results do not depend on
the order in which phenotype loops happen to execute.

## Scope boundary: repeated records are not random regression

The v1 proposal supports correlated coordinates of a fixed Gaussian vector that
are realized at different points in the simulation pipeline. For example,
on-test and off-test weight can be two distinct phenotypes with a fixed residual
covariance, even when the second is recorded after culling.

That is different from a longitudinal random-regression model. A typical random-
regression phenotype has a structure such as:

```text
y_it = x_it beta + z(t_it)' a_i + z(t_it)' p_i + e_it
```

where:

- `t_it` is the observation age, date, stage, or other continuous trajectory
  coordinate;
- `z(t_it)` is a basis evaluated at that coordinate;
- `a_i` is an animal-specific vector of genetic regression coefficients;
- `p_i` is an animal-specific vector of permanent-environment coefficients; and
- `e_it` is a record-specific residual whose variance or covariance may also
  depend on time.

Matching records by equal `pheno_number` does not express this model. It provides
only ordinal pairing between distinct phenotypes under the deliberately narrow
v1 contract.

The existing `weight_type` and `poly_order` columns provide potentially useful
metadata hooks, but they do not constitute a complete random-regression model.
In particular, the package currently lacks a durable observation coordinate,
basis-domain metadata, persisted coefficient vectors, and a covariance provider
that evaluates covariance as a function of time.

### Requirements for future random-regression support

A complete implementation will need explicit decisions for:

1. **Observation coordinate.** Every longitudinal record needs an age, date,
   stage, parity, test day, or user-defined numeric coordinate. This must be
   separate from `pheno_number`, which is only an ordinal record counter.
2. **Occasion identity.** When multiple phenotypes are measured at the same
   biological occasion, they need an explicit shared occasion key rather than an
   assumption that equal record numbers mean equal occasions.
3. **Basis definition.** Store basis family, polynomial order, centering/scaling
   range, knot placement where applicable, and extrapolation policy. Legendre
   coefficients are uninterpretable without the domain used to scale time.
4. **Coefficient covariance.** Genetic and permanent-environment regression
   coefficients require covariance matrices over basis coefficients and,
   potentially, across multiple phenotypes.
5. **Coefficient persistence.** Coefficient vectors should be sampled once per
   animal or effect level and reused when new records are generated. They should
   not be represented as record residuals.
6. **Time-dependent residuals.** Residual variance may be heterogeneous by age,
   and residual covariance may be a function of two observation times or their
   lag.
7. **Irregular schedules.** Individuals may be measured at different times and
   have different numbers of records. The model cannot depend on balanced record
   numbers.
8. **Boundary behavior.** Define whether observations outside the configured
   basis range error, warn, clamp, or extrapolate.
9. **Model versioning.** Once coefficients or residuals have been sampled, basis
   definitions and covariance parameters cannot be changed silently.
10. **Evaluation integration.** Generated records and exported evaluation models
    must use consistent basis construction and ordering.

### Recommended future record identity

Do not overload `pheno_number` with time semantics. A future additive extension
could introduce fields conceptually equivalent to:

```sql
occasion_id       VARCHAR
observation_time  DOUBLE
time_scale        VARCHAR
```

The exact schema should be designed with the random-regression feature rather
than added speculatively in v1. However, v1 documentation and code should avoid
assuming that `pheno_number` will eventually become the time coordinate.

## Future phenotype and dependence models

The fixed Gaussian resolver is a useful foundation, but it does not make every
phenotype family a multivariate Gaussian model. Future support should distinguish
the following cases explicitly.

### Threshold, binary, and ordinal phenotypes

For threshold traits, storing and conditioning on the liability-scale residual is
coherent. Category conversion occurs after construction of the latent liability.

Future work should nevertheless define:

- whether cross-trait covariance is among liabilities or observed categories;
- how fixed prevalence thresholds interact with heterogeneous residual variance;
- whether thresholds vary by condition or time; and
- how user-supplied observed categories can participate when their latent
  residual is unidentified.

Directly supplied categories must not be reverse-engineered into arbitrary latent
residuals.

### Count phenotypes

The current count path rounds and clips a Gaussian liability. Conditional Gaussian
residual sampling remains mechanically consistent with that existing model, but
it does not create a Poisson, negative-binomial, or generalized linear mixed
model.

Future true count-model support may require:

- a log or other link function;
- Poisson or negative-binomial conditional sampling;
- dispersion and exposure/offset metadata;
- zero-inflation or hurdle components; and
- a copula or shared latent-effect model for cross-trait dependence.

The v1 documentation should call the current behavior a transformed Gaussian
liability model rather than implying general correlated count support.

### Zero-inflated and hurdle phenotypes

These contain at least two stochastic processes: occurrence and positive
magnitude. A single `residual_value` cannot identify both components.

Future implementations will need component-specific latent draws and covariance
definitions, potentially allowing dependence between occurrence and magnitude
components across phenotypes.

### Survival and time-to-event phenotypes

Survival traits require event-history semantics rather than an ordinary residual
attached to one phenotype row. Important dimensions include:

- event time and entry time;
- right, left, or interval censoring;
- time-varying covariates;
- recurrent events;
- competing risks;
- frailty effects; and
- informative censoring or dropout.

The current proposal neither implements nor blocks such a model, but survival
draws should not be forced into the fixed residual-coordinate abstraction.

### Longitudinal state-space and autoregressive residuals

Some repeated phenotypes require residual dependence such as AR(1), continuous-
time exponential decay, Gaussian processes, or state-space innovations.

These differ from a fixed A/B residual covariance because covariance depends on
the pair of observation times. Efficient sequential sampling may exploit model-
specific state rather than repeatedly conditioning a growing dense MVN.

The future covariance interface should therefore permit specialized samplers and
state, not require every model to materialize a dense covariance matrix.

### Time-varying group and environmental effects

An animal can move between pens, herds, farms, management groups, or spatial
locations. A named effect keyed only by `level = "pen_1"` represents a permanent
pen realization, not a pen-period or exposure episode.

Future models may need composite identities such as:

```text
(pen_id, period_id)
(farm_id, year, season)
(location_id, observation_time)
```

The effect source must resolve the level applicable at the record's observation
occasion, rather than reading only the animal's current `ind_meta` value.

### Maternal, social, spatial, and dyadic effects

Not every stochastic entity is the focal individual or a scalar group level:

- maternal effects may vary by dam, parity, or litter;
- social effects may depend on contemporaneous group membership and group size;
- spatial effects may be correlated across locations or neighborhoods; and
- dyadic effects may belong to ordered or unordered pairs of individuals.

The shared resolver should therefore use an opaque, adapter-provided entity key.
It must not assume internally that every entity key is an `id_ind` or that every
named-effect level is time-invariant.

### Reaction norms and genotype-by-environment models

Reaction norms are closely related to random regression, with an environmental
gradient replacing time. They require:

- an explicit environmental coordinate;
- intercept/slope or richer basis coefficients;
- covariance among coefficients and traits;
- stable scaling of the environmental gradient; and
- a policy for environments outside the configured range.

The same future coefficient-storage and basis-definition system should serve both
random regression and reaction norms.

### Non-Gaussian multivariate dependence

Conditional MVN sampling applies only when the stored coordinates form a jointly
Gaussian vector. Future dependence among gamma, uniform, count, binary, or mixed
families requires an explicit model, such as:

- a Gaussian copula with marginal transformations;
- shared latent Gaussian effects;
- a multivariate generalized linear mixed model; or
- a purpose-built joint distribution.

The package should never infer one of these merely because a covariance matrix
was supplied. Until such an API exists, correlated non-normal named effects must
error.

## Additional edge cases and lifecycle rules

### Informative selection, culling, and missingness

Culling based on an already observed phenotype is naturally supported: later
records are generated conditionally for the surviving subset.

Different mechanisms require different models:

- missing completely at random needs no additional stochastic state;
- missingness based on observed history is a pipeline-selection rule;
- dropout based on an unobserved latent value is informative missingness and
  requires a joint dropout model; and
- censoring may preserve partial information about an unobserved liability rather
  than reveal its exact residual.

The v1 resolver conditions only on exact stored deviations, not interval or
selection information.

### Deletion and regeneration

Deleting an earlier phenotype record also deletes its residual under Option A,
but later coordinates may already have been sampled conditional on that value.
Their marginal distribution remains valid, yet the database no longer contains
the provenance needed to reconstruct the sequential conditional history.

A policy is required before advertising regeneration workflows. Possible future
choices include:

- prohibit deletion of a stochastic coordinate with realized dependents;
- track dependency/provenance edges;
- invalidate and regenerate the affected connected component; or
- define deletion as removing an observation while retaining an internal latent
  realization in separate storage.

For v1, document that deleting and resampling correlated phenotype records is not
guaranteed to reproduce a coherent conditional history.

### Covariance changes after realization

Covariance matrices, basis definitions, distribution families, and condition
mappings must not be silently changed after draws from the affected block exist.

The robust future solution is model-definition versioning: every stored draw is
associated with the definition version that generated it. A smaller v1 solution
is to reject covariance-definition updates when affected residuals or named draws
already exist.

### Duplicate and concurrent requests

Two sessions or nested pipeline operations could both observe a missing
coordinate and sample different replacements. Transactions alone may not prevent
this without logical uniqueness and suitable write conflict behavior.

The implementation should enforce the logical identities:

- residual: `(id_ind, phenotype_name, pheno_number)`; and
- named effect: `(phenotype_name, effect_name, level)`.

It should fail or deterministically reuse the winner rather than silently store
duplicates.

### Deterministic parallel sampling

Grouping entities by observation pattern can change RNG consumption when query
or batch order changes. Future parallel execution should use stable ordering or
entity/block-derived random streams so results do not depend on thread count,
database row order, or batch size.

V1 should at least sort blocks, patterns, coordinates, and entity identities
before drawing, and document reproducibility as conditional on the same call
sequence and model state.

### Empty, degenerate, and extreme blocks

Tests and validation should cover:

- no eligible entities after filtering;
- all requested coordinates already stored;
- one phenotype or one entity;
- zero residual variance;
- perfect positive or negative correlation where valid;
- very small or large variances;
- block members with explicit zero covariance;
- disconnected covariance components; and
- missing values in stored draws or entity keys.

## Future-compatible extension points

Two abstractions are more important than reserving speculative columns.

### Observation context

The resolver's adapters should be capable of receiving an observation context,
even if v1 uses only phenotype name, record number, and condition level:

```text
entity identity
phenotype coordinate
record identity
occasion identity
observation time or environmental coordinate
condition level
model-definition version
```

This does not require adding all fields to the v1 schema. It means avoiding helper
interfaces that accept only `(id_ind, phenotype_name)` and would later need to be
replaced wholesale.

### Covariance provider and sampler strategy

The fixed-matrix loader should sit behind a conceptual interface such as:

```r
resolve_stochastic_model(
  effect_name,
  phenotype_coordinates,
  entity_context,
  observation_context
)
```

For v1 it can return a fixed Gaussian covariance block. Future providers may
return:

- a covariance matrix evaluated at observation times;
- a random-regression coefficient model;
- a kernel or sparse precision representation;
- a state-space transition model;
- a copula plus marginal distributions; or
- a specialized survival/event sampler.

The conditional-MVN implementation should therefore be one sampler strategy,
not hard-coded as the permanent definition of all correlated stochastic effects.

## Proposed support roadmap

### Layer 1: fixed multivariate Gaussian blocks

This plan's implementation target:

- residual covariance across distinct phenotypes and pipeline stages;
- partial overlap and culling;
- liability-scale threshold traits;
- Gaussian named random effects; and
- fixed covariance matrices, optionally stratified by an unchanged condition.

### Layer 2: longitudinal Gaussian models

Add:

- observation times and shared occasion identity;
- random-regression and reaction-norm coefficients;
- permanent-environment trajectories;
- heterogeneous residual variance by time;
- covariance functions and AR/state-space residual models; and
- time-indexed group exposures.

### Layer 3: non-Gaussian and event models

Add explicit APIs for:

- generalized count traits;
- zero-inflated and hurdle processes;
- copula or latent-variable mixed-family dependence;
- survival, censoring, competing risks, and recurrent events; and
- informative dropout or observation processes.

Each layer should preserve the distinction between a realized latent stochastic
component and the observed phenotype value derived from it.

## Error and warning policy

Use errors when correctness cannot be guaranteed:

- incomplete declared covariance matrix;
- asymmetric or non-positive-semidefinite covariance matrix;
- non-normal distribution in a correlated named-effect block;
- incompatible source columns/tables within a named-effect block;
- prior and current residual condition levels differ;
- ambiguous duplicate stored draws; or
- incompatible stored coordinate metadata.

Warnings are appropriate for unusual but coherent behavior, not for silently
changing the probability model.

The current silent fallback to independent draws should be removed.

## Verification and testing priorities

### Core residual cases

1. Sequential A then B for identical individuals.
2. A followed by B for a culled subset.
3. A and B requested together for partially overlapping subsets.
4. B requested for a mixed group where only some individuals have stored A.
5. B requested first and A later, verifying order symmetry in distribution.
6. A/B/C covariance blocks with multiple observation patterns.
7. Individuals with no prior coordinate receive the correct marginal draw.
8. Individuals outside the later subset remain untouched.
9. Disconnected covariance components are resolved independently.
10. Explicit zero covariance within a declared block remains valid.

### Repeated records

11. Matching distinct phenotypes by equal `pheno_number`.
12. Repeated records of the same phenotype remain residual-independent.
13. Unequal record counts do not accidentally condition on the wrong record.
14. Irregular record schedules do not get mistaken for random regression.

### Named random effects

15. A-level draw stored first, followed by conditional B for the same levels.
16. Many levels with different partial-coordinate patterns.
17. New levels receive correct unconditional joint or marginal draws.
18. Stored level draws are reused exactly.
19. Incompatible source semantics are rejected.
20. Non-normal correlated named effects are rejected.
21. Time-varying group levels are not silently treated as permanent levels.

### Numerical behavior

22. One-dimensional conditional and unconditional draws.
23. Perfectly or nearly perfectly correlated covariance matrices.
24. Zero conditional variance.
25. Invalid negative-eigenvalue covariance rejection.
26. Tiny floating-point negative eigenvalues handled within tolerance.
27. Empty entity sets and fully resolved requests are RNG-neutral.
28. Extreme but valid variance scales remain numerically stable.

### Alternate input paths

29. `user_residual` is stored and conditions later draws.
30. `user_values` leaves `residual_value` null.
31. Derived phenotypes leave `residual_value` null.
32. Changed covariance condition level is rejected.
33. Threshold traits condition on latent residuals, not observed categories.
34. Count traits are documented and tested as transformed Gaussian liabilities.

### Reproducibility and integrity

35. Same seed and same call sequence produce identical output.
36. Joint and sequential paths reproduce the target covariance statistically.
37. Forced write failure rolls back residuals, named draws, and phenotype rows.
38. Existing databases migrate without losing phenotype records.
39. Stable entity ordering makes results independent of database row order.
40. Duplicate or concurrent requests cannot create duplicate logical draws.
41. Covariance updates are rejected once affected draws exist.
42. Deletion/regeneration behavior follows the documented v1 restriction.

Distributional tests should use enough entities for stable estimates and derive
tolerances from sampling uncertainty rather than arbitrary fixed margins.

## Recommended implementation sequence

### Phase 1: validation and observability

1. Make malformed or incomplete declared covariance matrices fail loudly.
2. Reject non-normal distributions in correlated named-effect blocks.
3. Add tests demonstrating every current silent-fallback defect.

### Phase 2: schema and storage

4. Add `residual_value` to the base `ind_phenotype` schema and migration path.
5. Add `residual_condition_level` if heterogeneous covariance remains supported.
6. Store ordinary and user-supplied residuals with phenotype records.

### Phase 3: shared conditional sampler

7. Implement covariance-block discovery.
8. Implement the robust conditional-Gaussian resolver.
9. Test it independently of database access across all coordinate patterns and
   numerical edge cases.

### Phase 4: residual integration

10. Resolve current-call residuals for all phenotype subsets before the
    per-phenotype write loop.
11. Condition on prior `ind_phenotype.residual_value` records.
12. Write residuals and phenotype records transactionally.

### Phase 5: named-effect integration

13. Replace the current correlated named-effect pre-draw path with the shared
    resolver and a `trait_random_effects` storage adapter.
14. Validate compatible effect-level sources across phenotype coordinates.
15. Make named-effect and phenotype writes atomic.

### Phase 6: documentation and performance

16. Document sequential sampling, repeated-record semantics, and condition-level
    restrictions.
17. Add culling and partial-overlap examples.
18. Benchmark large populations and optimize observation-pattern queries and
    batched writes without changing RNG semantics.

## Final recommendation

Adopt Option A and lazy conditional-MVN sampling, but do not implement the draft
pseudocode literally.

The durable design is:

- complete covariance-block discovery;
- explicit observed/requested/latent coordinate state;
- a numerically robust, database-independent Gaussian resolver;
- residual storage beside phenotype records;
- named-effect storage beside effect levels;
- strict validation instead of silent independence;
- explicit repeated-record and condition-level semantics; and
- transactional persistence.

This provides correct correlated effects across simulation stages without adding
a new user-facing concept, while leaving longitudinal covariance and
non-Gaussian dependence as explicit future extensions rather than accidental
semantics in the current implementation.

It should be presented specifically as **Layer 1: fixed multivariate Gaussian
blocks**. Random regression, covariance functions, true generalized phenotypes,
and event-history models require observation context and specialized stochastic
model providers. The v1 implementation should preserve those extension points,
but should neither reserve speculative biology in the schema nor claim those
future models are already supported.
