# Sampling correlated random effects at different points in simulated time

**Status**: Draft **v2** — for review. v2 folds in the review round (Gemini
proposal in `sample_correlated_effects_gemini.md`, Codex critique in
`sample_correlated_effects_codex_review.md`). Both reviewers independently
confirmed the three defects and converged on Option A + lazy conditional-MVN
sampling, so the storage decision is treated as settled and compressed below.
The substance of v2 is the **eight design corrections** Codex identified, two of
which reverse positions I took in v1. The four open questions raised in review
have since been **decided by the author** and are recorded with their implications
in §6; no open questions remain blocking implementation.

**Scope name (adopted from Codex)**: this is **Layer 1 — fixed multivariate
Gaussian blocks sampled across pipeline stages.** It is not a longitudinal,
random-regression, survival, or non-Gaussian dependence framework. Saying so up
front matters, because the failure mode of this feature is users assuming it
means more than it does.

**Author intent** (unchanged): a user who wrote
`define_residual_cov(c("on_test_wt", "off_test_wt"), R)` should get the covariance
they asked for whether the two phenotypes are recorded in one call or a hundred
simulated days apart, with culling in between — with zero new user-facing
concepts.

---

## 1. The motivating scenario

Pigs on test. On-test weight at entry; off-test weight ~100 days later, with
culling in between, so the off-test set is a **subset** of the on-test set.

```r
pop <- pop |>
  define_phenotype("on_test_wt",  type = "continuous", mean =  30) |>
  define_phenotype("off_test_wt", type = "continuous", mean = 120) |>
  define_residual_cov(
    c("on_test_wt", "off_test_wt"),
    matrix(c(16, 20, 20, 100), 2, 2,
           dimnames = list(c("on_test_wt", "off_test_wt"),
                           c("on_test_wt", "off_test_wt")))
  )

# ── simulated day 0 ───────────────────────────────────────────────
pop <- pop |> get_table("ind_meta") |> filter(gen == 1L) |>
  add_phenotype("on_test_wt")

# ... 100 simulated days, selection, culling ...

# ── simulated day 100 ─────────────────────────────────────────────
pop <- pop |> get_table("ind_meta") |> filter(gen == 1L, alive) |>
  add_phenotype("off_test_wt")
```

**Today this silently produces uncorrelated residuals.** No error, no warning.
The off-diagonal of 20 is never used.

There is currently **no way to express this scenario correctly**, including by
sampling both phenotypes in one call: the joint path requires *byte-identical*
`id_ind` sets ([add_phenotype.R:568-572](../R/add_phenotype.R#L568-L572)), and
culling guarantees they differ.

---

## 2. Confirmed defects

All three were independently confirmed by both reviewers.

### Defect 1 — residual deviations are never persisted

`sample_residuals()` ([phenotype_helpers.R:46-66](../R/phenotype_helpers.R#L46-L66))
returns an in-memory matrix, consumed at
[add_phenotype.R:766-788](../R/add_phenotype.R#L766-L788) to build the liability
and then discarded. A later call has nothing to condition on, so cross-call
correlation is not merely unimplemented — it is *unimplementable* without a
storage change. **This is a storage problem, not a missing sampling branch.**

### Defect 2 — the joint residual MVN is gated on same-call *and* identical ID sets

[add_phenotype.R:566-635](../R/add_phenotype.R#L566-L635) draws a joint
`MVN(0, R)` only when `length(phenos) >= 2`, `all_equal` (byte-identical sorted
`id_ind` vectors), and `R_unconditional` is non-`NULL`. Otherwise it falls
through **silently** to the independent `rnorm()` at
[add_phenotype.R:787](../R/add_phenotype.R#L787). The silence is the worst part —
a user cannot tell their covariance was ignored.

### Defect 3 — correlated *named* random effects lose correlation on partial re-draw

Live bug, independent of time separation. Draws are persisted per
`(phenotype_name, effect_name, level)`
([phenotype_helpers.R:207-255](../R/phenotype_helpers.R#L207-L255)), but the
correlated path at [add_phenotype.R:454-556](../R/add_phenotype.R#L454-L556) does:

```text
all_levels  = union of levels across the correlated phenotypes
new_levels  = union over phenotypes of (all_levels \ levels already stored for that phenotype)
draws_mat   = MASS::mvrnorm(length(new_levels), 0, R_eff)     # jointly drawn
for each phenotype et:
    write draws_mat[ setdiff(new_levels, existing_lvls_for_et), et ]
```

If pen `P1` has a stored draw for `ADG` but not `BW`, a fresh **joint**
`(ADG_P1, BW_P1)` pair is drawn, `BW_P1` is written, and the `ADG_P1` component
is **discarded** in favour of ADG's previously stored, independently drawn value.
The realized pair on disk is uncorrelated. This needs conditioning, not storage.

### Observation — `distribution` is ignored in the correlated path

`define_effect_random()` accepts `"normal" | "gamma" | "uniform"`
([define_effect_random.R:23](../R/define_effect_random.R#L23)) and the marginal
path honours it, but [add_phenotype.R:530](../R/add_phenotype.R#L530) calls
`MASS::mvrnorm()` unconditionally. A gamma effect inside a covariance block is
silently drawn normal.

---

## 3. Settled: storage and timing

Both reviewers reached the same conclusion, so this section is compressed. The
full three-option analysis is in the v1 history; the outcome:

**Timing — lazy conditional sampling, not eager pre-drawing.** Eager is a strict
subset: it requires knowing every correlated phenotype × every individual before
the first record, fails for animals born later or subsets touched later, and
writes deviations for records that may never exist. Lazy **subsumes today's
behaviour** — nothing observed ⇒ conditional mean 0 and conditional variance = full
`R` ⇒ the existing joint MVN is the degenerate case. §8.5 gets *simpler*.

```text
e_n | e_o  ~  N( R_no R_oo⁻¹ e_o ,  R_nn − R_no R_oo⁻¹ R_on )
```

**Storage — Option A**: `residual_value DOUBLE` on `ind_phenotype`, named-effect
draws stay in `trait_random_effects`, one shared mathematical resolver across
both.

Rejected: a unified `random_draw` table (Option B) — residual is per-**record**,
named effects are per-**level**, and unifying them needs a `record_number` column
that is always `1` for everything but residual. That is the same "one column
meaning two things" smell that drove the `chr_meta` → `chr_inheritance` +
`chr_recombination` split. It also erases `id_ind` FK typing (naming rule 4) and
shadows `ind_phenotype` row-for-row, creating a deletion-sync liability. Codex
put it well: **share the resolver, keep the schemas separate.**

Also rejected: a dedicated `ind_residual` table (Option C) — Option A plus a
table, justified only by pre-drawn lifetime trajectories, which are not on the
roadmap. A → C later is one `INSERT ... SELECT`.

---

## 4. Response to the review round

### From Gemini — two helper fixes, adopted

Gemini's proposal is substantively the same design. It contributes two real
corrections to my v1 helper sketch:

1. **PSD clamping** on the conditional covariance diagonal. Adopted, and
   generalized per Codex point 3 below (clamping the diagonal is not sufficient —
   it must be an eigenvalue-level operation).
2. **The undefined `n`** in my unconditional branch. Adopted — though note
   Gemini's own fix, `n <- if (!is.null(observed)) nrow(observed) else 1L`,
   introduces a different bug: an unconditional draw would return **one row**
   instead of one per entity. Codex's resolution is correct: take the entity
   count explicitly as an argument.

Gemini otherwise drops the heterogeneous-residual, incomplete-matrix,
non-Gaussian, and repeated-record decisions, so nothing else is folded in.

### From Codex — eight corrections

| # | Codex point | Disposition |
|---|---|---|
| 1 | Discover the full covariance block, not just requested phenotypes | **Adopted** — I missed this; it would have preserved the original bug |
| 2 | Model observed / requested / latent coordinate state explicitly | **Adopted** |
| 3 | Numerically robust conditioning (PSD, near-singular, edge cases) | **Adopted** — my `solve(Roo)` was naive |
| 4 | Error on condition-level change, don't warn | **Adopted — I was wrong in v1** |
| 5 | `pheno_number` is ordinal identity, not time | **Adopted** — confirms and sharpens v1 Open Decision 1 |
| 6 | Persist `user_residual`; leave `user_values`/`derived_formula` NULL | **Adopted** — good catch, I lumped all three together |
| 7 | Resolve in memory, then write draws + records atomically | **Adopted** |
| 8 | Distinguish absent / incomplete / explicit-zero covariance | **Adopted with a modification** — see Open Question 1 |

Plus, adopted: `residual_value` in the **base schema** rather than on-demand
`ALTER TABLE`; validation of source compatibility within a named-effect block;
stable sort ordering before every RNG-consuming step; rejection of covariance
redefinition once draws exist; resolve full current-call state before the write
loop.

**Two places I push back — both confirmed by the author:**

**(a) No migration path for existing databases, and existing migration code gets
deleted.** Codex recommends explicit migrations and lists "existing databases
migrate without losing phenotype records" as test 38. **Rejected.** CLAUDE.md is
unambiguous: pre-1.0.0 there is no backward-compatibility obligation of any kind,
and leftover compatibility code is technical debt, not harmless history. Nobody is
building on tidybreed.

`residual_value` and `residual_condition_level` go into the base `ind_phenotype`
`CREATE TABLE`. Databases created before the change are **regenerated, not
migrated** — no `ALTER TABLE ... IF NOT EXISTS`, no `restore_pop()` fixup, no
conditional column checks anywhere in the new code.

Beyond not *adding* legacy paths, this work should *remove* the ones already in
the phenotype layer — see **Phase 0** in §8.

**(b) Decline the speculative provider abstraction for v1.** Codex proposes a
`resolve_stochastic_model(effect_name, phenotype_coordinates, entity_context,
observation_context)` interface so that future providers can return
time-evaluated covariance, random-regression coefficient models, kernels,
state-space transitions, or copulas. That is architecture for models that do not
exist, and it runs into CLAUDE.md's schema-bias rule: reserve a dimension only
when altering later would be painful; do not implement future biology before it
is needed.

What I *do* adopt from that section is the cheap and genuinely valuable part:
**the resolver must be database-independent and must take an opaque entity key**,
never assuming the key is an `id_ind` or that a level is time-invariant. That
single discipline is what keeps maternal, social, spatial, dyadic, and
time-indexed group effects reachable later without a rewrite, and it costs
nothing now. The rest of Codex's future-models catalogue (random regression,
survival, zero-inflation, copulas, reaction norms, AR/state-space residuals) is
excellent and should stay in
`sample_correlated_effects_codex_review.md` as the reference for Layers 2–3
rather than being duplicated or half-built here.

**One correction to my own v1 text:** I wrote that detecting a covariance change
after realization "requires a write timestamp we do not currently keep, so this
may be docs-only." That was wrong. No timestamp is needed — a single query for
"does any non-NULL `residual_value` exist for any phenotype in this block"
answers it. Codex's recommendation is adopted as a real v1 check.

---

## 5. Revised design

### 5.1 Schema

`ind_phenotype` gains two columns in its **base** `CREATE TABLE`
([define_trait.R:209-217](../R/define_trait.R#L209-L217)):

```sql
residual_value           DOUBLE,   -- liability-scale realized residual
residual_condition_level VARCHAR   -- condition level of the R used; NULL = unconditional
```

**Scale matters.** `resid` enters the model on the *liability* scale
([add_phenotype.R:793](../R/add_phenotype.R#L793):
`liability <- pheno_mean + covariate_contrib + tbv + resid`), before any
categorical/count conversion. `residual_value` therefore stores the
liability-scale residual, which is the only scale on which conditioning is valid.
For a threshold trait, cross-phenotype covariance is among **liabilities**, never
among observed categories.

`residual_condition_level` exists to make Codex point 4 enforceable: it records
which `R` a realization was drawn under, so a later call can refuse to condition
across incompatible conditions instead of silently mixing them.

Population policy:

| Path | `residual_value` |
|---|---|
| Ordinary model-generated residual | written |
| `user_residual` (explicit liability-scale residual) | **written** — conditions later draws |
| `user_values` (final phenotype, decomposition unknown) | `NULL` |
| `derived_formula` (no residual sampled) | `NULL` |

Directly supplied categories or values must never be reverse-engineered into a
latent residual.

`trait_random_effects` is unchanged: logical identity stays
`(phenotype_name, effect_name, level)`, values stay in `draw_value`.

`phenotype_meta` gains one column, per **D2**, in its base `CREATE TABLE`
alongside the existing `missing_component_action`:

```sql
condition_change_action VARCHAR DEFAULT 'error'   -- 'error' | 'independent'
```

### 5.2 Covariance block discovery

The single most important thing missing from v1. When the user calls
`add_phenotype("off_test_wt")`, the requested vector contains one name; the
sampler must still find `on_test_wt`.

```r
find_covariance_block(pop, effect_name, target_phenotypes)
```

Returns the **connected component** of the covariance graph containing the
targets — if the user declared {A, B, C} together, the block for A is {A, B, C}
even where `Cov(A,C)` was written as an explicit `0`. Same rule for residual and
named-effect blocks. Restricting the lookup to phenotypes named in the current call
would preserve the exact bug we are fixing.

Per **D1**, the returned block is guaranteed complete and PSD — that is enforced at
`define_*` time — so the sampler never has to reason about a partial matrix. A
phenotype in no block keeps the existing fast marginal `rnorm()` path, and
`load_phenotype_cov()` returning `NULL` now means exactly that and nothing else.

### 5.3 Coordinate state model

For each **entity**, classify every coordinate in the block:

- **observed** — already stored from an earlier call; conditions the distribution.
- **requested** — missing and being generated in this call; sampled and stored.
- **latent** — in the block but not requested; **must not be eagerly realized.**

Entity identity, per adapter:

```text
Residual adapter        entity = (id_ind, pheno_number)   coordinate = phenotype_name
                        storage = ind_phenotype.residual_value

Named-effect adapter    entity = (effect_name, level)     coordinate = phenotype_name
                        storage = trait_random_effects.draw_value
```

Group entities by identical `(observed, requested)` pattern; compute the
conditional coefficients and covariance **once per pattern**. In practice there
are one to three patterns per call.

This state model is what makes all of these fall out of one mechanism: A and B
together on identical sets; A on more individuals than B; B later on a culled
subset; B on a mixture of individuals with and without stored A; partial
named-effect vectors; and blocks of more than two phenotypes.

### 5.4 The resolver

Database-independent, opaque entity keys, no writes:

```r
resolve_correlated_draws <- function(covariance,        # validated, complete, dimnamed
                                     requested,         # chr vector of coordinates
                                     entity_keys,       # opaque; sorted by caller
                                     stored,            # entity x coordinate, NA = not observed
                                     tolerance)
```

Returns a normalized entity × coordinate matrix. The caller persists.

Numerical requirements (Codex point 3):

- Cholesky solve when `R_oo` is positive-definite; tolerance-aware eigen or
  pseudoinverse path when it is only positive-semidefinite (perfect correlation
  and zero-variance coordinates are both **valid** inputs).
- Symmetrize the computed conditional covariance.
- Clamp only tiny negative **eigenvalues** from floating-point error — clamping
  `diag(Sc)` as Gemini's sketch does is not sufficient.
- Reject materially invalid conditional covariance.
- Explicit handling for: zero entities, one entity, one requested coordinate,
  zero observed coordinates, zero conditional variance.
- Entity count passed explicitly, never inferred from `nrow(observed)`.

Validation of `R` itself moves **upstream** into `define_residual_cov()` and
`define_effect_cov_matrix()`: matching unique dimnames, symmetry within tolerance,
finite entries, non-negative diagonal, PSD within tolerance, **connected-component
completeness (D1)**, and **no realized draws already exist for the block (D3)**.
Sampling code should never discover a malformed matrix deep inside phenotype
generation, and must never interpret one as an instruction to draw independently.

### 5.5 `add_phenotype()` control flow

```text
resolve phenos, subsets, TBVs                            (unchanged)

── resolve ALL stochastic state in memory, before any write ──
for each named-effect block (§7.5):
    block   <- find_covariance_block(pop, effect, phenos)
    validate: every coordinate normal; compatible (source_column, source_table)
    stored  <- SELECT level, phenotype_name, draw_value FROM trait_random_effects
    draws   <- resolve_correlated_draws(...)             # conditions on stored

block <- find_covariance_block(pop, "residual", phenos)
stored <- SELECT id_ind, phenotype_name, pheno_number,
                 residual_value, residual_condition_level
            FROM ind_phenotype
           WHERE phenotype_name IN block AND id_ind IN <subset>
             AND residual_value IS NOT NULL
check: does any entity mix condition levels?             (D2)
    condition_change_action = 'error'       → stop
    condition_change_action = 'independent' → drop those cross terms,
                                              warn with count + <=5 example IDs
resid <- resolve_correlated_draws(...)

── then write, in one transaction ──
BEGIN
  INSERT new trait_random_effects rows
  INSERT ind_phenotype rows incl. residual_value, residual_condition_level
COMMIT
```

Two consequences worth calling out:

**Same-call differing subsets must be resolved before the write loop.** If A is
requested for individuals 1–100 and B for 51–100, then 1–50 get a marginal A
draw, 51–100 get a joint A/B draw, and no B row is created for 1–50. Resolving
current-call state up front means the result does not depend on the order the
per-phenotype loop happens to execute in. The `all_equal` restriction at
[add_phenotype.R:568-572](../R/add_phenotype.R#L568-L572) is **deleted**.

**Stable ordering before every RNG-consuming step.** Sort blocks, patterns,
coordinates, and entity keys. DuckDB does not guarantee row order (CLAUDE.md
design principle 7), and this path consumes RNG, so unsorted iteration would make
seeded output depend on physical row order. This is a hard requirement, not an
optimization.

### 5.6 Named-effect block validation

Before treating two coordinates as members of the same named-effect vector,
validate that their `trait_effects` rows agree on `(source_column, source_table)`.
A pen identifier must not be paired with a herd identifier merely because both
effects were given the same `effect_name` string.

Also validate that every coordinate in a correlated block uses
`distribution = "normal"` — error otherwise, until an explicit copula or
multivariate non-Gaussian API exists. This fixes the silent gamma→normal
substitution.

### 5.7 Repeated records

Narrow, explicit v1 contract:

- Residual covariance applies **across distinct phenotype names** only.
- Prior residual lookup matches on the **same `pheno_number`** across those
  phenotypes (ordinal pairing).
- Repeated records of the *same* phenotype have **independent** residuals.
- Persistent within-animal covariance is the job of a permanent-environment named
  effect, `define_effect_random(..., source_column = "id_ind")` — which already
  exists and is already documented
  ([define_effect_random.R:41-45](../R/define_effect_random.R#L41-L45)).

Conflating p.e. with residual covariance would double-count. Both reviewers agree
with this reading.

**Documentation must state that `pheno_number` supplies ordinal pairing, not
simulated-time matching**, and code must not assume `pheno_number` will later
become the time coordinate. Irregular longitudinal data needs an explicit
occasion/time model (Layer 2).

### 5.8 Worked example — a pen effect across two time points

Most of this plan is written in terms of residuals, which makes it easy to miss
that **named random effects go through the identical machinery with a different
entity key.** This section walks one all the way through. Setup: a pen effect on
both average daily gain and backfat, negatively correlated.

```r
pop <- pop |>
  define_effect_random("ADG", "pen", source_column = "pen_id", variance = 150) |>
  define_effect_random("BF",  "pen", source_column = "pen_id", variance = 4) |>
  define_effect_cov_matrix("pen",
    matrix(c(150, -12, -12, 4), 2, 2,
           dimnames = list(c("ADG", "BF"), c("ADG", "BF"))))
```

#### What happens today

**Day 0 — `add_phenotype("ADG")`.** §7.5 is gated on `length(phenos) >= 2`
([add_phenotype.R:456](../R/add_phenotype.R#L456)), so a single-phenotype call
skips the correlated path entirely. The marginal path
([phenotype_helpers.R:207-255](../R/phenotype_helpers.R#L207-L255)) draws
`rnorm(sd = sqrt(150))` for each pen level and stores
`(ADG, pen, P1, +8.3)`, `(ADG, pen, P2, -4.1)`, …

**Day 100 — `add_phenotype("BF")`.** Same gate, skipped again. No `(BF, pen, P1)`
row exists, so it draws `rnorm(sd = sqrt(4))` — **independent of the `+8.3`
already on disk.** The `-12` is never used.

If instead both phenotypes are named in one call, §7.5 runs and hits **Defect 3**:
it draws a fresh joint `(ADG, BF)` pair for P1, writes the BF half, and discards
the ADG half in favour of the stored `+8.3`. The stored BF is then correlated with
an ADG realization that was thrown away.

#### What happens under this plan

Entity identity for a named effect is `(effect_name, level)` — for residuals it is
`(id_ind, pheno_number)`. Everything else is the same code: block discovery,
observed/requested/latent classification, pattern grouping, conditional draw.

**Day 0.** `find_covariance_block(pop, "pen", "ADG")` returns `{ADG, BF}`. For pen
`P1`: observed `{}`, requested `{ADG}`, latent `{BF}`. Zero observed coordinates ⇒
the conditional draw degenerates to a marginal draw from `R[ADG, ADG] = 150`.
Store `(ADG, pen, P1, +8.3)`. **`BF` stays latent — not drawn, not stored**, per
§5.3.

**Day 100.** Block discovery returns `{ADG, BF}` again, this time from `BF`. For
`P1`: observed `{ADG: +8.3}`, requested `{BF}`, latent `{}`.

```text
mean = R_BF,ADG · R_ADG,ADG⁻¹ · e_ADG  = -12 · (1/150) · 8.3   = -0.664
var  = R_BF,BF - R_BF,ADG · R_ADG,ADG⁻¹ · R_ADG,BF
     = 4 - (-12)(1/150)(-12)                                   =  3.04

BF_P1 ~ N(-0.664, 3.04)
```

Store `(BF, pen, P1, …)`. The realized pair `(8.3, BF_P1)` now has the correct
joint distribution across a 100-day gap and two separate R calls.

**Day 200 — `add_phenotype("ADG")` for a new batch of animals in `P1`.** Observed
`{ADG, BF}`, requested `{}`. Nothing is drawn; the stored `+8.3` is reused exactly,
so every animal that ever passes through `P1` receives the same pen shift. That
reuse is existing behaviour ([phenotype_helpers.R:213-222](../R/phenotype_helpers.R#L213-L222))
and this plan does not change it.

#### Mixed patterns

Pens appear and disappear over time, so a single call routinely spans several
patterns. If at day 100 pens `P1`–`P5` carry a stored ADG draw and `P6`–`P10` are
brand new, there are two groups:

| Pattern | Pens | Draw |
|---|---|---|
| observed `{ADG}`, requested `{BF}` | P1–P5 | conditional, `N(-0.08·e_ADG, 3.04)` |
| observed `{}`, requested `{BF}` | P6–P10 | marginal, `N(0, 4)` |

Two conditional-covariance computations, one per pattern. `P6`–`P10` keep `ADG`
latent, to be drawn conditionally on their `BF` value if ADG is ever recorded for
them — which is the general form of Defect 3, fixed.

#### How named effects differ from residuals

| | Residual | Named effect (pen, litter, herd, p.e.) |
|---|---|---|
| Entity | `(id_ind, pheno_number)` | `(effect_name, level)` |
| Storage | `ind_phenotype.residual_value` | `trait_random_effects.draw_value` |
| Reuse | never — one draw per record | always — one draw per level, reused by every animal in it |
| Declared by | `define_residual_cov()` | `define_effect_cov_matrix(effect_name, …)` |
| D1 complete-block rule | applies | applies |
| D3 redefinition lock | applies | applies |
| D2 `condition_change_action` | applies | **does not** — `phenotype_var_comp.condition_column` is residual-only |
| Distribution | Gaussian by construction | must be `"normal"` in a correlated block (§5.6) |
| Extra validation | none | `(source_column, source_table)` must agree (§5.6) |

The reuse row is the one real semantic difference, and it is pre-existing: a
residual is consumed by exactly one record, whereas a pen draw is a property of the
pen and is applied to every animal in it, in every call, forever.

**Permanent-environment effects get all of this for free.** A p.e. effect is a
named effect with `source_column = "id_ind"`, so its entity is `("pe", <id_ind>)`.
That means the mechanism §5.7 points to for across-time correlation of the *same*
repeated phenotype is the same mechanism described here — including correlated p.e.
effects across two phenotypes recorded months apart.

---

## 6. Resolved decisions

The four v2 open questions were decided by the author. Recorded here with their
implications, because three of them tighten the contract of existing exported
functions.

### D1 — A covariance block is declared in one call, as a complete matrix

**Decision: strict error** (Codex's position; my implicit-zero recommendation was
not taken).

If the connected component containing any touched phenotype is not fully
specified, `define_residual_cov()` / `define_effect_cov_matrix()` **error and roll
back**, naming the missing pairs. So this fails:

```r
define_residual_cov(c("A", "B"), R_ab)
define_residual_cov(c("B", "C"), R_bc)   # error: block {A,B,C} missing Cov(A,C)
```

and this is the supported form:

```r
define_residual_cov(c("A", "B", "C"), R_abc)   # complete 3x3, zeros written explicitly
```

**Rationale.** The rule states in one sentence: *a residual covariance block is
declared in one call, as a complete matrix.* It also matches the domain — breeders
get full `R` matrices out of REML, not pairwise fragments. Most importantly, an
explicit `0` in the matrix is now the **only** way to say "in the same block but
uncorrelated," so declared independence is always something the user typed, never
something the package inferred. That is the failure class this plan exists to
eliminate, and implicit-zero would have reintroduced it in a quieter form.

**Three consequences, all simplifying:**

1. Completeness and PSD validation live entirely in `define_*`, never in the
   sampler (Codex point 8, fully resolved).
2. `load_phenotype_cov()` returning `NULL`
   ([define_effect_cov_matrix.R:334](../R/define_effect_cov_matrix.R#L334)) now has
   exactly **one** meaning — *no block exists* — so the sampler routes to the
   marginal path without ambiguity. Its current double meaning is what made Codex
   point 8 necessary.
3. `resolve_correlated_draws()` may assume a complete, validated, PSD matrix as a
   precondition.

### D2 — Condition-level change: error by default, with a stored opt-out

**Decision: error, plus an explicit opt-out.**

Erroring is the mathematically honest default — if A was drawn under `R_a` and B
under `R_b`, nothing defines `Cov(A|a, B|b)`, and my v1 "warn and use the new `R`"
was simply wrong. But a hard error would make heterogeneous residual variance by a
*moving* condition (farm, management group) unrunnable rather than approximate,
which is too blunt.

**Where the opt-out lives matters more than that it exists.** It is **not** an
`add_phenotype()` argument. It is a column on `phenotype_meta`, set by
`define_phenotype()`:

```sql
condition_change_action VARCHAR DEFAULT 'error'   -- 'error' | 'independent'
```

This is deliberate, and follows three existing conventions:

- It mirrors `phenotype_meta.missing_component_action` exactly — same table, same
  shape, same defaulting, same "one unified field, never per-type arguments"
  design that was settled when that field was added.
- CLAUDE.md design principle: configuration belongs in a table, not in an R
  argument, so `restore_pop()` is complete without re-supplying options.
- A per-call argument would let the same phenotype behave differently on different
  calls, which is exactly the kind of hidden state that makes a simulation
  irreproducible.

`'independent'` proceeds with the cross-condition covariance dropped and **warns
with a count plus up to 5 example IDs** — matching how `missing_component_action =
"skip"` already reports. Sex as a condition column never triggers either path.

`residual_condition_level` is stored regardless of which action is chosen, so a
future `define_residual_cross_condition_cov()` can be added without a schema
change.

### D3 — Covariance redefinition is rejected once draws exist

**Decision: error by default, `force = TRUE` to override.**

`define_residual_cov()` / `define_effect_cov_matrix()` check whether any realized
draw from the affected block already exists — one query for a non-`NULL`
`residual_value` (or a `trait_random_effects` row) in the block — and error if so.
No timestamp is needed; my v1 claim that one was is withdrawn.

`force = TRUE` proceeds with a warning stating plainly that **existing
realizations become distributionally inconsistent with the new matrix**. Note what
`force` cannot do: the phenotype values were already computed from the old draws,
so clearing `residual_value` would leave `pheno_value` unexplainable. `force` means
"I accept the inconsistency," never "repair it." The error message should say so.

This tightens two existing exported functions. Interactive model-building is the
main casualty; the intended workflow is to settle `R` before generating records,
which is also the only workflow that produces a coherent database.

### D4 — Column names

**Decision: `residual_value` and `residual_condition_level`.**

Both follow naming rule 5 (no abbreviations). The `residual_` prefix on the second
column is kept rather than reusing the bare `condition_level` from
`phenotype_var_comp`: `ind_phenotype` is a wide mixed table, and a bare
`condition_level` would not say which effect's condition it records — ambiguity
that would bite if named effects ever gain conditions.

---

## 7. Verification and testing

Conditional draws are **order-dependent**. `add_phenotype("A")` then
`add_phenotype("B")` will not produce the same numbers under one seed as
`add_phenotype(c("A", "B"))`, even though both are correct draws from the same
joint distribution — the two paths consume RNG differently. Per CLAUDE.md this is
fine (same seed + same call sequence reproduces), but it constrains the tests:
**never assert numeric equality between the sequential and joint paths.**

Distributional tests should use enough entities for stable estimates and derive
tolerances from sampling uncertainty, not fixed arbitrary margins.

### Residual blocks

1. Sequential A → B, identical individuals.
2. A → B on a culled subset (the motivating scenario); realized correlation on the
   overlap matches target; the non-overlapping animals are untouched.
3. A and B together on partially overlapping subsets.
4. B for a mixed group where only some individuals have a stored A.
5. B first, A later — order symmetry in distribution.
6. A/B/C block with multiple observation patterns.
7. Individuals with no prior coordinate get the correct marginal draw.
8. Disconnected components resolve independently.
9. An explicit zero covariance inside a declared block stays valid — the two
   phenotypes are in one block, draw jointly, and come out uncorrelated.
10. **D1**: declaring `{A,B}` then `{B,C}` errors, names the missing `(A,C)` pair,
    and leaves `phenotype_var_comp` unchanged (rollback).
11. **D1**: declaring the complete `{A,B,C}` in one call succeeds.
12. **D1**: a non-PSD or asymmetric matrix is rejected at `define_*` time, not in
    the sampler.

### Repeated records

1. Distinct phenotypes match on equal `pheno_number`.
2. Repeated records of the same phenotype stay residual-independent.
3. Unequal record counts never condition on the wrong record.

### Named effects

1. Defect 3 directly: pen `P1` stored for ADG only, then a BW call; across many
   pens the realized `(ADG, BW)` pairs reproduce `R_eff`.
2. Many levels with different partial-coordinate patterns.
3. New levels get correct unconditional joint draws; stored levels are reused
   exactly.
4. Incompatible `(source_column, source_table)` within a block is rejected.
5. Non-normal distribution in a correlated block is rejected.

### Numerical

1. One-dimensional conditional and unconditional draws.
2. Perfect and near-perfect correlation; zero conditional variance.
3. Tiny negative eigenvalues absorbed within tolerance; materially negative
   rejected.
4. Empty entity sets and fully-resolved requests are **RNG-neutral** (consume no
   random numbers).
5. Extreme but valid variance scales stay stable.

### Alternate input paths

1. `user_residual` is stored and conditions later draws.
2. `user_values` and `derived_formula` leave `residual_value` NULL and do not
   constrain later draws.
3. **D2**: a changed condition level errors under the default
   `condition_change_action = 'error'`.
4. **D2**: under `'independent'`, the same case proceeds, drops the cross terms,
   and warns with a count plus example IDs; marginal variances are still correct.
5. **D2**: an immutable condition column (sex) never triggers either path.
6. Threshold traits condition on latent liabilities, not observed categories.

### Reproducibility and integrity

1. Same seed + same call sequence → identical output.
2. Joint and sequential paths reproduce the target covariance statistically.
3. Forced write failure rolls back residuals, named draws, and phenotype rows
   together.
4. Results are independent of database row order (stable-ordering requirement).
5. **D3**: redefining a block errors once any realized draw exists.
6. **D3**: `force = TRUE` proceeds and warns; the warning states that existing
   realizations are now distributionally inconsistent.

---

## 8. Implementation sequence

**Phase 0 — delete the legacy migration code in the phenotype layer.** Per §4(a),
and independent of everything else in this plan. `ensure_trait_tables()`
([define_trait.R:161-425](../R/define_trait.R#L161-L425)) documents itself as
running "any necessary schema migrations on older databases." All of it is dead
code on a fresh database, and all of it is exactly the technical debt CLAUDE.md
tells us to remove. Verified redundant:

| Location | Legacy code | Status |
|---|---|---|
| [define_trait.R:370-383](../R/define_trait.R#L370-L383) | `trait_effects`: rename `trait_name` → `phenotype_name`; ALTER-ADD `source_table`, `slope`, `center`, `poly_order`, `null_class_action` | All five columns already in the base DDL |
| [define_trait.R:385-391](../R/define_trait.R#L385-L391) | `trait_random_effects`: rename `trait_name` → `phenotype_name` | Base DDL already uses `phenotype_name` |
| [define_trait.R:393-399](../R/define_trait.R#L393-L399) | `ind_phenotype`: rename `trait_name` → `phenotype_name` | Base DDL already uses `phenotype_name` |
| [define_trait.R:400-407](../R/define_trait.R#L400-L407) | `phenotype_meta`: ALTER-ADD `missing_component_action` (v0.31.x → v0.32.0) | Already in base DDL |
| [define_trait.R:408-415](../R/define_trait.R#L408-L415) | `phenotype_meta`: ALTER-ADD `formula_tbv`, `formula` (v0.32.x → v0.34.0) | Both already in base DDL |
| [define_trait.R:420-423](../R/define_trait.R#L420-L423) | `register_schema_meta()` wrapped in `if ("_schema_meta" %in% ...)` — "guarded for legacy databases" | Guard removable once the table is unconditional |
| [add_phenotype.R:779-782](../R/add_phenotype.R#L779-L782) | "Backward-compat fallback for databases without `phenotype_var_comp` residual rows" → second `get_phenotype_var()` lookup | Verify it is redundant with `get_residual_cov()`, then delete; §5.5 replaces this branch anyway |

After Phase 0, `ensure_trait_tables()` should be pure idempotent DDL — create if
absent, nothing else — and the phrase "schema migrations on older databases"
should be gone from its roxygen. The same rule applies to everything this plan
adds: **no conditional column checks, no `IF NOT EXISTS` fixups, no version
branches.**

**Phase 1 — validation and observability.** At `define_*` time: dimnames,
symmetry, finiteness, non-negative diagonal, PSD, **connected-component
completeness (D1)**, and **no-existing-draws with `force` override (D3)**. Reject
non-normal distributions in correlated named-effect blocks. Add tests that
demonstrate each of the three current silent-fallback defects *before* fixing them.

**Phase 2 — schema.** Add `residual_value` and `residual_condition_level` to the
base `ind_phenotype` schema, and `condition_change_action` to the base
`phenotype_meta` schema with the `define_phenotype()` argument that sets it (D2).
Store ordinary and `user_residual` residuals.

**Phase 3 — resolver.** Implement `find_covariance_block()` and
`resolve_correlated_draws()`. Test both **without any database access**, across
every coordinate pattern and numerical edge case.

**Phase 4 — residual integration.** Resolve all current-call residuals before the
per-phenotype write loop; condition on prior `residual_value`; delete the
`all_equal` restriction; write transactionally.

**Phase 5 — named-effect integration.** Replace the §7.5 pre-draw path with the
shared resolver plus a `trait_random_effects` adapter. Validate source
compatibility. Make named-effect and phenotype writes atomic.

**Phase 6 — documentation and performance.** Document sequential sampling, ordinal
repeated-record pairing, the one-call-complete-block rule (D1),
`condition_change_action` (D2), and the redefinition lock (D3). Add a culling
example to `add_phenotype()`. Benchmark large populations under
`dev/benchmarks/`, optimizing the observation-pattern query and batched writes
**without changing RNG semantics**.

Rough size: ~250 lines net across `add_phenotype.R`, `phenotype_helpers.R`, and a
new `R/correlated_draws.R`, with §7.5 and §8.5 both shrinking as their special-case
branches collapse into the pattern grouping.

---

## 9. Known limitations to document in v1

**Deletion and regeneration.** Under Option A, deleting a phenotype record deletes
its residual — referentially clean, but *not* probabilistically coherent: later
coordinates may already have been drawn conditional on the deleted value. Their
marginal distributions remain valid, but the sequential conditional history can no
longer be reconstructed. v1 documents that deleting and resampling correlated
phenotype records is not guaranteed to produce a coherent conditional history.
(This softens a "free deletion consistency" claim I made in v1 — it is free at the
referential level only.)

**Informative missingness.** The v1 resolver conditions only on **exact stored
deviations**. Culling on an already-observed phenotype is naturally supported.
Dropout driven by an *unobserved* latent value is informative missingness and
requires a joint dropout model; interval censoring carries partial information
about a liability that v1 cannot use.

**Count traits are a transformed Gaussian liability model**, not a Poisson or
negative-binomial model. The docs must say so rather than implying general
correlated count support.

**Concurrency.** DuckDB is single-writer per file, so cross-session races are
largely moot, but the transaction in Phase 4 is still required for crash safety.
Logical identities to enforce: residual `(id_ind, phenotype_name, pheno_number)`;
named effect `(phenotype_name, effect_name, level)`.

---

## 10. Explicitly out of scope (Layers 2–3)

Not implemented and not blocked. See `sample_correlated_effects_codex_review.md`
for the full catalogue and the requirements each would carry.

- **Layer 2 — longitudinal Gaussian**: observation times and shared occasion
  identity (distinct from `pheno_number`), random regression and reaction norms,
  permanent-environment trajectories, time-heterogeneous residual variance,
  covariance functions and AR/state-space residuals, time-indexed group exposures.
- **Layer 3 — non-Gaussian and event models**: generalized count traits,
  zero-inflated/hurdle processes, copula or latent-variable mixed-family
  dependence, survival with censoring and competing risks, informative dropout.

The one v1 design constraint that keeps these reachable: **the resolver takes an
opaque entity key and a covariance matrix, is database-independent, and never
assumes the key is an `id_ind` or that a level is time-invariant.** That is the
whole extension-point budget for now — no speculative provider interface, no
reserved columns.
