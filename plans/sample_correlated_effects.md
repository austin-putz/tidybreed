# Sampling correlated random effects at different points in simulated time

**Status**: **v3.1, approved; implementation in progress.** Phases 0 and 1
shipped 2026-09-20 (`sample_correlated_effects_phase_0.md`,
`sample_correlated_effects_phase_1.md`); Phases 2–8 not started.
File:line references are refreshed after each phase. v3 is a re-baseline
against the codebase as of v0.70.0 (2026-09-20) plus the Codex review of v2
(`sample_correlated_effects_v2_review.md`). Every v2 design decision stands
except one: D3's `force = TRUE` escape hatch, which Codex showed to be
incoherent and which is now removed. **D3, D6 and D7 were decided by the author
on 2026-09-20, accepting the v3 recommendations** (§6). v3.1 adds six
implementation-level clarifications found in the final read-through (marked
*v3.1* inline); none changes a decision.

**What changed from v2** (details in §0 and §4.5):

1. Phase 0 (delete legacy migration code) **already landed** in v0.64.0 and
   v0.63.x. It is replaced by a much smaller list of six stragglers.
2. `trait_effects` → `phenotype_effects` and `trait_random_effects` →
   `phenotype_random_effects` (v0.64.0). All references updated.
3. **New Defect 4**: heterogeneous (conditional) residual variance is silently
   ignored on every single-phenotype call. Found while re-verifying; the existing
   test asserts only a row count.
4. **New bypass**: a third writer to `phenotype_var_comp`,
   `write_phenotype_var_diag()`, rewrites one diagonal in place and would slip
   past D1 and D3 as v2 wrote them. Closed by sharpening D1 and adding **D5**.
5. The repository's **RNG discipline** — `dbWriteTable()` advances R's RNG,
   register + `INSERT` does not — is now a hard rule across the package. The
   plan was silent on it; §5.5 now requires it.
6. `resolve_subset_ids()` (v0.70.0) returns sorted ids, which discharges half of
   the stable-ordering requirement for free.
7. The distribution check (§5.6) was placed at `add_phenotype()` in one section
   and at `define_*` time in another. Now placed at both, with the reason.
8. `liability_value` and `cat_name` are still added by on-demand `ALTER TABLE`,
   which the plan's own §4(a) rule forbids. Folded into the Phase 1 schema change.
9. All file:line references refreshed to v0.70.0.
10. **Codex's v2 review folded in** (§4.5): record planning before any draw
    (three-stage `add_phenotype()`), block identity by pair-row existence with
    an explicit replacement algorithm, a fourth **fixed** coordinate state for
    `user_residual`, support-consistency checking for singular `R_oo`,
    block-level agreement on `condition_change_action`, precise `'independent'`
    semantics, the RNG-on-error contract, thirteen more tests, and a
    re-sequenced implementation. One v2 decision is reopened: **D3's `force`**.

**What changed from v3 to v3.1** (final read-through, no decision changes):

- §5.1: the "no `ALTER TABLE` left in `add_phenotype.R`" claim is scoped to
  reserved columns — user `...` columns still go through
  `prepare_extra_cols()`, inside the Stage-3 transaction.
- §5.5: planned ids are registered as a temporary view and joined, never
  pasted into SQL text; the stratum lookup's contract for `condition_table`
  (one row per planned `id_ind`, `NULL` = no matching stratum) is stated; the
  resolver draws with `stats::rnorm()` through a factor of the conditional
  covariance, not `MASS::mvrnorm()`.
- §5.6: the `"normal"` requirement applies to blocks of **two or more**
  coordinates; a 1 × 1 non-normal named effect keeps the existing marginal
  gamma/uniform sampler.
- D3: the realization-lock predicate is defined (`residual_value IS NOT NULL`);
  `archive_replicate()` resets it, so between-replicate redefinition needs no
  `remove_rows()`.
- §7/§8: schema phase moved ahead of the writer phase (the residual lock
  predicate depends on it); "tests before the fix" made practical; `_schema_meta`
  column descriptions added to the schema phase; three tests added.

**Scope name**: this is **Layer 1 — fixed multivariate Gaussian blocks sampled
across pipeline stages.** It is not a longitudinal, random-regression, survival,
or non-Gaussian dependence framework. Saying so up front matters, because the
failure mode of this feature is users assuming it means more than it does.

**Author intent** (unchanged): a user who wrote
`define_residual_cov(c("on_test_wt", "off_test_wt"), R)` should get the covariance
they asked for whether the two phenotypes are recorded in one call or a hundred
simulated days apart, with culling in between — with zero new user-facing
concepts.

---

## 0. What changed in the codebase since v2

v2 was written against ~v0.60. Nine releases later:

| Change | Release | Effect on this plan |
|---|---|---|
| Migration section of `ensure_trait_tables()` and `.migrate_var_comp_tables()` deleted | v0.64.0 (`3fa68a8`) | v2 Phase 0 is done; see §8 for the stragglers that remain |
| Dead `dbExistsTable`/`has_meta` guards removed across `define_effect_*`, `schema()`, `add_phenotype()` | v0.63.x (`d882556`) | Same |
| `trait_effects` → `phenotype_effects`; `trait_random_effects` → `phenotype_random_effects` | v0.64.0 | Every reference in this plan renamed |
| `TABLE_RESERVED_COLS` gained entries for `phenotype_random_effects`, `phenotype_components`, `founder_haplotypes` | v0.64.0 | New `ind_phenotype` columns must be added to the reserved list (§8) |
| RNG discipline: `dbWriteTable()` advances R's RNG by a fixed amount (random temp-name generation); `duckdb_register()` + `INSERT` is RNG-neutral. Documented in `define_genome.R`, `founder_haplotype_helpers.R`, `add_offspring.R`, [define_effect_cov_matrix.R:129](../R/define_effect_cov_matrix.R#L129) | v0.5x–0.6x | New hard requirement in §5.5 |
| Transaction idiom standardized: `BEGIN TRANSACTION` + `on.exit(ROLLBACK unless committed)` + validator + `COMMIT` ([define_chromosome.R:219-227](../R/define_chromosome.R#L219-L227)) | v0.6x | §5.5 adopts it verbatim |
| `resolve_subset_ids()` — one SQL statement, ids returned **sorted** ([sql_utils.R:494](../R/sql_utils.R#L494)) | v0.70.0 | Residual entity keys arrive sorted; §5.5 ordering requirement narrows to blocks, patterns, coordinates, and named-effect levels |
| `remove_rows()` supports single-table deletion on every table except `_schema_meta`, keyed by `TABLE_ROW_KEYS` | v0.6x | The escape hatch that replaces D3's `force` (§6 D3) |
| `archive_replicate()` resets `ind_phenotype` and `phenotype_random_effects` per replicate | existing | Each replicate starts with **no** realizations, so nothing conditions across replicates. New columns travel with the table copy; no change needed |
| Genome-effects Phases B–E reshaped `add_phenotype()`'s TBV section | v0.65–0.68 | Line numbers shifted ~11 lines; no design impact |

`add_phenotype()` still has **no transaction** around its writes and still calls
`dbWriteTable()` at [add_phenotype.R:565](../R/add_phenotype.R#L565) (named-effect
draws) and [add_phenotype.R:869](../R/add_phenotype.R#L869) (records), both of
which sit *between* RNG draws. That interleaving is one of the things §5.5
removes.

---

## 1. The motivating scenario

Pigs on test. On-test weight at entry; off-test weight ~100–160 days later, with
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
`id_ind` sets ([add_phenotype.R:579-584](../R/add_phenotype.R#L579-L584)), and
culling guarantees they differ.

**The two ways to fix it, and which one this plan takes.** Either (a) draw both
residuals when the first phenotype is recorded and stash the second for later,
or (b) draw one, persist it, and draw the other *conditionally* when its turn
comes. This plan takes **(b)**. (a) requires knowing at day 0 every correlated
phenotype × every individual that will ever be recorded, fails for animals born
later, and writes residuals for records that never happen — culled pigs never get
an off-test weight. (b) subsumes today's behaviour: with nothing observed, the
conditional mean is 0 and the conditional variance is the full `R`, so the
current same-call joint MVN is just the degenerate case. See §3.

---

## 2. Confirmed defects

Defects 1–3 were independently confirmed by both v2 reviewers and re-verified
against v0.70.0. Defect 4 is new in v3.

### Defect 1 — residual deviations are never persisted

`sample_residuals()` ([phenotype_helpers.R:46-66](../R/phenotype_helpers.R#L46-L66))
returns an in-memory matrix, consumed at
[add_phenotype.R:775-798](../R/add_phenotype.R#L775-L798) to build the liability
and then discarded. A later call has nothing to condition on, so cross-call
correlation is not merely unimplemented — it is *unimplementable* without a
storage change. **This is a storage problem, not a missing sampling branch.**

### Defect 2 — the joint residual MVN is gated on same-call *and* identical ID sets

[add_phenotype.R:575-649](../R/add_phenotype.R#L575-L649) draws a joint
`MVN(0, R)` only when `length(phenos) >= 2`, `all_equal` (byte-identical sorted
`id_ind` vectors), and `R_unconditional` is non-`NULL`. Otherwise it falls
through **silently** to the independent `rnorm()` at
[add_phenotype.R:793](../R/add_phenotype.R#L793). The silence is the worst part —
a user cannot tell their covariance was ignored.

### Defect 3 — correlated *named* random effects lose correlation on partial re-draw

Live bug, independent of time separation. Draws are persisted per
`(phenotype_name, effect_name, level)`
([phenotype_helpers.R:208-250](../R/phenotype_helpers.R#L208-L250)), but the
correlated path at [add_phenotype.R:463-569](../R/add_phenotype.R#L463-L569) does:

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

### Defect 4 (new in v3) — heterogeneous residual variance is ignored on single-phenotype calls

`get_residual_cov()` ([phenotype_helpers.R:275-358](../R/phenotype_helpers.R#L275-L358))
builds `R_by_level` from the `condition_column` rows, but the only consumer is
the §8.5 branch, which is inside `if (length(phenos) >= 2 && all_equal)`. A
single-phenotype call — the overwhelmingly common case — reaches the independent
draw at [add_phenotype.R:787-793](../R/add_phenotype.R#L787-L793), which reads
`residual_var_unconditional[t]` and never looks at the level. So

```r
define_phenotype("BW", residual_var = 600)
define_residual_cov("BW", matrix(400), condition_column = "sex", condition_level = "M")
define_residual_cov("BW", matrix(800), condition_column = "sex", condition_level = "F")
pop |> get_table("ind_meta") |> add_phenotype("BW")
```

draws every animal from `N(0, 600)`. The existing test at
[test-phenotype_composite.R:198-226](../tests/testthat/test-phenotype_composite.R#L198-L226)
says *"variances should differ by sex"* in a comment and asserts only
`nrow(ph) == 200`. If the unconditional row is absent, the same call errors
"No residual variance found" even though two conditional rows exist.

This is not a separate feature gap. D2 (§6) assumes per-level `R` is actually
applied; the block-resolver design fixes it because a size-1 block with
conditions routes through the resolver like any other. It needs to be named,
tested before the fix, and tested after.

### Observation — `distribution` is ignored in the correlated path

`define_effect_random()` accepts `"normal" | "gamma" | "uniform"`
([define_effect_random.R:23](../R/define_effect_random.R#L23)) and the marginal
path honours it, but [add_phenotype.R:539](../R/add_phenotype.R#L539) calls
`MASS::mvrnorm()` unconditionally. A gamma effect inside a covariance block is
silently drawn normal.

### Observation (new in v3) — a third writer bypasses the block

`phenotype_var_comp` has **three** writers, not two:

| Writer | Called by | What it writes |
|---|---|---|
| `define_residual_cov()` | user | full `n × n` residual block, one condition slice |
| `define_effect_cov_matrix()` | user | full `n × n` block for a named effect |
| `write_phenotype_var_diag()` ([define_effect_cov_matrix.R:320-340](../R/define_effect_cov_matrix.R#L320-L340)) | `define_phenotype(residual_var = )`, `define_effect_random(variance = )` | **one diagonal cell**, in place, off-diagonals untouched |

Called after a block exists, the third writer can turn a valid `R` into a
non-PSD one, and it can do so after draws have been realized. v2's D1 and D3
named only the first two writers. Closed in §6 (D1 sharpened, D5 added).

### Observation (new in v3, from Codex) — individuals are excluded *after* the draw

Three exclusion points run inside the per-phenotype loop, **after** the §8.5
joint draw has already consumed RNG for them: `null_class_action = "skip"`
([add_phenotype.R:655-670](../R/add_phenotype.R#L655-L670)), formula-TBV
missing components ([add_phenotype.R:706-725](../R/add_phenotype.R#L706-L725)),
and composite-TBV missing components
([add_phenotype.R:735-749](../R/add_phenotype.R#L735-L749)). Today that only
wastes draws. Under Option A it is worse: a residual for an excluded animal has
no `ind_phenotype` row to live in, and seeded output depends on downstream
filtering details. This is why §5.5 plans records **before** drawing anything.

---

## 3. Settled: storage and timing

Both v2 reviewers reached the same conclusion, so this section is compressed.
The full three-option analysis is in the v1 history; the outcome:

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
draws stay in `phenotype_random_effects`, one shared mathematical resolver across
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

## 4. Response to the review rounds

### 4.1 From Gemini (v1 round) — two helper fixes, adopted

1. **PSD clamping** on the conditional covariance — adopted and generalized to
   an eigenvalue-level operation (Codex point 3); clamping `diag(Sc)` is not
   sufficient.
2. **The undefined `n`** in the unconditional branch — adopted, with Codex's
   resolution (entity count passed explicitly) rather than Gemini's
   (`nrow(observed) else 1L`, which returns one row for an unconditional draw).

### 4.2 From Codex (v1 round) — eight corrections

| # | Codex point | Disposition |
|---|---|---|
| 1 | Discover the full covariance block, not just requested phenotypes | **Adopted** |
| 2 | Model observed / requested / latent coordinate state explicitly | **Adopted**; extended to four states in v3 (§5.3) |
| 3 | Numerically robust conditioning (PSD, near-singular, edge cases) | **Adopted**; extended with support-consistency in v3 (§5.4) |
| 4 | Error on condition-level change, don't warn | **Adopted** (D2) |
| 5 | `pheno_number` is ordinal identity, not time | **Adopted** |
| 6 | Persist `user_residual`; leave `user_values`/`derived_formula` NULL | **Adopted** |
| 7 | Resolve in memory, then write draws + records atomically | **Adopted**; extended to plan → resolve → commit in v3 (§5.5) |
| 8 | Distinguish absent / incomplete / explicit-zero covariance | **Adopted via D1** |

Plus, adopted: `residual_value` in the **base schema** rather than on-demand
`ALTER TABLE`; validation of source compatibility within a named-effect block;
stable sort ordering before every RNG-consuming step; rejection of covariance
redefinition once draws exist; resolve full current-call state before the write
loop.

**Two pushbacks, both confirmed by the author:**

**(a) No migration path for existing databases.** CLAUDE.md is unambiguous:
pre-1.0.0 there is no backward-compatibility obligation of any kind. New columns
go into the base `CREATE TABLE`; databases created before the change are
regenerated, not migrated — no `ALTER TABLE ... IF NOT EXISTS`, no
`restore_pop()` fixup, no conditional column checks anywhere in the new code.
*(v3 note: the deletion of the pre-existing migration code that v2 attached to
this pushback has since landed independently — see §0 and §8 Phase 0.)*

**(b) No speculative provider abstraction for v1.** Codex's
`resolve_stochastic_model(...)` interface is architecture for models that do not
exist. What *is* adopted: **the resolver must be database-independent and must
take an opaque entity key**, never assuming the key is an `id_ind` or that a
level is time-invariant. That single discipline keeps maternal, social, spatial,
dyadic, and time-indexed group effects reachable later, and it costs nothing now.

### 4.5 Response to the Codex v2 review *(new in v3)*

Codex reviewed v2 (`sample_correlated_effects_v2_review.md`) in parallel with
the v3 re-baseline; several of its points had already been reached
independently (table names, Phase 0, the third writer, strict-subset
redeclaration). The rest are disposed of here. Items marked **decision** reverse
or extend an author decision and are written up in §6 with a recommendation.

| # | Codex point | Disposition | Where |
|---|---|---|---|
| B1 | Final record eligibility must be resolved **before** any RNG is consumed; today three exclusion points shrink `ids_t` *after* the §8.5 joint draw | **Adopted.** Verified against current code (§2, last observation). `add_phenotype()` becomes three explicit stages: **plan → resolve → commit**. This is a real restructuring, and Codex is right that ~250 net lines was optimistic | §5.5, §8 Phase 4 |
| B2 | Block identity by **pair-row existence** (an explicit `0` is still a pair row); exact replacement rule `N == U`; stratum = `(effect_name, condition_table, condition_column, condition_level)`; identical coordinate sets across levels; no cell-by-cell fallback | **Adopted.** Same rule as v3's sharpened D1, now stated as an algorithm. Added: a block has at most **one** `(condition_table, condition_column)`; only levels vary | §5.2, D1 |
| B3 | `force = TRUE` makes later conditional draws incoherent: `A` drawn under `R1`, block forced to `R2`, `B` then conditioned on `A` using `R2` — the pair has no declared joint distribution, and nothing on disk can tell the two eras apart | **Agree — decision.** Recommend removing `force` in v1. The escape hatch already exists: `remove_rows()` on the block's realizations, then redefine | D3 |
| B4 | Singular `R_oo` needs a **support-consistency** check `(I − R_oo R_oo⁺) e_o ≈ 0` before conditioning; scale-aware eigenvalue tolerance | **Adopted.** A pseudoinverse returns a number for an off-support observation, not a conditional draw | §5.4 |
| B5 | `user_residual` is a fourth coordinate state, **fixed**: current, caller-supplied, not drawn, but conditions same-call draws and is stored | **Adopted**, including a named list that supplies residuals for only a *subset* of the current phenotypes — the fixed/sample split makes that natural, and the current all-or-nothing branch goes | §5.3 |
| B6a | `condition_change_action` is per phenotype but the decision is per block — what if members disagree? | **Decision.** Recommend requiring agreement, enforced at sampling and, where the metadata already exist, at definition time | D6 |
| B6b | `'independent'` on blocks > 2: drop only stored coordinates whose stratum differs, keep conditioning on compatible ones; record the **selected** stratum (`NULL` for the unconditional fallback), never the raw column value; error when neither a matching slice nor a complete unconditional `R` exists — never `0` | **Adopted** | D2 |
| C1 | `trait_random_effects` → `phenotype_random_effects` throughout | Done in v3 | — |
| C2 | Phase 0 is stale; fold the `add_phenotype()` fallback into residual integration | Done in v3 as a short straggler list. Kept as a separate phase because most of the items are outside `add_phenotype()` | §8 |
| C3 | "No block" vs a 1 × 1 block is ambiguous | **Adopted.** A scalar variance *is* a complete 1 × 1 block; no rows at all is an error (already is: "No residual variance found"). There is one path — the resolver; the marginal `rnorm()` is what it degenerates to for one coordinate and zero observed | §5.2 |
| C4 | A stored pen draw persists for every animal ever in `P1`, across batches and seasons, unless occasion is part of the level | **Adopted** as a highlighted note beside the example and in the `define_effect_random()` docs | §5.8 |
| R | Resolver contract: `observed` = stored + fixed; return values for `sample_coordinates` only; nine adapter responsibilities | **Adopted** | §5.4 |
| T | Database transaction ≠ atomic operation: an R error after drawing advances `.Random.seed` even when the write rolls back. Pick a contract | **Decision.** Recommend *database atomic, RNG advances on failure* — consistent with every other RNG-consuming function in the package (nothing in `R/` touches `.Random.seed`) | D7 |
| S | Re-sequenced implementation with record planning as its own step | **Adopted.** New Phase 4; estimate revised | §8 |
| + | Thirteen additional tests | **Adopted** | §7 |

---

## 5. Revised design

### 5.1 Schema — ✅ shipped (Phase 1, 2026-09-20)

`ind_phenotype`'s base `CREATE TABLE`
([define_trait.R:204-214](../R/define_trait.R#L204-L214)) gained **four**
columns — the two this plan needs, plus the two that were previously added by
on-demand `ALTER TABLE` in violation of §4(a):

```sql
CREATE TABLE ind_phenotype (
  id_phenotype             INTEGER PRIMARY KEY,
  id_ind                   VARCHAR,
  phenotype_name           VARCHAR,
  pheno_value              DOUBLE,
  pheno_number             INTEGER,
  liability_value          DOUBLE,    -- was on-demand ALTER; NULL unless store_liability
  cat_name                 VARCHAR,   -- was on-demand ALTER; NULL unless categorical + cat_names
  residual_value           DOUBLE,    -- liability-scale realized residual
  residual_condition_level VARCHAR    -- the SELECTED stratum's level; NULL = unconditional R was used
)
```

There are now **no** `dbListFields()` checks or `ALTER TABLE` statements for
**reserved** columns in `add_phenotype.R`. The `store_liability` and `cat_names`
paths populate columns that already exist
([add_phenotype.R:855-866](../R/add_phenotype.R#L855-L866)). *(v3.1)* User `...`
columns are a different matter and are unchanged: they still go through
`prepare_extra_cols()` ([sql_utils.R:403](../R/sql_utils.R#L403)),
which issues `ALTER TABLE ADD COLUMN` for a column the table has not seen
before. Under §5.5 that call moves inside the Stage-3 transaction (DuckDB DDL is
transactional), ahead of the column-listed `INSERT`, and it is RNG-neutral. The
§7 "no `ALTER TABLE` is issued" assertion is therefore for a call **without**
`...`.

**Scale matters.** `resid` enters the model on the *liability* scale
([add_phenotype.R:798](../R/add_phenotype.R#L798):
`liability <- pheno_mean + covariate_contrib + tbv + resid`), before any
categorical/count conversion. `residual_value` therefore stores the
liability-scale residual, which is the only scale on which conditioning is valid.
For a threshold trait, cross-phenotype covariance is among **liabilities**, never
among observed categories. `liability_value` (the whole liability) and
`residual_value` (its residual component) are different quantities; storing both
is not redundant.

`residual_condition_level` records **which stratum's `R` was actually used**,
not the raw value of the condition column (Codex B6b). If an animal's level
matched no conditional slice and the unconditional `R` was used, the stored
value is `NULL`. This is what makes D2 enforceable: a later call can see exactly
which `R` a realization was drawn under and refuse to condition across
incompatible strata instead of silently mixing them. Because a block has at most
one `(condition_table, condition_column)` (D1), the level alone identifies the
stratum.

Population policy:

| Path | `residual_value` |
|---|---|
| Ordinary model-generated residual | written |
| `user_residual` (explicit liability-scale residual) | **written** — a *fixed* coordinate (§5.3); conditions same-call and later draws |
| `user_values` (final phenotype, decomposition unknown) | `NULL` |
| `derived_formula` (no residual sampled) | `NULL` |

Directly supplied categories or values must never be reverse-engineered into a
latent residual.

`phenotype_random_effects` is unchanged: logical identity stays
`(phenotype_name, effect_name, level)`, values stay in `draw_value`.

`phenotype_meta` ([open_pop.R:289-308](../R/open_pop.R#L289-L308)) gained
one column, per **D2**, alongside the existing `missing_component_action`, and
`define_phenotype()` the `condition_change_action = c("error", "independent")`
argument that sets it (validated with `match.arg()`, stored on the row,
replaced by `overwrite = TRUE` like every other field):

```sql
condition_change_action VARCHAR DEFAULT 'error'   -- 'error' | 'independent'
```

`TABLE_RESERVED_COLS` ([sql_utils.R:110-112](../R/sql_utils.R#L110-L112))
covers all nine `ind_phenotype` columns and `condition_change_action`; all
five new columns have `_schema_meta` descriptions. Until Phase 5 writes them,
`residual_value` and `residual_condition_level` are `NULL` on every row, which
the Phase 1 tests assert (`test-phenotype_schema.R`).

`archive_replicate()` already lists both realization tables under
`store_and_reset`, so every replicate begins with no stored coordinates; the new
columns are copied with the table. No change.

### 5.2 Covariance block discovery

The single most important thing missing from v1. When the user calls
`add_phenotype("off_test_wt")`, the requested vector contains one name; the
sampler must still find `on_test_wt`.

```r
find_covariance_block(pop, effect_name, target_phenotypes)
```

Returns the **connected component** of the covariance graph containing the
targets. **The graph is defined by the existence of a stored pair row, never by
`cov_value != 0`** (Codex B2): if the user declared {A, B, C} together, the block
for A is {A, B, C} even where `Cov(A,C)` was written as an explicit `0`. Same
rule for residual and named-effect blocks. Restricting the lookup to phenotypes
named in the current call would preserve the exact bug we are fixing.

**Strata.** For the residual effect, rows are partitioned into strata

```text
(effect_name, condition_table, condition_column, condition_level)
```

and the graph is computed per stratum. D1 requires (i) every stratum of a block
to name the same phenotype set and (ii) a block to have at most one
`(condition_table, condition_column)`, so block membership is well-defined
regardless of which stratum a given individual resolves to, and a single stored
level identifies the stratum. There is no cell-by-cell fallback between a
partial conditional matrix and the unconditional one — a stratum is complete or
it is rejected at definition time.

**A scalar variance is a complete 1 × 1 block** (Codex C3). There is no separate
"no block" marginal path: every phenotype with any residual row goes through
the resolver, and for one coordinate with nothing observed the resolver *is*
`rnorm(sd = sqrt(v))`. A phenotype with **no** residual row at all is an error,
as it is today ("No residual variance found"). `load_phenotype_cov()` returning
`NULL` ([define_effect_cov_matrix.R:275](../R/define_effect_cov_matrix.R#L275))
therefore has exactly one meaning: no rows exist.

Per **D1**, a discovered block is guaranteed complete and PSD — that is enforced
by the writers — so the sampler never has to reason about a partial matrix.

### 5.3 Coordinate state model

For each **entity**, classify every coordinate in the block into one of **four**
states (Codex B5 added the second):

- **stored** — realized in an earlier call; read from disk; conditions the draw.
- **fixed** — realized *in this call* by the caller (`user_residual`); not drawn;
  conditions the draw; written to disk.
- **sample** — being generated in this call; the only state that consumes RNG;
  written to disk.
- **latent** — in the block but neither requested nor realized; **must not be
  eagerly realized.**

The conditional observed set is `stored ∪ fixed`. If A and B are requested
together and the caller supplies A's residual, B is drawn conditionally on the
supplied A in the same call. If both are supplied, nothing is drawn and both are
stored to condition a later C.

**`user_residual` contract under this model.** A single phenotype takes a
numeric vector; multiple phenotypes take a named list keyed by `phenotype_name`,
and **the list may name only a subset of the current phenotypes** — the rest are
sampled. Each vector is positional over that phenotype's *planned* record list
(§5.5 Stage 1), which is the same "final per-phenotype individual list" the
current docs describe ([add_phenotype.R:57-64](../R/add_phenotype.R#L57-L64)),
now with a precise definition because planning precedes drawing.

Entity identity, per adapter:

```text
Residual adapter        entity = (id_ind, pheno_number)   coordinate = phenotype_name
                        storage = ind_phenotype.residual_value

Named-effect adapter    entity = (effect_name, level)     coordinate = phenotype_name
                        storage = phenotype_random_effects.draw_value
```

Group entities by identical `(stratum, observed coordinates, sample
coordinates)` pattern; compute the conditional coefficients and covariance
**once per pattern**. In practice there are one to three patterns per call.

This state model is what makes all of these fall out of one mechanism: A and B
together on identical sets; A on more individuals than B; B later on a culled
subset; B on a mixture of individuals with and without stored A; a supplied A
with a generated B; partial named-effect vectors; blocks of more than two
phenotypes; **and a size-1 block with heterogeneous variance (Defect 4)** — a
single sample coordinate, zero observed, drawn from the stratum matching the
entity's condition level.

### 5.4 The resolver

Database-independent, opaque entity keys, no writes, no knowledge of strata
(the adapter has already picked one matrix per pattern):

```r
resolve_correlated_draws <- function(covariance,          # one stratum's R: validated, complete, dimnamed
                                     sample_coordinates,  # chr; coordinates to draw
                                     entity_keys,         # opaque; sorted by caller
                                     observed,            # entity x coordinate; stored + fixed; NA = not observed
                                     tolerance)
# returns entity x sample_coordinates, in the caller's entity order
```

**Resolver responsibilities** (and nothing else):

- validate dimensions, dimnames, and finiteness;
- **support consistency for singular `R_oo`** (Codex B4): an observed vector
  must lie in the support of its Gaussian, `‖(I − R_oo R_oo⁺) e_o‖ ≈ 0` within a
  scale-aware tolerance; otherwise **error**. A nonzero supplied residual for a
  zero-variance coordinate has probability zero and defines no conditional
  distribution — a pseudoinverse would still return a number;
- Cholesky solve when `R_oo` is positive-definite; eigen/pseudoinverse path when
  it is only PSD (perfect correlation and zero-variance coordinates are both
  **valid** inputs);
- symmetrize the conditional covariance; project only numerically tiny negative
  eigenvalues to zero; reject a materially indefinite result;
- **scale-aware tolerance**: a function of dimension, the largest eigenvalue,
  and machine epsilon — never one absolute number reused across variances of
  different magnitude;
- explicit handling for zero entities, one entity, one sample coordinate, zero
  observed coordinates, zero conditional variance;
- entity count passed explicitly, never inferred from `nrow(observed)`;
- sample in stable order; **consume no RNG** when there is nothing to draw.

**Adapter responsibilities** (residual adapter and named-effect adapter; Codex R):

1. Discover the complete block by pair-row existence.
2. Select the exact stratum for every planned entity.
3. Produce stable entity and coordinate order.
4. Load stored values and check their stratum compatibility (D2).
5. Add fixed current values (`user_residual`).
6. Group entities by `(stratum, observed coordinates, sample coordinates)`.
7. Call the resolver once per pattern.
8. Merge fixed and sampled values into the planned rows.
9. Persist only inside the outer transaction.

Validation of `R` itself lives **upstream** in the three `phenotype_var_comp`
writers (§5.9): matching unique dimnames, symmetry within tolerance, finite
entries, non-negative diagonal, PSD within tolerance, **block replacement rule
(D1)**, and **no realized draws for the block (D3)**. Sampling code should never
discover a malformed matrix deep inside phenotype generation, and must never
interpret one as an instruction to draw independently.

### 5.5 `add_phenotype()` control flow — three stages

*(Restructured in v3 per Codex B1.)* Today the final population for each phenotype is only known deep
inside the per-phenotype loop, after §8.5 has already drawn for everyone. The new
flow separates planning, resolution, and commit so that **no random number is
consumed for a record that will not exist** and **no write happens before every
draw is done**.

```text
── Stage 1: PLAN — no RNG, no writes ──────────────────────────────────────────
resolve phenos, subsets (sorted ids via resolve_subset_ids()), TBVs   (as today)
for each phenotype t:
    apply expressed_sex, repeatable guard                              (as today, moved up)
    compute covariate contributions; apply null_class_action = "skip"  (moved up)
    evaluate formula_tbv / composite TBV; apply missing_component_action (moved up)
    classify path: model | user_values | derived_formula
    assign pheno_number for every planned record                       (next_pheno_numbers)
    look up each planned record's residual stratum (condition level)
    collect named-effect levels each planned record touches
→ one in-memory plan: for each phenotype, the final record list with
  id_ind, pheno_number, tbv, covariate_contrib, stratum, path

── Stage 2: RESOLVE — RNG, no writes ──────────────────────────────────────────
for each named-effect block (persistent per-level entities):
    block  <- find_covariance_block(pop, effect, phenos)
    validate: every coordinate normal; compatible (source_column, source_table)
    stored <- SELECT level, phenotype_name, draw_value FROM phenotype_random_effects
    draws  <- adapter → resolve_correlated_draws(...)  per pattern

block  <- find_covariance_block(pop, "residual", phenos)
stored <- SELECT id_ind, phenotype_name, pheno_number,
                 residual_value, residual_condition_level
            FROM ind_phenotype
           WHERE phenotype_name IN block AND id_ind IN <planned ids>
             AND residual_value IS NOT NULL
fixed  <- user_residual, mapped onto planned records
check condition_change_action agreement across the block              (D6)
per entity: stratum from plan; compare with stored residual_condition_level (D2)
    'error'       → stop
    'independent' → drop only the stored coordinates whose stratum differs,
                    warn with count + <=5 example IDs
resid  <- adapter → resolve_correlated_draws(...)  per pattern
liability / type conversion / record assembly, in memory              (as today, no writes)

── Stage 3: COMMIT — writes, no RNG ───────────────────────────────────────────
BEGIN TRANSACTION;  on.exit(ROLLBACK unless committed)
  register + INSERT new phenotype_random_effects rows
  register + INSERT ind_phenotype rows incl. residual_value, residual_condition_level
COMMIT
```

Five consequences worth calling out:

**Stage 1 is a real extraction, not a reordering.** The sex filter, repeatable
guard, covariate skip, and both TBV exclusions currently interleave with writes
inside one loop. Pulling them into a plan object is the bulk of the work in
Phase 4 (§8), and it is what makes the ~250-line estimate from v2 optimistic.
The payoff is that Stage 2 and Stage 3 become short and testable in isolation.

**Same-call differing subsets are resolved before anything is written.** If A is
planned for individuals 1–100 and B for 51–100, then 1–50 get a marginal A draw,
51–100 get a joint A/B draw, and no B row is created for 1–50. The `all_equal`
restriction at [add_phenotype.R:579-584](../R/add_phenotype.R#L579-L584) is
**deleted**.

**Stable ordering before every RNG-consuming step.** Sort blocks, patterns,
coordinates, and named-effect levels. DuckDB does not guarantee row order
(CLAUDE.md design principle 7), and this path consumes RNG, so unsorted iteration
would make seeded output depend on physical row order. Residual entity keys
already arrive sorted from `resolve_subset_ids()`; that is not an excuse to skip
the sort on the other four.

**All writes are register + `INSERT`, never `dbWriteTable()`.**
`dbWriteTable()` advances R's RNG by a fixed amount through its random temp-name
generation. Today `add_phenotype()` calls it at
[add_phenotype.R:565](../R/add_phenotype.R#L565) and
[add_phenotype.R:869](../R/add_phenotype.R#L869), interleaved with draws. Under
this plan every draw happens in Stage 2 and every write in Stage 3, and Stage 3
is RNG-neutral, so the stream a call consumes is a function of the model and the
plan only. Follow the idiom at
[define_chromosome.R:219-227](../R/define_chromosome.R#L219-L227). The same rule
applies to the marginal named-effect path at
[phenotype_helpers.R:247](../R/phenotype_helpers.R#L247), which this plan absorbs.

**One transaction covers named-effect draws and phenotype records together.**
A failed record write must not leave a pen draw on disk that conditions the next
call. What the transaction does **not** cover — the RNG state — is D7.

Three more, added in v3.1:

**Planned ids never appear in SQL text.** The Stage-2 pseudocode's
`WHERE id_ind IN <planned ids>` is shorthand. The plan's `(id_ind, pheno_number)`
list is registered with `duckdb_register()` as a temporary view and the stored
lookup is a `JOIN` against it — the same discipline as the genome-effects
evaluator, where individual identifiers never enter the statement. The current
code pastes id lists at [add_phenotype.R:605](../R/add_phenotype.R#L605); that
goes with it. `duckdb_register()` is RNG-neutral.

**Stratum lookup contract.** `condition_table` defaults to `ind_meta` and may be
any table with an `id_ind` column. Stage 1 reads `(id_ind, <condition_column>)`
from it for the planned ids by the same registered-view join and requires
**exactly one row per planned `id_ind`**; zero or several rows is an error
naming the table, the column, and up to 5 example ids (the current code at
[add_phenotype.R:601-612](../R/add_phenotype.R#L601-L612) silently takes the
first match). A `NULL` condition value, or a value matching no stratum, is "no
matching stratum" and follows D2's fallback rule: unconditional `R` if one
exists, stored as `residual_condition_level = NULL`; otherwise an error.

**The resolver draws with base R's RNG.** `stats::rnorm()` on a vector of
standard normals, multiplied through a Cholesky (PD) or eigen (PSD) factor of
the conditional covariance. Not `MASS::mvrnorm()`: it has its own
eigen-decomposition and sign conventions that the RNG-accounting test in §7
would have to replay. `MASS` stays in `Imports` for `define_additive_effects()`.

### 5.6 Named-effect block validation

Before treating two coordinates as members of the same named-effect vector,
validate that their `phenotype_effects` rows agree on
`(source_column, source_table)`. A pen identifier must not be paired with a herd
identifier merely because both effects were given the same `effect_name` string.

Also validate that every coordinate in a block of **two or more** coordinates
uses `distribution = "normal"` — error otherwise, until an explicit copula or
multivariate non-Gaussian API exists. This fixes the silent gamma→normal
substitution.

*(v3.1)* A **1 × 1 block is exempt.** `define_effect_random(distribution =
"gamma" | "uniform")` on a phenotype that shares its `effect_name` with no other
phenotype is a valid, supported model today, drawn by the marginal sampler at
[phenotype_helpers.R:230-239](../R/phenotype_helpers.R#L230-L239). Under this
plan the named-effect adapter dispatches: a 1 × 1 block whose coordinate is
non-normal keeps that marginal sampler (moved, not rewritten, and still
persisted per level in the Stage-3 transaction); every other block goes to the
resolver. The block-size check is what turns a legal gamma singleton into an
error the moment `define_effect_cov_matrix()` tries to join it to a second
phenotype — at the writer, with a message naming the non-normal coordinate.

**Where these checks run** *(clarified in v3; v2 placed them inconsistently)*:
in **both** writers and again in `add_phenotype()`.

- `define_effect_cov_matrix(effect_name, ...)` checks the `phenotype_effects`
  rows that already exist for `(effect_name, each phenotype)`.
- `define_effect_random(phenotype_name, effect_name, ...)` checks whether
  `effect_name` already has a block in `phenotype_var_comp` containing
  `phenotype_name` and, if so, that the new row is normal and source-compatible
  with the block's existing members.
- `add_phenotype()` re-checks in Stage 2 as a backstop, because the two
  `define_*` calls can arrive in either order and nothing prevents a third call
  from changing a row in between.

### 5.7 Repeated records

Narrow, explicit v1 contract:

- Residual covariance applies **across distinct phenotype names** only.
- Prior residual lookup matches on the **same `pheno_number`** across those
  phenotypes (ordinal pairing).
- Repeated records of the *same* phenotype have **independent** residuals.
- Persistent within-animal covariance is the job of a permanent-environment named
  effect, `define_effect_random(..., source_column = "id_ind")` — which already
  exists and is already documented
  ([define_effect_random.R:42-47](../R/define_effect_random.R#L42-L47)).

Conflating p.e. with residual covariance would double-count.

**Documentation must state that `pheno_number` supplies ordinal pairing, not
simulated-time matching**, and code must not assume `pheno_number` will later
become the time coordinate. Irregular longitudinal data needs an explicit
occasion/time model (Layer 2).

### 5.8 Worked example — a pen effect across two time points

Most of this plan is written in terms of residuals, which makes it easy to miss
that **named random effects go through the identical machinery with a different
entity key.** Setup: a pen effect on both average daily gain and backfat,
negatively correlated.

```r
pop <- pop |>
  define_effect_random("ADG", "pen", source_column = "pen_id", variance = 150) |>
  define_effect_random("BF",  "pen", source_column = "pen_id", variance = 4) |>
  define_effect_cov_matrix("pen",
    matrix(c(150, -12, -12, 4), 2, 2,
           dimnames = list(c("ADG", "BF"), c("ADG", "BF"))))
```

> **A pen draw is persistent.** (Codex C4.) The entity is `("pen", "P1")`, so
> the realized `(ADG, BF)` pair for `P1` applies to *every* animal that is ever
> in `P1`, across batches and seasons, in every call, forever. That is the
> model the user declared by using `pen_id` as the level. A pen effect that is
> re-realized per batch is a **different level** — `pen_batch_id`, or
> `interaction(pen_id, batch_id)` written to `ind_meta` — not a different
> feature. "Effects at different times" can mean *completing one persistent
> pen vector* (this plan) or *a new time-indexed pen realization* (Layer 2, not
> this plan). The `define_effect_random()` docs must say this next to the
> example.

#### What happens today

**Day 0 — `add_phenotype("ADG")`.** §7.5 is gated on `length(phenos) >= 2`
([add_phenotype.R:465](../R/add_phenotype.R#L465)), so a single-phenotype call
skips the correlated path entirely. The marginal path
([phenotype_helpers.R:208-250](../R/phenotype_helpers.R#L208-L250)) draws
`rnorm(sd = sqrt(150))` for each pen level and stores
`(ADG, pen, P1, +8.3)`, `(ADG, pen, P2, -4.1)`, …

**Day 100 — `add_phenotype("BF")`.** Same gate, skipped again. No `(BF, pen, P1)`
row exists, so it draws `rnorm(sd = sqrt(4))` — **independent of the `+8.3`
already on disk.** The `-12` is never used.

If instead both phenotypes are named in one call, §7.5 runs and hits **Defect 3**:
it draws a fresh joint `(ADG, BF)` pair for P1, writes the BF half, and discards
the ADG half in favour of the stored `+8.3`.

#### What happens under this plan

**Day 0.** `find_covariance_block(pop, "pen", "ADG")` returns `{ADG, BF}`. For pen
`P1`: stored `{}`, sample `{ADG}`, latent `{BF}`. Zero observed coordinates ⇒
marginal draw from `R[ADG, ADG] = 150`. Store `(ADG, pen, P1, +8.3)`. **`BF`
stays latent — not drawn, not stored**, per §5.3.

**Day 100.** Block discovery returns `{ADG, BF}` again, this time from `BF`. For
`P1`: stored `{ADG: +8.3}`, sample `{BF}`, latent `{}`.

```text
mean = R_BF,ADG · R_ADG,ADG⁻¹ · e_ADG  = -12 · (1/150) · 8.3   = -0.664
var  = R_BF,BF - R_BF,ADG · R_ADG,ADG⁻¹ · R_ADG,BF
     = 4 - (-12)(1/150)(-12)                                   =  3.04

BF_P1 ~ N(-0.664, 3.04)
```

Store `(BF, pen, P1, …)`. The realized pair `(8.3, BF_P1)` now has the correct
joint distribution across a 100-day gap and two separate R calls.

**Day 200 — `add_phenotype("ADG")` for a new batch of animals in `P1`.** Stored
`{ADG, BF}`, sample `{}`. Nothing is drawn; the stored `+8.3` is reused exactly.
See the persistence note above.

#### Mixed patterns

If at day 100 pens `P1`–`P5` carry a stored ADG draw and `P6`–`P10` are brand
new, there are two groups:

| Pattern | Pens | Draw |
|---|---|---|
| stored `{ADG}`, sample `{BF}` | P1–P5 | conditional, `N(-0.08·e_ADG, 3.04)` |
| stored `{}`, sample `{BF}` | P6–P10 | marginal, `N(0, 4)` |

Two conditional-covariance computations, one per pattern.

#### How named effects differ from residuals

| | Residual | Named effect (pen, litter, herd, p.e.) |
|---|---|---|
| Entity | `(id_ind, pheno_number)` | `(effect_name, level)` |
| Storage | `ind_phenotype.residual_value` | `phenotype_random_effects.draw_value` |
| Reuse | never — one draw per record | always — one draw per level, reused by every animal in it |
| Declared by | `define_residual_cov()` | `define_effect_cov_matrix(effect_name, …)` |
| Diagonal-only writer | `define_phenotype(residual_var = )` | `define_effect_random(variance = )` |
| Fixed coordinates | `user_residual` | none in v1 |
| D1 block rule | applies, per stratum | applies |
| D3 realization lock | applies | applies |
| D5 diagonal-writer rule | applies | applies |
| D2 / D6 `condition_change_action` | applies | **does not** — `phenotype_var_comp.condition_column` is residual-only |
| Distribution | Gaussian by construction | must be `"normal"` in a correlated block (§5.6) |
| Extra validation | none | `(source_column, source_table)` must agree (§5.6) |

**Permanent-environment effects get all of this for free.** A p.e. effect is a
named effect with `source_column = "id_ind"`, so its entity is `("pe", <id_ind>)`.

### 5.9 The three `phenotype_var_comp` writers

All validation in §5.4 is implemented once, in an internal
`validate_phenotype_cov_block(conn, effect_name, phenotype_names, stratum)`,
and called by all three writers **inside their write transaction, before
`COMMIT`** — the same shape as `validate_chr_inheritance()` and
`validate_genome_effects()`. None of the three writers currently opens a
transaction ([define_residual_cov.R:103](../R/define_residual_cov.R#L103),
[define_effect_cov_matrix.R:154](../R/define_effect_cov_matrix.R#L154),
[define_effect_cov_matrix.R:325](../R/define_effect_cov_matrix.R#L325) are all
bare `DELETE` then `INSERT`); each gets one, so a rejected block rolls back
instead of leaving a half-written matrix.

`write_phenotype_var_diag()` is not a user-facing function, but its two callers
are. Their behaviour under D1/D3/D5 is spelled out in §6.

---

## 6. Decisions

D1–D4 were decided by the author at v2. D1 is **sharpened** in v3. D3 is
**reopened** by Codex's v2 review with a recommendation to reverse it. D5–D7
are new; D5 records the recommendation the author asked for, D6 and D7 record
recommendations on questions Codex raised. **D3, D6 and D7 were decided by the
author on 2026-09-20, accepting each recommendation as written.**

### D1 — A covariance block is declared in one call, as a complete matrix

**Decision: strict error.** *(v3: sharpened and stated as an algorithm.)*

For a writer call over coordinate set `N` in stratum `s`:

1. Find every existing block in stratum `s` (and, for the residual effect, in
   every other stratum of the same block — see below) that intersects `N`.
2. Let `U = N ∪` every coordinate of those touched blocks.
3. **Require `N == U`**; otherwise error, naming the omitted coordinates and
   pairs, and roll back.
4. Validate the complete new `N × N` matrix (symmetry, finiteness, PSD).
5. Check the realization lock (D3).
6. Delete every row of the touched blocks in stratum `s`; insert the new rows.
   One transaction.

So this fails:

```r
define_residual_cov(c("A", "B"), R_ab)
define_residual_cov(c("B", "C"), R_bc)   # error: block {A,B,C} missing Cov(A,C)
```

this fails — the case v2's wording let through:

```r
define_residual_cov(c("A", "B", "C"), R_abc)
define_residual_cov(c("A", "B"), R_ab)   # error: {A,B} names a strict subset of block {A,B,C}
define_residual_cov("A", matrix(50))     # error: same rule, N = {A}
```

and this is the supported form:

```r
define_residual_cov(c("A", "B", "C"), R_abc)   # complete 3x3, zeros written explicitly
```

To make two phenotypes uncorrelated after declaring them together, redeclare the
block with an explicit `0` — distributionally identical to splitting it, and it
keeps the "declared independence is always something the user typed" property.
There is no operation that removes a phenotype from a block, and none is needed.

**Strata (residual only).** Completeness is checked per stratum. In addition:

- every stratum of a block must name the **same** phenotype set — declaring
  `{A, B}` unconditionally and then `{A}` alone for `condition_level = "M"` is an
  error, because the `M` stratum would have to be over `{A, B}` to be usable;
- a block has at most **one** `(condition_table, condition_column)`; only
  `condition_level` varies across its strata. A second condition column on the
  same block is an error;
- `define_residual_cov("BW", ..., condition_level = "M")` followed by the `"F"`
  stratum is fine — each is a complete 1 × 1 block over `{BW}`.

**Rationale.** The rule states in one sentence: *a residual covariance block is
declared in one call, as a complete matrix.* It matches the domain — breeders get
full `R` matrices out of REML, not pairwise fragments. Most importantly, an
explicit `0` in the matrix is now the **only** way to say "in the same block but
uncorrelated," so declared independence is always something the user typed, never
something the package inferred.

**Consequences, all simplifying:**

1. Completeness and PSD validation live entirely in the writers, never in the
   sampler.
2. `load_phenotype_cov()` returning `NULL` now has exactly **one** meaning —
   *no rows exist*.
3. `resolve_correlated_draws()` may assume a complete, validated, PSD matrix as a
   precondition.

### D2 — Condition-stratum change: error by default, with a stored opt-out

**Decision: error, plus an explicit opt-out.** *(v3: `'independent'` semantics
and stratum recording made precise per Codex B6b.)*

Erroring is the mathematically honest default — if A was drawn under `R_M` and
B is being drawn under `R_F`, nothing defines `Cov(A|M, B|F)`. But a hard error
would make heterogeneous residual variance by a *moving* condition (farm,
management group) unrunnable rather than approximate, which is too blunt.

**Where the opt-out lives.** It is **not** an `add_phenotype()` argument. It is
a column on `phenotype_meta`, set by `define_phenotype()`:

```sql
condition_change_action VARCHAR DEFAULT 'error'   -- 'error' | 'independent'
```

This mirrors `phenotype_meta.missing_component_action` exactly — same table, same
shape, same defaulting, same "one unified field, never per-type arguments"
design. Configuration belongs in a table, not in an R argument, so
`restore_pop()` is complete without re-supplying options, and a per-call argument
would let the same phenotype behave differently on different calls. How
disagreement between block members is handled is D6.

**Precise semantics of `'independent'`** in a block of any size: for each
entity, compare every stored coordinate's `residual_condition_level` with the
stratum the current record resolves to. Stored coordinates from a **different**
stratum are dropped from the observed set for this entity; stored coordinates
from the **same** stratum (and all fixed coordinates) continue to condition the
draw. The sample coordinates are then drawn from the current stratum's
conditional distribution. Warn once with a count of affected entities and up to
5 example IDs, naming the dropped coordinates — matching how
`missing_component_action = "skip"` reports.

**Stratum selection and fallback.** For each planned record, the adapter looks
up the entity's condition value and selects the matching stratum. If no stratum
matches and an unconditional stratum exists, the unconditional `R` is used and
`residual_condition_level` is stored as `NULL` — the *selected* stratum, not the
raw column value. If no stratum matches and there is no unconditional stratum,
**error** naming the unmatched level and count. The current code's
"residuals set to 0 (no unconditional fallback)" at
[add_phenotype.R:640-644](../R/add_phenotype.R#L640-L644) is not an acceptable
fallback and is deleted.

Sex as a condition column never triggers the change path, because an entity's
sex does not change between calls.

*(D2 presupposes that per-stratum `R` is applied at all, which Defect 4 shows is
not currently the case for single-phenotype calls. Under §5.5 every residual
draw — one phenotype or ten — selects its stratum by the entity's condition
level, so D2 becomes meaningful for the first time.)*

### D3 — Covariance redefinition is rejected once draws exist — **decided**

**v2 decision: error by default, `force = TRUE` to override.**
**v3 decision (author, 2026-09-20): error, no override in v1. `force` is removed.**

Codex's v2 review shows the override is not merely lossy but **incoherent**
(B3): suppose `A` was drawn under `R1`, the block is forced to `R2`, and `B` is
later requested. The sampler conditions the stored `A` using `R2`. The resulting
`(A, B)` pair has **no declared joint distribution** — it is not "an old record
that is distributionally inconsistent," it is a *new* record drawn from a
distribution nobody specified. Nothing on disk distinguishes the two eras:
named-effect draws have no covariance-version field, and `residual_condition_level`
records a stratum, not a covariance definition, so the sampler cannot detect or
isolate them. v2's own argument for D3 — "`force` means 'I accept the
inconsistency,' never 'repair it'" — assumed the inconsistency stayed in the
past. It does not; it propagates into every future conditional draw.

**Why not version the covariance instead.** The alternative that keeps `force`
is to stamp every realization with a covariance-definition id and forbid
cross-version conditioning. That is a new column on two tables, a new id
registry, and a new class of "cannot condition across versions" errors — all to
support a workflow (redefine `R` mid-simulation while keeping old records) that
produces a database whose records were generated under two different models.
Nothing on the roadmap needs it.

**The escape hatch already exists and is honest.** `remove_rows()` supports
single-table deletion on every table except `_schema_meta`. The D3 error names
it:

```text
Residual covariance block {on_test_wt, off_test_wt} has 1,240 realized draws in
ind_phenotype. A block cannot be redefined after realization. To redefine it,
remove the realizations first:

  pop |> get_table("ind_phenotype") |>
    filter(phenotype_name %in% c("on_test_wt", "off_test_wt")) |>
    remove_rows()

then call define_residual_cov() again.
```

For a named-effect block the message points at `phenotype_random_effects`
filtered by `effect_name`, and adds one sentence: phenotype records already
computed with those draws are **not** removed by that call and will no longer
be explainable from stored state — remove them too if the population must stay
coherent. Deleting realizations is a user's explicit, visible act on named rows;
that is the right shape for the sharp knife, not a flag on a `define_*` call.

**The lock predicate** *(v3.1)*. A residual block is *realized* when any row of
`ind_phenotype` with `phenotype_name` in the block has `residual_value IS NOT
NULL`. Rows written by `user_values` or `derived_formula` carry `NULL` and do
not lock the block — nothing was drawn under `R`, so nothing is made incoherent
by changing it. The `remove_rows()` recipe in the error text filters the same
way (`!is.na(residual_value)`), so it removes exactly the rows that hold the
lock. A named-effect block is realized when any `phenotype_random_effects` row
exists for `(effect_name, phenotype ∈ block)`. Because the residual predicate
reads a column that Phase 1 (schema) introduces, the schema phase precedes the
writer phase in §8.

**Replicates.** `archive_replicate()` moves and resets both realization tables,
so at the start of every replicate no block is locked. Redefining `R` between
replicates — the one legitimate mid-simulation redefinition workflow — needs no
`remove_rows()` call at all.

`force` appears nowhere in this plan's API. D5's messages point at
`remove_rows()` rather than at a `force` argument.

### D4 — Column names

**Decision: `residual_value` and `residual_condition_level`.**

Both follow naming rule 5 (no abbreviations). The `residual_` prefix on the second
column is kept rather than reusing the bare `condition_level` from
`phenotype_var_comp`: `ind_phenotype` is a wide mixed table, and a bare
`condition_level` would not say which effect's condition it records.

### D5 — The diagonal-only writers after a block exists *(recommendation)*

**Question.** `define_phenotype(residual_var = )` and
`define_effect_random(variance = )` both write a single diagonal cell through
`write_phenotype_var_diag()`. What should they do when the phenotype is already a
member of a multi-phenotype block, or when draws already exist?

**Recommendation: they are the D1 algorithm with `N = {phenotype}`. Error,
pointing at the matrix writer.**

| Situation | `define_phenotype(residual_var = )` / `define_effect_random(variance = )` |
|---|---|
| Phenotype in no block for this effect | write the 1 × 1 block (today's behaviour) |
| Phenotype is a **singleton** block, no draws exist | overwrite the 1 × 1 block (today's behaviour) |
| Phenotype is a singleton block, **draws exist** | **error** (D3). Message names the realizations and the `remove_rows()` call that clears them |
| Phenotype is in a **multi-member** block | **error** (D1: `N ≠ U`). Message: *"'BW' is in residual covariance block {BW, WW}; a block is redeclared as a whole with `define_residual_cov(c("BW", "WW"), R)`"* — regardless of whether draws exist |

**Why error rather than merge.** A diagonal-only rewrite of a joint block is
never what the user meant. Changing `Var(BW)` while holding `Cov(BW, WW)` fixed
changes the correlation implicitly, can break PSD, and does so from a call whose
signature (`define_phenotype`) gives no hint that a joint matrix is involved.
Making the user restate the block costs one line and produces a database whose
`R` was typed in full by someone who saw all of it.

**No `force` anywhere** (D3). `define_phenotype()` keeps `overwrite`, which
governs the `phenotype_meta` row only. On an overwrite call that supplies
`residual_var`, the residual write goes through the table above independently;
an overwrite with no `residual_var` leaves `phenotype_var_comp` untouched, as
today.

**Implementation.** `write_phenotype_var_diag()` becomes a thin call into the
same `validate_phenotype_cov_block()` that the matrix writers use, with a 1 × 1
matrix, so the four rows above are not a fourth code path.

### D6 — `condition_change_action` must agree across a residual block — **decided: option 1**

**Question (Codex B6a).** D2 stores the action per phenotype, but the decision
is taken for a block: when `B` is drawn conditionally on a stored `A` from a
different stratum, whose setting applies?

**Options.**

1. **Require agreement** — every phenotype in a residual block carries the same
   `condition_change_action`; a mismatch is an error naming the phenotypes.
2. **Strictest wins** — `'error'` beats `'independent'`; documented, no check.
3. **Move the setting to the block** — an argument of `define_residual_cov()`,
   stored on the `phenotype_var_comp` rows.

**Decision (author, 2026-09-20): option 1, require agreement.** It keeps D2's deliberate
placement (per phenotype, mirroring `missing_component_action`), and it fails
loudly at the right moment with a message that says what to fix. Option 2
produces the same *outcome* for the forgetful user (an error at sampling) but
with a message about stratum mixing rather than about the real cause, which is
the mismatch. Option 3 is schema-honest about *what* the setting governs but
moves a `define_phenotype()` argument onto a different function and would make
the same value repeat on every row of the block — a different agreement
invariant, not a smaller one.

**Where the check runs.** At sampling (Stage 2), always. At definition time
whenever the metadata exist: `define_residual_cov()` checks the
`phenotype_meta` rows of its block members that are already defined;
`define_phenotype()` checks whether the phenotype is already in a multi-member
block whose other members disagree. Because a block can be declared before its
phenotypes, the sampling-time check is the one that cannot be skipped.

### D7 — RNG state on failure — **decided: option 1**

**Question (Codex T).** Stage 3's transaction makes the *database* atomic. It
does not make the *operation* atomic: an R error after Stage 2 — including a
failed write — leaves `.Random.seed` advanced while the database rolls back.

**Options.**

1. **Database atomic; RNG advances on failure.** Documented and tested. A
   caller who retries after an error gets different draws.
2. **Capture `.Random.seed` before Stage 2 and restore it on any error before
   `COMMIT`.** Retrying after an error reproduces the draws.

**Decision (author, 2026-09-20): option 1.** Nothing in `R/` touches `.Random.seed` today;
every other RNG-consuming function in the package (`add_founders()`,
`add_offspring()`, `define_additive_effects()`, `define_founder_haplotypes()`)
advances the stream on failure, as does every base R function. Making
`add_phenotype()` alone restore the seed would be a global side effect that
holds for one function, interacts with `withr::with_seed()` and `dqrng` streams
in ways that would need their own tests, and buys retry-reproducibility for a
case — an error *after* all validation has passed — that should be rare enough
that the honest answer is "fix the cause and rerun from the last checkpoint."
If the package ever adopts seed restoration it should do so uniformly, as a
package-wide policy, not here.

The contract is stated in the `add_phenotype()` docs and tested directly (§7):
a forced Stage-3 failure leaves the database unchanged **and** `.Random.seed`
advanced by exactly the Stage-2 draws.

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
   phenotypes are in one block **because the pair row exists**, draw jointly,
   and come out uncorrelated.
10. **D1**: declaring `{A,B}` then `{B,C}` errors, names the missing `(A,C)` pair,
    and leaves `phenotype_var_comp` unchanged (rollback).
11. **D1**: declaring the complete `{A,B,C}` in one call succeeds.
12. **D1**: a non-PSD or asymmetric matrix is rejected at `define_*` time, not in
    the sampler.
13. **D1**: `define_residual_cov(c("A","B"), ...)` when the block is `{A,B,C}`
    errors and names `C`.
14. **D1**: `define_residual_cov("A", 1x1)` when `A ∈ {A, B}` errors and names
    the block.
15. **D1**: an unconditional `{A, B}` block followed by an `{A}`-only
    conditional stratum errors; `{A, B}` strata for each level succeed; a second
    `condition_column` on the same block errors.

### Defect 4 — heterogeneous residuals on single-phenotype calls

1. **Before the fix**: the existing composite test is extended to assert
   `var(resid | sex == "M") ≈ 400` and `≈ 800` for `F`; it must fail on current
   code, documenting the defect.
2. After the fix: same assertions pass on a single-phenotype call.
3. Conditional strata with **no** unconditional stratum: single-phenotype call
   succeeds (today it errors "No residual variance found").
4. Individuals whose level matches no stratum fall back to the unconditional `R`
   and store `residual_condition_level = NULL`; with no unconditional `R` they
   **error** rather than receiving `0`.

### Record planning (Stage 1)

1. An individual excluded by `null_class_action = "skip"` consumes **no**
   residual RNG and leaves no stochastic state: seeded output is identical
   whether or not the excluded animal is in the input subset.
2. Same for formula-TBV and composite-TBV exclusions.
3. `pheno_number` assigned in Stage 1 equals what is written in Stage 3.

### Fixed coordinates (`user_residual`)

1. A supplied `A` residual conditions a generated `B` residual **in the same
   call**; across many animals the realized pairs reproduce `R`.
2. Both supplied: nothing drawn (RNG-neutral), both stored, a later `C`
   conditions on them.
3. A named list naming only `A` out of `c("A", "B")` is accepted; `B` is
   sampled.
4. A supplied residual whose length does not match the *planned* record list
   errors with the planned count in the message.
5. A supplied value off the support of a singular `R` (e.g. nonzero for a
   zero-variance coordinate) errors in the resolver.

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
4. Reusing `P1` in a later batch reuses the same persistent draw; `P1:batch2`
   is a new level with a fresh draw.
5. Incompatible `(source_column, source_table)` within a block is rejected —
   by `define_effect_cov_matrix()`, by `define_effect_random()`, and by
   `add_phenotype()` when the rows were changed between the two.
6. Non-normal distribution in a block of two or more is rejected at all three
   sites.
7. *(v3.1)* A **1 × 1** `gamma` or `uniform` named effect still draws from that
   distribution (assert the sign/shape of the realized draws), is persisted per
   level, and is reused on the next call; joining it to a second phenotype with
   `define_effect_cov_matrix()` errors naming the coordinate.

### Numerical

1. One-dimensional conditional and unconditional draws.
2. Perfect and near-perfect correlation; zero conditional variance.
3. Tiny negative eigenvalues absorbed within tolerance; materially negative
   rejected.
4. **Support consistency**: an observed vector materially outside the support of
   a singular `R_oo` errors; one inside it (within scale-aware tolerance)
   proceeds.
5. Tolerance scales with the matrix: the same relative perturbation is accepted
   at variance `1e-6` and `1e6`.
6. Empty entity sets and fully-resolved requests are **RNG-neutral**.
7. Extreme but valid variance scales stay stable.

### Alternate input paths

1. `user_values` and `derived_formula` leave `residual_value` NULL and do not
   constrain later draws.
2. **D2**: a changed stratum errors under the default
   `condition_change_action = 'error'`.
3. **D2**: under `'independent'`, in a **three-coordinate** block, the stored
   coordinate from a different stratum is dropped while the compatible one still
   conditions the draw; the warning names the dropped coordinate; marginal
   variances are still correct.
4. **D2**: an immutable condition column (sex) never triggers either path.
5. **D6**: different `condition_change_action` values within one block are
   rejected at sampling; and at `define_residual_cov()` when the metadata
   already exist.
6. Threshold traits condition on latent liabilities, not observed categories.
7. `store_liability` and `cat_names` populate the now-base columns; for a call
   without `...`, no `ALTER TABLE` is issued (assert on `dbListFields()` before
   and after). A call **with** `...` still adds the user column, inside the
   Stage-3 transaction, and a forced Stage-3 failure rolls the column back too.
8. *(v3.1)* `condition_table` other than `ind_meta`: a planned `id_ind` with two
   rows in it errors naming the table and ids; a `NULL` condition value falls
   to the unconditional `R` and stores `residual_condition_level = NULL`.

### Diagonal writers (D5)

1. `define_phenotype(residual_var = )` on a phenotype in no block writes a
   singleton block.
2. Same call when the phenotype is in `{A, B}` errors and names the block; the
   error text names `define_residual_cov()`.
3. Same call on a singleton block with realized draws errors (D3) and the text
   contains the `remove_rows()` recipe.
4. `define_effect_random(variance = )` mirrors 1–3 against
   `define_effect_cov_matrix()`.
5. `define_phenotype(overwrite = TRUE)` without `residual_var` leaves
   `phenotype_var_comp` untouched.

### Reproducibility and integrity

1. Same seed + same call sequence → identical output.
2. Joint and sequential paths reproduce the target covariance statistically.
3. Forced Stage-3 write failure rolls back residuals, named draws, and phenotype
   rows together — **and** (D7) leaves `.Random.seed` advanced by exactly the
   Stage-2 draws. Both halves asserted.
4. Results are independent of database row order — shuffle `ind_meta` physical
   order between two seeded runs and compare.
5. **RNG accounting**: a seeded `add_phenotype()` call advances `.Random.seed`
   by exactly the draws it makes; no write in the call touches the RNG. Assert by
   comparing the post-call seed to one obtained by replaying the draws alone.
6. **D3**: redefining a block errors once any realized draw exists; **there is
   no override**; the message contains the `remove_rows()` recipe.
7. **D3**: after `remove_rows()` clears the realizations, redefinition succeeds.
8. Writer rollback: a rejected `define_residual_cov()` leaves the prior rows in
   place (today a rejection after the `DELETE` would lose them).

---

## 8. Implementation sequence

*(Re-sequenced in v3 per Codex S: covariance definitions first, record planning
extracted as its own phase.)*

**Phase 0 — remaining legacy code in the phenotype layer.** ✅ **Shipped
2026-09-20** — see `sample_correlated_effects_phase_0.md`. The v2 Phase 0 list
landed in v0.63.x–v0.64.0; the six stragglers listed in v3, plus five more of
the same two kinds found while removing them, are gone:

| What was removed | Kind |
|---|---|
| "Backward-compat fallback" second `get_phenotype_var()` lookup in `add_phenotype()`'s independent residual branch | dead branch — both lookups read the same unconditional diagonal |
| `dbListTables()` existence guards in `get_residual_cov()`, `load_phenotype_cov()`, `load_trait_cov()`, `get_phenotype_var()`, `get_trait_var()`, `add_phenotype()` (`phenotype_components`), `delete_existing_effect()` (`phenotype_random_effects`) | dead guard — `open_pop()` creates every one of these tables unconditionally |
| `ensure_phenotype_var_comp()` (3 call sites), `ensure_trait_var_comp()` (2 call sites), and their `man/` pages | duplicate DDL that could drift |
| The `phenotype_meta`, `phenotype_components`, `phenotype_var_comp` entries of `ensure_trait_tables()`'s DDL list — byte-identical copies of `open_pop()`'s | duplicate DDL that could drift |
| Warning text naming `phenotype_residual_cov` | stale string → `phenotype_var_comp` |
| `get_residual_cov(subset_df = )` and the argument `add_phenotype()` passed to it | unused parameter |

Net −182 lines; no behaviour change; full suite green. The phenotype-layer
DDL now has exactly one home per table: `open_pop.R` for `trait_var_comp`,
`phenotype_meta`, `phenotype_components`, `phenotype_var_comp`;
`ensure_trait_tables()` for `trait_meta`, `phenotype_effects`,
`phenotype_random_effects`, `ind_phenotype`, `ind_tbv`, `ind_tgv`, `ind_ebv`,
`ind_index`, `ind_true_index`. Phase 1 edits the `ind_phenotype` DDL in the
second and the `phenotype_meta` DDL in the first.

**Phase 1 — schema.** ✅ **Shipped 2026-09-20** — see
`sample_correlated_effects_phase_1.md`. *(v3.1: swapped with the writer phase —
the residual realization lock reads `residual_value`, so the column must exist
first.)* Base `ind_phenotype` has `liability_value`, `cat_name`,
`residual_value`, `residual_condition_level`; the two `ALTER TABLE` blocks in
`add_phenotype()` are gone and the two paths write the base columns directly;
`phenotype_meta` has `condition_change_action` and `define_phenotype()` the
argument that sets it (D2); `TABLE_RESERVED_COLS` and the `_schema_meta`
column descriptions updated; CLAUDE.md schema tables updated. New
`test-phenotype_schema.R` (21 expectations) covers the base column set, the
no-`ALTER TABLE` write, NULL-ness of the two residual columns before Phase 5,
reservation, descriptions, and the new argument. Full suite green.

**Phase 2 — centralize covariance definitions.** `validate_phenotype_cov_block()`
implementing the D1 algorithm (pair-row block discovery, `N == U`, per-stratum
completeness, one condition column per block, PSD) and the D3 realization lock
with the predicate defined in D3; the three writers wrapped in transactions and
calling it (§5.9); D5 for the two diagonal writers; distribution (blocks ≥ 2
only) and source checks in `define_effect_cov_matrix()` and
`define_effect_random()`; D6 agreement check at definition time.

**On "tests before the fix".** *(v3.1)* Defects 1–3 cannot be committed as
failing tests, and a characterization test that asserts the *wrong* behaviour
would have to be deleted two phases later. The practical form: each defect's
test is written in the phase that fixes it (Defect 3 and the Phase 6 named-effect
work; Defects 1, 2 and 4 with Phase 5), and the implementer runs the new test
file once against the pre-change commit to confirm it fails there — a manual
check noted in the commit message, not a committed artefact. The one exception
is the Defect 4 assertion strengthening on the existing composite test, which
lands with Phase 5 for the same reason.

**Phase 3 — pure resolver.** `find_covariance_block()` and
`resolve_correlated_draws()` in a new `R/correlated_draws.R`, including
support-consistency validation and the scale-aware tolerance. Test **without
any database access**, across every coordinate pattern and numerical edge case,
including RNG-neutrality of the nothing-to-draw cases.

**Phase 4 — extract record planning (Stage 1).** Pull sex expression, the
repeatable guard, covariate skip, formula/composite exclusion, path
classification, `pheno_number` assignment, stratum lookup, and named-effect
level collection out of the per-phenotype loop into an in-memory plan, with no
RNG and no writes. Behaviour-preserving for the non-correlated path; verified
by the existing suite plus the Stage-1 tests. This is the largest single step
and should land as its own commit.

**Phase 5 — residual integration (Stage 2 + 3 for residuals).** Residual
adapter over the plan: stratum per entity, stored + fixed observed set, D2/D6
checks, per-pattern resolver calls; delete the `all_equal` restriction, §8.5,
and the zero-residual fallback; rewrite `get_residual_cov()` around strata; one
transaction, register + `INSERT` only.

**Phase 6 — named-effect integration.** Replace §7.5 and the normal branch of
the marginal path in `compute_covariate_contribution()` with the named-effect
adapter over the same resolver, persistent per-level entity identity, source and
distribution checks as the `add_phenotype()` backstop. The gamma/uniform branch
of the marginal path is kept for 1 × 1 blocks (§5.6) and its `dbWriteTable()`
write moves into the Stage-3 transaction with everything else. Named-effect and
phenotype writes in the same transaction.

**Phase 7 — transaction/RNG boundary.** Implement D7 exactly as decided; the
two-part integrity test (database unchanged, seed advanced by the Stage-2 draws).

**Phase 8 — documentation, housekeeping, performance.**
- Roxygen: sequential sampling; ordinal repeated-record pairing; liability
  scale; D1 whole-block rule; D2 `condition_change_action` and `'independent'`
  semantics; D3 lock and the `remove_rows()` recipe; D5 on `define_phenotype()`
  and `define_effect_random()`; D6 agreement; D7 RNG contract; the persistent
  pen-identity note on `define_effect_random()`; a culling example on
  `add_phenotype()`; the `user_residual` subset-list contract.
- CLAUDE.md: `ind_phenotype` and `phenotype_meta` schema tables; the D1/D5
  rules in the `define_phenotype()` / `define_residual_cov()` sections; the
  `add_phenotype()` description of the three-stage flow and residual model;
  the `phenotype_var_comp` note that `condition_column` rows form strata of one
  block.
- `R/schema.R` column descriptions are added in Phase 1; confirm
  `test-schema-print.R` still passes (table lists are unchanged).
- `NEWS.md` under **0.71.0**, `DESCRIPTION` version bump.
- Benchmark under `dev/benchmarks/` for large populations, optimizing the
  observation-pattern query and batched writes **without changing RNG
  semantics**.

**Size.** v2's ~250 net lines covered Phases 5–6 only. With Phase 4's
extraction, the three writers' transactions and validator, and the resolver
with its numerical guards, a realistic estimate is **~600–800 lines net** across
`add_phenotype.R` (which should shrink), `phenotype_helpers.R`,
`define_effect_cov_matrix.R`, `define_residual_cov.R`, `define_phenotype.R`,
`define_effect_random.R`, and the new `R/correlated_draws.R`, plus tests.
§7.5 and §8.5 are deleted outright.

---

## 9. Known limitations to document in v1

**Deletion and regeneration.** Under Option A, deleting a phenotype record deletes
its residual — referentially clean, but *not* probabilistically coherent: later
coordinates may already have been drawn conditional on the deleted value. Their
marginal distributions remain valid, but the sequential conditional history can no
longer be reconstructed. v1 documents that deleting and resampling correlated
phenotype records is not guaranteed to produce a coherent conditional history.
The same applies, more sharply, to the D3 `remove_rows()` recipe: it is the
user's declaration that the old realizations no longer matter.

**Informative missingness.** The v1 resolver conditions only on **exact stored
deviations**. Culling on an already-observed phenotype is naturally supported.
Dropout driven by an *unobserved* latent value is informative missingness and
requires a joint dropout model; interval censoring carries partial information
about a liability that v1 cannot use.

**Count traits are a transformed Gaussian liability model**, not a Poisson or
negative-binomial model. The docs must say so rather than implying general
correlated count support.

**Persistent named-effect identity.** A level is forever (§5.8). Occasion-specific
realizations require occasion in the level. This is a documentation obligation,
not a limitation of the mechanism.

**Concurrency.** DuckDB is single-writer per file, so cross-session races are
largely moot, but the Stage-3 transaction is still required for crash safety.
Logical identities to enforce: residual `(id_ind, phenotype_name, pheno_number)`;
named effect `(phenotype_name, effect_name, level)`.

---

## 10. Explicitly out of scope (Layers 2–3)

Not implemented and not blocked. The original catalogue of these models, with
the requirements each would carry, is in the v1-round Codex review
(`sample_correlated_effects_codex_review.md`, at commit `ea693e4`; it has since
been removed from the working tree) and is summarized in
`sample_correlated_effects_v2_review.md`.

- **Layer 2 — longitudinal Gaussian**: observation times and shared occasion
  identity (distinct from `pheno_number`), random regression and reaction norms,
  permanent-environment trajectories, time-heterogeneous residual variance,
  covariance functions and AR/state-space residuals, time-indexed group
  exposures and occasion-specific pen realizations.
- **Layer 3 — non-Gaussian and event models**: generalized count traits,
  zero-inflated/hurdle processes, copula or latent-variable mixed-family
  dependence, survival with censoring and competing risks, informative dropout.

The one v1 design constraint that keeps these reachable: **the resolver takes an
opaque entity key and a covariance matrix, is database-independent, and never
assumes the key is an `id_ind` or that a level is time-invariant.** That is the
whole extension-point budget for now — no speculative provider interface, no
reserved columns.
