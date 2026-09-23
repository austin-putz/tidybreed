# Sampling correlated random effects at different points in simulated time

**Status**: **v3.8 — complete.** Phases 0–8 shipped (0–3 on 2026-09-20,
4–7 on 2026-09-21–22, 8 on 2026-09-22; `sample_correlated_effects_phase_0.md`
… `sample_correlated_effects_phase_8.md`). Every decision D1–D8 is decided,
implemented and tested.
File:line references are refreshed after each phase. v3 is a re-baseline
against the codebase as of v0.70.0 (2026-09-20) plus the Codex review of v2
(`sample_correlated_effects_v2_review.md`). Every v2 design decision stands
except one: D3's `force = TRUE` escape hatch, which Codex showed to be
incoherent and which is now removed. **D3, D6 and D7 were decided by the author
on 2026-09-20, accepting the v3 recommendations** (§6). v3.1 adds six
implementation-level clarifications found in the final read-through (marked
*v3.1* inline); none changes a decision.

**What changed from v3.7 to v3.8** (Phase 8, the closing phase):

1. **D8 is new and decided**: seeded output is **bit-identical**, not
   identical-within-tolerance. The v3.7 note that `add_tbv()`/`add_tgv()` were
   not bit-reproducible is resolved rather than documented-as-a-limitation —
   the cause was one parallel `SUM()` in the genome-effect evaluator, and an
   exact accumulator fixes it for free. See D8.
2. **D6 gains a mutability rule.** The agreement requirement made the value
   unchangeable once a block had two members; a block-scoped writer,
   `define_condition_change_action()`, is the answer. The v3.5 open question in
   §8 is closed.
3. The Phase 7 whole-database snapshot test no longer needs a tolerance on
   `ind_tbv`; it compares every table with `expect_identical()`.
4. §8 Phase 8 is marked shipped, with its list.

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

**What changed from v3.1 to v3.2** (Phase 2 implementation, no decision
changes; details in §5.9, D1, D5, D6):

- The validator runs **before** the `DELETE`, inside the writer's transaction —
  D1's steps 1–5 on the current table state, then step 6 — rather than as a
  post-write table check. Same guarantee (a rejected call changes nothing),
  better messages (the omitted coordinates are known before anything moves).
- D1 strata: **growing** a block that already has two or more strata
  (`{A, B}` unconditional + `{A, B}` for `M`, then `{A, B, C}`) is an error in
  every call order, because each stratum is checked against the others. The
  message gives the `remove_rows()` call on `phenotype_var_comp` that clears the
  block so every stratum can be redeclared. Blocks are declared once, from REML
  output; this is not a workflow worth a special case.
- D5: `define_phenotype(overwrite = TRUE)` **used to delete** the phenotype's
  unconditional residual rows (half of a multi-phenotype block's pair rows)
  whether or not `residual_var` was supplied — the plan's "as today" was wrong.
  It now leaves `phenotype_var_comp` untouched; only `residual_var` writes.
- D5: `define_effect_random(variance = )` **used to ignore** a supplied
  `variance` whenever a stored value existed. It now always writes through the
  block writer (singleton: overwrite; multi-member block: D1 error). Its D3 row
  is reached through `overwrite = TRUE`, which discards that phenotype's stored
  draws for the effect before the variance is rewritten — documented on the
  argument.
- `write_phenotype_var_diag()` is deleted rather than thinned; its only caller,
  `define_effect_random()`, calls the block writer directly inside its own
  transaction (which now also covers the `phenotype_effects` row and the
  overwrite deletes).
- §5.6 at the writers also refuses to join a **non-random** effect (a
  `fixed_class`/`fixed_cov` row sharing the `effect_name`) into a block.
- Unconditional rows store `condition_table = NULL` (previously `'ind_meta'`
  from `define_residual_cov()` and `NULL` from `define_effect_cov_matrix()`).
- All three writers are RNG-neutral (register + `INSERT`; `define_residual_cov()`
  used `dbWriteTable()`).

**What changed from v3.2 to v3.3** (Phase 3 implementation, no decision
changes; details in §5.2 and §5.4):

- The block loader is `find_covariance_blocks()` — **plural**, taking a
  connection — because the phenotypes of one `add_phenotype()` call can fall
  into several components. It returns one entry per component touching the
  targets, each with every stratum already assembled as a matrix
  (`unconditional`, `conditional[[level]]`); phenotypes with no rows are simply
  absent. It re-checks the D1 invariants on load (rows can be removed by hand
  with `remove_rows()`) and errors with the redeclaration call rather than
  guessing around a hole.
- The **pattern grouping lives inside the resolver**, not the adapter: the
  adapter passes one `observed` matrix per `(stratum, sample set)` with `NA`
  for "not observed", and the resolver groups by non-`NA` pattern itself.
  Adapter steps 6–7 in §5.4 collapse to "one call per `(stratum, sample set)`". The RNG contract is what makes this safe:
  a call consumes exactly `n × m` normals in entity order, whatever the
  patterns, so grouping is invisible to the stream.
- The resolver runs **every check before the first `rnorm()`** — a rejected
  call leaves `.Random.seed` untouched — and a coordinate with zero conditional
  variance still consumes its normal (multiplied by zero), so the accounting
  test in §7 is `n × m`, not a data-dependent count.
- Tolerance is **relative**: `tolerance × λ_max`, default
  `nrow(R) × sqrt(.Machine$double.eps)`. The PD/PSD decision for `R_oo` and
  for the conditional covariance uses it; Cholesky on the PD side, eigen with
  a **sign-normalized** `V` on the PSD side, so a seeded draw through the
  singular path does not depend on the LAPACK build.
- The writer now stores `(M + t(M)) / 2` (bit-identical for an already
  symmetric input), so stored `(i, j)` and `(j, i)` rows are exactly equal and
  the loader can require it.

**What changed from v3.3 to v3.4** (Phase 4 implementation, no decision
changes; details in §5.5 and `sample_correlated_effects_phase_4.md`):

- `add_phenotype()` now **is** the three-stage flow: `.ap_plan()` (Stage 1),
  `.ap_resolve()` (Stage 2), `.ap_commit()` (Stage 3) in the new
  `R/add_phenotype_stages.R`. The pre-Phase-4 §7.5 pre-draw, §8.5 joint
  residual draw, marginal named-effect draws and the independent residual
  draw still exist, but as Stage-2 steps over the plan, with their writes
  moved to Stage 3. Phases 5–6 replace their bodies; the stage boundaries
  do not move again.
- **Stage 3 is already one transaction, register + `INSERT` only.** The
  plan had this arriving in Phases 5–7; building the commit stage once was
  simpler than building it twice. Phase 7 is now the D7 integrity test plus
  whatever D7 says about the RNG on failure, not a transaction rewrite.
- **Random-effect levels are collected per planned record**, so a level
  touched only by an excluded individual is never drawn (the pre-Phase-4
  marginal path drew it). This is the §7 "no stochastic state" property,
  and it is what the pre-draw and marginal paths consume in Stage 2.
- **Stable ordering is now real**: the input subset is sorted by `id_ind`
  before planning; random-effect levels and `phenotype_effects` rows are
  sorted before any draw. Seeded output no longer depends on physical row
  order. `user_values` / `user_residual` positional matching is over this
  order, which the roxygen now states.
- The **stratum lookup contract** (§5.5, v3.1) is implemented in Stage 1
  for every planned record — exactly one row per planned id in the
  condition table, registered-view join, error otherwise. Until Phase 5 the
  value only feeds the retained §8.5 path; single-phenotype calls still
  ignore strata (Defect 4 stays open until Phase 5, as planned).
- The joint residual path (§8.5, retained) compares **planned** id sets,
  i.e. after exclusions, rather than pre-exclusion subsets. Two phenotypes
  whose exclusions differ now draw independently for that call instead of
  jointly over a superset that includes non-records. Phase 5 removes the
  equal-set restriction entirely, so this is a transient narrowing on a
  path that is being deleted, accepted for the RNG property above.
- `next_pheno_numbers()` and `.eval_derived_formula()` no longer write a
  temp table / paste ids into SQL; both use a registered view. A derived
  phenotype can read a feeder phenotype **planned in the same call** through
  the in-memory records, which is what lets Stage 3 write everything at
  once.
- Review-pass tightenings: named `user_values` must name planned
  individuals, each once (unknown ids used to be written as records);
  effects declared on a `derived_formula` phenotype are ignored rather than
  evaluated (and drawn) by accident; sorts before draws use byte order so
  seeded output is locale-independent; `sample_residuals()` is now
  `(n, R)` only, pending its Phase 5 deletion.

**What changed from v3.4 to v3.5** (Phase 5 implementation; details in
§5.3, §5.5, D2 and `sample_correlated_effects_phase_5.md`):

- The **residual adapter** exists: `.ap_resolve_residuals()` →
  `.ap_residual_block()` in `R/add_phenotype_stages.R`, over
  `find_covariance_blocks()` and `resolve_correlated_draws()`. The §8.5
  joint draw, its equal-planned-sets restriction, the independent residual
  branch, `sample_residuals()` and `get_residual_cov()` are deleted.
  Defects 1, 2 and 4 are closed; `residual_value` /
  `residual_condition_level` are written for every model-path record.
- **Stage-2 RNG order is now the plan's**: every named-effect draw (the
  pre-Phase-6 joint pre-draw, then the marginal draws for new levels, in
  plan order) precedes the residual adapter. Before Phase 5 the marginal
  named-effect draws were interleaved per phenotype with the independent
  residual draw. Phase 6 replaces the named-effect *bodies* without moving
  this boundary, so the order will not change again.
- **Stratum fallback reporting** (D2, refined): a record whose condition
  value is `NULL` falls back to the unconditional `R` silently — `NULL` is
  the documented "no group" state. A **non-`NULL` value matching no
  stratum** also falls back, but with a warning (count, the unmatched
  levels, ≤ 5 ids), since an unmodelled level is more likely a typo than a
  design. Without an unconditional stratum both are the D2 error, naming
  the levels and the count.
- **`user_residual` contract** (§5.3) as shipped: a plain vector is
  accepted only when exactly one phenotype of the call is on the model
  path (the pre-Phase-5 code applied one vector to every phenotype); a
  list must be named by `phenotype_name`, may name a subset, and may not
  name a `derived_formula` or `user_values` phenotype; named
  (per-`id_ind`) vectors and `user_values` + `user_residual` together are
  errors. Fixed values are stored with the stratum the record resolved to.
  A phenotype whose residuals are all supplied needs no residual block
  (nothing is drawn or conditioned for it); the check runs before any
  draw on every call.
- **Prevalence thresholds need an unconditional variance.** The
  categorical `prevalence` cut-point is `mean + z·√(V_A + V_E)` with `V_E`
  the *marginal* residual variance. When a phenotype has only conditional
  strata there is no single `V_E`; the pre-Phase-5 code silently used `0`.
  It is now an error pointing at `thresholds =` or an unconditional
  stratum. (Before Phase 5 such a phenotype could not be sampled at all —
  Defect 4 — so nothing that worked stops working.)
- **The D3 lock is now live** for residual blocks: the first model-path
  `add_phenotype()` call on a block member locks the block. The Phase 2
  tests that simulated a realized residual by `UPDATE` now use real
  records; a `user_values` record carries no residual and does not lock.
- **Open item found (not fixed here, for Phase 8):** the definition-time D6
  check refuses to flip `condition_change_action` on any member of a
  defined block, in either direction, so the value is immutable once every
  member is defined (`test-phenotype_cov_block.R` pins "flipping one member
  back is refused"). Changing a block from `'error'` to `'independent'`
  currently requires setting it on the phenotypes *before* the block is
  declared, or dropping and redefining the phenotypes. A same-call
  "set on all members" path is the obvious fix; noted in §8 Phase 8.

**What changed from v3.5 to v3.6** (Phase 6 implementation; details in
§5.6, §5.8 and `sample_correlated_effects_phase_6.md`):

- The **named-effect adapter** exists: `.ap_resolve_named_effects()` →
  `.ap_named_effect_block()` in `R/add_phenotype_stages.R`, the residual
  adapter's block loop with entity `(effect_name, level)`, storage
  `phenotype_random_effects`, no strata, no fixed coordinates, no D2. The
  §7.5 joint pre-draw (`.ap_predraw_named_effects()`, the last
  `MASS::mvrnorm()` on the phenotype path), the marginal path
  (`.ap_resolve_random_terms()`), `.ap_existing_draws()` and
  `load_phenotype_cov()` are deleted. **Defect 3 is closed**: a level with
  a stored `ADG` draw gets its `BF` draw conditional on it, whether the two
  phenotypes arrive in one call or a season apart, and a stored draw is
  never redrawn.
- **Stage-1 loads the named-effect blocks once**, one
  `find_covariance_blocks()` call per effect (byte-sorted), targeting only
  the phenotypes whose planned records carry a random term for it
  (`.ap_named_effect_targets()`). A block member with no such term is a
  latent coordinate: not drawn, but its stored draws condition.
- **Stage-2 RNG order, final**: effects in byte-sorted `effect_name` order
  → blocks in loader order → one resolver call per sample-set group in
  sorted group order (`n × m` normals, entity = level order) → then the
  residual adapter. A pen touched only by a `NULL`-level or excluded record
  is still never drawn (Phase 4's property, kept).
- **§5.6 backstop** is in: `validate_named_effect_block()` runs in Stage 2
  on every block of two or more, so a `phenotype_effects` row edited
  between the two `define_*` calls is refused at `add_phenotype()` before
  any draw. A named-effect block found with conditional strata (hand
  edit) is an integrity error, as is a random term whose variance rows are
  gone ("No variance stored for random effect").
- **1 × 1 `gamma` / `uniform` blocks** keep their marginal sampler (moved
  into the adapter, not rewritten); a 1 × 1 `normal` block goes through the
  resolver, which for one coordinate and nothing observed *is*
  `sqrt(v) · z` — the Phase 4 RNG-accounting test holds unchanged.
- **`load_phenotype_cov()` retired**: its one remaining caller,
  `write_renum_par()` (BLUPF90), now assembles the residual matrix from
  `find_covariance_blocks()` (`.blupf90_residual_cov()`): block-diagonal
  with explicit zeros between independent blocks, an error for a trait
  with only conditional strata (the old reader ignored `condition_column`
  and could return a conditional row's value).
- Docs: `define_effect_random()` carries the §5.8 persistence note ("a
  level's draw is persistent"; occasion goes in the level) and the
  conditional-sampling paragraph; `define_effect_cov_matrix()` and
  `add_phenotype()` describe the sequential named-effect draw.

**What changed from v3.6 to v3.7** (Phase 7 implementation, no decision
changes and **no code change**; details in D7, §7 and
`sample_correlated_effects_phase_7.md`):

- **D7 is stated and tested as decided.** Nothing in `R/` touches
  `.Random.seed` (verified); the Stage-3 transaction from Phase 4 is the
  whole of the database half. `?add_phenotype` and `?add_phenotype_stages`
  now carry the contract: any error — a Stage-1 rejection, a Stage-2 error
  after some draws, a failed write — leaves `ind_phenotype` and
  `phenotype_random_effects` exactly as they were, and `.Random.seed`
  advanced by exactly the draws made before the error; a retry draws
  different values unless it re-seeds.
- **The database half is stronger than "no rows"**: a failed call also rolls
  back the `ALTER TABLE ADD COLUMN` that `prepare_extra_cols()` may have
  issued earlier in the same transaction for a new `...` column (DuckDB's
  DDL is transactional). Both Rd pages say so. The one write a failed call
  *does* leave is the Stage-1 `add_tbv()` upsert, which costs no RNG and is
  rewritten by the retry; that is stated rather than hidden.
- **The integrity test is a whole-database comparison**, not a row count:
  `tests/testthat/test-add_phenotype_failure_contract.R` snapshots every
  base table before the failing call and asserts identity after, for a
  Stage-3 write failure (§7 item 3), a Stage-2 error in the named-effect
  adapter after the preceding effect drew, a Stage-2 error in the residual
  adapter (D2) after the named effect and the preceding residual block
  drew — proving a block resolved before the failing one is never written
  on its own — and three rejections before any draw (Stage 1; `user_residual`
  length; the resolver's support check). Each asserts the seed against a
  replay of the draws made.
- **One observation outside this plan's scope**: `ind_tbv` is the one table
  a failed call still touches (the Stage-1 `add_tbv()` upsert), and its
  values are *not* bit-reproducible across identical runs — DuckDB's
  multi-threaded `SUM()` in the genome-effects evaluator changes the
  floating-point summation order (differences ≈ 5e-16 on a 60-locus
  trait; identical with `SET threads = 1`). The test compares `ind_tbv`
  with tolerance and every other table for identity. This is the
  evaluator's concern (`R/genome_effects_eval.R`), not the phenotype
  layer's; it is recorded in §8 Phase 8 as an open item because it bears on
  the "same seed reproduces within the current code" contract at the last
  bit of `pheno_value`.

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

*(File:line references in §0–§2 describe the code as it was when the defects
were found, before Phase 4 restructured `add_phenotype()`; they are kept as
the record of what was wrong. §5 onward is refreshed.)*

v2 was written against ~v0.60. Nine releases later:

| Change | Release | Effect on this plan |
|---|---|---|
| Migration section of `ensure_trait_tables()` and `.migrate_var_comp_tables()` deleted | v0.64.0 (`3fa68a8`) | v2 Phase 0 is done; see §8 for the stragglers that remain |
| Dead `dbExistsTable`/`has_meta` guards removed across `define_effect_*`, `schema()`, `add_phenotype()` | v0.63.x (`d882556`) | Same |
| `trait_effects` → `phenotype_effects`; `trait_random_effects` → `phenotype_random_effects` | v0.64.0 | Every reference in this plan renamed |
| `TABLE_RESERVED_COLS` gained entries for `phenotype_random_effects`, `phenotype_components`, `founder_haplotypes` | v0.64.0 | New `ind_phenotype` columns must be added to the reserved list (§8) |
| RNG discipline: `dbWriteTable()` advances R's RNG by a fixed amount (random temp-name generation); `duckdb_register()` + `INSERT` is RNG-neutral. Documented in `define_genome.R`, `founder_haplotype_helpers.R`, `add_offspring.R`, [define_effect_cov_matrix.R:129](../R/define_effect_cov_matrix.R#L129) | v0.5x–0.6x | New hard requirement in §5.5; the three `phenotype_var_comp` writers comply since Phase 2 |
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

*(Found in the Phase 2 review, 2026-09-20.)* The same joint path crashes when
the multi-phenotype call has **zero** individuals left after the skips (e.g.
re-calling `add_phenotype(c("A", "B"))` on non-repeatable phenotypes everyone
already has): `sample_residuals()`
([phenotype_helpers.R:62](../R/phenotype_helpers.R#L62)) calls
`MASS::mvrnorm(n = 0, …)`, which fails with "non-conformable arguments" instead
of writing nothing. The single-phenotype path is fine (`rnorm(0)` is empty).
Not patched — Phase 5 deletes `sample_residuals()`, and the resolver's
empty-entity-set contract (§7 Numerical 6: RNG-neutral no-op) covers it.

### Defect 3 — correlated *named* random effects lose correlation on partial re-draw — ✅ closed (Phase 6, 2026-09-21)

*(v3.6: fixed by the named-effect adapter; `test-add_phenotype_named_effects.R`
"Defect 3 closed" pins the exact conditional draw and the untouched stored
value. The description below is of the pre-Phase-6 code.)*

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
([define_effect_random.R:31](../R/define_effect_random.R#L31)) and the marginal
path honours it, but [add_phenotype.R:539](../R/add_phenotype.R#L539) calls
`MASS::mvrnorm()` unconditionally. A gamma effect inside a covariance block is
silently drawn normal.

### Observation (new in v3) — a third writer bypasses the block

`phenotype_var_comp` has **three** writers, not two:

| Writer | Called by | What it writes |
|---|---|---|
| `define_residual_cov()` | user | full `n × n` residual block, one condition slice |
| `define_effect_cov_matrix()` | user | full `n × n` block for a named effect |
| `write_phenotype_var_diag()` (deleted in Phase 2) | `define_phenotype(residual_var = )`, `define_effect_random(variance = )` | **one diagonal cell**, in place, off-diagonals untouched |

Called after a block exists, the third writer can turn a valid `R` into a
non-PSD one, and it can do so after draws have been realized. v2's D1 and D3
named only the first two writers. Closed in §6 (D1 sharpened, D5 added).
*(Phase 2)* All three paths now go through one writer,
[phenotype_cov_block.R](../R/phenotype_cov_block.R).

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

### 5.2 Covariance block discovery — ✅ shipped (Phase 3, 2026-09-20)

The single most important thing missing from v1. When the user calls
`add_phenotype("off_test_wt")`, the requested vector contains one name; the
sampler must still find `on_test_wt`.

```r
find_covariance_blocks(conn, effect_name, phenotype_names)
# → list of blocks, one per connected component touching phenotype_names:
#   list(effect_name, phenotypes (sorted), condition_table, condition_column,
#        unconditional = R or NULL, conditional = list(<level> = R, ...))
```

*(v3.3: plural, and a connection rather than `pop`, because one call's
phenotypes may span several components and the adapters run on `pop$db_conn`.
Implemented in [correlated_draws.R](../R/correlated_draws.R) over
`.pvc_block_members()` from Phase 2; the loader re-checks completeness, one
condition column, exact symmetry and finiteness per stratum, since
`remove_rows()` can remove pair rows by hand, and errors with the
redeclaration call. Phenotypes with no rows are absent from the result:
`setdiff(targets, unlist(lapply(blocks, `[[`, "phenotypes")))`.)*

Each block is the **connected component** of the covariance graph containing
its targets. **The graph is defined by the existence of a stored pair row, never by
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
as it is today ("No residual variance found"). A phenotype absent from
`find_covariance_blocks()`'s result therefore has exactly one meaning: no
rows exist. *(v3.6: `find_covariance_blocks()` is now the only reader of
`phenotype_var_comp` matrices — `load_phenotype_cov()` is gone.)*

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
now with a precise definition because planning precedes drawing. *(v3.5,
shipped in `.ap_fixed_residuals()`: "single phenotype" means exactly one
phenotype of the call on the model path; the list may not name a
`derived_formula` / `user_values` phenotype; per-`id_ind` names are
refused; combining with `user_values` is an error.)*

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
*(v3.3: the adapter groups by stratum and sample set; the observed-pattern
grouping is done inside the resolver from the `NA` pattern of `observed`.)*

This state model is what makes all of these fall out of one mechanism: A and B
together on identical sets; A on more individuals than B; B later on a culled
subset; B on a mixture of individuals with and without stored A; a supplied A
with a generated B; partial named-effect vectors; blocks of more than two
phenotypes; **and a size-1 block with heterogeneous variance (Defect 4)** — a
single sample coordinate, zero observed, drawn from the stratum matching the
entity's condition level.

*(v3.5: every residual case in that list is a test in
`test-add_phenotype_residuals.R`, most of them exact replays of the
resolver's `n × m` stream. v3.6: the named-effect cases — partial
vectors, mixed stored/new levels, three-coordinate blocks — are the same
kind of test in `test-add_phenotype_named_effects.R`.)*

### 5.4 The resolver — ✅ shipped (Phase 3, 2026-09-20)

Database-independent, opaque entity keys, no writes, no knowledge of strata
(the adapter has already picked one matrix per stratum):

```r
resolve_correlated_draws <- function(covariance,          # one stratum's R: validated, complete, dimnamed
                                     sample_coordinates,  # chr; coordinates to draw (same set for every entity)
                                     entity_keys,         # opaque; vector/list or data frame; sorted by caller
                                     observed = NULL,     # entity x coordinate; stored + fixed; NA = not observed
                                     tolerance = NULL)    # relative; default nrow(R) * sqrt(eps)
# returns entity x sample_coordinates, in the caller's entity order
```

*(v3.3, as shipped in [correlated_draws.R](../R/correlated_draws.R).)* The
contract the adapters and the §7 tests rely on:

- **Every check precedes the first random number.** Validation, the PSD check
  on `R`, the support check on every entity, and the conditional-covariance
  factorization for every pattern all run first; a rejected call leaves
  `.Random.seed` untouched.
- **A successful call consumes exactly `n × m` standard normals** (`n`
  entities, `m` sample coordinates) from `stats::rnorm()`, in entity order and
  `sample_coordinates` order within an entity — whatever the observed
  patterns, and even for a coordinate with zero conditional variance (its
  normal is multiplied by zero and it is returned at its conditional mean
  exactly). Zero entities or zero sample coordinates consume nothing. The
  accounting test replays `rnorm(n * m)`.
- **Patterns are grouped inside.** Entities are grouped by the non-`NA`
  columns of `observed`; coefficients and the factor are computed once per
  pattern; the pre-drawn `z` rows are applied per pattern, so grouping cannot
  reorder the stream.
- **Tolerance is relative**: absolute `tolerance × λ_max(R)`. `R_oo` is
  inverted by Cholesky when its smallest eigenvalue exceeds it, else by an
  eigen pseudoinverse with the support check
  `‖V_null' e_o‖ ≤ tolerance × max(‖e_o‖, √λ_max)`. The conditional covariance
  is symmetrized, eigenvalues in `[−tol, 0)` are zeroed, anything below is an
  error; it is factored by Cholesky when PD and by `V √D` (eigenvector signs
  normalized so the largest-magnitude component is positive) otherwise.

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

1. Discover the complete block by pair-row existence
   (`find_covariance_blocks()`).
2. Select the exact stratum for every planned entity.
3. Produce stable entity and coordinate order.
4. Load stored values and check their stratum compatibility (D2).
5. Add fixed current values (`user_residual`).
6. Group entities by `(stratum, sample coordinates)` and build one `observed`
   matrix per group with `NA` for "not observed" *(v3.3: the observed-pattern
   split is the resolver's job)*.
7. Call the resolver once per group.
8. Merge fixed and sampled values into the planned rows.
9. Persist only inside the outer transaction.

Validation of `R` itself lives **upstream** in the three `phenotype_var_comp`
writers (§5.9): matching unique dimnames, symmetry within tolerance, finite
entries, non-negative diagonal, PSD within tolerance, **block replacement rule
(D1)**, and **no realized draws for the block (D3)**. Sampling code should never
discover a malformed matrix deep inside phenotype generation, and must never
interpret one as an instruction to draw independently.

### 5.5 `add_phenotype()` control flow — three stages — ✅ shipped in full (Phases 4–6, 2026-09-21)

*(v3.4: `.ap_plan()` / `.ap_resolve()` / `.ap_commit()` in
`R/add_phenotype_stages.R` implement this structure. Stage 1 is complete as
specified. Stage 3 is the single register + `INSERT` transaction. v3.5: the
residual part of Stage 2 is the adapter `.ap_resolve_residuals()` /
`.ap_residual_block()`. v3.6: the named-effect part is
`.ap_resolve_named_effects()` / `.ap_named_effect_block()`; nothing of the
pre-Phase-4 draw code remains.)*

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
    block  <- find_covariance_blocks(conn, effect, phenos)     # one per component
    validate: every coordinate normal; compatible (source_column, source_table)
    stored <- SELECT level, phenotype_name, draw_value FROM phenotype_random_effects
    draws  <- adapter → resolve_correlated_draws(...)  per (stratum, sample set)

blocks <- find_covariance_blocks(conn, "residual", phenos)   # one per component
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
resid  <- adapter → resolve_correlated_draws(...)  per (block, stratum, sample set)
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
restriction is **deleted** *(v3.5: done — `.ap_residual_block()` groups
entities by `(stratum, sample set)`, so 1–50 and 51–100 are two resolver
calls on the same block, in sorted group order)*.

**Stable ordering before every RNG-consuming step.** Sort blocks, patterns,
coordinates, and named-effect levels. DuckDB does not guarantee row order
(CLAUDE.md design principle 7), and this path consumes RNG, so unsorted iteration
would make seeded output depend on physical row order. Residual entity keys
already arrive sorted from `resolve_subset_ids()`; that is not an excuse to skip
the sort on the other four.

**All writes are register + `INSERT`, never `dbWriteTable()`.**
`dbWriteTable()` advances R's RNG by a fixed amount through its random temp-name
generation. Before Phase 4, `add_phenotype()` called it for named-effect draws
and for records, interleaved with draws. Under
this plan every draw happens in Stage 2 and every write in Stage 3, and Stage 3
is RNG-neutral, so the stream a call consumes is a function of the model and the
plan only. Follow the idiom at
[define_chromosome.R:219-227](../R/define_chromosome.R#L219-L227). *(v3.4:
done — `.ap_commit()` is that idiom, and the marginal named-effect path's
write now goes through it too.)*

**One transaction covers named-effect draws and phenotype records together.**
A failed record write must not leave a pen draw on disk that conditions the next
call. What the transaction does **not** cover — the RNG state — is D7.
*(v3.7: both halves tested in `test-add_phenotype_failure_contract.R`. The
transaction also covers the **schema**: an `ALTER TABLE ADD COLUMN` issued by
`prepare_extra_cols()` for a new `...` column rolls back with the rows.
The one write outside it is Stage 1's `add_tbv()` upsert, which is
RNG-neutral — verified — and idempotent.)*

Three more, added in v3.1:

**Planned ids never appear in SQL text.** The Stage-2 pseudocode's
`WHERE id_ind IN <planned ids>` is shorthand. The plan's `(id_ind, pheno_number)`
list is registered with `duckdb_register()` as a temporary view and the stored
lookup is a `JOIN` against it — the same discipline as the genome-effects
evaluator, where individual identifiers never enter the statement. *(v3.4:
Stage 1 already does this — `.ap_read_by_id()` is the one join helper for the
repeatable guard, the TBV read, effect source tables and the stratum lookup,
and `next_pheno_numbers()` / `.eval_derived_formula()` register their ids.
*v3.5:* the composite / formula contributor reads went the same way —
`R/contributor_tbv.R` — so nothing on the phenotype path pastes an id.)*
`duckdb_register()` is RNG-neutral.

**Stratum lookup contract.** `condition_table` defaults to `ind_meta` and may be
any table with an `id_ind` column. Stage 1 reads `(id_ind, <condition_column>)`
from it for the planned ids by the same registered-view join and requires
**exactly one row per planned `id_ind`**; zero or several rows is an error
naming the table, the column, and up to 5 example ids (the pre-Phase-4 code
silently took the first match; *v3.4:* `.ap_condition_values()` implements
the contract). A `NULL` condition value, or a value matching no stratum, is "no
matching stratum" and follows D2's fallback rule: unconditional `R` if one
exists, stored as `residual_condition_level = NULL`; otherwise an error.
*(v3.5: the stratum is selected per entity in `.ap_residual_block()`; a
non-`NULL` unmatched value warns on fallback, `NULL` does not — see the
v3.5 header.)*

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
phenotype is a valid, supported model today, drawn by the marginal sampler
(since Phase 4: `.ap_resolve_random_terms()`). Under this
plan the named-effect adapter dispatches: a 1 × 1 block whose coordinate is
non-normal keeps that marginal sampler (moved, not rewritten, and still
persisted per level in the Stage-3 transaction); every other block goes to the
resolver. *(v3.6: as shipped in `.ap_named_effect_block()`.)* The block-size check is what turns a legal gamma singleton into an
error the moment `define_effect_cov_matrix()` tries to join it to a second
phenotype — at the writer, with a message naming the non-normal coordinate.

**Where these checks run** *(clarified in v3; v2 placed them inconsistently)*:
in **both** writers and again in `add_phenotype()`. *(Phase 2)* The writer
sites are shipped: `validate_named_effect_block()` in
[phenotype_cov_block.R](../R/phenotype_cov_block.R), called by
`define_effect_cov_matrix()` (through the block validator) and by
`define_effect_random()` with its pending row. *(Phase 6)* The
`add_phenotype()` backstop is `.ap_named_effect_block()`'s first check,
calling the same function with `caller = "add_phenotype()"` on every block
of two or more, before any draw.

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
  ([define_effect_random.R:52-57](../R/define_effect_random.R#L52-L57)).

Conflating p.e. with residual covariance would double-count.

**Documentation must state that `pheno_number` supplies ordinal pairing, not
simulated-time matching**, and code must not assume `pheno_number` will later
become the time coordinate. Irregular longitudinal data needs an explicit
occasion/time model (Layer 2).

### 5.8 Worked example — a pen effect across two time points — ✅ this is what ships (Phase 6)

*(v3.6: the three days below are `test-add_phenotype_named_effects.R`'s
first test, exact to the resolver's stream; the mixed-pattern table is its
second. "What happens today" describes the pre-Phase-6 code.)*

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
> example. *(v3.6: they do — "A level's draw is persistent".)*

#### What happened before Phase 6

**Day 0 — `add_phenotype("ADG")`.** §7.5 was gated on `length(phenos) >= 2`
(Phases 4–5: `.ap_predraw_named_effects()`), so a single-phenotype call
skipped the correlated path entirely. The marginal path
(`.ap_resolve_random_terms()`, deleted in Phase 6) drew
`rnorm(sd = sqrt(150))` for each pen level and stores
`(ADG, pen, P1, +8.3)`, `(ADG, pen, P2, -4.1)`, …

**Day 100 — `add_phenotype("BF")`.** Same gate, skipped again. No `(BF, pen, P1)`
row exists, so it draws `rnorm(sd = sqrt(4))` — **independent of the `+8.3`
already on disk.** The `-12` is never used.

If instead both phenotypes are named in one call, §7.5 runs and hits **Defect 3**:
it draws a fresh joint `(ADG, BF)` pair for P1, writes the BF half, and discards
the ADG half in favour of the stored `+8.3`.

#### What happens under this plan

**Day 0.** `find_covariance_blocks(conn, "pen", "ADG")` returns the block `{ADG, BF}`. For pen
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

### 5.9 The three `phenotype_var_comp` writers — ✅ shipped (Phase 2, 2026-09-20)

All validation in §5.4 is implemented once, in
`validate_phenotype_cov_block(conn, effect_name, phenotype_names, cov_matrix,
condition_column, condition_table, condition_level, caller)` in
[phenotype_cov_block.R](../R/phenotype_cov_block.R), and every write to the
table goes through `write_phenotype_cov_block()` (own transaction) or
`.pvc_write_block()` (caller's transaction), which call the validator **before
the `DELETE`**. The validator, in order: matrix (dimnames, finite, symmetric,
diagonal ≥ 0, PSD with a scale-aware tolerance); block discovery by pair-row
existence across every stratum (`.pvc_block_members()`, a fixpoint over
`phenotype_name_1 IN (…) OR phenotype_name_2 IN (…)`); `N == U`; for the
residual, the strata rules and D6; for a named effect, §5.6; then the D3 lock.
The write is `DELETE` of the stratum's rows for the block, then
`duckdb_register()` + `INSERT … SELECT` — RNG-neutral, same shape as
`validate_chr_inheritance()` and `validate_genome_effects()`.

The three entry points:

| Entry point | Calls | Transaction |
|---|---|---|
| `define_residual_cov()` ([define_residual_cov.R](../R/define_residual_cov.R)) | `write_phenotype_cov_block("residual", …, stratum)` | its own |
| `define_effect_cov_matrix()` ([define_effect_cov_matrix.R:148-158](../R/define_effect_cov_matrix.R#L148-L158)) | `define_residual_cov()` for `"residual"`; `write_phenotype_cov_block(effect_name, …)` otherwise; the genetic route to `trait_var_comp` is unchanged | its own (phenotype routes) |
| `define_effect_random(variance = )` ([define_effect_random.R](../R/define_effect_random.R)) | `.pvc_write_block()` with a 1 × 1 matrix | one transaction around overwrite-delete, variance, §5.6 check and the `phenotype_effects` row |

`define_phenotype(residual_var = )` runs the validator once **before** the
`phenotype_meta` write (so a D1/D3 rejection leaves the metadata untouched)
and then calls `define_residual_cov()`, which validates again inside its
transaction. `write_phenotype_var_diag()` no longer exists. The
`get_phenotype_var()` reader drops its legacy `condition_column = ''`
predicate; writers never wrote `''`.

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
- *(v3.2)* Once a block has two or more strata it cannot be **grown**: with
  `{A, B}` unconditional and `{A, B}` for `M` on disk, `{A, B, C}` in any
  stratum errors (the other strata would be over a strict subset), and
  `{A, B, C}` cannot be written stratum-by-stratum either, because each write is
  checked against the strata that still name `{A, B}`. The error gives the
  `remove_rows()` call on `phenotype_var_comp` that clears the block; every
  stratum is then redeclared over `{A, B, C}`. A block with a single stratum
  grows freely (`N ⊋ U` passes `N == U` after the closure).

**Rationale.** The rule states in one sentence: *a residual covariance block is
declared in one call, as a complete matrix.* It matches the domain — breeders get
full `R` matrices out of REML, not pairwise fragments. Most importantly, an
explicit `0` in the matrix is now the **only** way to say "in the same block but
uncorrelated," so declared independence is always something the user typed, never
something the package inferred.

**Consequences, all simplifying:**

1. Completeness and PSD validation live entirely in the writers, never in the
   sampler.
2. A phenotype absent from `find_covariance_blocks()` has exactly **one**
   meaning — *no rows exist*.
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
"residuals set to 0 (no unconditional fallback)" is not an acceptable
fallback and is deleted. *(v3.5: the real rule is in `.ap_residual_block()`:
the error names the unmatched levels, the count, the strata stored and up
to 5 ids. The `'independent'` warning and the `'error'` message both list
up to 5 `id (coordinate: stored under X, now Y)` examples and, for
`'independent'`, the dropped coordinates.)*

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
an overwrite with no `residual_var` leaves `phenotype_var_comp` untouched.
*(v3.2 correction: before Phase 2 an overwrite deleted the phenotype's
unconditional residual rows unconditionally, which would have torn half the
pair rows out of a multi-phenotype block. That `DELETE` is gone.)*

**Implementation** *(as shipped in Phase 2)*. Both diagonal writers call the
block writer with a 1 × 1 matrix, so the four rows above are not a fourth code
path; `write_phenotype_var_diag()` is deleted. `define_phenotype()` runs the
validator before its `phenotype_meta` write and then calls
`define_residual_cov()`. `define_effect_random()` always writes a supplied
`variance` (it used to ignore it when a stored value existed) and, with
`variance = NULL`, requires a stored one; with `overwrite = TRUE` it first
discards that phenotype's stored draws for the effect, so the "singleton with
draws" row of the table is reached only through an explicit overwrite, which is
documented on the argument. Every rejection rolls back the whole call,
including the `phenotype_effects` delete an overwrite performs.

### D6 — `condition_change_action` must agree across a residual block — **decided: option 1** ✅ *(v3.8: mutability resolved, Phase 8)*

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
*(Phase 2)* Both definition-time sites are shipped, as
`.check_condition_change_agreement(conn, block, pending, caller)` in
[phenotype_cov_block.R](../R/phenotype_cov_block.R); `define_phenotype()`
passes its pending row so the check runs before the `phenotype_meta` write.
Stage 2 calls the same function without `pending` *(v3.5: shipped — run per
block whenever a stored coordinate exists, before D2 is applied)*.

**Mutability** *(v3.5 raised it; v3.8 decided and shipped it)*. Option 1 has a
consequence the decision did not spell out: agreement makes the value
**immutable** once a block has two or more defined members. Every route to
changing it goes through one member at a time, and each single-member change is
exactly the disagreeing state the check refuses — so the block is frozen at
whatever value its members happened to be registered with, and the error
message's advice ("set the same value on every phenotype") named an impossible
sequence.

The two candidates in §8 were a same-call "set on all members" form of
`define_phenotype()` and a relaxation letting the *last* member's change
through. Neither was taken:

- The relaxation is order-dependent (which member is "last" depends on the
  call history, not on the model) and it leaves the database in the
  disagreeing state between the two calls, where any `add_phenotype()` would
  fail. A mid-sequence invalid state is exactly what D1's whole-block rule
  exists to prevent.
- Overloading `define_phenotype()` looks smaller than it is: `overwrite = TRUE`
  replaces the whole `phenotype_meta` row, so flipping one flag means restating
  `type`, `mean`, `expressed_sex` and the rest, and forgetting one silently
  resets it to a default. A one-column change must not be spelled as a
  re-registration.

**Decision: the writer matches the scope of the property.**
`condition_change_action` is block-scoped, so
`define_condition_change_action(pop, phenotype_name, action)` writes every
member of the named phenotype's residual block in one transaction, touching
that column and nothing else. Agreement then holds by construction on this
path, and the D6 check remains the backstop for the case it was written for — a
`define_residual_cov()` call that *joins* phenotypes which already disagree,
where the user genuinely must choose. A phenotype in no block is a block of
one, so the same call works before a block exists and `define_phenotype()`
keeps setting the value at registration time.

The action is deliberately **not** subject to the D3 realization lock. D3 locks
the covariance *matrix*, because changing it would invalidate draws already
taken from it. The action changes nothing about a stored residual; it decides
how a *future* record treats one whose stratum no longer matches. Changing it
mid-simulation is a legitimate modelling choice — "from here on, animals that
moved farm draw independently rather than stopping the run" — and it leaves
`ind_phenotype` and `phenotype_random_effects` untouched.

### D7 — RNG state on failure — **decided: option 1** ✅ *(v3.7: shipped as decided, Phase 7)*

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
advanced by exactly the Stage-2 draws. *(v3.7: stated in `?add_phenotype`
and `?add_phenotype_stages`; tested in
`test-add_phenotype_failure_contract.R` for Stage-3, both Stage-2 adapters
mid-stream, and three pre-draw rejections, each against a whole-database
snapshot. The one write a failed call leaves is the Stage-1 `add_tbv()`
upsert, which is idempotent and RNG-independent.)* *(v3.8: that upsert is
bit-identical too, so the snapshot test no longer gives `ind_tbv` a
tolerance — see D8.)*

### D8 — "Same seed reproduces" means bit-identical — **decided, shipped (Phase 8)** ✅

**Question** *(raised in v3.7, from a Phase 7 finding)*. The Phase 7 integrity
test could not compare `ind_tbv` with `expect_identical()`: two runs of an
identical population differed by about 1e-15. CLAUDE.md's first reproducibility
contract says a seed "must produce identical output on repeated runs of the
current implementation". Does *identical* mean bit-identical, or identical
within tolerance?

**Cause** *(localized in Phase 8, not guessed)*. Every step of the evaluator's
one statement is bit-stable on its own except the last. `add_red` sums at most
two allele copies per group — floating-point addition **is** commutative; it is
associativity that fails, and with two summands there is no associativity
question. `state` sums integers, `string_agg` carries an explicit `ORDER BY`,
and `product()` runs over a fixed member count. The final
`SUM(genome_value * prod)` over a trait's terms is the only reduction with
enough summands for the order to matter, and DuckDB's parallel hash aggregate
combines partial sums in whatever order the threads finish. Measured: stable
under `SET threads = 1`, and the divergence grows with term count
(8.9e-16 at 60 loci, 2.2e-15 at 2000). The RNG stream was never involved —
residuals and random-effect draws were already bit-identical; `pheno_value`
inherited the wobble through its TBV term.

**Options.**

1. **Tolerance** — say in CLAUDE.md that evaluator-derived values reproduce
   within tolerance, and leave the evaluator alone.
2. **Ordered reduction** — `list_reduce(list(x ORDER BY id_genome_effect), ...)`.
   Deterministic, measured at 1.3× the plain `SUM()`.
3. **Single-threaded** — `SET threads = 1` around the statement.
4. **Exact accumulation** — sum in `DECIMAL`, which is integer arithmetic and
   therefore associative, then cast back.

**Decision: option 4, bit-identical via exact accumulation.**

Option 1 is the weaker contract and was tempting only while the wobble looked
cosmetic. It is not. A breeding simulation is *chaotic downstream of
selection*: once `select_parents()` exists (roadmap), a last-bit difference in
a TBV will occasionally flip a truncation-selection ranking between two nearly
tied animals, and from there the two runs produce different pedigrees, not
different last bits. That is a user-visible reproducibility failure, and it
would be reported as a bug in selection rather than traced back to a `SUM()`.

Option 3 serializes the expensive join along with the cheap aggregate. Option 2
is deterministic but changes the aggregate state from one value per group to
one per row, which regresses the memory profile of the hottest query in the
package and conflicts with design principle 2 ("enables simulations larger
than RAM").

Option 4 has neither cost: measured at 0.96–1.09× the plain `SUM()` across
500–5000 loci (the join dominates), with `O(groups)` state. Scale 18 is finer
than double precision for any contribution a breeding model produces, so the
per-term rounding sits below the inputs' own representation error.

**Consequences.**

- `GEV_ACC_TYPE` in [genome_effects_eval.R](../R/genome_effects_eval.R), with
  the reasoning inline so it is not "optimized" back to a plain `SUM()`.
- The accumulator's domain is finite at both ends. A single term below 1e-18
  rounds to zero before it is added — unreachable for a trait scaled to any
  sane variance (1e6 QTL on a unit-variance trait are still ~1e-3 each), so
  it is a statement about the representable domain rather than a practical
  limit. The ceiling is the end worth guarding: a model whose per-term
  contribution or running total exceeds 1e20 now errors. `.gev_accumulator_error()` rethrows
  DuckDB's bare conversion error as a tidybreed message naming the cause. The
  bound is **not** pre-checked: the only cheap bound
  (`|genome_value| × 2^n_members` summed over terms) is loose enough to reject
  legal high-order models.
- `tbv_value`, `tgv_value` and `pheno_value` move in the last bits relative to
  0.70.0. Pre-1.0.0 this is bookkeeping, and CLAUDE.md forbids testing against
  previous versions' output anyway.
- CLAUDE.md's reproducibility section now says "identical" means bit-identical
  and names the parallel-`SUM()` failure mode, because the next many-summand
  aggregate anyone adds will have the same problem.
- `test-genome-effects-determinism.R` pins it with `expect_identical()` — the
  only assertion that can catch a regression, since `expect_equal()` passes on
  the broken code. 16 of its 18 expectations fail on the pre-change evaluator.

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

*(Phase 2)* `tests/testthat/test-phenotype_cov_block.R` (110 expectations)
covers: Residual blocks 10–15; Named effects 5–6 at the two writer sites (the
`add_phenotype()` site is Phase 6) and the writer half of 7; Diagonal
writers 1–5; Reproducibility and integrity 6–8 (the residual half of 6–7 with a
`residual_value` set by SQL, since Phase 5 is what writes it; 8 by mocking
`next_int_id()` to fail after the `DELETE`); plus RNG-neutrality of the writers,
the strata rules, and D6 at both definition-time sites.

### Residual blocks — ✅ 1–9 covered (Phase 5, `test-add_phenotype_residuals.R`); 10–15 in Phase 2

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
9b. *(v3.8, added in the Phases 0–8 review)* **Unequal, non-unit variances.**
   Items 1–9 were all written against unit-diagonal matrices (`R_AB`, the 3 × 3
   with 1s on the diagonal, and the stratified `R_AB` / `4 * R_AB` pair). In
   every one of those the conditional slope `sigma_AB / var_A` equals the
   correlation `sigma_AB / sqrt(var_A var_B)`, so a resolver conditioning on
   the **correlation** would have passed all of them. The implementation was
   correct, but nothing proved it. Four tests now separate the two: a
   sequential and a joint call on `[[4, 1.8], [1.8, 9]]` (`beta = 0.45` versus
   a correlation of `0.30`, with an explicit assertion that the correlation
   form is *not* what came out), a three-phenotype block with three different
   variances recovered through the sequential path, and two strata carrying
   different unequal-variance matrices.
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

### Defect 4 — heterogeneous residuals on single-phenotype calls — ✅ covered (Phase 5)

*(v3.5: item 1's "must fail on current code" holds by construction rather
than by a re-run — the Phase 4 suite asserted `residual_value` and
`residual_condition_level` were `NULL` on every record, and the strengthened
composite assertions read both. Not a committed artefact, per §8 "On tests
before the fix".)*

1. **Before the fix**: the existing composite test is extended to assert
   `var(resid | sex == "M") ≈ 400` and `≈ 800` for `F`; it must fail on current
   code, documenting the defect.
2. After the fix: same assertions pass on a single-phenotype call.
3. Conditional strata with **no** unconditional stratum: single-phenotype call
   succeeds (today it errors "No residual variance found").
4. Individuals whose level matches no stratum fall back to the unconditional `R`
   and store `residual_condition_level = NULL`; with no unconditional `R` they
   **error** rather than receiving `0`.

### Record planning (Stage 1) — ✅ covered (Phase 4, `test-add_phenotype_stages.R`)

1. An individual excluded by `null_class_action = "skip"` consumes **no**
   residual RNG and leaves no stochastic state: seeded output is identical
   whether or not the excluded animal is in the input subset.
2. Same for formula-TBV and composite-TBV exclusions.
3. `pheno_number` assigned in Stage 1 equals what is written in Stage 3.

*(v3.4 adds: the repeatable guard behaves the same way; records are planned
and written in `id_ind` order; a call consumes exactly `n` residual normals
plus one per new random-effect level and nothing else; `user_values` calls
are RNG-neutral; Stage 3 rolls back both tables on failure; a derived
phenotype reads a feeder planned in the same call; the condition-table
exactly-one-row contract.)*

### Fixed coordinates (`user_residual`) — ✅ covered (Phase 5)

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

### Repeated records — ✅ covered (Phase 5)

1. Distinct phenotypes match on equal `pheno_number`.
2. Repeated records of the same phenotype stay residual-independent.
3. Unequal record counts never condition on the wrong record.

### Named effects — ✅ covered (Phase 6, `test-add_phenotype_named_effects.R`; 5–7 writer halves in Phase 2)

*(v3.6 mapping: 1 → "day 0 / day 100 / day 200" and "Defect 3 closed"
(exact draws; the covariance is exact by construction of the resolver, so
the statistical reproduction is asserted on the conditional coefficients
rather than a Monte-Carlo estimate); 2 → "mixed patterns" and "three-
phenotype block"; 3 → "Defect 3 closed"; 4 → "day 200" and the
permanent-environment test (a new level is a new entity by definition of
the entity key); 5–6 → "§5.6 backstop" (the `add_phenotype()` site, with
rows edited by SQL between the two `define_*` calls); 7 → "1 x 1 gamma or
uniform".)*

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

*(Phase 3)* `tests/testthat/test-correlated_draws.R` (115 expectations, no
database for the resolver half) covers 1–7 below, the mixed-pattern and latent
cases of §5.3, the `n × m` RNG accounting, RNG-neutrality of rejected calls,
the exact-mean return at zero conditional variance, a 20 000-entity
conditional-covariance check with sampling-error tolerances, input
validation, and `find_covariance_blocks()` (components, strata, sorted
members, absent phenotypes, hand-broken blocks, RNG-neutrality).

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

1. Same seed + same call sequence → identical output. ✅ *(v3.8: and
   **identical means bit-identical** — D8. The draws always were; the TBV
   term was not, so `pheno_value` was not either.
   `test-genome-effects-determinism.R` asserts it with `expect_identical()`,
   which is the only assertion that catches this class of regression.)*
2. Joint and sequential paths reproduce the target covariance statistically.
3. Forced Stage-3 write failure rolls back residuals, named draws, and phenotype
   rows together — **and** (D7) leaves `.Random.seed` advanced by exactly the
   Stage-2 draws. Both halves asserted. ✅ *(v3.7:
   `test-add_phenotype_failure_contract.R`, "a failed Stage-3 write";
   also a Stage-2 failure in each adapter after earlier draws, and the
   retry-draws-new-values consequence.)*
4. Results are independent of database row order — shuffle `ind_meta` physical
   order between two seeded runs and compare. ✅ *(v3.7: `test-add_phenotype_stages.R`,
   "seeded output is identical when ind_meta's physical row order is
   reversed" — both adapters, records and draws compared. Phase 4's
   neighbouring test asserts the **output** order only; the shuffle this
   item asks for was missing until the Phase 7 review.)*
5. **RNG accounting**: a seeded `add_phenotype()` call advances `.Random.seed`
   by exactly the draws it makes; no write in the call touches the RNG. Assert by
   comparing the post-call seed to one obtained by replaying the draws alone.
   ✅ *(Phase 4: `test-add_phenotype_stages.R`, "a call consumes exactly its
   draws"; Phases 5–7 add the per-block replays.)*
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

**Phase 2 — centralize covariance definitions.** ✅ **Shipped 2026-09-20** —
see `sample_correlated_effects_phase_2.md`. New `R/phenotype_cov_block.R`:
`validate_phenotype_cov_block()` (D1 pair-row block discovery, `N == U`,
per-stratum completeness, one condition column per block, symmetry/finite/PSD;
D3 lock with the D3 predicate; D6; §5.6) and the single writer
`write_phenotype_cov_block()` / `.pvc_write_block()` (validate → `DELETE` →
register + `INSERT`, one transaction). `define_residual_cov()`,
`define_effect_cov_matrix()` (phenotype routes) and `define_effect_random()`
all write through it; `define_phenotype(residual_var = )` pre-validates and
delegates; `write_phenotype_var_diag()` deleted. First user-visible behaviour
change of the plan (fragments, subsets, non-PSD matrices, locked blocks,
disagreeing `condition_change_action`, and incompatible named-effect rows are
now errors; `define_phenotype(overwrite = TRUE)` no longer deletes residual
rows; `define_effect_random(variance = )` no longer ignores a supplied value).
110 new expectations; full suite green. Details and the v3.2 clarifications in
the header.

**On "tests before the fix".** *(v3.1)* Defects 1–3 cannot be committed as
failing tests, and a characterization test that asserts the *wrong* behaviour
would have to be deleted two phases later. The practical form: each defect's
test is written in the phase that fixes it (Defect 3 and the Phase 6 named-effect
work; Defects 1, 2 and 4 with Phase 5), and the implementer runs the new test
file once against the pre-change commit to confirm it fails there — a manual
check noted in the commit message, not a committed artefact. The one exception
is the Defect 4 assertion strengthening on the existing composite test, which
lands with Phase 5 for the same reason.

**Phase 3 — pure resolver.** ✅ **Shipped 2026-09-20** — see
`sample_correlated_effects_phase_3.md`. New `R/correlated_draws.R`:
`find_covariance_blocks()` (one entry per component touching the targets,
every stratum as a matrix, D1 invariants re-checked on load) and
`resolve_correlated_draws()` (per-pattern conditional moments, Cholesky/eigen
with support consistency and relative tolerance, all checks before the first
`rnorm()`, exactly `n × m` normals consumed). No caller changed yet —
`add_phenotype()` still runs §7.5/§8.5 until Phases 5–6. The Phase 2 writer
now stores the symmetrized matrix. 115 new expectations; full suite green.
Details and the v3.3 clarifications in the header.

**Phase 4 — extract record planning (Stage 1).** ✅ **Shipped 2026-09-21** —
see `sample_correlated_effects_phase_4.md`. `add_phenotype()` is now
`.ap_plan()` → `.ap_resolve()` → `.ap_commit()` (`R/add_phenotype_stages.R`).
Stage 1 plans every record — sex expression, repeatable guard, fixed-effect
skip, formula/composite exclusion, path, `pheno_number`, residual condition
value, random-effect level per record — with no RNG and no writes (the
`add_tbv()` prerequisite aside). Stage 2 holds the pre-Phase-4 draw paths
over the plan; Stage 3 is one register + `INSERT` transaction. Seeded output
is independent of physical row order, and an excluded individual consumes
no RNG and leaves no random-effect draw behind. 57 new expectations in
`test-add_phenotype_stages.R`; full suite green. Phases 5–6 replace the
Stage-2 bodies; the transaction planned for Phases 5–7 already exists.

*(Original scope:)* Pull sex expression, the
repeatable guard, covariate skip, formula/composite exclusion, path
classification, `pheno_number` assignment, stratum lookup, and named-effect
level collection out of the per-phenotype loop into an in-memory plan, with no
RNG and no writes. Behaviour-preserving for the non-correlated path; verified
by the existing suite plus the Stage-1 tests. This is the largest single step
and should land as its own commit.

**Phase 5 — residual integration (Stage 2 + 3 for residuals).** ✅
**Shipped 2026-09-21** — see `sample_correlated_effects_phase_5.md`. The
residual adapter `.ap_resolve_residuals()` / `.ap_residual_block()` over
`find_covariance_blocks()` and `resolve_correlated_draws()`: stratum per
entity with the D2 fallback rule, stored coordinates of every block member
at the same `pheno_number` (registered-view join), `user_residual` as fixed
coordinates, D6 then D2 on the stored set, one resolver call per
`(stratum, sample set)` in sorted order, `residual_value` /
`residual_condition_level` written by Stage 3. `.ap_joint_residuals()`, the
independent-draw branch, `sample_residuals()` and `get_residual_cov()` are
deleted (there was nothing to "rewrite around strata" — the Phase 3 loader
already is that). Named-effect draws now all precede the residual adapter.
93 new expectations in `test-add_phenotype_residuals.R`; three Phase 1–2
tests that pinned the pre-Phase-5 `NULL` residuals updated; full suite
green.

*(Original scope:)* Residual
adapter over the plan: stratum per entity, stored + fixed observed set, D2/D6
checks, per-pattern resolver calls; delete the `all_equal` restriction, §8.5,
and the zero-residual fallback; rewrite `get_residual_cov()` around strata.

**Phase 6 — named-effect integration.** ✅ **Shipped 2026-09-21** — see
`sample_correlated_effects_phase_6.md`. The named-effect adapter
`.ap_resolve_named_effects()` / `.ap_named_effect_block()` over
`find_covariance_blocks()` and `resolve_correlated_draws()`: blocks loaded
once in Stage 1 per effect, entity = level, stored coordinates via a
registered view, §5.6 backstop per block, one resolver call per sample-set
group, 1 × 1 gamma/uniform marginal sampler kept; draws before the
residual adapter. `.ap_predraw_named_effects()`,
`.ap_resolve_random_terms()`, `.ap_existing_draws()`,
`.ap_append_random_effects()` and `load_phenotype_cov()` deleted. 74
expectations in `test-add_phenotype_named_effects.R` (the review pass
found and fixed one crash — every planned level `NULL` — and added its
test); full suite green.

*(Original scope:)* Replace §7.5 and the normal branch of
the marginal path (*v3.4:* `.ap_predraw_named_effects()` and the normal
branch of `.ap_resolve_random_terms()`) with the named-effect
adapter over the same resolver, persistent per-level entity identity, source and
distribution checks as the `add_phenotype()` backstop. The gamma/uniform branch
of the marginal path is kept for 1 × 1 blocks (§5.6); *(v3.4: its write
already goes through the Stage-3 transaction, and Stage 1 already supplies
the level per planned record via `.ap_covariate_terms()`.)* *(v3.5: the
residual adapter is the template — same block loop, entity =
`(effect_name, level)`, stored = `phenotype_random_effects`, no strata, no
fixed coordinates, no D2. Its draws must stay **before** the residual
adapter in `.ap_resolve()`, where the named-effect draws already are.)*

**Phase 7 — transaction/RNG boundary.** ✅ *(v3.7, shipped 2026-09-21; see
`sample_correlated_effects_phase_7.md`.)* Implement D7 exactly as decided; the
two-part integrity test (database unchanged, seed advanced by the Stage-2 draws).
*(v3.4: the Stage-3 transaction already exists since Phase 4 and
`test-add_phenotype_stages.R` has a rollback test; what remains is the
Stage-2-failure half of D7 and the seed-advanced assertion.)* *(v3.7: no
code change was needed — nothing in `R/` touches `.Random.seed`, and Stage
3 is the only writer. Shipped: the D7 paragraphs in `?add_phenotype` and
`?add_phenotype_stages`, the CLAUDE.md failure-contract note, and
`test-add_phenotype_failure_contract.R` — 38 expectations over whole-database
snapshots: Stage-3 write failure, Stage-2 failure in the named-effect
adapter after `herd` drew, Stage-2 failure in the residual adapter (D2)
after the pen level and residual block `{A}` drew, three pre-draw
rejections, the retry / re-seed consequence, the rollback of an
`ALTER TABLE` issued earlier in the failing transaction, and the one write
a failed call does leave (the Stage-1 `add_tbv()` upsert). The review also
added §7 item 4's missing physical-order shuffle to
`test-add_phenotype_stages.R`.)*

**Phase 8 — documentation, housekeeping, performance.** ✅ **Shipped
2026-09-22** — see `sample_correlated_effects_phase_8.md`. The closing phase.
Two decisions, the remaining docs, the benchmark, the release bookkeeping:

- **D8 (new): seeded output is bit-identical.** The v3.7 observation that
  `add_tbv()`/`add_tgv()` were not bit-reproducible is resolved, not
  documented as a limitation. Cause localized to one parallel `SUM()` in the
  genome-effect evaluator; fixed with an exact accumulator (`GEV_ACC_TYPE`)
  that measures free and keeps `O(groups)` state. New
  `.gev_accumulator_error()` for the finite range. New
  `test-genome-effects-determinism.R` (18 expectations; 16 fail on the
  pre-change evaluator). CLAUDE.md's reproducibility section now says
  bit-identical and names the failure mode. The Phase 7 snapshot test drops
  its `ind_tbv` tolerance.
- **D6 mutability (the v3.5 open item): decided and shipped.** New
  `define_condition_change_action()` — block-scoped, one transaction, one
  column. Neither §8 candidate was taken; the reasoning is in D6. The D6 error
  message now names a recipe that works. Four new tests in
  `test-phenotype_cov_block.R` (109 → 125 expectations).
- **Roxygen.** The list below was already complete except the culling example,
  which is now on `?add_phenotype`: record A on everyone, cull on the realized
  A, record B on the survivors conditional on their stored A residual. (The
  `user_residual` subset-list contract landed with Phase 5; D7 with Phase 7;
  the persistent pen-identity note with Phase 6; D1/D2/D3/D5/D6 with
  Phases 1–2.)
- **CLAUDE.md.** Confirmed the four listed sections current; added the
  `define_condition_change_action()` entry, the block-scoped note on the
  `phenotype_meta.condition_change_action` row, and the D8 paragraph.
- **Benchmark.** New `dev/benchmarks/benchmark_phenotype_scale.R`: PLAN /
  RESOLVE / COMMIT timed separately across four block shapes
  (`independent`, `correlated`, `conditional`, `named_effect`) and two passes,
  so the observation-pattern query reads apart from planning and from the
  batched write. Seeded, per CLAUDE.md. A `--check` mode prints seeded values
  as the guard that no future optimization changes RNG semantics. **Result: no
  optimization needed.** Per-individual cost *falls* with population size in
  every shape — roughly 3–4× from 1,000 to 16,000 individuals (independent:
  0.47 → 0.095 ms/ind). `commit` takes under 3× the time for 16× the rows, so
  the single transaction is batching; `plan` dominates the total everywhere;
  and the call-2-minus-call-1 gap in `resolve` — the observation-pattern
  query — stays in the tens of milliseconds rather than growing per entity.
- `R/schema.R` column descriptions: unchanged since Phase 1,
  `test-schema-print.R` green.
- `NEWS.md` under **0.71.0**; `DESCRIPTION` version.

*(Original scope, for the record:)*
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
- Make `condition_change_action` changeable on a defined block.
- Benchmark under `dev/benchmarks/` for large populations, optimizing the
  observation-pattern query and batched writes **without changing RNG
  semantics**.
- Decide whether "same seed reproduces within the current code" means
  bit-identical or identical within tolerance.

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
