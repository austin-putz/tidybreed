# Genome Effects v4.9 — term/member storage, executable contrasts, origin resolution

**Status:** design proposal, nothing implemented. **This file is canonical.**
**Revision:** v4.9 (2026-09-10). v4.2 answered the v4.1 pass of
`plans/update_genome_effects_v4_codex_review.md`; v4.3 renamed the result table;
v4.4 restored its component dimension and deferred the rest to
`plans/consolidate_genetic_values.md`; v4.5 closed the spec gaps found in a readiness
audit; v4.6 resolves imprinting; v4.7 fixes the three contracts that resolution
exposed and one pre-existing `CHECK` that never fired; v4.8 closed two
Phase-B blockers; **v4.9 specifies the two things every prior revision assumed —
the writer's input format and the evaluation strategy — and renames
`effect_set_name` to `effect_owner`.**

**Readiness: implementable.** No open question blocks any phase. The two remaining
are ergonomic and revisable after the fact.
★ **Phase A is complete** and required **no change to this schema** — results, hand
computations and four findings carried into B–D are in
`plans/update_genome_effects_phase_A.md`.
**Supersedes:** `update_genome_effects.md` (v1/v2, 16 tables) and
`update_genome_effects_v3_parsimony.md` (v3, 2 tables), retained as history.

---

## Revision history

**v1 → v2.** Generalized locus count but hard-coded per-locus state space; 16 tables.
**v2 → v3.** Indicator completeness removed the contrast registry and surface
families; 2 tables. **v3 → v4.** Storage-without-evaluation, hemizygous centering bug,
scalar `line_name`, no writer. **v4 → v4.1.** Origin precedence competed across
unrelated terms (a real bug); multi-locus fallback undefined; replacement had no
ownership boundary; dominance failed at evaluation. **v4.1 → v4.2** — four
correctness contracts, two of which were bugs I introduced in v4.1:

1. **Origin-row count is not a valid specificity order.** It ties a parent-qualified
   scope with its own generic fallback — which silently *disables the reciprocal
   feature v4.1 claimed to implement* — and it mis-ranks multi-row additive scopes
   whose AND/OR contract was never stated. **Replaced with predicate containment.**
2. **Dosage alone is not a complete genotype state.** At a variable-copy locus
   `dosage_value = 0` conflates "no copy", "one allele-0 copy", and "two allele-0
   copies". v4.1 directed users to indicators precisely at those loci, so **the
   escape hatch was broken where it was most needed.** One new column, no new table.
3. **`add_tbv()` cannot infer a breeding value from an arbitrary "additive" term.**
   Under the recommended functional `(a, d)` input the stored coefficient is `a`,
   while the average effect is `α = a + d(q − p)`. v4.1 said this about output labels
   and then failed to apply it to `add_tbv()`.
4. **The reserved and general writers shared `'default'`**, so rerunning
   `define_additive_effects()` in replace mode would delete a user's custom terms —
   the exact failure the effect-owner boundary was introduced to prevent.

**v4.2 → v4.3.** Maintainer correction, and it removes work rather than adding it.
A TBV is additive **by definition** — it does not become ambiguous when dominance
exists, so v4.2's Q1 ("what does `ind_tbv` mean once general effects exist?") was a
non-question. The genuinely new quantity is the **total genotypic value**, which needs
its own name, table, and function. v4.2's `ind_genetic_value` with a `component_name`
dimension is replaced by **`ind_tgv`**, structurally identical to `ind_tbv`. The
decomposition is deferred as a child table, per this plan's own decision rule.

**v4.3 → v4.4.** v4.3 dropped the component dimension from `ind_tgv` on the grounds
that `G − A` recovers the non-additive part. That was wrong for a **simulator**: the
evaluator already separates terms by family in order to sum them, so components are
free at the moment of computation and unreconstructable afterward. The dimension is
restored — and only the dimension. Adding `component_name` later would change the
unique key from `(id_ind, trait_name)` to `(id_ind, trait_name, component_name)`,
reshaping every row and reader, which is exactly what the decision rule forbids.
Everything else about consolidating genetic values — replacing `ind_tbv`, the
component vocabulary, surface classification, composite phenotypes, `ind_true_index` —
is rows, values, and migration, not shape, and moves to
**`plans/consolidate_genetic_values.md`**.

**v4.4 → v4.5.** Readiness audit. No design change; five spec gaps closed that would
each have stalled implementation: acceptance gates 1–25 existed only as a reference to
an overwritten revision and are now enumerated; the views Phases B and E depend on were
never defined; `ind_tgv`'s creation site was unstated; `remove_rows()` /
`archive_replicate()` interaction was unstated; and `add_phenotype()`'s status was
assumed rather than declared. One design hole found and **escalated, not closed** —
imprinting (Q1), which dropped out of the plan between v2 and v4 and which the origin
predicate now duplicates.

**v4.5 → v4.6.** Imprinting resolved (v4.5's Q1): `trait_meta.expressed_parent` is
**deleted** and imprinting becomes an origin predicate. Accepting that surfaced a
pre-existing bug — v4.5's containment lattice used `ANY` as a line value and its own
worked case `(ANY, 1)` was **unrepresentable in its own DDL**, because the `CHECK`
allowed only `'exact'` and `'unknown'` and `'any'` sat in Explicitly-out-of-scope.
`ANY` was reachable only as "zero origin rows", which cannot carry a `parent_origin`.
`line_match_type = 'any'` is now in scope, which both fixes that inconsistency and is
what makes "any line, from the sire" — imprinting — storable at all.

**v4.6 → v4.7.** The imprinting wrapper composes `line_name` with `parent_origin`
instead of stamping `'any'` unconditionally (it was erasing the line dimension);
`scale_to_target` becomes origin-aware (a parent-qualified term has one eligible copy,
not two, so imprinted models were landing at half the requested variance); multi-trait
`parent_origin` is specified as per-trait with mixed origins rejected under `G`; the
`'any'`-requires-a-parent rule moves from R into the `CHECK`; the no-eligible-copy
contract is stated; and the non-indicator `CHECK` branch gains
`center_value IS NOT NULL` — `BETWEEN` alone evaluates to UNKNOWN on NULL, which SQL
accepts, so the promised "every row-local invariant is a declared SQL constraint" was
false for the one column it mattered most on.

**v4.7 → v4.8** (2026-09-08). Final pre-implementation pass; no design change, two
blockers and three gaps. **(a)** The `locus_id` FK the plan has carried since v4.1 is
**not creatable**: `genome_meta` declares no `PRIMARY KEY` or `UNIQUE`
(`define_genome.R:272-275`), and DuckDB refuses an FK to an unconstrained column.
Creation order was necessary and not sufficient. **(b)** Adding the effect tables to
`GENOME_TABLES` while `open_pop()` still creates `genome_effects` makes
`define_genome()`'s preflight (`define_genome.R:146-153`) fire on **every fresh
population** — the two edits are one commit, not two. Plus: the `replace_scope`
behavior when only `parent_origin` changes, a `restore_pop()` contract for pre-change
database files, and view registration.

**v4.8 → v4.9** (2026-09-10). No new table, no new column, no change to the
resolution semantics. Four gaps, all of them in the layers *around* the schema:

1. **`terms` was never specified.** Every revision since v4.0 has taken a `terms`
   argument and described what the writer does *to* it, without once saying what a
   user types. That is the whole learnability surface of the feature. Specified in
   §Writer API with three worked examples.
2. **No evaluation strategy, and no performance gate among 50.** `add_tbv()` today is
   one set-based SQL statement per trait over `ind_haplotype` (`add_tbv.R:206-227`).
   The plan described per-individual tuple enumeration and never said how the two
   relate. Specified in §Evaluation strategy: **containment resolves once into a small
   variant map; evaluation stays set-based SQL**, and tuple enumeration is confined to
   scoped high-order terms. Gate 51.
3. **Phase D built `add_tbv()` twice.** It rewired `add_tbv()` to the reserved owner,
   and `plans/consolidate_genetic_values.md` then deletes that function. Phase D now
   builds **one** evaluator; `add_tbv()` becomes a thin filtered call into it, so
   nothing is thrown away when consolidation removes the wrapper. `ind_tbv` and
   `add_tbv()` survive this plan unchanged in meaning — non-additive genetic values in
   `ind_tbv` are the next plan's work, not this one's.
4. **`effect_set_name` → `effect_owner`.** The column is an ownership boundary for
   replacement; it is *in* the family signature, so owners always sum and can never be
   selected between. "Set" promised selectability the column does not have and the
   DDL comment had to disclaim. Renamed while pre-1.0 makes it free.

Plus: `component_name` values become explicit about order (`'order1_additive'`,
`'order1_dominance'`, `'order1_other'`, `'interaction'`) so they cannot be misread as
V_A / V_D / V_I; `family_key` is exposed on the term view so a user can *see* which
terms compete and which sum; and §Future limitations records six deferrals with the
shape of their eventual fix.

Size unchanged: **3 effect-definition tables + 1 result table.**

---

## The decision rule (unchanged)

> **Pre-build a dimension only when its absence would force existing rows to change
> shape. Everything else is a later child table.**

---

## Audit (verified against `4e9a72b`; ★ = added after v4.0)

1. **Multiallelic is blocked by the pipeline, not the datatype.** `add_dosage()`
   computes `CAST(SUM(h.allele) AS UTINYINT)` (`R/add_dosage.R:126-130`) — a count of
   allele 1, not a multiallelic genotype.
2. **No FK is possible in the current creation order.** `open_pop.R:286` creates
   `genome_effects`; `genome_meta` arrives at `define_genome.R:272`. All effect-table
   creation moves into `define_genome()`, beside `ind_haplotype`, `genome_map`, and
   `chr_inheritance` (`define_genome.R:300-371`).
   ★ **Order is necessary but not sufficient.** `genome_meta` is created with **no
   `PRIMARY KEY` and no `UNIQUE`** (`define_genome.R:272-275`), and DuckDB refuses
   `FOREIGN KEY … REFERENCES genome_meta(locus_id)` outright:
   *"Failed to create foreign key: there is no primary key or unique constraint for
   referenced table"* (verified, 1.5.5). `CLAUDE.md:174` documents `locus_id` as a
   primary key; the DDL never declared one. See §Schema.
3. ★ **DuckDB 1.5.5 enforces every proposed constraint — verified empirically.**
   Multi-branch `CHECK`, closed-set `CHECK`, range `CHECK`, single-column FK, and
   **composite `FOREIGN KEY (id, slot)`** all reject bad rows. Deletes do **not**
   cascade, so replacement must delete children before parents.
4. **`ind_tbv` is `UNIQUE (id_ind, trait_name)`** (`define_trait.R:222-228`).
5. ★ **`target_add_mean` is currently inert** — documented as "TBV centering mean for
   the base population" (`define_trait.R:20-23`), it appears nowhere outside
   `define_trait.R` and the two registries. `add_tbv()` never reads it.
6. **Per-copy additive semantics stay load-bearing.** `add_tbv.R:206-227` matches
   `e.line_name = h.line_origin` with a per-locus `NOT EXISTS` fallback. The
   `COALESCE(base_allele_freq, 0)` is a permissive bug.
7. ★ **Table budget:** 3 definition tables + 1 result table = **4 new tables, 24
   registry entries** (`sql_utils.R:85,135,161,252`; `schema.R:671`, `:100-117`).
8. **Blast radius:** 12 R files, 8 test files incl. `helper-parity.R`, 9 man pages,
   2 vignettes. ★ Plus one benchmark script under `dev/benchmarks/` for gate 51 and
   three views registered in the three schema lists.

---

## Indicator completeness — precise scope

Any function over a finite product of **modelled per-locus local states** is exactly a
sum of indicator terms. ★ A local state is `(realized_copy_count, dosage)`, not dosage
alone — see §Contrasts.

It does **not** cover phased haplotype combinations, multiallelic identity,
line-origin composition, sex, or environment. Those are separate dimensions, which is
why the origin table exists rather than being folded into indicators. ★ Line origin is
in scope here; the rest are recorded, with the shape of their eventual fix, in
§Future limitations.

---

## Schema

Created in `define_genome()`. Every row-local invariant is a declared SQL constraint
(audit item 3); cross-row rules are R-enforced.

★ **Prerequisite — `genome_meta` gains a primary key.** The `locus_id` FK below
cannot be declared without it (audit item 2):

```sql
-- define_genome.R:272-275, was: "locus_id INTEGER, locus_name VARCHAR, ..."
CREATE TABLE genome_meta (
  locus_id INTEGER PRIMARY KEY, locus_name VARCHAR, chr INTEGER,
  chr_name VARCHAR, pos_bp BIGINT
)
```

This makes `CLAUDE.md:174` true rather than aspirational, and is safe here where it
was not on `ind_haplotype`: `genome_meta` is `n_loci` rows written once, not
`n_ind × n_loci`, so the ART index the maintainers measured at ~53% of bulk-insert cost
(`define_genome.R:292`) is not a factor. Every existing writer reaches `genome_meta`
through `ALTER TABLE ADD COLUMN` + `UPDATE` — `define_chip.R:100`,
`founder_haplotype_helpers.R:219`, `mutate_table()` — and never rewrites the table, so
the constraint survives; verified that `ALTER`/`UPDATE` on a PK'd table works in 1.5.5.
`locus_id` is `seq_len(n_loci)`, so uniqueness is already an invariant, not a new one.

```sql
CREATE TABLE genome_effects (
  id_genome_effect INTEGER PRIMARY KEY,        -- next_int_id()
  trait_name       VARCHAR NOT NULL,           -- R-enforced FK to trait_meta
  effect_owner     VARCHAR NOT NULL,           -- who owns these rows for replacement; NOT a selectable scenario
  effect_name      VARCHAR,                    -- optional per-term label; no math meaning
  genome_value     DOUBLE  NOT NULL
);

CREATE TABLE genome_effect_members (
  id_genome_effect INTEGER  NOT NULL,
  member_slot      INTEGER  NOT NULL,          -- canonical: ascending locus_id
  locus_id         INTEGER  NOT NULL,
  contrast_name    VARCHAR  NOT NULL,
  copy_count_value UTINYINT,                   -- ★ realized copy count for indicator states
  dosage_value     UTINYINT,
  center_value     DOUBLE,
  PRIMARY KEY (id_genome_effect, member_slot),
  CHECK (contrast_name IN ('additive','dominance','indicator')),
  CHECK ( (contrast_name =  'indicator'
             AND copy_count_value IS NOT NULL
             AND dosage_value     IS NOT NULL
             AND dosage_value <= copy_count_value
             AND center_value     IS NULL)
       OR (contrast_name <> 'indicator'
             AND copy_count_value IS NULL
             AND dosage_value     IS NULL
             AND center_value     IS NOT NULL          -- ★ BETWEEN alone is UNKNOWN on NULL
             AND center_value BETWEEN 0 AND 1) ),
  FOREIGN KEY (id_genome_effect) REFERENCES genome_effects(id_genome_effect),
  FOREIGN KEY (locus_id)         REFERENCES genome_meta(locus_id)     -- ★ was promised, not declared
);

CREATE TABLE genome_effect_member_origins (
  id_genome_effect INTEGER  NOT NULL,
  member_slot      INTEGER  NOT NULL,
  origin_slot      INTEGER  NOT NULL,          -- canonical: sorted origin tuple
  line_match_type  VARCHAR  NOT NULL,
  line_name        VARCHAR,
  parent_origin    UTINYINT,                   -- NULL = either parent
  copy_count       INTEGER  NOT NULL,
  PRIMARY KEY (id_genome_effect, member_slot, origin_slot),
  CHECK (line_match_type IN ('exact','unknown','any')),
  CHECK (copy_count > 0),
  CHECK (parent_origin IS NULL OR parent_origin IN (1,2)),
  CHECK ( (line_match_type = 'exact'   AND line_name IS NOT NULL)
       OR (line_match_type = 'unknown' AND line_name IS NULL)
       OR (line_match_type = 'any'     AND line_name IS NULL
                                       AND parent_origin IS NOT NULL) ),  -- ★
  FOREIGN KEY (id_genome_effect, member_slot)
    REFERENCES genome_effect_members(id_genome_effect, member_slot)
);
```

`center_value` (renamed from `base_allele_freq` in v4.1) is a **per-copy centering
constant**: the base allele frequency `p` under Cockerham coding, `0.5` under
functional coding. Naming a stored `0.5` a "base allele frequency" would be false.

★ **`ind_tgv.id_ind` gets no SQL FK**, matching `ind_tbv`'s existing
precedent and keeping `remove_rows()` / `archive_replicate()` unblocked by referential
constraints. R-enforced, like `trait_name`.

### ★ Table lifecycle

| Table | Created in | Alongside |
|---|---|---|
| `genome_effects`, `genome_effect_members`, `genome_effect_member_origins` | `define_genome()` | `ind_haplotype`, `genome_map`, `chr_inheritance` (`define_genome.R:300-371`) — this is what makes the `locus_id` FK possible (audit item 2) |
| `ind_tgv` | `define_trait()`'s lazy DDL block | `ind_tbv` (`define_trait.R:221`), `ind_phenotype`, `ind_ebv` — both result tables then appear at the same lifecycle point |

★ **The `open_pop()` DDL and `GENOME_TABLES` move in one commit.**
`define_genome()`'s preflight refuses to run if **any** member of `GENOME_TABLES`
already exists, on the stated premise that "open_pop() creates none of these seven"
(`define_genome.R:141-153`). Adding the three effect tables to `GENOME_TABLES` while
`open_pop.R:286` still creates `genome_effects` makes that premise false and
`define_genome()` errors on **every fresh population**. So, atomically: delete the
`genome_effects` DDL from `.create_core_tables()`, add all three tables to
`GENOME_TABLES`, and update the preflight comment (seven → ten).

★ **`restore_pop()` and pre-change database files.** `genome_effects` keeps its name
with an incompatible schema, and `restore_pop()` does no structural validation
(`restore_pop.R:77-140`) — an old `.duckdb` would open cleanly and then fail deep inside
`add_tbv()` with an SQL column error. Pre-1.0 owes no migration, but it does owe a
legible failure: follow the existing precedent at `restore_pop.R:134-139`, which already
detects and reports a pre-v0.36.0 database. Detect `genome_effects` **without**
`genome_effect_members` and stop with "this population predates the term/member effect
schema; re-create it with `define_genome()`."

★ **`remove_rows()` and `archive_replicate()`.**

- `TABLE_ROW_KEYS` (`sql_utils.R:161`) gains: `genome_effects = "id_genome_effect"`;
  `genome_effect_members = c("id_genome_effect", "member_slot")`;
  `genome_effect_member_origins = c("id_genome_effect", "member_slot", "origin_slot")`;
  `ind_tgv = c("id_ind", "trait_name", "component_name")`.
- **`remove_rows()` on `genome_effects` must delete children first.** DuckDB refuses a
  parent delete with live children and does not cascade (audit item 3), so a naive
  single-table delete errors. Either teach `remove_rows()` the child order for these
  three tables or reject the call and point at `define_genome_effects(mode = ...)`.
  **Rejecting is preferred** — effect definitions are configuration and should be
  replaced through the writer, not row-deleted.
- `archive_replicate()` stamps `replicate` on **output** tables only. `ind_tgv` has the
  column and is stamped; the three effect-definition tables are configuration and are
  **not** stamped, exactly as `genome_effects` is treated today.

★ **`add_phenotype()` is unchanged by this plan.** It continues to call `add_tbv()`
internally and read `ind_tbv` for simple, composite, and SGE assembly
(`add_phenotype.R:754-763, 979-983, 1025-1029`). Switching contributor lookups to the
`ind_tgv` total is Task 3 of `plans/consolidate_genetic_values.md`, where it is a
no-op today and a correctness fix once dominance exists.

---

## Contrast definitions — executable specifications

`c` = the member's `center_value`; where `c = p`, `q = 1 − p`.

### `additive` — ploidy-aware

```
x_A(i, j) = Σ over eligible copies of individual i at locus j  of  (allele − c)
```

`dosage − 2c` at a diploid locus, `dosage − c` hemizygous, `0` when absent. Not
`dosage − 2p`.

### `dominance` — diploid loci only

HWE-orthogonal Cockerham coding, `g` = dosage of allele 1, `c = p`:
`x_D(0) = −2p²`, `x_D(1) = 2pq`, `x_D(2) = −2q²`.
Orthogonality (`E[x_D] = 0`, `E[x_A x_D] = 4p²q²(1 − p − q) = 0`) is a committed test.

★ **Rejected at write time** at any locus whose resolved `chr_inheritance` is not
`1,1` for every applicable offspring sex — **unless an exact origin multiset on that
same dominance member** has requested counts summing to 2. A diploid-proving scope on
a *different* member of the term proves nothing about this locus.

### ★ `indicator` — the local state is `(copy_count, dosage)`

```
x_I(i, j) = 1 if (realized copy count, dosage of allele 1) at locus j
              equals (copy_count_value, dosage_value), else 0
```

Dosage alone is not a state: at a variable-copy locus, dosage 0 is three different
biological situations. Both fields are **required** for indicators; the writer infers
`copy_count_value` for ordinary autosomal input, so most users never supply it.

★ **The evaluator must synthesize the zero-copy state.** An individual with no
`ind_haplotype` rows at a locus (Y in a female) has no row to join, so an inner join
would silently omit the `copy_count_value = 0` indicator. Use a left join against the
locus list and materialize the absent state.

### Aggregation for multi-locus terms

Each member reduces to **one value per individual** first; members are then
multiplied. Multiplying raw haplotype rows across members produces a Cartesian pairing
and is never correct.

★ Under an origin-scoped family the reduction grain is one value per
**(individual, label)** rather than per individual, because a different variant may be
selected per label — see §Evaluation strategy. The rule that raw haplotype rows are
never paired across members is unchanged, and the unscoped case is the one-label
special case of it.

### Functional coding needs no new contrast

Functional `x_A(g) = g − 1` is `additive` with `center_value = 0.5`; functional
`x_D(g) = 1[g = 1]` is `indicator` with `(copy_count_value, dosage_value) = (2, 1)`.

---

## Fallback families — precedence operates only within a family

Terms from different families always **sum**; only scope variants of the same
mathematical term compete. Family signature:

```
(trait_name, effect_owner,
 ordered [(locus_id, contrast_name, copy_count_value, dosage_value)])
```

- ★ `copy_count_value` is **included** — two indicator terms at one locus with the
  same dosage but different copy counts are different basis functions and must sum.
- `center_value` is **excluded** — a line-A variant legitimately carries a different
  frequency from its common fallback.
- Origin rows and `genome_value` are **excluded** — they distinguish variants *within*
  a family.
- `effect_owner` is **included**: because all owners sum, two owners defining the same
  common term must both fire rather than collide as a tie.

★ This signature is exposed to users as `family_key` on the `genome_effect_terms` view
(§Views). It is the only way to see, without reading this section, whether two
definitions compete or sum.

**Resolver contract:** partition into families → resolve specificity within each →
evaluate the selected variant(s) → **sum across families**. An additive main effect, a
dominance main effect, and an interaction sharing a locus are three families and never
suppress one another.

★ **Same-family, same-scope, different `center_value`.** Because centering is outside
the signature, a functional `additive`@0.5 and a Cockerham `additive`@`p` at one locus
under one owner collide as a duplicate logical term. They are genuinely combinable —
`a₁(g − c₁) + a₂(g − c₂) = (a₁ + a₂)(g − c′)` with
`c′ = (a₁c₁ + a₂c₂)/(a₁ + a₂)` — so the writer **rejects with the combined equivalent
in the message** rather than a bare duplicate error. See Q2.

---

## ★ Origin resolution — predicate containment, not row count

**This replaces v4.1's row-count rule, which was wrong.** Row count ties a
parent-qualified scope with its own generic fallback, so v4.1 would have *rejected*
the reciprocal models it claimed to support.

**Evaluation units.** `additive` → one unit per eligible allele copy.
`dominance` / `indicator` → one unit, the locus local state.
**Evaluation tuple** = one unit from each member (Cartesian product).

**Selection.** Among the variants of a family whose predicates match a tuple:

1. the common scope (zero origin rows) matches every tuple and is least specific;
2. if two matching predicates overlap and one is a **strict subset**, the subset wins;
3. if they overlap and **neither contains the other**, reject at write time;
4. equal predicates are a duplicate family/scope identity.

No match → the tuple contributes 0.

★ **…but a trait that matches *nothing* still errors.** `add_tbv()` stops today when an
individual has no haplotype row matching any additive effect (`add_tbv.R:230-247`), and
the message names both causes: an uninherited chromosome, and an imprinted trait
restricted to a parent the individual has no copies of. Per-tuple zero and that
diagnostic are not in conflict — they act at different granularities, and both are kept:

- a tuple whose predicate is unmatched contributes **0**; this is the definition of the
  predicate system and is what a *partial* miss does;
- an individual for whom **every** term of a trait matches no tuple is an **error**,
  preserving today's behavior message and all.

The guard is what stops a paternally-qualified X-linked trait from silently scoring every
male at 0 — the exact failure the current message was written to explain. It is
structural, not numerical: a matched tuple that evaluates to 0 is a value, not a miss.

### Containment is decidable over two small finite lattices

★ **At most one origin row per additive member per variant.** An A-row plus a B-row on
one additive member is incoherent: under OR it matches *more* copies than A-only while
row count calls it narrower; under AND no single copy can be both. "A or B" is
**expanded by the writer into separate variants**. `copy_count = 1` is also required
there, since additive matching is per copy.

**Additive predicate** = `(line, parent)` over one copy:

| Dimension | Values | Order |
|---|---|---|
| line | `ANY` (`line_match_type = 'any'`, or zero origin rows) · `exact(L)` · `unknown` | `ANY ⊐ exact(L)`, `ANY ⊐ unknown`; `exact(L)`, `exact(L')`, `unknown` pairwise **disjoint** |
| parent | `ANY` (`parent_origin IS NULL`) · `1` · `2` | `ANY ⊐ 1`, `ANY ⊐ 2`; `1`, `2` disjoint |

★ **`'any'` exists so that `(ANY, parent)` is storable.** `ANY` on the line dimension
was otherwise reachable only as *zero origin rows*, and a row is where `parent_origin`
lives — so "any line, from the sire" had no encoding. Two rules keep it from becoming
a second way to say something that already has one:

- an `'any'` row **must** carry a non-NULL `parent_origin`; `('any', NULL)` is the
  common scope written the long way and is **rejected** — ★ by the `CHECK`, not by the
  validator, since the rule is row-local;
- `'any'` is permitted **only on `additive` members**, where matching is per copy. On a
  genotype member the multiset must name lines and sum to the realized copy count, so
  an `'any'` entry constrains nothing and is **rejected**.

`P ⊆ Q` iff componentwise. Worked cases:

| P vs Q | Result |
|---|---|
| `(ANY,ANY)` vs `(exact A, ANY)` | strict subset → A wins on A copies (**today's behavior**) |
| `(exact A, ANY)` vs `(exact A, 1)` | strict subset → paternal-A wins (**the tie v4.1 got wrong**) |
| `(exact A, 1)` vs `(exact B, 2)` | disjoint → both apply, different copies, no conflict |
| `(exact A, ANY)` vs `(ANY, 1)` | overlap, neither contains → **rejected at write time** |

**Genotype predicate** = an exact multiset of `(line, parent, count)` demands whose
counts sum to the realized copy count. Different line multisets are **disjoint**;
within one line multiset, `P ⊆ Q` iff `P`'s parent assignments refine `Q`'s.

| P vs Q | Result |
|---|---|
| common vs `{A:1, B:1}` | strict subset |
| `{A:1, B:1}` parents `ANY` vs `{A:1@p1, B:1@p2}` | strict subset → **reciprocal override works** |
| `{A:1@p1, B:1@p2}` vs `{A:1@p2, B:1@p1}` | disjoint → both valid |
| `{A:2}` vs `{A:1, B:1}` | disjoint |

★ `origin_slot` is canonicalized by sorting
`(line_match_type, line_name, parent_origin, copy_count)`, NULLs last.

### Contribution and resource guard

```
family contribution = Σ over tuples of [ selected variant's genome_value
                                         × ∏ over members of (contrast value at that unit) ]
```

★ **Fast path — for a *common/unscoped* family only.** When no variant in a family
carries any origin row, the sum factors exactly into `∏ over members of (Σ over
units)` — plain per-member reduction. v4.1 wrongly called this the "additive-only"
path: **line-specific order-one additive effects do carry origin rows** and do not use
it. Their tuple cost is trivial (one tuple per copy), but the distinction matters for
correctness of the description.

★ **This is the semantic definition, not the execution plan.** Written per tuple it
looks like a per-individual loop; it is not one. Because a predicate reads only a
copy's `(line, parent)` **label**, tuples group by label-vector and the inner sum
factors inside each group — so the fast path above is the one-group special case of a
single set-based formulation. See **§Evaluation strategy**, which is binding on the
implementation.

★ **Tuple preflight.** Enumeration is over **label-vectors**, `∏_m |labels_m|` per
family, where a member's label alphabet is the small distinct set defined in
§Evaluation strategy — not `O(ploidy^order)` per individual. The evaluator estimates that count before running,
warns at a configurable threshold (default 10⁴ per family) and errors above a hard cap
(default 10⁶), so an accidental high-order scoped term fails loudly instead of
appearing to hang.

**`{A, NULL}` worked case.** One `line_origin = 'A'` copy and one `NULL` copy, with an
exact-A variant and a common variant: the A copy selects A, the NULL copy selects
common. No double count, no within-A masquerade.

★ **Imprinting worked case.** A paternally expressed additive effect is one `additive`
term carrying a single origin row `(line_match_type = 'any', parent_origin = 1,
copy_count = 1)`. Per-copy matching makes only paternal copies eligible, so
`x_A = Σ over paternal copies (allele − c)` — **numerically identical to today's
`expressed_parent = "parent_1"` filter** (`add_tbv.R:191-194`), but per locus and per
effect owner rather than per trait.

The incoherence v2 worried about resolves **structurally rather than by a validator
rule**: `'any'` is rejected on genotype members, so a `dominance` or `indicator`
contrast simply cannot be restricted to one parent's copy — which is right, since a
genotype needs both. Parent-*of-origin*-specific dominance remains expressible the
correct way, as an exact multiset with parent assignments (`{A:1@p1, B:1@p2}`) — the
reciprocal case already in the lattice above.

---

## ★ Evaluation strategy — resolve once, evaluate set-based

**This section did not exist before v4.9, and its absence was the largest
implementation risk in the plan.** `add_tbv()` today is *one* set-based SQL statement
per trait over `ind_haplotype` (`add_tbv.R:206-227`). Nothing above said how the tuple
language relates to that, and a literal reading — enumerate tuples per individual in R —
would replace a single join over millions of rows with a per-individual loop. The
package keeps three optimization plans and a `dev/benchmarks/` directory precisely
because this table is the hot one.

It does not have to be a loop. Two observations collapse the whole thing back into SQL.

**1. A predicate reads a *label*, never an individual.** Origin predicates are defined
over `(line_name, parent_origin)` and nothing else. So the winning variant is a function
of the label, not of who carries it. The label's shape follows the member's evaluation
unit, exactly as §Origin resolution defines it:

| Member contrast | Unit | Label | Alphabet |
|---|---|---|---|
| `additive` | one allele copy | `(line_origin, parent_origin)` | `SELECT DISTINCT line_origin, parent_origin FROM ind_haplotype` — bounded by `#lines × 2`, typically 2–8 rows |
| `dominance`, `indicator` | the locus local state | the **sorted multiset** of the copies' `(line_origin, parent_origin)` at that locus | the distinct multisets that actually occur — a `GROUP BY (id_ind, locus_id)` with a canonical string aggregate, then `DISTINCT` |

Both alphabets are computed once, in SQL, and both are tiny: a two-line cross has a
handful of genotype label multisets, not one per individual. `NULL` `line_origin` is a
label like any other and is what `line_match_type = 'unknown'` matches.

**2. Tuples group by label-vector, and the inner sum factors inside each group.**
Writing `L` for a vector of one label per member:

```text
family contribution = Σ over label-vectors L
                        [ variant(L).genome_value
                          × ∏ over members m ( Σ over units of m carrying label L_m
                                                 of the contrast value ) ]
```

`variant(L)` is constant within a group by construction, which is what makes the
factorization legal. The enumeration is therefore over **label-vectors, not over
individual allele copies**: `∏_m |labels_m|`, not `O(ploidy^order) per individual`. In
a purebred unscoped population every copy shares one label and there is exactly one
group — the fast path falls out as a special case rather than being a separate code
path.

### The three artifacts

| Artifact | Grain | Built |
|---|---|---|
| **Label alphabet** | one row per distinct label, per member kind (copy label or locus multiset) | once per evaluation, one `DISTINCT` each |
| **Resolved variant map** | one row per `(family_key, label-vector) → id_genome_effect` | once per evaluation, in R, from stored effect rows only — never per individual |
| **Member reduction** | one row per `(id_ind, id_genome_effect, member_slot, label)` | SQL, one statement per contrast type |

Containment search runs **only** while building the map, over a table whose size is
`#families × ∏_m |labels_m|`. It never runs during evaluation. The map is total and
unambiguous by construction: overlapping-but-incomparable predicates are rejected at
**write** time, so no tie can reach the evaluator.

### Order-one terms keep today's shape

An order-one family is one member, so a label-vector is one label, and the whole
evaluation is:

```text
ind_haplotype ⋈ resolved_variant_map (on line_origin, parent_origin)
              ⋈ genome_effect_members (on locus_id)
   → SUM(genome_value × contrast_value) GROUP BY id_ind
```

That is structurally the statement `add_tbv()` already runs, with the correlated
`NOT EXISTS` fallback replaced by a join against a precomputed map — strictly less work
than today, not more. **The line-fallback semantics move from a subquery into data.**

### Higher-order terms stay set-based

Each member reduces to one value per `(id_ind, label)` in SQL; members are then joined
on `id_ind` and multiplied. DuckDB's `product()` aggregate does the multiplication
directly, so an order-`k` term is `k` reductions plus `k − 1` joins — no row ever
carries a Cartesian pairing of raw haplotype rows, which §Aggregation already forbids.

★ **The zero-copy state is materialized, not joined away.** An individual with no
`ind_haplotype` row at a locus (Y in a female) has no row to reduce, so member
reduction left-joins against the `(individual × locus)` list and synthesizes
`copy_count = 0`.

★ **Open for Phase D — what does that left join join against?** Phase A showed the
zero-copy state is sharper than the Y-in-a-female framing suggests: a
`(copy_count 0, dosage 0)` indicator matches **every** individual with no row at that
locus, so it is a "carries no copy here" effect. The candidates are every locus in
`genome_meta`, or the loci the individual is expected to carry under
`chr_inheritance`. Recommend the former — it is what the Phase A fixtures do — with a
row that `chr_inheritance` says should exist but does not treated as a data error
surfaced elsewhere, not silently scored as absent.

### Tuple preflight, restated

The guard from §Contribution now has a computable subject: it estimates
`∏_m |labels_m|` **per family**, warns at 10⁴ and errors at 10⁶ (both configurable).
The earlier per-individual framing made the cap either meaningless or unreachable; this
one bounds the resolved map, which is the object that actually grows.

### The contract, stated as a gate

Evaluation of an order-one additive model **must not** degrade materially against the
current engine. Gate 51 measures it (`dev/benchmarks/`, following the precedent in
`benchmark_haplotype_scale.R`). If an implementation cannot hold that line, the
resolution artifacts are wrong — not the target.

---

## ★ Effect owners, TBV contract, and replacement

★ **Why `effect_owner`, not `effect_set_name` (renamed in v4.9).** The column is in
the family signature, so terms under different owners **always sum** — a user can never
define "model A" and "model B" and evaluate one of them. Its only job is to say which
writer owns a group of rows for replacement purposes, which is what the `replace_*`
modes below key on. "Set" advertised a selectable scenario the column does not
implement, and the v4.8 DDL comment had to carry an explicit `NOT selection`
disclaimer to say so. Pre-1.0 makes the rename free; scenario comparison stays what it
already is in this package — separate populations or separate replicates.

**Two distinct defaults.** `define_additive_effects()` owns the reserved owner
`generated_additive_tbv`; `define_genome_effects()` defaults to `custom`. v4.1 had
both on `'default'`, so rerunning the generator in replace mode would have deleted a
user's custom terms. The general writer **refuses to mutate reserved owner names**
without an explicit advanced override; reserved names are a closed package-owned list.

**`add_tbv()` reads only order-one `additive` variants from `generated_additive_tbv`.**
Filtering on `contrast_name = 'additive'` alone is not sufficient: under functional
`(a, d)` input the stored coefficient is `a` while the breeding-value coefficient in a
diploid HWE base is `α = a + d(q − p)`, and under epistasis average effects also
depend on other loci and LD. So:

- arbitrary terms written through `define_genome_effects()` contribute to
  `ind_tgv` but **never silently redefine TBV**;
- `add_tbv()` ignores additive members that appear inside interactions;
- `ind_tbv` keeps its exact current meaning and its committed oracle;
- deriving average effects from a general non-additive model is a separate
  calculation, documented as such. See Q1.

★ **Multi-line replacement needs scope-level granularity.** Common, A, and B additive
variants must share an effect owner *and* family so they compete rather than sum — so a
whole-owner replace on a per-line `define_additive_effects()` call would delete the
others. Putting each line under its own owner is **not** a fix: `effect_owner` is in
the family signature, so separate owners sum and give no common fallback. Of the four
modes below, `define_additive_effects()` uses `replace_scope`:

- `append` — insert into `(trait_name, effect_owner)`; reject duplicate
  family + scope identities.
- **`replace_scope`** — delete only variants whose origin predicate is **equal** to the
  supplied scope, within `(trait_name, effect_owner)`, then insert. Predicate
  equality is decidable from the containment lattices above. This is what
  `define_additive_effects()` uses, so successive common / A / B calls each replace
  only their own variant.
- `replace_owner` — delete and replace the whole `(trait_name, effect_owner)`,
  children before parents (audit item 3).
- `replace_trait` — clears every owner's terms for the trait. Never a default.

★ **Changing `parent_origin` on a re-run adds a variant; it does not replace one.**
`replace_scope` keys on origin-predicate *equality*, and `parent_origin` is part of the
predicate — so `define_additive_effects(line_name = "A")` followed by
`define_additive_effects(line_name = "A", parent_origin = 1)` leaves **both**
`(exact A, ANY)` and `(exact A, 1)` standing. They are in a containment relation, so the
result is a legal fallback pair: paternal copies take the imprinted value, maternal
copies fall back to the generic one. That is correct by the rules and almost certainly
not what a user re-running the call intended, and it is a behavior change from
`expressed_parent`, which was one trait-wide flag that could only be overwritten.

The mechanism stays as-is — it is the same mechanism gate 35 relies on for the
common/A/B line fallback — but the writer **warns** on exactly the confusable case: a
write that creates a containment pair whose two members differ **only in the parent
dimension** (`parent ANY` vs `parent 1`/`2` under the same line predicate). That never
fires for the line-fallback pattern, whose members differ in *line*, and never for a
reciprocal pair `(exact A, 1)` / `(exact A, 2)`, whose parents are disjoint rather than
nested. The message names both variants and points at `replace_owner`.

---

## Output contract — `ind_tgv`

**TGV = true genetic value.** The `t`/`e` prefix already carries meaning in this
package (`tbv` = true, `ebv` = estimated), and `ind_tgv` vs `ind_ebv` reads correctly
where `ind_gv` vs `ind_ebv` would lose the ground-truth distinction.

```sql
CREATE TABLE ind_tgv (
  id_tgv         INTEGER PRIMARY KEY,
  id_ind         VARCHAR NOT NULL,
  trait_name     VARCHAR NOT NULL,
  component_name VARCHAR NOT NULL,
  tgv_value      DOUBLE  NOT NULL,
  replicate      INTEGER,
  UNIQUE (id_ind, trait_name, component_name)
);
```

`replicate` mirrors `ind_tbv`'s column (`sql_utils.R:97`) so `archive_replicate()`
works identically. Written by **`add_tgv()`**, beside `add_tbv()` and `add_ebv()`.
Writes are idempotent — re-evaluation replaces an individual's rows for a trait in one
transaction.

**Long, not wide.** The component dimension is open: epistasis subdivides into A×A,
A×D, D×D and higher. Long format splits it as **rows**; wide columns would need a
schema change. Same precedent as `genome_effect_type` and the package's schema bias.

★ **`component_name` vocabulary — deliberately minimal, and explicitly about
declared model structure:**

| Value | Written when |
|---|---|
| `'order1_additive'` | order-1 term whose member contrast is `additive` |
| `'order1_dominance'` | order-1 term whose member contrast is `dominance` |
| `'order1_other'` | order-1 `indicator` term — a hand-entered surface mixes additive and dominance by construction, and separating them needs an orthogonal projection against base frequencies |
| `'interaction'` | any term with ≥ 2 members |

★ **Why the names carry the order.** These are **model-structure components, not
variance components.** A functional A×A term contributes to A, D *and* I in the
statistical sense; the label records how the term was *declared*, not an orthogonal
decomposition, and nothing here may be read as `V_A` / `V_D` / `V_I`. Bare
`'additive'` / `'dominance'` / `'epistatic'` invited exactly that misreading — they are
the names of the variance components — and the disclaimer that had to follow them was
a paragraph no user reads before filtering a table. `'order1_additive'` is ugly enough
to prompt a documentation lookup, and `'interaction'` says "≥ 2 members in one term",
which is the true criterion, where `'epistatic'` claimed a genetic interpretation the
value does not have.

**Total is derived by a view, never stored.** A stored `'total'` row would make every
`SUM(tgv_value)` double-count.

★ **Known, time-boxed redundancy.** `ind_tbv` survives this plan unchanged, so the
additive value appears in two places. They are not always the same number:
`ind_tbv.tbv_value` is the breeding value from the reserved `generated_additive_tbv`
owner, while `ind_tgv`'s `'order1_additive'` row sums **every** additive-structured term
in the model, including custom ones. Identical in the common case, subtly different with
custom additive terms — a user trap, and the reason consolidation is scheduled as the
very next plan rather than left open-ended. See `plans/consolidate_genetic_values.md`.

### Mean contract

- `ind_tgv` stores the **raw sum of stored terms. No mean is added.**
- A pure Cockerham model yields centered deviations; a raw `indicator` surface yields
  absolute genotypic values with a non-zero mean by construction, and is not silently
  re-centered.
- `trait_meta.target_add_mean` keeps its documented meaning and is **not** an
  effect-model intercept. Audit item 5 shows it is currently inert.
- ★ The `(a, d)` helper **reports** the implied mean of the *genetic component*,
  `μ = a(p − q) + 2pq·d`, and **writes it nowhere**. v4.1 told the user to place it in
  `phenotype_meta.mean`; that is **removed** — it would double-count once non-additive
  genetic values are integrated, because the raw genetic values already have
  expectation `μ`. When phenotype integration lands, the intercept must satisfy
  `intercept = M − E[genetic value] − other expected components`; for a pure raw
  functional model with desired mean `μ` the additional intercept is **zero**.
- ★ The reported "running total for a set" is valid only for **appended single-locus
  main-effect terms**. The expectation of a general epistatic term depends on joint
  genotype frequencies and LD, so no running sum is reported for those.
- No zero-member intercept term. Line- or cross-specific intercepts are out of scope.

---

## Writer API

```r
define_genome_effects(
  pop, trait_name, terms,
  effect_owner  = "custom",
  mode             = c("append", "replace_scope", "replace_owner", "replace_trait"),
  origin           = NULL,
  require_complete = FALSE
)
```

Resolves `locus_name` → `locus_id`; canonicalizes members by ascending `locus_id` and
origin rows by the documented tuple; expands multi-line additive scopes into separate
variants; infers `copy_count_value` for autosomal indicator input; assigns IDs via
`next_int_id()`; runs the full validator; writes in **one transaction**.

★ `require_complete = TRUE` validates coverage over all reachable
`(copy_count, dosage)` states — not dosage values alone.
★ `effect_owner`, `effect_name`, and line names take the package's normal
identifier/string validation.
`effect_name` is a per-term label supplied inside `terms`; there is no function-level
`effect_name` argument.

### ★ The `terms` input format

Every revision before v4.9 described what the writer does *to* `terms` without saying
what a user types. This is the entire learnability surface of the feature, so it is
specified before anything is built.

**`terms` is a long data frame (or tibble): one row per (term × member).** No nesting,
no S3 constructor, no surrogate IDs, no `member_slot` arithmetic — the same shape as
every other table the user already reads with `get_table()`.

| Column | Required | Meaning |
|---|---|---|
| `term_id` | yes | Groups rows into one term. **User-facing only; never stored.** Any atomic type; the writer replaces it with `next_int_id()` |
| `genome_value` | yes | The term's coefficient. Constant within a `term_id` — a varying value is an error, not a recycle |
| `effect_name` | no | Per-term label. Constant within a `term_id` |
| `locus_name` | yes | Resolved to `locus_id`; each locus at most once per term |
| `contrast_name` | yes | `"additive"`, `"dominance"`, or `"indicator"` |
| `center_value` | for non-indicator | `p` for Cockerham, `0.5` for functional |
| `copy_count_value`, `dosage_value` | for indicator | The local state. `copy_count_value` is inferred for ordinary autosomal loci and may be omitted |

Single-term calls may omit `term_id`; a one-row `terms` needs only
`locus_name`, `contrast_name`, `genome_value`, and a center.

**`origin` carries scope**, kept out of `terms` so the common case stays flat:

- `NULL` (default) — the common scope: zero origin rows, matches every tuple.
- **a named scalar list** — one scope applied to every member of every term, e.g.
  `list(line_name = "Duroc", parent_origin = 1)`. This is what
  `define_additive_effects()` uses.
- **a data frame** — per-member scopes, for exact multisets: columns `term_id`,
  `locus_name`, `line_match_type`, `line_name`, `parent_origin`, `copy_count`. Keyed by
  `locus_name`, not `member_slot`, so the user never touches canonical slot order.

The scalar-list form produces exactly one origin row per member, which is all an
`additive` member may have and is generally *not* enough for a genotype member, whose
multiset must sum to the realized copy count — so a scoped `dominance` or `indicator`
term normally needs the data-frame form. The validator says which, naming the locus.

#### Example 1 — one dominance term (Cockerham, `p = 0.3`)

```r
pop <- pop |> define_genome_effects(
  trait_name = "ADG",
  terms = data.frame(
    locus_name    = "Locus_10",
    contrast_name = "dominance",
    center_value  = 0.3,
    genome_value  = 0.8
  )
)
```

The `(a, d)` helper writes the additive partner in the same call; this is the shape it
expands to.

#### Example 2 — a hand-entered 3×3 A×A surface

Nine cells, nine terms, two members each. A surface is **rows, not a second
representation** — this is indicator completeness made concrete.

```r
cells <- expand.grid(g1 = 0:2, g2 = 0:2)
cells$value <- c(0, 0, 0,  0, 1.4, 2.1,  0, 2.1, 3.6)   # hand-computed

surface <- rbind(
  data.frame(term_id = seq_len(nrow(cells)), locus_name = "Locus_10",
             contrast_name = "indicator", dosage_value = cells$g1,
             genome_value  = cells$value),
  data.frame(term_id = seq_len(nrow(cells)), locus_name = "Locus_44",
             contrast_name = "indicator", dosage_value = cells$g2,
             genome_value  = cells$value)
)

pop <- pop |> define_genome_effects("ADG", surface[surface$genome_value != 0, ],
                                    effect_owner = "epistasis_AxA")
```

A cell you drop is a term you did not write, contributing zero — the only consistent
reading of a component-valued table. `require_complete = TRUE` is how a user asks for
the opposite: it **rejects** the filtered call above and demands every reachable
`(copy_count, dosage)` state on every member, including `copy_count_value = 0` where a
chromosome can be absent. Sparse by default, complete on request.

#### Example 3 — reciprocal dominance, per-member exact multiset

The `{A:1@p1, B:1@p2}` case from the containment lattice: an F1's dominance value
depends on which parent contributed which line.

```r
terms <- data.frame(
  term_id = 1L, locus_name = "Locus_10", contrast_name = "dominance",
  center_value = 0.3, genome_value = 1.2
)
origin <- data.frame(
  term_id         = 1L,
  locus_name      = "Locus_10",
  line_match_type = "exact",
  line_name       = c("Duroc", "Landrace"),
  parent_origin   = c(1L, 2L),
  copy_count      = c(1L, 1L)
)

pop <- pop |> define_genome_effects("ADG", terms, origin = origin,
                                    effect_owner = "reciprocal")
```

Writing the mirror image (`Duroc@2`, `Landrace@1`) as a second call gives the disjoint
reciprocal variant. The two never compete: disjoint predicates both apply, to different
individuals.

**Validation the format itself buys.** Constant `genome_value` per `term_id`; each
locus at most once per term; `contrast_name` ↔ state/center agreement checked before
any SQL is issued; unknown columns rejected rather than silently dropped. Error
messages name the `term_id` the user typed, never an `id_genome_effect` they have never
seen.

Helpers convert a genotype table into indicator terms and an `(a, d)` pair into
`additive` + `dominance` (or functional `additive`@0.5 + `indicator`@(2,1)) — no
additional tables. `define_additive_effects()` is reimplemented on the general writer,
owns `generated_additive_tbv`, and uses `replace_scope`.

★ **`define_additive_effects()` gains `parent_origin = NULL`** — the ergonomic
replacement for the deleted trait-wide `expressed_parent` flag: one argument reproduces
whole-genome imprinting, while the general writer still expresses the per-locus case the
flag never could. It **composes** `line_name` with `parent_origin` into the single origin
row an additive member is allowed. It must **not** stamp `'any'` unconditionally: that
erases the line dimension, so two line-specific imprinted calls
(`test-add_tbv.R:235-270` writes Duroc at 1.0 and Landrace at 4.0) would land on one
identical scope and the second would replace the first.

| `line_name` | `parent_origin` | Stored additive scope |
|---|---|---|
| `NULL` | `NULL` | zero origin rows (common) |
| `"A"`  | `NULL` | `('exact', 'A', parent NULL, copy_count = 1)` |
| `NULL` | `1` / `2` | `('any', NULL, parent, copy_count = 1)` |
| `"A"`  | `1` / `2` | `('exact', 'A', parent, copy_count = 1)` |

One row in every scoped case, consistent with the additive-member invariant. Only the
third combination needs `'any'`, and is the whole reason `'any'` exists.

★ **`parent_origin` is per trait.** `trait_name` accepts a vector and the deleted column
was per trait, so a scalar-only argument would lose exactly the expressiveness the
removal was meant to preserve. It takes a scalar (recycled) or a vector matching
`trait_name`, optionally named by trait. **A mixed-origin call is rejected when `G` is
supplied**: under random mating the paternal and maternal copies at a locus are
independent, so a paternal-only trait and a maternal-only trait have **zero** genetic
covariance however strongly their sampled coefficients are correlated — the requested
off-diagonal is unobtainable, not merely approximate. Uniform `parent_origin` across the
call is fine and scales as below.

★ **`scale_to_target` becomes origin-aware.** `rescale_effects_to_target()`
(`define_additive_effects.R:534-541`) uses `V_A = Σ 2 p q a²`, where the 2 is the number
of independent copies contributing to an unparented additive value. A parent-qualified
term contributes **one** copy, so its variance is `Σ p q a²` and an imprinted model asked
for `target_add_var = V` lands at `V/2`. The contract becomes

```text
V_A = Σ_j  n_eligible,j · p_j q_j a_j²      n_eligible = 2 unparented, 1 parent-qualified
```

`n_eligible ∈ {1, 2}` is a *complete* enumeration only because
`assert_qtl_autosomal()` (`chr_meta_helpers.R:369`) already refuses
`scale_to_target = TRUE` at any locus that is not `(1,1)` for both offspring sexes. That
guard stays, and it is what keeps the unparented factor at exactly 2 — without it a
hemizygous locus would need a third value and a sex ratio. Each variant keeps its own
line-specific `p_j`. This is a latent bug in the trait-wide implementation rather than a
regression introduced here; scaling only ever applied to the *sampled* branch
(`define_additive_effects.R:281-291`), never to manual `effects`, which is why no current
test catches it. The rewritten writer owns the scope and must not carry it forward.

---

## Validation and transactions

One transaction per public write: delete-then-insert-then-validate-then-commit.

**Declared in SQL** (verified in DuckDB 1.5.5): closed `contrast_name` /
`line_match_type` sets; `contrast_name` ↔ `copy_count_value` / `dosage_value` /
`center_value` agreement including `dosage_value <= copy_count_value`;
★ `center_value IS NOT NULL` **and** `BETWEEN 0 AND 1` for non-indicator members;
`copy_count > 0`; `parent_origin IN (1,2)` or NULL;
`line_match_type` ↔ `line_name` agreement (`'exact'` ⇒ named, otherwise NULL);
★ `'any'` ⇒ `parent_origin IS NOT NULL`;
`locus_id` FK; member → effect FK;
composite origin → member FK.

**Enforced in R:**

| Group | Checks |
|---|---|
| ★ Input shape | `genome_value` and `effect_name` constant within a `term_id`; no unknown `terms` columns; `origin` keys resolve to a member of the named term; every message names the user's `term_id` |
| Canonicalization | members by `locus_id`; origin rows by the documented tuple |
| Families | no duplicate family + scope identity; **no overlapping-but-incomparable predicates within a family** |
| Terms | ≥ 1 member; each locus at most once per term |
| Origins | ≤ 1 origin row and `copy_count = 1` on additive members; exact-multiset satisfiability on genotype members — ★ **by full bijection search, not greedy consumption**: a demand set mixing a parent-qualified row with an ANY-parent row can be satisfiable while greedy fails it (Phase A finding; ≤ 2 items at diploidy, so the cost is nil). The same search decides genotype containment; ★ `'any'` rejected on genotype members (cross-table — `contrast_name` lives in `genome_effect_members`) |
| Ploidy | `dominance` rejected at non-diploid loci **unless proven on that member** |
| Owners | reserved owner names not writable by the general writer without override |
| Cross-table | `trait_name` in `trait_meta`; children deleted before parents |

---

## ★ Views

Phase B builds these; Phase E moves existing readers onto them. The base tables must
not force a six-way join on anyone.

| View | Grain | Columns |
|---|---|---|
| `genome_effect_terms` | one row per term | `trait_name`, `effect_owner`, `effect_name`, `id_genome_effect`, `effect_order` (member count, **derived** — never stored), `contrast_signature` (ordered `contrast_name` list), ★ `family_key`, `scope_description`, `genome_value` |
| `genome_effect_loci` | one row per (term × locus) | `trait_name`, `effect_owner`, `id_genome_effect`, `member_slot`, `locus_id`, `locus_name` (joined from `genome_meta`), `contrast_name` |
| `ind_tgv_total` | one row per (individual × trait) | `id_ind`, `trait_name`, `tgv_total` = `SUM(tgv_value)`, `replicate` |

`locus_name` lives only in the view — v4.1 removed it from
`genome_effect_members` to kill the id/name agreement invariant.

★ **`family_key` makes the one invisible concept visible.** Whether two definitions
**compete** (specificity fallback) or **sum** is decided by the family signature — which
appears in no table, no column, and no error message a user sees before they are
surprised by it. `family_key` is that signature rendered as a stable hash or a readable
composite string, derived in the view and never stored. With it, *"why did my
line-A effect not override the common one?"* is one query against
`genome_effect_terms` instead of a re-read of §Fallback families. Terms sharing a
`family_key` compete; terms with different keys sum. That sentence is the whole
mental model, and the view is what lets a user check it.

★ **Register the views, or `schema()` diverges between a fresh and a restored
population.** `pop$tables` is a curated list everywhere except `restore_pop.R:77`, which
takes it from `DBI::dbListTables()` — and that **includes views** (verified). So a
restored population would show three extra undescribed rows that a freshly built one
does not. Add all three to `SYSTEM_TABLES`, `.schema_table_order()`, and
`.all_schema_descriptions()` (`test-schema-print.R:96-104` asserts those three agree),
and append them to `pop$tables` at creation so both paths list the same set.

### ★ Causal-locus terminology (this is what `tidybreed_pop.R:159` needs)

The print method currently reports `n_qtl` from a `genome_effects` row count, which
stops being meaningful once one coefficient spans several loci. Precise definitions,
all answerable from `genome_effect_loci`:

- a locus is **causal** for a trait if it appears as a member of any term;
- an **additive QTL** appears in an order-1 term whose contrast is `additive`;
- an **epistatic-only** locus is causal with no order-1 term at all.

`print.tidybreed_pop()` reports the **causal locus count**, since that is the one that
never undercounts. A breakdown by the three categories belongs in `describe_table()` or
a summary helper, not the print header.

---

## Acceptance gates

Every gate carries a **hand-computed expected value** and is evaluated, not merely
written. Storage representability gates Phase A; numerical evaluation gates Phase D.

**Core semantics (1–14).**

1. Order-one copy-additive values, including the per-copy line fallback — reproduces
   today's engine exactly.
2. Hemizygous locus — centering uses realized copy count, not `2p`.
3. Exact `dominance` values at dosages 0, 1, 2, plus both orthogonality identities
   (`E[x_D] = 0`, `E[x_A x_D] = 0`).
4. Pairwise A×D, with per-member reduction before multiplication.
5. Three-way term.
6. Complete **and** sparse 3×3 indicator surface.
7. Functional `(a, d)` vs Cockerham `(α, δ)` — including the intercept.
8. `{A, NULL}` origin handling — no within-A masquerade, no double count.
9. A/B and reciprocal origin cases.
10. Canonical equality for reversed input member order, compared as **canonical query
    results**, not byte-identical rows — surrogate IDs and physical row order are not
    semantic contracts.
11. Duplicate logical-term rejection.
12. Transaction rollback after any invalid member leaves no partial write.
13. Specificity-tie rejection at write time.
14. **An injected dominance or epistatic term changes a computed genetic value.**

**Families, owners, and output (15–25).**

15. An additive main effect, a dominance main effect, and an interaction sharing a
    locus all contribute; origin precedence never makes them suppress one another.
16. Common and specific variants compete **only inside the same fallback family**.
17. Two-locus common/specific partial-origin case, hand-computed, proving the
    tuple-level fallback semantics.
18. `copy_count = 1` vs `2` matching is unambiguous for A/A and A/B.
19. Replacing generated additive effects preserves terms under unrelated effect owners.
20. Replacing a set removes old loci absent from the new set.
21. Every accepted term type maps to exactly one `component_name`, and the components
    sum to the derived total.
22. `add_tgv()` is idempotent per `(id_ind, trait_name, component_name)`.
23. A validated `dominance` term can never fail at evaluation because a valid sex has
    a different realized copy count.
24. Deleting or replacing terms leaves no orphan member or origin rows.
25. Functional `(a, d)` input has an expected result for centered **and** raw modes,
    including repeated append calls.

**Origin containment and states (26–39).**

26. Parent-qualified A overrides generic A — no equal-row-count tie.
27. Reciprocal-specific A/B overrides generic A/B in each direction.
28. Overlapping origin predicates where neither contains the other are rejected.
29. Multiple origin rows on one additive member are expanded or rejected explicitly.
30. Indicator states distinguish copy counts 0, 1, 2 at dosage zero.
31. The evaluator produces an indicator value for the zero-copy state despite no
    haplotype row existing.
32. Custom functional dominance changes `tgv_value` without changing `tbv_value`.
33. `add_tbv()` excludes interactions and non-reserved custom basis terms.
34. General-writer defaults cannot be deleted by rerunning
    `define_additive_effects()`.
35. Successive common, A, and B additive definitions preserve every scope variant with
    correct fallback.
36. The reported functional genetic mean is not automatically added again as a
    phenotype mean.
37. Dominance ploidy proof occurs on the dominance member itself (mixed-locus fixture:
    one member diploid-scoped, the dominance member variable-copy).
38. ★ Label-vector preflight warns or stops at the configured threshold.
39. ★ A hand-entered order-1 surface is written as `component_name = 'order1_other'`
    and still sums into the derived total.

**★ Imprinting, origin scope, and generated-effect scaling (40–46).**

40. ★ **Imprinting parity — the full fixture, not a reduced one.** Preserve
    `test-add_tbv.R:235-270`: Duroc and Landrace line-specific additive variants with
    distinct coefficients and per-line centers, one `('exact', line, parent)` origin row
    each, an F1 offspring, paternal expression — value for value against the independent
    per-copy oracle, not against stored pre-change output. A population-wide
    `('any', 1)`-only fixture would pass with the line dimension silently dropped, which
    is precisely the wrapper bug this revision fixes.
41. ★ **`'any'` constraints.** `('any', parent_origin NULL)` rejected through the writer
    **and** by direct SQL insert (it is a `CHECK`); `'any'` on a genotype member rejected
    through the writer (cross-table, so R-only).
42. ★ **Wrapper scope matrix.** All four `(line_name, parent_origin)` combinations
    round-trip to the stored scope in the Writer-API mapping table, and `replace_scope`
    on one leaves the other three untouched.
43. ★ **Origin-aware variance target.** A generated unparented model and a generated
    parent-qualified model each realize `target_add_var` in the base population; the
    parent-qualified one lands at `V`, not `V/2`.
44. ★ **Multi-trait origin contract.** Per-trait `parent_origin` resolves per trait in
    both the scalar-recycled and vector forms; a mixed-origin call supplying `G` is
    rejected, and the message gives the covariance reason.
45. ★ **No eligible copy.** Sex-chromosome fixture. A paternally-qualified X-linked term
    contributes 0 where the same trait also has a matching autosomal term; the same trait
    with the X term alone **errors** for males, who have no paternal X copy.
46. ★ **Required center.** Direct insertion of an `additive` or `dominance` member with
    NULL `center_value` fails at the SQL constraint — the case a bare
    `CHECK (center_value BETWEEN 0 AND 1)` accepted.

**★ Structural prerequisites (47–50).**

47. ★ **`genome_meta` primary key.** A duplicate `locus_id` insert fails, and a
    `genome_effect_members` row naming a nonexistent `locus_id` is refused by the FK —
    the FK that could not be declared at all before the PK existed.
48. ★ **Fresh vs. restored parity.** `schema()` lists the same table set after
    `open_pop() |> define_genome()` as after `restore_pop()` of that same file, views
    included, each with a description.
49. ★ **Pre-change database.** `restore_pop()` on a file carrying the old
    `genome_effects` shape stops with the term/member message instead of failing later
    inside `add_tbv()`.
50. ★ **Parent-only re-run.** `define_additive_effects(line_name = "A")` then the same
    call with `parent_origin = 1` leaves both variants with correct per-copy fallback
    **and** warns; the common/A/B sequence of gate 35 and a reciprocal
    `(exact A, 1)` / `(exact A, 2)` pair both stay silent.

**★ Writer input, evaluation shape, and inspectability (51–54).**

51. ★ **Evaluation does not regress.** Order-one additive TBV over a benchmark
    population (`dev/benchmarks/`, following `benchmark_haplotype_scale.R`) stays within
    a stated factor of the current single-statement engine, and the statement count is
    **independent of the number of individuals** — the structural assertion that catches
    a per-individual loop even on a small fixture where wall-clock would not.
52. ★ **Resolution is precomputed, not per individual.** The resolved variant map built
    for 10 individuals is identical to the one built for 10,000 of the same lines, and a
    mixed-line multi-scope fixture agrees value-for-value with the per-tuple semantic
    definition in §Contribution — proving the label-vector factorization is not merely
    faster but the same function.
53. ★ **`terms` round-trips.** All three worked examples in §Writer API write and read
    back to the canonical stored form; a malformed call (varying `genome_value` within a
    `term_id`, a locus repeated in one term, an unknown column) is rejected with a
    message naming the **user's** `term_id`, never an `id_genome_effect`.
54. ★ **`family_key` is honest.** Two scope variants of one term share a `family_key`;
    an additive main effect and a dominance main effect at the same locus do not; two
    owners defining the same common term do not. The view answers "compete or sum?"
    without reading the plan.

The oracle at `tests/testthat/test-add_tbv.R:16-40` stays committed — it recomputes
the additive formula from first principles and is not pre-change golden output.

---

## Implementation order

| Phase | Work | Gate |
|---|---|---|
| **A** ✅ | **Complete — `plans/update_genome_effects_phase_A.md`.** Fixtures **with hand-computed expected values** before DDL, including the origin truth table: common vs A-specific additive · generic A vs paternal-A · common vs A/B dominance · generic A/B vs both reciprocals · overlapping-incomparable rejection · common vs origin-specific A×A · partial specificity at one member of two · two disjoint specific combinations that must both contribute · absent / hemizygous-allele-0 / diploid-dosage-0 indicator states | ✅ Met: 19 fixtures, **no schema change required**; both evaluators agree with the hand computations. Containment order and multi-locus fallback settled before DDL. ★ Writing each fixture in the `terms` format moves to Phase C as a **round-trip** against this registry (gate 53) — it only becomes a real check once a writer exists to canonicalize it |
| **B** | ★ `genome_meta` gains its `PRIMARY KEY` and the `open_pop.R:286` DDL is deleted **in the same commit** that adds the three tables to `GENOME_TABLES` — neither works alone; effect tables move into `define_genome()`; 4 tables + 24 registry entries; SQL constraints; containment checker; R validator; views, registered in all three schema lists | Registry, constraint, and validator tests pass; gates 46–48, ★ 54 |
| **C** | `define_genome_effects()`; `define_additive_effects()` rebuilt on it with `replace_scope` and `parent_origin`; `(a,d)` and genotype-table helpers; ★ origin-aware `scale_to_target`; ★ parent-only re-run warning; ★ **delete `trait_meta.expressed_parent`** (see below) | Writer round-trips every Phase-A fixture; gates 34–35, 41–44, 50, ★ 53 pass |
| **D** | ★ **One** evaluator, built to §Evaluation strategy: label alphabet, resolved variant map, member reduction (incl. synthesized zero-copy state), family partitioner, containment resolver, label-vector preflight, term evaluator, `add_tgv()` writing `ind_tgv`. ★ `add_tbv()` becomes a **thin filtered call into that same evaluator** — reserved owner, order-1 `additive` variants only — not a second implementation | Oracle agrees; gates 1–39, 40, 45, ★ 51–52 pass |
| **E** | Delete the old table shape; ★ add the `restore_pop()` guard for pre-change files; move QTL extraction (`tidybreed_pop.R:159`) and `extract_genotypes()`'s nominal `table_name` check (`:124-127`) to the locus view | No legacy columns remain; gate 49 |

★ **Phase D builds `add_tbv()` once, not twice.** Earlier revisions "rewired
`add_tbv()` to the reserved owner" as its own piece of work, which
`plans/consolidate_genetic_values.md` then deletes outright — throwaway code, and by
this package's own pre-1.0 policy a compatibility shim in everything but name. Phase D
instead builds a single evaluator and expresses `add_tbv()` as a filter over it:
reserved owner, order-1, `contrast_name = 'additive'`. Nothing is discarded when
consolidation later removes the wrapper, because the wrapper was never the substance.

`ind_tbv` and `add_tbv()` therefore **survive this plan with their current meaning
intact**. Carrying non-additive genetic values into them is the next plan's work,
alongside the simulation of non-additive QTL effects; this plan only has to stop
`add_tbv()` reading a table shape that no longer exists. The redundancy that leaves —
`ind_tbv.tbv_value` beside `ind_tgv`'s `'order1_additive'` — is knowingly accepted and
time-boxed; see §Output contract and Q1.

★ **Deleting `expressed_parent` (Phase C).** Per `CLAUDE.md` pre-1.0 policy the old
name goes completely — no shim, no alias.

| Site | Change |
|---|---|
| `define_trait.R:177` | drop the column from the `trait_meta` DDL |
| `define_trait.R:65,76,120` | drop the argument, `match.arg`, and the written value |
| `add_tbv.R:159,191-194` | drop the `SELECT` of the column and the `parent_filter` branch |
| `add_tbv.R:26,242` · `define_trait.R` roxygen | drop the imprinting prose; point at `parent_origin` |
| `sql_utils.R:100` · `schema.R:257` | drop from `TABLE_RESERVED_COLS` and the column description |
| `man/add_tbv.Rd`, `man/define_trait.Rd` | regenerate |
| `helper-parity.R`, `test-add_tbv.R`, `test-define_trait.R` | retarget to origin rows (gate 40) |
| `CLAUDE.md:130, 406, 921, 1027` | ★ **four** sites, not two: the `trait_meta` schema table, the `define_trait()` argument list in **Two-Layer Phenotype Design**, the `define_trait()` bullet under **Implemented Functions**, and the `add_tbv()` imprinting paragraph |
| ★ `vignettes/tidybreed-introduction.Rmd:112` | drop `expressed_parent` from the `define_trait()` argument table |
| ★ `vignettes/swine/swine-time-based-age-at-puberty-sex-semen.R:1406, 1439` | two **executable** `expressed_parent = "both"` calls that will error once the argument is gone |

Historical plan files and the historical `NEWS.md` entry stay as written — they describe
revisions that existed. Live documentation does not get that exemption.

The same phase carries the two writer contracts the imprinting move exposed:

| Site | Change |
|---|---|
| `define_additive_effects.R:152` | ★ add `parent_origin`, per trait; compose with `line_name` per the Writer-API mapping; reject mixed origins under `G` |
| `define_additive_effects.R:534-541` | ★ `rescale_effects_to_target()` takes `n_eligible` per locus; `assert_qtl_autosomal()` (`chr_meta_helpers.R:369`) stays and is what bounds it to `{1, 2}` |

It lands in Phase C rather than Phase E so that Phase D's `add_tbv()` rewire has no
legacy flag left to honour.

**Out of this plan, contract defined here:** wiring non-additive values into
`add_phenotype()` and activating the dead `phenotype_components.genome_effect_types`
(`open_pop.R:334`).

---

## Disposition of the v4.1 review

### Accepted (8 issues, 5 minor corrections, 13 gates)

| # | Finding | Resolution |
|---|---|---|
| 1 | **Origin-row count is not a valid specificity order** — silently disabled reciprocals | Predicate containment over two small finite lattices, fully specified with worked cases; ≤ 1 origin row per additive member |
| 2 | **Dosage alone is not a genotype state** — the indicator escape hatch was broken exactly where v4.1 sent users | `copy_count_value` column; state is `(copy_count, dosage)`; in the family signature; zero-copy state synthesized by the evaluator |
| 3 | **`add_tbv()` cannot infer TBV from arbitrary additive terms** | Reserved `generated_additive_tbv` owner; `add_tbv()` reads only order-1 additive variants from it |
| 4 | **Shared `'default'` permitted accidental deletion**, and per-line replace would delete sibling variants | Reserved vs `custom` defaults; reserved-name protection; **`replace_scope`** mode using predicate equality |
| 5 | Reported mean must not be copied into `phenotype_meta.mean` | Instruction removed; `μ` reported, never written; running total narrowed to single-locus main effects |
| 6 | Displayed DDL omitted the promised `locus_id` FK | Declared (**confirmed missing** — the DDL carried it only as a comment). ★ v4.8: declaring it also requires a `PRIMARY KEY` on `genome_meta`, which the table never had — see §Schema. `ind_tgv.id_ind` stays R-enforced, matching `ind_tbv` |
| 7 | Dominance ploidy proof must be member-specific | Stated explicitly; fixture 37 |
| 8 | Tuple enumeration needs a resource guard | Preflight estimate, configurable warn/hard-cap; fast path re-described as **common/unscoped**, not "additive-only" |

Minor: title carries the revision · identifier validation on owner/effect/line names ·
zero-row component behavior defined · `require_complete` over `(copy_count, dosage)`.

### Rejected (1)

**Minor correction 1 — "Phase E appears twice in the implementation table."** It does
not. `grep '^| \*\*[A-E]\*\*'` returns exactly five rows, A through E, each once. No
change made.

### Addition of my own

★ **The family signature cannot distinguish two additive terms with different
`center_value`.** A functional `additive`@0.5 and a Cockerham `additive`@`p` at one
locus under one owner have identical signatures and, if both common, collide as a duplicate
identity. They *are* combinable, so the writer rejects **with the combined equivalent
in the message** rather than a bare duplicate error. Raised as Q2.

### Settled by this review

Q2 (tuple-level resolution) and Q3
(functional default) from v4.1 are confirmed and closed. v4.1's Q1 (mechanical
components + derived total) is **withdrawn** — see the v4.3 revision note.

---

## Decided

### Imprinting — `trait_meta.expressed_parent` is deleted (v4.5 Q1, accepted 2026-09-06)

Imprinting becomes an **origin predicate**: an `additive` term with one
`('any', parent_origin, copy_count = 1)` origin row, matched per copy. Reasons:

1. **It matches the biology.** Imprinting is locus-specific — IGF2 in pigs is one
   imprinted locus, not a whole trait. A trait-wide flag cannot express that; an
   origin row can, per locus and per effect owner.
2. **It removes a mechanism instead of adding one.** Parent-of-origin was expressible
   in two places once `genome_effect_member_origins.parent_origin` existed — the same
   double-specification ambiguity v2 identified and closed for states-vs-scopes.
3. **The v2 incoherence resolves structurally.** `'any'` is rejected on genotype
   members, so a `dominance` contrast cannot be restricted to one parent's copy —
   correct, since a genotype needs both. No validator special case is required.
4. **Small blast radius:** functionally `add_tbv.R` and `define_trait.R`, plus two
   registries, two man pages, three tests, four `CLAUDE.md` sites, and two live
   vignettes — one argument table and two executable calls. Full inventory in Phase C.

**Consequence, and a bug this surfaced:** `line_match_type = 'any'` moves from
out-of-scope into the closed set. v4.5's lattice already used `ANY` as a line value and
its worked case `(ANY, 1)` was **unrepresentable in v4.5's own DDL**, since `ANY` was
reachable only as zero origin rows and a row is where `parent_origin` lives. Guarded by
two rules (`'any'` requires a non-NULL `parent_origin`; `'any'` is additive-members
only) so it cannot become a second spelling of the common scope.

★ **The wrapper composes, it does not stamp.** `define_additive_effects()` maps
`(line_name, parent_origin)` onto one origin row (Writer API); writing `'any'`
unconditionally would drop the line dimension and collapse the existing Duroc/Landrace
imprinted fixture onto a single scope. Numerical parity with the removed flag is gate 40,
against that full fixture rather than a population-wide reduction of it.

★ **One more thing the move owns.** Scope now lives on the effect rows, so the generated
writer — not `trait_meta` — is what knows how many copies a term is eligible for, and
`scale_to_target` must use that count. Carrying `2pq` forward would keep an imprinted
model at half its requested variance; see the Writer API and gate 43.

### ★ Scope of the effect model — settled in v4.9

Four decisions, recorded here so they are not re-litigated during implementation. None
adds a table, a column, or a resolution rule.

1. **Origin predicates stay copy-level.** `line` and `parent` describe an allele copy.
   Individual-level conditions — sex, environment, age, any covariate — do **not** go on
   the term, because that would make containment a product of copy-level and
   individual-level lattices. Sex-dimorphic and G×E genetics are modelled as **correlated
   component traits** composed at the phenotype layer, which the two-layer design already
   supports and which is what quantitative genetics does anyway. Future limitations 1
   and 4.
2. **Nonlinearity belongs to the phenotype layer.** The effect model is linear in its
   coefficients by construction; saturating or thresholded maps compose a genetic value
   here and transform it in `phenotype_meta.type = "derived_formula"`. Keeping it out is
   what preserves the set-based evaluator. Future limitations 5.
3. **`add_tbv()` is built once**, as a filter over the single evaluator — not rewired
   now and deleted next plan. `ind_tbv` keeps its current meaning through this plan;
   non-additive values in it are the next plan's work.
4. **`effect_set_name` → `effect_owner`.** The column never selected anything; it
   scoped replacement. See §Effect owners.

---

## Open questions

Neither blocks a phase; both are ergonomic and revisable after implementation.

### Q1 — Should `add_tbv()` warn when a trait has non-additive terms?

v4.2 asked what `ind_tbv` means once dominance exists. That was the wrong question:
a TBV is additive by definition and its meaning does not change. `add_tbv()` reads the
reserved additive set; `add_tgv()` computes `G`; the two never collide.

What survives is a narrower ergonomic question — and it now matters for longer, since v4.9 keeps `ind_tbv` alive through this plan rather than consolidating inside it. A user who defines a dominance model
and then selects on `ind_tbv` gets the additive component — correct, but possibly not
what they intended, since nothing forces them to call `add_tgv()`.

| Option | Notes |
|---|---|
| **Warn once from `add_tbv()`** ← recommended | When the trait has terms outside the reserved additive owner, note that `tbv_value` is the additive component only and point at `add_tgv()`. One message, no behaviour change |
| Silent | Defensible — TBV is doing exactly what it says. But the user has to already know they need a second function |
| Error unless acknowledged | Too aggressive; selecting on breeding value under a dominance model is a legitimate and common choice |

Low stakes. Flagging it only because it is the one place a user can get a correct
number that answers a question they did not mean to ask.

### Q2 — Mixed coding at one locus: reject, or allow and sum?

| Option | Notes |
|---|---|
| **Reject, with the combined equivalent in the message** ← recommended | Preserves the duplicate-term protection. The message tells the user exactly what single term to write instead. Costs nothing except that a user cannot stack codings even deliberately |
| Add `center_value` to the family signature | Different centerings become different families and sum. Mathematically fine, but the duplicate-term guard disappears at that locus and a user can silently stack two additive terms |
| Forbid mixing codings within a `(trait, owner)` entirely | Simplest to explain; blocks legitimate per-locus choices |

**Recommendation: reject with guidance.** Low stakes and easily revisited — flagging
it because it is a user-visible error message you will have to defend, and because the
review did not surface it.

---

## ★ Future limitations — recorded, not scheduled

Six entries: five things this design deliberately cannot express, and one thing it can
store but cannot yet help a user calibrate. Written down so each is recognized as a
known deferral rather than rediscovered as an architecture question. **None is scheduled
here and none blocks 1.0.0.** Each names the *shape* of the eventual fix, because that
is the part worth deciding once, cheaply, and in advance.

### 1. GxE and covariate-varying effect sizes — **revisit at 2.0.0 (plant breeding)**

A QTL whose coefficient depends on environment, year, location, or age has nowhere to
live: `genome_value` is one scalar per term.

**Not needed for the current livestock work.** Nucleus-vs-commercial GxE is already
modelled the way quantitative genetics conventionally models it — as **two correlated
traits** with a genetic correlation below 1, which the two-layer trait/phenotype design
handles today with no new machinery. Plant-breeding programs (multi-environment trials,
explicit E×QTL) are where the trait-pair device stops scaling, and that is the 2.0.0
conversation.

**Shape of the fix:** a term-level `weight_type` / covariate column on `genome_effects`,
mirroring `phenotype_components.weight_type` (`"fixed"`, `"covariate"`, `"legendre"`).
`ALTER TABLE` plus an evaluator branch. It does **not** enter the origin table — an
environment is a property of a record, not of an allele copy.

### 2. Multiallelic loci — **not needed; blocked upstream anyway**

The member state is `(copy_count, dosage-of-allele-1)`. Real multiallelic effects need
allele identity, and nothing in the package can produce a third allele: `ind_haplotype`
stores `allele UTINYINT` as 0/1 and `add_dosage()` computes `SUM(h.allele)`
(`add_dosage.R:126-130`), a count, not a genotype.

**Shape of the fix:** the state tuple widens — an `allele_id` column on
`genome_effect_members` with the state redefined as `(copy_count, allele_id, count)`, or
a per-member allele-count child table. This is the plan's **one acknowledged reshape**,
and it is correctly deferred: when `ind_haplotype` gains an allele dimension, effect
storage was going to change regardless, and that is a far larger project than this one.

### 3. Phase / cis-trans / haplotype-block effects — **noted; the expensive deferral**

"Both variants on the same strand" is unsayable. §Aggregation fixes that each member
reduces to **one value per individual** before members multiply, which is genotype-level
by construction; and origin predicates describe an allele's *line and parent*, never its
*strand*. Compound heterozygotes, cis-regulatory pairs, and MHC-style haplotype blocks
are all outside this.

**Shape of the fix, and why it is the costly one:** a phase constraint is a filter
**across** members of a term, so it fits neither the per-member origin table nor the
per-member reduction contract that makes set-based evaluation possible. It would need a
term-level constraint (a new small table or column pair) *and* an evaluator that keeps
tuples paired by strand instead of reducing members independently. Expect this to cost a
redesign of §Evaluation strategy rather than an `ALTER TABLE`. Not a reason to build it
now — a reason to know the price before promising it.

### 4. Sex-specific effect sizes — **out of scope; the answer is trait decomposition**

A QTL with a different coefficient in males and females has no home, and users will
look for one here because `genome_map` and `chr_inheritance` both carry a `sex`
dimension.

**This is a decision, not a gap.** Origin predicates are **copy-level** (line, parent);
sex is **individual-level**. Putting an individual-level dimension on the term would
make the containment lattice a product of copy-level and individual-level predicates —
mechanically an `ALTER`, conceptually a doubling of the resolution surface that §Origin
resolution spent four revisions getting right.

**The package already answers it:** sex-dimorphic genetics is two component traits
(`ADG_male`, `ADG_female`) carrying a genetic correlation in `trait_var_comp`, observed
through one phenotype per sex whose `phenotype_meta.expressed_sex` decides which
individuals get a record from which component. That is the same device as GxE above, it is already implemented, and it is what
quantitative genetics does anyway. If a future program needs sex-specific effects *on
one term* rather than as a trait pair, that is its own plan — **do not re-open the
origin lattice for it.**

### 5. Nonlinear genotype→value maps — **belongs to the phenotype layer, permanently**

Indicator completeness makes any *lookup* function representable, but a smooth or
saturating function (a cap, a ratio, a sigmoid of a multi-locus sum) is **dense** in the
indicator basis — representable in principle, useless in practice, since the row count
is the size of the genotype space.

**The escape hatch already exists and should stay where it is:** compose a linear
genetic value here, then transform it at the phenotype layer via
`phenotype_meta.type = "derived_formula"` and `formula_helpers.R`. Nonlinearity must not
leak into `genome_effects`; keeping it out is what preserves the set-based evaluator.

### 6. Realized variance is not reportable — **the first thing to build after this plan**

Variance targeting for non-additive terms is out of scope (below), which means this plan
ships *storage* for dominance and epistasis with no path from a user's `d²` to a set of
coefficients — and, worse, no way to find out what variance a model actually has.

**The cheap 80%:** not targeting, but a **diagnostic** — a helper that reports realized
`V_A` / `V_D` / `V_I` in the base population from stored effects plus founder allele
frequencies, stating its linkage-equilibrium assumption in the output. That turns
"storable" into "usable" without committing to inverting the variance function. Worth
scheduling immediately after `plans/consolidate_genetic_values.md`.

---

## Explicitly out of scope

- Polyploid and multiallelic **simulation** (see Future limitations 2).
- Lifting the `from_parent_1 + from_parent_2 <= 2` CHECK at `define_genome.R:368`.
- Phased haplotype-combination effects (Future limitations 3).
- ★ Sex-specific, environment-specific, or covariate-varying effect *sizes* on a term —
  individual-level scoping stays out of the origin lattice (Future limitations 1 and 4).
- ★ Nonlinear genotype→value maps; these belong to the phenotype layer
  (Future limitations 5).
- Mutation and any change to `make_gametes.cpp`.
- Estimating non-additive effects (`add_ebv()` / BLUPF90 remain additive).
- ★ **Simulating** non-additive QTL effects, and carrying non-additive values into
  `add_tbv()` / `ind_tbv` — in progress separately and scheduled after this plan; this
  plan owns storage and evaluation only.
- Wiring `target_add_mean` into an effect-model intercept; line- or cross-specific
  intercepts.
- Phenotype integration and `phenotype_components.genome_effect_types` activation.
- Variance targeting for non-additive terms — the LE assumption in
  `V_term = coefficient² × ∏ Var(contrast)` is load-bearing (but see Future
  limitations 6 for the diagnostic that should exist).
