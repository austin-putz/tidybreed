# Genome Effects v3 — the parsimonious design

**Status:** design proposal, nothing implemented.
**Revision:** v3 (2026-09-05), a parsimony review of `plans/update_genome_effects.md` (v2).
**Verdict on v2:** the problem statement is correct; the schema is roughly 8× larger
than the problem requires. v2 proposes 16 tables. This proposes **2** (net +1 over
today), and shows every v2 acceptance fixture is still representable.

---

## The disagreement in one paragraph

v2 is right that `genome_effects` cannot express epistasis, right that this must be
fixed pre-1.0, and right about most of what it rejects. Where it goes wrong is the
step from *"a coefficient must span multiple loci"* to *"therefore effect storage
needs a normalized state space, a contrast-function engine, an immutable reference
registry, an origin-matching scope language, and a parallel surface representation."*
Those are five independent generalizations. Only the first is forced. The other four
are pre-built for biology the rest of the package cannot currently simulate, and each
is addable later **without reshaping any existing row** — which is the only cost
`CLAUDE.md`'s Schema Design Bias actually tells us to avoid.

---

## The decision rule

> **Pre-build a dimension only when its absence would force existing rows to change
> shape. Everything else is a later child table, which is additive, not a rewrite.**

`CLAUDE.md` says: *"It is enough to reserve clean dimensions now when the schema
would be painful to alter later."* Painful is the operative word. Adding a nullable
column or a new child table keyed to an existing PK is not painful — this package
does it routinely (`define_founder_haplotypes()` adds `founder_allele_freq` by
`ALTER TABLE`; `genome_map` was added alongside `genome_meta` rather than inside it).

Applying the rule to v2's five generalizations:

| Generalization | Forces existing rows to change shape? | Verdict |
|---|---|---|
| Multiple loci per coefficient | **Yes.** Every row splits; every reader rewrites. | **Build now** |
| Arbitrary gene action / genotype tables | No — see "indicator completeness" below; it needs *zero* new structure | Build now, free |
| Origin-composition scopes | No. A child table on `(id_genome_effect, member_slot)` | Defer |
| Multiallelic / normalized states | No. And blocked upstream regardless | Defer |
| Reference-population provenance | No. One optional column, or a later table | Defer |

---

## Audit corrections to v2

Verified against the tree at `4e9a72b`.

**1. Multiallelic is blocked upstream, not by this table.**
`R/define_genome.R:302` — `ind_haplotype(… allele UTINYINT)`, and `genome_meta` has
no allele dimension at all. Nothing in the package can *produce* a third allele.
v2 spends `genome_states`, `genome_state_alleles`, `genome_effect_reference_alleles`,
and `genome_contrast_values`' state grain on multiallelic-readiness. That is a roof
over a foundation that does not exist. When `ind_haplotype` gains an allele
dimension, effect storage will need to change too — and it will be changing anyway,
because that is a far larger project.

**2. The per-allele-copy path is the one semantic that must survive.**
`R/add_tbv.R:206-227` sums `(h.allele − COALESCE(e.base_allele_freq,0)) * e.genome_value`
over `ind_haplotype`, matching `e.line_name = h.line_origin` with a per-locus
`NOT EXISTS` fallback to the common row. v2 identifies this correctly and it is the
single genuinely subtle thing in the current design. v3 preserves it verbatim — same
join, same fallback, same `line_name` column, same meaning.

**3. `COALESCE(base_allele_freq, 0)` is a real bug-shaped invariant.**
v2 is right. v3 fixes it by requiring `base_allele_freq` non-NULL for every
frequency-dependent contrast, validated at write time. This needs no new table.

**4. Registry cost is 6 edits per table.**
`R/sql_utils.R:85` (`TABLE_RESERVED_COLS`), `:135` (`TABLE_PRIMARY_KEYS`),
`:161` (`TABLE_ROW_KEYS`), `:252` (`SYSTEM_TABLES`); `R/schema.R:671`
(`.schema_table_order()`) and `:100-117` (`.genome_descriptions()`). Guarded by
`test-schema-registries.R` and `test-schema-print.R`. **16 tables ≈ 96 edits; 2
tables ≈ 6.** v2 acknowledges this and proposes a registration helper to make the
cost bearable. Building tooling to absorb the cost of unnecessary tables is the
wrong end of the problem.

**5. Blast radius (unchanged from v2, confirmed).** 12 R files, 8 test files
including `helper-parity.R`, 9 man pages, 2 vignettes.

---

## The insight v2 walked past: indicator completeness

Any function over a finite product of per-locus genotype state sets

```
f : S₁ × S₂ × … × Sₘ → ℝ
```

is *identically* a sum of indicator terms:

```
f(g₁,…,gₘ) = Σ_cells  value(cell) × ∏ⱼ 1[gⱼ = cell.stateⱼ]
```

with one term per non-zero cell. This is not an approximation and not a modelling
choice — it is what a lookup table *is*.

Consequences, each of which deletes tables from v2:

- **A "genotype surface" is not a second representation.** A hand-entered 3×3 table
  is 9 coefficient terms with 2 members each. v2's four surface tables
  (`genome_effect_surfaces`, `_members`, `_cells`, `_cell_states`), its
  `value_semantics` tagged-union discussion, and its `missing_state_action` policy
  all dissolve: a cell you did not enter is a term you did not write, contributing
  zero, which is the only consistent reading of a component-valued table.
- **A contrast-definition engine is redundant.** v2's `genome_contrasts` +
  `genome_contrast_values` exist to make storage self-describing — so a restored
  database doesn't depend on code remembering what `'additive'` meant. But
  `contrast = 'indicator'` plus a genotype value **is** self-describing, with no
  registry, no coding version, and no reference FK in the contrast key. Any
  user-defined basis, at any ploidy, is a linear combination of indicators.
- **Polyploid non-additive degrees of freedom are free.** v2 devotes a section to
  the fact that a locus with *k* copies has *k−1* non-additive DoF and that
  polyploid bases are basis-dependent. Indicators over `genotype ∈ 0..k` span every
  such basis without naming any of them.

The named contrasts `'additive'` and `'dominance'` survive only as *compressions* —
one row instead of three — and only because they let `rescale_effects_to_target()`
stay analytic. They are conveniences, not the mechanism.

---

## Schema

```sql
CREATE TABLE genome_effects (
  id_genome_effect INTEGER PRIMARY KEY,   -- next_int_id()
  trait_name       VARCHAR NOT NULL,      -- FK to trait_meta.trait_name
  model_name       VARCHAR NOT NULL DEFAULT 'default',
  line_name        VARCHAR,               -- NULL = common; else allele-copy line origin
  effect_class     VARCHAR NOT NULL,      -- display/filter only, never identity
  genome_value     DOUBLE  NOT NULL
);

CREATE TABLE genome_effect_members (
  id_genome_effect INTEGER  NOT NULL,
  member_slot      INTEGER  NOT NULL,     -- assigned by ascending locus_id, never input order
  locus_id         INTEGER  NOT NULL,     -- FK to genome_meta.locus_id
  locus_name       VARCHAR  NOT NULL,     -- denormalized, matches ind_haplotype convention
  contrast         VARCHAR  NOT NULL,     -- 'allele'|'indicator'|'additive'|'dominance'
  genotype         UTINYINT,              -- required iff contrast = 'indicator'
  base_allele_freq DOUBLE,                -- required for frequency-dependent contrasts
  PRIMARY KEY (id_genome_effect, member_slot)
);
```

```text
term contribution = genome_value × ∏ over members of (contrast value at that member)
```

### `contrast` — a closed set, two required values

| Value | Evaluation unit | Meaning | Needed for |
|---|---|---|---|
| `'allele'` | one `ind_haplotype` row | `allele − base_allele_freq` | **Required.** Today's additive path; the only unit that resolves per-copy `line_origin` for crossbreds |
| `'indicator'` | collapsed local genotype | `1` if genotype = `genotype`, else `0` | **Required.** Complete basis for all other gene action |
| `'additive'` | collapsed local genotype | `dosage − 2·base_allele_freq` | Optional 1-row compression |
| `'dominance'` | collapsed local genotype | HWE Cockerham `x_D` from `base_allele_freq` | Optional 1-row compression |

The two evaluation units are v2's `evaluation_unit` — kept, because the distinction
is real and load-bearing, but expressed as a column value rather than a table.

### Invariants (R-validated, one transaction per write)

1. A term has ≥ 1 member; `effect_class = 'intercept'` has exactly zero (empty
   product = 1).
2. A locus appears at most once per term.
3. `member_slot` is assigned **after** canonicalizing members by ascending
   `locus_id`. Supplying loci in either order stores byte-identical rows;
   swapping which locus carries dominance stores a different term.
4. Term identity is `(model_name, trait_name, line_name, ordered (locus_id, contrast,
   genotype) tuples)` — never `effect_class`.
5. `effect_class` is materialized for filtering and display, validated against the
   ordered member contrasts, and is **never** the source of a term's formula. In
   particular it does **not** promise that an `'additive'`-classed coefficient is a
   breeding value under epistasis (v2's §3.10 correction, adopted).
6. `genotype` is non-NULL iff `contrast = 'indicator'`.
7. `base_allele_freq` is non-NULL for `'allele'`, `'additive'`, `'dominance'`.
   No `COALESCE(…, 0)` anywhere.
8. No stored `effect_order`; member count is authoritative and a view exposes it.
9. `INTEGER` for `member_slot` — no artificial 255 ceiling.
10. Model replacement is explicit: `DELETE … WHERE trait_name = ? AND model_name = ?`
    inside the write transaction. Today's accidental "delete everything for this
    trait/line" becomes a named operation.

### Views

- `genome_effect_terms` — one row per term: trait, model, order, class, coefficient,
  human-readable member description.
- `genome_effect_loci` — one row per `(model, effect, locus)`; the causal-locus set.

QTL terminology becomes precise, changing `R/tidybreed_pop.R:159`: a locus is
**causal** if it appears in any term; an **additive QTL** appears in an order-1
`'allele'`/`'additive'` term; an **epistatic-only** locus is causal with no order-1
row.

`extract_genotypes()` must accept any filtered table carrying `locus_id`, replacing
the nominal `table_name == "genome_effects"` check at `R/extract_genotypes.R:124-127`.
(Adopted unchanged from v2.)

---

## Acceptance fixtures

Every fixture v2 lists, in v3 storage. These are the Phase A tests.

| Fixture | Rows |
|---|---|
| **Today's additive QTL** | 1 term (`effect_class='additive'`, `line_name=NULL`, value `a_j`) + 1 member (`contrast='allele'`, `base_allele_freq=p_j`). Summed over both haplotype rows gives `(dosage − 2p_j)·a_j` — **numerically identical to the current engine**, same SQL shape |
| **Line-A additive + common fallback** | Two terms, `line_name='A'` and `line_name=NULL`, each 1 member with its own `base_allele_freq`. Existing `NOT EXISTS` per-locus fallback unchanged |
| **Diploid dominance** | 1 term + 1 member `contrast='dominance'`; or 3 indicator terms if the user supplies raw values |
| **Pairwise A×D** | 1 term + 2 members, canonicalized by `locus_id`; lower-ID member `'additive'`, higher-ID `'dominance'` |
| **Three-way A×A×A** | 1 term + 3 members. No schema change from the pairwise case |
| **Arbitrary 2-locus 3×3 table** | ≤ 9 terms × 2 members, all `contrast='indicator'`. The day-one "type in your own gene action" feature |
| **Intercept** | 1 term, 0 members |
| **Autotetraploid locus** | indicator terms over `genotype ∈ 0..4`. No new columns, no reserved labels |
| **Model replacement** | Explicit DELETE by `(trait_name, model_name)` |
| **Canonicalization** | Loci supplied in either order → identical rows; dominance on the other locus → distinct term |

---

## What is deferred, and the retrofit path for each

This is the part that must be argued, not asserted. For each deferred capability:
what it costs to add later, and whether it reshapes anything.

### Origin-composition scopes (A/B-specific dominance, reciprocals)

**The genuine loss.** v3 keeps `line_name` on the term with today's meaning
("applies to allele copies of this line origin"). That expresses common, line-A, and
line-B effects. It does **not** express A/B-pair-specific dominance, reciprocal A/B
vs B/A, or per-member composition in an epistatic term.

**Retrofit:** one child table.

```sql
CREATE TABLE genome_effect_member_origins (
  id_genome_effect INTEGER  NOT NULL,
  member_slot      INTEGER  NOT NULL,
  origin_slot      INTEGER  NOT NULL,
  line_match_type  VARCHAR  NOT NULL,   -- 'exact'|'any'|'unknown'
  line_name        VARCHAR,
  parent_origin    UTINYINT,
  copy_count       INTEGER  NOT NULL
);
```

Keyed to `(id_genome_effect, member_slot)`, which already exists. Adding it changes
**zero existing rows** and **zero existing readers** (a term with no origin rows
behaves exactly as before). It is 6 registry edits and a validator extension.

**Why deferring is safe and deferring is correct:** v2 itself defers the matching and
precedence *resolver* (its Decision D3) — so v2 would ship scope tables that no code
can interpret. Shipping uninterpretable storage is worse than shipping none: it
invites users to write rows with no defined meaning. Build the tables when the
resolver is designed, in the same plan, with v2's truth table (common / A / B / A-B /
reciprocal / unknown / haploid / diploid / tetraploid) as its gate.

**Partial coverage available today:** `parent_origin`-specific effects are already
reachable through `trait_meta.expressed_parent` for the additive case.

### Multiallelic and normalized states

**Retrofit:** when `ind_haplotype.allele` gains a real allele dimension, add
`genome_effect_member_alleles` alongside the integer `genotype`, exactly as above.
But note that project will touch `genome_meta`, `ind_haplotype`, `ind_genotype`,
`add_dosage()`, `extract_genotypes()`, and `make_gametes.cpp`. Effect storage is a
rounding error in it. **Reserving for it now buys nothing and costs 4 tables.**

### Reference-population provenance

**Retrofit:** add `base_name VARCHAR` to `genome_effect_members` (one `ALTER TABLE`),
or a `genome_effect_references` table with a nullable FK. Neither reshapes a row.

**Why v2's version is over-built:** the immutability v2 wants is already achieved —
`base_allele_freq` is *copied into the member row at write time*, so it cannot drift
when the population changes. That is stronger immutability than an FK to a snapshot
table, not weaker. What the snapshot table adds is *provenance* (which generation,
which replicate, which filter). `CLAUDE.md` design principle 6 — "disdain and
intolerance for storing metadata" — argues directly against three tables for it.

### `value_semantics` (`'component'` vs `'total_genotypic'`)

Dissolved. Under indicator terms every row is a component summed into the model
total; "the whole genotypic value" is what you get by summing all of them. There is
no ambiguity to declare.

### `missing_state_action`

Dissolved. A cell you did not write is a term that does not exist and contributes
zero. If a user wants completeness enforced, that is a **write-time validation
option** on the writer call, not persisted state.

---

## Cost comparison

| | v2 | v3 |
|---|---|---|
| Tables | 16 (15 built, 1 reserved) | 2 (net **+1** over today) |
| Registry edits | ~96 | 6 |
| New registration tooling required first | Yes (v2's own recommendation) | No |
| Joins to read the additive path | 5–6 | 1 |
| Tables a user must understand to hand-enter a 2-locus table | 8 | 2 |
| Concepts documented | model, reference, reference-allele, state, state-allele, contrast, contrast-value, scope, scope-origin, term, member, surface, surface-member, cell, cell-state | term, member |
| Additive TBV numerically unchanged | Yes | Yes |
| Epistasis to any order | Yes | Yes |
| Arbitrary hand-entered gene action | Yes (4 tables) | Yes (0 extra tables) |
| Cross-specific dominance **storable** | Yes | Deferred (1 child table) |
| Cross-specific dominance **evaluable** | No (resolver deferred) | No (same) |

The last two rows are the whole trade. v2 buys the ability to *store* rows that no
code can interpret, at a cost of 14 tables and ~90 registry edits.

---

## Implementation order

| Phase | Work | Gate |
|---|---|---|
| **A** | **Fixtures before DDL.** Hand-build rows in a scratch DuckDB for every fixture in the table above | Every fixture representable with no serialized key and no special column. If one is not, stop |
| **B** | Create the two tables + 6 registry entries each; validator; the two views | `test-schema-registries.R`, `test-schema-print.R` pass |
| **C** | Write the new shape through `define_additive_effects()`: default model, `contrast='allele'`, `base_allele_freq` per member, order-1 canonical terms | Coefficients preserved exactly |
| **D** | Switch readers. `add_tbv.R` order-1 query joins `genome_effect_members` for `contrast='allele'`; verify against the existing independent oracle at `tests/testthat/test-add_tbv.R:16-40` | Oracle agrees |
| **E** | Delete the old table shape; move QTL extraction (`tidybreed_pop.R:159`) and `extract_genotypes()` (`:124-127`) to the locus view | No legacy columns remain |

No staging view is needed — the additive read path is a single join in both shapes,
so D is a query rewrite rather than a compatibility layer. This removes v2's Phase D
scaffold and the `CLAUDE.md` shim concern with it.

### Testing

Adopted from v2, minus the tests for tables that no longer exist:

- Each Phase A fixture becomes a test that writes rows and reads them through the views.
- **Step-C invariance is a development-time check, not a committed test.** Run the
  existing oracle on the same seed before and after. **Do not commit golden output**
  — `CLAUDE.md`'s forward-looking reproducibility contract forbids it.
- **Canonicalization:** either input locus order → byte-identical rows; dominance on
  the other locus → distinct term.
- **Indicator completeness:** a hand-computed 3×3 surface checked cell by cell, and
  its `(a, d)` coefficient-form equivalent, must agree.
- **Coding equivalence:** functional `(a, d)` at base `p` converted to Cockerham
  `(α, δ)` with `α = a + d(q − p)`, `δ = d` gives identical genotypic values at all
  three states.
- **NULL-origin regression:** an individual with `{A, NULL}` line origins must not be
  treated as within-A (v2's identified bug, still worth a test since `line_name`
  semantics are unchanged).
- **Skip path:** an additive-only model touches only the order-1 query — asserted on
  the branch, not by diffing output.

---

## Open questions

### D1. Is cross-specific dominance needed on day one?

The one capability v3 defers that a crossbreeding program might actually want. My
read: no, because v2 defers the resolver anyway, so v2 could not evaluate it either.
But this is the maintainer's call and it is the only cut worth reconsidering.

| Option | Notes |
|---|---|
| **Defer to its own plan, with the resolver** ← recommended | Storage and matching rules designed together. Retrofit is one child table, zero row changes |
| Add `genome_effect_member_origins` now | +1 table (3 total), still 13 fewer than v2. Storage with no interpreter — only worth it if you have a concrete near-term cross to parameterize |

### D2. Keep `'additive'` / `'dominance'` as named contrasts, or indicators only?

| Option | Notes |
|---|---|
| **Keep both, 4-value closed set** ← recommended | `'additive'`/`'dominance'` are 1 row instead of 3 and keep variance rescaling analytic. Cost is two branches in the evaluator |
| `'allele'` + `'indicator'` only | Absolutely minimal and fully general, but a 1000-QTL dominance model becomes 3000 rows and `rescale_effects_to_target()` must derive variance numerically |

### D3. Existing `.duckdb` files — migrate or regenerate?

Unchanged from v2's D4, and easier here: the v3 shape is close enough to today's
that migration is a single `INSERT … SELECT` per table. **Regenerate** is still
recommended per "Break Freely"; migration is now cheap enough to be a fallback rather
than a project.

### D4. Default contrast for `define_additive_effects()`

Unchanged from v2's D8. Today's engine is frequency-dependent (`allele − p`), so
`'allele'` with a stored `base_allele_freq` *is* the current default and no decision
is forced by this plan. The functional-vs-Cockerham choice becomes live only when a
genotype-level writer lands.

---

## Explicitly out of scope

Unchanged from v2, plus the deferred items above:

- Polyploid and multiallelic **simulation**.
- Lifting the `from_parent_1 + from_parent_2 <= 2` CHECK at `define_genome.R:368`.
- The scope resolver **and** its tables (v3 defers both together).
- Haplotype-level (phased-combination) effects.
- Mutation and any change to `make_gametes.cpp`.
- Estimating non-additive effects (`add_ebv()` / BLUPF90 remain additive).
- Non-additive evaluation, `ind_genetic_value`, variance-component renaming, and
  `phenotype_components.genome_effect_types` activation — each its own plan, as v2
  concluded.
