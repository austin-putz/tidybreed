# Genome Effects — A General Architecture for Additive, Dominance, and Epistasis

**Status:** design proposal, nothing implemented.
**Revision:** v2 (2026-09-04), after `plans/update_genome_effects_codex_review_v1.md`.
**Scope:** Part 1 (storage) is meant to be **locked** before any code. Part 2
(evaluation, outputs, variance) is a **sketch** that exists only to prove Part 1
is sufficient; it is explicitly not locked here.

## Revision history

**v1 → v2.** The Codex review was accepted on almost every substantive point.
The v1 "locked design" generalized the *number of loci* in a term but still
hard-coded the *state space at each locus* to one biallelic dosage, treated
free-text contrast names as if they defined contrast functions, and left
multi-line applicability on a single nullable `line_name`. Those three would each
have forced exactly the rewrite this project exists to prevent. v2 replaces them
with normalized genotype states, explicit contrast definitions, and origin-aware
scopes. A full accept/reject/decide disposition is in **§ Codex review v1 —
disposition** near the end; two items were rejected and four became decisions
for the maintainer.

---

## Why

`genome_effects` was built with a `genome_effect_type` dimension so dominance
could be added as **rows** rather than a schema change. That works for
dominance. It does **not** work for epistasis: an epistatic coefficient belongs
to a *set* of loci, and the current table has exactly one `locus_name` column,
so there is no way to say "this coefficient belongs to loci *j* **and** *k*".

Adding epistasis later on top of a one-locus-per-row table is precisely the
"another fundamental rewrite" that the Schema Design Bias section of `CLAUDE.md`
forbids. This plan reshapes effect storage **once**, pre-1.0.

Goals, in priority order:

1. Any number of loci per term, any mix of gene action, **now**.
2. A user can type in an arbitrary two-locus genotype table on day one, with no
   sampler in the package.
3. An effect's numeric meaning is recoverable **entirely from stored data** —
   `restore_pop()` must not depend on code remembering what a label meant.
4. Multi-line/crossbreeding applicability richer than one line name: A/B origin
   pairs, reciprocals, per-locus composition.
5. Ploidy-, multiallelic-, and hemizygous-ready state vocabulary, without
   implementing that biology now.
6. The order-1 additive engine keeps working, numerically unchanged.

---

## Where we stand today (audit, 2026-09-04)

### The table

`R/open_pop.R:285-295` — the only DDL:

```sql
CREATE TABLE genome_effects (
  id_genome_effect   INTEGER PRIMARY KEY,
  locus_name         VARCHAR NOT NULL,
  line_name          VARCHAR,
  trait_name         VARCHAR NOT NULL,
  genome_effect_type VARCHAR NOT NULL,
  genome_value       DOUBLE  NOT NULL,
  base_allele_freq   DOUBLE
)
```

No CHECK, no enum, no FK, no unique key. It is created in `open_pop.R` while
`genome_meta` is created later in `define_genome.R:272` — so a SQL foreign key
to a locus is not currently possible even in principle. Uniqueness of
`(trait_name, genome_effect_type, line_name)` is enforced only by the
DELETE-then-INSERT in `define_additive_effects()`.

### Every site that touches the type dimension

| Role | Sites |
|---|---|
| Writers | `define_additive_effects.R:307,318` (single-trait), `:502,513` (multi-trait), `:411` (`method="union"` membership read) |
| Readers | `add_tbv.R:179,215,221`; `add_phenotype.R:253`; `tidybreed_pop.R:159` (print `n_qtl`) |
| Registries | `sql_utils.R:85` (`TABLE_RESERVED_COLS`), `:135` (`TABLE_PRIMARY_KEYS`), `:161` (`TABLE_ROW_KEYS`), `:252` (`SYSTEM_TABLES`), `schema.R:100-117,671` |

`src/` has **zero** references. Tests reference it only at
`tests/testthat/test-add_tbv.R:25` (the independent R-side TBV oracle).

### Replacement scope is accidental, not designed

`define_additive_effects()` deletes **every** additive row for a trait/line
before inserting (`:297-322`, `:497-517`). So redefining a trait's architecture
silently replaces the whole locus set — a *model replacement*, expressed as a
DELETE predicate. v2 makes that an explicit, versioned model object.

### Dominance status

`"dominance"` / `"epistasis"` appear in exactly one functional line:

```r
# R/define_effect_cov_matrix.R:119
genetic_effects <- c("gen_add", "dominance", "epistasis")
```

Accurate statement (per the review's §1.5 correction, adopted):

> The covariance table can already store a matrix labeled `dominance`, but no
> functional dominance pipeline exists.

No effect writer, effect reader, simulation path, phenotype path, or result
table uses that matrix. This argues **against** coupling the new effect schema
to the current variance-component naming scheme — see Part 2.

Also stubbed and currently **dead**:
`phenotype_components.genome_effect_types VARCHAR DEFAULT 'additive'`
(`open_pop.R:334`, written at `define_phenotype.R:466,493`, read by nothing).

### The two structural blockers

**(a) The TBV query is per-allele-copy, and that is load-bearing.**
`add_tbv.R:206-227` sums `(h.allele − COALESCE(e.base_allele_freq, 0)) *
e.genome_value` over `ind_haplotype` rows, joining `genome_effects` on
`locus_name` with a per-locus `line_origin` precedence fallback. Consequences:

- Effects are **per allele copy**, not per collapsed genotype.
- A crossbred individual can use a **different coefficient and centering
  frequency for each inherited copy** at the same locus.
- `line_name` therefore does not mean "the individual belongs to this line". It
  means "this effect applies to an allele copy with this founding line origin".
  **The new schema must preserve that distinction explicitly** — it is why
  `evaluation_unit` and `match_unit = 'allele_copy'` exist below.
- `NULL base_allele_freq` is silently treated as zero. That is a bad invariant
  for any statistical coding; v2 makes frequency `NOT NULL` in the reference.

**(b) Collapsing to dosage destroys `line_origin`.** Dominance and epistasis
need the genotype; the `GROUP BY` that produces dosage erases which copy came
from which line. See the multi-line section.

### Two smaller findings

- `rescale_effects_to_target()` (`define_additive_effects.R:534-542`) computes
  `V_A = Σ 2p(1−p)a²` — the statistical-coding additive variance, which
  coincides with the biological one only because `d = 0`.
- `extract_genotypes()` validates `effects_tbl$table_name == "genome_effects"`
  (`R/extract_genotypes.R:124-127`). That **nominal** type check makes the schema
  hard to evolve; it should accept any filtered table carrying `locus_id`.

---

## The unifying abstraction

Every term in a general genetic model has one shape:

> **coefficient × ∏ over a set of loci of a per-locus contrast function**

| Term | Order | Contrast per locus |
|---|---|---|
| intercept | 0 | (empty product = 1) |
| additive at *j* | 1 | A at *j* |
| dominance at *j* | 1 | D at *j* |
| A×A(*j,k*) | 2 | A at *j*, A at *k* |
| A×D(*j,k*) | 2 | A at *j*, D at *k* |
| D×A(*j,k*) | 2 | D at *j*, A at *k* |
| A×A×A | 3 | A at each |

Two refinements v1 got wrong and v2 fixes:

**A contrast is data, not a name.** `x_A` and `x_D` are *functions from a local
genotype state to a number*. Storing the string `'additive'` records what a
contrast is called, not what it does. That breaks for multiallelic loci
(several additive and several dominance degrees of freedom), polyploid loci
(*k−1* non-additive contrasts), non-HWE bases, and any user-defined basis — and
it means a restored database is not self-describing. **Contrasts get their own
table with explicit per-state values.**

**Order alone does not identify a term.** A×D and D×A are distinct only
*relative to a canonical ordering of the loci*. Multiplication is commutative,
so user input order cannot define biology. **Members are canonicalized by
ascending `locus_id`.**

---

# Part 1 — Effect definition storage (lock this)

## Design invariants

These are the acceptance criteria. The schema below is one way to satisfy them;
the invariants themselves should not be traded away.

1. Any number of distinct loci per coefficient term; zero members = intercept.
2. Member order is canonical (`locus_id` ascending); the same product cannot be
   stored two ways.
3. A contrast's numeric meaning is recoverable entirely from stored rows.
4. Multiple contrasts of the same broad class may exist at one locus.
5. A local genotype state records copy number, an arbitrary allele-count vector,
   and optionally the allele-to-line/parent-origin binding.
6. No parsed compound keys anywhere. No serialized dosage strings.
7. Applicability supports common, per-copy line origin, cross-specific
   (A/B), reciprocal, and per-locus composition — without schema change.
8. Statistical effects reference an **immutable** base snapshot.
9. Coefficient payloads and surface payloads are never confused through a
   nullable column.
10. Hemizygous and polyploid local states are representable.
11. Missing surface cells have explicit, declared behavior.
12. Causal loci are discoverable through a view regardless of order or
    representation.
13. All writes and model replacements are single atomic transactions.

## Schema

Sixteen tables replace one. That cost is named honestly in **Decision D2**;
it buys invariants 1–13, and every one of them corresponds to a rewrite avoided.

All surrogate keys are assigned via `next_int_id()`, per package convention.
Every nullable-keyed table uses an **R-enforced NULL-normalized logical key**,
never a SQL `PRIMARY KEY` (DuckDB requires PK columns `NOT NULL`) — same
precedent as `genome_map`, `chr_inheritance`, and `founder_haplotypes`.

All the small VARCHAR discriminators below (`state_domain`, `evaluation_unit`,
`match_unit`, `line_match_type`, `missing_state_action`, `value_semantics`,
`source_type`, `contrast_class`) are **closed sets validated in R**, not open
labels. This is the distinction that makes them safe where a free-text
`locus_contrast` was not: they name *structure that code must branch on*, and
their valid set is part of the package, whereas `contrast_name` / `coding_name`
are open **provenance** whose semantics live in `genome_contrast_values`.

### 1. `genome_effect_models` — the model header

```sql
CREATE TABLE genome_effect_models (
  id_genome_effect_model INTEGER PRIMARY KEY,
  model_name             VARCHAR NOT NULL,
  trait_name             VARCHAR NOT NULL,
  model_version          INTEGER NOT NULL DEFAULT 1,
  notes                  VARCHAR
);
```

Logical key `(trait_name, model_name, model_version)`. Groups every term and
surface that forms one coherent architecture, turning today's accidental
"DELETE everything for this trait/line" into an explicit, atomic model
replacement. `model_name = 'default'` covers the whole existing API. No
`is_active` flag unless runtime selection genuinely needs persistent state.

### 2. Immutable reference populations

```sql
CREATE TABLE genome_effect_references (
  id_genome_effect_reference INTEGER PRIMARY KEY,
  reference_name             VARCHAR NOT NULL,
  source_type                VARCHAR NOT NULL,   -- 'founder_haplotypes'|'current_pop'|'legacy'|'user'
  source_description         VARCHAR,
  created_generation         INTEGER,
  created_replicate          INTEGER,
  notes                      VARCHAR
);

CREATE TABLE genome_effect_reference_alleles (
  id_genome_effect_reference INTEGER NOT NULL,
  locus_id                   INTEGER NOT NULL,
  allele_code                VARCHAR NOT NULL,
  line_name                  VARCHAR,
  allele_freq                DOUBLE  NOT NULL
);
```

Logical key `(id_genome_effect_reference, locus_id, allele_code, line_name)`,
NULL-normalized. Validate `0 <= allele_freq <= 1` and per-locus/context sums.

**Immutable once referenced.** A new snapshot is a new integer ID even if the
user-facing name repeats. This is what `base_name = 'current_pop'` could never
be: that string does not identify *which* current population, at what
generation, replicate, or filter.

`allele_code` makes the table multiallelic-ready; today every locus has codes
`'0'` and `'1'`, and `allele_freq` for `'1'` is exactly today's
`base_allele_freq`.

**Reserved, not built now:**

```sql
CREATE TABLE genome_effect_reference_states (
  id_genome_effect_reference INTEGER NOT NULL,
  id_genome_state            INTEGER NOT NULL,
  line_name                  VARCHAR,
  state_freq                 DOUBLE  NOT NULL
);
```

Needed only for bases that must represent departure from Hardy–Weinberg
proportions. Allele frequencies alone do not determine the genotype-state
distribution off HWE.

### 3. Normalized local genotype states

```sql
CREATE TABLE genome_states (
  id_genome_state INTEGER PRIMARY KEY,
  locus_id        INTEGER,            -- NULL = locus-independent template
  copy_number     INTEGER NOT NULL,
  state_domain    VARCHAR NOT NULL    -- 'allele'|'allele_line'|'allele_line_parent'
);

CREATE TABLE genome_state_alleles (
  id_genome_state_allele INTEGER PRIMARY KEY,
  id_genome_state        INTEGER NOT NULL,
  allele_code            VARCHAR NOT NULL,
  line_name              VARCHAR,
  parent_origin          UTINYINT,
  copy_count             INTEGER NOT NULL
);
```

Logical key for count rows: `(id_genome_state, allele_code, line_name,
parent_origin)`, NULL-normalized. The surrogate row ID lets one allele appear
in copies from multiple lines within one state.

`state_domain` declares whether the nullable origin columns participate in state
identity. **NULL is never an undocumented wildcard.**

**Amendment to the review (accepted with change):** `genome_states.locus_id` is
**nullable**. The review made it `NOT NULL`, which forces three state rows per
QTL even though the biallelic diploid states `{0:2}`, `{0:1,1:1}`, `{1:2}` are
structurally identical at every locus. Making it nullable means those three
rows are written **once for the whole genome**. A locus-specific state is needed
only when allele codes are locus-specific (multiallelic). Concretely, for a
1000-QTL statistical-coding model: 3 shared state rows instead of 3000, with no
loss of generality.

Validation: `copy_number >= 0`; `copy_count > 0`; `sum(copy_count) =
copy_number`; no duplicate normalized state at one locus (or among templates);
allele codes valid for the locus once `genome_meta` gains an allele dimension.

Examples: origin-agnostic biallelic diploid → `{0:2}`, `{0:1,1:1}`, `{1:2}`.
Hemizygous → `copy_number = 1`. Autotetraploid AABC → `A:2, B:1, C:1`.
Origin-aware → "allele 1 from line A, allele 0 from line B" is distinguishable
from the reverse. **No serialized dosage key anywhere.**

### 4. Explicit contrast definitions

```sql
CREATE TABLE genome_contrasts (
  id_genome_contrast         INTEGER PRIMARY KEY,
  contrast_name              VARCHAR NOT NULL,
  contrast_class             VARCHAR NOT NULL,   -- 'additive'|'dominance'|'non_additive'
  coding_name                VARCHAR NOT NULL,   -- provenance, e.g. 'functional_diploid_v1'
  coding_version             INTEGER NOT NULL DEFAULT 1,
  evaluation_unit            VARCHAR NOT NULL,   -- 'allele_copy'|'locus_genotype'
  locus_id                   INTEGER,            -- NULL = locus-independent
  id_genome_effect_reference INTEGER,            -- required for frequency-dependent codings
  notes                      VARCHAR
);

CREATE TABLE genome_contrast_values (
  id_genome_contrast INTEGER NOT NULL,
  id_genome_state    INTEGER NOT NULL,
  contrast_value     DOUBLE  NOT NULL,
  PRIMARY KEY (id_genome_contrast, id_genome_state)
);
```

Logical key on the header: `(contrast_name, coding_name, coding_version,
locus_id, id_genome_effect_reference)`, NULL-normalized. **The reference must be
part of contrast identity**: two traits can legitimately want a `cockerham_hwe_v1`
additive contrast at the same locus against *different* base populations (one
defined at the founders, one at generation 5), and those carry different `p` and
therefore different values. Omitting the reference from the key collides them.

**`genome_contrast_values` is authoritative.** `contrast_class` and
`coding_name` aid interpretation and filtering; they never define the math.
Future code evaluates a stored contrast without hard-coding its label, and a
restored database is fully self-describing.

**`evaluation_unit` preserves today's semantics.** `'allele_copy'` evaluates the
contrast once per `ind_haplotype` row against a `copy_number = 1` state — which
is exactly how additive TBV works today, and the only way per-copy `line_origin`
resolution survives. `'locus_genotype'` evaluates once on the collapsed local
genotype, which is what dominance requires. A genotype-level additive contrast
used inside an A×D term is simply a different contrast ID of the same broad
class.

**Values are always materialized, including for built-in codings** (Decision
D5). A built-in coding is a *generator* the writer runs once against the
reference; evaluation then reads rows, uniformly, for built-in and custom
contrasts alike. This is what makes invariant 3 hold.

Amendment to the review: the review's `locus_id NOT NULL` would force a separate
contrast row per locus even for a functional coding, where `x_A(g) = g − 1` is
the *same function everywhere*. Nullable `locus_id` collapses the whole
functional diploid additive basis to **1 contrast row + 3 value rows for the
entire genome**. Frequency-dependent codings still get one contrast per locus
(because `p_j` differs) but reuse the same three shared states.

### 5. General applicability scopes

```sql
CREATE TABLE genome_effect_scopes (
  id_genome_effect_scope INTEGER PRIMARY KEY,
  scope_name             VARCHAR NOT NULL,
  match_unit             VARCHAR NOT NULL,  -- 'common'|'allele_copy'|'locus_genotype'|'individual'
  notes                  VARCHAR
);

CREATE TABLE genome_effect_scope_origins (
  id_genome_effect_scope_origin INTEGER PRIMARY KEY,
  id_genome_effect_scope        INTEGER NOT NULL,
  member_slot                   INTEGER,           -- NULL = applies to all members
  origin_slot                   INTEGER NOT NULL,
  line_match_type               VARCHAR NOT NULL,  -- 'exact'|'any'|'unknown'
  line_name                     VARCHAR,
  parent_origin                 UTINYINT,
  copy_count                    INTEGER NOT NULL
);
```

Logical key `(id_genome_effect_scope, member_slot, origin_slot)`,
NULL-normalized.

`line_match_type` makes `exact`, `any`, and `unknown` **different stored
conditions**: `exact` requires `line_name`; `any` and `unknown` require it to be
NULL. Wildcard semantics are never inferred from a SQL NULL comparison.

Every term and surface carries a **non-null** scope ID; a distinguished common
scope (`match_unit = 'common'`, zero origin rows) replaces the old nullable
`line_name`. `match_unit = 'allele_copy'` reproduces today's per-copy additive
behavior exactly. An A/B dominance scope carries two origin rows at member slot
1. Reciprocals use `parent_origin`. Tetraploid compositions use `copy_count`
rather than four repeated rows.

**Amendment to the review — the origin double-specification rule.** Line origin
can now be expressed in *two* places: inside a `genome_states` row
(`state_domain = 'allele_line'`) and inside a scope origin row. Leaving that
unresolved would reintroduce exactly the ambiguity this schema exists to remove.
The rule:

> **Scope** answers *when does this effect apply*. **State** answers *what
> genotype is this*. Origin belongs in the state only when the contrast value
> depends on the **binding** of allele to origin (allele 1 from A vs. allele 1
> from B). Everything else is scope.

The validator **rejects** any term whose member simultaneously uses an
origin-aware state and an origin-carrying scope row at the same slot.

**Not yet designed:** the matching and precedence *resolver* — which scope wins
when several match. The review flags this and is right. The tables above fix the
**cardinality** (zero, one, or many origin conditions per member), which is the
part that cannot be added later. See Decision D3.

### 6. Coefficient terms

```sql
CREATE TABLE genome_effects (
  id_genome_effect       INTEGER PRIMARY KEY,
  id_genome_effect_model INTEGER NOT NULL,
  id_genome_effect_scope INTEGER NOT NULL,
  effect_class           VARCHAR NOT NULL,   -- materialized, validated, never identity
  genome_value           DOUBLE  NOT NULL
);

CREATE TABLE genome_effect_members (
  id_genome_effect   INTEGER NOT NULL,
  member_slot        INTEGER NOT NULL,
  locus_id           INTEGER NOT NULL,
  id_genome_contrast INTEGER NOT NULL,
  PRIMARY KEY (id_genome_effect, member_slot)
);
```

```text
term contribution = genome_value × ∏ (member contrast values)
```

Invariants:

- A coefficient term has ≥ 1 member; `effect_class = 'intercept'` has exactly
  zero (its empty product is 1). This is the zero-order term v1 omitted.
- A locus occurs at most once per term.
- **Members are canonicalized by ascending `locus_id`**; `member_slot` is
  assigned after canonicalization, never from input order.
- A contrast's `locus_id`, when non-NULL, equals the member's `locus_id`.
- **Term identity is `(model, scope, ordered contrast IDs)` — not the class
  label.** Reversing input locus order produces the *same* stored term;
  swapping which locus carries dominance produces a *different* one.
- `effect_class` is materialized for filtering/display and validated against the
  ordered broad contrast classes. It is never the term's identity and never the
  source of its formula.
- **No stored `effect_order`.** Member count is authoritative; a view exposes it.
- `INTEGER`, not `UTINYINT`, for `member_slot` — no artificial 255 ceiling on a
  small configuration table.

### 7. Direct genotype-value surfaces

```sql
CREATE TABLE genome_effect_surfaces (
  id_genome_effect_surface INTEGER PRIMARY KEY,
  id_genome_effect_model   INTEGER NOT NULL,
  id_genome_effect_scope   INTEGER NOT NULL,
  surface_name             VARCHAR NOT NULL,
  value_semantics          VARCHAR NOT NULL,  -- 'component' (reserved: 'total_genotypic')
  missing_state_action     VARCHAR NOT NULL,  -- 'error'|'zero'|'default'
  default_value            DOUBLE,
  notes                    VARCHAR
);

CREATE TABLE genome_effect_surface_members (
  id_genome_effect_surface INTEGER NOT NULL,
  member_slot              INTEGER NOT NULL,
  locus_id                 INTEGER NOT NULL,
  PRIMARY KEY (id_genome_effect_surface, member_slot)
);

CREATE TABLE genome_effect_surface_cells (
  id_genome_effect_cell    INTEGER PRIMARY KEY,
  id_genome_effect_surface INTEGER NOT NULL,
  genome_value             DOUBLE  NOT NULL
);

CREATE TABLE genome_effect_surface_cell_states (
  id_genome_effect_cell INTEGER NOT NULL,
  member_slot           INTEGER NOT NULL,
  id_genome_state       INTEGER NOT NULL,
  PRIMARY KEY (id_genome_effect_cell, member_slot)
);
```

Exactly one matching cell per surface contributes its `genome_value` to the
model sum. Separate tables mean `genome_value` is `NOT NULL` everywhere — v1's
tagged union (one nullable payload column meaning two different things) is gone.

`value_semantics` resolves the ambiguity v1 never stated: whether a hand-entered
3×3 table holds a **component to add** (`'component'`, the only value supported
now) or the **whole genotypic value** (`'total_genotypic'`, reserved). Those
produce different totals, so it must be declared, not assumed.

Validation: members canonicalized by `locus_id`; every cell has exactly one
state per surface member; each state belongs to the member's locus (or is a
template); cell state-tuples unique; `'default'` requires a non-NULL
`default_value`; completeness or deliberate sparsity agrees with
`missing_state_action`. **Missing cells are never silently interpreted.**

Same four tables hold 3×3 diploid, 5×5 tetraploid, mixed-copy-number, and
multiallelic tables.

### 8. Views (ergonomics — the base tables must not force six joins on users)

- `genome_effect_terms` — one row per term: trait, model, order, effect class,
  coefficient, human-readable member description.
- `genome_effect_loci` — one row per (model, effect, locus); the QTL-set view.
- `genome_effect_surfaces_long` — surface/cell/member states in tidy long form.

`extract_genotypes()` must accept **any filtered table carrying `locus_id`**,
replacing the nominal `table_name == "genome_effects"` check at
`R/extract_genotypes.R:124-127`.

QTL terminology becomes precise, which changes `tidybreed_pop.R:159`:

- a locus is **causal** for a model if it appears in any term or surface;
- an **additive QTL** appears in an order-1 additive-class term;
- an **epistatic-only** locus can be causal with no order-1 row at all.

## Worked examples (these are the acceptance fixtures)

**Existing additive effect.** Model `(ADG, 'default', 1)`. One contrast per QTL,
`evaluation_unit = 'allele_copy'`, `locus_id = j`, referencing the shared
`copy_number = 1` states `{0:1}` and `{1:1}` with values `−p_j` and `1 − p_j`.
Scope `match_unit = 'allele_copy'`, one `exact` origin row for line A. Term:
coefficient `a_j`, one member. Summing over both copies yields
`(dosage − 2p_j)·a_j` — **numerically identical to today's engine**. The common
fallback is a separate term with the common scope.

**Diploid dominance.** Contrast over states `{0:2}`, `{0:1,1:1}`, `{1:2}` with
values `0, 1, 0` for a functional heterozygote indicator. Common dominance uses
the common scope; A/B-specific dominance uses a `locus_genotype` scope with one
A-origin and one B-origin copy at member slot 1. **No meaning is inferred from
the word "dominance"** — the three stored values define it.

**Pairwise A×D.** Members canonicalized by `locus_id`; lower-ID locus takes an
additive contrast, higher-ID a dominance contrast; `effect_class = 'add_dom'` is
a validated display value. Reversing input order stores the same term; swapping
which locus is dominant stores the distinct `dom_add` term.

**Autotetraploid locus.** Define the five states and as many contrast IDs as the
chosen basis requires, storing every contrast's value at every state. The effect
tables need no new columns and no advance agreement on labels like `dom_2`,
`dom_3`, `dom_4` — because polyploid contrast systems are basis-dependent. **A
free VARCHAR reserves labels; explicit state values reserve semantics.**

**Arbitrary two-locus table.** One surface, two canonical member loci, one cell
per pair of state IDs, `missing_state_action = 'error'` for a complete manual
table. This is the day-one "type in your own gene action" feature.

**Reciprocal multi-line effect.** Two scopes — (parent 1/line A, parent 2/line
B) and (parent 1/line B, parent 2/line A) — with different terms attached. No
`line_name_1`, `line_name_2`, or future reciprocal-cross column.

## Validation and transaction rules

Constraints span parent/child tables, so they are R-enforced. **Every public
write is one DuckDB transaction running a single comprehensive validator before
commit**, mirroring `define_chromosome()`'s delete-then-insert-then-validate-then-commit
pattern.

| Group | Checks |
|---|---|
| Model / reference | logical keys unique; referenced traits/loci/scopes/references exist; references used by a contrast are immutable; numeric values finite |
| States | allele counts sum to `copy_number`; no duplicate normalized state; origin columns agree with `state_domain`; template vs. locus-specific used consistently |
| Contrasts | contrast states belong to the contrast's locus (or are templates); state coverage complete for built-in codings; frequency-dependent codings have a valid immutable reference; reference frequencies sum correctly per locus/context |
| Terms | non-intercept has members, intercept has none; member loci unique and canonical; contrast locus agrees with member locus; no duplicate term identity within model/scope; `effect_class` agrees with ordered contrast classes |
| Surfaces | members unique and canonical; cells carry exactly the surface's member slots; cell states agree with member loci; no duplicate state tuple; completeness agrees with `missing_state_action` |
| Scopes | common scopes have no origin rows; `copy_count > 0`; slot-specific rows reference existing member slots; valid `parent_origin`; normalized origin tuples unique; `line_match_type` and `line_name` agree; **no member uses an origin-aware state and an origin-carrying scope simultaneously** |

## Ploidy analysis — what was explicitly considered

**Current state.** `ploidy` must be `2` everywhere: `add_founders.R:113-118`
rejects anything else, `add_offspring.R:552-554` computes and asserts it, and
`assert_ploidy_2()` (`R/ploidy_helpers.R`) guards `add_dosage.R:83` and
`extract_genotypes.R:187`. The real blocker to polyploidy is elsewhere — the
`from_parent_1 + from_parent_2 <= 2` CHECK at `define_genome.R:368` — and is a
separate project. `ind_haplotype` is **already ploidy-general**: `strand` is
"copy within a parent's contribution; always 1 for diploids", so an
autotetraploid is `parent_origin ∈ {1,2} × strand ∈ {1,2}` = 4 rows per locus,
no schema change.

**What is ploidy-general in this design:**

- **`genome_states`** carries an explicit `copy_number` and an arbitrary
  allele-count vector. Hemizygous (`copy_number = 1`), diploid, and tetraploid
  states are all just rows. **This is the single biggest improvement over v1**,
  whose `dosage_key` string could not distinguish a haploid alternate genotype
  from a diploid heterozygote from one alternate copy in a tetraploid, and could
  not represent a tetraploid A=2/B=1/C=1 at all.
- **`genome_contrasts`.** A tetraploid locus gets as many contrast IDs as its
  basis requires. v1 tried to solve this with reserved labels `dom_2`/`dom_3`/
  `dom_4`; that reserves names, not math, and polyploid bases are
  basis-dependent.
- **The additive path.** `add_tbv()` never calls `assert_ploidy_2()` and does
  not need to: summing an `allele_copy` contrast over however many haplotype
  rows exist gives the right answer at any ploidy. This is a genuine strength of
  the existing per-copy formulation and the new schema preserves it via
  `evaluation_unit = 'allele_copy'`.
- **Scopes** use `copy_count`, so a tetraploid A/A/B/B origin composition is two
  rows, not four.
- **Surfaces.** A 5×5 tetraploid table is the same four tables as a 3×3 diploid
  one.

**What genuinely breaks at ploidy > 2:** dominance is not one thing. A locus
with *k* copies has *k+1* genotypes, hence 1 additive and **`k−1`** non-additive
degrees of freedom. Diploids get away with a single `x_D` because `k−1 = 1`.
Absorbed entirely by explicit contrasts, with no new columns.

**Considered and out of scope:** `strand` never enters a contrast under the
standard dosage/state-based decomposition (true haplotype effects — the value of
a specific phased combination — would be a different member semantics, reserved
not designed); mixed-ploidy crosses (`add_offspring.R:548-554` already computes
offspring ploidy as the sum of gamete contributions, so a triploid is
architecturally anticipated — storage fine, sampler absent); polysomic vs.
disomic inheritance and double reduction are *transmission*
(`make_gametes.cpp`, `chr_inheritance`), untouched here.

## Multi-line / crossbreeding analysis — what was explicitly considered

**What today's `line_name` actually means.** Not "the individual belongs to this
line" but "this effect applies to an allele copy with this founding line
origin". `add_tbv.R:215-226` resolves it **per allele copy**, per locus, with a
`NOT EXISTS` fallback to the common row. A crossbred individual legitimately
uses a different coefficient *and* a different centering frequency for each of
its two copies at one locus. The new schema preserves this exactly:
`evaluation_unit = 'allele_copy'` + `match_unit = 'allele_copy'`.

**The genuinely new problem.** Dominance and epistasis need the collapsed
genotype, and collapsing discards the origin composition. A diploid dosage of 1
could be A/A-origin, A/B-origin, B/B-origin, or include an unknown origin. There
is no principled way to pick one scalar line-specific dominance coefficient
after that information is gone.

**v1's answer was insufficient and partly buggy.** A single nullable `line_name`
can express common, line A, and line B — but not A/B-pair-specific dominance,
not reciprocal A/B vs. B/A, not "A/A at locus 1 interacting with A/B at locus
2", not a tetraploid A/A/B/B composition, and not one line condition per member
of a higher-order term. And v1's SQL sketch was wrong: `COUNT(DISTINCT
line_origin)` **ignores NULL**, so an origin set `{A, NULL}` reports one
distinct line with `MIN = A` and silently masquerades as a within-A genotype.

**v2's answer.** Storage must distinguish five concepts, and now does:

| Concept | Where it lives |
|---|---|
| 1. Allele-copy origin (today's additive `line_name`) | scope `match_unit = 'allele_copy'` + `exact` origin row |
| 2. Within-locus origin composition (A/A, A/B, B/B; ploidy multisets) | scope `match_unit = 'locus_genotype'`, several origin rows at one slot |
| 3. Per-member composition (different at each locus of an epistatic term) | `member_slot` on the origin row |
| 4. Parent-origin role (reciprocals, imprinting) | `parent_origin` on the origin row |
| 5. Common fallback | a deliberate `match_unit = 'common'` scope, **not** SQL NULL logic |

**On heterosis.** v1 claimed the "both copies same line, otherwise common" rule
made heterosis emergent. That overclaimed. The rule permits heterosis arising
from directional dominance plus divergent line frequencies, but with only one
`line_name` it **cannot represent an explicit A/B-specific dominance
deviation** — which is the thing a crossbreeding program most wants to
parameterize. v2's scope tables can. The corrected claim: *heterosis from
directional dominance is emergent; cross-specific dominance is now also
directly expressible.*

**Before implementing the resolver**, write a truth table covering common,
A-specific, B-specific, A/B-specific, reciprocal A/B, unknown origin, haploid,
diploid, and tetraploid cases. Accept the scope tables only if every row has an
unambiguous stored representation. (Decision D3.)

**Considered and out of scope:** the precedence resolver itself; line ×
environment beyond what scopes give; joint multilocus reference frequencies
(needed for exact epistatic orthogonality under LD — do not claim orthogonality
from marginal allele frequencies alone).

---

# Part 2 — Downstream (sketch only; not locked here)

The review is right that v1 mixed four projects and that mixing obscured what
must be locked. Part 2 is retained rather than deleted, for one reason: a
storage design cannot be accepted without evidence it can be evaluated. Nothing
below is a commitment, and each item gets its own plan.

**Evaluation.** Three query shapes, dispatched on `evaluation_unit` and member
count. (i) `allele_copy` order-1 — today's query, unchanged in shape, joining
contrast values instead of a scalar coefficient. (ii) `locus_genotype` order-1 —
derive the local state inline from `ind_haplotype` (never read `ind_genotype`,
which is an opt-in cache that may be empty), then join
`genome_contrast_values`. (iii) order ≥ 2 — the same state derivation, joined
once per member slot, with the contrast product formed per term; surfaces join
`genome_effect_surface_cell_states` on the full state tuple. A model with no
`locus_genotype` contrast and no order ≥ 2 term takes path (i) only, so an
additive-only simulation never builds any of this.

**Breeding value under epistasis.** v1 asserted `add_tbv()` should switch to
`α = a + d(q − p)`. That one-locus average substitution effect is correct for a
functional one-locus diploid model, but **once functional epistasis exists,
average effects also depend on interaction coefficients, frequencies at other
loci, and LD.** So the one-locus α does not generally convert a functional
epistatic genotypic model into a true breeding value. Storage implication, and
the reason this matters here: **do not encode the promise that an
`additive`-classified coefficient is the breeding value.** `effect_class` is a
filter, not a guarantee.

**Individual-level outputs.** A long `ind_genetic_value(id_ind, trait_name,
component, genetic_value)` table, and whether it carries an additive row — its
own plan.

**Variance targeting.** With orthogonal contrasts **under linkage equilibrium**,
`V_term = coefficient² × ∏_slots Var(contrast_slot)`, with `Var(x_A) = k·p·q` at
ploidy `k` and `Var(x_D) = (2pq)²` for the diploid HWE dominance contrast. This
generalizes `rescale_effects_to_target()`. The LE assumption is load-bearing and
must not be dropped silently.

**Contrast codings.** The diploid HWE formulas are:

| Coding | `x_A` | `x_D` |
|---|---|---|
| functional | `g − 1` | `1` if `g == 1` else `0` |
| Cockerham (HWE) | `g − 2p` | `x_D(0) = −2p²`, `x_D(1) = 2pq`, `x_D(2) = −2q²` |

**Correction adopted from the review:** v1 called the second one
"NOIA-Cockerham". That is misleading. General NOIA (Álvarez-Castro & Carlborg
2007, and the multiallelic extension) constructs design matrices from **observed
genotype frequencies** and supports departure from Hardy–Weinberg proportions;
allele frequency `p` alone is insufficient. The formulas above are the HWE
Cockerham special case. Both become `coding_name` values — `'functional_diploid_v1'`,
`'cockerham_hwe_v1'`, later `'noia_v1'` — and, crucially, the *stored state
values* are what evaluation reads, so a future NOIA basis needs no new schema.

**Variance-component naming.** v1 proposed renaming `trait_var_comp.effect_name`
to `gen_dom` / `gen_add_add` and routing by a `startsWith("gen_")` prefix. The
prefix rule is **rejected** (see disposition): `effect_name` is an open
user-facing label, and a user naming a random effect `gen_flock` would silently
route it into `trait_var_comp`. A lexical prefix is not a relational type. The
routing should use an explicit destination/class argument or a normalized
variance-component type — in its own plan, after Part 1 is stable.

**`phenotype_components.genome_effect_types`** activation: its own plan.

---

## Implementation order

| Phase | Work | Gate |
|---|---|---|
| **A** | **Fixtures before DDL.** Hand-build table rows for: common additive; line-A additive + common fallback; line-A and line-B additive with different base frequencies; an A/B dominance scope; an A×D pair; a tetraploid contrast; a 3×3 surface; a reciprocal line effect | Every fixture representable with no serialized key and no special column. If one is not, the schema is wrong — stop |
| **B** | Create definition/reference tables + all six registries each. Populate built-in biallelic template states and the built-in contrast generators | Registry tests pass |
| **C** | Write the new architecture through `define_additive_effects()`: default model, immutable reference from the current base-frequency computation, allele-copy scopes, order-1 canonical terms | Coefficients preserved exactly |
| **D** | Switch readers via an order-1 flattened **development** view; verify the existing independent TBV oracle against it; then move QTL extraction and the print summary to the locus view | Oracle agrees |
| **E** | Delete the old table shape and the staging view | No `genome_effects` legacy columns remain |

**The staging view in D must not survive to a release.** `CLAUDE.md` forbids
compatibility shims; it is a development scaffold, deleted in E.

Phase A is the cheapest possible way to find out the schema is wrong, and it
costs nothing but a scratch DuckDB file. It is the single most valuable item the
review contributed.

Blast radius (`grep -rl genome_effects`): 12 R files, 8 test files (including
`helper-parity.R`), 9 man pages, 2 vignettes. Non-additive computation, output
tables, phenotype inclusion, and samplers are separate projects after this one.

### Writer API sketch

`define_additive_effects()` stays as the ergonomic front door for the 95% case.
A general writer is needed but its **input shape** should be designed after
Phase A fixtures exist, not before — see Decision D6.

---

## Testing strategy

- **Phase-A fixtures as tests.** Each of the eight fixtures becomes a test that
  writes the rows and reads them back through the views.
- **Step-C invariance — a development-time check, not a committed test.** Run
  the existing `test-add_tbv.R` oracle (`tests/testthat/test-add_tbv.R:16-40`)
  on the same seed before and after, and confirm the TBVs match. **Do not commit
  a golden-output fixture**: `CLAUDE.md` ("The only reproducibility contract is
  forward-looking") forbids tests comparing against pre-change output. What ships
  is the existing oracle, which recomputes from first principles.
- **Independent R oracle** for each new compute path.
- **Canonicalization.** Writing an A×D term with the loci supplied in either
  order must produce byte-identical rows; supplying dominance on the *other*
  locus must produce a different term.
- **Self-describing storage.** Close and `restore_pop()` a database, then
  evaluate every stored contrast **without** consulting any `coding_name`
  branch. This is the direct test of invariant 3.
- **Coding equivalence.** Functional `(a, d)` at base `p` converted to Cockerham
  `(α, δ)` with `α = a + d(q − p)`, `δ = d` must give identical genotypic values
  at all three states.
- **Hand-computed 2-locus surface**, checked cell by cell, plus its
  coefficient-form equivalent.
- **Origin double-specification** must be rejected by the validator.
- **NULL-origin regression.** An individual with `{A, NULL}` line origins must
  not be treated as within-A — the specific v1 bug.
- **Skip path.** An additive-only model must take path (i) only — asserted on
  the branch, not by diffing old output.
- **Registry tests.** Every new table needs **six** registrations, each already
  test-guarded. In `R/sql_utils.R`: `TABLE_RESERVED_COLS` (`:77`),
  `TABLE_PRIMARY_KEYS` (`:130`), `TABLE_ROW_KEYS` (`:157`), `SYSTEM_TABLES`
  (`:250`) — `test-schema-registries.R:117` and `:157` fail otherwise. In
  `R/schema.R`: `.schema_table_order()` and the matching
  `.<group>_descriptions()` helper, asserted by `test-schema-print.R`. Sixteen
  tables × six registrations is a real cost; a helper that registers a table
  once, in one place, is worth building first.

---

## Options considered

**Term header + variable-length member list** — kept from v1; the review
confirms it as the correct hyperedge abstraction and the correct rejection of
wide `locus_name_2`, `locus_name_3` columns.

**Wide slots** — rejected. Hard-caps at two loci; three-way epistasis becomes
the fundamental rewrite the schema bias forbids.

**v1's flat design (`locus_contrast` label + `dosage_key` + single
`line_name`)** — rejected on review. It generalized the number of loci but not
the state space, the contrast semantics, or the origin cardinality, and would
have forced a rewrite for polyploids, multiallelic loci, or cross-specific
effects.

**Two representations rather than one nullable union** — kept and strengthened.
Contrast terms and genotype surfaces answer different questions:

| Property | Contrast term | Genotype surface |
|---|---|---|
| Compact additive architecture | Excellent | Wasteful |
| Named A/D/epistatic decomposition | Explicit | Not intrinsic |
| Arbitrary gene action | Only with a full basis | Exact |
| User-defined polyploid basis | Via contrast-state rows | Directly |
| Variance-component classification | Available as metadata | Requires decomposition |
| Missing genotype states | Contrast validation | Surface missing-cell policy |

A shared model header lets both coexist — compact additive main effects plus one
hand-entered two-locus surface — without a nullable payload column.

---

## Codex review v1 — disposition

Full detail is appended to `plans/update_genome_effects_codex_review_v1.md`.

**Accepted (18).** Contrast definition + value tables (§3.1); normalized
genotype states replacing `dosage_key` (§3.2); origin-aware scope tables
replacing scalar `line_name`, including the `COUNT(DISTINCT line_origin)` NULL
bug (§3.3); coding as versioned provenance at contrast grain, and the
NOIA-vs-Cockerham naming correction (§3.4); canonical member ordering by
`locus_id`, `effect_class` as validated non-identity, no stored `effect_order`,
`INTEGER` over `UTINYINT` (§3.5); immutable integer reference IDs with
allele-keyed frequencies (§3.6); separate strongly-typed surface tables with
`NOT NULL` payloads and declared value semantics (§3.7); zero-member intercept
terms (§3.8); `locus_id` as the member key (§3.9); the incomplete
breeding-value-under-epistasis correction (§3.10); the model/version header
(§5.1); `evaluation_unit` (§5.4); `line_match_type` as explicit stored
conditions (§5.5); ergonomic views, the `extract_genotypes()` nominal-type-check
removal, and precise QTL terminology (§10); the comprehensive validator and
one-transaction rule (§11); fixtures-before-DDL (§12 Phase A); the `dominance
is not end-to-end` wording (§1.5); and the acceptance criteria (§13).

**Accepted with amendment (3).**

1. **`genome_states.locus_id` made nullable** (review had `NOT NULL`). Without
   this, the three structurally identical biallelic diploid states must be
   written once per QTL — 3000 rows for a 1000-QTL model instead of 3. Locus-specific
   states are needed only for locus-specific allele codes.
2. **`genome_contrasts.locus_id` made nullable** (review had `NOT NULL`). A
   functional coding is the *same function at every locus*; nullable `locus_id`
   collapses the whole functional diploid additive basis to 1 contrast row + 3
   value rows genome-wide. Frequency-dependent codings still get one row per
   locus.
3. **Origin double-specification rule added.** The review introduces origin in
   *two* places — `state_domain = 'allele_line'` and scope origin rows — and
   does not say which is canonical. Added an explicit rule (state = genotype
   identity when the value depends on the allele↔origin binding; scope =
   applicability) plus a validator check that rejects using both at one member.

**Rejected (2).**

1. **Removing Part 2 entirely from this plan** (§4). The review is right that
   v1 mixed four projects and that the mixing obscured what must be locked; the
   *status* separation is accepted and implemented as Part 1 / Part 2. But
   deleting the downstream sketch is rejected: a storage design cannot be
   accepted without evidence it is evaluable, and the maintainer explicitly
   asked for the whole picture, not the DDL alone. Part 2 is demoted to
   non-binding, not removed. The specific items the review wanted moved
   (`ind_genetic_value`, the evaluator, `phenotype_components` activation, the
   line-fallback calculation, samplers, variance-component renaming) are all
   now marked "own plan".
2. **`state_domain` / `match_unit` / `line_match_type` / `evaluation_unit` /
   `missing_state_action` as free-text VARCHARs.** The review rejects a
   free-text `locus_contrast` because "a label does not define state values"
   (§3.1), then introduces seven free-text discriminators of its own. The
   principle is right but was applied unevenly. Resolution: these are **closed
   sets validated in R**, because they name structure the package must branch
   on; `contrast_name` and `coding_name` stay open because they are provenance
   whose semantics live in `genome_contrast_values`. Stated explicitly in the
   schema preamble.

**Escalated to maintainer decision (4).** Build-now scope (D2), the scope
resolver (D3), existing-database migration (D4), and materialized vs. computed
contrast values (D5) — all in Open Questions below.

---

## Explicitly out of scope

- Polyploid and multiallelic **simulation** (states and contrasts are
  representable; transmission and sampling are not implemented).
- Lifting the `from_parent_1 + from_parent_2 <= 2` CHECK at
  `define_genome.R:368` — the actual polyploidy blocker, a separate project.
- The scope precedence **resolver** (tables yes, matching rules no — D3).
- Haplotype-level (phased-combination) effects.
- Mutation, and any change to `make_gametes.cpp`.
- Estimating non-additive effects (`add_ebv()` / BLUPF90 remain additive).
- Joint multilocus reference frequencies for exact epistatic orthogonality
  under LD.

---

## Open questions

### D1. Plan scope — storage only, or storage plus a non-binding downstream sketch?

| Option | Notes |
|---|---|
| **Part 1 locked / Part 2 sketched** ← recommended, and what v2 does | Keeps the whole picture visible while making clear that only the schema is being decided. Every Part 2 item is labelled "own plan" |
| Storage only (the review's §4 position) | Cleaner decision boundary, but the schema's acceptance criteria become unfalsifiable — you cannot tell whether a storage design is sufficient without knowing what will read it |
| Keep v1's mixing | Rejected. It genuinely obscured which decisions were being locked |

### D2. How much to build now vs. reserve

This is the biggest question, because the design goes from 1 table to 16, each
needing six registry entries.

| Option | Tables | Delivers | Notes |
|---|---|---|---|
| **Tier 2 — everything except reference states** ← recommended | 15 | Additive, dominance, all-order epistasis, arbitrary genotype surfaces, full scopes | Surfaces are the *reason* this project exists ("users can input it now even if I can't sample it"). Postponing them postpones the payoff |
| Tier 1 — coefficient terms only | 11 | Everything except hand-entered genotype tables | The review's "safe subset". Acceptable **only if** the surface table shapes are agreed and reserved in this document now — which they are |
| Tier 0 — build v1's simpler schema, document v2 as the target | 5 | Additive + dominance | **Not recommended, and the review is right to call it unsafe.** `dosage_key` and single-`line_name` are exactly the two things that would force the next rewrite |

Recommendation: **Tier 2**, but build a registration helper first so that
sixteen tables cost one registration call each rather than six edits each.

### D3. The scope resolver — design it now or defer it?

The scope *tables* fix cardinality, which cannot be retrofitted. The *matching
and precedence rules* (which scope wins when several match) are undesigned, and
the review says so plainly.

| Option | Notes |
|---|---|
| **Defer, gated on a truth table** ← recommended | Write the common / A / B / A-B / reciprocal / unknown / haploid / diploid / tetraploid truth table during Phase A. Accept the scope tables only if every row has an unambiguous representation. Build the resolver in its own plan |
| Design it now, inline | Expands this plan substantially and delays the storage lock, for rules that only matter once non-additive evaluation exists |
| Ship only `common` + `allele_copy` scopes | Smallest resolver (reproduces today's behavior exactly), but you cannot then *store* an A/B dominance effect — losing the multi-line goal |

### D4. Existing populations — migrate or regenerate?

`CLAUDE.md` says no one is a downstream user before 1.0.0, so no *user* data
migration is owed. But you may have `.duckdb` files from your own runs.

| Option | Notes |
|---|---|
| **Regenerate** ← recommended if your populations are cheap to rebuild | No migration code at all. Simplest, and consistent with "Break Freely" |
| Migrate in place | The review's §12 Phase C. Needs care: if two current rows disagree on `base_allele_freq` for the same `(base context, locus)`, migration must **stop and report**, never pick one. Provenance is unknowable for existing rows, so migrated references get `source_type = 'legacy'` with values preserved and no invented provenance |

**I need your answer on this one** — it changes whether Phase C writes a
migrator or just a writer.

### D5. Contrast values — always materialized, or computed for built-in codings?

| Option | Notes |
|---|---|
| **Always materialized; built-in codings are generators run at write time** ← recommended | Evaluation is uniform for built-in and custom contrasts. Invariant 3 (self-describing storage) holds unconditionally. Cost is bounded: shared template states plus one contrast per locus only for frequency-dependent codings |
| Compute built-ins at read time from `coding_name` + reference | Fewer rows, but a restored database is only interpretable by a version that still knows the coding name — the exact failure mode §3.1 identifies |
| Hybrid, per contrast | Two evaluation paths to test and keep consistent, for no benefit at these row counts |

### D6. General writer API

| Option | Notes |
|---|---|
| **Design after Phase A fixtures exist** ← recommended | The fixtures reveal what the input shape must express. `define_additive_effects()` remains the front door for the common case regardless |
| `define_genome_effects()` now | v1's proposal. The review calls it "likely keep", but with normalized states/contrasts/scopes its argument list is no longer obvious |
| Split writers per gene action | Multiplies entry points and re-fragments what the schema just unified |

### D7. Imprinting × non-additive terms

`expressed_parent` restricts to one `parent_origin` (`add_tbv.R:191-196`), which
is incoherent for a term needing both copies. Note this is now *partly*
addressed by the schema: `parent_origin` exists on both state alleles and scope
origins, so parent-of-origin-specific dominance is **storable**.

| Option | Notes |
|---|---|
| **Error at the writer** ← recommended for v1 of the evaluator | Reject a non-additive term on a trait with `expressed_parent != "both"`, with a message pointing at the scope/state mechanism as the supported way to express it |
| Ignore | Silently applies non-additive terms across both copies while additive stays restricted. Cannot be right |
| Model it via scopes | The schema already permits it; needs the resolver (D3) first |

### D8. Default contrast coding

The review says not to lock a default here, since storage supports multiple
explicit bases. Partially agreed — the *storage* is basis-agnostic, but
`define_additive_effects()` still needs a default when the user says nothing.

| Option | Notes |
|---|---|
| **`'functional_diploid_v1'`** ← recommended | `a`/`d` are what breeders think and publish in; frequency-independent, so coefficients stay meaningful as the population drifts. Note the current engine's centered coding is *already* frequency-dependent, so this is a change in default, not just a label |
| `'cockerham_hwe_v1'` | Matches today's actual behavior most closely and partitions variance cleanly at the base; but coefficients drift out of orthogonality across generations |
| No default; require an explicit `coding_name` | Safest, least ergonomic. Reasonable if the two above are genuinely both common in your work |
