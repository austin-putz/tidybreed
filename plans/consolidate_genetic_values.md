# Consolidate genetic values — replace `ind_tbv` with `ind_tgv`

**Status:** design proposal, nothing implemented. **Scheduled immediately after**
`plans/update_genome_effects_v4.md` (v4.9) ships. **Updated 2026-09-10** for v4.9's
order-explicit `component_name` vocabulary and its `effect_owner` rename.
**Created:** 2026-09-06, splitting deferrable work out of the genome-effects plan.

---

## Why this is a separate plan

The genome-effects plan (v4.9) creates `ind_tgv` and writes components to it. It does
**not** touch `ind_tbv`. That split is deliberate and follows the decision rule both
plans run on:

> Pre-build a dimension only when its absence would force existing rows to change
> shape. Everything else is a later child table — or a later plan.

`component_name` had to exist in v4.9, because adding it afterwards would change
`ind_tgv`'s unique key and reshape every row. Everything in *this* plan is rows,
values, and migration. None of it reshapes anything v4.9 creates, and none of it is
needed for v4.9's effects to be correct or executable.

Keeping it separate also keeps v4.9 honest about scope: replacing `ind_tbv` touches
**10 R files, 15 test files, 10 man pages, and 2 vignettes**, none of which have
anything to do with whether a dominance coefficient evaluates correctly.

---

## The problem this plan closes

After v4.9 the additive genetic value lives in two places, and **they are not always
the same number**:

| | Meaning |
|---|---|
| `ind_tbv.tbv_value` | Breeding value from the reserved `generated_additive_tbv` effect owner |
| `ind_tgv` row `component_name = 'order1_additive'` | Sum of **every** additive-structured term in the model, including custom ones written through `define_genome_effects()` |

Identical whenever a user only ever calls `define_additive_effects()` — which is the
common case, and exactly why the discrepancy is dangerous. A user who adds one custom
additive term gets two tables that disagree, with no error and no warning beyond the
one v4.9 Q1 proposes.

`CLAUDE.md` is explicit that redundant paths are technical debt, not harmless history.
This plan removes the redundancy rather than documenting it.

---

## Target state

`ind_tbv` is **deleted**. `ind_tgv` is the single table of true genetic values.
`add_tbv()` is **deleted**. `add_tgv()` computes and writes every component in one
pass over the haplotypes and effects — which is also strictly less work than two
functions walking the same data.

```sql
-- unchanged from v4.9; no DDL change in this plan
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

**The breeding value is `component_name = 'order1_additive'`.** No separate table, no
separate function, no second definition of the same quantity.

Per `CLAUDE.md` pre-1.0 policy, the old names go completely: no `ind_tbv` view, no
`add_tbv()` alias, no deprecated wrapper. Man pages, roxygen examples, tests, and
vignettes move to the new path.

---

## Task 1 — Component vocabulary

v4.9 ships a deliberately minimal, order-explicit set: `'order1_additive'`,
`'order1_dominance'`, `'order1_other'`, `'interaction'`. This plan refines it.
**Every refinement is rows, not DDL.**

### 1a. Split `'interaction'` by interaction type

`'interaction'` collapses A×A, A×D, D×A, and D×D into one bucket. The member contrasts
already distinguish them, and members are canonically ordered by `locus_id`, so the
label is mechanically derivable:

`'add_add'` · `'add_dom'` · `'dom_dom'` · higher orders by the same rule

Note A×D and D×A are the *same* term under v4.9's canonical ordering — the label must
be read as "the lower-`locus_id` member is additive, the higher is dominant", not as a
biological direction. Decide whether third-order and above get explicit labels or fall
back to `'interaction'` with the order exposed in a view.

### 1b. Classify hand-entered surfaces

`'order1_other'` is v4.9's honest placeholder for an order-1 `indicator` term. A raw
genotype surface mixes additive and dominance by construction; separating them is an
orthogonal projection against base allele frequencies.

| Option | Notes |
|---|---|
| **Writer declares it** ← recommended | `define_genome_effects(..., component = "dominance")`. This is a simulator: the user is defining ground truth, so they are the authority on what they meant. Cheap, no statistics |
| Keep `'order1_other'` as the permanent answer | Honest but unhelpful — a user who typed a dominance surface sees no dominance component |
| Orthogonal projection at write time | Statistically correct and the only option that makes the label a real decomposition. Needs base frequencies, and under LD joint frequencies. **Its own project** — see "Related, not included" |

---

## Task 2 — Delete `ind_tbv` and `add_tbv()`

### Blast radius (verified at `4e9a72b`)

**R (10):** `add_tbv.R` (deleted), `add_phenotype.R`, `add_index.R`,
`archive_replicate.R`, `define_trait.R`, `formula_helpers.R`, `remove_rows.R`,
`schema.R`, `sql_utils.R`, `tidybreed_pop.R`

**Registries (6 per table):** `sql_utils.R:97,137,164,256,267`; `schema.R`
`.schema_table_order()` and `.ind_descriptions()`. Guarded by
`test-schema-registries.R` and `test-schema-print.R`.

**Tests (15):** incl. `helper-parity.R`, `test-add_tbv.R`, `test-add_tbv_index.R`,
`test-add_phenotype.R`, `test-phenotype_composite.R`, `test-remove_rows.R`

**Man (10) · Vignettes (2):** `tidybreed-introduction.Rmd`,
`swine/swine-time-based-age-at-puberty-sex-semen.R`

### 2a. `add_index()`

Auto-detection maps table → value column (`add_index.R:120`). `ind_tbv` → `"tbv_value"`
becomes `ind_tgv` → `"tgv_value"`, and callers must filter the component:

```r
pop |> get_table("ind_tgv") |> filter(component_name == "order1_additive") |> add_index("meat")
```

`add_index()` already errors when an individual has more than one value per trait, so
an unfiltered call fails loudly rather than silently summing components. **Confirm
that error fires on the component dimension** — it is the safety net for this change.

### 2b. `ind_true_index`

Computed by `add_tbv()` today from TBVs. Moves to `add_tgv()` with a component
parameter defaulting to `'order1_additive'` — selection indices are conventionally on breeding
values. Allowing `'total'` is a one-line generalization worth taking.

### 2c. The oracle

`tests/testthat/test-add_tbv.R:16-40` independently recomputes the additive formula
from first principles. It is **not** pre-change golden output and must survive, retargeted
at `ind_tgv` filtered to `'order1_additive'`. It is the only thing proving the migration
preserved values.

---

## Task 3 — Composite phenotypes should read the total, not the breeding value

`add_phenotype()` reads TBVs in four places to assemble composite and SGE phenotypes
(`add_phenotype.R:754-763, 979-983, 1025-1029`, plus `upsert_ind_tbv()` at `:1119`) —
`.fetch_contributor_tbvs()` and `.assemble_composite_tbv()`.

**This is a semantic bug the moment dominance exists.** A dam's realized maternal
ability includes her dominance deviation; a group-mate's realized social effect
includes theirs. Using only the additive component understates every contributor.

Change contributor lookups to the **derived total**, not `'order1_additive'`.

This is the one task in this plan that is genuinely coupled to something else — the
deferred phenotype-integration plan, which activates
`phenotype_components.genome_effect_types` (dead since `open_pop.R:334`). Two options:

| Option | Notes |
|---|---|
| **Do the swap here, phenotype integration later** ← recommended | The swap is correct independently: contributors should always have been total genetic value, and today total *equals* additive, so the change is a no-op until dominance exists. Making it now means phenotype integration inherits correct semantics instead of fixing them |
| Defer both together | Fewer moving parts, but leaves a known-wrong contributor lookup in the tree for a whole cycle |

---

## Acceptance gates

1. `ind_tbv` and `add_tbv()` do not appear anywhere in `R/`, `man/`, `tests/`, or
   `vignettes/`.
2. The retargeted first-principles oracle agrees with `add_tgv()`'s `'order1_additive'`
   component for a purely additive model.
3. For an additive-only model, `'order1_additive'` equals the derived total and no other
   component row is written.
4. Components sum to the derived total for a mixed additive + dominance + epistatic
   model, hand-computed.
5. A custom additive term written through `define_genome_effects()` appears in the
   `'order1_additive'` component — the discrepancy this plan exists to close.
6. `add_index()` on an unfiltered `ind_tgv` errors rather than summing components.
7. `ind_true_index` from `'order1_additive'` reproduces the pre-migration value.
8. A composite/maternal phenotype uses the contributor's total; with a dominance model
   present, it differs from the additive-only result by exactly the contributor's
   dominance deviation.
9. `archive_replicate()` and `remove_rows()` work on `ind_tgv` including `replicate`.
10. A×A, A×D, and D×D terms land in distinct component rows (Task 1a).
11. `schema()` and `describe_table()` render `ind_tgv` with all registries present.

---

## Open questions

### Q1 — Does `'total'` ever get stored?

Recommended: **no**, keep it derived. A stored total makes `SUM(tgv_value)`
double-count, which is a permanent footgun for every ad-hoc query a user writes.
Reconsider only if the view proves too slow on large populations, in which case a
materialized total is a performance decision with an explicit refresh contract, not a
schema decision.

### Q2 — Third-order and higher interaction labels

Recommended: **label explicitly up to order 3** (`add_add_add`, etc.), then fall back
to `'interaction'` with the order exposed in a view. Naming every combination at order 5
is a combinatorial vocabulary nobody reads. Low stakes.

### Q3 — Should `component_name` be constrained in SQL?

The set grows in this plan (Task 1a) and may grow again. A `CHECK` makes every
addition a DDL migration; R validation keeps it a one-line change. Recommended:
**R-validated closed set**, consistent with how v4.9 treats `contrast_name`'s
*membership* while still declaring the ones that are genuinely fixed.

---

## Related, not included

- **Orthogonal statistical decomposition** (`V_A` / `V_D` / `V_I`, NOIA-style
  projection of arbitrary surfaces onto a basis). This is the only thing that would
  make component labels a true variance decomposition rather than a record of how the
  model was declared. Needs base frequencies and, under LD, joint multilocus
  frequencies. **Its own plan**, and a prerequisite for any variance-targeting feature
  beyond additive.
- **Phenotype integration** — activating `phenotype_components.genome_effect_types`.
- **`ind_egv`** — an estimated counterpart to `ind_tgv`, if non-additive genomic
  prediction is ever added. `add_ebv()` and `ind_ebv` are untouched by this plan.
