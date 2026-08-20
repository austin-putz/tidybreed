# Plan: organize `schema()` output + rename the phenotype-effect tables

Status: draft for review. Two related pieces of work:

1. **Ordering** — `schema()` prints tables in an unstable, meaningless order.
2. **Naming** — `trait_effects` / `trait_random_effects` are misnamed; they are
   phenotype-layer tables. (Verified below — the rename is justified.)

They are separable, but doing the rename first means the ordering table below is
written once against final names.

---

## Part 1 — Why the current order looks random

`schema()` (`R/schema.R:549`) iterates `pop$tables` verbatim:

```r
user_tables <- setdiff(pop$tables, "_schema_meta")
```

`pop$tables` is populated two different ways:

| Entry point | Source | Resulting order |
|---|---|---|
| `open_pop()` | hard-coded `tables_created` vector (`R/open_pop.R:162`), plus appends from `define_genome()`, `define_trait()`, etc. | creation order — depends on which `define_*` calls were made, and in what order |
| `restore_pop()` | `DBI::dbListTables(db_conn)` (`R/restore_pop.R:80`) | DuckDB catalog order (effectively alphabetical) |

So **the same `.duckdb` file prints in a different order depending on whether you
just built it or restored it.** That instability is the real defect; the aesthetic
grouping is the improvement layered on top.

---

## Part 2 — Recommended ordering: group by pipeline stage

Order tables the way a user builds a population, and print section headers. This
mirrors the package's own `define_` (configuration) → `add_` (output) split, so
the printout doubles as a map of the workflow.

### Proposed groups and in-group order

In-group order is workflow order, **not** alphabetical.

| # | Group | Tables |
|---|---|---|
| 1 | **Genome** | `genome_meta`, `genome_map`, `chr_inheritance`, `chr_recombination` |
| 2 | **Founders** | `founder_haplotypes` |
| 3 | **Individuals** | `ind_meta`, `ind_haplotype`, `ind_genotype`, `ind_crossover` |
| 4 | **Genetic model** | `trait_meta`, `trait_var_comp`, `genome_effects` |
| 5 | **Observation model** | `phenotype_meta`, `phenotype_components`, `phenotype_var_comp`, `phenotype_effects`*, `phenotype_random_effects`* |
| 6 | **Selection** | `index_meta` |
| 7 | **Results** | `ind_tbv`, `ind_phenotype`, `ind_ebv`, `ind_index`, `ind_true_index` |
| 8 | **User tables** | anything unrecognized (e.g. from `define_table()`), alphabetical |
| 9 | **System** | `_schema_meta` — hidden unless `include_system = TRUE` |

\* post-rename names; see Part 5.

**Open judgment call:** `genome_effects` is genome-shaped by name but is trait
configuration by function (written by `define_additive_effects()`, keyed by
`trait_name`). Placing it in **Genetic model** next to `trait_meta` matches how
you reason about it; placing it in **Genome** matches the name. Recommendation:
Genetic model. Flag if you disagree — it is a one-line change.

Groups 4 and 5 are exactly the two-layer split from `CLAUDE.md`, which is the
main argument for this scheme over any lexical one: the print reinforces the
package's central design boundary.

### Sketch of the output

```
── Schema: MySim ───────────────────────────────────────────────────────
  Use describe_table(pop, "name") for column-level details.

  Genome
    genome_meta            1,000   7  Locus-level metadata. One row per locus…
    genome_map             1,000   7  Genetic map in long format. One row per…
    chr_inheritance           10   5  Per-chromosome copy counts, keyed by off…
    chr_recombination         10   4  Per-chromosome recombination, keyed by p…

  Individuals
    ind_meta                 550   8  Individual-level metadata. One row per i…
    ind_haplotype          1.10M   7  Phased haplotypes in long format. One ro…
    + 2 empty: ind_genotype, ind_crossover

  Genetic model
    trait_meta                 2   6  Genetic component trait definitions…
    trait_var_comp             4   5  Genetic variance component storage…
    genome_effects           120   7  QTL effect data in long format…
  …
```

---

## Part 3 — Where the ordering lives (pick one)

### Option A — hard-coded vector in R  ← recommended

An internal `.schema_table_order()` returning an ordered
`list(group_name = c(table, table, …))`, with unrecognized tables falling through
to **User tables**.

- Zero storage, single source of truth, trivially reviewable in a diff.
- A newly added table that someone forgets to register degrades *visibly*
  (shows up under "User tables") rather than silently misplacing itself.
- Cost: one more place to update when adding a table. Acceptable, and arguably a
  feature — it forces a deliberate decision about which layer the table belongs to.

### Option B — derive groups from the name prefix

Group by `genome_*`, `chr_*`, `ind_*`, `trait_*`, `phenotype_*`, `index_*`,
`founder_*`; alphabetical within prefix.

- Self-maintaining; new tables slot in with no code change.
- But the grouping is lexical, not conceptual: it splits `genome_effects` away
  from `trait_meta`, and files `ind_tbv` next to `ind_haplotype` (output next to
  raw genome data). It also cannot express workflow order within a group.
- Reasonable fallback *inside* Option A for unrecognized tables.

### Option C — `table_group` + `display_order` columns on `_schema_meta`

Store the grouping on the existing table-level `_schema_meta` rows.

- Most consistent with the "metadata lives in the database, not the R object"
  principle in `CLAUDE.md`, and it travels with the `.duckdb` file.
- Lets `define_schema_description()` gain a `table_group =` argument so users can
  place their own `define_table()` tables deliberately.
- Cost: a schema change to `_schema_meta` plus handling in `restore_pop()` for
  files written before the change.

**Recommendation:** A now, with C as a later additive enhancement if user-table
placement turns out to matter. A and C compose cleanly — the hard-coded vector is
the default, and a non-NULL `_schema_meta.table_group` overrides it.

---

## Part 4 — Print-method refinements to pair with the ordering

All in `print.tidybreed_schema()` (`R/schema.R:606`) unless noted.

1. **Add a `table_group` column to the returned tibble.** The grouping should
   exist as data, not only as printed text, so users can `filter()` / regroup.
   `schema()` returns a `tbl_df`; the new column is free to consumers who ignore it.

2. **Collapse empty tables within each group.** On a freshly built population most
   tables have zero rows and they bury the signal. Print
   `+ 3 empty: ind_ebv, ind_index, ind_crossover` as the group's last line. Add
   `show_empty = FALSE` as the argument to expand them.

3. **Abbreviate large counts** — `1.10M` rather than `1,100,000`. `ind_haplotype`
   alone currently dictates the width of the `Rows` column for every table.

4. **Suppress the group header when a group has no tables at all.** A population
   with no traits defined should not print an empty "Genetic model" heading.

5. **Keep `_schema_meta` hidden by default** (current behavior) but expose it via
   `include_system = TRUE` rather than making it permanently invisible.

---

## Part 5 — `schema(pop, order = ...)`

```r
schema(pop, order = c("pipeline", "name", "rows"), show_empty = TRUE,
       include_system = FALSE)
```

- `"pipeline"` (default) — grouped as above, with headers.
- `"name"` — flat alphabetical, no headers. Predictable lookup.
- `"rows"` — flat, descending row count. Answers "what is actually big here?",
  which is a genuine question for file-based DuckDB populations.

Headers print only for `"pipeline"`. `table_group` stays in the returned tibble
for all three.

---

## Part 6 — The rename: `trait_effects` → `phenotype_effects`

### Verification — yes, rename them

Both tables are **already keyed by `phenotype_name`**, including in their primary
keys (`R/define_trait.R:181-208`):

```sql
CREATE TABLE trait_effects (
  phenotype_name    VARCHAR NOT NULL,
  effect_name       VARCHAR NOT NULL,
  ...
  PRIMARY KEY (phenotype_name, effect_name)
)

CREATE TABLE trait_random_effects (
  phenotype_name VARCHAR NOT NULL,
  effect_name    VARCHAR NOT NULL,
  level          VARCHAR NOT NULL,
  ...
  PRIMARY KEY (phenotype_name, effect_name, level)
)
```

Everything else about them is observation-layer too:

- Written by `define_effect_fixed_class()`, `define_effect_fixed_cov()`,
  `define_effect_random()` — all phenotype-model configuration.
- Read by `add_phenotype()` / `phenotype_helpers.R` to build the fixed and random
  shifts added to the phenotype.
- `trait_random_effects` holds *sampled draws* per phenotype × effect × level —
  pure observation-layer noise, no genetic content whatsoever.
- Their own `_schema_meta` descriptions already describe them in phenotype terms
  (`R/schema.R:363-402`).

A column rename `trait_name` → `phenotype_name` already happened (there is a
migration for it at `R/define_trait.R:369-390`); the **table** names were simply
left behind. Under the naming-consistency rules in `CLAUDE.md`, the names are
wrong today, and the two-layer split makes `trait_*` actively misleading — it
implies `trait_meta` is the parent, when the FK is to `phenotype_meta`.

### Stale references found while verifying — FIXED in v0.63.1

These were fixed and committed separately from this plan; they are recorded here
because they are what proved the rename argument.

- ~~**`R/sql_utils.R:155` — real bug.**~~ `TABLE_ROW_KEYS$trait_effects` listed
  `trait_name` after the column became `phenotype_name`, so `remove_rows()` on
  this table always aborted with `"missing key column(s): trait_name"`.
- ~~**`R/sql_utils.R:98`**~~ — `TABLE_RESERVED_COLS$trait_effects` listed the same
  stale `trait_name` and omitted `poly_order` / `null_class_action`, so
  `mutate_table()` let a user overwrite the `phenotype_name` foreign key.
- ~~**`CLAUDE.md`**~~ — the `trait_effects` schema table documented `trait_name`
  and was missing both newer columns.
- `tests/testthat/test-schema-registries.R` now audits all three registries
  against the live schema, so this class of drift fails at test time.
- **`R/define_trait.R:369-390`** — the `trait_name` → `phenotype_name` migration
  block. Per the pre-1.0.0 policy ("no compatibility shims, no legacy paths"),
  this should be deleted rather than extended for the table rename.

### Newly found while auditing (not yet addressed)

The registry audit surfaced three tables absent from `TABLE_RESERVED_COLS`
entirely: `founder_haplotypes`, `phenotype_components`, and
`trait_random_effects`. When a table has no entry, `mutate_table()` blocks
nothing, so a user can overwrite e.g. `trait_random_effects.draw_value` or
`phenotype_components.contributor_type` in place. Each needs a judgment call
about whether user columns belong on that table at all, so they were left alone
rather than fixed blind.

### Rename scope

Roughly 60 non-plan references across the package:

| Area | Files |
|---|---|
| Schema creation | `R/define_trait.R` (~16 refs, incl. the migration block to delete) |
| Readers/writers | `R/add_phenotype.R`, `R/phenotype_helpers.R`, `R/effect_helpers.R`, `R/define_effect_random.R`, `R/define_effect_fixed_class.R`, `R/define_effect_fixed_cov.R` |
| Descriptions | `R/schema.R` (~22 refs in `.trait_layer_descriptions()`) |
| Constants | `R/sql_utils.R` (reserved cols, row keys, table list at `:216`) |
| Export | `R/blupf90_helpers.R` |
| Archiving | `R/archive_replicate.R:104` |
| Docs/tests | `man/*.Rd` (regenerate), `tests/testthat/`, `inst/swine/`, `CLAUDE.md` |

Mechanical, but touches the phenotype path broadly — worth its own commit,
separate from the printing change.

### Also worth considering in the same pass

`.trait_layer_descriptions()` in `R/schema.R:344` currently describes
`phenotype_meta`-adjacent tables *and* `ind_ebv` / `index_meta` / `ind_index`,
while `.core_layer_descriptions()` (`R/schema.R:175`) describes
`phenotype_meta`, `phenotype_components`, and `phenotype_var_comp`. The three
description functions do not partition along the same lines as the display groups
in Part 2. Re-splitting them to match (`.genome_*`, `.individual_*`,
`.genetic_model_*`, `.observation_model_*`, `.results_*`) would make the group
vector and the descriptions the same list read twice, which is harder to let drift.

Optional; not required for either piece of work above.

---

## Suggested sequencing

1. ~~Fix `R/sql_utils.R:98` and `:155`~~ — done in v0.63.1.
2. Rename the two tables; delete the `define_trait.R` migration block; update
   `CLAUDE.md`; regenerate `man/`.
3. Implement Part 2 + Part 3 Option A + Part 4.
4. Add the `order` / `show_empty` / `include_system` arguments (Part 5).
5. Optionally re-split the description helpers (Part 6, last section).

Version bump + `NEWS.md` entry per `CLAUDE.md` before commit.
