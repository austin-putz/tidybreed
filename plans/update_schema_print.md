# Plan: organize `schema()` output + rename the phenotype-effect tables

Status: **IMPLEMENTED at v0.64.0 (2026-08-31), uncommitted pending review.**
All six sequencing steps are done; see **Completion notes** at the end for what
landed, what was verified against the plan's own measurements, where the
implementation deviates from the plan and why, and what was deliberately left
out. Two related pieces of work:

1. **Ordering** — `schema()` prints tables in an unstable, meaningless order.
2. **Naming** — `trait_effects` / `trait_random_effects` are misnamed; they are
   phenotype-layer tables. (Verified below — the rename is justified.)

They are separable, but doing the rename first means the ordering table below is
written once against final names.

**Decisions (2026-08-20):**

- `genome_effects` is grouped under **Genome**, not Genetic model (Part 2).
- Ordering lives in a **hard-coded vector in R**, Option A (Part 3).

**Decisions (2026-08-31)** — closing everything left open at the first review:

- `show_empty` defaults to **`FALSE`** (empty tables collapse to one line);
  `show_empty = TRUE` expands them. Part 6's signature said `TRUE` and was wrong
  (Part 4, Part 6).
- The three tables missing from `TABLE_RESERVED_COLS` get **all columns
  reserved**, shipped as its own commit *before* the rename (Part 7).
- The rename gets **no `restore_pop()` migration**, and the legacy migration code
  that already exists is **deleted in the same commit** (Part 7).
- `sizes = TRUE` on an in-memory population **warns and omits the column**; only
  `order = "size"` without `sizes = TRUE` errors (Part 5d, Part 6).
- The description-helper re-split is **folded into the rename commit** rather than
  left as optional (Part 7, last section).

---

## Part 1 — Why the current order looks random

`schema()` (`R/schema.R:545`) iterates `pop$tables` verbatim:

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
| 1 | **Genome** | `genome_meta`, `genome_map`, `genome_effects`, `chr_inheritance`, `chr_recombination` |
| 2 | **Founders** | `founder_haplotypes` |
| 3 | **Individuals** | `ind_meta`, `ind_haplotype`, `ind_genotype`, `ind_crossover` |
| 4 | **Genetic model** | `trait_meta`, `trait_var_comp` |
| 5 | **Observation model** | `phenotype_meta`, `phenotype_components`, `phenotype_var_comp`, `phenotype_effects`*, `phenotype_random_effects`* |
| 6 | **Selection** | `index_meta` |
| 7 | **Results** | `ind_tbv`, `ind_phenotype`, `ind_ebv`, `ind_index`, `ind_true_index` |
| 8 | **User tables** | anything unrecognized (e.g. from `define_table()`), alphabetical |
| 9 | **System** | `_schema_meta` — hidden unless `include_system = TRUE` |

\* post-rename names; see Part 7.

**Coverage check (2026-08-31):** groups 1-9 above partition the 24 entries of
`SYSTEM_TABLES` (`R/sql_utils.R:209`) exactly — 5 + 1 + 4 + 2 + 5 + 1 + 5 + 1,
plus `_schema_meta`. No system table falls through to **User tables**, so a table
appearing there at runtime is a real signal that something new went unregistered.

**Decided:** `genome_effects` goes in **Genome**, not Genetic model. It is
locus-keyed (one row per locus x trait x effect type x line, FK on `locus_name`),
so it reads naturally next to `genome_meta` and `genome_map` — the three tables
that are all "one row per locus, plus dimensions." Genetic model is then exactly
the two trait-keyed configuration tables.

Groups 4 and 5 still carry the two-layer split from `CLAUDE.md`, which is the
main argument for this scheme over any lexical one: the print reinforces the
package's central design boundary. `genome_effects` sitting under Genome does not
weaken that — it is where the QTL *live*, and `trait_meta` is where the traits are
*defined*.

### Sketch of the output

```
── Schema: MySim ───────────────────────────────────────────────────────
  Use describe_table(pop, "name") for column-level details.

  Genome
    genome_meta            1,000   7  Locus-level metadata. One row per locus…
    genome_map             1,000   7  Genetic map in long format. One row per…
    genome_effects           120   7  QTL effect data in long format…
    chr_inheritance           10   5  Per-chromosome copy counts, keyed by off…
    chr_recombination         10   4  Per-chromosome recombination, keyed by p…

  Individuals
    ind_meta                 550   8  Individual-level metadata. One row per i…
    ind_haplotype          1.10M   7  Phased haplotypes in long format. One ro…
    + 2 empty: ind_genotype, ind_crossover

  Genetic model
    trait_meta                 2   6  Genetic component trait definitions…
    trait_var_comp             4   5  Genetic variance component storage…
  …
```

---

## Part 3 — Where the ordering lives

### Option A — hard-coded vector in R  ← DECIDED

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

**Decided: Option A.** Options B and C are kept above as the record of what was
weighed, not as live alternatives. C stays available as a purely additive
enhancement if user-table placement ever matters — the two compose cleanly, with
the hard-coded vector as the default and a non-NULL `_schema_meta.table_group`
overriding it. No `_schema_meta` schema change is needed for the work below.

---

## Part 4 — Print-method refinements to pair with the ordering

All in `print.tidybreed_schema()` (`R/schema.R:599`) unless noted.

1. **Add a `table_group` column to the returned tibble.** The grouping should
   exist as data, not only as printed text, so users can `filter()` / regroup.
   `schema()` returns a `tbl_df`; the new column is free to consumers who ignore it.

2. **Collapse empty tables within each group.** On a freshly built population most
   tables have zero rows and they bury the signal. Print
   `+ 3 empty: ind_ebv, ind_index, ind_crossover` as the group's last line.
   **Decided:** the argument is `show_empty`, defaulting to **`FALSE`**
   (collapsed); `show_empty = TRUE` prints every empty table as its own row. An
   earlier draft of the Part 6 signature had the default as `TRUE`, which
   contradicted the rationale above — collapsing *is* the feature, so it is on by
   default and the argument exists to turn it off.

3. **Abbreviate large counts** — `1.10M` rather than `1,100,000`. `ind_haplotype`
   alone currently dictates the width of the `Rows` column for every table.

4. **Suppress the group header when a group has no tables at all.** A population
   with no traits defined should not print an empty "Genetic model" heading.

5. **Keep `_schema_meta` hidden by default** (current behavior) but expose it via
   `include_system = TRUE` rather than making it permanently invisible.

---

## Part 5 — Table and database size in bytes

Row counts answer "how much data", but not "how much disk". For a file-based
package where the whole point is simulations larger than RAM, the second question
is the operational one. **Both are possible, but they are not equally cheap** —
the database total is free and exact, per-table is neither.

All figures below were measured (not assumed) against DuckDB via
`dev/`-style probes on a 60-founder / 2,000-locus population.

### 5a. Whole-database size — easy, exact, do it

`PRAGMA database_size` returns one row of already-formatted strings:

```r
ds <- DBI::dbGetQuery(conn, "PRAGMA database_size")
# database_size "4.0 MiB"   block_size 262144   total_blocks 16   used_blocks 16
# wal_size      "0 bytes"   memory_usage "2.1 MiB"   memory_limit "51.1 GiB"
```

- **File-backed pop**: report `database_size` + `wal_size`. Confirmed to match
  `file.size()` on the `.duckdb` (4.0 MiB vs 4,206,592 bytes).
- **In-memory pop** (`db_name = ":memory:"`): `database_size` is `"0 bytes"` and
  `block_size` is `0` — there is no file. Report `memory_usage` instead and label
  it as such. `pop$db_path` is the literal `":memory:"`, so branch on that.
- The WAL matters: it is nonzero until a `CHECKPOINT`, and a user looking at
  their file on disk will see both. Reporting only `database_size` under-reports.

Best home is the `schema()` header, where it costs one line and no per-table work:

```
── Schema: MySim ──────────────────── 24 tables · 4.0 MiB on disk ──
```

### 5b. Per-table size — possible, but qualified

DuckDB has **no direct per-table byte count**. Note the trap:
`duckdb_tables().estimated_size` sounds right but is estimated *cardinality*
(rows), not bytes — it read `189` for a 189-row `_schema_meta`.

The workable route is `PRAGMA storage_info('<table>')`, counting distinct
non-negative `block_id`s and multiplying by `block_size`. Measured results:

| Table | Rows | Blocks | Approx |
|---|---|---|---|
| `ind_haplotype` | 240,000 | 3 | 0.75 MB |
| `founder_haplotypes` | 120,000 | 2 | 0.50 MB |
| `genome_meta` | 2,000 | 1 | 0.25 MB |
| `chr_inheritance` | 5 | 1 | 0.25 MB |

Three caveats, all verified. **Decided: #2 and #3 are acceptable as long as the
printed output states them** — the numbers are useful even when coarse, so they
are a disclosure problem, not a blocker. See 5d for the exact footnote. #1 is
different in kind: it changes the database rather than misleading about it, so it
still governs the opt-in in 5c.

1. **It requires a `CHECKPOINT` first.** Before one, recently written data has
   `block_id = -1` and the table reports **0 bytes** — three of four tables read as
   zero on the first probe. `CHECKPOINT` is a *write* (it flushes the WAL), which
   is a surprising side effect for a read-only-looking `schema()` call.
   *Not* solved by a footnote — this one keeps `sizes` opt-in.
2. **Granularity is one block = 256 KB.** Every table with any data floors at
   0.25 MB, so `chr_inheritance` (5 rows) and `genome_meta` (2,000 rows) are
   indistinguishable. For a package with ~24 tables, most of which are small, the
   column would read as a wall of `0.25 MB`. → **footnoted**
3. **Per-table sizes do not sum to the file size** — 11 blocks attributed across
   tables vs 16 `used_blocks`, the remainder being catalog and header overhead.
   Blocks are *not* shared between tables (measured: 0 blocks claimed by more than
   one table), so attribution is at least sound as far as it goes. → **footnoted**

**The rejected alternative**: estimating `rows x sum(column type widths)` needs no
`CHECKPOINT` and has no quantization, but it ignores compression — and this schema
compresses extremely well (`ind_haplotype` stores `allele` as RLE and `id_ind` as
Dictionary; 240,000 x 7 columns lands in 0.75 MB, not the ~3.4 MB a naive width
sum predicts). A confidently-wrong number is worse than none.

### 5c. Recommendation

- **Do 5a unconditionally.** One line in the header, exact, free, no side effects.
- **Make 5b opt-in**: `schema(pop, sizes = FALSE)`. When `TRUE`, issue the
  `CHECKPOINT`, add a `size_bytes` column, and print the footnote from 5d.
  Never on by default — a silent `CHECKPOINT` inside a display function is the
  kind of hidden write this package should not have. This is caveat #1, and it is
  the *only* reason the argument is opt-in; #2 and #3 are handled by disclosure.
- **`order = "size"`** becomes meaningful only when `sizes = TRUE`; otherwise
  `order = "rows"` is the honest proxy and should stay the default answer to
  "what is big here?"
- Both belong in the returned tibble as data (`size_bytes`), formatted only at
  print time, consistent with the `table_group` decision in Part 4.

### 5d. The footnote (required whenever a Size column prints)

`print.tidybreed_schema()` must emit this whenever `sizes = TRUE`. It is not
optional and has no suppression argument — the column is only defensible because
the caveat travels with it.

```
── Schema: MySim ──────────────────────── 24 tables · 4.0 MiB on disk ──
  Use describe_table(pop, "name") for column-level details.

  Genome
    genome_meta            2,000   7   0.25 MB  Locus-level metadata…
    genome_map             2,000   7   0.25 MB  Genetic map in long format…
    chr_inheritance            5   5   0.25 MB  Per-chromosome copy counts…

  Individuals
    ind_meta                  60   6   0.25 MB  Individual-level metadata…
    ind_haplotype        240,000   7   0.75 MB  Phased haplotypes in long…
  …
  ─────────────────────────────────────────────────────────────────────
  Size is on-disk storage rounded up to whole 256 KB blocks, so any
  non-empty table reads as at least 0.25 MB and small tables are not
  distinguishable. Per-table sizes sum to less than the 4.0 MiB file —
  catalog and header blocks are not attributed to any table.
```

Requirements on that text:

- **State both caveats, not one.** The floor explains why small tables look
  identical; the shortfall explains why the column does not add up to the header.
  A reader who sees only one will treat the other as a bug.
- **Name the block size numerically** (`256 KB`, `0.25 MB`) rather than saying
  "approximate". A reader can then reason about which rows are real signal —
  anything at the floor is "small, unknown", anything above it is a true reading.
- **Reference the header total** so the shortfall is self-evidently expected
  rather than looking like missing data.
- Print it once as a block footer, not per row.
- Mirror it in `?schema` under a `Size reporting` heading, in the same words, so
  it is discoverable without running the function.

If `sizes = TRUE` is passed for an in-memory population, there are no blocks at
all (`block_size` is `0`). **Decided (2026-08-31): warn and omit the column, do
not error.** `schema()` is a display function and it still has something honest to
show in that case — 5a already prints `memory_usage` in the header for in-memory
populations — so refusing to print anything is a worse answer than printing
everything except a Size column. Never print a column of zeros under a footnote
that does not apply.

The harder line stays where it belongs: `order = "size"` without `sizes = TRUE`
**errors** (Part 6). The distinction is that there the requested ordering is
undefined, not merely unavailable — there is no honest fallback ordering to
silently substitute.

---

## Part 6 — `schema(pop, order = ...)`

```r
schema(pop, order = c("pipeline", "name", "rows", "size"), show_empty = FALSE,
       include_system = FALSE, sizes = FALSE)
```

- `"pipeline"` (default) — grouped as above, with headers.
- `"name"` — flat alphabetical, no headers. Predictable lookup.
- `"rows"` — flat, descending row count. Answers "what is actually big here?",
  which is a genuine question for file-based DuckDB populations.
- `"size"` — flat, descending on-disk bytes. Requires `sizes = TRUE` (Part 5b);
  error rather than silently falling back if it is not set.

Headers print only for `"pipeline"`. `table_group` stays in the returned tibble
for all four.

---

## Part 7 — The rename: `trait_effects` → `phenotype_effects`

### Verification — yes, rename them

Both tables are **already keyed by `phenotype_name`**, including in their primary
keys (`R/define_trait.R:182-208`):

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
- **`R/define_trait.R:342-416`** — the whole `── Migrations for databases created
  before v0.31.0 ──` section inside `ensure_trait_tables()`, of which the
  `trait_name` → `phenotype_name` column rename is one block. Per the pre-1.0.0
  policy ("no compatibility shims, no legacy paths"), it is deleted rather than
  extended for the table rename. The section is larger than this note originally
  implied and it is not the only such code — see **Legacy migration code** below.

### Newly found while auditing — DECIDED: reserve all columns

The registry audit surfaced three tables absent from `TABLE_RESERVED_COLS`
entirely: `founder_haplotypes`, `phenotype_components`, and
`trait_random_effects`. When a table has no entry, `mutate_table()` blocks
nothing, so a user can overwrite e.g. `trait_random_effects.draw_value` or
`phenotype_components.contributor_type` in place.

**Decided (2026-08-31): reserve every column on all three, and ship it as its own
commit before the rename.**

None of the three is a table a user annotates. `founder_haplotypes` is a sampling
pool consumed by `add_founders()`; `trait_random_effects` holds sampled draws
written by `define_effect_random()`; `phenotype_components` is configuration
generated from the `components =` data frame passed to `define_phenotype()`. The
tables that legitimately take user columns are the entity-shaped ones — `ind_meta`,
`genome_meta`, `ind_phenotype`, `index_meta` — which is also exactly where the
`...` forwarding in the `add_*` functions writes. `CLAUDE.md` already documents
"**Reserved**: all columns" for every table in this class (`genome_effects`,
`phenotype_meta`, `ind_true_index`); these three are omissions, not deliberate
openings.

Sequenced first and separately because it is a bug fix of the same kind as
v0.63.1, with no dependency on either the rename or the printing work. Add the
missing "**Reserved**" lines for `founder_haplotypes` and `phenotype_components`
to `CLAUDE.md` in the same change.

### Legacy migration code — delete it, do not extend it

**Decided (2026-08-31): the rename gets no `restore_pop()` migration, and the
migration code that already exists is deleted in the same commit.**

There is more of it than the bullet above implied, and there is a live precedent
pointing the other way that has to be resolved deliberately rather than by
omission:

| Location | What it does |
|---|---|
| `R/define_trait.R:342-416` | Eight guarded blocks inside `ensure_trait_tables()` under the `── Migrations for databases created before v0.31.0 ──` banner: drop `ind_tbv.date_calc`; drop `ind_ebv.date_calc` and add `eval_number`; add `index_meta.economic_weight`; the `trait_effects`, `trait_random_effects`, and `ind_phenotype` `trait_name` → `phenotype_name` column renames; add `phenotype_meta.missing_component_action`; add `phenotype_meta.formula_tbv` / `formula` |
| `.migrate_var_comp_tables()` (`R/open_pop.R:368`), called from `R/restore_pop.R:78` | `ALTER TABLE … RENAME` for two earlier **table** renames — `trait_effect_cov` → `trait_var_comp` and `phenotype_residual_cov` → `phenotype_var_comp` |

The second row is the precedent that matters: it is precisely the migration this
rename would be asked to imitate, for precisely the same kind of change. Keeping
it while deleting the column-rename block in `define_trait.R` would be incoherent
— same class of code, written for the same reason, both prohibited by the same
pre-1.0.0 policy. Both go, in the rename commit.

Consequence, and it must be **stated in `NEWS.md` rather than discovered**: a
`.duckdb` file written before this release will not restore usefully.
`restore_pop()` reads its table list from the DuckDB catalog (`R/restore_pop.R:80`),
so an old file's `trait_effects` would surface as an unrecognized **User table**
under the new grouping while every phenotype reader and writer looked for
`phenotype_effects` and failed. Pre-1.0.0 the answer is to rebuild the population,
not to migrate it.

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

### Re-split the description helpers — DECIDED: do it in the rename commit

`.trait_layer_descriptions()` in `R/schema.R:344` currently describes
`phenotype_meta`-adjacent tables *and* `ind_ebv` / `index_meta` / `ind_index`,
while `.core_layer_descriptions()` (`R/schema.R:175`) describes
`phenotype_meta`, `phenotype_components`, and `phenotype_var_comp`. The three
description functions do not partition along the same lines as the display groups
in Part 2. Re-splitting them to match (`.genome_*`, `.individual_*`,
`.genetic_model_*`, `.observation_model_*`, `.results_*`) would make the group
vector and the descriptions the same list read twice, which is harder to let drift.

**Decided (2026-08-31): fold this into the rename commit rather than leaving it
optional.** The rename already rewrites ~22 references inside
`.trait_layer_descriptions()` (`R/schema.R:344`), so that function is being edited
line by line regardless; the re-split costs little on top of that and is far more
expensive as a second pass over the same lines later. Doing it there also means
the Part 3 group vector and the description helpers land in the same shape at the
same time, instead of the group vector being written against a partitioning the
descriptions do not share.

---

## Suggested sequencing

1. ~~Fix `R/sql_utils.R:98` and `:155`~~ — done in v0.63.1.
2. ~~**Reserve all columns on `founder_haplotypes`, `phenotype_components`, and
   `trait_random_effects`**, plus the matching `CLAUDE.md` lines~~ — **DONE.**
   Three entries added to `TABLE_RESERVED_COLS`; the roxygen note above the list
   now states the "managed table reserves everything" rule so the next table
   added has to make the call deliberately. Two new tests in
   `test-schema-registries.R` (coverage + `mutate_table()` blocking).
   `CLAUDE.md` gained **Reserved** lines for `founder_haplotypes` and
   `phenotype_components`; `trait_random_effects` has no `CLAUDE.md` section at
   all, so one is written in step 4 under its final name.
3. ~~**Add the database-size header** (Part 5a)~~ — **DONE.** New internal
   `.schema_size_label()` in `R/schema.R`, attached to the `schema()` result as a
   `db_size` attribute and rendered on the right of the header rule. Verified
   against a live DuckDB: file-backed prints `2.0 MiB on disk`, in-memory prints
   `2.1 MiB in memory`, and an uncheckpointed database prints
   `0 bytes + 37.6 KiB WAL on disk` — which is why the WAL term is not optional
   (see the note added to Part 5a). `schema()` issues no `CHECKPOINT`; there is a
   test asserting the WAL is still non-empty after printing. Five tests in the
   new `tests/testthat/test-schema-print.R`.
4. ~~**Rename the two tables**~~ — **DONE.** (Details in the completion notes.)
5. ~~**Implement Part 2 + Part 3 Option A + Part 4**~~ — **DONE.** (Details in
   the completion notes at the end.)
6. ~~**Add the `order` / `show_empty` / `include_system` / `sizes` arguments**~~
   — **DONE.** (Details in the completion notes.)

Steps 2 and 3 are independently shippable in any order. Step 4 must precede step 5
so the group vector is written once. Version bump + `NEWS.md` entry per `CLAUDE.md`
before each commit — and step 4's entry must carry the "pre-existing `.duckdb`
files must be rebuilt" statement from **Legacy migration code**.

**Executed 2026-08-31 as a single uncommitted change set at v0.64.0**, not as five
commits: the work was done in plan order but the user asked to hold the commit
until they had reviewed it. The `NEWS.md` entry covers all five steps under one
version heading and carries the rebuild warning.


---

## Completion notes (2026-08-31)

Implemented at **v0.64.0**. Full suite: **0 failures**, 1 skip (`On CRAN`), 10
warnings — all ten pre-existing and unrelated (monomorphic-locus and
pooled-base-frequency advisories, the `formula_tbv` scalar-constant advisory).

`R CMD check --no-manual --no-vignettes --no-tests --no-examples`: 3 WARNINGs,
4 NOTEs, none attributable to this change. The non-ASCII WARNING lists 14 files;
`R/schema.R` was already among them at `HEAD` (14 non-ASCII lines before this
work), and every non-ASCII line this change *added* to it is a comment, which the
check explicitly permits — all new executable strings use `\uxxxx` escapes. The
"R code for possible problems" NOTE names no new function. The remaining findings
(vignettes not in `inst/doc`, invalid License DCF stub, hidden files) are
pre-existing and unrelated to R code.

### What landed, by step

**Step 2 — reserved columns.** `founder_haplotypes`, `phenotype_components`, and
`trait_random_effects` (renamed in step 4) added to `TABLE_RESERVED_COLS` with
every column listed. The roxygen block above the list now states the rule —
package-managed tables reserve everything, entity-shaped tables leave room — so
the next table added has to make the call deliberately. Two new tests in
`test-schema-registries.R`. `CLAUDE.md` gained the matching **Reserved** lines.

**Step 3 — database-size header.** `.schema_size_label()` reads
`PRAGMA database_size` and is attached to the `schema()` result as a `db_size`
attribute. Verified live: `2.0 MiB on disk`, `2.1 MiB in memory`, and
`0 bytes + 37.6 KiB WAL on disk` before a checkpoint.

**Step 4 — rename and deletions.** `trait_effects` → `phenotype_effects`,
`trait_random_effects` → `phenotype_random_effects` across 11 R files, 6 test
files, one vignette script, `README.md`, `CLAUDE.md`, and regenerated `man/`.
Deleted: the whole `R/define_trait.R` pre-v0.31.0 migration section (76 lines,
eight guarded blocks) and `.migrate_var_comp_tables()` plus its `restore_pop()`
call (44 lines). Description helpers re-split (below). `NEWS.md` carries the
rebuild warning.

**Step 5 — grouping.** `.schema_table_order()`, `.schema_group_of()`,
`.schema_format_rows()`; `schema()` gained `show_empty` / `include_system` and a
`table_group` factor column; the print method gained group headings, empty
collapse, and abbreviated counts.

**Step 6 — order and sizes.** `.schema_table_bytes()`, `.schema_format_size()`,
`.schema_print_size_footnote()`; `schema()` gained `order` and `sizes`.

### Verified against the plan's own claims

- **Part 1's defect is fixed and now has a test.** Built-then-restored orderings
  are `identical()`; before the change they differed. `test-schema-print.R`.
- **Part 2's coverage check holds in code, not just on paper.** A test asserts
  `.schema_table_order()`, `SYSTEM_TABLES`, and `.all_schema_descriptions()` name
  the same 24 tables, so the three lists cannot drift apart silently.
- **Part 5b's measurements reproduced exactly** on a 60-founder / 2,000-locus
  population: 4.0 MiB, 16 used blocks, `ind_haplotype` 3 blocks,
  `founder_haplotypes` 2, `genome_meta` 1, `chr_inheritance` 1, 10 of 16 blocks
  attributed. The `duckdb_tables().estimated_size` trap is real; `storage_info()`
  is the route used.
- **The `CHECKPOINT` caveat is enforced by test**, in both directions: the WAL is
  still non-empty after `schema()` prints, and empty after `schema(sizes = TRUE)`.

### Deviations from the plan, and why

1. **The 5d footnote text changed, because one of its claims is false.** The plan
   asserts "any non-empty table reads as at least 0.25 MB". Measured: a table
   small enough to live in the catalog reports **zero** blocks while holding rows
   (a 1-row table after `CHECKPOINT`). The footnote now says most non-empty
   tables read as 0.25 MiB and a catalog-resident one reads as 0.00 MiB. Both
   required caveats are still stated, the block size is still named numerically,
   and it still references the header total.
2. **`MiB`, not `MB`.** The header takes DuckDB's own formatting, which is `MiB`.
   Printing `0.25 MB` beside `4.0 MiB` in the same output would invite exactly
   the "these do not add up" confusion the footnote exists to prevent. One block
   is 262,144 bytes = 0.25 MiB exactly, so `MiB` is also the correct unit.
3. **`?schema`'s `Size reporting` section is expanded rather than verbatim.** It
   states the same caveats as the printed footnote plus the `CHECKPOINT` one,
   which belongs in the docs but not in a per-print footer.
4. **In-memory `order = "size"` errors even with `sizes = TRUE`.** The plan
   settled warn-and-omit for `sizes = TRUE` in memory, and error for
   `order = "size"` without `sizes`. Their intersection was unspecified; ordering
   by a column that cannot exist is undefined, so it takes the error branch.
5. **Descriptions are registered once by `open_pop()`.** The plan asked for the
   re-split but not for a call-site change. The two are inseparable: the old
   three functions were partitioned by *creation time* (`define_genome`,
   `open_pop`, `define_trait`) and the display groups cut across that — the
   Individuals group spans two creation times. Registering the whole set up front
   is what lets the functions be 1:1 with display groups. It also removes the
   window in which a later `define_*()` call clobbered a user's
   `define_schema_description()` override, and folds in the fourth registration
   site (`founder_haplotype_helpers.R`).
6. **Two tests changed, both testing deleted behavior.**
   `test-formula_phenotype.R`'s "old DB ... gets migrated by
   ensure_trait_tables" was deleted with the migration it exercised.
   `test-define_genome.R`'s rollback test dropped its `_schema_meta` assertion:
   descriptions are no longer written inside `define_genome()`'s transaction, so
   surviving its rollback is now correct. Its table and temp-view assertions are
   unchanged and still pass.
7. **Scope added beyond the plan, all small:** a `CLAUDE.md` section for
   `phenotype_random_effects` (it had none), the `phenotype_effects` section moved
   out of the genetic-layer block into the observation layer to match the new
   grouping, a `schema()` / `describe_table()` section in `CLAUDE.md` recording
   the two-list maintenance obligation, and removal of two now-dead locals
   (`is_new_table`, `avail`).

### Deliberately not done

All three are written up with suggested fixes in
**`plans/check_after_schema_update.md`**.

- **`TABLE_ROW_KEYS` / `TABLE_PRIMARY_KEYS` are still missing entries** for
  `phenotype_components` (which has an `id_phenotype_comp` integer PK),
  `founder_haplotypes`, and `phenotype_random_effects`. `remove_rows()` therefore
  cannot target rows in them. This is the same *class* of gap step 2 closed for
  `TABLE_RESERVED_COLS`, but it is a different registry with different
  consequences, and step 2's scope was explicitly the reserved columns. Worth its
  own decision.
- **`archive_replicate()`'s table lists were not revisited.** `store_once` names
  `phenotype_effects` but not `phenotype_random_effects`, `genome_map`,
  `chr_inheritance`, `chr_recombination`, or `founder_haplotypes`. Pre-existing,
  unrelated to this plan, and a real question for replicate archiving.
- **`phenotype_components` has duplicate `CREATE TABLE` DDL** in both
  `R/open_pop.R` and `R/define_trait.R`. Harmless today (both guarded by an
  existence check, and the definitions agree) but it is two places to keep in
  sync. Not touched.
