# Refactor: Wide → Long Haplotype & Genotype Tables

**Status**: Draft  
**Motivation**: Wide format blocks three critical capabilities:
1. Line-specific QTL effects in crossbreeding — you cannot join each allele to its line of origin
2. Novel mutations — new loci require `ALTER TABLE ADD COLUMN`, which does not scale and breaks existing queries
3. Correct heterozygous crossbred TBV — requires knowing which allele came from which genetic line

---

## Response to critique (`plans/refactor_haplotype_critique.md`)

This plan was reviewed and the review raised eight numbered issues plus a
staged-rollout recommendation. Decisions below; the plan body has been
edited to reflect each one.

| # | Critique point | Decision |
|---|-----------------|----------|
| Scope | Split into staged rollout instead of one refactor | **Agree.** Adopted as "Staged Implementation" section, replacing the old flat file-change table as the source of truth for sequencing. |
| 1 | `parent_origin`/`strand` semantics under-specified for meiosis | **Agree — real bug.** Fixed the `add_offspring()` description below: gamete generation must load *all* of a parent's own haplotype rows for a chromosome (both of that parent's own `parent_origin` values), not just one. |
| 2 | Dense vs. sparse haplotype semantics contradict each other | **Agree.** v1 is fully dense, no exceptions. Sparse mutation storage moved out of v1 scope entirely (Stage 5, see staged plan). |
| 3 | `line_origin` FK to `genome_effects.line_name` is wrong | **Agree.** Removed the FK; documented as unconstrained/convention-based. |
| 4 | Crossbreeding TBV SQL omits centering | **Agree.** Added explicit centered SQL for population-wide, line-specific-with-fallback, and imprinted cases. |
| 5 | Founder haplotypes lose `line_name` | **Agree — this is a real regression.** The current wide `founder_haplotypes` table has `line_name` (`R/define_founder_haplotypes.R`); the long schema proposed below had dropped it. Added back. |
| 6 | Table renaming (`genome_*` → `ind_*`) needs a compatibility decision | **Decided: rename to `ind_haplotype`/`ind_genotype`.** Any table keyed by `id_ind` uses the `ind_` prefix (`ind_meta`, `ind_phenotype`, `ind_tbv`, `ind_ebv`, `ind_index`); `genome_*` is reserved for non-individual data (`genome_meta`, `genome_effects`). No compatibility views — rename and break, per the "Table naming" note in Open Questions. |
| 7 | Sex chromosome design isn't complete enough to code safely | **Agree.** Marked out of v1 scope (Stage 4). Design kept as a draft for later, with known gaps flagged rather than resolved now. |
| 8 | Performance estimates need benchmarking before the schema is locked | **Agree.** Added a benchmark step as a Stage 1 exit criterion; numbers in this doc are relabeled as unverified estimates. |
| — | Rcpp section shouldn't be an acceptance criterion for v1 | **Agree.** Moved to an appendix marked "future optimization, out of scope for v1." |
| — | Migration via `pivot_longer()` won't scale | **Agree.** Replaced with a transactional SQL/`UNPIVOT` procedure. |
| — | `add_dosage()` vs. `add_genotypes()` naming confusion | **Agree the confusion is real, disagree on renaming to `materialize_dosage()`.** `add_dosage()` fits the project's own `add_*` convention (`add_tbv`/`add_ebv` also write derived, non-raw-simulation values under `add_`), so I kept the name and added an explicit doc contrast instead. |

---

## Non-goals for v1

This refactor bundles at least six related but separable design changes across five stages (see
"Staged Implementation" below). **"v1" means Stage 1 only** — the minimal long-format storage
rewrite with current behavior preserved exactly. Everything else ships as a later, separately
tested stage in the same roadmap — not as a "someday" idea — sequenced in this order:

1. **Wide-to-long storage** — v1 (Stage 1).
2. **Table renaming (`genome_*` → `ind_*`)** — v1 (Stage 1; same DDL/rewrite work either way):
   rename to `ind_haplotype`/`ind_genotype`.
3. **Populating `line_origin`** — v1 (Stage 1). `add_founders()`/`add_offspring()` set it from the
   start, even though nothing reads it yet — see the note on why in Stage 1 scope. This is
   plumbing, not the crossbreeding TBV feature itself.
4. **Out of v1, but next**: making `add_tbv()` actually *use* `line_origin` for crossbreeding TBV
   (Stage 2) — this is one of the primary motivations for doing this refactor at all (see the top
   of this document), just sequenced right after Stage 1 lands and passes its parity tests, not
   deferred indefinitely like items 6–7 below.
5. **Out of v1**: on-demand dosage caching (`add_dosage()`) (Stage 3) — needed for
   `extract_genotypes()`/BLUPF90 export, but that rewrite reads directly from `ind_haplotype` in
   Stage 1 regardless (see "Genotype Matrix Export for GBLUP"), so the cache itself can wait.
6. **Out of v1**: sex chromosome and polyploid inheritance rules (Stage 4).
7. **Out of v1**: sparse mutation storage (Stage 5) and Rcpp acceleration (appendix, future).

## Staged Implementation

This is the sequencing authority for the refactor — it supersedes any flat, unordered reading of
the "Files Requiring Changes" table further down. That table still exists as a per-file index of
what changes and why; use it to look up a given file's changes, but use this section to decide
*when* each change happens. Each stage should land, pass its exit criteria, and be usable on its
own before the next stage starts.

### Stage 1: Long diploid core, no new biology

Goal: preserve current behavior exactly while changing storage shape.

Scope:
- Diploid autosomes only; `chr_meta` created with default autosome rows, `define_chr()` not yet
  implemented.
- Keep current founder generation methods; `founder_haplotypes` gains `line_name` (see schema
  section) and stays otherwise equivalent.
- `add_founders()` and `add_offspring()` behavior equivalent to current tests, corrected per the
  "Meiosis semantics" section below.
- Fully dense `ind_haplotype`; no sparse mutation semantics.
- Rename `genome_haplotype`/`genome_genotype` → `ind_haplotype`/`ind_genotype` (decided — see
  Open Questions); update all ~13 referencing files and the migration procedure accordingly.
- Rewrite `add_tbv()` (with centering), `extract_genotypes()`, and BLUPF90 export against the long
  source of truth. `extract_genotypes()` also gains the new `loci_tbl` argument (see "Genotype
  Matrix Export for GBLUP" below) as a general locus filter independent of chips/QTL sets.
- `genome_meta.locus_name` gets a `UNIQUE` constraint or explicit validation.
- **Populate `line_origin` from the start** (unconstrained, no FK — see schema section):
  `add_founders()` sets it from the founder's `ind_meta.line_name`; `add_offspring()`'s
  recombination returns both an allele vector and a `line_origin` vector, inherited from whichever
  parental segment contributed each locus. This costs nothing extra — recombination already
  tracks which parental segment contributed each locus — and deferring it would be
  unrecoverable later: the crossover provenance isn't stored anywhere else, so any simulation run
  before this landed would have permanently lost `line_origin` data. TBV does not yet *use*
  `line_origin` in Stage 1 (see Stage 2); it's populated now purely so no data is lost.

Exit criteria:
- Current tests pass after being updated for the new table shape.
- A small deterministic simulation produces identical haplotypes, genotypes, TBVs, and exports
  before and after the refactor, same seed.
- `line_origin` is correctly populated for founders and offspring (verified against expected line
  ancestry in the deterministic test simulation above).
- The benchmark script (see "Storage, Memory, and Performance") exists and validates the long
  format's performance against wide at realistic scale.
- Migration procedure (see "Migration for Existing Databases") tested against a real pre-refactor
  database.

### Stage 2: Line-origin TBV

Goal: enable crossbreeding additive TBV by allele line of origin. `line_origin` itself is already
populated as of Stage 1 (see above) — this stage is entirely about `add_tbv()` learning to *use*
it.

Scope:
- Add the centered SQL for population-wide and line-specific additive effects (already drafted
  below).
- Tests for: F1, F2, and backcross animals; `line_origin IS NULL` fallback; a line with
  line-specific effects at some loci but not others; both a line-specific and population-wide row
  existing for the same locus/line; `base_allele_freq` differing by line.

Exit criteria:
- A hand-calculated F1 TBV test passes.
- An F2 recombination test proves `line_origin` follows the inherited segment (this can reuse the
  Stage 1 population data directly, since it's already populated).
- All edge-case tests above pass.

### Stage 3: Dosage cache and export hardening

Goal: make marker-assisted selection and external software export efficient without making cache
population a hidden requirement.

Scope:
- Add `add_dosage()` (see dedicated section below, including the `add_genotypes()` contrast and
  `overwrite_dosage` argument).
- Treat `ind_genotype` as a documented partial cache.
- Ensure `extract_genotypes()` and BLUPF90 export meet the acceptance criteria listed under
  "Genotype Matrix Export for GBLUP."

Exit criteria:
- Cache and direct extraction produce identical matrices.
- Repeated cache materialization is idempotent.
- Partial cache contents cannot silently change TBV or phenotype results.

### Stage 4: Non-diploid inheritance

Goal: add `chr_meta` rules beyond the default, `define_chr()`, sex chromosomes, organelles, and
later polyploid rules.

Scope:
- Resolve the `copies_M`/`copies_F` terminology clash (0/1/2 vs. 4/6/8) before writing any code.
- Start with sex chromosome and organellar cases.
- Defer autopolyploid/allopolyploid pairing until a clear meiosis model exists (see Open
  Questions).
- Make TBV centering and dosage export ploidy-aware.

Exit criteria:
- X/Y or Z/W hand tests pass.
- Absent chromosomes produce no haplotype rows and no accidental dosage.
- Export either handles sex-linked loci correctly or warns/errors when used for standard GBLUP.

### Stage 5: Mutation sparsity

Goal: support introduced loci without backfilling every existing individual.

Scope:
- Decide fully sparse vs. hybrid dense/sparse for introduced loci (not decided by this plan — see
  "Novel Mutation Support" below).
- Add query helpers that always coalesce missing alleles to zero, applied consistently across
  dosage, TBV, allele frequency, export, and recombination parent loading.

Exit criteria:
- A mutation in one individual can be inherited by descendants.
- Non-carriers export as dosage zero.
- TBV centering remains correct for carriers and non-carriers.

---

## Current Schema (Wide)

### `genome_haplotype` (wide — current name, being replaced)

```text
id_ind         VARCHAR    FK to ind_meta
parent_origin  INTEGER    1 = paternal (from id_parent_1), 2 = maternal (from id_parent_2)
locus_1        UTINYINT   0 or 1
locus_2        UTINYINT
...
locus_n        UTINYINT
```
Two rows per individual. Column count grows with `n_loci`.

### `genome_genotype` (wide — current name, being replaced)

```text
id_ind         VARCHAR    FK to ind_meta
locus_1        UTINYINT   0, 1, or 2
...
locus_n        UTINYINT
```
One row per individual. Derived as `locus_i_hap1 + locus_i_hap2`.

### Why wide fails

| Problem | Root cause |
|---------|-----------|
| Novel mutations | Adding locus requires `ALTER TABLE ADD COLUMN` across all rows |
| Line-specific QTL | No allele → line mapping; TBV JOIN has no `line_origin` to filter on |
| Crossbred TBV | Heterozygous crossbred individual has two alleles from different lines; wide format loses this |
| Panel subsets | Query must name specific columns (`locus_1, locus_42, ...`) rather than `WHERE locus_name IN (...)` |
| GWAS by locus | Requires expensive `UNPIVOT` before aggregation |

---

## Proposed Schema (Long)

### `ind_haplotype` (long)

One row per (individual × haplotype × locus). 

```text
id_ind         VARCHAR    FK to ind_meta; part of composite PK
parent_origin  UTINYINT   1 = from id_parent_1 (sire/male parent), 2 = from id_parent_2 (dam/female parent);
                          always 1 or 2 regardless of ploidy; part of composite PK
strand         UTINYINT   Which copy from that parent: always 1 for diploids; 1..copies/2 for polyploids
                          (e.g. 1 and 2 from each parent in a tetraploid); part of composite PK
line_origin    VARCHAR    Founding line this allele traces back to. Unconstrained (no FK) —
                          see "line_origin semantics" below. NULL = untracked/unknown origin.
locus_id       INTEGER    FK to genome_meta.locus_id; part of composite PK. Physical sort/PK key —
                          integer comparisons are cheaper than string comparisons for the hot
                          recombination and compression paths.
locus_name     VARCHAR    FK to genome_meta.locus_name; denormalized alongside locus_id (not part
                          of the PK) so TBV/GWAS queries can join straight to genome_effects
                          without routing through genome_meta.
allele         UTINYINT   0 or 1 (phased; one strand of the chromosome)
```

**Primary key**: `(id_ind, parent_origin, strand, locus_id)`, declared as an actual DuckDB
`PRIMARY KEY` (or `UNIQUE`) constraint in the DDL — not just a documented convention. This is
required for `INSERT OR REPLACE` semantics elsewhere (e.g. `ind_genotype`) to be meaningful.

**Indexes**:
- `(locus_name)` — for GWAS, per-locus lookups, TBV over QTL set (direct join to `genome_effects`)
- `(id_ind, parent_origin, strand)` — for loading a parent's stored chromosome copies during
  offspring generation

**`genome_meta.locus_name` uniqueness**: this key design assumes `locus_name` is unique in
`genome_meta`. That is not currently enforced by `define_genome()`/`initialize_genome()`. Add a
`UNIQUE` constraint (or an explicit validation check when custom `locus_names` are supplied) as
part of this refactor — a duplicated `locus_name` would silently corrupt joins throughout the
long schema.

**Both `locus_id` and `locus_name` are stored — this replaces the earlier "pick one" framing.**
Storing only `locus_name` meant every hot-table join, sort, and compression pass paid string-key
cost; storing only `locus_id` would have meant adding a redundant `locus_id` column to
`genome_effects` (or an extra join hop through `genome_meta`) just to compute TBV. Storing both
avoids the tradeoff entirely: `locus_id` is the PK/physical sort key (cheap integer comparisons
for the recombination hot path and RLE/zone-map compression), `locus_name` is a plain denormalized
column used only for the direct `genome_effects` join and user-facing filters — it is never part
of the primary key, so it adds one dictionary-encoded column's worth of storage (cheap, per the
compression table below) without touching the performance-critical key comparisons. The general
long-vs-wide and sort-order benchmark in "Staged Implementation" is still worth running, but it no
longer needs to adjudicate this specific fork.

**One consistency caveat**: because `locus_name` is now denormalized into every haplotype row,
`locus_name` must be treated as immutable once any `ind_haplotype` rows reference it — renaming a
locus in `genome_meta` after data exists would require cascading the update across every
referencing row (up to 500M+), which nothing in this plan supports and should not be attempted
casually. If locus renaming is ever needed, add an explicit function for it rather than a raw
`UPDATE genome_meta`.

**`parent_origin` and `strand`**: `parent_origin` is always 1 or 2 regardless of ploidy — for the
individual whose rows are being read, 1 means the strand(s) came from that individual's
`id_parent_1` (sire/male parent), and 2 means they came from that individual's `id_parent_2`
(dam/female parent). This preserves the biologically meaningful sire/dam distinction at any
ploidy level. For polyploids a single parent contributes more than one chromosome copy; `strand`
identifies which copy within that parent's contribution. For diploids `strand = 1` always. For a
tetraploid (4N) each parent contributes 2 copies, so `strand ∈ {1, 2}` per `parent_origin`,
giving 4 rows per individual per locus. For a hexaploid (6N) `strand ∈ {1, 2, 3}` per parent (6
rows total); for an octoploid (8N) `strand ∈ {1, 2, 3, 4}` (8 rows total). The number of strands
per parent is `copies / 2` where `copies` is `chr_meta.copies_M` or `chr_meta.copies_F` for the
individual's sex. `UTINYINT` (range 0–255) handles any realistic ploidy for both `parent_origin`
and `strand` without a type change. Imprinting queries filter on `parent_origin` directly — no
lookup into `chr_meta` needed at query time.

**`line_origin` semantics**:
- Founders: set from `ind_meta.line_name` at `add_founders()` time
- Offspring: inherited locus-by-locus from the contributing parental haplotype during recombination in `add_offspring()`
- NULL is valid only when line origin is genuinely untracked or unknown, not merely because the
  user has not defined line-specific effects

**v1 storage convention: fully dense, no exceptions.**  
Every individual has a row for every locus in `ind_haplotype`, including any locus added after
`initialize_genome()`. There is no sparse case in v1 — `add_mutation()` is **not** part of this
refactor (see "Staged Implementation" / Stage 5 below). An earlier draft of this plan proposed
sparse storage for novel mutations while keeping founding loci dense; that mixed convention is
rejected for v1 because it requires every downstream query (dosage, TBV, allele frequency,
`PIVOT` export, parent-haplotype loading for recombination) to consistently treat absent rows as
implicit zero, and none of the SQL in this plan does that today. If/when mutation support is
built (Stage 5), it gets its own design pass — including whether it's fully sparse, hybrid, or
solved a different way (e.g. backfilling zero rows in a background job) — rather than being
retrofitted onto this schema after the fact.

---

### `ind_genotype` (long)

One row per (individual × locus). Populated **on demand** by `add_dosage()` — never auto-populated by `add_founders()` or `add_offspring()`. Starts empty. Users call `add_dosage()` explicitly when they need dosage values for marker-assisted selection, allele frequency queries, or other downstream analysis.

```text
id_ind        VARCHAR    FK to ind_meta; part of composite PK
locus_id      INTEGER    FK to genome_meta.locus_id; part of composite PK; also drives stable
                         column ordering in the PIVOT export (see "Genotype Matrix Export")
locus_name    VARCHAR    FK to genome_meta.locus_name; denormalized alongside locus_id, same
                         rationale as ind_haplotype above
dosage_value  UTINYINT   Sum of allele across all haplotype strands. Range 0–2 for diploids,
                         0–4 for tetraploids, 0–8 for octoploids. UTINYINT (0–255) covers
                         any realistic ploidy without type change.
```

**Primary key**: `(id_ind, locus_id)`  
**Indexes**:
- `(locus_name)` — allele frequency computation, per-locus MAS queries
- `(id_ind)` — full individual dosage fetch

Derived from `ind_haplotype` on demand:
```sql
SELECT id_ind, locus_id, locus_name, SUM(allele) AS dosage_value
FROM ind_haplotype
WHERE locus_name IN (...)   -- filtered to chip or QTL loci
  AND id_ind IN (...)       -- filtered to candidate individuals
GROUP BY id_ind, locus_id, locus_name
```

INSERT uses `INSERT OR REPLACE` semantics — calling `add_dosage()` twice for the same individual × locus is idempotent (dosage cannot change once haplotypes are set).

---

### `add_dosage()` — new function

Follows the same pattern as `add_tbv()` and `add_phenotype()`: accepts a filtered `tidybreed_table` as its first argument, computes a SQL aggregation, and inserts results into a table. Returns `tidybreed_pop` invisibly.

```r
add_dosage(tbl, chip_name = NULL, locus_names = NULL, overwrite_dosage = FALSE)
```

- `tbl` — a `tidybreed_table` from `get_table()`, optionally filtered. Any table with an `id_ind` column is accepted. Consistent with the existing `add_tbv()`/`add_phenotype()` convention (see `get_table()` docs), the **distinct** `id_ind` values in the collected, filtered table form the candidate set — this matters explicitly for a table like `ind_phenotype` that can have multiple rows per individual; `add_dosage()` does not multiply work per extra row.
- `chip_name` — name of a chip previously defined via `define_chip()` (e.g. `"50K"`). `define_chip()` writes a `BOOLEAN` column named `is_{chip_name}` to `genome_meta`; `add_dosage()` filters loci using `WHERE genome_meta.is_{chip_name} = TRUE`. The `chip_name` value must exactly match a chip that already exists in `genome_meta` — if the column `is_{chip_name}` is not found, `add_dosage()` throws an error. If `NULL` and `locus_names` is also `NULL`, dosage is computed for all loci present in `ind_haplotype` for the target individuals.
- `locus_names` — explicit character vector of locus names. Takes precedence over `chip_name` if both are supplied.

**`add_dosage()` vs. `add_genotypes()` — do not confuse these.** `add_genotypes()` already
exists today and means "mark animals as physically genotyped on a chip" — it writes a `BOOLEAN`
`has_<chip_name>` column to `ind_meta` and touches no dosage data at all. `add_dosage()` is a
different operation: it **materializes simulated dosage values** (0/1/2, ground truth from
`ind_haplotype`) into `ind_genotype`. The names are close enough that user-facing docs must
contrast them explicitly (a short "see also" cross-reference in each function's roxygen block).
`add_dosage()` keeps the `add_` prefix rather than becoming `materialize_dosage()`: per this
project's naming convention (`add_*` = writes simulation *output* rows), `add_tbv()` and
`add_ebv()` already use `add_` for computed/derived individual-level values, not just raw
simulation draws, so `add_dosage()` is consistent with that existing precedent.

`INSERT OR REPLACE` requires the declared `PRIMARY KEY (id_ind, locus_id)` constraint on
`ind_genotype` from the schema section above — without it, idempotence here is only aspirational,
not enforced by DuckDB.

**Cache scope and `overwrite`**: no `overwrite` argument is needed for *value* correctness —
dosage is fully determined by the haplotypes already written, so calling `add_dosage()` twice for
the same individual × locus simply replaces the row with the same value via `INSERT OR REPLACE`.
But cache *scope* still matters across repeated calls over a simulation's lifetime (see Open
Questions): resolve that with an explicit `overwrite_dosage = FALSE` argument, named to match the
existing `overwrite_index` parameter on `add_index()`, that clears prior `ind_genotype` rows for
the candidate individual set before inserting when set to `TRUE`.

**Internal SQL** (diploid example):
```sql
INSERT OR REPLACE INTO ind_genotype (id_ind, locus_id, locus_name, dosage_value)
SELECT h.id_ind, h.locus_id, h.locus_name, CAST(SUM(h.allele) AS UTINYINT) AS dosage_value
FROM ind_haplotype h
WHERE h.id_ind IN (/* ids from tbl */)
  AND h.locus_name IN (/* chip or locus_names filter */)
GROUP BY h.id_ind, h.locus_id, h.locus_name
```

**Typical workflows:**

```r
# Compute dosage for all generation-5 selection candidates on the 50K chip
pop <- pop |>
  get_table("ind_meta") |>
  filter(gen == 5L) |>
  add_dosage(chip_name = "50K")

# Compute dosage at specific QTL only
pop <- pop |>
  get_table("ind_meta") |>
  filter(gen == 5L) |>
  add_dosage(locus_names = c("QTL_1", "QTL_42", "QTL_107"))

# Re-materialize dosage for this individual set, clearing any prior ind_genotype rows for them first
pop <- pop |>
  get_table("ind_meta") |>
  filter(gen == 5L) |>
  add_dosage(chip_name = "50K", overwrite_dosage = TRUE)

# Marker-assisted selection: keep only animals homozygous favorable at QTL_1
selected_ids <- pop |>
  get_table("ind_genotype") |>
  filter(locus_name == "QTL_1", dosage_value == 2L) |>
  pull(id_ind)

# MAS combined with index — pipe the MAS-filtered set straight to add_phenotype()
pop <- pop |>
  get_table("ind_genotype") |>
  filter(locus_name == "QTL_1", dosage_value == 2L) |>
  add_phenotype("ADG")
```

**Key design properties:**
- `add_founders()` and `add_offspring()` never write to `ind_genotype` — no hidden sync burden
- `ind_genotype` may contain rows for only a subset of generations, loci, or individuals — this is intentional and correct
- `extract_genotypes()` for GBLUP software computes from `ind_haplotype` directly (GROUP BY + PIVOT) and does not depend on `ind_genotype` being populated
- For dominance and epistasis TBV (future), `add_tbv()` will compute dosage inline within the TBV query rather than reading from `ind_genotype`

---

### `genome_meta` additions

Add one column to support mutation provenance:

```text
introduced_gen  INTEGER  DEFAULT NULL
```

- `NULL` = founding locus (present at `initialize_genome()`)
- Integer = generation at which a novel mutation was added via a future `add_mutation()` call

This lets users filter to founding QTL only (`WHERE introduced_gen IS NULL`) or identify segregating novel variants.

This is in addition to the `genome_meta.locus_name` uniqueness requirement described in the
`ind_haplotype` schema section.

---

### `founder_haplotypes` (long)

The `founder_haplotypes` table is a pool of real or simulated haplotypes created by `define_founder_haplotypes()` and sampled by `add_founders()`. It is currently wide format (one column per locus) and must be refactored to long format alongside `ind_haplotype`.

> **Fixed from the original draft**: the schema below restores `line_name`, which the current
> (wide) `founder_haplotypes` table already has (`R/define_founder_haplotypes.R`) and which
> `add_founders()` already relies on for line-specific founder pools (see its examples piping
> from `get_table("founder_haplotypes") |> filter(line_name == ...)`). The first draft of this
> long schema omitted it, which would have broken multi-line/crossbreeding founder setup entirely.

```text
line_name      VARCHAR    NULL = single/unlabeled pool; set for line-specific pools; logical key part
haplotype_id   INTEGER    Logical key part; sequential identifier for each haplotype in the pool
                          (unique within line_name, not necessarily globally)
locus_name     VARCHAR    Logical key part; FK to genome_meta.locus_name
allele         UTINYINT   0 or 1
```

**Logical row key**: `(line_name, haplotype_id, locus_name)` — `haplotype_id` is scoped per line,
matching the current wide table's behavior where each `define_founder_haplotypes(..., line_name =
...)` call manages its own pool. **This cannot be a literal DuckDB `PRIMARY KEY`/`UNIQUE`
constraint**, since `line_name` is nullable and PK/UNIQUE columns must be `NOT NULL`. This is the
same situation `index_meta` already has today (`index_name` is nullable, and per this project's
existing convention its uniqueness is "enforced in R code," not a DB constraint — see
`define_index()`). Apply the same pattern here: enforce `(line_name, haplotype_id, locus_name)`
uniqueness in R (e.g. in `define_founder_haplotypes()`/`import_founder_haplotypes()` before
insert), not via a DDL constraint.

`add_founders()` samples `haplotype_id` values from the pool matching the founder's `line_name` (with or without replacement depending on pool size vs. demand) and writes the selected alleles into `ind_haplotype` with the correct `id_ind`, `parent_origin`, `strand`, `line_origin` (set from the founder's `line_name`), `locus_name`, and `allele`. The `strand` value is determined by which `(parent_origin, strand)` slot within `chr_meta` the sampled haplotype fills.

**Why no `locus_id` here, unlike `ind_haplotype`/`ind_genotype`**: this is a deliberate asymmetry, not an oversight. `locus_id` was added to those two tables because they're scanned at up to 500M-row scale on every generation (recombination, TBV, export). `founder_haplotypes` is a much smaller, cold(er) table — a fixed pool sampled a handful of times during `add_founders()`, not scanned per-generation — so the join/sort cost `locus_id` optimizes for doesn't apply here with the same weight. `locus_name` alone is sufficient.

**`genome_meta.founder_allele_freq` stays as-is for this refactor.** Today it's a single
per-locus column populated "last founder pool wins" — already weak for multi-line simulations,
but that is a separate issue from the wide-to-long storage rewrite. Per the resolved decision in
"Open Questions," leave it in `genome_meta` for now and address line-aware base-frequency storage
as part of the upcoming QTL effect-sampling revision.

**`import_founder_haplotypes()` — new function**

Users will almost always have founder haplotype data in wide format — one row per haplotype, one column per locus — coming from PLINK `.raw`, VCF, or a simulation tool. A dedicated import function should accept either format and pivot internally before writing to the long `founder_haplotypes` table:

```r
# Wide input (one row per haplotype, columns are locus names)
pop <- pop |> import_founder_haplotypes(wide_df)

# Long input (already in long format)
pop <- pop |> import_founder_haplotypes(long_df, format = "long")
```

The function should:
1. Detect or accept `format = "wide"` (default) or `"long"`
2. For wide input: validate that column names match `genome_meta.locus_name`, pivot to long via DuckDB `UNPIVOT` or `tidyr::pivot_longer()`, then INSERT
3. For long input: validate that `locus_name` and `allele` columns are present, validate allele values are 0/1, then INSERT
4. Assign sequential `haplotype_id` values starting after any existing rows

This is the one place in the package where wide input is the expected norm and the function explicitly handles the conversion, keeping all internal storage in long format.

---

## Sex Chromosomes and Uniparental Inheritance

> **Status: Stage 4, out of v1 scope.** v1 implements diploid autosomes only, with `chr_meta`
> created and populated with default autosome rows (`copies_M = 2, copies_F = 2, hemi_parent =
> NULL, recombines = TRUE`) so the table exists and the column shape is stable, but
> `define_chr()` and any non-default inheritance rule is not implemented or tested in v1. This
> section is kept as a design draft for Stage 4 — several details below are known-incomplete
> (see callouts) and should not be treated as ready to implement as written. In particular:
> `copies_M`/`copies_F` are described as "0, 1, or 2" here but reused as "4, 6, 8" for polyploids
> later in this doc — that's an unresolved terminology clash, not a typo, and needs to be settled
> before `chr_meta` is finalized. Pseudoautosomal regions and per-segment (rather than
> per-chromosome) inheritance rules are also not representable by this schema and are not
> addressed here.

### The general problem

Some chromosomes break the standard diploid rule of two copies — one from each parent. The patterns that matter for breeding simulation:

| Chromosome | Species | Males | Females | Single copy always from |
|---|---|---|---|---|
| X | Mammals | 1 copy | 2 copies | dam (when male) |
| Y | Mammals | 1 copy | 0 copies | sire |
| Z | Birds | 2 copies | 1 copy | sire (when female) |
| W | Birds | 0 copies | 1 copy | dam |
| MT | Mammals | 1 copy | 1 copy | dam (always) |
| Chloroplast | Plants | 1 copy | 1 copy | sire or dam (species-specific) |

The long-format `ind_haplotype` handles all of these naturally — hemizygous chromosomes simply have **one row** per individual instead of two, and absent chromosomes (Y in females, W in males) have **zero rows**. No new columns are needed in `ind_haplotype` itself.

What is needed is a chromosome-level rule table so that `add_founders()` and `add_offspring()` know how many haplotypes to generate per individual per chromosome, and where each comes from.

---

### New table: `chr_meta`

One row per chromosome. Stores the inheritance rule for that chromosome. Autosome rows are written automatically by `initialize_genome()`. Users call `define_chr()` to set non-standard rules.

```text
chr_name      VARCHAR    PK; matches genome_meta.chr_name (e.g. "1", "X", "Y", "MT")
copies_M      UTINYINT   Copies carried by males: 0, 1, or 2 (default 2)
copies_F      UTINYINT   Copies carried by females: 0, 1, or 2 (default 2)
hemi_parent   VARCHAR    When an individual has exactly 1 copy: "parent_1" (sire) or
                         "parent_2" (dam); NULL for fully diploid chromosomes
recombines    BOOLEAN    TRUE = recombination occurs during gamete formation;
                         FALSE for Y, W, MT, most organelles (default TRUE)
```

**`hemi_parent` semantics**: this is the parent whose haplotype is passed to an offspring that carries exactly one copy of this chromosome. It is the same for all individuals that carry one copy:

- X (mammal males have 1 copy): `hemi_parent = "parent_2"` — the sire's single X came from his dam
- Y: `hemi_parent = "parent_1"` — inherited from the sire
- Z (bird females have 1 copy): `hemi_parent = "parent_1"` — the dam's single Z came from her sire
- W: `hemi_parent = "parent_2"` — inherited from the dam
- MT: `hemi_parent = "parent_2"` — inherited from the dam (egg) for all individuals

**Default rows** written by `initialize_genome()` for all chromosomes:
`copies_M = 2, copies_F = 2, hemi_parent = NULL, recombines = TRUE`

**Polyploidy**: `copies_M` and `copies_F` encode total chromosome copies — set them to 4 for a tetraploid, 6 for a hexaploid, 8 for an octoploid. `add_founders()` and `add_offspring()` compute `strands_per_parent = copies / 2` and write rows for `strand = 1 .. strands_per_parent` under each `parent_origin`. A tetraploid therefore gets `parent_origin ∈ {1, 2}` × `strand ∈ {1, 2}` = 4 rows per locus; a hexaploid gets `parent_origin ∈ {1, 2}` × `strand ∈ {1, 2, 3}` = 6 rows. The `hemi_parent` and `recombines` columns remain valid for polyploid sex chromosomes (e.g. a tetraploid plant with sex chromosomes) without modification.

---

### New function: `define_chr()`

```r
# Mammalian sex chromosomes
pop <- pop |>
  define_chr("X", copies_M = 1L, copies_F = 2L,
             hemi_parent = "parent_2", recombines = TRUE)  |>
  define_chr("Y", copies_M = 1L, copies_F = 0L,
             hemi_parent = "parent_1", recombines = FALSE)

# Avian sex chromosomes (ZW system — females are heterogametic)
pop <- pop |>
  define_chr("Z", copies_M = 2L, copies_F = 1L,
             hemi_parent = "parent_1", recombines = TRUE)  |>
  define_chr("W", copies_M = 0L, copies_F = 1L,
             hemi_parent = "parent_2", recombines = FALSE)

# Mitochondria (maternal in mammals)
pop <- pop |>
  define_chr("MT", copies_M = 1L, copies_F = 1L,
             hemi_parent = "parent_2", recombines = FALSE)

# Paternal plastid (some plant species)
pop <- pop |>
  define_chr("Chloroplast", copies_M = 1L, copies_F = 1L,
             hemi_parent = "parent_1", recombines = FALSE)
```

`define_chr()` upserts a row in `chr_meta`. Calling it before `add_founders()` is required — founding haplotypes are sampled according to these rules.

---

### How `add_founders()` changes

For each founder and each chromosome in `chr_meta`:
1. Look up `copies_M` or `copies_F` based on the founder's `sex`; compute `strands_per_parent = copies / 2`
2. If `copies = 0`: insert no haplotype rows for this chromosome
3. If `copies = 1`: sample one haplotype; assign `parent_origin` from `hemi_parent` (`"parent_1"` → 1, `"parent_2"` → 2) and `strand = 1`
4. If `copies >= 2`: for each `parent_origin ∈ {1, 2}` and each `strand ∈ {1 .. strands_per_parent}`, sample one independent haplotype — `copies` haplotypes total per founder per chromosome

Founders have no real parents, so `hemi_parent` just determines the `parent_origin` label used for consistency with offspring inheritance. For diploids `strand = 1` always and the resulting rows match the prior two-haplotype semantics exactly.

### Meiosis semantics — `parent_origin` means something different depending on whose row you're reading

**This is the highest-risk correctness issue in the plan and must be unambiguous before any code
is written.** `parent_origin` is not "the sire's contribution" or "the dam's contribution" as a
property of a stored row in general — it depends on whose individual you're looking at:

- **Within an existing individual's own `ind_haplotype` rows**, `parent_origin` records which of
  *that individual's* two parents each copy came from (grandparental origin, from the perspective
  of that individual's own offspring). It says nothing about which mating this individual is
  currently contributing a gamete to.
- **In a newly written offspring row**, `parent_origin` records which *mating* parent (sire =
  `id_parent_1` vs. dam = `id_parent_2`) contributed the gamete that produced that row.

These are different axes. **Gamete generation for a parent must load all of that parent's own
chromosome copies for the chromosome being recombined — i.e. both of that parent's own
`parent_origin ∈ {1, 2}` rows (and, for polyploids, all `strand` values within each) — not a
single `parent_origin` slice.** The earlier draft of step 4 below read as if you could recombine
"within one `parent_origin` group," which is wrong: for a diploid parent that group has exactly
one row per locus, so there would be nothing to recombine against.

### How `add_offspring()` changes

For each mating and each chromosome:
1. Determine offspring sex, look up `chr_meta` rule; compute `strands_per_parent = copies / 2`
2. **If `copies = 0`**: skip — no haplotype rows
3. **If `copies = 1`**: identify the contributing parent from `hemi_parent`. Load **all** of that
   parent's own stored haplotype rows for this chromosome (both of that parent's own
   `parent_origin` values, all `strand` values). If `recombines = FALSE`, pass one of that
   parent's stored copies through intact (no crossover). If `recombines = TRUE` (e.g. X from a
   dam who carries two X copies), recombine across all of the dam's stored copies for this
   chromosome and produce one gamete. Write the result to the offspring's `ind_haplotype` row
   with `parent_origin = 1` (if the contributing parent was the sire, `id_parent_1`) or
   `parent_origin = 2` (if the dam, `id_parent_2`), `strand = 1`. This offspring-row
   `parent_origin` value describes the mating, not the grandparental origin recorded on the
   parent's own rows.
4. **If `copies >= 2`**: for each mating parent (sire, then dam), load **all** of that parent's
   own stored haplotype rows for this chromosome — every `(parent_origin, strand)` combination
   that parent carries — and recombine across all of them to produce `strands_per_parent` gamete
   strands from that one parent. Write the sire's gamete strands to the offspring as
   `parent_origin = 1, strand = 1 .. strands_per_parent`, and the dam's as
   `parent_origin = 2, strand = 1 .. strands_per_parent`.

For diploids (`copies = 2`), `strands_per_parent = 1` — recombination happens across the parent's
own two stored rows (`parent_origin = 1` and `parent_origin = 2` *as recorded on that parent*),
producing one recombined gamete strand, identical to the prior two-haplotype model. The
offspring's `parent_origin` is then reassigned to reflect the mating (1 = from sire, 2 = from
dam), which is a relabeling, not a continuation of the parent's own `parent_origin` values.

**Polyploid pairing is explicitly unresolved and out of v1 scope** (Stage 4). Labeling homologs
as `parent_origin ∈ {1,2}` × `strand ∈ {1..n}` is a reasonable storage scheme, but which strands
are allowed to pair during meiosis (autopolyploid: any two; allopolyploid: only homeologous sets)
is a separate algorithmic question this plan does not answer. See Open Questions and Stage 4.

### Impact on `ind_genotype`

Hemizygous loci have `dosage_value ∈ {0, 1}` (one allele, not a sum of two). Absent loci (e.g. Y loci in a female) have no row in `ind_haplotype` and therefore no row in `ind_genotype` if `add_dosage()` is called. Code that computes the genomic relationship matrix or allele frequencies must be ploidy-aware for these loci. The `chr_meta` table provides the necessary copies-per-sex metadata for any such query.

---

## How Crossbreeding TBV Works with Long Format

> **Centering was missing from the original draft of this section and has been added below.**
> Current `add_tbv()` centers additive TBV using `base_allele_freq` (Falconer centering):
> non-imprinted is `genotype %*% a - 2 * sum(p_base * a)`; imprinted (`expressed_parent =
> "parent_1"`/`"parent_2"`) is `haplotype %*% a - sum(p_base * a)`. The long-format rewrite must
> reproduce both exactly, including when `base_allele_freq` differs by line.

### The crossbreeding problem

A Duroc × Landrace F1 offspring has:
- Paternal haplotype (origin 1) from a Duroc sire → `line_origin = "Duroc"` at each locus
- Maternal haplotype (origin 2) from a Landrace dam → `line_origin = "Landrace"` at each locus

Line-specific effects live in `genome_effects.line_name`. The correct centered TBV for this F1 is:

```text
TBV = Σ_loci [ (allele_hap1 − p_base_Duroc(locus))    × effect_Duroc(locus)
             + (allele_hap2 − p_base_Landrace(locus))  × effect_Landrace(locus) ]
```

### Centered TBV SQL with long format

**Key idea**: fold `(allele - base_allele_freq)` into the summed term itself, using whichever
`genome_effects` row (line-specific or population-wide fallback) matched that haplotype row. This
naturally generalizes Falconer centering to per-line `base_allele_freq` without a separate
centering constant, and collapses to the current formula exactly in the population-wide case (see
below).

**Population-wide effects, both parents expressed** (`expressed_parent = "both"`,
`e.line_name IS NULL` everywhere — the non-crossbreeding case):

```sql
SELECT h.id_ind,
       SUM((h.allele - e.base_allele_freq) * e.genome_value) AS tbv_value
FROM ind_haplotype h
JOIN genome_effects e
  ON  h.locus_name          = e.locus_name
  AND e.trait_name          = 'ADG'
  AND e.genome_effect_type  = 'additive'
  AND e.line_name           IS NULL
GROUP BY h.id_ind
```

Summing `(allele - p) * a` over both of an individual's haplotype rows at a locus is algebraically
identical to `(allele_hap1 + allele_hap2 - 2p) * a` = `genotype * a - 2p * a` — the current
formula, exactly.

**Line-specific effects with population-wide fallback** (the crossbreeding case):

```sql
SELECT h.id_ind,
       SUM((h.allele - e.base_allele_freq) * e.genome_value) AS tbv_value
FROM ind_haplotype h
JOIN genome_effects e
  ON  h.locus_name          = e.locus_name
  AND e.trait_name          = 'ADG'
  AND e.genome_effect_type  = 'additive'
  AND (
       e.line_name = h.line_origin              -- line-specific effect for this allele's line
    OR (e.line_name IS NULL AND NOT EXISTS (     -- fall back to pop-wide only if no line-specific
         SELECT 1 FROM genome_effects e2         -- row exists for this locus x line
         WHERE e2.locus_name        = h.locus_name
           AND e2.trait_name        = 'ADG'
           AND e2.genome_effect_type = 'additive'
           AND e2.line_name         = h.line_origin
       ))
  )
GROUP BY h.id_ind
```

This is impossible to express correctly in wide format without pivoting first. It also resolves
the edge cases flagged in review:
- **`h.line_origin IS NULL`** (untracked origin): `e.line_name = h.line_origin` is `NULL` (never
  true) for any line-specific row, and the `NOT EXISTS` subquery's `e2.line_name = h.line_origin`
  is likewise never true, so `NOT EXISTS` is always true — the row correctly falls back to the
  population-wide effect.
- **Both a line-specific and a population-wide row exist for the same locus/line/trait**: the
  `NOT EXISTS` guard excludes the population-wide row whenever a line-specific one matches, so
  exactly one effect row joins per haplotype row — line-specific wins, as intended.
- **A line has line-specific effects at some loci but not others**: the fallback is evaluated
  per-locus (the `NOT EXISTS` subquery is correlated on `locus_name`), so partial line-specific
  coverage falls back correctly locus-by-locus.
- **`base_allele_freq` differs by line**: each haplotype row's centering term uses the
  `base_allele_freq` from whichever row it actually joined to (line-specific or fallback), so this
  is correct by construction rather than needing separate handling.

**Imprinted expression** (`trait_meta.expressed_parent = "parent_1"` or `"parent_2"`): restrict
`h.parent_origin` to the expressed parent before joining, and center with a single
`base_allele_freq` (matching the current `haplotype %*% a - sum(p_base * a)` formula):

```sql
SELECT h.id_ind,
       SUM((h.allele - e.base_allele_freq) * e.genome_value) AS tbv_value
FROM ind_haplotype h
JOIN genome_effects e
  ON  h.locus_name          = e.locus_name
  AND e.trait_name          = 'ADG'
  AND e.genome_effect_type  = 'additive'
  AND e.line_name           IS NULL   -- or the line-specific/fallback OR-clause above
WHERE h.parent_origin = 1  -- 1 for expressed_parent = "parent_1", 2 for "parent_2"
GROUP BY h.id_ind
```

**What this design must still preserve** (carried over from current `add_tbv()` and not
mechanically solved by the SQL above — needs explicit test coverage in Stage 2):
- filtered `tidybreed_table` subsets (candidate individual set from `get_table()` + `filter()`);
- `trait_meta.expressed_parent` branch selection (which of the two SQL variants above runs);
- existing `ind_tbv` upsert semantics;
- `add_tbv(index_names = ...)` — unaffected by this change since it consumes the already-computed
  `tbv_value` and multiplies by index weights, same as today;
- formula/composite phenotype code paths in `add_phenotype()` that read from `ind_tbv`.

`R/define_additive_effects.R` (including its `base = "founder_haplotypes"` and
`base = "current_pop"` helpers, currently `founder_haplotype_helpers.R`) is a rewrite target here
too, not just `add_tbv.R` — both currently compute `base_allele_freq` from wide tables.

---

## Recombination in `add_offspring()` (Long Format)

> **How this relates to "How `add_offspring()` changes" under Sex Chromosomes above**: that
> section describes the per-chromosome *branching* logic driven by `chr_meta` (how many copies,
> whether it recombines, which parent contributes a hemizygous copy) — it's Stage 4 scoped. This
> section describes the underlying per-parent recombination *mechanics* each branch invokes (query
> haplotypes, build arrays, call `make_gamete()`, write results) — this is the Stage 1 diploid
> baseline that Stage 4 later generalizes. They are complementary, not duplicate or conflicting
> specs: Stage 4's branching logic decides how many times to invoke this mechanic and with which
> rows; this section is what that invocation actually does.

The current logic loads parent haplotypes as a matrix and calls `make_gamete()`. The long format version:

1. Query parent haplotypes including `id_ind`, `strand`, and `line_origin`:
   ```sql
   SELECT id_ind, parent_origin, strand, locus_id, locus_name, allele, line_origin
   FROM ind_haplotype
   WHERE id_ind IN ('Duroc_1', 'Landrace_5')
   ORDER BY id_ind, parent_origin, strand, locus_id
   ```

2. Build haplotype arrays in R **grouped by `id_ind`**, i.e. per parent — not per mating-side
   label. Within one parent's group, further split by that parent's own `parent_origin` value
   (their two grandparental homologs, for diploids) and `strand` (for polyploids). Recombination
   for a given parent happens **across that parent's own `parent_origin ∈ {1, 2}` rows** — see
   "Meiosis semantics" above. Do not filter to a single `parent_origin` before recombining; that
   would discard one of the parent's two homologs entirely.

3. `make_gamete(strands, chr_map)` recombines across a parent's homolog arrays (across that
   parent's own `parent_origin` values, and — for polyploids — across `strand` values within each)
   and returns one gamete allele vector plus the corresponding `line_origin` vector (inherited
   from whichever parental homolog contributed each locus after recombination).

4. Write offspring haplotypes as long rows with `(parent_origin, strand, locus_id, locus_name,
   allele, line_origin)`, where `parent_origin` is now reassigned to mean "from the sire" (1) or "from the
   dam" (2) for this mating — a relabeling relative to whatever `parent_origin` values were
   recorded on the parent's own rows. For diploids `strand = 1` always; for polyploids, one
   gamete strand is written per `strand` slot the parent contributes.

The recombination breakpoints are chromosome-position-based (same as today). Between breakpoints, all loci share the same contributing parental haplotype and therefore the same `line_origin`. This makes `line_origin` inheritance correct even for F2 and backcross animals.

---

## Novel Mutation Support

> **Status: Stage 5, out of v1 scope.** v1 is fully dense (see "Proposed Schema" above) — there
> is no sparse/implicit-zero convention to rely on, so `add_mutation()` cannot be implemented as
> sketched below without either (a) backfilling an explicit `allele = 0` row for every existing
> individual at the new locus, which is the O(n_ind) cost this section was trying to avoid, or
> (b) introducing sparse storage, which is a real design change deferred to its own stage. This
> section is kept only to record the shape of the future design; do not implement it against the
> v1 dense schema.

When `add_mutation()` is implemented (future, and pending a decision on dense vs. sparse storage
for introduced loci — see Stage 5):

1. Add new row to `genome_meta` (assign next `locus_id`, set `introduced_gen`)
2. Add rows to `ind_haplotype` for the mutant individual — one per `(parent_origin, strand)` combination; allele = 1 on the mutant strand, allele = 0 on all others — with `line_origin` inherited from that individual
3. If sparse storage is adopted for introduced loci: no rows added for any other individual (implicit allele 0). If dense storage is kept throughout: backfill `allele = 0` rows for every other individual at the new locus.
4. Do **not** write to `ind_genotype` — the user calls `add_dosage()` explicitly if they want dosage values for the new locus
5. Subsequent `add_offspring()` calls for that individual's children will pick up the new locus naturally — the locus appears in the parent's `ind_haplotype` rows and will be inherited with correct Mendelian probabilities

If sparse storage is adopted, every downstream query (dosage, TBV, allele frequency, `PIVOT`
export, parent-haplotype loading for recombination) must consistently use
`LEFT JOIN ... COALESCE(allele, 0)` — this needs to be true everywhere, not just in the queries
this plan happens to show, which was the core problem with the earlier hybrid dense/sparse draft.

---

## Genotype Matrix Export for GBLUP

The genomic relationship matrix and BLUPF90 both require a wide n×p matrix. `extract_genotypes()` always computes this on demand from `ind_haplotype` — it does not depend on `ind_genotype` being populated. DuckDB's native `PIVOT` makes this clean:

```sql
PIVOT (
  SELECT h.id_ind, h.locus_id, h.locus_name, CAST(SUM(h.allele) AS UTINYINT) AS dosage_value
  FROM ind_haplotype h
  WHERE h.locus_name IN (SELECT locus_name FROM genome_meta WHERE is_50K = TRUE)
    AND h.id_ind IN (/* genotyped individuals from tbl */)
  GROUP BY h.id_ind, h.locus_id, h.locus_name
) ON locus_name USING FIRST(dosage_value)
ORDER BY id_ind
```

Column names in the pivoted output are `locus_name` values (human-readable, matching the current
BLUPF90/GBLUP convention), but column *order* must be driven by `locus_id`, not incidental
`PIVOT` output order — DuckDB's `PIVOT` does not guarantee this on its own. Implement column
ordering explicitly (e.g. an ordered `IN (...)` locus list generated from
`SELECT locus_name FROM genome_meta WHERE is_50K = TRUE ORDER BY locus_id`, passed into the
`PIVOT ... ON locus_name IN (...)` form) rather than relying on default behavior.

The result is the standard 0/1/2 (or 0–ploidy) wide matrix required by external software. It is always authoritative because it reads directly from the haplotype source of truth rather than from a potentially partial `ind_genotype` cache.

**Acceptance criteria for the `extract_genotypes()`/BLUPF90 rewrite** (the `PIVOT` sketch above
does not yet guarantee all of these — add explicit tests for each):
- exactly one row for every requested individual, including those with all-reference genotypes;
- exactly one column for every requested locus;
- zero (not missing/NULL) for any locus a dense v1 individual has no non-reference allele at;
- stable column order, driven by `genome_meta.locus_id` (or equivalent), not join order;
- safe handling of arbitrary `locus_name` values in dynamically-generated SQL (no injection risk
  from user-supplied locus names reaching a `PIVOT`/column-name context — parameterize or
  validate against `genome_meta` first);
- ploidy-aware dosage for non-diploid/sex-linked loci is explicitly out of v1 scope (Stage 4) —
  until then, `extract_genotypes()` should assume/assert diploid autosomes only.

### New `loci_tbl` argument — a general locus filter, independent of chips or QTL sets

`extract_genotypes()` (already implemented today, `R/extract_genotypes.R`) currently supports two
loci-selection paths: `chip_name` (via `is_<chip>` flags in `genome_meta`, requires
`define_chip()` and physically marking animals genotyped via `add_genotypes()`) and `effects_tbl`
(a filtered `get_table("genome_effects")`, for QTL sets). Neither lets a user say "just give me
these loci" without first defining a chip or supplying QTL effects — there's no way today to
restrict to, say, autosomes only, without conflating that with genotyping-panel or QTL semantics.

Add a third, independent path:

```r
extract_genotypes(tbl, chip_name = NULL, effects_tbl = NULL, loci_tbl = NULL, col_name = NULL)
```

- `loci_tbl` — a `tidybreed_table` from `get_table("genome_meta")`, optionally filtered (e.g.
  `filter(chr_name %in% autosome_names)` or `filter(!chr_name %in% c("X", "Y", "MT"))`). The
  collected table's `locus_name` values become the locus set for this path. Mirrors how
  `effects_tbl` already works — a `tidybreed_table` from a specific source table, not a new
  convention.
- When multiple of `chip_name`/`effects_tbl`/`loci_tbl` are supplied, the locus sets are
  **unioned** (deduplicated, ordered by `locus_id`) — same documented behavior as the existing
  chip + effects_tbl combination.

```r
# Autosomes only, no chip or QTL definition required
geno <- pop |>
  get_table("ind_meta") |>
  extract_genotypes(
    loci_tbl = get_table(pop, "genome_meta") |> dplyr::filter(!chr_name %in% c("X", "Y", "MT"))
  )

# Chip loci restricted to autosomes (union then filter — or filter chip loci directly via loci_tbl)
geno <- pop |>
  get_table("ind_meta") |>
  extract_genotypes(chip_name = "50K", loci_tbl = get_table(pop, "genome_meta") |>
    dplyr::filter(!chr_name %in% c("X", "Y", "MT")))
```

**This de-risks Open Question on G-matrix ploidy for sex-linked loci (see "Open Questions"
below) without waiting for Stage 4.** Users get full manual control over which chromosomes feed
a GBLUP matrix today, via a plain `genome_meta` filter — independent of whether `chr_meta`/Stage 4
ploidy tracking exists yet. Stage 4 can later add a convenience (e.g. auto-deriving "autosomes
only" from `chr_meta` so users don't have to type the chromosome list themselves), but that
becomes a nice-to-have, not a blocking correctness requirement for v1.

---

## Storage, Memory, and Performance

> **These are unverified estimates, not benchmarked numbers.** Per critique, do not treat the
> figures below as settled — run the benchmark described at the end of this section (also a
> Stage 1 exit criterion) before locking the schema.

### Disk storage

Long format is naively larger than wide: every row repeats `id_ind`, `locus_name`, and `line_origin` alongside the allele. DuckDB's columnar compression recovers most of the overhead:

| Column | Compression mechanism | Why it works |
|--------|----------------------|--------------|
| `allele` | Bit-packing (1 bit/value) | Binary — only 0 or 1 |
| `locus_name` | Dictionary encoding | 50K unique values, referenced ~10K times each; the dictionary is tiny |
| `id_ind` | Dictionary + RLE | Each individual appears exactly `ploidy × n_loci` times in sorted order |
| `line_origin` | Dictionary encoding | Very few unique values (2–10 lines typical) |
| `parent_origin` | RLE | Always 1 or 2; long runs per individual in sorted order |
| `strand` | RLE / bit-packing | Always 1 for diploids; low cardinality (1–4) for polyploids |

Estimated on-disk footprint vs. wide format for 5,000 individuals × 50,000 loci:

| Table | Rows | Wide format | Long format (compressed) | Notes |
|-------|------|-------------|--------------------------|-------|
| `ind_haplotype` | 500M | ~60 MB | ~200–500 MB | always populated |
| `ind_genotype`  | 0–250M | ~30 MB | 0–~100 MB | on-demand only; often empty or partial |

Wide format wins on raw disk size because it stores only allele bits with no key overhead. Long format is roughly 3–5× larger even after compression. The gap narrows for simulations with many lines (more unique `line_origin` values compress less well in wide anyway via column proliferation).

**Sort order is the single most important storage tuning decision.** Inserting rows in `(id_ind, parent_origin, strand, locus_id)` order allows DuckDB's zone maps and RLE to compress maximally — all alleles for one individual's strand are contiguous, and `locus_id` (integer) sorts and compresses more cheaply than `locus_name` (string) would. If the dominant query is TBV computation (scan all individuals at QTL loci), sorting by `(locus_id, id_ind, parent_origin, strand)` is more cache-friendly. Choose `(id_ind, parent_origin, strand, locus_id)` as the primary physical sort since offspring generation (read one parent's strands at a time) is the hot path; TBV queries are a full-table scan regardless of sort order.

**For very large simulations (>50K individuals)**: restrict `ind_haplotype` to QTL and chip loci only; retain `pos_Mb` in `genome_meta` for the recombination map. Non-functional SNPs not on any panel or QTL list need not be stored at all.

---

### Memory (RAM) during `add_offspring()`

The memory bottleneck is loading parent haplotypes into R. For 200 matings with 100 unique parents × 50K loci × 2 haplotypes = 10M rows pulled into R memory (~500 MB as a data frame with character columns).

Mitigation strategies (in order of impact):

1. **Cache unique parents** — fetch each unique parent's haplotypes once, not once per mating. A sire used in 50 matings should be fetched once. The current code already does this partially; the long-format rewrite must preserve it explicitly.

2. **Process matings in chunks** — add an internal `chunk_size` parameter (default ~500 matings) so `add_offspring()` processes one chunk at a time, releasing parent haplotypes from R memory before loading the next chunk. Expose as `add_offspring(..., chunk_size = 500L)` for user tuning.

3. **Pull only needed loci** — when computing TBVs during offspring generation (if `add_tbv()` is called internally), load only QTL loci for the haplotype matrix, not all 50K loci. The full haplotype still needs to be stored, but the in-memory working set is smaller.

4. **Stream via Arrow** — DuckDB's Arrow integration allows reading query results as Arrow RecordBatches without materializing an R data frame. For the Rcpp path (see below) this is the right approach; even in pure R, `duckdb::duckdb_fetch_arrow()` reduces peak memory by processing batches.

---

### INSERT speed (`add_founders()` and `add_offspring()`)

For 200 offspring × 50K loci × 2 haplotypes = 20M rows inserted per generation. DuckDB INSERT throughput matters.

Insert methods ranked fastest to slowest:

| Method | Approx throughput | Notes |
|--------|-------------------|-------|
| DuckDB C++ `Appender` (via Rcpp) | 50–200M rows/s | Fastest; bypasses SQL parser entirely |
| `duckdb_register()` + `INSERT INTO ... SELECT` | 10–50M rows/s | Current approach for wide; works well for long too |
| `DBI::dbWriteTable(..., append = TRUE)` | 5–20M rows/s | Slower due to type conversion overhead |
| SQL `INSERT INTO ... VALUES (...)` | <1M rows/s | Never use for bulk data |

**Recommendation**: keep the `duckdb_register()` + `INSERT INTO ... SELECT` approach for the R implementation. When the Rcpp path is added later, switch to the C++ `Appender` for the hot loop.

Ensure the data frame is sorted by `(id_ind, parent_origin, strand, locus_id)` before registering — DuckDB can write sorted data more efficiently and subsequent reads will be faster.

---

### Benchmark plan (required before locking the schema — Stage 1 exit criterion)

Write a standalone DuckDB benchmark script, before implementing the full package refactor,
comparing:

- current wide format (baseline);
- long format with the `(locus_id, locus_name)` dual-key schema adopted above;
- `GROUP BY` TBV on QTL-only loci (the centered SQL above);
- `PIVOT` export for chip-sized loci sets, including the explicit `locus_id`-ordered column list.

Use at least one realistic scale (e.g. 5,000 individuals × 50,000 loci) and one CI-friendly small
scale for regression tests going forward. Storing both `locus_id` and `locus_name` (see schema
section above) means this benchmark no longer needs to adjudicate a "pick one key" decision — it
validates the long format's overall performance against wide, not a fork between key choices.

---

## Appendix: Rcpp / C++ Plan (future optimization — out of scope for v1)

> This section is preserved as future direction but is explicitly **not** an acceptance
> criterion for the first (Stage 1) refactor. Linking Rcpp directly against DuckDB and sharing a
> connection safely between R and C++ is a nontrivial undertaking on its own; the R
> implementation must be complete and correct first, with the same observable behavior the C++
> path would later target. Do not block Stage 1 on any of the below.

The recombination loop in `make_gamete()` is pure computation: given two parent haplotype vectors and a chromosome map, generate one offspring haplotype via crossover sampling. It is:
- Embarrassingly parallel across mating pairs
- Embarrassingly parallel across chromosomes within a mating
- A prime Rcpp + OpenMP target

**C++ can access DuckDB directly.** DuckDB ships a full C++ API (`libduckdb.hpp`). The `duckdb` R package includes the compiled shared library and headers, so Rcpp functions can link against it without any additional system dependency. Key C++ API features for this use case:

```cpp
#include <duckdb.hpp>

// Open the same .duckdb file the R session is using
duckdb::DuckDB db(db_path);
duckdb::Connection conn(db);

// Read parent haplotypes as Arrow (zero-copy into C++ arrays)
auto result = conn.Query(
  "SELECT id_ind, parent_origin, strand, locus_id, locus_name, allele, line_origin "
  "FROM ind_haplotype WHERE id_ind IN ('Duroc_1', 'Landrace_5') "
  "ORDER BY id_ind, parent_origin, strand, locus_id"
);

// Do recombination in C++ (OpenMP across matings)
#pragma omp parallel for
for (int i = 0; i < n_matings; i++) {
  make_gamete_cpp(parent_haps[sire[i]], parent_haps[dam[i]],
                  chr_info, offspring_haps[i]);
}

// Bulk-insert via Appender (fastest DuckDB write path)
duckdb::Appender appender(conn, "ind_haplotype");
for (auto& row : offspring_haplotypes) {
  appender.AppendRow(row.id_ind, row.parent_origin, row.strand,
                     row.locus_id, row.locus_name, row.allele, row.line_origin);
}
appender.Close();
```

**Important constraint**: DuckDB only allows one writer at a time. The Rcpp function must operate on the same connection the R session holds (passed in), or open the file in read-write mode when R has released it. The safest design is a single R-facing wrapper that calls the Rcpp function while temporarily suspending R's DBI connection, or — cleaner — using DuckDB's in-process connection sharing so both R and C++ use the same `duckdb::DuckDB` instance.

**What the Rcpp rewrite of `add_offspring()` would look like architecturally:**
1. R side: validate inputs, build `matings` data frame, call Rcpp function with `db_path` (or connection pointer)
2. C++ side: open connection, fetch unique parent haplotypes, run parallel recombination, bulk-insert offspring haplotypes into `ind_haplotype` via `Appender`, return offspring IDs to R
3. R side: write `ind_meta` rows (small, fast in R)

The long format works *better* for the C++ path than wide format would — you can stream rows without ever materialising a huge matrix, and the `Appender` row-by-row API maps naturally to the long format's one-row-per-allele structure.

**Do not design the R API around Rcpp.** The R implementation should be complete and correct first. The Rcpp path is a drop-in performance replacement for the same functions — same arguments, same return value, same observable behaviour. Add a `backend = c("R", "cpp")` argument to `add_offspring()` and `add_founders()` as the switchover point, defaulting to `"cpp"` once available.

---

## REJECTED FOR v1: Sparse Haplotype Storage

> **Status**: Considered and rejected. Dense storage is the correct choice for v1. This section is preserved so the reasoning is not re-litigated later. Revisit only if novel mutation support (`add_mutation()`) becomes a priority in a future version.

The current plan stores every allele for every individual at every locus in `ind_haplotype` (dense). An alternative is **sparse storage**: only store rows where `allele = 1`. Absence of a row means allele = 0 (reference). This section documents the tradeoffs so the decision can be made with full information.

### The core insight

`allele = 0` contributes nothing to any downstream computation:

- TBV: `genome_value × 0 = 0` — zero-allele rows add nothing to the sum
- Genotype dosage: `SUM(allele)` — zeros don't change the sum
- Allele frequency numerator: `SUM(allele)` — same

No computation requires the 0-rows to be physically present. This is the same principle behind sparse matrix formats in numerical computing.

### Storage savings

For founding loci, savings depend on allele frequency. At mean frequency 0.3 across the genome, roughly 30% of haplotype rows are allele = 1 and get stored; 70% are implicitly 0. Typical savings: **40–60% fewer rows** in `ind_haplotype` for a real breeding population after selection and drift.

For novel mutations the savings are dramatic: a new mutation in one individual requires **1–2 rows** instead of `2 × n_ind` rows. Without sparse storage, `add_mutation()` would need to backfill zero rows for every existing individual — expensive and wasteful.

### How each operation changes

**`add_founders()` and `add_offspring()`**: sample or recombine alleles as today, then filter to `allele = 1` before inserting. Simpler, not harder.

**`add_mutation()` (future)**: insert 1 haplotype row for the mutant haplotype. Nothing else. Every other individual implicitly has allele = 0. No backfill.

**TBV computation**: the Falconer-centered TBV splits cleanly into a sparse sum and a constant:

```text
TBV_i = Σ_k (genome_value_k × allele_ik)   ← sparse: only allele=1 rows contribute
       − Σ_k (2 × base_allele_freq_k × genome_value_k)  ← constant from genome_effects, no haplotype scan
```

SQL with sparse haplotype:
```sql
-- sparse sum (LEFT JOIN from ind_meta catches individuals with all-zero alleles)
SELECT i.id_ind, COALESCE(SUM(e.genome_value * h.allele), 0) AS raw_tbv
FROM ind_meta i
LEFT JOIN ind_haplotype h ON h.id_ind = i.id_ind
LEFT JOIN genome_effects e
  ON h.locus_name = e.locus_name AND e.trait_name = 'ADG' AND e.line_name IS NULL
GROUP BY i.id_ind

-- centering constant (no haplotype scan)
SELECT SUM(2.0 * base_allele_freq * genome_value) AS center
FROM genome_effects WHERE trait_name = 'ADG'
```

Final TBV = raw_tbv − center. Works correctly for all individuals including those with all reference alleles.

**Allele frequency tracking**: the denominator can no longer be derived from haplotype row counts — it must come from `ind_meta`:

```sql
SELECT locus_name,
       SUM(allele)::DOUBLE / (2.0 * (SELECT COUNT(*) FROM ind_meta WHERE gen = 5)) AS freq
FROM ind_haplotype
WHERE locus_name = 'Locus_42'
GROUP BY locus_name
```

Minor added complexity; easy to encapsulate in a helper.

### The hard case: GBLUP genotype matrix extraction

GBLUP and genomic relationship matrix computation require a **dense** n × p matrix — every individual at every chip locus, including zeros. With sparse `ind_haplotype` those zero cells have no rows to SELECT. The naive fix is a cross join:

```sql
SELECT i.id_ind, m.locus_name, COALESCE(SUM(h.allele), 0) AS dosage_value
FROM ind_meta i
CROSS JOIN (SELECT locus_name FROM genome_meta WHERE is_50K) m
LEFT JOIN ind_haplotype h
  ON h.id_ind = i.id_ind AND h.locus_name = m.locus_name
GROUP BY i.id_ind, m.locus_name
```

For 5,000 individuals × 50,000 chip loci this cross join generates 250M combinations. DuckDB handles this via hash join internally, but it is noticeably slower than selecting from a pre-materialised dense table.

**Proposed mitigation**: keep `ind_genotype` as a **dense materialized cache scoped to chip loci only**. When `define_chip()` is called it materialises dosage rows for that chip from the sparse haplotype table. `add_offspring()` inserts dosage rows (including zeros) for chip loci as part of offspring creation. Novel mutations only enter `ind_genotype` once assigned to a chip.

This means:
- `ind_haplotype` — fully sparse; ground truth for all loci including mutations
- `ind_genotype` — dense, but only for chip loci; fast cache for GBLUP extraction

Note: the adopted design (dense `ind_haplotype`, on-demand `add_dosage()`) already solves this problem without the complexity of a sparse haplotype + dense chip cache.

The maintenance burden (keeping `ind_genotype` in sync) is localised to `define_chip()` and `add_offspring()`.

### Complexity summary

| Operation | Dense (current plan) | Sparse (under consideration) |
|-----------|---------------------|------------------------------|
| `add_founders()` | baseline | simpler — filter to allele = 1 |
| `add_offspring()` | baseline | simpler — filter to allele = 1 |
| `add_mutation()` | O(n_ind) backfill rows | O(1) — one row only |
| TBV computation | baseline | equivalent; LEFT JOIN from ind_meta |
| Allele frequency | COUNT rows / 2 | denominator from ind_meta |
| GBLUP extraction | GROUP BY + PIVOT from ind_haplotype | cross join OR dense chip dosage cache |

### What needs to be decided before implementing

1. Is the storage saving worth the added complexity in `extract_genotypes()` and allele frequency queries?
2. Is the dense `ind_genotype` chip cache the right mitigation, or is the cross join fast enough in practice for the simulation scales we target?
3. Does sparse storage interact with the Rcpp path in any unexpected way? (Preliminary answer: no — the Appender writes whatever rows R/C++ hands it; filtering to allele = 1 before calling Appender is trivial.)

---

## Files Requiring Changes

| File | Change |
|------|--------|
| `R/schema.R` | Update table/column descriptions for `ind_haplotype`, `ind_genotype`, and new `chr_meta` |
| `R/define_genome.R` | New DDL for `ind_haplotype (id_ind, parent_origin, strand, line_origin, locus_id, locus_name, allele)`, `ind_genotype (id_ind, locus_id, locus_name, dosage_value)`, and `chr_meta (chr_name, copies_M, copies_F, hemi_parent, recombines)`; add `UNIQUE` constraint on `genome_meta.locus_name` |
| `R/initialize_genome.R` | Call updated DDL; add `introduced_gen` to `genome_meta` DDL; write default autosome rows to `chr_meta`; create `ind_genotype` empty (not auto-populated) |
| `R/open_pop.R` | Update schema registration; update `_schema_meta` descriptions |
| `R/define_chr.R` | New file — `define_chr()` upserts a row in `chr_meta` |
| `R/define_founder_haplotypes.R` | Refactor internal storage from wide to long `founder_haplotypes` table |
| `R/import_founder_haplotypes.R` | New file — accepts wide or long input, pivots wide internally, writes to long `founder_haplotypes` |
| `R/add_founders.R` | Read `chr_meta` per chromosome; build long haplotype DF respecting `copies_M`/`copies_F`, `hemi_parent`, and `strand`; does **not** write to `ind_genotype` |
| `R/add_offspring.R` | Read `chr_meta` per chromosome; compute `strands_per_parent = copies / 2`; route each chromosome through hemizygous or polyploid gamete logic; `make_gamete(strands, chr_map)` recombines across all strands and returns `(allele_vec, line_origin_vec)`; build long haplotype DF with `strand` column; does **not** write to `ind_genotype` |
| `R/add_dosage.R` | New file — `add_dosage(tbl, chip_name = NULL, locus_names = NULL, overwrite_dosage = FALSE)`; computes `SUM(allele) GROUP BY (id_ind, locus_id, locus_name)` from `ind_haplotype`; inserts into `ind_genotype` via `INSERT OR REPLACE` |
| `R/add_tbv.R` | Replace matrix multiply with the centered SQL `GROUP BY SUM((allele - base_allele_freq) * genome_value)` JOIN on `(locus_name, line_origin)` (see "Centered TBV SQL" above); dominance/epistasis TBV computes dosage inline — does not read from `ind_genotype` |
| `R/define_additive_effects.R` (incl. `founder_haplotype_helpers.R`) | Rewrite `base = "founder_haplotypes"` and `base = "current_pop"` helpers — both currently compute `base_allele_freq` from wide tables |
| `R/phenotype_helpers.R` | Update any code path that reads haplotype/genotype tables to assemble composite phenotypes |
| `R/extract_genotypes.R` | Compute wide matrix via GROUP BY + DuckDB `PIVOT` on `locus_name` from `ind_haplotype` directly — does not depend on `ind_genotype`; add the new `loci_tbl` argument (general locus filter via `get_table("genome_meta")`); must meet the acceptance criteria listed under "Genotype Matrix Export for GBLUP" above |
| `R/blupf90_helpers.R` | Same GROUP BY + PIVOT approach from `ind_haplotype` |
| `R/sql_utils.R` | Update `TABLE_ROW_KEYS` (and reserved-column lists) for `ind_haplotype`, `ind_genotype`, and `chr_meta` |
| `R/remove_rows.R` | DELETE by `id_ind` from `ind_haplotype` and `ind_genotype`; no column names to track |
| `R/summary_pop.R` | Update row count logic (hemizygous chromosomes have fewer rows per individual; `ind_genotype` may be partially populated) |
| `R/archive_replicate.R` | Update any hard-coded references to `genome_haplotype`/`genome_genotype` table names |
| `src/gamete.cpp` | Future Rcpp — `make_gamete_cpp()` with OpenMP parallel recombination across mating pairs |
| `src/add_offspring_cpp.cpp` | Future Rcpp — full offspring generation: DuckDB C++ read → recombination → Appender write to `ind_haplotype` only |

---

## Migration for Existing Databases

> **Rewritten per review.** `tidyr::pivot_longer()` on entire wide tables materializes the full
> wide table in R memory before writing — this does not scale to real database sizes (500M target
> rows). Migration must be SQL-first and transactional, not an R data-frame pivot.

`migrate_haplotype_to_long(pop)` should follow this transactional procedure:

1. **Create new long tables under temporary names** (`ind_haplotype_new`, `ind_genotype_new`)
   with the target schema, including the `PRIMARY KEY` constraints from this plan.
2. **Populate them via DuckDB `UNPIVOT` directly against the existing wide tables**, in chunks if
   needed for large genomes — never round-trip through an R data frame:
   ```sql
   INSERT INTO ind_haplotype_new (id_ind, parent_origin, strand, locus_id, locus_name, allele, line_origin)
   SELECT h.id_ind, h.parent_origin, 1 AS strand, gm.locus_id, gm.locus_name, u.allele,
          m.line_name AS line_origin
   FROM genome_haplotype h
   UNPIVOT (allele FOR locus_col IN (locus_1, locus_2, /* ... locus_n */)) AS u
   JOIN genome_meta gm ON u.locus_col = 'locus_' || gm.locus_id  -- locus_col -> locus_id/locus_name
   LEFT JOIN ind_meta m ON m.id_ind = h.id_ind                    -- line_origin from ind_meta.line_name
   ```
   (The old schema is always diploid, so `strand = 1` for every migrated row.) The old
   `genome_genotype` table migrates the same way into `ind_genotype_new`, renaming its dosage
   column to `dosage_value`. Note `genome_genotype` was auto-populated in the old schema, but
   `ind_genotype` is on-demand-only in the new one — migrating its existing contents is still
   correct (they're valid cached dosage values), just not required for correctness.
3. **Validate** row counts and a checksum (e.g. `SUM(allele)` per locus, or a hash of a sample of
   rows) between old and new tables before proceeding.
4. **Rename old tables to backup names** (`genome_haplotype_backup_<timestamp>`, etc.) — do not
   drop yet.
5. **Rename new tables into place** (`ind_haplotype_new` → `ind_haplotype`, etc.).
6. **Update bookkeeping**: `_schema_meta`, `pop$tables`, `TABLE_ROW_KEYS` assumptions in
   `R/sql_utils.R`, and any archive/replicate metadata that references the old table names.
7. **Drop backups only after explicit user confirmation** or a second successful validation pass
   — never automatically in the same transaction as the rename.

All of steps 1–6 should run inside a single DuckDB transaction where possible so a failure midway
leaves the original wide tables untouched.

**Migration must also cover**, not just the two haplotype/genotype tables:
- `founder_haplotypes` (wide → long, preserving `line_name` — see schema section above);
- chip `is_<chip_name>` columns in `genome_meta` (unaffected in shape, but referenced by the new
  `add_dosage()` chip filter — verify they still resolve correctly);
- `has_<chip_name>` columns in `ind_meta` (unaffected — `add_genotypes()` still writes these);
- existing `genome_effects` rows (unaffected in shape, but this is the first time they're queried
  against a `line_origin`-bearing table — verify joins resolve for pre-existing effect rows);
- archive/reset workflows (`R/archive_replicate.R`) that reference table names directly;
- older databases created before line-aware founder haplotypes existed (i.e. no `line_name`
  column on `founder_haplotypes` at all) — these need a distinct migration path or an explicit
  error directing the user to re-run `define_founder_haplotypes()`.

---

## Open Questions

**Genuine judgment calls — please decide these:**

1. **Polyploid meiosis pairing model** (Stage 4, not urgent now but worth flagging) — for diploids, gamete generation recombines across a parent's own two `parent_origin` rows (see "Meiosis semantics" above). For a tetraploid, which of the parent's four homologs are allowed to pair is a separate algorithmic question: in autopolyploids (AAAA) any two strands can pair; in allopolyploids (AABB, e.g. hexaploid wheat) only homeologous strand sets pair. The schema supports either by having the right number of `strand` rows, but no pairing algorithm is specified. Left for a future `define_ploidy_model()`; no decision needed until Stage 4 begins.

**Resolved during this review (recorded here for reference, not asking again):**

- **Table naming** — rename to `ind_haplotype`/`ind_genotype`. Decision: they store individual information, not metadata, so they follow the same `ind_` convention as `ind_meta`/`ind_phenotype`/`ind_tbv`/`ind_ebv`/`ind_index`. `genome_*` stays reserved for non-individual data (`genome_meta`, `genome_effects`). No compatibility views; rename and break, with the migration procedure above handling existing databases.
- **G matrix ploidy for sex-linked loci** — de-risked for v1 by the new `loci_tbl` argument on `extract_genotypes()` (see "New `loci_tbl` argument" under "Genotype Matrix Export for GBLUP" above): users can already restrict a GBLUP extraction to autosomes via a plain `get_table("genome_meta") |> filter(...)`, independent of `chr_meta`/Stage 4. Automatic ploidy-aware warnings/errors remain a Stage 4 convenience, not a blocking v1 requirement.
- **`locus_name` vs. `locus_id` as the internal key** — store both, rather than picking one. `locus_id INTEGER` is the physical sort/PK key (cheap comparisons for the recombination hot path and compression); `locus_name VARCHAR` is denormalized alongside it purely so TBV/GWAS queries can join straight to `genome_effects` without an extra hop through `genome_meta`. Never exposed differently to users either way — `locus_name` remains what users type in `filter()`, `define_chip()`, `add_dosage(locus_names = ...)`, etc. `locus_name` is treated as immutable once haplotype rows reference it (see caveat in the `ind_haplotype` schema section).
- **`founder_allele_freq`** — leave as-is in `genome_meta` for now ("last founder pool wins"). Decision: this is a separate, more complex issue than the wide→long storage change, and will be addressed as part of an upcoming revision to QTL effect sampling rather than bundled into this refactor.
- **`line_origin` when line is not tracked** — always set it in `add_founders()` when a line is known (cheap; skipping it silently breaks future crossbreeding). `NULL` is reserved for genuinely untracked origin, not for "population-wide effects mode."
- **`line_origin` in F2 / backcross** — remains accurate across generations; `make_gamete()` inherits `line_origin` from whichever parental homolog contributed each locus after recombination (see "Recombination in `add_offspring()`" above).
- **`NULL line_origin` in the TBV query** — the centered SQL's `NOT EXISTS` fallback handles this correctly by construction (see "Centered TBV SQL" above) — no special-casing needed.
- **`ind_genotype` cache scope across repeated `add_dosage()` calls** — add an `overwrite_dosage = FALSE` argument (named to match `add_index()`'s `overwrite_index`), not `delete_prior`.
- **Is mutation support part of this refactor?** — No. Deferred to Stage 5 as its own design pass (see "Novel Mutation Support" and "Staged Implementation" above); not blocking v1.

**Deferred as future design work, not blocking this refactor:**

- **Nucleotide identity (A/C/G/T)** — no reason to store actual nucleotides in `ind_haplotype`. It stores *which* allele (0 = reference, 1 = derived); *what* those alleles are would go in optional `genome_meta.ref_allele`/`alt_allele` columns if VCF export or transition/transversion-biased mutation models are ever needed. Mirrors VCF's `GT` (0/1 indices) + header (REF/ALT bases) split. Keep `allele UTINYINT` as 0/1; not worth revisiting now.
- **Multi-allelic loci** — the design assumes biallelic loci. `UTINYINT` technically allows values up to 255, but `ind_genotype.dosage_value = SUM(allele)` loses interpretability once alleles aren't binary (can't distinguish "two copies of allele 2" from "one copy each of alleles 1 and 3"). Would need a different representation (one row per allele copy, or a `genome_allele_count(id_ind, locus_name, allele_value)` table). Flagged for a future decision; out of scope here.
- **TBV performance at scale** — folded into the Stage 1 benchmark requirement above rather than left as a standalone question.
