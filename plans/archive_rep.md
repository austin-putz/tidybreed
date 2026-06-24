# archive_replicate() — Design Plan

## Problem

Users running many replicates (HPC job arrays or in-process loops) need a way
to accumulate results across replicates without ballooning the working DB or
re-querying with `WHERE replicate = ?` at every step. The working DB stays
lean; a separate archive DB accumulates all replicates.

---

## Global Options

### `tidybreed.replicate`

The package sets `options(tidybreed.replicate = 1L)` in `.onLoad()` as the
default starting replicate. Its only role is to tell `archive_replicate()` what
value to stamp on archived rows — no hidden filtering anywhere else in the
package.

**Loop users** never need to set it; `archive_replicate()` manages the counter,
incrementing `tidybreed.replicate + 1L` after each call.

**HPC/sbatch users** (one replicate per job) override it once at the top of the
script:

```r
options(tidybreed.replicate = as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID")))
```

This way users can launch all replicates at one time and at the end, save to a
common `.duckdb` file in a permanent folder specified as an option (shown below).

### `tidybreed.archive_path`

The package sets `options(tidybreed.archive_path = NULL)` in `.onLoad()`.
The directory where the archive DB will be written. HPC users typically set
this to their scratch or output directory once per session.

```r
options(tidybreed.archive_path = "~/user/my_simulations/results")
```

### `tidybreed.db_name_archive`

The package sets `options(tidybreed.db_name_archive = NULL)` in `.onLoad()`.
The filename of the archive DB. Combined with `tidybreed.archive_path` to form
the full path.

```r
options(tidybreed.db_name_archive = "scenario_1_all_replicates.duckdb")
```

**Resolution order** inside `archive_replicate()`:

1. Explicit `archive_path` argument — used as-is (full path), bypasses both options
2. `file.path(getOption("tidybreed.archive_path"), getOption("tidybreed.db_name_archive"))` — if `tidybreed.archive_path` is explicitly set
3. `file.path(dirname(pop$db_path), getOption("tidybreed.db_name_archive"))` — default: same directory as the working DB
4. `db_name_archive` is `NULL` — skip all archiving; only the reset phases run
   (`reset_only` deletes and `store_and_reset` tables are deleted without being
   copied to an archive)

---

## Function Signature

```r
archive_replicate(
  pop,
  replicate       = getOption("tidybreed.replicate"),
  archive_path    = NULL,
  store_and_reset = c("ind_meta", "ind_phenotype", "ind_tbv",
                      "ind_ebv", "ind_index", "ind_true_index"),
  store_once      = c("genome_meta", "genome_effects",
                      "trait_meta", "trait_effects", "trait_var_comp",
                      "phenotype_meta", "phenotype_components", "phenotype_var_comp",
                      "index_meta"),
  reset_only      = c("genome_haplotype", "genome_genotype")
)
```

Four-phase operation:

1. **store_and_reset** — ATTACH archive DB, INSERT each listed table with the
   current `replicate` value as a new column, DETACH. Then DELETE all rows from
   those tables in the working DB.
2. **store_once** — copy each listed table to the archive DB *without* a
   `replicate` column, but only if the table does not already exist in the
   archive. Detected by checking the archive schema — skipped on all subsequent
   calls. Working DB rows are left intact.
3. **reset_only** — DELETE all rows from each listed table in the working DB.
   Nothing is written to the archive. Intended for large genomic tables
   (`genome_haplotype`, `genome_genotype`) that are too large to archive by
   default and not needed for most downstream analyses.
4. **Increment** — `options(tidybreed.replicate = replicate + 1L)` so the next
   call in the loop automatically uses the next replicate number.

HPC users calling this once per job are unaffected by the increment — the
option bumps by 1 with no subsequent iteration to observe it.

Any table not listed in any argument is never touched.

---

## Table Categories (defaults)

### `store_and_reset` — archive with replicate, then delete

`ind_meta`, `ind_phenotype`, `ind_tbv`, `ind_ebv`, `ind_index`, `ind_true_index`

### `store_once` — copy to archive on first call only, leave working DB intact

On the first `archive_replicate()` call each of these is copied to the archive
DB without a `replicate` column. On subsequent calls the function detects the
table already exists in the archive and skips it. This makes the archive
self-contained for post-hoc analysis — trait names, phenotype definitions,
index weights, and variance components are all present alongside the
per-replicate results.

Default tables:
`genome_meta`, `genome_effects`,
`trait_meta`, `trait_effects`, `trait_var_comp`,
`phenotype_meta`, `phenotype_components`, `phenotype_var_comp`,
`index_meta`

Users with a **dynamic genome** (allele frequencies, effects, or loci that
change between replicates) should move `genome_meta` and/or `genome_effects` to
`store_and_reset` so each replicate's genetic architecture is preserved with its
replicate number.

### `reset_only` — delete rows, never archive

`genome_haplotype`, `genome_genotype`

Large genomic tables not needed for most downstream analyses. Users who need
per-replicate genomic data (e.g. to track inbreeding or run genomic prediction
post-hoc) move these to `store_and_reset`.

### Never touched (not listed in any argument)

`founder_haplotypes` — the base haplotype sampling pool. Large and not needed
for interpreting results. Users who need it can add it to `store_once`.

---

## The `replicate` Column

The column added to every `store_and_reset` table in the archive is named
`replicate`. This is unambiguous and self-documenting in post-hoc queries:
`filter(replicate == 5L)` reads clearly without any context.

### Collision guard

Before running any SQL, `archive_replicate()` checks whether `replicate`
already exists as a column in any table listed in `store_and_reset`. If it
does, the function stops immediately with a clear error:

> `"Table 'ind_meta' already contains a column named 'replicate'. Rename it before calling archive_replicate()."`

`replicate` is a reserved column name for all tables that may appear in
`store_and_reset` and should be added to the package's `TABLE_RESERVED_COLS`
guard so `mutate_table()` also blocks it. Users who need a replicate-like
column on `ind_meta` should name it something else (e.g. `sub_replicate`,
`within_replicate`).

---

## Implementation: DuckDB ATTACH

Use DuckDB's native `ATTACH` to avoid any R round-trip:

```sql
ATTACH 'results.duckdb' AS archive;
INSERT INTO archive.ind_meta SELECT *, 1 AS replicate FROM ind_meta;
-- repeat for each store_and_reset table
DETACH archive;
DELETE FROM ind_meta;
-- repeat DELETE for each store_and_reset table
```

No data is pulled into R memory. The cross-DB copy is a single SQL operation.

---

## Archive Path — Storage

The intended pattern is a user-created permanent results folder that accumulates
all archive DBs across scenarios and research questions, separate from the
temporary working DB directory where large genotype files, external software
outputs, and other disposable artifacts live.

```text
results/
  research_question_1/
    scenario_1_all_replicates.duckdb
    scenario_2_all_replicates.duckdb
  research_question_2/
    scenario_1_all_replicates.duckdb
```

Users set `tidybreed.archive_path` once to point at their results folder for
that research question, and `tidybreed.db_name_archive` to name the scenario.
The archive DB is permanent; the working DB directory is temporary and may be
wiped between runs to recover disk space.

Do not store `archive_path` on the `pop` object or pass it to
`initialize_genome()` — it is a workflow concern, not a genome concern. If a
convenience shortcut is added later, `open_pop()` would be the appropriate
place.

---

## Individual IDs and the Composite Key Problem

### Sequential loop (one working DB, replicates run one after another)

Auto-increment counters keep accumulating across replicates — replicate 1 ends
at `Line_500`, replicate 2 starts at `Line_501`. `id_ind` is globally unique in
the archive on its own. The `replicate` column is still present and useful for
filtering but is not strictly required for row identification.

Do **not** reset sequence counters between replicates. Let `id_phenotype`,
`id_tbv`, etc. continue accumulating.

### HPC parallel (separate working DB per job, replicates run simultaneously)

Each job starts from a fresh DB. Every job generates `Line_1`, `Line_2`, etc.
When all jobs write to the same archive DB, `id_ind` collides across replicates
— replicate 1 and replicate 5 both have a `Line_1`.

**Changing IDs is not the solution.** Prefixing or offsetting `id_ind` would
require transforming every FK column (`id_parent_1`, `id_parent_2`, and all
join columns) across all tables simultaneously — fragile and error-prone.

**The correct design: `(replicate, id_ind)` is the composite identifier in the
archive.** `id_ind` is only unique *within* a replicate. Across replicates,
`replicate` is always required to scope the result. This is standard practice
for replicated simulation studies.

**Rule for archive queries: always include `replicate` in joins and filters.**

```r
# correct
get_table(archive, "ind_phenotype") |>
  filter(replicate == 3L) |>
  left_join(get_table(archive, "ind_meta") |> filter(replicate == 3L), by = "id_ind")

# wrong — cross-joins replicates silently
get_table(archive, "ind_phenotype") |>
  left_join(get_table(archive, "ind_meta"), by = "id_ind")
```

---

## Intended Loop Pattern

```r
options(tidybreed.db_name_archive = "scenario_1_all_replicates.duckdb")
# options(tidybreed.replicate) starts at 1L automatically

pop <- initialize_genome(...)
# define traits, effects, phenotypes once — these persist across replicates

for (i in 1:100) {
  pop <- pop |>
    add_founders(...) |>
    # ... full replicate simulation
    archive_replicate()   # reads replicate and archive path from global options
}

# post-hoc analysis
archive <- open_archive("scenario_1_all_replicates.duckdb")
get_table(archive, "ind_phenotype") |> filter(replicate == 5L) |> collect()
```

## HPC Pattern (one replicate per job)

```r
options(tidybreed.replicate       = as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID")))
options(tidybreed.db_name_archive = "scenario_1_all_replicates.duckdb")

pop <- initialize_genome(...)
# define traits, effects, phenotypes...
pop <- pop |> add_founders(...) |> ...

archive_replicate(pop)   # stamps correct replicate, harmlessly increments option
```

---

## Open Questions

- Should `open_archive()` be a companion function that returns a read-only
  `tidybreed_pop`-like handle for querying the archive?
- Should `archive_replicate()` validate that all expected tables exist in the
  archive schema, or create them lazily on first call?
- For the HPC/sbatch case (one replicate per job, separate working DBs), users
  would use a separate post-processing step to merge results —
  `archive_replicate()` targets the in-process loop case primarily. Document
  both patterns.
