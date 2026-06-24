# Critique of `archive_rep()` Design Plan

## Overall Assessment

The current plan is much stronger than the original draft. It now keeps archive
configuration out of `initialize_genome()`, avoids storing archive metadata on
the `pop` object, separates `store_and_reset`, `store_once`, and `reset_only`
tables clearly, and documents that `tidybreed.rep` is only used to stamp
archived rows. DuckDB `ATTACH` remains the right implementation primitive for
copying rows without pulling data through R.

The remaining risks are mostly implementation-contract issues: archive schema
creation, schema evolution, transactional behavior, missing-table rules, and
the HPC workflow. These should be resolved before implementation because
`archive_rep()` is destructive: after success it deletes rows from the working
database.

## High-Priority Concerns

### 1. Archive schema creation is still unresolved

The plan's SQL sketch uses:

```sql
INSERT INTO archive.ind_meta SELECT *, 1 AS rep FROM ind_meta;
```

That only works if `archive.ind_meta` already exists with exactly the right
columns in exactly the right order. The open question about whether tables are
created lazily should be answered before implementation.

Recommended behavior:

- On first archive of a `store_and_reset` table, create the archive table from
  the current source schema plus a `rep` column.
- On first archive of a `store_once` table, create the archive table from the
  current source schema without `rep`.
- Insert with explicit column lists, never `SELECT *`.
- Do not copy working-table primary key constraints into archive tables unless
  the archive key includes `rep`.

Example shape:

```sql
INSERT INTO archive.ind_meta (id_ind, id_parent_1, id_parent_2, line_name, sex, rep)
SELECT id_ind, id_parent_1, id_parent_2, line_name, sex, 1 AS rep
FROM main.ind_meta;
```

### 2. Schema evolution needs an explicit policy

Tidybreed supports user-added columns through table mutation helpers and `...`
arguments. A later replicate may therefore have a source table schema that does
not exactly match the archive table created by the first replicate.

Recommended policy:

- Before any append, compare source and archive schemas.
- If the source gained columns, either add compatible nullable columns to the
  archive or error clearly. Pick one policy and document it.
- If a shared column has an incompatible type, error before modifying either
  database.
- If the archive has columns missing from the source, insert `NULL` only for
  columns that are nullable and not required by the archive contract.

Adding nullable columns is more user-friendly, but a strict error is simpler
and safer for the first implementation.

### 3. Copy/delete phases need failure semantics

The plan says to attach, insert, detach, then delete. If a failure happens
midway, the archive and working database can diverge. For example, rows may be
archived but not deleted, or some tables may be archived while later tables
fail.

Minimum safeguards:

- Validate `pop`, `rep`, archive path, table names, duplicate category
  assignment, source table existence, schemas, and archive writability before
  any insert or delete.
- Attach the archive once for the whole operation.
- Insert all archive tables before deleting any working rows.
- Delete only after all archive writes have succeeded.
- Increment `tidybreed.rep` only after all writes and deletes succeed.
- Use `on.exit()` cleanup for `DETACH`, without hiding the original error.

If DuckDB cannot guarantee a single transaction across the attached databases,
the plan should document the retry/failure model explicitly.

### 4. The HPC example still needs a safer default story

The plan now notes that `archive_rep()` primarily targets in-process loops and
that HPC users may need a separate post-processing merge. However, the HPC
example still shows every array job writing to the same
`scenario_1_all_reps.duckdb` archive.

That is likely to mislead users. DuckDB is not a high-concurrency write
database, and many array jobs appending to one archive file may fail or
serialize poorly.

Recommended documentation:

- In-process loops: one working DB, one archive DB, repeated `archive_rep()`
  calls.
- HPC arrays: each job writes a rep-local archive file, for example
  `scenario_1_rep_007.duckdb`.
- A later merge step combines rep-local archives into the final all-reps
  archive.

If same-file HPC writes are intentionally unsupported, say that directly.

### 5. `rep` option behavior should be precise

The current plan says `archive_rep()` always increments
`options(tidybreed.rep = rep + 1L)`. That is convenient for loops, but it also
means an explicit `rep = 10L` call mutates process-global state.

Recommended mitigation:

- Validate `rep` is a non-missing scalar integerish value.
- Decide whether explicit `rep` values should increment the option.
- Consider incrementing only when `rep` was taken from
  `getOption("tidybreed.rep")`.
- Add tests proving failed archives do not increment the option.

If the design intentionally keeps unconditional incrementing, document that
explicitly.

### 6. Missing-table behavior needs explicit rules

Not every population will contain every default table at every stage. Some
tables are created only after traits, phenotypes, or indexes are defined.
Custom tables may also exist.

Recommended behavior:

- Validate requested table names with the package's SQL identifier helper.
- Error for requested tables that do not exist unless an explicit
  `ignore_missing = TRUE` argument is added.
- Keep custom table support out of scope unless the reset/archive key semantics
  are documented.
- Error if the same table appears in more than one of `store_and_reset`,
  `store_once`, and `reset_only`.

## Medium-Priority Concerns

### Archive keys should be logical, not copied constraints

The plan says sequence counters should not be reset. That intent is sound, but
some IDs may still repeat after deleting rows because table-specific ID helpers
often derive the next value from the current table contents. After a table is
emptied, the next replicate may start again at `1`.

The archive should therefore treat `(rep, id_ind, ...)` or `(rep, id_*)` as the
logical identity for per-rep tables. It should not preserve source primary key
constraints that are unique only inside the working database.

### `id_ind` collisions across reps are expected

Repeated founder IDs such as `A_1`, `A_2`, etc. are normal across replicates.
The plan mentions this, but it should be promoted to a hard archive invariant:
per-individual archive queries must identify rows by `(rep, id_ind)`, not by
`id_ind` alone.

### Default table categories should be audited against dependencies

The defaults are sensible for common workflows:

- `store_and_reset`: individual-level results
- `store_once`: static interpretive metadata
- `reset_only`: large genomic state

Before implementation, audit every system table with an `id_ind`, individual
level, source table, random-effect level, or generated-individual dependency.
Any table that references reset individual rows should either be archived with
`rep` or reset alongside them.

The dynamic-genome note is useful. It should also mention dynamic trait,
phenotype, or index definitions: if those change by replicate, users must move
the relevant metadata/effect tables from `store_once` to `store_and_reset`.

### Archive provenance would make analysis and merging safer

A small archive metadata table would help downstream workflows and future
merge tooling. Useful fields include:

- tidybreed package version
- archive schema version
- archive creation timestamp
- population name or working DB path
- archived reps
- scenario label, if supplied

This is not required for a minimal `archive_rep()`, but it would reduce
ambiguity once archives are shared or merged.

### Path resolution should be tested

The current option design is straightforward:

1. explicit `archive_path` argument as a full path
2. `file.path(tidybreed.archive_path, tidybreed.db_name_archive)`
3. `db_name_archive = NULL` means reset without archive writes

The main remaining issue is testing and documenting relative paths. Because
`tidybreed.archive_path` defaults to `"."`, relative archives depend on the R
working directory at call time. That may be acceptable, but it should be
intentional and covered by tests.

## Implementation Recommendations

1. Add an internal archive table creator/synchronizer.

```r
ensure_archive_table(conn, table, archive_alias, add_rep)
```

It should inspect source/archive columns, create missing archive tables, and
reject unsupported schema drift.

2. Build all inserts from explicit quoted identifiers.

Avoid `SELECT *`. Quote table and column names through the package's existing
SQL helper style.

3. Do a validation pass before writes.

Check all arguments, table category overlaps, source table availability,
archive path resolution, schema compatibility, and archive writability before
copying or deleting rows.

4. Return `pop` for pipeline compatibility.

`archive_rep()` should return `pop`, probably invisibly, so this works:

```r
pop <- pop |> archive_rep()
```

If a summary is useful, attach it as an attribute, message it, or provide a
separate dry-run/verbose mode without breaking the pipeline contract.

5. Consider a dry-run or summary mode.

Because the function deletes working rows, a summary would be useful during
debugging:

```r
archive_rep(pop, dry_run = TRUE)
```

This can be deferred, but the implementation should at least have internal row
counts for tests.

## Tests To Add

- First call creates archive tables with `rep` columns for `store_and_reset`.
- First call creates `store_once` tables without `rep`.
- Second call appends per-rep rows without overwriting prior reps.
- `store_once` tables are copied once and skipped on later calls.
- `store_and_reset` and `reset_only` working tables are empty after success.
- `tidybreed.rep` increments only after a successful operation.
- Explicit `rep` behavior matches the documented option-increment policy.
- User-added columns are either preserved or rejected according to the chosen
  schema-evolution policy.
- Missing requested tables follow the documented missing-table policy.
- Duplicate table assignment across categories errors before writes.
- `db_name_archive = NULL` skips archive writes but still performs the reset
  phases.
- Relative and explicit archive path resolution behave as documented.
- Failed archive writes do not delete working rows.
- Failed deletes do not increment `tidybreed.rep`.
- Archive queries use `(rep, id_ind)` examples for per-individual tables.

## Suggested Plan Revisions

The plan is close to implementable. Before coding, resolve these choices:

1. Define lazy archive table creation for `store_and_reset` and `store_once`.
2. Choose strict-error versus nullable-column-add schema evolution.
3. Clarify transaction and retry behavior across the working and archive DBs.
4. Revise the HPC example to use per-job archives plus a merge step, or state
   clearly that concurrent writes to one archive are unsupported.
5. Decide whether explicit `rep` calls increment `tidybreed.rep`.
6. Specify missing-table behavior.
7. Decide whether archive provenance metadata is in scope for the first pass.

## Bottom Line

The revised design has addressed the earlier architectural concerns around the
deprecated initialization path and storing workflow metadata on `pop`. The
remaining work is to tighten destructive-operation guarantees and archive
schema contracts so `archive_rep()` is predictable, testable, and safe to use
in long replicate runs.
