# Follow-ups deferred from the v0.64.0 schema/rename work

Status: **reviewed and verified against the code (2026-09-04); ready to
implement.** Originally three gaps found while implementing
`plans/update_schema_print.md`, all deliberately left out of v0.64.0 because
each needed a decision the schema/print work did not settle.

The review confirmed the original diagnosis, **answered the one open question**,
and found **two additional defects** plus **three errors in the original plan**.
Revised order of value: **A > 1 > 3 > 2**.

| # | Gap | Blast radius | Suggested fix |
|---|-----|--------------|---------------|
| A | `phenotype_random_effects` never reset between replicates | **Correctness bug** — replicates are not independent | Add to `store_and_reset` |
| 1 | `TABLE_ROW_KEYS` missing entries **+ NULL-unsafe delete join** | `remove_rows()` hard-errors on 5 tables; silently deletes 0 rows on 2 more | Register the real keys; `IS NOT DISTINCT FROM` in the join |
| 2 | 4 tables have `CREATE TABLE` DDL in 2–3 files | Silent drift risk only | Single owner: delete the `ensure_*()` copies |
| 3 | `archive_replicate()` table lists are incomplete | Config tables not archived with the run | Extend the literal lists + completeness test |

---

## A. `phenotype_random_effects` carries across replicates (NEW — correctness bug)

### What is wrong

This was item 3's "check first" question. **The answer is: draws are reused.**

`R/phenotype_helpers.R:208-250` reads existing draws for the
`(phenotype_name, effect_name)` pair, computes
`new_lvls <- setdiff(unique_lvls, names(existing_map))`, and samples **only** for
levels not already stored:

```r
existing_draws <- DBI::dbGetQuery(pop$db_conn, "SELECT level, draw_value FROM ...")
new_lvls <- setdiff(unique_lvls, names(existing_map))
if (length(new_lvls) > 0) { ... rnorm/rgamma/runif ... }
per_ind <- unname(existing_map[as.character(group)])
```

That reuse is **correct within a replicate** — every animal in HYS level
`2020_Iowa` must receive the same shift across multiple `add_phenotype()` calls,
which is the entire point of a random effect. It is **wrong across replicates**,
and nothing clears the table: `R/archive_replicate.R:230` deletes only
`c(store_and_reset, reset_only)`, and `phenotype_random_effects` appears in
neither list (nor in `store_once`).

So replicate 2 inherits replicate 1's draws for every level name that repeats —
which, for `sex`, `line_name`, HYS, litter and pen grouping columns, is most of
them. The random effects are not re-drawn and the replicates are not
independent. This is a genuine correctness bug in the simulation, not an
archiving omission.

### Suggested fix

Add `"phenotype_random_effects"` to the `store_and_reset` default in
`archive_replicate()`. `store_and_reset` is exactly right: the draws are
per-replicate *output* worth keeping (stamped with `replicate`, so a run can be
audited or reproduced), and the working-DB rows are deleted afterwards, so the
next replicate re-draws from scratch.

`store_once` would be wrong — it archives once and leaves the working rows in
place, which is the current broken behaviour with extra steps.

`.ensure_archive_table(conn, tbl, add_replicate = TRUE)` adds the `replicate`
column to the **archive** copy only, so no working-DB schema change is needed and
the existing collision guard (`archive_replicate.R:145`) already covers it.

### Test to add

Two replicates with the same grouping levels must produce **different** draws:

```r
test_that("random effect draws are re-drawn after archive_replicate()", {
  # ... replicate 1, capture phenotype_random_effects
  pop <- archive_replicate(pop, archive_path = tempfile(fileext = ".duckdb"))
  expect_equal(nrow(collect(get_table(pop, "phenotype_random_effects"))), 0L)
  # ... replicate 2 with the same levels; draws must differ
})
```

### Scope

`R/archive_replicate.R` (one string + a Details note),
`tests/testthat/test-archive_replicate.R`. Small — release on its own, ahead of
the tidy-ups.

---

## 1. `TABLE_ROW_KEYS` gaps **and** a NULL-unsafe delete join

### What is wrong

**1a. Missing registry entries (as originally described — confirmed).**
v0.64.0 closed this class of gap for `TABLE_RESERVED_COLS` (every one of the 24
`SYSTEM_TABLES` now has an entry). The other two registries in `R/sql_utils.R`
were left as they were:

```
no TABLE_ROW_KEYS:      _schema_meta, founder_haplotypes,
                        phenotype_random_effects, phenotype_meta,
                        phenotype_components
no TABLE_PRIMARY_KEYS:  _schema_meta, ind_haplotype, ind_genotype,
                        chr_inheritance, chr_recombination, founder_haplotypes,
                        phenotype_effects, phenotype_random_effects,
                        phenotype_meta, phenotype_components
```

`TABLE_ROW_KEYS` is the one that bites. `remove_rows()` in single-table mode
(`R/remove_rows.R:216`) does:

```r
key_cols <- TABLE_ROW_KEYS[[table_name]]
if (is.null(key_cols)) {
  stop("Cannot delete from '", table_name,
       "': table is not registered in TABLE_ROW_KEYS. ...")
}
```

Verified live:

```
> pop |> get_table("phenotype_meta") |> filter(phenotype_name == "ADG") |> remove_rows()
Error : Cannot delete from 'phenotype_meta': table is not registered in
TABLE_ROW_KEYS. ...
```

`phenotype_meta` and `phenotype_components` both have a perfectly good integer
PK, so this is a pure omission. Same shape as the v0.63.1 bug
(`TABLE_ROW_KEYS$trait_effects` listed a renamed column, so every delete
aborted), except it fails at "not registered" rather than "missing key column".

**1b. `delete_exact_rows()` cannot delete rows with NULL key values (NEW).**
`R/remove_rows.R:6-46` builds its join as:

```r
join_sql <- paste(paste0("t.", key_cols, " = f.", key_cols), collapse = " AND ")
```

`NULL = NULL` is `NULL`, not `TRUE`, so any row whose key column is `NULL` never
matches. Two **already-registered** tables have nullable key columns —
`chr_inheritance` (`offspring_sex`, `line_name`) and `chr_recombination`
(`parent_sex`, `line_name`) — and `define_genome()` seeds exactly those rows with
`NULL` in both. Verified live:

```
> pop |> get_table("chr_inheritance") |> filter(chr_name == "1") |> remove_rows()
Deleted 0 rows from `chr_inheritance`      # <- and the row is still there
```

This is worse than 1a: it reports **success** while doing nothing. It also
blocks registering `founder_haplotypes`, whose `line_name` is `NULL` for the
shared pool.

`define_chromosome()` already uses the right idiom (`IS NOT DISTINCT FROM`) for
its delete-then-insert upsert, so the fix is to match it.

`TABLE_PRIMARY_KEYS` is lower-stakes: it drives vector updates and filtered
updates in `mutate_table()`. Several tables in that missing list genuinely have
no single-column PK (`ind_haplotype`, `chr_inheritance`, …), so a missing entry
there is often *correct*.

### Suggested fix

**1a. Register the four tables where row deletion is well-defined.** In
`R/sql_utils.R`:

```r
TABLE_ROW_KEYS <- list(
  ...
  phenotype_meta           = "id_phenotype_meta",
  phenotype_components     = "id_phenotype_comp",
  phenotype_random_effects = c("phenotype_name", "effect_name", "level"),
  founder_haplotypes       = c("line_name", "haplotype_id", "locus_name"),
  ...
)

TABLE_PRIMARY_KEYS <- list(
  ...
  phenotype_meta       = "id_phenotype_meta",
  phenotype_components = "id_phenotype_comp",
  ...
)
```

`phenotype_random_effects` and `founder_haplotypes` get a row key but no primary
key: their PKs are composite, which is exactly what `TABLE_ROW_KEYS` is for and
exactly what `TABLE_PRIMARY_KEYS` is not.

**Correction to the original plan: do not blanket-refuse `founder_haplotypes`.**
The original 1b proposed a `TABLE_NO_ROW_DELETE` entry for it, on the grounds
that deleting one row leaves a non-rectangular pool. That reasoning holds for a
single `(haplotype_id, locus_name)` row, but not for the operation a user
actually wants: dropping **one line's entire pool** (`filter(line_name == "B")`)
is well-defined, leaves the remaining lines rectangular, and is a plausible
crossbreeding-setup operation. A flat refusal forecloses it to prevent a misuse
that a `nrow`-level check cannot distinguish anyway. Register the composite key
and let the operation work.

**1b. Make the delete join NULL-safe.** In `delete_exact_rows()`:

```r
join_sql <- paste(
  paste0("t.", key_cols, " IS NOT DISTINCT FROM f.", key_cols),
  collapse = " AND "
)
```

`IS NOT DISTINCT FROM` is DuckDB's NULL-safe equality and is already the idiom
used in `define_chromosome()`. It behaves identically to `=` for non-NULL values,
so no currently-working delete changes behaviour.

**1c. Keep `TABLE_NO_ROW_DELETE`, but only for `_schema_meta`.** An omission and
a decision still look identical from the outside, and the registry should say
which it is:

```r
#' Tables where single-row deletion is not a meaningful operation
#'
#' A missing TABLE_ROW_KEYS entry is ambiguous: it may be an oversight. These
#' tables are listed so the refusal is a decision with a reason attached.
#' @keywords internal
TABLE_NO_ROW_DELETE <- c(
  `_schema_meta` = paste(
    "schema descriptions are package-managed and rebuilt by open_pop();",
    "editing them by hand would be silently overwritten."
  )
)
```

Branch on it in `remove_rows()` before the generic "not registered" error, so the
user gets a reason and an alternative rather than a dead end.

**1d. Extend the existing registry audit.** `tests/testthat/test-schema-registries.R`
already asserts that every *listed* column exists. Add the converse, so a new
table cannot be added without answering the question:

```r
test_that("every system table either has a row key or is explicitly excluded", {
  expect_setequal(
    c(names(TABLE_ROW_KEYS), names(TABLE_NO_ROW_DELETE)),
    SYSTEM_TABLES
  )
})
```

Plus a regression test that a NULL-keyed row actually deletes.

This is the same visible-degradation principle `.schema_table_order()` uses: the
next person to add a table is forced to make the call, and forgetting fails a
test instead of surfacing months later as "why can't I delete from this table".

### Scope

`R/sql_utils.R`, `R/remove_rows.R`, `tests/testthat/test-schema-registries.R`,
`tests/testthat/test-remove_rows.R`, plus a `remove_rows()` roxygen note.

---

## 2. Duplicate `CREATE TABLE` DDL (four tables, up to three sites each)

### What is wrong

**Correction to the original plan: the duplication is wider than three tables.**
The original count missed `R/define_effect_cov_matrix.R`, which contains
`ensure_trait_var_comp()` and `ensure_phenotype_var_comp()`. Full enumeration of
every `CREATE TABLE` across the package:

| Table | DDL sites |
|---|---|
| `phenotype_var_comp` | `open_pop.R:341`, `define_trait.R:322`, `define_effect_cov_matrix.R:217` — **three** |
| `trait_var_comp` | `open_pop.R:276`, `define_effect_cov_matrix.R:195` — two |
| `phenotype_meta` | `open_pop.R:298`, `define_trait.R:279` — two |
| `phenotype_components` | `open_pop.R:319`, `define_trait.R:300` — two |

Single-site tables: `_schema_meta`, `genome_effects`, `ind_meta` (`open_pop.R`);
`ind_ebv`, `ind_index`, `ind_phenotype`, `ind_tbv`, `ind_true_index`,
`index_meta`, `phenotype_effects`, `phenotype_random_effects`, `trait_meta`
(`define_trait.R`); the genome and per-chromosome tables (`define_genome.R`).
`founder_haplotypes` is created implicitly by `dbWriteTable()` in
`R/founder_haplotype_helpers.R:274`, which is a fourth pattern again.

All copies were diffed column-for-column during the review and **currently agree
exactly**, so nothing is broken today. The problem is four tables with two or
three definitions each and nothing enforcing that they stay in sync. Add a column
to one copy and the behaviour depends on which function ran first — exactly the
kind of drift the v0.63.1 registry bug came from.

The split is also not principled: there is no property that makes
`phenotype_meta` belong in both files while `phenotype_effects` belongs in only
one. It reads like accretion, not design.

### Suggested fix

**Give these tables one owner: `open_pop()`.** It already creates
`phenotype_meta`, `phenotype_var_comp`, `trait_var_comp` and the rest of the
observation-layer tables unconditionally, so `phenotype_components` sitting there
is consistent; the `ensure_*()` helpers are the odd ones out.

1. Delete the `phenotype_components`, `phenotype_meta` and `phenotype_var_comp`
   entries from the `ddl` list in `R/define_trait.R`, and delete
   `ensure_trait_var_comp()` / `ensure_phenotype_var_comp()` from
   `R/define_effect_cov_matrix.R` (replacing their call sites with a plain
   existence assertion, or nothing at all if `open_pop()` guarantees the table).
2. Watch `pop$tables`: `ensure_trait_tables()` updates it with `names(ddl)`, and
   the `ensure_*_var_comp()` helpers append their own name, so dropping them also
   drops those tables from that union. They must still reach `pop$tables` via
   `open_pop()`'s own `tables_created` vector, or `schema()` will stop listing
   them.
3. Confirm no path reaches those tables without going through `open_pop()` — it
   should not, since `open_pop()` is the only entry point that creates a
   database, but `restore_pop()` is worth a look.

**Then make the drift impossible rather than merely absent.** The better version,
if it is worth the effort: move every `CREATE TABLE` string into one internal
`.tidybreed_ddl()` list keyed by table name, have `open_pop()`,
`ensure_trait_tables()` and `define_genome()` all create from it, and test that
its names equal `SYSTEM_TABLES`. Each call site keeps deciding *which* tables it
creates and when; only the column definitions are single-sourced.

That would also let the registry audit get much stronger: with the DDL as data,
`TABLE_RESERVED_COLS` could be checked against it directly instead of against a
live database, and `test-schema-registries.R` would stop needing to build a
population just to list columns.

(The original plan sketched a `grep`-over-source-files test as the cheap option.
It is file-path-bound, needs one copy per duplicated table, and would not have
caught the `define_effect_cov_matrix.R` sites it was written against — skip it in
favour of the DDL-as-data version or nothing.)

### Scope

`R/define_trait.R`, `R/define_effect_cov_matrix.R`, `R/open_pop.R`, one test.
Small if you take the delete-only option; a day if you take the
single-DDL-registry option, which also touches `R/define_genome.R` and
`R/founder_haplotype_helpers.R`.

---

## 3. `archive_replicate()`'s table lists are incomplete

### What is wrong

`R/archive_replicate.R:101` hard-codes three lists of table names:

```r
store_and_reset = c("ind_meta", "ind_phenotype", "ind_tbv",
                    "ind_ebv", "ind_index", "ind_true_index"),
store_once      = c("genome_meta", "genome_effects",
                    "trait_meta", "phenotype_effects", "trait_var_comp",
                    "phenotype_meta", "phenotype_components",
                    "phenotype_var_comp", "index_meta"),
reset_only      = c("ind_haplotype", "ind_genotype", "ind_crossover")
```

That accounts for 18 of the 24 `SYSTEM_TABLES`. Missing entirely:

| Table | Should be | Consequence of the omission |
|---|---|---|
| `phenotype_random_effects` | `store_and_reset` | **See item A** — promoted to its own fix; replicates are not independent |
| `genome_map` | `store_once` | The genetic map is not archived with the run that used it |
| `chr_inheritance` | `store_once` | Per-chromosome inheritance rules not archived |
| `chr_recombination` | `store_once` | Per-chromosome recombination rules not archived |
| `founder_haplotypes` | `store_once` | The founder pool a replicate was drawn from is not recorded |
| `_schema_meta` | neither | Correctly excluded — it is system metadata |

With item A landed separately, what remains here is a pure archiving-completeness
gap: four configuration tables that describe how a replicate was generated are
not copied into the archive, so an archived run is not self-describing.

### Suggested fix

**Correction to the original plan: do not derive the lists from the display
groups.** The original proposal mapped `.schema_table_order()` groups to archive
categories with a single `phenotype_random_effects` override. That mapping is
wrong: it places **Individuals** → `reset_only`, but `ind_meta` is in the
Individuals group and is currently `store_and_reset` — and deliberately so; the
roxygen at `archive_replicate.R:45` documents that archiving `ind_meta` is what
makes `(replicate, id_ind)` a valid composite key. Deriving from groups would
silently stop archiving pedigree.

That is two overrides out of eight groups, which is the tell: display grouping
answers *where a table sits in the pipeline*, archiving answers *whether its
contents are per-replicate output*. Those are different axes and they disagree on
`ind_meta` and `phenotype_random_effects` — the two most important tables in the
mapping. A derivation that needs an override for its own headline cases is not a
derivation.

**Instead: keep the three literal lists, extend them, and add a completeness
test.** The lists are the right shape; they were merely never updated when
`genome_map`, `chr_inheritance`, `chr_recombination` and `founder_haplotypes`
were added.

```r
store_once = c("genome_meta", "genome_map", "genome_effects",
               "chr_inheritance", "chr_recombination", "founder_haplotypes",
               "trait_meta", "phenotype_effects", "trait_var_comp",
               "phenotype_meta", "phenotype_components",
               "phenotype_var_comp", "index_meta")
```

**Test to add:** every `SYSTEM_TABLES` entry appears in exactly one of the three
default lists, or in an explicit `ARCHIVE_EXCLUDED` vector (`_schema_meta`) —
the same visible-degradation pattern as item 1. This gives the "add a table"
checklist its fourth item, but makes forgetting it a test failure rather than a
silently incomplete archive.

### Scope

`R/archive_replicate.R`, `tests/testthat/test-archive_replicate.R`. Small once
item A is separated out.

---

## Suggested sequencing

1. **Item A** — the `phenotype_random_effects` reset. Correctness bug, own
   commit, own release.
2. **Item 1**, as its own commit. Self-contained; removes a hard error and a
   silent no-op that users can hit today.
3. **Item 3**, once 1 is in — it is a short list edit plus the completeness test.
4. **Item 2** whenever convenient; it protects against a future problem rather
   than fixing a present one.

Version bump + `NEWS.md` entry per `CLAUDE.md` before each commit.
