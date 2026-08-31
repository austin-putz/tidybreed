# Follow-ups deferred from the v0.64.0 schema/rename work

Status: **not started — for review.** Three gaps found while implementing
`plans/update_schema_print.md`. All three were deliberately left out of v0.64.0
because each needs a decision that the schema/print work did not settle, and
guessing would have buried a design choice inside a mechanical change.

They are independent of each other. Rough order of value: **1 > 3 > 2**.

| # | Gap | Blast radius | Suggested fix |
|---|-----|--------------|---------------|
| 1 | `TABLE_ROW_KEYS` / `TABLE_PRIMARY_KEYS` missing entries | `remove_rows()` hard-errors on 5 tables | Register the 3 that have real keys; make the other 2 error deliberately |
| 2 | 3 tables have `CREATE TABLE` DDL in two files | Silent drift risk only | Single owner: delete the `ensure_trait_tables()` copies |
| 3 | `archive_replicate()` table lists are incomplete | Replicates silently lose config tables | Derive from `.schema_table_order()` groups, not a hand list |

---

## 1. `TABLE_ROW_KEYS` / `TABLE_PRIMARY_KEYS` have missing entries

### What is wrong

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

So `remove_rows()` cannot delete a row from `phenotype_meta`,
`phenotype_components`, `founder_haplotypes`, or `phenotype_random_effects` —
even though the first two have a perfectly good integer PK. This is the same
shape as the v0.63.1 bug (`TABLE_ROW_KEYS$trait_effects` listed a column that had
been renamed, so every delete aborted), except it fails at "not registered"
rather than "missing key column".

`TABLE_PRIMARY_KEYS` is lower-stakes: it drives vector updates and filtered
updates in `mutate_table()`. Several tables in that missing list genuinely have
no single-column PK (`ind_haplotype`, `chr_inheritance`, …), so a missing entry
there is often *correct*.

### Why it was not fixed in v0.64.0

Two different questions are tangled together, and only one is mechanical:

- Tables with a real single-column PK (`phenotype_meta.id_phenotype_meta`,
  `phenotype_components.id_phenotype_comp`) are pure omissions — register them.
- Tables with only a composite logical key (`founder_haplotypes`,
  `phenotype_random_effects`) raise an actual design question: *should*
  `remove_rows()` delete from them at all? Deleting one row of a founder
  haplotype pool leaves a pool that is no longer rectangular, and every
  downstream sampler assumes it is. That is a decision, not a lookup.

### Suggested fix

**Split by whether deletion is meaningful, and make both answers explicit.**

**1a. Register the three tables where row deletion is well-defined.** In
`R/sql_utils.R`:

```r
TABLE_ROW_KEYS <- list(
  ...
  phenotype_meta       = "id_phenotype_meta",
  phenotype_components = "id_phenotype_comp",
  phenotype_random_effects = c("phenotype_name", "effect_name", "level"),
  ...
)

TABLE_PRIMARY_KEYS <- list(
  ...
  phenotype_meta       = "id_phenotype_meta",
  phenotype_components = "id_phenotype_comp",
  ...
)
```

`phenotype_random_effects` gets a row key but no primary key: its PK is the
composite `(phenotype_name, effect_name, level)`, which is exactly what
`TABLE_ROW_KEYS` is for and exactly what `TABLE_PRIMARY_KEYS` is not.

**1b. Refuse deliberately for `founder_haplotypes`, do not just omit it.** An
omission and a decision currently look identical from the outside. Add a small
registry next to the others:

```r
#' Tables where single-row deletion is not a meaningful operation
#'
#' A missing TABLE_ROW_KEYS entry is ambiguous: it may be an oversight. These
#' tables are listed so the refusal is a decision with a reason attached.
#' @keywords internal
TABLE_NO_ROW_DELETE <- c(
  founder_haplotypes =
    paste("the pool is rectangular (every haplotype_id has a row at every",
          "locus) and add_founders() assumes it. Drop and rebuild the pool",
          "with define_founder_haplotypes() instead."),
  "_schema_meta" =
    "descriptions are managed by define_schema_description()."
)
```

and branch on it in `remove_rows()` before the generic "not registered" error,
so the user gets the reason and the alternative rather than a dead end.

**1c. Extend the existing registry audit.** `tests/testthat/test-schema-registries.R`
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

This is the same visible-degradation principle `.schema_table_order()` uses: the
next person to add a table is forced to make the call, and forgetting fails a
test instead of surfacing months later as "why can't I delete from this table".

### Scope

`R/sql_utils.R`, `R/remove_rows.R`, `tests/testthat/test-schema-registries.R`,
plus a `remove_rows()` roxygen note. Small — an afternoon, most of it the
`founder_haplotypes` decision.

---

## 2. Three tables have duplicate `CREATE TABLE` DDL

### What is wrong

**Three** tables are defined twice, not one. Enumerating every `CREATE TABLE` in
the two files that create tables:

| Only in `R/open_pop.R` | Only in `R/define_trait.R` (`ensure_trait_tables()`) | **In both** |
|---|---|---|
| `_schema_meta`, `genome_effects`, `ind_meta`, `trait_var_comp` | `ind_ebv`, `ind_index`, `ind_phenotype`, `ind_tbv`, `ind_true_index`, `index_meta`, `phenotype_effects`, `phenotype_random_effects`, `trait_meta` | `phenotype_components`, `phenotype_meta`, `phenotype_var_comp` |

(`define_genome()` creates the remaining genome and per-chromosome tables;
`founder_haplotypes` is created implicitly by `dbWriteTable()` in
`R/founder_haplotype_helpers.R`, which is a fourth pattern again.)

Both copies are guarded by an existence check and the definitions currently agree
column-for-column, so nothing is broken today. The problem is three tables with
two definitions each and nothing enforcing that they stay in sync. Add a column
to one copy and the behaviour depends on which function ran first — exactly the
kind of drift the v0.63.1 registry bug came from.

The split is also not principled: there is no property that makes
`phenotype_meta` belong in both files while `phenotype_effects` belongs in only
one. It reads like accretion, not design.

### Why it was not fixed in v0.64.0

It is unrelated to the rename and to the printing work, and touching table
creation order is the sort of change that wants its own test run rather than
being buried in a 44-file diff.

### Suggested fix

**Give these tables one owner: `open_pop()`.** It already creates
`phenotype_meta`, `phenotype_var_comp` and the rest of the observation-layer
tables unconditionally, so `phenotype_components` sitting there is consistent;
`ensure_trait_tables()` is the odd one out.

1. Delete the `phenotype_components`, `phenotype_meta` and `phenotype_var_comp`
   entries from the `ddl` list in `R/define_trait.R`. Diff each pair first —
   they are believed identical, but confirm rather than assume.
2. Watch `pop$tables`: `ensure_trait_tables()` updates it with `names(ddl)`, so
   dropping the three entries also drops them from that union. They must still
   reach `pop$tables` via `open_pop()`'s own `tables_created` vector, or
   `schema()` will stop listing them.
3. Confirm no path reaches those tables without going through `open_pop()` — it
   should not, since `open_pop()` is the only entry point that creates a
   database, but `restore_pop()` is worth a look.

**Then make the drift impossible rather than merely absent.** A test that
compares the live column set of every table against a single declared source
would catch this class permanently. The cheapest version reuses what already
exists:

```r
test_that("no table is created by more than one DDL site", {
  ddl_sites <- c(
    length(grep("CREATE TABLE phenotype_meta", readLines("../../R/open_pop.R"))),
    length(grep("CREATE TABLE phenotype_meta", readLines("../../R/define_trait.R")))
  )
  expect_equal(sum(ddl_sites), 1L)
})
```

That is ugly, file-path-bound, and needs one copy per duplicated table. A better
version, if it is worth the effort: move every `CREATE TABLE` string into one
internal `.tidybreed_ddl()` list keyed by table name, have `open_pop()`,
`ensure_trait_tables()` and `define_genome()` all create from it, and test that
its names equal `SYSTEM_TABLES`. Each call site keeps deciding *which* tables it
creates and when; only the column definitions are single-sourced.

That would also let the registry audit get much stronger: with the DDL as data,
`TABLE_RESERVED_COLS` could be checked against it directly instead of against a
live database, and `test-schema-registries.R` would stop needing to build a
population just to list columns.

### Scope

`R/define_trait.R`, `R/open_pop.R`, one test. Small if you take the delete-only
option; a day if you take the single-DDL-registry option, which also touches
`R/define_genome.R` and `R/founder_haplotype_helpers.R`.

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

| Table | Should probably be | Consequence of the omission |
|---|---|---|
| `phenotype_random_effects` | `store_and_reset` | Sampled random-effect draws are neither archived nor cleared between replicates, so replicate 2 silently inherits replicate 1's draws |
| `genome_map` | `store_once` | The genetic map is not archived with the run that used it |
| `chr_inheritance` | `store_once` | Per-chromosome inheritance rules not archived |
| `chr_recombination` | `store_once` | Per-chromosome recombination rules not archived |
| `founder_haplotypes` | `store_once` (or `reset_only`) | The founder pool a replicate was drawn from is not recorded |
| `_schema_meta` | neither | Correctly excluded — it is system metadata |

`phenotype_random_effects` is the one that looks like a genuine correctness bug
rather than an archiving omission: those are *sampled values*, not configuration.
Carrying them across replicates means the random effects are not re-drawn, which
is probably not what a replicate is supposed to mean. This needs confirming
against how `define_effect_random()` and `add_phenotype()` treat existing draws
before it is called a bug.

### Why it was not fixed in v0.64.0

Every row of that table is a semantic question about what a replicate *is*, and
the rename only touched one name in one of the three lists. Answering
"should the founder pool be archived per replicate?" inside a rename commit would
have been the wrong place for it.

### Suggested fix

**Stop hand-maintaining the lists; derive them from the display groups.**
v0.64.0 established `.schema_table_order()` as the single place that knows every
table and what layer it belongs to, and the archiving categories line up with
those groups almost exactly:

| Display group | Archive category | Reasoning |
|---|---|---|
| Genome, Founders | `store_once` | Fixed across replicates of one design |
| Genetic model, Observation model, Selection | `store_once` | Model configuration |
| Individuals | `reset_only` | Regenerated each replicate; too large to archive |
| Results | `store_and_reset` | The per-replicate output — the point of archiving |
| System | excluded | |

The one place the mapping is not mechanical is `phenotype_random_effects`: it
lives in the **Observation model** group (it is keyed by phenotype and written by
`define_effect_random()`) but behaves like **Results** (it holds sampled values).
That is the real finding here — the display grouping and the archiving semantics
disagree about exactly one table, and that disagreement is worth resolving
explicitly rather than papering over.

Suggested shape:

```r
#' Archive category per display group, with per-table overrides
#' @keywords internal
.archive_category <- function() {
  by_group <- c(
    "Genome" = "store_once", "Founders" = "store_once",
    "Genetic model" = "store_once", "Observation model" = "store_once",
    "Selection" = "store_once",
    "Individuals" = "reset_only", "Results" = "store_and_reset",
    "System" = NA_character_
  )
  cat <- by_group[as.character(.schema_group_of(SYSTEM_TABLES)$table_group)]
  names(cat) <- SYSTEM_TABLES

  # Override: sampled draws, not configuration. Lives in the Observation model
  # group because it is phenotype-keyed, but must be re-drawn each replicate.
  cat["phenotype_random_effects"] <- "store_and_reset"
  cat
}
```

Keep the three arguments on `archive_replicate()` so a caller can still override,
but default them from this rather than from literals. Then a new table is
archived correctly the moment it is registered in `.schema_table_order()`, and
the "add a table" checklist does not grow a fourth item.

**Test to add:** every `SYSTEM_TABLES` entry is either assigned a category or
explicitly `NA` — the same visible-degradation pattern as items 1 and 3.

**Check first, before writing any of this:** whether `add_phenotype()` re-draws
`phenotype_random_effects` when rows already exist for a phenotype × effect. If
it re-draws, the omission is cosmetic (stale rows get overwritten) and this is
just an archiving gap. If it reuses existing draws, replicates are not
independent and that is a correctness bug that should be fixed and released on
its own, ahead of the tidy-up.

### Scope

`R/archive_replicate.R`, `R/schema.R` (one small helper), `tests/testthat/test-archive_replicate.R`.
Medium — the code is easy; the `phenotype_random_effects` question is the work.

---

## Suggested sequencing

1. **Answer the `phenotype_random_effects` re-draw question first** (item 3's
   "check first"). It is a yes/no read of `add_phenotype()` and it decides
   whether item 3 is a tidy-up or a bug fix.
2. **Item 1**, as its own commit. Self-contained, and it removes a hard error
   users can hit today.
3. **Item 3**, once 1 is in — it reuses `.schema_group_of()` and reads better
   after the registries are consistent.
4. **Item 2** whenever convenient; it protects against a future problem rather
   than fixing a present one.

Version bump + `NEWS.md` entry per `CLAUDE.md` before each commit.
