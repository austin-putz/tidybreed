# Plan: individual ids out of SQL text, and id numbering that survives `_` in line names

Status: revised after a Codex review (`plans/individual_id_sql_codex_review.md`,
2026-11-03); its accepted points are folded in below. Nothing implemented.
Decisions taken 2026-11-03: D1(a) per-line numbering, D2(a) focused fix plus
a follow-up plan, D3 no migration. **Implemented in 0.76.14** (`.next_ind_number()` in `R/sql_utils.R`, `tests/testthat/test-individual_ids.R`;
follow-up: `plans/individual_id_sql_followup.md`). One change from §2a found
by test 5: ids are built as `(start - 1L) + k`, since `start + k - 1L`
overflows at the limit before the subtraction. Test 9 checks individual ids in
every statement, and line names only in `ind_meta` statements: the map and
rule helpers quote line names with `dbQuoteLiteral()`, which the rule allows.

Target version for the implementation: 0.76.14 (this plan ships with the
0.76.13 docs/test work).

## 1. The problems

Three problems, all in how `add_founders()` and `add_offspring()` talk to
`ind_meta`. Problems 1 and 2 are bugs that stop real crossbreeding designs.
Problem 3 is a breach of the CLAUDE.md hard rule "individual ids never appear
in SQL text".

All three were reproduced on 0.76.13 (2026-11-03), and again independently by
the review.

### Problem 1 — a line whose name contains `_` can be bred only once

Line names may contain underscores: `add_founders()` and `add_offspring()` both
accept `^[a-zA-Z][a-zA-Z0-9_]*$`. Ids are `{line_name}_{n}`. The next `n` is
found by (`R/add_founders.R:197-201`, `R/add_offspring.R:462-470`):

```sql
SELECT MAX(CAST(SUBSTRING(id_ind FROM POSITION('_' IN id_ind) + 1) AS INTEGER))
FROM ind_meta WHERE LEFT(id_ind, <nchar(ln) + 1>) = '<ln>_'
```

`POSITION('_' IN id_ind)` is the **first** underscore. For line `Duroc_F1` the
first call works (no rows yet, so `MAX` is `NULL` and numbering starts at 1),
writing `Duroc_F1_1`, `Duroc_F1_2`. Every later call fails:

```
Conversion Error: Could not convert string 'F1_1' to INT32
```

Same failure for a second `add_founders()` call into such a line.

Why tests miss it: lines with `_` appear in `test-add_founders.R`,
`test-extract_allele_freq.R`, `test-add_tgv_breeding_value.R`, but each line
receives animals only once.

### Problem 2 — a line name that is a prefix of another line's name

The rows are found by the id prefix `'<ln>_'`, not by `line_name`. Line `A`'s
prefix `A_` also matches line `A_B`'s ids. Once `A_B_1` exists, any
`add_founders()` / `add_offspring()` into line `A` fails with
`Could not convert string 'B_1' to INT32`. If the cast had succeeded (say a line
`A_7` with ids `A_7_1`), line `A` would silently skip numbers instead.

### Problem 3 — parent ids pasted into SQL in `add_offspring()`

`R/add_offspring.R:234`:

```r
parent_id_list <- paste(paste0("'", unique_parents, "'"), collapse = ", ")
```

used at line 238 (`ind_meta` lookup) and line 367 (`ind_haplotype` load).
`unique_parents` comes straight from the user's `matings`, unescaped:

- a parent id containing `'` gives a raw DuckDB parser error
  (`syntax error at or near "A_3"`) instead of "Parent ID(s) not found";
- it is an injection path in principle (arbitrary user strings in SQL text);
- the statement text grows with the number of parents. This plan makes no
  performance claim; any claim would first need a benchmark under
  `dev/benchmarks/` (CLAUDE.md).

The two numbering queries also put the line name in SQL text. That one is safe
in practice (validated by the regex, which allows no quote) but is the same
pattern.

### Problem 4 (found in review) — the counter has no defined numeric range

Not a likely population size, but undefined behaviour that the new helper
should not inherit:

- `TRY_CAST(... AS INTEGER)` turns a suffix above 2,147,483,647 into `NULL`,
  so a line holding only `Huge_2147483648` would restart at 1 (reproduced).
- At `max = 2147483647`, R's `max + 1L` is `NA` with a warning, and
  `add_offspring()` increments its counter after assigning each id
  (`R/add_offspring.R:494`), so the increment after the last permitted id
  overflows even when that id itself is fine.
- `add_founders()` builds ids with `seq(start_id, start_id + n_founders - 1)`
  (`R/add_founders.R:222`); `n_founders` can be a double, so this can write an
  id beyond the int32 range today.

### The wider pattern (not in scope unless D2 says so)

Ids reach SQL text in other modules too. None of these take ids from user input
directly — they come from `resolve_subset_ids()` or a table, so they are
tidybreed-generated `{line}_{n}` strings with no quote — so they are rule
breaches, not live bugs:

| File | Lines | Quoting |
|---|---|---|
| `R/add_dosage.R` | 100–111 | `sql_in_list()` (escapes `'`) |
| `R/add_genotypes.R` | 79–82 | `sql_in_list()` |
| `R/extract_genotypes.R` | 163–167, 251 | `sql_in_list()` |
| `R/ploidy_helpers.R` | 20 (`assert_ploidy_2()`, called from `add_dosage.R:77`, `extract_genotypes.R:189`) | `sql_in_list()` |
| `R/add_ebv.R` | 291–317 | **none** |
| `R/blupf90_helpers.R` | 34–49, 118–127, 156–160, 214, 236 | **none** |

The pattern already used where the rule is honoured: register the ids as a view
and join (`.ap_read_by_id()` in `R/add_phenotype_stages.R:514`,
`resolve_subset_ids()`, `R/contributor_tgv.R:93`, the `tmp` joins in
`R/add_tgv.R:479, 544`, `remove_rows()`).

## 2. Proposed fix

### 2a. One id-allocation helper, by `line_name`

New internal helper in `R/sql_utils.R`, used by both writers:

```r
# First id number to use per line for ids "{line_name}_{n}", given how many
# new ids each line needs. A registered-view query: line names never appear
# in SQL text. Errors (naming the line) on malformed stored ids or when the
# allocation would leave the supported range.
.next_ind_number <- function(conn, n_new)   # n_new: integer, named by line
```

**The id format** this helper reads and the writers produce: `id_ind` is
exactly `line_name || '_' || n`, where `n` is a positive decimal integer with
no leading zero (`[1-9][0-9]*`) and `1 <= n <= 2147483647`
(`.Machine$integer.max`, the int32 range the counter uses). Nothing in R writes
any other form.

**The query.** Register `data.frame(line_name = names(n_new))` as
`__next_id_lines` and read, per line, the stored ids of that line:

```sql
WITH r AS (
  SELECT l.line_name, m.id_ind,
         substr(m.id_ind, length(l.line_name) + 2) AS suffix,
         starts_with(m.id_ind, l.line_name || '_')
           AND regexp_full_match(substr(m.id_ind, length(l.line_name) + 2),
                                 '[1-9][0-9]*') AS well_formed
  FROM __next_id_lines l
  LEFT JOIN ind_meta m ON m.line_name = l.line_name
)
SELECT line_name,
       COUNT(*) FILTER (WHERE id_ind IS NOT NULL AND NOT well_formed) AS n_malformed,
       COUNT(*) FILTER (WHERE well_formed AND length(suffix) > 10)    AS n_too_long,
       MAX(CASE WHEN well_formed AND length(suffix) <= 10
                THEN TRY_CAST(suffix AS BIGINT) END)                AS max_num
FROM r
GROUP BY line_name
```

plus up to five malformed example ids per line for the error message. Reading
the **whole remainder** after the exact `line_name || '_'` prefix (rather than
the digits after the last `_`) is what makes `Bad_extra_99` in line `Bad`
malformed instead of `n = 99`. `BIGINT` holds every suffix of up to 10 digits,
so one above the int32 range is seen as a number, not dropped; a longer one is
counted in `n_too_long`. The cast must stay **inside** the `CASE`: DuckDB
evaluates a plain `CAST` before an aggregate `FILTER`, so
`MAX(CAST(suffix AS BIGINT)) FILTER (WHERE well_formed)` fails on a malformed
row (checked on DuckDB 1.5.5 while revising this plan; the query above was run
on lines `A`, `A_B`, `A_7`, `A_`, `Duroc_F1`, `Bad`, `Huge`, an 11-digit suffix
and an unseen line, and gave the expected result for each).

**Policy, in R, before any RNG use or write:**

- **Malformed ids are an error**, naming the line, the count and the examples
  (`Bad_0` and `Bad_007` are malformed too: no zero, no leading zero). Only a
  hand-edited table can produce one (`line_name` is reserved against
  `mutate_table()`), and silently skipping corrupt identity data could hide a
  damaged population.
- **Range.** `next = max + 1` (or `1` for a line with no rows) and the whole
  allocation are checked in double arithmetic, which cannot overflow at these
  sizes: `next + n_new - 1 <= .Machine$integer.max`, else an error naming the
  line and the supported range. A stored suffix already above the range
  (`max_num > .Machine$integer.max`, or `n_too_long > 0`) fails the same way.
- **Mapping.** The result is built from the returned `line_name` column and
  indexed by `names(n_new)`; `GROUP BY` promises no order.

**The writers** get the start per line and build every id in one vectorised
step, `paste0(line, "_", start + seq_len(k) - 1L)` in integer arithmetic that
the range check has already proven safe. There is no increment after the last
id (removes problem 4's post-increment). `add_offspring()` assigns ids in
mating-row order exactly as today: row order decides which offspring gets which
id and which per-gamete stream.

- Rows are selected by the `line_name` **column**, so line `A` never sees
  `A_B`'s animals (problem 2).
- `n` is the remainder after the line's own prefix, so `Duroc_F1_3` reads as 3
  (problem 1), `A_7_1` in line `A_7` as 1, and `A__1` in line `A_` as 1.
- Neither writer's RNG use changes: the helper only registers a view and reads
  (no `dbWriteTable()`), and it runs before every draw as today.

**What uniqueness means.** "Current maximum" is over the rows `ind_meta`
holds now. Within a population the new id is unique by construction
(`{line}_{n}` splits uniquely at the line's own prefix, and `n` exceeds every
stored `n` of that line), with `ind_meta.id_ind`'s primary key
(`R/open_pop.R:267`) as the backstop. It is **not** a lifetime counter: deleting
a line's highest ids lets their numbers be reused, and `archive_replicate()`
deliberately restarts numbering every replicate (`R/archive_replicate.R:43-48`).
That design stays; the docs must not imply global or gap-free numbering.

### 2b. Parent ids through a registered view in `add_offspring()`

**Validation first.** After alias normalisation (`id_sire`/`id_dam`) and the
required-column check, and before the ids are combined or registered, check
`id_parent_1` and `id_parent_2` each: a factor is converted with
`as.character()`; any other non-character type, an `NA`, or an empty string is
an error naming the column. The generated-id regex is **not** applied to these
lookup values: an unknown id of any shape must reach the ordinary "Parent ID(s)
not found in ind_meta" error. Today an `NA` parent becomes the string `'NA'`
and is reported as a missing parent.

**Then one registered view** of `unique_parents`, used by both reads, with
function-scoped cleanup so it is removed on success and on error:

```r
duckdb::duckdb_register(conn, "__ao_parents",
                        data.frame(id_ind = unique_parents, stringsAsFactors = FALSE))
on.exit(try(duckdb::duckdb_unregister(conn, "__ao_parents"), silent = TRUE), add = TRUE)
```

- `ind_meta` read: `SELECT m.id_ind, m.sex, m.ploidy, m.line_name FROM ind_meta m
  JOIN __ao_parents USING (id_ind)`.
- `ind_haplotype` read: same join, keeping
  `ORDER BY id_ind, parent_origin, locus_id` exactly (the per-parent split
  depends on it).

The numbering helper's `__next_id_lines` view gets the same function-scoped
cleanup. Fixed internal names follow repository practice; no general helper
framework is needed for this patch.

`duckdb_register()` does not advance the RNG (that is why the
`add_phenotype()` path uses it). Section 4 asserts this directly.

## 3. Decisions

**D1 — numbering scope.** (a) Per line, as today: `n` counts within each line,
so `A_1` and `B_1` both exist. (b) One counter across the population, ids
still prefixed (`A_1`, `B_2`, `A_3`). Recommendation: **(a)**. The
`line_name` column defines the allocation scope, it keeps the current id format
and every current id, and an id is an identifier, not an ordering (row order is
never implied). Review: agrees. **Decided: (a).**

**D2 — scope of the id-in-SQL cleanup.** (a) This change: `add_founders()` and
`add_offspring()` only, plus a follow-up plan for the six files in the table
above, `assert_ploidy_2()` and its callers included. (b) All of them now,
through one shared `.with_id_view()` helper generalising `.ap_read_by_id()`.
Recommendation: **(a)**. Problems 1–4 are small and testable on their own, and
the wider change deserves its own review. That follow-up is more tractable than
CLAUDE.md's coverage figures for `add_ebv.R` / `blupf90_helpers.R` (29% / 22%)
suggest: `ebv_parent_avg()`, pedigree collection and input-file preparation
can be tested on small database fixtures without the BLUPF90 binary. Review:
agrees. **Decided: (a).**

**D3 — existing databases.** No migration code. This is a query and
validation change with no schema change. Every id the old code wrote in
practice reads correctly with the new helper (a second call into a problem line
failed instead of writing). The one exception is theoretical: an id above the
int32 range, which `add_founders()`'s double `seq()` could write (problem 4).
The new helper reports it as an out-of-range error rather than misreading it;
that is validation, not a compatibility layer. Review: agrees (and pointed out
that the earlier wording, "every old id reads correctly", was too strong).
**Decided: no migration.**

## 4. Tests

New `tests/testthat/test-individual_ids.R`. Every comparison of collected
tables sorts by explicit keys first (`id_ind`; `id_ind, parent_origin,
locus_id`).

**Numbering**

1. **Underscore line, both writers, both directions.** Into `Duroc_F1`:
   founders, founders, offspring, offspring, then founders again. Ids run
   `Duroc_F1_1 ...` continuously with no gap; parent ids are the ones given.
2. **Prefix lines, interleaved.** Lines `A`, `A_B`, `A_7` and `A_` (trailing
   underscore, ids `A__1`) created in mixed order, several calls each, and one
   `add_offspring()` call whose mating rows interleave all four. Each line
   numbers from 1 with no gap and no line sees another's ids; the id each
   mating row receives follows row order.
3. **Helper mapping.** `.next_ind_number()` called with an unsorted request
   that includes a never-seen line returns the right start for each name.
   (The writers pass unique names; duplicate names in the request are an
   error.)
4. **Malformed stored id.** A tiny `ind_meta` fixture with line `Bad` and id
   `Bad_extra_99` (also `Bad_0`, `Bad_007`): the writer errors naming the line
   and examples, and the database is unchanged.
5. **Range boundary.** Tiny fixtures, no large population:
   - stored `Lim_2147483646`, one new id → `Lim_2147483647`, no warning;
   - stored `Lim_2147483647`, one new id → range error;
   - stored `Lim_2147483640`, ten new ids (an allocation that crosses the
     limit part-way) → range error;
   - stored `Huge_2147483648` → range error, never a restart at 1;
   - `add_founders()` with `n_males` / `n_females` given as doubles near the
     limit → range error.
   Every error names the line and the supported range, and leaves `ind_meta`
   and `ind_haplotype` unchanged.
6. **Reopen.** A small file-backed population: add animals to `Duroc_F1`,
   `close_pop()`, `restore_pop()`, add again. Numbering continues: it comes
   from stored data only.

**Parent ids**

7. **SQL-shaped unknown parent.** `id_parent_1 = "A_1'); --"` gives "Parent
   ID(s) not found in ind_meta", not a parser error, and the population is
   unchanged.
8. **Parent validation.** For `id_parent_1` and `id_parent_2` separately: `NA`,
   `""` and a numeric column each give the named-column error; a factor column
   works like the character one.
9. **No id literals in executed SQL.** Wrap `DBI::dbGetQuery()`,
   `DBI::dbExecute()` and `DBI::dbSendQuery()` with
   `testthat::local_mocked_bindings(..., .package = "DBI")` so every statement is
   recorded and forwarded to the real connection. Run one successful
   `add_offspring()` (which reaches both the `ind_meta` read and the
   `ind_haplotype` read) and one successful `add_founders()`, using parent ids
   and line names that are distinctive strings. Assert that no recorded
   statement contains any of them. This checks the rule itself; escaping would
   pass test 7 but fail this one. (No source-file grep: it misses `= '<id>'`
   forms, rejects the legitimate `IN (SELECT id_ind FROM <view>)`, and reads
   repository source, which the installed-package coverage run does not have.)
10. **View cleanup.** After a successful call and after a failed lookup, no
    `__ao_parents` / `__next_id_lines` view remains, and a valid call on the
    same connection then succeeds.

**RNG**

11. **RNG neutrality, asserted.** `.Random.seed` is identical before and after
    `.next_ind_number()`, after a parent-validation error, and after a
    missing-parent error (all of which happen before any draw). Test code may
    read `.Random.seed`; production code never restores it.
12. **Reproducibility.** For lines `A` and `B`, a seeded `add_founders()` +
    `add_offspring()` writes the literal id sequence expected, and two runs with
    the same seed give `expect_identical()` `ind_meta` and `ind_haplotype`.
    Per CLAUDE.md this is within-current-code reproducibility. Matching the
    previous version's seeded output is **not** required and not pinned.

Then run `test-add_founders.R`, `test-add_offspring*.R`,
`test-recombination*.R`, `test-make_gametes_parity.R`, `test-genome_map.R`,
`test-archive*.R`, and the full suite.

## 5. Docs

- `?add_founders` / `?add_offspring`: ids are `{line_name}_{n}`, numbered per
  line from the largest `n` the line currently holds in `ind_meta`; line names
  may contain `_`; the supported range of `n`; numbers can be reused after
  deletion and restart after `archive_replicate()`.
- `tidybreed-schema` skill: the `ind_meta.id_ind` row states the format, range
  and current-maximum semantics.
- `tidybreed-api` skill: the `add_founders()` section ("ID format:
  `{line_name}_{n}`") gains the same, plus "allocated by `.next_ind_number()`".
- `NEWS.md` 0.76.14: the two bugs (with the error text users saw), the range
  errors, and the parent-id change.
- If D2 = (a): `plans/individual_id_sql_followup.md` listing the six files and
  `assert_ploidy_2()`'s callers.

## 6. Out of scope

- Changing the id format or widening the counter beyond int32.
- A persisted counter (it would be stored metadata; CLAUDE.md).
- The C++ kernel and recombination (verified separately in 0.76.13).
- `define_chip.R:102` and `extract_genetic_variance.R:633` paste **locus** ids
  (integers). Not individual ids, not user strings; leave unless D2 = (b).
- Performance work.
