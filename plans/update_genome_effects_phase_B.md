# Genome Effects — Phase B results

**Spec:** `plans/update_genome_effects_v4.md` (v4.9). **Branch:** `feat/genome-effects-v49`.
**Follows:** `plans/update_genome_effects_phase_A.md`.
**Status:** complete. Gates 46–48 and 54 met.
**Date:** 2026-09-10.

Phase B commits the schema. It creates the three effect tables and `ind_tgv`,
gives `genome_meta` the primary key that makes the `locus_id` foreign key
declarable, registers everything, and builds the containment checker and
cross-row validator. **It writes no effect rows and computes no value** — the
writer is Phase C and the evaluator is Phase D.

**This is the first phase that leaves the package red.** `genome_effects` keeps
its name with an incompatible shape, so `add_tbv()` and
`define_additive_effects()` — which still `SELECT locus_name, line_name,
genome_effect_type, base_allele_freq` from it — fail until Phase C and D. That is
the expected consequence of not writing a compatibility shim; see
**Expected red** below for the exact list.

---

## What shipped

| File | Change |
|---|---|
| `R/define_genome.R` | `genome_meta` gains `PRIMARY KEY (locus_id)`; the three effect tables and two views are created here, inside the existing transaction; `GENOME_TABLES` grows to ten and a new `GENOME_VIEWS` holds the two views; the preflight and `pop$tables` cover both |
| `R/open_pop.R` | the old `genome_effects` DDL is **deleted**, and the table drops out of `tables_created` |
| `R/define_trait.R` | `ind_tgv` joins the lazy DDL block beside `ind_tbv`; the `ind_tgv_total` view is created alongside it and registered |
| `R/genome_effects_helpers.R` | **new** — view SQL, the family signature, the containment lattices, and `validate_genome_effects()` |
| `R/sql_utils.R` | six registries updated (see below) |
| `R/schema.R` | descriptions for four tables and three views; `.schema_table_order()` placements |
| `R/archive_replicate.R` | `ind_tgv` joins `store_and_reset`; the two new effect tables join `store_once` |
| `R/tidybreed_pop.R` | `print()` counts **causal loci** from `genome_effect_loci` instead of QTL rows from a column that no longer exists (finding 10) |
| `tests/testthat/test-genome-effects-schema.R` | **new** — 30 tests |
| `tests/testthat/test-schema-registries.R` | `ind_tgv$replicate` joins `DEFERRED_COLS` |
| `tests/testthat/test-print-pop.R` | asserts "causal loci" rather than "QTL" |
| `tests/testthat/test-open_pop.R` | asserts the effect tables are **absent** until `define_genome()` |
| `CLAUDE.md`, `NEWS.md`, `DESCRIPTION` | schema sections rewritten; 0.65.0 |

### Registry entries

| Registry | Added |
|---|---|
| `TABLE_RESERVED_COLS` | `genome_effects` (rewritten), `genome_effect_members`, `genome_effect_member_origins`, `ind_tgv`, and all three views (every column of a view is derived, so every column is reserved) |
| `TABLE_PRIMARY_KEYS` | `ind_tgv` |
| `TABLE_ROW_KEYS` | `ind_tgv`; `genome_effects` **removed** — see finding 4 |
| `TABLE_NO_ROW_DELETE` | the three effect tables, with a reason naming `define_genome_effects(mode = ...)`; the three views, with a derived-view reason |
| `SYSTEM_TABLES` | three tables, `ind_tgv`, three views |
| `IND_TABLE_ID_IND_COLS` | `ind_tgv` |
| `.schema_table_order()` | Genome group + Results group |
| `.all_schema_descriptions()` | four tables and three views, every column described |

`test-schema-print.R` asserts `SYSTEM_TABLES`, `.schema_table_order()` and
`.all_schema_descriptions()` name the same objects; all three agree.

### Constraints, verified empirically in DuckDB 1.5.5 before being written

Every row-local invariant the plan promises is a declared SQL constraint, and
each was probed against a real insert:

- closed sets on `contrast_name` and `line_match_type`;
- the `contrast_name` ↔ state/centre branch, including `dosage_value <=
  copy_count_value` and `center_value IS NOT NULL` (gate 46);
- `copy_count > 0`; `parent_origin IN (1,2)` or NULL;
- `line_match_type` ↔ `line_name` agreement, and `'any'` ⇒ a parent (gate 41's
  row-local half);
- `genome_meta` PK, the `locus_id` FK, the member → term FK, and the **composite**
  `(id_genome_effect, member_slot)` FK (gate 47).

Deletes do not cascade, so a parent delete with live children is refused — which
is why row deletion on these tables is refused outright (finding 4).

---

## Errors and fixes found reviewing the code against the plan

Ten, of which two were errors **in the plan** that would have shipped as bugs.

### 1. `ind_tgv` must not carry a `replicate` column — plan error

v4.9's DDL declares `replicate INTEGER`, justified as mirroring "`ind_tbv`'s
column (`sql_utils.R:97`)". But `sql_utils.R:97` is the *reserved-columns
registry*, not DDL: `ind_tbv`'s actual `CREATE TABLE` has no such column, and
neither do any of the other result tables. The column exists only in the
**archive** copy, added by `.ensure_archive_table(add_replicate = TRUE)`.

Shipping the plan's DDL verbatim would have broken `archive_replicate()` the
first time it ran on a population with true genetic values:
`archive_replicate.R:153-158` is a collision guard that **refuses to stamp a
table that already contains a `replicate` column**. A column added to be helpful
would have made the function that uses it error.

**Fixed.** `ind_tgv` has no `replicate` column; `replicate` stays in
`TABLE_RESERVED_COLS` (so a user cannot create a conflicting one) and joins
`DEFERRED_COLS` in `test-schema-registries.R`, exactly as `ind_tbv` does.

### 2. `ind_tgv_total` cannot expose `replicate` — consequence of 1

The plan's view spec lists `replicate` among its columns. The working table never
has that column, so the view is `(id_ind, trait_name, tgv_total)`. Nothing is
lost: archived rows live in the archive database, where the column exists on the
archived table itself.

### 3. `'unknown'` inside a genotype multiset was incomparable with itself — implementation bug

My first cut of `.ge_predicate()` carried the raw `line_name` as the line
dimension. For an `'unknown'` row that is `NULL`, and the bijection matcher
requires `!is.na(line)`, so **any predicate containing an `unknown` demand failed
to compare with anything, including an identical one** — a duplicate scope would
have passed validation and then been ambiguous at evaluation.

`'unknown'` is legal on a genotype member: the plan forbids only `'any'` there,
and correctly so — `'any'` constrains nothing, while `'unknown'` matches copies
whose founding line really is NULL. Fixed with `.ge_line_token()`, which encodes
the dimension as `exact:<line>` / `unknown` / `any`. That also closes a latent
aliasing hole on **both** lattices: a line literally named `"unknown"` would
otherwise have masqueraded as the unknown-line scope. Both cases are now tested.

### 4. The plan contradicted itself on `remove_rows()` — plan error

§Table lifecycle both lists the three effect tables in `TABLE_ROW_KEYS` *and*
says "**Rejecting is preferred** — effect definitions are configuration and
should be replaced through the writer, not row-deleted." Those cannot both hold:
`remove_rows()` checks `TABLE_NO_ROW_DELETE` first (`remove_rows.R:228`), so
`TABLE_ROW_KEYS` entries for the same tables would be unreachable.

**Implemented the stated preference.** The three tables are in
`TABLE_NO_ROW_DELETE` with a reason that names the parent/child structure, the
absence of cascade, and `define_genome_effects(mode = ...)`. They are not in
`TABLE_ROW_KEYS`.

*Known cosmetic wrinkle:* `remove_rows()`'s no-filter guard fires **before** the
refusal, so an unfiltered call on `genome_effects` reports the generic "no filter
applied" message rather than the specific one. A filtered call gets the right
message. Left as is — reordering the guards is a `remove_rows()` change with its
own blast radius, and the user is stopped either way.

### 5. "The multiset must sum to the realized copy count" is checkable at write time

The plan states this rule (§Origin resolution) but assigns it no enforcement, and
my first cut only bounded the sum at ≤ 2. It is exactly mechanical:

- a `dominance` member's state is the diploid genotype, so its multiset must sum
  to **2**;
- an `indicator` member's state declares its own `copy_count_value`, so its
  multiset must sum to **that**.

Both are now enforced and tested. This is what makes a one-copy scope on a
dominance member — meaningless, since a genotype needs both copies — a write-time
error rather than an evaluation-time surprise.

### 6. The dominance ploidy rule needs the database

The plan lists "dominance rejected at non-diploid loci unless proven on that
member" under *Enforced in R*, but its inputs are `chr_inheritance` and
`genome_meta`, not the candidate rows. Implemented as
`.ge_validate_dominance_ploidy(conn, ...)`, separate from the frame-level rules,
honouring the "proven on **that same member**" escape — a two-copy multiset on a
*different* member of the term proves nothing, and there is a test for exactly
that.

**Limitation recorded:** resolution is called with `line_name = NULL`, so a
line-specific karyotype rule is not consulted. Every current fixture and every
`define_chromosome()` example is line-agnostic. Revisit if line-specific
inheritance rules are ever written.

### 7. Orphan member rows were skipped when there were no terms

`.ge_validate_frames()` returned early on `nrow(terms) == 0`, so member rows with
no term at all — the very thing the orphan check exists for — passed silently.
The early return now requires all three frames to be empty. Only reachable
through the frame-level entry point (the database's foreign keys make it
impossible in a live population), which is precisely why it needed a test:
Phase C validates candidate frames *before* issuing any SQL, where no FK protects
it.

### 8. A roxygen block was silently re-attached to the wrong object

Inserting `.GE_NO_DELETE_REASON` between the `TABLE_NO_ROW_DELETE` roxygen block
and its assignment moved the documentation onto the new constant.
`roxygen2::roxygenise()` then **deleted `man/TABLE_NO_ROW_DELETE.Rd`** and wrote
`man/dot-GE_NO_DELETE_REASON.Rd` in its place — no error, no warning, and no test
covers it. Caught by regenerating the docs and reading `git status`. Fixed by
defining the constant **above** the block.

Worth remembering for Phase C, which inserts into several existing files: an
in-place insert immediately after a `#'` block silently steals that block.

### 9. The views were in `SYSTEM_TABLES` but in no row-key registry — regression

`test-schema-registries.R:160` asserts
`names(TABLE_ROW_KEYS) ∪ names(TABLE_NO_ROW_DELETE) == SYSTEM_TABLES`, a guard
added in 0.64.2 precisely so that registering a new object **forces a decision**:
give it a row key, or say why deletion is refused. The three views had neither,
so Phase B broke the guard that exists to catch Phase B.

Fixed by adding all three to `TABLE_NO_ROW_DELETE` with a view-specific reason
("this is a derived view, not a table: it has no rows of its own"), and by giving
them `TABLE_RESERVED_COLS` entries so `mutate_table()` reports *reserved* rather
than letting `ALTER TABLE` fail on a view with a raw DuckDB message. A test now
asserts each view lists exactly the columns it reserves.

The same guard-ordering wrinkle as finding 4 applies: `remove_rows()` checks for
an unfiltered or empty selection **before** consulting `TABLE_NO_ROW_DELETE`, so
the specific refusal is reached only by a filtered call that matches rows.

### 10. `print(pop)` broke in Phase B, but its fix was scheduled for Phase E

`print.tidybreed_pop()` reported `n_qtl` from
`SELECT COUNT(DISTINCT locus_name) FROM genome_effects WHERE genome_effect_type
= 'additive'` (`tidybreed_pop.R:159`). Both of those columns disappear in Phase
B, so **printing any population** — the most basic operation there is — would
have errored for the whole of Phases B, C and D. The plan schedules that line for
Phase E, but Phase E is about the *terminology*; the mechanical break lands where
the column does.

Pulled forward and fixed: the count is now
`SELECT COUNT(DISTINCT locus_id) FROM genome_effect_loci`, and the label reads
**"causal loci"** rather than "QTL" — which is the plan's own terminology
decision, and the honest word once one coefficient can span several loci. The
additive-QTL / epistatic-only breakdown stays out of the print header, as the
plan says.

`extract_genotypes()`'s `effects_tbl` path (`extract_genotypes.R:118-127`) has
the same problem and **stays in Phase E**: it is an opt-in argument, not a basic
operation, and moving it properly means moving it onto `genome_effect_loci`.

---

## Deviations from the plan, and why

| Plan says | Built | Reason |
|---|---|---|
| `ind_tgv` has `replicate INTEGER` | no such column | finding 1 — it would break `archive_replicate()` |
| `ind_tgv_total` exposes `replicate` | `(id_ind, trait_name, tgv_total)` | finding 2 |
| `TABLE_ROW_KEYS` gains the three effect tables | `TABLE_NO_ROW_DELETE` instead | finding 4 — the plan's own stated preference |
| (validator table lists the ploidy rule with no home) | `.ge_validate_dominance_ploidy(conn, ...)` | finding 6 |
| (canonicalization listed, unspecified) | `member_slot` must be `1..n` in ascending `locus_id` order | makes "canonicalized by ascending `locus_id`" mechanical, and the view's `family_key` depends on it |
| preflight covers ten tables | ten tables **and** the two views | a leftover view would break `CREATE VIEW` just as a leftover table breaks `CREATE TABLE` |
| `tidybreed_pop.R:159` moves in Phase **E** | moved in Phase **B** | finding 10 — the column it reads disappears here, and `print(pop)` is not something to leave broken for three phases |

---

## Gate status

| Gate | Status |
|---|---|
| **46** — direct insert of an `additive`/`dominance` member with NULL `center_value` fails at the SQL constraint | ✅ |
| **47** — duplicate `locus_id` refused; member naming a nonexistent `locus_id` refused by the FK that could not be declared before the PK existed | ✅ |
| **48** — `schema()` lists the same objects fresh as restored, views included, each with a description | ✅ |
| **54** — `family_key` is honest: scope variants share it, different contrasts and different member sets do not; the view's key equals the R signature | ✅ |
| **41** (row-local half, early) — `('any', parent NULL)` refused by the `CHECK`, not the validator | ✅ |

---

## Expected red

The suite is **not** green at the end of Phase B, by design. Every failure below
is a reader of the old `genome_effects` shape, and each is scheduled:

| File | Cause | Fixed in |
|---|---|---|
| `test-define_additive_effects.R` | `define_additive_effects()` writes `locus_name`, `line_name`, `genome_effect_type`, `base_allele_freq` | **C** |
| `test-add_tbv.R`, `test-add_tbv_index.R` | `add_tbv()` reads the same columns | **D** |
| `test-add_phenotype.R`, `test-phenotype_composite.R`, `test-formula_phenotype.R` | call `add_tbv()` internally | **D** |
| `test-parity.R`, `helper-parity.R` | seeded end-to-end run through both | **D** |
| `test-schema-registries.R`, `test-print-pop.R`, `test-summary_pop.R`, `test-mutate_table_defaults.R`, `test-mutate_derived.R`, `test-sql_injection_hardening.R`, `test-define_effect_fixed_cov.R`, `test-define_effect_intercept.R`, `test-extract_genotypes.R`, `test-dosage_extract_parity.R` | fixtures call `define_additive_effects()` | **C** |
| `test-archive_replicate.R`, `test-remove_rows.R` | same fixtures | **C** |

`test-open_pop.R` was **not** in this category — see Verification.

`test-genome-effects-fixtures.R` (Phase A) and `test-genome-effects-schema.R`
(Phase B) are green and stay green.

---

## Verification

```
test_file("tests/testthat/test-genome-effects-schema.R")   96 pass, 0 fail
test_file("tests/testthat/test-genome-effects-fixtures.R") 248 pass, 0 fail   (Phase A, still green)
```

Before any DDL was written into the package, the whole schema was probed against
a scratch DuckDB 1.5.5 database — every `CHECK`, every foreign key including the
composite one, `ALTER`/`UPDATE` on a primary-keyed `genome_meta`, and all three
views. Nothing in the plan's schema section was taken on trust.

**Full suite: 1,887 pass, 3 fail, 180 error.** Every one of those 183 traces to a
single cause — `Binder Error: Referenced column "genome_effect_type" not found` —
i.e. `define_additive_effects()` (Phase C) or `add_tbv()` (Phase D) reading the
old column shape, directly or through a fixture. **Nothing failed for any other
reason.**

One test was a genuine Phase-B obligation rather than expected red and was fixed
here: `test-open_pop.R` asserted that `open_pop()` creates `genome_effects`,
which Phase B deliberately reverses. It now asserts the opposite, plus that
`ind_tbv` / `ind_tgv` / `ind_tgv_total` *are* present — the trait/phenotype block
still runs in `open_pop()`.

---

## Next

**Phase C** — `define_genome_effects()` with the `terms` data-frame format;
`define_additive_effects()` rebuilt on it with `replace_scope` and
`parent_origin`; the `(a, d)` and genotype-table helpers; origin-aware
`scale_to_target`; the parent-only re-run warning; and the complete deletion of
`trait_meta.expressed_parent` (fourteen sites, including two executable vignette
calls). Gates 34–35, 41–44, 50 and 53.

Phase C's validator work is already done: `.ge_validate_frames()` takes candidate
data frames precisely so the writer can validate before issuing any SQL, and
finding 7 exists because that path has no foreign keys behind it.
