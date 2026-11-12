# Plan (follow-up): individual ids out of SQL text in the remaining modules

Status: revised after two Codex review rounds
(`plans/individual_id_sql_followup_codex_review.md`, 2026-11-03 and 2026-11-04);
the second approved the design, and its test clarifications are folded in
below, along with a final review against the code (2026-11-04). Implemented in 0.76.15; see `plans/individual_id_sql_followup_report.md`. Follows `plans/individual_id_sql.md`
(D2 = (a)), which fixed `add_founders()` and `add_offspring()` in 0.76.14.

## Why

CLAUDE.md: "Individual ids never appear in SQL text — register them as a view."
These sites still paste ids into `IN (...)` lists. None takes ids from user
input directly (they come from `resolve_subset_ids()` or a table, so they are
tidybreed-generated `{line}_{n}` strings with no quote), so they are rule
breaches, not live bugs. Two files do not even escape.

**Scope:** this removes manual interpolation of *resolved individual-id
vectors* in the six modules below. It does not make that hard rule true
repository-wide: `resolve_subset_ids()` renders pending user filters with
`dbplyr::sql_render()` (`R/sql_utils.R`), so a user's
`filter(id_ind %in% c("A_1", "A_2"))` still puts those strings in SQL. Rewriting
dbplyr predicates would need its own design and is out of scope. The hard rule's
wording is narrowed to match (see Docs).

## Inventory (re-checked on the 0.76.14 working tree)

| File | Lines | What | Quoting |
|---|---|---|---|
| `R/add_dosage.R` | 100–101, 109–111 | read `ind_haplotype`; `DELETE FROM ind_genotype` | `sql_in_list()` |
| `R/add_genotypes.R` | 79–82 | `UPDATE ind_meta SET <chip> = TRUE` | `sql_in_list()` |
| `R/extract_genotypes.R` | 163–167, 244–251 | genotype reads | `sql_in_list()` |
| `R/ploidy_helpers.R` | 20 | `assert_ploidy_2()`; callers `add_dosage.R:77`, `extract_genotypes.R:189` | `sql_in_list()` |
| `R/add_ebv.R` | 291–295, 302–317 | subset read; `ebv_parent_avg()` parent read | **none** |
| `R/blupf90_helpers.R` | 34–38, 45–49 (`trace_pedigree()`), 118–127 (`build_data_file()`), 156–160 (fixed-effect source tables), 210–214, 230–236 (`write_geno_file()`) | pedigree, phenotypes, fixed-effect source tables, genotyped animals | **none** |

A repository-wide search for `id_ind IN (`, `id_ind = '`, `id_ind NOT IN (`
and `id_parent_[12] IN (` (excluding `IN (SELECT ...)`) finds exactly these
sites and nothing else. Every id list comes from `resolve_subset_ids()` or
from a table read, so the values are ids tidybreed wrote; since 0.76.14 the
only writers are `add_founders()` / `add_offspring()` through
`.next_ind_number()`, so they are `{line_name}_{n}` with no quote (both writers
hold `line_name` to `^[a-zA-Z][a-zA-Z0-9_]*$`).

Not individual ids, out of scope here: trait / phenotype names in `add_ebv.R`
(206, 238–240, 303) and `blupf90_helpers.R` (81–88), and the `model` label in
`add_ebv.R` (315, 322) — validated by `validate_sql_identifier()` in `add_ebv()`
(line 156–157) and could move to `DBI::dbQuoteLiteral()`; the table and column
identifiers pasted at `blupf90_helpers.R:159`; integer `id_phenotype` and
`locus_id` lists.

## Proposed change

One internal helper in `R/sql_utils.R`, so genotype and EBV code do not depend
on phenotype-stage internals:

```r
# Run `fun(view_name)` with `ids` registered as a one-column view `id_ind`,
# unregistered on exit (success or error), and return fun()'s result.
.with_id_view <- function(conn, ids, fun)
```

Contract:

- **Set semantics.** It registers `unique(ids)`, as `.ap_read_by_id()` does. A
  join against repeated ids multiplies source rows, which changes a dosage
  `SUM()` rather than only duplicating output (DuckDB 1.5.5: 2 through the
  join, 1 through `IN`).
- **Input.** `ids` is a character vector without `NA`, keeping what
  `sql_in_list()` enforces today, and `NULL` is an error. There is no id-format
  check: an unknown string is data, not misuse.
- **Empty.** `character(0)` means nobody and never everyone: `fun()` runs
  against an empty `VARCHAR` view. `data.frame(id_ind = NULL)` has no column,
  so build the frame with `character()`. Public guards are unchanged:
  `add_dosage()`, `add_genotypes()` and `add_ebv()` warn and return,
  `extract_genotypes()` errors. `assert_ploidy_2(pop, character(0))` changes
  from an error (via `sql_in_list()`) to passing on nobody; both callers guard
  first, so nothing public moves.
- **Fresh names.** `__id_view_<k>`, from a counter in a private environment
  (namespace bindings are locked; no RNG). Leave `overwrite` at its default, so
  `duckdb_register()` errors on a live name instead of replacing it. A fixed
  name would therefore fail when calls nest, not silently swap the view. Today's
  callers run sequentially and nothing nests; fresh names keep a shared callback
  helper safe if one later holds its view open across another helper.
- **Cleanup.** Install `on.exit()` only after registration succeeds, so a failed
  registration cannot unregister a view the helper does not own. Wrap the
  unregister in `try(..., silent = TRUE)`, the repository's pattern, so a cleanup
  failure cannot replace the original error. RNG-neutral: register and
  unregister only.
- **No open transaction.** The helper is not for use inside a transaction: in an
  aborted transaction `duckdb_unregister()` errors and the view survives the
  rollback. None of the six modules opens one.
- **Eager.** `fun()` finishes all database work before returning. A lazy table
  over the view is not a valid return value.

Each site becomes a `JOIN <view> USING (id_ind)` or
`id_ind IN (SELECT id_ind FROM <view>)`. `UPDATE` / `DELETE` use the
`IN (SELECT ...)` form. Where `resolve_subset_ids()` returns `NULL` (all
individuals), keep the existing no-filter branch; do not call the helper.

**One view per id set, not per function.** A view holds the frame it was
registered with, so these sites need care:

| Site | View lifetime |
|---|---|
| `add_dosage()` | One view spans the `DELETE` and the `INSERT`. Keep the concrete `resolve_subset_ids()` snapshot: a selection from `ind_genotype` itself must not become a live subquery, or the `DELETE` empties the selection before the `INSERT`. |
| `trace_pedigree()` | `all_ids` grows each generation: register per iteration, then the final vector for the final read. One view from before the loop would miss later generations. |
| `extract_genotypes()` | The chip path reads `subset_ids`; the genotype read uses `has_ids`, which that read narrows. Two id sets. |
| `ebv_parent_avg()` | The candidate read uses `subset_ids`; the parent EBV read uses `parent_ids`. Two id sets. |
| `build_data_file()` | Phenotypes read `subset_ids`; fixed-effect sources read `pheno_wide$id_ind`. Two id sets. |
| `write_geno_file()` | The genotyped-animal read uses `all_ped_ids`; the dosage read uses `geno_ids`. Reusing the broader view returns animals absent from the matrix. |

**Row order becomes explicit.** Neither form preserves row order: on a
1M-row table with 300k ids at 8 threads, two identical runs of the same `JOIN`,
and likewise of the same `IN (SELECT ...)`, returned different orders. That
would break the bit-identical contract for output built from the result, so add
`ORDER BY id_ind` (matching `resolve_subset_ids()`) to:

- the `has_ids` reads in `extract_genotypes()` (matrix rows), including the
  unrestricted `SELECT DISTINCT id_ind FROM ind_haplotype`;
- the final read in `trace_pedigree()` (`pedigree.txt`);
- the `geno_ids` read in `write_geno_file()` (`genotype.txt`).

The other reads are matched back by key (`match()`, `merge()`, per-id lookup)
and need no ordering, provided the key is unique in the rows read. Every one of
them has a unique key except the fixed-effect source read, covered next.
Pre-1.0 policy does not require keeping the old incidental order.

**Refuse duplicate source rows.** `build_data_file()` reads a fixed effect from
any `source_table` and takes the first `match()`. With several rows per animal,
"first" depends on row order, so the value written to `data.txt` is not stable:
on 300,000 ids × 2 rows at 8 threads, a standalone DuckDB run of the planned
query forms changed the first match for up to 191,520 animals with
`IN (SELECT ...)` and 68,640 with `JOIN`. `ORDER BY id_ind` cannot settle it, and
a generic source table has no known column to break the tie.

`add_phenotype()` already refuses this (`R/phenotype_helpers.R:123`), so give
`build_data_file()` the same refusal, naming the effect and source table, after
the source read and before any `match()`. This is an existing gap, not one the
views introduce, and it is narrow: it needs duplicates outside the animals
phenotyped at the time, or rows added later. It is a behavior change, so it gets
its own NEWS line. One row per animal is unaffected.

`sql_in_list()` stays: after this change it still has many non-id callers
(trait, locus, chromosome, component, index and owner names in
`add_tgv.R`, `genome_effects_eval.R`, `chr_meta_helpers.R`,
`extract_genetic_variance.R`, `define_genome_effect_terms.R`, and the locus
names in `add_dosage.R` / `extract_genotypes.R`). Its `what = "individual ID"`
uses all go, along with the `"individual ID"` example in its roxygen
(`R/sql_utils.R:550`) and the `id_ind` mention in the header of
`tests/testthat/test-sql_injection_hardening.R`.

Leave 0.76.14's `add_founders()` / `add_offspring()` / `.ap_read_by_id()`
code as it is; migrating them onto the new helper is optional cleanup, not
part of this fix.

## Same statement, separate defect: the repeated-record mean

`build_data_file()` computes `AVG(pheno_value)` across repeated records. That
is a parallel floating aggregate, which CLAUDE.md's bit-identical rule forbids,
and the views do not fix it. The real `build_data_file()` on 0.76.14 shows it,
as does the planned `JOIN` form with plain `AVG()`. Records are
`rnorm(n, 10, 2)`, 5 per animal, stored interleaved by record position; the
counts compare 1 thread against 4 or 8:

| Animals × records | Rows | Means that differ | `data.txt` lines that differ |
|---|---|---|---|
| 20,000 × 5 | 100,000 | 226 | 23 |
| 100,000 × 5 | 500,000 | 24,863–31,834 | 2,519–3,162 |

So `write.table()`'s 15 significant digits do not absorb the drift: about a
tenth of the differing means change the file, and `data.txt` is the solver's
input. (No solver was run, so the effect on stored EBVs is inferred.) The counts
are one seed on one machine, not a threshold. Whether drift appears depends on:

- **Value precision.** R's `runif(n, 0, 100)` never drifted, even at 500,000
  rows, because it carries about 32 random bits: sums of a few such values are
  exact in a double, so order cannot matter. The same values divided by 3 fill
  the mantissa and drift like `rnorm()`. Simulated phenotypes are full-precision.
- **Physical layout.** Contiguous records barely drift (0–2 means at 500,000
  rows); interleaved ones do, because each animal's rows reach different threads.
- **Size.** With interleaved `rnorm()` values, 50,000 rows showed nothing and
  100,000 did.

Fix it in the same rewrite, as its own step and NEWS line, using the fixed-order
pattern `ind_tgv_total` already uses:
`list_avg(list(pheno_value ORDER BY id_phenotype))`. `id_phenotype` is the
primary key, so the order is total. `list_avg()` skips `NULL` and returns `NULL`
for an all-`NULL` group, as `AVG()` does, and was bit-identical across 1, 4 and
8 threads in the same probe. Not `GEV_ACC_TYPE`: its `DOUBLE → DECIMAL → DOUBLE`
round trip moves a single-record value (`R/genome_effects_helpers.R`), and an
accumulator range limit would be a new failure mode here.

## Tests

**Statement capture.** Reuse the pattern from `test-individual_ids.R` ("no
individual id appears in the SQL the writers execute": mocked `DBI`
`dbGetQuery` / `dbExecute` / `dbSendQuery` forwarding to the real connection,
ids from a distinctive line name such as `Qzline_F1`, asserting no statement
matches `Qzline_F1_[0-9]`).

- Select through non-id predicates (`sex`, `gen`, a fixture column). An explicit
  `filter(id_ind %in% ...)` renders the ids and fails by design (see Scope).
  Internal helpers can take distinctive id vectors directly.
- Keep fixture setup outside the captured operation, so test inserts do not
  contaminate the statement list.
- Assert the intended read or mutation actually ran, as the existing test does
  for both parent reads. No id in SQL proves nothing if the query never ran.

**Results.** Capture alone cannot catch the wrong animals being selected. Use
small fixtures with deliberately different dosages, EBVs and fixed-effect
values:

| Path | Required coverage and result assertion |
|---|---|
| `add_dosage()` | Prefill the cache; overwrite a proper subset with fewer loci. Selected animals lose their old cache rows, others keep theirs, new values equal the haplotype sums. Exercise unrestricted overwrite too. Regression: select from `ind_genotype` itself with `overwrite_dosage = TRUE`. |
| `add_genotypes()` | Mark one subset, then another. Only the requested animals change; existing `TRUE` flags stay. |
| `extract_genotypes()` | Chip/subset intersection, and the non-chip `loci_tbl` or `effects_tbl` path. Exact ids, locus columns, integer dosages, and row order. |
| `assert_ploidy_2()` | A bad-ploidy animal outside the set is ignored; inside it, the existing error. Cover the unrestricted branch. |
| `ebv_parent_avg()` | Offspring with two known parents and stored parental EBVs, including several evaluations and another model; assert the latest-within-model average. A founders-only fixture never reaches the parent-EBV query. |
| `trace_pedigree()` | Three generations and a shared ancestor. Depths 0, 1 and 2, no duplicates, `NULL` parent to `"0"`, and row order. |
| `build_data_file()` | A fixed class and a fixed covariate, so the source-table read runs. Selected animals, effect values, column mapping, and a nonempty `pheno_ids` restriction. A `source_table` with two rows for a selected animal is refused, naming the effect and table; duplicates only among unselected animals are not. |
| `write_geno_file()` | Only some pedigree animals genotyped, and a proper chip locus subset. Fixed-width ids, exact dosage strings and row order in the file; also the no-genotyped-animal return. |

No solver is needed: `ebv_parent_avg()` runs through
`add_ebv(..., parent_avg = TRUE)`, and the three file builders only read the
database and write to a temporary `eval_dir`.

Row-order tests assert the deliberate order directly (ids equal their
`sort(method = "radix")`). A small fixture does not reorder (one order across 20
repeats of a 50-id join), so comparing two runs would pass on the broken code.

**Helper** (`.with_id_view()`):

- Repeated ids against a real aggregate, not just the frame's row count; empty
  character input; `NULL` and `NA` refused.
- Lifetime: outer call reads its ids; inner call with different ids reads its
  own; the outer view still works after return. Repeat with an inner error
  caught by the outer callback: the inner view is gone, the outer usable. Let an
  error escape: neither view remains, and a later call on the same connection
  succeeds. An unrelated registered view survives all of it.
- `.Random.seed` unchanged after success and after failure.
- An error from `fun()` reaches the caller unchanged, not replaced by a cleanup
  error.

**Mean determinism.** Build `ind_phenotype` by direct insert and run
`build_data_file()` under `SET threads` of 1 and 4. Start from the fixture
measured above: 100,000 animals × 5 records, seeded `rnorm(n, 10, 2)`, stored
interleaved by record position (`rep(ids, times = 5)`), about 1.4 s per build.
Do not use raw `runif()`, which cannot fail, or records contiguous by animal.
Row count alone proves nothing, so the fixture is a starting point: the test is
valid only once it has been seen to fail with plain `AVG()` restored, per
`test-genome-effects-determinism.R`. That check is about the algorithm's
determinism. Do not compare against stored pre-change output.

Assert `expect_identical()` on the returned `data` and on the `data.txt` lines
across thread counts, so the solver's input is covered as well as the in-memory
result.

**Mean edge cases.** A small separate fixture, through `build_data_file()`:

- a singleton `0.30000000000000004` stays exact, asserted on `data`, because
  `data.txt` rounds it to `0.3`;
- a mixed `NULL` / value group averages only its values (`NULL, 2, 4` gives 3);
- an all-`NULL` group is `NULL` from SQL and becomes the builder's missing value
  0 (`NA` to 0 after the merge).

**Regression set.** `test-individual_ids.R`, `test-resolve_subset_ids.R`,
`test-add_dosage.R`, `test-add_genotypes.R`, `test-extract_genotypes.R`,
`test-ploidy_helpers.R`, `test-add_ebv.R`, `test-dosage_extract_parity.R`,
`test-sql_injection_hardening.R`, plus the new tests; then the full suite, since
this crosses modules. Existing parity and reproducibility tests stay unchanged
(within current code). A `covr` run is optional.

## Docs

- `tidybreed-api` skill: the `.with_id_view()` contract, and that
  `build_data_file()`'s repeated-record mean is deliberately fixed-order, so it
  is not simplified back to `AVG()`, and that it refuses duplicate source rows.
- CLAUDE.md hard rule **IDs**: name `.with_id_view()` and state the user-filter
  exception, so the rule stops claiming more than the code does.
- Regenerate `man/`: `man/sql_in_list.Rd` carries the `"individual ID"` example,
  and the `build_data_file()` page should state the duplicate refusal.
- No schema-skill change: no table, column or invariant moves.
- `NEWS.md` and `DESCRIPTION` before commit and push.

## Out of scope

Performance claims (would need `dev/benchmarks/` first); changing the id format
or numbering; schema migration; user-filter SQL rendering (see Scope).
