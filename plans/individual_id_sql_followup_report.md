# Report: individual ids out of SQL text in the remaining modules

Implements `plans/individual_id_sql_followup.md`. Version 0.76.15, 2026-11-04.
Not committed.

## Outcome

All six modules now reach SQL through the new `.with_id_view()`
(`R/sql_utils.R`). No resolved id vector is pasted into SQL in `add_dosage()`,
`add_genotypes()`, `extract_genotypes()`, `assert_ploidy_2()`, `add_ebv()` or
the BLUPF90 builders, and `sql_in_list()` has no individual-id caller left.

Three behavior changes come with it, each with its own NEWS line:

- `extract_genotypes()` rows, `pedigree.txt` and `genotype.txt` are sorted by
  `id_ind`.
- `build_data_file()`'s repeated-record mean is
  `list_avg(list(pheno_value ORDER BY id_phenotype))`, not `AVG()`.
- `build_data_file()` refuses duplicate `source_table` rows for animals it reads.

## Final plan check

The design held; nothing in it changed. Checked against the code:

- **Inventory:** complete. The plan's regex was only one way to look, so the
  check also covered `IN (` builders with variable columns, hand-rolled quoted
  lists, `sql_in_list()`, `format_sql_value()` / `dbQuoteLiteral()` callers and
  `dplyr::filter()` in package code. The one extra hit, `R/remove_rows.R:144`,
  is a doc comment; that code uses a temporary table.
- **Transactions:** none of the six modules opens one.
- **Unique keys:** `parent_ebv` is unique through
  `UNIQUE (id_ind, trait_name, model, eval_number)`, and `id_phenotype` is the
  primary key, so the plan's "matched back by key" claim holds.
- **DuckDB 1.5.5 probes:** an empty frame registers as a `VARCHAR` view,
  re-registering a live name errors, and `list_avg()` gives the planned `NULL`,
  mixed and singleton results.

## What differed from the plan

| Plan | Done | Why |
|---|---|---|
| Reuse the mocked-`DBI` capture from `test-individual_ids.R` | Used `record_sql()` from `tests/testthat/helper-sql.R` | It already exists and also sees dbplyr's statements. The plan did not mention it. The tests assert `DELETE`, `UPDATE` and reads were captured. |
| `JOIN` or `IN (SELECT ...)` per site | `IN (SELECT ...)` at every site | One form, closest to the old statements. `JOIN` remains documented for later callers. |
| (not in plan) | Added a test that a failed registration leaves a view it does not own alone | The contract states it; the mutation below shows it is otherwise unpinned. |

`man/` regeneration changed exactly the three expected pages
(`sql_in_list`, `build_data_file`, `assert_ploidy_2`). `DESCRIPTION` has no
`RoxygenNote`; roxygen was run in a scratch copy and those three files copied
back.

## Verification

- **Full suite:** 1,223 tests, 5,205 expectations, 0 failed, 0 errors, 0
  skipped. It ran on the final code and tests. Only docs, `man/`, `NEWS.md` and
  the version changed after it.
- **New file:** `tests/testthat/test-individual_id_views.R`, 123 expectations.
- **Mean determinism:** the plan's fixture (100,000 animals x 5 interleaved
  `rnorm()` records) failed with `AVG()` restored, on both `data` and the
  `data.txt` lines, and passes with `list_avg()`. This was one seed on one
  machine. The plan's drift counts were not re-measured.

**Mutation check.** Each defect was put into a scratch copy, one at a time, and
the new test file was run. All 14 were caught:

| Defect | Caught by |
|---|---|
| `AVG()` restored | mean determinism; fixed-effects read |
| No `unique()` in the helper | set semantics |
| Duplicate refusal removed | duplicate refusal |
| No `ORDER BY`: `trace_pedigree()` final read | pedigree |
| No `ORDER BY`: `geno_ids` | `write_geno_file()` |
| No `ORDER BY`: filtered chip read | `extract_genotypes()` chip and filter |
| No `ORDER BY`: unfiltered chip read | `extract_genotypes()` unfiltered paths |
| Pedigree view not refreshed per generation | pedigree |
| `add_dosage()` selection as a live subquery | `ind_genotype` self-selection |
| Literal ids in `add_genotypes()` | capture |
| Literal ids in `assert_ploidy_2()` | capture |
| `parent_ids` read with the wrong id set | parent average |
| Fixed-effect read with the wrong id set | `build_data_file()` |
| `on.exit()` installed before registration | failed registration |

The unfiltered `SELECT DISTINCT id_ind FROM ind_haplotype` ordering was not
mutated; a test asserts the statement text and the result order.

**Not run:** BLUPF90 itself (no solver; the builders were exercised directly),
`R CMD check`, and `covr` (optional in the plan).

## Added after review: empty `phenotype` filter, per-trait counts

Found during the check and present before this work: in `build_data_file()`,
`pheno_ids = integer(0)` added no restriction, so
`add_ebv(phenotype = <filter matching nothing>)` warned that the "data file may
be empty" and then evaluated on **every** record. Decided with the user (option
D) and fixed in the same 0.76.15:

- `add_ebv()` now stops during validation: "`phenotype` filter matched no
  records in ind_phenotype; nothing to evaluate." Nothing is written.
- `build_data_file()` treats a non-`NULL` empty `pheno_ids` as no records (the
  existing "No phenotypic records found" error), as a backstop for any caller.
- `build_data_file()` prints `Records [<trait>]: <n> animals` per trait,
  counting animals with a non-missing value. A trait with 0 (e.g. litter size
  early in a simulation) is allowed and labelled: in a multi-trait model its
  EBVs come through the genetic correlations. This is a `message()`, not a
  `warning()`, because it is expected, not a likely mistake. Every trait at 0
  is the "No phenotypic records found" error, so an all-missing data file no
  longer reaches the solver.

Five new tests (14 expectations) in `test-individual_id_views.R`. Mutation
check: restoring the warning, the `length(pheno_ids) > 0` condition, a row
count instead of a non-missing count, or dropping the every-trait-zero stop
each fails at least one of them.

The BLUPF90 missing-value coding found along the way (`OPTION missing 0`, which
covers responses only) is left for the `add_ebv()` redesign:
`plans/blupf90_missing_values.md`.

## Smaller notes

- `ebv_parent_avg()`'s "multiple evaluations" warnings follow the row order of
  `parent_ebv`, so with several traits their order can vary. Values are matched
  by key and do not vary. Left as is.
- `eval_number` counts across models within a trait: stored evaluations 1-3
  under two models make the next one 4, whichever model it uses. Existing
  behavior; the parent-average test records it.
- `skip_on_cran()` was tried on the determinism test and removed: it skipped the
  test whenever `NOT_CRAN` was unset. The test adds a few seconds per run.
