# Plan (follow-up): individual ids out of SQL text in the remaining modules

Status: draft, not implemented. Follows `plans/individual_id_sql.md` (D2 = (a)),
which fixed `add_founders()` and `add_offspring()` in 0.76.14.

## Why

CLAUDE.md: "Individual ids never appear in SQL text — register them as a view."
These sites still paste ids into `IN (...)` lists. None takes ids from user
input directly (they come from `resolve_subset_ids()` or a table, so they are
tidybreed-generated `{line}_{n}` strings with no quote), so they are rule
breaches, not live bugs. Two files do not even escape.

## Inventory (re-checked on the 0.76.14 working tree)

| File | Lines | What | Quoting |
|---|---|---|---|
| `R/add_dosage.R` | 100–101, 109–111 | read `ind_haplotype`; `DELETE FROM ind_genotype` | `sql_in_list()` |
| `R/add_genotypes.R` | 79–82 | `UPDATE ind_meta SET <chip> = TRUE` | `sql_in_list()` |
| `R/extract_genotypes.R` | 163–167, 244–251 | genotype reads | `sql_in_list()` |
| `R/ploidy_helpers.R` | 20 | `assert_ploidy_2()`; callers `add_dosage.R:77`, `extract_genotypes.R:189` | `sql_in_list()` |
| `R/add_ebv.R` | 291–295, 302–317 | subset read; `ebv_parent_avg()` parent read | **none** |
| `R/blupf90_helpers.R` | 34–38, 45–49, 118–127, 156–160, ~210–236 | pedigree, phenotypes, fixed-effect source tables, genotyped animals | **none** |

Not individual ids, out of scope here: trait / phenotype names in `add_ebv.R`
(206, 238–240, 303) and `blupf90_helpers.R` (81–88) — these are validated
names and could move to `DBI::dbQuoteLiteral()`; integer `id_phenotype` and
`locus_id` lists.

## Proposed change

One internal helper generalising `.ap_read_by_id()`
(`R/add_phenotype_stages.R`):

```r
# Run `fun(view_name)` with `ids` registered as a one-column view `id_ind`,
# unregistered on exit (success or error). RNG-neutral (register only).
.with_id_view <- function(conn, ids, fun, view = "__id_view")
```

Each site becomes a `JOIN <view> USING (id_ind)` or
`id_ind IN (SELECT id_ind FROM <view>)`. `UPDATE` / `DELETE` use the
`IN (SELECT ...)` form. Distinct view names where calls nest
(`assert_ploidy_2()` inside `add_dosage()`).

`sql_in_list()` may then have no individual-id callers left; keep it only if a
non-id caller remains (locus names), else remove it.

## Tests

- Extend the statement-capture test from `test-individual_ids.R` (mocked DBI
  bindings forwarding to the real connection) to `add_dosage()`,
  `add_genotypes()`, `extract_genotypes()`: no individual id in any statement.
- `ebv_parent_avg()`, pedigree collection and BLUPF90 input-file preparation
  run on small database fixtures without the BLUPF90 binary; cover them with
  the same capture (coverage of these files is low only because the solver
  call is not run).
- View cleanup after success and error; `.Random.seed` unchanged.
- Existing parity / reproducibility tests unchanged (within current code).

## Out of scope

Performance claims (would need `dev/benchmarks/` first); changing the id format.
