# Stage 3: write the resolved call in one transaction

`duckdb_register()` + `INSERT` only — never `dbWriteTable()`, which
advances R's RNG. Any failure rolls back both tables. The per-phenotype
"Wrote ..." messages are emitted after the commit.

## Usage

``` r
.ap_commit(pop, plan, resolved, extra_cols = list())
```

## Arguments

- plan:

  The Stage-1 plan (for the path of each phenotype).

- resolved:

  The Stage-2 result.

- extra_cols:

  Scalar extra columns for `ind_phenotype`.
