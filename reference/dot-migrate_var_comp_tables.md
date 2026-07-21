# Migrate old variance-component table names to the v0.42.0 naming convention

Renames `trait_effect_cov` → `trait_var_comp` and
`phenotype_residual_cov` → `phenotype_var_comp` in an existing database,
adding the `effect_name` column (backfilled as 'residual') when needed.
Called automatically by
[`restore_pop()`](https://austin-putz.github.io/tidybreed/reference/restore_pop.md).

## Usage

``` r
.migrate_var_comp_tables(con)
```

## Arguments

- con:

  An active DuckDB connection.
