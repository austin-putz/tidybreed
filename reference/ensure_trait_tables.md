# Ensure the trait-layer tables exist in the database

Creates all core trait / phenotype / individual tables if they do not
already exist. Idempotent. Pre-1.0.0 there is no cross-version
compatibility contract, so this creates the current schema only and
never migrates an older database.

## Usage

``` r
ensure_trait_tables(pop)
```

## Arguments

- pop:

  A `tidybreed_pop` object.

## Value

The `tidybreed_pop` object with `$tables` updated.
