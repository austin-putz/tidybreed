# Ensure the trait-layer tables exist in the database

Creates all core trait / phenotype / individual tables if they do not
already exist, and runs any necessary schema migrations on older
databases. Idempotent.

## Usage

``` r
ensure_trait_tables(pop)
```

## Arguments

- pop:

  A `tidybreed_pop` object.

## Value

The `tidybreed_pop` object with `$tables` updated.
