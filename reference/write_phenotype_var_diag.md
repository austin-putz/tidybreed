# Write a single diagonal entry to phenotype_var_comp

Used internally by define_effect_random() to write a per-phenotype
variance. Uses dbExecute() to avoid consuming R's RNG.

## Usage

``` r
write_phenotype_var_diag(pop, effect_name, phenotype_name, variance)
```

## Arguments

- pop:

  A `tidybreed_pop` object.

- effect_name:

  Character.

- phenotype_name:

  Character.

- variance:

  Numeric scalar.

## Value

The modified `tidybreed_pop` (invisibly).
