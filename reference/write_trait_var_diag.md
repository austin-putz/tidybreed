# Write a single diagonal entry to trait_var_comp

Used internally by define_trait() to write a per-trait variance as a 1x1
diagonal entry. Uses dbExecute() to avoid consuming R's RNG.

## Usage

``` r
write_trait_var_diag(pop, effect_name, trait_name, variance)
```

## Arguments

- pop:

  A `tidybreed_pop` object.

- effect_name:

  Character.

- trait_name:

  Character.

- variance:

  Numeric scalar.

## Value

The modified `tidybreed_pop` (invisibly).
