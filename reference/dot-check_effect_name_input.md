# Refuse a reserved name where a user-named effect is expected

`genetic_ok = TRUE` (only
[`define_effect_cov_matrix()`](https://austin-putz.github.io/tidybreed/reference/define_effect_cov_matrix.md))
lets the genetic names through, because that function routes them to
`trait_var_comp`. Every phenotype-layer definer passes `FALSE`: a random
or fixed effect named `"additive"` would collide with the genetic
vocabulary.

## Usage

``` r
.check_effect_name_input(effect_name, genetic_ok = FALSE)
```
