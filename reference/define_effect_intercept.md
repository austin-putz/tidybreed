# Set the intercept (population mean) for a phenotype

Updates `mean` in `phenotype_meta` for the named phenotype. This value
is the overall population mean added to every individual's liability
before residuals and fixed/random shifts. It can also be set at
definition time via the `mean` argument in
[`define_phenotype()`](https://austin-putz.github.io/tidybreed/reference/define_phenotype.md);
`define_effect_intercept()` lets you update it without redefining the
whole phenotype.

## Usage

``` r
define_effect_intercept(pop, phenotype_name, mean)
```

## Arguments

- pop:

  A `tidybreed_pop` object.

- phenotype_name:

  Character. Name of an existing phenotype in `phenotype_meta`.

- mean:

  Numeric scalar. The new population mean / intercept value.

## Value

The modified `tidybreed_pop` (invisibly).

## See also

[`define_phenotype()`](https://austin-putz.github.io/tidybreed/reference/define_phenotype.md),
[`define_effect_fixed_class()`](https://austin-putz.github.io/tidybreed/reference/define_effect_fixed_class.md),
[`define_effect_fixed_cov()`](https://austin-putz.github.io/tidybreed/reference/define_effect_fixed_cov.md),
[`define_effect_random()`](https://austin-putz.github.io/tidybreed/reference/define_effect_random.md)

## Examples

``` r
if (FALSE) { # \dontrun{
pop <- pop |> define_effect_intercept("ADG", mean = 850)
} # }
```
