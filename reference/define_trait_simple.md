# Define a trait with QTL and sampled effects in one call

Convenience wrapper that chains
[`define_trait()`](https://austin-putz.github.io/tidybreed/reference/define_trait.md),
[`define_additive_effects()`](https://austin-putz.github.io/tidybreed/reference/define_additive_effects.md),
and
[`define_phenotype()`](https://austin-putz.github.io/tidybreed/reference/define_phenotype.md)
for a single uncorrelated trait using random QTL placement. For
correlated multi-trait simulations or custom QTL placement, use the
functions individually with
`get_table("genome_meta") |> filter(...) |> define_additive_effects()`.

## Usage

``` r
define_trait_simple(
  pop,
  trait_name,
  n_qtl,
  target_add_var,
  mean = 0,
  residual_var,
  type = "continuous",
  expressed_sex = "both",
  repeatable = FALSE,
  effect_distribution = "normal",
  scale_to_target = TRUE,
  seed = NULL,
  ...
)
```

## Arguments

- pop:

  A `tidybreed_pop` object.

- trait_name:

  Character. Trait name; must be a valid SQL identifier.

- n_qtl:

  Integer. Number of QTL to randomly select from `genome_meta`.

- target_add_var:

  Numeric. Additive genetic variance target. Passed to
  [`define_trait()`](https://austin-putz.github.io/tidybreed/reference/define_trait.md).

- mean:

  Numeric. Phenotypic population mean. Passed to
  [`define_phenotype()`](https://austin-putz.github.io/tidybreed/reference/define_phenotype.md)
  as the `mean` argument. Default `0`.

- residual_var:

  Numeric. Residual variance. Passed to
  [`define_phenotype()`](https://austin-putz.github.io/tidybreed/reference/define_phenotype.md).

- type:

  Character. One of `"continuous"` (default), `"count"`,
  `"categorical"`. Passed to
  [`define_phenotype()`](https://austin-putz.github.io/tidybreed/reference/define_phenotype.md).

- expressed_sex:

  Character. `"both"` (default), `"M"`, or `"F"`. Passed to
  [`define_phenotype()`](https://austin-putz.github.io/tidybreed/reference/define_phenotype.md).

- repeatable:

  Logical. Whether the trait allows repeated records. Passed to
  [`define_phenotype()`](https://austin-putz.github.io/tidybreed/reference/define_phenotype.md).
  Default `FALSE`.

- effect_distribution:

  Character. `"normal"` (default) or `"gamma"`. Passed to
  [`define_additive_effects()`](https://austin-putz.github.io/tidybreed/reference/define_additive_effects.md).

- scale_to_target:

  Logical. Passed to
  [`define_additive_effects()`](https://austin-putz.github.io/tidybreed/reference/define_additive_effects.md).
  Default `TRUE`.

- seed:

  Optional integer for reproducibility (applied before QTL draw and
  effect sampling).

- ...:

  Additional arguments forwarded to
  [`define_phenotype()`](https://austin-putz.github.io/tidybreed/reference/define_phenotype.md)
  (e.g. `prevalence`, `thresholds`, `cat_values`, `cat_names`,
  `min_value`, `max_value`, `store_liability`).

## Value

The modified `tidybreed_pop` (invisibly).

## See also

[`define_trait()`](https://austin-putz.github.io/tidybreed/reference/define_trait.md),
[`define_additive_effects()`](https://austin-putz.github.io/tidybreed/reference/define_additive_effects.md),
[`define_phenotype()`](https://austin-putz.github.io/tidybreed/reference/define_phenotype.md)

## Examples

``` r
if (FALSE) { # \dontrun{
pop <- pop |>
  define_trait_simple(
    trait_name     = "ADG",
    n_qtl          = 100,
    target_add_var = 100,
    mean           = 850,
    residual_var   = 120
  )
} # }
```
