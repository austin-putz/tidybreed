# Define a random group effect in a phenotype model

Inserts a row into `trait_effects` for a random effect. One value is
drawn per distinct level of `source_column`; all individuals sharing
that level receive the same shift. Drawn values are stored in
`trait_random_effects` so they are reproducible across repeated calls to
[`add_phenotype()`](https://austin-putz.github.io/tidybreed/reference/add_phenotype.md)
without requiring a fixed `seed`.

To correlate this effect across multiple phenotypes (e.g. the same herd
affects both ADG and BW), call
[`define_effect_cov_matrix()`](https://austin-putz.github.io/tidybreed/reference/define_effect_cov_matrix.md)
with the appropriate `effect_name` — either before or after this call.

## Usage

``` r
define_effect_random(
  pop,
  phenotype_name,
  effect_name,
  source_column,
  variance = NULL,
  distribution = c("normal", "gamma", "uniform"),
  source_table = "ind_meta",
  overwrite = FALSE
)
```

## Arguments

- pop:

  A `tidybreed_pop` object.

- phenotype_name:

  Character. Name of an existing phenotype in `phenotype_meta`.

- effect_name:

  Character. Unique label for this effect within the phenotype.

- source_column:

  Character. Column in `source_table` whose distinct values define the
  groups (e.g. `"herd_id"`, `"litter"`, `"id_ind"` for PE).

- variance:

  Numeric scalar. Variance of the random effect. Optional if already
  stored via
  [`define_effect_cov_matrix()`](https://austin-putz.github.io/tidybreed/reference/define_effect_cov_matrix.md);
  required otherwise.

- distribution:

  Character. Sampling distribution: `"normal"` (default), `"gamma"`, or
  `"uniform"`.

- source_table:

  Character. Table containing `source_column`. Default `"ind_meta"`.

- overwrite:

  Logical. Replace an existing effect with the same name.

## Value

The modified `tidybreed_pop` (invisibly).

## See also

[`define_effect_cov_matrix()`](https://austin-putz.github.io/tidybreed/reference/define_effect_cov_matrix.md),
[`define_effect_fixed_class()`](https://austin-putz.github.io/tidybreed/reference/define_effect_fixed_class.md),
[`define_effect_fixed_cov()`](https://austin-putz.github.io/tidybreed/reference/define_effect_fixed_cov.md),
[`add_phenotype()`](https://austin-putz.github.io/tidybreed/reference/add_phenotype.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# Herd random effect
pop <- pop |>
  define_effect_random("ADG", "herd",
    source_column = "herd_id",
    variance = 150)

# Permanent environment (PE) for repeatability — one draw per animal
pop <- pop |>
  define_effect_random("litter_size", "pe",
    source_column = "id_ind",
    variance = 0.3)
} # }
```
