# Define a genetic component trait

Creates one row in `trait_meta` describing a **genetic component
trait**: a quantity with QTL effects in `genome_effects`, TBVs in
`ind_tbv`, and additive genetic variance in `trait_var_comp`. Contains
no phenotype-level information.

To register the **observed phenotype** that individuals receive records
for, call
[`define_phenotype()`](https://austin-putz.github.io/tidybreed/reference/define_phenotype.md)
after this function.

The trait's genetic **targets** are not set here. Store them with
[`define_effect_cov_matrix()`](https://austin-putz.github.io/tidybreed/reference/define_effect_cov_matrix.md)
(`effect_name = "additive"`), or pass `G` to
[`define_additive_effects()`](https://austin-putz.github.io/tidybreed/reference/define_additive_effects.md),
which writes the target together with the effects calibrated to it. The
usual chain is `define_trait()` → `define_additive_effects(G = )` →
[`define_phenotype()`](https://austin-putz.github.io/tidybreed/reference/define_phenotype.md).

## Usage

``` r
define_trait(
  pop,
  trait_name,
  description = NULL,
  units = NULL,
  overwrite = FALSE
)
```

## Arguments

- pop:

  A `tidybreed_pop` object.

- trait_name:

  Character. Unique identifier for this genetic component trait. Must be
  a valid SQL identifier.

- description:

  Character. Free-text description of the trait.

- units:

  Character. Measurement units, e.g. `"kg"`, `"count"`.

- overwrite:

  Logical. If `TRUE` and a trait with the same name already exists,
  replace its `trait_meta` row and clear associated `phenotype_effects`
  rows. Default `FALSE` errors if the trait already exists.

## Value

The modified `tidybreed_pop` (invisibly).

## See also

[`define_phenotype()`](https://austin-putz.github.io/tidybreed/reference/define_phenotype.md),
[`define_additive_effects()`](https://austin-putz.github.io/tidybreed/reference/define_additive_effects.md),
[`define_effect_fixed_class()`](https://austin-putz.github.io/tidybreed/reference/define_effect_fixed_class.md),
[`define_effect_fixed_cov()`](https://austin-putz.github.io/tidybreed/reference/define_effect_fixed_cov.md),
[`define_effect_random()`](https://austin-putz.github.io/tidybreed/reference/define_effect_random.md),
[`add_phenotype()`](https://austin-putz.github.io/tidybreed/reference/add_phenotype.md),
[`define_effect_cov_matrix()`](https://austin-putz.github.io/tidybreed/reference/define_effect_cov_matrix.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# Simple genetic component trait, its target and its effects:
pop <- pop |>
  define_trait("ADG", units = "g/day")
pop <- get_table(pop, "genome_meta") |>
  define_additive_effects("ADG", G = 100)

# Maternal component traits (no define_phenotype call needed for WWD/WWM):
pop <- pop |>
  define_trait("WWD") |>
  define_trait("WWM")
pop <- get_table(pop, "genome_meta") |>
  define_additive_effects(c("WWD", "WWM"),
    G = matrix(c(200, -40, -40, 80), 2, 2))

# Then define the observed composite phenotype:
pop <- pop |>
  define_phenotype("WW", type = "continuous", mean = 230,
    residual_var = 180,
    components = tibble::tribble(
      ~source_trait_name, ~contributor_type,
      "WWD", "self",
      "WWM", "dam"
    ))
} # }
```
