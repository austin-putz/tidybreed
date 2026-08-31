# Define a genetic component trait

Creates one row in `trait_meta` describing a **genetic component
trait**: a quantity with QTL effects in `genome_effects`, TBVs in
`ind_tbv`, and additive genetic variance in `trait_var_comp`. Contains
no phenotype-level information.

To register the **observed phenotype** that individuals receive records
for, call
[`define_phenotype()`](https://austin-putz.github.io/tidybreed/reference/define_phenotype.md)
after this function. For the common one-off case use
[`define_trait_simple()`](https://austin-putz.github.io/tidybreed/reference/define_trait_simple.md),
which chains both steps together.

## Usage

``` r
define_trait(
  pop,
  trait_name,
  target_add_var = NULL,
  target_add_mean = 0,
  expressed_parent = c("both", "parent_1", "parent_2"),
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

- target_add_var:

  Numeric. Target additive genetic variance. Written to `trait_var_comp`
  as a diagonal entry under `effect_name = "gen_add"`. Used by
  [`define_additive_effects()`](https://austin-putz.github.io/tidybreed/reference/define_additive_effects.md)
  to rescale effects. If already set via
  [`define_effect_cov_matrix()`](https://austin-putz.github.io/tidybreed/reference/define_effect_cov_matrix.md),
  leave `NULL`.

- target_add_mean:

  Numeric. TBV centering mean for the base population. Default `0`;
  `E[TBV] = 0` when TBVs are centered on base allele frequencies. The
  phenotypic population mean (intercept) is set separately in
  [`define_phenotype()`](https://austin-putz.github.io/tidybreed/reference/define_phenotype.md).

- expressed_parent:

  Character. Parent-of-origin expression: `"both"` (default),
  `"parent_1"` (paternal), or `"parent_2"` (maternal). Imprinted traits
  use only the haplotype from the specified parent when computing TBVs.

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
[`define_trait_simple()`](https://austin-putz.github.io/tidybreed/reference/define_trait_simple.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# Simple genetic component trait:
pop <- pop |>
  define_trait("ADG", target_add_var = 100, units = "g/day")

# Maternal component traits (no define_phenotype call needed for WWD/WWM):
pop <- pop |>
  define_trait("WWD", target_add_var = 200) |>
  define_trait("WWM", target_add_var = 80)

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
