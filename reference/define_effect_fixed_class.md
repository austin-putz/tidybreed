# Define a discrete fixed-class effect in a phenotype model

Inserts a row into `phenotype_effects` representing a discrete fixed
effect: each level of `source_column` maps to a numeric shift on the
phenotype scale. The `null_class_action` argument controls what happens
when an individual has a `NULL` value for `source_column` at phenotyping
time.

## Usage

``` r
define_effect_fixed_class(
  pop,
  phenotype_name,
  effect_name,
  source_column,
  levels,
  source_table = "ind_meta",
  null_class_action = c("skip", "error", "zero"),
  overwrite = FALSE
)
```

## Arguments

- pop:

  A `tidybreed_pop` object.

- phenotype_name:

  Character. Name of an existing phenotype in `phenotype_meta`.

- effect_name:

  Character. Unique label for this effect within the phenotype. Must be
  a valid SQL identifier.

- source_column:

  Character. Column in `source_table` holding the group membership. Must
  be a valid SQL identifier.

- levels:

  Named numeric vector mapping level strings to shifts, e.g.
  `c(M = 30, F = 0)`. Every level present in the data at phenotyping
  time must appear here.

- source_table:

  Character. Database table containing `source_column`. Default
  `"ind_meta"`.

- null_class_action:

  Character. Behaviour when `source_column` is `NULL` for an individual:

  - `"skip"` (default) — individual is excluded from `ind_phenotype` for
    this phenotype; a warning lists the count.

  - `"error"` — hard stop when any `NULL` is encountered. Use for
    effects like `sex` that should never be missing.

  - `"zero"` — treat `NULL` as a zero-shift contribution (reference
    level).

- overwrite:

  Logical. Replace an existing effect with the same name.

## Value

The modified `tidybreed_pop` (invisibly).

## See also

[`define_effect_intercept()`](https://austin-putz.github.io/tidybreed/reference/define_effect_intercept.md),
[`define_effect_fixed_cov()`](https://austin-putz.github.io/tidybreed/reference/define_effect_fixed_cov.md),
[`define_effect_random()`](https://austin-putz.github.io/tidybreed/reference/define_effect_random.md),
[`add_phenotype()`](https://austin-putz.github.io/tidybreed/reference/add_phenotype.md)

## Examples

``` r
if (FALSE) { # \dontrun{
pop <- pop |>
  define_effect_fixed_class("ADG", "sex",
    source_column = "sex",
    levels = c(M = 30, F = 0))
} # }
```
