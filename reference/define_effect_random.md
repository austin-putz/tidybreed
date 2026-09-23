# Define a random group effect in a phenotype model

Inserts a row into `phenotype_effects` for a random effect. One value is
drawn per distinct level of `source_column`; all individuals sharing
that level receive the same shift. Drawn values are stored in
`phenotype_random_effects` so they are reproducible across repeated
calls to
[`add_phenotype()`](https://austin-putz.github.io/tidybreed/reference/add_phenotype.md)
without requiring a fixed `seed`.

**A level's draw is persistent.** The realized value for pen `P1`
applies to *every* individual that is ever in `P1` — across batches and
seasons, in every later
[`add_phenotype()`](https://austin-putz.github.io/tidybreed/reference/add_phenotype.md)
call, forever. That is the model declared by using `pen_id` as the
level. An effect that should be re-realized per batch is a **different
level** — `pen_batch_id`, or an `interaction(pen_id, batch_id)` column
written to `ind_meta` — not a different feature. A level is drawn the
first time a planned record touches it; a level touched only by
individuals that end up without a record is never drawn.

To correlate this effect across multiple phenotypes (e.g. the same herd
affects both ADG and BW), call
[`define_effect_cov_matrix()`](https://austin-putz.github.io/tidybreed/reference/define_effect_cov_matrix.md)
with the appropriate `effect_name` — either before or after this call.
Once the phenotype belongs to a block of two or more phenotypes for
`effect_name`, this call must use `distribution = "normal"` and the same
`(source_column, source_table)` as the block's other members, and
`variance` can no longer be set here — the block is redeclared as a
whole with
[`define_effect_cov_matrix()`](https://austin-putz.github.io/tidybreed/reference/define_effect_cov_matrix.md).
Within such a block a level's draw for one phenotype is conditional on
the draws it already has stored for the block's other phenotypes,
whichever phenotype was generated first and however many calls apart:
`add_phenotype("ADG")` today and `add_phenotype("BF")` next season gives
pen `P1` a `(ADG, BF)` pair with the declared covariance. A block member
that is not in a call, or that has no random term for the effect, is
simply not drawn — its coordinate stays latent until it is needed. A
`"gamma"` or `"uniform"` effect is supported only while its phenotype is
alone in its block.

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

  Numeric scalar or `NULL`. Variance of the random effect. `NULL`
  (default) uses the value already stored in `phenotype_var_comp` via
  [`define_effect_cov_matrix()`](https://austin-putz.github.io/tidybreed/reference/define_effect_cov_matrix.md)
  and errors if there is none. A number writes (or overwrites) a 1 × 1
  block for this phenotype; it is an error when the phenotype is already
  in a multi-phenotype block for `effect_name`.

- distribution:

  Character. Sampling distribution: `"normal"` (default), `"gamma"`
  (shape 1, rate `1 / sqrt(variance)`), or `"uniform"` (on
  `± sqrt(3 * variance)`). The last two are marginal samplers for a
  phenotype alone in its block; a block of two or more requires
  `"normal"`.

- source_table:

  Character. Table containing `source_column`. Default `"ind_meta"`.

- overwrite:

  Logical. Replace an existing effect with the same name. The stored
  draws of that effect for this phenotype in `phenotype_random_effects`
  are discarded with it.

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
