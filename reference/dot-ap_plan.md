# Stage 1: plan every record of an `add_phenotype()` call

Stage 1: plan every record of an
[`add_phenotype()`](https://austin-putz.github.io/tidybreed/reference/add_phenotype.md)
call

## Usage

``` r
.ap_plan(tbl, phenos, user_values = NULL)
```

## Arguments

- tbl:

  The `tidybreed_table` passed to
  [`add_phenotype()`](https://austin-putz.github.io/tidybreed/reference/add_phenotype.md).

- phenos:

  Character vector of validated phenotype names.

- user_values:

  The `user_values` argument, or `NULL`.

## Value

`NULL` when no individual matched (a warning is issued), otherwise a
list with `pop` (after
[`add_tbv()`](https://austin-putz.github.io/tidybreed/reference/add_tbv.md)),
`phenos` (in evaluation order), `pheno_meta` (rows in that order) and
`entries`: one list per phenotype, in the same order, each with

- `phenotype_name`:

- `path`:

  `"model"`, `"derived_formula"` or `"user_values"`.

- `id_ind`, `pheno_number`:

  The planned records, `id_ind`-ordered. `pheno_number` is the value
  Stage 3 writes.

- `tbv`:

  Numeric per record (`"model"` path only).

- `fixed`:

  Fixed-effect contribution per record (`"model"` path; zeros for
  `"derived_formula"`, which has no model terms).

- `random`:

  Random-effect terms: a list of
  `list(effect_name, distribution, level)` where `level` is the
  character level of the grouping column per record (`NA` = none).

- `condition_table`, `condition_column`, `condition_value`:

  The residual stratum lookup: `NULL` when the phenotype's residual
  block is unconditional or absent, otherwise the raw condition value
  per record as character (`NA` = `NULL` in the table).

- `user_values`:

  Numeric per record (`"user_values"` path).

- `formula`:

  The derived formula string (`"derived_formula"` path).

and the covariance blocks Stage 2 draws through, from
[`find_covariance_blocks()`](https://austin-putz.github.io/tidybreed/reference/find_covariance_blocks.md):
`residual_blocks` (the residual blocks touching the call's phenotypes),
`named_targets` (per random effect, the model-path phenotypes whose
planned records carry it; see
[`.ap_named_effect_targets()`](https://austin-putz.github.io/tidybreed/reference/dot-ap_named_effect_targets.md))
and `named_blocks` (per effect, the blocks touching those phenotypes).
