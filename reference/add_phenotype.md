# Generate phenotype records for a subset of individuals

Simulates phenotype values for one or more phenotypes and writes them to
`ind_phenotype`. Also computes and stores the underlying true breeding
value (TBV) per individual per trait in `ind_tbv`.

**Model** (per phenotype, on the liability / continuous scale):


      y_i = mean + sum(fixed_shifts) + sum(random_shifts) + TBV_i + e_i

- `mean` comes from `phenotype_meta.mean`.

- Fixed and random shifts come from `phenotype_effects` rows.

- For **simple** phenotypes (`phenotype_name == trait_name`), `TBV_i` is
  the standard additive TBV from `genome_effects` (computed via
  [`add_tbv()`](https://austin-putz.github.io/tidybreed/reference/add_tbv.md),
  which this function calls internally for every source trait it needs).

- For **composite** phenotypes (rows in `phenotype_components`, written
  by `define_phenotype(..., components = ...)`), `TBV_i` is the weighted
  sum of contributor TBVs (self, dam, sire, or group) — see
  [`.assemble_composite_tbv()`](https://austin-putz.github.io/tidybreed/reference/dot-assemble_composite_tbv.md).

- For **`formula_tbv`** composite phenotypes
  (`phenotype_meta.formula_tbv` set, written by
  `define_phenotype(..., formula_tbv = ...)`), `TBV_i` is evaluated from
  a small DSL expression referencing self/dam/sire/group TBVs instead of
  a `phenotype_components` data frame.

- `e_i` is residual: drawn from `MVN(0, R)` across phenotypes when a
  covariance matrix is stored in `phenotype_var_comp` and multiple
  phenotypes share the same subset; otherwise drawn independently.

**`derived_formula` phenotypes are the one exception to the model
above.** When `phenotype_meta.type == "derived_formula"`
(`phenotype_meta.formula` set), the phenotype value is computed directly
as an arithmetic expression over other individuals' already-written
`ind_phenotype` records — there is no TBV, no mean/fixed/random
contribution, and no residual draw for that phenotype. When a call mixes
`derived_formula` phenotypes with others that feed them, the phenotypes
are topologically sorted first so dependencies are written before the
formulas that consume them.

**Subset selection**: pipe a `tidybreed_table` (from
[`get_table()`](https://austin-putz.github.io/tidybreed/reference/get_table.md)
and optionally
[`dplyr::filter()`](https://dplyr.tidyverse.org/reference/filter.html))
as the first argument.

**Escape hatches**:

- `user_values`: skip model computation and write these values as
  phenotype records for the subset.

- `user_residual`: supply a numeric vector (or named list) to override
  the residual draw.

## Usage

``` r
add_phenotype(
  tbl,
  phenotype_name = NULL,
  user_residual = NULL,
  user_values = NULL,
  seed = NULL,
  ...
)
```

## Arguments

- tbl:

  A `tidybreed_table` object from
  [`get_table()`](https://austin-putz.github.io/tidybreed/reference/get_table.md)
  (optionally piped through
  [`dplyr::filter()`](https://dplyr.tidyverse.org/reference/filter.html)).
  The table must contain an `id_ind` column.

- phenotype_name:

  Character vector of phenotype name(s). When `NULL` (default), all
  phenotypes in `phenotype_meta` are used in `id_phenotype_meta` order.

- user_residual:

  Optional override for residual draws (skips the `MVN`/independent
  residual sampling step described above; the mean, covariate, and TBV
  contributions are still computed and added). For a single
  `phenotype_name`, a plain numeric vector matched **by position** to
  the final per-phenotype individual list — i.e. after sex-expression
  filtering (`expressed_sex`), repeatable-record exclusion, and any
  missing-component exclusion, so its length must equal the resulting
  subset size, which may be smaller than the original filtered `tbl`.
  For multiple phenotypes, a named list keyed by `phenotype_name`, each
  element following the same positional rule for that phenotype's own
  subset. Unlike `user_values` below, named (per-`id_ind`) vectors are
  **not** supported here.

- user_values:

  Optional override for the full phenotype value — skips the model
  entirely (mean, covariates, and residual are not evaluated), though
  TBVs are still computed and stored in `ind_tbv`. For a single
  `phenotype_name`: a plain numeric vector matching the filtered subset
  by position (`id_ind` order), or a named numeric vector (e.g.
  `c(id_1 = 555, id_2 = 560)`) to match by `id_ind` regardless of order.
  For multiple phenotypes: a named list keyed by `phenotype_name`, each
  element following the same (positional-or-named) rule.

- seed:

  Optional integer for reproducibility.

- ...:

  Optional scalar extra columns written to `ind_phenotype` (broadcast to
  all records). Supply per-record vectors with
  [`mutate_table()`](https://austin-putz.github.io/tidybreed/reference/mutate_table.md)
  after the call.

## Value

The modified `tidybreed_pop` (invisibly).

## See also

[`define_phenotype()`](https://austin-putz.github.io/tidybreed/reference/define_phenotype.md),
[`define_trait()`](https://austin-putz.github.io/tidybreed/reference/define_trait.md),
[`define_additive_effects()`](https://austin-putz.github.io/tidybreed/reference/define_additive_effects.md),
[`define_residual_cov()`](https://austin-putz.github.io/tidybreed/reference/define_residual_cov.md),
[`define_effect_cov_matrix()`](https://austin-putz.github.io/tidybreed/reference/define_effect_cov_matrix.md),
[`define_effect_fixed_class()`](https://austin-putz.github.io/tidybreed/reference/define_effect_fixed_class.md),
[`define_effect_fixed_cov()`](https://austin-putz.github.io/tidybreed/reference/define_effect_fixed_cov.md),
[`define_effect_random()`](https://austin-putz.github.io/tidybreed/reference/define_effect_random.md),
[`add_tbv()`](https://austin-putz.github.io/tidybreed/reference/add_tbv.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# All individuals — all phenotypes
pop <- pop |> get_table("ind_meta") |> add_phenotype()

# Named phenotype, filtered subset
pop <- pop |>
  get_table("ind_meta") |>
  dplyr::filter(sex == "F", gen == 1L) |>
  add_phenotype("ADG")

# Composite (maternal) phenotype: WW = direct (self) + maternal (dam) TBV,
# registered once via define_phenotype(components = ...)
pop <- pop |>
  define_phenotype("WW", type = "continuous", mean = 230, residual_var = 180,
    components = tibble::tribble(
      ~source_trait_name, ~contributor_type,
      "WWD",              "self",
      "WWM",              "dam"
    ))
pop <- pop |>
  get_table("ind_meta") |>
  dplyr::filter(gen == 1L) |>
  add_phenotype("WW")

# Escape hatch: supply phenotype values directly (skips the model, but
# still computes and stores TBVs); named vector matches by id_ind
pop <- pop |>
  get_table("ind_meta") |>
  dplyr::filter(gen == 0L, sex == "M") |>
  add_phenotype("ADG", user_values = c(A_1 = 555, A_2 = 560))
} # }
```
