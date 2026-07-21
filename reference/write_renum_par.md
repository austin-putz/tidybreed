# Write the renumf90 parameter file (renum.par)

Write the renumf90 parameter file (renum.par)

## Usage

``` r
write_renum_par(
  eval_dir,
  col_map,
  distinct_effects,
  effects_df,
  trait,
  pop,
  chip_name,
  estimate_var,
  alpha_size
)
```

## Arguments

- eval_dir:

  path to evaluation folder; writes renum.par there

- col_map:

  named integer vector mapping mu/id_ind/effect_name/trait_name to their
  data.txt column numbers (as returned by
  [`build_data_file()`](https://austin-putz.github.io/tidybreed/reference/build_data_file.md))

- distinct_effects:

  data.frame of fixed effects (one row per effect_name), as returned by
  [`build_data_file()`](https://austin-putz.github.io/tidybreed/reference/build_data_file.md)

- effects_df:

  data.frame of (trait_name x effect_name) fixed-effect rows from
  `trait_effects`, as returned by
  [`build_data_file()`](https://austin-putz.github.io/tidybreed/reference/build_data_file.md)

- trait:

  character vector of trait names (in model order)

- pop:

  tidybreed_pop; used to look up residual and additive genetic
  (co)variance matrices via
  [`load_phenotype_cov()`](https://austin-putz.github.io/tidybreed/reference/load_phenotype_cov.md)
  and
  [`load_trait_cov()`](https://austin-putz.github.io/tidybreed/reference/load_trait_cov.md)

- chip_name:

  character or NULL; when non-NULL, adds a `SNP_FILE` line for
  single-step GBLUP

- estimate_var:

  logical; if TRUE, sets `OPTION method VCE` instead of BLUP-only

- alpha_size:

  numeric; minimum size hint for `OPTION alpha_size` (rounded up to a
  multiple of 5, floor 20)

## Value

integer animal_effect_num (effect number to extract EBVs from solutions)
