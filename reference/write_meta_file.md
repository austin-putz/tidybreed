# Write the human-readable meta.txt for the evaluation folder

Write the human-readable meta.txt for the evaluation folder

## Usage

``` r
write_meta_file(
  eval_dir,
  eval_id,
  col_map,
  distinct_effects,
  trait,
  effects_df,
  chip_name,
  n_loci,
  id_width,
  animal_effect_num
)
```

## Arguments

- eval_dir:

  path to evaluation folder; writes meta.txt there

- eval_id:

  character scalar; evaluation label used in the file header

- col_map:

  named integer vector mapping mu/id_ind/effect_name/trait_name to their
  data.txt column numbers (as returned by
  [`build_data_file()`](https://austin-putz.github.io/tidybreed/reference/build_data_file.md))

- distinct_effects:

  data.frame of fixed effects (one row per effect_name), as returned by
  [`build_data_file()`](https://austin-putz.github.io/tidybreed/reference/build_data_file.md)

- trait:

  character vector of trait names (in model order)

- effects_df:

  data.frame of (trait_name x effect_name) fixed-effect rows from
  `trait_effects`, as returned by
  [`build_data_file()`](https://austin-putz.github.io/tidybreed/reference/build_data_file.md)

- chip_name:

  character or NULL; chip name to report in the genotype file section
  (omitted entirely when NULL)

- n_loci:

  integer; number of loci written to the genotype file

- id_width:

  integer; fixed field width used for animal IDs in the genotype file

- animal_effect_num:

  integer; effect number for the animal random effect, as returned by
  [`write_renum_par()`](https://austin-putz.github.io/tidybreed/reference/write_renum_par.md)

## Value

`NULL` invisibly; writes meta.txt as a side effect
