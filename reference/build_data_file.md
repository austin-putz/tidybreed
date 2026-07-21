# Build the BLUPF90 data file

Column layout: mu(1), id_ind(2), fixed_class effects, fixed_cov effects,
trait observations (missing → 0).

## Usage

``` r
build_data_file(pop, subset_ids, trait_name, eval_dir, pheno_ids = NULL)
```

## Arguments

- pop:

  tidybreed_pop

- subset_ids:

  character vector of animal IDs whose phenotypes to include

- trait_name:

  character vector of trait names

- eval_dir:

  path to evaluation folder; writes data.txt there

- pheno_ids:

  optional integer vector of `ind_phenotype.id_phenotype` values to
  restrict which phenotype records are included (e.g. to exclude future
  records). `NULL` (default) includes all matching records.

## Value

list with data, col_map, distinct_effects, effects_df, n_fixed_effects,
trait_cols
