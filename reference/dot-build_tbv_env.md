# Pre-fetch every TBV vector a `formula_tbv` expression needs

One contributor lookup per trait reference (see
[`?contributor_tbv`](https://austin-putz.github.io/tidybreed/reference/contributor_tbv.md)),
returned as a named list ready to be the
[`eval()`](https://rdrr.io/r/base/eval.html) environment.

## Usage

``` r
.build_tbv_env(pop, trait_refs, subset_df, phenotype_name)
```

## Arguments

- trait_refs:

  List from `.walk_formula_tbv_ast()$trait_refs`.

- subset_df:

  The planned `ind_meta` rows.

## Value

Named list: placeholder -\> numeric vector (`NA` = missing piece).
