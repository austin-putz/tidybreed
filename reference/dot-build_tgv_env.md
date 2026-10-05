# Pre-fetch every genetic-value vector a `formula_tgv` expression needs

One contributor lookup per reference (see
[`?contributor_tgv`](https://austin-putz.github.io/tidybreed/reference/contributor_tgv.md)),
reading the reference's `component` (`"total"` by default), returned as
a named list ready to be the
[`eval()`](https://rdrr.io/r/base/eval.html) environment.

## Usage

``` r
.build_tgv_env(pop, trait_refs, subset_df, phenotype_name)
```

## Arguments

- trait_refs:

  List from `.walk_formula_tgv_ast()$trait_refs`.

- subset_df:

  The planned `ind_meta` rows.

## Value

Named list: placeholder -\> numeric vector (`NA` = missing piece).
