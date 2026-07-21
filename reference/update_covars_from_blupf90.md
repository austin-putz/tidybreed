# Stub for VCE writeback — parse blupf90.out and update trait_var_comp

Not yet implemented; called by
[`add_ebv()`](https://austin-putz.github.io/tidybreed/reference/add_ebv.md)
when `estimate_var = TRUE` and `update_covars = TRUE`.

## Usage

``` r
update_covars_from_blupf90(pop, eval_dir, trait_name)
```

## Arguments

- pop:

  tidybreed_pop

- eval_dir:

  path to evaluation folder containing blupf90.out

- trait_name:

  character vector of trait names

## Value

`NULL` invisibly
