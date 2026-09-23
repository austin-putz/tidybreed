# Validate `user_residual` against the plan

A single model-path phenotype takes a numeric vector; otherwise a named
list keyed by `phenotype_name` that may name any subset of the
model-path phenotypes. Each vector is positional over that phenotype's
planned records.

## Usage

``` r
.ap_fixed_residuals(entries, user_residual)
```

## Arguments

- entries:

  The plan's entries.

- user_residual:

  The `user_residual` argument, or `NULL`.

## Value

A list named by phenotype of finite numeric vectors (possibly empty).
