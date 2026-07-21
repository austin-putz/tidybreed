# Run renumf90 in the eval directory

Run renumf90 in the eval directory

## Usage

``` r
run_renumf90(eval_dir, renumf90_path)
```

## Arguments

- eval_dir:

  path to evaluation folder containing renum.par

- renumf90_path:

  path to the renumf90 binary (from
  [`find_blupf90_binary()`](https://austin-putz.github.io/tidybreed/reference/find_blupf90_binary.md))

## Value

`NULL` invisibly; errors if renumf90 exits non-zero
