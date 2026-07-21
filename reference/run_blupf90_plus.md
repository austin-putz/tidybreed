# Run blupf90+ in the eval directory (reads renf90.par written by renumf90)

Run blupf90+ in the eval directory (reads renf90.par written by
renumf90)

## Usage

``` r
run_blupf90_plus(eval_dir, blupf90_path)
```

## Arguments

- eval_dir:

  path to evaluation folder containing renf90.par

- blupf90_path:

  path to the blupf90+ binary (from
  [`find_blupf90_binary()`](https://austin-putz.github.io/tidybreed/reference/find_blupf90_binary.md))

## Value

`NULL` invisibly; errors if blupf90+ exits non-zero
