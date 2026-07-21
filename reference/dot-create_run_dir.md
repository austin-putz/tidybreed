# Create a timestamped layer-5 run directory inside a tool dir

Called at the start of every external tool wrapper (e.g.
`run_blupf90()`, `run_plink()`). Creates a unique subdirectory inside
`pop$run_dirs[[tool]]` for the current invocation and returns the path.

The tool must have been declared via the `tools` argument of
[`open_pop()`](https://austin-putz.github.io/tidybreed/reference/open_pop.md).
If `tool` is not in `pop$run_dirs`, the function errors with a clear
message listing the registered tools.

For v0.40.0 the `keep` argument is accepted for API completeness but
cleanup is not yet automatically registered. The caller is responsible
for deleting the run directory on success if `keep = FALSE` or
`keep = "on_error"`.

## Usage

``` r
.create_run_dir(pop, tool, keep = "on_error")
```

## Arguments

- pop:

  A `tidybreed_pop` object with `run_dirs` populated.

- tool:

  Character scalar. Which tool dir to use (must be a key in
  `pop$run_dirs`).

- keep:

  `FALSE`, `TRUE`, or `"on_error"`. Cleanup semantics (not yet enforced
  automatically; reserved for a future release).

## Value

Path to the newly created run directory (invisibly).
