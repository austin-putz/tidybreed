# Reconnect to an existing tidybreed database

Re-opens a `.duckdb` file and reconstructs the `tidybreed_pop` R object
so that simulation can continue. No metadata is stored on the object;
all genome statistics are derived from live DB queries when needed.

This is the intended way to resume a simulation after
[`close_pop()`](https://austin-putz.github.io/tidybreed/reference/close_pop.md),
after a crash, or in a new R session. A common pattern is to call
[`open_pop()`](https://austin-putz.github.io/tidybreed/reference/open_pop.md)
once, configure all traits and covariances, then start each replicate
with `restore_pop()` followed by
[`add_founders()`](https://austin-putz.github.io/tidybreed/reference/add_founders.md).

`run_dirs` are reconstructed from the database location and the `tools`
argument. Set `options(tidybreed.tools = ...)` before calling
`restore_pop()` to get the same tool dirs back automatically.

## Usage

``` r
restore_pop(
  db_path,
  pop_name = NULL,
  tools = getOption("tidybreed.tools", NULL)
)
```

## Arguments

- db_path:

  Character. Path to an existing `.duckdb` file.

- pop_name:

  Character or `NULL`. Population name to assign to the restored object.
  When `NULL` (default) the name is inferred from the filename by
  stripping a trailing `_tidybreed.duckdb` or `.duckdb` suffix.

- tools:

  Character vector or `NULL`. Tool subdirectory names to reconstruct in
  `pop$run_dirs`. Defaults to `getOption("tidybreed.tools", NULL)`.

## Value

A fully operational `tidybreed_pop` object.

## Examples

``` r
if (FALSE) { # \dontrun{
# One-time setup. db_name is placed inside open_pop()'s layer-2/3 folder
# structure (see ?open_pop), so capture pop$db_path rather than assuming
# the file lives at "cattle_sim.duckdb" in the working directory.
pop <- open_pop(pop_name = "cattle", db_name = "cattle_sim.duckdb") |>
  define_genome(n_loci = 50000, n_chr = 29, chr_len_Mb = 100) |>
  define_founder_haplotypes(n_haplotypes = 100, line_name = "A")
pop <- define_trait(pop, trait_name = "milk", target_add_var = 10)
db_path <- pop$db_path
close_pop(pop)

# Each replicate
pop <- restore_pop(db_path)
pop <- pop |>
  get_table("founder_haplotypes") |>
  add_founders(n_males = 50, n_females = 200, line_name = "A")
} # }
```
