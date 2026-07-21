# Define a selection index

Registers a named selection index in the `index_meta` table by
specifying which traits are included and their weighting coefficients.
Call this once per index (before calling
[`add_index()`](https://austin-putz.github.io/tidybreed/reference/add_index.md)).
Re-calling with the same `index_name` and `trait_name` pairs updates the
weights in place.

Not all traits in the population need to appear in an index — you may
have five traits but index only three. Traits not listed in
`trait_names` are ignored when
[`add_index()`](https://austin-putz.github.io/tidybreed/reference/add_index.md)
computes the index value.

## Usage

``` r
define_index(
  pop,
  index_name,
  trait_names,
  index_wts,
  economic_wts = NULL,
  overwrite = FALSE,
  ...
)
```

## Arguments

- pop:

  A `tidybreed_pop` object.

- index_name:

  Character scalar. Name of the index (e.g. `"terminal"`, `"maternal"`).
  Must be a valid SQL identifier.

- trait_names:

  Character vector. Trait names to include in the index. Each must be a
  valid SQL identifier. Order must match `index_wts`.

- index_wts:

  Numeric vector. Selection index weights, one per trait in
  `trait_names`. Positive weights favour higher trait values; negative
  weights favour lower values.

- economic_wts:

  Numeric vector (optional). Economic value per unit for each trait, in
  the same order as `trait_names`. Some values may be `0`. When `NULL`
  (default), the `economic_weight` column is not written.

- overwrite:

  Logical. If `FALSE` (default), re-calling for an existing
  `(index_name, trait_name)` pair is a no-op — existing weights and
  economic values are preserved. If `TRUE`, all non-key columns
  (`index_weight`, `economic_weight`, and any `...` columns) are updated
  in place.

- ...:

  Optional extra columns to add to `index_meta`. Scalars are broadcast
  to all rows (one per trait). Vectors must have length equal to
  `length(trait_names)`. Types are inferred via
  [`infer_duckdb_type()`](https://austin-putz.github.io/tidybreed/reference/infer_duckdb_type.md).

## Value

The `tidybreed_pop` (invisibly). Assign the result back.

## See also

[`add_index()`](https://austin-putz.github.io/tidybreed/reference/add_index.md),
[`define_trait()`](https://austin-putz.github.io/tidybreed/reference/define_trait.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# Index weights alone
pop <- define_index(pop, "terminal",
                    trait_names = c("ADG", "FCR"),
                    index_wts   = c(1.2, -0.8))

# Also record the underlying economic weights (per-unit $ value), stored
# alongside index_weight in index_meta$economic_weight — e.g. for later use
# with add_tbv(index_names = ..., type = "economic")
pop <- define_index(pop, "maternal",
                    trait_names  = c("ADG", "FCR"),
                    index_wts    = c(0.5, -0.3),
                    economic_wts = c(2.5, -1.0))

# With an extra metadata column
pop <- define_index(pop, "growth_only",
                    trait_names = c("ADG"),
                    index_wts   = c(1.0),
                    source = "genetic_team_v2")

# Re-defining an index: overwrite = FALSE (default) is a no-op if the
# (index_name, trait_name) pair already exists — the call below does NOT
# change "terminal"'s weights since overwrite defaults to FALSE
pop <- define_index(pop, "terminal",
                    trait_names = c("ADG", "FCR"),
                    index_wts   = c(1.5, -0.5))

# overwrite = TRUE updates the weights in place
pop <- define_index(pop, "terminal",
                    trait_names = c("ADG", "FCR"),
                    index_wts   = c(1.5, -0.5),
                    overwrite   = TRUE)
} # }
```
