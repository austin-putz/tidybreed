# Assign sequential integer group IDs to filtered rows

Randomly shuffles the filtered rows and assigns them to consecutive
integer group labels (pen 1, 2, 3 …). The column is always `INTEGER`.
Useful for pen assignment, contemporary group numbering, and any
situation where each group is identified by a unique integer.

Pipe `get_table("ind_meta")` — optionally through
[`dplyr::filter()`](https://dplyr.tidyverse.org/reference/filter.html) —
into this function. Filtered primary keys are sorted ascending before
shuffling, giving a stable, reproducible starting point.

## Usage

``` r
mutate_group_seq(
  tbl,
  col_name,
  n_per_group = NULL,
  n_groups = NULL,
  start = NULL,
  include_leftover = FALSE,
  seed = NULL,
  overwrite = FALSE,
  quiet = FALSE
)
```

## Arguments

- tbl:

  A `tidybreed_table` from
  [`get_table()`](https://austin-putz.github.io/tidybreed/reference/get_table.md)
  (optionally filtered).

- col_name:

  Character scalar. Name of the group column to create or update. Must
  be a valid SQL identifier and not a reserved column.

- n_per_group:

  Integer scalar or length-2 integer vector `c(min, max)`. Mutually
  exclusive with `n_groups`.

- n_groups:

  Positive integer scalar. Split filtered rows into this many groups.
  Mutually exclusive with `n_per_group`.

- start:

  Integer scalar. First group label. Default: `MAX(col_name) + 1` if the
  column exists, else `1`.

- include_leftover:

  Logical. If `TRUE`, remainder animals are assigned to a final
  (smaller) group. If `FALSE` (default), they are left `NULL`.

- seed:

  Integer or `NULL`. If supplied, sets the RNG seed before shuffling so
  results are reproducible. Validation happens before the seed is set,
  so a failed call does not consume RNG state.

- overwrite:

  Logical. If `FALSE` (default), error when any filtered row already has
  a non-`NULL` value in `col_name`.

- quiet:

  Logical. If `TRUE`, suppress informational messages.

## Value

The modified `tidybreed_pop` object (invisibly).

## Group-size modes

Exactly one of `n_per_group` or `n_groups` must be supplied.

- **`n_per_group` scalar**: every group has exactly `n_per_group`
  animals. Remainder animals are left `NULL` unless
  `include_leftover = TRUE`.

- **`n_per_group` length-2 vector `c(min, max)`**: each group is sized
  by drawing uniformly from `min:max`. When remaining animals fall below
  `min`, they are left `NULL` (or included if
  `include_leftover = TRUE`).

- **`n_groups`**: splits the filtered set into `n_groups` groups of size
  `floor(n / n_groups)`. When `include_leftover = TRUE`, the first
  `n %% n_groups` groups each receive one extra animal (balanced split).
  When `include_leftover = FALSE`, the remainder animals are left
  `NULL`.

## Start numbering

When `start = NULL` (default), the next label is `MAX(col_name) + 1` if
the column already exists, or `1` for a new column. Supplying
`start = 1L` resets numbering — useful when assigning pens independently
per generation. Setting `start` does **not** require `overwrite = TRUE`;
overwrite only applies when **filtered rows** already hold non-`NULL`
values.

## See also

[`mutate_group_named()`](https://austin-putz.github.io/tidybreed/reference/mutate_group_named.md),
[`mutate_group_concatenate()`](https://austin-putz.github.io/tidybreed/reference/mutate_group_concatenate.md),
[`mutate_table()`](https://austin-putz.github.io/tidybreed/reference/mutate_table.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# Assign males to pens of 10 (remainder left NULL)
pop <- pop |>
  get_table("ind_meta") |>
  dplyr::filter(sex == "M") |>
  mutate_group_seq(col_name = "pen", n_per_group = 10L)

# Assign all individuals to 5 balanced groups
pop <- pop |>
  get_table("ind_meta") |>
  mutate_group_seq(col_name = "pen", n_groups = 5L, include_leftover = TRUE)

# Restart pen numbering from 1 for generation 2 (no overwrite = TRUE needed)
pop <- pop |>
  get_table("ind_meta") |>
  dplyr::filter(gen == 2L) |>
  mutate_group_seq(col_name = "pen", n_per_group = 10L, start = 1L)

# Reproducible assignment
pop <- pop |>
  get_table("ind_meta") |>
  mutate_group_seq(col_name = "pen", n_per_group = 10L, seed = 42L)
} # }
```
