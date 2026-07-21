# Assign named group labels to filtered rows by count or proportion

Randomly assigns each filtered row to one of the supplied `group_names`.
Group membership is determined either by explicit `counts` (how many
animals per group) or by `proportions` (relative fractions). The column
is always `VARCHAR`. Useful for farm assignment, treatment groups, line
splitting, and any situation where groups have user-defined string
labels.

Pipe `get_table("ind_meta")` — optionally through
[`dplyr::filter()`](https://dplyr.tidyverse.org/reference/filter.html) —
into this function. Filtered primary keys are sorted ascending before
shuffling, giving a stable, reproducible starting point.

## Usage

``` r
mutate_group_named(
  tbl,
  col_name,
  group_names,
  counts = NULL,
  proportions = NULL,
  prop_method = "exact",
  include_leftover = FALSE,
  underfill_action = "error",
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

- group_names:

  Character vector. Group labels. Must be the same length as `counts` or
  `proportions`.

- counts:

  Integer vector of animals per group. Mutually exclusive with
  `proportions`.

- proportions:

  Numeric vector of relative group fractions. Mutually exclusive with
  `counts`.

- prop_method:

  `"exact"` (default) uses the largest-remainder method for
  deterministic integer rounding; `"random"` draws each label
  independently.

- include_leftover:

  Logical. When `sum(counts) < n_filtered` (or when proportions sum to
  \< n due to rounding), add the unassigned rows to the last group.
  Default `FALSE` leaves them `NULL`.

- underfill_action:

  `"error"` (default), `"proportional"`, or `"sequential"`. Only applies
  when `counts` is supplied and `sum(counts) > n_filtered`.

- seed:

  Integer or `NULL`. If supplied, sets the RNG seed before shuffling.
  Validation happens before the seed is set.

- overwrite:

  Logical. If `FALSE` (default), error when any filtered row already has
  a non-`NULL` value in `col_name`.

- quiet:

  Logical. If `TRUE`, suppress informational messages.

## Value

The modified `tidybreed_pop` object (invisibly).

## Counts vs proportions

Exactly one of `counts` or `proportions` must be supplied.

- **`counts`**: integer vector, same length as `group_names`. When
  `sum(counts) > n_filtered`, `underfill_action` decides what happens.
  When `sum(counts) < n_filtered`, remaining rows are left `NULL` unless
  `include_leftover = TRUE` (which adds them to the last group).

- **`proportions`**: numeric vector, same length as `group_names`.
  Values must be positive and are normalized to sum to 1 within
  `sqrt(.Machine$double.eps)` tolerance. With `prop_method = "exact"`
  (default), the largest-remainder method converts proportions to
  integer counts, tie-breaking left-to-right. With
  `prop_method = "random"`, each animal is independently drawn from the
  distribution.

## Underfill behavior (counts mode only)

When `sum(counts) > n_filtered` and `underfill_action`:

- `"error"` (default): stop before sampling.

- `"proportional"`: rescale counts to fit `n_filtered` using the
  largest-remainder method.

- `"sequential"`: fill groups in order until animals run out; later
  groups may be partially or fully empty.

## See also

[`mutate_group_seq()`](https://austin-putz.github.io/tidybreed/reference/mutate_group_seq.md),
[`mutate_group_concatenate()`](https://austin-putz.github.io/tidybreed/reference/mutate_group_concatenate.md),
[`mutate_table()`](https://austin-putz.github.io/tidybreed/reference/mutate_table.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# Assign females to two farms by count
pop <- pop |>
  get_table("ind_meta") |>
  dplyr::filter(sex == "F") |>
  mutate_group_named(
    col_name    = "farm",
    group_names = c("farm_A", "farm_B"),
    counts      = c(60L, 40L)
  )

# Assign all individuals by proportion (60 / 40 split)
pop <- pop |>
  get_table("ind_meta") |>
  mutate_group_named(
    col_name     = "farm",
    group_names  = c("farm_A", "farm_B"),
    proportions  = c(0.6, 0.4)
  )

# If available animals are fewer than requested counts, rescale proportionally
pop <- pop |>
  get_table("ind_meta") |>
  dplyr::filter(gen == 1L) |>
  mutate_group_named(
    col_name         = "trt",
    group_names      = c("control", "high", "low"),
    counts           = c(50L, 50L, 50L),
    underfill_action = "proportional"
  )
} # }
```
