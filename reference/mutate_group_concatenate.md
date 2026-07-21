# Create a composite group column by concatenating existing columns

Deterministically builds a `VARCHAR` group label for each filtered row
by concatenating the values of one or more existing columns in the same
table, separated by `sep`. Useful for contemporary groups
(`sex_line_gen`), sib groups (`sire_id_dam_id`), and any label derived
from existing row data.

Pipe `get_table("ind_meta")` — optionally through
[`dplyr::filter()`](https://dplyr.tidyverse.org/reference/filter.html) —
into this function. All `source_columns` must exist in the same table as
the target column. No randomness is involved; results are fully
deterministic.

## Usage

``` r
mutate_group_concatenate(
  tbl,
  col_name,
  source_columns,
  sep = "_",
  null_action = "propagate",
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

- source_columns:

  Character vector of one or more column names in the same table to
  concatenate. Order is preserved.

- sep:

  Character scalar. Separator string. Default `"_"`.

- null_action:

  One of `"propagate"` (default), `"skip"`, or `"literal"`. See section
  NULL handling above.

- overwrite:

  Logical. If `FALSE` (default), error when any filtered row already has
  a non-`NULL` value in `col_name`.

- quiet:

  Logical. If `TRUE`, suppress informational messages.

## Value

The modified `tidybreed_pop` object (invisibly).

## NULL handling

The `null_action` argument controls what happens when any source column
value is `NULL` (`NA`) for a given row:

- `"propagate"` (default): any `NULL` source produces a `NULL` result
  for that row. Safe for sib-group labels where a missing parent ID
  should produce no group, not a partial label.

- `"skip"`: `NULL` sources are silently omitted; remaining non-`NULL`
  values are concatenated. A row where all sources are `NULL` still
  yields a `NULL` result.

- `"literal"`: `NULL` sources are replaced with the string `"NA"` before
  concatenation.

## See also

[`mutate_group_seq()`](https://austin-putz.github.io/tidybreed/reference/mutate_group_seq.md),
[`mutate_group_named()`](https://austin-putz.github.io/tidybreed/reference/mutate_group_named.md),
[`mutate_table()`](https://austin-putz.github.io/tidybreed/reference/mutate_table.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# Contemporary group: sex × line_name
pop <- pop |>
  get_table("ind_meta") |>
  mutate_group_concatenate(
    col_name       = "cg",
    source_columns = c("sex", "line_name")
  )

# Sib group from parent IDs — NULL propagates when a parent is unknown
pop <- pop |>
  get_table("ind_meta") |>
  mutate_group_concatenate(
    col_name       = "sib_group",
    source_columns = c("id_parent_1", "id_parent_2"),
    null_action    = "propagate"
  )

# Generation × sex label, NULLs as literal "NA"
pop <- pop |>
  get_table("ind_meta") |>
  dplyr::filter(gen == 1L) |>
  mutate_group_concatenate(
    col_name       = "gen_sex",
    source_columns = c("gen", "sex"),
    sep            = ":",
    null_action    = "literal"
  )
} # }
```
