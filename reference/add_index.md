# Compute a selection index from a tidybreed table

Calculates a named selection index by multiplying each individual's
values by the index weights defined in
[`define_index()`](https://austin-putz.github.io/tidybreed/reference/define_index.md)
and summing the products. Results are appended to the `ind_index` table
with an auto-incrementing `index_number` so that successive runs are
preserved.

The first argument must be a `tidybreed_table` obtained from
[`get_table()`](https://austin-putz.github.io/tidybreed/reference/get_table.md)
(optionally filtered). The table must contain `id_ind`, a
trait/phenotype name column, and a numeric value column (`value_col`).
For `ind_phenotype` the name column is `phenotype_name`; for every other
table (`ind_ebv`, `ind_tbv`, or a user-defined table) it is `trait_name`
— this is detected automatically from the table's columns. Works with
`ind_ebv` (EBVs), `ind_phenotype` (phenotypes), `ind_tbv` (TBVs), or any
user-defined table with the same structure.

### Uniqueness requirement

There must be exactly one row per `(id_ind, trait_name)` after any
filter is applied. If duplicates remain, an error is thrown — filter the
table down to a single model, evaluation, or phenotype record before
calling `add_index()`.

### Completeness requirement

Every individual must have a value for **every** trait in the index. If
any individual is missing a value for any required trait, an error is
thrown before any rows are written.

## Usage

``` r
add_index(
  tbl,
  index_name,
  value_col = NULL,
  overwrite_index = FALSE,
  delete_all = FALSE,
  ...
)
```

## Arguments

- tbl:

  A `tidybreed_table` from
  [`get_table()`](https://austin-putz.github.io/tidybreed/reference/get_table.md)
  (optionally filtered). Must contain `id_ind`, a name column
  (`trait_name`, or `phenotype_name` for `ind_phenotype`), and the
  column specified by `value_col`. Any table with this structure is
  accepted: `ind_ebv`, `ind_phenotype`, `ind_tbv`, or a custom table.

- index_name:

  Character scalar. Name of the index to compute; must already exist in
  `index_meta` (created via
  [`define_index()`](https://austin-putz.github.io/tidybreed/reference/define_index.md)).

- value_col:

  Character scalar or `NULL`. The column in `tbl` that holds the numeric
  value to weight. When `NULL` (default), auto-detected from the table
  name: `ind_ebv` → `"ebv_value"`, `ind_phenotype` → `"pheno_value"`,
  `ind_tbv` → `"tbv_value"`. An error is thrown for unknown tables if
  `value_col` is not supplied.

- overwrite_index:

  Logical. If `TRUE`, all existing `ind_index` rows for this
  `index_name` are deleted before inserting; new rows receive
  `index_number = 1`. Default `FALSE`.

- delete_all:

  Logical. If `TRUE`, **all** rows in `ind_index` are deleted before
  inserting; new rows receive `index_number = 1`. Takes precedence over
  `overwrite_index`. Default `FALSE`.

- ...:

  Optional extra columns to add to `ind_index`. Scalar values only
  (broadcast to all inserted rows). Types are inferred via
  [`infer_duckdb_type()`](https://austin-putz.github.io/tidybreed/reference/infer_duckdb_type.md).

## Value

The `tidybreed_pop` (invisibly). Assign the result back.

## See also

[`define_index()`](https://austin-putz.github.io/tidybreed/reference/define_index.md),
[`get_table()`](https://austin-putz.github.io/tidybreed/reference/get_table.md),
[`add_ebv()`](https://austin-putz.github.io/tidybreed/reference/add_ebv.md),
[`add_phenotype()`](https://austin-putz.github.io/tidybreed/reference/add_phenotype.md),
[`add_tbv()`](https://austin-putz.github.io/tidybreed/reference/add_tbv.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# Compute index from EBVs — filter to a specific model and eval_number
pop <- pop |>
  get_table("ind_ebv") |>
  dplyr::filter(model == "blup_v1", eval_number == 1L) |>
  add_index("terminal")

# Compute phenotypic index — filter to first record per individual per trait
pop <- pop |>
  get_table("ind_phenotype") |>
  dplyr::filter(pheno_number == 1L) |>
  add_index("terminal")

# Compute true-value index from TBVs (value_col auto-detected)
pop <- pop |>
  get_table("ind_tbv") |>
  add_index("terminal")

# Explicit value_col for a user-defined table
pop <- pop |>
  get_table("my_scores") |>
  add_index("terminal", value_col = "composite_score")

# Replace previous run with a fresh one
pop <- pop |>
  get_table("ind_ebv") |>
  dplyr::filter(model == "blup_v2", eval_number == 1L) |>
  add_index("terminal", overwrite_index = TRUE)
} # }
```
