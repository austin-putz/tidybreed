# Insert a batch of crossover events into `ind_crossover`

Assigns `id_crossover` R-side via the BIGINT-safe
[`next_row_id()`](https://austin-putz.github.io/tidybreed/reference/next_row_id.md)
(the table can accumulate far more rows than any INTEGER-keyed table),
then registers the frame and does one `INSERT`. RNG-neutral (register +
INSERT). Called only when `add_offspring(store_crossovers = TRUE)`;
within a transaction,
[`next_row_id()`](https://austin-putz.github.io/tidybreed/reference/next_row_id.md)
sees prior batches' rows, so ids stay sequential across batches.

## Usage

``` r
.write_crossovers(conn, xframe)
```

## Arguments

- conn:

  DuckDB connection.

- xframe:

  data.frame with columns `id_ind`, `parent_origin`, `chr`, `chr_name`,
  `pos_cM` (one row per crossover event).

## Value

Invisibly `NULL`.
