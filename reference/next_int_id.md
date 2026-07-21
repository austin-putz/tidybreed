# Get the next integer ID for a table's primary-key sequence.

Returns MAX(id_col) + 1, or 1 if the table is empty. For multi-row
inserts, pre-generate a range:
`seq.int(next_int_id(conn, tbl, col), length.out = n)`

## Usage

``` r
next_int_id(conn, table, id_col)
```

## Arguments

- conn:

  A DBI connection.

- table:

  Character. Table name.

- id_col:

  Character. Integer PK column name.

## Value

Integer scalar.
