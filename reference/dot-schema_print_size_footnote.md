# The Size-column footnote, required whenever a Size column prints

The column is only defensible because the caveat travels with it, so
this has no suppression argument. It must state *both* caveats: the
block quantization (why small tables are indistinguishable) and the
shortfall against the header total (why the column does not add up). A
reader shown only one will read the other as a bug.

## Usage

``` r
.schema_print_size_footnote(width, db_size)
```

## Arguments

- width:

  Integer console width.

- db_size:

  Character. The header's whole-database label, referenced so the
  shortfall is self-evidently expected rather than looking like loss.

## Value

`NULL`, invisibly; called for the side effect of printing.
