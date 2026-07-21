# Insert a long haplotype frame into `ind_haplotype` (no UNPIVOT)

The single shared write core for
[`add_offspring()`](https://austin-putz.github.io/tidybreed/reference/add_offspring.md)
and
[`add_founders()`](https://austin-putz.github.io/tidybreed/reference/add_founders.md).
Registers `long_frame` and does one `INSERT ... SELECT` that joins
`genome_meta` on `locus_id` to attach `locus_name`. `strand` is always
`1L` in this ploidy-2 version. RNG-neutral (register + INSERT, no
`dbWriteTable`).

## Usage

``` r
.write_long_haplotypes(conn, long_frame)
```

## Arguments

- conn:

  DuckDB connection.

- long_frame:

  data.frame with columns `id_ind`, `parent_origin`, `strand`,
  `line_origin`, `locus_id`, `allele`.

## Value

Invisibly `NULL`.

## Details

Write mechanism (Stage-1 Appender vs register+INSERT decision): we use
`duckdb_register()` + `INSERT ... SELECT`, not the DuckDB Appender API.
It is RNG-neutral (unlike `dbWriteTable`, it advances no R RNG — the
founder/offspring RNG discipline relies on this), zero-copy over the R
frame, and already the fast path the pre-Stage-1 UNPIVOT write used. The
Appender API would save only the per-batch view registration (marginal,
and batching already bounds that), at the cost of a C++ boundary; not
worth it here. Revisit if profiling ever shows registration overhead
dominating.
