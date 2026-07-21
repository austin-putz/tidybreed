# Fold (offspring, parent_origin, stream_kind) into a single dqrng stream id

Fold (offspring, parent_origin, stream_kind) into a single dqrng stream
id

## Usage

``` r
.gamete_stream_id(o, parent_origin, kind)
```

## Arguments

- o:

  Integer offspring index (the original `matings` row index; global,
  never batch-local — so any `batch_size` yields identical output).

- parent_origin:

  Integer 1 (parent_1/sire) or 2 (parent_2/dam).

- kind:

  Integer 1 (autosome) or 2 (special chromosome).

## Value

Integer stream id `sid`. Fits signed 32-bit for `o < ~5e8`.
