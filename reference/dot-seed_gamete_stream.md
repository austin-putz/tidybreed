# Seed the global dqrng stream for one gamete (Scheme A)

Seed the global dqrng stream for one gamete (Scheme A)

## Usage

``` r
.seed_gamete_stream(base_seed, sid)
```

## Arguments

- base_seed:

  Integer base seed for the whole
  [`add_offspring()`](https://austin-putz.github.io/tidybreed/reference/add_offspring.md)
  call.

- sid:

  Integer stream id from
  [`.gamete_stream_id()`](https://austin-putz.github.io/tidybreed/reference/dot-gamete_stream_id.md).

## Value

Invisibly `NULL`. Sets the dqrng generator kind + seed as a side effect;
subsequent
[`dqrng::dqrunif()`](https://daqana.github.io/dqrng/reference/dqrng-functions.html)
draws consume this gamete's stream.
