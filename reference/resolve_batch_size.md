# Resolve the offspring/founders-per-batch size `B`

Precedence: explicit `batch_size` \> explicit `max_batch_mem` \>
RAM-aware auto-pick (target `target_frac` of detected available memory)
\> conservative fixed fallback (`fallback_bytes`) when detection fails.
Always clamped to `[1, n_total]`. Any `B` yields byte-identical output —
this only bounds peak memory, never changes values.

## Usage

``` r
resolve_batch_size(
  n_total,
  n_loci,
  batch_size = NULL,
  max_batch_mem = NULL,
  bytes_per_offspring_row = 200,
  target_frac = 0.25,
  fallback_bytes = 512 * 1024^2
)
```

## Arguments

- n_total:

  Integer. Total rows to process (offspring or founders).

- n_loci:

  Integer. Loci count (drives per-row memory).

- batch_size:

  Optional explicit batch size (wins when set).

- max_batch_mem:

  Optional explicit per-batch memory budget (bytes or `"512MB"`-style
  string).

- bytes_per_offspring_row:

  Heuristic peak bytes per (offspring x locus) during the in-R long
  build + register. Calibrated (~200) from Phase-0/Stage-1 measurements:
  each offspring emits `2 x n_loci` long rows (two haplotypes), each ~90
  bytes at peak including the `rbind`/register transient copy — so the
  per-(offspring x locus) constant is ~2 x 90 ~ 190, rounded up for
  headroom.

- target_frac:

  Fraction of available memory to target for one batch.

- fallback_bytes:

  Conservative per-batch budget when detection fails.

## Value

Integer `B` in `[1, n_total]` (or `0L` when `n_total == 0`).
