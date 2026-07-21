# Trace pedigree n generations back from a set of candidate IDs

Trace pedigree n generations back from a set of candidate IDs

## Usage

``` r
trace_pedigree(pop, subset_ids, n_gen)
```

## Arguments

- pop:

  tidybreed_pop

- subset_ids:

  character vector of starting animal IDs

- n_gen:

  integer number of generations to trace back

## Value

data.frame with columns id_ind, id_parent_1, id_parent_2 (unknown
parents replaced with "0")
