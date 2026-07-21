# Pass a chromosome copy through to a gamete without recombination

Used for non-recombining / single-copy special chromosomes (Y, W, MT,
most organelles). Draws from the **current, already-seeded** dqrng
stream (the caller seeds the gamete's `"special"` stream once, then
processes all of that gamete's special chromosomes in sequence). When
`k == 2`, one dqrng uniform chooses the homolog (`u < 0.5`); when
`k == 1` there is nothing to choose, so strictly hemizygous inheritance
consumes **no** RNG.

## Usage

``` r
pass_through_gamete(hap_matrix, lo_matrix)
```

## Arguments

- hap_matrix:

  `k x n_chr_loci` integer matrix, `k` in `{1, 2}`.

- lo_matrix:

  `k x n_chr_loci` character `line_origin` matrix, parallel.

## Value

List with `allele` and `line_origin`, each length `n_chr_loci`.
