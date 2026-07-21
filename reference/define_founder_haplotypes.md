# Define founder haplotypes for a tidybreed population

Generates a pool of founder haplotypes and stores them in the
`founder_haplotypes` table. This table is required by
[`add_founders()`](https://austin-putz.github.io/tidybreed/reference/add_founders.md)
to assign phased alleles to new individuals.

Call this function after
[`define_genome()`](https://austin-putz.github.io/tidybreed/reference/define_genome.md)
and before
[`add_founders()`](https://austin-putz.github.io/tidybreed/reference/add_founders.md).
You may call it **multiple times** with different `line_name` values to
create line-specific haplotype pools (e.g., different LD structures per
line).
[`add_founders()`](https://austin-putz.github.io/tidybreed/reference/add_founders.md)
will then sample only from the pool matching its own `line_name`
argument.

Also writes a `founder_allele_freq` column to `genome_meta` recording
the per-locus allele frequency of the most recently generated pool (the
empirical column mean for LD methods). For multi-line setups needing
accurate per-line Falconer centering, use `base = "current_pop"` in
[`define_additive_effects()`](https://austin-putz.github.io/tidybreed/reference/define_additive_effects.md)
instead.

## Usage

``` r
define_founder_haplotypes(
  pop,
  n_haplotypes,
  method = "uniform",
  allele_freq = NULL,
  min_allele_freq = NULL,
  max_allele_freq = NULL,
  beta_shape1 = NULL,
  beta_shape2 = NULL,
  fst = NULL,
  mean_freq = NULL,
  n_templates = NULL,
  switch_rate = NULL,
  decay_rate = NULL,
  line_name = NULL
)
```

## Arguments

- pop:

  A `tidybreed_pop` object returned by
  [`define_genome()`](https://austin-putz.github.io/tidybreed/reference/define_genome.md).

- n_haplotypes:

  Positive integer. Number of founder haplotypes to generate.

- method:

  Character scalar. One of `"uniform"` (default), `"fixed"`, `"beta"`,
  `"balding_nichols"`, `"mosaic"`, or `"gaussian_copula"`. Passing an
  argument that belongs to a different method raises an error.

- allele_freq:

  *(method = "fixed" only)* Numeric scalar in \[0, 1\]. Allele frequency
  applied to every locus. Default `0.5`.

- min_allele_freq:

  *(method = "uniform" only)* Lower bound. Default `0.01`.

- max_allele_freq:

  *(method = "uniform" only)* Upper bound. Default `0.99`.

- beta_shape1:

  *(method = "beta" only)* First Beta shape parameter (α). Must be \> 0.
  Default `0.5`.

- beta_shape2:

  *(method = "beta" only)* Second Beta shape parameter (β). Must be
  \> 0. Default `0.5`.

- fst:

  *(method = "balding_nichols" only)* Wright's F_ST. Scalar in (0, 1).
  Controls spread of frequencies around `mean_freq`; larger values
  produce more extreme frequencies. Default `0.1`.

- mean_freq:

  *(method = "balding_nichols" only)* Ancestral mean allele frequency.
  Scalar in (0, 1). Default `0.5`.

- n_templates:

  *(method = "mosaic" only)* Number of template haplotypes. Must be
  \>= 2. Default `max(2, ceiling(sqrt(n_haplotypes)))`.

- switch_rate:

  *(method = "mosaic" only)* Expected template switches per cM (genetic
  distance). Higher values create shorter LD blocks. Default `1.0`. (At
  the default `cM_per_Mb = 1`, per-cM equals per-Mb.)

- decay_rate:

  *(method = "gaussian_copula" only)* LD decay rate λ in ρ = exp(−λ ×
  d_cM). Higher values → faster LD decay. Default `1.0` gives ρ ≈ 0.37
  at 1 cM.

- line_name:

  Optional character scalar. A label for this haplotype pool (e.g.
  `"LineA"`). Must start with a letter and contain only letters,
  numbers, or underscores. When supplied, `founder_haplotypes` rows are
  tagged with this label (`haplotype_id` is sequential within the pool).
  [`add_founders()`](https://austin-putz.github.io/tidybreed/reference/add_founders.md)
  with the same `line_name` will sample exclusively from these rows.
  When `NULL` (default), rows are stored with `line_name = NA` and
  [`add_founders()`](https://austin-putz.github.io/tidybreed/reference/add_founders.md)
  falls back to them when no named pool exists for the requested line.

## Value

The `tidybreed_pop` object (invisibly), with `founder_haplotypes`
registered in `pop$tables` and `founder_allele_freq` added to
`genome_meta`.

## Methods

- `"uniform"` (default):

  Per-locus allele frequency sampled from Uniform(`min_allele_freq`,
  `max_allele_freq`). No LD structure.

- `"fixed"`:

  Every locus gets the same `allele_freq` (default 0.5). No LD
  structure. Useful for quick sanity checks.

- `"beta"`:

  Per-locus allele frequency sampled from Beta(`beta_shape1`,
  `beta_shape2`). Default `shape1 = shape2 = 0.5` (Jeffreys prior) gives
  a U-shaped MAF distribution with many rare and common alleles —
  biologically realistic. No LD structure.

- `"balding_nichols"`:

  Balding-Nichols model: per-locus frequency from Beta(α, β) where α =
  `mean_freq` × (1 − `fst`) / `fst` and β = (1 − `mean_freq`) × (1 −
  `fst`) / `fst`. Models allele frequency drift around a mean. No LD
  structure.

- `"mosaic"`:

  Quick LD via haplotype block copying. Generates `n_templates` template
  haplotypes, then builds each new haplotype as a mosaic: copies from a
  template and switches templates with probability
  `1 - exp(-switch_rate × d_cM)` at each locus. Creates realistic LD
  blocks without external software. Uses the resolved default genetic
  map (`genome_map`) for adjacent-locus genetic distances. Performance
  is O(n_haplotypes × n_loci); use `"gaussian_copula"` for very large
  simulations.

- `"gaussian_copula"`:

  Fast LD via AR(1) latent normal. A latent Gaussian AR(1) process runs
  along each chromosome; adjacent-locus correlation decays as ρ =
  exp(−`decay_rate` × d_cM). Haplotypes are generated by thresholding
  latent values at each locus's target allele frequency. Fully
  vectorised over haplotypes. Uses the resolved default genetic map
  (`genome_map`).

## See also

[`define_genome()`](https://austin-putz.github.io/tidybreed/reference/define_genome.md),
[`add_founders()`](https://austin-putz.github.io/tidybreed/reference/add_founders.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# Each method demo below starts from a fresh pop — define_founder_haplotypes()
# errors if called again with line_name = NULL once a NULL-line pool already
# exists, so these are alternative choices, not a sequential recipe.
new_pop <- function() {
  open_pop(pop_name = "A", db_name = ":memory:") |>
    define_genome(n_loci = 1000, n_chr = 10, chr_len_Mb = 100)
}

# Uniform allele frequencies (default)
pop <- new_pop() |> define_founder_haplotypes(n_haplotypes = 200)

# Fixed allele frequency of 0.5
pop <- new_pop() |> define_founder_haplotypes(n_haplotypes = 200, method = "fixed")

# Beta(0.5, 0.5) — biologically realistic MAF distribution
pop <- new_pop() |> define_founder_haplotypes(n_haplotypes = 200, method = "beta")

# Balding-Nichols with FST = 0.05
pop <- new_pop() |> define_founder_haplotypes(n_haplotypes = 200,
                                 method = "balding_nichols", fst = 0.05)

# Quick LD via mosaic block copying
pop <- new_pop() |> define_founder_haplotypes(n_haplotypes = 200, method = "mosaic")

# Fast LD via Gaussian copula AR(1)
pop <- new_pop() |> define_founder_haplotypes(n_haplotypes = 200,
                                 method = "gaussian_copula", decay_rate = 0.5)

# Two lines with different LD structures (same pop, distinct line_name)
pop <- new_pop()
pop <- pop |> define_founder_haplotypes(n_haplotypes = 200, line_name = "A",
                                 method = "mosaic")
pop <- pop |> define_founder_haplotypes(n_haplotypes = 200, line_name = "B",
                                 method = "gaussian_copula", decay_rate = 0.5)
} # }
```
