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
empirical column mean for LD methods and for `exact_freq = TRUE`; the
sampled target otherwise). This column is **informational only** — no
other tidybreed function reads it, and because each call rewrites it for
every locus, in a multi-line setup it describes only the pool written
last. For multi-line setups needing accurate per-line Falconer
centering, use `base = "current_pop"` in
[`define_additive_effects()`](https://austin-putz.github.io/tidybreed/reference/define_additive_effects.md)
instead: with `base = "founder_haplotypes"` the base frequency is
recomputed by pooling **all** lines' haplotypes together.

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
  mean_allele_freq = NULL,
  n_templates = NULL,
  template_switch_rate = NULL,
  ld_decay_rate = NULL,
  exact_freq = NULL,
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
  realized *exactly* at every locus. Default `0.5`. Frequencies live on
  a `1 / n_haplotypes` grid; if `allele_freq × n_haplotypes` is not a
  whole number it is rounded, with a warning naming the frequency
  actually used.

- min_allele_freq:

  *(method = "uniform" only)* Lower bound. Default `0.01`.

- max_allele_freq:

  *(method = "uniform" only)* Upper bound. Default `0.99`.

- beta_shape1:

  *(method = "beta" only)* First Beta shape parameter (α). Must be
  finite and \> 0. Default `0.5`.

- beta_shape2:

  *(method = "beta" only)* Second Beta shape parameter (β). Must be
  finite and \> 0. Default `0.5`.

- fst:

  *(method = "balding_nichols" only)* Wright's F_ST. Scalar in (0, 1).
  Controls spread of frequencies around `mean_allele_freq`; larger
  values produce more extreme frequencies. Default `0.1`.

- mean_allele_freq:

  *(method = "balding_nichols" only)* Ancestral mean allele frequency.
  Scalar in (0, 1). Default `0.5`.

- n_templates:

  *(method = "mosaic" only)* Number of template haplotypes. Must be a
  whole number in `[2, n_haplotypes]`. Also controls the MAF spectrum —
  see the `"mosaic"` entry above. Default
  `max(2, ceiling(sqrt(n_haplotypes)))`.

- template_switch_rate:

  *(method = "mosaic" only)* Template re-draw rate per cM (genetic
  distance). Higher values create shorter LD blocks. `0` means never
  switch (complete LD within a chromosome). Note that observable
  template *changes* occur at
  `template_switch_rate × (n_templates − 1) / n_templates` — see the
  `"mosaic"` entry above. Default `1.0`. (At the default
  `cM_per_Mb = 1`, per-cM equals per-Mb.)

- ld_decay_rate:

  *(method = "gaussian_copula" only)* LD decay rate λ in ρ = exp(−λ ×
  d_cM). Higher values → faster LD decay. `0` means no decay (complete
  LD within a chromosome). Default `1.0` gives ρ ≈ 0.37 at 1 cM.

- exact_freq:

  *(methods "fixed", "uniform", "beta", "balding_nichols")* Logical.
  When `TRUE`, each locus receives exactly `round(p × n_haplotypes)`
  copies of the 1-allele on an independently drawn random subset of
  haplotypes, so the **realized** pool frequency equals the per-locus
  target `p` with no binomial fluctuation (no LD is induced — each locus
  draws its own subset). When `FALSE`, alleles are drawn as independent
  `Bernoulli(p)` trials, so realized frequencies scatter around `p` with
  sd `sqrt(p(1-p)/n_haplotypes)`.

  Defaults to `TRUE` for `method = "fixed"` (a "fixed" frequency that
  drifts is not fixed) and `FALSE` for the distribution-based methods,
  where drawing a frequency and then sampling alleles binomially is the
  correct generative model. Set `TRUE` on those methods when you want
  `founder_allele_freq` to be the exact base frequency for Falconer
  centering. Frequencies are snapped to the `1 / n_haplotypes` grid; for
  `method = "fixed"` a warning names the frequency actually used when
  the requested one is off-grid (for the distribution methods a
  continuous draw is essentially never on the grid, so rounding is
  silent).

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

  Every locus gets the same `allele_freq` (default 0.5), **exactly**.
  Rather than drawing each allele as a Bernoulli(`allele_freq`) trial —
  which leaves the realized pool frequency fluctuating binomially around
  the target — exactly `round(allele_freq × n_haplotypes)` of the
  haplotypes carry the 1-allele at each locus, on an independently drawn
  random subset per locus (so no LD is induced). No LD structure.

- `"beta"`:

  Per-locus allele frequency sampled from Beta(`beta_shape1`,
  `beta_shape2`). Default `shape1 = shape2 = 0.5` (Jeffreys prior) gives
  a U-shaped MAF distribution with many rare and common alleles —
  biologically realistic. No LD structure.

- `"balding_nichols"`:

  Balding-Nichols model: per-locus frequency from Beta(α, β) where α =
  `mean_allele_freq` × (1 − `fst`) / `fst` and β = (1 −
  `mean_allele_freq`) × (1 − `fst`) / `fst`. Models allele frequency
  drift around a mean. No LD structure.

- `"mosaic"`:

  Quick LD via haplotype block copying. Generates `n_templates` template
  haplotypes, then builds each new haplotype as a mosaic: copies from a
  template and re-draws a template with probability
  `1 - exp(-template_switch_rate × d_cM)` at each locus. Creates
  realistic LD blocks without external software. Uses the genetic map
  (`genome_map`) resolved for this pool's `line_name` for adjacent-locus
  distances. Performance is O(n_haplotypes × n_loci); use
  `"gaussian_copula"` for very large simulations.

  Two properties of this model are easy to miss. First, the re-draw is
  uniform over **all** templates including the current one, so the rate
  of *observable* template changes is
  `template_switch_rate × (n_templates − 1) / n_templates`. This is
  deliberate — it is the standard Li-Stephens copying kernel, and it is
  what makes realized LD invariant to marker density. Second, templates
  are the only source of allelic variation, so a locus where all
  `n_templates` templates agree is monomorphic **regardless of
  `n_haplotypes`**; roughly `2 / (n_templates + 1)` of loci are
  monomorphic, and the MAF spectrum is quantized to multiples of
  `1 / n_templates`. Raise `n_templates` if you need rare variants or a
  dense MAF spectrum.

- `"gaussian_copula"`:

  Fast LD via AR(1) latent normal. A latent Gaussian AR(1) process runs
  along each chromosome; adjacent-locus correlation decays as ρ =
  exp(−`ld_decay_rate` × d_cM). Haplotypes are generated by thresholding
  latent values at each locus's target allele frequency. Fully
  vectorised over haplotypes, and unlike `"mosaic"` it produces an
  unquantized MAF spectrum. Uses the genetic map (`genome_map`) resolved
  for this pool's `line_name`.

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

# Beta frequencies, realized exactly (no binomial scatter around the draw)
pop <- new_pop() |> define_founder_haplotypes(n_haplotypes = 200,
                                 method = "beta", exact_freq = TRUE)

# Quick LD via mosaic block copying
pop <- new_pop() |> define_founder_haplotypes(n_haplotypes = 200, method = "mosaic")

# Fast LD via Gaussian copula AR(1)
pop <- new_pop() |> define_founder_haplotypes(n_haplotypes = 200,
                                 method = "gaussian_copula", ld_decay_rate = 0.5)

# Two lines with different LD structures (same pop, distinct line_name)
pop <- new_pop()
pop <- pop |> define_founder_haplotypes(n_haplotypes = 200, line_name = "A",
                                 method = "mosaic")
pop <- pop |> define_founder_haplotypes(n_haplotypes = 200, line_name = "B",
                                 method = "gaussian_copula", ld_decay_rate = 0.5)
} # }
```
