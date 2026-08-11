# Add offspring via recombination

Produces offspring by simulating chromosomal crossovers from explicit
parent pairs. Each row of `matings` defines exactly one offspring,
giving full control over the mating design.

## Usage

``` r
add_offspring(
  pop,
  matings,
  seed = NULL,
  store_crossovers = FALSE,
  batch_size = NULL,
  max_batch_mem = NULL
)
```

## Arguments

- pop:

  A `tidybreed_pop` object

- matings:

  A tibble or data.frame where each row = one offspring.

  **Required columns** (canonical names match `ind_meta`):

  - `id_parent_1` — sire / parent 1 ID; must exist in `ind_meta`. Alias:
    `id_sire`

  - `id_parent_2` — dam / parent 2 ID; must exist in `ind_meta`. Alias:
    `id_dam`

  - `sex` — `"M"` or `"F"`

  - `line_name` — line name for offspring IDs (same format as
    [`add_founders()`](https://austin-putz.github.io/tidybreed/reference/add_founders.md))

  **Optional extra columns** (e.g. `gen = 2L`, `farm = "Iowa"`) are
  validated and written directly to `ind_meta`. If a column does not yet
  exist in `ind_meta` it is added automatically. Scalar values in a
  tibble are recycled to all rows before being passed to this function.

- seed:

  Optional integer base seed for the per-gamete recombination streams.
  Each gamete draws from its own deterministic `dqrng` sub-stream keyed
  on its global offspring index, so output is independent of
  `batch_size` and (later) thread count. `seed = NULL` (default) draws
  one base seed from base-R's RNG, so an upstream
  [`set.seed()`](https://rdrr.io/r/base/Random.html) still fully
  reproduces the call; pass an integer to fix the gamete streams
  independently of surrounding base-R state. The resolved base seed is
  reported in the message and attached to the returned object as
  `attr(pop, "base_seed")`; it is not written to any table.

- store_crossovers:

  Logical (default `FALSE`). When `TRUE`, every crossover drawn during
  meiosis is recorded as one row in `ind_crossover` (`id_ind`,
  `parent_origin`, `chr`, `chr_name`, `pos_cM`). Absence of a row for a
  `(id_ind, parent_origin, chr)` means that gamete's chromosome did not
  recombine. Includes crossovers on recombining special chromosomes.
  `FALSE` costs nothing (no buffer emitted, no `ind_crossover` write).

- batch_size:

  Optional integer. Offspring generated + written per batch. Bounds peak
  memory at roughly `batch_size x n_loci` long rows regardless of the
  total number of offspring (a single huge mating streams to disk batch
  by batch). Any value — including `batch_size = nrow(matings)` —
  produces byte-identical output for a fixed seed; it only trades peak
  memory against per-batch overhead. Overrides `max_batch_mem` when both
  are set.

- max_batch_mem:

  Optional per-batch memory budget as bytes (e.g. `512e6`) or a string
  (e.g. `"512MB"`, `"2GB"`). Used to derive `batch_size` when
  `batch_size` is `NULL`. When both are `NULL` (the default), the batch
  size is auto-picked from detected available system memory (falling
  back to a conservative fixed budget if detection is unavailable).

## Value

The modified `tidybreed_pop` object (invisibly). Assign the result back:
`pop <- add_offspring(pop, matings)`

## Details

**Mating design flexibility:**

- Multiple offspring per pair: repeat the pair row in `matings`

- Multiple sires per dam (pooled semen / polyspermy): include each sire
  as a separate row with the same `id_parent_2`

- Multiple dams per sire: include each dam as a separate row with the
  same `id_parent_1`

- Cross-line matings: use any `line_name` value in the `line_name`
  column

**Column aliases:** `id_sire` is accepted as an alias for `id_parent_1`,
and `id_dam` for `id_parent_2`. Both naming styles produce identical
results.

**Recombination model:** Gametes are simulated using the Haldane map
function. Each gamete draws from its own deterministic `dqrng`
sub-stream keyed on `(base seed, offspring index, parent role)`, so
output is identical for any `batch_size` (and, in later versions, any
thread count) and is set by `seed` (see above). The crossover count per
chromosome ~ Poisson(chr_len_cM\[i\] / 100), i.e. Poisson(genetic length
in Morgans) since 1 Morgan = 100 cM, drawn by a uniform-only inversion
sampler. Genetic positions come from the resolved genetic map
(`genome_map`, per the gamete-producing parent's sex/line); crossover
positions are uniform in cM within each chromosome, and the starting
haplotype is chosen at random. This applies to plain autosomes (the
seeded default: `chr_inheritance` `from_parent_1 = from_parent_2 = 1`
for both sexes and `chr_recombination` `recombines = TRUE` for both
parent sexes). Sex chromosomes and organelles configured via
[`define_chromosome()`](https://austin-putz.github.io/tidybreed/reference/define_chromosome.md)
follow their own inheritance rule instead: the contributing parent's
copy recombines via the same Haldane model only when that parent carries
two copies of the chromosome AND the producing parent's
`chr_recombination.recombines` resolves `TRUE` for its sex; otherwise
the parent's single stored copy is passed straight through unchanged (no
recombination, no additional RNG draws).

**Offspring IDs:** IDs follow the same `"{line_name}_{n}"` format as
[`add_founders()`](https://austin-putz.github.io/tidybreed/reference/add_founders.md).
Numbering continues from the current maximum for each line.

## Examples

``` r
if (FALSE) { # \dontrun{
pop <- open_pop(pop_name = "test", db_name = ":memory:") |>
  define_genome(n_loci = 1000, n_chr = 5, chr_len_Mb = 100) |>
  define_founder_haplotypes(n_haplotypes = 200, line_name = "A")
pop <- pop |>
  get_table("founder_haplotypes") |>
  add_founders(n_males = 5, n_females = 5, line_name = "A", gen = 1L)

# One offspring per mating, extra metadata column (gen)
matings <- tibble::tibble(
  id_parent_1 = rep("A_1", 5),
  id_parent_2 = paste0("A_", 6:10),
  sex         = c("M", "F", "M", "F", "M"),
  line_name   = "A",
  gen         = 2L
)
pop <- pop |> add_offspring(matings)

# Animal-breeder-style aliases
matings2 <- tibble::tibble(
  id_sire   = rep("A_1", 3),
  id_dam    = paste0("A_", 6:8),
  sex       = c("M", "F", "M"),
  line_name = "A",
  gen       = 2L
)
pop <- pop |> add_offspring(matings2)
} # }
```
