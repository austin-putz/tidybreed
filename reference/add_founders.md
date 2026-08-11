# Add founder individuals to population

Creates founder individuals by sampling haplotypes from the
`founder_haplotypes` table. Each founder receives two randomly sampled
haplotypes (with replacement), which are used to populate the long
`ind_haplotype` table (`ind_genotype` dosage is materialized on demand
via
[`add_dosage()`](https://austin-putz.github.io/tidybreed/reference/add_dosage.md)).

The `ind_meta` table is created (if it doesn't exist) or appended to
with the new founders. Founder individuals have `NULL` for both parent
IDs.

## Usage

``` r
add_founders(
  tbl,
  n_males,
  n_females,
  line_name,
  ploidy = 2L,
  ...,
  batch_size = NULL,
  max_batch_mem = NULL
)
```

## Arguments

- tbl:

  A `tidybreed_table` from `get_table("founder_haplotypes")` (optionally
  piped through
  [`dplyr::filter()`](https://dplyr.tidyverse.org/reference/filter.html)).
  The filtered rows supply the haplotype pool for sampling; use the
  filter to select a line-specific or custom subset.

- n_males:

  Integer. Number of male founders to create

- n_females:

  Integer. Number of female founders to create

- line_name:

  Character. Line identifier used for individual IDs. IDs are formatted
  as `"{line_name}_{number}"` (e.g., "A_1", "A_2")

- ploidy:

  Integer scalar. Genome ploidy for these founders. Must be `2` in this
  version of tidybreed (real polyploidy is not yet supported).

- ...:

  Optional named arguments for custom `ind_meta` columns, e.g.
  `gen = 0L`, `farm = "Iowa"`. Scalar values are broadcast to all new
  founders; vectors must have length `n_males + n_females`. Column types
  are inferred from the R type: use `0L` for INTEGER, `0` for DOUBLE,
  `"text"` for VARCHAR, `TRUE`/`FALSE` for BOOLEAN. Reserved column
  names (`id_ind`, `sex`, `line_name`, `ploidy`, etc.) are blocked.

- batch_size:

  Optional integer. Founders materialized + written per batch. Bounds
  peak memory at roughly `batch_size x n_loci` long rows regardless of
  the total number of founders. Any value produces byte-identical output
  for a fixed seed (only the *write* is batched, never the sampling).
  Overrides `max_batch_mem` when both are set.

- max_batch_mem:

  Optional per-batch memory budget as bytes (e.g. `512e6`) or a string
  (e.g. `"512MB"`). Used to derive `batch_size` when it is `NULL`. When
  both are `NULL` (the default), the batch size is auto-picked from
  detected available system memory (conservative fallback if
  unavailable).

## Value

The modified `tidybreed_pop` object (invisibly). **Important:** Assign
the result back to update your object: `pop <- add_founders(pop, ...)`

## Details

**Requirements:**

- The `founder_haplotypes` table must exist. Create it by calling
  [`define_founder_haplotypes()`](https://austin-putz.github.io/tidybreed/reference/define_founder_haplotypes.md)
  (after
  [`open_pop()`](https://austin-putz.github.io/tidybreed/reference/open_pop.md)
  and
  [`define_genome()`](https://austin-putz.github.io/tidybreed/reference/define_genome.md)).

**What it does:**

1.  Samples 2 haplotypes per founder from `founder_haplotypes` (with
    replacement)

2.  Creates/updates `ind_meta` table with founder metadata

3.  Populates `ind_haplotype` (long format; `line_origin` set to the
    founder's line and `strand = 1`. Row count per chromosome follows
    the resolved `chr_inheritance` `from_parent_1`/`from_parent_2` for
    the founder's sex — 2 rows/locus for a plain autosome (`1, 1`, the
    default), 1 for a hemizygous sex chromosome (e.g. `0, 1`), 0 for an
    absent chromosome (`0, 0`); see
    [`define_chromosome()`](https://austin-putz.github.io/tidybreed/reference/define_chromosome.md))

**ID Format:**

- Individual IDs: `"{line_name}_{number}"` (e.g., "A_1", "A_2", "B_1")

- Numbers are sequential within each line

- If founders already exist for a line, numbering continues from max ID

**Multiple Lines:**

- Can be called multiple times to add different lines to same database

- Each line has independent ID numbering

## Examples

``` r
if (FALSE) { # \dontrun{
pop <- open_pop(pop_name = "test", db_name = ":memory:") |>
  define_genome(n_loci = 1000, n_chr = 10, chr_len_Mb = 100) |>
  define_founder_haplotypes(n_haplotypes = 100, line_name = "A")

# Simple — all haplotypes in the table
pop <- pop |>
  get_table("founder_haplotypes") |>
  add_founders(n_males = 10, n_females = 100, line_name = "A")

# Filtered by line (line-specific pools need their own
# define_founder_haplotypes() call before they can be filtered to)
pop <- pop |>
  define_founder_haplotypes(n_haplotypes = 100, line_name = "Yorkshire")
pop <- pop |>
  get_table("founder_haplotypes") |>
  dplyr::filter(line_name == "Yorkshire") |>
  add_founders(n_males = 10, n_females = 50, line_name = "Yorkshire")

# With custom ind_meta columns via ...
pop <- pop |>
  get_table("founder_haplotypes") |>
  dplyr::filter(line_name == "A") |>
  add_founders(n_males = 10, n_females = 100, line_name = "A",
               gen = 0L, farm = "FarmA")

# Add a second line
pop <- pop |>
  define_founder_haplotypes(n_haplotypes = 50, line_name = "B")
pop <- pop |>
  get_table("founder_haplotypes") |>
  dplyr::filter(line_name == "B") |>
  add_founders(n_males = 5, n_females = 50, line_name = "B")

# View founders
get_table(pop, "ind_meta") |> dplyr::collect()
} # }
```
