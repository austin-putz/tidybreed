# Materialize genotype dosage values into `ind_genotype`

Computes 0/1/2 genotype dosages (sum of alleles across an individual's
haplotype strands) from the long `ind_haplotype` table and writes them
to the on-demand `ind_genotype` cache. `ind_genotype` is **never**
auto-populated by
[`add_founders()`](https://austin-putz.github.io/tidybreed/reference/add_founders.md)
or
[`add_offspring()`](https://austin-putz.github.io/tidybreed/reference/add_offspring.md)
— call `add_dosage()` explicitly when you need dosage values for
marker-assisted selection, allele-frequency queries, or other downstream
analysis.

Pipe a `tidybreed_table` (from
[`get_table()`](https://austin-putz.github.io/tidybreed/reference/get_table.md)
and optionally
[`dplyr::filter()`](https://dplyr.tidyverse.org/reference/filter.html))
as the first argument to select individuals. As with
[`add_tbv()`](https://austin-putz.github.io/tidybreed/reference/add_tbv.md),
the **distinct** `id_ind` values in the collected table form the
candidate set, so a table with multiple rows per individual (e.g.
`ind_phenotype`) does not multiply work.

## Usage

``` r
add_dosage(tbl, chip_name = NULL, locus_names = NULL, overwrite_dosage = FALSE)
```

## Arguments

- tbl:

  A `tidybreed_table` from
  [`get_table()`](https://austin-putz.github.io/tidybreed/reference/get_table.md)
  (optionally filtered). Any table with an `id_ind` column is accepted.

- chip_name:

  Character or `NULL`. Name of a chip defined via
  [`define_chip()`](https://austin-putz.github.io/tidybreed/reference/define_chip.md)
  (which writes `is_<chip_name>` to `genome_meta`). Loci are restricted
  to `is_<chip_name> = TRUE`. Errors if the column is missing.

- locus_names:

  Character vector or `NULL`. Explicit loci to materialize. Takes
  precedence over `chip_name` when both are supplied.

- overwrite_dosage:

  Logical. When `TRUE`, existing `ind_genotype` rows for the candidate
  individuals are deleted before inserting (cache-scope reset). When
  `FALSE` (default), rows are upserted via `INSERT OR REPLACE` — dosage
  is fully determined by the haplotypes, so re-running is idempotent.

## Value

The modified `tidybreed_pop` (invisibly).

## `add_dosage()` vs. [`add_genotypes()`](https://austin-putz.github.io/tidybreed/reference/add_genotypes.md)

These are different operations despite similar names.
[`add_genotypes()`](https://austin-putz.github.io/tidybreed/reference/add_genotypes.md)
marks animals as *physically genotyped* on a chip by writing a BOOLEAN
`has_<chip>` column to `ind_meta`; it touches no dosage data.
`add_dosage()` **materializes simulated dosage values** (ground truth
from `ind_haplotype`) into `ind_genotype`.

## See also

[`add_genotypes()`](https://austin-putz.github.io/tidybreed/reference/add_genotypes.md),
[`extract_genotypes()`](https://austin-putz.github.io/tidybreed/reference/extract_genotypes.md),
[`define_chip()`](https://austin-putz.github.io/tidybreed/reference/define_chip.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# Dosage for all generation-5 candidates on the 50K chip
pop <- pop |>
  get_table("ind_meta") |>
  dplyr::filter(gen == 5L) |>
  add_dosage(chip_name = "50K")

# Dosage at specific QTL only
pop <- pop |>
  get_table("ind_meta") |>
  add_dosage(locus_names = c("Locus_1", "Locus_42"))

# Marker-assisted selection: homozygous-favorable at a QTL
ids <- pop |>
  get_table("ind_genotype") |>
  dplyr::filter(locus_name == "Locus_1", dosage_value == 2L) |>
  dplyr::pull(id_ind)
} # }
```
