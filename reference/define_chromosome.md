# Define a chromosome's inheritance or recombination rule

Sets a non-default rule for one chromosome — sex chromosomes (X/Y, Z/W,
X0/Z0) and organelles (mitochondria, plastids).
[`define_genome()`](https://austin-putz.github.io/tidybreed/reference/define_genome.md)
seeds every chromosome with the default diploid-autosome rule (one copy
from each parent, recombining in both sexes); call `define_chromosome()`
before
[`add_founders()`](https://austin-putz.github.io/tidybreed/reference/add_founders.md)
to override that default.

The rule is stored across two explicit long tables, and **each call sets
exactly one concern**:

- **Inheritance** (copy count): supply `from_parent_1` and
  `from_parent_2` — the number of copies of this chromosome an offspring
  of `offspring_sex` inherits from parent_1 (sire) and parent_2 (dam).
  Writes `chr_inheritance`.

- **Recombination**: supply `recombines` — whether a parent of
  `parent_sex` recombines this chromosome when making gametes. Writes
  `chr_recombination`.

Copy counts are **absolute** (correct at ploidy 2, which this version
enforces): an autosome is `from_parent_1 = 1, from_parent_2 = 1`; a
male's X is `0, 1`; a male's Y is `1, 0`; a female's Y is `0, 0`;
maternal mitochondria are `0, 1`. The two concerns are kept in separate
calls because "which sex" means the **offspring** for copy count but the
**producing parent** for recombination — mixing them in one call would
make that ambiguous.

## Usage

``` r
define_chromosome(
  pop,
  chr_name,
  offspring_sex = NULL,
  parent_sex = NULL,
  line_name = NULL,
  from_parent_1 = NULL,
  from_parent_2 = NULL,
  recombines = NULL,
  overwrite = TRUE
)
```

## Arguments

- pop:

  A `tidybreed_pop` object.

- chr_name:

  Character scalar. Must already exist as a `genome_meta.chr_name` value
  (defined via
  [`define_genome()`](https://austin-putz.github.io/tidybreed/reference/define_genome.md)'s
  `chr_names` argument).

- offspring_sex:

  `NULL` (default; the sex-agnostic default row), `"M"`, or `"F"`. The
  **inheritance** row key — the offspring sex whose copy counts this
  call sets. Only valid on an inheritance call.

- parent_sex:

  `NULL` (default; both parent sexes), `"M"`, or `"F"`. The
  **recombination** row key — the producing-parent sex whose
  recombination this call sets. Only valid on a recombination call.

- line_name:

  `NULL` (default; all lines) or a line name. Reserved for line-specific
  rules (crossbreeding); resolves with the same precedence as the sex
  dimension.

- from_parent_1, from_parent_2:

  Non-negative integers supplied **together** on an inheritance call.
  Copies inherited from parent_1 (sire) and parent_2 (dam); their sum
  must be `<= 2` in this diploid release.

- recombines:

  Single logical supplied on a recombination call. `TRUE` if the
  chromosome recombines during gamete formation, `FALSE` for
  non-recombining chromosomes (Y, W, MT, most organelles). For a
  genome-wide achiasmatic sex, set `recombines_M`/`recombines_F` on
  [`define_genome()`](https://austin-putz.github.io/tidybreed/reference/define_genome.md)
  instead of one call per chromosome.

- overwrite:

  Logical. `TRUE` (default) upserts the row for this logical key.
  `FALSE` errors if the exact `(chr_name, sex, line_name)` key already
  exists (note:
  [`define_genome()`](https://austin-putz.github.io/tidybreed/reference/define_genome.md)
  pre-seeds the `(chr, NULL, NULL)` default row, so `overwrite = FALSE`
  on a default-keyed call always collides — it is meant for guarding
  *new* sex/line-specific rows).

## Value

The `tidybreed_pop` (invisibly). Assign the result back.

## See also

[`define_genome()`](https://austin-putz.github.io/tidybreed/reference/define_genome.md)

## Examples

``` r
if (FALSE) { # \dontrun{
pop <- open_pop(pop_name = "test", db_name = ":memory:") |>
  define_genome(n_loci = 1000, n_chr = 3, chr_len_Mb = 100,
                chr_names = c("1", "X", "Y"))

# Mammalian sex chromosomes — inheritance (override only the deviating sexes)
pop <- pop |>
  define_chromosome("X", offspring_sex = "M", from_parent_1 = 0, from_parent_2 = 1) |>
  define_chromosome("Y", offspring_sex = "M", from_parent_1 = 1, from_parent_2 = 0) |>
  define_chromosome("Y", offspring_sex = "F", from_parent_1 = 0, from_parent_2 = 0) |>
  # Y is non-recombining (both parent sexes)
  define_chromosome("Y", recombines = FALSE)

# Avian Z/W — females (ZW) are the heterogametic sex, males are ZZ
pop_zw <- open_pop(pop_name = "chicken", db_name = ":memory:") |>
  define_genome(n_loci = 1000, n_chr = 3, chr_len_Mb = 100,
                chr_names = c("1", "Z", "W"))
pop_zw <- pop_zw |>
  define_chromosome("Z", offspring_sex = "F", from_parent_1 = 1, from_parent_2 = 0) |>
  define_chromosome("W", offspring_sex = "F", from_parent_1 = 0, from_parent_2 = 1) |>
  define_chromosome("W", offspring_sex = "M", from_parent_1 = 0, from_parent_2 = 0) |>
  # W is non-recombining (both parent sexes)
  define_chromosome("W", recombines = FALSE)

# Mitochondria — maternal inheritance + non-recombining (separate calls)
pop2 <- open_pop(pop_name = "test2", db_name = ":memory:") |>
  define_genome(n_loci = 1000, n_chr = 2, chr_len_Mb = 100,
                chr_names = c("1", "MT"))
pop2 <- pop2 |>
  define_chromosome("MT", from_parent_1 = 0, from_parent_2 = 1) |>
  define_chromosome("MT", recombines = FALSE)
} # }
```
