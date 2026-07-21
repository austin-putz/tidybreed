# Define a chromosome's inheritance rule

Sets a non-default inheritance rule for one chromosome in `chr_meta` —
sex chromosomes (X/Y, Z/W, X0/Z0) and organelles (mitochondria,
plastids).
[`define_genome()`](https://austin-putz.github.io/tidybreed/reference/define_genome.md)
seeds every chromosome with the default autosome rule
(`copy_mode_M = copy_mode_F = "full"`,
`recombines_M = recombines_F = TRUE`); call `define_chr()` before
[`add_founders()`](https://austin-putz.github.io/tidybreed/reference/add_founders.md)
to override that default for a specific chromosome.

`copy_mode_M`/`copy_mode_F` are **relative to an individual's own
ploidy**, not an absolute copy count: `"full"` means the same as that
individual's ploidy for this chromosome, `"half"` means hemizygous (half
of that ploidy), and `"none"` means absent. In this version of tidybreed
every individual's ploidy is `2`, so `"full"` = 2 copies and `"half"` =
1 copy.

## Usage

``` r
define_chr(
  pop,
  chr_name,
  copy_mode_M = "full",
  copy_mode_F = "full",
  hemi_parent = NULL,
  recombines = TRUE,
  recombines_M = NULL,
  recombines_F = NULL,
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

- copy_mode_M:

  Character scalar: `"full"` (default), `"half"`, or `"none"`. Copy
  count for males, relative to ploidy.

- copy_mode_F:

  Character scalar: `"full"` (default), `"half"`, or `"none"`. Copy
  count for females, relative to ploidy.

- hemi_parent:

  `NULL`, `"parent_1"`, or `"parent_2"`. Which parent supplies the
  reduced copy when either sex's copy_mode is not `"full"`. Must be
  `NULL` when both copy_modes are `"full"`; must be
  `"parent_1"`/`"parent_2"` otherwise.

- recombines:

  Logical. Primary shorthand that sets **both** `recombines_M` and
  `recombines_F`. `TRUE` (default) if recombination occurs during gamete
  formation; `FALSE` for non-recombining chromosomes (Y, W, MT, most
  organelles).

- recombines_M, recombines_F:

  Logical or `NULL`. Per-sex recombination switches (matching the
  `copy_mode_M`/`copy_mode_F` pattern). `NULL` (default) means "use
  `recombines`"; pass one explicitly to set achiasmy in a single sex
  (e.g. *Drosophila* males: `recombines_M = FALSE`,
  `recombines_F = TRUE`).

- overwrite:

  Logical. If `TRUE` (default), re-calling `define_chr()` for the same
  `chr_name` updates the existing row in place — this is a normal edit
  workflow (e.g. fixing a typo'd `hemi_parent`), not append-only
  simulation output, so the default favors correction over silent no-op.

## Value

The `tidybreed_pop` (invisibly). Assign the result back.

## See also

[`define_genome()`](https://austin-putz.github.io/tidybreed/reference/define_genome.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# chr_names must exist in genome_meta before define_chr() can target them
pop <- open_pop(pop_name = "test", db_name = ":memory:") |>
  define_genome(n_loci = 1000, n_chr = 3, chr_len_Mb = 100,
                chr_names = c("1", "X", "Y"))

# Mammalian sex chromosomes
pop <- pop |>
  define_chr("X", copy_mode_M = "half", copy_mode_F = "full",
             hemi_parent = "parent_2", recombines = TRUE) |>
  define_chr("Y", copy_mode_M = "half", copy_mode_F = "none",
             hemi_parent = "parent_1", recombines = FALSE)

# Mitochondria (maternal, non-recombining) — needs its own chr_name
pop2 <- open_pop(pop_name = "test2", db_name = ":memory:") |>
  define_genome(n_loci = 1000, n_chr = 2, chr_len_Mb = 100,
                chr_names = c("1", "MT"))
pop2 <- pop2 |>
  define_chr("MT", copy_mode_M = "half", copy_mode_F = "half",
             hemi_parent = "parent_2", recombines = FALSE)
} # }
```
