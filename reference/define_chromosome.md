# Define a chromosome's inheritance or recombination rule

Changes how one chromosome is passed on, for sex chromosomes (X/Y, Z/W,
X0, Z0) and organelles (mitochondria, plastids).
[`define_genome()`](https://austin-putz.github.io/tidybreed/reference/define_genome.md)
starts every chromosome as a diploid autosome: one copy from each
parent, recombining in both sexes. Call `define_chromosome()` after
[`define_genome()`](https://austin-putz.github.io/tidybreed/reference/define_genome.md)
and **before
[`add_founders()`](https://austin-putz.github.io/tidybreed/reference/add_founders.md)**,
once per rule that differs from that default.

Each call sets **one** of two things:

- **Inheritance**: give `from_parent_1` and `from_parent_2`, the number
  of copies an offspring inherits from its sire (parent 1) and its dam
  (parent 2), optionally for one `offspring_sex`. Writes
  `chr_inheritance`.

- **Recombination**: give `recombines`, whether a parent recombines this
  chromosome when it makes gametes, optionally for one `parent_sex`.
  Writes `chr_recombination`.

They are separate calls because "which sex" means different animals: the
**offspring** for inheritance, the **parent** for recombination.

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

  A `tidybreed_pop` with a genome from
  [`define_genome()`](https://austin-putz.github.io/tidybreed/reference/define_genome.md).

- chr_name:

  Character scalar: a chromosome name from
  [`define_genome()`](https://austin-putz.github.io/tidybreed/reference/define_genome.md)'s
  `chr_names` (e.g. `"X"`).

- offspring_sex:

  Inheritance calls only: `"M"` or `"F"` for a rule that applies to
  offspring of that sex, or `NULL` (default) for both.

- parent_sex:

  Recombination calls only: `"M"` or `"F"` for a rule that applies when
  a parent of that sex makes gametes, or `NULL` (default) for both.

- line_name:

  Optional line name: the rule applies only to animals of that line (the
  offspring's line for inheritance, the parent's for recombination), for
  crossbreeding designs whose lines differ. `NULL` (default) applies to
  every line.

- from_parent_1, from_parent_2:

  Inheritance calls: whole numbers, always given together, of copies
  from the sire (parent 1) and the dam (parent 2). Their sum is 0, 1 or
  2; see the table in *Details*.

- recombines:

  Recombination calls: `TRUE` if the chromosome recombines when a parent
  makes gametes, `FALSE` if it is passed on whole (a Y, a W). To stop
  one sex recombining on every chromosome, use `recombines_M` /
  `recombines_F` in
  [`define_genome()`](https://austin-putz.github.io/tidybreed/reference/define_genome.md)
  instead.

- overwrite:

  Logical. `TRUE` (default) replaces an existing rule with the same
  chromosome, sex and line. `FALSE` refuses instead, to guard against
  redefining a rule by accident. Every chromosome already has a default
  rule (both sexes, all lines) from
  [`define_genome()`](https://austin-putz.github.io/tidybreed/reference/define_genome.md),
  so `overwrite = FALSE` is useful only for sex- or line-specific rules.

## Value

The `tidybreed_pop`, invisibly. Assign it back.

## Details

**Copy counts** are absolute numbers of copies, at most 2 in total:

|  |  |  |  |
|----|----|----|----|
| Chromosome | `offspring_sex` | `from_parent_1` (sire) | `from_parent_2` (dam) |
| Autosome (the default) | either | 1 | 1 |
| Mammalian X in males | `"M"` | 0 | 1 |
| Mammalian Y in males | `"M"` | 1 | 0 |
| Mammalian Y in females | `"F"` | 0 | 0 |
| Avian Z in females | `"F"` | 1 | 0 |
| Avian W in females | `"F"` | 0 | 1 |
| Avian W in males | `"M"` | 0 | 0 |
| Maternal mitochondria | either (`NULL`) | 0 | 1 |

Only the sexes that differ from the default need a rule: a mammalian X
in females keeps the default `1, 1`. A chromosome a parent carries in
one copy is passed on unchanged; one it carries in two copies recombines
as `chr_recombination` says.

**Which rule applies.** For an animal of sex S in line L, the most
specific stored rule wins: `(S, L)`, then `(S, any line)`, then
`(any sex, L)`, then the default `(any sex, any line)`. Inheritance
rules are looked up with the offspring's sex and line, recombination
rules with the parent's.

**Timing.** Founders and offspring are written with the rules in force
at the time, so set every inheritance rule before
[`add_founders()`](https://austin-putz.github.io/tidybreed/reference/add_founders.md).
Once individuals exist, an inheritance call that would change the copies
any of them carry is refused (nothing is written): their stored
haplotypes would no longer match, and breeding from them would fail. A
call that changes nothing for existing individuals, such as a rule for a
line with no individuals yet, is still allowed. Recombination rules can
be changed at any time; they apply to gametes made afterwards.

Every write is checked in the same transaction (one rule per key, both
sexes resolvable, at most 2 copies); a rejected call leaves the previous
rule in place. See the rules with `get_table(pop, "chr_inheritance")`
and `get_table(pop, "chr_recombination")`.

## See also

[`define_genome()`](https://austin-putz.github.io/tidybreed/reference/define_genome.md),
[`add_founders()`](https://austin-putz.github.io/tidybreed/reference/add_founders.md),
[`add_offspring()`](https://austin-putz.github.io/tidybreed/reference/add_offspring.md).

## Examples

``` r
# Mammalian X and Y: inheritance for the sexes that differ, and a
# non-recombining Y
pop <- open_pop(pop_name = "mammal", db_name = ":memory:") |>
  define_genome(n_loci = 300, n_chr = 3, chr_len_Mb = 100,
                chr_names = c("1", "X", "Y"))
#> duckdb keeps downloaded extensions and secrets in a temporary directory:
#> ℹ /tmp/RtmplKq43k/duckdb
#> This is removed when the R session ends.
#> • Extensions are re-downloaded each session.
#> • Secrets are lost.
#> ℹ Run duckdb(shared_home = TRUE) (or create ~/.duckdb) to keep them (suitable for most users).
#> ℹ Run duckdb(shared_home = FALSE) to accept the temporary directory (and silence this message).
#> ℹ See ?duckdb_storage for details and alternatives.
#> Opened in-memory population 'mammal'
#> Defined genome: 300 loci across 3 chromosomes | chr lengths (Mb): all equal to 100 Mb
#>   Tables written: genome_meta, genome_map, ind_haplotype, ind_genotype, ind_crossover, chr_inheritance, chr_recombination, genome_effects, genome_effect_members, genome_effect_member_origins
#>   Views created: genome_effect_terms, genome_effect_loci
pop <- pop |>
  define_chromosome("X", offspring_sex = "M", from_parent_1 = 0, from_parent_2 = 1) |>
  define_chromosome("Y", offspring_sex = "M", from_parent_1 = 1, from_parent_2 = 0) |>
  define_chromosome("Y", offspring_sex = "F", from_parent_1 = 0, from_parent_2 = 0) |>
  define_chromosome("Y", recombines = FALSE)
#> Defined chr_inheritance rule for 'X' (offspring_sex=M, line_name=NULL): from_parent_1=0, from_parent_2=1
#> Defined chr_inheritance rule for 'Y' (offspring_sex=M, line_name=NULL): from_parent_1=1, from_parent_2=0
#> Defined chr_inheritance rule for 'Y' (offspring_sex=F, line_name=NULL): from_parent_1=0, from_parent_2=0
#> Defined chr_recombination rule for 'Y' (parent_sex=NULL, line_name=NULL): recombines=FALSE
get_table(pop, "chr_inheritance")
#> <tidybreed_table: chr_inheritance>  [6 rows × 5 fields]
#> # A tibble: 6 × 5
#>   chr_name offspring_sex line_name from_parent_1 from_parent_2
#>   <chr>    <chr>         <chr>             <int>         <int>
#> 1 1        NA            NA                    1             1
#> 2 X        NA            NA                    1             1
#> 3 Y        NA            NA                    1             1
#> 4 X        M             NA                    0             1
#> 5 Y        M             NA                    1             0
#> 6 Y        F             NA                    0             0

# Founders are written with those rules: a male has one X copy and one Y
set.seed(1)
pop <- pop |> define_founder_haplotypes(n_haplotypes = 40)
#> Generating 40 founder haplotypes (method = "uniform")...
#>   Created founder_haplotypes (40 haplotypes x 300 loci)
pop <- pop |>
  get_table("founder_haplotypes") |>
  add_founders(n_males = 2, n_females = 2, line_name = "A")
#> Added 4 founders (2 males, 2 females) to line 'A'
close_pop(pop)

# Avian Z and W: females (ZW) are the heterogametic sex
pop <- open_pop(pop_name = "chicken", db_name = ":memory:") |>
  define_genome(n_loci = 300, n_chr = 3, chr_len_Mb = 100,
                chr_names = c("1", "Z", "W"))
#> duckdb keeps downloaded extensions and secrets in a temporary directory:
#> ℹ /tmp/RtmplKq43k/duckdb
#> This is removed when the R session ends.
#> • Extensions are re-downloaded each session.
#> • Secrets are lost.
#> ℹ Run duckdb(shared_home = TRUE) (or create ~/.duckdb) to keep them (suitable for most users).
#> ℹ Run duckdb(shared_home = FALSE) to accept the temporary directory (and silence this message).
#> ℹ See ?duckdb_storage for details and alternatives.
#> Opened in-memory population 'chicken'
#> Defined genome: 300 loci across 3 chromosomes | chr lengths (Mb): all equal to 100 Mb
#>   Tables written: genome_meta, genome_map, ind_haplotype, ind_genotype, ind_crossover, chr_inheritance, chr_recombination, genome_effects, genome_effect_members, genome_effect_member_origins
#>   Views created: genome_effect_terms, genome_effect_loci
pop <- pop |>
  define_chromosome("Z", offspring_sex = "F", from_parent_1 = 1, from_parent_2 = 0) |>
  define_chromosome("W", offspring_sex = "F", from_parent_1 = 0, from_parent_2 = 1) |>
  define_chromosome("W", offspring_sex = "M", from_parent_1 = 0, from_parent_2 = 0) |>
  define_chromosome("W", recombines = FALSE)
#> Defined chr_inheritance rule for 'Z' (offspring_sex=F, line_name=NULL): from_parent_1=1, from_parent_2=0
#> Defined chr_inheritance rule for 'W' (offspring_sex=F, line_name=NULL): from_parent_1=0, from_parent_2=1
#> Defined chr_inheritance rule for 'W' (offspring_sex=M, line_name=NULL): from_parent_1=0, from_parent_2=0
#> Defined chr_recombination rule for 'W' (parent_sex=NULL, line_name=NULL): recombines=FALSE
close_pop(pop)

# Mitochondria: one maternal copy in both sexes, never recombining
pop <- open_pop(pop_name = "mito", db_name = ":memory:") |>
  define_genome(n_loci = 200, n_chr = 2, chr_len_Mb = c(100, 0.02),
                chr_names = c("1", "MT"))
#> duckdb keeps downloaded extensions and secrets in a temporary directory:
#> ℹ /tmp/RtmplKq43k/duckdb
#> This is removed when the R session ends.
#> • Extensions are re-downloaded each session.
#> • Secrets are lost.
#> ℹ Run duckdb(shared_home = TRUE) (or create ~/.duckdb) to keep them (suitable for most users).
#> ℹ Run duckdb(shared_home = FALSE) to accept the temporary directory (and silence this message).
#> ℹ See ?duckdb_storage for details and alternatives.
#> Opened in-memory population 'mito'
#> Defined genome: 200 loci across 2 chromosomes | chr lengths (Mb): 1:100, MT:0.02
#>   Tables written: genome_meta, genome_map, ind_haplotype, ind_genotype, ind_crossover, chr_inheritance, chr_recombination, genome_effects, genome_effect_members, genome_effect_member_origins
#>   Views created: genome_effect_terms, genome_effect_loci
pop <- pop |>
  define_chromosome("MT", from_parent_1 = 0, from_parent_2 = 1) |>
  define_chromosome("MT", recombines = FALSE)
#> Defined chr_inheritance rule for 'MT' (offspring_sex=NULL, line_name=NULL): from_parent_1=0, from_parent_2=1
#> Defined chr_recombination rule for 'MT' (parent_sex=NULL, line_name=NULL): recombines=FALSE
close_pop(pop)
```
