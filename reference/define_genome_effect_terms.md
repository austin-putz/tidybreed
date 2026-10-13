# Define genome effects from coefficients you supply

Writes genome-effect terms whose coefficients you already know (GWAS
estimates, a published QTL map, a hand-built test model) for one trait.
Each term is one coefficient (`genome_effects`) over one or more loci
(`genome_effect_members`), each locus optionally restricted to allele
copies from a given line and/or parent (`genome_effect_member_origins`).
[`add_tgv()`](https://austin-putz.github.io/tidybreed/reference/add_tgv.md)
then evaluates every term into `ind_tgv`.

Nothing is sampled or scaled: the coefficients are written as given, so
the trait has whatever variance they imply. To sample effects calibrated
to a target variance use
[`define_additive_effects()`](https://austin-putz.github.io/tidybreed/reference/define_additive_effects.md)
or
[`define_genome_effects()`](https://austin-putz.github.io/tidybreed/reference/define_genome_effects.md).
Those generators write the reserved owner `"generated"` through the same
internal engine; this function refuses that owner, so your terms and
generated ones never replace each other.

Most models are easiest to build with the helpers, which return `terms`
rows that [`rbind()`](https://rdrr.io/r/base/cbind.html) together:

- [`ad_terms()`](https://austin-putz.github.io/tidybreed/reference/ad_terms.md):
  additive and dominance coefficients per locus;

- [`aa_terms()`](https://austin-putz.github.io/tidybreed/reference/aa_terms.md):
  additive-by-additive coefficients per locus pair;

- [`genotype_terms()`](https://austin-putz.github.io/tidybreed/reference/genotype_terms.md):
  any genotypic value table over one or more loci.

## Usage

``` r
define_genome_effect_terms(
  pop,
  trait_name,
  terms,
  effect_owner = "custom",
  mode = c("append", "replace_scope", "replace_owner", "replace_trait"),
  origin = NULL,
  base_tbl = NULL,
  require_complete = FALSE
)
```

## Arguments

- pop:

  A `tidybreed_pop`.

- trait_name:

  Character scalar: an existing trait (from
  [`define_trait()`](https://austin-putz.github.io/tidybreed/reference/define_trait.md)).

- terms:

  A data frame of terms; see *The `terms` data frame*. Usually built
  with
  [`ad_terms()`](https://austin-putz.github.io/tidybreed/reference/ad_terms.md),
  [`aa_terms()`](https://austin-putz.github.io/tidybreed/reference/aa_terms.md)
  and
  [`genotype_terms()`](https://austin-putz.github.io/tidybreed/reference/genotype_terms.md).

- effect_owner:

  Character scalar labelling who owns these terms, so a later call can
  replace them without touching anyone else's. Owners never compete: all
  owners' terms add. Default `"custom"`; `"generated"` is reserved for
  the generators and refused.

- mode:

  `"append"` (default), `"replace_scope"`, `"replace_owner"` or
  `"replace_trait"`; see *Replacement modes*.

- origin:

  `NULL` (default, every allele copy), a named list of scalars or a
  per-member data frame; see *Scope (`origin`)*.

- base_tbl:

  Optional `tidybreed_table` from
  [`get_table()`](https://austin-putz.github.io/tidybreed/reference/get_table.md)
  (optionally filtered) whose allele frequencies fill any missing
  `center_value` on `additive` and `dominance` members:
  `founder_haplotypes`, `ind_haplotype`, or any table with an `id_ind`
  column, from the same `pop` (see
  [`extract_allele_freq()`](https://austin-putz.github.io/tidybreed/reference/extract_allele_freq.md)).
  Given centres are never overwritten, and `indicator` members have
  none. `NULL` (default) fills nothing, so a missing centre is an error.
  The fill is Cockerham `p`; functional coding writes `0.5` itself (as
  [`ad_terms()`](https://austin-putz.github.io/tidybreed/reference/ad_terms.md)
  does). One call gives one `p` per locus, so write line-specific terms
  one line per call.

- require_complete:

  Logical. `TRUE` requires an `indicator` table to cover every reachable
  genotype state of its loci (including `copy_count_value = 0` where a
  chromosome can be absent). Default `FALSE`: a state you leave out
  contributes zero.

## Value

The `tidybreed_pop`, invisibly. Run
[`add_tgv()`](https://austin-putz.github.io/tidybreed/reference/add_tgv.md)
to evaluate the terms; `get_table(pop, "genome_effect_terms")` lists
them.

## The `terms` data frame

One row per (term x locus), the long shape every
[`get_table()`](https://austin-putz.github.io/tidybreed/reference/get_table.md)
table has.

|  |  |  |
|----|----|----|
| Column | Required | Meaning |
| `term_id` | yes, unless one term | Groups rows into one term. **Input only**: never stored; the writer assigns an `id_genome_effect`. Any atomic type |
| `genome_value` | yes | The term's coefficient. Constant within a `term_id` |
| `effect_name` | no | A free label per term (constant within a `term_id`); no mathematical meaning |
| `locus_name` | yes | A locus in `genome_meta`; each locus at most once per term |
| `contrast_name` | yes | `"additive"`: each counted allele-1 copy minus `center_value` (`g - 2 * center_value` over both copies). `"dominance"`: the HWE dominance contrast at `center_value`. `"indicator"`: 1 when the locus is in the state below, else 0 |
| `center_value` | `additive`, `dominance` | The centre: the allele-1 frequency `p` for Cockerham (statistical) coding, `0.5` for functional coding. May be left `NA` with `base_tbl` |
| `copy_count_value`, `dosage_value` | `indicator` | The genotype state: copies present and copies of allele 1. `copy_count_value` may be left `NA` at ordinary diploid autosomal loci (it is 2) |

A term with several members is their product, so a two-locus `additive`
term is an A x A interaction and a two-locus `indicator` term is one
cell of a two-locus genotype table.

## Scope (`origin`)

- `NULL` (default): the common scope, every allele copy.

- A **named list of scalars** applied to every member of every term,
  with names among `line_name`, `parent_origin` (`1` sire, `2` dam),
  `line_match_type` and `copy_count`. `list(line_name = "Duroc")` is
  copies from line Duroc (`line_match_type = "exact"`);
  `list(parent_origin = 1)` is copies from the sire, any line
  (`line_match_type = "any"`).

- A **data frame** with columns `term_id`, `locus_name`,
  `line_match_type` (`"exact"`, `"unknown"` or `"any"`), `line_name`,
  `parent_origin`, `copy_count`: one scope per member, needed for the
  copy combinations a `dominance` or `indicator` member takes, such as
  reciprocal dominance.

An `additive` member takes at most one origin row. A `dominance` or
`indicator` member takes one row per kind of copy, whose `copy_count`s
sum to the state's copy count (2 for `dominance`), so a scalar list is
usually not enough for one; the error names the locus.

Where a copy matches both a line-specific term and a common one in the
same family (same trait, owner, loci and contrasts), the more specific
term applies to it, not both; terms of different owners always add.

## Replacement modes

- `"append"` (default): insert; a term with the same family and scope as
  an existing one is refused.

- `"replace_scope"`: delete this owner's terms for the trait whose scope
  equals `origin`, then insert. Needs `origin` to be `NULL` or a scalar
  list.

- `"replace_owner"`: delete all of this owner's terms for the trait,
  then insert.

- `"replace_trait"`: delete every owner's terms for the trait, then
  insert. Refused while the trait has generated terms: those are
  replaced only by re-running a generator, or removed with
  [`remove_generated_effects()`](https://austin-putz.github.io/tidybreed/reference/remove_generated_effects.md).

Rows are never deleted from the genome-effect tables any other way.

## See also

[`ad_terms()`](https://austin-putz.github.io/tidybreed/reference/ad_terms.md),
[`aa_terms()`](https://austin-putz.github.io/tidybreed/reference/aa_terms.md)
and
[`genotype_terms()`](https://austin-putz.github.io/tidybreed/reference/genotype_terms.md)
build `terms`;
[`define_additive_effects()`](https://austin-putz.github.io/tidybreed/reference/define_additive_effects.md)
and
[`define_genome_effects()`](https://austin-putz.github.io/tidybreed/reference/define_genome_effects.md)
for sampled, calibrated effects;
[`add_tgv()`](https://austin-putz.github.io/tidybreed/reference/add_tgv.md).

## Examples

``` r
set.seed(1)
pop <- open_pop(pop_name = "terms_example", db_name = ":memory:") |>
  define_genome(n_loci = 60, n_chr = 2, chr_len_Mb = 100) |>
  define_founder_haplotypes(n_haplotypes = 100, line_name = "Duroc")
#> duckdb keeps downloaded extensions and secrets in a temporary directory:
#> ℹ /tmp/RtmplKq43k/duckdb
#> This is removed when the R session ends.
#> • Extensions are re-downloaded each session.
#> • Secrets are lost.
#> ℹ Run duckdb(shared_home = TRUE) (or create ~/.duckdb) to keep them (suitable for most users).
#> ℹ Run duckdb(shared_home = FALSE) to accept the temporary directory (and silence this message).
#> ℹ See ?duckdb_storage for details and alternatives.
#> Opened in-memory population 'terms_example'
#> Defined genome: 60 loci across 2 chromosomes | chr lengths (Mb): all equal to 100 Mb
#>   Tables written: genome_meta, genome_map, ind_haplotype, ind_genotype, ind_crossover, chr_inheritance, chr_recombination, genome_effects, genome_effect_members, genome_effect_member_origins
#>   Views created: genome_effect_terms, genome_effect_loci
#> Generating 100 founder haplotypes (method = "uniform")...
#>   Created founder_haplotypes for line 'Duroc' (100 haplotypes x 60 loci)
pop <- pop |>
  get_table("founder_haplotypes") |>
  add_founders(n_males = 10, n_females = 10, line_name = "Duroc")
#> Added 20 founders (10 males, 10 females) to line 'Duroc'
pop <- define_trait(pop, "ADG")
#> Added trait 'ADG'.

# Additive and dominance effects of two loci, functional coding
# (a = half the homozygote difference, d = heterozygote deviation)
pop <- pop |>
  define_genome_effect_terms("ADG",
    ad_terms(c("Locus_10", "Locus_20"), a = c(0.4, -0.2), d = c(0.1, 0),
             p = c(0.3, 0.5)),
    effect_owner = "qtl_map")
#> Implied genetic mean (functional coding), reported only — written to no table: Locus_10 mu = -0.118; Locus_20 mu =  0.000. Running total over these single-locus main effects: -0.118.
#> Wrote 3 genome-effect terms (3 member rows) for trait 'ADG' under owner 'qtl_map' [append].

# An A x A interaction between them, added under its own owner
pop <- pop |>
  define_genome_effect_terms("ADG",
    aa_terms("Locus_10", "Locus_20", e = 0.05, p_1 = 0.3, p_2 = 0.5),
    effect_owner = "epistasis")
#> Implied genetic mean (functional coding), reported only — written to no table: (Locus_10, Locus_20) mu = 0. Running total over these pairs: 0.
#> Wrote 1 genome-effect term (2 member rows) for trait 'ADG' under owner 'epistasis' [append].

# The same kind of rows written directly: one Cockerham additive term
# whose centre is filled from the founders' allele frequency
pop <- pop |>
  define_genome_effect_terms("ADG",
    data.frame(locus_name = "Locus_30", contrast_name = "additive",
               genome_value = 0.25),
    effect_owner = "candidate_gene",
    base_tbl = get_table(pop, "founder_haplotypes"))
#> Wrote 1 genome-effect term (1 member rows) for trait 'ADG' under owner 'candidate_gene' [append].

# An effect only on copies inherited from the sire (imprinting)
pop <- pop |>
  define_genome_effect_terms("ADG",
    ad_terms("Locus_40", a = 0.3, p = 0.5),
    effect_owner = "imprinted", origin = list(parent_origin = 1))
#> Implied genetic mean (functional coding), reported only — written to no table: Locus_40 mu = 0. Running total over these single-locus main effects: 0.
#> Wrote 1 genome-effect term (1 member rows) for trait 'ADG' under owner 'imprinted' [append].

# Reciprocal dominance: the heterozygote value depends on which parent
# gave the Duroc copy and which the Landrace copy
pop <- pop |>
  define_genome_effect_terms("ADG",
    terms = data.frame(term_id = 1L, locus_name = "Locus_50",
                       contrast_name = "dominance", center_value = 0.5,
                       genome_value = 1.2),
    origin = data.frame(term_id = 1L, locus_name = "Locus_50",
                        line_match_type = "exact",
                        line_name = c("Duroc", "Landrace"),
                        parent_origin = c(1L, 2L), copy_count = 1L),
    effect_owner = "reciprocal")
#> Wrote 1 genome-effect term (1 member rows) for trait 'ADG' under owner 'reciprocal' [append].

# Change your mind: replace everything under one owner
pop <- pop |>
  define_genome_effect_terms("ADG",
    ad_terms(c("Locus_10", "Locus_20"), a = c(0.5, -0.1), p = c(0.3, 0.5)),
    effect_owner = "qtl_map", mode = "replace_owner")
#> Implied genetic mean (functional coding), reported only — written to no table: Locus_10 mu = -0.2; Locus_20 mu =  0.0. Running total over these single-locus main effects: -0.2.
#> Wrote 2 genome-effect terms (2 member rows) for trait 'ADG' under owner 'qtl_map' [replace_owner].

get_table(pop, "genome_effect_terms")
#> <tidybreed_table: genome_effect_terms>  [6 rows × 9 fields]
#> # A tibble: 6 × 9
#>   trait_name effect_owner   effect_name id_genome_effect effect_order
#>   <chr>      <chr>          <chr>                  <int>        <dbl>
#> 1 ADG        imprinted      NA                         6            1
#> 2 ADG        reciprocal     NA                         7            1
#> 3 ADG        epistasis      NA                         4            2
#> 4 ADG        candidate_gene NA                         5            1
#> 5 ADG        qtl_map        NA                         8            1
#> 6 ADG        qtl_map        NA                         9            1
#> # ℹ 4 more variables: contrast_signature <chr>, family_key <chr>,
#> #   scope_description <chr>, genome_value <dbl>
pop <- pop |> get_table("ind_meta") |> add_tgv("ADG")
#> Computed TGV for 20 individuals on trait 'ADG' (additive, interaction).
close_pop(pop)
```
