# Genetic models: writing, generating and measuring QTL effects

A trait’s genetics in `tidybreed` is a set of **genome-effect terms**:
additive, dominance and additive-by-additive (A x A) coefficients on
QTL, stored in the `genome_effects` tables and evaluated into `ind_tgv`
by
[`add_tgv()`](https://austin-putz.github.io/tidybreed/reference/add_tgv.md).
There are four ways in, one per job:

| Job | Function |
|----|----|
| You know the coefficients | [`define_genome_effect_terms()`](https://austin-putz.github.io/tidybreed/reference/define_genome_effect_terms.md) with [`ad_terms()`](https://austin-putz.github.io/tidybreed/reference/ad_terms.md), [`aa_terms()`](https://austin-putz.github.io/tidybreed/reference/aa_terms.md), [`genotype_terms()`](https://austin-putz.github.io/tidybreed/reference/genotype_terms.md) |
| Sample additive effects that hit a target `G` exactly | [`define_additive_effects()`](https://austin-putz.github.io/tidybreed/reference/define_additive_effects.md) |
| Sample additive + dominance + A x A effects that hit `G_A`, `G_D`, `G_AA` exactly | [`define_genome_effects()`](https://austin-putz.github.io/tidybreed/reference/define_genome_effects.md) |
| Measure what a population actually has | [`extract_genetic_variance()`](https://austin-putz.github.io/tidybreed/reference/extract_genetic_variance.md) |

The two generators sample **and** calibrate: their effects are written
under the reserved owner `"generated"`, which means “scaled to the
stored target”. Known coefficients always go through the writer, under
your own owner.

## A small population

Everything below runs on a 60-locus genome held in memory.

``` r

library(tidybreed)
library(dplyr)

set.seed(1)
pop <- open_pop(pop_name = "models", db_name = ":memory:") |>
  define_genome(n_loci = 60, n_chr = 2, chr_len_Mb = 100) |>
  define_founder_haplotypes(n_haplotypes = 400)
pop <- pop |> get_table("founder_haplotypes") |>
  add_founders(n_males = 100, n_females = 100, line_name = "A")
for (t in c("ADG", "BF", "WW", "FI", "LS")) pop <- define_trait(pop, t)
```

## 1. Known coefficients: the writer

With functional coding, `a` is half the homozygote difference, `d` the
heterozygote’s deviation and `e` the A x A coefficient on the centred
dosages. `p` is required: it is the basis of the reported mean. The
builders only build rows;
[`define_genome_effect_terms()`](https://austin-putz.github.io/tidybreed/reference/define_genome_effect_terms.md)
writes them.

``` r

pop <- pop |>
  define_genome_effect_terms(
    trait_name = "ADG",
    terms = rbind(
      ad_terms(locus_name = c("Locus_10", "Locus_44"),
               a = c(0.30, -0.12), d = c(0.10, 0.05), p = c(0.35, 0.60),
               coding = "functional"),
      aa_terms(locus_name_1 = "Locus_10", locus_name_2 = "Locus_44",
               e = 0.08, p_1 = 0.35, p_2 = 0.60, coding = "functional")))
#> Implied genetic mean (functional coding), reported only — written to no table: Locus_10 mu = -0.0445; Locus_44 mu =  0.0000. Running total over these single-locus main effects: -0.0445.
#> Implied genetic mean (functional coding), reported only — written to no table: (Locus_10, Locus_44) mu = -0.0048. Running total over these pairs: -0.0048.
#> Wrote 5 genome-effect terms (6 member rows) for trait 'ADG' under owner 'custom' [append].

pop <- pop |> get_table("ind_meta") |> add_tgv("ADG")
#> Computed TGV for 200 individuals on trait 'ADG' (additive, indicator, interaction).
get_table(pop, "ind_tgv") |> collect() |> arrange(id_ind) |> head(6)
#> # A tibble: 6 × 5
#>   id_tgv id_ind trait_name component_name tgv_value
#>    <int> <chr>  <chr>      <chr>              <dbl>
#> 1    130 A_1    ADG        additive            0.42
#> 2    446 A_1    ADG        interaction        -0.08
#> 3    379 A_1    ADG        indicator           0   
#> 4     12 A_10   ADG        additive            0.12
#> 5    557 A_10   ADG        interaction         0   
#> 6    410 A_10   ADG        indicator           0.1
```

Each individual gets one value per component; the `ind_tgv_total` view
adds them. Functional coding stores the heterozygote deviation as an
`indicator` term on the heterozygous genotype, so the components are the
functional ones (`additive`, `indicator`, `interaction`), not breeding
values: the breeding value of this model also carries part of `d` and
`e`. Nothing here is calibrated to a variance: a hand-written model has
the variance its coefficients give.

## 2. Additive effects to a target: `define_additive_effects()`

Pass a target `G` and the sampled effects are rescaled so the additive
covariance is `G` **exactly** under the anchor: `"genic"` (the HWE
expectation at the base allele frequencies, default) or `"realised"`
(the base individuals’ own genotypes). The target is stored in
`trait_var_comp` with the effects, in one transaction.

``` r

G <- matrix(c(1.0, 0.4,
              0.4, 2.0), 2, dimnames = list(c("BF", "WW"), c("BF", "WW")))
set.seed(2)
pop <- get_table(pop, "genome_meta") |>
  filter(locus_id <= 40) |>
  define_additive_effects(c("BF", "WW"), G = G)
#> The founder pool (pool expectation under random pairing) sees relative spectrum [0.851, 1.01] of the target (sampling LD of 400 haplotypes).
#> Set correlated additive effects for traits: BF, WW (method: shared; base: founder_haplotypes; scope: all lines, both parents' copies): exact under the genic anchor; delivered [BF,BF = 1; BF,WW = 0.4; WW,WW = 2].
```

The message reports what other populations would see (here the founder
pool, with the sampling LD of 400 haplotypes). `line_name =` writes a
line-specific variant calibrated to that line’s own target, for
crossbreeding designs.

## 3. Additive, dominance and epistasis: `define_genome_effects()`

The same idea with three targets. Dominance degrees are drawn around
`dominance_degree_mean` (0.19 by default), A x A effects sit on random
pairs of QTL (or the `pairs` you supply), and all three blocks are
calibrated together, highest order first.

``` r

set.seed(3)
pop <- get_table(pop, "genome_meta") |>
  define_genome_effects("FI", G_A = 1, G_D = 0.3, G_AA = 0.2)
#> Drew 30 random A x A pairs: every QTL paired once (floor(60 / 2)), as AlphaSimR does. Pass n_pairs for fewer pairs, or pairs for a chosen design (a locus may then appear in several pairs).
#> Set genome effects for 60 QTL and 30 pairs on trait FI (base: founder_haplotypes): exact under the genic anchor; delivered additive variance 1; dominance variance 0.3; additive-by-additive variance 0.2.
#> Additive floor for this sampled architecture (the least additive variance its dominance and A x A effects allow; another draw gives another floor): FI = 0.3706. G_A - floor has smallest eigenvalue 0.6294 on the correlation scale.
#> Inbreeding depression (drop in the mean per unit of F, positive = depression; sum 2pq d under the genic anchor's frequencies): implied FI = 3.114.
#> Functional effects (sampling: dominance_degree_mean = 0.19, dominance_degree_sd = 0.097): FI: 60 QTL with a != 0, 60 with d != 0, 30 pairs with e != 0.
#> The founder pool (pool expectation under random pairing; additive block only, dominance and additive-by-additive not compared) sees relative spectrum [0.992, 0.992] of the target (sampling LD of 400 haplotypes).
```

The messages report the additive floor of this draw (below), the
inbreeding depression the dominance effects imply
(`inbreeding_depression =` sets it instead), and the founder pool’s
view. A re-run replaces the trait’s `"generated"` terms and never
touches terms you wrote yourself.

**The additive floor.** Dominance and epistasis induce additive effects
of their own (an allele substitution changes the dominance and pair
values too), so the additive variance cannot be made arbitrarily small.
A `G_A` below what the drawn architecture can cancel is refused, with
the floor named:

``` r

set.seed(4)
pop <- get_table(pop, "genome_meta") |>
  filter(locus_id <= 10) |>
  define_genome_effects("LS", G_A = 0.001, G_D = 1)
#> Error:
#> ! `G_A` is below the additive floor for this sampled architecture under the "genic" anchor. The dominance and additive-by-additive effects already induce additive variance that the drawn additive architecture cannot cancel: at least LS = 0.451916 (the floor's diagonal; G_A - floor has smallest eigenvalue -450.9 on the correlation scale). Raise `G_A` or lower `G_D` / `G_AA`. The floor belongs to this draw: another seed gives another floor. Nothing was written.
```

The floor belongs to the sampled architecture: another seed, or more
loci, gives another.

## 4. Measuring: `extract_genetic_variance()`

The extractor decomposes the stored model on any group of individuals,
under the same two anchors, and writes nothing. Its rows use
`trait_var_comp`’s column names, so a target comparison is a join.

``` r

measured <- get_table(pop, "ind_meta") |>
  extract_genetic_variance(c("BF", "WW"), anchor = "realised")
#> Realised genetic covariances of the 200 selected individuals (frequencies and regressions are the cohort's own).
measured
#> # A tibble: 12 × 7
#>    effect_name    trait_name_1 trait_name_2 cov_value n_ind decomposition anchor
#>    <chr>          <chr>        <chr>            <dbl> <int> <chr>         <chr> 
#>  1 additive       BF           BF               0.991   200 full          reali…
#>  2 additive       BF           WW               0.274   200 full          reali…
#>  3 additive       WW           BF               0.274   200 full          reali…
#>  4 additive       WW           WW               1.65    200 full          reali…
#>  5 between_compo… BF           BF               0       200 full          reali…
#>  6 between_compo… BF           WW               0       200 full          reali…
#>  7 between_compo… WW           BF               0       200 full          reali…
#>  8 between_compo… WW           WW               0       200 full          reali…
#>  9 total          BF           BF               0.991   200 full          reali…
#> 10 total          BF           WW               0.274   200 full          reali…
#> 11 total          WW           BF               0.274   200 full          reali…
#> 12 total          WW           WW               1.65    200 full          reali…

targets <- get_table(pop, "trait_var_comp") |> collect()
inner_join(measured, targets,
           by = c("effect_name", "trait_name_1", "trait_name_2"),
           suffix = c("_measured", "_target")) |>
  select(effect_name, trait_name_1, trait_name_2,
         cov_value_measured, cov_value_target)
#> # A tibble: 4 × 5
#>   effect_name trait_name_1 trait_name_2 cov_value_measured cov_value_target
#>   <chr>       <chr>        <chr>                     <dbl>            <dbl>
#> 1 additive    BF           BF                        0.991              1  
#> 2 additive    BF           WW                        0.274              0.4
#> 3 additive    WW           BF                        0.274              0.4
#> 4 additive    WW           WW                        1.65               2
```

The generation targets were calibrated under the genic anchor at the
founder pool’s frequencies; these 200 individuals are a finite sample of
it, so the realised values differ a little.
`anti_join(targets, measured, ...)` lists targets with no measured block
(a trait without terms, say).

## What “exact” promises

- A generator’s result is exact for a feasible target **under the named
  anchor and the selected QTL**: other populations, later generations
  and other anchors see what their frequencies and LD give. Every call
  reports one such comparison.
- A line-scoped call calibrates its own variant only.
- `method = "union"` (each trait its own QTL set) delivers the variances
  exactly but the covariances only approximately, with a warning that
  gives them; `method = "shared"` is exact.
- One call writes one `parent_origin` scope.
- The extractor measures the population **you select**. Compare like
  with like: the generation anchor, base and scope.
