# Define additive QTL effects for one or more traits

Selects QTL from a filtered `genome_meta` table, samples an effect
architecture, **calibrates** it to the stored additive target, and
writes one order-one `additive` term per locus and trait through the
same engine as
[`define_genome_effect_terms()`](https://austin-putz.github.io/tidybreed/reference/define_genome_effect_terms.md),
under the reserved effect owner `"generated"`. Only generators write
that owner, so effects written here and effects a user writes with
[`define_genome_effect_terms()`](https://austin-putz.github.io/tidybreed/reference/define_genome_effect_terms.md)
can never be confused for one another.
[`add_tgv()`](https://austin-putz.github.io/tidybreed/reference/add_tgv.md)
evaluates both; its `additive` component is the breeding value.

## Usage

``` r
define_additive_effects(
  tbl,
  trait_name,
  distribution = c("normal", "gamma"),
  G = NULL,
  trait_var_comp_tbl = NULL,
  anchor = c("genic", "realised"),
  method = c("shared", "union"),
  base_tbl = NULL,
  line_name = NULL,
  parent_origin = NULL,
  warn_bounds = c(0.8, 1.25),
  seed = NULL
)
```

## Arguments

- tbl:

  A `tidybreed_table` from
  [`get_table()`](https://austin-putz.github.io/tidybreed/reference/get_table.md)`("genome_meta")`
  (with an optional
  [`dplyr::filter()`](https://dplyr.tidyverse.org/reference/filter.html)).
  The filtered rows determine which loci are QTL.

- trait_name:

  Character scalar **or** vector. Name(s) of existing traits (from
  [`define_trait()`](https://austin-putz.github.io/tidybreed/reference/define_trait.md)).
  With two or more, the architecture is drawn jointly from `MVN(0, G)`,
  `G` is a `k x k` matrix across them, and `method` applies.

- distribution:

  Character. The shape of a single trait's sampled architecture before
  calibration: `"normal"` (default; normal effects) or `"gamma"`
  (gamma-distributed sizes with random signs: a few large QTL and many
  small ones). Calibration then scales the effects to the target either
  way. Ignored with several traits, which are always drawn from
  `MVN(0, G)`.

- G:

  Optional additive-genetic (co)variance target: a `k x k` matrix (named
  in `trait_name` order, or unnamed), or a single number for one trait.
  Written to `trait_var_comp` (with the call's `line_name`) in the same
  transaction as the effects, and never over a stored block. `NULL`
  reads the stored target (see *Targets*).

- trait_var_comp_tbl:

  Optional filtered `get_table(pop, "trait_var_comp")`: the stored rows
  to calibrate to. Use it to pick a block explicitly, e.g. to leave a
  stored non-additive block out or to calibrate one trait of a stored
  block alone. Not with `G`.

- anchor:

  Character. `"genic"` (default) or `"realised"`: the reference
  covariance the calibration is exact for. See *Details*.

- method:

  Character, multi-trait only. `"shared"` (default): every trait uses
  all the filtered loci as one shared QTL set, and the whole `G` is
  delivered exactly. `"union"`: each trait keeps the QTL set of its
  existing generated terms at this scope (restricted to the filtered
  loci), so the traits must already have generated effects there; each
  trait is scaled to its own variance, and the covariances are
  approximate (a warning gives the delivered values).

- base_tbl:

  Optional `tidybreed_table` from
  [`get_table()`](https://austin-putz.github.io/tidybreed/reference/get_table.md)
  (optionally filtered) selecting the allele copies that define base
  allele frequencies: `founder_haplotypes`, `ind_haplotype`, or any
  table with an `id_ind` column. Must come from the same `pop` as `tbl`.
  `NULL` (default) resolves to the founder pool of the line the effect
  applies to — see *Which population centers the effects*. Under
  `anchor = "realised"` it is required and names the individuals whose
  `Cov(X)` is the anchor.

- line_name:

  Optional character. When set, effects are scoped to allele copies of
  this genetic line: a copy whose `line_origin` matches takes these
  values, and falls back per copy to the common variant where no
  line-specific one exists. Also selects the default `base_tbl` and the
  line's own target block. `NULL` (default) means the common scope.

- parent_origin:

  Optional `1` (sire / parent_1) or `2` (dam / parent_2) — imprinting,
  restricting the term to copies inherited from that parent. `NULL`
  (default) means both parents' copies. **Per trait**: a scalar is
  recycled, a vector must match `trait_name` positionally, or name its
  entries by trait; each value must be exactly `1` or `2`. One call must
  use one origin for every trait: it calibrates against one reference
  covariance, defined for one set of inherited copies (a mixed-scope
  anchor is not supported yet). Define differently-scoped traits in
  separate calls. For imprinting that varies locus by locus, write the
  terms with
  [`define_genome_effect_terms()`](https://austin-putz.github.io/tidybreed/reference/define_genome_effect_terms.md).

- warn_bounds:

  Numeric length 2, `c(lower, upper)` with `0 < lower <= upper`, or
  `NULL` to turn the comparison off. Outside the bounds, an observed or
  genic-limit comparison warns; a founder-pool comparison adds the
  realised-anchor hint to its message. Default `c(0.8, 1.25)` (±25%,
  multiplicative).

- seed:

  Optional integer, applied with
  [`set.seed()`](https://rdrr.io/r/base/Random.html) immediately before
  the draw, after every input check and after the anchor's feasibility
  check (`rank(G) <= rank(M)`), so those refusals never touch the RNG. A
  failure that depends on the draw itself – the drawn architecture's
  rank, or a calibration that fails verification – comes after the seed
  and has consumed RNG draws. `NULL` (default) draws from the current
  stream, so a preceding
  [`set.seed()`](https://rdrr.io/r/base/Random.html) reproduces the call
  too.

## Value

The modified `tidybreed_pop` (invisibly). The effects are written to the
`genome_effects` tables under the owner `"generated"`, and a passed `G`
to `trait_var_comp`. Run
[`add_tgv()`](https://austin-putz.github.io/tidybreed/reference/add_tgv.md)
to evaluate them.

## Details

**When the result is exact.** The requested covariance `G` is delivered
exactly, `B' M B = G` to machine precision, for one trait or for several
with `method = "shared"`, when the rank is feasible:
`rank(G) <= rank(M)` (the anchor has enough independent segregating
directions at the selected loci) and `rank(G) <= rank(B0' M B0)` (the
drawn architecture does too). Each infeasibility is its own error. Rank
and positive semidefiniteness are judged on the target's **correlation**
scale, so a trait recorded in small units is never truncated away.
"Exact" is checked, not assumed: the delivered `B' M B` is compared with
`G` as stored, entry by entry, to a relative tolerance of `1e-8` on the
correlation scale; a calibration that misses it (a numerically
ill-conditioned architecture) is an error before anything is written.
The closing message says "exact" or "approximate" and gives the
delivered covariance under the anchor.

**How.** Effects are drawn first (one draw per QTL for one trait; joint
`MVN(0, G)` rows for several), and that draw is only the *architecture*
`B0`. It is then right-multiplied by a `k x k` matrix `A` so that
`B = B0 A` satisfies `B' M B = G` exactly. For one trait this is exactly
the scalar rescale `b0 * sqrt(G / sum(w b0^2))`. For two or more it also
fixes the genetic **correlations**, which a per-trait rescale cannot: at
200 QTL and a target correlation of 0.4 a scalar rescale delivers
anything from about 0.18 to 0.60. The same seed gives the same `B0`.

**The anchor `M`** is the reference-population genotype covariance the
calibration is exact for:

- `anchor = "genic"` (default): `M = diag(n_eligible * p * q)` at the
  base allele frequencies, the random-mating (HWE + linkage-equilibrium)
  limit. Use it for multi-generation studies: under random mating the
  realised covariance converges to it. With `base_tbl` the weights use
  that population's frequencies instead, e.g. one reference line.

- `anchor = "realised"`: `M = Cov(X)` of the individuals `base_tbl`
  selects, linkage disequilibrium included. The single-generation /
  clonal option. `base_tbl` must select **individuals** (a table with
  `id_ind`, not `founder_haplotypes` or `ind_haplotype`) with complete
  genotypes at the QTL, and the call must use the common scope (no
  `line_name`, no `parent_origin`). The genotype matrix is collected
  into memory, so there is a size limit; above it the call errors and
  suggests `"genic"`.

Under `method = "union"` each trait keeps its own QTL set and is scaled
by its own variance only, so the variances are exact and the covariances
are **approximate** – including a zero target covariance, which
overlapping QTL sets do not deliver. A warning gives the delivered
covariance and correlation and names `method = "shared"` as the exact
option. A trait with a positive target variance and no QTL in the call
is an error.

After calibration the delivered covariance is compared with what another
population sees: the **pool expectation** `2 Cov(H)` when the base is
the founder pool, the **observed** `Cov(X)` when it selects individuals,
and the genic limit under `anchor = "realised"`. The founder pool's
comparison is always a message: a small pool's departure is its sampling
LD, not a mistake in the call. To get the target exactly in the
founders, add them first and calibrate with `anchor = "realised"` on
them. The observed and genic-limit comparisons warn when the relative
spectrum leaves `warn_bounds`. Nothing is stored.

## Targets

The target is the population-wide (or line) `additive` block of
`trait_var_comp`, the single source of generation targets:

- pass `G` (a `k x k` matrix, or a number for one trait) to write it
  **and** calibrate to it, in the same transaction as the effects. If a
  block is already stored for any of the traits at that `line_name` the
  call errors, even for an identical matrix, and gives the
  [`remove_rows()`](https://austin-putz.github.io/tidybreed/reference/remove_rows.md)
  call;

- or leave `G = NULL` to use the stored rows, optionally chosen with
  `trait_var_comp_tbl = get_table(pop, "trait_var_comp") |> filter(...)`.
  The rows for the call's traits must form one complete symmetric block.
  A stored block that pairs one of the traits with a trait outside the
  call is an error (calibrating one trait alone would break the stored
  covariance), as is a stored `dominance` or `additive_by_additive`
  block for the traits: this generator calibrates the additive block
  only and never silently ignores a stored target. Filter them away with
  `trait_var_comp_tbl` to say so explicitly.

With `line_name = "C"` the default reads line C's block when one exists
and otherwise the population-wide one.

A new `G` must not leave other generated terms describing the old
target. A call replaces only its own scope, so `G` is refused when
generated additive terms of a line with no target of its own fell back
to the target it would write; the error gives the route (give that line
its own `G` first). Terms of another `parent_origin` at the same target
scope cannot get their own target (targets are per line), so the call
writes and warns, naming the scopes to re-run without `G`.

## Generated means calibrated

The generator always samples **and** calibrates: every term it writes
under the `"generated"` owner delivers its stored target. That is what
lets
[`define_phenotype()`](https://austin-putz.github.io/tidybreed/reference/define_phenotype.md)`(prevalence = )`
trust the stored target. It has no option to write fixed or unscaled
effects. Exact coefficients (GWAS estimates, a published QTL map, a
hand-built test) go through
[`define_genome_effect_terms()`](https://austin-putz.github.io/tidybreed/reference/define_genome_effect_terms.md)
with
[`ad_terms()`](https://austin-putz.github.io/tidybreed/reference/ad_terms.md)
under a user owner;
[`add_tgv()`](https://austin-putz.github.io/tidybreed/reference/add_tgv.md)
evaluates them like any other term, and a liability phenotype on such a
trait then needs `thresholds =` instead of `prevalence =`.

## Which population centers the effects

Base allele frequencies center the true breeding value (the Falconer
`allele - p` term) and set the genic weights `n_eligible * p q`. They
come from `base_tbl`, a filtered `tidybreed_table` whose identity says
*what kind of thing* is selected and whose
[`dplyr::filter()`](https://dplyr.tidyverse.org/reference/filter.html)
says *which* (see
[`extract_allele_freq()`](https://austin-putz.github.io/tidybreed/reference/extract_allele_freq.md)
for the three accepted shapes: the founder pool, allele copies in
`ind_haplotype`, or individuals from any table with `id_ind`).

When `base_tbl = NULL` the base is **the population the effect applies
to**, resolved with the same `line -> NULL` precedence as
`resolve_genome_map()`: a line-specific effect (`line_name = "A"`)
centers on line A's own founder pool, or on the shared
(`line_name = NULL`) pool when no named pool exists; a population-wide
effect (`line_name = NULL`) centers on the whole founder table. Only
that last case warns when the founder table holds more than one pool:
pooling divergent lines overstates within-line heterozygosity — the
Wahlund effect — so the calibration **under**-scales the effects and the
realised within-line additive variance falls short of the target. An
explicit `base_tbl` is an intentional selection and never warns — pass
`base_tbl = get_table(pop, "founder_haplotypes")` to pool on purpose,
which is how the common fallback variant of a crossbreeding model is
defined.

A selected QTL locus with no allele copies in the base is an error,
never silently centered at `p = 0`.

The centering constant is stored per member as
`genome_effect_members.center_value` and travels with its
`genome_value`, so evaluation applies each allele copy's own line's
centering — a crossbred animal's line-A alleles are centered on line A
and its line-B alleles on line B.

## Lines and line means

Line-specific effects (`line_name = "A"`) are centred on line A's own
base, so they add **no** difference between line means. Differences
between lines come from allele-frequency differences at QTL whose
effects are shared: use common effects centred on one reference line,
e.g.
`base_tbl = get_table(pop, "founder_haplotypes") |> filter(line_name == "Terminal")`.
The calibration then hits the target within that line, and the other
lines get whatever variance their frequencies give.

## Scope, and what a re-run replaces

`line_name` and `parent_origin` compose into the single origin row an
`additive` member is allowed:

|             |                 |                                               |
|-------------|-----------------|-----------------------------------------------|
| `line_name` | `parent_origin` | Stored scope                                  |
| `NULL`      | `NULL`          | no origin rows (the common scope)             |
| `"A"`       | `NULL`          | `('exact', 'A', parent NULL, copy_count = 1)` |
| `NULL`      | `1` / `2`       | `('any', NULL, parent, copy_count = 1)`       |
| `"A"`       | `1` / `2`       | `('exact', 'A', parent, copy_count = 1)`      |

Re-running replaces **only the variant at the same scope**
(`mode = "replace_scope"`), so successive common / line-A / line-B calls
each keep the others: the per-copy fallback that makes crossbred
breeding values correct depends on all of them standing. It also means
changing `parent_origin` on a re-run **adds** a variant rather than
replacing one — the two are in a containment relation and both apply, to
different copies. That is legal and rarely intended, so the function
warns on exactly that case.

## See also

[`define_genome_effects()`](https://austin-putz.github.io/tidybreed/reference/define_genome_effects.md)
for dominance and A x A effects;
[`extract_genetic_variance()`](https://austin-putz.github.io/tidybreed/reference/extract_genetic_variance.md)
to measure what a population has;
[`define_trait()`](https://austin-putz.github.io/tidybreed/reference/define_trait.md),
[`define_effect_cov_matrix()`](https://austin-putz.github.io/tidybreed/reference/define_effect_cov_matrix.md),
[`define_genome_effect_terms()`](https://austin-putz.github.io/tidybreed/reference/define_genome_effect_terms.md)
for known coefficients,
[`remove_generated_effects()`](https://austin-putz.github.io/tidybreed/reference/remove_generated_effects.md),
[`extract_allele_freq()`](https://austin-putz.github.io/tidybreed/reference/extract_allele_freq.md).

## Examples

``` r
# A small in-memory population with two founder lines
set.seed(1)
pop <- open_pop(pop_name = "dae_example", db_name = ":memory:") |>
  define_genome(n_loci = 200, n_chr = 2, chr_len_Mb = 100) |>
  define_founder_haplotypes(n_haplotypes = 100, line_name = "Duroc") |>
  define_founder_haplotypes(n_haplotypes = 100, line_name = "Landrace")
#> duckdb keeps downloaded extensions and secrets in a temporary directory:
#> ℹ /tmp/Rtmp7yuhM6/duckdb
#> This is removed when the R session ends.
#> • Extensions are re-downloaded each session.
#> • Secrets are lost.
#> ℹ Run duckdb(shared_home = TRUE) (or create ~/.duckdb) to keep them (suitable for most users).
#> ℹ Run duckdb(shared_home = FALSE) to accept the temporary directory (and silence this message).
#> ℹ See ?duckdb_storage for details and alternatives.
#> Opened in-memory population 'dae_example'
#> Defined genome: 200 loci across 2 chromosomes | chr lengths (Mb): all equal to 100 Mb
#>   Tables written: genome_meta, genome_map, ind_haplotype, ind_genotype, ind_crossover, chr_inheritance, chr_recombination, genome_effects, genome_effect_members, genome_effect_member_origins
#>   Views created: genome_effect_terms, genome_effect_loci
#> Generating 100 founder haplotypes (method = "uniform")...
#>   Created founder_haplotypes for line 'Duroc' (100 haplotypes x 200 loci)
#> Generating 100 founder haplotypes (method = "uniform")...
#>   Created founder_haplotypes for line 'Landrace' (100 haplotypes x 200 loci)
pop <- pop |>
  get_table("founder_haplotypes") |>
  dplyr::filter(line_name == "Duroc") |>
  add_founders(n_males = 20, n_females = 40, line_name = "Duroc")
#> Added 60 founders (20 males, 40 females) to line 'Duroc'
for (t in c("ADG", "BW", "FCR", "LMA", "IMF", "BF")) pop <- define_trait(pop, t)
#> Added trait 'ADG'.
#> Added trait 'BW'.
#> Added trait 'FCR'.
#> Added trait 'LMA'.
#> Added trait 'IMF'.
#> Added trait 'BF'.

# The Duroc founder pool as the base (the default would pool both lines)
duroc <- get_table(pop, "founder_haplotypes") |>
  dplyr::filter(line_name == "Duroc")

# One trait: every locus on chromosome 1 is a QTL. G writes the target to
# trait_var_comp and calibrates to it in one transaction.
pop <- pop |>
  get_table("genome_meta") |>
  dplyr::filter(chr_name == "1") |>
  define_additive_effects("ADG", G = 0.25, base_tbl = duroc, seed = 11)
#> The founder pool (pool expectation under random pairing) sees relative spectrum [0.826, 0.826] of the target (sampling LD of 100 haplotypes).
#> Set additive effects for 100 QTL on trait 'ADG' (base: founder_haplotypes [1 filter]; scope: all lines, both parents' copies): exact under the genic anchor; delivered variance 0.25.

# Two correlated traits on one shared QTL set: both variances and the
# correlation are delivered exactly
G <- matrix(c(0.25, 0.10, 0.10, 0.30), 2, 2,
            dimnames = list(c("BW", "FCR"), c("BW", "FCR")))
pop <- pop |>
  get_table("genome_meta") |>
  define_additive_effects(c("BW", "FCR"), G = G, base_tbl = duroc, seed = 12)
#> The founder pool (pool expectation under random pairing) sees relative spectrum [1.01, 1.12] of the target (sampling LD of 100 haplotypes).
#> Set correlated additive effects for traits: BW, FCR (method: shared; base: founder_haplotypes [1 filter]; scope: all lines, both parents' copies): exact under the genic anchor; delivered [BW,BW = 0.25; BW,FCR = 0.1; FCR,FCR = 0.3].

# A target stored beforehand is read back when G is omitted
pop <- pop |> define_effect_cov_matrix("additive", 0.4, trait_name = "LMA")
#> Stored 'additive' covariance matrix for: LMA.
pop <- pop |>
  get_table("genome_meta") |>
  define_additive_effects("LMA", base_tbl = duroc, seed = 13)
#> The founder pool (pool expectation under random pairing) sees relative spectrum [1.21, 1.21] of the target (sampling LD of 100 haplotypes).
#> Set additive effects for 200 QTL on trait 'LMA' (base: founder_haplotypes [1 filter]; scope: all lines, both parents' copies): exact under the genic anchor; delivered variance 0.4.

# Exact in the founders themselves, linkage disequilibrium included
pop <- pop |>
  get_table("genome_meta") |>
  define_additive_effects("IMF", G = 0.1, anchor = "realised",
    base_tbl = get_table(pop, "ind_meta"), seed = 14)
#> Set additive effects for 200 QTL on trait 'IMF' (base: ind_meta; scope: all lines, both parents' copies): exact under the realised anchor; delivered variance 0.1.
pop |>
  get_table("ind_meta") |>
  extract_genetic_variance(trait_name = "IMF")
#> Realised genetic covariances of the 60 selected individuals (frequencies and regressions are the cohort's own).
#> # A tibble: 3 × 7
#>   effect_name     trait_name_1 trait_name_2 cov_value n_ind decomposition anchor
#>   <chr>           <chr>        <chr>            <dbl> <int> <chr>         <chr> 
#> 1 additive        IMF          IMF              0.100    60 full          reali…
#> 2 between_compon… IMF          IMF              0        60 full          reali…
#> 3 total           IMF          IMF              0.100    60 full          reali…

# Crossbreeding: a common variant centred on the pooled founder table
# (named on purpose, so no Wahlund warning), plus one variant per line,
# each centred on that line's own pool and reading the common target
gm <- pop |> get_table("genome_meta") |> dplyr::filter(chr_name == "2")
pop <- gm |> define_additive_effects("BF", G = 0.2,
  base_tbl = get_table(pop, "founder_haplotypes"), seed = 15)
#> The founder pool (pool expectation under random pairing) sees relative spectrum [0.889, 0.889] of the target (sampling LD of 200 haplotypes).
#> Set additive effects for 100 QTL on trait 'BF' (base: founder_haplotypes; scope: all lines, both parents' copies): exact under the genic anchor; delivered variance 0.2.
pop <- gm |> define_additive_effects("BF", line_name = "Duroc", seed = 16)
#> Set additive effects for 100 QTL on trait 'BF' (base: founder_haplotypes [1 filter]; scope: line Duroc, both parents' copies): exact under the genic anchor; delivered variance 0.2.
#> Line-specific effects are centred on line Duroc's own base, so they add no difference between line means. For differences between lines, use common effects centred on one reference line (see ?define_additive_effects).
pop <- gm |> define_additive_effects("BF", line_name = "Landrace", seed = 17)
#> Set additive effects for 100 QTL on trait 'BF' (base: founder_haplotypes [1 filter]; scope: line Landrace, both parents' copies): exact under the genic anchor; delivered variance 0.2.
#> Line-specific effects are centred on line Landrace's own base, so they add no difference between line means. For differences between lines, use common effects centred on one reference line (see ?define_additive_effects).

# Evaluate the breeding values
pop <- pop |> get_table("ind_meta") |> add_tgv(c("ADG", "BW", "FCR", "BF"))
#> Computed TGV for 60 individuals on trait 'ADG' (additive).
#> Computed TGV for 60 individuals on trait 'BW' (additive).
#> Computed TGV for 60 individuals on trait 'FCR' (additive).
#> Computed TGV for 60 individuals on trait 'BF' (additive).
close_pop(pop)
```
