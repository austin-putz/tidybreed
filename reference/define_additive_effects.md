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

  Character scalar **or** vector. Name(s) of existing traits in
  `trait_meta`. When length \>= 2, the architecture is drawn jointly
  from `MVN(0, G)` and `method` becomes active.

- distribution:

  Character. `"normal"` (default) or `"gamma"`, the single-trait
  architecture. Ignored for multi-trait (always MVN).

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

  Character. `"shared"` (default) or `"union"`. Multi-trait only.
  `"shared"` — all listed traits use the filtered loci as their shared
  QTL set; exact. `"union"` — per-trait QTL sets are read from existing
  generated terms at this scope, restricted to the filtered loci;
  per-trait scaling, approximate for a non-zero covariance.

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
  and has consumed RNG draws.

## Value

The modified `tidybreed_pop` (invisibly).

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

**How.** Effects are drawn as today (one draw per QTL for one trait;
joint `MVN(0, G)` rows for several), and that draw is only the
*architecture* `B0`. It is then right-multiplied by a `k x k` matrix `A`
so that `B = B0 A` satisfies `B' M B = G` exactly (the congruence of
Proposition 2 in the source method). For one trait this is exactly the
scalar rescale `b0 * sqrt(G / sum(w b0^2))`. For two or more it also
fixes the genetic **correlations**, which a per-trait rescale cannot: at
200 QTL and a target correlation of 0.4 a scalar rescale delivers
anything from about 0.18 to 0.60. The same seed gives the same `B0`.

**The anchor `M`** is the reference-population genotype covariance the
calibration is exact for:

- `anchor = "genic"` (default): `M = diag(n_eligible * p * q)` at the
  base allele frequencies, the random-mating (HWE + linkage-equilibrium)
  limit. Use it for multi-generation studies: under random mating the
  realised covariance converges to it. The manuscript's "reference"
  anchor is this plus a non-default `base_tbl`, whose frequencies the
  weights use.

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
population sees (§7.4 of the plan): the **pool expectation** `2 Cov(H)`
when the base is the founder pool, the **observed** `Cov(X)` when it
selects individuals, and the genic limit under `anchor = "realised"`.
The founder pool's comparison is always a message: a small pool's
departure is its sampling LD, not a mistake in the call. To get the
target exactly in the founders, add them first and calibrate with
`anchor = "realised"` on them. The observed and genic-limit comparisons
warn when the relative spectrum leaves `warn_bounds`. Nothing is stored.

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

[`define_trait()`](https://austin-putz.github.io/tidybreed/reference/define_trait.md),
[`define_effect_cov_matrix()`](https://austin-putz.github.io/tidybreed/reference/define_effect_cov_matrix.md),
[`define_genome_effect_terms()`](https://austin-putz.github.io/tidybreed/reference/define_genome_effect_terms.md),
[`extract_allele_freq()`](https://austin-putz.github.io/tidybreed/reference/extract_allele_freq.md).

## Examples

``` r
if (FALSE) { # \dontrun{
# Single trait — all loci on chr 1-5 become QTL; write the target and
# calibrate to it in one call
pop <- pop |> define_trait("ADG")
pop <- pop |>
  get_table("genome_meta") |>
  dplyr::filter(chr_name %in% as.character(1:5)) |>
  define_additive_effects("ADG", G = 0.25)

# Correlated traits — shared QTL set, exact G (variances and correlation)
pop <- pop |> define_trait("BW")
G <- matrix(c(0.25, 0.10, 0.10, 0.30), 2, 2,
            dimnames = list(c("ADG", "BW"), c("ADG", "BW")))
pop <- pop |>
  get_table("genome_meta") |>
  define_additive_effects(c("ADG", "BW"), G = G)

# A target stored beforehand is read back: no G
pop <- pop |> define_trait("FCR") |>
  define_effect_cov_matrix("additive", 0.1, trait_name = "FCR")
pop <- pop |> get_table("genome_meta") |> define_additive_effects("FCR")

# Exact in the generation-0 individuals themselves, LD included
pop <- pop |>
  get_table("genome_meta") |>
  define_additive_effects("FCR", anchor = "realised",
    base_tbl = get_table(pop, "ind_meta") |> dplyr::filter(gen == 0L))

# Crossbreeding: three variants. The common fallback names its base to say
# "yes, pool"; each line's variant centers on its own founder pool by default.
gm <- pop |> get_table("genome_meta") |> dplyr::filter(chr_name == "1")
pop <- gm |> define_additive_effects("ADG",
               base_tbl = get_table(pop, "founder_haplotypes"))
pop <- gm |> define_additive_effects("ADG", line_name = "Duroc")
pop <- gm |> define_additive_effects("ADG", line_name = "Landrace")
} # }
```
