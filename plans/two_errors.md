# Two outstanding errors from the `define_founder_haplotypes()` audit

Status: **not fixed.** Both were found during the v0.60.x audit of
`define_founder_haplotypes()`. Neither is fixed in that release — the first is a
real correctness bug deliberately deferred because the fix changes
`define_additive_effects()`'s API; the second is currently only documented and
warned about, not structurally addressed.

Everything below was reproduced empirically, not inferred.

---

## Error 1 — `base = "founder_haplotypes"` pools all lines (Wahlund)

**Severity: real bug, silent, affects every multi-line simulation.**

### What happens

`compute_base_allele_freq()` in `R/define_additive_effects.R` computes base
allele frequencies for Falconer centering with this query:

```sql
SELECT gm.locus_id, AVG(CAST(fh.allele AS DOUBLE)) AS f
FROM founder_haplotypes fh
JOIN genome_meta gm ON fh.locus_name = gm.locus_name
GROUP BY gm.locus_id ORDER BY gm.locus_id
```

There is **no `line_name` filter**. With more than one founder pool, `p_base` is
the pool-size-weighted average across all lines.

Verified on a two-line population (line A mosaic, line B gaussian_copula,
200 loci): `cor(p_base, line_B_freq) = 0.68`, while
`cor(p_base, (f_A + f_B)/2) = 1` and `max|p_base − (f_A+f_B)/2| = 1.1e-16`.

### Why it matters

`2p(1−p)` evaluated at a pooled frequency **overstates** within-line
heterozygosity — this is the Wahlund effect. Two divergent lines fixed for
opposite alleles (`p_A = 0`, `p_B = 1`) each have zero within-line variance, but
pool to `p = 0.5` and an apparent `2pq = 0.5`.

`define_additive_effects(scale_to_target = TRUE)` divides by
`sum(2pq·a²)`, so an inflated `2pq` **under-scales** the effects and the realized
additive variance within each line falls short of `target_add_var`. Nothing in
the output signals this — the run looks clean and the numbers look plausible.

### What was done in v0.60.x

- Corrected the roxygen at `R/define_additive_effects.R:29`, which claimed the
  function reads `genome_meta.founder_allele_freq`. It does not — that column is
  write-only (its only production reference was that wrong doc line).
- Added a warning in `compute_base_allele_freq()` when `founder_haplotypes`
  contains more than one distinct `line_name`.

Neither makes the number correct. The warning tells you the result is wrong; it
does not give you a right one.

### Proposed fix

Thread a line into the base-frequency computation:

1. Add a `base_line_name` argument to `define_additive_effects()` (default
   `NULL` = current pooled behavior, so nothing silently changes).
2. When set, add `WHERE fh.line_name = ?` to the query in
   `compute_base_allele_freq()` — parameterized, matching the hardening already
   done in `define_founder_haplotypes()`.
3. Decide the interaction with the existing `line_name` argument (which tags the
   *effects* as line-specific in `genome_effects`). Most likely
   `base_line_name` should default to `line_name` when that is set, since
   line-specific effects almost certainly want line-specific centering — but
   that is a semantic decision worth making explicitly rather than by default.
4. Downgrade the multi-line warning to fire only when neither is set.

**Open question to resolve first:** for a crossbred animal whose alleles trace to
two lines via `ind_haplotype.line_origin`, which base frequency should center its
TBV? `add_tbv()` already joins `genome_effects` on `(locus_name, line_origin)`
and prefers a line-specific effect row over the population-wide one, so the
machinery for per-line centering partly exists — the gap is that
`base_allele_freq` is computed once, pooled, at effect-definition time. Fixing
`compute_base_allele_freq()` without checking this against `add_tbv()`'s join
risks making the two layers disagree.

### Interim workaround

Use `base = "current_pop"` with a line-filtered `base_tbl`:

```r
line_a <- get_table(pop, "ind_meta") |> dplyr::filter(line_name == "A")
pop |> get_table("genome_meta") |> dplyr::filter(...) |>
  define_additive_effects("ADG", base = "current_pop", base_tbl = line_a)
```

### Verification when fixed

Two lines fixed for opposite alleles at a set of QTL. Assert that per-line
`var(TBV)` matches `target_add_var` within sampling tolerance, which it does not
today. Also assert the pooled path still reproduces the old numbers when
`base_line_name` is unset.

---

## Error 2 — `"mosaic"` quantizes the MAF spectrum and fixes monomorphic loci

**Severity: statistical hazard. The code is correct; the model has a property
that is easy to walk into unaware.**

### What happens

In `.gen_haplotypes_mosaic()` (`R/founder_haplotype_helpers.R`), the
`n_templates` template haplotypes are the **only** source of allelic variation.
Every output haplotype is a mosaic of templates, so a locus where all templates
happen to agree is monomorphic in the pool **no matter how many haplotypes are
generated**. With template frequencies drawn `Uniform(0.01, 0.99)`:

> P(locus monomorphic) = E[p^K + (1−p)^K] ≈ 2 / (K + 1)

Measured at 600 loci / 100 haplotypes:

| `n_templates` | monomorphic | predicted `2/(K+1)` |
|---|---|---|
| 2 | 66.5% | 66.7% |
| 5 | 34.5% | 33.3% |
| 10 | 17.2% | 18.2% |
| 15 | 12.3% | 12.5% |
| `gaussian_copula` (any) | 0.5% | — |

The default is `n_templates = max(2, ceiling(sqrt(n_haplotypes)))`, so a typical
100-haplotype mosaic pool has **~17% of loci dead**. The surviving loci have a
MAF spectrum quantized to multiples of `1/K` — there are no genuinely rare
variants, which is often the whole point of simulating a realistic founder pool.

### Why it matters

It compounds with Error 1. QTL placed on monomorphic loci contribute nothing to
`sum(2pq·a²)`, so `scale_to_target = TRUE` inflates the rescale factor and
concentrates the intended genetic variance into the surviving fraction of QTL.
At `n_templates = 2`, two-thirds of QTL are silent.

### What was done in v0.60.x

- Documented both properties in the `"mosaic"` method entry.
- Added a warning when more than 10% of loci come out monomorphic, naming
  `n_templates` and the `2/(n_templates + 1)` relationship.

The model is unchanged — a user who reads the warning and proceeds still gets a
quantized spectrum.

### Options (not yet chosen)

1. **Do nothing more.** Defensible: this is inherent to Li-Stephens copying, and
   `gaussian_copula` already exists for an unquantized spectrum. The docs and
   warning may be sufficient.
2. **Raise the default `n_templates`.** `sqrt(n_haplotypes)` is a weak default —
   it yields K=10 at 100 haplotypes (17% dead). Something like
   `max(20, sqrt(n_haplotypes))` would cut that to ~10%. Cheap, but only shifts
   the problem.
3. **Draw template frequencies from a realistic MAF distribution** rather than
   `Uniform(0.01, 0.99)` — e.g. Beta(0.5, 0.5), matching the `"beta"` method.
   This changes which loci are monomorphic but does not remove quantization.
4. **Guarantee polymorphism structurally** — resample or force a minor allele at
   loci where all templates agree. This breaks the clean Li-Stephens
   interpretation and should probably be opt-in if done at all.

**Recommendation: 2 as a cheap default improvement, and leave 3/4 alone unless a
concrete simulation need appears.** Option 1 is a legitimate outcome.

### Do NOT "fix" the related switch-rate discrepancy

While auditing this, note that `sample.int(n_templates, 1L)` can redraw the
*current* template, so observable template changes occur at
`template_switch_rate × (K−1)/K`, not `template_switch_rate`. Verified at
r = 1, d = 0.01 cM, 2e5 steps: K=2 → 0.5025, K=5 → 0.789, K=15 → 0.9345.

This is **correct and deliberate**. With
`T(d) = e^(−rd)·I + (1 − e^(−rd))·J/K` the kernel composes exactly
(`T(d1)T(d2) = T(d1+d2)`), which is what makes realized LD invariant to marker
density. Excluding self would break that and make LD depend on how finely the
chromosome is marked up. It is documented in v0.60.x; it should not be
"corrected" in code.

---

## Confirmed correct — do not re-audit

The v0.60.x audit verified these are sound, with empirical checks:

- Gaussian-copula thresholding direction and exact `N(0,1)` marginals
  (`P(z < qnorm(p)) = p`; measured 0.3007 at p = 0.3, n = 2e5).
- Per-chromosome AR(1) restart: `rho = 0` zeroes the carry-over exactly
  (within-chromosome adjacent |r| = 0.282 vs 0.003–0.023 across boundaries).
- Step-vs-locus indexing in both LD helpers (`traversal_order` is a permutation;
  `eps`/`U` columns are i.i.d. and exchangeable).
- No locus-id contiguity assumption; disagreeing chr/locus_id orders are caught
  by `resolve_genome_map()`'s monotonicity check.
- `allele_freqs` ↔ `genome_meta` alignment (both sides `ORDER BY locus_id`;
  measured max discrepancy 0 over 600 loci), now also guarded by a length check.
