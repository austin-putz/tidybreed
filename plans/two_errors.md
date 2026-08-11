# Two outstanding errors from the `define_founder_haplotypes()` audit

Status: **both resolved.** Found during the v0.60.x audit of
`define_founder_haplotypes()`, fixed in v0.61.0. The measurements and reasoning
below are kept because they justify the design; see *Resolution* under each.

Everything below was reproduced empirically, not inferred.

## Resolution summary

- **Error 1 — fixed.** `define_additive_effects()` gained `base_line_name`,
  which **defaults to `line_name`**, so a line-specific effect is centered on
  its own founder pool and a population-wide effect on the pooled base. The
  open question below resolved *in favor* of the fix: `add_tbv()`'s centering is
  already per-haplotype-row, so per-line base frequencies make the two layers
  agree rather than diverge.
- **Error 2 — addressed by the chosen option (raise the default).** The mosaic
  `n_templates` default went from `max(2, sqrt(n_hap))` to
  `max(20, sqrt(n_hap))`, clamped to `n_haplotypes`. Measured monomorphic
  fraction at 100 haplotypes: **17% → 7.2%**. The Li-Stephens model itself is
  unchanged, so quantization remains — just finer.

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

**Open question — RESOLVED, in favor of the fix.** The worry was that per-line
base frequencies might disagree with `add_tbv()`'s per-`line_origin` join. Reading
the actual SQL (`R/add_tbv.R:205-226`) settles it: `base_allele_freq` is selected
from the **same joined row** that supplies `genome_value`, and the `NOT EXISTS`
fallback is correlated on both `h.locus_name` *and* `h.line_origin`. So precedence
is per haplotype row, and the consumption layer was **already** fully per-line —
only the producer was pooled. An F1's line-A alleles are centered on line A's
stored frequency and its line-B alleles on line B's, automatically. Existing test
`tests/testthat/test-add_tbv.R:184` already asserted this.

Residual asymmetry, inherent to the fallback rather than to the fix: where a locus
has a line-specific row for A but only a population-wide row for B, the B allele
is centered on the population-wide frequency.

### Resolution (v0.61.0)

`base_line_name` added to `define_additive_effects()`, defaulting to `line_name`.
`compute_base_allele_freq()` gained a `line_name` parameter and applies
`WHERE fh.line_name = ?` as a **bound** parameter, with the name checked by the
existing `validate_sql_identifier()`. A named line with no `founder_haplotypes`
rows is a hard error, not a silent zero — `out` is zero-initialised, so a typo
would otherwise center every allele at 0 and contribute nothing to `V_A` while
still passing a `[0,1]` range assertion. The multi-line warning now fires only for
a population-wide effect, where the ambiguity is real.

The former workaround (`base = "current_pop"` with a line-filtered `base_tbl`)
still works and remains right when you want the *realized* founder-sample
frequency rather than the founder-pool frequency.

### Verification (implemented)

`tests/testthat/test-define_additive_effects.R` — two lines fixed for opposite
alleles: line A stores `p = 0`, line B `p = 1`, a population-wide effect `p = 0.5`.
A second test with A at `p ≈ 0.1` and B at `p ≈ 0.9` shows the payoff: with
per-line centering, within-line realized variance hits `target_add_var = 2`
exactly; with `base_line_name = NULL` forcing the old pooled behaviour, the
Falconer bookkeeping still "hits" 2 while the variance actually realized within
line A is **below 1.0**. Plus: unknown line errors, `base_line_name` with
`base = "current_pop"` errors, an injection payload is rejected, and a new test
covers the previously untested guarantee that per-line calls don't clobber each
other's rows.

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

### Resolution (v0.61.0) — option 2 chosen

Default is now `max(20, ceiling(sqrt(n_haplotypes)))`, clamped to `n_haplotypes`
so a small pool cannot trip the function's own `n_templates <= n_haplotypes`
check with a default the user never chose. Measured monomorphic fraction at
600 loci:

| `n_haplotypes` | K before → after | monomorphic after |
|---|---|---|
| 5 | 3 → 5 (clamped) | 47.5% |
| 30 | 6 → 20 | 9.5% |
| 100 | 10 → 20 | 7.2% |
| 400 | 20 → 20 | 7.2% |
| 1000 | 32 → 32 (unchanged) | 4.3% |

Above ~400 haplotypes `sqrt` dominates and nothing changes. Options 3 and 4 were
**not** taken: quantization to multiples of `1/K` remains, and `gaussian_copula`
is still the method to reach for when an unquantized MAF spectrum matters.
**This changes seeded `"mosaic"` output.**

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
