# Import QTL-effect methods — Step 4 results

Plan: [import_qtl_effect_methods_phase_4_plan.md](import_qtl_effect_methods_phase_4_plan.md),
revised after the [Codex review](import_qtl_effect_methods_phase_4_codex_review.md).
Run in one pass (both sub-steps, no review pause between them, decided by the user
2026-10-05); each sub-step is its own commit with the full suite green.

## 4a — Builders and the NOIA conversion (0.75.0)

**Built**
- `aa_terms()` (exported), in `R/genome_effect_terms_builders.R`, as 4a.2 specifies.
  Inputs are validated before `e = 0` pairs are dropped (D7); non-finite coefficients and
  frequencies are refused (also now in `ad_terms()`).
- One fixed column set for `ad_terms()`, `aa_terms()`, `genotype_terms()` through the internal
  `.terms_frame()`; `.ad_rbind_fill()` is deleted.
- Collision-free `term_id`s through one encoder, `.term_id(builder, loci, suffix)`:
  `"ad:2:L1#a"`, `"aa:1:A|3:BxC"`, `"geno:1:A|1:B#4"` (D5 as revised).
- `.stored_to_functional()`, `.noia_to_stored()` and `.noia_terms()` (the writer-frame step
  4a.3 left implicit), with the two sign conventions stated in their headers.
- The writer needed no change: the typed frame with `NA` in irrelevant columns writes as is.

**Gates** — `tests/testthat/test-genome-effect-terms-builders.R`, 14 tests, all pass.
- C16: identical names and types for every builder and coding, `rbind()` in every order;
  the collision gates (`c("A","B")` vs `"AxB"`, `(A, BxC)` vs `(AxB, C)`, names containing
  `:`, `|`, `#`, one locus in every builder) assert the written term count, member orders
  and that the evaluated total equals the sum of the separately written parts; two surfaces
  over the same loci are refused by the writer.
- N1: each conversion row in its own trait, evaluated by `add_tgv()` on a 200-individual
  panel where all single-locus and all nine pair states occur, at centres ≠ 0.5; plus one
  two-owner model with every shape. 1e-12.
- N2: random `(a, d, e, p)` → `.noia_terms()` written → read back through
  `.gev_read_model()` → `.stored_to_functional()` returns `(a, d, e)` and `kappa = −mu`
  (1e-12); the functional-coding model differs from it by `mu` alone (1e-10); `mu` equals
  the §3 formula computed in the test; `ad_terms()`' reported `mu` equals the no-pair `mu`.
- N3: Cockerham `aa_terms()` writes centres `p_k`, `p_l` and evaluates
  `e (g_k − 2p_k)(g_l − 2p_l)`.

**Mutation check** — flipping the dominance row's `a −= v(1 − 2c)` fails N1 (both tests)
and N2.

**Suite** — full suite with `NOT_CRAN=true`: 0 failures, 0 errors, 4,080 expectations
passed, 6 warnings (the 0.74.5 count; the same five files as before: `add_founders`,
`add_phenotype`, `genome_map`, `parity`, `phenotype_composite`). `devtools::document()` and
`pkgdown::check_pkgdown()` clean.

**Census of `term_id` assertions (risk in the plan).** No existing test compared builder
`term_id` values; the only builder assertions were counts of distinct ids, which are
unchanged.

## 4b — `extract_genetic_variance()` (0.75.1)

**Built** — `R/extract_genetic_variance.R`, exported, read-only, as 4b.1–4b.5 of the revised
plan specify:
- trait default restricted to traits with terms (`.egv_traits()`; the shared resolver is
  unchanged), whole-genotype default base frequencies through a registered id view, an
  explicit `base_tbl` through `extract_allele_freq()`, missing base frequencies an error;
- classification per family (`.egv_classify()`; owners summed only after it), the three
  cases, block availability on the canonical coefficients, the induced `e·c` of a fixed
  pair member kept;
- the source's realised algebra without its dense anchors; A×A accumulated in
  deterministic pair chunks (`.egv_aa_values()`, `QTL_REALISED_MAX_CELLS / n` pairs per
  chunk); the dosage guard before any evaluation; `between_components` computed directly
  from the value matrices; the internal sum-to-total assertion per trait (mixed tolerance
  `1e-10 + 1e-8 · max|g − ḡ|`);
- the `anchor` column and the population message (D1).
- The dosage collector is the shared `.collect_dosages()` + `.dosage_guard()` in
  `R/genome_effects_helpers.R`; `.dae_collect_dosages()` is a thin caller with its old
  messages (A21 still passes). Its row names now follow the SQL order instead of R's
  `sort()`, which could have mislabelled rows under a collation that orders ids differently
  from DuckDB.

**Gates** — `tests/testthat/test-extract_genetic_variance.R`, 24 tests (B1–B18), all pass.
The oracle is `tests/testthat/helper-nonadd-oracle.R`, a verbatim copy of the source's
`nonadd_covariates()` / `nonadd_decompose()` (unchanged between `8f8a97c` and the source's
current `318e54f`). Fixture notes:
- B2 writes the model as Cockerham terms through `.noia_terms(.noia_to_stored(...))` and
  compares every block, `total` and `between_components` (including Cov(D, A×A), both
  orientations) with the oracle at 1e-10 on a 12-haplotype LD panel.
- B3 uses a 6-haplotype panel with an asymmetric two-trait model, asserting the blocks
  alone miss `total` by more than 1e-6; a second test adds a multi-locus indicator surface
  so `unpartitioned` cross terms are part of `between_components`.
- B4/B6/B12/B13 build F1s of two lines with `add_offspring()`; the scoped variants are
  custom-owner terms (`with_additive_terms()`), which exercise the same scope machinery as
  generated ones without needing per-line targets.
- B12 fixes a locus by rewriting its `ind_haplotype` alleles (dosage 0, 1, 2) and compares
  with the reduced model carrying `e·c`.
- B17's cross-chunk comparison is `expect_equal(1e-12)`, not identical: different chunk
  sizes group the floating-point sums differently. A given chunk size is bit-identical on
  repeat, and the chunk size depends only on `n`, so the output is a function of the inputs.
- B18's repeated records come from `add_phenotype()` twice on a repeatable phenotype.

**Mutation checks** (scripted, file restored after each):

| Mutation | Fails |
|---|---|
| Classify per term, not per family | B13 |
| Omit the D–A×A cross covariance | B2, B3 |
| Only one orientation of each cross pair | B2, B3 (both), B13 |
| Omit `unpartitioned` cross terms | B3 (surface test), B13 |
| Drop the induced `e·c` at monomorphic loci | B2, B3, B9, B12, B17 |
| Availability from stored contrast names | B16 |

(The plan's "total minus blocks" mutation is not used: it is algebraically equal, Codex
finding 8.)

**Benchmark** — `dev/benchmarks/benchmark_extract_genetic_variance.R`: 2,000 individuals,
500 loci (5 chromosomes, 400 founder haplotypes). Two runs on the same machine; the second
ran right after the full suite, so its times are higher.

| Case | Run 1 | Run 2 | Peak R heap (run 2) |
|---|---:|---:|---:|
| typical (2 traits, A + D + 1,000 pairs), realised | 8.2 s | 12.9 s | 304 MB |
| typical, genic | 1.1 s | 1.8 s | 227 MB |
| all pairs of 200 loci (19,900), realised | 65.0 s | 113.2 s | 849 MB |
| all pairs of 200 loci, genic | 13.1 s | 13.1 s | 587 MB |

- A literal port of `nonadd_covariates()` at the all-pairs size would allocate a
  2,000 × 19,900 pair matrix (~318 MB) and a 19,900² pair covariance (~3.2 GB). The extractor
  ran it in two chunks of 10,000 pairs; its peak is one chunk's matrix plus temporaries,
  bounded by the cell budget whatever the pair count.
- The realised time is dominated by the evaluator (`.gev_evaluate()`, the same cost as
  `add_tgv()`): the genic run, which evaluates nothing, takes 13 s at this size, spent in the
  R-side classification and conversion loops over 20,000 terms.
- **Not the planned all 124,750 pairs of 500 loci.** Writing that model through
  `define_genome_effect_terms()` takes most of an hour: the writer's per-term R loop in
  `.ge_build()` scales slightly worse than linearly (1,000 / 2,000 / 4,000 pairs: 2.2 / 5.0 /
  11.5 s; the 19,900-pair write took 109–196 s). **Risk for step 5**: the
  `define_genome_effects()` generator writes through the same engine, so a model with tens of
  thousands of A×A pairs will spend minutes in the writer. Worth a writer benchmark before
  step 5 relies on large pair sets.

**Suite** — full suite with `NOT_CRAN=true`: 0 failures, 0 errors, 4,223 expectations
passed (4,080 at 0.75.0), 6 warnings (the same five files as 0.74.5). A21
(`define_additive_effects()` realised size guard, through the shared collector) passes.
`devtools::document()` and `pkgdown::check_pkgdown()` clean.

## Implementation-review follow-up (0.75.2)

[Codex reviewed the built phase](import_qtl_effect_methods_phase_4_codex_review.md) and
found four issues. All four were accepted and fixed:

| # | Finding | Change | Gate |
|---|---|---|---|
| 1 | Under LD the realised `additive` row is not the cohort's joint least-squares additive projection | Documented (roxygen section "What the realised blocks are, under LD", API skill, main plan §4.3.3 and §8): the rows are NOIA contrast components; `full` means supported shapes. The estimator is unchanged — it is the source's | Codex's 32-individual HWE-margin LD panel: `additive = 0`, `between_components = 0`, while `lm()` explains 0.021 |
| 2 | Floating-point residue in an accumulated `d` (or summed `e`) decided which rows exist | `.stored_to_functional()` sets a coefficient to 0 when `|x| ≤ n·eps·Σ|contributions|` (`.cancelled()`); a single contribution is never residue | B16: a linear surface in all six row orders, and `ad_terms(d = 0)`, have no dominance row under both anchors; `d = 1e-9` and `d = 1e-20` keep theirs; three owners summing a pair to 0 have no A×A row |
| 3 | `.noia_to_stored()` gave `NA` alpha for pair-only loci and `.noia_terms()` then silently wrote only the pair | The converter is sparse: coefficients span every locus the model names, missing main effects 0. `.noia_terms()` refuses misaligned or incomplete input | Pair-only at `p = (0.3, 0.6)` (alpha `0.2, −0.4`, `mu = −0.08`), one endpoint missing, permuted/disjoint `a`/`d` names: functional − converted = `mu` through the evaluator; refusal cases |
| 4 | `genotype_terms()` truncated `copy_count = 2.9` to 2 (and `−0.5` to 0) | Refused before `as.integer()` and before zero rows drop | `2.9`, `−0.5`, a fractional count on a zero-valued row, `NA` refused; 0, 1, 2 accepted |

Every new gate fails on the 0.75.1 code (checked by restoring the old builder file) and
passes on 0.75.2. Also fixed: a roxygen inline-code false positive (`` `r x k` `` in an
internal comment) that made `document()` print an error.
