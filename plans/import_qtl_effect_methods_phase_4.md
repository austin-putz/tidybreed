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
