# Import QTL-effect methods — Step 5 results

Plan: [import_qtl_effect_methods_phase_5_plan.md](import_qtl_effect_methods_phase_5_plan.md),
revised twice after the [Codex review](import_qtl_effect_methods_phase_5_codex_review.md)
and its re-review. Decisions D1–D8 were made by the user 2026-10-06 to 2026-10-09, each as
recommended; D5 adds a rank `message()` for singular targets (built in 5b). Three commits
(D1); the user reviews between them.

## Step 5 at a glance — complete (0.75.3 → 0.76.3)

| Commit | Version | What |
|---|---|---|
| 5a | 0.75.3 | The genome-effect writer and evaluator preparation made linear in model size (all 124,750 pairs of 500 loci written in 7.1 s; 16,000 pairs went from 71.8 s to 0.7 s). Internal only. |
| 5a review | 0.75.4 | Exact pair keys above 2^53. |
| 5b | 0.76.0 | `define_genome_effects()`: additive, dominance and A×A effects sampled and calibrated to `G_A`, `G_D`, `G_AA` exactly under `"genic"` or `"realised"`; the additive floor; inbreeding depression (exact for one trait); additive-only targets take `define_additive_effects()`'s path, with identical rows. |
| 5b review | 0.76.1 | The stored additive coefficients are the verified ones (they lost digits when `G_A` ≪ `G_D`); realised designs only for non-zero blocks, so the size guard is exact; the extractor's precision limit documented. Re-verified by Codex. |
| 5c | 0.76.2 | End-to-end gates (phenotypes, prevalence, removal, both anchors measured back), the generator benchmark, the "Genetic models" vignette. |
| 5c review | 0.76.3 | `?remove_generated_effects` qualified (kept targets still count for generation; an empty model cannot be re-evaluated, so old `ind_tgv` values stay); "No genome effects found" names `define_genome_effects()`; a test pins it. |

Every Part C gate of the main plan (§11 C1–C20) and the phase plan's G1–G8 is in the suite;
the full suite after 5c is 1,170 tests, 4,855 expectations, no failures. Open after step
5 (not defects): `extract_genetic_variance()` loses relative precision on a block many
orders of magnitude below its coupling (documented, decided with the user), and the
out-of-scope list of the phase plan stands.

## 5a — Writer speed (0.75.3)

Internal only: no API, schema, stored-value or message change.

### Order of work

1. **Contract test first.** `tests/testthat/test-genome-effects-writer-order.R` (6 tests,
   9 expectations) was written against the unmodified 0.75.2 code and passed there before
   any refactor. It is unchanged through 5a and passes after it. It pins:
   - two terms broken by the same rule: the first **typed** (not the first sorted) is named;
   - a later rule failing on an earlier term beats an earlier rule on a later term (repeated
     locus on term 1 vs two coefficients on term 2; `effect_name` on term 1 vs coefficient
     on term 2), and rule order inside one term;
   - the whole `Invalid 'terms':` list (term order, then `locus_id` inside a term, then the
     two per-member rules) and the whole `Invalid 'origin':` list;
   - `.ge_validate_frames()`'s whole violation vector over every rule at once: orphan
     members and origins first, then per term in **row** order (not id order), then the
     family clashes; labelled and unlabelled ids;
   - a whole-table conflict with a stored row names the user's `term_id` and rolls back.
2. **Baseline** with the new `dev/benchmarks/benchmark_genome_effect_writer.R` on 0.75.2.
3. **Vectorisation** (below), then the old-vs-new equivalence check, the tests, the
   benchmark again, a profile, and the full suite.

### Changes

Every change replaces a per-term or per-id subset of a data frame (`x[x$id == id, ]` inside
a loop, `O(terms × rows)`) with one grouping pass. Messages are built for all rows at once
and put back in the order the loop produced them.

| Function | File | Change |
|---|---|---|
| `.ge_build()` | `R/define_genome_effect_terms.R` | Per-term scalar rules decided for every term at once (`.ge_any_by()`); the error is the first term with any violation, then its first rule. Members built by one `order(term, locus_id)` and `sequence()` for `member_slot`. |
| `.ge_infer_copy_counts()` | same | Vectorised lookup; the error names the first member (in member order) on a chromosome with more than one copy count. |
| `.ge_check_member_fields()`, `.ge_check_origin_fields()` | same | Each rule a vector over rows; messages interleaved row by row, as the loop did. |
| `.ge_resolve_deletes()` (`replace_scope`) | same | `.ge_terms_at_scope()`: members and origins split once (`.ge_rows_by_id()`). |
| `.ge_validate_frames()` | `R/genome_effects_helpers.R` | The per-term phase is `.ge_frame_term_messages()`: every rule vectorised, messages ordered by (term row, rule, member row, member rule). Family predicates read pre-split rows. |
| `.ge_validate_dominance_ploidy()` | same | Origins are read only for members off a diploid chromosome. |
| `.stored_to_functional()` | `R/genome_effect_terms_builders.R` | Contributions per shape gathered as vectors and summed with `rowsum()` in term order (`.s2f_accumulate()`). `rowsum()` adds in input order in double precision, so `a`, `d`, the cancellation bound and `kappa` are **bit-identical** to the loop. |
| `.gev_target_kind()`, `.gev_term_line()`, `.gev_term_parent()` | `R/genome_effects_eval.R` | One grouping pass (`.gev_term_origin_values()`). Not in the plan's list; see "Deviations". |
| `.dae_existing_loci()`, `.dae_target_dependents()`, `.dae_warn_parent_only()` | `R/define_additive_effects.R` | Pre-split rows (`.ge_rows_by_id()`), computed once per call. Not in the plan's list. |
| `.gev_variant_map()`, `.gev_preflight()` | `R/genome_effects_eval.R` | Families and their first term reached by position, not by name. Not in the plan's list; found by the extractor benchmark. |

Helpers added: `.ge_any_by()`, `.ge_pair_key()` (exact numeric key over two integer
columns, with a pasted-string fallback for `NA`), `.ge_rows_by_id()`, `.ge_terms_at_scope()`,
`.ge_frame_term_messages()`, `.gev_term_origin_values()`, `.s2f_accumulate()`.

`validate_genome_effects()` still validates the **whole** stored table before every
`COMMIT` (CLAUDE.md); it was not narrowed to the new rows. `.ge_commit()` and the
transaction are unchanged.

### Equivalence (development check, not committed)

A scratch script sourced the 0.75.2 versions of the five changed files into their own
environment and compared old and new with `identical()` (CLAUDE.md forbids committed
golden-from-old tests):

- `.ge_build()` on 60 random term sets: additive, dominance, indicator surfaces with and
  without inferred copy counts, pairs and triples, rows interleaved across terms, character
  and numeric `term_id`s, origins (scalar list), centres filled from a base;
- 18 malformed `terms` (each rule, several rules at once, NA/NaN coefficients, `effect_name`
  with NA, X-chromosome copy-count inference, unknown columns/contrasts/loci, fractional
  dosage) and two malformed `origin` frames: identical error messages;
- `.ge_validate_frames()` on 200 random frames with corruptions (bad slots, missing
  members, unsorted loci, orphans, duplicated terms causing family clashes, random origin
  rows, shuffled rows, partial labels);
- `.stored_to_functional()` on 80 random covered models, half with a second owner whose
  coefficients exactly cancel, plus large locus ids and both internal errors;
- `.gev_target_kind()` / `.gev_term_line()` / `.gev_term_parent()`, `.ge_resolve_deletes()`
  for every mode × scope × owner, and `validate_genome_effects()` on a stored multi-owner,
  line- and parent-scoped model;
- locus names with spaces, `:`, `|`, `#`, a literal `"NA"`, non-ASCII and all digits,
  through `.ge_build()` (valid and failing), `validate_genome_effects()` and
  `.stored_to_functional()`;
- the new `.ge_build()` against itself (determinism on repeat).

**Result: 455 comparisons identical (357 of them non-trivial: errors, violation vectors or
result lists), 0 mismatches**; 462 (360 non-trivial) with the evaluator and unusual-name comparisons
added later (see "Deviations"). `.stored_to_functional()` agrees bit for bit, not just within
the plan's 1e-12.

### Benchmark

`dev/benchmarks/benchmark_genome_effect_writer.R`; R 4.5.3, macOS (Darwin 24.6.0), x86_64,
16 cores, DuckDB 1.5.5 with 16 threads; one trait; 500 loci on 5 chromosomes.

| A×A pairs | 0.75.2 write | 0.75.3 write | 0.75.3 setup (`aa_terms()`) | 0.75.3 per pair | 0.75.3 peak R memory |
|---:|---:|---:|---:|---:|---:|
| 1,000 | 2.2 s | 0.24 s | 0.03 s | 237 µs | 165 Mb |
| 4,000 | 9.6 s | 0.46 s | 0.15 s | 115 µs | 179 Mb |
| 16,000 | 71.8 s | 0.70 s | 0.42 s | 44 µs | 194 Mb |
| 64,000 | not run | 3.5 s | 1.6 s | 54 µs | 243 Mb |
| 124,750 (all pairs) | not run | 7.1 s | 3.4 s | 57 µs | 390 Mb |

- **Targets met:** time per pair from 4,000 to 64,000 pairs grows **0.47×** (target at most
  2×; the fixed per-call cost dominates small writes), and all 124,750 pairs take **7.1 s**
  (target under 60 s).
- **Replacement into a populated model** (16,000 pairs of one trait replaced, `replace_owner`,
  in a database of 20,200 terms: three traits, three owners, 100 line-scoped variants):
  **74.6 s → 1.4 s**, peak 280 Mb.
- **`.stored_to_functional()`** on 16,000 terms: **8.56 s → 0.17 s**.
- The baseline's peak-memory column read the wrong `gc()` column (R 4.5 adds `limit (Mb)`), so
  0.75.2 memory is not recorded; the script now finds the column by name.

**Profile after** (124,750 pairs, 7.2 s): `.ge_commit()` 91%, of which `.ge_validate_frames()`
71% (it runs twice: on the candidate frames and, inside `validate_genome_effects()`, on the
whole table) and `.ge_family_keys()` 33% (string keys, `tapply`). Everything left is linear
(string pasting, `factor()`, `order()`, DuckDB I/O). Not optimised further: the targets are met
by an order of magnitude, and the family key must keep agreeing with the SQL view's
`family_key`.

**Step-4 extractor benchmark at the planned size** (`benchmark_extract_genetic_variance.R`,
now all 124,750 pairs of 500 loci; 2,000 individuals):

| Case | Step 4 (0.75.1) | 0.75.3 | Peak R heap (0.75.3) |
|---|---:|---:|---:|
| writing the all-pairs model | 109–196 s (19,900 pairs) | 10.6 s (124,750 pairs) | — |
| typical (2 traits, A + D + 1,000 pairs), realised | 8.2 s | 6.4 s | 319 MB |
| typical, genic | 1.1 s | 0.38 s | 228 MB |
| all pairs, realised | 65.0 s (19,900) | 295 s (124,750) | 1,057 MB |
| all pairs, genic | 13.1 s (19,900) | 7.2 s (124,750) | 528 MB |

- The genic run evaluates nothing. Its 13 s at 19,900 pairs was the per-term R loops
  (`.stored_to_functional()` and the target lookups); with those now one pass it takes 7.2 s
  at 6× the pairs.
- **The realised run was not linear in the first attempt.** It ran 14 minutes before failing
  (below). Profiling realised extraction at 5,000 and 20,000 pairs showed it 10× slower for
  4× the pairs: the evaluator's R-side preprocessing, `.gev_variant_map()` (0.9 → 10.3 s)
  and `.gev_preflight()` (0.2 → 2.3 s), looked up per-term lists **by name** inside a
  per-family loop, which is quadratic. Both now index by position (see "Deviations"). After
  the change: 5,000 pairs 10.4 s, 20,000 pairs **25.5 s (was 124 s)**, 124,750 pairs 295 s
  (2.4 ms per pair, 3.3 ms at step 4's 19,900).
- What remains in the realised run is DuckDB's one evaluation statement (`.gev_sql()`: the
  same cost as `add_tgv()`), whose CPU time is linear in the pair count. At 124,750 pairs ×
  2,000 individuals it spills to disk. **Risk for 5c:** `add_tgv()` / `add_phenotype()` on a
  generated model with tens of thousands of pairs costs minutes per evaluation; the 5c
  fixtures should stay at hundreds of pairs, and a faster evaluator is its own change.
- **In-memory spill failure (environment finding, not a 5a change).** An in-memory
  (`":memory:"`) database spills to `tempdir()/duckdb/temp`, and DuckDB does not create the
  missing `duckdb/` parent: the query fails with `IO Error: Failed to create directory`.
  File-backed databases (tidybreed's default) spill next to their file and are not
  affected. The benchmark creates the directory first. If in-memory populations are meant to
  handle spilling queries, `open_pop(db_name = ":memory:")` could create it; not done here.

### Deviations from the plan

- **More functions vectorised than listed.** 5a.2 named `.ge_build()`, `.ge_validate_frames()`,
  `validate_genome_effects()`, `.ge_resolve_deletes()` and `.stored_to_functional()`. The same
  per-id pattern was in `.gev_target_kind()` / `.gev_term_line()` / `.gev_term_parent()`
  (target refusals, the prevalence threshold), `.dae_existing_loci()`,
  `.dae_target_dependents()` and `.dae_warn_parent_only()`; 5b's generator reaches all of
  them with large models, so they were fixed here under the same equivalence check.
- **The evaluator's R-side preprocessing** (`.gev_variant_map()`, `.gev_preflight()` in
  `R/genome_effects_eval.R`) was quadratic in the term count for the same reason (a
  lookup by name per family). Found by the extractor benchmark the plan asks for; fixed by
  indexing families and their first term by position. The evaluation SQL is untouched. The
  equivalence script compares the whole old and new `.gev_evaluate()` on a model with
  common, line-scoped and parent-scoped variants, dominance and pairs (identical, bit for
  bit, after sorting rows: the statement has no `ORDER BY`), and `.gev_preflight()` with and
  without its warning firing. Final count, with the unusual-name cases: **462 comparisons
  identical (360 non-trivial), 0 mismatches.**
- **The contract test was committed with the refactor, not before it.** D1 fixes three
  commits, so the test lands in the 5a commit; it was run and passed against the unmodified
  0.75.2 code first (step 1 above), which is what the plan's "committed before" protects.
- **A pre-existing roxygen slip fixed:** an `@noRd` comment in
  `R/extract_genetic_variance.R` held the inline code `` `r x 2` ``, which roxygen tried to
  evaluate as R (`devtools::document()` printed a parse failure on 0.75.2 too). Reworded.
- **`validate_genome_effects()` itself needed no separate change**: after the frame
  validator, its remaining cost is linear (profile above).

### Suite

Full suite with `NOT_CRAN=true`: **0 failures, 0 errors, 4,295 expectations passed**, 6
warnings (the same five files as 0.74.5: `add_founders`, `add_phenotype`, `genome_map`,
`parity`, `phenotype_composite`). No existing test file was edited. `devtools::document()`
(no `man/` or `NAMESPACE` change) and `pkgdown::check_pkgdown()` clean.

Note on the run: a `git stash` / `stash pop` used to confirm the roxygen slip predates 5a ran
for under a minute while the suite was running. The suite loads the package once at start,
the new test file is untracked (untouched by the stash), and the only test file that reads
files from disk (`test-schema.R`) was rerun on its own afterwards: 34 passed.

### Implementation-review follow-up (0.75.4)

[Codex reviewed the built 5a](import_qtl_effect_methods_phase_5_codex_review.md) and
recommended proceeding with 5b. It reran both writer gates (0.48×, 7.2 s), the targeted test
files, and its own 7,020 old-vs-new comparisons (all identical). One finding, accepted:

- **`.ge_pair_key()` precision (low).** The key `x · (max(y) + 1) + y` was not checked for
  exactness: at radix 2^31 and `x = 2^22`, the pairs `(2^22, 3)` and `(2^22, 4)` rounded to one
  double, which would report a false repeated locus. The header's `max(x) · max(y)` bound was
  also wrong (it left out the radix's `+ 1` and the final `+ y`). Fixed: the numeric key is
  used only when `(max(x) + 1) · (max(y) + 1) ≤ 2^53`, which bounds every key, with the same
  non-negative-whole check on `x` as on `y`; otherwise string keys. A hand-derived test in
  `test-genome-effects-writer-order.R` checks the boundary pair, an ordinary small-key case,
  and equal pairs and `NA` on both paths; it fails on the 0.75.3 helper. The equivalence
  script is unchanged (462 identical); no realistic model reaches the boundary.

Two qualifications from the review are added to the profile above, which said "everything
left is linear":

- the evidence is for growing numbers of **common** A×A families. Family validation still
  compares every pair of variants **within** a family, and the scoped paths keep some
  per-term work (`.ge_term_at_scope()`, the cache misses of `.gev_variant_map()`). A model
  with very many scoped variants of one term, or of very high order, is not shown linear;
- large realised evaluation stays expensive (the DuckDB statement, above); 5b/5c fixtures
  stay small, and the forward `.noia_to_stored()` is profiled in 5b as planned.

## 5b — `define_genome_effects()` (0.76.0)

The generator, its calibration internals, the owner rules elsewhere, and gates C1–C12,
C14, C15, C17, C18, C20 and G1–G8. Everything the plan lists for 5b is built; the
deviations are below.

### What was built

| Piece | File | Notes |
|---|---|---|
| `define_genome_effects()` | `R/define_genome_effects.R` (new) | Exported. Validation in the plan's 5b.2 order (no RNG, no write), the D2 draw (`.dge_draw()`), the two routes (D5), storage, one `.ge_commit()` with the passed targets in `before_commit`, the eight messages of 5b.7, the D6 diagnostics (`.dge_diagnostics()`). |
| Calibration | `R/genome_effects_calibration.R` (new) | `.na_anchors()` / `.na_aa_anchor()` (anchor objects, never dense `M`), `.na_coupling()` (one `rowsum()`, hubs exact), `.na_calibrate()` (pure), `.na_additive_stage()` (correlation scale, floor from the residual coupling, budget `1e-10 max(1, ||floor_s||)`), `.na_solve_dd_mean()` (robust solver). |
| Anchor objects | `R/qtl_congruence.R` | `cross()` on both kinds; the design anchor caches its SVD; `.qtl_anchor_rank_check()` can name the block; new `.qtl_rank_note()`. |
| Part A, shared | `R/define_additive_effects.R` | `.dae_anchor_at()`, `.dae_calibrate_shared()`, `.dae_build_traits()` (the additive-only path both generators call, C4 (b)); `.tvc_block_from_rows()`, `.tvc_collect_explicit()`, `.tvc_hint_filter()` (the per-block target resolver both share; Part A's messages unchanged for `additive`); `.dae_refuse_nonadditive_model()` (5b.6); the A16 third fix; the rank note. |
| Owner rules elsewhere | `R/define_genome_effect_terms.R`, `R/define_effect_cov_matrix.R`, `R/remove_generated_effects.R` | `replace_trait` refusal names both generators; `define_effect_cov_matrix()` emits the rank note and names `define_genome_effects()` as the regeneration route for non-additive blocks; removal docs. |
| Tests | `test-define_genome_effects.R` (29 tests), `test-genome-effects-calibration.R` (31 tests), `helper-nonadd-generator-oracle.R` | Below. |

### Gates

**Internals** (`test-genome-effects-calibration.R`, no database). The source suite on a
300 × 120 LD panel with 30 % inbred rows, architectures drawn in the test and supplied:
source tests 1–6, 9–15, 17, 18 as the plan's table says (7, 8 not ported; 16 is public).
Test 1 checks each delivered block three ways: `.na_calibrate()`'s own anchor, the
source's `nonadd_decompose()` on our functional `(a, d, e)` (its `real_*` or `genic_*`
blocks, and its identity `g = const + BV + DD + AA`), and the source generator on the
same supplied architectures. C6 against `zeng_appendix_A()` at 1e-12. C4 (a): with no
non-zero D / A×A block the result is `identical()` to `.qtl_calibrate()`. G4: the nine
solver cases (degenerate `ρ = ±1`; `A = 0` with the linear root `μ = 0.1`, built from a
symmetric pair of loci; Codex's `B = 0` case `μ = ±0.1`; the preferred mean on `ID = 0`
giving 0.29 / 0.09; a repeated root at `x = 0`; `ρ = 0` with zero variance at the root;
`ρ = 0`, negative, just inside and just outside the bound; an unattainable one-locus
request; the `sd` sweep 1e-8 … 1e8 with the same branch and `x`). G6: per-trait units
`S = diag(1e-6, 1e6)`, the floor computed independently from the residual coupling
(on / just above / below it), a small-unit trait below its own floor beside a large-unit
trait far above it, and the `G_A = 4` hand case (floor 1, standardised residual 0.75).
The oracle is isolated (parent `baseenv()`, every function's environment checked) and
runs on its own.

**Public** (`test-define_genome_effects.R`). C1, C3, C9 from the stored rows and the
evaluator; C2 over all 27 genotypes of three loci and one pair, both anchors, at 1e-12
(plus the realised round trip, below); C4 (b) for `G = 0.7`, a 2 × 2 PD target, `0`,
`diag(1, 0)` and a rank-one target, both anchors: `identical()` rows, `trait_var_comp`
**and post-call `.Random.seed`**; explicit zero `G_D` / `G_AA`; replacement with
custom owners. C5 / C11 (floor error; injected failure on the second target write
restores all four tables, scoped rows and custom owners included). C7, C12 (exact
counts `n(p² + Fpq, …)` at `F = 0.5`, mean dominance value `−F Σ 2pq d` at 1e-12, and
individual values). C8 (a)–(d), C10, C14, C15 (a)–(d), C17 with the D2 half on
`.dge_draw()` directly, C18. C20: every knowable refusal checked with and without an
existing `.Random.seed` (`expect_refusal()`), database unchanged. G1 (D5 and the
three rank-note reasons from all three entry points), G2 (threads 1 / 8, and a file
database after `restore_pop()`), G3, G5, G7 (mocked `QTL_REALISED_MAX_CELLS = 4000`),
G8, D3 (a) with a locus fixed in the base pool, and source test 16 as D6.

### Mutation checks

Scripted (`scratchpad/mutate5b.R`, files restored in `finally`, restoration checked
with `cmp`). Every mutation fails its gate:

| Mutation | Fails |
|---|---|
| Drop `E_c` from `C` | C1/C3/C9, C2, G3, C15 (a) |
| Observed `b` in storage under `"realised"` | C2 (the realised round trip) |
| Floor test on the raw scale | G6 (on / above / below the floor) |
| Floor as `R − QᵀP⁻¹Q` with the source's tolerance | G6 |
| Verbatim `solve_dd_mean()` | G4 (six of the cases) |
| Additive-only storage through `.noia_terms()` | C4 (b) |
| Solve unrequested traits' degree means too | G5 |
| D5 check without the zero-block exception | G1 |
| Draw pairs before dominance | C17 / D2 |
| Additive-only through the quadratic | C4 (b) (all three), G1 |
| Skip the D5 check | G1, C20 |
| `define_additive_effects()` content refusal removed | C8 (b) |

### Forward conversion

Profiled as the plan asked. `.noia_to_stored()` + `.noia_terms()`: 1,000 QTL with 500
pairs, 0.03 + 0.03 s; 1,000 QTL with 19,900 hub pairs on 200 loci, 0.79 + 0.54 s. Not
vectorised: it does not show. The generator benchmark itself is 5c.

### Deviations from the plan

- **The repeated root is taken at the vertex.** The plan's `q = 0` branch covers only
  `B = 0`, `disc = 0` exactly. Probing showed that any discriminant within rounding of 0
  breaks the stable formula: `√(residue)` moves both roots by about `√eps` (so `ρ = 0`
  was refused on a clean case), and `C / q` of two residues is noise. A discriminant
  within the same budget is now one root, `−B / (2A)`, which also covers the plan's
  `q = 0` case.
- **G4's bound cases use the attainable side.** The line `d = x u + v` reaches the
  Cauchy–Schwarz bound on one side only (where the maximising direction puts positive
  weight on `v`), so the fixture's bound carries that sign. The repeated-root test asserts
  `|x| < 1e-3` and the exact ratio: a double root's location is ill-conditioned.
- **C2's mutation guard.** Recovering functional coefficients from the stored ones is
  self-consistent under the "observed `b` in storage" mutation, so the identity alone
  cannot catch it. The C2 test adds the realised round trip
  (`extract_genetic_variance(anchor = "realised")` on the base returns the targets at
  1e-10), which does; that is the core of 5c's C20, run here on a 54-individual fixture.
- **One-locus inbreeding depression.** With one locus `ID / √V_D = sign(d)`, so the public
  one-locus case requests `±√G_D` (the solver's degenerate branch), not an arbitrary
  value.
- **D6 under `"genic"` with individuals** collects the dosages again for the comparison
  (Part A does the same); above the size limit it is skipped with a message.
- **Messages.** Calibration errors from the dominance and A×A stages now name the block
  and the anchor (the plan asks every error to name the anchor). Small delivered
  off-diagonals print as round-off (e.g. `-3.7e-17`), as Part A's do.

### Suite

Full suite with `NOT_CRAN=true` (`testthat::test_dir()`, 78 files): **1,159 tests, 4,779
expectations, 0 failures, 0 errors**, 6 warnings, the 0.75.2 baseline in the same five
files (`add_founders`, `add_phenotype`, `genome_map`, `parity` ×2, `phenotype_composite`).
`devtools::document()` and `pkgdown::check_pkgdown()` clean. The 484 new expectations are
the two new test files.

### Implementation-review follow-up (0.76.1)

[Codex reviewed the built 5b](import_qtl_effect_methods_phase_5_codex_review.md): the planned
gates and the full suite pass, and its independent scientific checks (oracle integrity, a
public trait-unit sweep to `u = 1e-10`, three traits with hubs, a fixed partner and an
inbred cohort under both anchors, individual values against `add_tgv()`) found no error at
ordinary scale. Two medium findings, both accepted:

- **F1: the stored additive coefficients were not the verified ones.** `.na_calibrate()`
  verified `B_alpha`, then returned the functional `a = B_alpha - C`; storage
  (`.noia_to_stored()`) added the coupling back. When `B_alpha` is tiny beside `C` that
  cancels: `G_A = 1e-24`, `G_D = 1`, one QTL at `p = 0.3`, seed 2 stored a variance off by
  −2.0e-4 behind a message reporting it exact (reproduced through the public path, other
  seeds at the same ratio +8.5e-5). Fixed by `.na_store_alpha()`, which never removes and
  re-adds the coupling: under `"genic"` the anchor's `b`, `c` are the HWE ones at the stored
  centre, so `B_alpha` is stored as is (now exact, relative error 0); under `"realised"` the
  stored coefficient is `B_alpha + Delta`, `Delta` the coupling (linear in `b`, `c`) of the
  HWE-minus-observed `b` and `c`, and the realised coefficient the stored model implies
  (stored − `Delta`) is verified against the target again; a miss is refused ("cannot be
  stored exactly", e.g. an out-of-HWE locus at `G_A = 1e-24`). `.dge_build_nonadditive()`
  writes these in place of `.noia_to_stored()`'s recovered `alpha`; the delivered covariance
  in the message is measured from them.
- **F2: the realised size guard undercounted kept designs.** The guard counts dominance and
  pair columns only for non-zero blocks, but `.na_anchors("realised")` always built the
  dominance design and a zero `G_AA` still built the pair design (A + A×A, or D + zero A×A,
  kept 800 cells where the guard counted 480 / 640). Fixed on the allocation side, so the
  documented contract stands: `.na_anchors(dominance =)` builds `D` only for a non-zero
  block, the A×A anchor is built only for a non-zero block, and `.na_calibrate()` gives a
  zero block zero coefficients and a zero delivered covariance without an anchor (its
  draws are still consumed, so C17 holds). The genic-with-individuals diagnostic builds the
  dominance design only when it compares that block.

**Tests.** Public (`test-define_genome_effects.R`): G3 at `G_A = 1e-24` / `G_D = 1` checks
the stored genic variance to 1e-12; G7 captures the anchors handed to `.na_calibrate()`
for A + A×A, D + zero A×A and zero D + A×A (exactly `n m`, `n m_D`, `n r` cells), admits
A + A×A at `n (m + r)` and refuses D + zero A×A one cell over `n (m + m_D)`, with database
and RNG unchanged. Calibration (`test-genome-effects-calibration.R`): the short genic
reproduction, the realised `Delta` identity and its refusal, and zero blocks without
anchors. Mutations: storing `.noia_to_stored()`'s `alpha` again fails the G3 test; always
building the dominance design fails G7.

**Not changed: the extractor at such ratios.** `extract_genetic_variance()` canonicalises
stored terms to functional effects and re-projects (`.egv_coefficients()`, `.egv_alpha()`),
the same subtraction and re-addition. On the fixed model above it reports the additive
block as `1.000085e-24` (the stored model is exact). Its precision is absolute, about
`eps` times the coupling, which is ordinary for a decomposition whose own consistency check
is relative to the total genetic value; an additive block twenty-four orders below the
dominance block is below that resolution. Decided with the user: documented, not changed
(a new "Precision" section in `?extract_genetic_variance`).

**Suite.** Full suite with `NOT_CRAN=true` (78 files): **1,162 tests, 4,805 expectations, 0
failures, 0 errors**, the same 6 baseline warnings in the same five files. No exported
signature or roxygen changed.

## 5c — End to end and the vignette (0.76.2)

Tests, a benchmark and a vignette over the code of 5a and 5b; no function's behaviour
changed. Codex's re-verification of 0.76.1 (both findings closed, "proceed to 5c") is
committed with it.

### Gates — `tests/testthat/test-define_genome_effects-integration.R` (8 tests)

**Fixture.** An exact HWE + LE cohort, built as the full factorial of per-locus genotype
counts (1, 2, 1) at `p = 1/2` and (9, 6, 1) at `p = 1/4`: three loci, 256 individuals. Every
function of one locus is uncorrelated with every function of the others and each locus is
in exact HWE, so every realised block is exactly `n/(n − 1)` times its genic value. The
`p = 1/4` locus keeps the coupling `c = 2p − 1` non-zero. The other fixtures are the
200- or 400-founder population of the 5b tests (sampled from 300 haplotypes, so in LD), and
an inbred version of it.

| Gate | What is checked |
|---|---|
| **C13** | A + D (`G_A = 1`, `G_D = 2`, realised, 400 individuals): the record minus its residual minus the mean is `ind_tgv_total` per individual (1e-12); `var(ind_tgv_total)` is the extractor's realised `total` (1e-10); the blocks are the targets. Coarse: the phenotypic minus the residual variance is in (2, 4), i.e. near `V_A + V_D = 3`, not `V_A = 1`. |
| **Prevalence (a)** | A + D + A×A on the factorial fixture: the stored centres are the exact frequencies; the genic blocks are the targets (1e-10); realised `between_components` is 0 and `total` is `n/(n − 1) × 1.5`; the threshold's sum is 1.5, i.e. `(n − 1)/n` times the realised total. `add_phenotype()` on a `prevalence` phenotype records category 2 exactly when the stored liability exceeds `qnorm(0.8) √(1.5 + 1)`. |
| **Prevalence (b)** | On the LD founders: blocks + `between_components` = `total` (1e-10), with `|between_components| > 1e-4` and a block more than 1e-3 off its target, while the threshold still sums the targets (1.5). |
| **Prevalence (c)** | `define_additive_effects(parent_origin = 1)` beside a generated D / A×A model is refused (the 5b.6 refusal), nothing written; the common-scope model gives the threshold 1.5. |
| **Removal** | `remove_generated_effects()` on a real A + D + A×A model (term kinds `additive`, `dominance`, `interaction`) leaves only the custom term and the targets unchanged; `define_additive_effects()` then works and writes additive terms only. |
| **C19** | Two traits, A + D, genic, on the factorial fixture: genic blocks and cross-covariances are the targets (1e-10); `cov` of `ind_tgv`'s `additive` component equals the realised additive block (1e-10) and is `n/(n − 1)` times `G_A`. Off the fixture (half the cohort fully inbred) the two additive measures differ by more than 1e-3. |
| **C20** | Two traits, A + D + A×A, realised, generated and measured on the same 200 individuals: all nine block entries are the targets (1e-10), and the 12 rows `inner_join()` the stored targets one to one. |

### Benchmark — `dev/benchmarks/benchmark_define_genome_effects.R`

2,000 individuals, 1,000 QTL, k = 2, `G_A`, `G_D`, `G_AA` all 2 × 2. Each call timed whole
(resolution, draw, calibration, diagnostics, write, messages). Peak is R's `gc()` "max
used", DuckDB excluded.

| Run | Time | Peak R memory |
|---|---:|---:|
| genic, 500 random pairs (founder-pool base) | 2.7 s | 612 Mb |
| realised, 500 random pairs (designs 2,000 × 2,500 = 5e6 cells) | 10.5 s | 533 Mb |
| genic, 20,000 supplied hub pairs (21 hubs; 44,000 stored terms) | 7.7 s | 693 Mb |

All three are well within interactive use, so nothing was optimised. The realised run was
not profiled.

### Vignette — `vignettes/genetic-models.Rmd`

The four paths in order: known coefficients (§9.3's worked A + D + A×A through
`define_genome_effect_terms()`), `define_additive_effects()` with a two-trait `G`,
`define_genome_effects()` with the messages it prints, and the floor refusal shown live,
then `extract_genetic_variance()` joined to `trait_var_comp`. It closes with the scope
promises. A 60-locus `:memory:` population; it renders in about 5 s.

### Deviations from the plan

- **C20's `warn_bounds` check** is not repeated: the 5b test "16 / D6" already fires it on
  an inbred panel under both anchors and silences it with `NULL`. A comment points there.
- **Removal, then `define_additive_effects()`** needs `trait_var_comp_tbl` limited to the
  additive rows: the stored D and A×A targets outlive the terms (as documented), and the
  additive generator refuses to ignore a stored target silently. What the plan asked for
  holds: the 5b.6 refusal of a non-additive *model* no longer applies.
- **`_pkgdown.yml`** has no `articles:` index (the introduction vignette was never listed
  either); pkgdown lists every vignette on its own, and `check_pkgdown()` is clean.
- **The vignette** says functional `d` is stored as an `indicator` term, so the hand-written
  model's `ind_tgv` components are `additive`, `indicator` and `interaction`. The first
  draft said `dominance`; rendering caught it.
- **No mutation checks** for these gates. They measure code whose own gates were
  mutation-checked in 4 and 5b; each asserts an exact identity, not a loose bound, except
  C13's single coarse check.

### Suite

Full suite with `NOT_CRAN=true` (79 files): **1,170 tests, 4,855 expectations, 0 failures,
0 errors**, the six baseline warnings in the same five files. `pkgdown::check_pkgdown()`
clean; the vignette renders against `load_all()`. `R CMD check` was not run.

### Implementation-review follow-up (0.76.3)

[Codex reviewed 5c](import_qtl_effect_methods_phase_5_codex_review.md): every gate passes,
its independent checks (two traits at exact HWE + LE under both anchors; three traits on an
inbred base with hubs and a fixed partner; replacement and repeated records; a liability
on the cutpoint) found no scientific or implementation defect, and the benchmark and the
vignette ran. One low finding, accepted: **`?remove_generated_effects` overpromised
recovery.**

- It said a target without terms "blocks nothing" and that `define_additive_effects()`
  "then accepts the trait again". Both are true of the prevalence threshold and of the
  non-additive-model refusal, but a stored D or A×A target still makes
  `define_additive_effects()` refuse unless `trait_var_comp_tbl` selects the additive
  rows (5c's own test does this; the deviation note above explains it).
- It said `ind_tgv` values are stale until `add_tgv()` re-evaluates. When the removal
  leaves **no terms at all**, `add_tgv()` refuses ("No genome effects found") and so does
  `add_phenotype()`, so the old values stay in `ind_tgv` until a new model is written.
  Nothing is ever recorded from them; a reader of `get_table("ind_tgv")` still sees them.

Fixed in the documentation (now two bullets, targets and `ind_tgv`; the paragraph about a
`define_genome_effects()` model moved out of `@return`, where it had been filed by mistake)
and in the API skill. The three "No genome effects found" errors (`add_tgv()`,
`add_phenotype()`, `extract_genotypes()`) now name `define_genome_effects()` too. A new
integration test removes a trait's whole generated model and checks that `add_tgv()` refuses
naming the generators, `ind_tgv` is unchanged, an unfiltered `define_additive_effects()` is
refused naming the stored targets, and `define_genome_effects()` to the kept targets then
refreshes every value. The empty-model refusal itself is unchanged, as Codex recommended.

Codex's C13 note is right and already how the test reads: under LD the realised total
includes covariances between components, so the coarse bound is a sanity check on that
fixture, not an equality.

**Suite.** Full suite with `NOT_CRAN=true` (79 files): **1,171 tests, 4,862 expectations, 0
failures, 0 errors**, the six baseline warnings. `devtools::document()` regenerated
`man/remove_generated_effects.Rd`.
