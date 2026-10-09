# Import QTL-effect methods — Step 5 results

Plan: [import_qtl_effect_methods_phase_5_plan.md](import_qtl_effect_methods_phase_5_plan.md),
revised twice after the [Codex review](import_qtl_effect_methods_phase_5_codex_review.md)
and its re-review. Decisions D1–D8 were made by the user 2026-10-06 to 2026-10-09, each as
recommended; D5 adds a rank `message()` for singular targets (built in 5b). Three commits
(D1); the user reviews between them.

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
