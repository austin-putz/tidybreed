# Review of `update_genome_effects_base_tbl.md`

**Review date:** 2026-09-20  
**Reviewed proposal:** `plans/update_genome_effects_base_tbl.md` (through v4)  
**Status:** approved for implementation.

## Update after v4 — the two-function boundary

V4 resolves all three amendments from the v3 re-review and the three non-blocking
wording notes:

| Re-review item | v4 result |
|---|---|
| Obsolete `founder_haplotypes.line_name` projection requirement | Removed in §2.3 and §3.1 |
| Omitted `center_value` failed to trigger a fill | Corrected in §3.4; both omitted and explicit-`NA` cases are tested |
| Realized line-origin copies were compared with a founder-pool expectation | Replaced by a deterministic, realized-copy comparison in test 17 |
| Partial versus wholly empty base contract | Stated explicitly in §2.2 |
| Fallback query count | Corrected in Q10 |
| Explicit pooled call called a default call | Corrected in test 16 |

I also reviewed the new §2.5 decision about whether to merge
`define_additive_effects()` and `define_genome_effects()`.

### Verdict: keep both functions

I agree with the split. They are two user-facing operations at different levels,
not duplicate spellings of one operation:

- `define_genome_effects()` is a deterministic general writer. The caller owns the
  terms, coefficients, owner, origin predicates, and replacement mode.
- `define_additive_effects()` is a model generator. It selects QTL, samples or
  accepts additive coefficients, optionally realizes a multi-trait covariance,
  scales to target variance, fixes the reserved TBV owner, and replaces one scope.

Combining them would require dispatching on the presence of `terms` versus a
filtered `genome_meta` pipe subject, expose mutually exclusive sampling and writing
arguments, make RNG conditional on argument combinations, and make reserved-owner
and replacement rules mode-dependent. That would be harder to explain and easier to
misuse than two verbs with one shared engine.

The important seam is internal, and v4 puts it in the right place: both paths use
the same `.ge_build()` / resolution / commit machinery and the same
`extract_allele_freq()` interpretation of `base_tbl`. The surface distinction is
also clear:

- use `define_additive_effects()` when asking tidybreed to generate an additive
  architecture; and
- use `define_genome_effects()` when supplying an architecture yourself.

This remains coherent even though `base_tbl = NULL` differs intentionally: the
generator has a domain default and resolves the appropriate founder pool, while the
general writer has no basis for inventing a population and therefore requires
explicit centres or an explicit `base_tbl`.

### Acceptance-test clarification

Test 22 is the right drift guard, but its writer half should **omit
`center_value`** (or set it to `NA`) when reconstructing the additive terms. If it
copies the stored centres from `genome_effect_loci`, the supplied `base_tbl` is
validated but deliberately not queried, so the test proves storage-path parity but
not the claimed shared base semantics. The strongest version asserts both:

1. captured generated `genome_value`s + omitted centres + the same `base_tbl`
   reproduce the three effect tables (surrogate ids aside); and
2. `add_tbv()` produces identical values for common and composed
   `line_name`/`parent_origin` scopes.

This is a test-specification clarification, not a design blocker. With it, v4 is
approved for implementation. I do not recommend merging the two public functions.

## Update after the maintainer / Claude response (v3)

The response accepted and incorporated all five required revisions from the first
review. I rechecked the revised validator, SQL shape, `.ge_build()` integration,
call-site inventory, acceptance tests, and the new multi-line walkthrough against
the current code. The architecture is ready; the remaining issues are local plan
corrections, not objections to `base_tbl` or to the new default behavior.

### Disposition of the original review

| Original point | v3 disposition |
|---|---|
| Shared validation, including connection provenance and projected columns | Accepted in §3.1 and applied to both writers |
| Fill centres inside `.ge_build()` before member validation | Accepted in §3.4 |
| Correct the `ad_terms()` workflow | Accepted in §2.4 / Q6 |
| Tighten SQL, empty-result language, and founder diagnostics | Accepted in §3.2–3.3 |
| Correct the migration blast radius and expand tests | Accepted in §3.6–3.8 |

The response also correctly replaced the rationale for retaining
`genome_meta.founder_allele_freq`: it remains only to keep this change scoped, not
because it is a reliable multi-line snapshot.

### New Q10 — shared-pool fallback

**Agree with (a): a line-scoped implicit default should resolve its named founder
pool first and the shared (`line_name IS NULL`) pool second.** This matches the
meaning already assigned to the shared pool by `add_founders()`, makes the common
one-pool/multiple-lines setup work, and remains loud when neither pool exists. The
named pool must win whenever present, as §3.3 specifies.

This is pool-level fallback, not per-locus stitching: once a named pool exists, an
incomplete named pool must not be silently completed from the shared pool. The
plan's later QTL-level missing-copy check provides the correct failure behavior for
that malformed case.

### New Q11 — Wahlund warning boundary

**Agree with (b): warn only for an implicit, population-wide founder default that
contains multiple pools.** An explicit `base_tbl` is an intentional population
selection, so warning on every explicit pooled call would be noise. Keeping
`extract_allele_freq()` free of modelling warnings also gives the exported helper a
clean analysis contract. Counting the unnamed pool as a distinct pool is correct.

Moving the warning has one consequence for §3.1: `line_name` is no longer a
required projected column for an explicit `founder_haplotypes` base. See amendment
1.

### Required amendments to v3

#### 1. Remove `line_name` from the public helper's projected-column requirement

For `founder_haplotypes`, `extract_allele_freq()` needs only `locus_name` and
`allele`. In v3, neither remaining line-aware operation reads `line_name` from the
rendered subquery:

- the Wahlund query runs only in `.dae_default_base()` and reads the physical
  `founder_haplotypes` table; and
- the empty-result helper lists available pools from the physical table as well.

Therefore this valid call should not be rejected:

```r
get_table(pop, "founder_haplotypes") |>
  dplyr::select(locus_name, allele) |>
  extract_allele_freq()
```

Change the §2.3 / §3.1 requirement to `c("locus_name", "allele")`. The implicit
default still receives the unprojected table, so Q10/Q11 behavior is unaffected.
This also aligns the response with the original qualification: require `line_name`
only if the helper itself retains the Wahlund diagnostic; v3 deliberately moved it
out.

#### 2. Fix missing-column detection before deciding whether to query

The §2.4 contract and example allow hand-written terms to omit `center_value`, but
the §3.4 pseudocode currently does this:

```r
any(is.na(terms$center_value) &
    terms$contrast_name %in% c("additive", "dominance"))
```

When `center_value` is absent, `terms$center_value` is `NULL`, the expression has
length zero, and `any()` returns `FALSE`. No frequencies are queried; `.ge_build()`
then creates the missing column as `NA` and raises the old error. In other words,
the plan's showcased hand-written dominance example would fail.

Normalize for this decision exactly as `.ge_build()` will:

```r
centres <- if ("center_value" %in% names(terms)) {
  terms$center_value
} else {
  rep(NA_real_, nrow(terms))
}
needs_fill <- !is.null(base_tbl) &&
  "contrast_name" %in% names(terms) &&
  any(is.na(centres) &
      terms$contrast_name %in% c("additive", "dominance"))
```

The cleanest implementation may instead extract the small terms-normalization
stage so malformed terms are rejected before any database query. Either shape is
fine; the acceptance test must cover both an explicitly `NA` centre and an omitted
`center_value` column.

#### 3. Replace the non-deterministic assertion in acceptance test 17

`ind_haplotype |> filter(line_origin == "Duroc")` in a realized population does
not generally equal the Duroc founder-pool frequency. Those copies are a finite,
transmitted sample of the pool; equality holds only in expectation over repeated
simulations, which is not a stable unit-test assertion.

Instead, hand-compute `AVG(allele)` from the selected realized Duroc-origin copies
and compare `extract_allele_freq()` with that result. Separately prove that an
`ind_meta |> filter(line_name == "F1")` base includes both origins and matches the
hand-computed frequency over all copies of those F1 animals. A deliberately fixed
fixture can then assert the two selections differ.

### Minor wording cleanups (non-blocking)

- The §2.2 public contract should explicitly say that a partly covered base returns
  per-locus `NA`, while a base with no usable matches at any locus errors. The code
  and tests already specify this; the headline contract currently mentions only
  the first half.
- §3.3 says the fallback "adds one `EXISTS` query"; the missing-named-pool path
  executes two (`has_line`, then `has_shared`). This has no design consequence.
- Acceptance test 16 calls all three writes "default-base calls", although the
  common fallback write deliberately supplies an explicit pooled `base_tbl`.

After the three amendments above, I approve implementation. Q10 and Q11 should be
recorded as accepted, and no further schema or evaluator review is needed for this
change.

## Original executive verdict (review of v1)

The section below is retained as the audit trail for the first review. Its five
required revisions are resolved by v3 as summarized above.

The main simplification is sound: replace the coupled `base`, `base_tbl`, and
`base_line_name` controls with one table-shaped `base_tbl`, while retaining a
line-aware default. This fits the package's existing filtered-table idiom, removes a
fragile `missing()`-based sentinel, permits allele-copy selection that the current API
cannot express, and eliminates the current collect-IDs-then-build-an-`IN`-list path.

I agree with recommendations Q1-Q7 and Q9, subject to the qualifications below. I
agree that Q8 should stay out of this change, but I do **not** agree with the stated
rationale for keeping `genome_meta.founder_allele_freq`: in a multi-line population it
is a last-call-wins value, not an unambiguous snapshot of the complete founder pool.

The plan is close, but it should be revised before implementation in five areas:

1. specify one shared validation contract for `base_tbl`, including connection
   provenance and required projected columns;
2. show the actual `define_genome_effects()` / `.ge_build()` refactor needed to fill
   centres before validation;
3. correct the claims about preserved empty-line diagnostics and the call-site blast
   radius;
4. make the SQL composition and empty-result behavior more precise;
5. expand acceptance tests for repeated individual rows, selected-away columns, and
   both writers' cross-population rejection.

No change is needed to the proposed storage or evaluation model.

## What the proposal gets right

- `base_tbl` makes the population selection explicit without inventing another enum.
- Keeping `genome_meta` as the pipe subject is the right choice. Locus selection and
  base-population selection are independent inputs.
- The three accepted shapes represent genuinely different useful selections:
  founder-pool haplotypes, arbitrary allele copies, and individuals.
- `SELECT DISTINCT id_ind` in the individual-table path gives correct set semantics
  even for tables such as `ind_phenotype`, where one individual may have many rows.
- Returning `NA` for a locus with no selected copies is much safer than the current
  zero-initialized vector.
- Rejecting a selected QTL with no base copies is preferable to dropping it or
  inventing `p = 0.5`.
- `extract_allele_freq()` is a useful public analysis helper and avoids making
  `ad_terms()` database-aware.
- Explicit `center_value` winning over an optional fill is a clear precedence rule.
- Keeping the `add_ebv()` cleanup and `founder_allele_freq` schema question out of
  this implementation prevents avoidable scope expansion.

## Required revisions

### 1. Define and reuse a complete `base_tbl` validator

The sketch validates connection identity only in `define_additive_effects()`. The
same requirement applies to `define_genome_effects()`, because its `base_tbl` is also
rendered and executed against the writer population. Put this in one internal helper,
used by both writers before calling `extract_allele_freq()`.

The helper should validate:

- `inherits(base_tbl, "tidybreed_table")`;
- `validate_tidybreed_pop(base_tbl$pop)`;
- `identical(base_tbl$pop$db_conn, pop$db_conn)` for either writer;
- required columns after any user `select()`:
  - `founder_haplotypes`: `locus_name`, `allele`, and `line_name` if the Wahlund
    diagnostic is retained;
  - `ind_haplotype`: `locus_id` and `allele`;
  - every other table: `id_ind`;
- a useful error naming `base_tbl$table_name` and the missing columns.

This matters because `select.tidybreed_table()` preserves the wrapper while changing
the lazy projection. Dispatching only on `table_name` is therefore insufficient: a
valid `tidybreed_table` can still lack columns needed by the generated SQL.

The default line filter should also use unambiguous tidy evaluation, for example
`.data$line_name == .env$line_name`, rather than relying on two meanings of
`line_name` in the same expression.

### 2. Show the real `.ge_build()` integration point

The proposed fill cannot be added outside `.ge_build()` as the current code is
structured. `.ge_build()` resolves `locus_name -> locus_id`, creates `members`, and
immediately calls `.ge_check_member_fields()`, which rejects missing centres before
returning `built`.

The implementation plan should explicitly choose one of these shapes:

```r
base_freq <- if (is.null(base_tbl)) NULL else extract_allele_freq(base_tbl)
built <- .ge_build(conn, trait_name, terms, origin, effect_owner,
                   base_freq = base_freq)
```

with the fill inside `.ge_build()` after `members` are created and before
`.ge_check_member_fields()`, or a small extracted normalization stage that returns
unvalidated members for the writer to fill and then validate. The first option is the
smaller change.

Compute frequencies only when at least one additive/dominance input row has a missing
`center_value`; an explicit `base_tbl` should not cause an unnecessary full-genome
query when every centre is already supplied.

The resulting error should retain the existing `term_id` and locus context and append
that the supplied base contains no allele copies at that locus. Do not replace the
current row-specific validation with a generic list of loci.

### 3. Correct the `ad_terms()` example and contract discussion

The Q5 discussion briefly suggests `ad_terms(..., p = NA)` and then retracts it. The
current builder rejects any `NA` through `.ad_recycle()`, so the final plan should
simply state:

- use `extract_allele_freq()` when building terms through `ad_terms()`;
- use a hand-written terms data frame with missing `center_value` when intentionally
  asking `define_genome_effects(base_tbl = ...)` to fill centres.

This keeps `ad_terms()` pure and avoids implying an unsupported workflow.

### 4. Tighten the SQL and empty-result specification

Prefer `dbplyr::sql_render(tbl$tbl)` for rendering the lazy query unless the pinned
dbplyr version specifically requires `remote_query()`. `sql_render()` communicates
the public operation intended here more directly.

The individual path should be described as a semi-join, even if implemented with
`IN (SELECT DISTINCT ...)`. A join against a distinct ID subquery is also clear:

```sql
FROM ind_haplotype h
JOIN (SELECT DISTINCT id_ind FROM (<filtered base>) b) ids USING (id_ind)
```

The final left join to `genome_meta` should remain: it is what guarantees one result
row per locus and distinguishes "no copies" from an observed frequency of zero.

`all(is.na(out$allele_freq))` establishes that no usable locus matched; it does not
prove that the filtered source itself had zero rows. The error should therefore say
"the filtered base contains no usable allele copies matching `genome_meta`", or the
query should also return enough counts to distinguish an empty selection from broken
locus keys.

The plan currently says that the existing line-specific error listing available
founder lines is preserved by the generic empty check. It is not. A generic all-`NA`
check cannot produce the old "Available: ..." diagnostic. Either retain a specialized
founder diagnostic on failure or explicitly accept the simpler new error.

Also state that the frequency query is one SQL statement, while the optional Wahlund
diagnostic may be a second statement. That is still a substantial improvement over
collecting IDs into R.

### 5. Correct the blast-radius inventory

The claimed "26 lines" is too small for the current tree. A repository search finds
roughly twice that many live references across implementation, roxygen, tests, and
the swine vignette. In particular, the vignette has multiple executable
`base = "current_pop"` calls, not one.

This distinction matters because removing `base = "current_pop"` without supplying
an explicit `base_tbl` silently changes those calls to the founder-pool default. The
migration inventory should be generated with a repository-wide search and recorded
by file/call, not summarized as a small fixed line count.

The migration rule itself is correct:

```r
base = "current_pop"
```

becomes an explicit whole-population selection such as:

```r
base_tbl = get_table(pop, "ind_meta")
```

Do not merely delete the old argument in these cases.

## Additional implementation improvements

### Preserve selection semantics, not incidental row multiplicity

For a table carrying `id_ind`, frequencies must depend only on the distinct selected
individuals. Add a test proving that `ind_phenotype` with repeated records produces
the same frequencies as `ind_meta` for the same IDs. This guards against a future
rewrite from accidentally weighting individuals by their record counts.

### Keep the output order and type contractual

The public helper should promise:

- exactly one row per `genome_meta` locus;
- ascending `locus_id` order;
- `locus_id` as integer, `locus_name` as character, and `allele_freq` as double;
- `NA_real_` for loci with no selected copies;
- no database writes and no population mutation.

Validate that every non-missing returned frequency is finite and in `[0, 1]`. The
core haplotype tables should guarantee this, but the public helper should fail clearly
if database corruption violates that invariant.

### Avoid overclaiming parity

The default line-aware founder behavior is numerically unchanged for valid, complete
founder pools. Diagnostics for malformed or empty selections are changing unless the
specialized old messages are deliberately recreated. Describe this as behavioral
parity on valid inputs, not exact semantic parity in every failure mode.

### Document the warning boundary

The Wahlund warning is proposed only for `founder_haplotypes`, where distinct
`line_name` values are directly available. Explicit `ind_haplotype` or individual
selections can also pool origins or lines, but there is no uniformly valid line field
for all accepted table shapes. State that explicit non-founder bases are treated as
intentional selections and are not diagnosed for pooling.

## Open-question verdicts

### Q1. Required base or line-aware default?

**Agree with (a): keep `base_tbl = NULL` and the line-aware founder default.** It is a
domain default, not merely typing convenience. Requiring the common correct choice at
every call would add noise without improving safety. Keep the pooled-founder warning
when `line_name = NULL`.

### Q2. Argument name?

**Agree with `base_tbl`.** It is concise, retains the existing name, and makes the
expected object class clearer than `base`. `base_pop_tbl` adds little information.

### Q3. Which table kinds?

**Agree with all three shapes**, provided the required-column and same-connection
checks above are part of the contract. Restricting the API to only IDs or only direct
haplotype tables would remove useful and distinct selections.

### Q4. Missing copies at a selected QTL?

**Agree: error and name the affected loci.** Dropping loci changes the requested
model, while `p = 0.5` fabricates data. The check should be limited to loci that will
actually be written for each trait, especially in the multi-trait path.

### Q5. Add `base_tbl` to `define_genome_effects()`?

**Agree with (a), with the `.ge_build()` revision above.** Supplying `base_tbl` is an
explicit opt-in that gives missing `center_value` a clear contextual meaning. Without
`base_tbl`, missing centres must retain their current error. Query only if a fill is
needed, preserve explicit values, and never fill indicators.

### Q6. Add database access to `ad_terms()`?

**Agree with (a): no.** Keep the builder deterministic and population-independent.
The helper-plus-`p` workflow is explicit, and the writer fill covers concise
hand-written terms.

### Q7. Refactor `add_ebv()` now?

**Agree with (a): out of scope.** Record the naming and collection issue in
`plans/TODO.md`, but do not couple an uncovered BLUPF90 refactor to this change.

### Q8. Drop `genome_meta.founder_allele_freq` now?

**Agree that it should be decided separately; disagree with the keep rationale.**
`.write_founder_haplotypes()` documents that the column is last-call-wins for
multi-line populations. It therefore is not a reliable snapshot of the combined pool
and cannot identify which line it describes. Keep it during this change solely to
avoid scope creep, then separately choose among dropping it, making it explicitly
line-keyed, or documenting its limited meaning.

### Q9. Cross-population tables pointing at the same file?

**Agree with strict connection identity.** Same-file/different-connection objects are
not a supported composition boundary, and accepting them complicates ownership and
connection-lifetime reasoning. The plan overstates that rendered SQL is necessarily
wrong on the second connection, but the strict provenance rule is still the safer
public contract. Apply and test it in both writers.

## Revised acceptance checklist

In addition to the proposed tests, add or make explicit:

1. `select()` removing a required column produces a package error naming the table
   and missing column, not a raw DuckDB error.
2. `ind_meta`, repeated-row `ind_phenotype`, and direct `ind_haplotype` selections
   agree when they select the same allele copies.
3. Both `define_additive_effects()` and `define_genome_effects()` reject a
   `base_tbl` from another population/connection.
4. `define_genome_effects()` does not query the base when all centres are explicit.
5. Missing-centre errors retain `term_id` and locus names after a failed fill.
6. A base with some, but not all, loci returns `NA` only at absent loci; a writer
   errors only when an absent locus is actually being written.
7. An observed frequency of exactly `0` or `1` is retained and is not confused with
   no copies.
8. Founder Wahlund warnings use the filtered subquery and the warning boundary for
   non-founder base tables is documented.
9. The public helper is read-only and returns stable locus ordering and column types.
10. Every executable `base = "current_pop"` call site is migrated to an explicit
    whole-population `base_tbl` rather than inheriting the new founder default.

## Recommended implementation order

1. Add the shared validator and `extract_allele_freq()` with focused SQL/contract
   tests.
2. Rewire `define_additive_effects()` and migrate all current-population call sites.
3. Refactor `.ge_build()` to accept optional base frequencies, then expose
   `base_tbl` on `define_genome_effects()`.
4. Update documentation, namespace, NEWS, and version only after both writer paths
   pass the full test suite.
5. Log, but do not implement here, the `add_ebv()` table-query cleanup and the
   `founder_allele_freq` redesign/removal decision.

With these revisions, the proposal is a clear improvement over the current API and
is ready to implement without changing the broader genome-effects architecture.
