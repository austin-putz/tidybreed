# Review of `individual_id_sql.md`

Reviewed 2026-11-03 against the current working tree (`DESCRIPTION` reports
0.76.13), `CLAUDE.md`, and the project schema/API skills.

**Recommendation: proceed with the core design, with the refinements below.**
The underscore and prefix failures are real, and registered parent-ID views are
the right fix. I agree with D1(a), D2(a), and D3's recommendation against a
migration. I would revise the numeric-boundary behavior and regression-test
plan before implementation. The malformed-ID policy also needs clarification,
although it does not affect populations produced normally within the supported
numeric range.

This is a proposal review. No production code or original plan was changed.

## Findings to address

### 1. Define numeric limits; do not silently discard an overflowing suffix

**Priority: medium. Plan §2a and §4.**

The proposed expression treats two different situations as `NULL`:

```sql
TRY_CAST(regexp_extract(m.id_ind, '_([0-9]+)$', 1) AS INTEGER)
```

One is an unparseable ID; the other is a perfectly numeric suffix outside
DuckDB's signed 32-bit `INTEGER` range. In an isolated fixture, a line containing
only `Huge_2147483648` produces `max_num = NULL`, so the proposed helper would
restart that line at 1. At the other boundary, `Limit_2147483647` produces the
correct maximum, but R's `2147483647L + 1L` produces `NA` with an overflow warning.
The existing offspring counter also increments with `+ 1L`
(`R/add_offspring.R:499–505`).

This is not a likely population size in ordinary development, and some boundary
failures predate the proposal. It is still worth defining the behavior while
centralizing numbering. Notably, the existing founder ID expression can cross
the boundary: `n_founders` can be an R double, and
`seq(2147483647L, 2147483647L + 2 - 1)` produces both 2147483647 and 2147483648.
Thus §3's statement that every old ID necessarily reads correctly with the new
helper is too strong.

Recommended implementation:

- Keep an explicit supported suffix range. An int32 limit is reasonable for
  this focused change; widening the whole counter path is unnecessary.
- Detect a numeric suffix outside that range rather than ignoring it through
  `TRY_CAST`. Parsing into a wider type can help detect this, but does not by
  itself make the R integer counter safe.
- Check the entire requested allocation per line before RNG or writes:
  `current_max + n_new_for_line` must fit. Perform the check in arithmetic that
  cannot overflow before the comparison. Checking only the first new ID misses
  a batch that crosses the limit.
- Avoid an unnecessary overflowing increment after assigning the final
  permitted ID, or reject conservatively with a documented limit.

Use tiny database fixtures at the boundary; these tests need no large
population. Errors should name the line and supported range and leave the
database unchanged. No persisted counter or migration is needed.

### 2. The claimed malformed-ID guard does not validate the complete format

**Priority: low; clarification/correction. Plan §2a.**

This combination checks a prefix and a suffix independently:

```sql
starts_with(m.id_ind, l.line_name || '_')
regexp_extract(m.id_ind, '_([0-9]+)$', 1)
```

For a row with `line_name = 'Bad'` and `id_ind = 'Bad_extra_99'`, the proposed
query returns 99. The ID does not have the promised `{line_name}_{n}` format for
that row's line. I reproduced this with the exact query from the plan.

Normal writers do not create this mismatch, and `line_name` is reserved against
`mutate_table()`. Therefore this is not another ordinary crossbreeding bug.
However, if the plan deliberately promises protection against hand-edited rows,
the guard should deliver it.

Extract the entire remainder after the exact `line_name || '_'` prefix and
validate that remainder as a positive decimal integer. Then explicitly choose
whether malformed rows are rejected or ignored. I prefer a clear error for
corrupt identity data; silently ignoring it can hide a damaged population. If
that validation is outside this fix, remove the stronger protection claim.
Consider zero and leading-zero suffixes when specifying what “format” means.

### 3. Replace the source-grep test with checks of actual behavior and statements

**Priority: medium. Plan §4, test 7.**

Searching two source files for `id_ind IN (` and `LEFT(id_ind` does not enforce
the rule:

- `WHERE id_ind = '<value>'`, a differently spaced `IN`, or interpolation in a
  helper would pass it.
- A valid `IN (SELECT id_ind FROM registered_view)` would fail it; existing
  `add_tgv()` code uses exactly this legitimate pattern.
- Reading repository R source needs special handling under installed-package
  testing, including the installed-load coverage workflow in `CLAUDE.md`.

Keep a source lint as an optional development check if useful, but make the
regression tests exercise the contract. Capture the executed statements around
the relevant reads, forwarding them to the real database, and verify that test
parent values are not interpolated into SQL. Exercise both the metadata read
and the successful haplotype read. The proposed quoted-*unknown*-parent test
stops at the first read and cannot cover the second one.

Also retain the behavioral unknown-ID test and include an SQL-shaped unknown ID,
such as one containing a quote and comment marker. It must remain a lookup
value, report a missing parent, and leave the population unchanged. Escaping
values would pass that behavioral test too, which is why statement inspection
adds useful coverage for the stronger no-ID-literals rule.

## Implementation and test refinements

These do not require expanding the feature's scope.

1. **Preserve line-to-counter mapping.** `GROUP BY` does not promise the
   order of the helper's result. Construct the named vector from the returned
   `line_name` column and index it by the requested lines. Test an unsorted
   multi-line request with duplicate names and a previously unseen line. Keep
   the existing mating-row order: it determines which offspring gets each ID
   and each per-gamete stream.

2. **Specify parent validation precisely.** Validate both canonical parent
   columns after alias normalization and required-column checks, before
   combining or registering them. Test each column independently for `NA`,
   empty strings, and non-character values. If factors are deliberately
   rejected by the character-only contract, test that error too. Do not impose
   the generated-ID regex on parent lookup values: the quoted-ID test should
   reach the ordinary missing-parent error as proposed.

3. **Test registered-view cleanup on success and error.** Give both new views
   function-scoped cleanup, including the numbering helper's view, and verify
   a failed lookup can be followed by a valid call on the same connection.
   Assert that the internal views do not linger after either path. Fixed
   internal names follow existing repository practice; a new general helper
   framework is not required for this patch.

4. **Make RNG neutrality an assertion.** Two repeated seeded simulations can
   both consume an unwanted extra random draw and still be identical. Add a
   direct RNG-state assertion for the numbering helper and for early parent
   validation/missing-parent failures, which happen before the writer's draws.
   Test code may inspect `.Random.seed`; production code must not restore it.
   Retain the full within-current-code reproducibility tests, sorting collected
   `ind_meta` and `ind_haplotype` by explicit keys before `expect_identical()`.

5. **Test the shared founder/offspring namespace in both directions.** The
   proposed founder → founder → offspring → offspring case is valuable. Add
   founders again after offspring, and test a mixed mating table with
   interleaved `A`, `A_B`, and `A_7` rows. Assert IDs, parent relationships, and
   per-line continuity, rather than just row counts. Include a trailing-
   underscore line such as `A_`; it is accepted by the current regex and writes
   IDs such as `A__1`.

6. **Include a reopen test.** Use a small file-backed population, close and
   restore it, then append into an underscore line. This verifies that the
   numbering derives entirely from stored data, without session-only state.

7. **Bound the uniqueness claim.** “Current maximum” means the maximum among
   retained `ind_meta` rows. Deleting the highest IDs can allow their numbers
   to be reused; `archive_replicate()` deliberately resets numbering between
   replicates (`R/archive_replicate.R:43–48`). Preserve that design and avoid
   implying lifetime/global uniqueness or gap-free numbering after deletions.

## Decisions and scope

| Decision | Review |
|---|---|
| D1(a): per-line numbering | Agree. The line column defines the allocation scope, and the text ID remains an identifier rather than a population-wide ordering mechanism. |
| D2(a): focused fix, explicit follow-up | Agree. Fix the two writers now and inventory the remaining six files in the follow-up. The follow-up must include `assert_ploidy_2()` as well as its callers. |
| D3: no migration | Agree. This is a query/validation change with no necessary schema change. Detect unsupported or malformed identity data explicitly if that policy is adopted; detection is not a compatibility layer. |

One qualification to the D2 rationale: low coverage for `add_ebv.R` and
`blupf90_helpers.R` does not mean all their SQL paths require the external
solver. `ebv_parent_avg()`, pedigree collection, and input-file preparation can
be tested with small database fixtures without launching BLUPF90. That makes
the follow-up more tractable than the aggregate coverage percentages suggest.
The six-file inventory is consistent with my search for the direct ID-list
construction patterns.

I would keep integer locus-ID lists outside this patch's scope. The proposed
one-time before/after check can provide useful evidence, but preserving
previous-version seeded output should remain optional: current-code
determinism and R/C++ parity are the requirements. Passing the same existing
tests before and after a change does not itself demonstrate identical output
between versions. Direct RNG-neutrality assertions and repeated current-code
runs provide the relevant regression protection.

The plan mentions SQL parsing/scale benefits. If implementation work measures
or claims a performance improvement, start with the small deterministic script
under `dev/benchmarks/` required by `CLAUDE.md`, varying parent and line counts.
The correctness fix itself does not need to become a performance project.

For implementation documentation, update the schema skill's `ind_meta` ID
invariant as well as the API skill and roxygen/generated help. Describe the
range and current-maximum semantics consistently. Call §2a a registered-view
query rather than a parameterized query unless it actually uses bound query
parameters.

## Validation performed for this review

- Read the plan, repository rules, relevant project references, both writers,
  SQL/subset utilities, existing registered-ID read pattern, and the listed
  follow-up call sites and tests.
- Ran isolated in-memory probes with DuckDB R 1.5.5 using the plan's exact
  proposed numbering query. `A`, `A_B`, `A_7`, and `Duroc_F1` were isolated and
  parsed correctly; an unseen line returned `NULL` as intended. The malformed
  and overflowing examples behaved as described above.
- Confirmed registration/query/unregistration left R's RNG state identical in
  that probe.
- Loaded the current package with `pkgload::load_all()` and reproduced the
  second `add_founders()` failure for `Duroc_F1` and the parser error for a
  quoted unknown parent. Separately reproduced the prefix-line conversion
  failure with the current SQL expression.

The proposed helper and joins are not implemented, so these probes validate
the diagnosis and proposed SQL behavior, not a completed patch. I did not run
the full package suite for this documentation-only review.
