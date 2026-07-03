# Stage 3 — Dosage Cache & Export Hardening (COMPLETE)

**Status**: Done (tidybreed 0.48.1, 2026-07-02)
**Parent plan**: [refactor_haplotype.md](refactor_haplotype.md) (Stage 3 of 5)
**Precedes**: Stage 4 (non-diploid inheritance)

---

## Outcome

- **Full test suite: 0 failures, 0 errors, 1055 passing** (up from 1024 at Stage 2) — 31 new tests
  across 3 new files and 2 extended existing files.
- **Parity preserved**: `tests/testthat/test-parity.R` still compares against its Stage 1 golden
  files (no recapture) — this stage touches no haplotype/TBV/phenotype-writing code path.
- As the analysis below found, this stage's real deliverable ended up being **hardening + tests**,
  not new features — `add_dosage()` and `extract_genotypes(loci_tbl = ...)` had already shipped in
  Stage 1.

## What shipped

- **`sql_in_list(values, what)`** (`R/sql_utils.R`) — vectorized IN-list builder: escapes embedded
  single quotes (quote-doubling, matching `format_sql_value()`'s VARCHAR branch), rejects empty
  input and `NA`. Replaces every `paste0("'", x, "'", collapse = ", ")` call site in
  `R/add_dosage.R` and `R/extract_genotypes.R` (locus names, QTL/loci-filter names, individual IDs)
  with no change to surrounding SQL-building code.
- **`assert_diploid_only(pop)`** (new `R/ploidy_helpers.R`) — queries `chr_meta` and errors
  (mentioning Stage 4 / `define_chr()`) if any chromosome has non-default `copies_M`/`copies_F`.
  Called at the top of both `add_dosage()` and `extract_genotypes()`, before any dosage/export SQL
  runs.
- **`chip_name` identifier hardening**, with one deviation from the plan's original fix (see
  below): both functions now validate the *constructed* `chip_col` (`is_<chip_name>`) as a SQL
  identifier, not the raw `chip_name` value.
- **New tests**: `tests/testthat/test-ploidy_helpers.R` (guard unit tests),
  `tests/testthat/test-sql_injection_hardening.R` (`sql_in_list()` unit tests, `chip_name`
  rejection, and a realistic end-to-end check that a real `locus_name` containing a quote — e.g.
  from a careless import — still resolves correctly rather than corrupting the query),
  `tests/testthat/test-dosage_extract_parity.R` (cache-vs-direct parity; partial `ind_genotype`
  population never affects `add_tbv()`/`add_phenotype()`, using two identically-seeded pops).
  Extended `test-add_dosage.R` and `test-extract_genotypes.R` with diploid-guard integration tests
  and (for `extract_genotypes()`) four `loci_tbl` tests that didn't exist despite the argument
  already shipping in Stage 1.

## Deviations from this plan (intentional)

- **`chip_name` validation targets the constructed `chip_col`, not the raw value** — the plan's
  draft fix ("call `validate_sql_identifier(chip_name, ...)`") would have broken the
  industry-standard, digit-leading chip-name convention (`"50k"`, `"770k"`) used throughout
  `test-extract_genotypes.R` and the package's own roxygen `@examples`, since
  `SQL_IDENTIFIER_RE` requires starting with a letter. Validating `chip_col <- paste0("is_", chip_name)`
  instead closes the same injection vector (any unsafe character still fails the regex on the
  combined string) while `"50k"` → `"is_50k"` passes unchanged — zero collateral changes to test
  fixtures, docs, or the shared regex used by every other caller. Verified via the full existing
  test suite (all `"50k"`/`"mychip"`/`"HD"` usages still pass) plus a dedicated regression test.
- **Escaping, not parameterized queries** — the plan left this open ("scope the final choice during
  implementation review"). The codebase has zero existing usage of `dbBind`/`params=` anywhere;
  escaping via `sql_in_list()` matches the project's real, existing convention (already used inline
  at `add_phenotype.R:146`) and fixes every call site with a mechanical substitution.

---

## Why this stage looks different from the parent plan's Stage 3 section

The parent plan's Stage 3 scope (lines 122–137) reads as three build items: add `add_dosage()`,
document `ind_genotype` as a partial cache, and make `extract_genotypes()`/BLUPF90 export meet the
six acceptance criteria in "Genotype Matrix Export for GBLUP" (lines 849–931).

**Both `add_dosage()` (`R/add_dosage.R`) and `extract_genotypes()`'s `loci_tbl` argument
(`R/extract_genotypes.R`) already shipped in Stage 1**, ahead of the plan's own sequencing (see
[refactor_haplotype_stage_1.md](refactor_haplotype_stage_1.md), "What shipped" section). Reviewing
the current code against the six acceptance criteria:

| # | Criterion | Current state |
|---|-----------|----------------|
| 1 | One row per requested individual, incl. all-reference | **Met by construction** — `extract_genotypes()` builds an R matrix over `unique(has_ids)` initialized to `0L`, not a SQL `PIVOT`. |
| 2 | One column per requested locus | **Met by construction** — same matrix, columns = `locus_<id>` for every resolved `locus_id`. |
| 3 | Zero (not NULL) for no-non-reference-allele loci | **Met by construction** — matrix starts as `0L`, only cells with a matching `(id_ind, locus_id)` aggregate row get overwritten. |
| 4 | Stable column order by `locus_id` | **Met** — `locus_ids <- sort(locus_ids)` before building `locus_cols`. |
| 5 | No SQL-injection risk from `locus_name`/`chip_name` values | **Not met** — see "Gap 1" below. |
| 6 | Ploidy explicitly out of scope; assume/assert diploid autosomes | **Not met** — no check exists; behavior is silently wrong (not an error) if `chr_meta` ever has non-default rows before Stage 4 ships. |

So Stage 3 does not need to build `add_dosage()` or `loci_tbl` — it needs to **close the two real
gaps (5 and 6)** and add the exit-criteria tests the plan calls for that don't exist yet. This is a
hardening/test pass, not new-feature work.

---

## Gap 1 — SQL-injection hardening (criterion 5)

Both `add_dosage()` and `extract_genotypes()` interpolate user-supplied strings into SQL via
`paste0("'", x, "'", collapse = ", ")` with no escaping and, for `chip_name`, no identifier
validation at all:

- **`add_dosage(tbl, chip_name = ...)`** (`R/add_dosage.R:89-98`): builds
  `chip_col <- paste0("is_", chip_name)` and interpolates it directly as a bare SQL identifier
  (`... WHERE is_<chip_name> = TRUE`). `chip_name` is never passed through
  `validate_sql_identifier()` in this function.
- **`extract_genotypes(tbl, chip_name = ...)`** (`R/extract_genotypes.R:107`): same pattern —
  `chip_col <- paste0("is_", chip_name)` used as a bare identifier, no validation on `chip_name`
  itself (only the separate `col_name` argument is validated, and only when the caller doesn't
  supply their own `col_name`).
- **`add_dosage(locus_names = ...)`** and **`extract_genotypes()`'s `loci_tbl`/`effects_tbl`
  paths**: `locus_names`/`qtl_names`/`sel_names` values are wrapped in single quotes and
  concatenated into `IN (...)` clauses with no quote-escaping. A `locus_name` containing an
  embedded `'` (whether from a malicious caller or just a careless `import_founder_haplotypes()`
  load) breaks the query or changes its meaning.
- Same pattern for `subset_ids`/`has_ids` (`id_ind` values) in both functions.

This mirrors the project's own existing defense pattern (`validate_sql_identifier()` /
`SQL_IDENTIFIER_RE` in `R/sql_utils.R`), already used by `add_tbv()`, `add_phenotype()`,
`add_ebv()`, `define_chip()`, etc. for values interpolated into SQL. `add_dosage()` and
`extract_genotypes()` are the two places in the codebase that skip it for `chip_name`.

**Fix**:
- Call `validate_sql_identifier(chip_name, what = "chip name")` in both functions before building
  `chip_col`.
- For `locus_names` / `id_ind` values (which are legitimate free-text-ish values, not identifiers,
  so `validate_sql_identifier()`'s letter/digit/underscore whitelist is too strict) — reject any
  value containing a single quote with a clear error, or switch to parameterized queries
  (`dbGetQuery(conn, "... IN (?)", params = list(...))` / `glue::glue_sql()`). Parameterization is
  preferable since it fixes the whole class of interpolation call sites at once rather than adding
  a per-value blacklist check; scope the final choice during implementation review.

---

## Gap 2 — Diploid-only guard (criterion 6)

Per the user's decision, add an **explicit runtime check**, not just documentation: query
`chr_meta` and error if any row has `copies_M != 2 OR copies_F != 2` before computing dosage or
extracting a genotype matrix. This makes the "diploid autosomes only until Stage 4" assumption
fail loudly instead of silently producing wrong (summed-across-mismatched-ploidy) dosage once
`chr_meta` ever gains non-default rows.

- Add a small shared helper (e.g. `assert_diploid_only(pop)` in `R/sql_utils.R` or a new
  `R/ploidy_helpers.R`) that runs:
  ```sql
  SELECT COUNT(*) AS n FROM chr_meta WHERE copies_M != 2 OR copies_F != 2
  ```
  and errors with a clear message (mentioning Stage 4 / `define_chr()`) if `n > 0`.
- Call it from `add_dosage()` and `extract_genotypes()` before the dosage/export SQL runs.
- **Not** added to `add_tbv()`/`add_phenotype()` in this stage — those are out of this stage's
  scope; ploidy-awareness across the whole codebase is Stage 4's job. Scoping the guard to the two
  functions this stage actually touches avoids scope creep.
- Since `chr_meta` today only ever has the default autosome rows (no `define_chr()` exists yet),
  this check is unreachable in current usage — it's forward defense for the gap between "Stage 3
  ships" and "Stage 4 ships," not a fix for an active bug.

---

## Remaining exit-criteria tests to add

The parent plan's Stage 3 exit criteria and the acceptance-criteria table above imply tests that
don't exist in `tests/testthat/test-add_dosage.R` / `tests/testthat/test-extract_genotypes.R`
today:

1. **Cache vs. direct extraction parity** — for the same individual/locus set, `add_dosage()` into
   `ind_genotype` followed by reading it back must equal `extract_genotypes()`'s direct-from-
   `ind_haplotype` matrix.
2. **Partial-cache safety** — populate `ind_genotype` for only a subset of individuals/loci via
   `add_dosage()`, then confirm `add_tbv()` and `add_phenotype()` results are unaffected and
   identical to a run where `ind_genotype` was never populated at all (i.e. confirm by test, not
   just by code inspection, that TBV/phenotype never read `ind_genotype`).
3. **`loci_tbl` coverage** — no existing test exercises `extract_genotypes(loci_tbl = ...)` alone
   or unioned with `chip_name`/`effects_tbl`; add tests mirroring the two usage examples in the
   parent plan (autosomes-only extraction; chip ∩ autosomes).
4. **Injection-safety regression tests** — a `chip_name`/`locus_names`/`id_ind`-adjacent value
   containing a quote or SQL keyword is rejected with a clear error, not executed.
5. **Diploid guard** — a `chr_meta` row with non-default `copies_M`/`copies_F` causes `add_dosage()`
   and `extract_genotypes()` to error before touching data (test by directly mutating `chr_meta`,
   since `define_chr()` doesn't exist yet).

---

## Scope for Stage 3

- Add `validate_sql_identifier()` calls for `chip_name` in `add_dosage()` and `extract_genotypes()`.
- Decide and implement a fix for quoted-value interpolation (`locus_names`, `id_ind` sets) —
  parameterized queries preferred; single-quote rejection as a fallback if parameterization proves
  too invasive for this stage.
- Add `assert_diploid_only()` (or equivalent) and call it from `add_dosage()` and
  `extract_genotypes()`.
- Add the five test groups above.
- Update `ind_genotype`'s roxygen/schema docs (already partially done in `CLAUDE.md`) to explicitly
  state partial-cache semantics if not already clear at the function-doc level.
- No changes to `add_tbv()`, `add_phenotype()`, `PIVOT`-based BLUPF90 export helpers, or new
  functions — those are unaffected by this stage's findings.

## Out of scope

- Ploidy-aware dosage/export math itself (Stage 4).
- Any change to `add_tbv()`/`add_phenotype()` beyond the parity/safety tests above (they already
  read only from `ind_haplotype`, per Stage 2).
- Revisiting `PIVOT` usage — current `extract_genotypes()` doesn't use it, so the parent plan's
  `PIVOT`-specific caveats (column-order stability) don't apply to the shipped implementation.

## Exit criteria

- All five new test groups pass.
- `grep` for un-validated `chip_name` interpolation and unparameterized quoted-value interpolation
  in `R/add_dosage.R` / `R/extract_genotypes.R` comes back clean.
- Full test suite: 0 failures, 0 errors.
- Parity golden (`tests/testthat/test-parity.R`) unaffected — this stage touches no TBV/phenotype/
  haplotype-writing code path, so no re-baseline expected.

## Deviations from the parent plan (intentional, recorded here for Stage 4+ readers)

- `add_dosage()` and `extract_genotypes(loci_tbl = ...)` were built in Stage 1, not Stage 3 — see
  "Why this stage looks different" above. This doc supersedes the parent plan's Stage 3 section as
  the sequencing authority for what remains to be done.
- See "Deviations from this plan" at the top of this doc for the `chip_name`-validation and
  escaping-vs-parameterization decisions made during implementation.

---

## What's next (not done in Stage 3)

- **Stage 4** — `define_chr()`, sex chromosomes, organelles, polyploidy. `assert_diploid_only()`
  (this stage) will need to be relaxed/replaced once `chr_meta` legitimately gains non-default rows.
- **Stage 5** — sparse mutation storage (`add_mutation()`).
- The related-but-out-of-scope observation from planning: `add_tbv()` and `add_ebv()` validate
  `trait_name`/`model` as SQL identifiers but then also interpolate the same value unescaped inside
  a quoted string literal elsewhere. Not touched here (out of this stage's scope), but worth a
  dedicated look if injection hardening is ever revisited more broadly.
