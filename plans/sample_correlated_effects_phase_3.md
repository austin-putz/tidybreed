# Correlated effects — Phase 3 results

**Spec:** `plans/sample_correlated_effects.md` (v3.3), §5.2, §5.3, §5.4, §7
"Numerical" and §8 Phase 3.
**Branch:** `feat/genome-effects-v49`.
**Follows:** `plans/sample_correlated_effects_phase_2.md`.
**Status:** complete. Full suite green.
**Date:** 2026-09-20.

Phase 3 is the mathematical core of the feature, delivered with **no caller
change**: a loader that turns the rows Phase 2 guarantees into matrices, and a
pure resolver that draws the coordinates of a Gaussian block conditional on
whatever each entity has already realized. `add_phenotype()` still runs its old
§7.5/§8.5 sampling; Phases 5 and 6 swap it for adapters over these two
functions. Because nothing calls them yet, everything here was testable in
isolation — the resolver without a database at all — which was the point of
sequencing it before the record-planning extraction.

---

## What shipped

| File | Change |
|---|---|
| `R/correlated_draws.R` | **new** (~430 lines, half of it roxygen). `find_covariance_blocks(conn, effect_name, phenotype_names)` — components over pair-row existence (union–find over `.pvc_block_members()` + `.pvc_block_rows()` from Phase 2), one matrix per stratum, D1 invariants re-checked on load. `resolve_correlated_draws(covariance, sample_coordinates, entity_keys, observed = NULL, tolerance = NULL)` — per-pattern conditional moments, Cholesky / eigen-pseudoinverse with support consistency, relative tolerance, all checks before the first `rnorm()`, exactly `n × m` normals consumed. Helper `.cd_factor()` (Cholesky when PD, sign-normalized `V √D` when PSD) and `.cd_build_block()` |
| `R/phenotype_cov_block.R` | `.pvc_write_block()` stores `(M + t(M)) / 2` — bit-identical for an already-symmetric input — so the stored `(i, j)` and `(j, i)` rows are exactly equal and the loader can demand it |
| `tests/testthat/test-correlated_draws.R` | **new**, 115 expectations; see below |
| `man/` | `correlated_draws.Rd`, `find_covariance_blocks.Rd`, `resolve_correlated_draws.Rd` (all `@keywords internal`) |
| `plans/sample_correlated_effects.md` | v3.3 header, §5.2 and §5.4 marked shipped with the as-built contract, adapter steps 6–7 revised, Stage 2 pseudocode and worked example renamed, §7 coverage note, §8 Phase 3 |
| `CLAUDE.md`, `NEWS.md` | `phenotype_var_comp` covariance-block paragraph gains the loader/resolver pointer and the exact-symmetry note; 0.71.0 "Changed" entry |

`NAMESPACE` is unchanged: both functions are internal.

---

## What a user sees differently

Nothing. `add_phenotype()`, the three writers' accepted inputs, and every
seeded output are exactly as at the end of Phase 2. The one observable
difference is inside `phenotype_var_comp`: a matrix that was symmetric only
within the writer's tolerance (say `R[A,B] = 0.3 + 1e-12`, `R[B,A] = 0.3`) is
now stored as its symmetric part in both rows.

---

## Why the shape is what it is

**The loader is plural and takes a connection.** The plan wrote
`find_covariance_block(pop, effect_name, target_phenotypes)` returning "the
connected component containing the targets". One `add_phenotype(c("A", "B",
"C"))` call routinely spans components — `{A, B}` correlated, `C` on its own —
so the adapter needs the partition, not one component. It returns a list, one
entry per component touching the targets, ordered by first member; phenotypes
with no rows are simply absent (the caller's `setdiff` finds them, and for the
residual effect that is still "No residual variance found"). It takes `conn`
because its siblings in `phenotype_cov_block.R` do and the adapters will hold
`pop$db_conn`.

**Every stratum comes back as a matrix, already checked.** Phase 5's residual
adapter needs "the matrix for stratum *s*" per entity; making it assemble
matrices from rows would put the D1 completeness logic in a second place.
The loader builds `unconditional` and `conditional[[level]]` in sorted member
order and refuses anything that violates what the writers guarantee: a stratum
with fewer than `n²` rows (names the missing pairs), duplicate rows, a
non-finite `cov_value`, an asymmetric pair, or two condition columns on one
block. These states are reachable — the D1 "growing" recipe and the D3 lock
recipe both send users to `remove_rows()` on `phenotype_var_comp`, and a
partial `filter()` leaves exactly such a hole — so the error carries the
redeclaration call rather than guessing around it.

**Pattern grouping moved into the resolver.** The plan split the work as
"adapter groups by `(stratum, observed, sample)`; resolver called once per
pattern". Grouping by observed pattern requires looking at the `NA` structure
of a matrix the adapter has already built, so the resolver does it: the adapter
passes one `observed` matrix per `(stratum, sample set)` with `NA` for "not
observed", and the resolver splits rows by non-`NA` pattern, computes
coefficients and the factor once per pattern, and applies them to the rows.
What makes this safe is the RNG contract below — the stream cannot depend on
how the rows were grouped.

**The RNG contract is `n × m`, unconditionally.** A successful call draws
`rnorm(n * m)` once, reshaped by row so entity *i*'s normals are consecutive,
and each pattern multiplies its rows by its own factor. Two consequences were
chosen deliberately:

- A coordinate with **zero conditional variance** (perfect correlation, a
  zero-variance member) still consumes its normal; the factor's column is
  zero, so the value returned is the conditional mean exactly. The alternative
  — skipping degenerate coordinates — makes the count data-dependent and the
  §7 accounting test a re-implementation of the resolver.
- **Every check runs before the draw**: validation, the PSD check on `R`, the
  support check for every entity, and the factorization of every pattern's
  conditional covariance. A rejected call leaves `.Random.seed` untouched. (D7
  is about the *outer* `add_phenotype()` call, where a Stage-3 write failure
  legitimately follows Stage-2 draws; within the resolver there is no reason
  to spend the stream on a call that fails.)

**Tolerance is relative to `λ_max`.** The plan asked for "a function of
dimension, the largest eigenvalue, and machine epsilon". The default relative
tolerance is `nrow(R) × sqrt(.Machine$double.eps)` (≈ `1.5e-8 × n`) and the
absolute threshold is that times the largest eigenvalue of `R`. It is used
three times: to decide whether `R_oo` is PD (Cholesky) or only PSD
(pseudoinverse); to decide whether an observed vector is on the support of a
singular `R_oo` — `‖V_null' e_o‖ ≤ tolerance × max(‖e_o‖, √λ_max)`, a
sd-scale comparison; and to project eigenvalues of the conditional covariance
in `[−tol, 0)` to zero while rejecting anything below. The tests run the same
relative perturbation at variance `1e-6`, `1`, and `1e6` and get the same
accept/reject decision. `sqrt(eps)` rather than `eps` because the Schur
complement of an ill-conditioned `R_oo` carries rounding error of order
`eps × κ(R_oo) × λ_max`, and near-perfect correlation is a legitimate input.

**Cholesky when PD, eigen when PSD — with fixed signs.** Cholesky is
platform-deterministic; eigenvector signs are not (they depend on the LAPACK
build), which would make a seeded draw through the singular path differ across
machines. `.cd_factor()` therefore flips each eigenvector so its
largest-magnitude component is positive before forming `V √D`. Not a
cross-version promise (CLAUDE.md), but it removes one avoidable source of
same-code, different-machine divergence. `MASS::mvrnorm()` is not used: it
has its own eigen path and sign convention that the accounting test would have
to replay.

**Entity keys may be a data frame.** The residual entity is
`(id_ind, pheno_number)`; the natural R object for a list of those is a data
frame, whose `length()` is its column count. The resolver counts rows for a
data frame and elements otherwise, and uses nothing else about the keys.

---

## Tests (`test-correlated_draws.R`, 115 expectations)

Reference values come from the textbook conditional formulas with
`MASS::ginv()` for the pseudoinverse and an explicit replay of the normals the
resolver consumes (`matrix(rnorm(n * m), byrow = TRUE)`), so each assertion is
"the resolver equals `μ + z F` for the `z` it must have drawn", not "the
resolver equals itself".

**One-dimensional.** 1 × 1 unconditional (`√v z`); one sample given one
observed (closed-form `μ`, `σ²`); one entity, two sample coordinates given one
observed (Cholesky factor).

**Mixed patterns.** Five entities with four different observed patterns in
one call, each checked against its own conditional moments with its own `z`
row; an all-`NA` column equals not supplying it; a latent third coordinate
gives the same draw as the 2 × 2 block; a data frame for `observed`; a data
frame for `entity_keys` (rows, not columns).

**RNG contract.** Exactly `n × m` normals with mixed patterns and with none;
zero entities (`character(0)`, `NULL`, empty data frame), zero sample
coordinates, and both, are RNG-neutral and correctly shaped; a rejected call
(off-support, unknown coordinate) leaves the seed untouched; zero conditional
variance consumes its normals and returns the mean exactly.

**Singular and near-singular.** Perfect correlation determines the sample
coordinate; `ρ = 1 − 1e-6` gives `√(1 − ρ²)`; a zero-variance member is valid
both sampled (exactly 0) and observed at 0 (conditions nothing); support
consistency — off-support errors naming the count and an example entity,
`1e-12` noise on a zero-variance coordinate passes, two perfectly correlated
observed coordinates that disagree error while the same values with one of
them latent pass; the same relative perturbation accepted/rejected identically
at variance `1e-6`, `1`, `1e6`; a rank-1 matrix perturbed at `1e-13` (really
slightly indefinite) is accepted and gives the rank-1 answer, perturbed at
`0.5` is rejected as not PSD; variance scales `1e-12` and `1e12` match the
reference to `1e-10` relative.

**Distributional.** 20 000 entities, `{B, C} | A`: empirical conditional
covariance and mean within 4 sampling standard errors derived from the target
(`√((σ_ii σ_jj + σ_ij²) / n)`), not a fixed margin.

**Validation.** Missing/unequal/duplicate dimnames, non-finite and asymmetric
`R`, duplicate or unknown sample coordinates, non-character
`sample_coordinates`, row-count mismatch, unknown or overlapping `observed`
columns, `Inf` in `observed`, unnamed `observed` columns, character
`observed`, negative tolerance.

**`find_covariance_blocks()`** (a bare `open_pop()` is enough — no genome, no
phenotypes): two residual components plus an absent phenotype from one call,
members sorted although declared `{B, A}`, strata `F`/`M` as named matrices
with `condition_table`/`condition_column`, a singleton with `NULL` condition
fields and an empty `conditional` list, one member returns the whole
component, a different `effect_name` is a different graph, empty input; a
block with only conditional strata; exact symmetry of a stored matrix whose
input was symmetric only within tolerance; hand-removed rows — one pair row
(names the pair and the redeclaration call), a whole stratum (the other stratum
still loads), one phenotype's rows from a stratum (names both missing pairs
and the stratum); a forged second condition column; RNG-neutrality; and the
loader's output fed straight into the resolver equals the in-memory matrix.

---

## Verification

- `testthat::test_file("tests/testthat/test-correlated_draws.R")`: 115
  expectations, 0 failures.
- Full suite (`testthat::test_dir("tests/testthat")` after `load_all()`):
  **3167 passed, 0 failed, 0 errors, 1 skipped.** The Phase 2 total was 3042;
  114 of the +125 are this file (115 after one more loader assertion added in
  the final review). The other 11 are not from any edited test —
  `git status` shows the new file as the only test change — but from
  expectation counts elsewhere that vary between runs (assertions inside loops
  over sampled data).
- `roxygen2::roxygenise(roclets = "rd")` clean; three new internal Rd pages,
  `NAMESPACE` unchanged.

---

## Review pass (after the suite)

- **Bug found and fixed before the suite:** the nothing-observed pattern has
  key `""`, and `groups[[""]]` is `NULL` in R, so those entities were never
  filled (the first smoke test showed `NA` rows). The pattern loop now
  iterates by index. The 1 × 1 and zero-variance tests would have caught it;
  they did.
- **Wrong test, not wrong code:** an early "rejected call is RNG-neutral"
  test fed a single observed coordinate with positive variance and expected a
  support error. A single non-degenerate coordinate is always on its support;
  the support check only bites when `R_oo` itself is singular. The test now
  observes two perfectly correlated coordinates that disagree.
- Confirmed by reading: the PD branch solves `R_oo X = R_os` with two
  triangular solves and transposes (`W = R_so R_oo⁻¹`); the conditional
  covariance is `R_ss − W R_os`; `z %*% F` with `t(F) F = C` for both the
  Cholesky (`F = U`) and eigen (`F = t(V √D)`) factors; the pseudoinverse and
  the null-space projection are sign-invariant, so sign normalization is only
  needed in `.cd_factor()`.
- `lambda_max` is first `max|R_ij|` (for the symmetry check, before `R` is
  known to be PSD) and then the largest eigenvalue; for a PSD matrix the
  latter is ≥ the former, so the symmetry threshold is never looser than the
  rest.
- The all-zero matrix is handled: `tol_abs = 0`, no eigen call, every observed
  value must be exactly 0, every draw is exactly 0 (and still consumes its
  normals).
- Final review added one loader check: a conditional row whose
  `condition_level` (or `condition_table`) is `NULL` is an error rather than
  being silently excluded from every stratum by the `==` level filter.
- No caller of `load_phenotype_cov()` / `get_residual_cov()` was touched;
  `grep` confirms `find_covariance_blocks` and `resolve_correlated_draws`
  appear only in `R/correlated_draws.R`, its tests, and the docs.

---

## Plan bookkeeping

- Status header → v3.3, Phases 0–3 shipped; new "What changed from v3.2 to
  v3.3" (five bullets: plural loader, grouping inside the resolver, checks
  before RNG / `n × m`, relative tolerance + sign-normalized eigen path,
  writer symmetrization).
- §5.2 and §5.4 retitled "✅ shipped (Phase 3)", signature blocks replaced
  with the as-built ones, the resolver's contract spelled out; adapter steps
  6–7 revised; §5.3 grouping sentence annotated.
- §5.5 Stage 2 pseudocode, §5.8 worked example, and D1 consequence 2 renamed
  to `find_covariance_blocks()`; `load_phenotype_cov()`'s `NULL` meaning kept
  "until Phase 6".
- §7 "Numerical" gains the coverage note; §8 Phase 3 marked shipped.

---

## Next

**Phase 4 — extract record planning (Stage 1).** The largest single step:
pull sex expression, the repeatable guard, covariate skip, formula/composite
exclusion, path classification, `pheno_number` assignment, stratum lookup and
named-effect level collection out of `add_phenotype()`'s per-phenotype loop
into an in-memory plan with no RNG and no writes. Behaviour-preserving for the
non-correlated path; the existing suite is the check. Phase 5 then writes the
residual adapter over that plan and these two functions, and deletes
`sample_residuals()` (and with it the recorded `MASS::mvrnorm(n = 0)` crash).
