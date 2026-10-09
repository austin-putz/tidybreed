# Import QTL-effect methods — Phase 5 Codex review

Re-reviewed 2026-10-06 against the revised `plans/import_qtl_effect_methods_phase_5_plan.md` (926 lines), tidybreed 0.75.2, the main specification, `CLAUDE.md`, relevant production code/tests, and source commit `318e54f`.

**Assessment: the revision addresses the original eight findings at the design level. Keep the revised design; fix the remaining scalar-solver cases, settle the floor's coordinates, and tighten the validation/test details before 5b.** The writer work in 5a now has a substantially stronger contract. None of the remaining calibration findings requires redesigning it.

These are findings about a plan, not defects verified in a Phase 5 implementation. The proposed generator and new solver do not exist yet.

## Original findings: status after the revision

| Original finding | Status in the revised plan |
|---|---|
| 1. Singular additive targets are not mathematically incompatible with dominance | Addressed in D5: implementation restriction, zero-block exception, and G1's feasible counterexample. |
| 2. The source dominance-mean solver cannot be ported verbatim | Addressed in direction: anchor-based forms, linear/degenerate branches, signed-ratio verification and G4. The newly specified arithmetic still has finding 1 below. |
| 3. Additive-only identity needs shared storage, including zeros | Addressed in 5b.4–5b.5 and C4: shared draw/calibration/build path, zero rows retained. |
| 4. Depression subsets and pre-draw zero-dominance checks need a contract | Addressed in D4, 5b.2–5b.4 and G5: align by name, solve requested columns only, supply `z`, no multivariate closeness guarantee. |
| 5. Floor tolerance, independent trait units and rank convention need gates | Substantially addressed: residual floor, operand-relative rounding budget, G6, cached SVD rank. Clarify coordinates as in finding 2 below. |
| 6. Exact HWE/LE counts do not eliminate the sample-variance denominator | Addressed in 5c.1 with the `n/(n−1)` factor. |
| 7. Resource limits must follow the active path; knowable pair shortages precede draws | Addressed in principle: additive-only parity, retained-cell counts, chunked diagnostics, structural shortages. Remaining validation/G7 details are below. |
| 8. Writer error order needs a full rule and explicit tests | Addressed in D7 and 5a.3: term/rule/member order, whole violation vector, and a contract test before refactoring. |

The smaller recommendations were also incorporated: per-block target precedence (D8), mid-transaction target-write injection, populated-model benchmarks, forward-conversion profiling, sampling-parameter summaries, and oracle isolation. The oracle dependency list still needs completion.

## 1. High — The new dominance-mean formula still rejects feasible cases

**Plan locations:** 5b.4, `.na_solve_dd_mean()`, especially lines 495–506; G4.

### The stable quadratic formula needs a nonzero sign at B = 0

The plan specifies:

```r
q = -(B + sign(B) * sqrt(disc)) / 2
roots = c(q / A, C / q)
```

In R, `sign(0) = 0`. A quadratic with `B = 0` and two nonzero roots therefore gets `q = 0`, losing both roots.

This is a feasible **two-locus genic** example:

```text
p = (0.5, 0.5)
u = (1, 1), v = (1, -1), w = (0.5, 0.5)
M_D = 0.25 I, sd = 0.1, rho = 1

A = -0.5, B = 0, C = 0.005, discriminant = 0.01
true roots = ±0.1
```

At `mu = 0.1`, `d = (0.2, 0)`, so `ID = 0.1`, `V_D = 0.01`, and the requested ratio is exactly 1. The plan's formula returns `0, -Inf`. Final verification should prevent wrong output, but would reject an attainable request. I reproduced these values in R.

**Revision:** use a sign convention with `sign_plus(0) = +1`, or an explicit `B = 0` branch. Handle `q = 0` separately for a repeated zero root; do not divide `C/q` blindly. Retain finite-root filtering and signed-ratio verification.

### Reflection does not move a preferred mean already at ID = 0

The all-zero-polynomial branch chooses `dominance_degree_mean` if its sign fits, otherwise reflects it across the root of `ID(mu) = 0`. If the preferred mean is that root, reflection leaves it unchanged.

For one genic locus, take `u = 1`, `v = -1.9`, `sd = 0.1`, `w = 0.5`, `M_D = 0.25`, `rho = 1`, and preferred mean `0.19`. All polynomial coefficients are zero. The preferred mean and its reflection are both `0.19`, giving **`ID = V_D = 0`**, where the ratio is undefined. A finite mean greater than `0.19` attains the request; a negative request needs a mean below the root. I reproduced the zero-variance candidate in R.

**Revision:** every accepted candidate needs finite **positive** `V_D` as well as the signed ratio. If reflection leaves the candidate on the zero, choose a finite displacement on the required side and verify it. Apply this positive-variance condition to every root branch, including `rho = 0`; zero ID alone does not make a zero architecture calibratable to positive `G_D`.

**Numerical clarification:** the stated coefficient scale `s` omits the `sd` and `sd²` factors in `B` and `C`. Very small or large valid `sd` can change branch classification for the same dimensionless problem. Since depression requests require `sd > 0`, solving in `x = mu/sd` removes those factors before classification. Define the coefficient/discriminant rounding budgets explicitly.

**Gates:** add both hand-derived examples to G4, a repeated-zero-root case with positive remaining variance, and a small/large-`sd` sweep of a supplied architecture. Assert finite positive dominance variance and the signed ratio. Preserve the existing one-locus ±1 and feasible-linear probes.

## 2. Medium — Make the additive floor's coordinates unambiguous

**Plan location:** 5b.4, `.na_additive_stage()`, lines 517–530.

The first bullet moves both architectures to correlation coordinates: `B_a S⁻¹` and `C S⁻¹`. The residual-floor bullet defines `floor = M.cov(C_res)`, but the following bullet subtracts `S⁻¹ floor S⁻¹` from `R_A`.

Those definitions are consistent only if the residual floor uses the **original-unit** coefficients. If computed from the standardised coefficients, it is already on the correlation scale and must not be divided by `S` again. The current text leaves both readings available.

For `G_A = 4` and an original-unit floor of 1, the correct standardised residual target is `1 − 1/4 = 0.75`. Standardising the floor twice gives `1 − 1/16 = 0.9375`. Final verification catches the resulting miss, but a port of that reading would fail feasible calls.

**Revision:** name the coordinates explicitly, for example:

```text
B_s = B_a S⁻¹; C_s = C S⁻¹
P_s = M.cov(B_s); Q_s = M.cross(B_s, C_s)
C_res_s = C_s − B_s solve(P_s, Q_s)
floor_s = M.cov(C_res_s)
G_tilde_s = R_A − floor_s
floor_original = S floor_s S       # messages only
```

Apply the rounding budget to `floor_s` and rescale final coefficients once. Computing the residual in original units is also valid; explicitly choose that alternative if intended.

**Gate:** extend G6 with a non-unit scalar variance and supplied coupling with a hand-computed nonzero floor. Assert the reported floor in original units and final covariance, alongside the existing per-trait-unit tests.

## 3. Medium — Finish pre-draw validation and its dependencies

**Plan location:** 5b.2 steps 1, 3, 5 and 6; C20/G7.

The plan allows `dominance_degree_mean = 0` and `dominance_degree_sd = 0`. With nonzero resolved `G_D` and no depression request, this makes `B_d = 0` for **every** additive draw. The eventual `.qtl_calibrate()` architecture-rank error is knowable from the inputs; consuming A/D draws first would violate the “every knowable refusal” contract.

**Revision:** before drawing, refuse nonzero dominance targets when both degree parameters are zero. Keep that combination valid for absent/exactly zero dominance targets. Explain it as an architecture-parameter restriction.

The numbered ordering also has an unresolved dependency: step 3 uses `m`, the selected locus set, and under the realised anchor `n`, while locus/base resolution supplying them is step 5.

**Revision:** resolve sorted locus keys/count and lightweight base information before pair/count checks that need them. Resolve the realised individual's count and restrictions before `rank(G_c) ≤ n−1`; apply the cell guard before collecting dosage/design arrays. Separate early pair syntax checks from checks needing resolved keys/counts. Keep supplied-pair anchor-rank checks after design construction.

**Gates:** positive `G_D` with both degree parameters zero refuses without changing database/RNG, with and without an existing seed. Explicit zero `G_D` with those parameters succeeds via the additive-only route. Retain pair-shortage and cohort-rank pre-draw gates.

## 4. Medium — G7 cannot require successful sampling to preserve RNG

**Plan location:** 5b.9, G7, lines 764–769.

G7 requires an additive-only realised call just below Part A's limit to **succeed**, then says “all of these preserve `.Random.seed`.” Successful positive-target sampling advances RNG. This contradicts 5b.3, additive-only parity, and the intended failure contract. A focused positive-target additive calibration probe confirms the seed changes.

**Revision:** scope preservation to G7's deterministic **refusals**. For the successful boundary case, assert the same rows and post-call RNG state as `define_additive_effects()` under the same starting seed/inputs. For validation starting without a seed, assert that `.Random.seed` remains absent; do not read an undefined variable.

The boundary fixture can temporarily lower the cell-limit constant. There is no need to allocate almost 20 million cells to test dispatch/guard parity. Keep large memory measurements in the benchmark.

## 5. Medium — The source-oracle list omits required dependencies

**Plan location:** 5b.9, “Oracle isolation,” lines 674–685.

The list includes `sim_qtl_effects_nonadd()` but omits `nonadd_decompose()`, which it calls unconditionally for founder diagnostics. It also needs source `.na_*` validation/target/spectrum helpers. Copying only the six listed definitions leaves the oracle unusable, or permits lookup outside the intended environment.

**Revision:** copy the transitive dependencies needed by the entry points exercised, including `nonadd_decompose()` and relevant source helpers. Keep them in an isolated environment with an explicit parent providing base/recommended-package functions, rather than the package namespace or global environment. Alternatively, exercise fewer source entry points and copy their complete dependency set.

**Gate:** smoke-test the isolated oracle on a small A + D + A×A call without sourcing the checkout into `.GlobalEnv`. Confirm every source-owned function resolves to the oracle environment. The supplied-architecture and hand-derived cases remain the substantive correctness gates.

## Smaller improvements and clarifications

- **Pin D8 with explicit new tests.** Main-spec C17 does not cover all new combinations. Add passed D + explicitly selected stored A; passed A + stored D/A×A; a passed block duplicated in the explicit selection; an overwrite hidden by filtering; and an intentionally filtered-out stored block. Preserve the current resolver's scope check: explicit line rows cannot describe newly generated common-scope terms. `.dae_resolve_target()` already enforces this; the per-block extraction must retain it for A, D and A×A.
- **Correct 5b.5's nonempty-model explanation.** Positive-definite `G_A` guarantees nonzero calibrated `B_alpha`, but stored HWE `stat$alpha` differs under the realised anchor and need not itself be nonzero. The model is still nonempty: if all stored α, d and e were zero, the realised additive component would also be zero. Use that argument rather than identity of observed/HWE α.
- **Describe the floor as conditional on the sampled architecture.** Its diagonal is the minimum reachable within the chosen additive span and calibrated non-additive effects, not a universal biological minimum determined by `G_D`/`G_AA` alone. Retain this qualification in errors/messages/roxygen. With `C = 0`, the floor can be zero even with nonzero dominance/epistasis.
- **Keep D3 gates specific to the relevant block.** A fixed locus's centred contrast has zero base variance; a pair containing it can still induce a partner's additive effect through `e*c`. The revision correctly retains that coupling. Do not require the entire model's variance to vanish because one locus/pair is fixed.

## Decisions D1–D8

| Decision | Review recommendation |
|---|---|
| D1 — three commits | Agree; measure and verify 5a independently. |
| D2 — A, D, pairs, A×A draw order | Agree; distinguish draws from final coefficients. |
| D3 — retain effects at fixed loci | Agree with (a), including induced partner effects. |
| D4 — exact single-trait depression; report multivariate delivery | Agree with (a), subject to scalar-solver corrections. |
| D5 — positive-definite A for nonzero non-additive targets | Accept the release limitation and zero-block exception. |
| D6 — diagnostics by block/base kind | Agree, including chunking, pre-commit computation and post-commit reporting. |
| D7 — preserve writer rows/messages/error order | Agree with explicit contract tests and full-table validation. |
| D8 — passed/stored sources per block | Agree with (a); add mixed-source and scope gates. |

These are review recommendations, not user approval of the pending decisions.

## Verification performed and limits

I reread the revised plan and compared its changed contracts with the previous review, main specification, `CLAUDE.md`, current additive target resolver, congruence/anchor helpers, conversion/builders, and source calibration implementation/tests. The package remains at 0.75.2; no Phase 5 generator implementation is present.

I ran independent **pure-R probes of the newly specified arithmetic**. They reproduce the `B = 0` root failure, the degenerate reflection's zero-variance candidate, the floor-coordinate distinction, and successful sampling's RNG advancement. I also inspected source-oracle dependencies. Script: `/tmp/tidybreed_phase5_rereview_probes.R`; output: `/tmp/tidybreed_phase5_rereview_probes.log`. These do not claim to test an implemented `.na_solve_dd_mean()`.

The initial source/current-helper probes remain in `/tmp/tidybreed_phase5_codex_review_probes.R` and its `.log`; their findings are now marked addressed in the plan rather than repeated as outstanding bugs.

This is a documentation-only re-review. I did not run the full package suite, rebuild docs, or certify performance targets. Only this review file changed; the implementation plan was left untouched.
