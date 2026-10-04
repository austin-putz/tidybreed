# Import the QTL-effect methods — exact genetic covariances, dominance, A×A epistasis, and true genetic values

**Status:** design proposal, nothing implemented. **For review.**
**Created:** 2026-09-22.
**Baseline:** v0.71.1 on `main` (`feat/genome-effects-v49` is merged, so §10 step 0 is done).
Line references to `R/` below were taken at 0.71.0. Re-check them before editing.
**Source project:** `~/Claude/simulate_qtl_effects/` — the GSE manuscript method
(`R/qtl_effects.R`) and its non-additive extension
(`non-additive/R/qtl_effects_nonadd.R`). The status note this plan grew out of is
`docs/notes/note_tidybreed_readiness.qmd` there. It has render-time checks of every
numerical claim below.
**Relates to:**
- `plans/update_genome_effects_v4.md` (v4.9): the term/member/origin storage this plan writes into. **No DDL change to the `genome_effect*` tables.** The DDL this plan does change: `trait_var_comp.line_name` added and `trait_meta.target_add_mean` dropped (step 2, §6C); `ind_tbv` dropped, `phenotype_meta.formula_tbv` → `formula_tgv`, and the `phenotype_components.component_names` default → `'total'` (step 3, §6A, Q18).
- `plans/update_genome_effects_base_tbl.md` §2.5: the writer/generator split this plan follows.
- `plans/consolidate_genetic_values.md`: `ind_tbv` → `ind_tgv`. **Integrated into §6 of this plan** (2026-10-02); that file is now a pointer here.
- `plans/update_genome_effects_v4_overall_summary.md` §3.2: `extract_genetic_variance()`, which this plan specifies (§8).

**Review 2026-09-22 (in-repo integration pass).** Changes from the draft:
1. Added precondition **P2** (§6A): the phenotype layer must read the *total* genetic
   value for `self`, not only `additive`. Without it Part C writes dominance and
   A×A effects that never reach a phenotype. Consolidation Task 3 covers only
   dam/sire/group contributors.
2. §7.4 / Q7: diagnostics are **not stored** (not in `trait_var_comp`, not as a pop
   attribute). `extract_genetic_variance()` can measure any population that can still be
   named. That is not necessarily the exact generation-time population (§7.4). See
   CLAUDE.md principle 6 and "never in the R object".
3. §8: `reference_tbl` → `base_tbl` (same concept, same name). Output columns match
   `trait_var_comp` so target and realised values join directly. Purely additive models
   of any scope are decomposed exactly, so they do not fall into the unassigned remainder.
4. §7.1: `anchor = "realised"` is restricted to individual-selecting `base_tbl` and the
   common scope.
5. §9: the ported calibrator takes `p`, not `X`, on the genic path. Effects are normally
   defined **before** `add_founders()`, so there are no individuals yet. Argument
   names are spelled out, and `pairs` uses `locus_name_1` / `locus_name_2`.
6. §11: gate A1 rewritten. It compared against v0.71.0 output, which CLAUDE.md forbids.
   Added determinism (`expect_identical`) and phenotype gates.
7. New open questions Q10–Q14.
8. **Decided 2026-09-23 — naming (§6B).** One vocabulary of lowercase full words:
   `additive`, `dominance`, `additive_by_additive` for variance components, and
   `additive`, `dominance`, `indicator`, `interaction` for `ind_tgv` value components.
   `gen_add`, `epistasis` and the `order1_` prefix go. The `interaction` row stays as one
   row (consolidation Task 1a dropped).

**Review 2026-09-23 (Codex critique, `plans/import_qtl_effect_methods_codex_review.md`).**
Its code claims were checked against 0.71.x. They hold: `define_additive_effects()` writes a
passed `G` (`R/define_additive_effects.R:328`) before its `.ge_commit()` (`:485`);
`define_effect_cov_matrix()` checks only symmetry and non-negative diagonals (`:115-120`);
manual `effects` is a positional numeric vector (`:256-262`). Changes made:
1. §3, §4.3.4, §6A: the two codings agree on the total **up to a trait-specific
   constant**, not exactly. The constant is written out. `phenotype_meta.mean` is an
   intercept, defined as such (§6A).
2. §8 rewritten. It has three supported cases, a genic and a realised definition, and
   cross-block covariances. `per_ind` is removed; the generation-t breeding value
   moves to a later `extract_breeding_values()`. There is a size guard.
3. §7.1: manual `effects` becomes a named table *(superseded 2026-10-03: removed, Q21)*. Supplied values are never scaled. The
   exactness promise is stated at the call. `warn_bounds` is exposed. Targets are
   validated as PSD. Target and term writes are atomic.
4. §9.1 / §9.3: `aa_terms()` is specified. There is a worked manual A+D+A×A example,
   and pair validation rules.
5. §6B rule 5: the reserved genetic `effect_name` set is a finite list, not `*_by_*`.
6. §11: A2, A7, A8, A9, B1–B8, C2, C4, C11, C12 revised; A10–A12, B9, C15, C16 and
   P1–P4 added.
7. New open questions Q15–Q18. **Four points are left for the author**: target source
   and persistence (Q15), how cross-block covariances are reported (Q16), a manual
   Cockerham path for A×A models (Q17), and the new name of `formula_tbv` (Q18).
8. **2026-09-24/25:** `additive_by_dominance` / `dominance_by_dominance` are reserved and
   refused as "not yet supported" (§6B). New §7.6 covers lines: line means come only
   from DNA via divergent selection; there is a message for line-scoped centring
   (gate A13). Per-line covariance targets are Q19. Frequency-targeting and
   cross-database seeding are deferred to `plans/TODO.md`. Q20 separates generation
   targets from evaluation parameters (`add_ebv()`).

**Review 2026-10-01 (readiness pass).** Each Codex finding was checked against the plan
again. All seven findings and every "other correction" are answered in the text. The open
ones are open questions the review asked the author to settle (Q15, Q16, Q18), not gaps.
This pass found internal inconsistencies, which are now fixed:
1. §6 amendment C5 (now K5) contradicted §4.3.1. Under `"realised"` the extractor re-projects with the observed
   within-locus regression, while the stored `additive` uses $b=q-p$. The two agree only on
   an HWE base. C5 is restated.
2. §8 case 1 did not include `ad_terms(coding = "functional")`'s heterozygote `indicator`
   term, so the manual example of §9.3 would have fallen into case 3. That contradicted Q17.
   Case 1 now names the shapes it accepts.
3. §8 `anchor = "genic"` has no definition for cases 2 and 3, because the evaluator gives
   only realised values. It is now restricted, and the role of `base_tbl` under
   `"realised"` is stated.
4. §9.3: `ad_terms()` requires `p` under both codings, because it reports $\mu$. The worked
   example omitted `p` and would have errored. `aa_terms()` now mirrors this: `p_1` / `p_2`
   are required, and the pair's $\mu$ share is reported. The pair order inside a term is
   canonical.
5. §9.1: `inbreeding_depression` needs a dominance block and is a call argument, not a
   stored target. A duplicated bullet was removed.
6. §11: a gate for target semantics (C17, written once Q15 is decided), a dimnames gate
   (C18), and a reopen gate for the Part A realised anchor (A14) were added. C13/C14 were
   moved back under Part C.
7. Q19 was moved from "schedule later" to a decision Part A needs: under Q15 (a), Part A
   would persist a line-scoped `G` into a table with no `line_name`.
8. **Decided 2026-10-01 — targets (new §6C, closes Q15 and Q19).** `trait_var_comp` is
   the single source of truth. A passed matrix (or a number, as 1×1) is written, and an
   existing block is a hard error, even an identical one. `trait_var_comp_tbl` picks stored
   rows, and the surviving `effect_name`s decide the model. `define_trait()` loses
   `target_add_var` and `target_add_mean`. `trait_var_comp` gains `line_name`, with fallback
   to `NULL`. The variance half of the §6B rename moves ahead of Part A (step 1, §10).
9. **Decided 2026-10-01 — function names (§0A, closes Q1).** The exported low-level
   writer `define_genome_effects()` becomes **`define_genome_effect_terms()`**. The freed
   name goes to the Part C generator, which was drafted as `define_nonadditive_effects()`.
   That name suggested the generator writes *only* non-additive effects, but it writes the
   additive block too, and it now also accepts an additive-only target. The writer rename
   is its **own first release** (step 1, §10), with no change in behaviour, so the name is
   never ambiguous inside one version. Below, `define_genome_effect_terms()` always means
   the writer and `define_genome_effects()` always means the new generator.
10. **Decided 2026-10-01 — Q16 and Q3.** Realised cross-block covariances go into one
    `between_components` row per trait pair (§8). Both generators share one reserved owner,
    `generated`, with the replacement rules of §5. No mutual exclusion, and no clearing
    call is needed.
11. **Decided 2026-10-02 — Q18.** `formula_tbv` → `formula_tgv` (argument and column),
    defaulting to total genetic values. The DSL gains a named-only `component =`, and its
    existing `table =` becomes named-only, documented, validated and tested.
12. **Decided 2026-10-02 — Q8.** Random A×A pairs are a random matching (each locus at
    most once). `n_pairs = NULL` means $\lfloor m/2\rfloor$, with a message. Supplied
    `pairs` may reuse loci.
13. **Final review 2026-10-02** (four parallel reviews: phase order, code claims, the
    consolidation plan, internal consistency). §10 is rewritten as a per-step checklist
    with a fixed order 1 → 2 → 3 → 4 → 5; steps 3 and 4 no longer swap. Changes:
    - **Code facts:** `trait_var_comp` stores 7 significant digits (fixed in step 2, A18);
      `ind_tgv_total` is a parallel `SUM()` (exact in step 3, PH5); the simple-phenotype
      precheck and the prevalence threshold change under P2 (PH6, PH7);
      `load_trait_cov()` / `get_trait_var()` get line resolution.
    - **Decisions:** `G` with manual `effects` or `scale_to_target = FALSE` errors (A19);
      `G_A` / `G_D` / `G_AA` kept as a stated exception; `define_trait_simple()` removed;
      `define_effect_cov_matrix(trait_names =)` → `trait_name =` in step 1.
    - **Moves and renames:** `aa_terms()` and the NOIA conversion move to step 4. The §6
      amendments are relabelled K1–K8 and the phenotype gates PH1–PH8, to stop label
      clashes. Amendment K5 becomes gate C19, with corrected conditions. A checklist for
      the consolidation plan's stale text is in §6. Gate C8 (b) is corrected (the filter is
      still needed). Gates A18–A22, B11, C19, C20 and PH5–PH8 were added.
14. **2026-10-02 — consolidation integrated.** `plans/consolidate_genetic_values.md` is
    merged into §6 (precondition P1), so this plan is the single source for all five
    steps. The K1–K8 amendment labels are gone. Its gates are now T1–T9 (§11), and the
    old file is a pointer to §6.
15. **2026-10-03 — Codex review of steps 0b and 1**
    (`plans/import_qtl_effect_methods_phase_0b_1_review.md`, response at its end).
    - Finding 2 fixed in 0.72.3: `restore_pop()` refuses the pre-0.72.0 stored strings.
    - Finding 3 fixed in 0.72.3: a `phenotype_components` group-contributor
      determinism test, checked to fail on the old floating `SUM()`.
    - Finding 1 is a design gap, not a step-0b defect. A stored `additive` target does
      not prove the active terms were calibrated to it. The §6A active-block rule
      could not tell the difference either. Recorded as **Q21**, and **decided the same
      day: (a), the owner is the provenance.** `define_additive_effects()` loses
      `effects` and `scale_to_target` in step 3. Exact coefficients go through
      `define_genome_effect_terms()`. The §7.1 named-table `effects` decision and gate
      A12 are withdrawn. §0A, §6A, §7.1, Q2, steps 2–3 and gates A19/PH7 are updated.
    - A full re-read afterwards (same day) fixed the leftovers before step 2:
      - §6C's option list;
      - A19's message;
      - target resolution on the surviving manual paths;
      - one internal target insert, so the Q21 refusal never blocks `G =`;
      - `trait_var_comp.line_name` goes in the base DDL, and `restore_pop()` keys on it;
      - `seed` moves after validation (A22);
      - `_pkgdown.yml`, the roxygen example and the stale error strings;
      - the Q21 refusal is scoped by line;
      - the `ad_terms()` coding for step 3.
16. **2026-10-03 — step 2 built (0.73.0).** Deviations in Step 2's "As built" paragraph
    and `_phase_2.md`. New open question Q22 (default `warn_bounds` on small pools).
17. **2026-10-03 — Q22 decided (b), built 0.73.1.** The founder-pool comparison is a
    message; only observed individuals and the genic limit warn (§7.4, Q22).

---

## 0. Summary

Three things are imported, after a rename-only release (§0A, §10 step 1):

| Part | What | Where in tidybreed | Schema change |
|---|---|---|---|
| **A** | Exact multi-trait additive covariance (`anchor =`); exact for sampled effects under `method = "shared"`, approximate under `"union"` (§7.3) | `define_additive_effects()`, `define_effect_cov_matrix()`, `define_trait()` (§6C) | `trait_var_comp.line_name` added; `trait_meta.target_add_mean` dropped (§6C) |
| **B** | Measuring realised genetic variance components | new `extract_genetic_variance()` | none |
| **C** | Calibrated additive + dominance + A×A epistasis (additive-only allowed) | new generator `define_genome_effects()` (the name is freed in step 1) | none |

Plus two **preconditions**, both before Part C:
- **P1**: consolidation, §6 (delete `ind_tbv`; the breeding value lives in `ind_tgv`).
- **P2**: the phenotype's own genetic value is the total, not the additive component (§6A).

The decision this plan asks for, and the reason it exists as more than a port:

> **Generated effects are always written in statistical (Cockerham) coding at the base
> allele frequencies.** The total genetic value differs between the two codings only by
> a trait-specific constant, so every difference between individuals is the same. Only
> under this coding, though, is `ind_tgv`'s `additive` component the true breeding
> value, and `dominance` / `interaction` the dominance and A×A values. Under
> functional coding `additive` correlates only 0.81–0.87 with the breeding value in
> a realistic example (§4.2).

`ind_tbv` is **not** widened to hold dominance and epistasis. A breeding value is additive
by definition. `ind_tgv` is already long by `component_name` for exactly this purpose.

---

## 0A. Function names *(decided 2026-10-01)*

| Function | Role | Status |
|---|---|---|
| `define_genome_effect_terms()` | **Writer.** Writes exact, user-supplied terms (`ad_terms()`, `aa_terms()`, `genotype_terms()`) unchanged: no sampling, no calibration. Any owner, scope or contrast | renamed from `define_genome_effects()` in step 1 |
| `define_additive_effects()` | **Generator**, additive only. Samples effects and calibrates them to the `additive` target. Has the additive-only options: `line_name`, `parent_origin`, `method = "union"`, `distribution`. Always calibrates: `effects` and `scale_to_target` are removed in step 3 (Q21) | exists; Part A extends it |
| `define_genome_effects()` | **Generator**, additive + dominance + A×A, common scope. Samples and calibrates to whichever of the three blocks the targets contain. An additive-only target is accepted and gives the same coefficients as `define_additive_effects()` with the same seed and options (gate C4) | new in Part C |
| `define_effect_cov_matrix()` | Stores targets in `trait_var_comp` (§6C) | exists |

Why these names:
- **`define_genome_effect_terms()`.** The rows it writes appear in the existing
  `genome_effect_terms` view. Its main argument is already `terms`. Its builders are the
  `*_terms()` family in `R/genome_effect_terms_builders.R`. The name says what sets it
  apart: you supply the terms.
- **`define_genome_effects()`.** It writes into the `genome_effect*` tables like the other
  two, and it is the general generator, with `define_additive_effects()` as its additive
  special case. It does not yet do line or origin scoping, A×D or D×D; the roxygen says so,
  and Q14 adds scoping later.
- Both generators call the writer's **internal** engine (`.ge_build()` / `.ge_commit()`), not
  the exported writer. The rename does not touch that path. The generator gets a new file
  in step 5, `R/define_genome_effects.R`; at 0.71.1 that path holds the writer, which step 1
  moves.
- **Target arguments `G`, `G_A`, `G_D`, `G_AA`** (decided 2026-10-02). They keep the
  matrix notation of the literature and the source method, and pair with
  `define_additive_effects(G =)`. This is a deliberate exception: §6B's full-word rule governs
  values stored in rows, not argument names. The roxygen spells each out ("`G_D`: the
  dominance covariance target, `effect_name = "dominance"`").
- **`define_effect_cov_matrix(trait_names =)` becomes `trait_name =`** in step 1, matching
  the generators and naming rule 1 (the columns are `trait_name_1` / `trait_name_2`).

**The rename rule (pre-1.0, CLAUDE.md).** In step 1 the old name goes **completely**:
the R file becomes `R/define_genome_effect_terms.R`, plus `NAMESPACE`, the man page, every
roxygen cross-reference (`[define_genome_effects()]` in other functions' docs), error and
message strings (e.g. `R/define_additive_effects.R:705`), tests, `dev/` scripts,
`package_summary.md`, `CLAUDE.md` and the `tidybreed-api` / `tidybreed-schema` skills. There is
no alias. Past `NEWS.md` entries are history and stay as written; the new entry records the
rename. Closed plans in `plans/` are history too. The active ones
(`TODO.md` and this plan) use the new name. The name stays
**unused** from step 1 until Part C (step 5). In between, an old call fails with "could not
find function", and in Part C an old-style call (`pop`, `terms`) fails on its arguments.
Neither fails silently. Gate R1 (§11) checks the rename with `grep`.

---

## 1. Background — what is being imported

### 1.1 The additive method (manuscript v4, frozen)

Given founder genotypes and a target $k\times k$ additive covariance $\mathbf G$, find the
$m\times k$ effect matrix $\mathbf B$ with

$$\mathbf B^\top \mathbf M \mathbf B = \mathbf G$$

**exactly**. Here $\mathbf M$ is an explicitly named reference-population genotype
covariance, the **anchor**:

| `anchor` | $\mathbf M$ | Names | Use |
|---|---|---|---|
| `"genic"` (default) | $\mathrm{diag}(n_{\text{eligible},j}\,p_jq_j)$ | the panmictic (HWE + LE) limit at base frequencies | multi-generation studies. Theorem 2: under random mating the realised covariance converges to this |
| `"realised"` | $\mathbf S_0=\mathrm{Cov}(\mathbf X)$ of the base individuals | the base cohort exactly as supplied, LD included | single-generation / clonal studies |
| `"reference"` | user-named weights or genotype panel | whatever the user names | a separate reference population |

Algorithm (Proposition 2): sample any architecture $\mathbf B_0$ ($m\times k$), then
$\mathbf B=\mathbf B_0\mathbf A$ with
$\mathbf A = \mathbf U_c\mathbf D_c^{-1/2}\mathbf Q\,\boldsymbol\Gamma^{1/2}\mathbf V^\top$
($\mathbf C=\mathbf B_0^\top\mathbf M\mathbf B_0=\mathbf U_c\mathbf D_c\mathbf U_c^\top$,
$\mathbf G=\mathbf V\boldsymbol\Gamma\mathbf V^\top$, $\mathbf Q$ the polar factor of a thin
SVD of $\mathbf U_c^\top\mathbf V$). Cost $O(nmk)$, never an $m\times m$ matrix. $\mathbf Q$ is
what makes rank-deficient architectures correct. Feasible **iff**
$\mathrm{rank}(\mathbf G)\le\mathrm{rank}(\mathbf M)$ (Theorem 1). That is an existence
claim for *some* effects. On a fixed drawn $\mathbf B_0$ the congruence also needs
$\mathrm{rank}(\mathbf C)\ge\mathrm{rank}(\mathbf G)$, since $\mathbf B^\top\mathbf M\mathbf B=\mathbf A^\top\mathbf C\mathbf A$.
A sampled $\mathbf B_0$ meets this almost surely when the anchor does, but a supplied or
degenerate architecture may not. The two failures get distinct errors (gate A5). For $k=1$ this is exactly
today's scalar rescale.

### 1.2 Why tidybreed needs it

`define_additive_effects()` draws multi-trait rows from `MASS::mvrnorm(0, G)` and rescales
each trait by a scalar. A diagonal rescale cannot change a correlation, so the delivered
genetic correlation is whatever the finite draw gave. At 200 QTL, target $r_g=0.4$, 500
draws: mean 0.399, **sd 0.065, range [0.177, 0.598]**. The congruence delivers 0.400
exactly at any $m\ge k$. The miss shrinks only as $1/\sqrt m$.

### 1.3 The non-additive method

Functional model (AlphaSimR's, Chu et al. 2024's), $g$ = dosage of allele 1:

$$g_i=\sum_j a_j(x_{ij}-1)+\sum_j d_j\,\mathbb 1[x_{ij}=1]+\sum_{(k,l)}e_{kl}(x_{ik}-1)(x_{il}-1).$$

Its exact statistical (NOIA) re-expression is

$$g-\text{const}=\mathbf Z_A\boldsymbol\alpha+\mathbf Z_D\mathbf d+\mathbf Z_{AA}\mathbf e,\qquad
\alpha_j=a_j+b_jd_j+\sum_l e_{jl}c_l,\quad c_l=2p_l-1.$$

Under the genic anchor $b_j=q_j-p_j$. The re-expression turns three functional effect
matrices into three **targets** $(\mathbf G_A,\mathbf G_D,\mathbf G_{AA})$, each with its own
anchor:
$\mathbf M_A=\mathrm{diag}(2pq)$, $\mathbf M_D=\mathrm{diag}((2pq)^2)$,
$\mathbf M_{AA}=\mathrm{diag}(4p_kq_kp_lq_l)$.

Construction, highest order first, all closed form:

1. $\mathbf B_{aa}\leftarrow$ congruence to $\mathbf G_{AA}$;
2. $\mathbf B_d\leftarrow$ congruence to $\mathbf G_D$ (dominance degrees $h\sim N(0.19, 0.097)$ scale $|\mathbf B_a^*|$ for the architecture);
3. with $\mathbf C=\mathrm{diag}(\mathbf b)\mathbf B_d+\mathbf E_c$ now fixed, solve the matrix quadratic
   $(\mathbf B_a^*\mathbf T+\mathbf C)^\top\mathbf M_A(\mathbf B_a^*\mathbf T+\mathbf C)=\mathbf G_A$ by completing the square.

Stage 3 is feasible **iff** $\mathbf G_A\succeq\mathbf R-\mathbf Q^\top\mathbf P^{-1}\mathbf Q$. In
words: dominance and epistasis already imply a floor on the additive covariance, and the
function refuses a $\mathbf G_A$ below it. For $k=1$ this is Zeng et al. (2013) Appendix A.
An inbreeding-depression target is a third quadratic, exact for $k=1$ or diagonal
$\mathbf G_D$ (the closed form of AlphaSimR's `altAddTraitAD()` optimiser). With no
$\mathbf G_D$, $\mathbf G_{AA}$ it reduces exactly to 1.1.

**Verification in the source project** (suites listed in §1A): additive — 12 core tests + 13 v2 tests, 2026-08-28
bug audit closed. Non-additive — 18 tests / 89 checks, exact to $10^{-13}$ under both anchors;
validation/diagnostics parity with the additive file (2026-09-19,
`non-additive/verification_nonadd_implementation.html`).

---

## 1A. Source code and test suites — where everything lives

Everything to be ported lives in **`/Users/austinputz/Claude/simulate_qtl_effects/`**
(below, `$SRC`). It is a git repository, not a package: plain R files `source()`d by
plain-script test files. Each test file finds the repo root by walking up to `AGENTS.md`,
so it runs from any working directory. Suite results below are from re-running each suite
on 2026-09-22 (all pass).

### Implementation files

| File | Role | Port? |
|---|---|---|
| `$SRC/R/qtl_effects.R` | **The additive method, current version.** `sim_qtl_effects(X, G, architecture, anchor, reference, ploidy, rel_tol, warn_bounds, warn_on_mismatch, marginal, marginal_iters)`. Congruence (the `.congruence` closure inside `sim_qtl_effects()`), anchors (`genic`, `marginal_observed`, `realised`, `dual`, `reference`), validation helpers `.qtl_validate_*`, `.qtl_psd_eigen()`, diagnostics `.qtl_relative_spectrum()` / `.qtl_spectrum_text()`, dual basis `.qtl_dual_basis()`, samplers `arch_gaussian()` … `arch_maf()`, `marginal_*()` and the `marginal =` iteration, `print.qtl_effects()` | **Yes — Part A** (congruence, anchors, validation, diagnostics). Samplers and `marginal =` later (Q5) |
| `$SRC/R/qtl_effects_paper.R` | The same algebra frozen at manuscript v4. The paper's numbers are built from it | **No.** Reference only; its test suite is ported, pointed at the current code |
| `$SRC/non-additive/R/qtl_effects_nonadd.R` | **The non-additive method.** `sim_qtl_effects_nonadd(X, G_A, G_D, G_AA, pairs, anchor, B_a, B_d, B_aa, dd_mean, dd_sd, inbr_depr, rel_tol, warn_bounds, warn_on_mismatch)`. Stages: `congruence()` (thin-SVD form, fixed 2026-09-19), `solve_additive_stage()` (matrix quadratic + floor), `solve_dd_mean()` (inbreeding depression). Measurement: `nonadd_covariates()` (NOIA design + per-component anchors), `nonadd_decompose()`. Comparators: `scale_package_style()`, `ss_noia()`, `zeng_appendix_A()`. Helpers `.na_*`, `print.qtl_effects_nonadd()` | **Yes — Parts B and C.** `scale_package_style()` / `ss_noia()` are comparison code, not ported. `zeng_appendix_A()` becomes a test oracle |

The two files carry **two copies of the same congruence algebra** (`.congruence` in the
additive file, `congruence()` in the non-additive one). Port it **once**, to
`R/qtl_congruence.R`, and have both generators call it.

### Test suites

| File | Covers | Size | Run (from `$SRC`) |
|---|---|---|---|
| `$SRC/tests/test_qtl_effects_paper.R` | **The core additive surface** (v1): realised exactness, all-heterozygote F1 case (Codex F1), rank-1 target / rank-deficient architecture (Codex F2), matrix diagnostics (F4), dual anchor, custom `reference` (weights and genotypes), singular target under dual, input validation (F5), inbred panel with repulsion LD, generation-1 miss, non-HWE vs HWE founders, linked transmission | 12 tests, ~2 s | `Rscript tests/test_qtl_effects_paper.R` |
| `$SRC/tests/test_qtl_effects.R` | **Only the v2 additions** (samplers, `marginal =`). Test 0 checks v2 reproduces v1 bit-for-bit | 13 tests (0–12), ~2 s | `Rscript tests/test_qtl_effects.R` |
| `$SRC/non-additive/tests/test_qtl_effects_nonadd.R` | Exactness under both anchors, reduction to additive (test 4), Zeng et al. Appendix A, feasibility floor, reduced-rank congruence (test 9), rank-deficient `B_a` message (test 10), mismatch warning (test 16), supplied architectures (test 18), validation parity | 18 numbered tests, 89 PASS checks, ~20 s | `Rscript non-additive/tests/test_qtl_effects_nonadd.R` |

**Trap:** `test_qtl_effects.R` is *not* the additive method's main suite. The congruence,
anchor and validation tests Part A needs are in `test_qtl_effects_paper.R`, which
`source()`s the frozen paper file. When ported they must call the tidybreed port instead.

### How they become testthat files

The source suites are plain scripts: `stop()` on failure, `ok(name, cond)` / `cat("PASS")`
reporting, shared fixtures such as `mk_panel()` defined inline. For each ported test, keep
its **number and name in the `test_that()` description**
(`"paper-3: rank-1 target, rank-deficient architecture (Codex F2)"`,
`"nonadd-9: reduced-rank congruence"`) so each testthat case traces back to its source. Move the
inline fixtures to `tests/testthat/helper-qtl-effects.R`.

| Source suite | tidybreed target | Part |
|---|---|---|
| `test_qtl_effects_paper.R` tests 1–4, 6, 8, 9 (on the ported internals) | `tests/testthat/test-qtl-congruence.R` | A |
| same, tests 5, 7, 11 (dual) | same file, `skip_on_cran()`; internals only (Q6) | A |
| same, tests 10, 12 (generational) | `tests/testthat/test-define_additive_effects-anchor.R`, using `add_offspring()` instead of the script's own transmission code | A |
| *(as built, 0.73.0)* | Test 10 is a property of the **dual** anchor (internal, Q6), so it is an algebraic check inside `paper-9` in `test-qtl-congruence.R`. Test 11 (dual) is there too, with the source's own transmission code. Only test 12 goes through `add_offspring()` | A |
| `test_qtl_effects.R` | not ported with Part A (Q5) | — |
| `test_qtl_effects_nonadd.R` construction tests | `tests/testthat/test-genome-effects-calibration.R` | C |
| same, decomposition / measurement tests | `tests/testthat/test-extract_genetic_variance.R` | B |
| gates in §11 (tidybreed-level: storage, `ind_tgv`, owners) | `test-define_genome_effects.R`, `test-define_additive_effects-anchor.R` | A, C |

### Supporting documents in `$SRC` (read, not ported)

| Path | What it is |
|---|---|
| `$SRC/paper/gse_manuscript.tex` (v4) | The additive method's mathematics — Propositions 1–2, Theorems 1–4 |
| `$SRC/docs/notes/note_tidybreed_readiness.qmd` / `.html` | Status note; render-time checks of §1.2 and §4.2 of this plan |
| `$SRC/non-additive/verification_nonadd_implementation.qmd` / `.html` | Review of the non-additive file after the 2026-09-19 parity pass |
| `$SRC/non-additive/README.md` | The non-additive derivation and simulator comparisons |
| `$SRC/AGENTS.md` | Orientation for the whole source project |

---

## 2. Current state in tidybreed (audit, 0.71.0)

- **Storage** (`R/define_genome.R:446-535`): `genome_effects` (one row per term:
  `trait_name`, `effect_owner`, `effect_name`, `genome_value`), `genome_effect_members`
  (per term × locus: `contrast_name` ∈ {`additive`, `dominance`, `indicator`},
  `center_value` or `(copy_count_value, dosage_value)`), and `genome_effect_member_origins`
  (scope). Views: `genome_effect_terms`, `genome_effect_loci`, `ind_tgv_total`
  (`R/genome_effects_helpers.R:31,71,89`).
- **Contrasts** (v4.9 §Contrast definitions). `additive` = $g-2c$ (ploidy-aware).
  `dominance` = Cockerham $(-2p^2, 2pq, -2q^2)$ at $c=p$, **diploid loci only** (the writer
  rejects it elsewhere). `indicator` = $\mathbb 1[\text{state}]$. A multi-member term is the
  product of member values.
- **Components** (`R/genome_effects_eval.R:80`): today `order1_additive`,
  `order1_dominance`, `order1_other` (order-1 indicator), `interaction` (≥ 2 members).
  These describe *declared structure*, not variance components. §6B renames them to
  `additive`, `dominance`, `indicator`, `interaction`, and the rest of this plan uses the
  new names.
- **Writer**: `define_genome_effect_terms()` (`R/define_genome_effect_terms.R`; named
  `define_genome_effects()` until step 1, §0A), with `ad_terms(coding = "functional" |
  "cockerham")` and `genotype_terms()`. Everywhere else in this plan,
  `define_genome_effects()` means the **new Part C generator**.
- **Generator**: `define_additive_effects()` (`R/define_additive_effects.R:176`), reserved
  owner `GE_GENERATED_OWNER = "generated"` (`:507`; renamed in step 1, §5). Frequencies come from
  `base_tbl` via `extract_allele_freq()`. The scalar rescale is `rescale_effects_to_target()`
  (`:741`). The multi-trait path is `MASS::mvrnorm` (`:417, :425`) then the same rescale
  (`:443-456`).
- **Targets**: `trait_var_comp`, with `define_effect_cov_matrix()` already accepting
  `"additive"`, `"dominance"`, `"additive_by_additive"` (`R/define_effect_cov_matrix.R:124`;
  renamed from `gen_add` / `epistasis` in step 1, §6B). `schema.R:331` documents the
  latter two as reserved, with no generator yet.
- **Results**: `ind_tbv` (`add_tbv()`, reserved additive owner only) and `ind_tgv`
  (`add_tgv()`, every term).

---

## 3. The mapping — no DDL change to the `genome_effect*` tables

Every block of the non-additive output is an existing contrast. Two codings are possible:

| Block | Functional coding | **Statistical coding at base $p$ (proposed)** |
|---|---|---|
| additive, locus $j$ | `additive`, $c=0.5$, value $a_j$ | `additive`, $c=p_j$, value $\alpha_j=a_j+(q_j-p_j)d_j+\sum_l e_{jl}(2p_l-1)$ |
| dominance, locus $j$ | `indicator` $(2,1)$, value $d_j$ | `dominance`, $c=p_j$, value $d_j$ |
| A×A, pair $(k,l)$ | **one** term, two `additive` members at $c=0.5$, value $e_{kl}$ | **one** term, two `additive` members at $c=p_k$, $c=p_l$, value $e_{kl}$ |
| `ind_tgv` rows | `additive`, `indicator`, `interaction` | `additive`, `dominance`, `interaction` |
| base mean of the genetic value | $\mu\ne0$ | 0 in expectation under HWE + LE at base $p$ |

The two codings agree on the total **up to a trait-specific constant**, for every
genotype:

$$g_{\text{functional}}-g_{\text{stored}}=\mu=\sum_j\big[a_j(2p_j-1)+2p_jq_jd_j\big]+\sum_{(k,l)}e_{kl}(2p_k-1)(2p_l-1).$$

The identities used are
$\mathbb 1[g=1]=2pq+(q-p)(g-2p)+x_D(g;p)$ and
$(g_k-1)(g_l-1)=(u_k+c_k)(u_l+c_l)$, with $u=g-2p$. Both are algebraic identities in $g$ at
any fixed $p$, not population statements. So contrasts between individuals are identical,
and raw totals differ by $\mu$. The ~1e-14 in §4.2 was measured **after centring**, which
removes $\mu$. The stored model's mean is 0 only in expectation. A finite, selected,
non-HWE or LD base can have a nonzero sample mean of `dominance` and `interaction`. The
generator never adds $\mu$ anywhere (§6A).

An A×A pair is **one two-member term**, not a 3×3 `genotype_terms()` surface. The surface
would be up to nine indicator terms for one coefficient, all labelled `interaction`, with
no link back to $e_{kl}$.

Worked `terms` for one locus with $(\alpha, d)$ and one pair, as the generator would build
them per trait:

```r
data.frame(
  term_id       = c("L10_a",    "L10_d",     "L10xL44",  "L10xL44"),
  locus_name    = c("Locus_10", "Locus_10",  "Locus_10", "Locus_44"),
  contrast_name = c("additive", "dominance", "additive", "additive"),
  center_value  = c(p10,        p10,         p10,        p44),
  genome_value  = c(alpha10,    d10,         e_10_44,    e_10_44))
```

---

## 4. Why the coding decides what the "TBV" is

### 4.1 The claim

Under statistical coding at the anchor's base frequencies, `add_tgv()`'s components **are**
the NOIA components, each up to a constant:

| `component_name` | equals | genic covariance (by construction) |
|---|---|---|
| `additive` | $\mathbf Z_A\boldsymbol\alpha$ — the breeding value | $\mathbf G_A$ |
| `dominance` | $\mathbf Z_D\mathbf d$ — the dominance deviation | $\mathbf G_D$ |
| `interaction` | $\mathbf Z_{AA}\mathbf e$ — the A×A value | $\mathbf G_{AA}$ |

Under functional coding `additive` = $\sum_j a_j(g_j-1)$. That is the functional
additive effect, which is not the breeding value once $d\ne0$ or $e\ne0$.

### 4.2 The evidence (source project, rendered in the status note)

Two traits, 300 loci, 100 A×A pairs, $\mathbf G_A$ diag 4 / 3, $\mathbf G_D$ diag 1 / 0.6,
$\mathbf G_{AA}$ diag 0.5 / 0.4, $n=1500$:

| Anchor | total, statistical − functional (centred) | cor(`additive`, BV), statistical | cor(`additive`, BV), functional | genic var of stored $\alpha$ |
|---|---|---|---|---|
| genic | 1.3e-14 | 1.000 / 1.000 | 0.807 / 0.874 | 4.000 / 3.000 |
| realised | 1.3e-14 | 1.000 / 1.000 | 0.828 / 0.869 | 4.134 / 3.049 |

### 4.3 Four consequences to write into the docs

1. **Realised anchor.** The method's realised $\boldsymbol\alpha$ uses the observed within-locus
   regression $b_j=\mathrm{Cov}(w_j,x_j)/\mathrm{Var}(x_j)$. tidybreed's `dominance` contrast
   is the HWE form, which fixes $b_j=q_j-p_j$. So under `anchor = "realised"` the generator
   **re-derives** $\alpha_j=a_j+(q_j-p_j)d_j+\sum_l e_{jl}c_l$ for storage rather than storing
   the method's `B_alpha`. The total is still exact, but the additive/dominance split in
   `ind_tgv` is then HWE-referenced (4.134 vs target 4 above). Document it: the realised
   decomposition is a *measurement* (`extract_genetic_variance()`, §8), not a storage
   coding. A realised `dominance` contrast would need three genotype frequencies per locus
   rather than one `center_value`. That is a DDL change nothing needs yet.
2. **One centre per locus.** All three blocks take $p$ from one `base_tbl`. That rules out
   the mismatched-centre trap (v4.9 review §3.1) for generated effects by construction.
   The §3.1 write-time warning is still wanted for hand-written terms.
3. **The breeding value is referenced to the base.** $\alpha_j$ depends on $p$ once
   $d\ne0$ or $e\ne0$. As frequencies drift, `additive` stays the additive value
   relative to the **base** population. The total stays correct. The breeding value *in
   generation t* is a re-projection at generation-t frequencies. It belongs in a later
   `extract_breeding_values()` (§8), not in `ind_tgv`. It is defined for the common-scope
   diploid class only. This is the manuscript's "the anchor names a population" point, one
   level up.
4. **Inbreeding depression is correct under either coding.** Under inbreeding $F$,
   $E[x_D]=-2Fpq$, which equals the functional $2pq(1-F)d$ minus its base value. The mean
   contract (v4.9 §Mean contract: raw sum, no mean added) is unchanged. Under statistical
   coding the base mean is 0 in HWE/LE expectation (§3), and inbreeding shows up as a
   negative shift, which is what users expect.

---

## 5. Owners *(decided 2026-10-01: one reserved owner, Q3)*

**One reserved owner, `generated`, for both generators.** `define_additive_effects()` and
`define_genome_effects()` write every term they build under `GE_GENERATED_OWNER =
"generated"`, the only entry in `GE_RESERVED_OWNERS`. After consolidation the owner has two
jobs. It **protects** generated terms: the writer `define_genome_effect_terms()` refuses
the reserved owner with no exported override (Q23, 0.73.2; only the internal
`.ge_write_terms()` engine can write it). And it is the
**replacement boundary**: a generator re-run replaces `generated` terms and never touches
`custom` or other user owners, which still sum, per v4.9. Which function wrote a model is
not recorded. What the model *is* can be read from its terms (principle 6).

**Replacement rules** (Part C, step 5):
- **`define_genome_effects()` replaces the trait's whole `generated` model**
  (`mode = "replace_owner"`). It always builds a fresh, common-scope model, so line-scoped
  additive variants from earlier `define_additive_effects()` calls are deleted with the
  rest. Its `message()` says how many terms it replaced and how many of them were
  line-scoped, so a crossbred model is never wiped silently. The roxygen says so too.
- **`define_additive_effects()` keeps `mode = "replace_scope"`** within `generated`, so
  common / line-A / line-B variants still coexist. It **refuses** when the trait's
  `generated` terms include anything other than order-one `additive` terms, i.e. a
  `dominance` or interaction term from `define_genome_effects()`. Replacing $\boldsymbol\alpha$
  would leave $\mathbf d$, $\mathbf e$ calibrated against the old one, which can silently
  put $\mathbf G_A$ below its floor. The error names the fix, which exists without any new
  clearing call: re-run `define_genome_effects()` with an additive-only target
  (`trait_var_comp_tbl = get_table(pop, "trait_var_comp") |> filter(effect_name ==
  "additive")`). That replaces the model with additive-only terms, which
  `define_additive_effects()` can then extend. The stored `dominance` /
  `additive_by_additive` targets are still in `trait_var_comp`, though, so that later
  `define_additive_effects()` call needs the same filter, or the rows removed with
  `remove_rows()` (§6C).
- Both checks run **before any write**, including target writes (§7.1).
- An additive-only model is **row-identical** whichever generator wrote it, owner included
  (gate C4).

**Rename timing.** The string `generated_additive_tbv` → `generated` and the constant
`GE_ADDITIVE_OWNER` → `GE_GENERATED_OWNER` are part of the rename-only release (step 1).
`add_tbv()` keeps reading that owner until consolidation removes `ind_tbv` (step 3); only
`define_additive_effects()` writes it until then. Messages that say the reserved owner is
"owned by `define_additive_effects()`" (e.g. `.ge_resolve_deletes()`) are reworded to name
both generators in step 5.

**Why consolidation must come first.** `add_tbv()` reads only the reserved owner's
order-one additive terms. A trait defined by `define_genome_effects()` would have an
`ind_tbv` holding $\boldsymbol\alpha$ only, and `add_phenotype()` reads `ind_tbv`, so
dominance and A×A would never reach a phenotype.

---

## 6. Precondition P1 — consolidation: `ind_tbv` is replaced by `ind_tgv` (step 3)

*Integrated 2026-10-02 from `plans/consolidate_genetic_values.md` (created 2026-09-06),
which is now a pointer to this section. Its old text survives in git history. What
changed in the merge: its Task 1a / 1b and their gate and question were dropped, its
vocabulary, owner and writer names follow this plan, and its gates were corrected (see
"Corrections" below).*

### 6.1 The problem

`ind_tbv.tbv_value` is built by `add_tbv()` from the reserved owner's order-one additive
terms only. `ind_tgv`'s additive row sums **every** additive term, custom owners
included. They agree whenever a user only calls `define_additive_effects()`, which is the
common case and exactly why a disagreement is dangerous. One custom additive term, and the
two tables differ with no error. It is worse for non-additive terms: custom-owner
`dominance` / `indicator` terms written with the writer **today** reach `ind_tgv` but
never `ind_tbv` or a phenotype. Two definitions of one quantity are technical debt
(CLAUDE.md), and this removes one.

### 6.2 Target state

- **`ind_tbv` and `add_tbv()` are deleted.** No view, alias or wrapper is kept (pre-1.0).
  `ind_tgv` is the single table of true genetic values. `add_tgv()` computes and writes
  every component in one pass.
- **The breeding value is `component_name = 'additive'`**, for effects written in
  statistical coding at a single base `p` per locus, which every generator guarantees
  (§4). For hand-written functional terms (`ad_terms(coding = "functional")`, indicator
  surfaces), `additive` is the functional additive effect, not the breeding value. That
  sentence goes in the roxygen of `add_tgv()`, `define_genome_effect_terms()` and
  `ad_terms()`. It replaces the stale-coefficient warning, which is deleted with
  `add_tbv()`.
- **Value components** are §6B's: `additive`, `dominance`, `indicator` (a one-locus
  term's `contrast_name`), and one `interaction` row for any term over two or more loci.
  The interaction row is **not** split by type (A×A / A×D / D×D). Per-type variances come
  from `extract_genetic_variance()`, and per-type values can be computed from
  `genome_effect_terms`. Writer-declared component labels for hand-entered surfaces were
  considered and dropped. For generated effects the components already **are** the
  orthogonal (genic) decomposition (§4). For hand-written surfaces that decomposition
  stays a separate future project.
- **The total is never stored.** It is the `ind_tgv_total` view (CLAUDE.md: no `'total'`
  row), which sums exactly (§6A). A stored total would make `SUM(tgv_value)` double-count
  in every ad-hoc query.
- **`component_name` is an R-validated closed set**, not an SQL `CHECK`, so adding a
  name later is a one-line change.
- `ind_tgv`'s DDL is unchanged. It has **no** `replicate` column (CLAUDE.md;
  `R/define_trait.R:227-241`); only archive copies carry one.

### 6.3 Tasks

1. **Value-component rename** (§6B): `order1_additive` / `order1_dominance` /
   `order1_other` → `additive` / `dominance` / `indicator`, at `R/genome_effects_eval.R:80-83`
   and `R/schema.R:423,545`, and in `add_tgv.R`, `man/add_tgv.Rd`,
   `test-genome-effects-eval.R` and `test-genome-effects-schema.R`.
2. **Delete `ind_tbv` and `add_tbv()`.** Also remove:
   - the table registries (`R/sql_utils.R:113,155,182,299,310`; `R/schema.R`
     `.schema_table_order()`, `.ind_descriptions()`, and lines 322, 524-534, 612 and 780;
     guarded by `test-schema-registries.R` and `test-schema-print.R`);
   - the evaluator helpers that exist only for `add_tbv()`: `.gev_reserved_additive()`,
     `.gev_require_terms(tbv = TRUE)` and `.gev_warn_tbv_stale()`
     (`R/genome_effects_eval.R:170-185, 808-890`).

   `restore_pop()` refuses a file that still has `ind_tbv` (pattern at
   `R/restore_pop.R:118-142`).
3. **What `add_tgv()` inherits from `add_tbv()`:**
   - `index_names`, and `type` renamed `weight_type` (naming rule 1);
   - `...` custom-field forwarding (`R/add_tbv.R:151`, `R/sql_utils.R:362`);
   - **`ind_true_index`**, computed from a `component_name` argument (naming rule 1;
     decided 2026-10-04) that defaults to `"additive"` (selection indices are on
     breeding values), with `"total"` allowed. `ind_true_index` gains a
     `component_name` column, so the two coexist.
4. **`add_index()`.** Its table map (`R/add_index.R:117-121`) gains `ind_tgv` →
   `tgv_value`. Callers filter the component:
   `get_table(pop, "ind_tgv") |> filter(component_name == "additive") |> add_index("meat")`.
   Its existing "more than one value per individual × trait" error is the safety net
   against silently summing components; gate T6 confirms it fires on the component
   dimension.
5. **Phenotypes read the total**: precondition P2 (§6A), for every contributor and path,
   including `self`, simple phenotypes and the formula DSL (`formula_tgv`, Q18).
6. **Keep the first-principles oracle.** `tests/testthat/test-add_tbv.R:16-40` recomputes
   the additive value from first principles. It is not golden output. It is retargeted at
   `ind_tgv` filtered to `additive` and must survive (gate T2).
7. **Docs:** CLAUDE.md's hard rules (the D7 sentence, "One evaluator … `add_tbv()` reads
   only reserved-owner…"), design principle 4's action list, the Roadmap line, both
   skills, the 2 vignettes (`tidybreed-introduction.Rmd`,
   `swine/swine-time-based-age-at-puberty-sex-semen.R`), and `package_summary.md`.

**Blast radius** at 0.71.1 (recount when step 3 starts): `ind_tbv` / `add_tbv` appear in 17
`R/` files, 25 test files (including `helper-parity.R`, `test-add_tbv.R`,
`test-add_tbv_index.R`, `test-add_phenotype.R`, `test-phenotype_composite.R`,
`test-remove_rows.R`), 3 `dev/` benchmarks, 2 vignettes, and `package_summary.md`.

**Corrections made in the merge** (cross-check 2026-10-02):
- Its `ind_tgv` DDL listed a `replicate` column. Removed.
- Its gate "`ind_true_index` reproduces the pre-migration value" compared against old
  output, which CLAUDE.md forbids. Now computed in the test (T7).
- Its "unfiltered `add_index()` errors" gate passes vacuously on an additive-only model.
  It now runs on a multi-component model (T6).
- `genome_effect_types` was renamed `component_names` in v0.68.0.
- Its contributor task chose "do the swap here". That is now P2, extended to `self`,
  simple phenotypes and the DSL.

### 6.4 Why it must come before Part C

`add_tbv()` reads only the reserved owner's order-one additive terms, and `add_phenotype()`
reads `ind_tbv`. A trait defined by `define_genome_effects()` would therefore reach its
phenotypes with $\boldsymbol\alpha$ only, and the calibrated $\mathbf G_D$ and
$\mathbf G_{AA}$ would never be expressed (§5, §6A).

**Not included:** an estimated counterpart `ind_egv`, if non-additive genomic prediction is
ever added. `add_ebv()` and `ind_ebv` are untouched.

---

## 6A. Precondition P2 — phenotypes must see the total genetic value

**The gap.** Today `add_phenotype()` reads only the additive value for every
contributor, `self` included. `phenotype_components.component_names` defaults to
`"order1_additive"` (`"additive"` after §6B) and is documented as reserved (`R/define_phenotype.R:68`,
`R/open_pop.R:324`). The original consolidation plan moved only **dam/sire/group**
lookups to the total. It did not move `self`, or simple phenotypes, which have no
`phenotype_components` row at all. After Part C, a trait written by
`define_genome_effects()` would have its $\mathbf G_D$ and $\mathbf G_{AA}$ calibrated
exactly, and its phenotypes would ignore them. The acceptance tests would all pass,
because they measure `ind_tgv`, and nobody would notice.

**The change** (rows, and one DDL default):

- The genetic value that enters a phenotype is `ind_tgv_total`, for every
  `contributor_type`, including `self` and simple phenotypes.
- `phenotype_components.component_names` default becomes `"total"`. That is a DDL default
  in `CREATE TABLE phenotype_components` (`R/open_pop.R:324`), also set in
  `R/define_phenotype.R:538`. Listing specific components (e.g. `"additive"`) stays
  possible, and that is now the only thing the column does, so it stops being "reserved".
- `add_phenotype()`'s prerequisite call becomes `add_tgv()` (it is `add_tbv()` today).
- **`ind_tgv_total` sums deterministically.** The view was a plain `SUM(tgv_value)`.
  Once phenotypes read it, its result must be bit-identical whatever DuckDB's thread
  count (CLAUDE.md, "Identical means bit-identical"). With up to four components per
  individual × trait, a parallel `SUM()` is not guaranteed to be. *As built (0.74.0):*
  the view adds the components in a fixed order,
  `list_sum(list(tgv_value ORDER BY component_name))`. The plan said `GEV_ACC_TYPE`,
  but a `DOUBLE → DECIMAL(38, 18) → DOUBLE` round trip is not exact: it moved about
  10% of one-component totals by one ulp, so an additive-only trait's total would not
  have been its breeding value bit for bit (gate T3). The ordered sum returns a single
  row unchanged and is still a function of the stored rows alone. Listed components
  (`component_names`) are summed the same way.
- **Group-contributor sums accumulate exactly too.** `group_sum()` / `group_mean()` run a
  second reduction over mates in `.group_mate_tbv()` (`R/contributor_tbv.R`). Since
  0.71.2 (step 0b, B-1) that sum goes through `GEV_ACC_TYPE`. When step 3 moves its
  input from `ind_tbv` to `ind_tgv_total` (or a listed component's rows), it keeps that
  accumulation, and the 0.71.2 test `test-group-contributor-determinism.R` is retargeted.
  The `mean` divides the exact sum by the integer `COUNT`. Gate PH5 covers it.
- **The simple-phenotype precheck changes.** `R/add_phenotype_stages.R:148-163` errors
  unless the reserved owner has order-one additive terms ("call
  define_additive_effects() first"). After P2 a trait with **any** terms is valid, including
  a model written only with `define_genome_effect_terms()`. The check becomes "the trait has
  at least one term", and its message stops naming `formula_tbv`.
- **Prevalence thresholds use the total genetic variance.** `R/add_phenotype_stages.R:1146`
  sets a binary trait's liability threshold from `get_trait_var(pop, "additive", t)` alone.
  Once the liability carries the total genetic value, the threshold uses the variance of
  the **active** model *(decided 2026-10-02, Codex review finding 2)*: the sum of the
  stored diagonals (`line_name IS NULL` rows, §6C) of `additive`, `dominance` and
  `additive_by_additive`, counting a block **only if the trait's model has terms of that
  kind** (any owner). A stored target the generator was told to leave out (§6C,
  "available targets") therefore does not enter the threshold. The call **errors**, naming
  `define_phenotype(thresholds = )`, when the threshold cannot be derived from targets:
  - the model has terms of a component with no stored target, or
    terms outside the three blocks (indicator surfaces, higher-order interactions);
  - the trait has **any active term not owned by `generated`** (hand-written
    `define_genome_effect_terms()` terms), even when a target is stored. Only generators
    write `generated` terms, and from step 3 they always calibrate (Q21), so
    "`generated`" proves the terms were calibrated to the stored target. A user-owner
    term proves nothing about it;
  - the phenotype is composite (it has `phenotype_components` or `formula_tgv`): its
    genetic liability is a combination of several traits and contributors, which no
    stored diagonal describes.
  The result is still the *target* at the reference population, as today: an
  approximation the roxygen states. Rejected: an empirical threshold from each batch's
  realised liability, which makes prevalence exact per call but different between
  batches. The composite case was a silent bug, B-2 (§10, step 0b). Since 0.71.2
  `define_phenotype()` already refuses `prevalence` with `components` / `formula_tbv`,
  and `.ap_check_prevalence()` (PLAN stage) refuses a missing target. Step 3 extends
  that check to the active-block rule rather than adding a new one.
  The owner rule is Q21's decision (a). To keep it sound after generation,
  `define_effect_cov_matrix()` refuses to write a genetic block whose traits already have
  `generated` terms of that kind **at the block's scope**: a `line_name = NULL` block
  against common-scope terms, and a `line_name = "C"` block against line-C-scoped terms.
  So writing a line-C target and then generating line-C effects still works (A17).
  Without this refusal, `remove_rows()` plus a new target would leave the old terms
  calibrated to a target that is gone. To change a target: `remove_rows()` the old block,
  then regenerate with `define_additive_effects(G = )`, which writes target and terms in
  one transaction through the shared internal insert (§7.1). Step 3 rewords step 2's
  "already stored" error so that it gives this sequence, not `remove_rows()` alone.
  `method = "union"` and line-scoped effects are calibrated *toward* the target, not
  exactly (§7.3, §7.6). The threshold is an approximation in those cases, as the roxygen
  already says of the reference-population target.

Under an additive-only model the total equals `additive`, so this is a no-op for
every existing test. It belongs in the consolidation release (step 3, §10), not in
Part C.

**Every phenotype path.** Stage 1 today calls `add_tbv()` for simple phenotypes, for
explicit `phenotype_components`, and for `formula_tbv` (renamed `formula_tgv`, Q18), and `R/contributor_tbv.R` reads
`ind_tbv` for dam/sire/group. P2 moves **all** of these: `self`, dam, sire, group, simple,
composite and the formula DSL. Gates PH1–PH3 (§11) exercise each one. The D7 failure
contract is unchanged in substance: its one surviving write becomes the `add_tgv()` upsert.
CLAUDE.md's D7 sentence and `test-add_phenotype_failure_contract.R` are retargeted.

**The formula DSL** is specified in Q18 (`formula_tgv`, named-only `component =` and
`table =`). `component = "total"` is accepted and is the default, matching
`component_names`.

**What `component_names` means.** `"total"` reads `ind_tgv_total`. A listed component
(`"additive"`, `"dominance"`, …) reads that row. A component the trait's model has no
terms for contributes 0, because the model has none. A name outside §6B's value vocabulary
errors in `define_phenotype()`. An individual with **no** `ind_tgv` row for a component
trait is missing, and `missing_component_action` applies as today.

**Mean: `phenotype_meta.mean` is an intercept.** The phenotype is
`mean + (genetic value as stored) + random effects + residual`, per the v4.9 mean contract
(raw sum, nothing added). It is **not** a target for the realised base mean. Under
generated (statistical) coding the genetic value has mean 0 in HWE/LE expectation at base
$p$ (§3). So `mean` equals the base phenotypic mean in expectation, and the realised mean of
a finite, selected or non-HWE base differs by that sample's genetic mean. For
hand-written functional terms the offset is $\mu$ (§3). A user who wants a particular
realised base mean sets `mean` to the target minus the base's mean `ind_tgv_total`,
measured after `add_founders()`. Say this in the `define_phenotype()` roxygen. Do not add
$\mu$ anywhere; this is `ad_terms()`'s existing rule.

---

## 6B. Naming decision — genetic effect vocabulary *(decided 2026-09-23)*

Two concepts, two columns, **one set of words**:

- a **variance component** (a share of $\mathrm{Var}(g)$): `trait_var_comp.effect_name`
  for targets, `extract_genetic_variance()`'s `effect_name` for measurements;
- a **value component** (the part of an individual's genetic value that came from
  terms of one kind): `ind_tgv.component_name`.

| Concept | Variance: `effect_name` | Value: `component_name` | Replaces |
|---|---|---|---|
| additive | `additive` | `additive` | `gen_add` / `order1_additive` |
| dominance | `dominance` | `dominance` | `dominance` (unchanged) / `order1_dominance` |
| additive-by-additive | `additive_by_additive` | `interaction` | `epistasis` / `interaction` |
| additive-by-dominance (future) | `additive_by_dominance` | `interaction` | none |
| dominance-by-dominance (future) | `dominance_by_dominance` | `interaction` | none |
| hand-entered order-1 genotype surface | none | `indicator` | `order1_other` |
| total | `total` (extract output only) | `ind_tgv_total` view, never a row | none |
| variance of terms outside the supported decomposition | `unpartitioned` (extract output only) | none | draft `residual` |
| sum of cross-block covariances (realised, §8) | `between_components` (extract output only) | none | none |

**Rules behind the choice**

1. **Lowercase snake_case full words.** No letters (`A`, `D`, `AA`), no symbols (`AxA`),
   no abbreviations (`gen_dom`). `E` in particular is never used: in quantitative
   genetics it means the environment ($P = G + E$).
2. **Name the specific component, not the umbrella.** `epistasis` is the family name;
   the calibrator targets A×A specifically. `"by"` is the literature's own wording
   ("additive-by-additive epistasis").
3. **Same word in both columns on purpose.** For effects generated with
   `anchor = "genic"` (Cockerham coding at the base $p$, §4), the genic covariance of the
   stored `additive` values at the generation base **is** the `additive` target, so a
   join across the two tables says something true. (Variance is a property of a
   population, not of one individual's value.) Under `anchor = "realised"` the stored
   values are HWE-coded and their covariance can differ from the target (§4.3.1, 4.134
   vs 4); compare the target with `extract_genetic_variance(anchor = "realised")`
   instead. For
   hand-written functional terms it is not; that is documented where those terms are
   written (`ad_terms(coding = "functional")`, `genotype_terms()`), not encoded in the
   name.
4. **The value-component rule is one sentence:** *the `contrast_name` of a one-locus
   term, or `interaction` for a term over two or more loci.* No `order1_` prefix. The
   interaction row is **not** split by type (§6.2).
5. **The genetic `effect_name` values are reserved words.** `define_effect_cov_matrix()`
   routes on the string. A **finite** constant, `GENETIC_EFFECT_NAMES = c("additive",
   "dominance", "additive_by_additive")`, goes to `trait_var_comp`. `residual` and anything
   else goes to `phenotype_var_comp`. Names are added to the constant when the package can
   calibrate them. A pattern such as `*_by_*` would capture user-chosen random-effect names
   and targets nothing can realise (review 2026-09-23). A user-named random effect that uses
   a reserved name is an error. The derived output names `total`, `unpartitioned` and
   `between_components` (§8) are also refused as input.
   **Reserved but not yet supported** (decided 2026-09-24): `additive_by_dominance` and
   `dominance_by_dominance`, in a second constant `GENETIC_EFFECT_NAMES_FUTURE`.
   `define_effect_cov_matrix()` refuses them with a "not yet supported" error that names
   `define_genome_effect_terms()` as the way to write such effects by hand. Without this,
   today's routing would silently store them in `phenotype_var_comp` as a random effect.
   A name moves to `GENETIC_EFFECT_NAMES` when a generator can calibrate it.
6. **Sex codes stay `M` / `F`.** They are the one deliberate exception to rule 1 (an
   external standard code).

**What else keeps its name.** `ind_tgv` / `tgv_value` stay (TBV / EBV / TGV is standard
shorthand and pairs with `ind_ebv`). `trait_var_comp` stays (renaming it touches
everything for no gain). `ind_ebv` stays a separate table: it stacks repeated
evaluations (`eval_number`, plus user columns such as `eval_date`) for accuracy and
bias studies such as the LR method, while `ind_tgv` holds one current value per
individual × trait × component.

**Blast radius** (as of 0.71.0; recount at step 1): `gen_add` appears in 6 R files
(`define_additive_effects.R` 9×, `define_effect_cov_matrix.R` 4×, `define_trait.R`,
`schema.R`, `blupf90_helpers.R`, `add_phenotype_stages.R`) and 15 files across
`R/`, `tests/`, `man/`, `vignettes/`. `order1_*` appears in 10 files. Pre-1.0 rule: the
old strings go completely. No alias, no accepting `"gen_add"` as input.

**When.** Two halves, two releases. The **variance** strings (`gen_add` → `additive`,
`epistasis` → `additive_by_additive`) go in the rename-only release (step 1, §10), next to
the writer rename. It is a pure string change, and Part A's §6C filters and error messages
must name the new values from the start. The **behaviour** that goes with them (the
reserved constants `GENETIC_EFFECT_NAMES` / `GENETIC_EFFECT_NAMES_FUTURE`, refusing future
and derived names) is Part A (step 2). The **value** names (`order1_*` → `additive` /
`dominance` / `indicator`) go in the consolidation release (step 3), which already rewrites
every reader of `ind_tbv` / `ind_tgv`. None of them waits for Part C.

Proposed addition to CLAUDE.md "Naming Consistency Rules": *6. Categorical values stored
in rows are lowercase snake_case full words: no single letters, no symbols, no
abbreviations. Standard external codes (`M`/`F`) are the exception. Values that share a
namespace with user-chosen names (`effect_name`) are reserved and listed in one
constant.*

---

## 6C. Generation targets — one rule for both generators *(decided 2026-10-01; was Q15)*

`trait_var_comp` is the **single source of truth** for generation targets (Q15 option (a)):
every target a generator calibrates to comes from it, and a generator never overwrites it.
It holds the **available** targets, not the active model *(decided 2026-10-02, Codex
review finding 1)*. A call may calibrate to a subset of the stored blocks
(`trait_var_comp_tbl`, below), and the excluded rows stay in the table. So a stored
block does not prove the model has that component. What the model *is* is read from its
terms (principle 6), and the target-vs-delivered join (§8) compares only blocks the model
has. Rejected: treating the table as the active specification, where filtering a block
out would delete or deactivate it. That reverses Q15's "nothing deleted".

**One entry path.** `define_trait()` loses `target_add_var` and `target_add_mean`
(breaking, no alias). `target_add_var` wrote 1×1 pieces with no covariances, which cannot
be combined into a multi-trait block (`load_trait_cov()` returns `NULL` when an off-diagonal
is missing), and it would bypass the refusal below. `target_add_mean` is written to
`trait_meta` and read nowhere in `R/`. After P2 the mean is the intercept in
`phenotype_meta.mean` (§6A). The `trait_meta` column goes too. Targets enter only through
`define_effect_cov_matrix()` or a generator's matrix argument. The CLAUDE.md "Two-Layer
Phenotype Design" paragraph and the `tidybreed-api` / `tidybreed-schema` skills are
updated in the same change.

**`define_trait_simple()` is removed** (decided 2026-10-02). It is an exported wrapper that
forwards `target_add_var` to `define_trait()` (`R/define_trait_simple.R:53,70`), so it would
break here anyway. Users chain the three calls it wrapped:
`define_trait()` → `define_additive_effects(G = …)` → `define_phenotype()`. Its R file,
`NAMESPACE` entry and man page go, and its callers are updated: 7 non-plan files, plus `man/`, the `tidybreed-api` skill,
`README.md` and `package_summary.md` (recount at step 2). The internal
`write_trait_var_diag()` (`R/define_effect_cov_matrix.R:246`) served only `define_trait()`'s
`target_add_var` and goes too.

**Targets are stored at full precision.** Today `define_effect_cov_matrix()` writes
values with `format(x, scientific = FALSE)` (`R/define_effect_cov_matrix.R:143, 258`),
which keeps 7 significant digits: 1/3 is stored as `0.3333333`. A calibrated model then
misses its own stored target at ~1e-7, and every 1e-10 gate fails. Part A writes targets
through a parameterised insert (or `sprintf("%.17g")`), which round-trips a double exactly,
inside one transaction with the rest of the call. This uses `dbExecute()`, never
`dbWriteTable()` (CLAUDE.md, RNG).

**Passing a matrix writes it, and never overwrites.**
- `G` (`define_additive_effects()`), and `G_A` / `G_D` / `G_AA` (`define_genome_effects()`), accept a
  matrix, or a single number when there is one trait (taken as 1×1).
- A passed block is validated (§7.1: PSD, dimnames) and then written to `trait_var_comp`
  in the same transaction as the terms.
- If **any** row already exists for that block's key, meaning `effect_name` × any of the
  call's traits × `line_name` (NULL-safe), the call is a **hard error**. This holds even
  when the passed matrix is identical to the stored one. The check runs on the whole
  table, not on `trait_var_comp_tbl`, so filtering the old rows away does not let a second
  matrix in next to them. The error names the stored block and gives the exact removal
  call, which works today because `trait_var_comp` has a row key:

  ```r
  get_table(pop, "trait_var_comp") |>
    filter(effect_name == "dominance", trait_name_1 %in% c("ADG", "BF")) |>
    remove_rows()
  ```
- `define_effect_cov_matrix()` gets the same refusal. Today it silently deletes and
  rewrites (`R/define_effect_cov_matrix.R:127-132`), which would make the generators'
  check pointless. This mirrors `phenotype_var_comp`, whose blocks are declared whole and
  then locked (CLAUDE.md, covariance blocks).

**Reading stored targets: `trait_var_comp_tbl`.**
- The argument takes `get_table(pop, "trait_var_comp") |> filter(...)`. The name says
  which table to filter. Any other table is an error.
- `NULL` (the default) means the stored rows for the call's traits and line (below).
- One argument serves all blocks, because `effect_name` is a column. **The `effect_name`
  values that survive the filter decide which blocks are in the model.** "Additive +
  dominance only, ignoring a stored A×A target" is `filter(effect_name !=
  "additive_by_additive")`, with nothing deleted. That settles absent vs stored without an
  extra argument.
- **A stored block is used whole, or explicitly not at all** *(decided 2026-10-02, Codex
  review finding 1)*. With `trait_var_comp_tbl = NULL`, if a resolved row pairs one of the
  call's traits with a trait **outside** the call (e.g. a stored `additive` block for
  `ADG` and `BF`, and a call with `trait_name = "ADG"`), the call errors. Reading the
  `(ADG, ADG)` row alone as a 1×1 block would recalibrate `ADG` and silently break the
  stored `ADG`–`BF` covariance. The error names the block's traits and gives the two fixes:
  pass all of them in `trait_name`, or pass `trait_var_comp_tbl` filtered so that no row
  leaves the call's traits (an explicit choice to calibrate `ADG` alone). The block is
  found from the rows themselves (the traits linked by off-diagonal rows within one
  `effect_name` × `line_name`); no block id is stored. The check runs before any write
  or draw.
- After filtering, the rows for the call's traits must form exactly **one complete,
  symmetric k×k block per `effect_name`**. A partial block (e.g. `filter(trait_name_1 ==
  "ADG")` drops the `(BF, ADG)` rows) is an error. So are two candidate sets for one
  `effect_name`; the error asks the user to filter to one.
- The same block passed as a matrix **and** present in the table is the "already stored"
  error above.
- A zero matrix is an exact constraint, stored as zeros. A block with no rows is absent.

**Which blocks each generator accepts.**
- `define_additive_effects()` uses `additive` only. If the resolved rows also hold
  `dominance` or `additive_by_additive` for the call's traits, it errors. A stored
  non-additive target is never silently ignored by the default call. **From step 2** the
  error names `trait_var_comp_tbl` filtered to `effect_name == "additive"`, and the
  `remove_rows()` call. `define_genome_effects()` does not exist until step 5, and no
  message before then names it. **Step 5** adds it as a third fix.
- `define_genome_effects()` needs an `additive` block; `dominance` and
  `additive_by_additive` are optional. With no `additive` block it errors: dominance and
  epistasis imply a floor on the additive covariance (§1.3), so there is no model without
  it. **An additive-only target is accepted** and gives the same rows as
  `define_additive_effects()` with the same seed, anchor and `base_tbl` (gate C4). The two
  differ in their replacement rule (§5) and in the options only
  `define_additive_effects()` has (`line_name`, `parent_origin`, `method = "union"`,
  `distribution`). The roxygen says which to use: `define_additive_effects()` for scoped or
  crossbred additive models, `define_genome_effects()` when dominance or A×A may be added.

**Lines (Q19 (a), decided 2026-10-01).** `trait_var_comp` gains `line_name VARCHAR`
(NULL = all lines), added with `ALTER TABLE` in Part A. The block key becomes
`(effect_name, line_name, trait_name_1, trait_name_2)`. A generator's passed matrix is
written with the call's `line_name`. Default resolution with `trait_var_comp_tbl = NULL`
follows the package's shared/default pattern:
- a common-scope call reads `line_name IS NULL` rows;
- a `line_name = "C"` call reads line C's rows if that block exists, and otherwise falls
  back to the `NULL` rows. The fallback is decided **per `effect_name` block**: line C
  can have its own `additive` block and share the population-wide `dominance` block. A
  *partial* line-C block (some of the call's traits, or one triangle) is an error, not a
  fallback. The non-additive refusal above looks at the same resolved rows.

So lines A and B can share the population-wide set while line C has its own, with no
duplicated rows. `trait_var_comp_tbl` overrides the default (e.g. line C deliberately
using the shared set). `define_effect_cov_matrix()` gains `line_name = NULL`.
The internal readers `load_trait_cov()` and `get_trait_var()`
(`R/define_effect_cov_matrix.R:178, 216`) gain the same line resolution. Without it they
would mix rows from different lines. Their other callers read the population-wide rows
explicitly (`line_name IS NULL`): `add_ebv()`'s BLUPF90 parameter file
(`R/blupf90_helpers.R:326`; Q20 later separates evaluation parameters) and the
prevalence threshold (§6A). The column is registered in `R/schema.R` (`.sm_col()`) and in
`R/sql_utils.R`'s column registry. `archive_replicate()` copies it with the table.
`extract_genetic_variance()` gets no `line_name` column. It measures whatever individuals
`tbl` selects, and that selection has no line of its own. To compare with a line's target,
filter `trait_var_comp` to that `line_name` first; the §8 join is then unchanged.

The biology limits what this is for (§7.6). A per-line target is meaningful only for
**line-specific** effects. One set of common effects can be calibrated to one population's
covariance only. Other lines get whatever their frequencies give, and that is measured,
not targeted.

---

## 7. Part A — `anchor =` on `define_additive_effects()`

### 7.1 API

Signature from step 3. Step 2 still has `effects` and `scale_to_target` (see below).

```r
define_additive_effects(tbl, trait_name,
  distribution    = c("normal", "gamma"),               # see Q5 on `architecture`
  G               = NULL,                               # matrix, or a number for one trait; written, never overwrites (§6C)
  trait_var_comp_tbl = NULL,                            # NEW; filtered trait_var_comp; NULL = stored rows (§6C)
  anchor          = c("genic", "realised"),             # NEW; "genic" = today's behaviour
  method          = c("shared", "union"),
  base_tbl        = NULL,                               # also names the realised population
  line_name       = NULL, parent_origin = NULL,
  warn_bounds     = c(0.8, 1.25),                       # NEW; NULL turns the §7.4 warnings off
  seed            = NULL)
```

- **When the result is exact.** The requested covariance is delivered exactly **only**
  with `method = "shared"` and a feasible
  rank ($\mathrm{rank}(\mathbf G)\le\mathrm{rank}(\mathbf M)$ and
  $\le\mathrm{rank}(\mathbf B_0^\top\mathbf M\mathbf B_0)$, otherwise one of two distinct errors, §1.1, A5). This
  sentence goes at the top of the roxygen `@details`. The closing `message()` says
  "exact" or "approximate", and gives the delivered covariance under the anchor.
- **No manual `effects`, no `scale_to_target`** *(decided 2026-10-03, Q21 (a); removed
  in step 3)*. The generator always samples and calibrates, like `define_genome_effects()`
  (§9.1). Exact, user-chosen coefficients (GWAS estimates, a published QTL map) are
  written with `define_genome_effect_terms()` and `ad_terms()`, under a user owner. That
  makes "`generated`" mean "calibrated to the stored target", which the prevalence
  threshold relies on (§6A). The earlier plan to make `effects` a named table is
  withdrawn. Until step 3 the old arguments stay as they are, because `add_tbv()` reads
  only `generated` terms before consolidation (§5). Step 2 adds only the refusal of `G`
  together with `effects` or `scale_to_target = FALSE` (A19): nothing would be
  calibrated, so `trait_var_comp` would record a target the model does not deliver.
  **Target resolution on the two surviving paths (step 2 only):**
  - with `effects`, no target is resolved or required, and the §6C block checks
    (partial block, non-additive block, missing target) are skipped;
  - with `scale_to_target = FALSE`, a stored target is read only where today's code reads
    it, as the `Sigma` of the k ≥ 2 joint draw (`R/define_additive_effects.R:417`), and
    no congruence is applied. A missing block errors only on that k ≥ 2 path.
- **Targets are validated before any write.** A passed or stored `G` must be finite,
  symmetric and positive semidefinite (the source's `.qtl_psd_eigen()` relative tolerance).
  An indefinite `G` errors, naming the block and the most negative eigenvalue. The same
  check goes into `define_effect_cov_matrix()`, which today checks only symmetry and
  diagonals. A zero-rank target and a rank-deficient but feasible target are distinct
  tested cases (A11).
- **Dimnames are checked, never overwritten.** Today `define_additive_effects()`
  (`:327`) and `define_effect_cov_matrix()` (`:113`, always: with no `trait_names` it copies the
  row names onto the columns) assign dimnames over whatever the matrix carries. So a `G` named `c("BF", "ADG")` passed with
  `trait_name = c("ADG", "BF")` is stored with its rows the wrong way round. A named matrix
  must match `trait_name` in order, or the call errors. An unnamed matrix is taken in
  `trait_name` order. The same rule applies to all three Part C targets (C18).
- **Atomic.** Today a passed `G` is written to `trait_var_comp` (line 328) before effect
  validation and `.ge_commit()` (line 485), so a failed call can change the target and keep
  the old coefficients. That is an existing bug. Fix it in Part A: resolve and validate
  every target, build every trait's terms, and only then commit targets and terms
  together. (The §5 owner-content replacement check joins this sequence in step 5, with
  `define_genome_effects()`.) Target resolution, the write of a passed `G`, and the
  "already stored" refusal follow §6C.
- **One internal target insert.** A passed `G` is no longer written by calling the
  exported `define_effect_cov_matrix()` (today `:328`). Both writers call one internal,
  transactional, full-precision insert (§6C), which carries the "already stored" refusal.
  Step 3's Q21 refusal (no new block under existing `generated` terms, §6A) goes only in
  the exported `define_effect_cov_matrix()`, never in that shared insert. Otherwise
  `define_additive_effects(G = )`, the documented way to set a new target and regenerate,
  would refuse itself.
- **`seed` after validation.** Today `set.seed(seed)` runs before any validation
  (`:216`), so a refused call with `seed =` changes `.Random.seed`. Move it after every
  check, immediately before the first draw (A22).

- `"reference"` needs no argument of its own. `base_tbl` **already is** the reference
  selection, and the genic weights are computed from whatever it names. The manuscript's
  `reference` anchor is therefore `anchor = "genic"` plus a non-default `base_tbl`. Say so
  in the roxygen.
- `anchor = "genic"`: weights $w_j=n_{\text{eligible},j}\,p_jq_j$ (the same quantity
  `rescale_effects_to_target()` uses today). For $k\ge2$ the per-trait scalar is replaced
  by the $k\times k$ congruence. **For $k=1$ the result is the scalar rescale** (gate A1, checked within the current code).
- `anchor = "realised"`: $\mathbf M=\mathrm{Cov}(\mathbf X)$ of the individuals `base_tbl`
  selects, via `extract_genotypes()`. The congruence uses the centred design, $O(nmk)$.
  Document it as the single-generation / clonal option (manuscript Recommendation 3).
  First-release restrictions, each an error that gives the reason:
  - `base_tbl` must select **individuals**, meaning a table with `id_ind` other than
    `ind_haplotype`. The `ind_haplotype` shape means "these allele **copies**". A filter
    such as `line_origin == "Duroc"` leaves crossbreds with partial genotypes, and
    $\mathrm{Cov}(\mathbf X)$ of partial genotypes is not a population covariance.
    `founder_haplotypes` has no individuals (see Q11 for a pool-based alternative).
  - Common scope only (`line_name = NULL`, `parent_origin = NULL`). A scoped
    realised anchor needs the covariance of the *eligible* copies (paternal only,
    line-L only), not of dosages. That is well defined, but no one needs it yet.
  - $\mathbf M$ is formed in R from the collected, `id_ind`-then-`locus_id`-ordered
    matrix. Never form it with a SQL `SUM()` over doubles (CLAUDE.md, bit-identical
    rule). Allele-frequency counts from `extract_allele_freq()` are integer sums and
    are fine.
  - **Size guard.** Collecting an $n\times m$ genotype matrix defeats the larger-than-RAM
    design. The first release errors above an internal cell limit. The error gives $n$,
    $m$ and the limit, and suggests `anchor = "genic"` or a smaller `base_tbl`. A later
    bounded-memory path accumulates $\mathbf X^\top\mathbf X$ block by block in a fixed
    `id_ind` order, which keeps it deterministic. A test checks that the same filter gives
    the same row and locus order after `restore_pop()`.
- `"marginal_observed"` and `"dual"` stay out of the user API. Keep them in the ported
  internals, untested from the public surface (Q6).

### 7.2 Algorithm inside the generator

1. Draw $\mathbf B_0$ exactly as today (`mvrnorm` rows, or the single-trait draw). The draw
   is now only the *architecture*. Part A moves it into one internal helper,
   `.draw_additive_architecture()`, which `define_genome_effects()` calls too, first,
   before any dominance or pair draw (§9.2). That keeps an additive-only
   `define_genome_effects()` call identical to `define_additive_effects()` (gate C4).
2. Build $\mathbf M$ from the anchor.
3. $\mathbf B\leftarrow\mathbf B_0\mathbf A$ by the congruence (`.qtl_congruence()`, ported).
4. Write through the existing path (`.dae_build()` → writer engine), with target writes
   in the same commit (§7.1, Atomic).
5. Report diagnostics (§7.4).

**RNG.** For $k=1$ with `"genic"`, no extra RNG draw, and the values are the scalar
rescale's (gate A1, checked inside the current code, not against an old version). For $k\ge2$ the
draw is unchanged too, so the same seed gives the same $\mathbf B_0$. Only the rescale
changes, so values change by design. `NEWS.md` must say so.

### 7.3 `method = "union"` — the one real conflict

The $k\times k$ right factor mixes columns. A trait-specific zero pattern (a locus that is a
QTL for trait 1 only) is **destroyed**: after $\mathbf B_0\mathbf A$ every locus in the union
affects every trait. Options in Q4. Recommended for the first release: under `"union"`,
keep today's per-trait scalar, which preserves the distinct QTL sets. When the target has
a nonzero off-diagonal, issue a **warning** (not a message) that says "approximate" and
gives the delivered covariance and correlation next to the target, and names
`method = "shared"` as the exact option. Rejecting non-diagonal targets under `"union"` was
considered. It would make an approximation users accept today impossible, for no gain in
safety once the warning is explicit.

### 7.4 Diagnostics

Port the additive file's per-contrast relative spectrum and `warn_bounds`, exposed on
both generators (default `c(0.8, 1.25)`, `NULL` = off). The source's separate
`warn_on_mismatch` flag only gates the same warning, so it is not ported. Warn when the
covariance the **other** population sees departs from the target by more than those
bounds. That is the manuscript's Recommendation 6: the anchor you choose is exact, and the
warning tells you what the other population sees. Each warning names which reference was
**observed** (a genotype covariance of real individuals) and which is an **expectation**
(the genic limit, or the pool expectation below).

**The comparison population at generator time.** The usual order is
`define_founder_haplotypes()` → `define_*_effects()` → `add_founders()`, so there are
usually no individuals yet. When `base_tbl` is `founder_haplotypes`, the comparison is the
**pool expectation under random pairing** of haplotypes,
$\mathbf M_{\text{pool}} = 2\,\mathrm{Cov}(\mathbf H)$, which keeps the pool's LD. The
covariance uses the **population divisor** $n_h$: `add_founders()` draws each haplotype
independently and with replacement, so a founder's two copies are iid draws from the
empirical pool (0.73.2; the sample divisor $n_h - 1$ overstated it by $n_h/(n_h-1)$). It is
**not** the covariance of the founders that `add_founders()` will draw. A finite sample,
or a structured pairing, gives something else. Label it "pool expectation" in messages,
never "founders". When `base_tbl` selects individuals, use their $\mathrm{Cov}(\mathbf X)$
and label it "observed".

*(Q22, decided (b), 0.73.1: the pool comparison is always a `message()`, with a hint
to calibrate `anchor = "realised"` on the founders when it leaves `warn_bounds`; only the
observed and genic-limit comparisons `warning()`.)*

**Storage: none.** Diagnostics are printed (`message()` / `warning()`) and not stored: not
as an attribute on the pop (CLAUDE.md: never in the R object), and not as rows in
`trait_var_comp` (a table of *targets* that generators read back, principle 6). This is a
choice, and it has a cost: not every diagnostic can be rebuilt later. The stored terms
keep coefficients and centres. They do **not** keep the arbitrary `base_tbl` filter, and
they do not keep the founder pool's state at generation time. `extract_genetic_variance()`
(§8) measures the **current** named selection. A user who wants to compare against the
generation-time population later must save that population's selection (e.g. its
`id_ind` set). No shadow metadata table is added for this.

### 7.5 Files

Paths relative to `$SRC`; full inventory in §1A.

| Source (`simulate_qtl_effects/`) | tidybreed |
|---|---|
| `R/qtl_effects.R`: `.qtl_psd_eigen`, `.qtl_relative_spectrum`, the `.congruence` closure in `sim_qtl_effects()`, validation helpers | new `R/qtl_congruence.R` (internal, `@noRd`) |
| `tests/test_qtl_effects_paper.R` (not `test_qtl_effects.R` — see §1A): congruence, rank-deficiency (Codex F2), Theorem 1 error, validation | `tests/testthat/test-qtl-congruence.R` |
| `arch_*()` samplers, `marginal =` | **not in Part A** (Q5) |

Port from `R/qtl_effects.R`, **never** `R/qtl_effects_paper.R`, which is frozen to the
manuscript.

### 7.6 Lines, the reference population, and line means *(discussed 2026-09-24/25)*

**Line differences live in the DNA.** A difference in mean genetic value between lines
comes from allele-frequency differences at QTL whose effects are shared by every line.
Under additive effects centred on reference line R,
$\bar g_L-\bar g_R=\sum_j 2\alpha_j(p_{Lj}-p_{Rj})$, and an F1 averages its parent lines.
There is **no** line-mean input (a fixed breed-composition mean was considered and
rejected: it cannot respond to selection and does not segregate). The documented
workflow:

1. **Burn-in.** Start from one founder pool (p = 0.5 is fine here; it is only a
   starting point) and **common** effects (`line_name = NULL`). This phase can run in
   tidybreed or in a faster tool such as AlphaSimR.
2. **Divergence.** Split the pool into lines and select each on its own index over many
   generations. For example, maternal lines go up for litter size and the terminal line
   goes for growth, with restricted indexes computed outside tidybreed to hold traits that
   should not move. Selection is **by filtering tables** and calling `add_offspring()`.
   There is no `select_parents()`, by design. The lines end with realistic allele
   frequencies, LD, and their own (smaller, selected) variance components.
3. **Main simulation.** Carry the generation-*t* animals, their haplotypes and the
   **same** effects (with their `center_value`s) into the main run. Choose a reference
   line (e.g. the terminal line) and set `phenotype_meta.mean` to its target phenotypic
   mean minus its mean `ind_tgv_total` (§6A: `mean` is an intercept). The other lines'
   means then follow from their genotypes.

Each line's covariance matrices **emerge** from selection and drift. They are not inputs
to the effect generator. They are measured by `extract_genetic_variance()` on that
line, and a breeder's *estimates* of them are evaluation parameters (Q20), not generation
targets. Tools that set line frequencies directly to target mean differences, and a
supported import of one database's end state as another's founders, are deferred
(`plans/TODO.md`).

**Common effects centred on a reference line** work today:
`define_additive_effects(tbl, ..., base_tbl = get_table(pop, "founder_haplotypes") |>
filter(line_name == "Terminal"))`. The calibration then hits the target **within that
line**. The other lines get whatever variance their frequencies give.

**Line-specific effects** (`line_name = "A"`) are a different model: the same allele acts
differently in each genetic background. By default each variant is centred on its own
line's base, so those terms contribute **no** mean difference between lines. Nothing says
so today. Add a `message()` on every line-scoped call: *"Line-specific effects are
centred on line A's own base, so they add no difference between line means. For
differences between lines, use common effects centred on one reference line (see
?define_additive_effects)."* Put the same text in the roxygen. It is a message, not a
warning, because the model is legitimate. Gate A13.

**Per-line covariance targets** are decided in §6C (was Q19).

---

## 8. Part B — `extract_genetic_variance()`

v4.9 review §3.2 named this "the first thing to build after this plan" and put it ahead
of consolidation. It is **not** a straight port of `nonadd_decompose()`. That function
takes one diploid dosage matrix and dense functional $\mathbf B_a,\mathbf B_d,\mathbf B_{aa}$.
tidybreed's writer allows several owners, line and parent-origin scopes, indicator
surfaces and higher-order terms. `.gev_evaluate()` evaluates all of them, but it cannot
give a unique functional $(a, d, e)$ for every one. The function therefore declares which
models it can decompose.

```r
extract_genetic_variance(tbl, trait_name = NULL,
  base_tbl = NULL,             # genic only: frequencies; NULL = the selected individuals
  anchor   = c("realised", "genic"))
```

- `tbl`: a `tidybreed_table` selecting individuals (as `add_tgv()` takes).
- `base_tbl`, not `reference_tbl`: like the generators' `base_tbl`, a selection that
  defines $p$ via `extract_allele_freq()`. Two differences are stated in both roxygen
  blocks:
  - Its `NULL` default is **the selected individuals**, not the founder pool. That default
    does **not** reproduce a generation target once frequencies have drifted. To compare
    with a target, pass the generation base (e.g. `founder_haplotypes`).
  - It is **genic only**. Under `"realised"` the projection uses the frequencies and
    regressions of the individuals `tbl` selects, so a non-`NULL` `base_tbl` is an error,
    not silently ignored.
- **Read-only.** Evaluates through `.gev_evaluate()`
  (`R/genome_effects_eval.R:633`), which returns values without writing `ind_tgv`. An
  `extract_*` function must not change simulation state.
- **No `per_ind`.** A switch that sometimes returns covariance rows and sometimes
  individual values is two functions. The generation-t breeding value (§4.3.3) becomes a
  separate, later `extract_breeding_values()`, defined for supported case 1 only. A
  generation-t breeding value for a scoped model needs its own reference-population
  projection, which nobody has specified.
- **Size guard** as in §7.1: the realised path collects genotypes, so it has the same
  cell limit and the same error.
- **Cohort.** Every row of the output is computed on the same individuals: those `tbl`
  selects. `.gev_evaluate()` omits an individual that contributes nothing to a trait
  (`R/genome_effects_eval.R:629`). If any selected individual is absent for any of the
  call's traits, the call errors, giving the count per trait and asking for a narrower
  `tbl`. It never silently switches to per-trait or pairwise-complete subsets. So `n_ind`
  is one number for the whole result. Fewer than two individuals is an error (a sample
  covariance needs two). A locus with zero observed variance in the cohort contributes
  nothing to any realised block: its regression $b_j$ is set to 0 rather than divided
  by zero, and the same holds for a pair with a monomorphic member.

**Two anchors, two definitions.**

- `anchor = "genic"`: the HWE + LE orthogonal decomposition at `base_tbl`'s
  frequencies. The blocks are orthogonal by construction, so for a fully covered model
  the block covariances sum to the genic `total`.
- `anchor = "realised"`: the observed covariances of the selected individuals. Under LD or
  departures from HWE the blocks correlate, and
  $\mathrm{Var}(g)=\sum_b\mathrm{Var}(b)+2\sum_{b<b'}\mathrm{Cov}(b,b')$. The cross-block
  terms are **reported**, never folded into a remainder. The source's
  `nonadd_decompose()` returns `cov_A_D` and `cov_A_AA` but **not** a D–A×A covariance,
  which is also nonzero under LD or non-HWE. The port computes every block pair,
  including $\mathrm{Cov}(D, AA)$, in both orientations ($\mathrm{Cov}(b_{t_1}, b'_{t_2})$
  is not symmetric in the traits). They go into **one
  `between_components` row per trait pair** (Q16, decided): for traits $t_1, t_2$ its value is
  $\sum_{b\ne b'}\mathrm{Cov}(b_{t_1}, b'_{t_2})$ over every ordered pair of different blocks,
  `unpartitioned` included. So the block rows plus `between_components` sum exactly to the
  observed `total` for every trait pair, and a test asserts it (B3). Under `"genic"` the
  value is 0 by construction and the row is omitted.

**Supported cases** (dispatch on the stored model per trait, before projecting):

1. **Common-scope diploid A + D + A×A**, generated or hand-written, no line or origin
   scope. Accepted term shapes, at any centre: order-one `additive`; order-one
   `dominance`; order-one `indicator` on the diploid heterozygote
   (`copy_count_value = 2`, `dosage_value = 1`, which is what
   `ad_terms(coding = "functional")` writes for $d$); two-member `additive × additive`.
   Several owners are summed first. These canonicalise to functional $(a, d, e)$ (the
   inverse of §3, Q13) and are re-projected on NOIA (realised: observed $b$; genic:
   $q-p$ at `base_tbl`'s frequencies). Both codings of the same model therefore give the
   same report.
2. **Scoped additive-only** (all terms order-one `additive`, at least one with a line or
   origin scope: line-specific, parent-of-origin, crossbred fallback). A common-scope
   additive-only model is case 1, which is checked first. Then $g$ is `additive` plus a constant, and the
   function reports the covariance of the evaluated `additive`, labelled as the
   **evaluated** additive variance. This is **not** a NOIA re-projection at the new
   frequencies: a line-scoped crossbred can carry different coefficients and centres per
   inherited copy. It covers the main crossbreeding model today.
3. **Anything else** (indicator surfaces, A×D, D×D, scoped non-additive terms, several
   owners whose terms mix these classes): `total` is always reported, because the
   evaluator gives it. The covered blocks are reported where they can be separated.
   Terms outside the supported set go into `unpartitioned`, which is the variance of
   **those terms' value only**. Their covariance with the covered blocks goes with the
   other cross-block terms in `between_components`. A `decomposition` column says `"full"`, `"additive_only"`
   or `"partial"`. Do not refuse.

**`"genic"` is case 1 only.** The genic report is an expectation over HWE + LE genotypes at
`base_tbl`'s frequencies. The evaluator gives realised values, not that expectation. For a
scoped model the expectation also needs per-line and per-origin frequencies, which nobody has
specified. So `anchor = "genic"` on a case 2 or case 3 model errors. The error says that
`"realised"` is available, and that a common-scope additive-only model is case 1 and works.

**Output.** A tibble shaped like `trait_var_comp`:
`(effect_name, trait_name_1, trait_name_2, cov_value, n_ind, decomposition)`. The
`effect_name` values come from the **same vocabulary** as `trait_var_comp` (§6B), plus
`total`, `unpartitioned` and `between_components` (realised only, §8 above).
`inner_join(targets, realised, by = c("effect_name", "trait_name_1", "trait_name_2"))`
is then the whole target-vs-delivered check, which is the test loop this work will be
verified with. **A block row is returned only when the model has terms of that kind**, so
a stored target the generator was told to leave out (§6C, "available targets") drops out
of the `inner_join` instead of being compared with a block the model never had; an
`anti_join()` lists such targets. The comparison is a like-for-like check only when the
trait's model is the `generated` owner alone. Custom-owner terms sum with generated ones
(§5), and the extractor measures the whole model, so with custom terms present a mismatch
is expected and is not an error. Say so in the roxygen; no owner filter is added. (`unpartitioned` is not `residual`: in this package that word means
environmental noise in `phenotype_var_comp`.)

It is the measurement AlphaSimR gets wrong under epistasis: `calcGenParamE` has a
locus-index typo, found in the source project and not yet reported upstream.

Files: `non-additive/R/qtl_effects_nonadd.R` `nonadd_covariates()`, `nonadd_decompose()` →
`R/extract_genetic_variance.R`, plus the dispatch and canonicalisation.

---

## 9. Part C — `define_genome_effects()`

### 9.1 API

```r
define_genome_effects(tbl, trait_name,            # tbl = get_table(pop, "genome_meta") |> filter(...)
  G_A = NULL, G_D = NULL, G_AA = NULL,                 # matrices (or numbers for one trait); written, never overwrite (§6C)
  trait_var_comp_tbl = NULL,                           # filtered trait_var_comp; surviving effect_names = the model (§6C)
  pairs      = NULL,                                   # data frame (locus_name_1, locus_name_2); NULL with G_AA => random pairing (AlphaSimR)
  n_pairs    = NULL,                                   # random pairs; NULL = floor(m/2), every QTL once (Q8)
  anchor     = c("genic", "realised"),
  dominance_degree_mean = 0.19,                        # dominance degrees for the B_d architecture
  dominance_degree_sd   = 0.097,
  inbreeding_depression = NULL,                        # per-trait; needs a dominance block; exact for k = 1 or diagonal G_D
  base_tbl   = NULL,
  warn_bounds = c(0.8, 1.25))                          # as §7.4; NULL = off
```

- **Samples and calibrates only.** There is no `effects` input. Users with exact
  coefficients use `define_genome_effect_terms()` with `ad_terms()` / `aa_terms()` (§9.3, with a
  worked example). "Write these coefficients unchanged" and "calibrate this architecture"
  are different operations, and one argument must not mean both. A supplied architecture
  ($\mathbf B_a,\mathbf B_d,\mathbf B_{aa}$ as starting points) waits for Q5.
- **Any subset of blocks that includes `additive`** (§6C). With an additive-only target the
  calibrator reduces exactly to Part A's congruence (§1.3), and the call gives the same
  coefficients as `define_additive_effects()` under the matching settings (gate C4). So a
  user can start additive-only and later add a dominance target, then re-run the same
  function. Re-running replaces the trait's whole `generated` model (`mode =
  "replace_owner"`, §5).
- **Pairs** (Q8, decided 2026-10-02). Two routes, both used only when an A×A block is
  present. `pairs` or `n_pairs` without one is an error, and so is passing both.
  - **Supplied `pairs`** (`locus_name_1`, `locus_name_2`): the user's explicit design.
    **A locus may appear in several pairs**, e.g. a hub gene interacting with five
    others. The method accumulates them exactly (§1.3; `qtl_effects_nonadd.R:245-246`),
    and the genic anchor stays diagonal, since pairs sharing a locus are uncorrelated under
    HWE + LE. Unknown loci, loci outside the filtered QTL set, self-pairs and repeated pairs
    in either order are errors that name the keys.
  - **Random pairs** (`pairs = NULL`): a **random matching**, each locus in at most one
    pair (AlphaSimR's structure). The QTL are sorted by `locus_id`, permuted with one
    `sample()` call (§9.2 fixes where that falls in the RNG order), paired off
    consecutively, and the first `n_pairs` pairs kept. Within a pair the loci are put in
    `aa_terms()`'s canonical order (§9.3). The pairs are then sorted by
    `(locus_id_1, locus_id_2)` of that order, so the written order does not depend on the
    draw, or on the platform. The same rule applies to supplied `pairs`. With an odd $m$, one locus stays unpaired.
  - **`n_pairs = NULL` means $\lfloor m/2\rfloor$**: every QTL paired once, as AlphaSimR
    does. A `message()` says so, so the default is never silent: *"Drew 150 random A×A
    pairs: every QTL paired once (floor(300 / 2)), as AlphaSimR does. Pass n_pairs for
    fewer pairs, or pairs for a chosen design (a locus may then appear in several
    pairs)."* `n_pairs` must be an integer in $1..\lfloor m/2\rfloor$. A larger value errors,
    gives the maximum, and points to `pairs` for designs where loci repeat.
- **Targets** follow §6C. A `NULL` matrix argument means "not passed": the block comes
  from `trait_var_comp_tbl`'s rows (default: the stored `line_name IS NULL` rows), and is
  absent if there are none. Matrix dimnames must equal `trait_name`, in order. Mismatched
  dimnames are an error, never renamed. A zero matrix is an explicit exact constraint and
  is different from an absent block.

- Arguments are spelled out (`dd_mean` → `dominance_degree_mean`, `inbr_depr` →
  `inbreeding_depression`), following the no-abbreviation rule. `pairs` columns follow
  the `trait_name_1` / `trait_name_2` pattern. The pair refusals above match the source
  (`qtl_effects_nonadd.R:169`). A repeated pair would also be two terms with the same
  family signature and scope, which the writer refuses.
- **`inbreeding_depression`** is a named numeric vector keyed by `trait_name`. With no
  dominance block it is an error. It is a calibration argument for this call and is **not**
  stored: it has no `trait_var_comp` row, and the implied depression is reported (§9.2 step 6).
  The stored $d$ keep the model reproducible after `restore_pop()`.
- No `seed` argument (Q12): the caller uses `set.seed()`.
- **The port takes `p`, not `X`, on the genic path.** `sim_qtl_effects_nonadd()`
  computes `p <- colMeans(X)/2`. In tidybreed the default `base_tbl` is the founder
  pool, and there are no individuals when effects are defined. So `.na_sim()` takes
  `p` (from `extract_allele_freq()`) for `"genic"` and `X` only for `"realised"`. The
  same split applies to `.qtl_congruence()`'s callers in Part A.
- **Crossbreeding.** Common scope still produces heterosis. With one `d_j` per locus,
  F1 dominance comes from allele-frequency differences between the parent lines,
  which is the classical model. Default `base_tbl` resolution and the Wahlund warning
  follow `define_additive_effects()`: the pooled founder table is centred on purpose,
  and the targets then hold at the **pooled** frequencies, not within each line. Say
  so in the roxygen.

- Same pipe subject, `base_tbl` semantics and default base as `define_additive_effects()`.
  A second verb in the generator family, per `update_genome_effects_base_tbl.md` §2.5.
- Common scope only in the first release. `line_name` / `parent_origin` are not arguments,
  and the roxygen says so (and points to Q14).
- Refuses QTL at non-diploid or non-autosomal loci (reuse `assert_qtl_autosomal()`). The
  source method's `X ∈ [0,2]` validation already enforces diploidy, and the writer rejects
  `dominance` at non-diploid loci anyway.

### 9.2 Steps

1. Resolve loci, traits and targets (§6C), and validate any supplied `pairs`. Validate
   each target as a PSD matrix (§7.1) and count the `generated` terms the call will
   replace, for the §5 message. Nothing is written
   and no random number is drawn yet.
2. $p$ from `extract_allele_freq(base_tbl)`. For `"realised"` only, the genotype matrix
   from `extract_genotypes()` (same restrictions as §7.1).
3. Draw, in this fixed order: the additive architecture with
   `.draw_additive_architecture()` (shared with `define_additive_effects()`, §7.2), then
   random pairs (only if `pairs = NULL` and A×A is present), then the dominance degrees
   and the A×A architecture (only for the blocks present). Then call the ported
   `.na_sim(p, X = NULL, ...)` (the body of `sim_qtl_effects_nonadd()`, refactored to take
   `p` and the drawn architectures). The order is what makes an additive-only call match
   `define_additive_effects()` (C4), and adding a block never changes the additive draw.
4. **Convert to statistical coding** (§3). Under `"realised"`, re-derive $\alpha$ with
   $b=q-p$ (§4.3.1).
5. Build `terms` per trait and write through the engine with
   `effect_owner = GE_GENERATED_OWNER, mode = "replace_owner"`, all traits **and any
   passed-target writes** (§6C) in **one transaction**. A half-written multi-trait model is a
   different model, and a target that changed while its terms did not is a lie in
   `trait_var_comp`.
6. **Report, do not store** (§7.4): delivered covariances, the floor
   $\mathbf R-\mathbf Q^\top\mathbf P^{-1}\mathbf Q$, the implied inbreeding depression and
   the functional $(a,d,e)$ summary go to `message()`. The targets are already in
   `trait_var_comp`. Delivered values are recomputable from the stored terms. The floor
   only matters when it is violated, and then it is in the error.
7. Errors: $\mathbf G_A$ below the floor (with the discriminant message, naming the minimum
   feasible diagonal); rank infeasibility (Theorem 1 per component); rank-deficient
   $\mathbf B_a$ (source test 10's message).

### 9.3 New builder

`aa_terms()`: a pure builder alongside `ad_terms()`, returning the two-member rows of §3.
The generator uses it, and users writing hand-calibrated epistasis get the correct shape
instead of reaching for `genotype_terms()`. Exported.

- **Input**, mirroring `ad_terms(locus_name, a, d, p, coding, effect_name, report)`:
  `aa_terms(locus_name_1, locus_name_2, e, p_1, p_2, coding = c("functional",
  "cockerham"), effect_name = NULL, report = TRUE)`. `e` is one coefficient per pair.
  `p_1` / `p_2` are **required under both codings**, exactly as `ad_terms()` requires `p`:
  they are the centres under Cockerham and the basis of the reported mean under
  functional. All are recycled to the pairs and aligned with the locus names they sit
  next to, as in `ad_terms()`, so every value is keyed by a name.
- **Output**: two rows per pair, both `contrast_name = "additive"`, the same coefficient on
  each, centres $0.5$ (functional) or $p_k$, $p_l$ (Cockerham). Within a pair the two
  loci are put in a canonical order (by `locus_name`, with
  `order(..., method = "radix")`, i.e. C-locale byte order, so it is the same on every
  platform; each `p` is swapped along with its locus), and `term_id` is built from them. So `(L10, L44)` and `(L44, L10)`
  give the same rows, and the output is stable under input reordering.
- **Report**: under functional coding, a `message()` with each pair's share of $\mu$,
  $e_{kl}(2p_k-1)(2p_l-1)$ (§3), as `ad_terms()` does for its terms. It is written to no
  table.
- **Pure**: it never changes or rescales $e$. Pair validation is as in §9.1.
- **One caveat, stated in its roxygen.** `ad_terms(coding = "cockerham")` takes $\alpha$,
  not the functional $a$: $\alpha_j=a_j+(q_j-p_j)d_j$ already whenever $d\ne0$. With pairs,
  $\alpha_j$ also needs $\sum_l e_{jl}c_l$, which `ad_terms()` never sees. So functional
  $a$ fed to the Cockerham branch is **not** the statistical form of the model. The total is still right up to a constant; only
  the `additive` / `interaction` split changes. A manual Cockerham path is Q17.

**Worked manual A + D + A×A** (the §9 roxygen and the vignette use this). It uses
functional coding, so the total is exact, and the components are the functional ones, not
breeding values:

```r
pop |>
  define_genome_effect_terms(
    trait_name = "ADG",
    terms = rbind(
      ad_terms(locus_name = c("Locus_10", "Locus_44"),
               a = c(0.30, -0.12), d = c(0.10, 0.05), p = c(0.35, 0.60),
               coding = "functional"),
      aa_terms(locus_name_1 = "Locus_10", locus_name_2 = "Locus_44",
               e = 0.08, p_1 = 0.35, p_2 = 0.60, coding = "functional")))
```

(Signatures checked against 0.71.1, where the writer was still named
`define_genome_effects()`: `(pop, trait_name, terms, effect_owner = "custom", ...)`, and `ad_terms(locus_name, a, d = 0, p, coding, effect_name, report)`.
Base `rbind()` needs identical columns. Functional `ad_terms()` output carries
`copy_count_value` / `dosage_value` for its `indicator` rows, and a Cockerham one does not.
So every term builder (`ad_terms()`, `aa_terms()`, `genotype_terms()`) returns **one fixed
column set**, with `NA` where a column does not apply, and C16 checks that the three
outputs `rbind()` in any combination.)

### 9.4 Files

Paths relative to `$SRC`; full inventory in §1A.

| Source | tidybreed |
|---|---|
| `non-additive/R/qtl_effects_nonadd.R`: `sim_qtl_effects_nonadd()`, `congruence()` (thin-SVD form, fixed 2026-09-19), `solve_additive_stage()`, `solve_dd_mean()`, validation/diagnostics helpers | `R/genome_effects_calibration.R` (internal) |
| `nonadd_covariates()`, `nonadd_decompose()` | `R/extract_genetic_variance.R` (Part B) |
| `zeng_appendix_A()` | test oracle only |
| `non-additive/tests/test_qtl_effects_nonadd.R` (18 tests / 89 checks) | `tests/testthat/test-genome-effects-calibration.R` |

Share one congruence implementation with Part A. The source project has two with
identical algebra (`.congruence` and `congruence()`). Port it once.

---

## 10. Implementation order

**Fixed order: 1 → 2 → 3 → 4 → 5.** No swaps. Step 4 reads the `ind_tgv` value names
that step 3 introduces (`additive`, not `order1_additive`), and its gates use `aa_terms()`,
which moves to step 4 (below). Each step ships on its own with the **full** suite passing.
There is no compatibility shim between steps (CLAUDE.md, pre-1.0).

| Step | Content | Version | Depends on |
|---|---|---|---|
| 0 | Merge `feat/genome-effects-v49` to `main` | — | **done** (v0.71.1) |
| 0b | Two live bug fixes (below) | 0.71.2 | **done** (`_phase_0b.md`) |
| 1 | Rename only | 0.72.0 | **done** (`_phase_1.md`) |
| 2 | Part A + §6C targets | 0.73.0 | **done** (`_phase_2.md`) |
| 3 | Consolidation + P2 + Q18, in three commits (3a / 3b / 3c) | 0.74.0 / 0.74.1 / 0.74.2 | 2 (the `line_name` readers, §6C); 3a **done** |
| 4 | Part B | 0.75.0 | 2 (genotype collection, size guard, PSD helper in `R/qtl_congruence.R`) and 3 (value names) |
| 5 | Part C | 0.76.0 | 2, 3, 4 |

**Every step**, before its commit: bump `DESCRIPTION` `Version:` and add a `NEWS.md` entry
(CLAUDE.md). The entry says that databases written by earlier versions are not readable
(pre-1.0, no migration) whenever the step changes stored strings or DDL (steps 1, 2, 3).
Update `CLAUDE.md`, both skills (`tidybreed-api`, `tidybreed-schema`), `package_summary.md`
and any vignette the step touches. Regenerate `man/` and `NAMESPACE` with
`devtools::document()`, never by hand. Run the full suite.

**After every step, write a phase summary** at
`plans/import_qtl_effect_methods_phase_<step>.md` (`_phase_0b.md`, `_phase_1.md`, …,
`_phase_5.md`), like the `sample_correlated_effects_phase_*.md` files. It records what was
built, the files and functions changed, schema/DDL and stored-string changes, tests and
gates added (with their results), deviations from this plan and why, bugs found along the
way, and anything deferred to a later step.

**Keep this plan current while implementing.** Implementation often finds things the plan
got wrong or missed. When that happens, fix the code **and** update the affected sections
and gates of this plan in the same step, so the plan always describes what will be (or
was) built. The phase summary lists each such plan change.

### Step 0b — live bug fixes (0.71.2) *(found 2026-10-02, Codex review findings 2–3; **done** 2026-10-02, see `import_qtl_effect_methods_phase_0b.md`)*

Two bugs in today's code, independent of this plan. Fix them first, as a patch release,
so each lands with its own regression test instead of inside a large step. Steps 3 and 5
then retarget the same code (§6A), keeping the tests.

- **B-1. Group-contributor sums are not bit-identical.** `.group_mate_tbv()`
  (`R/contributor_tbv.R:136`) runs `SUM(t.tbv_value)` over mates. With three or more mates,
  DuckDB's parallel aggregation can add them in a different order per run, so
  `group_sum()` / `group_mean()` and the phenotypes built from them can differ in the last
  bits between thread counts. That breaks CLAUDE.md's bit-identical rule. Fix: accumulate
  through `GEV_ACC_TYPE` (`CAST(SUM(CAST(t.tbv_value AS DECIMAL(38, 18))) AS DOUBLE)`), and
  divide that sum by the integer `COUNT` for the mean. Test: a composite phenotype with
  `group_sum()` and `group_mean()` over groups of at least five mates is
  `expect_identical()` under `SET threads = 1` and `SET threads = 8` (the step 3 form of
  this test is PH5).
- **B-2. A prevalence threshold silently uses zero genetic variance.**
  `R/add_phenotype_stages.R:1146` calls `get_trait_var(pop, "additive", t)` with the
  **phenotype** name `t`. For a composite phenotype (`phenotype_components` or
  `formula_tbv`) no `trait_var_comp` row has that name, and for a simple trait with no
  stored additive target there is no row either. Both return `NA`, which the next line
  turns into 0, so the threshold is computed from the residual variance alone and the
  realised prevalence is wrong with no message. Fix: error in both cases, naming
  `define_phenotype(thresholds = )`. Test: a composite phenotype with `prevalence` errors;
  a simple trait with a stored target is unchanged; explicit `thresholds` works for both.
  Step 3 extends the same rule to the total genetic variance and the active-block rule
  (§6A, PH7).

Both are bug fixes. B-2 is a behaviour change: calls that used to succeed with a wrong
threshold now error. `NEWS.md` lists them under 0.71.2.

**As built (deviations).** B-2's composite refusal is in `define_phenotype()` as well as
in `add_phenotype()`, so it fails at definition time. B-1's regression test needed pens of
200: at pens of 10 the old `SUM()` never diverged, so a small fixture would pass on the
broken code.

### Step 1 — rename only (0.72.0) *(**done** 2026-10-02, see `import_qtl_effect_methods_phase_1.md`)*

- Exported writer `define_genome_effects()` → `define_genome_effect_terms()` (§0A).
- `trait_var_comp` strings `gen_add` → `additive`, `epistasis` → `additive_by_additive`
  (§6B). This includes `define_effect_cov_matrix()`'s routing vector, `get_trait_var()` /
  `load_trait_cov()` callers, `blupf90_helpers.R:326,331`, `add_phenotype_stages.R:1146`,
  and the "future" wording at `R/schema.R:331`.
- Reserved owner `generated_additive_tbv` → `generated`, `GE_ADDITIVE_OWNER` →
  `GE_GENERATED_OWNER` (§5). This includes `R/define_genome_effects.R:12`,
  `R/genome_effects_eval.R:182,812,849,881`, `R/define_additive_effects.R:290,480,507,608,640,674`,
  `R/schema.R:108`, `R/add_tbv.R:18`, and the test literals (`helper-genome-effects-db.R`,
  `test-add_tbv.R`, `test-genome-effects-writer.R`).
- `define_effect_cov_matrix(trait_names =)` → `trait_name =` (§0A).
- The active plans name the writer by its new name: `plans/TODO.md:20`.

**Step 1 is deliberately mechanical.** It changes names and nothing else, so the whole
existing suite must pass after it with only renamed calls and strings in the tests. That
makes any later failure a behaviour change, not a missed rename. How to do it safely:
1. `grep -rn "define_genome_effects" .`, `grep -rn "gen_add\|\"epistasis\"\|'epistasis'" .`,
   `grep -rn "generated_additive_tbv\|GE_ADDITIVE_OWNER" .` and
   `grep -rn "trait_names" R/define_effect_cov_matrix.R`. List every hit before editing. At
   0.71.1:
   - The writer name is in 11 `R/` files (30 hits), 4 test files (93) and 7 man pages
     (33), plus `NAMESPACE`, `CLAUDE.md`, both skills, `NEWS.md`, `package_summary.md`
     and 2 `dev/` scripts.
   - Quoted `gen_add` / `epistasis` are in 19 files: 6 `R/`, 4 tests, 3 man pages, 2
     vignettes, 2 skills, `package_summary.md` and `dev/package_summary/package_summary.html`.
     The bare word "epistasis" in prose (25 files) is left alone unless it names the
     `effect_name`.
   - Watch for **substrings**. `define_genome_effects` is not a substring of
     `define_genome_effect_terms`, so a word-bounded replace is safe. Check by eye anything
     ending in `_effects_`.
2. Rename the file with `git mv` (history kept), edit roxygen, then run
   `devtools::document()`.
3. Update every error and message string that names an old function or string.
4. Run the full suite, and gates R1 and R2.

**As built (deviations).**
- The 0.71.2 census found four files the list above missed:
  - `_pkgdown.yml` (the reference index; `pkgdown::check_pkgdown()` fails without it);
  - `README.md`;
  - `vignettes/swine/swine-time-based-age-at-puberty-sex-semen.R`;
  - `dev/benchmarks/benchmark_tgv_scale.R`.
- The test helper view `gen_add_flat` (`helper-genome-effects-db.R`) is renamed `additive_flat`, so R1's grep comes back clean.
- The `trait_var_comp` descriptions in `R/schema.R` say "reserved, with no generator yet" instead of "future".
- **Left alone, as outside this step:**
  - `define_index(trait_names =)`, a different function. It breaks naming rule 1 the same way; this is a possible follow-up.
  - Internal vector arguments named `trait_names`: `load_trait_cov()`, `.gev_read_model()`, `.dae_warn_parent_only()`.
  - §2's audit, a 0.71.0 snapshot that still names the old writer.
- **Review (0.72.1):** realigned continuation lines after the rename, and reworded the `trait_var_comp` schema description. Whitespace and wording only (`_phase_1.md`, "Review pass").
### Step 2 — Part A and the §6C target rules (0.73.0) *(**done** 2026-10-03, see `import_qtl_effect_methods_phase_2.md`)*

- `R/qtl_congruence.R` (congruence, PSD validation, diagnostics), `anchor =`, and
  `.draw_additive_architecture()` (§7.1–§7.5).
- `G` with `effects` or `scale_to_target = FALSE` is refused (§7.1, A19). The arguments
  themselves stay until step 3 (Q21). There is no named-table `effects` work.
- Validation: PSD and dimnames checks in `define_additive_effects()` and
  `define_effect_cov_matrix()`. Atomic target + term writes.
- `"union"`: the "approximate" warning (§7.3). `warn_bounds`, "observed" vs "pool
  expectation" labels (§7.4). The realised-anchor restrictions and size guard (§7.1). The
  line-scoped centring message (§7.6).
- §6C:
  - `trait_var_comp_tbl`, and no-overwrite in both target writers (its error names the
    filter and `remove_rows()` only).
  - `trait_var_comp.line_name` in the base `CREATE TABLE` (`R/open_pop.R`). There is no
    migration, so an `ALTER TABLE` path would have nothing to upgrade. The line resolution
    goes in `load_trait_cov()` / `get_trait_var()` and
    `define_effect_cov_matrix(line_name =)`.
  - Full-precision, transactional target inserts.
  - `target_add_var` / `target_add_mean` (and the `trait_meta` column) removed from
    `define_trait()`. `define_trait_simple()` and `write_trait_var_diag()` deleted.
  - `R/schema.R` `.sm_col()` and `R/sql_utils.R` registries updated (`line_name` added,
    `target_add_mean` dropped). `restore_pop()` refuses a file whose `trait_var_comp` has
    no `line_name`. It also refuses a file that still has `trait_meta.target_add_mean`.
    The first check is the one that always fires: `trait_var_comp` is created by
    `open_pop()`, while `trait_meta` is created lazily (`ensure_trait_tables()`), so a
    pre-0.73 file may have no `trait_meta` at all.
- Reserved `effect_name` constants (§6B rule 5): future names refused with "not yet
  supported". Derived names and reserved names refused as user random effects.
- CLAUDE.md: the Two-Layer paragraph loses `target_add_var` / `target_add_mean`, and
  naming rule 6 (§6B) is added.
- **Test churn to count first:** `target_add_var` has 179 hits across `tests/` and
  `vignettes/`. Every test that re-declares a target block, or re-runs
  `define_additive_effects(G =)` on a stored block, must now pass `G` once and later rely
  on the stored rows. Tests switch to `define_effect_cov_matrix()` or `G =` in the same
  commit. Write fixtures as "target first, then generate". A test that calls
  `define_effect_cov_matrix()` after generating effects breaks again in step 3 (Q21).
  Recount first: the count was 185 at 0.72.4.
- **Small items:**
  - remove `define_trait_simple` from `_pkgdown.yml`;
  - fix the `define_additive_effects()` roxygen example, which stores `G` with
    `define_effect_cov_matrix()` and then passes `G = G` again (refused by no-overwrite);
  - reword error strings that name `define_trait(target_add_var = )`
    (`R/define_additive_effects.R:274`, `R/add_phenotype_stages.R:1165`);
  - scalar `G` for one trait (A15) is **new** code: today's single-trait path ignores `G`
    entirely, and it is read only in the multi-trait path (`:321`).

**As built (deviations), 0.73.0.** Full list in `_phase_2.md`.
- **k = 1 goes through the congruence too.** `rescale_effects_to_target()` is deleted;
  one calibration path serves every k (A1 holds to ~1e-16). As a consequence, a QTL set
  with no segregating locus is now the "anchor cannot carry the target" error instead of
  the old "Falconer V_A is zero" warning that wrote the draw unscaled.
- **`G` and `trait_var_comp_tbl` together** are refused ("not both"). §6C only said the
  same block in both is the "already stored" error; refusing the combination outright is
  simpler and covers it.
- **`anchor = "realised"` also refuses** manual `effects` and `scale_to_target = FALSE`
  (an anchor only means something for a calibration).
- **`define_effect_cov_matrix()` accepts a single number** when one name is given, like
  the generators' `G` (used by the test helper `with_additive_target()`).
- **Reserved names** are refused for fixed effects too (`define_effect_fixed_class()`,
  `define_effect_fixed_cov()`), not only random effects (A20), and in
  `write_phenotype_cov_block()`.
- **Diagnostics (§7.4).** An `ind_haplotype` base is skipped (a copy selection has no
  pairing to compare with); a comparison above `QTL_REALISED_MAX_CELLS`, or individuals
  whose genotypes cannot be collected, is skipped with a `message()`. Diagnostics run
  for the common scope only. **Observed: the default `warn_bounds` fires on most small
  founder pools** (116 of the suite's 127 warnings, pools of 20–100 haplotypes), because
  the pool's sampling LD moves `2 Cov(H)` well away from the genic limit. Q22 decided
  (b) in 0.73.1: the pool comparison is a message.
- **`seed`** must be an integer scalar (validated with the other arguments).
- Paper tests 10 and 11 stay on the internals (dual anchor); test 12 uses
  `add_offspring()` (§1A).

**Corrections, 0.73.2** (Codex review `_phase_2_codex_review.md`, findings 1–3, 5–9;
response at the end of that file). Gates R1–R9 in
`test-define_additive_effects-anchor.R`, plus `test-genome-effects-writer.R` (gate 44,
`parent_origin`) and the union gap test in `test-define_additive_effects.R`.
- Target rank and PSD are decided on the correlation scale (`.qtl_target_std()`); `B0`'s
  columns are normalised to unit anchor variance before the congruence. "Exact" is a
  verified claim: `.qtl_calibrate()` checks `B' M B` against the stored `G` at
  `QTL_CALIBRATION_TOL = 1e-8` (correlation scale) and errors before any write.
- A passed `G` is refused next to a stored non-additive target (same rule as `G = NULL`).
- `"union"`: a positive-variance trait with no QTL is an error; "exact" is computed, so
  overlapping sets with a zero target covariance warn "approximate".
- Pool expectation divisor `n_h` (§7.4). Diagnostics computed before the commit, on the
  un-projected filtered pool.
- `parent_origin` validated before coercion. Anchor-rank check before `set.seed()`.
- The mixed-`parent_origin` message states the real reason (one anchor per call).

### Step 3 — consolidation, P2 and Q18 (0.74.0–0.74.2) *(planned 2026-10-04, see `import_qtl_effect_methods_phase_3_plan.md`)*

**Decided while planning (2026-10-04):**
- `add_tgv()`'s true-index argument is `component_name` (naming rule 1), not `component`
  as §6.3 task 3 said. The DSL keeps `component =` (Q18).
- `ind_true_index` gains `component_name`; row key `(id_ind, index_name, weight_type,
  component_name)`, so an additive and a total index can coexist.
- Three commits, each with the full suite green and a review pause: **3a** (0.74.0)
  consolidation + P2 readers + value names + the active-block prevalence rule; **3b**
  (0.74.1) Q21 removal, owner rule, `define_effect_cov_matrix()` refusal; **3c** (0.74.2)
  Q18. The results go in one `_phase_3.md`, a section per sub-step as it lands.

**As built, 3a (0.74.0).** Gates T3–T9, PH1, PH3–PH6 and the step-3a part of PH7
in `test-tgv-consolidation.R` and `test-phenotype-total-genetic-value.R`; PH5's
group form in `test-group-contributor-determinism.R`.
- **`ind_tgv_total` is an ordered floating sum, not `GEV_ACC_TYPE`** (§6A,
  rewritten): the DECIMAL round trip moved ~10% of one-component totals by one ulp.
  Group-mate sums keep `GEV_ACC_TYPE` (many summands, step 0b).
- **`add_tgv()` re-evaluation upserts** surviving components and deletes only the
  ones the model no longer produces. The old delete-then-insert would have wiped
  custom columns on every `add_phenotype()` call, now that `ind_tgv` is the table
  phenotypes materialise (the old breeding-value table upserted).
- **Bug found:** `define_phenotype()` validated `components` after writing
  `phenotype_meta` (and after deleting the old rows under `overwrite = TRUE`). All
  component checks now run first; PH3 pins it.
- `.gev_component()` asserts the closed set `TGV_COMPONENT_NAMES`.
- The readers `.tbv_by_id()` / `.group_mate_tbv()` became `.tgv_by_id()` /
  `.group_mate_tgv()` here (they changed body anyway); the rest of the Q18 internal
  renames stay in 3c. One shared reader, `.tgv_read()` (`R/add_tgv.R`), serves
  phenotypes and true indices.
- Test files renamed: `test-add_tbv.R` → `test-add_tgv_breeding_value.R`,
  `test-add_tbv_index.R` → `test-add_tgv_index.R`. The `.gev_warn_tbv_stale()`
  tests were replaced by gate T5 (a custom additive term is in `additive`).
- PH8's grep must except the refusal code and tests that name the removed table
  (`restore_pop()`'s pre-0.74 check, T8, `test-open_pop.R`'s absence check).


- Consolidation, as specified in §6.3 (tasks 1–7).
- P2 (§6A):
  - Every phenotype path reads `ind_tgv_total` (or a listed component).
  - The `component_names` DDL default becomes `'total'`.
  - `add_phenotype()` calls `add_tgv()`.
  - `ind_tgv_total` sums exactly.
  - The simple-phenotype precheck and the prevalence threshold change, including Q21's
    owner rule (§6A).
- Q21 (a):
  - `define_additive_effects()` loses `effects` and `scale_to_target`.
  - Their call sites move to `define_genome_effect_terms()` + `ad_terms()`. Estimates
    range from 76 to 93 in `tests/`, plus vignettes and examples, so count them first.
    To keep the same genetic values, use `ad_terms(coding = "cockerham", p = <base p>)`,
    because the generator centres at `p_base`. Functional coding shifts TBVs by a
    constant.
  - `define_effect_cov_matrix()` refuses a genetic block whose traits already have
    `generated` terms of that kind.
  - The `define_phenotype(prevalence = )` roxygen's 0.72.3 caveat is replaced by the
    owner rule.
  - `define_phenotype()` roxygen states that `mean` is an intercept.
  - Q23 is closed (0.73.2): `allow_reserved_owner` is gone from the exported writer.
    PH7 keeps a regression that the exported writer cannot write `generated`, not only
    ordinary generator calls.
  - A passed `G` next to a stored non-additive block stays refused, with the
    store-then-`trait_var_comp_tbl` route (decided 2026-10-04: no new argument).
    `define_genome_effects()` follows the same rule in step 5.
- The `ind_tgv` half of §6B: `order1_*` → `additive` / `dominance` / `indicator`.
- The §6.2 breeding-value roxygen sentence on `add_tgv()`, `define_genome_effect_terms()` and `ad_terms()`.
- Q18: `formula_tbv` → `formula_tgv` (argument, `phenotype_meta` column, internals); the
  DSL's named-only `component =` / `table =`, validated.
- `restore_pop()` refuses `ind_tbv` / `phenotype_meta.formula_tbv`.
- CLAUDE.md hard rules (D7, "one evaluator"), design principle 4's action list, and the
  skills drop `add_tbv()` / `ind_tbv`.

### Step 4 — Part B (0.75.0)

- Report the reference population and interpretation of every estimate in the output
  (Codex review, item 3): `anchor = "genic"` on a selected cohort is a projection at
  `base_tbl` frequencies, `"realised"` measures that cohort, and `base_tbl = NULL` drifts
  with the cohort (not a check of the generation target). Label case 2 "evaluated
  additive variance". Add a worked target-vs-measured example filtering
  `trait_var_comp` by `line_name`.

- `extract_genetic_variance()` (§8), with `between_components` (Q16).
- The NOIA conversion pair `.noia_to_stored()` / `.stored_to_functional()` (Q13), which
  case 1 needs.
- **`aa_terms()`**, and the fixed column set for `ad_terms()` / `genotype_terms()` (§9.3,
  gate C16). They move here from Part C because B2 and B10 build A + D + A×A fixtures by
  hand through `define_genome_effect_terms()`.

### Step 5 — Part C (0.76.0)

- Generator `define_genome_effects()` (the name freed in step 1) in a new
  `R/define_genome_effects.R`. Calibration internals go in `R/genome_effects_calibration.R`.
- The §5 replacement rules. The `.ge_resolve_deletes()` message names both generators.
- §6C's non-additive refusal gains `define_genome_effects()` as a third fix.
- Gate C19.
- One short vignette with the four paths of the Codex review's "recommended first-release
  contract" (writer, additive generator, genome generator, extractor).
- Apply 0.73.2's target rules to `G_A`, `G_D`, `G_AA` and the additive floor
  (`.qtl_target_std()` + `.qtl_calibrate()` verification), and test non-additive values
  in **phenotypes**, not only `ind_tgv`.
- The vignette states the scope promise (Codex review, "Changes to the remaining plan"
  item 5 and its table "What 'target this G in that population' means"): exact for a
  feasible `G` under the named anchor and QTL set; a line call calibrates its own
  variant only; `"union"` does not hit non-zero off-diagonals; one call has one
  `parent_origin` scope; the extractor measures the population the user means.

---

## 11. Acceptance gates

**Rename release (step 1)**

- R1. After step 1, `grep -rnw "define_genome_effects"` over `R/`, `tests/`, `man/`, `vignettes/`, `dev/`, `NAMESPACE`, `CLAUDE.md`, `package_summary.md`, `README.md`, `_pkgdown.yml`, `tools/`, `plans/` (except the rename records: this plan, `_phase_1.md`, the Codex review) and `.claude/skills/` returns nothing. The same holds for `gen_add` (including the test view `gen_add_flat`, renamed `additive_flat`), for `epistasis` used as an `effect_name`, for `generated_additive_tbv` / `GE_ADDITIVE_OWNER`, and for `define_effect_cov_matrix(trait_names =`. `NEWS.md` (past entries) and closed plans are excluded. `define_genome_effect_terms` is exported, documented, and listed in `NAMESPACE`. In step 5 this gate's `define_genome_effects` check is replaced by "no message or roxygen written in steps 1–4 names `define_genome_effects()`" (checked by `grep` before step 5 starts), since the name then returns as the generator.
- R2. The full suite passes with only renamed calls and strings changed in `tests/`: no assertion, tolerance or fixture changes. A `git diff --stat` of `tests/` shows renames only.

**Part A**

- A0. Determinism: the same `set.seed()` twice gives `expect_identical()` `genome_effects` / `genome_effect_members` rows, for both anchors and $k\in\{1,2\}$. Not `expect_equal()` (CLAUDE.md).
- A1. $k=1$, `anchor = "genic"`: the congruence output equals the scalar rescale $\mathbf b_0\sqrt{G/\sum_j w_jb_{0j}^2}$ computed **in the test** from the same draw (1e-12). No comparison with v0.71.0 output: CLAUDE.md forbids golden-from-old tests.
- A2. $k=2$, `"genic"`: $\mathbf B^\top\mathrm{diag}(n_{\text{eligible}}pq)\mathbf B=\mathbf G$ to 1e-10, including off-diagonals. This is an algebraic identity, so a few cases suffice: three seeds × $\{m=k,\ m\gg k\}$ × {full-rank, rank-1 $\mathbf G$}, each also checked against an independent dense-matrix oracle computed in the test ($\mathbf A$ from `eigen()` of $\mathbf C$ and $\mathbf G$ directly). No 50-seed sweep.
- A3. `"realised"`: $\mathrm{Cov}(\mathbf X\mathbf B)=\mathbf G$ on the base individuals to 1e-10.
- A4. Rank-deficient architecture (two proportional columns) with a rank-1 target in a different direction: exact. This is the Codex F2 / 2026-09-19 regression.
- A5. Two distinct rank errors, each naming the ranks it compared. (a) $\mathrm{rank}(\mathbf G)>\mathrm{rank}(\mathbf M)$: the anchor cannot carry the target, for any effects. (b) $\mathrm{rank}(\mathbf M)\ge\mathrm{rank}(\mathbf G)$ but $\mathrm{rank}(\mathbf B_0^\top\mathbf M\mathbf B_0)<\mathrm{rank}(\mathbf G)$ (e.g. two proportional architecture columns and a full-rank target): the drawn architecture cannot, and the message says so rather than blaming the anchor.
- A6. `"realised"` with a `founder_haplotypes` or `ind_haplotype` `base_tbl`, or with `line_name` / `parent_origin`, errors with the reason.
- A7. `method = "union"` with $k\ge2$ and a non-diagonal target: per-trait QTL sets are preserved (no locus gains a trait), and a **warning** containing "approximate" and the delivered correlation is raised. The same target under `"shared"` is exact and raises no such warning. The roxygen `@details` states the exactness conditions of §7.1.
- A8. `warn_bounds` fires when the comparison population is far from the genic limit (an inbred panel) and not on HWE/LE data. `warn_bounds = NULL` silences it. The warning text says "observed" or "pool expectation" as appropriate.
- A9. Parent-origin scoping: $n_{\text{eligible}}=1$ enters the weights, i.e. $\sum_j w_j b_j^2 = V_A$ with $w_j = p_jq_j$ for a `parent_origin = 1` trait (computed in the test).
- A10. Atomicity: a passed `G` for traits with no stored block, and trait 2 invalid after trait 1 → `trait_var_comp` still has no block for them, and all three `genome_effect*` tables are unchanged. The same holds for an infeasible rank. (A passed `G` *over* a stored block never gets this far; that is A15.)
- A11. Target validation: an indefinite `G` (passed or stored) errors, naming the most negative eigenvalue, with no write. A zero-rank target and a rank-deficient but feasible target are separate cases with separate expected outcomes. `define_effect_cov_matrix()` refuses an indefinite matrix, and refuses `additive_by_dominance` / `dominance_by_dominance` with "not yet supported", writing nothing to either var-comp table. A `G` whose dimnames disagree with `trait_name` (or with `trait_name` in `define_effect_cov_matrix()`) errors instead of being relabelled.
- A12. *(Withdrawn 2026-10-03, Q21 (a): manual `effects` is removed instead of becoming a named table.)*
- A13. A line-scoped `define_additive_effects()` call emits the §7.6 centring message. A common-scope call does not.
- A14. `anchor = "realised"` after `restore_pop()`: the same `base_tbl` filter gives the same row and locus order, and `expect_identical()` stored terms for the same seed.
- A15. No overwrite (§6C): passing `G` when an `additive` block exists for any of the traits (same `line_name`) errors, **including an identical matrix**, and the error contains a working `remove_rows()` call. After that call, the same `G` succeeds. A filtered `trait_var_comp_tbl` that hides the stored block does not let the write through. `define_effect_cov_matrix()` refuses an existing block the same way. A scalar `G` with one trait is stored as a 1×1 block.
- A16. `trait_var_comp_tbl` (§6C): a table other than `trait_var_comp` errors. A partial block (one triangle filtered away) errors. Two candidate sets for one `effect_name` error and ask for a filter. A stored `dominance` target with the default `NULL` errors in `define_additive_effects()`, naming the filter and the `remove_rows()` call (and, from step 5, `define_genome_effects()`), and the same call with `filter(effect_name == "additive")` succeeds. Partial trait set: with a stored `additive` block for `ADG` and `BF`, `trait_name = "ADG"` with the default `NULL` errors, naming both traits and the two fixes, before any write and without consuming RNG. The same call with `trait_var_comp_tbl` filtered to the `(ADG, ADG)` row succeeds, and with `trait_name = c("ADG", "BF")` it succeeds.
- A17. Lines (§6C): a `line_name = "C"` call reads line C's block when it exists and falls back to the `NULL` block otherwise. A passed `G` is written with the call's `line_name`. A line-C block and a population-wide block for the same traits coexist. Fallback is per `effect_name`: line C's own `additive` block plus the shared `dominance` block resolve together. A partial line-C block errors. `load_trait_cov()` / `get_trait_var()` never mix lines, and `add_ebv()`'s parameter file reads the `NULL` rows. `define_effect_cov_matrix(line_name = "C")` writes line rows. `define_trait()` no longer accepts `target_add_var` / `target_add_mean`, `trait_meta` has no `target_add_mean` column, `define_trait_simple()` is not exported, and `restore_pop()` refuses a file whose `trait_var_comp` lacks `line_name`, or that still has `trait_meta.target_add_mean`.
- A18. Precision (§6C): `define_effect_cov_matrix()` with `cov_matrix = matrix(c(1/3, 1/7, 1/7, 2/3), 2)` stores values that read back `expect_identical()` to `G`, and a generator calibrated to it hits it to 1e-12. A failure halfway through a target write leaves `trait_var_comp` unchanged.
- A19. *(Step 2 only; moot from step 3, when both arguments are removed.)* `G` with manual `effects`, and `G` with `scale_to_target = FALSE`, each error before any write. The message says that `G` is a calibration target for sampled effects only, and that manual or unscaled effects take no target. It must **not** suggest storing one with `define_effect_cov_matrix()`: a target next to uncalibrated terms is the Q21 hole.
- A20. Reserved names (§6B rule 5): `define_effect_cov_matrix()` refuses `total`, `unpartitioned` and `between_components` as `effect_name`. A user random effect named `additive` (or any reserved name) errors in the phenotype layer. Nothing is written to either var-comp table.
- A21. Part A size guard: `anchor = "realised"` above the cell limit errors with $n$, $m$ and the limit, before any write.
- A22. RNG on refusal (CLAUDE.md): every validation error in `define_additive_effects()` (bad target, overwrite, `G` + `effects`; owner content from step 5) leaves `.Random.seed` unchanged, **including with `seed =`** (§7.1, "`seed` after validation"). No draw happens before validation passes.

**Part B**

- B1. `anchor = "genic"`, fully covered A + D + A×A model: the block covariance matrices sum to the genic `total` to 1e-10, and `additive` equals $\sum 2pq\alpha^2$ computed in the test.
- B2. On an A + D + A×A model written by hand through `define_genome_effect_terms()` from the source project's calibrated $\mathbf B_\alpha, \mathbf B_d, \mathbf B_{aa}$ (Cockerham terms via `ad_terms()` / `aa_terms()`), evaluated on the base individuals, realised blocks **and** cross-block covariances equal the source project's `nonadd_decompose()` (`real_*`, `cov_A_D`, `cov_A_AA`) to 1e-10. The source returns no D–A×A covariance, so the test computes $\mathrm{Cov}(D, AA)$ itself from the source's `DD` and `AA` value matrices, both orientations for the off-diagonal trait pair, and the port matches it to 1e-10.
- B3. Accounting under LD / non-HWE: on a structured cohort (an inbred or repulsion-LD panel) with every term covered, block variances do **not** sum to `total`. Blocks plus the `between_components` row do, to 1e-10, for every trait pair including the off-diagonal ones ($k=2$), and `between_components` equals the sum of **every** ordered cross-block covariance computed in the test ($A$–$D$, $A$–$AA$ and $D$–$AA$, both orientations), not only the two the source returns. There is no `unpartitioned` row. Under `"genic"` on the same model no `between_components` row is returned.
- B4. Dispatch: a model of two custom owners on the same loci (additive + dominance, common scope) is case 1, not two separate models. A line-scoped F1 additive model is case 2, reported as "evaluated additive variance" (`decomposition = "additive_only"`), and is not re-projected as one unscoped coefficient matrix.
- B5. A hand-written indicator surface gives a correct `total` (equal to `var()` of `ind_tgv_total` computed in the test), an `unpartitioned` row holding that surface's own variance only, and `decomposition = "partial"`, and no error.
- B6. A line-scoped additive crossbreeding model (common + line-A + line-B variants) on F1s: `additive` equals the variance of the evaluated `additive` exactly, with no `unpartitioned` row.
- B7. Read-only: `ind_tgv` and every other persistent table are unchanged after the call.
- B8. Output joins to `trait_var_comp` on `(effect_name, trait_name_1, trait_name_2)` with no renaming. A model with no `dominance` terms returns no `dominance` row, so a stored `dominance` target that the generator was told to leave out (§6C) is absent from the `inner_join` and present in the `anti_join()`.
- B9. Determinism: two calls, and a call after `restore_pop()`, give `expect_identical()` output. The realised path above the size limit errors with $n$, $m$ and the limit.
- B10. Coding and anchor scope: the same A + D + A×A model written once in functional coding (§9.3 example) and once in Cockerham coding gives the same report under both anchors (1e-10). `anchor = "genic"` on a case 2 or case 3 model errors and names `"realised"`. A non-`NULL` `base_tbl` with `anchor = "realised"` errors.
- B11. `decomposition` is `"full"` for a common-scope A + D + A×A model (case 1), `"additive_only"` for a scoped additive model (case 2), and `"partial"` for one with an indicator surface (case 3).
- B12. Cohort (§8): a `tbl` that selects an individual with no value for one of the call's traits errors, giving the count. A `tbl` selecting one individual errors. Every output row carries the same `n_ind`. A locus that is monomorphic in the cohort gives finite output, equal to the same model without that locus.

**Part C**

- C1. Genic: stored terms give the genic covariances of the blocks present ($\mathbf G_A$, and $\mathbf G_D$ / $\mathbf G_{AA}$ when present) to 1e-10 (from the stored $\alpha$, $d$, $e$ and $p$ alone).
- C2. For **every** diploid genotype combination at a small model (3 loci, one pair: all $3^3$ genotypes), functional $g$ minus `add_tgv()`'s total is the same constant, equal to $\mu$ of §3, to 1e-12, under both anchors. Test the per-genotype difference, not centred totals.
- C3. `additive` equals $\mathbf Z_A\mathbf B_\alpha$ under `"genic"`, to 1e-10.
- C4. Additive-only reduction, two levels. (a) **Internals**: with no dominance or A×A block and a **supplied** architecture $\mathbf B_0$, the calibrator's output equals `.qtl_congruence()`'s on the same $\mathbf B_0$ and $\mathbf M$, to 1e-12 (the source method's reduction test). (b) **Public** (decided 2026-10-01, §0A): with the same `set.seed()`, traits, target, `anchor` and `base_tbl`, `define_genome_effects()` with an additive-only target writes `expect_identical()` term, member and origin rows to `define_additive_effects()` with its defaults (`distribution = "normal"`, `method = "shared"`, `seed = NULL`), `effect_owner = "generated"` included (§5). Row ids are compared after dropping the `id_*` key columns, which `next_int_id()` assigns per database. This holds because both draw through `.draw_additive_architecture()` first (§7.2, §9.2), so it also constrains Q5: a future sampler must be added to that helper, not to one generator. Checked for $k=1$ and $k=2$ and both anchors.
- C5. $\mathbf G_A$ below the floor errors, naming the floor.
- C6. $k=1$ matches `zeng_appendix_A()` to 1e-12.
- C7. The inbreeding-depression target is exact for $k=1$ and for diagonal $\mathbf G_D$ with $k\ge2$, and reported (not enforced) for non-diagonal $\mathbf G_D$.
- C8. One owner (§5): (a) `define_genome_effects()` on a trait with common, line-A and line-B `generated` additive variants replaces all of them, and its message reports the replaced count and the line-scoped count. (b) On a trait whose `generated` model has a `dominance` or interaction term, `define_additive_effects(line_name = "A", G = …)` (no line-A block is stored, so the A15 overwrite check does not fire) errors before any write, on the §5 content check. The error names the additive-only `define_genome_effects()` re-run, and `trait_var_comp` has no line-A block afterwards. After that re-run, `define_additive_effects()` **without `G`**, with `trait_var_comp_tbl` filtered to `effect_name == "additive"`, succeeds. Without the filter it still errors on the stored `dominance` target (§6C). (c) Custom-owner terms on the same trait are untouched by every call. (d) The writer refuses `effect_owner = "generated"` without `allow_reserved_owner = TRUE`.
- C9. A×A terms have two `additive` members and land in `interaction`.
- C10. Non-diploid / non-autosomal QTL are refused with the reason.
- C11. Multi-trait write is one transaction: an error on trait 2 (invalid target, infeasible floor) leaves trait 1's prior model **and** `trait_var_comp` intact.
- C12. Inbreeding depression, two parts. (a) **Exact, deterministic**: build a genotype fixture whose genotype counts equal $n\,(p^2+Fpq,\ 2pq(1-F),\ q^2+Fpq)$ exactly at each locus (loci independent, $n$ chosen so counts are integers). The mean of the stored `dominance` value over it is $-F\sum_j 2p_jq_jd_j$ to 1e-12. (b) **Individual level**: for a handful of offspring, `ind_tgv` equals the value computed in the test from their genotypes and the stored terms. A stochastic full-sib cohort check is optional. If kept, it needs enough replicates and a tolerance derived from the sampling variance, and it must not rely on pedigree $F$ holding in a small sample.
- C13. (Needs precondition P2) A single-trait `define_genome_effects()` model with $G_D$ large relative to a feasible $G_A$ and no A×A: `add_phenotype()` phenotypic variance minus residual variance tracks $V_A+V_D$, not $V_A$. This is the end-to-end check that the effects reach the observation layer.
- C14. Determinism: the same `set.seed()` twice gives `expect_identical()` stored terms, including random pairing.
- C15. Pairs (Q8): (a) supplied `pairs`: unknown loci, loci outside the filter, self-pairs and reversed duplicates each error naming the keys. A hub design (one locus in three pairs) is accepted, and its genic $\mathbf G_{AA}$ and $\mathbf G_A$ are exact (C1). (b) Random matching: no locus appears in two pairs; the pairs are written in sorted order; the same seed gives `expect_identical()` pairs. (c) `n_pairs = NULL` gives $\lfloor m/2\rfloor$ pairs (odd $m$: one locus unpaired) with the message. $n_{\text{pairs}} > \lfloor m/2\rfloor$ errors naming the maximum and `pairs`. (d) `pairs` and `n_pairs` together, or either without an A×A block, error.
- C16. `aa_terms()` output is identical under reordering of its input pairs and under swapping the two loci of a pair, and never rescales `e`. `aa_terms()`, `ad_terms()` (both codings) and `genotype_terms()` outputs `rbind()` with each other. Missing `p_1` / `p_2` errors under both codings. Under functional coding `aa_terms()` reports each pair's $\mu$ share, $e(2p_1-1)(2p_2-1)$, and writes it nowhere. (Ships in step 4.)
- C17. Target semantics (§6C): a block absent from the resolved rows is absent from the model; a stored block is used; a zero matrix is calibrated as an exact zero and stored as zeros; a passed matrix is written. `trait_var_comp_tbl |> filter(effect_name != "additive_by_additive")` gives an A + D model with the stored A×A row untouched. Passing `G_D` while a `dominance` block is stored errors. No `additive` block errors. An additive-only target is accepted (C4). Adding a `dominance` target later and re-running replaces the trait's `generated` model, and the additive architecture draw is unchanged by the new block (same seed). `inbreeding_depression` with no dominance block errors and writes nothing.
- C18. Dimnames: a `G_A` / `G_D` / `G_AA` whose dimnames differ from `trait_name` (other names, or the same names in another order) errors and writes nothing. Unnamed matrices are taken in `trait_name` order.
- C19. (Moved here from the consolidation gates: it needs steps 4 and 5.) For a model **generated with `anchor = "genic"`**, `extract_genetic_variance(anchor = "genic", base_tbl = <the generation base>)`'s blocks equal the targets to 1e-10. On a genotype fixture with **exact** HWE genotype counts at frequencies equal to the stored centres, and LE at every paired locus pair, the covariance of `ind_tgv` `additive` also equals the `"realised"` `additive` block to 1e-10. Away from those conditions they differ, because the stored contrast fixes $b=q-p$ and the realised projection uses the observed $b$ (§4.3.1). The test shows one such difference.
- C20. Generator, realised anchor: for a model generated with `anchor = "realised"`, `extract_genetic_variance(anchor = "realised")` on the same base individuals returns $\mathbf G_A$, $\mathbf G_D$, $\mathbf G_{AA}$ to 1e-10. The §7.1 restrictions (individuals only, size guard) error as in A6 and A21. `warn_bounds` fires on an inbred panel, and `NULL` silences it. A validation error leaves `.Random.seed` unchanged (as A22).

**Consolidation (precondition P1, step 3, §6)**

- T1. `ind_tbv`, `add_tbv()` and `tbv_value` appear nowhere outside history (part of PH8's `grep`).
- T2. The retargeted first-principles oracle (§6.3 task 6) agrees with `add_tgv()`'s `additive` row for a purely additive model.
- T3. For an additive-only model, `additive` is the only component row written, and it equals `ind_tgv_total`.
- T4. For a hand-computed mixed model (additive + dominance + an A×A term + an indicator term), the `additive`, `dominance`, `indicator` and `interaction` rows each equal the test's own computation, and they sum to `ind_tgv_total`.
- T5. A custom-owner additive term written with `define_genome_effect_terms()` appears in the `additive` row: the discrepancy §6.1 exists to close.
- T6. On a model with dominance, `add_index()` on an unfiltered `ind_tgv` errors (more than one value per individual × trait) instead of summing. Filtered to `additive` it succeeds.
- T7. `ind_true_index` equals the `ind_tgv` `additive` values × `index_meta` weights computed in the test, and with `component = "total"` it equals the total × weights. `add_tgv()` accepts `index_names`, `weight_type` and custom fields via `...`.
- T8. `archive_replicate()` stamps `replicate` on the archive copy of `ind_tgv`. `remove_rows()` works on `ind_tgv`. `restore_pop()` refuses a file that still has `ind_tbv`.
- T9. `schema()` and `describe_table()` render `ind_tgv` with every registry entry present, and no `ind_tbv` entry remains.

**Phenotype integration (precondition P2, consolidation release, step 3)**

- PH1. `self` for a simple phenotype, and dam / sire / group contributors for a composite, read `ind_tgv_total` by default. On an additive-only model the output is `expect_identical()` to reading `additive`. On a model with dominance it is not.
- PH2. The formula DSL (`formula_tgv`, Q18) reads the total by default: on a model with dominance, `"WWD + dam(WWM)"` equals the calf's `ind_tgv_total` for `WWD` plus the dam's for `WWM`. `dam(WWM, component = "additive")` reads the dam's `additive` row only. `component = "bogus"` errors in `define_phenotype()`. `group_sum(trait, col, table = other_table)` reads the named table. An unknown table or column, a non-identifier name, and a positional third argument each error in `define_phenotype()` before anything is written. `formula_tbv` is gone from the API and from `phenotype_meta` (grep, as R1).
- PH3. `component_names`: `"total"`, a listed component, a component the model has no terms for (contributes 0), and a name outside the vocabulary (errors at `define_phenotype()`) each behave as §6A says. A failing `add_phenotype()` leaves `ind_phenotype` and `phenotype_random_effects` unchanged (D7).
- PH4. Mean: on a non-HWE / LD base with dominance, the realised base phenotypic mean minus `phenotype_meta.mean` equals the base mean of `ind_tgv_total` (plus the sample mean of residuals), showing that `mean` is an intercept.
- PH5. Exact total (§6A): `ind_tgv_total` for a model with all four components is `expect_identical()` under `SET threads = 1` and `SET threads = 8`, and so is the phenotype built from it. The same holds for a composite phenotype with a `group_sum()` and a `group_mean()` contributor, read both as `"total"` and as one listed component. Use groups large enough for DuckDB to parallelise the aggregate (pens of 200, as in `test-group-contributor-determinism.R`): at pens of 10 the old floating `SUM()` never diverged, so a small fixture passes on broken code.
- PH6. Precheck (§6A): a simple phenotype whose trait has **only** custom-owner terms (written with `define_genome_effect_terms()`) records phenotypes, with no "call define_additive_effects() first" error. A trait with no terms still errors.
- PH7. Prevalence (§6A): with stored `additive` and `dominance` targets and a model with both kinds of terms, the binary threshold uses their sum. With `additive` only it is unchanged from the additive-only formula. A stored A + D + A×A target with the model generated from a filter that left A×A out: the threshold uses A + D only. A model with terms but no stored target (custom terms only) errors, naming `thresholds =`, as does a composite phenotype with `prevalence`. A stored target plus **any** non-`generated` active term errors too (Q21). Fixture: target 1, `ad_terms()` effects of 10 at ten loci (genic variance about 550), and `prevalence = 0.1`. This is Codex finding 1's reproduction, which passes silently at 0.72.3. After generating effects, `define_effect_cov_matrix()` on that trait's `additive` block errors, including after `remove_rows()` of the old target. Explicit `thresholds` works in every one of these cases. Each error leaves `ind_phenotype` and `phenotype_random_effects` unchanged (D7).
- PH8. Step-3 grep gate: `ind_tbv`, `add_tbv`, `tbv_value`, `formula_tbv` and `order1_` appear nowhere in `R/`, `tests/`, `man/`, `vignettes/`, `dev/`, `NAMESPACE`, `CLAUDE.md` or the skills (past `NEWS.md` entries and closed plans excepted). `test-add_phenotype_failure_contract.R` asserts the D7 contract against `add_tgv()`.

---

## 12. Open questions

### Q1 — Function names *(decided 2026-10-01 → §0A)*

The Part C generator is **`define_genome_effects()`**. The current low-level writer of that
name becomes **`define_genome_effect_terms()`** in a rename-only first release. Rejected:
`define_nonadditive_effects()` (the draft name: it suggests the function writes *only*
non-additive effects, but it writes the additive block too and accepts additive-only
targets); `define_genetic_effects()` (close to the chosen name, but drops the link to the
`genome_effect*` tables); `define_dominance_effects()` (too narrow once A×A is in).
Rejected for the writer: `define_effect_terms()` (loses the table link) and
`define_custom_effects()` (the writer is not limited to the `custom` owner).

### Q2 — Should `define_additive_effects()` become a thin wrapper of Part C? *(decided 2026-10-01: no, keep both)*

**Decision:** two generators, each keeping its own arguments. `define_genome_effects()`
accepts any set of blocks that includes `additive`, so it covers additive-only models too,
with the same coefficients as `define_additive_effects()` (C4 (b)). `define_additive_effects()`
keeps the options it alone has: `line_name`, `parent_origin` and the crossbred per-copy
fallback, `method = "union"`, `distribution`. Neither becomes a wrapper of
the other. Both share `.qtl_congruence()` and `.draw_additive_architecture()`, so there is
one implementation of the shared algebra.

A single merged generator was considered (2026-10-01): since §6C, the rows of
`trait_var_comp` already decide which blocks a model has. It was not taken because the two
argument sets barely overlap. A merged signature would have about 16 arguments, roughly
half conditional on what is stored. If Part C later gains line scoping (Q14), most of
those arguments would apply to both, and a merge can be reconsidered then.

### Q3 — One reserved owner or two? *(decided 2026-10-01: one, `generated` → §5)*

**Decision:** one reserved owner, `generated`, shared by both generators. The string and
constant are renamed in step 1. The replacement rules (§5) land with Part C.

Why. Consolidation removes the owner's third job (`add_tbv()` defining the breeding value
from one owner's terms). That leaves protection and the replacement boundary, which one
owner expresses exactly. The alternative, two owners (`generated_additive`,
`generated_genome`) with mutual exclusion per trait, was rejected:
- Switching generators for a trait would need a "clear the other owner" call, and **none
  exists**: the writer rejects empty `terms` (`R/define_genome_effect_terms.R:246`), and row
  deletion from the `genome_effect*` tables is refused by design.
- An additive-only model would differ by owner depending on which generator wrote it.
- The owner would record provenance (which function), which nothing reads.

The cost of one owner is that `define_additive_effects()` needs a content check (refuse
when non-additive `generated` terms exist), and `define_genome_effects()` deletes
line-scoped variants on purpose. Both are stated in §5 and pinned by C8.

### Q4 — `method = "union"` under the congruence *(decided: the recommended option → §7.3, gate A7)*

| Option | Notes |
|---|---|
| **Keep per-trait scalar for `"union"`, report delivered $r_g$** ← recommended | No surprise loss of trait-specific QTL sets. Exactness is available via `"shared"` |
| Congruence anyway, warn that supports merge | Exact $\mathbf G$, but every union locus affects every trait, which is a different architecture from the one requested |
| Block-structured solve (shared loci carry the covariance, private loci only variance) | Preserves supports. Feasible only if the shared block can carry $\mathbf G$'s off-diagonal. New mathematics, not in the manuscript |

### Q5 — Architecture samplers (`arch_*()`) and `marginal =` *(deferred to a later change; constrained by C4 (b): new samplers go into `.draw_additive_architecture()`)*

The source has seven samplers (`gaussian`, `mvt`, `gamma`, `laplace`, `sparse`, `major`,
`maf`) and a rank-restoring iteration for non-elliptical marginals at $k>1$.
Recommended: **a separate, later change.** Replacing `distribution` with `architecture` is
an argument-surface change on its own. Part A delivers exactness with today's draws.

### Q6 — Expose `marginal_observed` / `dual`? *(decided: no → §7.1)*

Recommended: no. `dual` (hit founders and the limit at once) is feasible only under
Theorem 3's condition and gives up the choice of architecture. `marginal_observed` is a
diagnostic. Port the internals so the diagnostics can use them.

### Q7 — Diagnostics storage and the `effect_name` vocabulary *(both decided)*

**Storage: none** (§7.4). The draft proposed `gen_add_founder`, `gen_add_equilibrium`,
`dominance_delivered` and `gen_add_floor` rows in `trait_var_comp`. That would mix
derived measurements into the table generators read their targets from
(`trait_var_comp_tbl = NULL` reads the stored rows, §6C). Each can be measured again on any population that can
still be named (principle 6). §7.4 says what cannot be rebuilt.

**Vocabulary: §6B.** The draft options (`A`/`D`/`AxA`, `gen_add`/`gen_dom`/`gen_add_add`,
full words) were weighed on 2026-09-23; full words won.

### Q8 — Pair selection *(decided 2026-10-02 → §9.1)*

**Decision.** Random pairs (`pairs = NULL`) are a **random matching**: each locus in at most
one pair, AlphaSimR's structure. `n_pairs = NULL` defaults to $\lfloor m/2\rfloor$ (every
QTL once), announced by a message that also says how to change it. **Supplied `pairs` may
reuse a locus** across pairs (hub designs). The math supports both: the method accumulates
a locus's pairs, and pairs sharing a locus stay uncorrelated under HWE + LE, so the genic
anchor is still diagonal.

Rejected: free sampling over all $m(m-1)/2$ pairs (option (b)). It gives random accidental
hub loci and departs from AlphaSimR. "All pairs among a subset" (MoBPS-like) is not
offered; a supplied `pairs` table covers it. If someone asks for random designs where loci
repeat, a `pairing =` argument can add option (b) later without changing the default.

### Q9 — Report functional $(a, d, e)$ anywhere persistent?

Users from AlphaSimR think in functional effects, and they are recoverable from the stored
statistical ones plus $p$. Recommended: a helper
`extract_functional_effects(pop, trait_name)`, not storage. Low priority.

### Q10 — Default of `phenotype_components.component_names` under P2 *(decided: `"total"` → §6A, PH3)*

§6A proposes `"total"`. Alternative: keep `"additive"` and make users opt into
non-additive expression per phenotype. Recommended: **`"total"`**. A phenotype that
silently ignores calibrated dominance is the worse surprise. The additive-only case is
unchanged either way.

### Q11 — A pool-based `"realised"` anchor *(decided: diagnostics only for now → §7.4, A6)*

$\mathbf M_{\text{pool}} = 2\,\mathrm{Cov}(\mathbf H)$ over `founder_haplotypes` is the
expected founder genotype covariance under random pairing, LD included, and it is
available **before** `add_founders()`. §7.4 already uses it for diagnostics. It could
also serve as `anchor = "realised"` when `base_tbl` is the pool, instead of the A6
error. Recommended: **diagnostics only in Part A**, and revisit if users ask for an
LD-aware anchor without individuals. It is the same congruence with a different
$\mathbf M$, so it is cheap to add later.

### Q12 — `seed` arguments *(decided for the new function: none → §9.1; the package-wide question stays separate)*

`define_additive_effects(seed =)` calls `set.seed()` (`R/define_additive_effects.R:216`),
which silently resets the caller's RNG stream partway through a script.
(`add_offspring(seed =)` is different and fine: it seeds its own `dqrng` sub-stream and
leaves base-R state alone.) Recommended: **no `seed` on the new function**, and a separate
package-wide pre-1.0 decision on removing it from the existing ones. Do not mix it into
this plan.

### Q13 — Where the NOIA → storage conversion lives *(decided: the recommendation, in step 4 → §10)*

The conversion $(a,d,e) \to (\alpha, d, e)$ at $p$ (§3) and its inverse (§8) will be
needed by the generator, by `extract_genetic_variance()`, by `aa_terms()` and by any
`extract_functional_effects()` (Q9). Recommended: **one internal pair**
(`.noia_to_stored()` / `.stored_to_functional()`) in `R/genome_effect_terms_builders.R`
next to `ad_terms()`. Today `ad_terms()`'s Cockerham branch takes $\alpha$ as given and does
no $(a, d) \to \alpha$ conversion; only its functional branch computes $\mu$. So the
conversion is new code. `ad_terms()` calls it rather than duplicating it.

### Q14 — Line-specific non-additive effects *(deferred: a follow-up release after Part C)*

The first release of Part C writes common-scope terms only. Heterosis still appears,
because F1 dominance follows from allele-frequency differences between lines.
Line-specific *additive* effects already exist (`define_additive_effects(line_name =)`).
A line-scoped Part C is a moderate extension, with no DDL:
- call once per line with `line_name = L` and `base_tbl` defaulting to L's founder pool,
  so the targets hold **within** line L;
- store the `additive` members with an `exact L` origin, the `dominance` members with
  the multiset {`exact L`, `exact L`}, and the A×A members with `exact L` on each
  member;
- crossbred copies fall back to a **common-scope** variant, exactly as for additive
  effects. A common variant is therefore **required** once any line variant exists.
  Otherwise an L×M heterozygote matches no dominance term, and heterosis disappears.

Recommended: **a follow-up release after Part C**, with the "common variant required"
rule enforced at write time. Sex-specific genetic effects are **not** a
`genome_effects` dimension. They are modelled the standard way, as two traits (e.g.
`BW_M`, `BW_F`) with a genetic correlation below 1 in `G`, which Part A already
calibrates exactly. `parent_origin` is something else: which parent an allele came from
(imprinting), not the individual's sex.

### Q15 — Where targets come from, and whether generators persist them *(decided 2026-10-01 → §6C)*

**Decision:** option (a), refined. `trait_var_comp` is the single source of truth. A passed
matrix is written, and is a hard error if the block already exists (even an identical one).
Stored rows are chosen with `trait_var_comp_tbl` (`NULL` = the stored rows), and the
surviving `effect_name`s decide which blocks are in the model, so "additive only" needs no
deletion. `define_trait()` loses `target_add_var` / `target_add_mean`. A number is taken as
a 1×1 matrix. The analysis that led there is kept below.


Review findings 5 and 6. Today `define_additive_effects(G =)` persists a passed `G`, and
`G = NULL` reads the stored `additive` target. In the source method `G_D = NULL` means
"no dominance". The draft §9.1 used the tidybreed meaning for all three blocks. That leaves
no way to say "additive only" for a trait that has a stored dominance target, and no way
to tell "absent" from "zero".

| Option | `NULL` means | Passed matrix | Consequence |
|---|---|---|---|
| **(a) The table is the single source of truth** ← recommended | the stored row, or absent if there is no row (all three blocks) | persisted, in the same transaction as the terms | every calibrated target is in `trait_var_comp` (*refined 2026-10-02: the table holds available targets, and a call may use a subset, §6C; §8's join compares only blocks the model has*), and `restore_pop()` has everything. Needs a way to delete a stored block (*superseded: `remove_rows()` already works on `trait_var_comp`, §6C*). `define_additive_effects()` on a trait with a stored `dominance` / `additive_by_additive` target errors, naming `define_genome_effects()` (from step 5, §6C) and the removal call |
| (b) Arguments decide, nothing persisted | `G_A`: stored; `G_D`, `G_AA`: absent | used for this call only; `define_effect_cov_matrix()` alone persists | No implicit configuration writes (the review's preference). But a passed target that differs from the stored one leaves `trait_var_comp` disagreeing with the model. Stored non-additive targets are never read by anything |
| (c) Explicit mode, `targets = c("stored", "arguments")` | per mode | per mode | Most explicit, one more argument, and two behaviours to test |

Under every option a zero matrix is an exact constraint and is stored as zeros, not
treated as absence. Recommendation **(a)**: it keeps one rule for all three blocks,
and it is the only option under which every calibrated target can be compared with what was delivered (§8).
It also matches what `define_additive_effects(G =)` and `define_trait(target_add_var =)`
already do.

### Q16 — How realised cross-block covariances are reported *(decided 2026-10-01: (a))*

**Decision: (a).** One `between_components` row per trait pair, as specified in §8. The
output keeps `trait_var_comp`'s shape, so the target-vs-delivered check stays one
`inner_join`. The reasoning: cross-block covariances are 0 in expectation under random
mating and noticeable mainly under inbreeding or strong recent selection, so the breakdown
is diagnostic detail. If someone needs the individual pairs later, that is a **separate
function** (e.g. `extract_genetic_covariance()`, option (b)'s shape), not an argument here:
an argument that changes the output shape is two functions (§8, `per_ind`). Adding it would
not change (a)'s rows. The analysis follows.

Review finding 2. Under LD or non-HWE, $\mathrm{Cov}(A,D)$, $\mathrm{Cov}(A,AA)$,
$\mathrm{Cov}(D,AA)$ (and covariances with `unpartitioned`) are nonzero, and
`(effect_name, trait_name_1, trait_name_2)` has no slot for a pair of effects.

| Option | Shape | Notes |
|---|---|---|
| **(a) One `between_components` row per trait pair** ← recommended | $\sum_{b\ne b'}\mathrm{Cov}(b_{t_1}, b'_{t_2})$ | Rows sum to `total` exactly. The output keeps `trait_var_comp`'s shape, so the join works unchanged. It is 0 under `"genic"` and is omitted there. The individual cross terms are not shown |
| (b) `effect_name_1` / `effect_name_2` columns | full block × block covariance, long | Complete, but the join to targets needs `filter(effect_name_1 == effect_name_2)` first, and the output is wider than users need |
| (c) Omit cross terms | blocks + `total` only | Simplest. Blocks visibly fail to sum, and users will ask why |

(a) was chosen. (b) is deferred to a possible separate function, never an argument on
`extract_genetic_variance()`.

### Q17 — A manual *Cockerham* path for models with A×A *(decided for the first release: (a), functional coding; (b) when someone needs it)*

Review finding 4 / §9.3 caveat. With pairs, the statistical $\alpha$ needs $e$, which
`ad_terms()` never sees. Options: (a) manual A×A models use functional coding only,
documented, with the total exact and the components functional ← **recommended for the
first release**; (b) one exported builder `nonadditive_terms(locus_name, a, d, pairs, p)`
that runs `.noia_to_stored()` (Q13) and returns the whole statistical term set; (c) make
`ad_terms()` accept pairs. (a) loses nothing measurable: §8 canonicalises either coding
(case 1), so `extract_genetic_variance()` gives the same answer. Only the `ind_tgv`
component split differs. Add (b) when a user needs Cockerham components from hand-set
coefficients. (c) overloads a single-locus builder.

### Q18 — New name for `formula_tbv`, and the DSL's optional arguments *(decided 2026-10-02 → §6A)*

**Name: `formula_tgv`.** After P2 the DSL reads `ind_tgv` (the total, or one of its
components) and `ind_tbv` no longer exists, so `formula_tbv` describes neither the table
nor the value. `formula_tgv` names the table it reads, stays in the TBV/TGV/EBV family §6B
keeps, and pairs with `formula` (which evaluates recorded phenotypes). Rejected:
`formula_genetic` / `genetic_formula` (more descriptive, but they lose the link to
`ind_tgv`), and renaming `formula` as well (churn for no gain). The argument and the
`phenotype_meta` column are renamed together (naming rule 1), with no alias, in the
consolidation release. The internal names go in the same change: `.FORMULA_TBV_DSL_FUNS`,
`.eval_formula_tbv()`, `.validate_formula_tbv()`, `.walk_formula_tbv_ast()`,
`.substitute_tbv_ast()`, `.build_tbv_env()`, `.ap_materialize_tbvs()`,
`.assemble_composite_tbv()`, `R/contributor_tbv.R` and its `.tbv_*` placeholders, and the
roxygen wording "composite TBV". At 0.71.1: 159 lines (183 occurrences, 22 files) across `R/`, `tests/`, `man/`
and one vignette (`vignettes/swine/`).

**DSL arguments** (the syntax gate PH2 promised but the plan had not specified):
- **Default is the total genetic value.** `"WWD + dam(WWM)"` is the calf's total `WWD` plus
  the dam's total `WWM`, additive and non-additive. A phenotype is the sum of all effects,
  so this is the normal case.
- **`component =`, named-only**, on `self()`, `dam()`, `sire()`, `group_sum()`,
  `group_mean()`: `dam(WWM, component = "additive")`. It is an opt-in for unusual cases
  (ignoring non-additive maternal effects, or comparing with an additive-only model).
  Values come from §6B's value vocabulary (`additive`, `dominance`, `indicator`,
  `interaction`). Anything else errors in `define_phenotype()`. A component the model
  has no terms for contributes 0 (§6A).
- **`table =`, named-only**, on `group_sum()` / `group_mean()`:
  `group_sum(ADG_social, pen_id, table = pen_history)`, default `ind_meta`. It **already
  exists** in the parser (`R/formula_helpers.R:201-216`), named or as the third positional
  argument. But it is undocumented (the roxygen says `ind_meta` only), untested, and
  unvalidated: `col` and `table` reach SQL text (`R/contributor_tbv.R:92-93`) without the
  `validate_sql_identifier()` the `components` route applies (`R/define_phenotype.R:505,
  511`). Fix all three. The positional form is removed, so no optional argument depends
  on position. `define_phenotype()` also checks that the table exists and has the column.
  `table.col` notation was rejected: R parses `pen_history.pen_id` as one symbol, and it
  would differ from the `components` route's separate `group_table` / `group_column`.

### Q19 — Per-line additive covariance targets *(decided 2026-10-01 → §6C)*

`define_additive_effects(line_name = "A", G = G_A)` can already be called once per line,
and after Part A each line-specific effect set would hit its own `G` exactly within its
line. What clashes:

1. **`trait_var_comp` has no `line_name`.** Each call overwrites the single stored `G` for
   the trait, so after three calls the table holds line C's target, labelled as the
   population's. This is an existing defect whenever `G` is passed per line.
2. **It needs line-specific effects**, which is a different allele effect in each line.
   That breaks "an allele means the same thing across populations" (§7.6), and those
   terms add no line mean differences.
3. **One common effect set cannot hit three lines' `G` at once.** Hitting several anchor
   populations simultaneously is the two-anchor problem (GSE Theorem 3, §13 out of
   scope). It is feasible only under restrictive conditions.

Options:

| Option | Notes |
|---|---|
| **(a) No per-line targets for common effects; add a `line_name` column to `trait_var_comp` (NULL = population-wide)** ← recommended | Fixes defect 1, and follows the CLAUDE.md schema-bias rule (reserve the line dimension). A line-scoped generator call reads and writes its line's row. Common effects keep one target, at the named reference. Lines' own covariances emerge from divergence (§7.6), and `extract_genetic_variance()` on each line measures them |
| (b) Also calibrate common effects to several lines at once | Theorem 3 machinery. Rarely feasible, and nobody needs it yet |
| (c) Leave as is, and document "last call wins" | Keeps a table that can misreport its own targets |

(a) is a small DDL change (pre-1.0, no alias).

**Decided 2026-10-01: (a), in Part A** (§6C). Under §6C a passed `G` is persisted, so
without the column every per-line call would collide with the population-wide block. A
line-scoped call reads its line's block and falls back to the `NULL` block, the package's
shared/default pattern.

### Q20 — Generation targets vs evaluation parameters *(open; mostly outside this plan)*

**Three different things are called "variance components" today:**

| Role | What it is | Stored? | Who reads it |
|---|---|---|---|
| **Generation target** | the covariance the effect generators calibrate to, at a named reference population | yes, in `trait_var_comp` (plus `phenotype_var_comp` for residual / random-effect noise that `add_phenotype()` **simulates**) | `define_*_effects()`, `add_phenotype()` |
| **Truth now** | the covariance a given line or generation actually has after selection and drift | **no**; derived (principle 6) | `extract_genetic_variance()` measures it |
| **Evaluation parameters** | what the breeder's genetic evaluation *assumes*, per line and over time: user-supplied, or REML estimates | **not yet**; there is no place for them | `add_ebv()` |

**Current defect.** `add_ebv()` builds the BLUPF90 parameter file from the generation
targets: `load_trait_cov(pop, "additive")` and the residual blocks of
`phenotype_var_comp` (`R/blupf90_helpers.R:325-326`). So every evaluation is told the
base-population truth, whichever line it runs on and however far selection has moved
the variances. Worse, `add_ebv(update_covars = TRUE)` is documented to write REML
estimates **back into `trait_var_comp`** (`R/add_ebv.R:61-66`, not implemented yet). That
would overwrite the generating targets. Writing into `phenotype_var_comp` would change
the residual noise `add_phenotype()` simulates. Estimates must never be written into
either table.

**Recommended:**

- A separate table for evaluation parameters, e.g. `eval_var_comp`, with columns
  `eval_model_name`, `effect_name` (`animal`, `maternal`, `residual`, random effects),
  `phenotype_name_1`, `phenotype_name_2`, `line_name` (NULL = all lines) and `cov_value`.
  It is keyed on **phenotypes**, because an evaluation models observed records, not the
  genetic components in `trait_var_comp`.
- One entry point, per the "one general function" preference:
  `define_effect_cov_matrix(..., eval_model_name = "maternal_2027")` writes an evaluation
  set. Without it, the call writes a generation target, as today.
- `add_ebv(eval_model_name = ...)` reads the named set. A REML run writes its estimates
  as a named set (new, or replacing one), never into the target tables. Keeping every
  REML round as a history, like `ind_ebv` stacks evaluations, is a sub-decision.
- Default when no set is named: use the generation targets as today, and say so in a
  message, so simple simulations keep working.

This is an `add_ebv()` / schema change, not part of Parts A–C. It constrains Q15, though:
under §6C, `trait_var_comp` holds **generation targets only**, and nothing else
writes to it.

### Q21 — How does the prevalence threshold know a target describes the active terms? *(decided 2026-10-03: (a) → §6A, §7.1, step 3)*

**The gap** (Codex review 2026-10-03, finding 1). The `prevalence` threshold is
`mean + qnorm(1 − prev)·sqrt(target + Ve)`. It is right only if the trait's terms
deliver `target`. Today, and under §6A's active-block rule as written, the code checks
only that a target row **exists**. Some calls write terms that ignore a stored target:
- `define_additive_effects(effects = )`;
- `scale_to_target = FALSE`;
- `define_genome_effect_terms()` terms next to generated ones.

Such a call leaves an old target row that the model does not deliver. A reproduction
has target 1, a delivered variance of about 553, `prevalence = 0.1`, and a realised
45.5%, with no message. §7.1's pre-decision wording made this deliberate ("a stored block
neither blocks them nor is checked"). Until step 3, the workaround is explicit `thresholds =`, now stated in
the `define_phenotype()` roxygen (0.72.3).

**Options:**

| Option | How calibration is proven | Cost |
|---|---|---|
| **(a) The owner is the provenance** ← **decided** | Generators write only calibrated terms. Manual `effects` and `scale_to_target` leave `define_additive_effects()`, as §9.1 already does for `define_genome_effects()` ("samples and calibrates only"). Exact coefficients go through `define_genome_effect_terms()` (`ad_terms()`) under a user owner. The threshold uses the targets only when **every** active term of the trait is `generated`. A user-owner term errors, naming `thresholds =`. One remaining hole: a target removed (`remove_rows()`) and rewritten after generation. Closed by `define_effect_cov_matrix()` refusing a genetic block whose traits already have `generated` terms of that kind. The fix is to regenerate with `G =`, which writes target and terms atomically (§7.1) | Breaking: reverses §7.1's named-table `effects` decision and removes `scale_to_target`. About 76 test call sites move to `define_genome_effect_terms()`. Must land in step 3 or later, because before consolidation `add_tbv()` reads only `generated` terms (§5) |
| (b) Check against delivered variance | At PLAN, recompute the genic variance of the active terms from their stored `center_value` and origin scope. Compare it to the target, and error on a mismatch | No API change. Wrong under `anchor = "realised"` (step 2): a correctly calibrated model fails a genic check. Dominance and A×A need step 4's Part B machinery, which comes after step 3. Needs a tolerance |
| (c) Store calibration provenance | A column on `genome_effects` (e.g. the target the term was calibrated to) | Stored derived metadata (CLAUDE.md principle 6), and a DDL change to the `genome_effect*` tables, which this plan otherwise avoids |

(a) gives "generated ⇔ calibrated to the stored target" by construction, with no
tolerance and no new metadata. It also makes `define_additive_effects()` and
`define_genome_effects()` the same kind of function.

**Decision (2026-10-03): (a).** Skipping calibration is only ever wanted for the user's
own numbers, and those belong in the writer. Applied in §0A, §6A, §7.1, Q2, steps 2–3 and
gates A12 (withdrawn), A19 and PH7.

---

### Q22 — Default of `warn_bounds` on small founder pools *(decided (b), 2026-10-03; built 0.73.1)*

§7.4 compares the calibrated covariance with what another population sees. With a
founder-pool base the comparison is the pool expectation `2 Cov(H)`, which carries the
pool's **sampling** LD. On pools of 20–100 haplotypes it routinely leaves `[0.8, 1.25]`
(a 100-haplotype pool with 200 QTL gave a spectrum of `[0.72, 1.05]`), so the default
call warns in most small simulations: 116 of the suite's 127 warnings at 0.73.0. The
warning is true (the founders `add_founders()` draws will see roughly that covariance),
but a warning that almost always fires gets ignored.

| Option | Effect |
|---|---|
| (a) Keep as built | Honest; noisy for small pools. Users pass `warn_bounds = NULL` |
| (b) Pool comparison as a `message()`, `warning()` only for an **observed** base (real individuals) or the realised anchor's genic limit | Keeps the information; warns only when real individuals are off |
| (c) Widen the default for the pool comparison by its expected sampling error (~`1/sqrt(n_haplotypes)`) | Warns only beyond sampling noise; more machinery |

Simulated fire rate of the ±25% default, 200 unlinked QTL (pure sampling LD):

| Haplotypes | 1 trait | 2 traits, r_g 0.4 |
|---|---|---|
| 20 | 48% | 89% |
| 50 | 29% | 68% |
| 100 | 10% | 45% |
| 200 | 4% | 10% |
| 500+ | 0% | 0% |

*(0.73.2: the rates above were measured with the sample divisor `n_h − 1`, which
overstates the pool expectation by `n_h/(n_h − 1)` — 5% at 20 haplotypes, 0.5% at 200.
The divisor is now `n_h`. The decision does not depend on the exact rates; not re-run.)*

**Decision (b).** The founder-pool comparison is always a `message()` giving the
relative spectrum and the pool size; outside `warn_bounds` it adds the fix ("add the
founders first and calibrate with `anchor = "realised"` on them", AlphaSimR's default
behaviour). The observed (`base_tbl` selects individuals) and genic-limit (realised
anchor) comparisons keep the `warning()`: there the user chose the animals, and the fix
is one argument. (c) stays possible later as a structure detector on top of (b). Step 5
reuses the same rule.

### Q23 — The reserved owner can still be written by hand *(decided 2026-10-04: (a); built 0.73.2)*

Q21 (a) makes "every active term is `generated`" the proof that the stored target
describes the model. Two public paths break that premise; one claimed path does not:

- **Open:** `define_genome_effect_terms(effect_owner = "generated",
  allow_reserved_owner = TRUE)` writes any coefficients under the reserved owner. The
  flag exists for the "generator == writer" test.
- **Open, already planned:** a target removed with `remove_rows()` and rewritten with
  `define_effect_cov_matrix()` — closed by the step-3 refusal in Q21.
- **Not a path:** removing a subset of generated *terms*. `remove_rows()` refuses all
  three `genome_effect*` tables (pinned in `test-genome-effects-schema.R`).
- **Not a prevalence issue:** `"union"` misses off-diagonals only; the threshold reads
  a diagonal, which `"union"` delivers exactly (verified since 0.73.2).

| Option | Effect |
|---|---|
| (a) Remove `allow_reserved_owner` from the exported writer; the test calls an internal entry point ← **recommended** | The owner becomes unforgeable through exported functions. Direct SQL remains possible, as for any table |
| (b) Keep the flag, and at PLAN verify the active model's genic variance against the target | Q21 (b)'s problems: wrong under `anchor = "realised"`, needs step 4 for non-additive terms, needs a tolerance |
| (c) Keep the flag, and document that it voids the prevalence guarantee | No code; relies on users reading it |

**Decision (a).** `define_genome_effect_terms()` loses `allow_reserved_owner`; it calls the
internal `.ge_write_terms()` with `FALSE`. The "generator == writer" test calls the
engine directly. `replace_trait` still refuses to delete generated terms, now with no
override: they are replaced only by re-running a generator.

## 13. Explicitly out of scope

- A×D, D×D, third-order epistasis. Storable today; the calibrator's extra NOIA blocks are "mechanical but unwritten".
- Two anchors at once (Chu's pPG; purebred-for-crossbred): GSE Theorem 3 per component.
- Line- or parent-of-origin-scoped **non-additive** generation. The writer can store it; the calibrator assumes one base per locus.
- A realised-coding `dominance` contrast (DDL).
- Enforced inbreeding-depression targets for non-diagonal $\mathbf G_D$ (reported instead).
- Non-diploid dominance.

---

## 14. External dependencies before anything is *cited*

The source project has three unreported upstream findings: AlphaSimR `calcGenParamE`
locus-index typo; MoBPS `var.*.l` are SDs; XSim scalar-`vg` `TypeError`. None blocks the
implementation. All three should be reported to maintainers before tidybreed docs or a
vignette cite comparisons with those packages.
