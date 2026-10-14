# List of TODO Items

Last reviewed **2026-10-12** against the v0.76.4 working tree, after the QTL-effect
import plan (`plans/import_qtl_effect_methods.md`, steps 0–5) was finished. This is
the active backlog; completed implementation/review details stay in their original
plans. Proposed APIs below are proposals, not existing functions.

## Recommended next steps

0. **Write `plans/check_policy.md`** (section below): the check policy and the
   inventory of refusals, then the chromosome / genome-effect checks as its first use.
1. **Finish the small documentation and release checks below.** The new effect
   generators and variance extractor are implemented; another genetics expansion
   is not needed to start documenting the package.
2. **Plan evaluation parameters separately from generation targets (Q20).** This is
   the largest remaining design issue. Include the `add_ebv()` argument cleanup and
   tests of evaluation input files; finish it before teaching configurable genetic
   evaluation as a stable workflow.
3. **Decide the remaining public-API questions before writing those recipes.** In
   particular, generator RNG policy and scope/replacement vocabulary. A decision
   to retain the current vocabulary is sufficient; unification is not mandatory.
4. **Build the first runnable Quarto chapters now.** Start with a small complete
   breeding cycle, genetic models, and target-versus-measured variance. Expand into
   the full cookbook and methods book using [the documentation proposal](documentation.md).

## Check policy: refuse by default, `allow =` to bypass (next plan)

Agreed with the user on 2026-10-13, while reviewing
`plans/line_specific_chr_rules_add_offspring.md` (its §6). The next plan is
**`plans/check_policy.md`**: a written policy plus an inventory of every refusal in
`R/`, reviewed by the user and Codex before any code.

**Decisions already made:**

- **Users may model anything.** Unusual genomes, chromosomes and custom effects are
  allowed. By default the package refuses or warns loudly when an assumption breaks.
  The user can relax a named check, and from then on the package computes what is
  well-defined and accepts what they did.
- **Three kinds of check:**
  1. **Integrity and safety. Never bypassable.** DB consistency and orphans, id rules,
     determinism, the failure/rollback contract, the closed formula grammar (`eval()`
     security). The 0.76.10 `define_chromosome()` guard is here too: if bypassed,
     stored haplotypes would disagree with the rules.
  2. **Not implemented. Not a check.** For example, `2,0` transmission in
     `add_offspring()`. Implement it, or define a behavior before offering a bypass.
  3. **Modelling assumptions. Bypassable.**
- **A check is bypassable only when every later function has a defined, documented
  behavior for the state it lets through.** Otherwise the bypass only moves the crash.
- **Mechanism.** One argument, named the same across functions: `allow =`, a character
  vector of named relaxations (e.g. `allow = "non_diploid_dominance"`). Every refusal
  message ends with the exact value to pass. Not a blanket `check = FALSE`, and not a
  global option.
- **Default level.**
  - **Refuse** when the result would be silently wrong: calibration, a missing genotype
    state, a variance target.
  - **Warn** when the result is literally what was asked for, only unusual (e.g.
    additive effects on a hemizygous chromosome, counted per copy, which is correct).
- **Later code that trusts an earlier check re-checks the condition itself** and needs
  its own `allow =`. A relaxation is not stored in the database. Example:
  `define_phenotype(prevalence = )` trusts calibration of `"generated"` terms (Q21).
  If a generator ran with the diploid check relaxed, prevalence re-checks and refuses
  without its own `allow =`.

**First use: the chromosome / genome-effect checks.** They read only the line-agnostic
`chr_inheritance` rules today:

- `assert_qtl_autosomal()` (generators)
- `.ge_validate_dominance_ploidy()` (every genome-effect write)
- `.ge_reachable_copy_counts()` (`copy_count_value` inference, `require_complete`)

With a line rule that makes a chromosome hemizygous in one line, these approve effects
that are then wrong:

- dominance fails late in `add_tgv()`;
- generated-QTL calibration misses its target, silently;
- indicator members get no effect, silently;
- `require_complete` passes an incomplete surface.

Separately, `define_chromosome()` never re-runs these checks after it writes a rule.
So "effects first, then `define_chromosome()`" gets through even for line-agnostic
X/Y rules.

The follow-up covers:

- resolve copy counts over every (sex, line) context with rules;
- diploid means "two copies in every context" (recombination is irrelevant);
- reachable copy counts are the union over contexts;
- dominance keeps its exact two-copy-origin exception;
- `define_chromosome()` re-validates after a write;
- `allow =` values for each check;
- `add_tgv()` / `add_phenotype()` `allow =` to skip dominance for animals with a copy
  count other than 2, with a message giving the counts;
- tests for both configuration orders.

## Documentation and release cleanup

- **Audit the README against the current API, then shorten it to an entry point.**
  Concrete drift: the badge says 0.71.1; the founder-frequency note still teaches
  removed `base = "current_pop"` / `base = "founder_haplotypes"` arguments
  (`README.md`, near lines 490–501). Its three-trait additive block is later generated
  in separate one-/two-trait calls, which the current complete-block validation
  rejects. Replace this with one executed minimal example and links to the
  introduction, genetic-models article, book, and API reference.
- **Fix `inst/CITATION`.** (Its version now comes from `meta$Version`, 0.76.4.) Its
  second entry points to the nonexistent
  `vignettes/qmd/tidybreed-quick-start.qmd` and describes it as distributed with
  the package. Point citations to actual maintained documentation; add a book
  citation once the book is published. Do not promise a DOI before one exists.
- **Verify the tightened CI gate on GitHub.** v0.76.4 fixed the non-ASCII
  warning and the build notes (license stub, imports, partial match, ignored
  files) and switched `.github/workflows/R-CMD-check.yaml` to `error-on:
  '"warning"'`. Confirm the Linux/macOS matrix passes on the next push.
- **Keep generated/internal documentation out of the public site.**
  `package_summary.md` is a generated snapshot (currently stamped 0.74.0) and the
  pkgdown workflow strips only `CLAUDE`. Decide where snapshots belong and remove
  internal pages from deployment output. Keep `tools/quarto/legacy/` out of the
  new book: its quick start is v0.27.0 and teaches removed functions/tables.
- **Create maintained documentation with executable examples.** Keep pkgdown for
  function reference and short package vignettes; add `book/` for the Quarto
  cookbook, methods, validation, and bibliography. Add it to `.Rbuildignore`,
  check its examples on PRs, and coordinate publishing with the existing pkgdown
  deployment. Details and chapter order: [documentation.md](documentation.md).
- **Preserve schema descriptions in archives, and fix the restore message.**
  A new archive made with the default `archive_replicate()` settings omits
  `_schema_meta`. `restore_pop()` then warns that it "predates schema metadata
  (pre-v0.36.0)", even though it was just written by v0.76.4. Reproduced in the
  file-backed selection/restore/archive smoke run. Decide how descriptions travel
  to the archive (including user columns), and add a focused round-trip check.
  Explain that the default archive is for analysis: haplotypes/genotypes are
  reset-only, so reopening it is not the same as resuming the working population.
- **Run the release checks on the final candidate.** Full test suite, package build
  and `R CMD check`, `pkgdown::check_pkgdown()`, executed vignettes/book, and a
  small file-backed multi-generation example including close/restore and archive.
  Include a documented BLUPF90 smoke run separately from binary-free checks.
  Windows is not in the current Linux/macOS CI matrix; decide whether to add it
  before claiming tested Windows support.
- **Clean up test warnings and explicitly run the non-CRAN checks before release.**
  This review's full suite passed 4,843 assertions with no failures/errors, but
  emitted six warnings. Investigate the uninitialised `hap_id` column warning in
  `test-add_founders.R`; explicitly expect the intended repeat-phenotype,
  monomorphic-pool and pooled-line warnings in the other tests. Five tests were
  CRAN-gated skips: one founder performance check and four internal dual-anchor
  congruence checks. Run those separately with `NOT_CRAN=true` when validating the
  release or methods evidence; their presence does not expose a public dual-anchor API.

## Evaluation parameters and `add_ebv()`

- **Separate evaluation parameters from generation targets** (import plan Q20 — the most
  important open design issue). `add_ebv()` builds the BLUPF90 parameter file from the
  generation targets (`load_trait_cov(pop, "additive", trait)` and the residual blocks of
  `phenotype_var_comp`, `R/blupf90_helpers.R`), so the default evaluation uses
  population-wide generation targets whatever line or generation it runs on. `update_covars` is
  documented to write REML estimates back into `trait_var_comp`; writeback is not
  implemented, but the message it prints today tells the user to update
  `trait_var_comp` by hand with `define_effect_cov_matrix()`, which contradicts the
  design (targets are generation-only and never overwritten; that function refuses a
  stored block anyway). Proposed: an `eval_var_comp` table keyed by `eval_model_name`,
  `effect_name`, `phenotype_name_1/2` and `line_name`; `define_effect_cov_matrix(...,
  eval_model_name =)` to write a set; `add_ebv(eval_model_name =)` to read one; REML
  estimates written as a named set, never into the target tables; default to the
  generation targets with a message. Schema change: needs its own plan and review.
- **`add_ebv()`: rename `phenotype` → `phenotype_tbl`** for `*_tbl` consistency with every
  other second-table argument (`base_tbl`, `effects_tbl`, `loci_tbl`), and replace its
  `collect()`-the-whole-filtered-table-to-pull-ids (`R/add_ebv.R`, "Validate phenotype
  filter") with the rendered-subquery approach `extract_allele_freq()` uses. Untouched
  because full BLUPF90 solver runs have no CI coverage. Routing, data-file building
  and solution parsing do have tests; add binary-free parameter/filter tests and
  a separately run solver smoke test. Natural to do with Q20.
  (From `plans/update_genome_effects_base_tbl.md` §3.5, v0.69.0.)

## Genome effects — future extensions, not prerequisites for the book

- **Line-specific non-additive effects** (import plan Q14). `define_genome_effects(line_name
  =)`, calibrated within each line on that line's founder pool, with additive members
  `exact L`, dominance members {`exact L`, `exact L`} and A×A members `exact L` on each.
  Crossbred copies fall back to the common variant, so a **common variant is required**
  once any line variant exists (enforced at write time), or an L×M heterozygote matches
  no dominance term and heterosis disappears. No DDL. Parent-of-origin non-additive
  effects belong with it.
- **`extract_breeding_values()`** (import plan §8). Generation-t breeding values per
  individual, for the supported common-scope case only, naming its estimator under LD
  (the NOIA contrast component vs the joint least-squares projection `lm(g ~ dosages)`).
  A scoped model needs a reference-population projection nobody has specified.
- **`extract_functional_effects(pop, trait_name)`** (Q9, low priority). Report the
  functional `(a, d, e)` AlphaSimR users think in, recovered from the stored terms and
  their centres. Nothing stored.
- **Architecture samplers and `marginal =`** (Q5). The source's seven samplers
  (`gaussian`, `mvt`, `gamma`, `laplace`, `sparse`, `major`, `maf`) and its rank-restoring
  iteration for non-elliptical marginals at k > 1. New samplers go into
  `.draw_additive_architecture()` so both generators keep identical additive draws
  (gate C4 (b)). Replacing `distribution` with an `architecture` argument is an API change
  of its own. Supplied architectures at the public level belong here too.
- **A whole-model Cockerham builder** (Q17 (b)): `nonadditive_terms(locus_name, a, d,
  pairs, p)` running `.noia_to_stored()` and returning the full statistical term set, for
  users who want Cockerham components from hand-set coefficients with A×A. Today manual
  A×A models use functional coding (documented). Only when someone needs it.
- **A pool-based `"realised"` anchor** (Q11). `M_pool = 2 Cov(H)` over
  `founder_haplotypes` (LD-aware, available before `add_founders()`) as `anchor =
  "realised"` when `base_tbl` is the pool, instead of the current error. Diagnostics use it
  already; cheap to add if asked for.
- **`extract_genetic_covariance()`** (Q16 (b)). The full block × block covariance (A with
  D, A with A×A, …) per trait pair, as a separate function; `extract_genetic_variance()`
  keeps its single `between_components` row.
- **Out of scope for now** (import plan §13), listed so they are not forgotten: A×D, D×D
  and third-order epistasis (storable today; the calibrator's extra NOIA blocks are
  unwritten); two anchors at once (Chu's pPG, purebred-for-crossbred); a realised-coding
  `dominance` contrast (DDL); enforced inbreeding-depression targets for non-diagonal
  `G_D`; non-diploid dominance.

## Genome-effect writers and schema

- **Settle generator RNG policy before the reproducibility chapter** (Q12).
  Recommended: remove `seed =` from `define_additive_effects()` and teach
  `set.seed()` explicitly, as with `define_genome_effects()`. The current argument
  calls `set.seed()` and resets the caller's RNG stream. `add_offspring(seed =)`
  seeds its own `dqrng` sub-stream and is a different contract.
- **Unify the two genome-effect writers' vocabulary** (base_tbl plan §2.5). Scope is
  `line_name` + `parent_origin` on `define_additive_effects()`, but `origin = list(...)`
  on `define_genome_effect_terms()`. `define_genome_effects()` currently has neither
  scope argument and replaces the traits' whole generated model; the additive
  generator replaces its scope, while the writer exposes `mode =`. Preserve those
  different replacement meanings even if the vocabulary is unified. One scope spelling —
  either the additive generator accepts `origin =` with `line_name` / `parent_origin` kept as
  sugar that composes into one origin row, or the writer's `origin` accepts the short form
  `list(line_name = , parent_origin = )` — and `mode =` exposed on the generators if
  wanted. The "generator == writer" test in `tests/testthat/test-genome-effects-writer.R`
  guards whichever direction is taken.
- **`genome_meta.founder_allele_freq`**: decide whether to drop this unused column
  before 1.0 or retain it explicitly as informational. Its last-call-wins meaning
  is already documented in the README (the adjacent base-selection advice needs
  the cleanup above). Written by `define_founder_haplotypes()`
  (`R/founder_haplotype_helpers.R`) and read by nothing (all base frequencies come from
  `extract_allele_freq()`); in a multi-line population it describes only the pool
  written last. (base_tbl plan Q8.)

## Lines and populations

Deferred on purpose in import plan §7.6. The documented route to line differences today
is divergent selection from a common founder pool, continued in the same database.

- **Target line means through founder frequencies.** A tool that derives a line's founder
  pool from a reference pool plus target mean differences, for one or several traits: the
  minimum-norm frequency shift $\delta_j \propto \alpha_j$ solving
  $\sum_j 2\alpha_j\delta_j = \Delta$ (a $k$-equation system for several traits; nonlinear
  once dominance is present), clipped to $[0, 1]$, with an error when the target needs
  impossible frequencies. Neutral loci get background differentiation (Fst) too, so the
  QTL do not stand out. It must run after the effects are defined, which reverses the
  usual call order. Open: a new function or a new `define_founder_haplotypes()` method.
- **Seed a new population from another simulation's end state.** Copy the genome, the
  genome-effect tables (with their `center_value`s) and the haplotypes of chosen
  individuals at generation *t* into a new database as its founder pool. No supported
  import path today. Interim routes: (a) copy the `.duckdb` file and continue in the copy
  (everything stays consistent, pedigree included), or (b) export the tables and inject
  them with DBI, which bypasses the writers' validation and so needs identical
  `genome_meta` / `locus_id`s, ids from `next_int_id()`, the three `genome_effect*` tables
  copied together (then `validate_genome_effects()`), and the selected animals'
  `ind_haplotype` rows becoming `founder_haplotypes` rows with their `line_name`.
- **Mutation, and per-line mutation rates** over generations. There is no mutation in the
  package yet.
- **Per-gamete fast-path routing for line-specific chromosome rules (option C).**
  - **Current rule** (0.76.12, `plans/line_specific_chr_rules_add_offspring.md`, D1 =
    option A): a chromosome that *any* line's rule makes special goes through the
    per-offspring R path for *every* line.
  - **Two costs of that rule:**
    - Every line runs slower. In one benchmark (10 chromosomes, 10k loci, 2,000
      offspring), the first special chromosome added ~10 s to 15.7 s, and each
      further one ~3.5 s.
    - A rule for line B changes line A's seeded offspring from then on, so paired
      scenario comparisons lose common random numbers.
  - **Option C:** choose the path per gamete from its own context (offspring sex/line,
    parent sex/line), so a rule affects only the animals it names. The C++ kernel
    probably stays as it is, since it already takes chromosome bounds per gamete group.
    The parent matrices, per-gamete row counts and batch assembly in `add_offspring()`
    need rework, and hemizygous parents need placeholder rows.
  - **When to do it:** when a real cross-strain or cross-species use case needs the
    speed or isolated seeds.

## Accepted limitations to explain in the book

These are documented decisions, not unfinished implementation tasks. Reopen them
only when a concrete use case needs different behaviour.

- Removing a trait's final terms leaves its old `ind_tgv` rows. `add_tgv()` and
  `add_phenotype()` refuse an empty model; direct queries can still see those old
  rows. Explain this in the replacement/removal recipe (5c review, Option A).
- `extract_genetic_variance()` can lose relative precision for an additive block
  many orders smaller than dominance/epistatic coupling (5b review). Explain the
  numerical limit; reading stored `alpha` directly at the matching projection
  frequency is a possible later refinement.
- Under LD, the extractor's additive NOIA contrast component is not the joint
  least-squares breeding-value projection. Name the measured quantity precisely;
  `extract_breeding_values()` remains a separate proposed feature.
