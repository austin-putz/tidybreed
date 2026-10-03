# List of TODO Items

## From `plans/update_genome_effects_base_tbl.md` (v0.69.0, 2026-09-20)

Logged there as §3.5; deliberately not implemented in that change.

- **`add_ebv()`: rename `phenotype` → `phenotype_tbl`** for `*_tbl` consistency
  with every other second-table argument (`base_tbl`, `effects_tbl`,
  `loci_tbl`), and replace its `collect()`-the-whole-filtered-table-to-pull-ids
  (`R/add_ebv.R`, "Validate phenotype filter") with the same rendered-subquery
  approach `extract_allele_freq()` uses. Untouched because the BLUPF90 paths
  have no CI coverage.
- **`genome_meta.founder_allele_freq`**: decide between dropping it, making it
  line-keyed, or documenting its last-call-wins meaning. It is written by
  `define_founder_haplotypes()` and read by no writer (all base frequencies
  now come from `extract_allele_freq()`); in a multi-line population it
  describes only the pool written last. See the plan's Q8.
- **Unify the two genome-effect writers' vocabulary** (plan §2.5). Scope is
  `line_name` + `parent_origin` on `define_additive_effects()` but
  `origin = list(...)` on `define_genome_effect_terms()`; replacement is an
  implicit `replace_scope` on the generator but `mode =` on the writer. One
  scope spelling — either the generator accepts `origin =` with
  `line_name`/`parent_origin` kept as sugar that composes into one origin row,
  or the writer's `origin` accepts the short form
  `list(line_name = , parent_origin = )` — and `mode =` exposed on the
  generator. Separate change; it touches argument surfaces the `base_tbl` plan
  did not. The "generator == writer" test in
  `tests/testthat/test-genome-effects-writer.R` is the guard for whichever
  direction is taken.

## From `plans/import_qtl_effect_methods.md` §7.6 (2026-09-25)

Deferred on purpose. The documented route to line differences today is divergent
selection from a common founder pool, continued in the same database.

- **Target line means through founder frequencies.** A tool that derives a line's
  founder pool from a reference pool plus target mean differences, for one or several
  traits: the minimum-norm frequency shift $\delta_j \propto \alpha_j$ solving
  $\sum_j 2\alpha_j\delta_j = \Delta$ (a $k$-equation system for several traits;
  nonlinear once dominance is present), clipped to $[0, 1]$, with an error when the
  target needs impossible frequencies. Neutral loci get background differentiation
  (Fst) too, so the QTL do not stand out. It must run after the effects are defined,
  which reverses the usual call order. Open: a new function or a new
  `define_founder_haplotypes()` method.
- **Seed a new population from another simulation's end state.** Copy the genome, the
  genome-effect tables (with their `center_value`s) and the haplotypes of chosen
  individuals at generation *t* into a new database as its founder pool. Today there is
  no supported import path. The interim routes are (a) copy the `.duckdb` file and
  continue in the copy (everything stays consistent, pedigree included), or (b) export
  the tables and inject them with DBI. Route (b) bypasses the writers' validation, so
  it needs: identical `genome_meta` / `locus_id`s, ids from `next_int_id()`, the three
  `genome_effect*` tables copied together (then `validate_genome_effects()`), and the
  selected animals' `ind_haplotype` rows becoming `founder_haplotypes` rows with their
  `line_name`.
- **Per-line mutation** over generations (rates may differ by line).
- **Separate evaluation parameters from generation targets** (plan Q20). `add_ebv()`
  feeds BLUPF90 the generating `trait_var_comp` / `phenotype_var_comp` values. Its
  `update_covars` is documented to write REML estimates back into `trait_var_comp`,
  which would overwrite the simulation's targets. Proposed: an `eval_var_comp` table
  keyed by `eval_model_name` and `line_name`, and `add_ebv(eval_model_name =)`.
