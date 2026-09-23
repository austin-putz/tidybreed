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
  `origin = list(...)` on `define_genome_effects()`; replacement is an
  implicit `replace_scope` on the generator but `mode =` on the writer. One
  scope spelling — either the generator accepts `origin =` with
  `line_name`/`parent_origin` kept as sugar that composes into one origin row,
  or the writer's `origin` accepts the short form
  `list(line_name = , parent_origin = )` — and `mode =` exposed on the
  generator. Separate change; it touches argument surfaces the `base_tbl` plan
  did not. The "generator == writer" test in
  `tests/testthat/test-genome-effects-writer.R` is the guard for whichever
  direction is taken.
