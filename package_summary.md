# tidybreed — Package Summary

**Version:** 0.68.3 · **Snapshot date:** 2026-09-16 · **Branch:** `feat/genome-effects-v49` (`3917890`)

A database-first (DuckDB) breeding-program simulator in R. All genomic and
individual data lives in a file-based DuckDB database rather than R memory,
enabling simulations larger than RAM, resumable runs, and replicate archiving.
A small Rcpp kernel handles meiosis/recombination.

## Highlights

| | |
|---|---|
| Exported functions | 43 (+17 S3 methods; 275 functions total incl. internals) |
| R source | 59 files, ~20,700 lines (~5,600 of which are roxygen docs) |
| C++ (Rcpp) | 2 files, ~220 lines (gamete/recombination kernel) |
| Documentation | 180 man pages, 1 vignette, 1,150-line README |
| Tests | 55 testthat files, ~17,300 lines, 794 tests, ~1,860 assertions |
| Test : source ratio | 0.83 : 1 |
| History | 198 commits, 128 released versions in NEWS.md |
| Dependencies | duckdb, DBI, dbplyr, cli, dplyr, tibble, rlang, dqrng, MASS, Rcpp |

## Detailed Counts

| Metric | Count |
|---|---:|
| Exported functions (`NAMESPACE`) | 43 |
| S3 methods registered | 17 |
| Total R function definitions | 275 |
| R source files (`R/`) | 59 |
| R lines of code | 20,671 |
| Roxygen doc lines (`#'`) in `R/` | 5,631 |
| C++ source files (`src/`) | 2 |
| C++ lines of code | 217 |
| Man pages (`man/*.Rd`) | 180 |
| testthat test files | 55 |
| testthat lines of code | 17,260 (+ 4 helper files, 999 lines) |
| `test_that()` blocks | 794 |
| `expect_*()` assertions | 1,864 |
| Vignettes | 1 (`tidybreed-introduction.Rmd`, 649 lines) |
| `NEWS.md` | 3,595 lines, 128 version headings |
| `README.md` | 1,149 lines |
| Design docs (`plans/`) | 53 files, 22,798 lines |
| Git commits | 198 |

## Exported API (43 functions)

| Prefix | Functions |
|---|---|
| `open_` / `restore_` / `close_` | `open_pop`, `restore_pop`, `close_pop` |
| `define_` | `define_genome`, `define_chromosome`, `define_founder_haplotypes`, `define_chip`, `define_trait`, `define_trait_simple`, `define_additive_effects`, `define_genome_effects`, `define_phenotype`, `define_residual_cov`, `define_effect_cov_matrix`, `define_effect_random`, `define_effect_fixed_class`, `define_effect_fixed_cov`, `define_effect_intercept`, `define_index`, `define_table`, `define_schema_description` |
| `add_` | `add_founders`, `add_offspring`, `add_phenotype`, `add_tbv`, `add_tgv`, `add_ebv`, `add_index`, `add_dosage`, `add_genotypes` |
| `mutate_` | `mutate_table`, `mutate_derived`, `mutate_group_seq`, `mutate_group_named`, `mutate_group_concatenate` |
| `extract_` / `remove_` / `archive_` | `extract_genotypes`, `remove_rows`, `archive_replicate` |
| Term builders | `ad_terms`, `genotype_terms` |
| Inspection | `get_table`, `schema`, `describe_table` |

## Largest Source Files

| File | LOC | Purpose |
|---|---:|---|
| `schema.R` | 1,520 | Table registry, descriptions, `schema()` / `describe_table()` |
| `add_phenotype.R` | 1,146 | Phenotype simulation (composite, SGE, fixed/random effects, residuals) |
| `define_genome_effects.R` | 917 | General genome-effect writer: terms, members, origins, replace modes |
| `define_additive_effects.R` | 912 | QTL effect sampling, Falconer rescale, multi-trait MVN, line/parent-origin scope |
| `genome_effects_eval.R` | 874 | The one evaluator behind `add_tbv()` / `add_tgv()` |
| `add_offspring.R` | 823 | Mating, gamete formation, offspring haplotype writes |

## Largest Test Files

| Test file | LOC |
|---|---:|
| `test-genome-effects-writer.R` | 991 |
| `test-genome-effects-eval.R` | 861 |
| `test-define_founder_haplotypes.R` | 803 |
| `test-mutate_group.R` | 763 |
| `test-genome-effects-schema.R` | 674 |
| `test-define_index.R` | 644 |

## Dependencies

- **Imports:** cli, dplyr, tibble, rlang, duckdb, DBI, dbplyr, dqrng, MASS, Rcpp, stats
- **LinkingTo:** Rcpp, dqrng, BH, sitmo
- **Suggests:** testthat (≥ 3.2.0), withr, knitr, rmarkdown, covr
- **SystemRequirements:** C++17
- **License:** MIT
