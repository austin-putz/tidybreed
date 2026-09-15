# tidybreed — Package Summary

**Version:** 0.65.0 · **Snapshot date:** 2026-09-10 · **Branch:** `feat/genome-effects-v49`

A database-first (DuckDB) breeding-program simulator in R. All genomic and
individual data lives in a file-based DuckDB database rather than R memory,
enabling simulations larger than RAM, resumable runs, and replicate archiving.
A small Rcpp kernel handles meiosis/recombination.

## Highlights

| | |
|---|---|
| Exported functions | 39 (+17 S3 methods; 213 functions total incl. internals) |
| R source | 55 files, ~18,200 lines (~4,800 of which are roxygen docs) |
| C++ (Rcpp) | 2 files, ~220 lines (gamete/recombination kernel) |
| Documentation | 176 man pages, 1 vignette, 1,100-line README |
| Tests | 53 testthat files, ~15,100 lines, 715 tests, ~1,600 assertions |
| Test : source ratio | 0.83 : 1 |
| History | 190 commits, 122 released versions in NEWS.md |
| Dependencies | duckdb, DBI, dbplyr, dplyr, tibble, rlang, cli, dqrng, MASS, Rcpp |

## Detailed Counts

| Metric | Count |
|---|---:|
| Exported functions (`NAMESPACE`) | 39 |
| S3 methods registered | 17 |
| Total R function definitions | 213 |
| R source files (`R/`) | 55 |
| R lines of code | 18,152 |
| Roxygen doc lines (`#'`) in `R/` | 4,817 |
| C++ source files (`src/`) | 2 |
| C++ lines of code | 217 |
| Man pages (`man/*.Rd`) | 176 |
| testthat test files | 53 |
| testthat lines of code | 15,115 (+ 3 helper files, 969 lines) |
| `test_that()` blocks | 715 |
| `expect_*()` assertions | 1,637 |
| Vignettes | 1 (`tidybreed-introduction.Rmd`, 645 lines) |
| `NEWS.md` | 3,268 lines, 122 version headings |
| `README.md` | 1,145 lines |
| Design docs (`plans/`) | 49 files, 21,269 lines |
| Git commits | 190 |

## Exported API (39 functions)

| Prefix | Functions |
|---|---|
| `open_` / `restore_` / `close_` | `open_pop`, `restore_pop`, `close_pop` |
| `define_` | `define_genome`, `define_chromosome`, `define_founder_haplotypes`, `define_chip`, `define_trait`, `define_trait_simple`, `define_additive_effects`, `define_phenotype`, `define_residual_cov`, `define_effect_cov_matrix`, `define_effect_random`, `define_effect_fixed_class`, `define_effect_fixed_cov`, `define_effect_intercept`, `define_index`, `define_table`, `define_schema_description` |
| `add_` | `add_founders`, `add_offspring`, `add_phenotype`, `add_tbv`, `add_ebv`, `add_index`, `add_dosage`, `add_genotypes` |
| `mutate_` | `mutate_table`, `mutate_derived`, `mutate_group_seq`, `mutate_group_named`, `mutate_group_concatenate` |
| `extract_` / `remove_` / `archive_` | `extract_genotypes`, `remove_rows`, `archive_replicate` |
| Inspection | `get_table`, `schema`, `describe_table` |

## Largest Source Files

| File | LOC | Purpose |
|---|---:|---|
| `schema.R` | 1,520 | Table registry, descriptions, `schema()` / `describe_table()` |
| `add_phenotype.R` | 1,146 | Phenotype simulation (composite, SGE, fixed/random effects, residuals) |
| `add_offspring.R` | 823 | Mating, gamete formation, offspring haplotype writes |
| `mutate_table.R` | 753 | Generic typed column add/update on any table |
| `formula_helpers.R` | 703 | Formula-based derived phenotypes |
| `define_additive_effects.R` | 645 | QTL effect sampling, Falconer rescale, multi-trait MVN |

## Largest Test Files

| Test file | LOC |
|---|---:|
| `test-define_founder_haplotypes.R` | 803 |
| `test-mutate_group.R` | 763 |
| `test-define_index.R` | 644 |
| `test-formula_phenotype.R` | 634 |
| `test-genome-effects-schema.R` | 621 |
| `test-phenotype_composite.R` | 590 |

## Dependencies

- **Imports:** cli, dplyr, tibble, rlang, duckdb, DBI, dbplyr, dqrng, MASS, Rcpp, stats
- **LinkingTo:** Rcpp, dqrng
- **Suggests:** testthat (≥ 3.2.0), withr, knitr, rmarkdown, covr
- **License:** MIT
