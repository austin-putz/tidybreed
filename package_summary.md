# tidybreed — Package Summary

**Version:** 0.68.4 · **Snapshot date:** 2026-09-16 · **Branch:** `feat/genome-effects-v49` (`a190bf4`)

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
| History | 199 commits, 129 released versions in NEWS.md |
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
| `NEWS.md` | 3,607 lines, 129 version headings |
| `README.md` | 1,149 lines |
| Design docs (`plans/`) | 53 files, 22,798 lines |
| Git commits | 199 |

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

## Database Tables

27 tables and 3 views in 8 groups, 201 columns in total, 201 of them described in `_schema_meta` (`schema()`, `describe_table()`).

Column types: VARCHAR ×111, INTEGER ×40, DOUBLE ×32, UTINYINT ×11, BIGINT ×3, BOOLEAN ×3, DATE ×1.

Keys: 22 tables declare a SQL `PRIMARY KEY`; 5 use a logical key enforced in R (`TABLE_ROW_KEYS`), by design where a DuckDB constraint would block bulk inserts or transactional replacement. `Archive` is how `archive_replicate()` treats the table: copied and stamped *per replicate*, copied *once*, *reset only*, or *kept* in the working database.

### Genome (9)

| Table | Kind | Cols | Key | Created by | Archive | Description |
|---|---|---:|---|---|---|---|
| `genome_meta` | table | 6 | `locus_id` | `define_genome()` | once | Locus-level metadata. One row per locus. Stores chromosome assignment, physical position (pos_bp), and chip membership flags added by define_chip(). Genetic-map positions live in the separate genome_map table. |
| `genome_map` | table | 7 | `id_genome_map` (logical) | `define_genome()` | kept | Genetic map in long format. One row per (locus x sex x line x map) with a defined genetic position (pos_cM). sex NULL = both sexes; line_name NULL = all lines. define_genome() writes a single default map (map_name 'default'); sex/line/version-specific maps are added as rows. Logical key (locus_id, sex, line_name, map_name). |
| `genome_effects` | table | 5 | `id_genome_effect` | `define_genome()` | once | Genome-effect terms. One row per term: a single coefficient over one or more loci. The loci are in genome_effect_members; scope (line, parent of origin) is in genome_effect_member_origins. A locus is causal for a trait if it appears as a member of any term. |
| `genome_effect_members` | table | 7 | `id_genome_effect`, `member_slot` | `define_genome()` | once | The loci a term spans, one row per (term x locus), canonicalized by ascending locus_id. contrast_name says which basis function the locus contributes: 'additive' (per allele copy), 'dominance' (Cockerham, diploid loci), or 'indicator' (one genotype state). |
| `genome_effect_member_origins` | table | 7 | `id_genome_effect`, `member_slot`, `origin_slot` | `define_genome()` | once | Scope of a member, as a predicate over allele copies. No rows = the common scope, which matches every copy. An additive member takes at most one row; a genotype member takes an exact multiset whose copy counts sum to the realized copy count. Among matching variants of one family the most specific wins; overlapping but incomparable scopes are refused at write time. |
| `genome_effect_terms` | view | 9 | — | `define_genome()` | kept | View: one row per term, with effect_order (member count) and family_key derived rather than stored. Terms sharing a family_key are scope variants of one mathematical term and compete under specificity fallback; terms with different keys sum. |
| `genome_effect_loci` | view | 8 | — | `define_genome()` | kept | View: one row per (term x locus), with locus_name joined from genome_meta. The place to answer which loci are causal for a trait; locus_name lives only here, so there is no id/name agreement invariant to keep in the base tables. |
| `chr_inheritance` | table | 5 | `chr_name`, `offspring_sex`, `line_name` (logical) | `define_genome()` | kept | Per-chromosome copy counts, keyed by offspring sex. One row per (chr_name, offspring_sex, line_name). Seeded default is a diploid autosome (from_parent_1 = from_parent_2 = 1). Non-default rules (sex chromosomes, organelles) are set via define_chromosome(). |
| `chr_recombination` | table | 4 | `chr_name`, `parent_sex`, `line_name` (logical) | `define_genome()` | kept | Per-chromosome recombination, keyed by producing-parent sex. One row per (chr_name, parent_sex, line_name). Seeded from define_genome()'s genome-wide recombines_M/recombines_F defaults. Non-default rules (Y, W, achiasmy) are set via define_chromosome(). |

### Founders (1)

| Table | Kind | Cols | Key | Created by | Archive | Description |
|---|---|---:|---|---|---|---|
| `founder_haplotypes` | table | 4 | `line_name`, `haplotype_id`, `locus_name` (logical) | `define_founder_haplotypes()` | kept | Pool of founder haplotypes in long format (one row per haplotype x locus) sampled by add_founders() to assign phased alleles. |

### Individuals (4)

| Table | Kind | Cols | Key | Created by | Archive | Description |
|---|---|---:|---|---|---|---|
| `ind_meta` | table | 6 | `id_ind` | `open_pop()` | per replicate | Individual-level metadata. One row per individual. Core columns are managed by the system; user-defined columns can be added via mutate_table() or the ... arguments of add_founders() and add_offspring(). |
| `ind_haplotype` | table | 7 | `id_ind`, `parent_origin`, `strand`, `locus_id` (logical) | `define_genome()` | reset only | Phased haplotypes in long format. One row per individual x haplotype x locus. parent_origin (1/2) and strand (1 for diploids) identify the copy; line_origin traces the allele's founding line. No DB primary key (dropped for insert speed); (id_ind, parent_origin, strand, locus_id) is unique by construction, guaranteed R-side. |
| `ind_genotype` | table | 4 | `id_ind`, `locus_id` | `define_genome()` | reset only | Genotype dosage cache in long format. One row per individual x locus. NOT auto-populated; filled on demand by add_dosage() from ind_haplotype. May be empty or partial. PRIMARY KEY (id_ind, locus_id). |
| `ind_crossover` | table | 6 | `id_crossover` | `define_genome()` | reset only | Crossover events in long format, one row per crossover drawn during meiosis. Created empty by define_genome(); populated only when add_offspring(store_crossovers = TRUE) (row writes land with the Stage-2 kernel). Absence of a row for a (id_ind, parent_origin, chr) means that gamete's chromosome did not recombine. |

### Genetic model (2)

| Table | Kind | Cols | Key | Created by | Archive | Description |
|---|---|---:|---|---|---|---|
| `trait_meta` | table | 5 | `id_trait` | `define_trait()` | once | Genetic component trait definitions. One row per trait. Genetic layer only — no observation-layer metadata. Populated by define_trait(). |
| `trait_var_comp` | table | 5 | `id_trait_var_comp` | `open_pop()` | once | Genetic variance component storage. One row per (effect_name, trait_name_1, trait_name_2); both (i,j) and (j,i) stored. Reserved effect_name values: 'gen_add', 'dominance' (future), 'epistasis' (future). |

### Observation model (5)

| Table | Kind | Cols | Key | Created by | Archive | Description |
|---|---|---:|---|---|---|---|
| `phenotype_meta` | table | 16 | `id_phenotype_meta` | `open_pop()` | once | Observed phenotype definitions. One row per phenotype. Manages the observation layer: type, mean, sex expression, and distributional parameters. Populated by define_phenotype(). |
| `phenotype_components` | table | 17 | `id_phenotype_comp` | `open_pop()` | once | Component definitions for composite phenotypes. One row per phenotype x component. Enables maternal effects, social genetic effects, and multi-contributor phenotypes. Populated by define_phenotype(..., components = ...). |
| `phenotype_var_comp` | table | 10 | `id_phenotype_var_comp` | `open_pop()` | once | Phenotype-level variance component storage. One row per (effect_name, phenotype pair, optional condition). Stores residual covariances (effect_name = 'residual') and named random effects (hys, litter, pen, etc.). Populated by define_phenotype(), define_residual_cov(), and define_effect_random(). |
| `phenotype_effects` | table | 12 | `phenotype_name`, `effect_name` | `define_trait()` | once | Fixed and random effect configurations for phenotype models. One row per phenotype x effect. Populated by define_effect_fixed_class(), define_effect_fixed_cov(), define_effect_random(). |
| `phenotype_random_effects` | table | 5 | `phenotype_name`, `effect_name`, `level` | `define_trait()` | per replicate | Sampled random effect levels. One row per phenotype x effect x level. Populated by add_phenotype() on first use; subsequent calls reuse the stored value for consistency. |

### Selection (1)

| Table | Kind | Cols | Key | Created by | Archive | Description |
|---|---|---:|---|---|---|---|
| `index_meta` | table | 5 | `id_index_name` | `define_trait()` | once | Selection index definitions. One row per index x trait. A special row with index_name = NULL holds the global economic weight per trait written by define_trait(). Named indices hold selection weights from define_index(). |

### Results (7)

| Table | Kind | Cols | Key | Created by | Archive | Description |
|---|---|---:|---|---|---|---|
| `ind_tbv` | table | 4 | `id_tbv` | `define_trait()` | per replicate | True breeding values (simulation ground truth). One row per individual x trait. Populated by add_phenotype() and add_tbv(). Values computed from genome effects in genome_effects. |
| `ind_tgv` | table | 5 | `id_tgv` | `define_trait()` | per replicate | True genetic values (simulation ground truth): the total genotypic value, split by declared model structure. One row per individual x trait x component. Populated by add_tgv(). The total is the derived view ind_tgv_total, never a stored row, so SUM(tgv_value) cannot double-count. |
| `ind_tgv_total` | view | 3 | — | `define_trait()` | kept | View: the total genetic value per individual x trait, summing every component of ind_tgv. Derived rather than stored so it can never disagree with its parts. |
| `ind_phenotype` | table | 5 | `id_phenotype` | `define_trait()` | per replicate | Phenotype records in long format. One row per individual x phenotype x record number. Populated by add_phenotype(). User-defined columns can be added via the ... argument of add_phenotype(). |
| `ind_ebv` | table | 8 | `id_ebv` | `define_trait()` | per replicate | Estimated breeding values from external BLUP or GBLUP analyses. One row per individual x trait x model x evaluation number. Populated by add_ebv(). |
| `ind_index` | table | 5 | `id_index` | `define_trait()` | per replicate | Computed selection index values. One row per individual x index x run. Multiple runs are distinguished by index_number. Populated by add_index(). |
| `ind_true_index` | table | 5 | `id_true_index` | `define_trait()` | per replicate | True selection index values computed from TBVs. One row per individual x index x weight type. Populated by add_tbv() when index_names is supplied. |

### System (1)

| Table | Kind | Cols | Key | Created by | Archive | Description |
|---|---|---:|---|---|---|---|
| `_schema_meta` | table | 6 | `id_schema_meta` | `open_pop()` | kept | System table storing table and column descriptions for all tidybreed database objects. |

## Dependencies

- **Imports:** cli, dplyr, tibble, rlang, duckdb, DBI, dbplyr, dqrng, MASS, Rcpp, stats
- **LinkingTo:** Rcpp, dqrng, BH, sitmo
- **Suggests:** testthat (≥ 3.2.0), withr, knitr, rmarkdown, covr
- **SystemRequirements:** C++17
- **License:** MIT
