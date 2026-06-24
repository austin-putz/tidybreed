# Scenario Management System — Design Plan

**Package:** `tidybreed`  
**Version bump:** `0.9.4` → `0.10.0`  
**Status:** Draft for team review

---

## What This Adds

A YAML-based scenario system that lets users:

1. Define breeding program parameters in `scenarios/*.yaml` files
2. Write one `main.R` simulation script that reads params from a YAML
3. Launch `main.R` for every YAML in `scenarios/` — locally via `Rscript`, or via `sbatch` on SLURM HPC clusters

**Design model: parameter file, not declarative automation.**  
The YAML is a _parameter store_. The user's `main.R` contains all simulation logic (mating loops, selection, custom columns, etc.) and reads parameters via `load_scenario()`. The package does not generate R code from YAML — the R script stays in control.

---

## Files to Create / Modify

| File | Action | Description |
|------|--------|-------------|
| `R/scenario.R` | New | `load_scenario`, `print.tidybreed_scenario`, `list_scenarios`, `initialize_from_scenario` |
| `R/run_scenarios.R` | New | `run_scenarios`, `use_scenario_template`, `use_main_template`, `.build_sbatch_script` |
| `inst/templates/scenario_template.yaml` | New | Fully commented example YAML |
| `inst/templates/main_template.R` | New | Template `main.R` showing canonical usage |
| `inst/templates/slurm_template.sh` | New | Standalone sbatch template for manual HPC submission |
| `R/tidybreed-package.R` | Modify | Add `@importFrom yaml read_yaml` and `@importFrom tools file_path_sans_ext` |
| `DESCRIPTION` | Modify | Add `yaml` to Imports, `parallel` to Suggests, bump version to `0.10.0` |
| `NEWS.md` | Modify | Add `0.10.0` entry |
| `NAMESPACE` | Regenerate | Via `devtools::document()` |

---

## YAML Schema

Seven top-level sections. Only `genome` and `founders` are required.

```yaml
# ============================================================
# tidybreed scenario: high-performance pig selection program
# ============================================================

scenario:
  name: "high_sel_pig"          # used for {scenario_name} interpolation below
  description: "High-intensity 2-trait selection with BLUP index"
  version: 1                    # user-defined integer for tracking

# ---- Reproducibility / Output ----------------------------------------
output:
  db_path: "results/{scenario_name}.duckdb"   # {scenario_name} is auto-filled
  results_dir: "results/{scenario_name}"
  seed: 42

# ---- Genome Architecture ---------------------------------------------
genome:
  pop_name: "PigSim"
  n_loci: 50000
  n_chr: 18
  chr_len_Mb: 120.0             # scalar (all chromosomes same) OR list of n_chr values
  n_haplotypes: 200
  min_allele_freq: 0.05         # use this pair -OR- fixed_allele_freq below (not both)
  max_allele_freq: 0.95
  # fixed_allele_freq: 0.5      # alternative: single fixed frequency for all loci

# ---- Founding Lines --------------------------------------------------
founders:
  - line_name: "A"
    n_males: 20
    n_females: 200
    gen: 0                      # any extra keys become ... args to add_founders()
    farm: "Iowa"
  - line_name: "B"
    n_males: 10
    n_females: 100
    gen: 0
    farm: "Kansas"

# ---- Traits ----------------------------------------------------------
traits:
  - trait_name: "ADG"
    description: "Average daily gain"
    units: "g/day"
    trait_type: "continuous"        # continuous | count | binary | categorical
    expressed_sex: "both"           # both | M | F
    target_add_mean: 850.0
    target_add_var: 2500.0          # option 1: explicit variances
    residual_var: 7500.0
    # h2: 0.25                      # option 2 (alternative): heritability
    # pheno_var: 10000.0            #   loader computes add_var = h2 * pheno_var
    n_qtl: 500
    qtl_method: "chromosome_even"   # random | even | chromosome_even
    effect_distribution: "normal"   # normal | gamma
    # min_value: 0                  # for count traits only
    # max_value: ~                  # ~ = null = no upper cap
    # prevalence: 0.3               # for binary traits only
    index_weight: 1.0
    economic_value: 0.05

  - trait_name: "BF"
    description: "Backfat thickness"
    units: "mm"
    trait_type: "continuous"
    expressed_sex: "both"
    target_add_mean: 12.0
    target_add_var: 4.0
    residual_var: 12.0
    n_qtl: 300
    qtl_method: "chromosome_even"
    effect_distribution: "normal"
    index_weight: -0.5
    economic_value: -0.10

# ---- Covariance Matrices ---------------------------------------------
# Omit entirely for single-trait simulations.
# effect_name must be "gen_add", "residual", or a named random effect.
covariances:
  - effect_name: "gen_add"
    traits: ["ADG", "BF"]
    matrix:                         # rows listed top-to-bottom, same order as traits
      - [2500.0, -500.0]
      - [-500.0,    4.0]

  - effect_name: "residual"
    traits: ["ADG", "BF"]
    matrix:
      - [7500.0, -200.0]
      - [-200.0,   12.0]

# ---- Fixed Effects ---------------------------------------------------
fixed_effects:
  - trait: "ADG"
    effect_name: "sex"
    source_column: "sex"            # column in ind_meta (or source_table)
    levels:
      M: 30.0
      F:  0.0

  - trait: "BF"
    effect_name: "sex"
    source_column: "sex"
    levels:
      M: -1.5
      F:  0.0

# ---- Random Effects --------------------------------------------------
random_effects:
  - trait: "ADG"
    effect_name: "herd"
    source_column: "farm"           # column in ind_meta defining groups
    variance: 500.0
    distribution: "normal"          # normal | gamma | uniform

# ---- User Simulation Parameters (not validated by package) -----------
# Access in main.R via:  p <- scenario$simulation
simulation:
  n_sires: 10
  n_dams: 100
  n_generations: 5
  selection_method: "index"         # user-defined strings; not parsed by package
  genotyping_rate: 0.5
  chip_name: "50K"
```

---

## Exported Functions

### `load_scenario(path)`

Reads a YAML file, validates all sections, interpolates `{scenario_name}` in output paths, converts `h2 + pheno_var` to canonical `target_add_var + residual_var`, and returns a `tidybreed_scenario` S3 object (a named list with class).

```r
scenario <- load_scenario("scenarios/high_sel_pig.yaml")
print(scenario)
# <tidybreed scenario: high_sel_pig>
# Description: High-intensity 2-trait selection with BLUP index
# DB path:     results/high_sel_pig.duckdb
# Genome:      50000 loci, 18 chr, 200 haplotypes
# Founders:    2 lines (460 total)
# Traits:      ADG, BF [h2: 0.25, 0.25]
# Simulation:  6 user params
```

**Validations performed:**
- File exists and parses as valid YAML
- Required `genome` keys all present (`pop_name`, `n_loci`, `n_chr`, `chr_len_Mb`, `n_haplotypes`)
- Allele freq mode: exactly one of `fixed_allele_freq` OR `(min_allele_freq + max_allele_freq)`
- `founders`: list with ≥1 entry, each with `line_name`, `n_males`, `n_females`
- Per trait: `trait_name`, `trait_type` (valid enum), `n_qtl` required; variance spec is mutually exclusive (`add_var + residual_var` vs `h2 + pheno_var`); `h2` must be in (0, 1)
- Covariance matrix dimensions match `length(traits)` and matrix is square
- Fixed/random effect entries have all required keys

---

### `list_scenarios(dir = "scenarios")`

Scans a directory for `*.yaml` / `*.yml` files and returns a data frame.

```r
list_scenarios("scenarios/")
#   name            path                          description
#   high_sel_pig    scenarios/high_sel_pig.yaml   High-intensity 2-trait selection...
#   low_sel_pig     scenarios/low_sel_pig.yaml    Low-intensity baseline program
```

Files that fail to parse get a `warning()` and `description = "[parse error]"` — one bad YAML does not block the rest.

---

### `initialize_from_scenario(scenario, db_path = NULL)`

Convenience wrapper that calls the tidybreed setup functions in the correct order from YAML params. Returns a fully-initialized `tidybreed_pop`.

**Order of operations:**
1. `initialize_genome()` from `scenario$genome`
2. `add_founders()` per entry in `scenario$founders` (extra YAML keys → `...` args)
3. `add_trait_simple()` per trait in `scenario$traits`
4. `add_effect_cov_matrix()` per entry in `scenario$covariances`
5. `add_effect_fixed_class()` per entry in `scenario$fixed_effects`
6. `add_effect_random()` per entry in `scenario$random_effects`

> **Note on multi-trait correlated QTL effects:** `add_trait_simple()` handles single-trait QTL effect sampling. For correlated multi-trait effects (using `set_qtl_effects_multi()`), call it manually in `main.R` after `initialize_from_scenario()`.

---

### `run_scenarios(script, scenarios_dir = "scenarios", method = c("local", "sbatch"), parallel = FALSE, dry_run = FALSE, slurm_opts = list(), ...)`

Launches `main.R` for every YAML in `scenarios_dir`. Returns a data frame of results.

```r
# Preview what would run (no execution)
run_scenarios("main.R", "scenarios/", method = "local", dry_run = TRUE)

# Run sequentially, locally
run_scenarios("main.R", "scenarios/", method = "local")

# Run in parallel locally (Unix/Mac only)
run_scenarios("main.R", "scenarios/", method = "local", parallel = TRUE)

# Submit to SLURM
run_scenarios("main.R", "scenarios/",
              method     = "sbatch",
              slurm_opts = list(mem = "16G", time = "4:00:00",
                                partition = "normal", account = "mylab"))
```

**`slurm_opts` recognized keys:**

| Key | Default | SLURM flag |
|-----|---------|------------|
| `mem` | `"8G"` | `--mem` |
| `time` | `"01:00:00"` | `--time` |
| `job_name` | scenario name | `--job-name` |
| `output` | `"logs/{scenario_name}_%j.out"` | `--output` |
| `error` | `"logs/{scenario_name}_%j.err"` | `--error` |
| `ntasks` | `1` | `--ntasks` |
| `partition` | (omitted if NULL) | `--partition` |
| `account` | (omitted if NULL) | `--account` |
| `modules` | (omitted if NULL) | `module load ...` lines |

`{scenario_name}` in output/error paths is replaced on the R side before submission; `%j` (SLURM job ID) is left as-is.

**dry_run for sbatch** prints the _full generated `.sh` script content_ — not just the `sbatch` command — so HPC users can inspect exactly what will be submitted.

---

### `use_scenario_template(path = "scenarios/default.yaml")`
### `use_main_template(path = "main.R")`

Scaffold helper functions (similar to `usethis::use_*`). Copy the template from `inst/templates/` to the specified path, create parent directories if needed, warn without overwriting if the file already exists.

```r
use_scenario_template("scenarios/program_A.yaml")
# Created 'scenarios/program_A.yaml'.

use_main_template("main.R")
# Created 'main.R'.
```

---

## Template `main.R` (canonical structure)

```r
#!/usr/bin/env Rscript
# Usage: Rscript main.R scenarios/program_A.yaml

library(tidybreed)

# 1. Parse scenario file from command-line argument
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) stop("Usage: Rscript main.R <scenario.yaml>")

scenario <- load_scenario(args[1])
print(scenario)

# 2. Set seed and create output directory
set.seed(scenario$output$seed)
if (!dir.exists(scenario$output$results_dir))
  dir.create(scenario$output$results_dir, recursive = TRUE)

# 3. Initialize population from YAML parameters
pop <- initialize_from_scenario(scenario)

# OR: call functions manually for full control
# pop <- initialize_genome(
#   pop_name     = scenario$genome$pop_name,
#   n_loci       = scenario$genome$n_loci,
#   n_chr        = scenario$genome$n_chr,
#   chr_len_Mb   = scenario$genome$chr_len_Mb,
#   n_haplotypes = scenario$genome$n_haplotypes,
#   db_path      = scenario$output$db_path
# )

# 4. Access user-defined simulation parameters
p <- scenario$simulation     # user-defined section; no package validation
n_sires <- p$n_sires
n_dams  <- p$n_dams

# 5. Simulation loop (user writes this section)
for (gen in seq_len(p$n_generations)) {

  # --- Phenotyping ---
  pop <- pop |>
    get_table("ind_meta") |>
    dplyr::filter(gen == !!gen) |>
    add_phenotype(c("ADG", "BF"))

  # --- Selection (example: top TBV) ---
  sires <- get_table(pop, "ind_tbv") |>
    dplyr::filter(trait_name == "ADG") |>
    dplyr::slice_max(tbv, n = n_sires) |>
    dplyr::pull(id_ind)

  # --- Mating ---
  matings <- ... # build mating tibble from selected parents

  pop <- pop |>
    add_offspring(matings)
}

# 6. Save / export results (population already persists in DuckDB)
results <- get_table(pop, "ind_tbv") |> dplyr::collect()
saveRDS(results, file.path(scenario$output$results_dir, "tbv_final.rds"))

close_pop(pop)
message("Done: ", scenario$scenario_name)
```

---

## Dependencies

| Package | Placement | Reason |
|---------|-----------|--------|
| `yaml` | `Imports` | Required at runtime for `load_scenario()` |
| `parallel` | `Suggests` | Optional; used only when `run_scenarios(parallel = TRUE)`; warn + fall back to serial on Windows |
| `tools` | `Imports` (already base R) | `file_path_sans_ext()` for scenario name derivation |

---

## Questions for Team Review

1. **YAML schema scope** — The `simulation:` section is completely free-form and not validated. Is this the right call, or should we define a fixed set of required keys (e.g., `n_generations`) that the package enforces?

2. **`initialize_from_scenario()` scope** — For multi-trait correlated QTL effects, users need to call `set_qtl_effects_multi()` themselves in `main.R`. Should we attempt to auto-detect this case and call it automatically, or leave it to the user?

3. **Scenario versioning** — Should `load_scenario()` validate a `schema_version:` field so we can evolve the YAML format without silently breaking old files?

4. **HPC: job arrays** — SLURM job arrays (`sbatch --array=1-10 job.sh`) are more efficient than N separate submissions for large scenario sets. Should `run_scenarios()` support this as an option, or is N-separate-jobs sufficient for the initial release?

5. **Error handling in launcher** — If one scenario's `Rscript` call fails, should `run_scenarios()` stop immediately, continue and collect errors, or re-try? Current plan: continue and report exit codes in the returned data frame.

6. **Naming** — `use_scenario_template()` / `use_main_template()` follow the `usethis` naming convention. Is this the right naming, or would `scaffold_scenario()` / `scaffold_main()` feel more natural in the tidybreed context?

---

## Verification Checklist

- [ ] `use_scenario_template("scenarios/test.yaml")` creates the file
- [ ] `load_scenario("scenarios/test.yaml")` returns a `tidybreed_scenario`; `print()` shows clean summary
- [ ] `list_scenarios("scenarios/")` returns a data frame with the test scenario
- [ ] `initialize_from_scenario(scenario)` creates a valid `tidybreed_pop`; `get_table(pop, "ind_meta") |> collect()` has rows
- [ ] `run_scenarios("main.R", "scenarios/", method = "local", dry_run = TRUE)` prints the `Rscript` command without executing
- [ ] `run_scenarios("main.R", "scenarios/", method = "sbatch", dry_run = TRUE)` prints the full sbatch script
- [ ] `run_scenarios("main.R", "scenarios/", method = "local")` executes and returns exit code 0
- [ ] `R CMD check` passes with no new warnings
