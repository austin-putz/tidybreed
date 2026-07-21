# Open a new breeding population

The single entry point to start a simulation. Creates the DuckDB
database, all core metadata tables, and the per-scenario folder
structure (layers 2–4). Genome tables are added separately by
[`define_genome()`](https://austin-putz.github.io/tidybreed/reference/define_genome.md).

Reads global options for default argument values so simulations can be
parameterized via [`options()`](https://rdrr.io/r/base/options.html) or
a YAML config without changing the script body.

### Folder layout

    <base_dir>/
      <output_dir>/                   # layer 2 (e.g. "tidybreed_output")
        <scenario_dir>/               # layer 3 (e.g. "scenario_A" or auto timestamp)
          <db_name>                   # the DuckDB file
          blupf90/                    # layer 4 (if "blupf90" in tools)
          plink/                      # layer 4 (if "plink" in tools)

When `db_name = ":memory:"` or an absolute/UNC path, folder creation is
skipped entirely (including layer-4 `tools` dirs) and `pop$run_dirs` is
`character(0)`.

## Usage

``` r
open_pop(
  pop_name = getOption("tidybreed.pop_name", "sim"),
  base_dir = getOption("tidybreed.base_dir", getwd()),
  output_dir = getOption("tidybreed.output", "tidybreed_output"),
  scenario_dir = getOption("tidybreed.scenario", NULL),
  tools = getOption("tidybreed.tools", NULL),
  db_name = getOption("tidybreed.db_name", "sim.duckdb"),
  clean = TRUE
)
```

## Arguments

- pop_name:

  Character scalar. Optional label stored on the pop object. Default:
  `getOption("tidybreed.pop_name", "sim")`.

- base_dir:

  Character scalar. Root location (layer 1). Default:
  `getOption("tidybreed.base_dir", getwd())`.

- output_dir:

  Character scalar. Output folder name (layer 2). Default:
  `getOption("tidybreed.output", "tidybreed_output")`.

- scenario_dir:

  Character scalar or `NULL`. Scenario name (layer 3). When `NULL`
  (default), a `YYYYMMDD_HHMMSS` timestamp folder is auto-generated so
  runs are always isolated. Default:
  `getOption("tidybreed.scenario", NULL)`.

- tools:

  Character vector or `NULL`. Tool subdirectory names to create at layer
  4 (e.g. `c("blupf90", "plink")`). These become the keys of
  `pop$run_dirs`. Default: `getOption("tidybreed.tools", NULL)`. Ignored
  (no tool dirs are created; `pop$run_dirs` is `character(0)`) whenever
  `db_name` is `":memory:"` or an absolute/UNC path — see `db_name`
  below.

- db_name:

  Character scalar. DuckDB filename placed inside the layer-3 folder.
  Use `":memory:"` for an in-memory database, or supply an absolute/UNC
  path to use that exact file location directly — both skip all
  layer-2/3/4 folder creation entirely (including any `tools` dirs;
  `pop$run_dirs` will be `character(0)` even if `tools` is non-`NULL`).
  Default: `getOption("tidybreed.db_name", "sim.duckdb")`.

- clean:

  Logical. If `TRUE` (default), deletes the existing database file and
  all timestamped layer-5 run subdirectories inside each tool dir before
  recreating. Does **not** delete the tool dirs or the scenario dir.
  Silently ignored for in-memory databases.

## Value

A `tidybreed_pop` object with all core tables created and `pop$run_dirs`
populated (empty `character(0)` for in-memory databases).

## See also

[`define_genome()`](https://austin-putz.github.io/tidybreed/reference/define_genome.md)
to add genome tables,
[`restore_pop()`](https://austin-putz.github.io/tidybreed/reference/restore_pop.md)
to reconnect to an existing database.

## Examples

``` r
if (FALSE) { # \dontrun{
# In-memory (common in tests and development)
pop <- open_pop(pop_name = "test", db_name = ":memory:") |>
  define_genome(n_loci = 1000, n_chr = 10, chr_len_Mb = 100)

# File-based with folder structure and tool dirs
options(
  tidybreed.output   = "tidybreed_output",
  tidybreed.scenario = "scenario_A",
  tidybreed.tools    = c("blupf90", "plink"),
  tidybreed.db_name  = "sim.duckdb"
)
pop <- open_pop() |>
  define_genome(n_loci = 50000, n_chr = 18, chr_len_Mb = 100)
# creates: tidybreed_output/scenario_A/blupf90/
#          tidybreed_output/scenario_A/plink/
#          tidybreed_output/scenario_A/sim.duckdb

# HPC array job: scenario encoded in option from YAML
options(tidybreed.scenario = paste0("scenario_A_rep_", Sys.getenv("SLURM_ARRAY_TASK_ID")))
pop <- open_pop()
} # }
```
