# Plan: External Tool Run Directory Structure

## Problem

External software wrappers (BLUPF90, PLINK, etc.) currently create many
individual working directories directly in the user's working directory, which
quickly becomes unmanageable — especially across many scenarios and replicates.

Users need access to individual runs to verify outputs and need full isolation
between scenarios and replicates launched as separate processes (HPC, background jobs).

## Goals

- All simulation output nested under one clearly named output folder
- One subfolder per scenario (possibly from YAML) inside the output folder
- Per-tool subdirectories inside each scenario folder
- Per-run timestamped subdirectories inside each tool folder — deleted when done
- No folder name collisions across parallel runs
- All stable settings controlled via global options

---

## Folder Layout

### With scenario name set

```text
<base_dir>/
  tidybreed_output/         ← Layer 2: output folder (tidybreed.output option)
    scenario_A/             ← Layer 3: scenario name (tidybreed.scenario option)
      blupf90/              ← Layer 4: tool dir (tidybreed.tools option)
        20260608_143022_a3f9b1/    ← Layer 5: runtime temp dir, deleted when done
        20260608_150311_c7d2e8/
      plink/
        20260608_143055_9f1a2c/
      sim.duckdb
    scenario_B/
      blupf90/
      plink/
      sim.duckdb
```

### Without scenario name (`tidybreed.scenario = NULL`)

When no scenario name is provided, `open_pop()` auto-generates a `YYYYMMDD_HHMMSS`
folder at layer 3 so the DuckDB and tool dirs are always isolated — even across
multiple sequential runs without a scenario name.

```text
<base_dir>/
  tidybreed_output/
    20260608_143022/        ← auto-generated timestamp folder
      blupf90/
      plink/
      sim.duckdb
    20260608_150311/        ← second run, isolated automatically
      blupf90/
      plink/
      sim.duckdb
```

### HPC replicates (replicate encoded in scenario name)

For HPC array jobs, encode the replicate number in the scenario name. Each job
reads its own YAML file containing a unique `scenario_name`. No separate rep
layer is needed — the scenario folder provides full isolation.

```text
tidybreed_output/
  scenario_A_rep_01/
    blupf90/
    plink/
    sim.duckdb
  scenario_A_rep_02/
    blupf90/
    plink/
    sim.duckdb
```

### Layer descriptions

**Layer 1 — Base directory** (`tidybreed.base_dir`): physical root, defaults to
`getwd()`. Constant across all runs.

**Layer 2 — Output directory** (`tidybreed.output`): collects all tidybreed
output under one named folder. Defaults to `"tidybreed_output"` to avoid
colliding with existing user folders. Users can override if desired.

**Layer 3 — Scenario directory** (`tidybreed.scenario`): one folder per scenario
or replicate. Set from a YAML `scenario_name` field or directly via `options()`.
When `NULL`, `open_pop()` auto-generates a `YYYYMMDD_HHMMSS` folder name so runs
are always isolated.

**Layer 4 — Tool directories** (`tidybreed.tools`): one subfolder per external
tool (`blupf90/`, `plink/`, etc.). Declared once via global option. Created
eagerly by `open_pop()`.

**Layer 5 — Run directories**: timestamped + 6-char random hex, created at the
start of each individual tool call. Deleted on completion unless `keep = TRUE`
or `keep = "on_error"`. Format: `YYYYMMDD_HHMMSS_rrrrrr`. Six hex chars
(`16^6` ≈ 16M combinations per second) provides near-zero collision risk even
when 50+ jobs launch simultaneously.

`pop$run_dirs$base` always points to the layer-3 folder (named or auto-generated),
so `.create_run_dir()` and all wrappers work identically regardless of which
options are set.

---

## Scenarios Folder Convention

Users are encouraged to keep a `scenarios/` folder in their project root
alongside `tidybreed_output/`:

```text
<project_root>/
  scenarios/
    scenario_A.yaml
    scenario_B.yaml
  tidybreed_output/
    scenario_A/
    scenario_B/
```

The `scenarios/` folder is user-managed — tidybreed does not create it.
The `scenario_name` field inside each YAML file is what gets passed to
`options(tidybreed.scenario = ...)`.

For HPC replicates, use one YAML file per replicate with a distinct
`scenario_name` (e.g., `scenario_A_rep_01`):

```text
<project_root>/
  scenarios/
    scenario_A.yaml
    scenario_B.yaml
  tidybreed_output/
    scenario_A_rep_01/
    scenario_A_rep_02/
    scenario_B_rep_01/
    scenario_B_rep_02/
```

---

## Global Options

| Option | Default | Description |
|--------|---------|-------------|
| `tidybreed.base_dir` | `getwd()` | Root location — layer 1 |
| `tidybreed.output` | `"tidybreed_output"` | Output folder name — layer 2 |
| `tidybreed.scenario` | `NULL` | Scenario name; `NULL` triggers auto-timestamp folder — layer 3 |
| `tidybreed.tools` | `NULL` | Tool subdir names (e.g. `c("blupf90", "plink")`) — layer 4 |
| `tidybreed.db_name` | `"sim.duckdb"` | DuckDB filename placed inside the layer-3 folder |

---

## Functions

### `open_pop(pop_name, tools, base_dir, output_dir, scenario_dir, db_name, clean)`

Public. The single entry point to start a simulation. Reads all global options,
creates the full folder structure (layers 2–4), creates the DuckDB file, creates
all empty core tables, and returns a `tidybreed_pop`. Replaces the old
`initialize_genome()` for database and folder setup — genome structure is then
defined separately with `define_genome()`.

| Argument       | Type      | Default                            | Description |
|----------------|-----------|------------------------------------|-------------|
| `pop_name`     | character | `NULL`                             | Optional label stored on the pop object |
| `base_dir`     | character | `getOption("tidybreed.base_dir")`  | Root location (layer 1) |
| `output_dir`   | character | `getOption("tidybreed.output")`    | Output folder name (layer 2) |
| `scenario_dir` | character | `getOption("tidybreed.scenario")`  | Scenario name; `NULL` auto-generates `YYYYMMDD_HHMMSS` |
| `tools`        | character | `getOption("tidybreed.tools")`     | Tool subdirs to create (layer 4) |
| `db_name`      | character | `getOption("tidybreed.db_name")`   | DuckDB filename |
| `clean`        | logical   | `TRUE`                             | If `TRUE`, delete existing layer-5 run subdirs inside each tool dir |

**What `open_pop()` does internally:**
1. Resolves the layer-3 folder: uses `scenario` if provided, otherwise generates `YYYYMMDD_HHMMSS`
2. Creates layers 2–4 (`tidybreed_output/`, scenario folder, tool subdirs)
3. Creates the DuckDB at `layer3/db_name`
4. Creates all empty core tables (genome, individual, trait tables)
5. Stores resolved paths in `pop$run_dirs`
6. Returns `tidybreed_pop`

**`clean` semantics**: deletes timestamped layer-5 run subdirs inside each tool
dir. Does not delete the tool dirs or the DuckDB. Idempotent — existing dirs are
reused without error. Also, delete the `db_name` from layer 3 
folder before it's recreated. 

Example — single scenario:

```r
options(
  tidybreed.base_dir = getwd(),
  tidybreed.output   = "tidybreed_output",
  tidybreed.scenario = "scenario_A",
  tidybreed.tools    = c("blupf90", "plink"),
  tidybreed.db_name  = "sim.duckdb"
)

pop <- open_pop() |>
  define_genome(n_loci = 50000, n_chr = 18, chr_len_Mb = 100) |>
  define_founder_haplotypes(...) |>
  add_founders(n_males = 50, n_females = 500, line_name = "A")
# creates: tidybreed_output/scenario_A/blupf90/
#          tidybreed_output/scenario_A/plink/
#          tidybreed_output/scenario_A/sim.duckdb
```

Example — no scenario name (auto-timestamp):

```r
options(
  tidybreed.tools = c("blupf90", "plink")
  # tidybreed.scenario left NULL
)

pop <- open_pop()
# creates: tidybreed_output/20260608_143022/blupf90/
#          tidybreed_output/20260608_143022/plink/
#          tidybreed_output/20260608_143022/sim.duckdb
```

Example — HPC array job:

```r
config <- read_sim_config(Sys.getenv("SIM_YAML"))
# sets options(tidybreed.scenario = "scenario_A_rep_01")

pop <- open_pop() |>
  define_genome(n_loci = config$n_loci, ...) |>
  ...
# creates: tidybreed_output/scenario_A_rep_01/blupf90/
#          tidybreed_output/scenario_A_rep_01/sim.duckdb
```

Returns `pop$run_dirs`:

```r
pop$run_dirs
# $base      "/proj/tidybreed_output/scenario_A"
# $blupf90   "/proj/tidybreed_output/scenario_A/blupf90"
# $plink     "/proj/tidybreed_output/scenario_A/plink"
```

---

### `restore_pop(path)`

Opens an existing DuckDB and reconnects the pop object. Also recreates tool dirs
from the current global options (`tidybreed.tools`) since tool dirs are
session-level and not stored in the database. Users call `options()` to set
tools before calling `restore_pop()`.

---

### `.create_run_dir(pop, tool, keep)` — internal

Called at the start of every external tool wrapper. Creates a timestamped
layer-5 subfolder inside the named tool dir and registers cleanup via `on.exit()`
in the caller's frame.

| Argument | Type          | Default      | Description |
|----------|---------------|--------------|-------------|
| `pop`    | tidybreed_pop | —            | Source of `pop$run_dirs` |
| `tool`   | character     | —            | Which tool dir to use (must exist in `pop$run_dirs`) |
| `keep`   | character     | `"on_error"` | `FALSE` = always delete, `TRUE` = always keep, `"on_error"` = keep only if call errored |

Returns the path to the new run directory. Errors with a clear message if
`pop$run_dirs` is `NULL` or if `tool` is not present in `pop$run_dirs`.

---

## Storage on `pop`

`open_pop()` attaches resolved paths to `pop$run_dirs` as a named list.
Session-level only — not stored in the database. Users who restart a session
call `restore_pop()` which also recreates tool dirs.

---

## Wrapper Integration Pattern

```r
run_plink <- function(pop, ..., keep = "on_error") {
  run_dir <- .create_run_dir(pop, tool = "plink", keep = keep)
  # write input files to run_dir
  # call plink via system2()
  # read output files from run_dir
  # return pop  (cleanup handled by on.exit registered in .create_run_dir)
}
```

---

## Simulation Workflow (full example)

```r
# 1. Set global options once at top of script
options(
  tidybreed.base_dir = getwd(),
  tidybreed.output   = "tidybreed_output",
  tidybreed.scenario = "baseline_v1",
  tidybreed.tools    = c("blupf90", "plink"),
  tidybreed.db_name  = "sim.duckdb"
)

# 2. Open population — creates folders + DuckDB + empty tables
pop <- open_pop(pop_name = "baseline") |>

# 3. Define genome architecture
  define_genome(n_loci = 50000, n_chr = 18, chr_len_Mb = 100) |>

# 4. Build founder haplotype pool
  define_founder_haplotypes(...) |>

# 5. Define traits
  define_trait("ADG", target_add_var = 0.25) |>
  define_additive_effects("ADG") |>

# 6. Add founders — simulation output begins here
  add_founders(n_males = 50, n_females = 500, line_name = "A")
```

---

## Open Questions

1. Should `pop$run_dirs` be a named list or named character vector?
   Named list allows attaching metadata later; named character vector is simpler.
2. Should tools not in `pop$run_dirs` be created on-demand by `.create_run_dir()`
   or cause an error? Current thinking: error — forces upfront declaration via
   `tidybreed.tools` option.

