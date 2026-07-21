# Archive a completed replicate and reset the working database

Stamps per-replicate tables with a `replicate` column and appends them
to a separate archive DuckDB file, copies static metadata tables to the
archive on the **first call only** (no `replicate` column), then deletes
the working-DB rows so the next replicate starts clean. This is the only
function that writes to an archive DB; the working DB stays lean
throughout the loop.

## Usage

``` r
archive_replicate(
  pop,
  replicate = getOption("tidybreed.replicate"),
  archive_path = NULL,
  store_and_reset = c("ind_meta", "ind_phenotype", "ind_tbv", "ind_ebv", "ind_index",
    "ind_true_index"),
  store_once = c("genome_meta", "genome_effects", "trait_meta", "trait_effects",
    "trait_var_comp", "phenotype_meta", "phenotype_components", "phenotype_var_comp",
    "index_meta"),
  reset_only = c("ind_haplotype", "ind_genotype", "ind_crossover")
)
```

## Arguments

- pop:

  A `tidybreed_pop` object with an open DuckDB connection.

- replicate:

  Integer scalar. Replicate number stamped on every archived row in
  `store_and_reset` tables. Defaults to
  `getOption("tidybreed.replicate")` (starts at `1L`). Always increments
  the option by 1 after a successful archive so the next loop iteration
  uses the next number automatically.

- archive_path:

  Character or `NULL`. Full path to the archive DuckDB file. When `NULL`
  (default) the path is resolved from global options — see Details.
  Passing `NULL` and leaving both archive options unset means only the
  reset phases run (no archive is written).

- store_and_reset:

  Character vector. Tables to archive with a `replicate` column and then
  DELETE from the working DB after archiving succeeds.

- store_once:

  Character vector. Tables copied to the archive on the first call only
  (detected by whether the table already exists in the archive schema).
  Working-DB rows are left intact. Use for static configuration that
  does not change between replicates (trait definitions, genome
  structure, variance components, index weights).

- reset_only:

  Character vector. Tables whose rows are deleted from the working DB
  without being archived. Intended for large genomic tables
  (`ind_haplotype`, `ind_genotype`, `ind_crossover`) that are
  unnecessary for most downstream analyses.

## Value

`pop` invisibly (pipeline-compatible).

## Details

**Archive path resolution** — first non-`NULL` result wins:

1.  Explicit `archive_path` argument (full file path)

2.  `file.path(tidybreed.archive_path, tidybreed.db_name_archive)` when
    both options are set

3.  `file.path(dirname(pop$db_path), tidybreed.db_name_archive)` when
    only `tidybreed.db_name_archive` is set (archive lands next to the
    working DB)

4.  `tidybreed.db_name_archive = NULL` → `archive_path` resolves to
    `NULL` → skip all archive writes; only the reset phases run

**Individual-ID invariant**: `id_ind` is **not** globally unique across
replicates, even in a single-process sequential loop.
`archive_replicate()` deletes `ind_meta` (a `store_and_reset` table)
after archiving, so
[`add_founders()`](https://austin-putz.github.io/tidybreed/reference/add_founders.md)
and
[`add_offspring()`](https://austin-putz.github.io/tidybreed/reference/add_offspring.md)
see an empty table on the next replicate and restart numbering from
`{line_name}_1`. Every replicate (and every HPC parallel job, which
independently restarts numbering too) reuses the same `id_ind` values —
always use `(replicate, id_ind)` as the composite key in archive
queries.

**HPC note**: concurrent writes from multiple jobs to one shared archive
file are **not supported**. Use per-job archive files
(`scenario_rep_007.duckdb`) and merge in a post-processing step.

## Usage patterns

**In-process loop**:

    options(tidybreed.replicate = 1L)
    options(tidybreed.db_name_archive = "scenario_1_all_reps.duckdb")

    pop <- open_pop(...) |> define_genome(...) |>
      define_trait(...) |> define_phenotype(...)

    for (i in seq_len(100)) {
      pop <- pop |>
        get_table("founder_haplotypes") |>
        add_founders(n_males = 10, n_females = 90, line_name = "A") |>
        get_table("ind_meta") |>
        add_phenotype("ADG") |>
        archive_replicate()
    }

    # post-hoc analysis. When only 'tidybreed.db_name_archive' is set (as
    # above), the archive file is written next to the working DB, which may
    # itself be nested under open_pop()'s layer-2/3 folders (see ?open_pop) —
    # so resolve the path from the pop object rather than assuming a bare
    # filename in the current working directory.
    archive_path <- file.path(dirname(pop$db_path), "scenario_1_all_reps.duckdb")
    archive <- restore_pop(archive_path)
    get_table(archive, "ind_phenotype") |>
      dplyr::filter(replicate == 5L) |>
      dplyr::collect()

**HPC / SLURM array** (one job per replicate):

    rep_id <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID"))
    options(tidybreed.replicate = rep_id)
    options(tidybreed.db_name_archive = paste0("scenario_1_rep_", rep_id, ".duckdb"))
    pop <- open_pop(...) |> ... |> archive_replicate()
