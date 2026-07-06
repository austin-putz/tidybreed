#' Open a new breeding population
#'
#' @description
#' The single entry point to start a simulation. Creates the DuckDB database,
#' all core metadata tables, and the per-scenario folder structure (layers 2–4).
#' Genome tables are added separately by [define_genome()].
#'
#' Reads global options for default argument values so simulations can be
#' parameterized via `options()` or a YAML config without changing the script
#' body.
#'
#' ## Folder layout
#'
#' ```
#' <base_dir>/
#'   <output_dir>/                   # layer 2 (e.g. "tidybreed_output")
#'     <scenario_dir>/               # layer 3 (e.g. "scenario_A" or auto timestamp)
#'       <db_name>                   # the DuckDB file
#'       blupf90/                    # layer 4 (if "blupf90" in tools)
#'       plink/                      # layer 4 (if "plink" in tools)
#' ```
#'
#' When `db_name = ":memory:"` or an absolute/UNC path, folder creation is
#' skipped entirely (including layer-4 `tools` dirs) and `pop$run_dirs` is
#' `character(0)`.
#'
#' @param pop_name Character scalar. Optional label stored on the pop object.
#'   Default: `getOption("tidybreed.pop_name", "sim")`.
#' @param base_dir Character scalar. Root location (layer 1). Default:
#'   `getOption("tidybreed.base_dir", getwd())`.
#' @param output_dir Character scalar. Output folder name (layer 2). Default:
#'   `getOption("tidybreed.output", "tidybreed_output")`.
#' @param scenario_dir Character scalar or `NULL`. Scenario name (layer 3).
#'   When `NULL` (default), a `YYYYMMDD_HHMMSS` timestamp folder is
#'   auto-generated so runs are always isolated.
#'   Default: `getOption("tidybreed.scenario", NULL)`.
#' @param tools Character vector or `NULL`. Tool subdirectory names to create
#'   at layer 4 (e.g. `c("blupf90", "plink")`). These become the keys of
#'   `pop$run_dirs`. Default: `getOption("tidybreed.tools", NULL)`. Ignored
#'   (no tool dirs are created; `pop$run_dirs` is `character(0)`) whenever
#'   `db_name` is `":memory:"` or an absolute/UNC path — see `db_name` below.
#' @param db_name Character scalar. DuckDB filename placed inside the layer-3
#'   folder. Use `":memory:"` for an in-memory database, or supply an
#'   absolute/UNC path to use that exact file location directly — both skip
#'   all layer-2/3/4 folder creation entirely (including any `tools` dirs;
#'   `pop$run_dirs` will be `character(0)` even if `tools` is non-`NULL`).
#'   Default: `getOption("tidybreed.db_name", "sim.duckdb")`.
#' @param clean Logical. If `TRUE` (default), deletes the existing database
#'   file and all timestamped layer-5 run subdirectories inside each tool dir
#'   before recreating. Does **not** delete the tool dirs or the scenario dir.
#'   Silently ignored for in-memory databases.
#'
#' @return A `tidybreed_pop` object with all core tables created and
#'   `pop$run_dirs` populated (empty `character(0)` for in-memory databases).
#'
#' @seealso [define_genome()] to add genome tables, [restore_pop()] to
#'   reconnect to an existing database.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # In-memory (common in tests and development)
#' pop <- open_pop(pop_name = "test", db_name = ":memory:") |>
#'   define_genome(n_loci = 1000, n_chr = 10, chr_len_Mb = 100)
#'
#' # File-based with folder structure and tool dirs
#' options(
#'   tidybreed.output   = "tidybreed_output",
#'   tidybreed.scenario = "scenario_A",
#'   tidybreed.tools    = c("blupf90", "plink"),
#'   tidybreed.db_name  = "sim.duckdb"
#' )
#' pop <- open_pop() |>
#'   define_genome(n_loci = 50000, n_chr = 18, chr_len_Mb = 100)
#' # creates: tidybreed_output/scenario_A/blupf90/
#' #          tidybreed_output/scenario_A/plink/
#' #          tidybreed_output/scenario_A/sim.duckdb
#'
#' # HPC array job: scenario encoded in option from YAML
#' options(tidybreed.scenario = paste0("scenario_A_rep_", Sys.getenv("SLURM_ARRAY_TASK_ID")))
#' pop <- open_pop()
#' }
open_pop <- function(pop_name     = getOption("tidybreed.pop_name",  "sim"),
                     base_dir     = getOption("tidybreed.base_dir",   getwd()),
                     output_dir   = getOption("tidybreed.output",     "tidybreed_output"),
                     scenario_dir = getOption("tidybreed.scenario",   NULL),
                     tools        = getOption("tidybreed.tools",      NULL),
                     db_name      = getOption("tidybreed.db_name",    "sim.duckdb"),
                     clean        = TRUE) {

  stopifnot(is.character(pop_name), length(pop_name) == 1)
  stopifnot(is.logical(clean), length(clean) == 1)

  # Detect absolute path or :memory: — bypass folder creation
  in_memory <- identical(db_name, ":memory:")
  is_abs_path <- !in_memory && (
    grepl("^/", db_name) ||                  # Unix absolute path
    grepl("^[A-Za-z]:[/\\\\]", db_name) ||  # Windows absolute path
    grepl("^\\\\\\\\", db_name)              # UNC path
  )

  if (in_memory || is_abs_path) {
    # No folder structure — use db_name directly as db_path
    db_path  <- db_name
    run_dirs <- character(0)
    # Still apply clean semantics: delete the existing file if clean = TRUE
    if (!in_memory && clean && file.exists(db_path))
      file.remove(db_path)
  } else {
    # Resolve layer-3 (scenario) directory
    if (is.null(scenario_dir))
      scenario_dir <- format(Sys.time(), "%Y%m%d_%H%M%S")

    layer3 <- file.path(base_dir, output_dir, scenario_dir)

    if (clean) {
      # Delete existing DB file (will be recreated)
      db_full <- file.path(layer3, db_name)
      if (file.exists(db_full)) file.remove(db_full)
      # Delete layer-5 run subdirs inside each declared tool dir
      if (!is.null(tools)) {
        for (tool in tools) {
          tool_dir <- file.path(layer3, tool)
          if (dir.exists(tool_dir)) {
            run_subdirs <- list.dirs(tool_dir, recursive = FALSE, full.names = TRUE)
            if (length(run_subdirs) > 0)
              unlink(run_subdirs, recursive = TRUE)
          }
        }
      }
    }

    # Create layer 2 + layer 3
    dir.create(layer3, recursive = TRUE, showWarnings = FALSE)

    # Create layer-4 tool dirs and build run_dirs vector
    if (!is.null(tools) && length(tools) > 0) {
      tool_paths <- vapply(tools, function(t) {
        p <- file.path(layer3, t)
        dir.create(p, recursive = TRUE, showWarnings = FALSE)
        p
      }, character(1))
      run_dirs <- c(base = layer3, setNames(tool_paths, tools))
    } else {
      run_dirs <- c(base = layer3)
    }

    db_path <- file.path(layer3, db_name)
  }

  # Create DuckDB connection
  drv     <- duckdb::duckdb()
  db_conn <- DBI::dbConnect(drv, dbdir = db_path)

  # Create all core tables (all except genome tables)
  .create_core_tables(db_conn)

  # Register schema descriptions for core tables
  register_schema_meta(db_conn, .core_layer_descriptions())

  tables_created <- c(
    "_schema_meta", "ind_meta", "trait_var_comp", "genome_effects",
    "phenotype_meta", "phenotype_components", "phenotype_var_comp"
  )

  pop <- new_tidybreed_pop(
    db_conn  = db_conn,
    pop_name = pop_name,
    db_path  = db_path,
    tables   = tables_created,
    run_dirs = run_dirs
  )

  # Create remaining trait/phenotype/index tables (idempotent)
  pop <- ensure_trait_tables(pop)

  validate_tidybreed_pop(pop)

  if (in_memory) {
    message("Opened in-memory population '", pop_name, "'")
  } else if (!is_abs_path) {
    message("Opened population '", pop_name, "' at: ", db_path)
  } else {
    message("Opened population '", pop_name, "' at: ", db_path)
  }

  pop
}


#' Create a timestamped layer-5 run directory inside a tool dir
#'
#' @description
#' Called at the start of every external tool wrapper (e.g. `run_blupf90()`,
#' `run_plink()`). Creates a unique subdirectory inside `pop$run_dirs[[tool]]`
#' for the current invocation and returns the path.
#'
#' The tool must have been declared via the `tools` argument of [open_pop()].
#' If `tool` is not in `pop$run_dirs`, the function errors with a clear message
#' listing the registered tools.
#'
#' For v0.40.0 the `keep` argument is accepted for API completeness but cleanup
#' is not yet automatically registered. The caller is responsible for deleting
#' the run directory on success if `keep = FALSE` or `keep = "on_error"`.
#'
#' @param pop A `tidybreed_pop` object with `run_dirs` populated.
#' @param tool Character scalar. Which tool dir to use (must be a key in
#'   `pop$run_dirs`).
#' @param keep `FALSE`, `TRUE`, or `"on_error"`. Cleanup semantics (not yet
#'   enforced automatically; reserved for a future release).
#'
#' @return Path to the newly created run directory (invisibly).
#' @keywords internal
.create_run_dir <- function(pop, tool, keep = "on_error") {

  if (!inherits(pop, "tidybreed_pop"))
    stop("pop must be a tidybreed_pop object", call. = FALSE)

  registered <- names(pop$run_dirs)
  if (!tool %in% registered) {
    stop(
      "Tool '", tool, "' is not registered in pop$run_dirs. ",
      "Declare it via options(tidybreed.tools = c(\"", tool, "\", ...)) ",
      "before calling open_pop(). ",
      "Registered: ", paste(registered, collapse = ", "),
      call. = FALSE
    )
  }

  stamp   <- format(Sys.time(), "%Y%m%d_%H%M%S_")
  suffix  <- paste0(sample(c(letters, 0:9), 6, replace = TRUE), collapse = "")
  run_dir <- file.path(pop$run_dirs[[tool]], paste0(stamp, suffix))
  dir.create(run_dir, recursive = TRUE)

  invisible(run_dir)
}


#' Create all core (non-genome) database tables
#'
#' @description
#' Internal helper called by [open_pop()] (and the deprecated
#' `initialize_genome()` wrapper). Creates the `_schema_meta` system table and
#' all non-genome metadata tables. Genome tables (`genome_meta`,
#' `ind_haplotype`, `ind_genotype`, `chr_meta`) are created by [define_genome()].
#'
#' @param db_conn An active DuckDB connection.
#' @keywords internal
.create_core_tables <- function(db_conn) {

  # System metadata table — must exist before any register_schema_meta() calls
  DBI::dbExecute(db_conn, "
    CREATE TABLE _schema_meta (
      id_schema_meta INTEGER PRIMARY KEY,
      object_type    VARCHAR NOT NULL,
      table_name     VARCHAR NOT NULL,
      column_name    VARCHAR,
      description    VARCHAR NOT NULL,
      notes          VARCHAR
    )
  ")

  DBI::dbExecute(db_conn, "
    CREATE TABLE ind_meta (
      id_ind      VARCHAR PRIMARY KEY,
      id_parent_1 VARCHAR,
      id_parent_2 VARCHAR,
      line_name   VARCHAR,
      sex         VARCHAR,
      ploidy      UTINYINT NOT NULL DEFAULT 2
    )
  ")

  DBI::dbExecute(db_conn, "
    CREATE TABLE trait_var_comp (
      id_trait_var_comp INTEGER PRIMARY KEY,
      effect_name       VARCHAR,
      trait_name_1      VARCHAR,
      trait_name_2      VARCHAR,
      cov_value         DOUBLE
    )
  ")

  DBI::dbExecute(db_conn, "
    CREATE TABLE genome_effects (
      id_genome_effect   INTEGER PRIMARY KEY,
      locus_name         VARCHAR NOT NULL,
      line_name          VARCHAR,
      trait_name         VARCHAR NOT NULL,
      genome_effect_type VARCHAR NOT NULL,
      genome_value       DOUBLE  NOT NULL,
      base_allele_freq   DOUBLE
    )
  ")

  DBI::dbExecute(db_conn, "
    CREATE TABLE phenotype_meta (
      id_phenotype_meta        INTEGER PRIMARY KEY,
      phenotype_name           VARCHAR UNIQUE NOT NULL,
      type                     VARCHAR,
      mean                     DOUBLE DEFAULT 0,
      expressed_sex            VARCHAR DEFAULT 'both',
      repeatable               BOOLEAN DEFAULT FALSE,
      min_value                DOUBLE,
      max_value                DOUBLE,
      prevalence               DOUBLE,
      thresholds               VARCHAR,
      cat_values               VARCHAR,
      cat_names                VARCHAR,
      store_liability          BOOLEAN DEFAULT FALSE,
      missing_component_action VARCHAR DEFAULT 'skip',
      formula_tbv              VARCHAR,
      formula                  VARCHAR
    )
  ")

  DBI::dbExecute(db_conn, "
    CREATE TABLE phenotype_components (
      id_phenotype_comp   INTEGER PRIMARY KEY,
      phenotype_name      VARCHAR NOT NULL,
      source_trait_name   VARCHAR NOT NULL,
      contributor_type    VARCHAR NOT NULL,
      group_column        VARCHAR,
      group_table         VARCHAR DEFAULT 'ind_meta',
      aggregation         VARCHAR DEFAULT 'sum',
      weight              DOUBLE  DEFAULT 1.0,
      weight_type         VARCHAR DEFAULT 'fixed',
      covariate_name      VARCHAR,
      covariate_table     VARCHAR,
      poly_order          INTEGER,
      poly_scale_min      DOUBLE,
      poly_scale_max      DOUBLE,
      genome_effect_types VARCHAR DEFAULT 'additive',
      missing_action      VARCHAR DEFAULT 'skip',
      contributor_filter  VARCHAR
    )
  ")

  DBI::dbExecute(db_conn, "
    CREATE TABLE phenotype_var_comp (
      id_phenotype_var_comp INTEGER PRIMARY KEY,
      effect_name           VARCHAR NOT NULL DEFAULT 'residual',
      phenotype_name_1      VARCHAR NOT NULL,
      phenotype_name_2      VARCHAR NOT NULL,
      cov_value             DOUBLE  NOT NULL,
      condition_column      VARCHAR,
      condition_table       VARCHAR DEFAULT 'ind_meta',
      condition_level       VARCHAR,
      weight_type           VARCHAR DEFAULT 'fixed',
      poly_order            INTEGER
    )
  ")

  invisible(NULL)
}


#' Migrate old variance-component table names to the v0.42.0 naming convention
#'
#' Renames `trait_effect_cov` → `trait_var_comp` and
#' `phenotype_residual_cov` → `phenotype_var_comp` in an existing database,
#' adding the `effect_name` column (backfilled as 'residual') when needed.
#' Called automatically by [restore_pop()].
#'
#' @param con An active DuckDB connection.
#' @keywords internal
.migrate_var_comp_tables <- function(con) {
  existing <- DBI::dbListTables(con)

  if ("trait_effect_cov" %in% existing && !"trait_var_comp" %in% existing) {
    DBI::dbExecute(con, "ALTER TABLE trait_effect_cov RENAME TO trait_var_comp")
    tryCatch(
      DBI::dbExecute(con,
        "ALTER TABLE trait_var_comp RENAME id_trait_effect_cov TO id_trait_var_comp"),
      error = function(e) NULL  # column rename not supported in older DuckDB; harmless
    )
    message("[tidybreed] Migrated trait_effect_cov -> trait_var_comp")
  }

  if ("phenotype_residual_cov" %in% existing && !"phenotype_var_comp" %in% existing) {
    DBI::dbExecute(con, "ALTER TABLE phenotype_residual_cov RENAME TO phenotype_var_comp")
    tryCatch(
      DBI::dbExecute(con,
        "ALTER TABLE phenotype_var_comp RENAME id_residual_cov TO id_phenotype_var_comp"),
      error = function(e) NULL
    )
    # Add effect_name column if not yet present; backfill existing rows as 'residual'
    col_info <- DBI::dbGetQuery(con, "PRAGMA table_info('phenotype_var_comp')")
    if (!"effect_name" %in% col_info$name) {
      DBI::dbExecute(con,
        "ALTER TABLE phenotype_var_comp ADD COLUMN effect_name VARCHAR DEFAULT 'residual'")
      DBI::dbExecute(con,
        "UPDATE phenotype_var_comp SET effect_name = 'residual' WHERE effect_name IS NULL")
    }
    message("[tidybreed] Migrated phenotype_residual_cov -> phenotype_var_comp")
  }

  invisible(NULL)
}
