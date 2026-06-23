# Tests for add_ebv() routing logic (managed vs. fallback directory mode).
# Full BLUPF90 end-to-end tests are omitted here because they require the
# renumf90 / blupf90+ binaries to be on PATH.

# ---------------------------------------------------------------------------
# Helper: minimal file-based pop with one trait in trait_meta
# ---------------------------------------------------------------------------
make_managed_pop <- function(tmp, tools) {
  open_pop(
    pop_name     = "t",
    base_dir     = tmp,
    output_dir   = "out",
    scenario_dir = "sc1",
    tools        = tools,
    db_name      = "sim.duckdb"
  ) |>
    define_genome(n_loci = 20, n_chr = 1, chr_len_Mb = 50) |>
    define_founder_haplotypes(n_haplotypes = 10) |>
    get_table("founder_haplotypes") |>
    add_founders(n_males = 2, n_females = 2, line_name = "A") |>
    define_trait("ADG", target_add_var = 0.25)
}

# ---------------------------------------------------------------------------
# Managed mode: blupf90 not registered → clear error before binary lookup
# ---------------------------------------------------------------------------

test_that("add_ebv() errors with clear message when blupf90 not in pop$run_dirs", {
  tmp <- tempfile()
  dir.create(tmp)
  on.exit(unlink(tmp, recursive = TRUE))

  pop <- make_managed_pop(tmp, tools = "plink")
  on.exit(close_pop(pop), add = TRUE)

  expect_error(
    pop |> get_table("ind_meta") |>
      add_ebv("ADG", software = "blupf90"),
    regexp = "blupf90.*not set up"
  )
})

test_that("add_ebv() error message mentions how to fix (add blupf90 to open_pop)", {
  tmp <- tempfile()
  dir.create(tmp)
  on.exit(unlink(tmp, recursive = TRUE))

  pop <- make_managed_pop(tmp, tools = "plink")
  on.exit(close_pop(pop), add = TRUE)

  err <- tryCatch(
    pop |> get_table("ind_meta") |>
      add_ebv("ADG", software = "blupf90"),
    error = function(e) conditionMessage(e)
  )
  expect_match(err, "open_pop")
  expect_match(err, "blupf90")
})

test_that("add_ebv() error lists currently registered tools", {
  tmp <- tempfile()
  dir.create(tmp)
  on.exit(unlink(tmp, recursive = TRUE))

  pop <- make_managed_pop(tmp, tools = "plink")
  on.exit(close_pop(pop), add = TRUE)

  err <- tryCatch(
    pop |> get_table("ind_meta") |>
      add_ebv("ADG", software = "blupf90"),
    error = function(e) conditionMessage(e)
  )
  expect_match(err, "plink")
})

test_that("add_ebv() error fires even when no tools are registered (only base)", {
  tmp <- tempfile()
  dir.create(tmp)
  on.exit(unlink(tmp, recursive = TRUE))

  # open_pop with no tools — run_dirs has only "base"
  pop <- open_pop(
    pop_name     = "t",
    base_dir     = tmp,
    output_dir   = "out",
    scenario_dir = "sc1",
    tools        = NULL,
    db_name      = "sim.duckdb"
  ) |>
    define_genome(n_loci = 20, n_chr = 1, chr_len_Mb = 50) |>
    define_founder_haplotypes(n_haplotypes = 10) |>
    get_table("founder_haplotypes") |>
    add_founders(n_males = 2, n_females = 2, line_name = "A") |>
    define_trait("ADG", target_add_var = 0.25)
  on.exit(close_pop(pop), add = TRUE)

  err <- tryCatch(
    pop |> get_table("ind_meta") |>
      add_ebv("ADG", software = "blupf90"),
    error = function(e) conditionMessage(e)
  )
  expect_match(err, "No tools are currently registered")
})

# ---------------------------------------------------------------------------
# Managed mode: blupf90 registered → warnings for ignored run_dir / eval_id
# (The binary lookup fires after warnings, so we swallow the downstream error)
# ---------------------------------------------------------------------------

test_that("add_ebv() warns when eval_id is passed in managed mode", {
  tmp <- tempfile()
  dir.create(tmp)
  on.exit(unlink(tmp, recursive = TRUE))

  pop <- make_managed_pop(tmp, tools = "blupf90")
  on.exit(close_pop(pop), add = TRUE)

  expect_warning(
    tryCatch(
      pop |> get_table("ind_meta") |>
        add_ebv("ADG", software = "blupf90", eval_id = "my_run"),
      error = function(e) NULL
    ),
    regexp = "eval_id.*ignored"
  )
})

test_that("add_ebv() warns when non-default run_dir is passed in managed mode", {
  tmp <- tempfile()
  dir.create(tmp)
  on.exit(unlink(tmp, recursive = TRUE))

  pop <- make_managed_pop(tmp, tools = "blupf90")
  on.exit(close_pop(pop), add = TRUE)

  expect_warning(
    tryCatch(
      pop |> get_table("ind_meta") |>
        add_ebv("ADG", software = "blupf90", run_dir = "/tmp/mydir"),
      error = function(e) NULL
    ),
    regexp = "run_dir.*ignored"
  )
})

test_that("add_ebv() does NOT warn for default run_dir in managed mode", {
  tmp <- tempfile()
  dir.create(tmp)
  on.exit(unlink(tmp, recursive = TRUE))

  pop <- make_managed_pop(tmp, tools = "blupf90")
  on.exit(close_pop(pop), add = TRUE)

  # Default run_dir = "." should produce no warning about run_dir
  # (there may be other warnings from the binary not being found; we only check
  # that the specific "run_dir.*ignored" warning does NOT appear)
  warnings_caught <- character(0)
  withCallingHandlers(
    tryCatch(
      pop |> get_table("ind_meta") |>
        add_ebv("ADG", software = "blupf90"),
      error = function(e) NULL
    ),
    warning = function(w) {
      warnings_caught <<- c(warnings_caught, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  expect_false(any(grepl("run_dir.*ignored", warnings_caught)))
})

# ---------------------------------------------------------------------------
# Fallback mode: in-memory pop → run_dir / eval_id behaviour unchanged
# ---------------------------------------------------------------------------

test_that("add_ebv() fallback mode creates eval dir under run_dir (not managed)", {
  tmp_run <- tempfile()
  dir.create(tmp_run)
  on.exit(unlink(tmp_run, recursive = TRUE))

  pop <- make_test_pop(n_loci = 20, n_chr = 1, n_males = 2, n_females = 2,
                       n_haplotypes = 10) |>
    define_trait("ADG", target_add_var = 0.25)
  on.exit(close_pop(pop), add = TRUE)

  # pop$run_dirs is empty (in-memory) → fallback path
  expect_equal(length(pop$run_dirs), 0L)

  # We can't run blupf90 in CI, but we can verify the dir is created
  # by calling the code up to the point it fails on the binary.
  # Capture what directory was created under tmp_run.
  tryCatch(
    pop |> get_table("ind_meta") |>
      add_ebv("ADG", software = "blupf90",
              run_dir = tmp_run, eval_id = "test_run"),
    error = function(e) NULL
  )
  created <- list.dirs(tmp_run, recursive = FALSE)
  expect_true(any(grepl("eval_test_run", basename(created))))
})
