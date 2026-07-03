# Parity regression for the wide -> long haplotype/genotype refactor (Stage 1).
#
# On the FIRST run (against the current wide code) this captures golden files
# under tests/testthat/parity_golden/. On EVERY later run (against the long code)
# it compares the freshly-simulated artifacts against those goldens. Any semantic
# drift in haplotypes, dosages, TBVs, or the exported genotype matrix fails here.
#
# To intentionally re-capture goldens (e.g. after an approved behavior change),
# delete tests/testthat/parity_golden/ and re-run the suite.

# Capture-or-compare a single artifact. Returns TRUE if it compared against an
# existing golden, FALSE if it captured a new one.
parity_check_artifact <- function(name, value, tolerance = 1e-8) {
  dir <- parity_golden_dir()
  if (!dir.exists(dir)) dir.create(dir, recursive = TRUE)
  path <- file.path(dir, paste0(name, ".rds"))

  if (!file.exists(path)) {
    saveRDS(value, path)
    return(FALSE)
  }
  golden <- readRDS(path)
  testthat::expect_equal(value, golden, tolerance = tolerance)
  TRUE
}

test_that("parity simulation runs and matches (or captures) golden artifacts", {
  artifacts <- run_parity_sim()

  captured_any <- FALSE
  for (nm in names(artifacts)) {
    compared <- parity_check_artifact(nm, artifacts[[nm]])
    if (!compared) captured_any <- TRUE
  }

  if (captured_any) {
    testthat::skip(
      "Captured new parity golden file(s) under tests/testthat/parity_golden/. Re-run to compare."
    )
  }

  # Sanity: the sim produced the expected individuals (16 founders + 4 F1 + 2 F2).
  expect_setequal(
    unique(artifacts$haplotypes$id_ind),
    c(paste0("A_", 1:8), paste0("B_", 1:8), paste0("F1_", 1:4), paste0("F2_", 1:2))
  )
})
