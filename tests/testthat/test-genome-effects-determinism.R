# CLAUDE.md's first reproducibility contract -- "same seed reproduces within
# the current code" -- read literally, not within a tolerance. The RNG stream
# was never the problem; the genetic values were. DuckDB's parallel hash
# aggregate combined the per-term partial sums in whatever order the threads
# finished, so two runs of an identical population differed in the last bits
# of tbv_value (~1e-15 at 2000 loci), and pheno_value inherited it through the
# TBV term.
#
# The evaluator now accumulates that one sum exactly (GEV_ACC_TYPE), so the
# result is a function of the stored model alone. These tests fail on the
# pre-change code: the margin grows with the number of terms, which is why
# they use a QTL count where the old path diverged within a handful of runs.
#
# See plans/sample_correlated_effects.md D8 and the header of
# R/genome_effects_eval.R.

det_pop <- function(name, n_loci = 500L, n_ind = 40L) {
  set.seed(4)
  pop <- suppressMessages(make_test_pop(name, n_loci = n_loci, n_chr = 1,
                                        n_males = n_ind / 2,
                                        n_females = n_ind / 2))
  pop <- suppressMessages(define_trait(pop, "A", target_add_var = 1))
  suppressMessages(pop |> get_table("genome_meta") |> define_additive_effects("A"))
}

# Repeated evaluation of one unchanged population. `expect_identical()` on
# doubles is the whole point -- `expect_equal()` would pass on the old code.
test_that("add_tbv() is bit-identical across repeated runs of the same population", {
  pop <- det_pop("det_tbv")
  on.exit(close_pop(pop), add = TRUE)

  once <- function() {
    suppressMessages(pop |> get_table("ind_meta") |> add_tbv("A"))
    DBI::dbGetQuery(pop$db_conn,
      "SELECT tbv_value FROM ind_tbv ORDER BY id_ind")$tbv_value
  }
  first <- once()
  expect_length(first, 40)
  for (i in 1:7) expect_identical(once(), first)
})

test_that("add_tgv() is bit-identical across repeated runs of the same population", {
  pop <- det_pop("det_tgv")
  on.exit(close_pop(pop), add = TRUE)

  once <- function() {
    suppressMessages(pop |> get_table("ind_meta") |> add_tgv("A"))
    DBI::dbGetQuery(pop$db_conn,
      "SELECT tgv_value FROM ind_tgv ORDER BY id_ind, component_name")$tgv_value
  }
  first <- once()
  for (i in 1:7) expect_identical(once(), first)
})

# The evaluator feeds add_phenotype()'s TBV term, so a wobble there reached
# pheno_value even though every draw was already reproducible. Two separately
# built populations, same build seed and same phenotype seed.
test_that("a seeded add_phenotype() is bit-identical between two identical populations", {
  build <- function(name) {
    pop <- det_pop(name)
    pop <- suppressMessages(define_phenotype(
      pop, "A", type = "continuous", mean = 10, residual_var = 1))
    suppressMessages(pop |> get_table("ind_meta") |> add_phenotype("A", seed = 9))
    pop
  }
  p1 <- build("det_ph_1"); on.exit(close_pop(p1), add = TRUE)
  p2 <- build("det_ph_2"); on.exit(close_pop(p2), add = TRUE)

  vals <- function(pop) DBI::dbGetQuery(pop$db_conn,
    "SELECT pheno_value, residual_value FROM ind_phenotype ORDER BY id_ind")

  a <- vals(p1); b <- vals(p2)
  expect_identical(a$residual_value, b$residual_value)  # never was the problem
  expect_identical(a$pheno_value,    b$pheno_value)     # this one was
})

# The exact accumulator has a finite range, so a model whose values fall
# outside it must fail as a tidybreed error naming the cause, not as a bare
# DuckDB conversion error.
test_that("a genetic value too large for the exact accumulator errors with a tidybreed message", {
  pop <- suppressMessages(make_test_pop("det_big", n_loci = 4, n_chr = 1,
                                        n_males = 1, n_females = 1))
  on.exit(close_pop(pop), add = TRUE)
  pop <- suppressMessages(define_trait(pop, "A", target_add_var = 1))
  pop <- suppressMessages(define_genome_effects(pop, "A", data.frame(
    locus_name = "Locus_1", contrast_name = "additive",
    center_value = 0.5, genome_value = 1e30)))

  expect_error(
    suppressMessages(pop |> get_table("ind_meta") |> add_tgv("A")),
    "too large to evaluate exactly")
})
