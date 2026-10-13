# CLAUDE.md's first reproducibility contract -- "same seed reproduces within
# the current code" -- read literally, not within a tolerance. The RNG stream
# was never the problem; the genetic values were. DuckDB's parallel hash
# aggregate combined the per-term partial sums in whatever order the threads
# finished, so two runs of an identical population differed in the last bits
# of the breeding value (~1e-15 at 2000 loci), and pheno_value inherited it through the
# genetic-value term.
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
  pop <- suppressMessages(with_additive_target(pop, "A", 1))
  suppressMessages(pop |> get_table("genome_meta") |> define_additive_effects("A"))
}

# Repeated evaluation of one unchanged population. `expect_identical()` on
# doubles is the whole point -- `expect_equal()` would pass on the old code.
test_that("add_tgv() is bit-identical across repeated runs of the same population", {
  pop <- det_pop("det_tbv")
  on.exit(close_pop(pop), add = TRUE)

  once <- function() {
    suppressMessages(pop |> get_table("ind_meta") |> add_tgv("A"))
    DBI::dbGetQuery(pop$db_conn,
      "SELECT tgv_value FROM ind_tgv WHERE component_name = 'additive' ORDER BY id_ind")$tgv_value
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
  pop <- suppressMessages(with_additive_target(pop, "A", 1))
  pop <- suppressMessages(define_genome_effect_terms(pop, "A", data.frame(
    locus_name = "Locus_1", contrast_name = "additive",
    center_value = 0.5, genome_value = 1e30)))

  expect_error(
    suppressMessages(pop |> get_table("ind_meta") |> add_tgv("A")),
    "too large to evaluate exactly")
})

# The ids are written values too. The evaluator's GROUP BY returns rows in
# whatever order DuckDB's threads finish, and `id_tgv` is numbered in row
# order, so two identical populations got the same genetic values under
# different `id_tgv`s until the evaluator sorted its result.
test_that("add_tgv() numbers id_tgv identically in two identical populations", {
  run <- function(name) {
    set.seed(5)
    pop <- suppressMessages(make_test_pop(name, n_loci = 60, n_chr = 2,
                                          n_males = 100, n_females = 100))
    on.exit(close_pop(pop), add = TRUE)
    for (t in c("A", "B")) {
      pop <- suppressMessages(define_trait(pop, t))
      pop <- suppressMessages(define_genome_effect_terms(pop, t, rbind(
        ad_terms(c("Locus_10", "Locus_44"), a = c(0.3, -0.1), d = c(0.1, 0.05),
                 p = c(0.35, 0.6), report = FALSE),
        aa_terms("Locus_10", "Locus_44", e = 0.08, p_1 = 0.35, p_2 = 0.6,
                 report = FALSE))))
    }
    suppressMessages(pop |> get_table("ind_meta") |> add_tgv(c("A", "B")))
    DBI::dbGetQuery(pop$db_conn, paste(
      "SELECT id_tgv, id_ind, trait_name, component_name, tgv_value",
      "FROM ind_tgv ORDER BY id_tgv"))
  }
  first <- run("det_id_1")
  expect_equal(nrow(first), 200L * 2L * 3L)
  for (i in 2:4) expect_identical(run(paste0("det_id_", i)), first)
})
