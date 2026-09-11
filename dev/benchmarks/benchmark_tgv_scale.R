#!/usr/bin/env Rscript
# ------------------------------------------------------------------------
# Gate 51 of plans/update_genome_effects_v4.md: evaluation must not regress.
#
# The structural half of that gate -- statement count independent of the
# number of individuals -- is asserted in
# tests/testthat/test-genome-effects-eval.R, where it runs on every suite
# execution and catches a per-individual loop even on an 8-individual
# fixture. This script is the wall-clock half, which is too slow for CI.
#
# It times order-one additive TBV through the term/member evaluator at
# increasing population size and reports the per-individual cost. The claim
# under test is that the cost per individual is flat: the resolved variant
# map is built once from the effect model and the label alphabet, so growing
# the population adds rows to one GROUP BY rather than work to R.
#
# Three model shapes, because they exercise different parts of the evaluator:
#   common      one unscoped variant per locus  -- the fast path
#   lines       common + two line-specific variants -- containment search
#   dominance   common additive + a dominance term -- the genotype grain,
#               including the materialized zero-copy state
#
# Run manually:
#   Rscript dev/benchmarks/benchmark_tgv_scale.R
#   TIDYBREED_BENCH_LARGE=1 Rscript dev/benchmarks/benchmark_tgv_scale.R
# ------------------------------------------------------------------------

devtools::load_all(quiet = TRUE)
suppressMessages(library(dplyr))

run_large <- identical(Sys.getenv("TIDYBREED_BENCH_LARGE"), "1")
sizes  <- if (run_large) c(500L, 2000L, 8000L) else c(200L, 800L)
n_loci <- 2000L
n_qtl  <- 500L

timed <- function(expr) {
  t <- proc.time()
  value <- force(expr)
  list(value = value, elapsed = (proc.time() - t)[["elapsed"]])
}

build_pop <- function(n_ind, shape) {
  pop <- open_pop(pop_name = paste0("bench_", shape, "_", n_ind),
                  db_name = ":memory:") |>
    define_genome(n_loci = n_loci, n_chr = 10L, chr_len_Mb = 100)
  pop <- define_founder_haplotypes(pop, n_haplotypes = 100L, line_name = "A",
                                   method = "fixed", allele_freq = 0.5)
  pop <- define_founder_haplotypes(pop, n_haplotypes = 100L, line_name = "B",
                                   method = "fixed", allele_freq = 0.5)
  for (ln in c("A", "B")) {
    pop <- pop |> get_table("founder_haplotypes") |>
      filter(line_name == ln) |>
      add_founders(n_males = n_ind / 4, n_females = n_ind / 4, line_name = ln)
  }
  pop <- define_trait(pop, "ADG", target_add_var = 1.0)

  qtl <- pop |> get_table("genome_meta") |> collect() |>
    arrange(locus_id) |> slice_head(n = n_qtl) |> pull(locus_name)
  qtl_tbl <- function() get_table(pop, "genome_meta") |> filter(locus_name %in% !!qtl)

  pop <- suppressWarnings(
    qtl_tbl() |> define_additive_effects("ADG", effects = rep(1.0, n_qtl)))
  if (shape == "lines") {
    pop <- qtl_tbl() |> define_additive_effects("ADG", effects = rep(2.0, n_qtl),
                                                line_name = "A")
    pop <- qtl_tbl() |> define_additive_effects("ADG", effects = rep(3.0, n_qtl),
                                                line_name = "B")
  }
  if (shape == "dominance") {
    pop <- define_genome_effects(pop, "ADG", data.frame(
      locus_name = qtl[seq_len(50)], contrast_name = "dominance",
      center_value = 0.5, genome_value = 0.4), effect_owner = "dom")
  }
  pop
}

cat("tidybreed TGV/TBV evaluation scale benchmark\n")
cat("loci:", n_loci, " QTL:", n_qtl, " sizes:", paste(sizes, collapse = ", "), "\n\n")

for (shape in c("common", "lines", "dominance")) {
  cat("== model shape:", shape, "==\n")
  for (n_ind in sizes) {
    pop <- suppressMessages(build_pop(n_ind, shape))
    ids <- get_table(pop, "ind_meta")

    tbv <- timed(suppressMessages(add_tbv(ids, "ADG")))
    tgv <- timed(suppressMessages(add_tgv(ids, "ADG")))
    cat(sprintf("  n = %6d   add_tbv %7.3fs (%6.2f ms/ind)   add_tgv %7.3fs (%6.2f ms/ind)\n",
                n_ind, tbv$elapsed, 1000 * tbv$elapsed / n_ind,
                tgv$elapsed, 1000 * tgv$elapsed / n_ind))
    close_pop(pop)
  }
  cat("\n")
}
cat("Per-individual cost should be flat or falling across sizes. A rising\n")
cat("ms/ind means work is being done per individual somewhere it should not.\n")
