#!/usr/bin/env Rscript
# ------------------------------------------------------------------------
# Long-format storage scale benchmark (deferred Stage 1 exit criterion).
#
# plans/refactor_haplotype.md's "Benchmark plan" section originally asked for
# a wide-vs-long comparison before locking the schema. Stage 1 (tidybreed
# 0.47.0) already shipped as a greenfield rewrite with no compatibility layer
# -- the old wide tables (genome_haplotype/genome_genotype) no longer exist
# anywhere in the codebase and there is no pre-refactor git tag to check out.
# Re-litigating "wide vs. long" is moot; that decision already shipped.
#
# This script instead validates the CURRENT long-format implementation at
# scale: insert throughput for add_founders()/add_offspring(), the Stage 2
# line-origin-aware add_tbv() query (both the population-wide-only and
# line-specific-with-fallback branches), and the PIVOT-based
# extract_genotypes() export. See plans/refactor_haplotype_stage_2.md.
#
# Not part of R CMD check or the testthat suite -- multi-thousand-individual
# runs are too slow for CI. Run manually:
#
#   Rscript dev/benchmarks/benchmark_haplotype_scale.R
#   TIDYBREED_BENCH_LARGE=1 Rscript dev/benchmarks/benchmark_haplotype_scale.R
# ------------------------------------------------------------------------

devtools::load_all(quiet = TRUE)

run_large <- identical(Sys.getenv("TIDYBREED_BENCH_LARGE"), "1")

scales <- list(
  small = list(n_loci = 5000L,  n_chr = 10L, n_males = 250L, n_females = 250L)
)
if (run_large) {
  scales$large <- list(n_loci = 50000L, n_chr = 20L, n_males = 2500L, n_females = 2500L)
}

timed <- function(expr) {
  t <- proc.time()
  value <- force(expr)
  elapsed <- (proc.time() - t)[["elapsed"]]
  list(value = value, elapsed = elapsed)
}

run_scale <- function(scale_name, cfg) {
  cat("\n==", scale_name, "scale: n_ind =", cfg$n_males + cfg$n_females,
      " n_loci =", cfg$n_loci, "==\n")

  results <- list()

  pop <- open_pop(pop_name = paste0("bench_", scale_name), db_name = ":memory:") |>
    define_genome(n_loci = cfg$n_loci, n_chr = cfg$n_chr, chr_len_Mb = 100)

  set.seed(1)
  pop <- define_founder_haplotypes(pop, n_haplotypes = cfg$n_males + cfg$n_females,
                                   line_name = "A", method = "uniform")
  set.seed(2)
  pop <- define_founder_haplotypes(pop, n_haplotypes = cfg$n_males + cfg$n_females,
                                   line_name = "B", method = "uniform")

  # -- add_founders() insert throughput --------------------------------------
  t_founders_a <- timed({
    pop <<- pop |>
      get_table("founder_haplotypes") |>
      dplyr::filter(line_name == "A") |>
      add_founders(n_males = cfg$n_males, n_females = cfg$n_females,
                   line_name = "A", gen = 0L)
  })
  t_founders_b <- timed({
    pop <<- pop |>
      get_table("founder_haplotypes") |>
      dplyr::filter(line_name == "B") |>
      add_founders(n_males = cfg$n_males, n_females = cfg$n_females,
                   line_name = "B", gen = 0L)
  })
  n_founder_rows <- 2L * cfg$n_loci * (cfg$n_males + cfg$n_females) * 2L
  founders_elapsed <- t_founders_a$elapsed + t_founders_b$elapsed
  results$add_founders <- list(
    rows = n_founder_rows, elapsed = founders_elapsed,
    rows_per_sec = n_founder_rows / founders_elapsed
  )

  # -- add_offspring() insert throughput (one F1 generation) ------------------
  sires <- paste0("A_", seq_len(cfg$n_males))
  dams  <- paste0("B_", cfg$n_males + seq_len(cfg$n_females))
  matings_f1 <- tibble::tibble(
    id_parent_1 = sires,
    id_parent_2 = dams,
    sex         = rep(c("M", "F"), length.out = length(sires)),
    line_name   = "F1",
    gen         = 1L
  )
  t_offspring <- timed({
    pop <<- add_offspring(pop, matings_f1)
  })
  n_offspring_rows <- 2L * cfg$n_loci * nrow(matings_f1)
  results$add_offspring <- list(
    rows = n_offspring_rows, elapsed = t_offspring$elapsed,
    rows_per_sec = n_offspring_rows / t_offspring$elapsed
  )

  # -- add_tbv(): population-wide-only effects ---------------------------------
  pop <- define_trait(pop, "ADG_popwide", target_add_var = 1.0)
  set.seed(3)
  pop <- pop |>
    get_table("genome_meta") |>
    define_additive_effects("ADG_popwide", distribution = "normal")

  t_tbv_popwide <- timed({
    pop <<- pop |>
      get_table("ind_meta") |>
      add_tbv("ADG_popwide")
  })
  n_ind_total <- DBI::dbGetQuery(pop$db_conn, "SELECT COUNT(*) AS n FROM ind_meta")$n
  results$add_tbv_popwide <- list(
    n_ind = n_ind_total, elapsed = t_tbv_popwide$elapsed,
    ind_per_sec = n_ind_total / t_tbv_popwide$elapsed
  )

  # -- add_tbv(): line-specific effects with fallback (Stage 2 query) ---------
  pop <- define_trait(pop, "ADG_lineSpec", target_add_var = 1.0)
  set.seed(4)
  qtl_tbl <- pop |> get_table("genome_meta") |> dplyr::filter(locus_id %% 10L == 0L)
  pop <- define_additive_effects(qtl_tbl, "ADG_lineSpec",
                                 distribution = "normal", line_name = "A")
  set.seed(5)
  pop <- define_additive_effects(qtl_tbl, "ADG_lineSpec",
                                 distribution = "normal", line_name = "B")

  t_tbv_linespec <- timed({
    pop <<- pop |>
      get_table("ind_meta") |>
      add_tbv("ADG_lineSpec")
  })
  results$add_tbv_linespec <- list(
    n_ind = n_ind_total, elapsed = t_tbv_linespec$elapsed,
    ind_per_sec = n_ind_total / t_tbv_linespec$elapsed
  )

  # -- extract_genotypes(): PIVOT export for a chip-sized locus subset ---------
  chip_loci <- pop |>
    get_table("genome_meta") |>
    dplyr::filter(locus_id %% 10L == 0L)
  pop <- define_chip(chip_loci, "bench_chip")
  pop <- pop |> get_table("ind_meta") |> add_genotypes("bench_chip")

  t_export <- timed({
    export_df <<- pop |>
      get_table("ind_meta") |>
      extract_genotypes(chip_name = "bench_chip")
  })
  results$extract_genotypes <- list(
    n_ind = nrow(export_df), n_loci = ncol(export_df) - 1L,
    elapsed = t_export$elapsed
  )

  close_pop(pop)
  results
}

all_results <- lapply(names(scales), function(nm) run_scale(nm, scales[[nm]]))
names(all_results) <- names(scales)

cat("\n\n==== Summary ====\n")
for (scale_name in names(all_results)) {
  cat("\n--", scale_name, "--\n")
  res <- all_results[[scale_name]]
  cat(sprintf("  add_founders():        %10d rows in %7.2fs  (%.0f rows/sec)\n",
              res$add_founders$rows, res$add_founders$elapsed, res$add_founders$rows_per_sec))
  cat(sprintf("  add_offspring():       %10d rows in %7.2fs  (%.0f rows/sec)\n",
              res$add_offspring$rows, res$add_offspring$elapsed, res$add_offspring$rows_per_sec))
  cat(sprintf("  add_tbv() pop-wide:    %10d ind  in %7.2fs  (%.0f ind/sec)\n",
              res$add_tbv_popwide$n_ind, res$add_tbv_popwide$elapsed, res$add_tbv_popwide$ind_per_sec))
  cat(sprintf("  add_tbv() line-spec:   %10d ind  in %7.2fs  (%.0f ind/sec)\n",
              res$add_tbv_linespec$n_ind, res$add_tbv_linespec$elapsed, res$add_tbv_linespec$ind_per_sec))
  cat(sprintf("  extract_genotypes():   %10d ind x %d loci in %7.2fs\n",
              res$extract_genotypes$n_ind, res$extract_genotypes$n_loci, res$extract_genotypes$elapsed))
}

if (!run_large) {
  cat("\n(Set TIDYBREED_BENCH_LARGE=1 to also run the ~5,000 x 50,000 realistic scale.)\n")
}
