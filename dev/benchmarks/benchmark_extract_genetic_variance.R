#!/usr/bin/env Rscript
# ------------------------------------------------------------------------
# extract_genetic_variance() at a realistic size, and with every pair of loci
# (plans/import_qtl_effect_methods_phase_4_plan.md, "Risks and things to
# watch"; Codex review finding 4).
#
# Case "typical": 2,000 individuals, 500 loci, 2 traits, A + D + A x A with
#   1,000 pairs, realised and genic.
# Case "all_pairs": the same cohort, every pair of the first 200 loci (19,900
#   pairs). A literal port of nonadd_covariates() would build a 2,000 x 19,900
#   double matrix (~318 MB) plus a 19,900^2 pair covariance (~3.2 GB); the
#   extractor accumulates the A x A values in pair chunks (two here) instead.
#   All 124,750 pairs of the 500 loci are not used: writing them through
#   define_genome_effect_terms() takes most of an hour (the writer's per-term
#   R loop is the bottleneck, ~2.9 ms per pair at 4,000 pairs), which is a
#   writer cost, not the extractor's. The write is timed separately below.
#
# Deterministic (fixed seeds). Not part of R CMD check. Run manually:
#
#   Rscript dev/benchmarks/benchmark_extract_genetic_variance.R
# ------------------------------------------------------------------------

devtools::load_all(quiet = TRUE)

n_ind  <- 2000L
n_loci <- 500L
loci   <- sprintf("L%03d", seq_len(n_loci))

set.seed(2026)
pop <- suppressMessages(
  open_pop(pop_name = "bench_egv", db_name = ":memory:") |>
    define_genome(n_loci = n_loci, n_chr = 5, chr_len_Mb = 100,
                  locus_names = loci) |>
    define_founder_haplotypes(n_haplotypes = 400) |>
    get_table("founder_haplotypes") |>
    add_founders(n_males = n_ind / 2, n_females = n_ind / 2, line_name = "A"))

write_model <- function(pop, trait, pairs, seed) {
  set.seed(seed)
  terms <- rbind(
    suppressMessages(ad_terms(loci, a = stats::rnorm(n_loci, sd = 0.1),
                              d = stats::rnorm(n_loci, sd = 0.05), p = 0.5)),
    aa_terms(loci[pairs[, 1]], loci[pairs[, 2]],
             e = stats::rnorm(nrow(pairs), sd = 0.02),
             p_1 = 0.5, p_2 = 0.5, report = FALSE))
  pop <- suppressMessages(define_trait(pop, trait))
  suppressMessages(define_genome_effect_terms(pop, trait, terms))
}

set.seed(1)
every <- t(utils::combn(n_loci, 2))
some_pairs <- every[sort(sample(nrow(every), 1000L)), , drop = FALSE]
all_pairs <- t(utils::combn(200L, 2))

pop <- write_model(pop, "T1", some_pairs, 11)
pop <- write_model(pop, "T2", some_pairs, 12)
t_write <- system.time(pop <- write_model(pop, "ALL", all_pairs, 13))[["elapsed"]]
cat(sprintf("writing the all_pairs model (%d pairs): %.1f s\n", nrow(all_pairs),
            t_write))

time_it <- function(label, expr) {
  gc(reset = TRUE)
  t <- system.time(res <- suppressMessages(expr))[["elapsed"]]
  g <- gc()
  # Peak R heap since the reset: max-used cons cells (56 bytes) and vector
  # cells (8 bytes).
  mem <- (g[1, "max used"] * 56 + g[2, "max used"] * 8) / 2^20
  cat(sprintf("%-34s %8.2f s   peak R memory %8.1f Mb   %d rows\n",
              label, t, mem, nrow(res)))
  invisible(res)
}

cat("n_ind =", n_ind, " n_loci =", n_loci, "\n")
cat("typical: 2 traits, 1,000 pairs | all_pairs: 1 trait,",
    format(nrow(all_pairs), big.mark = ","), "pairs of 200 loci\n\n")
cohort <- get_table(pop, "ind_meta")
time_it("typical, realised",  extract_genetic_variance(cohort, c("T1", "T2")))
time_it("typical, genic",     extract_genetic_variance(cohort, c("T1", "T2"),
                                                       anchor = "genic"))
time_it("all_pairs, realised", extract_genetic_variance(cohort, "ALL"))
time_it("all_pairs, genic",    extract_genetic_variance(cohort, "ALL",
                                                        anchor = "genic"))
close_pop(pop)
