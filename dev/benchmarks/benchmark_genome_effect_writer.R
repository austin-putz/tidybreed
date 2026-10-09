#!/usr/bin/env Rscript
# ------------------------------------------------------------------------
# The genome-effect writer at scale (step 5a of
# plans/import_qtl_effect_methods_phase_5_plan.md).
#
# 1. Fresh writes of 1,000 / 4,000 / 16,000 / 64,000 A x A pairs, and all
#    124,750 pairs of 500 loci, through define_genome_effect_terms(). Time
#    per pair should grow at most 2x from 4,000 to 64,000 (near linear), and
#    the 124,750-pair write should take under 60 s.
# 2. Replacement into a populated model: re-write one trait's 16,000 pairs
#    while the database holds two other traits' models plus line-scoped and
#    multi-owner terms. The writer validates the whole table before COMMIT,
#    so its cost depends on what is already stored.
# 3. .stored_to_functional() on the largest model (the extractor's
#    conversion).
#
# Setup (building the terms frame) is timed apart from the write. Peak R
# memory is gc()'s "max used" since the last reset.
#
# Deterministic (fixed seeds). Not part of R CMD check. Run manually:
#
#   Rscript dev/benchmarks/benchmark_genome_effect_writer.R [max_pairs]
#
# `max_pairs` (default all) skips the larger sizes, for timing the slow
# pre-5a writer without waiting an hour.
# ------------------------------------------------------------------------

devtools::load_all(quiet = TRUE)

args      <- commandArgs(trailingOnly = TRUE)
max_pairs <- if (length(args)) as.numeric(args[1]) else Inf

n_loci <- 500L
loci   <- sprintf("L%03d", seq_len(n_loci))
every  <- t(utils::combn(n_loci, 2))          # 124,750 pairs

new_pop <- function() {
  suppressMessages(
    open_pop(pop_name = "bench_gew", db_name = ":memory:") |>
      define_genome(n_loci = n_loci, n_chr = 5, chr_len_Mb = 100,
                    locus_names = loci))
}

pair_terms <- function(pairs, seed) {
  set.seed(seed)
  aa_terms(loci[pairs[, 1]], loci[pairs[, 2]],
           e = stats::rnorm(nrow(pairs), sd = 0.02),
           p_1 = 0.5, p_2 = 0.5, report = FALSE)
}

# "max used" in Mb, Ncells + Vcells: the column after "max used" (gc() may
# also report a "limit (Mb)" column, so never index it by position).
peak_mb <- function() {
  g <- gc()
  sum(g[, which(colnames(g) == "max used") + 1L])
}

hw <- Sys.info()
cat(sprintf("R %s | %s %s | %s | %d cores\n", getRversion(), hw[["sysname"]],
            hw[["release"]], hw[["machine"]], parallel::detectCores()))
probe <- new_pop()
cat(sprintf("DuckDB %s, threads = %s\n\n", utils::packageVersion("duckdb"),
            DBI::dbGetQuery(probe$db_conn,
                            "SELECT current_setting('threads') AS t")$t))
close_pop(probe)

# ── 1. Fresh writes ─────────────────────────────────────────────────────────

sizes <- c(1000, 4000, 16000, 64000, nrow(every))
sizes <- sizes[sizes <= max_pairs]
set.seed(1)
res <- data.frame()
for (n in sizes) {
  pop <- new_pop()
  pop <- suppressMessages(define_trait(pop, "T"))
  pairs <- if (n == nrow(every)) every else
    every[sort(sample(nrow(every), n)), , drop = FALSE]
  invisible(gc(reset = TRUE))
  t_setup <- system.time(terms <- pair_terms(pairs, 11))[["elapsed"]]
  t_write <- system.time(
    suppressMessages(define_genome_effect_terms(pop, "T", terms)))[["elapsed"]]
  res <- rbind(res, data.frame(pairs = n, setup_s = t_setup, write_s = t_write,
                               us_per_pair = 1e6 * t_write / n,
                               peak_mb = peak_mb()))
  close_pop(pop)
}
cat("Fresh writes (one trait, A x A pairs only)\n")
print(res, row.names = FALSE, digits = 4)
if (all(c(4000, 64000) %in% res$pairs)) {
  g <- res$us_per_pair[res$pairs == 64000] / res$us_per_pair[res$pairs == 4000]
  cat(sprintf("growth in time per pair, 4,000 -> 64,000: %.2fx (target <= 2)\n", g))
}
if (nrow(every) %in% res$pairs) {
  cat(sprintf("all %d pairs: %.1f s (target < 60 s)\n", nrow(every),
              res$write_s[res$pairs == nrow(every)]))
}

# ── 2. Replacement into a populated model ───────────────────────────────────

if (16000 <= max_pairs) {
  pop <- new_pop()
  set.seed(2)
  other <- every[sort(sample(nrow(every), 1000L)), , drop = FALSE]
  for (tr in c("O1", "O2")) {
    pop <- suppressMessages(define_trait(pop, tr))
    set.seed(match(tr, c("O1", "O2")))
    pop <- suppressMessages(define_genome_effect_terms(pop, tr, rbind(
      suppressMessages(ad_terms(loci, a = stats::rnorm(n_loci, sd = 0.1),
                                d = stats::rnorm(n_loci, sd = 0.05), p = 0.5)),
      pair_terms(other, 20))))
  }
  pop <- suppressMessages(define_trait(pop, "T"))
  # Line-scoped additive variants and a second owner on the replaced trait.
  scoped <- suppressMessages(ad_terms(loci[1:100], a = 0.01, p = 0.5))
  pop <- suppressMessages(define_genome_effect_terms(
    pop, "T", scoped, origin = list(line_name = "A"), effect_owner = "lineA"))
  pop <- suppressMessages(define_genome_effect_terms(
    pop, "T", suppressMessages(ad_terms(loci[101:200], a = 0.02, p = 0.5)),
    effect_owner = "other_owner"))
  pairs16 <- every[sort(sample(nrow(every), 16000L)), , drop = FALSE]
  pop <- suppressMessages(define_genome_effect_terms(pop, "T",
                                                     pair_terms(pairs16, 30)))
  stored <- DBI::dbGetQuery(pop$db_conn,
                            "SELECT COUNT(*) AS n FROM genome_effects")$n
  terms <- pair_terms(pairs16, 31)
  invisible(gc(reset = TRUE))
  t_rep <- system.time(suppressMessages(define_genome_effect_terms(
    pop, "T", terms, mode = "replace_owner")))[["elapsed"]]
  cat(sprintf(paste0("\nReplacement: 16,000 pairs of trait T (owner 'custom') ",
                     "into %d stored terms (3 traits, 3 owners, 100 line-",
                     "scoped): %.1f s, peak %.0f Mb\n"), stored, t_rep, peak_mb()))

  # ── 3. .stored_to_functional() ────────────────────────────────────────────
  model <- .gev_read_model(pop$db_conn, "T")
  cov_t <- model$terms[model$terms$effect_owner == "custom", , drop = FALSE]
  cov_m <- model$members[model$members$id_genome_effect %in%
                           cov_t$id_genome_effect, , drop = FALSE]
  t_s2f <- system.time(.stored_to_functional(cov_t, cov_m))[["elapsed"]]
  cat(sprintf(".stored_to_functional(): %d terms in %.2f s\n", nrow(cov_t),
              t_s2f))
  close_pop(pop)
}
