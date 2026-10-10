#!/usr/bin/env Rscript
# ------------------------------------------------------------------------
# define_genome_effects() at scale (step 5c of
# plans/import_qtl_effect_methods_phase_5_plan.md).
#
# 1. 2,000 base individuals, 1,000 QTL, two traits, A + D + A x A with the
#    default floor(m / 2) = 500 random pairs, under anchor = "genic" (base:
#    the founder haplotypes) and anchor = "realised" (base: the individuals;
#    designs of 2,000 x (1,000 + 1,000 + 500) cells, under the 2e7 guard).
# 2. A hub design: 20,000 supplied pairs (21 hub loci, each paired with the
#    979 non-hub loci, the first 20,000), anchor = "genic".
#
# Each call is timed whole: target resolution, draw, calibration,
# diagnostics, the write and the messages. Peak R memory is gc()'s "max
# used" since the last reset (DuckDB's own memory is not included).
#
# Deterministic (fixed seeds). Not part of R CMD check. Run manually:
#
#   Rscript dev/benchmarks/benchmark_define_genome_effects.R
# ------------------------------------------------------------------------

devtools::load_all(quiet = TRUE)

n_loci <- 1000L
n_ind  <- 2000L
traits <- c("T1", "T2")
nm <- function(M) { dimnames(M) <- list(traits, traits); M }
G_A  <- nm(matrix(c(1, 0.3, 0.3, 2), 2))
G_D  <- nm(matrix(c(0.3, 0.1, 0.1, 0.5), 2))
G_AA <- nm(matrix(c(0.2, 0.05, 0.05, 0.3), 2))

# "max used" in Mb, Ncells + Vcells: the column after "max used".
peak_mb <- function() {
  g <- gc()
  sum(g[, which(colnames(g) == "max used") + 1L])
}

new_pop <- function() {
  set.seed(1)
  suppressMessages({
    p <- open_pop(pop_name = "bench_dge", db_name = ":memory:") |>
      define_genome(n_loci = n_loci, n_chr = 10, chr_len_Mb = 100) |>
      define_founder_haplotypes(n_haplotypes = 2 * n_ind)
    p <- p |> get_table("founder_haplotypes") |>
      add_founders(n_males = n_ind / 2, n_females = n_ind / 2, line_name = "A")
    for (t in traits) p <- define_trait(p, t)
    p
  })
}

run <- function(label, f) {
  invisible(gc(reset = TRUE))
  t <- system.time(suppressWarnings(suppressMessages(f())))[["elapsed"]]
  cat(sprintf("%-46s %7.1f s   peak %6.0f Mb\n", label, t, peak_mb()))
}

cat(sprintf("%d individuals, %d QTL, k = %d\n\n", n_ind, n_loci, length(traits)))

pop <- new_pop()
gm  <- get_table(pop, "genome_meta")
run("genic, A + D + A x A, 500 random pairs", function() {
  set.seed(11)
  define_genome_effects(gm, traits, G_A = G_A, G_D = G_D, G_AA = G_AA)
})
close_pop(pop)

pop <- new_pop()
gm  <- get_table(pop, "genome_meta")
run("realised, A + D + A x A, 500 random pairs", function() {
  set.seed(12)
  define_genome_effects(gm, traits, G_A = G_A, G_D = G_D, G_AA = G_AA,
                        anchor = "realised",
                        base_tbl = get_table(pop, "ind_meta"))
})
close_pop(pop)

pop <- new_pop()
gm  <- get_table(pop, "genome_meta")
loci <- DBI::dbGetQuery(pop$db_conn,
  "SELECT locus_name FROM genome_meta ORDER BY locus_id")$locus_name
hubs <- loci[1:21]
rest <- loci[-(1:21)]
pairs <- data.frame(locus_name_1 = rep(hubs, each = length(rest)),
                    locus_name_2 = rep(rest, times = length(hubs)),
                    stringsAsFactors = FALSE)[1:20000, ]
run("genic, A + D + A x A, 20,000 hub pairs", function() {
  set.seed(13)
  define_genome_effects(gm, traits, G_A = G_A, G_D = G_D, G_AA = G_AA,
                        pairs = pairs)
})
n_terms <- DBI::dbGetQuery(pop$db_conn,
  "SELECT COUNT(*) AS n FROM genome_effects")$n
cat(sprintf("\n(hub model: %d stored terms)\n", n_terms))
close_pop(pop)
