#!/usr/bin/env Rscript
# ------------------------------------------------------------------------
# Phase 0 per-step profiler for the haplotype write paths.
#
#   Rscript dev/benchmarks/profile_haplotype_paths.R
#   TIDYBREED_BENCH_LARGE=1 Rscript dev/benchmarks/profile_haplotype_paths.R
#
# This is the Phase-0 deliverable of plans/optimize_add_offspring.md: a
# per-STEP time + peak-memory breakdown of define_founder_haplotypes(),
# add_founders(), and add_offspring() across realistic mating designs, so
# Stages 1-3 are ordered by measured ROI rather than assumption. It is
# distinct from dev/benchmarks/benchmark_haplotype_scale.R (broad insert
# throughput) -- this one attributes wall time to named internal steps.
#
# Method: base-R Rprof(line.profiling=TRUE) attributes self-time to source
# lines, which are bucketed into named steps below. peak memory = max Vcells
# used since gc(reset=TRUE) (a true peak, not just total allocations).
# profvis (flame-graph HTML) and bench (mem_alloc) are OPTIONAL enrichment --
# requireNamespace()-guarded and skipped with a message if absent; the core
# table needs neither. Neither is a package Suggests (dev-only).
#
# NOT part of R CMD check or the testthat suite.
# ------------------------------------------------------------------------

# Line profiling of package internals needs source refs retained at load.
options(keep.source = TRUE, keep.source.pkgs = TRUE)
suppressMessages(pkgload::load_all(".", quiet = TRUE))

has_profvis <- requireNamespace("profvis", quietly = TRUE)
has_bench   <- requireNamespace("bench",   quietly = TRUE)
if (!has_profvis) message("[skip] profvis not installed -- no flame-graph HTML.")
if (!has_bench)   message("[skip] bench not installed -- no mem_alloc column.")

run_large <- identical(Sys.getenv("TIDYBREED_BENCH_LARGE"), "1")

# --- Named-step buckets: file -> list(c(from_line, to_line, step_label)) ------
# Line ranges are for the CURRENT (pre-Stage-0) working tree. add_offspring.R
# and recombination_helpers.R both feed the "gamete-generation" step so they
# aggregate under one label. The two founder-pool files are one step each.
BUCKETS <- list(
  "add_offspring.R" = list(
    c(91,  211, "validate + parent-meta"),
    c(213, 294, "map-resolve + cache"),
    c(296, 309, "parent-load SQL"),
    c(310, 377, "parent-matrix build (O(P*N) scan)"),
    c(379, 406, "offspring-ids"),
    c(408, 502, "gamete-generation"),
    c(504, 522, "wide-build"),
    c(524, 553, "ind_meta-build"),
    c(555, 603, "DB-write (UNPIVOT+INSERT)")
  ),
  "recombination_helpers.R" = list(
    c(1, 200, "gamete-generation")
  ),
  "add_founders.R" = list(
    c(95,  149, "validate"),
    c(151, 180, "read-pool + rebuild-matrix"),
    c(182, 191, "id-seq"),
    c(193, 205, "sample()"),
    c(207, 232, "ind_meta-build"),
    c(234, 257, "wide-build"),
    c(259, 308, "DB-write (per-chr UNPIVOT+INSERT)")
  ),
  "define_founder_haplotypes.R" = list(
    c(1, 5000, "pool-generation (wide matrix)")
  ),
  "founder_haplotype_helpers.R" = list(
    c(1, 5000, "pool long-materialize + write")
  )
)

bucket_of <- function(file, line) {
  specs <- BUCKETS[[file]]
  if (is.null(specs)) return(NA_character_)
  for (s in specs) if (line >= as.integer(s[1]) && line <= as.integer(s[2])) return(s[3])
  NA_character_
}

# Profile a single thunk: returns elapsed, peak MB, and a per-step self-time
# table derived from Rprof line profiling.
profile_call <- function(thunk) {
  tmp <- tempfile(fileext = ".Rprof")
  gc(reset = TRUE)
  t0 <- proc.time()[["elapsed"]]
  Rprof(tmp, line.profiling = TRUE, interval = 0.004, memory.profiling = FALSE)
  invisible(thunk())
  Rprof(NULL)
  elapsed <- proc.time()[["elapsed"]] - t0
  g <- gc()
  peak_mb <- sum(g[, "max used"] * c(NA, 8)[seq_len(nrow(g))], na.rm = TRUE) # fallback
  # Vcells row is the meaningful one for numeric/integer data:
  peak_mb <- unname(g["Vcells", "max used"]) * 8 / 1e6

  sr <- tryCatch(summaryRprof(tmp, lines = "show"), error = function(e) NULL)
  steps <- data.frame(step = character(0), self_s = numeric(0))
  if (!is.null(sr) && !is.null(sr$by.line) && nrow(sr$by.line) > 0) {
    bl <- sr$by.line
    keys <- rownames(bl)
    parts <- strsplit(keys, "#", fixed = TRUE)
    file <- vapply(parts, function(p) p[1], character(1))
    line <- suppressWarnings(as.integer(vapply(parts,
                     function(p) if (length(p) > 1) p[2] else NA, character(1))))
    lbl <- mapply(bucket_of, file, line)
    self <- bl[["self.time"]]
    keep <- !is.na(lbl)
    if (any(keep)) {
      steps <- aggregate(self[keep], by = list(step = lbl[keep]), FUN = sum)
      names(steps)[2] <- "self_s"
      steps <- steps[order(-steps$self_s), , drop = FALSE]
    }
    other <- sum(self[!keep])
    if (other > 0) steps <- rbind(steps, data.frame(step = "(other/unbucketed)", self_s = other))
  }
  unlink(tmp)
  list(elapsed = elapsed, peak_mb = peak_mb, steps = steps)
}

print_result <- function(title, res) {
  cat(sprintf("\n-- %s --  elapsed %.3fs | peak Vcells %.0f MB\n",
              title, res$elapsed, res$peak_mb))
  if (nrow(res$steps) == 0) {
    cat("   (no line-profile samples -- run longer / larger, or srcrefs missing)\n")
    return(invisible())
  }
  tot <- sum(res$steps$self_s)
  for (i in seq_len(nrow(res$steps))) {
    s <- res$steps[i, ]
    cat(sprintf("   %-38s %7.3fs  %5.1f%%\n", s$step, s$self_s,
                100 * s$self_s / max(tot, 1e-9)))
  }
}

# ---- scenario configs -------------------------------------------------------
scales <- list(small = list(n_loci = 5000L, n_chr = 10L, n_hap = 600L,
                            n_m = 250L, n_f = 250L, n_off = 500L))
if (run_large) {
  scales$large <- list(n_loci = 50000L, n_chr = 30L, n_hap = 5000L,
                       n_m = 2500L, n_f = 2500L, n_off = 10000L)
}

profile_scale <- function(nm, cfg) {
  cat("\n========================================================\n")
  cat(sprintf("SCALE %s: n_loci=%d n_chr=%d founders=%d offspring=%d\n",
              nm, cfg$n_loci, cfg$n_chr, cfg$n_m + cfg$n_f, cfg$n_off))
  cat("========================================================\n")

  pop <- open_pop(pop_name = paste0("prof_", nm), db_name = ":memory:") |>
    define_genome(n_loci = cfg$n_loci, n_chr = cfg$n_chr, chr_len_Mb = 100)

  # -- FOUNDER PATH -----------------------------------------------------------
  set.seed(1)
  r_defA <- profile_call(function() {
    pop <<- define_founder_haplotypes(pop, n_haplotypes = cfg$n_hap,
                                      line_name = "A", method = "uniform")
  })
  print_result("define_founder_haplotypes() line A", r_defA)

  set.seed(2)
  pop <- define_founder_haplotypes(pop, n_haplotypes = cfg$n_hap,
                                   line_name = "B", method = "uniform")

  set.seed(11)
  r_fndA <- profile_call(function() {
    pop <<- pop |> get_table("founder_haplotypes") |>
      dplyr::filter(line_name == "A") |>
      add_founders(n_males = cfg$n_m, n_females = cfg$n_f, line_name = "A", gen = 0L)
  })
  print_result("add_founders() line A", r_fndA)

  set.seed(12)
  pop <- pop |> get_table("founder_haplotypes") |>
    dplyr::filter(line_name == "B") |>
    add_founders(n_males = cfg$n_m, n_females = cfg$n_f, line_name = "B", gen = 0L)

  # -- OFFSPRING PATH: design 1 -- many parents, 1 offspring each (litter/hatch)
  sires <- paste0("A_", seq_len(cfg$n_m))
  dams  <- paste0("B_", cfg$n_m + seq_len(cfg$n_f))
  npair <- min(length(sires), length(dams))
  matings_litter <- tibble::tibble(
    id_parent_1 = rep(sires[seq_len(npair)], length.out = cfg$n_off),
    id_parent_2 = rep(dams[seq_len(npair)],  length.out = cfg$n_off),
    sex = rep(c("M", "F"), length.out = cfg$n_off), line_name = "F1", gen = 1L)
  set.seed(21)
  r_off1 <- profile_call(function() pop <<- add_offspring(pop, matings_litter))
  print_result(sprintf("add_offspring() litter (%d parents x %d off)",
                       npair, cfg$n_off), r_off1)

  # -- OFFSPRING PATH: design 2 -- 1 sire x many offspring (per-parent-loop stress)
  matings_1xN <- tibble::tibble(
    id_parent_1 = sires[1], id_parent_2 = rep(dams, length.out = cfg$n_off),
    sex = rep(c("M", "F"), length.out = cfg$n_off), line_name = "S1", gen = 1L)
  set.seed(22)
  r_off2 <- profile_call(function() pop <<- add_offspring(pop, matings_1xN))
  print_result(sprintf("add_offspring() 1-sire x %d off", cfg$n_off), r_off2)

  if (has_bench) {
    mm <- bench::mark(add_offspring(pop, matings_litter[1:min(200, cfg$n_off), ]),
                      iterations = 1, check = FALSE, filter_gc = FALSE)
    cat(sprintf("\n   [bench] add_offspring(200) mem_alloc = %s\n",
                format(mm$mem_alloc)))
  }
  close_pop(pop)
  invisible(NULL)
}

for (nm in names(scales)) profile_scale(nm, scales[[nm]])

# -- Special-chromosome scenario (exercises the R special branch) -------------
cat("\n========================================================\n")
cat("SPECIAL-CHR scenario (define_chr X/Y): small scale\n")
cat("========================================================\n")
pop <- open_pop(pop_name = "prof_sexchr", db_name = ":memory:") |>
  define_genome(n_loci = 2000L, n_chr = 4L, chr_len_Mb = 100) |>
  define_chr("4", copy_mode_M = "half", copy_mode_F = "full",
             hemi_parent = "parent_2", recombines = TRUE)
set.seed(31); pop <- define_founder_haplotypes(pop, n_haplotypes = 200L,
                                               line_name = "A", method = "uniform")
set.seed(32)
pop <- pop |> get_table("founder_haplotypes") |>
  add_founders(n_males = 50, n_females = 50, line_name = "A", gen = 0L)
matings_sx <- tibble::tibble(
  id_parent_1 = paste0("A_", 1:50),
  id_parent_2 = paste0("A_", 51:100),
  sex = rep(c("M", "F"), 25), line_name = "F1", gen = 1L)
set.seed(33)
r_sx <- profile_call(function() pop <<- add_offspring(pop, matings_sx))
print_result("add_offspring() with sex chromosome (chr 4 half/full)", r_sx)
close_pop(pop)

cat("\nDone. See the committed baseline table + conclusion at the bottom of this file.\n")
if (!run_large) cat("(Set TIDYBREED_BENCH_LARGE=1 for the 50k-loci / 10k-offspring scale.)\n")

# ============================================================================
# COMMITTED BASELINE  (pre-Stage-0 tree; small scale; macOS, R 4.5.3, in-memory
# DuckDB. Absolute seconds are machine-specific -- read the % split, not the s.)
# ============================================================================
# small = 5000 loci, 10 chr, 600-hap pool, 250M+250F founders/line, 500 offspring.
#
# FOUNDER PATH
#   define_founder_haplotypes() A   0.80s  peak 106 MB
#     pool long-materialize + write ............ 95.8%   <- .write_founder_haplotypes (rep()+dbWriteTable)
#     pool-generation (wide matrix) ............  1.8%
#   add_founders() A               14.23s  peak 301 MB
#     DB-write (per-chr UNPIVOT+INSERT) ........ 68.8%   <- wide->UNPIVOT (Stage 1)
#     read-pool + rebuild-matrix ............... 27.9%   <- dense hap_data_matrix rebuild from long pool
#
# OFFSPRING PATH
#   add_offspring() litter (250 parents x 500)  40.05s  peak 402 MB
#     parent-matrix build (O(P*N) scan) ........ 58.6%   <- hot spot #2; Stage 0 split() fix
#     DB-write (UNPIVOT+INSERT) ................. 34.3%   <- wide->UNPIVOT (Stage 1)
#     parent-load SQL ...........................  4.8%
#     gamete-generation .........................  1.7%   <- the C++-kernel target (Stages 2-3)
#   add_offspring() 1-sire x 500                18.70s  peak 289 MB
#     DB-write (UNPIVOT+INSERT) ................. 60.5%
#     parent-matrix build (O(P*N) scan) ........ 29.8%
#     gamete-generation .........................  3.2%
#   add_offspring() + sex chr (100 off)          1.33s  peak 288 MB
#     DB-write 50.4% | parent-matrix build 35.1% | gamete-generation 1.9%
#
#   bench: add_offspring(200 offspring) mem_alloc = 43.2 GB  <- dense wide-matrix
#          allocation churn; the memory case for Stage 1's long/batched write.
#
# CONCLUSION (orders the next increment):
#   1. Gamete generation is 1.7-3.2% of add_offspring -- NOT the bottleneck at
#      this scale. The dqrng/C++/parallel gamete kernel (Stages 2-4) optimizes a
#      <5% slice. It should be RE-PRIORITIZED / deferred behind the two wins below;
#      the RNG-parity spike is done so Stage 2/3 remain unblocked when wanted.
#   2. The O(P*N) parent-matrix scan is the #1 add_offspring hot spot in
#      many-parent designs (58.6%). Stage 0's split() fix (C2) targets it directly
#      and should collapse it to a single O(N) pass.
#   3. The wide->UNPIVOT DB write dominates everywhere else (34-69%: add_founders
#      68.8%, add_offspring 34-60%, founder pool write 95.8%) and drives 43 GB of
#      allocation churn for 200 offspring. Stage 1 (direct long write + batching,
#      shared across all three paths) is the highest-ROI change after C2.
#   4. Net re-order: Stage 0 (split + preallocate) -> Stage 1 (long/batched write,
#      all three paths) capture the large majority of the win; the gamete kernel
#      (2-4) is a later, smaller optimization, not the priority the design implied.
#
# large (50000 loci, 30 chr, 5000-hap, 2500+2500, 10000 off) is gated behind
# TIDYBREED_BENCH_LARGE=1. Not run for this baseline: at target scale the dense
# wide matrices (~2 GB each x4) are expected to thrash/OOM on the UNMODIFIED tree
# -- itself the point of Stage 1. Re-run large AFTER Stage 1 to confirm the ceiling.
#
# ----------------------------------------------------------------------------
# AFTER STAGE 0 (seam extraction + index-split parent scan + preallocated
# special_rows). Same small config; absolute s vary with machine load, so the
# load-independent metrics are the % split and mem_alloc.
#   add_offspring() litter (250 parents x 500):
#     parent-matrix build (O(P*N) scan)   58.6% (19.06s)  ->  ~10-11% (1.5-2.2s)   ~12x on that step
#     total                               40.05s          ->  12.5s (clean run)    ~3.2x
#     DB-write is now the dominant step (~70%) -> the Stage 1 target.
#   add_offspring(200) mem_alloc          43.2 GB         ->  1.55 GB              ~28x less churn
#   peak Vcells (litter)                  402 MB          ->  418 MB (index-split; a
#     transient data.frame-split variant hit 659 MB and was reverted to splitting
#     ROW INDICES so only one parent's subset materializes at a time).
#   Output-neutral: test-parity.R golden (autosomes) green; make_gametes_batch()
#     proven == inline make_gamete() for the same seed (special path); full suite
#     1229 pass / 0 fail. gamete-generation stays ~2-6% (unchanged; the C++-kernel
#     slice), confirming Stage 0 is a pure refactor + the two structural wins.
