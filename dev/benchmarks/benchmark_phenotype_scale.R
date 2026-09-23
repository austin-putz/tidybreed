#!/usr/bin/env Rscript
# ------------------------------------------------------------------------
# Phase 8 of plans/sample_correlated_effects.md: the correlated-effects
# phenotype path must scale, and the two places it can stop scaling are
# named in the plan -- the observation-pattern query (Stage 2, reading back
# the residuals a block member has already realized) and the writes
# (Stage 3).
#
# add_phenotype() runs in three stages (?add_phenotype_stages):
#   PLAN     no RNG, no writes  -- one query per planning concern
#   RESOLVE  RNG, no writes     -- the stored-coordinate lookup and the
#                                 conditional-MVN draws
#   COMMIT   writes, no RNG     -- one transaction, register + INSERT
# The script times each stage separately, because a regression in one says
# something completely different from a regression in another.
#
# **Nothing here may change RNG semantics.** The draw order is part of the
# contract (tests/testthat/test-add_phenotype_*.R pin it), so an
# optimization is only admissible if the seeded output is byte-identical
# before and after. The `--check` mode below is the guard: it records the
# phenotype values for a fixed seed so a candidate optimization can be
# diffed against them.
#
# Four shapes, because they exercise different parts of Stage 2:
#   independent   two 1x1 residual blocks     -- no conditioning at all
#   correlated    one 2x2 residual block      -- the stored-coordinate join
#   conditional   2x2 block, strata by farm   -- per-record stratum select
#   named_effect  2x2 'pen' covariance block  -- the named-effect adapter
#
# The conditioning cost only appears on the *second* phenotype of a block,
# when the first one's residuals are already on disk, so every shape is
# timed twice: call 1 writes A, call 2 writes B conditional on stored A.
#
# Run manually:
#   Rscript dev/benchmarks/benchmark_phenotype_scale.R
#   TIDYBREED_BENCH_LARGE=1 Rscript dev/benchmarks/benchmark_phenotype_scale.R
#   Rscript dev/benchmarks/benchmark_phenotype_scale.R --check
# ------------------------------------------------------------------------

devtools::load_all(quiet = TRUE)
suppressMessages(library(dplyr))

run_large  <- identical(Sys.getenv("TIDYBREED_BENCH_LARGE"), "1")
check_mode <- "--check" %in% commandArgs(trailingOnly = TRUE)
sizes  <- if (check_mode) 200L else if (run_large) c(1000L, 4000L, 16000L) else c(250L, 1000L)
n_loci <- 1000L
n_qtl  <- 200L

sym <- function(nm, v) matrix(v, length(nm), length(nm), byrow = TRUE,
                              dimnames = list(nm, nm))

# ── Stage timing ─────────────────────────────────────────────────────────
# The three stages are internal, so they are wrapped in the namespace for
# the duration of the run. Wrapping rather than tracing keeps the return
# value untouched: these are pure timers, and the RNG never sees them.
.bench_time <- new.env(parent = emptyenv())

with_stage_timers <- function(expr) {
  ns  <- asNamespace("tidybreed")
  nms <- c(".ap_plan", ".ap_resolve", ".ap_commit")
  orig <- mget(nms, envir = ns)
  for (nm in nms) .bench_time[[nm]] <- 0
  on.exit({
    for (nm in nms) utils::assignInNamespace(nm, orig[[nm]], ns = "tidybreed")
  }, add = TRUE)
  for (nm in nms) {
    local({
      key <- nm
      fn  <- orig[[key]]
      utils::assignInNamespace(key, function(...) {
        t <- proc.time()
        on.exit(.bench_time[[key]] <-
                  .bench_time[[key]] + (proc.time() - t)[["elapsed"]], add = TRUE)
        fn(...)
      }, ns = "tidybreed")
    })
  }
  total <- proc.time()
  value <- force(expr)
  list(value   = value,
       total   = (proc.time() - total)[["elapsed"]],
       plan    = .bench_time[[".ap_plan"]],
       resolve = .bench_time[[".ap_resolve"]],
       commit  = .bench_time[[".ap_commit"]])
}

# ── Population ───────────────────────────────────────────────────────────
build_pop <- function(n_ind, shape) {
  # Seeded: the founder haplotypes and the QTL effects are sampled, so an
  # unseeded build would give every run a different genome and make the
  # --check comparison meaningless. CLAUDE.md requires deterministic
  # benchmarks for the same reason.
  set.seed(20260922L)
  pop <- open_pop(pop_name = paste0("bench_ph_", shape, "_", n_ind),
                  db_name = ":memory:") |>
    define_genome(n_loci = n_loci, n_chr = 10L, chr_len_Mb = 100) |>
    define_founder_haplotypes(n_haplotypes = 100L, method = "fixed",
                              allele_freq = 0.5) |>
    get_table("founder_haplotypes") |>
    add_founders(n_males = n_ind / 2, n_females = n_ind / 2, line_name = "A")

  qtl <- pop |> get_table("genome_meta") |> collect() |>
    arrange(locus_id) |> slice_head(n = n_qtl) |> pull(locus_name)
  for (t in c("A", "B")) {
    pop <- define_trait(pop, t, target_add_var = 1.0)
    pop <- suppressWarnings(
      pop |> get_table("genome_meta") |> filter(locus_name %in% !!qtl) |>
        define_additive_effects(t, effects = rep(1.0, n_qtl)))
  }

  # A condition column and a pen column, both with a realistic number of
  # levels: the stratum select and the named-effect adapter both scale with
  # level count, not with individuals.
  ids <- pop |> get_table("ind_meta") |> collect() |> arrange(id_ind) |>
    pull(id_ind)
  farms <- paste0("F", seq_len(4))
  pens  <- paste0("P", seq_len(max(2L, n_ind %/% 20L)))
  for (f in farms) {
    who <- ids[seq_along(ids) %% length(farms) == match(f, farms) - 1L]
    if (length(who)) pop <- pop |> get_table("ind_meta") |>
      filter(id_ind %in% !!who) |> mutate_table(farm = f)
  }
  for (p in pens) {
    who <- ids[seq_along(ids) %% length(pens) == match(p, pens) - 1L]
    if (length(who)) pop <- pop |> get_table("ind_meta") |>
      filter(id_ind %in% !!who) |> mutate_table(pen = p)
  }

  if (shape == "independent") {
    for (t in c("A", "B")) pop <- define_phenotype(
      pop, t, type = "continuous", mean = 10, residual_var = 1)
  } else {
    for (t in c("A", "B")) pop <- define_phenotype(
      pop, t, type = "continuous", mean = 10)
  }

  if (shape == "correlated" || shape == "named_effect") {
    pop <- define_residual_cov(pop, c("A", "B"), sym(c("A", "B"), c(2, .8, .8, 1)))
  }
  if (shape == "conditional") {
    for (f in farms) pop <- define_residual_cov(
      pop, c("A", "B"), sym(c("A", "B"), c(2, .8, .8, 1)),
      condition_column = "farm", condition_level = f)
  }
  if (shape == "named_effect") {
    for (t in c("A", "B")) pop <- define_effect_random(
      pop, t, "pen", source_column = "pen", variance = 1)
    pop <- define_effect_cov_matrix(pop, "pen", sym(c("A", "B"), c(1, .5, .5, 1)))
  }
  pop
}

# ── Run ──────────────────────────────────────────────────────────────────
shapes <- c("independent", "correlated", "conditional", "named_effect")

if (check_mode) {
  # RNG-semantics guard. Any optimization to Stage 2 or Stage 3 must leave
  # these numbers byte-identical; diff two runs of this mode across the
  # change. Values, not a hash, so a diff says *where* it moved.
  cat("tidybreed phenotype RNG-semantics check (seeds 42/43, n = ", sizes,
      ")\n\n", sep = "")
  for (shape in shapes) {
    pop <- suppressMessages(suppressWarnings(build_pop(sizes, shape)))
    suppressMessages(suppressWarnings(
      pop |> get_table("ind_meta") |> add_phenotype("A", seed = 42)))
    suppressMessages(suppressWarnings(
      pop |> get_table("ind_meta") |> add_phenotype("B", seed = 43)))
    v <- DBI::dbGetQuery(pop$db_conn,
      "SELECT phenotype_name, id_ind, pheno_value, residual_value
       FROM ind_phenotype ORDER BY phenotype_name, id_ind")
    cat("== ", shape, " ==\n", sep = "")
    print(utils::head(v, 6), digits = 17)
    cat("  checksum(pheno_value)   ", format(sum(v$pheno_value), digits = 17), "\n")
    cat("  checksum(residual_value)",
        format(sum(v$residual_value, na.rm = TRUE), digits = 17), "\n\n")
    close_pop(pop)
  }
  quit(save = "no")
}

cat("tidybreed add_phenotype() scale benchmark\n")
cat("loci:", n_loci, " QTL:", n_qtl, " sizes:", paste(sizes, collapse = ", "),
    "\n\n")

for (shape in shapes) {
  cat("== shape:", shape, "==\n")
  for (n_ind in sizes) {
    pop <- suppressMessages(suppressWarnings(build_pop(n_ind, shape)))
    # TBVs first: add_phenotype() would otherwise materialize them inside
    # Stage 1 and the evaluator's cost would be charged to planning.
    suppressMessages(pop |> get_table("ind_meta") |> add_tbv(c("A", "B")))

    for (pass in 1:2) {
      t <- with_stage_timers(suppressMessages(suppressWarnings(
        pop |> get_table("ind_meta") |>
          add_phenotype(c("A", "B")[pass], seed = 40 + pass))))
      cat(sprintf(
        "  n = %6d  call %d (%s)  total %7.3fs (%6.3f ms/ind)  plan %6.3f  resolve %6.3f  commit %6.3f\n",
        n_ind, pass, c("no stored residual", "conditional on stored")[pass],
        t$total, 1000 * t$total / n_ind, t$plan, t$resolve, t$commit))
    }
    close_pop(pop)
  }
  cat("\n")
}

cat("Reading the output:\n")
cat("  * ms/ind should be flat across sizes in every stage.\n")
cat("  * call 2 - call 1 in `resolve` is the observation-pattern query: the\n")
cat("    cost of reading back the residuals already realized for the block.\n")
cat("    It should be a constant number of statements, not one per entity.\n")
cat("  * `commit` is one transaction; a rise with n that outpaces the row\n")
cat("    count means the INSERT stopped being batched.\n")
cat("  * `independent` is the floor -- no conditioning, no stratum select.\n")
cat("    The gap to `correlated` is what the block machinery costs.\n")
