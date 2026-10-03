# group_sum() / group_mean() contributors must be bit-identical whatever
# DuckDB's thread count (CLAUDE.md, "Identical means bit-identical"). The mate
# sum in .group_mate_tbv() was a plain floating SUM() over the mates, which a
# parallel aggregate may add in any order once a group has three or more
# mates. It now accumulates exactly through GEV_ACC_TYPE, as the genome-effect
# evaluator does (tidybreed 0.71.2, plans/import_qtl_effect_methods.md step 0b,
# B-1). Only expect_identical() can catch a regression here.
#
# Size matters: DuckDB only splits the aggregate across threads when the join
# is large enough, and whether it does also depends on the query plan. With
# pens of 10 the old plain SUM() never diverged. On this fixture (two pens of
# 200) the end-to-end phenotype test below fails on the old code; the direct
# helper comparison happened not to, and stays for the sum/mean semantics.
# Verified 2026-10-02 by reverting the fix.

grp_pop <- function(name) {
  set.seed(2101)
  pop <- open_pop(pop_name = name, db_name = ":memory:") |>
    define_genome(n_loci = 100, n_chr = 1, chr_len_Mb = 100) |>
    define_founder_haplotypes(n_haplotypes = 80, method = "fixed")
  pop <- pop |> get_table("founder_haplotypes") |>
    add_founders(n_males = 200L, n_females = 200L, line_name = "A")
  ind_ids <- sort(dplyr::collect(get_table(pop, "ind_meta"))$id_ind)
  pens <- tibble::tibble(id_ind = ind_ids,
                         pen_id = paste0("pen", rep(1:2, each = 200)))
  pop <- get_table(pop, "ind_meta") |> mutate_table(pen_id = pens)
  pop <- define_trait(pop, "D")
  pop <- define_trait(pop, "S")
  pop <- suppressMessages(pop |> get_table("genome_meta") |>
    define_additive_effects(c("D", "S"),
      G = matrix(c(1, 0.1, 0.1, 0.3), 2, 2,
                 dimnames = list(c("D", "S"), c("D", "S")))))
  suppressMessages(pop |> get_table("ind_meta") |> add_tbv(c("D", "S")))
}

test_that("group-mate sums and means are bit-identical across thread counts", {
  pop <- grp_pop("grp_det_helper")
  on.exit(close_pop(pop), add = TRUE)
  conn <- pop$db_conn
  ids  <- sort(dplyr::collect(get_table(pop, "ind_meta"))$id_ind)

  at_threads <- function(n, aggregation) {
    DBI::dbExecute(conn, paste0("SET threads = ", n))
    .group_mate_tbv(conn, "S", ids, "pen_id", "ind_meta", aggregation,
                    what = "test")
  }
  for (agg in c("sum", "mean")) {
    one <- at_threads(1L, agg)
    expect_length(one, 400L)
    for (i in 1:5) expect_identical(at_threads(8L, agg), one)
    expect_identical(at_threads(1L, agg), one)
  }

  # Semantics unchanged: the sum over the *other* pen members' TBVs.
  tbv <- dplyr::collect(get_table(pop, "ind_tbv")) |>
    dplyr::filter(trait_name == "S")
  pens <- dplyr::collect(get_table(pop, "ind_meta"))[, c("id_ind", "pen_id")]
  v    <- tbv$tbv_value[match(ids, tbv$id_ind)]
  pen  <- pens$pen_id[match(ids, pens$id_ind)]
  expected <- vapply(seq_along(ids), function(i)
    sum(v[pen == pen[[i]] & ids != ids[[i]]]), numeric(1))
  DBI::dbExecute(conn, "SET threads = 1")
  expect_equal(.group_mate_tbv(conn, "S", ids, "pen_id", "ind_meta", "sum",
                               what = "test"), expected, tolerance = 1e-12)
  expect_equal(.group_mate_tbv(conn, "S", ids, "pen_id", "ind_meta", "mean",
                               what = "test"), expected / 199, tolerance = 1e-12)
})

test_that("a seeded SGE phenotype is bit-identical at 1 and 8 threads", {
  run <- function(name, threads) {
    pop <- grp_pop(name)
    on.exit(close_pop(pop), add = TRUE)
    DBI::dbExecute(pop$db_conn, paste0("SET threads = ", threads))
    pop <- define_phenotype(pop, "W", residual_var = 1,
                            formula_tbv = "D + group_sum(S, pen_id) + group_mean(D, pen_id)")
    set.seed(77)
    pop <- pop |> get_table("ind_meta") |> add_phenotype("W")
    DBI::dbGetQuery(pop$db_conn,
      "SELECT pheno_value FROM ind_phenotype ORDER BY id_ind")$pheno_value
  }
  one <- run("grp_det_t1", 1L)
  expect_length(one, 400L)
  expect_identical(run("grp_det_t8", 8L), one)
})

test_that("a phenotype_components group contributor is bit-identical at 1 and 8 threads", {
  # The second caller of .group_mate_tbv(): the components route, not the
  # formula DSL. Step 3 retargets the helper to ind_tgv_total, and both
  # routes must keep the exact sum. Same contributors as the formula test
  # above, which is what makes the old floating SUM() diverge here.
  run <- function(name, threads) {
    pop <- grp_pop(name)
    on.exit(close_pop(pop), add = TRUE)
    DBI::dbExecute(pop$db_conn, paste0("SET threads = ", threads))
    pop <- define_phenotype(pop, "W", residual_var = 1,
      components = tibble::tribble(
        ~source_trait_name, ~contributor_type, ~group_column, ~aggregation,
        "D",                "self",            NA_character_, "sum",
        "S",                "group",           "pen_id",      "sum",
        "D",                "group",           "pen_id",      "mean"))
    set.seed(77)
    pop <- pop |> get_table("ind_meta") |> add_phenotype("W")
    DBI::dbGetQuery(pop$db_conn,
      "SELECT pheno_value FROM ind_phenotype ORDER BY id_ind")$pheno_value
  }
  one <- run("grp_cmp_t1", 1L)
  expect_length(one, 400L)
  expect_identical(run("grp_cmp_t8", 8L), one)
})
