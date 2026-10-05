# Part A gates (plans/import_qtl_effect_methods.md §11, A0-A22) at the
# generator level: define_additive_effects() calibrates its architecture by the
# congruence and resolves its target from trait_var_comp (§6C). The gates on
# the internals (A1, A2's oracle, A4, A5) are in test-qtl-congruence.R.

# Founders from a uniform pool: many loci, HWE and (near) linkage equilibrium.
anchor_pop <- function(name, n_loci = 60, n_hap = 400, n_ind = 400,
                       db_name = ":memory:") {
  set.seed(7001)
  pop <- open_pop(pop_name = name, db_name = db_name) |>
    define_genome(n_loci = n_loci, n_chr = 2, chr_len_Mb = 100) |>
    define_founder_haplotypes(n_haplotypes = n_hap)
  pop <- pop |> get_table("founder_haplotypes") |>
    add_founders(n_males = n_ind / 2, n_females = n_ind / 2, line_name = "A")
  pop <- get_table(pop, "ind_meta") |> mutate_table(gen = 0L)
  suppressMessages(pop |> define_trait("T1") |> define_trait("T2"))
}

# Stored generated effects as an n_qtl x k matrix with the centre p, locus order.
stored_B <- function(pop, traits) {
  m <- DBI::dbGetQuery(pop$db_conn, paste0(
    "SELECT g.trait_name, m.locus_id, g.genome_value, m.center_value ",
    "FROM genome_effects g JOIN genome_effect_members m USING (id_genome_effect) ",
    "WHERE g.effect_owner = 'generated' ORDER BY m.locus_id"))
  loci <- sort(unique(m$locus_id))
  B <- matrix(0, length(loci), length(traits), dimnames = list(loci, traits))
  for (t in traits) {
    r <- m[m$trait_name == t, ]
    B[match(r$locus_id, loci), t] <- r$genome_value
  }
  list(B = B, locus_id = loci,
       p = m$center_value[match(loci, m$locus_id)])
}

ge_rows <- function(pop) {
  list(
    terms   = DBI::dbGetQuery(pop$db_conn,
      "SELECT * FROM genome_effects ORDER BY id_genome_effect"),
    members = DBI::dbGetQuery(pop$db_conn,
      "SELECT * FROM genome_effect_members ORDER BY id_genome_effect, locus_id"),
    origins = DBI::dbGetQuery(pop$db_conn,
      "SELECT * FROM genome_effect_member_origins ORDER BY id_genome_effect"))
}
tvc_rows <- function(pop) DBI::dbGetQuery(pop$db_conn,
  "SELECT * FROM trait_var_comp ORDER BY id_trait_var_comp")

G2 <- matrix(c(1, 0.6, 0.6, 2), 2, 2, dimnames = list(c("T1", "T2"), c("T1", "T2")))

gen0 <- function(pop) get_table(pop, "ind_meta") |> dplyr::filter(gen == 0L)

quiet <- function(expr) suppressWarnings(suppressMessages(expr))


test_that("A0: the same seed gives identical stored terms, both anchors, k = 1 and 2", {
  run <- function(anchor, traits) {
    pop <- anchor_pop("a0")
    on.exit(close_pop(pop))
    G <- if (length(traits) == 1L) 0.5 else G2
    quiet(get_table(pop, "genome_meta") |> dplyr::filter(locus_id <= 40) |>
      define_additive_effects(traits, G = G, anchor = anchor,
        base_tbl = if (anchor == "realised") gen0(pop), seed = 11))
    ge_rows(pop)
  }
  for (anchor in c("genic", "realised")) for (traits in list("T1", c("T1", "T2"))) {
    expect_identical(run(anchor, traits), run(anchor, traits))
  }
})

test_that("A2: k = 2 genic delivers G exactly, off-diagonals included", {
  pop <- anchor_pop("a2")
  on.exit(close_pop(pop))
  quiet(get_table(pop, "genome_meta") |> dplyr::filter(locus_id <= 40) |>
    define_additive_effects(c("T1", "T2"), G = G2, seed = 3))
  s <- stored_B(pop, c("T1", "T2"))
  D <- crossprod(s$B * sqrt(2 * s$p * (1 - s$p)))
  expect_lt(max(abs(D - G2)), 1e-10)
})

test_that("A3: 'realised' delivers Cov(X B) = G on the base individuals", {
  pop <- anchor_pop("a3")
  on.exit(close_pop(pop))
  quiet(get_table(pop, "genome_meta") |> dplyr::filter(locus_id <= 40) |>
    define_additive_effects(c("T1", "T2"), G = G2, anchor = "realised",
                            base_tbl = gen0(pop), seed = 4))
  s <- stored_B(pop, c("T1", "T2"))
  X <- .dae_collect_dosages(pop, gen0(pop), s$locus_id)$X
  expect_lt(max(abs(stats::cov(X %*% s$B) - G2)), 1e-10)

  # k = 1 too, and the delivered variance is the TBV variance of the base.
  quiet(get_table(pop, "trait_var_comp") |> remove_rows(confirm_all = TRUE))
  quiet(get_table(pop, "genome_meta") |> dplyr::filter(locus_id <= 40) |>
    define_additive_effects("T1", G = 0.7, anchor = "realised",
                            base_tbl = gen0(pop), seed = 5))
  quiet(get_table(pop, "ind_meta") |> add_tgv("T1"))
  tbv <- tgv_additive(pop)
  expect_equal(stats::var(tbv$tgv_value), 0.7, tolerance = 1e-10)
})

test_that("A5: an anchor that cannot carry G is its own error, through the generator", {
  pop <- anchor_pop("a5")
  on.exit(close_pop(pop))
  expect_error(get_table(pop, "genome_meta") |> dplyr::filter(locus_id == 1L) |>
                 define_additive_effects(c("T1", "T2"), G = G2),
               "anchor cannot carry the target.*rank\\(G\\) = 2.*rank 1")
  expect_equal(nrow(tvc_rows(pop)), 0L)
})

test_that("A6: 'realised' refuses pools, copies, scoped calls and no base", {
  pop <- anchor_pop("a6")
  on.exit(close_pop(pop))
  gm <- get_table(pop, "genome_meta") |> dplyr::filter(locus_id <= 10)
  expect_error(define_additive_effects(gm, "T1", G = 1, anchor = "realised"),
               "needs `base_tbl`")
  expect_error(define_additive_effects(gm, "T1", G = 1, anchor = "realised",
                 base_tbl = get_table(pop, "founder_haplotypes")),
               "select individuals.*no individuals")
  expect_error(define_additive_effects(gm, "T1", G = 1, anchor = "realised",
                 base_tbl = get_table(pop, "ind_haplotype")),
               "select individuals.*partial genotypes")
  expect_error(define_additive_effects(gm, "T1", G = 1, anchor = "realised",
                 base_tbl = gen0(pop), line_name = "A"), "common scope only")
  expect_error(define_additive_effects(gm, "T1", G = 1, anchor = "realised",
                 base_tbl = gen0(pop), parent_origin = 1), "common scope only")
  expect_equal(nrow(tvc_rows(pop)), 0L)
  expect_equal(nrow(ge_rows(pop)$terms), 0L)
})

test_that("A7: 'union' keeps per-trait QTL sets and warns 'approximate'", {
  pop <- anchor_pop("a7")
  on.exit(close_pop(pop))
  # Per-trait sets, planted uncalibrated (test-only; the generator has no
  # unscaled mode since 0.74.1): T1 on 1-30, T2 on 20-50.
  set.seed(1)
  quiet(get_table(pop, "genome_meta") |> dplyr::filter(locus_id <= 30) |>
    plant_generated_additive("T1", effects = stats::rnorm(30)))
  quiet(get_table(pop, "genome_meta") |>
    dplyr::filter(locus_id >= 20, locus_id <= 50) |>
    plant_generated_additive("T2", effects = stats::rnorm(31)))
  w <- character()
  withCallingHandlers(
    suppressMessages(get_table(pop, "genome_meta") |> dplyr::filter(locus_id <= 50) |>
      define_additive_effects(c("T1", "T2"), G = G2, method = "union",
                              warn_bounds = NULL, seed = 3)),
    warning = function(cnd) { w <<- c(w, conditionMessage(cnd))
                              invokeRestart("muffleWarning") })
  expect_true(any(grepl("approximate.*correlation T1,T2 = ", w)))
  s <- stored_B(pop, c("T1", "T2"))
  expect_true(all(s$B[s$locus_id > 30, "T1"] == 0))
  expect_true(all(s$B[s$locus_id < 20, "T2"] == 0))
  # Each trait's own variance is still exact.
  D <- crossprod(s$B * sqrt(2 * s$p * (1 - s$p)))
  expect_equal(diag(D), diag(G2), tolerance = 1e-10, ignore_attr = TRUE)

  # The same target under "shared" is exact and raises no such warning.
  pop2 <- anchor_pop("a7b")
  on.exit(close_pop(pop2), add = TRUE)
  expect_no_warning(suppressMessages(
    get_table(pop2, "genome_meta") |> dplyr::filter(locus_id <= 50) |>
      define_additive_effects(c("T1", "T2"), G = G2, warn_bounds = NULL)))
})

test_that("A8: warn_bounds fires on an inbred base, not on HWE/LE data; a pool only messages", {
  pop <- anchor_pop("a8", n_loci = 10, n_hap = 4000, n_ind = 2000)
  on.exit(close_pop(pop))
  gm <- get_table(pop, "genome_meta")
  # HWE / LE individuals: observed Cov(X) is close to the genic limit.
  expect_no_warning(suppressMessages(
    define_additive_effects(gm, "T1", G = 1, base_tbl = gen0(pop), seed = 1)))
  # Fully inbred individuals: both copies identical, Cov(X) ~ 2 x genic.
  DBI::dbExecute(pop$db_conn, paste0(
    "UPDATE ind_haplotype AS h SET allele = s.allele FROM ind_haplotype AS s ",
    "WHERE h.parent_origin = 2 AND s.parent_origin = 1 ",
    "AND s.id_ind = h.id_ind AND s.locus_id = h.locus_id"))
  quiet(get_table(pop, "trait_var_comp") |> remove_rows(confirm_all = TRUE))
  expect_warning(suppressMessages(
    define_additive_effects(gm, "T1", G = 1, base_tbl = gen0(pop), seed = 1)),
    "base individuals \\(observed\\).*outside warn_bounds")
  quiet(get_table(pop, "trait_var_comp") |> remove_rows(confirm_all = TRUE))
  expect_no_warning(suppressMessages(
    define_additive_effects(gm, "T1", G = 1, base_tbl = gen0(pop), seed = 1,
                            warn_bounds = NULL)))
  # A founder-pool base is the pool expectation, reported as a message with
  # the realised-anchor hint, never a warning (Q22).
  pop2 <- anchor_pop("a8b", n_loci = 10, n_hap = 6, n_ind = 10)
  on.exit(close_pop(pop2), add = TRUE)
  msgs <- character()
  expect_no_warning(withCallingHandlers(
    define_additive_effects(get_table(pop2, "genome_meta"), "T1", G = 1,
                            warn_bounds = c(0.999, 1.001), seed = 1),
    message = function(cnd) {
      msgs <<- c(msgs, conditionMessage(cnd))
      invokeRestart("muffleMessage")
    }))
  pool_msg <- grep("pool expectation", msgs, value = TRUE)
  expect_length(pool_msg, 1L)
  expect_match(pool_msg, "sampling LD of 6 haplotypes")
  expect_match(pool_msg, 'anchor = "realised"', fixed = TRUE)
  # Inside the bounds the message stays, without the hint.
  quiet(get_table(pop2, "trait_var_comp") |> remove_rows(confirm_all = TRUE))
  msgs <- character()
  withCallingHandlers(
    define_additive_effects(get_table(pop2, "genome_meta"), "T1", G = 1,
                            warn_bounds = c(1e-6, 1e6), seed = 1),
    message = function(cnd) {
      msgs <<- c(msgs, conditionMessage(cnd))
      invokeRestart("muffleMessage")
    })
  pool_msg <- grep("pool expectation", msgs, value = TRUE)
  expect_length(pool_msg, 1L)
  expect_no_match(pool_msg, "realised")
})

test_that("A9: parent_origin uses n_eligible = 1 in the genic weights", {
  pop <- anchor_pop("a9")
  on.exit(close_pop(pop))
  quiet(get_table(pop, "genome_meta") |> dplyr::filter(locus_id <= 30) |>
    define_additive_effects("T1", G = 0.8, parent_origin = 1, seed = 9))
  s <- stored_B(pop, "T1")
  expect_equal(sum(s$p * (1 - s$p) * s$B[, 1]^2), 0.8, tolerance = 1e-12)
})

test_that("A10: target and terms commit together or not at all", {
  pop <- anchor_pop("a10")
  on.exit(close_pop(pop))
  before <- ge_rows(pop)
  # A failure inside the transaction, after the terms were inserted.
  local_mocked_bindings(.tvc_write_block = function(conn, ...) {
    DBI::dbExecute(conn, paste0(
      "INSERT INTO trait_var_comp VALUES (999, 'additive', NULL, 'T1', 'T1', 1)"))
    stop("boom")
  })
  expect_error(suppressMessages(get_table(pop, "genome_meta") |>
    dplyr::filter(locus_id <= 30) |>
    define_additive_effects(c("T1", "T2"), G = G2)), "boom")
  expect_equal(nrow(tvc_rows(pop)), 0L)
  expect_identical(ge_rows(pop), before)
})

test_that("A10: an infeasible rank writes neither target nor terms", {
  pop <- anchor_pop("a10b")
  on.exit(close_pop(pop))
  before <- ge_rows(pop)
  expect_error(get_table(pop, "genome_meta") |> dplyr::filter(locus_id == 1L) |>
                 define_additive_effects(c("T1", "T2"), G = G2), "rank")
  expect_equal(nrow(tvc_rows(pop)), 0L)
  expect_identical(ge_rows(pop), before)
})

test_that("A11: targets are validated before any write; names never relabelled", {
  pop <- anchor_pop("a11")
  on.exit(close_pop(pop))
  gm <- get_table(pop, "genome_meta") |> dplyr::filter(locus_id <= 30)
  bad <- matrix(c(1, 2, 2, 1), 2, 2)
  expect_error(define_additive_effects(gm, c("T1", "T2"), G = bad),
               "positive semi-definite; the smallest eigenvalue of its correlation matrix is -1")
  expect_error(define_effect_cov_matrix(pop, "additive", bad,
                                        trait_name = c("T1", "T2")),
               "positive semi-definite")
  # A stored indefinite block (written behind the writer's back) is refused too.
  DBI::dbExecute(pop$db_conn, paste0(
    "INSERT INTO trait_var_comp VALUES (1, 'additive', NULL, 'T1', 'T1', 1), ",
    "(2, 'additive', NULL, 'T1', 'T2', 2), (3, 'additive', NULL, 'T2', 'T1', 2), ",
    "(4, 'additive', NULL, 'T2', 'T2', 1)"))
  expect_error(define_additive_effects(gm, c("T1", "T2")), "positive semi-definite")
  DBI::dbExecute(pop$db_conn, "DELETE FROM trait_var_comp")
  expect_equal(nrow(ge_rows(pop)$terms), 0L)

  # Zero-rank target: all-zero effects. Rank-deficient but feasible: exact.
  quiet(define_additive_effects(gm, c("T1", "T2"), G = matrix(0, 2, 2), seed = 1))
  expect_true(all(stored_B(pop, c("T1", "T2"))$B == 0))
  quiet(get_table(pop, "trait_var_comp") |> remove_rows(confirm_all = TRUE))
  G1 <- matrix(c(1, 2, 2, 4), 2, 2)
  quiet(define_additive_effects(gm, c("T1", "T2"), G = G1, seed = 1))
  s <- stored_B(pop, c("T1", "T2"))
  expect_lt(max(abs(crossprod(s$B * sqrt(2 * s$p * (1 - s$p))) - G1)), 1e-10)

  # Future names are refused, writing nothing to either table.
  n_pvc <- nrow(DBI::dbGetQuery(pop$db_conn, "SELECT * FROM phenotype_var_comp"))
  expect_error(define_effect_cov_matrix(pop, "additive_by_dominance", 1,
                                        trait_name = "T1"), "not yet supported")
  expect_error(define_effect_cov_matrix(pop, "dominance_by_dominance", 1,
                                        trait_name = "T1"), "not yet supported")
  expect_equal(nrow(DBI::dbGetQuery(pop$db_conn,
    "SELECT * FROM phenotype_var_comp")), n_pvc)

  # Dimnames that disagree with trait_name are an error, never a relabel.
  named <- matrix(c(2, 0.5, 0.5, 1), 2, 2, dimnames = list(c("T2", "T1"), c("T2", "T1")))
  quiet(get_table(pop, "trait_var_comp") |> remove_rows(confirm_all = TRUE))
  expect_error(define_additive_effects(gm, c("T1", "T2"), G = named),
               "never relabelled")
  expect_error(define_effect_cov_matrix(pop, "additive", named,
                                        trait_name = c("T1", "T2")),
               "never relabelled")
})

test_that("A13: a line-scoped call says line effects add no mean difference", {
  pop <- anchor_pop("a13")
  on.exit(close_pop(pop))
  gm <- get_table(pop, "genome_meta") |> dplyr::filter(locus_id <= 20)
  msgs <- character()
  withCallingHandlers(suppressWarnings(
    define_additive_effects(gm, "T1", G = 1, line_name = "A")),
    message = function(cnd) { msgs <<- c(msgs, conditionMessage(cnd))
                              invokeRestart("muffleMessage") })
  expect_true(any(grepl("add no difference between line means", msgs)))
  msgs <- character()
  withCallingHandlers(suppressWarnings(
    define_additive_effects(gm, "T2", G = 1)),
    message = function(cnd) { msgs <<- c(msgs, conditionMessage(cnd))
                              invokeRestart("muffleMessage") })
  expect_false(any(grepl("add no difference between line means", msgs)))
})

test_that("A14: 'realised' after restore_pop() gives the same order and terms", {
  tmp <- tempfile(fileext = ".duckdb")
  on.exit(unlink(tmp), add = TRUE)
  pop <- anchor_pop("a14", db_name = tmp)
  quiet(get_table(pop, "genome_meta") |> dplyr::filter(locus_id <= 30) |>
    define_additive_effects("T1", G = 0.5, anchor = "realised",
                            base_tbl = gen0(pop), seed = 21))
  first <- stored_B(pop, "T1")
  close_pop(pop)
  pop <- quiet(restore_pop(tmp))
  on.exit(close_pop(pop), add = TRUE)
  quiet(get_table(pop, "genome_meta") |> dplyr::filter(locus_id <= 30) |>
    define_additive_effects("T2", G = 0.5, anchor = "realised",
                            base_tbl = gen0(pop), seed = 21))
  second <- stored_B(pop, "T2")
  expect_identical(unname(second$B[, "T2"]), unname(first$B[, "T1"]))
})

test_that("A15: a stored block is never overwritten; the error's remove_rows() call works", {
  pop <- anchor_pop("a15")
  on.exit(close_pop(pop))
  gm <- get_table(pop, "genome_meta") |> dplyr::filter(locus_id <= 30)
  quiet(define_additive_effects(gm, c("T1", "T2"), G = G2, seed = 1))
  err <- tryCatch(define_additive_effects(gm, c("T1", "T2"), G = G2),
                  error = function(e) conditionMessage(e))
  expect_match(err, "already stored for T1, T2")
  expect_match(err, "even by an identical")
  call <- sub("^.*remove it first:\n", "", err)
  pop <- eval(parse(text = call))
  expect_equal(nrow(tvc_rows(pop)), 0L)
  quiet(define_additive_effects(gm, c("T1", "T2"), G = G2, seed = 1))
  expect_identical(load_trait_cov(pop, "additive", c("T1", "T2")), G2)

  # A filtered table cannot smuggle a second matrix in.
  expect_error(define_additive_effects(gm, c("T1", "T2"), G = G2,
    trait_var_comp_tbl = get_table(pop, "trait_var_comp") |>
      dplyr::filter(effect_name == "none")), "not both")
  # define_effect_cov_matrix() refuses too: here first because the traits
  # have generated terms calibrated to the stored block (0.74.1, Q21).
  expect_error(define_effect_cov_matrix(pop, "additive", G2),
               "already have generated 'additive' terms")

  # A scalar G with one trait is a 1 x 1 block.
  suppressMessages(define_trait(pop, "T3"))
  quiet(define_additive_effects(gm, "T3", G = 0.25, seed = 1))
  expect_identical(get_trait_var(pop, "additive", "T3"), 0.25)
})

test_that("A16: trait_var_comp_tbl chooses the block, and refusals name the fixes", {
  pop <- anchor_pop("a16")
  on.exit(close_pop(pop))
  gm <- get_table(pop, "genome_meta") |> dplyr::filter(locus_id <= 30)
  quiet(define_effect_cov_matrix(pop, "additive", G2))
  tvc <- function() get_table(pop, "trait_var_comp")

  expect_error(define_additive_effects(gm, c("T1", "T2"),
                 trait_var_comp_tbl = get_table(pop, "ind_meta")),
               "must be get_table\\(pop, \"trait_var_comp\"\\)")
  expect_error(define_additive_effects(gm, c("T1", "T2"),
                 trait_var_comp_tbl = tvc() |>
                   dplyr::filter(!(trait_name_1 == "T2" & trait_name_2 == "T1"))),
               "complete k x k block: missing \\(T2, T1\\)")

  # Partial trait set: refused with both fixes, before any write or draw.
  seed_before <- .Random.seed
  expect_error(define_additive_effects(gm, "T1"),
               "links T1 with T2.*trait_name = c\\(\"T1\", \"T2\"\\).*trait_var_comp_tbl")
  expect_identical(.Random.seed, seed_before)
  expect_equal(nrow(ge_rows(pop)$terms), 0L)
  quiet(define_additive_effects(gm, "T1", seed = 1,
    trait_var_comp_tbl = tvc() |>
      dplyr::filter(trait_name_1 == "T1", trait_name_2 == "T1")))
  quiet(define_additive_effects(gm, c("T1", "T2"), seed = 1))

  # A stored dominance target is never silently ignored ...
  quiet(define_effect_cov_matrix(pop, "dominance", 0.1, trait_name = "T1"))
  expect_error(define_additive_effects(gm, "T1", seed = 1,
                 trait_var_comp_tbl = tvc() |>
                   dplyr::filter(trait_name_1 == "T1", trait_name_2 == "T1")),
               "stored 'dominance' target.*effect_name == \"additive\".*remove_rows")
  # ... and filtering it away says so explicitly.
  quiet(define_additive_effects(gm, "T1", seed = 1,
    trait_var_comp_tbl = tvc() |>
      dplyr::filter(effect_name == "additive", trait_name_1 == "T1",
                    trait_name_2 == "T1")))

  # Two candidate sets for one effect_name.
  quiet(define_effect_cov_matrix(pop, "additive", 0.3, trait_name = "T1",
                                 line_name = "C"))
  expect_error(define_additive_effects(gm, "T1",
                 trait_var_comp_tbl = tvc() |>
                   dplyr::filter(effect_name == "additive", trait_name_1 == "T1",
                                 trait_name_2 == "T1")),
               "two candidate 'additive' blocks")
})

test_that("A17: line targets resolve per line with fallback, and never mix", {
  pop <- anchor_pop("a17")
  on.exit(close_pop(pop))
  gm <- get_table(pop, "genome_meta") |> dplyr::filter(locus_id <= 30)
  quiet(define_effect_cov_matrix(pop, "additive", 1, trait_name = "T1"))
  # No line-C block: falls back to the population-wide one.
  expect_identical(get_trait_var(pop, "additive", "T1", line_name = "C"), 1)
  quiet(define_effect_cov_matrix(pop, "additive", 0.4, trait_name = "T1",
                                 line_name = "C"))
  expect_identical(get_trait_var(pop, "additive", "T1", line_name = "C"), 0.4)
  expect_identical(get_trait_var(pop, "additive", "T1"), 1)
  expect_identical(load_trait_cov(pop, "additive", "T1"),
                   matrix(1, 1, 1, dimnames = list("T1", "T1")))
  rows <- tvc_rows(pop)
  expect_setequal(rows$line_name, c(NA, "C"))

  # A passed G is written with the call's line_name.
  quiet(define_additive_effects(gm, "T2", G = 0.7, line_name = "D", seed = 1))
  expect_identical(get_trait_var(pop, "additive", "T2", line_name = "D"), 0.7)
  expect_true(is.na(get_trait_var(pop, "additive", "T2")))

  # Fallback is per effect_name: line C's own additive, the shared dominance.
  quiet(define_effect_cov_matrix(pop, "dominance", 0.1, trait_name = "T1"))
  r <- .dae_default_target_rows(pop$db_conn, "T1", "C")
  expect_identical(r$line_name[r$effect_name == "additive"], "C")
  expect_true(is.na(r$line_name[r$effect_name == "dominance"]))

  # line_name is for genetic effects only.
  expect_error(define_effect_cov_matrix(pop, "pen", 1, trait_name = "T1",
                                        line_name = "C"),
               "applies to genetic effects")
})

test_that("A18: targets are stored at full precision and hit to 1e-12", {
  pop <- anchor_pop("a18")
  on.exit(close_pop(pop))
  G <- matrix(c(1/3, 1/7, 1/7, 2/3), 2, 2, dimnames = list(c("T1", "T2"), c("T1", "T2")))
  quiet(define_effect_cov_matrix(pop, "additive", G))
  expect_identical(load_trait_cov(pop, "additive", c("T1", "T2")), G)
  quiet(get_table(pop, "genome_meta") |> dplyr::filter(locus_id <= 40) |>
    define_additive_effects(c("T1", "T2"), seed = 2))
  s <- stored_B(pop, c("T1", "T2"))
  expect_lt(max(abs(crossprod(s$B * sqrt(2 * s$p * (1 - s$p))) - G)), 1e-12)

  # A failure halfway through a target write leaves trait_var_comp unchanged.
  before <- tvc_rows(pop)
  local_mocked_bindings(.tvc_write_block = function(conn, ...) {
    DBI::dbExecute(conn, paste0(
      "INSERT INTO trait_var_comp VALUES (999, 'additive', NULL, 'X', 'X', 1)"))
    stop("halfway")
  })
  expect_error(define_effect_cov_matrix(pop, "additive", 1, trait_name = "X"),
               "halfway")
  expect_identical(tvc_rows(pop), before)
})

# A19 (G with manual or unscaled effects refused) is moot since 0.74.1: the
# generator has no manual or unscaled mode (Q21). test-define_additive_effects.R
# pins that `effects =` and `scale_to_target =` are unused arguments.

test_that("A20: reserved effect names are refused as input", {
  pop <- anchor_pop("a20")
  on.exit(close_pop(pop))
  for (nm in c("total", "unpartitioned", "between_components")) {
    expect_error(define_effect_cov_matrix(pop, nm, 1, trait_name = "T1"),
                 "reserved for derived output")
  }
  suppressMessages(define_phenotype(pop, "T1", residual_var = 1))
  expect_error(define_effect_random(pop, "T1", effect_name = "additive",
                                    source_column = "gen", variance = 1),
               "reserved genetic variance component")
  expect_error(define_effect_random(pop, "T1", effect_name = "dominance_by_dominance",
                                    source_column = "gen", variance = 1),
               "not yet supported")
  expect_equal(nrow(tvc_rows(pop)), 0L)
  expect_equal(nrow(DBI::dbGetQuery(pop$db_conn, paste0(
    "SELECT * FROM phenotype_var_comp WHERE effect_name <> 'residual'"))), 0L)
})

test_that("A21: the realised size guard errors before any write", {
  pop <- anchor_pop("a21")
  on.exit(close_pop(pop))
  local_mocked_bindings(QTL_REALISED_MAX_CELLS = 100)
  expect_error(get_table(pop, "genome_meta") |> dplyr::filter(locus_id <= 10) |>
                 define_additive_effects("T1", G = 1, anchor = "realised",
                                         base_tbl = gen0(pop)),
               "400 x 10 genotype matrix.*above the limit of 100.*anchor = \"genic\"")
  expect_equal(nrow(tvc_rows(pop)), 0L)
  expect_equal(nrow(ge_rows(pop)$terms), 0L)
})

test_that("A22: no refusal touches the RNG, with or without seed =", {
  pop <- anchor_pop("a22")
  on.exit(close_pop(pop))
  gm <- get_table(pop, "genome_meta") |> dplyr::filter(locus_id <= 10)
  quiet(define_effect_cov_matrix(pop, "additive", G2))
  set.seed(99)
  for (s in list(NULL, 5)) {
    before <- .Random.seed
    expect_error(define_additive_effects(gm, c("T1", "T2"),
                   G = matrix(c(1, 2, 2, 1), 2, 2), seed = s))          # bad target / overwrite
    expect_error(define_additive_effects(gm, "T1", seed = s))           # partial block
    expect_error(define_additive_effects(gm, c("T1", "T2"), G = G2, seed = s),
                 "already stored")                                      # overwrite
    expect_error(define_additive_effects(gm, c("T1", "T2"), anchor = "realised",
                                         seed = s))                     # no base
    expect_identical(.Random.seed, before)
  }
})

# ── Step 2 corrections (Codex review, 0.73.2) ──────────────────────────────

test_that("R1: a trait in small units keeps its variance; rank is unit-free", {
  # Unit level: diag(1, 1e-11) is full rank on the correlation scale.
  std <- .qtl_target_std(diag(c(1, 1e-11)))
  expect_identical(std$rank, 2L)
  cal <- .qtl_calibrate(diag(2), std, .qtl_anchor_diag(c(1, 1)))
  expect_equal(diag(cal$delivered), c(1, 1e-11), tolerance = 1e-12)
  # Rescaling the traits never changes the rank.
  G <- matrix(c(1, 0.5, 0.5, 1), 2)
  for (u in c(1e-12, 1, 1e12)) {
    expect_identical(.qtl_target_std(G * c(1, u) %o% c(1, u))$rank, 2L)
  }
  # Through the generator: the stored target is what is delivered.
  pop <- anchor_pop("r1")
  on.exit(close_pop(pop))
  Gs <- matrix(c(1, 0, 0, 1e-11), 2, dimnames = list(c("T1", "T2"), c("T1", "T2")))
  quiet(get_table(pop, "genome_meta") |> dplyr::filter(locus_id <= 40) |>
    define_additive_effects(c("T1", "T2"), G = Gs, seed = 2))
  s <- stored_B(pop, c("T1", "T2"))
  D <- crossprod(s$B * sqrt(2 * s$p * (1 - s$p)))
  expect_equal(D[2, 2], 1e-11, tolerance = 1e-8)
  expect_lt(abs(D[1, 2]) / sqrt(D[1, 1] * D[2, 2]), 1e-8)
})

test_that("R1: a calibration that misses the target errors before any write", {
  pop <- anchor_pop("r1b")
  on.exit(close_pop(pop))
  before <- ge_rows(pop)
  local_mocked_bindings(.qtl_congruence = function(B0, ...) list(B = B0))
  expect_error(suppressMessages(get_table(pop, "genome_meta") |>
    dplyr::filter(locus_id <= 30) |>
    define_additive_effects(c("T1", "T2"), G = G2, seed = 1)),
    "did not reach the target.*Nothing was written")
  expect_equal(nrow(tvc_rows(pop)), 0L)
  expect_identical(ge_rows(pop), before)
})

test_that("R2: the pool expectation uses the with-replacement divisor n_h", {
  # Pool {0, 1} at one QTL. A founder's dosage is Binomial(2, 1/2), variance
  # 1/2 = the genic target at p = 1/2, so the relative spectrum is exactly 1
  # (the n_h - 1 divisor reported 2).
  pop <- open_pop(pop_name = "r2", db_name = ":memory:") |>
    define_genome(n_loci = 2, n_chr = 1, chr_len_Mb = 10) |>
    define_founder_haplotypes(n_haplotypes = 2)
  on.exit(close_pop(pop))
  DBI::dbExecute(pop$db_conn, paste0(
    "UPDATE founder_haplotypes SET allele = CASE WHEN haplotype_id = ",
    "(SELECT MIN(haplotype_id) FROM founder_haplotypes) THEN 0 ELSE 1 END"))
  quiet(define_trait(pop, "T1"))
  msgs <- character()
  withCallingHandlers(
    get_table(pop, "genome_meta") |> dplyr::filter(locus_id == 1L) |>
      define_additive_effects("T1", G = 0.5, seed = 1),
    message = function(cnd) { msgs <<- c(msgs, conditionMessage(cnd))
                               invokeRestart("muffleMessage") })
  expect_match(grep("pool expectation", msgs, value = TRUE),
               "relative spectrum \\[1, 1\\].*sampling LD of 2 haplotypes")
  # Exact enumeration of the founder dosage distribution, effect b:
  b <- stored_B(pop, "T1")$B[1, 1]
  dos <- c(0, 1, 1, 2); expect_equal(mean((dos * b)^2) - mean(dos * b)^2, 0.5)
})

test_that("R3: a passed G never bypasses a stored non-additive target", {
  for (e in c("dominance", "additive_by_additive")) {
    pop <- anchor_pop("r3")
    quiet(define_effect_cov_matrix(pop, e, 0.2, trait_name = "T1"))
    before <- ge_rows(pop)
    expect_error(get_table(pop, "genome_meta") |> dplyr::filter(locus_id <= 30) |>
                   define_additive_effects("T1", G = 1, seed = 1),
                 paste0("stored '", e, "' target exists for T1.*",
                        "store `G` first and select it explicitly"))
    expect_equal(nrow(tvc_rows(pop)), 1L)
    expect_identical(ge_rows(pop), before)
    # The route the error gives works.
    quiet(define_effect_cov_matrix(pop, "additive", 1, trait_name = "T1"))
    quiet(get_table(pop, "genome_meta") |> dplyr::filter(locus_id <= 30) |>
      define_additive_effects("T1", seed = 1,
        trait_var_comp_tbl = get_table(pop, "trait_var_comp") |>
          dplyr::filter(effect_name == "additive", is.na(line_name))))
    expect_gt(nrow(ge_rows(pop)$terms), 0L)
    close_pop(pop)
  }
})

test_that("R5: union never calls a zero target covariance exact where QTL overlap", {
  pop <- anchor_pop("r5")
  on.exit(close_pop(pop))
  set.seed(1)
  quiet(get_table(pop, "genome_meta") |> dplyr::filter(locus_id <= 30) |>
    plant_generated_additive("T1", effects = stats::rnorm(30)))
  quiet(get_table(pop, "genome_meta") |> dplyr::filter(locus_id >= 20, locus_id <= 50) |>
    plant_generated_additive("T2", effects = stats::rnorm(31)))
  G0 <- diag(2); dimnames(G0) <- list(c("T1", "T2"), c("T1", "T2"))
  expect_warning(suppressMessages(
    get_table(pop, "genome_meta") |> dplyr::filter(locus_id <= 50) |>
      define_additive_effects(c("T1", "T2"), G = G0, method = "union",
                              warn_bounds = NULL, seed = 3)),
    "union.*approximate")
})

test_that("R5: union with no QTL for a zero-variance trait says so and stores the target", {
  pop <- anchor_pop("r5b")
  on.exit(close_pop(pop))
  set.seed(1)
  quiet(get_table(pop, "genome_meta") |> dplyr::filter(locus_id <= 30) |>
    plant_generated_additive("T1", effects = stats::rnorm(30)))
  Gz <- matrix(c(1, 0, 0, 0), 2, dimnames = list(c("T1", "T2"), c("T1", "T2")))
  expect_message(
    get_table(pop, "genome_meta") |> dplyr::filter(locus_id <= 50) |>
      define_additive_effects(c("T1", "T2"), G = Gz, method = "union",
                              warn_bounds = NULL, seed = 3),
    "Trait 'T2' has no QTL in this call; its target variance is 0")
  expect_equal(nrow(tvc_rows(pop)), 4L)
})

test_that("R7: a projected founder selection works; a diagnostic failure writes nothing", {
  pop <- anchor_pop("r7")
  on.exit(close_pop(pop))
  proj <- get_table(pop, "founder_haplotypes") |> dplyr::select(locus_name, allele)
  expect_message(
    get_table(pop, "genome_meta") |> dplyr::filter(locus_id <= 30) |>
      define_additive_effects("T1", G = 1, base_tbl = proj, seed = 1),
    "pool expectation")
  expect_equal(nrow(tvc_rows(pop)), 1L)

  pop2 <- anchor_pop("r7b")
  on.exit(close_pop(pop2), add = TRUE)
  before <- ge_rows(pop2)
  local_mocked_bindings(.dae_diagnostics = function(...) stop("diag boom"))
  expect_error(suppressMessages(get_table(pop2, "genome_meta") |>
    dplyr::filter(locus_id <= 30) |> define_additive_effects("T1", G = 1)),
    "diag boom")
  expect_equal(nrow(tvc_rows(pop2)), 0L)
  expect_identical(ge_rows(pop2), before)
})

test_that("R9: an anchor-rank refusal does not touch the RNG; seed or not", {
  pop <- anchor_pop("r9")
  on.exit(close_pop(pop))
  gm <- get_table(pop, "genome_meta") |> dplyr::filter(locus_id == 1L)
  set.seed(99)
  for (s in list(NULL, 5)) {
    before <- .Random.seed
    expect_error(define_additive_effects(gm, c("T1", "T2"), G = G2, seed = s),
                 "anchor cannot carry the target")
    expect_identical(.Random.seed, before)
  }
})

test_that("paper-12: genic is the random-mating limit, approached slowly under linkage", {
  # Source paper-12, through tidybreed's own transmission (add_offspring())
  # instead of the script's crossover code. Founders drawn from a pool in
  # strong LD (two backgrounds, 5% switches) start far from the genic target;
  # random mating moves the TBV variance toward it, but not all the way in
  # six generations on one short chromosome.
  set.seed(61)
  pop <- open_pop(pop_name = "p12", db_name = ":memory:") |>
    define_genome(n_loci = 20, n_chr = 1, chr_len_Mb = 20) |>
    define_founder_haplotypes(n_haplotypes = 400)
  on.exit(close_pop(pop))
  fh <- DBI::dbGetQuery(pop$db_conn, "SELECT * FROM founder_haplotypes")
  fh <- fh[order(fh$haplotype_id, fh$locus_name), ]
  bg <- sample(0:1, 400, TRUE)
  fh$allele <- as.integer((bg[match(fh$haplotype_id, sort(unique(fh$haplotype_id)))] +
                             stats::rbinom(nrow(fh), 1, 0.05)) %% 2)
  DBI::dbExecute(pop$db_conn, "DELETE FROM founder_haplotypes")
  DBI::dbAppendTable(pop$db_conn, "founder_haplotypes", fh)
  pop <- pop |> get_table("founder_haplotypes") |>
    add_founders(n_males = 200, n_females = 200, line_name = "A", gen = 0L)
  pop <- suppressMessages(define_trait(pop, "T"))
  pop <- suppressMessages(get_table(pop, "genome_meta") |>
    define_additive_effects("T", G = 4, seed = 62, warn_bounds = NULL))

  v <- numeric(7)
  for (g in 0:6) {
    pop <- suppressMessages(get_table(pop, "ind_meta") |>
      dplyr::filter(gen == !!g) |> add_tgv("T"))
    v[g + 1] <- stats::var(DBI::dbGetQuery(pop$db_conn, paste0(
      "SELECT t.tgv_value FROM ind_tgv t JOIN ind_meta i USING (id_ind) ",
      "WHERE i.gen = ", g))$tgv_value)
    if (g == 6) break
    par <- DBI::dbGetQuery(pop$db_conn, paste0(
      "SELECT id_ind, sex FROM ind_meta WHERE gen = ", g, " ORDER BY id_ind"))
    mat <- tibble::tibble(
      id_parent_1 = sample(par$id_ind[par$sex == "M"], 400, TRUE),
      id_parent_2 = sample(par$id_ind[par$sex == "F"], 400, TRUE),
      sex = rep(c("M", "F"), 200), line_name = "A", gen = g + 1L)
    pop <- suppressMessages(add_offspring(pop, mat))
  }
  expect_lt(abs(v[7] - 4), abs(v[1] - 4))
  expect_gt(abs(v[7] - 4), 0.05)
})
