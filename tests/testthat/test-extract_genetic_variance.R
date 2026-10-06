# extract_genetic_variance(): gates B1-B18 of plans/import_qtl_effect_methods.md
# §11 and plans/import_qtl_effect_methods_phase_4_plan.md 4b.6.
#
# Expectations are computed here, from the individuals' dosages and the
# coefficients the test wrote, or from the source project's own measuring
# instrument (helper-nonadd-oracle.R, decision D4). Neither shares code with
# the package's conversion or projection.

egv_pop <- function(loci = paste0("L", 1:6), n_ind = 80, seed = 1,
                    n_hap = 30, pop_name = "egv", db_name = ":memory:") {
  set.seed(seed)
  suppressMessages(
    open_pop(pop_name = pop_name, db_name = db_name) |>
      define_genome(n_loci = length(loci), n_chr = 1, chr_len_Mb = 10,
                    locus_names = loci) |>
      define_founder_haplotypes(n_haplotypes = n_hap) |>
      get_table("founder_haplotypes") |>
      add_founders(n_males = n_ind / 2, n_females = n_ind / 2,
                   line_name = "A"))
}

egv_write <- function(pop, trait, terms, owner = "custom", origin = NULL) {
  if (!trait %in% DBI::dbGetQuery(pop$db_conn,
                                  "SELECT trait_name FROM trait_meta")$trait_name) {
    pop <- suppressMessages(define_trait(pop, trait))
  }
  suppressMessages(define_genome_effect_terms(pop, trait, terms,
                                              effect_owner = owner,
                                              origin = origin))
}

egv_ids <- function(pop) {
  DBI::dbGetQuery(pop$db_conn, "SELECT id_ind FROM ind_meta ORDER BY id_ind")$id_ind
}

# Dosages, rows in `ids` order, columns in locus_id order.
egv_dosage <- function(pop, ids = egv_ids(pop)) {
  d <- DBI::dbGetQuery(pop$db_conn, paste0(
    "SELECT id_ind, locus_id, CAST(SUM(allele) AS DOUBLE) AS g ",
    "FROM ind_haplotype GROUP BY 1, 2 ORDER BY 1, 2"))
  loci <- sort(unique(d$locus_id))
  X <- matrix(NA_real_, length(ids), length(loci))
  X[cbind(match(d$id_ind, ids), match(d$locus_id, loci))] <- d$g
  X
}

# Evaluated totals via add_tgv() and the ind_tgv_total view.
egv_total <- function(pop, trait, ids = egv_ids(pop)) {
  suppressMessages(pop |> get_table("ind_meta") |> add_tgv(trait))
  v <- DBI::dbGetQuery(pop$db_conn, paste0(
    "SELECT id_ind, tgv_total FROM ind_tgv_total WHERE trait_name = '", trait,
    "'"))
  v$tgv_total[match(ids, v$id_ind)]
}

egv_get <- function(r, effect, t1, t2 = t1) {
  v <- r$cov_value[r$effect_name == effect & r$trait_name_1 == t1 &
                     r$trait_name_2 == t2]
  if (length(v) == 0L) NA_real_ else v
}

egv <- function(...) suppressMessages(extract_genetic_variance(...))

# Functional coefficients -> one writer frame (functional coding).
egv_functional_terms <- function(loci, a, d, pairs = NULL, e = NULL) {
  out <- list(suppressMessages(ad_terms(loci, a = a, d = d, p = 0.5)))
  if (!is.null(pairs)) {
    out[[2]] <- aa_terms(loci[pairs[, 1]], loci[pairs[, 2]], e = e,
                         p_1 = 0.5, p_2 = 0.5, report = FALSE)
  }
  do.call(rbind, out)
}

# Sum over every ordered pair of different blocks of Cov(b_t1, b'_t2).
egv_between <- function(V, t1, t2) {
  s <- 0
  for (b in names(V)) for (b2 in names(V)) {
    if (b != b2) s <- s + stats::cov(V[[b]][, t1], V[[b2]][, t2])
  }
  s
}


# -- B1: genic closed forms ------------------------------------------------------

test_that("B1: genic blocks sum to total and additive is sum 2pq alpha^2", {
  pop <- egv_pop()
  on.exit(close_pop(pop), add = TRUE)
  loci <- paste0("L", 1:6)
  a <- c(0.4, -0.2, 0.3, 0.1, -0.5, 0.25); d <- c(0.2, 0.1, -0.3, 0, 0.4, 0.15)
  pairs <- rbind(c(1, 2), c(3, 5)); e <- c(0.6, -0.4)
  pop <- egv_write(pop, "T", egv_functional_terms(loci, a, d, pairs, e))
  r <- egv(get_table(pop, "ind_meta"), anchor = "genic")

  X <- egv_dosage(pop)
  p <- colMeans(X) / 2; q <- 1 - p; cc <- 2 * p - 1
  alpha <- a + (q - p) * d
  alpha[1] <- alpha[1] + e[1] * cc[2]; alpha[2] <- alpha[2] + e[1] * cc[1]
  alpha[3] <- alpha[3] + e[2] * cc[5]; alpha[5] <- alpha[5] + e[2] * cc[3]
  vA  <- sum(2 * p * q * alpha^2)
  vD  <- sum((2 * p * q)^2 * d^2)
  vAA <- sum(e^2 * 4 * p[pairs[, 1]] * q[pairs[, 1]] * p[pairs[, 2]] *
               q[pairs[, 2]])
  expect_equal(egv_get(r, "additive", "T"), vA, tolerance = 1e-12)
  expect_equal(egv_get(r, "dominance", "T"), vD, tolerance = 1e-12)
  expect_equal(egv_get(r, "additive_by_additive", "T"), vAA, tolerance = 1e-12)
  expect_equal(egv_get(r, "total", "T"), vA + vD + vAA, tolerance = 1e-10)
  expect_false("between_components" %in% r$effect_name)
})


# -- B2: the source oracle ---------------------------------------------------

test_that("B2: realised blocks and cross-block covariances equal nonadd_decompose()", {
  pop <- egv_pop(n_ind = 120, n_hap = 12, seed = 21)   # few haplotypes: LD
  on.exit(close_pop(pop), add = TRUE)
  loci <- paste0("L", 1:6)
  set.seed(5)
  Ba <- matrix(stats::rnorm(12), 6, 2); Bd <- matrix(stats::rnorm(12), 6, 2)
  pairs <- rbind(c(1, 2), c(3, 5), c(2, 6)); Baa <- matrix(stats::rnorm(6), 3, 2)
  p0 <- stats::setNames(stats::runif(6, 0.2, 0.8), loci)
  for (j in 1:2) {
    stat <- .noia_to_stored(
      stats::setNames(Ba[, j], loci), stats::setNames(Bd[, j], loci),
      data.frame(locus_1 = loci[pairs[, 1]], locus_2 = loci[pairs[, 2]],
                 e = Baa[, j]), p0)
    pop <- egv_write(pop, paste0("T", j), .noia_terms(stat))   # Cockerham
  }
  r <- egv(get_table(pop, "ind_meta"))
  o <- nonadd_decompose(egv_dosage(pop), Ba, Bd, Baa, pairs)
  T <- c("T1", "T2")
  for (i in 1:2) for (j in 1:2) {
    expect_equal(egv_get(r, "additive", T[i], T[j]), o$real_A[i, j], tolerance = 1e-10)
    expect_equal(egv_get(r, "dominance", T[i], T[j]), o$real_D[i, j], tolerance = 1e-10)
    expect_equal(egv_get(r, "additive_by_additive", T[i], T[j]), o$real_AA[i, j],
                 tolerance = 1e-10)
    expect_equal(egv_get(r, "total", T[i], T[j]), o$real_G[i, j], tolerance = 1e-10)
    # The source returns Cov(A, D) and Cov(A, AA); Cov(D, AA) from its values.
    V <- list(A = o$BV, D = o$DD, AA = o$AA)
    expect_equal(egv_get(r, "between_components", T[i], T[j]),
                 egv_between(V, i, j), tolerance = 1e-10)
  }
  expect_gt(abs(o$cov_A_D[1, 2]), 1e-4)            # the panel has LD
  expect_equal(egv_between(list(A = o$BV, D = o$DD), 1, 2),
               o$cov_A_D[1, 2] + o$cov_A_D[2, 1], tolerance = 1e-12)
  expect_lt(o$identity_error, 1e-10)
})


# -- B3: accounting under LD --------------------------------------------------

test_that("B3: blocks plus between_components sum to total; blocks alone do not", {
  pop <- egv_pop(n_ind = 100, n_hap = 6, seed = 8)     # strong LD
  on.exit(close_pop(pop), add = TRUE)
  loci <- paste0("L", 1:6)
  pairs <- rbind(c(1, 2), c(4, 6))
  A1 <- c(0.5, -0.3, 0.2, 0.4, 0, -0.2); D1 <- c(0.6, 0.1, 0, -0.4, 0.3, 0.2)
  A2 <- c(-0.1, 0.4, 0.3, 0, 0.5, 0.2);  D2 <- c(0, 0.5, -0.3, 0.2, 0.1, -0.6)
  pop <- egv_write(pop, "T1", egv_functional_terms(loci, A1, D1, pairs, c(0.8, -0.5)))
  pop <- egv_write(pop, "T2", egv_functional_terms(loci, A2, D2, pairs, c(-0.3, 0.9)))
  r <- egv(get_table(pop, "ind_meta"))
  o <- nonadd_decompose(egv_dosage(pop), cbind(A1, A2), cbind(D1, D2),
                        rbind(c(0.8, -0.3), c(-0.5, 0.9)), pairs)
  V <- list(A = o$BV, D = o$DD, AA = o$AA)
  # Asymmetric: the two orientations of a cross term differ.
  expect_gt(abs(stats::cov(o$BV[, 1], o$DD[, 2]) - stats::cov(o$DD[, 1], o$BV[, 2])),
            1e-6)
  expect_gt(abs(stats::cov(o$DD[, 1], o$AA[, 2])), 1e-6)
  T <- c("T1", "T2")
  for (i in 1:2) for (j in 1:2) {
    blocks <- egv_get(r, "additive", T[i], T[j]) + egv_get(r, "dominance", T[i], T[j]) +
      egv_get(r, "additive_by_additive", T[i], T[j])
    tot <- egv_get(r, "total", T[i], T[j])
    btw <- egv_get(r, "between_components", T[i], T[j])
    expect_gt(abs(tot - blocks), 1e-6)
    expect_equal(blocks + btw, tot, tolerance = 1e-10)
    expect_equal(btw, egv_between(V, i, j), tolerance = 1e-10)
  }
  expect_false("unpartitioned" %in% r$effect_name)
  expect_false("between_components" %in%
                 egv(get_table(pop, "ind_meta"), anchor = "genic")$effect_name)
})

test_that("B3: unpartitioned cross terms are in between_components too", {
  pop <- egv_pop(n_ind = 100, n_hap = 6, seed = 8)
  on.exit(close_pop(pop), add = TRUE)
  loci <- paste0("L", 1:6)
  surf <- genotype_terms(expand.grid(L2 = 0:2, L3 = 0:2),
                         c(0, 0.5, 1, -0.4, 0.8, 0.3, 1.1, -0.2, 0.6))
  pop <- egv_write(pop, "T1", rbind(
    egv_functional_terms(loci[1:4], c(0.5, -0.3, 0.2, 0.4), c(0.6, 0.1, 0, -0.4)),
    surf))
  pop <- egv_write(pop, "T2", egv_functional_terms(loci[1:4], c(-0.1, 0.4, 0.3, 0),
                                                   c(0, 0.5, -0.3, 0.2)))
  r <- egv(get_table(pop, "ind_meta"))
  X <- egv_dosage(pop)
  o <- nonadd_decompose(X[, 1:4], cbind(c(0.5, -0.3, 0.2, 0.4), c(-0.1, 0.4, 0.3, 0)),
                        cbind(c(0.6, 0.1, 0, -0.4), c(0, 0.5, -0.3, 0.2)))
  pop <- egv_write(pop, "S", surf)
  U <- cbind(egv_total(pop, "S"), 0)
  V <- list(A = o$BV, D = o$DD, U = U)
  T <- c("T1", "T2")
  for (i in 1:2) for (j in 1:2) {
    expect_equal(egv_get(r, "between_components", T[i], T[j]),
                 egv_between(V, i, j), tolerance = 1e-10)
  }
  expect_equal(egv_get(r, "unpartitioned", "T1"), stats::var(U[, 1]),
               tolerance = 1e-10)
  expect_true(is.na(egv_get(r, "unpartitioned", "T2")))
  expect_true(is.na(egv_get(r, "unpartitioned", "T1", "T2")))
})


# -- B4 / B6 / B13: dispatch and the family rule --------------------------------

# Two lines, then F1s of A sires x B dams.
egv_cross_pop <- function(seed = 3, n_f1 = 60) {
  set.seed(seed)
  pop <- suppressMessages(
    open_pop(pop_name = "egvx", db_name = ":memory:") |>
      define_genome(n_loci = 6, n_chr = 1, chr_len_Mb = 10,
                    locus_names = paste0("L", 1:6)))
  pop <- suppressMessages(define_founder_haplotypes(pop, n_haplotypes = 20,
                                                    line_name = "LA"))
  pop <- suppressMessages(define_founder_haplotypes(pop, n_haplotypes = 20,
                                                    line_name = "LB"))
  pop <- suppressMessages(pop |> get_table("founder_haplotypes") |>
    dplyr::filter(line_name == "LA") |>
    add_founders(n_males = 10, n_females = 10, line_name = "LA", gen = 0L))
  pop <- suppressMessages(pop |> get_table("founder_haplotypes") |>
    dplyr::filter(line_name == "LB") |>
    add_founders(n_males = 10, n_females = 10, line_name = "LB", gen = 0L))
  matings <- tibble::tibble(
    id_parent_1 = paste0("LA_", rep(1:10, length.out = n_f1)),
    id_parent_2 = paste0("LB_", rep(11:20, length.out = n_f1)),
    sex = rep(c("M", "F"), length.out = n_f1), line_name = "F1", gen = 1L)
  suppressMessages(add_offspring(pop, matings))
}

egv_line_additive <- function(pop, trait, effects, line_name = NULL) {
  suppressMessages(pop |> get_table("genome_meta") |>
    with_additive_terms(trait, effects = effects, line_name = line_name,
                        base_tbl = get_table(pop, "founder_haplotypes")))
}

test_that("B4: two owners on the same loci are one case-1 model", {
  pop <- egv_pop()
  on.exit(close_pop(pop), add = TRUE)
  loci <- paste0("L", 1:3)
  pop <- egv_write(pop, "T", suppressMessages(ad_terms(loci, a = c(0.3, 0.2, -0.4),
                                                       d = 0, p = 0.5)), "one")
  pop <- egv_write(pop, "T", genotype_terms(data.frame(L1 = 1L), 0.7), "two")
  pop <- egv_write(pop, "T", data.frame(locus_name = "L2", contrast_name = "dominance",
                                        center_value = 0.4, genome_value = 0.5), "two")
  r <- egv(get_table(pop, "ind_meta"))
  expect_true(all(r$decomposition == "full"))
  expect_true(all(c("additive", "dominance") %in% r$effect_name))
  expect_false("unpartitioned" %in% r$effect_name)
})

test_that("B4/B6: a line-scoped F1 additive model is case 2, the evaluated additive variance", {
  pop <- egv_cross_pop()
  on.exit(close_pop(pop), add = TRUE)
  pop <- suppressMessages(define_trait(pop, "T"))
  pop <- egv_line_additive(pop, "T", c(0.2, -0.1, 0.3, 0.05, -0.2, 0.15))
  pop <- egv_line_additive(pop, "T", c(0.5, 0.1, -0.3, 0.2, 0.1, -0.4), "LA")
  pop <- egv_line_additive(pop, "T", c(-0.2, 0.4, 0.1, -0.3, 0.2, 0.3), "LB")
  f1 <- get_table(pop, "ind_meta") |> dplyr::filter(line_name == "F1")
  expect_message(r <- extract_genetic_variance(f1), "Evaluated additive variance")
  expect_true(all(r$decomposition == "additive_only"))
  expect_setequal(unique(r$effect_name), c("additive", "between_components", "total"))
  ids <- DBI::dbGetQuery(pop$db_conn,
    "SELECT id_ind FROM ind_meta WHERE line_name = 'F1' ORDER BY id_ind")$id_ind
  suppressMessages(f1 |> add_tgv("T"))
  a <- DBI::dbGetQuery(pop$db_conn, paste0(
    "SELECT id_ind, tgv_value FROM ind_tgv WHERE trait_name = 'T' AND ",
    "component_name = 'additive'"))
  v <- a$tgv_value[match(ids, a$id_ind)]
  expect_equal(egv_get(r, "additive", "T"), stats::var(v), tolerance = 1e-12)
  expect_equal(egv_get(r, "total", "T"), stats::var(v), tolerance = 1e-12)
  expect_equal(egv_get(r, "between_components", "T"), 0)
  expect_error(suppressMessages(extract_genetic_variance(f1, anchor = "genic")),
               "realised")
})

test_that("B13: a family with a scoped variant is uncovered as a whole", {
  pop <- egv_cross_pop()
  on.exit(close_pop(pop), add = TRUE)
  pop <- suppressMessages(define_trait(pop, "T"))
  pop <- egv_line_additive(pop, "T", c(0.2, -0.1, 0.3, 0.05, -0.2, 0.15))
  pop <- egv_line_additive(pop, "T", c(0.5, 0.1, -0.3, 0.2, 0.1, -0.4), "LA")
  pop <- egv_write(pop, "T", data.frame(locus_name = "L2", contrast_name = "dominance",
                                        center_value = 0.5, genome_value = 0.8), "dom")
  cohort <- get_table(pop, "ind_meta")
  r <- egv(cohort)
  expect_true(all(r$decomposition == "partial"))
  ids <- egv_ids(pop)
  g <- egv_total(pop, "T", ids)
  expect_equal(egv_get(r, "total", "T"), stats::var(g), tolerance = 1e-10)
  # The additive families (common + line variants) are unpartitioned: their
  # evaluated value alone.
  pop <- suppressMessages(define_trait(pop, "ADD"))
  pop <- egv_line_additive(pop, "ADD", c(0.2, -0.1, 0.3, 0.05, -0.2, 0.15))
  pop <- egv_line_additive(pop, "ADD", c(0.5, 0.1, -0.3, 0.2, 0.1, -0.4), "LA")
  expect_equal(egv_get(r, "unpartitioned", "T"), stats::var(egv_total(pop, "ADD", ids)),
               tolerance = 1e-10)
  expect_equal(egv_get(r, "additive", "T") + egv_get(r, "dominance", "T") +
                 egv_get(r, "unpartitioned", "T") +
                 egv_get(r, "between_components", "T"),
               egv_get(r, "total", "T"), tolerance = 1e-10)
})


# -- B5 / B11: partial models and the labels ------------------------------------

test_that("B5: a multi-locus indicator surface is partial and goes to unpartitioned", {
  pop <- egv_pop()
  on.exit(close_pop(pop), add = TRUE)
  surf <- genotype_terms(expand.grid(L1 = 0:2, L4 = 0:2),
                         c(0, 0, 0, 0, 1.4, 2.1, 0, 2.1, 3.6))
  pop <- egv_write(pop, "T", surf)
  r <- egv(get_table(pop, "ind_meta"))
  g <- egv_total(pop, "T")
  expect_equal(egv_get(r, "total", "T"), stats::var(g), tolerance = 1e-10)
  expect_equal(egv_get(r, "unpartitioned", "T"), stats::var(g), tolerance = 1e-10)
  expect_true(all(r$decomposition == "partial"))
  expect_false("additive" %in% r$effect_name)   # nothing decomposed
})

test_that("B11: decomposition labels per case and on mixed off-diagonals", {
  pop <- egv_pop()
  on.exit(close_pop(pop), add = TRUE)
  pop <- egv_write(pop, "F", egv_functional_terms(paste0("L", 1:2), c(0.3, 0.1),
                                                  c(0.2, 0)))
  pop <- egv_write(pop, "P", genotype_terms(expand.grid(L1 = 0:2, L2 = 0:2),
                                            seq(0.1, 0.9, by = 0.1)))
  r <- egv(get_table(pop, "ind_meta"))
  expect_equal(unique(r$decomposition[r$trait_name_1 == "F" & r$trait_name_2 == "F"]),
               "full")
  expect_equal(unique(r$decomposition[r$trait_name_1 == "P" & r$trait_name_2 == "P"]),
               "partial")
  expect_equal(unique(r$decomposition[r$trait_name_1 == "F" & r$trait_name_2 == "P"]),
               "partial")
})


# -- B7: read-only ------------------------------------------------------------

test_that("B7: the call writes nothing", {
  pop <- egv_pop()
  on.exit(close_pop(pop), add = TRUE)
  pop <- egv_write(pop, "T", egv_functional_terms(paste0("L", 1:4), c(0.3, 0.1, 0, 0.2),
                                                  c(0.2, 0, 0.4, 0), rbind(c(1, 2)), 0.5))
  suppressMessages(get_table(pop, "ind_meta") |> add_tgv("T"))
  snap <- function() {
    tabs <- DBI::dbListTables(pop$db_conn)
    # Row order of an unordered SELECT is not part of the state: sort.
    stats::setNames(lapply(tabs, function(t) {
      d <- DBI::dbGetQuery(pop$db_conn, paste0("SELECT * FROM \"", t, "\""))
      if (nrow(d) > 0L) d <- d[do.call(order, unname(as.list(d))), , drop = FALSE]
      rownames(d) <- NULL
      d
    }), tabs)
  }
  before <- snap()
  egv(get_table(pop, "ind_meta"))
  egv(get_table(pop, "ind_meta"), anchor = "genic")
  expect_identical(snap(), before)
})


# -- B8: the join to trait_var_comp ---------------------------------------------

test_that("B8: rows join trait_var_comp; a target with no block is in the anti_join", {
  pop <- egv_pop()
  on.exit(close_pop(pop), add = TRUE)
  pop <- with_additive_target(pop, "T", 1)
  pop <- suppressMessages(define_effect_cov_matrix(pop, "dominance", 0.3,
                                                   trait_name = "T"))
  pop <- egv_write(pop, "T", suppressMessages(
    ad_terms(paste0("L", 1:4), a = c(0.3, 0.1, -0.2, 0.2), p = 0.5)))
  r <- egv(get_table(pop, "ind_meta"))
  targets <- get_table(pop, "trait_var_comp") |> dplyr::collect() |>
    dplyr::filter(is.na(line_name)) |>
    dplyr::select(effect_name, trait_name_1, trait_name_2, target = cov_value)
  by <- c("effect_name", "trait_name_1", "trait_name_2")
  expect_equal(dplyr::inner_join(targets, r, by = by)$effect_name, "additive")
  expect_equal(dplyr::anti_join(targets, r, by = by)$effect_name, "dominance")
})


# -- B9: determinism and the size guard -------------------------------------------

test_that("B9: identical output on repeat, across threads and after restore_pop()", {
  f <- tempfile(fileext = ".duckdb")
  pop <- egv_pop(n_ind = 120, pop_name = "egvdet", db_name = f)
  loci <- paste0("L", 1:6)
  pop <- egv_write(pop, "T1", egv_functional_terms(loci, c(0.3, 0.1, 0, 0.2, -0.4, 0.1),
                                                   c(0.2, 0, 0.4, 0, 0.1, -0.2),
                                                   rbind(c(1, 2), c(3, 6)), c(0.5, -0.7)))
  pop <- egv_write(pop, "T2", genotype_terms(expand.grid(L2 = 0:2, L5 = 0:2),
                                             seq(-0.4, 0.4, by = 0.1)))
  one <- egv(get_table(pop, "ind_meta"))
  expect_identical(egv(get_table(pop, "ind_meta")), one)
  for (th in c(1, 4)) {
    DBI::dbExecute(pop$db_conn, paste0("SET threads = ", th))
    expect_identical(egv(get_table(pop, "ind_meta")), one)
  }
  close_pop(pop)
  pop <- suppressMessages(restore_pop(f))
  on.exit({ close_pop(pop); unlink(f) }, add = TRUE)
  expect_identical(egv(get_table(pop, "ind_meta")), one)
})

test_that("B9/B17: the dosage guard fires before any evaluation, naming n, m and the limit", {
  pop <- egv_pop()
  on.exit(close_pop(pop), add = TRUE)
  pop <- egv_write(pop, "T", egv_functional_terms(paste0("L", 1:6), rep(0.1, 6),
                                                  rep(0.1, 6)))
  local_mocked_bindings(QTL_REALISED_MAX_CELLS = 100)
  local_mocked_bindings(.gev_evaluate = function(...) stop("evaluated"))
  expect_error(egv(get_table(pop, "ind_meta")),
               "80 x 6 genotype matrix.*above the limit of 100")
})


# -- B10: coding invariance and anchor scope -----------------------------------

test_that("B10: functional and Cockerham codings of one model give the same report", {
  pop <- egv_pop(n_ind = 100, seed = 4)
  on.exit(close_pop(pop), add = TRUE)
  loci <- paste0("L", 1:5)
  a <- stats::setNames(c(0.3, -0.2, 0.1, 0.4, -0.1), loci)
  d <- stats::setNames(c(0.2, 0.3, -0.1, 0, 0.25), loci)
  pr <- data.frame(locus_1 = c("L1", "L3"), locus_2 = c("L4", "L5"), e = c(0.5, -0.6))
  pop <- egv_write(pop, "FUN", egv_functional_terms(loci, a, d,
                                                    rbind(c(1, 4), c(3, 5)), pr$e))
  pop <- egv_write(pop, "COC", .noia_terms(.noia_to_stored(
    a, d, pr, stats::setNames(c(0.3, 0.45, 0.6, 0.7, 0.2), loci))))
  for (anc in c("realised", "genic")) {
    r1 <- egv(get_table(pop, "ind_meta"), trait_name = "FUN", anchor = anc)
    r2 <- egv(get_table(pop, "ind_meta"), trait_name = "COC", anchor = anc)
    expect_identical(r1$effect_name, r2$effect_name)
    expect_equal(r1$cov_value, r2$cov_value, tolerance = 1e-10)
  }
  expect_error(egv(get_table(pop, "ind_meta"), anchor = "realised",
                   base_tbl = get_table(pop, "founder_haplotypes")),
               "genic")
})


# -- B12: the cohort -----------------------------------------------------------

test_that("B12: fewer than two individuals, or one without a value, is an error", {
  pop <- egv_cross_pop()
  on.exit(close_pop(pop), add = TRUE)
  pop <- suppressMessages(define_trait(pop, "T"))
  # Only LB-scoped terms: LA purebreds carry no matching copy.
  pop <- egv_line_additive(pop, "T", c(0.2, -0.1, 0.3, 0.05, -0.2, 0.15), "LB")
  expect_error(egv(get_table(pop, "ind_meta") |> dplyr::filter(line_name == "LA")),
               "T \\(20\\)")
  expect_error(egv(get_table(pop, "ind_meta") |> dplyr::filter(id_ind == "LB_11")),
               "at least 2 individuals")
  r <- egv(get_table(pop, "ind_meta") |> dplyr::filter(line_name == "F1"))
  expect_true(all(r$n_ind == 60L))
})

# A pair-only functional model with locus L2 fixed at `dose` in the cohort.
egv_fixed_pair <- function(dose) {
  pop <- egv_pop(n_ind = 60, seed = 9)
  DBI::dbExecute(pop$db_conn, paste0(
    "UPDATE ind_haplotype SET allele = CASE WHEN parent_origin = 1 THEN ",
    as.integer(dose >= 1), " ELSE ", as.integer(dose == 2), " END ",
    "WHERE locus_id = (SELECT locus_id FROM genome_meta WHERE locus_name = 'L2')"))
  pop
}

test_that("B12: a fixed pair member keeps its induced additive effect (fixed 0, 1, 2)", {
  for (dose in 0:2) {
    pop <- egv_fixed_pair(dose)
    e <- 0.9
    pop <- egv_write(pop, "PAIR", aa_terms("L1", "L2", e = e, p_1 = 0.5, p_2 = 0.5,
                                           report = FALSE))
    # Reduced model: drop the pair, carry e * (g2 - 1) onto L1, drop the constant.
    if (dose != 1) {
      pop <- egv_write(pop, "RED", suppressMessages(
        ad_terms("L1", a = e * (dose - 1), p = 0.5)))
    }
    r <- egv(get_table(pop, "ind_meta"))
    expect_true(all(is.finite(r$cov_value)), info = dose)
    expect_equal(egv_get(r, "additive_by_additive", "PAIR"), 0, info = dose)
    if (dose == 1) {
      expect_equal(egv_get(r, "additive", "PAIR"), 0, info = dose)
      expect_equal(egv_get(r, "total", "PAIR"), 0, info = dose)
    } else {
      expect_gt(egv_get(r, "additive", "PAIR"), 0.01)
      expect_equal(egv_get(r, "additive", "PAIR"), egv_get(r, "additive", "RED"),
                   tolerance = 1e-12, info = dose)
      expect_equal(egv_get(r, "total", "PAIR"), egv_get(r, "total", "RED"),
                   tolerance = 1e-12, info = dose)
    }
    close_pop(pop)
  }
})


# -- B14 / B16: block availability ----------------------------------------------

test_that("B14: a one-locus genotype surface is case 1 and equals ad_terms()", {
  pop <- egv_pop()
  on.exit(close_pop(pop), add = TRUE)
  # values g0, g1, g2 = 0.1, 0.9, 0.5: a = (g2 - g0)/2, d = g1 - (g0 + g2)/2
  pop <- egv_write(pop, "SURF", genotype_terms(data.frame(L3 = 0:2), c(0.1, 0.9, 0.5)))
  pop <- egv_write(pop, "AD", suppressMessages(ad_terms("L3", a = 0.2, d = 0.6, p = 0.5)))
  for (anc in c("realised", "genic")) {
    s <- egv(get_table(pop, "ind_meta"), trait_name = "SURF", anchor = anc)
    a <- egv(get_table(pop, "ind_meta"), trait_name = "AD", anchor = anc)
    expect_true(all(s$decomposition == "full"))
    expect_identical(s$effect_name, a$effect_name)
    expect_equal(s$cov_value, a$cov_value, tolerance = 1e-10)
  }
})

test_that("B16: rows follow the canonical decomposition, not stored term kinds", {
  pop <- egv_pop(n_ind = 100, seed = 12)
  on.exit(close_pop(pop), add = TRUE)
  # Dosage-only indicator surface: all additive, no dominance row.
  pop <- egv_write(pop, "DOSE", genotype_terms(data.frame(L1 = 0:2), c(0, 1, 2)))
  # Heterozygote only: additive (q - p) d and dominance.
  pop <- egv_write(pop, "HET", genotype_terms(data.frame(L2 = 1L), 1))
  # Functional pair only: additive (through e c) and A x A.
  pop <- egv_write(pop, "FPAIR", aa_terms("L3", "L4", e = 1, p_1 = 0.5, p_2 = 0.5,
                                          report = FALSE))
  r <- egv(get_table(pop, "ind_meta"))
  has <- function(eff, t1, t2 = t1) !is.na(egv_get(r, eff, t1, t2))
  g_dose <- egv_total(pop, "DOSE")
  expect_equal(egv_get(r, "additive", "DOSE"), stats::var(g_dose), tolerance = 1e-12)
  expect_equal(egv_get(r, "total", "DOSE"), stats::var(g_dose), tolerance = 1e-12)
  expect_false(has("dominance", "DOSE"))
  expect_false(has("additive_by_additive", "DOSE"))
  expect_true(has("additive", "HET") && has("dominance", "HET"))
  expect_gt(egv_get(r, "additive", "HET"), 1e-4)       # unequal genotype frequencies
  expect_true(has("additive", "FPAIR") && has("additive_by_additive", "FPAIR"))
  expect_false(has("dominance", "FPAIR"))
  # Different block support between traits: no off-diagonal dominance row.
  expect_false(has("dominance", "HET", "DOSE"))
  expect_true(has("additive", "HET", "DOSE"))
  expect_true(has("between_components", "HET", "DOSE"))
  expect_true(all(r$decomposition == "full"))
})

test_that("B16: a Cockerham pair-only model measured after drift reports additive", {
  pop <- egv_pop(n_ind = 100, seed = 12)
  on.exit(close_pop(pop), add = TRUE)
  # Centred at frequencies the cohort does not have.
  pop <- egv_write(pop, "CPAIR", aa_terms("L3", "L4", e = 1, p_1 = 0.1, p_2 = 0.9,
                                          coding = "cockerham"))
  r <- egv(get_table(pop, "ind_meta"))
  expect_gt(egv_get(r, "additive", "CPAIR"), 1e-4)
  expect_false(is.na(egv_get(r, "additive_by_additive", "CPAIR")))
  expect_equal(egv_get(r, "total", "CPAIR"), stats::var(egv_total(pop, "CPAIR")),
               tolerance = 1e-10)
})


# -- B15: anchor and message -------------------------------------------------------

test_that("B15: every row carries the anchor; the message names the population", {
  pop <- egv_pop()
  on.exit(close_pop(pop), add = TRUE)
  pop <- egv_write(pop, "T", egv_functional_terms(paste0("L", 1:3), c(0.3, 0.1, 0.2),
                                                  c(0.1, 0, 0)))
  expect_message(r1 <- extract_genetic_variance(get_table(pop, "ind_meta")),
                 "Realised .* 80 selected individuals")
  expect_true(all(r1$anchor == "realised"))
  expect_message(r2 <- extract_genetic_variance(get_table(pop, "ind_meta"),
                                                anchor = "genic"),
                 "whole-genotype allele frequencies of the 80 selected")
  expect_true(all(r2$anchor == "genic"))
  expect_message(extract_genetic_variance(
    get_table(pop, "ind_meta"), anchor = "genic",
    base_tbl = get_table(pop, "founder_haplotypes")),
    "'founder_haplotypes', an explicit base \\(allele copies\\)")
})


# -- B17: chunked pairs -------------------------------------------------------------

test_that("B17: A x A values do not depend on the pair chunk size", {
  set.seed(2)
  X <- matrix(stats::rbinom(50 * 20, 2, 0.4), 50, 20)
  Z <- sweep(X, 2L, colMeans(X), "-")
  pairs <- t(utils::combn(20, 2))
  E <- matrix(stats::rnorm(nrow(pairs) * 2), nrow(pairs), 2)
  whole <- .egv_aa_values(Z, pairs, E, chunk = nrow(pairs))
  # The source's dense form.
  ZAA <- Z[, pairs[, 1]] * Z[, pairs[, 2]]
  ZAA <- sweep(ZAA, 2L, colMeans(ZAA), "-")
  expect_equal(unname(whole), unname(ZAA %*% E), tolerance = 1e-12)
  for (ch in c(1, 7, 64)) {
    expect_equal(.egv_aa_values(Z, pairs, E, chunk = ch), whole, tolerance = 1e-12)
  }
  expect_identical(.egv_aa_values(Z, pairs, E, chunk = 7),
                   .egv_aa_values(Z, pairs, E, chunk = 7))
  expect_equal(.egv_pair_chunk(2000), 10000L)
})

test_that("B17: a many-pairs model runs in chunks and matches the oracle", {
  pop <- egv_pop(loci = paste0("L", 1:12), n_ind = 60, seed = 6)
  on.exit(close_pop(pop), add = TRUE)
  loci <- paste0("L", 1:12)
  pairs <- t(utils::combn(12, 2))
  set.seed(1)
  e <- stats::rnorm(nrow(pairs), sd = 0.2)
  pop <- egv_write(pop, "T", egv_functional_terms(loci, rep(0.1, 12), rep(0, 12),
                                                  pairs, e))
  local_mocked_bindings(QTL_REALISED_MAX_CELLS = 720)   # 12 pairs per chunk
  r <- egv(get_table(pop, "ind_meta"))
  o <- nonadd_decompose(egv_dosage(pop), matrix(rep(0.1, 12)), NULL,
                        matrix(e), pairs)
  expect_equal(egv_get(r, "additive_by_additive", "T"), o$real_AA[1, 1],
               tolerance = 1e-10)
  expect_equal(egv_get(r, "total", "T"), o$real_G[1, 1], tolerance = 1e-10)
})


# -- B18: the default base and the trait default ----------------------------------

test_that("B18: base_tbl = NULL uses the cohort's whole genotypes, whatever table selects it", {
  pop <- egv_pop(n_ind = 40)
  on.exit(close_pop(pop), add = TRUE)
  pop <- egv_write(pop, "T", egv_functional_terms(paste0("L", 1:4), c(0.3, 0.1, -0.2, 0.2),
                                                  c(0.2, 0.1, 0, 0.3)))
  pop <- suppressMessages(define_phenotype(pop, "T", residual_var = 1,
                                           repeatable = TRUE))
  set.seed(1)
  suppressMessages(suppressWarnings(get_table(pop, "ind_meta") |> add_phenotype("T")))
  suppressMessages(suppressWarnings(get_table(pop, "ind_meta") |> add_phenotype("T")))
  g <- function(tbl) egv(tbl, anchor = "genic")
  ref <- g(get_table(pop, "ind_meta"))
  expect_identical(g(get_table(pop, "ind_haplotype") |>
                       dplyr::filter(parent_origin == 1L)), ref)
  expect_identical(g(get_table(pop, "ind_phenotype")), ref)

  copies <- egv(get_table(pop, "ind_meta"), anchor = "genic",
                base_tbl = get_table(pop, "ind_haplotype") |>
                  dplyr::filter(parent_origin == 1L))
  expect_false(isTRUE(all.equal(copies$cov_value, ref$cov_value)))
  expect_error(egv(get_table(pop, "ind_meta"), anchor = "genic",
                   base_tbl = get_table(pop, "ind_haplotype") |>
                     dplyr::filter(locus_id != 2L)),
               "no allele copies at decomposed locus.*'L2'")
})

test_that("B18: trait_name = NULL skips traits without terms; naming one errors", {
  pop <- egv_pop()
  on.exit(close_pop(pop), add = TRUE)
  pop <- suppressMessages(define_trait(pop, "Empty"))
  expect_error(egv(get_table(pop, "ind_meta")), "No trait has genome-effect terms")
  pop <- egv_write(pop, "T", egv_functional_terms("L1", 0.3, 0))
  r <- egv(get_table(pop, "ind_meta"))
  expect_equal(unique(r$trait_name_1), "T")
  expect_error(egv(get_table(pop, "ind_meta"), trait_name = c("T", "Empty")),
               "No genome-effect terms for trait\\(s\\) 'Empty'")
})


# -- Codex implementation review: cancellation and the LD meaning ---------------

test_that("B16: equivalent linear surfaces have the same blocks in any row order", {
  pop <- egv_pop(n_ind = 100, seed = 12)
  on.exit(close_pop(pop), add = TRUE)
  # One linear surface on L1 (g0, g1, g2 = 0.3, 0.2, 0.1): a = -0.1, d = 0.
  # Summing its rows' dominance contributions leaves a 1e-17 residue in some
  # orders and exactly 0 in others.
  vals <- c(0.3, 0.2, 0.1)
  perms <- list(1:3, c(1, 3, 2), c(2, 1, 3), c(2, 3, 1), c(3, 1, 2), 3:1)
  for (i in seq_along(perms)) {
    o <- perms[[i]]
    pop <- egv_write(pop, paste0("S", i),
                     genotype_terms(data.frame(L1 = (0:2)[o]), vals[o]))
  }
  pop <- egv_write(pop, "AD", ad_terms("L1", a = -0.1, d = 0, p = 0.5,
                                       report = FALSE))
  # Genuinely small, uncancelled dominance must survive.
  pop <- egv_write(pop, "SMALL", genotype_terms(data.frame(L1 = 0:2),
                                                c(0.3, 0.2 + 1e-9, 0.1)))
  pop <- egv_write(pop, "TINY", ad_terms("L1", a = 0, d = 1e-20, p = 0.5,
                                         report = FALSE))
  for (anc in c("realised", "genic")) {
    for (t in c(paste0("S", seq_along(perms)), "AD")) {
      r <- egv(get_table(pop, "ind_meta"), trait_name = t, anchor = anc)
      expect_false("dominance" %in% r$effect_name, info = paste(anc, t))
      expect_true("additive" %in% r$effect_name, info = paste(anc, t))
    }
    for (t in c("SMALL", "TINY")) {
      r <- egv(get_table(pop, "ind_meta"), trait_name = t, anchor = anc)
      expect_true("dominance" %in% r$effect_name, info = paste(anc, t))
    }
  }
  # Owners summing a pair to residue: no additive_by_additive row.
  pop <- egv_write(pop, "PC", aa_terms("L2", "L3", e = 0.1, p_1 = 0.5, p_2 = 0.5,
                                       report = FALSE), owner = "o1")
  pop <- egv_write(pop, "PC", aa_terms("L2", "L3", e = 0.2, p_1 = 0.5, p_2 = 0.5,
                                       report = FALSE), owner = "o2")
  pop <- egv_write(pop, "PC", aa_terms("L2", "L3", e = -0.3, p_1 = 0.5, p_2 = 0.5,
                                       report = FALSE), owner = "o3")
  r <- egv(get_table(pop, "ind_meta"), trait_name = "PC")
  expect_false("additive_by_additive" %in% r$effect_name)
})

test_that("the realised additive block is the contrast component, not lm() under LD", {
  # Codex implementation review finding 1. 32 individuals, both loci with
  # counts (8, 16, 8) -- exact HWE margins, p = 0.5 -- but in LD.
  pop <- egv_pop(loci = c("L1", "L2"), n_ind = 32, seed = 4)
  on.exit(close_pop(pop), add = TRUE)
  cells <- expand.grid(L1 = 0:2, L2 = 0:2)
  X <- cells[rep(1:9, c(3, 3, 2, 2, 9, 5, 3, 4, 1)), ]
  ids <- egv_ids(pop)
  set_g <- data.frame(id_ind = rep(ids, 2),
                      locus_name = rep(c("L1", "L2"), each = 32),
                      g = c(X$L1, X$L2))
  duckdb::duckdb_register(pop$db_conn, "egv_set_g", set_g)
  DBI::dbExecute(pop$db_conn, paste0(
    "UPDATE ind_haplotype h SET allele = CASE WHEN h.parent_origin = 1 ",
    "THEN CAST(s.g >= 1 AS INTEGER) ELSE CAST(s.g = 2 AS INTEGER) END ",
    "FROM egv_set_g s JOIN genome_meta m USING (locus_name) ",
    "WHERE h.id_ind = s.id_ind AND h.locus_id = m.locus_id"))
  duckdb::duckdb_unregister(pop$db_conn, "egv_set_g")
  G <- egv_dosage(pop, ids)
  expect_equal(colMeans(G) / 2, c(0.5, 0.5))
  expect_equal(as.numeric(table(G[, 1])), c(8, 16, 8))

  pop <- egv_write(pop, "E", aa_terms("L1", "L2", e = 1, p_1 = 0.5, p_2 = 0.5,
                                      report = FALSE))
  r <- egv(get_table(pop, "ind_meta"))
  g <- (G[, 1] - 1) * (G[, 2] - 1)
  expect_equal(egv_get(r, "additive", "E"), 0)
  expect_equal(egv_get(r, "additive_by_additive", "E"), stats::var(g),
               tolerance = 1e-12)
  expect_equal(egv_get(r, "between_components", "E"), 0, tolerance = 1e-12)
  expect_equal(egv_get(r, "total", "E"), stats::var(g), tolerance = 1e-12)
  expect_true(all(r$decomposition == "full"))
  # The joint additive regression is a different estimator: it explains part
  # of g, which the documented contrast component does not claim to.
  fit <- stats::lm(g ~ G[, 1] + G[, 2])
  expect_equal(stats::var(stats::fitted(fit)), 0.02099937, tolerance = 1e-6)
})
