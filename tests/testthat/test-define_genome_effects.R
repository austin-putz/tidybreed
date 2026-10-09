# define_genome_effects(): the public gates of step 5b
# (plans/import_qtl_effect_methods.md §11 C1-C18, C20;
# plans/import_qtl_effect_methods_phase_5_plan.md 5b.9 G1-G8). The
# calibration internals and the source suite are in
# test-genome-effects-calibration.R. Expectations are computed here from the
# stored rows and the individuals' dosages.

dge_pop <- function(name = "dge", n_loci = 40, n_ind = 200, n_hap = 300,
                    db_name = ":memory:", traits = c("T1", "T2", "T3")) {
  set.seed(8001)
  pop <- suppressMessages(suppressWarnings({
    p <- open_pop(pop_name = name, db_name = db_name) |>
      define_genome(n_loci = n_loci, n_chr = 2, chr_len_Mb = 100) |>
      define_founder_haplotypes(n_haplotypes = n_hap)
    p <- p |> get_table("founder_haplotypes") |>
      add_founders(n_males = n_ind / 2, n_females = n_ind / 2, line_name = "A")
    p <- get_table(p, "ind_meta") |> mutate_table(gen = 0L)
    for (t in traits) p <- define_trait(p, t)
    p
  }))
  pop
}

quiet <- function(expr) suppressWarnings(suppressMessages(expr))
gm_at <- function(pop, n = 20) get_table(pop, "genome_meta") |> dplyr::filter(locus_id <= n)
gen0  <- function(pop) get_table(pop, "ind_meta") |> dplyr::filter(gen == 0L)
nm2   <- function(M, t = c("T1", "T2")) { dimnames(M) <- list(t, t); M }

GA2  <- nm2(matrix(c(1, 0.3, 0.3, 2), 2))
GD2  <- nm2(matrix(c(0.3, 0.1, 0.1, 0.5), 2))
GAA2 <- nm2(matrix(c(0.2, 0.05, 0.05, 0.3), 2))

# Every row of the three genome-effect tables, keyed without the id columns
# next_int_id() assigns, in (trait_name, term content, locus) order.
rows_noid <- function(pop) {
  DBI::dbGetQuery(pop$db_conn, paste0(
    "SELECT e.trait_name, e.effect_owner, e.effect_name, e.genome_value, ",
    "       m.member_slot, m.locus_id, m.contrast_name, m.copy_count_value, ",
    "       m.dosage_value, m.center_value, ",
    "       (SELECT COUNT(*) FROM genome_effect_member_origins o ",
    "         WHERE o.id_genome_effect = e.id_genome_effect) AS n_origins ",
    "FROM genome_effects e JOIN genome_effect_members m USING (id_genome_effect) ",
    "ORDER BY e.trait_name, m.locus_id, m.contrast_name, e.genome_value, m.member_slot"))
}
tvc <- function(pop) DBI::dbGetQuery(pop$db_conn, paste0(
  "SELECT effect_name, line_name, trait_name_1, trait_name_2, cov_value ",
  "FROM trait_var_comp ORDER BY effect_name, trait_name_1, trait_name_2"))
db_state <- function(pop) list(
  t = DBI::dbGetQuery(pop$db_conn, "SELECT * FROM genome_effects ORDER BY id_genome_effect"),
  m = DBI::dbGetQuery(pop$db_conn, "SELECT * FROM genome_effect_members ORDER BY id_genome_effect, member_slot"),
  o = DBI::dbGetQuery(pop$db_conn, "SELECT * FROM genome_effect_member_origins ORDER BY id_genome_effect, member_slot, origin_slot"),
  v = tvc(pop))

# The stored common-scope model as coefficient matrices: alpha and d per
# (locus, trait), e per (pair, trait), with the stored centres.
stored_model <- function(pop, traits) {
  r <- DBI::dbGetQuery(pop$db_conn, paste0(
    "SELECT e.id_genome_effect AS id, e.trait_name, e.effect_owner, e.genome_value, ",
    "       m.member_slot, m.locus_id, m.contrast_name, m.center_value ",
    "FROM genome_effects e JOIN genome_effect_members m USING (id_genome_effect) ",
    "WHERE e.effect_owner = 'generated' AND e.trait_name IN (",
    paste0("'", traits, "'", collapse = ", "), ") ",
    "ORDER BY e.id_genome_effect, m.member_slot"))
  n_mem <- table(r$id)
  r$n <- as.integer(n_mem[as.character(r$id)])
  loci <- sort(unique(r$locus_id))
  m <- length(loci); k <- length(traits)
  al <- dd <- matrix(0, m, k, dimnames = list(loci, traits))
  p <- rep(NA_real_, m)
  one <- r[r$n == 1L, ]
  p[match(one$locus_id, loci)] <- one$center_value
  a1 <- one[one$contrast_name == "additive", ]
  al[cbind(match(a1$locus_id, loci), match(a1$trait_name, traits))] <- a1$genome_value
  d1 <- one[one$contrast_name == "dominance", ]
  dd[cbind(match(d1$locus_id, loci), match(d1$trait_name, traits))] <- d1$genome_value
  two <- r[r$n == 2L, ]
  ids <- unique(two$id)
  first  <- two[two$member_slot == 1L, ][match(ids, two$id[two$member_slot == 1L]), ]
  second <- two[two$member_slot == 2L, ][match(ids, two$id[two$member_slot == 2L]), ]
  key <- paste(pmin(first$locus_id, second$locus_id), pmax(first$locus_id, second$locus_id))
  uk <- unique(key[order(pmin(first$locus_id, second$locus_id), pmax(first$locus_id, second$locus_id))])
  E <- matrix(0, length(uk), k, dimnames = list(uk, traits))
  E[cbind(match(key, uk), match(first$trait_name, traits))] <- first$genome_value
  pr <- do.call(rbind, lapply(strsplit(uk, " "), as.integer))
  pc <- if (length(uk)) {
    cbind(first$center_value[match(uk, key)], second$center_value[match(uk, key)])
  }
  pi <- if (length(uk)) cbind(match(pr[, 1], loci), match(pr[, 2], loci))
  list(alpha = al, d = dd, E = E, pairs = pi, pair_loci = pr, pair_centres = pc,
       p = p, loci = loci, n_members = r$n, rows = r)
}

genic_blocks <- function(sm) {
  w <- 2 * sm$p * (1 - sm$p)
  out <- list(A = crossprod(sm$alpha, w * sm$alpha), D = crossprod(sm$d, w^2 * sm$d))
  if (nrow(sm$E)) {
    wAA <- w[sm$pairs[, 1]] * w[sm$pairs[, 2]]
    out$AA <- crossprod(sm$E, wAA * sm$E)
  }
  out
}

dosages <- function(pop, ids, locus_ids) {
  d <- DBI::dbGetQuery(pop$db_conn, paste0(
    "SELECT id_ind, locus_id, CAST(SUM(allele) AS DOUBLE) AS g FROM ind_haplotype ",
    "GROUP BY 1, 2"))
  X <- matrix(NA_real_, length(ids), length(locus_ids))
  ok <- d$locus_id %in% locus_ids & d$id_ind %in% ids
  X[cbind(match(d$id_ind[ok], ids), match(d$locus_id[ok], locus_ids))] <- d$g[ok]
  X
}

ind_ids <- function(pop) DBI::dbGetQuery(pop$db_conn, "SELECT id_ind FROM ind_meta ORDER BY id_ind")$id_ind

tgv_component <- function(pop, trait, component, ids) {
  v <- DBI::dbGetQuery(pop$db_conn, paste0(
    "SELECT id_ind, tgv_value FROM ind_tgv WHERE trait_name = '", trait,
    "' AND component_name = '", component, "'"))
  out <- numeric(length(ids))
  out[match(v$id_ind, ids)] <- v$tgv_value
  out
}
tgv_total <- function(pop, trait, ids) {
  v <- DBI::dbGetQuery(pop$db_conn, paste0(
    "SELECT id_ind, tgv_total FROM ind_tgv_total WHERE trait_name = '", trait, "'"))
  v$tgv_total[match(ids, v$id_ind)]
}

# Assign genotypes: `G` is n x m dosages for the individuals `ids` at `loci`.
set_genotypes <- function(pop, ids, loci, G) {
  df <- data.frame(id_ind = rep(ids, length(loci)),
                   locus_name = rep(loci, each = length(ids)),
                   g = as.vector(G), stringsAsFactors = FALSE)
  duckdb::duckdb_register(pop$db_conn, "dge_set_g", df)
  on.exit(duckdb::duckdb_unregister(pop$db_conn, "dge_set_g"))
  DBI::dbExecute(pop$db_conn, paste0(
    "UPDATE ind_haplotype h SET allele = CASE WHEN h.parent_origin = 1 ",
    "THEN CAST(s.g >= 1 AS INTEGER) ELSE CAST(s.g = 2 AS INTEGER) END ",
    "FROM dge_set_g s JOIN genome_meta m USING (locus_name) ",
    "WHERE h.id_ind = s.id_ind AND h.locus_id = m.locus_id"))
  invisible(pop)
}

seed_state <- function() {
  if (exists(".Random.seed", envir = globalenv(), inherits = FALSE))
    get(".Random.seed", envir = globalenv()) else NULL
}

# A knowable refusal: an error matching `pattern`, the database unchanged, and
# .Random.seed unchanged (with a seed) or still absent (without one).
expect_refusal <- function(pop, call, pattern) {
  before <- db_state(pop)
  set.seed(99)
  s0 <- seed_state()
  expect_error(force(call()), pattern)
  expect_identical(seed_state(), s0)
  if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
    rm(".Random.seed", envir = globalenv())
  }
  expect_error(force(call()), pattern)
  expect_false(exists(".Random.seed", envir = globalenv(), inherits = FALSE))
  expect_identical(db_state(pop), before)
  set.seed(NULL)
}


# -- C1, C3, C9: the stored terms deliver the targets ------------------------

test_that("C1 / C3 / C9: genic stored terms give G_A, G_D, G_AA; A x A lands in interaction", {
  pop <- dge_pop("c1")
  on.exit(close_pop(pop))
  set.seed(1)
  quiet(define_genome_effects(gm_at(pop), c("T1", "T2"), G_A = GA2, G_D = GD2,
                              G_AA = GAA2))
  sm <- stored_model(pop, c("T1", "T2"))
  gb <- genic_blocks(sm)
  expect_lt(max(abs(gb$A - GA2)), 1e-10)
  expect_lt(max(abs(gb$D - GD2)), 1e-10)
  expect_lt(max(abs(gb$AA - GAA2)), 1e-10)
  # C9: every A x A term has two additive members.
  two <- sm$rows[sm$rows$n == 2L, ]
  expect_true(nrow(two) > 0 && all(two$contrast_name == "additive"))
  expect_equal(nrow(sm$E), 10L)
  # The passed targets were written.
  v <- tvc(pop)
  expect_setequal(unique(v$effect_name), c("additive", "dominance", "additive_by_additive"))
  # C3 and C9 through the evaluator.
  quiet(get_table(pop, "ind_meta") |> add_tgv(c("T1", "T2")))
  ids <- ind_ids(pop)
  X <- dosages(pop, ids, sm$loci)
  Z <- sweep(X, 2, 2 * sm$p)
  for (j in 1:2) {
    t <- c("T1", "T2")[j]
    expect_lt(max(abs(tgv_component(pop, t, "additive", ids) - Z %*% sm$alpha[, j])), 1e-10)
    inter <- (Z[, sm$pairs[, 1]] * Z[, sm$pairs[, 2]]) %*% sm$E[, j]
    expect_lt(max(abs(tgv_component(pop, t, "interaction", ids) - inter)), 1e-10)
  }
})

test_that("C2 / G3: for every genotype, functional g minus the evaluated total is mu", {
  for (anc in c("genic", "realised")) {
    pop <- dge_pop(paste0("c2", anc), n_loci = 3, n_ind = 54, n_hap = 40,
                   traits = "T1")
    ids <- ind_ids(pop)
    loci <- c("Locus_1", "Locus_2", "Locus_3")
    all27 <- as.matrix(expand.grid(0:2, 0:2, 0:2))
    set.seed(3)
    extra <- sapply(c(0.3, 0.6, 0.8), function(p) stats::rbinom(27, 2, p))
    G <- rbind(all27, extra)
    set_genotypes(pop, ids, loci, G)
    set.seed(4)
    quiet(define_genome_effects(get_table(pop, "genome_meta"), "T1", G_A = 1,
                                G_D = 0.3, G_AA = 0.2,
                                pairs = data.frame(locus_name_1 = "Locus_1",
                                                   locus_name_2 = "Locus_3"),
                                anchor = anc, base_tbl = get_table(pop, "ind_meta")))
    sm <- stored_model(pop, "T1")
    p <- sm$p
    # Functional coefficients from the stored statistical ones (plan §3).
    e <- sm$E[, 1]; pr <- sm$pairs
    a <- sm$alpha[, 1] - (1 - 2 * p) * sm$d[, 1]
    a[pr[, 1]] <- a[pr[, 1]] - e * (2 * p[pr[, 2]] - 1)
    a[pr[, 2]] <- a[pr[, 2]] - e * (2 * p[pr[, 1]] - 1)
    g <- (G - 1) %*% a + (G == 1) %*% sm$d[, 1] +
      ((G[, pr[, 1]] - 1) * (G[, pr[, 2]] - 1)) * e
    mu <- sum(a * (2 * p - 1) + 2 * p * (1 - p) * sm$d[, 1]) +
      e * (2 * p[pr[, 1]] - 1) * (2 * p[pr[, 2]] - 1)
    quiet(get_table(pop, "ind_meta") |> add_tgv("T1"))
    diff <- as.vector(g) - tgv_total(pop, "T1", ids)
    expect_lt(max(abs(diff - mu)), 1e-12)
    # The calibration's own guarantee, measured independently of the storage
    # choice: under "realised" the extractor recovers the targets on the base
    # (it would not if storage used the observed regression b).
    if (anc == "realised") {
      r <- suppressMessages(extract_genetic_variance(get_table(pop, "ind_meta"),
                                                     anchor = "realised"))
      got <- r$cov_value[match(c("additive", "dominance", "additive_by_additive"),
                               r$effect_name)]
      expect_lt(max(abs(got - c(1, 0.3, 0.2))), 1e-10)
    } else {
      gb <- genic_blocks(sm)
      expect_lt(max(abs(c(gb$A, gb$D, gb$AA) - c(1, 0.3, 0.2))), 1e-10)
    }
    close_pop(pop)
  }
})

test_that("G3: under genic the stored alpha is the calibrated B_alpha", {
  set.seed(5)
  m <- 30; k <- 2
  p <- stats::runif(m, 0.1, 0.9)
  an <- .na_aa_anchor(.na_anchors("genic", p = p), cbind(1:10, 11:20))
  cal <- .na_calibrate(an, GA2, GD2, GAA2, B_a = matrix(stats::rnorm(m * k), m),
                       z = matrix(stats::rnorm(m * k), m),
                       B_aa = matrix(stats::rnorm(10 * k), 10), pairs = cbind(1:10, 11:20))
  nm <- paste0("L", 1:m)
  for (j in 1:k) {
    st <- .noia_to_stored(stats::setNames(cal$B_a[, j], nm),
                          stats::setNames(cal$B_d[, j], nm),
                          data.frame(locus_1 = nm[1:10], locus_2 = nm[11:20],
                                     e = cal$B_aa[, j]),
                          stats::setNames(p, nm))
    expect_lt(max(abs(st$alpha - cal$B_alpha[, j])), 1e-12)
  }
})

test_that("G3: an additive target far below the dominance one is stored exactly (review r6 F1)", {
  # One QTL at p = 0.3, G_A / G_D = 1e-24: recovering alpha from the
  # functional a by cancellation stored a variance off by 2e-4.
  pop <- dge_pop("g3x")
  on.exit(close_pop(pop))
  set_genotypes(pop, ind_ids(pop), "Locus_1",
                matrix(rep(c(0, 0, 0, 1, 2), 40), ncol = 1))
  set.seed(2)
  quiet(define_genome_effects(gm_at(pop, 1), "T1", G_A = 1e-24, G_D = 1,
                              base_tbl = gen0(pop), warn_bounds = NULL))
  sm <- stored_model(pop, "T1")
  expect_equal(sm$p, 0.3)
  g <- genic_blocks(sm)
  expect_lt(abs(g$A[1, 1] / 1e-24 - 1), 1e-12)
  expect_lt(abs(g$D[1, 1] - 1), 1e-12)
})


# -- C4 (b): the additive-only route is define_additive_effects() ------------

test_that("C4 (b): additive-only targets write define_additive_effects()'s rows, k = 1, 2, both anchors", {
  targets <- list(0.7, GA2, 0, nm2(diag(c(1, 0))), nm2(matrix(1, 2, 2)))
  for (anc in c("genic", "realised")) for (G in targets) {
    traits <- if (length(G) == 1L) "T1" else c("T1", "T2")
    p1 <- dge_pop("c4a"); p2 <- dge_pop("c4b")
    set.seed(21)
    quiet(define_additive_effects(gm_at(p1), traits, G = G, anchor = anc,
                                  base_tbl = if (anc == "realised") gen0(p1)))
    s1 <- seed_state()
    set.seed(21)
    quiet(define_genome_effects(gm_at(p2), traits, G_A = G, anchor = anc,
                                base_tbl = if (anc == "realised") gen0(p2)))
    expect_identical(rows_noid(p2), rows_noid(p1))
    expect_identical(tvc(p2), tvc(p1))
    expect_identical(seed_state(), s1)
    close_pop(p1); close_pop(p2)
  }
})

test_that("C4 (b) / D5: explicit zero G_D and G_AA give the same additive rows and no other terms", {
  p1 <- dge_pop("c4z1"); p2 <- dge_pop("c4z2")
  on.exit({ close_pop(p1); close_pop(p2) })
  z2 <- nm2(matrix(0, 2, 2))
  set.seed(22)
  quiet(define_additive_effects(gm_at(p1), c("T1", "T2"), G = GA2))
  set.seed(22)
  quiet(define_genome_effects(gm_at(p2), c("T1", "T2"), G_A = GA2, G_D = z2, G_AA = z2))
  expect_identical(rows_noid(p2), rows_noid(p1))
  v <- tvc(p2)
  expect_true(all(v$cov_value[v$effect_name != "additive"] == 0))
  expect_setequal(unique(v$effect_name), c("additive", "dominance", "additive_by_additive"))
})

test_that("C4 (b): replacement, target persistence and custom owners", {
  p1 <- dge_pop("c4r1"); p2 <- dge_pop("c4r2")
  on.exit({ close_pop(p1); close_pop(p2) })
  for (p in list(p1, p2)) {
    quiet(define_genome_effect_terms(p, "T1", ad_terms("Locus_33", a = 0.5, d = 0.1, p = 0.5),
                                     effect_owner = "mine"))
  }
  set.seed(23)
  quiet(define_additive_effects(gm_at(p1), "T1", G = 0.5))
  set.seed(23)
  quiet(define_genome_effects(gm_at(p2), "T1", G_A = 0.5))
  # Re-run with the stored target: both replace their generated model.
  set.seed(24)
  quiet(define_additive_effects(gm_at(p1, 30), "T1"))
  set.seed(24)
  quiet(define_genome_effects(gm_at(p2, 30), "T1"))
  expect_identical(rows_noid(p2), rows_noid(p1))
  expect_identical(tvc(p2), tvc(p1))
  r <- rows_noid(p2)
  expect_equal(sum(r$effect_owner == "mine"), 2L)
  expect_equal(sum(r$effect_owner == "generated"), 30L)
})


# -- C5, C11: the floor, and atomicity ----------------------------------------

test_that("C5 / C11: below the floor errors naming it, and nothing is written", {
  pop <- dge_pop("c5")
  on.exit(close_pop(pop))
  set.seed(31)
  quiet(define_genome_effects(gm_at(pop), "T1", G_A = 1, G_D = 0.2))
  before <- db_state(pop)
  set.seed(32)
  err <- tryCatch(quiet(define_genome_effects(
    gm_at(pop), c("T2", "T3"), G_A = nm2(diag(c(1, 1e-6)), c("T2", "T3")),
    G_D = nm2(diag(c(0.2, 5)), c("T2", "T3")))), error = identity)
  expect_match(conditionMessage(err), "below the additive floor for this sampled architecture")
  expect_match(conditionMessage(err), "T3 = [0-9.e+-]+")
  expect_identical(db_state(pop), before)
})

test_that("C11: a failure inside the transaction restores every table", {
  pop <- dge_pop("c11")
  on.exit(close_pop(pop))
  quiet(define_genome_effect_terms(pop, "T1", ad_terms("Locus_33", a = 0.5, p = 0.5),
                                   effect_owner = "mine"))
  quiet(define_additive_effects(gm_at(pop), "T1", G = 0.5, line_name = "A", seed = 1))
  set.seed(41)
  quiet(define_genome_effects(gm_at(pop), "T2", G_A = 1, G_D = 0.2))
  before <- db_state(pop)
  calls <- 0L
  real <- .tvc_write_block
  testthat::local_mocked_bindings(.tvc_write_block = function(conn, ...) {
    calls <<- calls + 1L
    if (calls == 2L) stop("injected failure")
    real(conn, ...)
  }, .package = "tidybreed")
  set.seed(42)
  expect_error(quiet(define_genome_effects(gm_at(pop), "T3", G_A = 1, G_D = 0.2)),
               "injected failure")
  expect_equal(calls, 2L)
  expect_identical(db_state(pop), before)
})


# -- C7, C12: inbreeding depression -------------------------------------------

test_that("C7: one trait's inbreeding depression is exact, and two traits' is reported", {
  pop <- dge_pop("c7")
  on.exit(close_pop(pop))
  set.seed(51)
  quiet(define_genome_effects(gm_at(pop, 30), "T1", G_A = 1, G_D = 0.4,
                              inbreeding_depression = c(T1 = 1.5)))
  sm <- stored_model(pop, "T1")
  expect_lt(abs(sum(2 * sm$p * (1 - sm$p) * sm$d[, 1]) - 1.5), 1e-10)
  # One locus with A + D: ID / sqrt(V_D) is sign(d), so the attainable
  # depression is +/- sqrt(G_D) (the solver's degenerate branch).
  set.seed(52)
  quiet(define_genome_effects(get_table(pop, "genome_meta") |> dplyr::filter(locus_id == 5L),
                              "T2", G_A = 1, G_D = 0.1, inbreeding_depression = c(T2 = -sqrt(0.1))))
  sm2 <- stored_model(pop, "T2")
  w <- 2 * sm2$p * (1 - sm2$p)
  expect_lt(abs(sum(w * sm2$d[, 1]) + sqrt(0.1)), 1e-10)
  expect_lt(abs(sum(w^2 * sm2$d[, 1]^2) - 0.1), 1e-10)
  # Two traits: requested and delivered by name, approximate.
  pop2 <- dge_pop("c7b")
  on.exit(close_pop(pop2), add = TRUE)
  set.seed(53)
  expect_message(define_genome_effects(gm_at(pop2, 30), c("T1", "T3"),
                   G_A = nm2(diag(2), c("T1", "T3")),
                   G_D = nm2(diag(c(0.3, 0.4)), c("T1", "T3")),
                   inbreeding_depression = c(T3 = 0.3, T1 = 0.5)) |> suppressWarnings(),
                 "T1 requested 0.5, delivered [0-9.e-]+; T3 requested 0.3, delivered .*no closeness guarantee")
})

test_that("C12: the mean dominance value under exact inbred genotype counts is -F sum 2pq d", {
  pop <- dge_pop("c12", n_loci = 3, n_ind = 32, n_hap = 40, traits = "T1")
  on.exit(close_pop(pop))
  ids <- ind_ids(pop)
  # Counts n (p^2 + Fpq, 2pq(1 - F), q^2 + Fpq) at F = 0.5, n = 32.
  counts <- list(c(21, 6, 5), c(12, 8, 12), c(5, 6, 21))   # dosage 0, 1, 2
  set.seed(6)
  G <- sapply(counts, function(cn) sample(rep(0:2, cn)))
  set_genotypes(pop, ids, c("Locus_1", "Locus_2", "Locus_3"), G)
  set.seed(7)
  quiet(define_genome_effects(get_table(pop, "genome_meta"), "T1", G_A = 1, G_D = 0.5,
                              base_tbl = get_table(pop, "ind_meta"),
                              inbreeding_depression = c(T1 = 0.8)))
  quiet(get_table(pop, "ind_meta") |> add_tgv("T1"))
  sm <- stored_model(pop, "T1")
  expect_equal(sm$p, c(0.25, 0.5, 0.75))
  dom <- tgv_component(pop, "T1", "dominance", ids)
  expect_lt(abs(mean(dom) - (-0.5 * 0.8)), 1e-12)
  # (b) individual values from the genotypes and the stored terms.
  p <- sm$p; q <- 1 - p
  # The dominance contrast at dosage x of allele 1: -2p^2 (x = 0), 2pq (1), -2q^2 (2).
  Dc <- sapply(1:3, function(j) c(-2 * p[j]^2, 2 * p[j] * q[j], -2 * q[j]^2)[G[, j] + 1])
  own <- sweep(G, 2, 2 * p) %*% sm$alpha[, 1] + Dc %*% sm$d[, 1]
  expect_lt(max(abs(tgv_total(pop, "T1", ids) - own)), 1e-12)
})


# -- C8: one owner ------------------------------------------------------------

test_that("C8 (a): a re-run replaces common and line-scoped variants and counts them", {
  pop <- dge_pop("c8a")
  on.exit(close_pop(pop))
  quiet(define_additive_effects(gm_at(pop), "T1", G = 0.5, seed = 1,
                                base_tbl = get_table(pop, "founder_haplotypes")))
  quiet(define_additive_effects(gm_at(pop), "T1", seed = 2, line_name = "A"))
  quiet(define_additive_effects(gm_at(pop), "T1", seed = 3, line_name = "B",
                                base_tbl = get_table(pop, "founder_haplotypes")))
  set.seed(61)
  expect_message(define_genome_effects(gm_at(pop), "T1", G_D = 0.1) |> suppressWarnings(),
                 "Replaced 60 generated terms of trait 'T1' \\(40 of them line-scoped")
  r <- rows_noid(pop)
  expect_equal(sum(r$n_origins > 0), 0L)
})

test_that("C8 (b) / (c) / (d): define_additive_effects() refuses a non-additive generated model", {
  pop <- dge_pop("c8b")
  on.exit(close_pop(pop))
  quiet(define_genome_effect_terms(pop, "T1", ad_terms("Locus_33", a = 0.5, p = 0.5),
                                   effect_owner = "mine"))
  set.seed(62)
  quiet(define_genome_effects(gm_at(pop), "T1", G_A = 1, G_D = 0.2))
  before <- db_state(pop)
  expect_error(define_additive_effects(gm_at(pop), "T1", G = 1, line_name = "A"),
               "define_genome_effects\\(c\\(\"T1\"\\), trait_var_comp_tbl")
  expect_identical(db_state(pop), before)
  expect_false(any(!is.na(tvc(pop)$line_name)))
  # The additive-only re-run, then define_additive_effects() with the filter.
  set.seed(63)
  quiet(define_genome_effects(gm_at(pop), "T1", trait_var_comp_tbl =
    get_table(pop, "trait_var_comp") |> dplyr::filter(effect_name == "additive", is.na(line_name))))
  flt <- get_table(pop, "trait_var_comp") |> dplyr::filter(effect_name == "additive")
  expect_no_error(quiet(define_additive_effects(gm_at(pop), "T1", trait_var_comp_tbl = flt, seed = 4)))
  expect_error(define_additive_effects(gm_at(pop), "T1", seed = 4), "stored 'dominance' target")
  # (c) the custom owner is untouched throughout; (d) the writer refuses the owner.
  expect_equal(sum(rows_noid(pop)$effect_owner == "mine"), 1L)
  expect_error(define_genome_effect_terms(pop, "T1", ad_terms("Locus_1", a = 1, p = 0.5),
                                          effect_owner = "generated"),
               "reserved effect owner")
})


# -- C10, C14, C15 ------------------------------------------------------------

test_that("C10: non-autosomal QTL are refused", {
  pop <- open_pop(pop_name = "c10", db_name = ":memory:") |>
    define_genome(n_loci = 20, n_chr = 2, chr_names = c("1", "X"), chr_len_Mb = 100) |>
    define_chromosome("X", offspring_sex = "M", from_parent_1 = 0, from_parent_2 = 1) |>
    define_founder_haplotypes(n_haplotypes = 20, method = "fixed")
  on.exit(close_pop(pop))
  pop <- quiet(define_trait(pop, "T1"))
  expect_error(get_table(pop, "genome_meta") |> dplyr::filter(chr_name == "X") |>
                 define_genome_effects("T1", G_A = 1, G_D = 0.1),
               "assume diploid/autosomal QTL")
})

test_that("C14 / C15 (b): the same seed gives identical terms and a random matching", {
  run <- function() {
    pop <- dge_pop("c14")
    on.exit(close_pop(pop))
    set.seed(71)
    quiet(define_genome_effects(gm_at(pop, 21), c("T1", "T2"), G_A = GA2, G_D = GD2,
                                G_AA = GAA2))
    rows_noid(pop)
  }
  a <- run(); b <- run()
  expect_identical(a, b)
  pr <- a[a$trait_name == "T1" & a$contrast_name == "additive", ]
  pop <- dge_pop("c15b")
  on.exit(close_pop(pop))
  set.seed(71)
  expect_message(define_genome_effects(gm_at(pop, 21), c("T1", "T2"), G_A = GA2,
                                       G_D = GD2, G_AA = GAA2) |> suppressWarnings(),
                 "Drew 10 random A x A pairs: every QTL paired once \\(floor\\(21 / 2\\)\\)")
  sm <- stored_model(pop, c("T1", "T2"))
  expect_equal(anyDuplicated(as.vector(sm$pair_loci)), 0L)
  expect_identical(sm$pair_loci[order(sm$pair_loci[, 1], sm$pair_loci[, 2]), ], sm$pair_loci)
})

test_that("C15 (a): supplied pairs: bad keys error by name; a hub design is exact", {
  pop <- dge_pop("c15a")
  on.exit(close_pop(pop))
  pr <- function(a, b) data.frame(locus_name_1 = a, locus_name_2 = b)
  call <- function(pairs) function() define_genome_effects(gm_at(pop), "T1", G_A = 1,
                                                           G_AA = 0.2, pairs = pairs)
  expect_refusal(pop, call(pr("Locus_1", "Nope")), "not in genome_meta: 'Nope'")
  expect_refusal(pop, call(pr("Locus_1", "Locus_35")), "outside the filtered QTL set: 'Locus_35'")
  expect_refusal(pop, call(pr("Locus_2", "Locus_2")), "pairs a locus with itself: 'Locus_2'")
  expect_refusal(pop, call(pr(c("Locus_1", "Locus_2"), c("Locus_2", "Locus_1"))),
                 "repeats a pair \\(in either order\\): \\(Locus_2, Locus_1\\)")
  hub <- pr("Locus_3", c("Locus_4", "Locus_5", "Locus_6"))
  set.seed(72)
  quiet(define_genome_effects(gm_at(pop), "T1", G_A = 1, G_D = 0.1, G_AA = 0.2, pairs = hub))
  sm <- stored_model(pop, "T1")
  gb <- genic_blocks(sm)
  expect_equal(nrow(sm$E), 3L)
  expect_lt(abs(gb$A[1, 1] - 1), 1e-10)
  expect_lt(abs(gb$AA[1, 1] - 0.2), 1e-10)
})

test_that("C15 (c) / (d): n_pairs default, odd m, the maximum, and misuse", {
  pop <- dge_pop("c15c")
  on.exit(close_pop(pop))
  set.seed(73)
  expect_message(define_genome_effects(gm_at(pop, 7), "T1", G_A = 1, G_AA = 0.2) |>
                   suppressWarnings(), "Drew 3 random A x A pairs")
  sm <- stored_model(pop, "T1")
  expect_equal(length(unique(as.vector(sm$pair_loci))), 6L)
  expect_refusal(pop, function() define_genome_effects(gm_at(pop, 7), "T2", G_A = 1,
                                                       G_AA = 0.2, n_pairs = 4),
                 "`n_pairs` = 4 exceeds the 3 pairs.*pass `pairs`")
  expect_refusal(pop, function() define_genome_effects(gm_at(pop, 7), "T2", G_A = 1,
                   G_AA = 0.2, n_pairs = 2,
                   pairs = data.frame(locus_name_1 = "Locus_1", locus_name_2 = "Locus_2")),
                 "not both")
  expect_refusal(pop, function() define_genome_effects(gm_at(pop, 7), "T2", G_A = 1,
                                                       n_pairs = 2),
                 "`n_pairs` is used only with an 'additive_by_additive' target")
})


# -- C17, C18, G8: targets ----------------------------------------------------

test_that("C17: absent, stored, zero and passed blocks", {
  pop <- dge_pop("c17")
  on.exit(close_pop(pop))
  quiet(pop |> define_effect_cov_matrix("additive", 1, trait_name = "T1") |>
          define_effect_cov_matrix("dominance", 0.2, trait_name = "T1") |>
          define_effect_cov_matrix("additive_by_additive", 0.1, trait_name = "T1"))
  # Stored blocks are used.
  set.seed(81)
  quiet(define_genome_effects(gm_at(pop), "T1"))
  sm <- stored_model(pop, "T1")
  expect_true(nrow(sm$E) > 0 && any(sm$d != 0))
  # Filtering A x A out gives an A + D model; the stored A x A row stays.
  set.seed(81)
  quiet(define_genome_effects(gm_at(pop), "T1", trait_var_comp_tbl =
    get_table(pop, "trait_var_comp") |> dplyr::filter(effect_name != "additive_by_additive")))
  sm <- stored_model(pop, "T1")
  expect_equal(nrow(sm$E), 0L)
  expect_true(any(sm$d != 0))
  expect_true("additive_by_additive" %in% tvc(pop)$effect_name)
  # A passed block over a stored one, and a missing additive block.
  expect_refusal(pop, function() define_genome_effects(gm_at(pop), "T1", G_D = 0.2),
                 "'dominance' block is already stored for T1")
  expect_refusal(pop, function() define_genome_effects(gm_at(pop), "T2", G_D = 0.2),
                 "No 'additive' target for T2")
  expect_refusal(pop, function() define_genome_effects(gm_at(pop), "T2", G_A = 1,
                                                       inbreeding_depression = c(T2 = 1)),
                 "needs a 'dominance' target")
})

test_that("C17 / D2: adding a block never changes the earlier draws", {
  set.seed(1)
  mask <- matrix(TRUE, 30, 2, dimnames = list(NULL, c("T1", "T2")))
  nm <- paste0("L", 1:30); id <- 1:30
  draw <- function(has_d, has_aa, pairs = NULL, n_pairs = 5L) {
    set.seed(91)
    .dge_draw(mask, GA2, has_d, has_aa, pairs, n_pairs, nm, id)
  }
  a  <- draw(FALSE, FALSE)
  ad <- draw(TRUE, FALSE)
  adaa <- draw(TRUE, TRUE)
  sup <- draw(TRUE, TRUE, pairs = cbind(1:3, 4:6))
  expect_identical(ad$B_a, a$B_a)
  expect_identical(adaa$B_a, a$B_a)
  expect_identical(adaa$z, ad$z)
  expect_identical(sup$z, ad$z)
  expect_identical(sup$B_a, a$B_a)
})

test_that("C18: dimnames must equal trait_name; unnamed matrices are taken in order", {
  pop <- dge_pop("c18")
  on.exit(close_pop(pop))
  bad <- nm2(GA2, c("T2", "T1"))
  expect_refusal(pop, function() define_genome_effects(gm_at(pop), c("T1", "T2"), G_A = bad),
                 "is named \\(T2, T1\\) but `trait_name` is \\(T1, T2\\)")
  expect_refusal(pop, function() define_genome_effects(gm_at(pop), c("T1", "T2"), G_A = GA2,
                                                       G_D = nm2(GD2, c("A", "B"))),
                 "`G_D` is named")
  set.seed(92)
  quiet(define_genome_effects(gm_at(pop), c("T1", "T2"), G_A = unname(GA2)))
  v <- tvc(pop)
  expect_equal(v$cov_value[v$trait_name_1 == "T2" & v$trait_name_2 == "T2"], 2)
})

test_that("G8: one source per block, and the scope check", {
  pop <- dge_pop("g8")
  on.exit(close_pop(pop))
  quiet(pop |> define_effect_cov_matrix("additive", 1, trait_name = "T1") |>
          define_effect_cov_matrix("additive_by_additive", 0.1, trait_name = "T1"))
  add_only <- get_table(pop, "trait_var_comp") |> dplyr::filter(effect_name == "additive")
  # Passed G_D with the explicit additive rows: allowed.
  set.seed(93)
  expect_no_error(quiet(define_genome_effects(gm_at(pop), "T1", G_D = 0.2,
                                              trait_var_comp_tbl = add_only)))
  expect_true("dominance" %in% tvc(pop)$effect_name)
  # A passed block also in the explicit selection.
  expect_refusal(pop, function() define_genome_effects(gm_at(pop), "T1", G_AA = 0.3,
                   trait_var_comp_tbl = get_table(pop, "trait_var_comp")),
                 "two sources")
  # A passed block already stored but filtered out: the whole-table check.
  expect_refusal(pop, function() define_genome_effects(gm_at(pop), "T1", G_AA = 0.3,
                   trait_var_comp_tbl = add_only),
                 "'additive_by_additive' block is already stored for T1")
  # Explicit line rows for any block are refused by the scope check.
  quiet(define_effect_cov_matrix(pop, "dominance", 0.3, trait_name = "T2", line_name = "A"))
  quiet(define_effect_cov_matrix(pop, "additive", 1, trait_name = "T2"))
  expect_refusal(pop, function() define_genome_effects(gm_at(pop), "T2",
                   trait_var_comp_tbl = get_table(pop, "trait_var_comp") |>
                     dplyr::filter(trait_name_1 == "T2")),
                 "selects the line 'A' 'dominance' block")
  # Passed G_A with the stored D / A x A (default NULL table).
  quiet(define_effect_cov_matrix(pop, "dominance", 0.1, trait_name = "T3"))
  set.seed(94)
  expect_no_error(quiet(define_genome_effects(gm_at(pop), "T3", G_A = 1)))
  expect_true(any(stored_model(pop, "T3")$d != 0))
})


# -- C20: knowable refusals leave the seed and the database alone -------------

test_that("C20: every knowable refusal preserves .Random.seed and the database", {
  pop <- dge_pop("c20")
  on.exit(close_pop(pop))
  gm <- gm_at(pop)
  ref <- function(pattern, ...) {
    args <- list(...)
    expect_refusal(pop, function() do.call(define_genome_effects, c(list(gm), args)), pattern)
  }
  ref("dominance_degree_sd", "T1", G_A = 1, G_D = 0.1, dominance_degree_sd = -1)
  ref("warn_bounds", "T1", G_A = 1, warn_bounds = c(2, 1))
  ref("needs `dominance_degree_sd` > 0", "T1", G_A = 1, G_D = 0.1,
      dominance_degree_sd = 0, inbreeding_depression = c(T1 = 1))
  ref("names trait\\(s\\) not in this call: T9", "T1", G_A = 1, G_D = 0.1,
      inbreeding_depression = c(T9 = 1))
  ref("must name each of its traits once", "T1", G_A = 1, G_D = 0.1,
      inbreeding_depression = c(T1 = 1, T1 = 2))
  ref("must be positive semi-definite", "T1", G_A = -1)
  ref("whose dominance variance is 0", c("T1", "T2"), G_A = GA2,
      G_D = nm2(diag(c(0.4, 0))), inbreeding_depression = c(T2 = 1))
  ref("both 0, so every dominance effect drawn is 0", "T1", G_A = 1, G_D = 0.1,
      dominance_degree_mean = 0, dominance_degree_sd = 0)
  ref("positive-definite `G_A`", c("T1", "T2"), G_A = nm2(matrix(1, 2, 2)), G_D = GD2)
  # Explicit zero G_D with zero degree parameters is valid.
  set.seed(95)
  expect_no_error(quiet(define_genome_effects(gm, "T1", G_A = 1, G_D = 0,
                                              dominance_degree_mean = 0,
                                              dominance_degree_sd = 0)))
})


# -- G1: D5 and the rank note --------------------------------------------------

test_that("G1: a singular G_A is fine additive-only and refused, as a limit, with non-zero D", {
  pop <- dge_pop("g1")
  on.exit(close_pop(pop))
  one <- nm2(matrix(1, 2, 2))
  z2  <- nm2(matrix(0, 2, 2))
  set.seed(101)
  expect_no_error(quiet(define_genome_effects(gm_at(pop), c("T1", "T2"), G_A = one)))
  pop2 <- dge_pop("g1b")
  on.exit(close_pop(pop2), add = TRUE)
  set.seed(102)
  expect_no_error(quiet(define_genome_effects(gm_at(pop2), c("T1", "T2"), G_A = one,
                                              G_D = z2, G_AA = z2)))
  pop3 <- dge_pop("g1c")
  on.exit(close_pop(pop3), add = TRUE)
  err <- tryCatch(define_genome_effects(gm_at(pop3), c("T1", "T2"), G_A = one, G_D = GD2),
                  error = identity)
  expect_match(conditionMessage(err), "This release calibrates dominance or epistasis only")
  expect_false(grepl("infeasible|impossible", conditionMessage(err)))
  # Codex's hand-derived case is attainable (p = 0.5, B_a = [1 1; -1 -1],
  # B_d = I gives G_A = [1 1; 1 1], G_D = diag(1/16)): refused on purpose, as
  # this release's limit, not as an impossible target.
  expect_refusal(pop3, function() define_genome_effects(gm_at(pop3, 2), c("T1", "T2"),
                   G_A = one, G_D = nm2(diag(c(1, 1) / 16))),
                 "positive-definite `G_A`")
})

test_that("G1: the rank note names the reason, from all three entry points", {
  pop <- dge_pop("g1n")
  on.exit(close_pop(pop))
  set.seed(103)
  expect_message(define_genome_effects(gm_at(pop), c("T1", "T2"), G_A = GA2,
                                       G_D = nm2(diag(c(0.4, 0)))) |> suppressWarnings(),
                 "The 'dominance' target for T1, T2 is singular \\(rank 1 of 2\\): T2 has zero dominance variance")
  pop2 <- dge_pop("g1n2")
  on.exit(close_pop(pop2), add = TRUE)
  set.seed(104)
  expect_message(define_genome_effects(gm_at(pop2), c("T1", "T2"),
                                       G_A = nm2(matrix(c(1, 2, 2, 4), 2))) |> suppressWarnings(),
                 "genetic correlation \\+1 between T1 and T2")
  G3 <- matrix(c(1, 0.2, 1.2, 0.2, 1, 1.2, 1.2, 1.2, 2.4), 3,
               dimnames = list(c("T1", "T2", "T3"), c("T1", "T2", "T3")))
  expect_message(.qtl_rank_note(G3, "additive"),
                 "a linear dependency among traits T1, T2, T3")
  expect_silent(.qtl_rank_note(GA2, "additive"))
  pop3 <- dge_pop("g1n3")
  on.exit(close_pop(pop3), add = TRUE)
  expect_message(define_additive_effects(gm_at(pop3), c("T1", "T2"),
                                         G = nm2(matrix(1, 2, 2)), seed = 1) |> suppressWarnings(),
                 "genetic correlation \\+1 between T1 and T2")
  expect_message(define_effect_cov_matrix(pop3, "dominance", nm2(diag(c(0.4, 0)), c("T2", "T3"))),
                 "T3 has zero dominance variance")
})


# -- G2: determinism ------------------------------------------------------------

test_that("G2: identical rows across seeded runs, thread counts and restore_pop()", {
  run <- function(threads, name, db_name = ":memory:") {
    pop <- dge_pop(name, db_name = db_name)
    DBI::dbExecute(pop$db_conn, paste0("SET threads = ", threads))
    set.seed(111)
    quiet(define_genome_effects(gm_at(pop), c("T1", "T2"), G_A = GA2, G_D = GD2,
                                G_AA = GAA2, anchor = "realised", base_tbl = gen0(pop)))
    out <- rows_noid(pop)
    list(pop = pop, rows = out)
  }
  a <- run(1, "g2a"); b <- run(8, "g2b")
  expect_identical(a$rows, b$rows)
  close_pop(a$pop); close_pop(b$pop)
  tmp <- tempfile(fileext = ".duckdb")
  on.exit(unlink(tmp), add = TRUE)
  pop <- dge_pop("g2c", db_name = tmp)
  path <- tmp
  close_pop(pop)
  pop <- suppressMessages(restore_pop(path))
  set.seed(111)
  quiet(define_genome_effects(gm_at(pop), c("T1", "T2"), G_A = GA2, G_D = GD2,
                              G_AA = GAA2, anchor = "realised", base_tbl = gen0(pop)))
  expect_identical(rows_noid(pop), a$rows)
  close_pop(pop)
})


# -- G5, G7 ---------------------------------------------------------------------

test_that("G5: reversed name order gives identical rows; bad names are refused", {
  run <- function(ib) {
    pop <- dge_pop("g5")
    on.exit(close_pop(pop))
    set.seed(121)
    quiet(define_genome_effects(gm_at(pop, 30), c("T1", "T2"), G_A = GA2,
                                G_D = nm2(diag(c(0.3, 0.4))), inbreeding_depression = ib))
    rows_noid(pop)
  }
  expect_identical(run(c(T1 = 0.5, T2 = 0.2)), run(c(T2 = 0.2, T1 = 0.5)))
})

test_that("G7: the realised size guard counts the kept designs; knowable shortages refuse early", {
  pop <- dge_pop("g7", traits = c("T1", "T2", "T3", "T4"))
  on.exit(close_pop(pop))
  testthat::local_mocked_bindings(QTL_REALISED_MAX_CELLS = 4000, .package = "tidybreed")
  # Additive-only, n x m = 200 x 20 = 4000: at the limit, accepted, and equal
  # to define_additive_effects() including the RNG state afterwards.
  p1 <- dge_pop("g7a")
  on.exit(close_pop(p1), add = TRUE)
  set.seed(131)
  quiet(define_additive_effects(gm_at(p1), "T1", G = 1, anchor = "realised", base_tbl = gen0(p1)))
  s1 <- seed_state()
  set.seed(131)
  quiet(define_genome_effects(gm_at(pop), "T1", G_A = 1, anchor = "realised", base_tbl = gen0(pop)))
  expect_identical(rows_noid(pop), rows_noid(p1))
  expect_identical(seed_state(), s1)
  expect_refusal(pop, function() define_genome_effects(gm_at(pop, 20), "T2", G_A = 1,
                   G_D = 0.1, G_AA = 0.1, anchor = "realised", base_tbl = gen0(pop)),
                 "200 individuals x \\(20 additive \\+ 20 dominance \\+ 10 pair columns\\)")
  expect_refusal(pop, function() define_genome_effects(gm_at(pop, 1), "T2", G_A = 1,
                                                       G_AA = 0.1),
                 "Random A x A pairs need at least 2 QTL")
  expect_refusal(pop, function() define_genome_effects(gm_at(pop), c("T2", "T3"),
                   G_A = nm2(diag(2), c("T2", "T3")), G_AA = nm2(diag(2), c("T2", "T3")),
                   pairs = data.frame(locus_name_1 = "Locus_1", locus_name_2 = "Locus_2")),
                 "has rank 2 but there is 1 pair")
  # Only non-zero blocks keep a design (review r6 F2): the anchors handed to
  # the calibration are exactly the guard's count, for a zero or absent
  # dominance or A x A block.
  kept <- function(...) {
    got <- NULL
    testthat::local_mocked_bindings(.na_calibrate = function(anchors, ...) {
      got <<- anchors
      stop("captured")
    }, QTL_REALISED_MAX_CELLS = 1e6, .package = "tidybreed")
    expect_error(define_genome_effects(gm_at(pop, 20), "T3", G_A = 1,
                                       anchor = "realised", base_tbl = gen0(pop), ...),
                 "captured")
    cells <- function(a) if (is.null(a)) 0 else length(a$design)
    c(A = cells(got$A), D = cells(got$D), AA = cells(got$AA))
  }
  expect_identical(kept(G_AA = 0.1), c(A = 4000, D = 0, AA = 2000))
  expect_identical(kept(G_D = 0.1, G_AA = 0), c(A = 4000, D = 4000, AA = 0))
  expect_identical(kept(G_D = 0, G_AA = 0.1), c(A = 4000, D = 0, AA = 2000))
  # So the guard admits A + A x A at 200 x (20 + 10) = 6000 cells, where it
  # previously kept 8000 (the dominance design too).
  testthat::local_mocked_bindings(QTL_REALISED_MAX_CELLS = 6000, .package = "tidybreed")
  set.seed(132)
  expect_no_error(quiet(define_genome_effects(gm_at(pop, 20), "T4", G_A = 1,
                    G_AA = 0.1, anchor = "realised", base_tbl = gen0(pop))))
  expect_refusal(pop, function() define_genome_effects(gm_at(pop, 20), "T2", G_A = 1,
                   G_D = 0.1, G_AA = 0, anchor = "realised", base_tbl = gen0(pop),
                   n_pairs = 10),
                 "200 individuals x \\(20 additive \\+ 20 dominance \\+ 0 pair columns\\)")
  small <- get_table(pop, "ind_meta") |> dplyr::filter(id_ind %in% !!ind_ids(pop)[1:2])
  expect_refusal(pop, function() define_genome_effects(gm_at(pop, 5), c("T2", "T3"),
                   G_A = nm2(diag(2), c("T2", "T3")), anchor = "realised", base_tbl = small),
                 "at most 1 independent directions")
})


# -- Source test 16 / D6: the comparison with another population ---------------

test_that("16 / D6: an inbred base warns per block; warn_bounds = NULL silences it", {
  pop <- dge_pop("d6")
  on.exit(close_pop(pop))
  # Half the individuals fully inbred: the realised cohort is far from HWE.
  ids <- ind_ids(pop)
  inb <- ids[seq(1, length(ids), by = 2)]
  duckdb::duckdb_register(pop$db_conn, "d6_ids", data.frame(id_ind = inb))
  DBI::dbExecute(pop$db_conn, paste0(
    "UPDATE ind_haplotype AS h SET allele = s.allele FROM ind_haplotype AS s ",
    "WHERE h.id_ind = s.id_ind AND h.locus_id = s.locus_id AND h.parent_origin = 2 ",
    "AND s.parent_origin = 1 AND h.id_ind IN (SELECT id_ind FROM d6_ids)"))
  duckdb::duckdb_unregister(pop$db_conn, "d6_ids")
  set.seed(151)
  expect_warning(quiet_msg <- suppressMessages(define_genome_effects(
    gm_at(pop), "T1", G_A = 1, G_D = 0.3, anchor = "realised", base_tbl = gen0(pop))),
    "genic limit \\(expectation\\) sees covariances departing from the targets: .*dominance")
  set.seed(152)
  expect_warning(suppressMessages(define_genome_effects(
    gm_at(pop), "T2", G_A = 1, G_D = 0.3, base_tbl = gen0(pop))),
    "base individuals \\(observed\\) sees covariances")
  set.seed(153)
  expect_no_warning(suppressMessages(define_genome_effects(
    gm_at(pop), "T3", G_A = 1, G_D = 0.3, anchor = "realised", base_tbl = gen0(pop),
    warn_bounds = NULL)))
  # The founder pool compares the additive block only, as a message.
  pop2 <- dge_pop("d6b")
  on.exit(close_pop(pop2), add = TRUE)
  set.seed(154)
  expect_message(define_genome_effects(gm_at(pop2), "T1", G_A = 1, G_D = 0.3) |> suppressWarnings(),
                 "additive block only, dominance and additive-by-additive not compared")
})


# -- D3 (a): a QTL fixed in the base keeps its effects -------------------------

test_that("D3 (a): a QTL fixed in the base keeps non-zero effects that other populations see", {
  pop <- dge_pop("d3")
  on.exit(close_pop(pop))
  # Fix Locus_1 in the founder pool used as the base; individuals still segregate.
  DBI::dbExecute(pop$db_conn, paste0(
    "UPDATE founder_haplotypes SET allele = 0 WHERE locus_name = 'Locus_1'"))
  set.seed(141)
  quiet(define_genome_effects(gm_at(pop), "T1", G_A = 1, G_D = 0.2))
  sm <- stored_model(pop, "T1")
  expect_equal(sm$p[1], 0)
  expect_true(sm$alpha[1, 1] != 0 && sm$d[1, 1] != 0)
  # Zero variance at the base for its own contrasts; real in the individuals.
  expect_equal(2 * sm$p[1] * (1 - sm$p[1]), 0)
  quiet(get_table(pop, "ind_meta") |> add_tgv("T1"))
  ids <- ind_ids(pop)
  X <- dosages(pop, ids, sm$loci)
  expect_gt(stats::var(X[, 1]), 0)
  others <- sweep(X[, -1], 2, 2 * sm$p[-1]) %*% sm$alpha[-1, 1]
  expect_gt(stats::var(tgv_component(pop, "T1", "additive", ids) - others), 0)
})
