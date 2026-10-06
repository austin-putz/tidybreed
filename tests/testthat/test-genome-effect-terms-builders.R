# Step 4a of plans/import_qtl_effect_methods_phase_4_plan.md: the fixed builder
# column set and collision-free term ids (gate C16), aa_terms(), and the NOIA
# conversion pair (gates N1-N3).
#
# Every equivalence here is checked against the real evaluator: terms are
# written through define_genome_effect_terms() and evaluated by add_tgv(), and
# the test computes its expectation from the individuals' dosages itself.

# A small panel with arbitrary locus names, so the id-collision gates can use
# names that contain the encoder's own delimiters.
tb_pop <- function(loci = c("A", "B", "AxB", "C", "BxC", "x:1", "y|2", "z#3"),
                   n_ind = 120, seed = 7) {
  set.seed(seed)
  open_pop(pop_name = "tb", db_name = ":memory:") |>
    define_genome(n_loci = length(loci), n_chr = 1, chr_len_Mb = 10,
                  locus_names = loci) |>
    define_founder_haplotypes(n_haplotypes = 60) |>
    get_table("founder_haplotypes") |>
    add_founders(n_males = n_ind / 2, n_females = n_ind / 2, line_name = "A")
}

# Dosage of allele 1 per individual (rows) and locus name (columns).
tb_dosage <- function(pop) {
  d <- DBI::dbGetQuery(pop$db_conn, paste0(
    "SELECT h.id_ind, g.locus_name, CAST(SUM(h.allele) AS INTEGER) AS g ",
    "FROM ind_haplotype h JOIN genome_meta g USING (locus_id) ",
    "GROUP BY 1, 2 ORDER BY 1, 2"))
  x <- stats::xtabs(g ~ id_ind + locus_name, data = d)
  as.matrix(as.data.frame.matrix(x))
}

# Write `terms` for a fresh trait and return its evaluated total, by id_ind.
tb_total <- function(pop, trait, terms, owner = "custom") {
  pop <- define_trait(pop, trait)
  pop <- suppressMessages(
    define_genome_effect_terms(pop, trait, terms, effect_owner = owner))
  suppressMessages(pop |> get_table("ind_meta") |> add_tgv(trait))
  v <- DBI::dbGetQuery(pop$db_conn, paste0(
    "SELECT id_ind, tgv_total AS tgv_value FROM ind_tgv_total WHERE trait_name = '", trait,
    "' ORDER BY id_ind"))
  stats::setNames(v$tgv_value, v$id_ind)
}

tb_n_terms <- function(pop, trait) {
  DBI::dbGetQuery(pop$db_conn, paste0(
    "SELECT COUNT(*) AS n FROM genome_effects WHERE trait_name = '", trait,
    "'"))$n
}

BUILDER_COLS <- c("term_id", "locus_name", "contrast_name", "center_value",
                  "copy_count_value", "dosage_value", "genome_value",
                  "effect_name")
BUILDER_TYPES <- c("character", "character", "character", "numeric",
                   "integer", "integer", "numeric", "character")


# -- C16: one column set ------------------------------------------------------

test_that("every builder returns the same columns, order and types (C16)", {
  outs <- list(
    ad_f  = suppressMessages(ad_terms(c("A", "B"), a = c(1, 2), d = c(0.5, 0),
                                      p = c(0.3, 0.6))),
    ad_c  = suppressMessages(ad_terms(c("A", "B"), a = c(1, 2), d = c(0.5, 0),
                                      p = c(0.3, 0.6), coding = "cockerham")),
    aa_f  = suppressMessages(aa_terms("A", "B", e = 1, p_1 = 0.3, p_2 = 0.6)),
    aa_c  = aa_terms("A", "B", e = 1, p_1 = 0.3, p_2 = 0.6,
                     coding = "cockerham"),
    geno  = genotype_terms(data.frame(A = 0:2), c(0, 1, 2)))
  for (nm in names(outs)) {
    expect_identical(names(outs[[nm]]), BUILDER_COLS, info = nm)
    expect_identical(unname(vapply(outs[[nm]], function(x) class(x)[1], "")),
                     BUILDER_TYPES, info = nm)
  }
  # rbind() in any order gives the same column set and types.
  for (i in seq_along(outs)) for (j in seq_along(outs)) {
    b <- rbind(outs[[i]], outs[[j]])
    expect_identical(names(b), BUILDER_COLS)
    expect_identical(unname(vapply(b, function(x) class(x)[1], "")),
                     BUILDER_TYPES)
  }
})


# -- C16: collision-free term ids (D5) ----------------------------------------

test_that("surfaces on c('A','B') and on 'AxB' stay separate terms (D5)", {
  pop <- tb_pop()
  on.exit(close_pop(pop), add = TRUE)
  s1 <- genotype_terms(expand.grid(A = 0:2, B = 0:2), rep(1, 9))
  s2 <- genotype_terms(data.frame(AxB = 0:2), rep(1, 3))
  expect_length(intersect(s1$term_id, s2$term_id), 0L)

  t1 <- tb_total(pop, "S1", s1)
  t2 <- tb_total(pop, "S2", s2)
  t3 <- tb_total(pop, "S3", rbind(s1, s2))
  # 9 two-member terms + 3 one-member terms, not merged into 3-member ones.
  expect_equal(tb_n_terms(pop, "S3"), 12)
  ord <- DBI::dbGetQuery(pop$db_conn, paste0(
    "SELECT effect_order, COUNT(*) AS n FROM genome_effect_terms ",
    "WHERE trait_name = 'S3' GROUP BY 1 ORDER BY 1"))
  expect_equal(ord$n, c(3, 9))
  expect_equal(t3, t1 + t2, tolerance = 1e-12)
})

test_that("aa_terms() pairs (A, BxC) and (AxB, C) stay separate terms (D5)", {
  pop <- tb_pop()
  on.exit(close_pop(pop), add = TRUE)
  p1 <- aa_terms("A", "BxC", e = 1, p_1 = 0.5, p_2 = 0.5, report = FALSE)
  p2 <- aa_terms("AxB", "C", e = 1, p_1 = 0.5, p_2 = 0.5, report = FALSE)
  expect_length(intersect(p1$term_id, p2$term_id), 0L)
  t1 <- tb_total(pop, "P1", p1)
  t2 <- tb_total(pop, "P2", p2)
  t3 <- tb_total(pop, "P3", rbind(p1, p2))
  expect_equal(tb_n_terms(pop, "P3"), 2)
  expect_equal(t3, t1 + t2, tolerance = 1e-12)
})

test_that("locus names containing the encoder's delimiters never collide (D5)", {
  pop <- tb_pop()
  on.exit(close_pop(pop), add = TRUE)
  odd <- c("x:1", "y|2", "z#3")
  tt <- rbind(
    suppressMessages(ad_terms(odd, a = c(0.2, -0.3, 0.4), d = c(0.1, 0.2, 0),
                              p = 0.5)),
    aa_terms("x:1", "y|2", e = 0.5, p_1 = 0.5, p_2 = 0.5, report = FALSE),
    aa_terms("y|2", "z#3", e = -0.4, p_1 = 0.5, p_2 = 0.5, report = FALSE),
    genotype_terms(stats::setNames(expand.grid(0:2, 0:2), c("x:1", "z#3")),
                   seq(0.1, 0.9, by = 0.1)),
    genotype_terms(stats::setNames(data.frame(0:2), "y|2"), c(1, 0, 3)))
  # Distinct terms have distinct ids; equal ids only within one term.
  per_id <- tapply(tt$locus_name, tt$term_id, function(l) !anyDuplicated(l))
  expect_true(all(per_id))
  n_expected <- 5 + 2 + 9 + 2   # the y|2 surface skips the het state ad_terms writes
  expect_equal(length(unique(tt$term_id)), n_expected)

  parts <- list(
    suppressMessages(ad_terms(odd, a = c(0.2, -0.3, 0.4), d = c(0.1, 0.2, 0),
                              p = 0.5)),
    aa_terms("x:1", "y|2", e = 0.5, p_1 = 0.5, p_2 = 0.5, report = FALSE),
    aa_terms("y|2", "z#3", e = -0.4, p_1 = 0.5, p_2 = 0.5, report = FALSE),
    genotype_terms(stats::setNames(expand.grid(0:2, 0:2), c("x:1", "z#3")),
                   seq(0.1, 0.9, by = 0.1)),
    genotype_terms(stats::setNames(data.frame(0:2), "y|2"), c(1, 0, 3)))
  sep <- Reduce(`+`, lapply(seq_along(parts), function(i)
    tb_total(pop, paste0("Q", i), parts[[i]])))
  all_t <- tb_total(pop, "QALL", tt)
  expect_equal(tb_n_terms(pop, "QALL"), n_expected)
  expect_equal(all_t, sep, tolerance = 1e-12)
})

test_that("one locus used by every builder gives distinct ids per term (D5)", {
  ids <- c(
    suppressMessages(ad_terms("A", a = 1, d = 1, p = 0.5))$term_id,
    aa_terms("A", "B", e = 1, p_1 = 0.5, p_2 = 0.5, report = FALSE)$term_id[1],
    genotype_terms(data.frame(A = 0:2), c(1, 2, 3))$term_id)
  expect_equal(length(unique(ids)), length(ids))
})

test_that("two surfaces over the same loci are refused by the writer, not merged", {
  pop <- tb_pop()
  on.exit(close_pop(pop), add = TRUE)
  pop <- define_trait(pop, "R")
  s1 <- genotype_terms(expand.grid(A = 0:2, B = 0:2), rep(1, 9))
  s2 <- genotype_terms(expand.grid(B = 0:2, A = 0:2), rep(2, 9))
  expect_error(suppressMessages(
    define_genome_effect_terms(pop, "R", rbind(s1, s2))))
  expect_equal(tb_n_terms(pop, "R"), 0)
})


# -- aa_terms() ---------------------------------------------------------------

test_that("aa_terms() canonicalises pair order and moves p with its locus", {
  a <- aa_terms("B", "A", e = 0.5, p_1 = 0.2, p_2 = 0.7, coding = "cockerham")
  b <- aa_terms("A", "B", e = 0.5, p_1 = 0.7, p_2 = 0.2, coding = "cockerham")
  expect_identical(a, b)
  expect_equal(a$locus_name, c("A", "B"))
  expect_equal(a$center_value, c(0.7, 0.2))
  expect_true(all(a$contrast_name == "additive"))
  f <- suppressMessages(aa_terms("B", "A", e = 0.5, p_1 = 0.2, p_2 = 0.7))
  expect_equal(f$center_value, c(0.5, 0.5))
})

test_that("aa_terms() validates before dropping zero pairs (D7)", {
  expect_error(aa_terms("A", "A", e = 1, p_1 = 0.5, p_2 = 0.5),
               "same locus twice")
  expect_error(aa_terms(c("A", "B"), c("B", "A"), e = 1, p_1 = 0.5, p_2 = 0.5),
               "repeated")
  expect_error(aa_terms("A", "B", e = 1), "required")
  expect_error(aa_terms("A", c("B", "C"), e = 1, p_1 = 0.5, p_2 = 0.5),
               "same length")
  expect_error(aa_terms("A", "B", e = 1, p_1 = 1.5, p_2 = 0.5), "between 0")
  expect_error(aa_terms("A", "B", e = 0, p_1 = 0.5, p_2 = 0.5),
               "no term to write")
  # A malformed input is reported even when its pair has e = 0.
  expect_error(aa_terms(c("A", "C"), c("B", "D"), e = c(1, 0),
                        p_1 = c(0.5, NA), p_2 = 0.5), "free of NA")
  expect_error(aa_terms(c("A", "C"), c("B", "C"), e = c(1, 0),
                        p_1 = 0.5, p_2 = 0.5), "same locus twice")
  expect_error(aa_terms("A", "B", e = Inf, p_1 = 0.5, p_2 = 0.5), "finite")
  kept <- aa_terms(c("A", "C"), c("B", "D"), e = c(1, 0), p_1 = 0.5,
                   p_2 = 0.5, report = FALSE)
  expect_equal(unique(kept$locus_name), c("A", "B"))
})

test_that("aa_terms() reports each pair's share of mu under functional coding only", {
  mu <- 0.8 * (2 * 0.3 - 1) * (2 * 0.6 - 1)
  expect_message(aa_terms("A", "B", e = 0.8, p_1 = 0.3, p_2 = 0.6),
                 format(round(mu, 6)), fixed = TRUE)
  expect_message(aa_terms("A", "B", e = 0.8, p_1 = 0.3, p_2 = 0.6),
                 "written to no table")
  expect_silent(aa_terms("A", "B", e = 0.8, p_1 = 0.3, p_2 = 0.6,
                         coding = "cockerham"))
})


# -- N1: every conversion row against the real evaluator ----------------------

# The functional value of a .stored_to_functional() result, plus kappa, from
# the dosages (locus names as columns).
tb_functional <- function(pop, f, G) {
  nm <- DBI::dbGetQuery(pop$db_conn,
    "SELECT locus_id, locus_name FROM genome_meta")
  name_of <- function(id) nm$locus_name[match(as.integer(id), nm$locus_id)]
  v <- rep(f$kappa, nrow(G))
  for (l in names(f$a)) {
    g <- G[, name_of(l)]
    v <- v + f$a[[l]] * (g - 1) + f$d[[l]] * (g == 1)
  }
  for (r in seq_len(nrow(f$pairs))) {
    gk <- G[, name_of(f$pairs$locus_id_1[r])]
    gl <- G[, name_of(f$pairs$locus_id_2[r])]
    v <- v + f$pairs$e[r] * (gk - 1) * (gl - 1)
  }
  stats::setNames(as.numeric(v), rownames(G))
}

test_that("each stored term shape equals its functional form plus kappa (N1)", {
  pop <- tb_pop(loci = c("A", "B", "C", "D"), n_ind = 200, seed = 3)
  on.exit(close_pop(pop), add = TRUE)
  G <- tb_dosage(pop)
  # Every single-locus state and all nine pair states occur.
  expect_true(all(0:2 %in% G[, "A"]))
  expect_equal(nrow(unique(G[, c("A", "B")])), 9L)

  shapes <- list(
    additive  = data.frame(locus_name = "A", contrast_name = "additive",
                           center_value = 0.37, genome_value = 0.8),
    dominance = data.frame(locus_name = "A", contrast_name = "dominance",
                           center_value = 0.23, genome_value = -1.3),
    ind_het   = data.frame(locus_name = "A", contrast_name = "indicator",
                           copy_count_value = 2L, dosage_value = 1L,
                           genome_value = 0.6),
    ind_hom2  = data.frame(locus_name = "A", contrast_name = "indicator",
                           copy_count_value = 2L, dosage_value = 2L,
                           genome_value = 0.9),
    ind_hom0  = data.frame(locus_name = "A", contrast_name = "indicator",
                           copy_count_value = 2L, dosage_value = 0L,
                           genome_value = -0.7),
    aa        = data.frame(term_id = 1L, locus_name = c("A", "B"),
                           contrast_name = "additive",
                           center_value = c(0.31, 0.72), genome_value = 1.1))
  for (nm in names(shapes)) {
    tr <- paste0("N1_", nm)
    got <- tb_total(pop, tr, shapes[[nm]])
    m <- .gev_read_model(pop$db_conn, tr)
    f <- .stored_to_functional(m$terms, m$members)
    want <- tb_functional(pop, f, G)
    expect_equal(got, want[names(got)], tolerance = 1e-12, info = nm)
  }
})

test_that("a two-owner model with every shape converts as a whole (N1)", {
  pop <- tb_pop(loci = c("A", "B", "C", "D"), n_ind = 200, seed = 3)
  on.exit(close_pop(pop), add = TRUE)
  G <- tb_dosage(pop)
  pop <- define_trait(pop, "MIX")
  pop <- suppressMessages(define_genome_effect_terms(pop, "MIX", rbind(
    suppressMessages(ad_terms(c("A", "B"), a = c(0.4, -0.2), d = c(0.3, 0.1),
                              p = c(0.2, 0.65))),
    aa_terms("A", "C", e = 0.7, p_1 = 0.4, p_2 = 0.45, coding = "cockerham")),
    effect_owner = "one"))
  pop <- suppressMessages(define_genome_effect_terms(pop, "MIX", rbind(
    suppressMessages(ad_terms(c("A", "D"), a = c(0.1, 0.5), d = c(-0.2, 0.25),
                              p = c(0.33, 0.58), coding = "cockerham")),
    aa_terms("B", "D", e = -0.6, p_1 = 0.5, p_2 = 0.5, report = FALSE),
    genotype_terms(data.frame(C = 0:2), c(0.2, -0.4, 0.9))),
    effect_owner = "two"))
  suppressMessages(pop |> get_table("ind_meta") |> add_tgv("MIX"))
  got <- DBI::dbGetQuery(pop$db_conn, paste0(
    "SELECT id_ind, tgv_total AS tgv_value FROM ind_tgv_total WHERE trait_name = 'MIX' ",
    "ORDER BY id_ind"))
  m <- .gev_read_model(pop$db_conn, "MIX")
  f <- .stored_to_functional(m$terms, m$members)
  want <- tb_functional(pop, f, G)
  expect_equal(got$tgv_value, unname(want[got$id_ind]), tolerance = 1e-12)
})


# -- N2: round trip -----------------------------------------------------------

test_that(".noia_to_stored() and .stored_to_functional() round-trip (N2)", {
  pop <- tb_pop(loci = c("A", "B", "C", "D"), n_ind = 150, seed = 11)
  on.exit(close_pop(pop), add = TRUE)
  set.seed(101)
  loci <- c("A", "B", "C", "D")
  a <- stats::setNames(stats::rnorm(4), loci)
  d <- stats::setNames(stats::rnorm(4), loci)
  p <- stats::setNames(stats::runif(4, 0.1, 0.9), loci)
  pairs <- data.frame(locus_1 = c("A", "A", "B"), locus_2 = c("B", "C", "D"),
                      e = stats::rnorm(3))

  stat <- .noia_to_stored(a, d, pairs, p)
  # mu is the functional model's HWE + LE mean (plan §3), computed here.
  mu_src <- sum(a * (2 * p - 1) + 2 * p * (1 - p) * d) +
    sum(pairs$e * (2 * p[pairs$locus_1] - 1) * (2 * p[pairs$locus_2] - 1))
  expect_equal(stat$mu, unname(mu_src), tolerance = 1e-12)

  s_tot <- tb_total(pop, "STAT", .noia_terms(stat))
  m <- .gev_read_model(pop$db_conn, "STAT")
  f <- .stored_to_functional(m$terms, m$members)
  ids <- DBI::dbGetQuery(pop$db_conn,
    "SELECT locus_id, locus_name FROM genome_meta")
  nm <- function(id) ids$locus_name[match(as.integer(id), ids$locus_id)]
  expect_equal(unname(f$a[ids$locus_id[match(loci, ids$locus_name)]]),
               unname(a), tolerance = 1e-12)
  expect_equal(unname(f$d[ids$locus_id[match(loci, ids$locus_name)]]),
               unname(d), tolerance = 1e-12)
  got_pairs <- paste(nm(f$pairs$locus_id_1), nm(f$pairs$locus_id_2))
  expect_equal(f$pairs$e[match(paste(pairs$locus_1, pairs$locus_2), got_pairs)],
               pairs$e, tolerance = 1e-12)
  expect_equal(f$kappa, -stat$mu, tolerance = 1e-12)

  # The same model in functional coding differs by mu alone.
  fun <- rbind(suppressMessages(ad_terms(loci, a = unname(a), d = unname(d),
                                         p = unname(p))),
               aa_terms(pairs$locus_1, pairs$locus_2, e = pairs$e,
                        p_1 = unname(p[pairs$locus_1]),
                        p_2 = unname(p[pairs$locus_2]), report = FALSE))
  f_tot <- tb_total(pop, "FUN", fun)
  expect_equal(unname(f_tot - s_tot), rep(stat$mu, length(f_tot)),
               tolerance = 1e-10)
})

test_that("ad_terms()' reported mu equals .noia_to_stored()'s with no pairs", {
  a <- c(A = 0.4, B = -0.3); d <- c(A = 0.2, B = 0.5); p <- c(A = 0.3, B = 0.8)
  stat <- .noia_to_stored(a, d, data.frame(locus_1 = character(0),
                                           locus_2 = character(0),
                                           e = numeric(0)), p)
  mu_ad <- sum(a * (2 * p - 1) + 2 * p * (1 - p) * d)
  expect_equal(stat$mu, mu_ad, tolerance = 1e-15)
  expect_message(ad_terms(names(a), a = unname(a), d = unname(d),
                          p = unname(p)),
                 format(round(mu_ad, 6)), fixed = TRUE)
})


# -- N3: Cockerham aa_terms() -------------------------------------------------

test_that("Cockerham aa_terms() writes centres p and evaluates e(g_k-2p_k)(g_l-2p_l) (N3)", {
  pop <- tb_pop(loci = c("A", "B", "C", "D"), n_ind = 120, seed = 5)
  on.exit(close_pop(pop), add = TRUE)
  G <- tb_dosage(pop)
  tt <- aa_terms("B", "A", e = 0.9, p_1 = 0.62, p_2 = 0.27,
                 coding = "cockerham")
  got <- tb_total(pop, "N3", tt)
  cen <- DBI::dbGetQuery(pop$db_conn, paste0(
    "SELECT l.locus_name, m.center_value FROM genome_effect_members m ",
    "JOIN genome_effect_loci l USING (id_genome_effect, member_slot) ",
    "JOIN genome_effects e USING (id_genome_effect) ",
    "WHERE e.trait_name = 'N3' ORDER BY l.locus_name"))
  expect_equal(cen$center_value, c(0.27, 0.62))
  want <- 0.9 * (G[, "A"] - 2 * 0.27) * (G[, "B"] - 2 * 0.62)
  expect_equal(got, as.numeric(want[names(got)]) |>
                 stats::setNames(names(got)), tolerance = 1e-12)
})
