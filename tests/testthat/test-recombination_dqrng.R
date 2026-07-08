# Stage 2: dqrng per-gamete-stream recombination kernel.
# Locks (a) the log-accumulation inversion Poisson sampler, (b) the seed API,
# (c) crossover storage + its consistency with the haplotype it produced, and
# (d) that store_crossovers has no effect on the allele stream.

# ---- inversion Poisson sampler --------------------------------------------

test_that(".draw_poisson_dqrng matches Poisson mean/variance and handles edges", {
  dqrng::dqRNGkind("Xoroshiro128++")

  # lambda = 0 -> always 0, draws nothing.
  dqrng::dqset.seed(c(1L, 1L))
  expect_true(all(vapply(1:100, function(i) .draw_poisson_dqrng(0), integer(1)) == 0L))

  for (lam in c(0.3, 1, 5)) {
    dqrng::dqset.seed(c(11L, lam * 1000))
    x <- vapply(1:20000, function(i) .draw_poisson_dqrng(lam), integer(1))
    expect_equal(mean(x), lam, tolerance = 0.05)
    expect_equal(stats::var(x), lam, tolerance = 0.10)
    expect_true(min(x) >= 0L)
  }

  # Above the documented ceiling it errors rather than spinning O(lambda).
  expect_error(.draw_poisson_dqrng(31), "ceiling")

  # Reproducible for the same stream.
  dqrng::dqset.seed(c(3L, 9L)); a <- vapply(1:50, function(i) .draw_poisson_dqrng(1), integer(1))
  dqrng::dqset.seed(c(3L, 9L)); b <- vapply(1:50, function(i) .draw_poisson_dqrng(1), integer(1))
  expect_identical(a, b)
})

# ---- seed API --------------------------------------------------------------

# Build a fixed pop with founders seeded deterministically; return the pop plus a
# fresh mating design. Founders use base-R sample() (seeded here), so the pop is
# identical every call regardless of surrounding RNG state.
build_seed_pop <- function() {
  set.seed(2024)
  pop <- open_pop(pop_name = "sd", db_name = ":memory:") |>
    define_genome(n_loci = 150, n_chr = 3, chr_len_Mb = 100) |>
    define_founder_haplotypes(n_haplotypes = 30, line_name = "A", method = "uniform")
  pop |> get_table("founder_haplotypes") |>
    add_founders(n_males = 5, n_females = 5, line_name = "A")
}
seed_matings <- function() tibble::tibble(
  id_parent_1 = rep(paste0("A_", 1:5), length.out = 12),
  id_parent_2 = rep(paste0("A_", 6:10), length.out = 12),
  sex = rep(c("M", "F"), length.out = 12), line_name = "A")

off_haps <- function(pop) {
  h <- DBI::dbGetQuery(pop$db_conn,
    "SELECT id_ind, parent_origin, locus_id, allele FROM ind_haplotype
     WHERE id_ind LIKE 'A\\_1%' ESCAPE '\\' OR id_ind LIKE 'A\\_2%' ESCAPE '\\'
     ORDER BY id_ind, parent_origin, locus_id")
  h
}

test_that("explicit seed is reproducible and independent of surrounding base-R state", {
  p1 <- build_seed_pop(); set.seed(1)
  p1 <- suppressMessages(add_offspring(p1, seed_matings(), seed = 42L)); h1 <- off_haps(p1); close_pop(p1)

  p2 <- build_seed_pop(); set.seed(999)   # different outer base-R state
  p2 <- suppressMessages(add_offspring(p2, seed_matings(), seed = 42L)); h2 <- off_haps(p2); close_pop(p2)

  expect_equal(h1, h2)                     # explicit seed bypasses base-R

  # A different explicit seed changes the offspring stream.
  p3 <- build_seed_pop()
  p3 <- suppressMessages(add_offspring(p3, seed_matings(), seed = 7L)); h3 <- off_haps(p3); close_pop(p3)
  expect_false(isTRUE(all.equal(h1, h3)))
})

test_that("seed = NULL is reproducible under an upstream set.seed()", {
  run <- function() {
    set.seed(555)
    pop <- open_pop(pop_name = "sn", db_name = ":memory:") |>
      define_genome(n_loci = 120, n_chr = 2, chr_len_Mb = 100) |>
      define_founder_haplotypes(n_haplotypes = 24, line_name = "A", method = "uniform")
    pop <- pop |> get_table("founder_haplotypes") |>
      add_founders(n_males = 4, n_females = 4, line_name = "A")
    pop <- suppressMessages(add_offspring(pop, tibble::tibble(
      id_parent_1 = rep(paste0("A_", 1:4), 2), id_parent_2 = rep(paste0("A_", 5:8), 2),
      sex = rep(c("M", "F"), 4), line_name = "A")))
    bs <- attr(pop, "base_seed")
    h  <- DBI::dbGetQuery(pop$db_conn, "SELECT * FROM ind_haplotype ORDER BY id_ind,parent_origin,locus_id")
    close_pop(pop); list(h = h, bs = bs)
  }
  a <- run(); b <- run()
  expect_equal(a$h, b$h)
  expect_equal(a$bs, b$bs)                 # same base seed drawn from the same base-R stream
})

# ---- crossover storage -----------------------------------------------------

test_that("store_crossovers populates ind_crossover without perturbing ind_haplotype", {
  run <- function(sx) {
    set.seed(88)
    pop <- open_pop(pop_name = "cx", db_name = ":memory:") |>
      define_genome(n_loci = 200, n_chr = 3, chr_len_Mb = 200) |>
      define_founder_haplotypes(n_haplotypes = 30, line_name = "A", method = "uniform")
    pop <- pop |> get_table("founder_haplotypes") |>
      add_founders(n_males = 5, n_females = 5, line_name = "A")
    pop <- suppressMessages(add_offspring(pop, tibble::tibble(
      id_parent_1 = rep(paste0("A_", 1:5), 4), id_parent_2 = rep(paste0("A_", 6:10), 4),
      sex = rep(c("M", "F"), 10), line_name = "A"), seed = 314L, store_crossovers = sx))
    h  <- DBI::dbGetQuery(pop$db_conn, "SELECT * FROM ind_haplotype ORDER BY id_ind,parent_origin,locus_id")
    xc <- DBI::dbGetQuery(pop$db_conn, "SELECT * FROM ind_crossover")
    close_pop(pop); list(h = h, xc = xc)
  }
  off <- run(FALSE); on <- run(TRUE)

  expect_equal(off$h, on$h)                # no side effect on the allele stream
  expect_equal(nrow(off$xc), 0L)           # off => empty
  expect_gt(nrow(on$xc), 0L)
  expect_setequal(names(on$xc),
    c("id_crossover", "id_ind", "parent_origin", "chr", "chr_name", "pos_cM"))
  expect_false(any(duplicated(on$xc$id_crossover)))          # unique BIGINT ids
  expect_true(all(on$xc$parent_origin %in% c(1L, 2L)))
  # mean crossovers per gamete-chromosome ~ Poisson(len_M); len = 200 cM -> lambda 2.
  n_gam_chr <- 20 * 2 * 3                   # offspring x parent_origin x chr
  expect_equal(nrow(on$xc) / n_gam_chr, 2, tolerance = 0.25)
  expect_true(all(on$xc$pos_cM >= 0 & on$xc$pos_cM <= 200))
})

test_that("ind_crossover EVENTS are batch-invariant (surrogate id excluded)", {
  # Crossover events (id_ind, parent_origin, chr, chr_name, pos_cM) must be
  # identical across batch sizes for a fixed seed; only id_crossover — a surrogate
  # assigned in write order — may differ, so it is dropped from the comparison.
  run <- function(bs) {
    set.seed(51)
    pop <- open_pop(pop_name = "xb", db_name = ":memory:") |>
      define_genome(n_loci = 150, n_chr = 3, chr_len_Mb = 150) |>
      define_founder_haplotypes(n_haplotypes = 30, line_name = "A", method = "uniform")
    pop <- pop |> get_table("founder_haplotypes") |>
      add_founders(n_males = 5, n_females = 5, line_name = "A")
    pop <- suppressMessages(add_offspring(pop, tibble::tibble(
      id_parent_1 = rep(paste0("A_", 1:5), 5), id_parent_2 = rep(paste0("A_", 6:10), 5),
      sex = rep(c("M", "F"), length.out = 25), line_name = "A"),
      seed = 707L, store_crossovers = TRUE, batch_size = bs))
    xc <- DBI::dbGetQuery(pop$db_conn,
      "SELECT id_ind, parent_origin, chr, chr_name, pos_cM FROM ind_crossover
       ORDER BY id_ind, parent_origin, chr, pos_cM")
    close_pop(pop); xc
  }
  small <- run(1L); big <- run(1000L)
  expect_gt(nrow(small), 0L)
  expect_equal(small, big)
})

test_that("recorded crossovers reproduce the homolog switches in the gamete", {
  # Controlled fixture: the two homologs differ at every locus (row1 all 0, row2
  # all 1), so allele == homolog - 1 and every homolog switch is visible. Assert
  # the returned crossover positions delimit exactly the observed allele switches.
  L    <- 400L
  pos  <- seq(0.5, 299.5, length.out = L)      # single 300 cM chromosome (lambda 3)
  ci   <- list("1" = list(locus_idx = seq_len(L), pos_cM = pos, chr_len = max(pos)))
  hapm <- rbind(rep(0L, L), rep(1L, L))

  for (sid in 1:25) {
    g  <- make_gamete(hapm, ci, base_seed = 2718L, sid = sid, store_crossovers = TRUE)
    xp <- g$xover_pos_cM
    obs_switch <- which(diff(g$allele) != 0L)   # adjacent loci with different allele
    exp_switch <- which(vapply(seq_len(L - 1L),
      function(i) sum(xp > pos[i] & xp < pos[i + 1L]) %% 2L == 1L, logical(1)))
    expect_identical(obs_switch, exp_switch)
    expect_true(all(diff(xp) >= 0))             # positions returned sorted
  }
})
