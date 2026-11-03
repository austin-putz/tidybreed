# Recombination delivers the genetic map: realised recombination fractions
# match the Haldane map function r = (1 - exp(-2d/100)) / 2 for the cM
# distances in genome_map, crossover counts are Poisson(L / 100), and the map
# a gamete uses is its producing parent's. Every other recombination test
# checks reproducibility (seed, batch size, R <-> C++); these check that the
# model is the right one. All runs are seeded, so the outcome is fixed; the
# |z| < 4 bound is wide enough to hold for any seed of a correct model and
# narrow enough to catch a wrong distance, rate or map.

haldane <- function(d_cM) 0.5 * (1 - exp(-2 * d_cM / 100))

# z-scores of the observed recombination fraction of every locus pair in `W`
# (gametes x loci, TRUE = first homolog / line A) against `r_exp(i, j)`.
pair_z <- function(W, r_exp) {
  G <- nrow(W)
  z <- c()
  for (i in 1:(ncol(W) - 1)) for (j in (i + 1):ncol(W)) {
    e <- r_exp(i, j)
    z <- c(z, (mean(W[, i] != W[, j]) - e) / sqrt(e * (1 - e) / G))
  }
  z
}

# ---- the kernel alone ------------------------------------------------------

test_that("the gamete kernel delivers Haldane recombination on an irregular map", {
  skip_if_not(exists("make_gametes_batch_cpp", mode = "function"),
              "compiled kernel unavailable")

  # Irregular spacing; chromosome 1 starts at 7 cM, chromosome 3 is short
  pos <- c(7, 7.5, 9, 15, 30, 31, 60, 95, 140, 150,   # chr 1, 150 cM
           0.1, 1, 2, 5, 10, 20, 40, 50,              # chr 2,  50 cM
           3, 12)                                     # chr 3,  12 cM
  chr <- rep(1:3, c(10, 8, 2))
  L   <- length(pos)
  chr_start <- c(1L, 11L, 19L)
  chr_end   <- c(10L, 18L, 20L)
  chr_len   <- c(150, 50, 12)

  # One parent whose homologs differ everywhere: allele == homolog - 1
  pa <- rbind(rep(0L, L), rep(1L, L))
  G  <- 50000L
  out <- make_gametes_batch(pa, pa + 1L, rep(1L, G), seq_len(G), rep(1L, G),
                            chr_start, chr_end, pos, chr_len, 20261103L, TRUE)
  W <- matrix(out$allele == 0L, nrow = G, byrow = TRUE)

  # Same chromosome: Haldane at the cM distance. Different chromosomes: 0.5.
  z <- pair_z(W, function(i, j)
    if (chr[i] == chr[j]) haldane(pos[j] - pos[i]) else 0.5)
  expect_lt(max(abs(z)), 4)

  # No transmission bias at any locus
  expect_lt(max(abs(colMeans(W) - 0.5)) / sqrt(0.25 / G), 4)

  # No interference: double recombinants at the product of the two rates
  r12 <- W[, 1] != W[, 5]
  r23 <- W[, 5] != W[, 8]
  expect_equal(mean(r12 & r23) / (mean(r12) * mean(r23)), 1, tolerance = 0.05)

  # Crossovers per chromosome: Poisson(L / 100), uniform on [0, L]
  for (c_ in 1:3) {
    n   <- tabulate(out$xover_gamete_o[out$xover_chr == c_], nbins = G)
    lam <- chr_len[c_] / 100
    expect_lt(abs(mean(n) - lam) / sqrt(lam / G), 4)
    expect_equal(stats::var(n), lam, tolerance = 0.05)
    xp <- out$xover_pos_cM[out$xover_chr == c_]
    expect_gt(suppressWarnings(
      stats::ks.test(xp, "punif", 0, chr_len[c_])$p.value), 1e-4)
  }
})

# ---- end to end through add_offspring() ------------------------------------

test_that("add_offspring() recombines at the producing parent's map rates", {
  # Line A sires x line B dams give F1 that carry one A and one B copy at every
  # locus, so in an F1 gamete every line_origin switch is a recombination,
  # whatever the alleles. Chromosome 1 has a male map twice the default;
  # chromosome 2 runs at 0.5 cM/Mb; chromosome 3 is an X.
  set.seed(42)
  pop <- open_pop(pop_name = "rr", db_name = ":memory:") |>
    define_genome(n_loci = 45, n_chr = 3, chr_len_Mb = c(150, 100, 80),
                  cM_per_Mb = c(1, 0.5, 1)) |>
    define_chromosome("3", offspring_sex = "M",
                      from_parent_1 = 0, from_parent_2 = 1) |>
    define_founder_haplotypes(n_haplotypes = 10, line_name = "A") |>
    define_founder_haplotypes(n_haplotypes = 10, line_name = "B")
  on.exit(close_pop(pop), add = TRUE)

  base <- resolve_genome_map(pop$db_conn)
  m1   <- base[base$chr == 1, ]
  male <- data.frame(
    id_genome_map = seq.int(next_int_id(pop$db_conn, "genome_map", "id_genome_map"),
                            length.out = nrow(m1)),
    locus_id = m1$locus_id, locus_name = m1$locus_name, sex = "M",
    line_name = NA_character_, map_name = "default", pos_cM = m1$pos_cM * 2)
  duckdb::duckdb_register(pop$db_conn, "__male_map", male)
  DBI::dbExecute(pop$db_conn, "INSERT INTO genome_map SELECT * FROM __male_map")
  duckdb::duckdb_unregister(pop$db_conn, "__male_map")
  validate_genome_map(pop$db_conn)

  pop <- pop |> get_table("founder_haplotypes") |>
    dplyr::filter(line_name == "A") |>
    add_founders(n_males = 5, n_females = 0, line_name = "A")
  pop <- pop |> get_table("founder_haplotypes") |>
    dplyr::filter(line_name == "B") |>
    add_founders(n_males = 0, n_females = 5, line_name = "B")
  pop <- suppressMessages(add_offspring(pop, data.frame(
    id_parent_1 = paste0("A_", rep(1:5, 4)), id_parent_2 = paste0("B_", rep(1:5, 4)),
    sex = rep(c("M", "F"), 10), line_name = "F1"), seed = 1))

  N <- 2000L
  pop <- suppressMessages(add_offspring(pop, data.frame(
    id_parent_1 = paste0("F1_", rep(c(1, 3, 5, 7, 9), length.out = N)),
    id_parent_2 = paste0("F1_", rep(c(2, 4, 6, 8, 10), length.out = N)),
    sex = rep(c("M", "F"), N / 2), line_name = "F2"),
    seed = 2, store_crossovers = TRUE))

  hap <- DBI::dbGetQuery(pop$db_conn, "
    SELECT h.id_ind, m.sex AS off_sex, h.parent_origin, h.locus_id, h.line_origin
    FROM ind_haplotype h JOIN ind_meta m USING (id_ind)
    WHERE m.line_name = 'F2' ORDER BY h.id_ind, h.parent_origin, h.locus_id")

  # Gametes of parent_origin `po` at chromosome `chr_id` against the map of
  # parent sex `map_sex`
  check_map <- function(po, chr_id, map_sex, off_sex = c("M", "F")) {
    mp <- resolve_genome_map(pop$db_conn, sex = map_sex)
    mp <- mp[mp$chr == chr_id, ]
    h  <- hap[hap$parent_origin == po & hap$locus_id %in% mp$locus_id &
              hap$off_sex %in% off_sex, ]
    W  <- matrix(h$line_origin == "A", ncol = nrow(mp), byrow = TRUE)
    pair_z(W, function(i, j) haldane(mp$pos_cM[j] - mp$pos_cM[i]))
  }
  expect_lt(max(abs(check_map(1L, 1L, "M"))), 4)   # sire: the 2x male map
  expect_lt(max(abs(check_map(2L, 1L, "F"))), 4)   # dam: the default map
  expect_lt(max(abs(check_map(1L, 2L, "M"))), 4)   # 0.5 cM/Mb chromosome
  expect_lt(max(abs(check_map(2L, 3L, "F"))), 4)   # X in a female recombines

  # A sire read against the default map is rejected: the test can tell maps apart
  expect_gt(max(abs(check_map(1L, 1L, "F"))), 10)

  # The hemizygous F1 male passes his X (from his B dam) unchanged to
  # daughters, and nothing to sons
  x_ids <- base$locus_id[base$chr == 3]
  pat_x <- hap[hap$parent_origin == 1L & hap$locus_id %in% x_ids, ]
  expect_true(all(pat_x$off_sex == "F"))
  expect_true(all(pat_x$line_origin == "B"))

  # Crossover counts: Poisson(map length / 100) for the producing parent's map
  xo <- DBI::dbGetQuery(pop$db_conn, "
    SELECT c.parent_origin, c.chr, COUNT(*) AS n
    FROM ind_crossover c JOIN ind_meta m USING (id_ind)
    WHERE m.line_name = 'F2' GROUP BY ALL")
  len <- function(sex, chr_id) {
    mp <- resolve_genome_map(pop$db_conn, sex = sex)
    max(mp$pos_cM[mp$chr == chr_id])
  }
  for (k in list(list(1L, 1L, "M"), list(2L, 1L, "F"), list(2L, 3L, "F"))) {
    lam <- len(k[[3]], k[[2]]) / 100
    obs <- xo$n[xo$parent_origin == k[[1]] & xo$chr == k[[2]]] / N
    expect_lt(abs(obs - lam) / sqrt(lam / N), 4)
  }
  expect_false(any(xo$parent_origin == 1L & xo$chr == 3L))
})
