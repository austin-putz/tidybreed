# Line-specific chr_inheritance / chr_recombination rules in add_offspring().
#
# A chromosome that is a plain autosome under the line-agnostic rules used to
# take the gamete kernel, which never reads a rule, so a line rule on it was
# ignored. A chromosome is now a plain autosome only when it is plain for every
# line with rules (plans/line_specific_chr_rules_add_offspring.md, D1/D2).

# Population with chromosome rules applied before any founder exists. `rules`
# is a function(pop) -> pop; `founders` a list of list(line, n_males, n_females).
lcr_pop <- function(rules = identity, n_chr = 3, n_loci = 60,
                    founders = list(list("A", 4, 4), list("B", 4, 4))) {
  set.seed(11)
  pop <- suppressMessages(
    open_pop(pop_name = "lcr", db_name = ":memory:") |>
      define_genome(n_loci = n_loci, n_chr = n_chr, chr_len_Mb = 100))
  pop <- suppressMessages(rules(pop))
  pop <- suppressWarnings(suppressMessages(
    define_founder_haplotypes(pop, n_haplotypes = 30)))
  for (f in founders) {
    pop <- suppressMessages(pop |> get_table("founder_haplotypes") |>
      add_founders(n_males = f[[2]], n_females = f[[3]], line_name = f[[1]]))
  }
  pop
}

lcr_hap <- function(pop, ids = NULL) {
  where <- if (is.null(ids)) "" else paste0(
    "WHERE h.id_ind IN (", paste0("'", ids, "'", collapse = ", "), ") ")
  DBI::dbGetQuery(pop$db_conn, paste0(
    "SELECT h.id_ind, h.parent_origin, g.chr_name, h.locus_id, h.allele, ",
    "h.line_origin FROM ind_haplotype h JOIN genome_meta g USING (locus_id) ",
    where, "ORDER BY h.id_ind, h.parent_origin, h.locus_id"))
}

lcr_xover <- function(pop, ids = NULL) {
  where <- if (is.null(ids)) "" else paste0(
    "WHERE id_ind IN (", paste0("'", ids, "'", collapse = ", "), ") ")
  DBI::dbGetQuery(pop$db_conn, paste0(
    "SELECT id_ind, parent_origin, chr, chr_name, pos_cM FROM ind_crossover ",
    where, "ORDER BY id_ind, parent_origin, chr, pos_cM"))
}

# Parent origins an individual carries on a chromosome.
lcr_origins <- function(h, id, chr) {
  sort(unique(h$parent_origin[h$id_ind == id & h$chr_name == chr]))
}

lcr_mating <- function(p1, p2, sex, line) {
  data.frame(id_parent_1 = p1, id_parent_2 = p2, sex = sex, line_name = line,
             stringsAsFactors = FALSE)
}

# chr 1: F1 males inherit only the dam's copy; line-A sires do not recombine.
lcr_mixed_rules <- function(pop) {
  pop |>
    define_chromosome("1", offspring_sex = "M", line_name = "F1",
                      from_parent_1 = 0, from_parent_2 = 1) |>
    define_chromosome("1", parent_sex = "M", line_name = "A",
                      recombines = FALSE)
}

lcr_mixed_matings <- function() {
  rbind(
    lcr_mating(c("A_1", "A_2", "A_3", "A_4"), c("B_5", "B_6", "B_7", "B_8"),
               c("M", "F", "M", "F"), "F1"),
    lcr_mating(c("B_1", "B_2", "B_3", "B_4"), c("A_5", "A_6", "A_7", "A_8"),
               c("M", "F", "M", "F"), "C"),
    lcr_mating(c("A_1", "A_2"), c("A_5", "A_6"), c("M", "F"), "A"))
}


test_that("a line inheritance rule on a default autosome is honored", {
  pop <- lcr_pop(function(pop) define_chromosome(
    pop, "1", offspring_sex = "M", line_name = "F1",
    from_parent_1 = 0, from_parent_2 = 1))
  on.exit(close_pop(pop))
  mat <- rbind(
    lcr_mating(c("A_1", "A_2", "A_3"), c("A_5", "A_6", "A_7"),
               c("M", "M", "F"), "F1"),
    lcr_mating(c("A_1", "A_2"), c("A_5", "A_6"), c("M", "F"), "A"))
  pop <- suppressMessages(add_offspring(pop, mat, seed = 3L))

  h <- lcr_hap(pop)
  for (id in c("F1_1", "F1_2")) {
    expect_identical(lcr_origins(h, id, "1"), 2L)
    expect_identical(lcr_origins(h, id, "2"), c(1L, 2L))
    expect_identical(lcr_origins(h, id, "3"), c(1L, 2L))
  }
  for (id in c("F1_3", "A_9", "A_10")) {
    for (chr in c("1", "2", "3")) {
      expect_identical(lcr_origins(h, id, chr), c(1L, 2L))
    }
  }
})


test_that("a parent made hemizygous by its founder line's rule can breed", {
  pop <- lcr_pop(
    function(pop) define_chromosome(pop, "1", offspring_sex = "M",
                                    line_name = "H", from_parent_1 = 0,
                                    from_parent_2 = 1),
    founders = list(list("H", 2, 0), list("A", 0, 4)))
  on.exit(close_pop(pop))

  h0 <- lcr_hap(pop, "H_1")
  expect_identical(lcr_origins(h0, "H_1", "1"), 2L)   # founder: one copy

  mat <- lcr_mating(c("H_1", "H_1"), c("A_1", "A_2"), c("M", "F"), "A")
  pop <- suppressMessages(add_offspring(pop, mat, seed = 5L))

  h <- lcr_hap(pop)
  sire_chr1 <- h0$allele[h0$chr_name == "1"]
  for (id in c("A_5", "A_6")) {
    expect_identical(lcr_origins(h, id, "1"), c(1L, 2L))
    # The sire's single copy is passed through unchanged.
    expect_identical(
      h$allele[h$id_ind == id & h$chr_name == "1" & h$parent_origin == 1L],
      sire_chr1)
  }
})


test_that("a line recombination rule is resolved at the producing parent's line", {
  pop <- lcr_pop()
  pop <- suppressMessages(define_chromosome(
    pop, "1", parent_sex = "M", line_name = "A", recombines = FALSE))
  on.exit(close_pop(pop))

  # Offspring are line C, so resolving at the offspring's line would miss the rule.
  sires_a <- rep(c("A_1", "A_2", "A_3", "A_4"), 5)
  sires_b <- rep(c("B_1", "B_2", "B_3", "B_4"), 5)
  mat <- rbind(
    lcr_mating(sires_a, rep(c("B_5", "B_6", "B_7", "B_8"), 5), "F", "C"),
    lcr_mating(sires_b, rep(c("A_5", "A_6", "A_7", "A_8"), 5), "F", "C"))
  pop <- suppressMessages(add_offspring(pop, mat, seed = 7L,
                                        store_crossovers = TRUE))

  ids_a <- paste0("C_", 1:20)
  ids_b <- paste0("C_", 21:40)
  xo <- lcr_xover(pop)
  chr1 <- xo[xo$chr_name == "1", ]

  # Line-A sires: no crossover on chr 1, and the paternal copy is one of the
  # sire's own two homologs, unchanged.
  expect_false(any(chr1$id_ind %in% ids_a & chr1$parent_origin == 1L))
  h <- lcr_hap(pop)
  for (k in seq_along(ids_a)) {
    got <- h$allele[h$id_ind == ids_a[k] & h$chr_name == "1" &
                      h$parent_origin == 1L]
    s   <- h[h$id_ind == sires_a[k] & h$chr_name == "1", ]
    expect_true(identical(got, s$allele[s$parent_origin == 1L]) ||
                  identical(got, s$allele[s$parent_origin == 2L]))
  }
  # Line-B sires and both lines' dams still recombine on chr 1.
  expect_true(any(chr1$id_ind %in% ids_b & chr1$parent_origin == 1L))
  expect_true(any(chr1$id_ind %in% ids_a & chr1$parent_origin == 2L))
  expect_true(any(chr1$id_ind %in% ids_b & chr1$parent_origin == 2L))
})


test_that("line rules: the same seed gives bit-identical offspring", {
  run <- function() {
    pop <- lcr_pop(lcr_mixed_rules)
    on.exit(close_pop(pop))
    pop <- suppressMessages(add_offspring(pop, lcr_mixed_matings(), seed = 9L,
                                          store_crossovers = TRUE))
    list(h = lcr_hap(pop), x = lcr_xover(pop))
  }
  expect_identical(run(), run())
})


test_that("classification depends on stored rules, not on the matings in the call", {
  # A rule for line B2 only. Line-A offspring of rows 1..n must not depend on
  # whether B2 matings follow them in the same call (streams are keyed on row).
  rules <- function(pop) define_chromosome(
    pop, "1", offspring_sex = "M", line_name = "B2",
    from_parent_1 = 0, from_parent_2 = 1)
  mat_a  <- lcr_mating(c("A_1", "A_2", "A_3"), c("A_5", "A_6", "A_7"),
                       c("M", "F", "M"), "A")
  mat_b2 <- lcr_mating(c("B_1", "B_2"), c("B_5", "B_6"), c("M", "F"), "B2")
  run <- function(mat) {
    pop <- lcr_pop(rules)
    on.exit(close_pop(pop))
    pop <- suppressMessages(add_offspring(pop, mat, seed = 13L,
                                          store_crossovers = TRUE))
    ids <- c("A_9", "A_10", "A_11")
    list(h = lcr_hap(pop, ids), x = lcr_xover(pop, ids))
  }
  expect_identical(run(mat_a), run(rbind(mat_a, mat_b2)))
})


test_that(".plain_autosome_chrs() matches the line-agnostic classification without line rules", {
  pop <- lcr_pop(function(pop) {
    pop |>
      define_chromosome("3", offspring_sex = "M",
                        from_parent_1 = 0, from_parent_2 = 1) |>
      define_chromosome("2", recombines = FALSE)
  })
  on.exit(close_pop(pop))
  rules <- get_chr_rules_map(pop$db_conn)
  old   <- names(rules)[vapply(rules, is_plain_autosome, logical(1))]
  expect_identical(.plain_autosome_chrs(pop$db_conn), old)
  expect_identical(old, "1")
})


test_that("a genome with no fast-path chromosome still breeds", {
  pop <- lcr_pop(function(pop) define_chromosome(
    pop, "1", offspring_sex = "M", line_name = "F1",
    from_parent_1 = 0, from_parent_2 = 1), n_chr = 1, n_loci = 30)
  on.exit(close_pop(pop))
  expect_identical(.plain_autosome_chrs(pop$db_conn), character(0))

  mat <- rbind(lcr_mating(c("A_1", "A_2"), c("A_5", "A_6"), c("M", "F"), "F1"),
               lcr_mating("A_1", "A_5", "M", "A"))
  pop <- suppressMessages(add_offspring(pop, mat, seed = 17L))
  h <- lcr_hap(pop)
  expect_identical(lcr_origins(h, "F1_1", "1"), 2L)
  expect_identical(lcr_origins(h, "F1_2", "1"), c(1L, 2L))
  expect_identical(lcr_origins(h, "A_9",  "1"), c(1L, 2L))
  expect_identical(sum(h$id_ind == "A_9"), 60L)
})


test_that("a line rule of 0,0 leaves the chromosome out for those offspring", {
  pop <- lcr_pop(function(pop) define_chromosome(
    pop, "2", offspring_sex = "M", line_name = "F1",
    from_parent_1 = 0, from_parent_2 = 0))
  on.exit(close_pop(pop))
  mat <- rbind(lcr_mating(c("A_1", "A_2"), c("A_5", "A_6"), c("M", "F"), "F1"),
               lcr_mating("A_1", "A_5", "M", "A"))
  pop <- suppressMessages(add_offspring(pop, mat, seed = 19L))
  h <- lcr_hap(pop)
  expect_identical(lcr_origins(h, "F1_1", "2"), integer(0))
  expect_identical(lcr_origins(h, "F1_1", "1"), c(1L, 2L))
  expect_identical(lcr_origins(h, "F1_2", "2"), c(1L, 2L))
  expect_identical(lcr_origins(h, "A_9",  "2"), c(1L, 2L))
})


test_that("an unsupported line rule fails and rolls back every offspring table", {
  pop <- lcr_pop(function(pop) define_chromosome(
    pop, "1", offspring_sex = "M", line_name = "F1",
    from_parent_1 = 2, from_parent_2 = 0))
  on.exit(close_pop(pop))
  snap <- function() lapply(
    c(ind_meta = "SELECT * FROM ind_meta ORDER BY id_ind",
      ind_haplotype = "SELECT * FROM ind_haplotype ORDER BY id_ind, parent_origin, locus_id",
      ind_crossover = "SELECT * FROM ind_crossover ORDER BY id_crossover"),
    function(q) DBI::dbGetQuery(pop$db_conn, q))
  before <- snap()

  # The valid line-A row is written in batch 1 before batch 2 fails.
  mat <- rbind(lcr_mating("A_1", "A_5", "F", "A"),
               lcr_mating("A_2", "A_6", "M", "F1"))
  expect_error(
    suppressMessages(add_offspring(pop, mat, seed = 23L, batch_size = 1L,
                                   store_crossovers = TRUE)),
    "uniparental")
  expect_identical(snap(), before)
})


test_that("a fully shadowed line rule keeps the chromosome on the fast path", {
  shadowed <- function(pop) {
    pop |>
      define_chromosome("1", line_name = "L",
                        from_parent_1 = 0, from_parent_2 = 1) |>
      define_chromosome("1", offspring_sex = "M", line_name = "L",
                        from_parent_1 = 1, from_parent_2 = 1) |>
      define_chromosome("1", offspring_sex = "F", line_name = "L",
                        from_parent_1 = 1, from_parent_2 = 1)
  }
  run <- function(rules) {
    pop <- lcr_pop(rules)
    on.exit(close_pop(pop))
    plain <- .plain_autosome_chrs(pop$db_conn)
    pop <- suppressMessages(add_offspring(pop, lcr_mixed_matings()[9:10, ],
                                          seed = 29L, store_crossovers = TRUE))
    list(plain = plain, h = lcr_hap(pop), x = lcr_xover(pop))
  }
  with_rows <- run(shadowed)
  expect_identical(with_rows$plain, c("1", "2", "3"))
  expect_identical(with_rows, run(identity))
})


test_that("line rules: R and C++ kernels and any batch size agree", {
  run <- function(kernel, batch_size) {
    withr::local_options(tidybreed.kernel = kernel)
    pop <- lcr_pop(lcr_mixed_rules)
    on.exit(close_pop(pop))
    mat <- lcr_mixed_matings()
    pop <- suppressMessages(add_offspring(pop, mat, seed = 31L,
                                          store_crossovers = TRUE,
                                          batch_size = batch_size))
    list(h = lcr_hap(pop), x = lcr_xover(pop))
  }
  ref <- run("auto", 1L)
  expect_identical(run("auto", 10L), ref)
  expect_identical(run("r", 1L), ref)
  expect_identical(run("r", 10L), ref)
})
