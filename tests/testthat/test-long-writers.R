# Checkpoints for the long-format writers. These verify ind_haplotype structure,
# the line_origin / strand columns, and line_origin inheritance across
# generations. (Whole-simulation value parity is covered by test-parity.R.)

test_that("add_founders writes long ind_haplotype with line_origin/strand", {
  pop <- make_test_pop(n_loci = 40, n_chr = 2, n_males = 3, n_females = 3,
                       n_haplotypes = 20)
  on.exit(close_pop(pop), add = TRUE)

  ih <- DBI::dbGetQuery(pop$db_conn, "SELECT * FROM ind_haplotype")
  expect_true(all(ih$allele %in% c(0L, 1L)))
  # strand always 1 for diploids; line_origin = founder line ('A').
  expect_true(all(ih$strand == 1L))
  expect_true(all(ih$line_origin == "A"))
  # 2 haplotype rows x 40 loci x 6 founders = 480 rows.
  expect_equal(nrow(ih), 2L * 40L * 6L)

  # ind_genotype is on-demand only — add_founders must NOT populate it.
  expect_equal(DBI::dbGetQuery(pop$db_conn,
    "SELECT COUNT(*) AS n FROM ind_genotype")$n, 0)
})

test_that("add_founders line_origin tracks the line for multiple lines", {
  pop <- open_pop(pop_name = "ml", db_name = ":memory:") |>
    define_genome(n_loci = 20, n_chr = 2, chr_len_Mb = 100) |>
    define_founder_haplotypes(n_haplotypes = 15, line_name = "A") |>
    define_founder_haplotypes(n_haplotypes = 15, line_name = "B")
  on.exit(close_pop(pop), add = TRUE)

  pop <- pop |> get_table("founder_haplotypes") |>
    dplyr::filter(line_name == "A") |>
    add_founders(n_males = 2, n_females = 2, line_name = "A")
  pop <- pop |> get_table("founder_haplotypes") |>
    dplyr::filter(line_name == "B") |>
    add_founders(n_males = 2, n_females = 2, line_name = "B")

  lo <- DBI::dbGetQuery(pop$db_conn,
    "SELECT DISTINCT id_ind, line_origin FROM ind_haplotype")
  expect_true(all(lo$line_origin == sub("_[0-9]+$", "", lo$id_ind)))
})

test_that("add_offspring writes long ind_haplotype and inherits line_origin (F1/F2)", {
  set.seed(123)
  pop <- open_pop(pop_name = "cross", db_name = ":memory:") |>
    define_genome(n_loci = 60, n_chr = 3, chr_len_Mb = 100) |>
    define_founder_haplotypes(n_haplotypes = 20, line_name = "A") |>
    define_founder_haplotypes(n_haplotypes = 20, line_name = "B")
  on.exit(close_pop(pop), add = TRUE)

  pop <- pop |> get_table("founder_haplotypes") |>
    dplyr::filter(line_name == "A") |>
    add_founders(n_males = 3, n_females = 3, line_name = "A")
  pop <- pop |> get_table("founder_haplotypes") |>
    dplyr::filter(line_name == "B") |>
    add_founders(n_males = 3, n_females = 3, line_name = "B")

  # F1 = A sire x B dam.
  f1 <- tibble::tibble(id_parent_1 = c("A_1", "A_2"),
                       id_parent_2 = c("B_4", "B_5"),
                       sex = c("M", "F"), line_name = "F1")
  pop <- add_offspring(pop, f1)

  # F1: parent_origin 1 (from A sire, a pure-A founder) is all line A;
  # parent_origin 2 (from B dam) is all line B.
  f1_lo <- DBI::dbGetQuery(pop$db_conn,
    "SELECT DISTINCT parent_origin, line_origin FROM ind_haplotype WHERE id_ind LIKE 'F1_%' ORDER BY parent_origin")
  expect_equal(f1_lo$line_origin[f1_lo$parent_origin == 1], "A")
  expect_equal(f1_lo$line_origin[f1_lo$parent_origin == 2], "B")

  # F2 = F1 x F1: gametes recombine across each F1's A and B homologs, so F2
  # haplotypes are a mosaic containing BOTH founding lines.
  f2 <- tibble::tibble(id_parent_1 = "F1_1", id_parent_2 = "F1_2",
                       sex = "M", line_name = "F2")
  pop <- add_offspring(pop, f2)

  f2_lo <- DBI::dbGetQuery(pop$db_conn,
    "SELECT DISTINCT line_origin FROM ind_haplotype WHERE id_ind = 'F2_1'")$line_origin
  expect_setequal(f2_lo, c("A", "B"))
  # No allele's origin is untracked.
  expect_false(any(is.na(f2_lo)))
})
