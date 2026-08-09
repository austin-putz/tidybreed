# Living record of the cross-species coverage claims in plans/update_chr_meta.md:
# build each configuration and assert the RESOLVED per-parent copy counts and
# per-parent-sex recombination.

.inh <- function(conn, sex, chr) {
  r <- resolve_chr_inheritance(conn, sex)
  unname(as.integer(r[r$chr_name == chr, c("from_parent_1", "from_parent_2")]))
}
.rec <- function(conn, sex, chr) {
  r <- resolve_chr_recombination(conn, sex)
  r$recombines[r$chr_name == chr]
}

test_that("ZW (avian/silkworm): heterogametic females carry the single-copy W", {
  suppressMessages({
    pop <- open_pop(pop_name = "cov_zw", db_name = ":memory:") |>
      define_genome(n_loci = 9, n_chr = 3, chr_names = c("1", "Z", "W"),
                    chr_len_Mb = 50) |>
      define_chromosome("Z", offspring_sex = "F", from_parent_1 = 1, from_parent_2 = 0) |>
      define_chromosome("W", offspring_sex = "M", from_parent_1 = 0, from_parent_2 = 0) |>
      define_chromosome("W", offspring_sex = "F", from_parent_1 = 0, from_parent_2 = 1) |>
      define_chromosome("W", recombines = FALSE)
  })
  conn <- pop$db_conn
  on.exit(close_pop(pop), add = TRUE)
  expect_equal(.inh(conn, "M", "Z"), c(1L, 1L))   # ZZ son
  expect_equal(.inh(conn, "F", "Z"), c(1L, 0L))   # Z from sire
  expect_equal(.inh(conn, "M", "W"), c(0L, 0L))   # no W
  expect_equal(.inh(conn, "F", "W"), c(0L, 1L))   # W from dam
  expect_false(.rec(conn, "M", "W"))
})

test_that("X0: hemizygous X with no partner chromosome", {
  suppressMessages({
    pop <- open_pop(pop_name = "cov_x0", db_name = ":memory:") |>
      define_genome(n_loci = 9, n_chr = 3, chr_names = c("1", "2", "X"),
                    chr_len_Mb = 50) |>
      define_chromosome("X", offspring_sex = "M", from_parent_1 = 0, from_parent_2 = 1)
  })
  conn <- pop$db_conn
  on.exit(close_pop(pop), add = TRUE)
  expect_equal(.inh(conn, "M", "X"), c(0L, 1L))
  expect_equal(.inh(conn, "F", "X"), c(1L, 1L))
})

test_that("Conifer organelles: paternal chloroplast + maternal mitochondria", {
  suppressMessages({
    pop <- open_pop(pop_name = "cov_conifer", db_name = ":memory:") |>
      define_genome(n_loci = 9, n_chr = 3, chr_names = c("1", "CP", "MT"),
                    chr_len_Mb = 50) |>
      define_chromosome("CP", from_parent_1 = 1, from_parent_2 = 0) |>
      define_chromosome("CP", recombines = FALSE) |>
      define_chromosome("MT", from_parent_1 = 0, from_parent_2 = 1) |>
      define_chromosome("MT", recombines = FALSE)
  })
  conn <- pop$db_conn
  on.exit(close_pop(pop), add = TRUE)
  # Both sexes share the organelle default (offspring_sex = NULL).
  expect_equal(.inh(conn, "M", "CP"), c(1L, 0L))   # paternal chloroplast
  expect_equal(.inh(conn, "F", "CP"), c(1L, 0L))
  expect_equal(.inh(conn, "M", "MT"), c(0L, 1L))   # maternal mitochondria
  expect_equal(.inh(conn, "F", "MT"), c(0L, 1L))
  expect_false(.rec(conn, "M", "CP"))
  expect_false(.rec(conn, "F", "MT"))
})

test_that("Drosophila: genome-wide male achiasmy + non-recombining dot chr4", {
  suppressMessages({
    pop <- open_pop(pop_name = "cov_dmel", db_name = ":memory:") |>
      define_genome(n_loci = 12, n_chr = 4, chr_names = c("2", "3", "4", "X"),
                    chr_len_Mb = 50, recombines_M = FALSE) |>
      # Genome-wide male achiasmy already seeds per-parent-sex rows (M FALSE, F
      # TRUE) for every chromosome. The dot chr4 also does not recombine in
      # FEMALES, so override only the deviating (female) sex.
      define_chromosome("4", parent_sex = "F", recombines = FALSE) |>
      define_chromosome("X", offspring_sex = "M", from_parent_1 = 0, from_parent_2 = 1)
  })
  conn <- pop$db_conn
  on.exit(close_pop(pop), add = TRUE)

  # Male parent recombines nowhere (genome-wide achiasmy seeded at define_genome).
  rec_M <- resolve_chr_recombination(conn, "M")
  expect_true(all(!rec_M$recombines))
  # Female parent recombines everywhere EXCEPT the dot chr4.
  rec_F <- resolve_chr_recombination(conn, "F")
  expect_false(rec_F$recombines[rec_F$chr_name == "4"])
  expect_true(all(rec_F$recombines[rec_F$chr_name != "4"]))
  # X copy counts still correct.
  expect_equal(.inh(conn, "M", "X"), c(0L, 1L))
})
