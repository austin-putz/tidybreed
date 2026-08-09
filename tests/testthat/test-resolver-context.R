# v5 review #4: inheritance resolves by the OFFSPRING's sex+line; recombination
# resolves by the PRODUCING PARENT's sex+line. The two tables resolve
# independently, so a sex-specific copy rule and a line-specific recombination
# rule can never interact.

test_that("inheritance and recombination are keyed by different, independent lines", {
  pop <- make_pop_base(n_loci = 6, n_chr = 1, chr_len_Mb = 50)
  on.exit(close_pop(pop), add = TRUE)
  conn <- pop$db_conn

  # Line-specific INHERITANCE: line A offspring get 0,1; line B offspring get 1,0.
  suppressMessages({
    pop <- define_chromosome(pop, "1", line_name = "A", from_parent_1 = 0, from_parent_2 = 1)
    pop <- define_chromosome(pop, "1", line_name = "B", from_parent_1 = 1, from_parent_2 = 0)
    # Line-specific RECOMBINATION: parent line A recombines, parent line B does not.
    pop <- define_chromosome(pop, "1", parent_sex = NULL, line_name = "A", recombines = TRUE)
    pop <- define_chromosome(pop, "1", parent_sex = NULL, line_name = "B", recombines = FALSE)
  })

  # Inheritance follows the OFFSPRING line.
  inh_A <- resolve_chr_inheritance(conn, "M", "A")
  inh_B <- resolve_chr_inheritance(conn, "M", "B")
  expect_equal(c(inh_A$from_parent_1, inh_A$from_parent_2), c(0L, 1L))
  expect_equal(c(inh_B$from_parent_1, inh_B$from_parent_2), c(1L, 0L))

  # Recombination follows the PRODUCING-PARENT line, independently.
  expect_true(resolve_chr_recombination(conn, "M", "A")$recombines)
  expect_false(resolve_chr_recombination(conn, "M", "B")$recombines)

  # Cross-independence: the line-B inheritance rule (1,0) does not change line-A
  # recombination, and the line-B recombination rule (FALSE) does not change
  # line-A inheritance.
  expect_true(resolve_chr_recombination(conn, "M", "A")$recombines)
  expect_equal(c(resolve_chr_inheritance(conn, "F", "A")$from_parent_1,
                 resolve_chr_inheritance(conn, "F", "A")$from_parent_2), c(0L, 1L))
})

test_that("get_chr_rules_map() resolves at the requested line, defaulting when absent", {
  pop <- make_pop_base(n_loci = 6, n_chr = 1, chr_len_Mb = 50)
  on.exit(close_pop(pop), add = TRUE)
  conn <- pop$db_conn
  suppressMessages(
    pop <- define_chromosome(pop, "1", line_name = "A", from_parent_1 = 0, from_parent_2 = 1))

  # Line A sees its own rule; an unknown line C falls back to the (NULL,NULL) seed.
  mapA <- get_chr_rules_map(conn, "A")
  mapC <- get_chr_rules_map(conn, "C")
  expect_equal(unname(mapA[["1"]]$from_parent_1[["M"]]), 0L)
  expect_equal(unname(mapA[["1"]]$from_parent_2[["M"]]), 1L)
  expect_equal(unname(mapC[["1"]]$from_parent_1[["M"]]), 1L)  # seed default 1,1
  expect_equal(unname(mapC[["1"]]$from_parent_2[["M"]]), 1L)
})
