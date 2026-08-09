# resolve_chr_inheritance()/resolve_chr_recombination() precedence and the
# validate_chr_*() structural checks (mirror of genome_map's resolver/validator).

.ins_inh <- function(conn, chr, sex, line, f1, f2) {
  DBI::dbExecute(conn, sprintf(
    "INSERT INTO chr_inheritance (chr_name, offspring_sex, line_name, from_parent_1, from_parent_2) VALUES ('%s', %s, %s, %d, %d)",
    chr, if (is.na(sex)) "NULL" else sprintf("'%s'", sex),
    if (is.na(line)) "NULL" else sprintf("'%s'", line), f1, f2))
}

test_that("resolve_chr_inheritance() applies the four-tier precedence", {
  pop <- make_pop_base(n_loci = 6, n_chr = 1, chr_len_Mb = 50)
  on.exit(close_pop(pop), add = TRUE)
  conn <- pop$db_conn
  # Seed already inserted (1, NULL, NULL) = (1,1). Add the three higher tiers.
  .ins_inh(conn, "1", "M", NA,   0, 1)   # (sex=M, line=NULL)  prio 2
  .ins_inh(conn, "1", NA,  "LA", 2, 0)   # (sex=NULL, line=LA) prio 1
  .ins_inh(conn, "1", "M", "LA", 1, 0)   # (sex=M, line=LA)    prio 3

  r_m_la <- resolve_chr_inheritance(conn, "M", "LA")   # -> (M,LA)   = 1,0
  r_m_lb <- resolve_chr_inheritance(conn, "M", "LB")   # -> (M,NULL) = 0,1
  r_f_la <- resolve_chr_inheritance(conn, "F", "LA")   # -> (NULL,LA)= 2,0
  r_f_lb <- resolve_chr_inheritance(conn, "F", "LB")   # -> (NULL,NULL)=1,1

  expect_equal(c(r_m_la$from_parent_1, r_m_la$from_parent_2), c(1L, 0L))
  expect_equal(c(r_m_lb$from_parent_1, r_m_lb$from_parent_2), c(0L, 1L))
  expect_equal(c(r_f_la$from_parent_1, r_f_la$from_parent_2), c(2L, 0L))
  expect_equal(c(r_f_lb$from_parent_1, r_f_lb$from_parent_2), c(1L, 1L))
})

test_that("resolvers error when a chromosome has no resolved row", {
  pop <- make_pop_base(n_loci = 6, n_chr = 2, chr_len_Mb = 50)
  on.exit(close_pop(pop), add = TRUE)
  conn <- pop$db_conn
  DBI::dbExecute(conn, "DELETE FROM chr_inheritance WHERE chr_name = '2'")
  expect_error(resolve_chr_inheritance(conn, "M"), "no chr_inheritance row")

  DBI::dbExecute(conn, "DELETE FROM chr_recombination WHERE chr_name = '2'")
  expect_error(resolve_chr_recombination(conn, "M"), "no chr_recombination row")
})

test_that("validate_chr_inheritance() flags a duplicate NULL-normalized logical key", {
  pop <- make_pop_base(n_loci = 6, n_chr = 1, chr_len_Mb = 50)
  on.exit(close_pop(pop), add = TRUE)
  conn <- pop$db_conn
  # A second (1, NULL, NULL) row — a plain UNIQUE could not catch two NULLs.
  .ins_inh(conn, "1", NA, NA, 0, 1)
  expect_error(validate_chr_inheritance(conn), "duplicate logical key")
})

test_that("shadowing check errors on conflict and passes on agreement", {
  pop <- make_pop_base(n_loci = 6, n_chr = 1, chr_len_Mb = 50)
  on.exit(close_pop(pop), add = TRUE)
  conn <- pop$db_conn

  # (M, NULL) = 0,1 and (NULL, LA) = 1,0 differ, no explicit (M, LA) -> conflict.
  .ins_inh(conn, "1", "M", NA,   0, 1)
  .ins_inh(conn, "1", NA,  "LA", 1, 0)
  expect_error(validate_chr_inheritance(conn), "shadowing conflict")

  # Add the explicit (M, LA) row -> resolves the ambiguity, no error.
  .ins_inh(conn, "1", "M", "LA", 0, 1)
  expect_silent(validate_chr_inheritance(conn))
})

test_that("shadowing check does NOT fire when the two rules agree", {
  pop <- make_pop_base(n_loci = 6, n_chr = 1, chr_len_Mb = 50)
  on.exit(close_pop(pop), add = TRUE)
  conn <- pop$db_conn
  # (M, NULL) = 0,1 and (NULL, LA) = 0,1 agree -> no spurious failure.
  .ins_inh(conn, "1", "M", NA,   0, 1)
  .ins_inh(conn, "1", NA,  "LA", 0, 1)
  expect_silent(validate_chr_inheritance(conn))
})

test_that("validate_chr_recombination() checks its own key + shadowing", {
  pop <- make_pop_base(n_loci = 6, n_chr = 1, chr_len_Mb = 50)
  on.exit(close_pop(pop), add = TRUE)
  conn <- pop$db_conn
  DBI::dbExecute(conn,
    "INSERT INTO chr_recombination (chr_name, parent_sex, line_name, recombines) VALUES ('1', 'M', NULL, FALSE)")
  DBI::dbExecute(conn,
    "INSERT INTO chr_recombination (chr_name, parent_sex, line_name, recombines) VALUES ('1', NULL, 'LA', TRUE)")
  expect_error(validate_chr_recombination(conn), "shadowing conflict")
})
