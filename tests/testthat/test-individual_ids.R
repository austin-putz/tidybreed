# Individual ids are "{line_name}_{n}", numbered per line from the largest n
# the line holds in ind_meta (.next_ind_number()). Line names may contain "_",
# and one line's name may be a prefix of another's. Parent ids reach DuckDB
# only as rows of a registered view. See plans/individual_id_sql.md.

ids_pop <- function(db_name = ":memory:") {
  set.seed(11)
  open_pop(pop_name = "ids", db_name = db_name) |>
    define_genome(n_loci = 20, n_chr = 2, chr_len_Mb = 50) |>
    define_founder_haplotypes(n_haplotypes = 10, line_name = "A")
}

founders <- function(pop, line, n_males = 1, n_females = 1) {
  suppressMessages(get_table(pop, "founder_haplotypes") |>
                     add_founders(n_males, n_females, line))
}

offspring <- function(pop, sire, dam, line, seed = 1, sex = "M") {
  suppressMessages(add_offspring(pop, data.frame(
    id_parent_1 = sire, id_parent_2 = dam, sex = sex, line_name = line,
    stringsAsFactors = FALSE), seed = seed))
}

line_ids <- function(pop, line) {
  DBI::dbGetQuery(pop$db_conn, "SELECT id_ind, line_name FROM ind_meta") |>
    (\(d) d$id_ind[d$line_name == line])()
}

n_of <- function(ids, line) sort(as.integer(substring(ids, nchar(line) + 2L)))

# A hand-made ind_meta row: the fixtures for malformed and boundary ids
put_ind <- function(pop, id, line, sex = "M") {
  DBI::dbExecute(pop$db_conn,
    "INSERT INTO ind_meta (id_ind, line_name, sex) VALUES (?, ?, ?)",
    params = list(id, line, sex))
}

row_counts <- function(pop) {
  c(meta = DBI::dbGetQuery(pop$db_conn, "SELECT COUNT(*) AS n FROM ind_meta")$n,
    hap  = DBI::dbGetQuery(pop$db_conn, "SELECT COUNT(*) AS n FROM ind_haplotype")$n)
}

# ---- numbering -------------------------------------------------------------

test_that("a line whose name contains '_' keeps numbering across both writers", {
  pop <- ids_pop()
  on.exit(close_pop(pop), add = TRUE)

  pop <- founders(pop, "Duroc_F1", 2, 2)                 # 1-4
  pop <- founders(pop, "Duroc_F1", 1, 1)                 # 5-6 (failed before)
  pop <- offspring(pop, "Duroc_F1_1", "Duroc_F1_3", "Duroc_F1")      # 7
  pop <- offspring(pop, c("Duroc_F1_2", "Duroc_F1_7"),
                   c("Duroc_F1_4", "Duroc_F1_5"), "Duroc_F1", seed = 2)  # 8-9
  pop <- founders(pop, "Duroc_F1", 1, 0)                 # 10

  ids <- line_ids(pop, "Duroc_F1")
  expect_identical(n_of(ids, "Duroc_F1"), 1:10)
  meta <- DBI::dbGetQuery(pop$db_conn,
    "SELECT id_ind, id_parent_1, id_parent_2 FROM ind_meta ORDER BY id_ind")
  expect_identical(meta$id_parent_1[meta$id_ind == "Duroc_F1_9"], "Duroc_F1_7")
  expect_identical(meta$id_parent_2[meta$id_ind == "Duroc_F1_9"], "Duroc_F1_5")
})

test_that("lines whose names are prefixes of each other number independently", {
  pop <- ids_pop()
  on.exit(close_pop(pop), add = TRUE)

  pop <- founders(pop, "A_B")
  pop <- founders(pop, "A")            # failed before: 'B_1' read as A's number
  pop <- founders(pop, "A_7")
  pop <- founders(pop, "A_")           # ids A__1, A__2
  pop <- founders(pop, "A")
  pop <- founders(pop, "A_B")

  # One call whose rows interleave all four lines: ids follow row order
  lines <- c("A_7", "A", "A_", "A_B", "A", "A_7", "A_")
  pop <- suppressMessages(add_offspring(pop, data.frame(
    id_parent_1 = "A_1", id_parent_2 = "A_2", sex = "F", line_name = lines),
    seed = 3))
  new <- DBI::dbGetQuery(pop$db_conn,
    "SELECT id_ind, line_name FROM ind_meta WHERE id_parent_1 IS NOT NULL")
  expect_setequal(new$id_ind, c("A_7_3", "A_5", "A__3", "A_B_5", "A_6",
                                "A_7_4", "A__4"))

  expect_identical(n_of(line_ids(pop, "A"),   "A"),   1:6)
  expect_identical(n_of(line_ids(pop, "A_B"), "A_B"), 1:5)
  expect_identical(n_of(line_ids(pop, "A_7"), "A_7"), 1:4)
  expect_identical(n_of(line_ids(pop, "A_"),  "A_"),  1:4)
  expect_true(all(startsWith(line_ids(pop, "A_"), "A__")))

  # Row order decides which offspring gets which number
  rank_ids <- DBI::dbGetQuery(pop$db_conn,
    "SELECT id_ind FROM ind_meta WHERE line_name = 'A_7' AND id_parent_1 IS NOT NULL")
  expect_setequal(rank_ids$id_ind, c("A_7_3", "A_7_4"))
})

test_that(".next_ind_number() maps results by line name", {
  pop <- ids_pop()
  on.exit(close_pop(pop), add = TRUE)
  pop <- founders(pop, "Zeta", 2, 1)        # 1-3
  pop <- founders(pop, "Alpha", 1, 0)       # 1
  pop <- founders(pop, "Mid_X", 3, 3)       # 1-6

  got <- .next_ind_number(pop$db_conn,
                          c(Zeta = 2, NewLine = 1, Mid_X = 4, Alpha = 1))
  expect_identical(got, c(Zeta = 4L, NewLine = 1L, Mid_X = 7L, Alpha = 2L))

  expect_error(.next_ind_number(pop$db_conn, c(Zeta = 1, Zeta = 2)), "unique")
  expect_error(.next_ind_number(pop$db_conn, c(Zeta = 0)), "positive")
})

test_that("a malformed stored id is an error naming the line and examples", {
  pop <- ids_pop()
  on.exit(close_pop(pop), add = TRUE)
  pop <- founders(pop, "Bad")
  for (id in c("Bad_extra_99", "Bad_0", "Bad_007")) put_ind(pop, id, "Bad")
  before <- row_counts(pop)

  err <- expect_error(founders(pop, "Bad"), "line 'Bad': 3")
  expect_match(conditionMessage(err), "Bad_0, Bad_007, Bad_extra_99", fixed = TRUE)
  expect_error(offspring(pop, "Bad_1", "Bad_2", "Bad"), "line 'Bad'")
  expect_identical(row_counts(pop), before)

  # Another line is unaffected
  expect_no_error(founders(pop, "Good"))
})

test_that("id numbers stay within 1..int32 max for the whole allocation", {
  pop <- ids_pop()
  on.exit(close_pop(pop), add = TRUE)
  lim <- .Machine$integer.max

  put_ind(pop, paste0("Lim_", lim - 1L), "Lim")
  expect_no_warning(pop <- founders(pop, "Lim", 1, 0))
  expect_true(paste0("Lim_", lim) %in% line_ids(pop, "Lim"))

  before <- row_counts(pop)
  expect_error(founders(pop, "Lim", 1, 0), "already holds an id number at or above")

  # An allocation that crosses the limit part-way
  put_ind(pop, paste0("Cross_", lim - 7L), "Cross")
  expect_error(founders(pop, "Cross", 5, 5), "line 'Cross' would need numbers")

  # A stored number above the range is never read as "no rows"
  put_ind(pop, "Huge_2147483648", "Huge")
  expect_error(founders(pop, "Huge", 1, 0), "line 'Huge' already holds")
  put_ind(pop, "Long_12345678901", "Long")
  expect_error(founders(pop, "Long", 1, 0), "line 'Long' already holds")

  # Doubles near the limit (the old seq() wrote past it)
  put_ind(pop, paste0("Dbl_", lim - 2L), "Dbl")
  expect_error(founders(pop, "Dbl", 2, 1), "line 'Dbl' would need numbers")

  # add_offspring() checks the same way, before writing anything
  pop <- founders(pop, "P", 1, 1)
  put_ind(pop, paste0("Off_", lim), "Off")
  before <- row_counts(pop)
  expect_error(offspring(pop, "P_1", "P_2", "Off"), "line 'Off' already holds")
  expect_identical(row_counts(pop), before)

  expect_error(founders(pop, "P", 1.5, 0), "whole numbers")
})

test_that("numbering continues after close_pop() and restore_pop()", {
  tmp <- tempfile(fileext = ".duckdb")
  on.exit(unlink(tmp), add = TRUE)
  pop <- ids_pop(tmp)
  pop <- founders(pop, "Duroc_F1", 1, 1)
  pop <- offspring(pop, "Duroc_F1_1", "Duroc_F1_2", "Duroc_F1")
  close_pop(pop)

  pop <- restore_pop(tmp, pop_name = "ids")
  pop <- founders(pop, "Duroc_F1", 1, 0)
  pop <- offspring(pop, "Duroc_F1_1", "Duroc_F1_2", "Duroc_F1", seed = 2)
  expect_identical(n_of(line_ids(pop, "Duroc_F1"), "Duroc_F1"), 1:5)
  close_pop(pop)
})

# ---- parent ids ------------------------------------------------------------

test_that("an SQL-shaped unknown parent id is reported as not found", {
  pop <- ids_pop()
  on.exit(close_pop(pop), add = TRUE)
  pop <- founders(pop, "A")
  before <- row_counts(pop)

  expect_error(offspring(pop, "A_1'); --", "A_2", "A"),
               "Parent ID(s) not found in ind_meta: A_1'); --", fixed = TRUE)
  expect_identical(row_counts(pop), before)
})

test_that("parent id columns must be non-missing character (factors converted)", {
  pop <- ids_pop()
  on.exit(close_pop(pop), add = TRUE)
  pop <- founders(pop, "A")

  base <- data.frame(id_parent_1 = "A_1", id_parent_2 = "A_2", sex = "M",
                     line_name = "A", stringsAsFactors = FALSE)
  for (col in c("id_parent_1", "id_parent_2")) {
    m <- base; m[[col]] <- NA_character_
    expect_error(add_offspring(pop, m, seed = 1), paste0("matings\\$", col, " has missing"))
    m <- base; m[[col]] <- ""
    expect_error(add_offspring(pop, m, seed = 1), paste0("matings\\$", col, " has missing"))
    m <- base; m[[col]] <- 1
    expect_error(add_offspring(pop, m, seed = 1),
                 paste0("matings\\$", col, " must be character"))
  }

  # A factor matings table gives the same offspring as a character one
  as_fac <- base
  as_fac[] <- lapply(as_fac, factor)
  pop <- suppressMessages(add_offspring(pop, as_fac, seed = 4))
  pop <- suppressMessages(add_offspring(pop, base, seed = 4))
  hap <- DBI::dbGetQuery(pop$db_conn, "
    SELECT id_ind, parent_origin, locus_id, allele FROM ind_haplotype
    WHERE id_ind IN (SELECT id_ind FROM ind_meta WHERE id_parent_1 IS NOT NULL)
    ORDER BY id_ind, parent_origin, locus_id")
  expect_identical(hap$allele[hap$id_ind == "A_3"], hap$allele[hap$id_ind == "A_4"])
})

test_that("no individual id appears in the SQL the writers execute", {
  pop <- ids_pop()
  on.exit(close_pop(pop), add = TRUE)

  # Record every statement sent through DBI, forwarding it unchanged
  stmts <- character()
  rec <- function(f) {
    force(f)
    function(conn, statement, ...) {
      stmts <<- c(stmts, paste(as.character(statement), collapse = "\n"))
      f(conn, statement, ...)
    }
  }
  testthat::local_mocked_bindings(
    dbGetQuery  = rec(DBI::dbGetQuery),
    dbExecute   = rec(DBI::dbExecute),
    dbSendQuery = rec(DBI::dbSendQuery),
    .package = "DBI")

  pop <- founders(pop, "Qzline_F1", 2, 2)
  pop <- founders(pop, "Qzline_F1", 1, 0)
  pop <- suppressMessages(add_offspring(pop, data.frame(
    id_parent_1 = c("Qzline_F1_1", "Qzline_F1_2"),
    id_parent_2 = c("Qzline_F1_3", "Qzline_F1_4"),
    sex = "F", line_name = "Qzline_F1"), seed = 5))
  expect_error(offspring(pop, "Qzline_F1_99", "Qzline_F1_3", "Qzline_F1"),
               "not found")

  expect_gt(length(stmts), 0L)
  # Both reads of a successful add_offspring() ran: ind_meta and ind_haplotype
  expect_true(any(grepl("JOIN __ao_parents", stmts) & grepl("ind_haplotype", stmts)))
  expect_true(any(grepl("JOIN __ao_parents", stmts) & grepl("ind_meta", stmts)))

  # No individual id, given or generated, in any statement
  expect_false(any(grepl("Qzline_F1_[0-9]", stmts)))
  # The line name reaches ind_meta only through the numbering view
  expect_false(any(grepl("Qzline", stmts) & grepl("ind_meta", stmts)))
})

test_that("internal views are removed after success and after failure", {
  pop <- ids_pop()
  on.exit(close_pop(pop), add = TRUE)
  views <- function() DBI::dbGetQuery(pop$db_conn,
    "SELECT view_name FROM duckdb_views() WHERE NOT internal")$view_name
  leftover <- function() intersect(c("__ao_parents", "__next_id_lines"),
                                   c(views(), DBI::dbListTables(pop$db_conn)))

  pop <- founders(pop, "A")
  pop <- offspring(pop, "A_1", "A_2", "A")
  expect_length(leftover(), 0L)

  expect_error(offspring(pop, "A_404", "A_2", "A"), "not found")
  expect_length(leftover(), 0L)
  put_ind(pop, "Bad_x", "Bad")
  expect_error(founders(pop, "Bad"), "line 'Bad'")
  expect_length(leftover(), 0L)

  expect_no_error(pop <- offspring(pop, "A_1", "A_3", "A", seed = 2))
  expect_identical(n_of(line_ids(pop, "A"), "A"), 1:4)
})

# ---- RNG -------------------------------------------------------------------

test_that("numbering and early add_offspring() errors draw no random numbers", {
  pop <- ids_pop()
  on.exit(close_pop(pop), add = TRUE)
  pop <- founders(pop, "A")

  set.seed(99)
  s <- .Random.seed
  .next_ind_number(pop$db_conn, c(A = 3, B = 1))
  expect_identical(.Random.seed, s)

  bad <- data.frame(id_parent_1 = NA_character_, id_parent_2 = "A_2",
                    sex = "M", line_name = "A")
  expect_error(add_offspring(pop, bad), "has missing")
  expect_identical(.Random.seed, s)

  expect_error(add_offspring(pop, transform(bad, id_parent_1 = "A_404")),
               "not found")
  expect_identical(.Random.seed, s)
})

test_that("seeded founders and offspring reproduce, with the expected ids", {
  run <- function() {
    pop <- ids_pop()
    pop <- founders(pop, "A", 2, 2)
    pop <- founders(pop, "B_x", 1, 1)
    pop <- suppressMessages(add_offspring(pop, data.frame(
      id_parent_1 = c("A_1", "B_x_1", "A_2"), id_parent_2 = c("B_x_2", "A_3", "A_4"),
      sex = c("M", "F", "M"), line_name = c("AB", "B_x", "AB")), seed = 6))
    out <- list(
      meta = DBI::dbGetQuery(pop$db_conn, "SELECT * FROM ind_meta ORDER BY id_ind"),
      hap  = DBI::dbGetQuery(pop$db_conn,
        "SELECT * FROM ind_haplotype ORDER BY id_ind, parent_origin, locus_id"))
    close_pop(pop)
    out
  }
  a <- run()
  b <- run()
  expect_identical(a, b)
  expect_identical(a$meta$id_ind,
                   c("AB_1", "AB_2", "A_1", "A_2", "A_3", "A_4",
                     "B_x_1", "B_x_2", "B_x_3"))
  expect_identical(a$meta$id_parent_1[a$meta$id_ind == "AB_2"], "A_2")
})
