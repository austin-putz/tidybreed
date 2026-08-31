# Covers the schema() header, pipeline grouping, and the print-method
# refinements from plans/update_schema_print.md (Parts 2-6).

# Header line of a printed schema object.
schema_header <- function(pop, ...) {
  capture.output(print(schema(pop, ...)))[1]
}


# ── Part 5a: whole-database size in the header ────────────────────────────────

test_that("header reports on-disk size for a file-backed population", {
  dir  <- withr::local_tempdir()
  pop  <- open_pop(pop_name = "Disk", db_name = file.path(dir, "sz.duckdb")) |>
    define_genome(n_loci = 50, n_chr = 1, chr_len_Mb = 10)
  on.exit(close_pop(pop), add = TRUE)

  hdr <- schema_header(pop)
  expect_match(hdr, "Schema: Disk")
  expect_match(hdr, "on disk")
  expect_false(grepl("in memory", hdr, fixed = TRUE))
})

test_that("header reports memory usage for an in-memory population", {
  # PRAGMA database_size returns database_size "0 bytes" and block_size 0 when
  # there is no file, so reporting it as an on-disk figure would read as zero.
  pop <- open_pop(pop_name = "Mem", db_name = ":memory:") |>
    define_genome(n_loci = 50, n_chr = 1, chr_len_Mb = 10)
  on.exit(close_pop(pop), add = TRUE)

  hdr <- schema_header(pop)
  expect_match(hdr, "in memory")
  expect_false(grepl("on disk", hdr, fixed = TRUE))
})

test_that("header counts the tables actually shown", {
  pop <- make_test_pop(n_males = 2, n_females = 2, n_loci = 30, n_chr = 1)
  on.exit(close_pop(pop), add = TRUE)

  s <- schema(pop)
  expect_match(capture.output(print(s))[1], paste0(nrow(s), " tables"))
})

test_that("uncheckpointed WAL is reported alongside the file size", {
  # Before a CHECKPOINT the whole population can sit in the WAL with
  # database_size reading "0 bytes"; printing only that would be wrong.
  dir <- withr::local_tempdir()
  pop <- open_pop(pop_name = "Wal", db_name = file.path(dir, "wal.duckdb")) |>
    define_genome(n_loci = 200, n_chr = 2, chr_len_Mb = 100)
  on.exit(close_pop(pop), add = TRUE)

  expect_match(schema_header(pop), "WAL")

  DBI::dbExecute(pop$db_conn, "CHECKPOINT")
  hdr <- schema_header(pop)
  expect_match(hdr, "on disk")
  expect_false(grepl("WAL", hdr, fixed = TRUE))
})

test_that("schema() itself does not checkpoint", {
  # A silent CHECKPOINT inside a display function is a hidden write. The WAL
  # must still be non-empty after schema() has been called and printed.
  dir <- withr::local_tempdir()
  pop <- open_pop(pop_name = "NoCkpt", db_name = file.path(dir, "nc.duckdb")) |>
    define_genome(n_loci = 200, n_chr = 2, chr_len_Mb = 100)
  on.exit(close_pop(pop), add = TRUE)

  invisible(capture.output(print(schema(pop))))
  wal <- DBI::dbGetQuery(pop$db_conn, "PRAGMA database_size")$wal_size[1]
  expect_false(grepl("^0 +bytes$", wal))
})


# ── Parts 2-4: pipeline grouping and print refinements ────────────────────────

test_that("ordering is identical whether the population was built or restored", {
  # The real defect Part 1 describes: open_pop() fills pop$tables in creation
  # order and restore_pop() in DuckDB catalog order, so the same .duckdb file
  # printed in a different order depending on how it was opened.
  dir <- withr::local_tempdir()
  f   <- file.path(dir, "ord.duckdb")

  pop <- open_pop(pop_name = "Ord", db_name = f) |>
    define_genome(n_loci = 40, n_chr = 2, chr_len_Mb = 50)
  built <- schema(pop, show_empty = TRUE)$table_name
  close_pop(pop)

  pop2 <- restore_pop(f)
  on.exit(close_pop(pop2), add = TRUE)
  restored <- schema(pop2, show_empty = TRUE)$table_name

  expect_identical(built, restored)
})

test_that("every system table is placed in a group, none fall through", {
  # If a new table is added to SYSTEM_TABLES without registering it in
  # .schema_table_order(), it shows up under "User tables" and this fails.
  placed <- unlist(.schema_table_order(), use.names = FALSE)
  expect_setequal(placed, SYSTEM_TABLES)
})

test_that("the group vector and the description helpers name the same tables", {
  described <- unique(.all_schema_descriptions()$table_name)
  expect_setequal(described, unlist(.schema_table_order(), use.names = FALSE))
})

test_that("tables are returned in pipeline-group order with a table_group column", {
  pop <- make_test_pop(n_males = 2, n_females = 2, n_loci = 30, n_chr = 1)
  on.exit(close_pop(pop), add = TRUE)

  s <- schema(pop, show_empty = TRUE)
  expect_true("table_group" %in% names(s))
  expect_s3_class(s$table_group, "factor")

  # Group blocks are contiguous and in the declared display order.
  grp <- as.character(s$table_group)
  expect_identical(grp, rep(rle(grp)$values, rle(grp)$lengths))
  expect_identical(rle(grp)$values,
                   intersect(names(.schema_table_order()), grp))

  # Within Genome, in-group order is workflow order, not alphabetical.
  genome <- s$table_name[s$table_group == "Genome"]
  expect_identical(genome[1:2], c("genome_meta", "genome_map"))
})

test_that("unregistered tables land under User tables, sorted by name", {
  pop <- make_test_pop(n_males = 2, n_females = 2, n_loci = 30, n_chr = 1)
  on.exit(close_pop(pop), add = TRUE)

  pop <- define_table(pop, "zzz_notes",   columns = c(note = "VARCHAR"))
  pop <- define_table(pop, "aaa_scratch", columns = c(x = "INTEGER"))

  s    <- schema(pop, show_empty = TRUE)
  user <- s$table_name[s$table_group == "User tables"]
  expect_identical(user, c("aaa_scratch", "zzz_notes"))
  expect_identical(tail(as.character(s$table_group), 1L), "User tables")
})

test_that("_schema_meta is hidden by default and shown by include_system", {
  pop <- make_test_pop(n_males = 2, n_females = 2, n_loci = 30, n_chr = 1)
  on.exit(close_pop(pop), add = TRUE)

  expect_false("_schema_meta" %in% schema(pop)$table_name)

  s <- schema(pop, include_system = TRUE)
  expect_true("_schema_meta" %in% s$table_name)
  expect_identical(as.character(s$table_group[s$table_name == "_schema_meta"]),
                   "System")
  expect_identical(tail(as.character(s$table_group), 1L), "System")
})

test_that("empty tables collapse by default and expand with show_empty", {
  pop <- make_test_pop(n_males = 2, n_females = 2, n_loci = 30, n_chr = 1)
  on.exit(close_pop(pop), add = TRUE)

  collapsed <- capture.output(print(schema(pop)))
  expect_true(any(grepl("^\\s+\\+ [0-9]+ empty: ", collapsed)))
  expect_false(any(grepl("^\\s+ind_genotype\\s", collapsed)))

  expanded <- capture.output(print(schema(pop, show_empty = TRUE)))
  expect_false(any(grepl("empty:", expanded)))
  expect_true(any(grepl("^\\s+ind_genotype\\s", expanded)))

  # The tibble always carries every table; collapsing is print-time only.
  expect_identical(nrow(schema(pop)), nrow(schema(pop, show_empty = TRUE)))
})

test_that("group headings print only for groups that have tables", {
  pop <- make_test_pop(n_males = 2, n_females = 2, n_loci = 30, n_chr = 1)
  on.exit(close_pop(pop), add = TRUE)

  out <- capture.output(print(schema(pop)))
  expect_true(any(grepl("^  Genome$", out)))
  # No traits defined and no user tables created in this population.
  expect_false(any(grepl("^  User tables$", out)))
})

test_that("row counts above a million are abbreviated", {
  expect_identical(.schema_format_rows(c(0L, 999L, 12345L)),
                   c("0", "999", "12,345"))
  expect_identical(.schema_format_rows(1100000), "1.10M")
  expect_identical(.schema_format_rows(2500000000), "2.50B")
  expect_identical(.schema_format_rows(NA_integer_), "?")
})


# ── Parts 5b-6: order and sizes ───────────────────────────────────────────────

test_that("order = 'name' and 'rows' produce flat orderings without headings", {
  pop <- make_test_pop(n_males = 2, n_females = 2, n_loci = 30, n_chr = 1)
  on.exit(close_pop(pop), add = TRUE)

  by_name <- schema(pop, order = "name", show_empty = TRUE)
  expect_identical(by_name$table_name, sort(by_name$table_name))

  by_rows <- schema(pop, order = "rows", show_empty = TRUE)
  expect_false(is.unsorted(rev(by_rows$n_rows)))

  # table_group survives in the tibble even when headings are not printed.
  expect_true("table_group" %in% names(by_name))
  out <- capture.output(print(by_name))
  expect_false(any(grepl("^  Genome$", out)))
  expect_true(any(grepl("^\\s+Table\\s+Rows\\s+Cols", out)))
})

test_that("order = 'size' errors unless sizes = TRUE", {
  pop <- make_test_pop(n_males = 2, n_females = 2, n_loci = 30, n_chr = 1)
  on.exit(close_pop(pop), add = TRUE)

  # Undefined, not merely unavailable — so this errors rather than falling back.
  expect_error(schema(pop, order = "size"), "requires sizes = TRUE")
})

test_that("sizes = TRUE on an in-memory population warns and omits the column", {
  pop <- make_test_pop(n_males = 2, n_females = 2, n_loci = 30, n_chr = 1)
  on.exit(close_pop(pop), add = TRUE)

  expect_warning(s <- schema(pop, sizes = TRUE), "in memory")
  expect_false("size_bytes" %in% names(s))

  # Ordering by size is undefined there, so that combination still errors.
  expect_error(schema(pop, order = "size", sizes = TRUE), "in-memory")
})

test_that("sizes = TRUE adds size_bytes and prints the caveat footnote", {
  dir <- withr::local_tempdir()
  pop <- open_pop(pop_name = "Sz", db_name = file.path(dir, "sz.duckdb")) |>
    define_genome(n_loci = 500, n_chr = 2, chr_len_Mb = 100)
  on.exit(close_pop(pop), add = TRUE)

  s <- schema(pop, sizes = TRUE)
  expect_true("size_bytes" %in% names(s))
  expect_type(s$size_bytes, "double")
  # Attribution is in whole blocks.
  block <- DBI::dbGetQuery(pop$db_conn, "PRAGMA database_size")$block_size[1]
  expect_true(all(s$size_bytes %% block == 0))

  out <- paste(capture.output(print(s)), collapse = " ")
  expect_match(out, "256 KiB blocks")            # the quantization caveat
  expect_match(out, "sum to less than")          # the shortfall caveat
  expect_match(out, "Size")                      # the column itself
})

test_that("order = 'size' sorts by bytes descending", {
  # Big enough that ind_haplotype spans more than one block: at ~256 KiB
  # granularity a small population leaves most tables tied at one block, and a
  # test that passed on the tie-break would not be testing the ordering.
  dir <- withr::local_tempdir()
  pop <- open_pop(pop_name = "Sz2", db_name = file.path(dir, "sz2.duckdb")) |>
    define_genome(n_loci = 1000, n_chr = 2, chr_len_Mb = 100)
  pop <- pop |> define_founder_haplotypes(n_haplotypes = 40) |>
    get_table("founder_haplotypes") |>
    add_founders(n_males = 30, n_females = 30, line_name = "A")
  on.exit(close_pop(pop), add = TRUE)

  s <- schema(pop, order = "size", sizes = TRUE, show_empty = TRUE)
  expect_false(is.unsorted(rev(s$size_bytes)))
  expect_identical(s$size_bytes[1], max(s$size_bytes))

  block <- DBI::dbGetQuery(pop$db_conn, "PRAGMA database_size")$block_size[1]
  hap   <- s$size_bytes[s$table_name == "ind_haplotype"]
  expect_gt(hap, block)                     # genuinely larger, not a tie-break
  expect_identical(s$table_name[1], "ind_haplotype")
})

test_that("sizes = FALSE issues no CHECKPOINT", {
  dir <- withr::local_tempdir()
  pop <- open_pop(pop_name = "NoCk2", db_name = file.path(dir, "nc2.duckdb")) |>
    define_genome(n_loci = 200, n_chr = 2, chr_len_Mb = 100)
  on.exit(close_pop(pop), add = TRUE)

  invisible(schema(pop))
  wal <- DBI::dbGetQuery(pop$db_conn, "PRAGMA database_size")$wal_size[1]
  expect_false(grepl("^0 +bytes$", wal))

  # ...and sizes = TRUE does, which is why it is opt-in.
  invisible(schema(pop, sizes = TRUE))
  wal2 <- DBI::dbGetQuery(pop$db_conn, "PRAGMA database_size")$wal_size[1]
  expect_match(wal2, "^0 +bytes$")
})
