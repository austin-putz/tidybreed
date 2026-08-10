test_that("define_genome() adds genome tables to pop$tables", {
  pop <- open_pop(db_name = ":memory:") |>
    define_genome(n_loci = 100, n_chr = 2, chr_len_Mb = 100)
  on.exit(close_pop(pop))

  expect_true("genome_meta"   %in% pop$tables)
  expect_true("genome_map"    %in% pop$tables)
  expect_true("ind_haplotype" %in% pop$tables)
  expect_true("ind_genotype"  %in% pop$tables)
  expect_true("chr_inheritance"   %in% pop$tables)
  expect_true("chr_recombination" %in% pop$tables)
})

test_that("define_genome() genome_meta has correct row count and columns", {
  n <- 200
  pop <- open_pop(db_name = ":memory:") |>
    define_genome(n_loci = n, n_chr = 4, chr_len_Mb = 50)
  on.exit(close_pop(pop))

  gm <- get_table(pop, "genome_meta") |> dplyr::collect()
  expect_equal(nrow(gm), n)
  # Physical coordinate is pos_bp (BIGINT); genetic pos_cM lives in genome_map.
  expect_true(all(c("locus_id", "locus_name", "chr", "chr_name", "pos_bp") %in%
                    names(gm)))
  # Old columns are gone.
  expect_false("pos_Mb"         %in% names(gm))
  expect_false("pos_cM"         %in% names(gm))
  expect_false("introduced_gen" %in% names(gm))

  # pos_bp is a genuine BIGINT column (not DOUBLE/INTEGER).
  pt <- DBI::dbGetQuery(pop$db_conn,
    "SELECT data_type FROM information_schema.columns
     WHERE table_name = 'genome_meta' AND column_name = 'pos_bp'")$data_type
  expect_equal(pt, "BIGINT")
})

test_that("define_genome() distributes loci evenly across chromosomes", {
  pop <- open_pop(db_name = ":memory:") |>
    define_genome(n_loci = 100, n_chr = 4, chr_len_Mb = 100)
  on.exit(close_pop(pop))

  gm <- get_table(pop, "genome_meta") |> dplyr::collect()
  chr_counts <- table(gm$chr)
  expect_equal(length(chr_counts), 4)
  # Each chromosome should have approximately 25 loci
  expect_true(all(chr_counts >= 24 & chr_counts <= 26))
})

test_that("define_genome() scalar chr_len_Mb applies to all chromosomes", {
  pop <- open_pop(db_name = ":memory:") |>
    define_genome(n_loci = 100, n_chr = 3, chr_len_Mb = 150)
  on.exit(close_pop(pop))

  gm <- get_table(pop, "genome_meta") |> dplyr::collect()
  for (c in 1:3) {
    chr_pos <- gm$pos_bp[gm$chr == c]        # base pairs now
    expect_true(max(chr_pos) <= 150 * 1e6)
    expect_true(min(chr_pos) > 0)            # 1-based, strictly positive
    expect_false(anyDuplicated(chr_pos) > 0) # strictly increasing, no dup
  }
})

test_that("define_genome() vector chr_len_Mb respected per chromosome", {
  lengths <- c(50, 100, 200)
  pop <- open_pop(db_name = ":memory:") |>
    define_genome(n_loci = 90, n_chr = 3, chr_len_Mb = lengths)
  on.exit(close_pop(pop))

  gm <- get_table(pop, "genome_meta") |> dplyr::collect()
  for (i in seq_along(lengths)) {
    chr_pos <- gm$pos_bp[gm$chr == i]        # base pairs
    expect_true(max(chr_pos) <= lengths[i] * 1e6)
  }
})

test_that("define_genome() custom locus_names are stored correctly", {
  names_vec <- paste0("SNP_", 1:50)
  pop <- open_pop(db_name = ":memory:") |>
    define_genome(n_loci = 50, n_chr = 2, chr_len_Mb = 100,
                  locus_names = names_vec)
  on.exit(close_pop(pop))

  gm <- get_table(pop, "genome_meta") |> dplyr::collect()
  expect_equal(sort(gm$locus_name), sort(names_vec))
})

test_that("define_genome() custom chr_names are stored correctly", {
  chr_names_vec <- c("chr1", "chr2", "chrX")
  pop <- open_pop(db_name = ":memory:") |>
    define_genome(n_loci = 60, n_chr = 3, chr_len_Mb = 100,
                  chr_names = chr_names_vec)
  on.exit(close_pop(pop))

  gm <- get_table(pop, "genome_meta") |> dplyr::collect()
  expect_equal(sort(unique(gm$chr_name)), sort(chr_names_vec))
})

test_that("define_genome() ind_haplotype and ind_genotype are empty tables", {
  pop <- open_pop(db_name = ":memory:") |>
    define_genome(n_loci = 50, n_chr = 2, chr_len_Mb = 100)
  on.exit(close_pop(pop))

  n_hap  <- DBI::dbGetQuery(pop$db_conn,
                             "SELECT COUNT(*) AS n FROM ind_haplotype")$n
  n_geno <- DBI::dbGetQuery(pop$db_conn,
                             "SELECT COUNT(*) AS n FROM ind_genotype")$n
  expect_equal(n_hap,  0L)
  expect_equal(n_geno, 0L)
})

test_that("define_genome() registers genome table descriptions in _schema_meta", {
  pop <- open_pop(db_name = ":memory:") |>
    define_genome(n_loci = 50, n_chr = 2, chr_len_Mb = 100)
  on.exit(close_pop(pop))

  sm <- DBI::dbGetQuery(
    pop$db_conn,
    "SELECT DISTINCT table_name FROM _schema_meta WHERE object_type = 'table'"
  )$table_name
  expect_true("genome_meta"   %in% sm)
  expect_true("ind_haplotype" %in% sm)
  expect_true("ind_genotype"  %in% sm)
  expect_true("chr_inheritance"   %in% sm)
  expect_true("chr_recombination" %in% sm)
})

test_that("define_genome() is pipe-friendly and returns pop", {
  pop <- open_pop(db_name = ":memory:") |>
    define_genome(n_loci = 50, n_chr = 2, chr_len_Mb = 100)
  on.exit(close_pop(pop))

  expect_s3_class(pop, "tidybreed_pop")
})

test_that("define_genome() rejects duplicate/NA/empty chr_names before writing genome_meta", {
  pop <- open_pop(db_name = ":memory:")
  on.exit(close_pop(pop))

  # Duplicate names error, and the failure is EARLY — genome_meta is never
  # created, so a corrected re-run on the same pop succeeds (no half-genome that
  # trips the "genome already defined" guard).
  expect_error(
    define_genome(pop, n_loci = 30, n_chr = 3, chr_len_Mb = 10,
                  chr_names = c("1", "1", "X")),
    "chr_names must be unique")
  expect_false(DBI::dbExistsTable(pop$db_conn, "genome_meta"))

  expect_error(
    define_genome(pop, n_loci = 20, n_chr = 2, chr_len_Mb = 10,
                  chr_names = c("1", "")),
    "non-missing, non-empty")
  expect_error(
    define_genome(pop, n_loci = 20, n_chr = 2, chr_len_Mb = 10,
                  chr_names = c("1", NA_character_)),
    "non-missing, non-empty")

  # Corrected re-run on the SAME pop works — nothing was left behind.
  expect_no_error(suppressMessages(
    define_genome(pop, n_loci = 30, n_chr = 3, chr_len_Mb = 10,
                  chr_names = c("1", "2", "X"))))
})

# ---------------------------------------------------------------------------
# Argument validation — every rejection must happen BEFORE any table is written
# ---------------------------------------------------------------------------

# The seven tables define_genome() owns. A failed call must leave none of them.
genome_tables_present <- function(conn) {
  intersect(
    c("genome_meta", "genome_map", "ind_haplotype", "ind_genotype",
      "ind_crossover", "chr_inheritance", "chr_recombination"),
    DBI::dbListTables(conn)
  )
}

test_that("define_genome() rejects non-whole / NA / Inf / non-numeric n_loci and n_chr", {
  pop <- open_pop(db_name = ":memory:")
  on.exit(close_pop(pop))

  # n_loci = 100.5 used to survive this check and blow up ~170 lines later inside
  # tibble() — seq.int(length.out = 100.5) gives 101 rows but seq_len(100.5) gives
  # 100 — after genome_meta had already been written.
  expect_error(define_genome(pop, n_loci = 100.5, n_chr = 2, chr_len_Mb = 100),
               "`n_loci` must be a whole number")
  expect_error(define_genome(pop, n_loci = NA, n_chr = 2, chr_len_Mb = 100),
               "`n_loci` must be numeric")
  expect_error(define_genome(pop, n_loci = NA_real_, n_chr = 2, chr_len_Mb = 100),
               "`n_loci` must not be NA")
  expect_error(define_genome(pop, n_loci = Inf, n_chr = 2, chr_len_Mb = 100),
               "`n_loci` must be finite")
  expect_error(define_genome(pop, n_loci = "10", n_chr = 2, chr_len_Mb = 100),
               "`n_loci` must be numeric")
  expect_error(define_genome(pop, n_loci = 0, n_chr = 1, chr_len_Mb = 100),
               "`n_loci` must be at least 1")

  expect_error(define_genome(pop, n_loci = 10, n_chr = 2.5, chr_len_Mb = 100),
               "`n_chr` must be a whole number")
  expect_error(define_genome(pop, n_loci = 10, n_chr = c(1, 2), chr_len_Mb = 100),
               "`n_chr` must be a single value")

  expect_error(define_genome(pop, n_loci = 10, n_chr = 1, chr_len_Mb = 10,
                             recombines_M = NA),
               "`recombines_M` must be TRUE or FALSE")
  expect_error(define_genome(pop, n_loci = 10, n_chr = 1, chr_len_Mb = 10,
                             recombines_F = "yes"),
               "`recombines_F` must be TRUE or FALSE")

  expect_length(genome_tables_present(pop$db_conn), 0L)
})

test_that("define_genome() requires n_loci >= n_chr", {
  pop <- open_pop(db_name = ":memory:")
  on.exit(close_pop(pop))

  # Fewer loci than chromosomes leaves some chromosomes empty, which made the
  # pos_bp index `[2:(n_chr_loci + 1)]` collapse to the reversed `[2:1]` and
  # reported a nonsense "physical span ... too short to place 0 loci" error.
  err <- expect_error(define_genome(pop, n_loci = 3, n_chr = 5, chr_len_Mb = 100))
  expect_match(conditionMessage(err),
               "n_loci \\(3\\) must be at least n_chr \\(5\\)")
  expect_no_match(conditionMessage(err), "too short")

  expect_length(genome_tables_present(pop$db_conn), 0L)

  # n_loci == n_chr is the boundary and must be allowed (one locus per chromosome).
  expect_no_error(suppressMessages(
    define_genome(pop, n_loci = 5, n_chr = 5, chr_len_Mb = 100)))
  expect_equal(
    DBI::dbGetQuery(pop$db_conn, "SELECT COUNT(*) AS n FROM genome_meta")$n, 5L)
})

test_that("define_genome() rejects NA / empty / wrong-length locus_names before writing", {
  pop <- open_pop(db_name = ":memory:")
  on.exit(close_pop(pop))

  # locus_name is the denormalized join key into genome_effects / ind_haplotype /
  # genome_map. NA and "" used to be written straight through and silently drop
  # those loci from every downstream join.
  expect_error(
    define_genome(pop, n_loci = 4, n_chr = 1, chr_len_Mb = 100,
                  locus_names = c("L1", NA_character_, "L3", "L4")),
    "locus_names must be non-missing, non-empty")
  expect_error(
    define_genome(pop, n_loci = 4, n_chr = 1, chr_len_Mb = 100,
                  locus_names = c("L1", "", "L3", "L4")),
    "locus_names must be non-missing, non-empty")
  expect_error(
    define_genome(pop, n_loci = 4, n_chr = 1, chr_len_Mb = 100,
                  locus_names = c("L1", "L2")),
    "locus_names must have length n_loci \\(4\\), got 2")

  expect_length(genome_tables_present(pop$db_conn), 0L)

  # Corrected re-run on the SAME pop works — nothing was left behind.
  expect_no_error(suppressMessages(
    define_genome(pop, n_loci = 4, n_chr = 1, chr_len_Mb = 100,
                  locus_names = c("L1", "L2", "L3", "L4"))))
})

test_that("define_genome() preflight covers all seven genome tables", {
  # An empty-but-existing genome table used to pass the old "genome_meta is
  # non-empty" preflight, and the plain CREATE then failed mid-flight leaving a
  # populated genome_meta behind.
  for (tbl in c("genome_meta", "ind_haplotype", "chr_recombination")) {
    pop <- open_pop(db_name = ":memory:")
    DBI::dbExecute(pop$db_conn, paste0("CREATE TABLE ", tbl, " (x INTEGER)"))

    expect_error(define_genome(pop, n_loci = 10, n_chr = 1, chr_len_Mb = 10),
                 paste0("Genome already defined.*", tbl))
    # Nothing new was created alongside the planted table.
    expect_equal(genome_tables_present(pop$db_conn), tbl)

    close_pop(pop)
  }
})

test_that("define_genome() rejects a pop whose connection is closed", {
  pop <- open_pop(db_name = ":memory:")
  close_pop(pop)
  expect_error(define_genome(pop, n_loci = 10, n_chr = 1, chr_len_Mb = 10),
               "connection is no longer valid")
})

# ---------------------------------------------------------------------------
# Atomicity
# ---------------------------------------------------------------------------

test_that("define_genome() rolls back completely when a validator fails mid-flight", {
  pop <- open_pop(db_name = ":memory:")
  on.exit(close_pop(pop))

  # Argument validation is now tight enough that no ordinary input reaches a
  # mid-transaction failure, so force one. This exercises the ROLLBACK path and
  # the duckdb_unregister() cleanup (registered views are session-level and are
  # NOT undone by ROLLBACK).
  testthat::local_mocked_bindings(
    validate_genome_map = function(...) stop("simulated validator failure"))

  expect_error(define_genome(pop, n_loci = 10, n_chr = 2, chr_len_Mb = 100),
               "simulated validator failure")

  # All seven tables gone, the _schema_meta descriptions rolled back with them,
  # and no orphaned temp views left registered on the connection.
  expect_length(genome_tables_present(pop$db_conn), 0L)
  expect_equal(
    DBI::dbGetQuery(pop$db_conn, paste0(
      "SELECT COUNT(*) AS n FROM _schema_meta ",
      "WHERE table_name IN ('genome_meta', 'genome_map', 'ind_crossover')"))$n,
    0L)
  expect_length(
    grep("^__tmp",
         DBI::dbGetQuery(pop$db_conn, "SELECT view_name FROM duckdb_views()")$view_name),
    0L)
})

test_that("define_genome() succeeds on the same pop after a rolled-back attempt", {
  pop <- open_pop(db_name = ":memory:")
  on.exit(close_pop(pop))

  local({
    testthat::local_mocked_bindings(
      validate_genome_map = function(...) stop("simulated validator failure"))
    expect_error(define_genome(pop, n_loci = 10, n_chr = 2, chr_len_Mb = 100))
  })

  # The whole point of the transaction: the population is still usable.
  pop <- suppressMessages(define_genome(pop, n_loci = 10, n_chr = 2, chr_len_Mb = 100))

  expect_length(genome_tables_present(pop$db_conn), 7L)
  expect_equal(
    DBI::dbGetQuery(pop$db_conn, "SELECT COUNT(*) AS n FROM genome_meta")$n, 10L)
  expect_equal(
    DBI::dbGetQuery(pop$db_conn, "SELECT COUNT(*) AS n FROM genome_map")$n, 10L)
})

test_that("define_genome() registers each genome table in pop$tables exactly once", {
  pop <- open_pop(db_name = ":memory:") |>
    define_genome(n_loci = 20, n_chr = 2, chr_len_Mb = 100)
  on.exit(close_pop(pop))

  expect_equal(anyDuplicated(pop$tables), 0L)
  expect_true(all(genome_tables_present(pop$db_conn) %in% pop$tables))
  # pop$tables must not drift from the database.
  expect_length(setdiff(pop$tables, DBI::dbListTables(pop$db_conn)), 0L)
})
