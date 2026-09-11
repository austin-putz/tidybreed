test_that("restore_pop() produces a valid pop object from a freshly initialized DB", {
  tmp <- tempfile(fileext = ".duckdb")
  pop <- open_pop(pop_name = "test_pop", db_name = tmp) |>
    define_genome(n_loci = 100, n_chr = 2, chr_len_Mb = 50) |>
    define_founder_haplotypes(n_haplotypes = 30)
  close_pop(pop)

  pop2 <- restore_pop(tmp, pop_name = "test_pop")

  expect_s3_class(pop2, "tidybreed_pop")
  expect_equal(pop2$pop_name, "test_pop")
  expect_equal(pop2$db_path,  tmp)

  expect_equal(
    as.integer(DBI::dbGetQuery(pop2$db_conn, "SELECT COUNT(*) AS n FROM genome_meta")$n),
    100L
  )
  expect_equal(
    as.integer(DBI::dbGetQuery(pop2$db_conn,
      "SELECT COUNT(DISTINCT chr) AS n FROM genome_meta")$n),
    2L
  )
  expect_equal(
    DBI::dbGetQuery(pop2$db_conn,
      "SELECT chr_name FROM genome_meta GROUP BY chr_name ORDER BY MIN(chr)")$chr_name,
    c("1", "2")
  )
  expect_equal(
    as.integer(DBI::dbGetQuery(pop2$db_conn, "SELECT COUNT(*) AS n FROM ind_meta")$n),
    0L
  )

  close_pop(pop2)
  unlink(tmp)
})


test_that("restore_pop() returns pop with correct individual count in ind_meta", {
  tmp <- tempfile(fileext = ".duckdb")
  pop <- open_pop(pop_name = "count_test", db_name = tmp) |>
    define_genome(n_loci = 50, n_chr = 2, chr_len_Mb = 50) |>
    define_founder_haplotypes(n_haplotypes = 20)
  pop <- pop |>
    get_table("founder_haplotypes") |>
    add_founders(n_males = 10, n_females = 20, line_name = "A")
  close_pop(pop)

  pop2 <- restore_pop(tmp)

  expect_equal(
    as.integer(DBI::dbGetQuery(pop2$db_conn, "SELECT COUNT(*) AS n FROM ind_meta")$n),
    30L
  )

  close_pop(pop2)
  unlink(tmp)
})


test_that("restore_pop() handles heterogeneous chromosome lengths", {
  tmp <- tempfile(fileext = ".duckdb")
  pop <- open_pop(pop_name = "cattle", db_name = tmp) |>
    define_genome(n_loci = 90, n_chr = 3, chr_len_Mb = c(158, 137, 121)) |>
    define_founder_haplotypes(n_haplotypes = 20)
  close_pop(pop)

  pop2 <- restore_pop(tmp)

  chr_lengths <- DBI::dbGetQuery(pop2$db_conn,
    "SELECT chr, MAX(pos_bp) AS max_pos FROM genome_meta GROUP BY chr ORDER BY chr")

  expect_equal(nrow(chr_lengths), 3L)
  # Loci distributed proportional to chromosome length, so ordering is preserved
  expect_true(chr_lengths$max_pos[1] >= chr_lengths$max_pos[2])
  expect_true(chr_lengths$max_pos[2] >= chr_lengths$max_pos[3])

  close_pop(pop2)
  unlink(tmp)
})


test_that("add_founders() works after restore_pop() on empty DB", {
  tmp <- tempfile(fileext = ".duckdb")
  pop <- open_pop(pop_name = "fresh", db_name = tmp) |>
    define_genome(n_loci = 50, n_chr = 2, chr_len_Mb = 50) |>
    define_founder_haplotypes(n_haplotypes = 20)
  close_pop(pop)

  pop2 <- restore_pop(tmp)

  expect_no_error({
    pop2 <- pop2 |>
      get_table("founder_haplotypes") |>
      add_founders(n_males = 5, n_females = 5, line_name = "A")
  })

  ind_meta <- get_table(pop2, "ind_meta") |> dplyr::collect()
  expect_equal(nrow(ind_meta), 10L)

  close_pop(pop2)
  unlink(tmp)
})


test_that("add_offspring() works after restore; IDs continue without collision", {
  tmp <- tempfile(fileext = ".duckdb")
  pop <- open_pop(pop_name = "offspring_test", db_name = tmp) |>
    define_genome(n_loci = 50, n_chr = 2, chr_len_Mb = 50) |>
    define_founder_haplotypes(n_haplotypes = 20)
  pop <- pop |>
    get_table("founder_haplotypes") |>
    add_founders(n_males = 5, n_females = 5, line_name = "A")
  close_pop(pop)

  pop2 <- restore_pop(tmp)

  matings <- tibble::tibble(
    id_parent_1 = "A_1",
    id_parent_2 = "A_6",
    sex         = "M",
    line_name   = "A"
  )

  expect_no_error({
    pop2 <- add_offspring(pop2, matings)
  })

  ind_meta <- get_table(pop2, "ind_meta") |> dplyr::collect()
  expect_equal(nrow(ind_meta), 11L)
  expect_equal(length(unique(ind_meta$id_ind)), 11L)
  expect_true("A_11" %in% ind_meta$id_ind)

  close_pop(pop2)
  unlink(tmp)
})


test_that("pop_name inferred from _tidybreed.duckdb suffix", {
  tmp_dir  <- tempdir()
  db_path  <- file.path(tmp_dir, "mysim_tidybreed.duckdb")

  pop <- open_pop(pop_name = "mysim", db_name = db_path) |>
    define_genome(n_loci = 50, n_chr = 2, chr_len_Mb = 50) |>
    define_founder_haplotypes(n_haplotypes = 10)
  close_pop(pop)

  pop2 <- restore_pop(db_path)
  expect_equal(pop2$pop_name, "mysim")

  close_pop(pop2)
  unlink(db_path)
})


test_that("pop_name override parameter works", {
  tmp <- tempfile(fileext = ".duckdb")
  pop <- open_pop(pop_name = "original", db_name = tmp) |>
    define_genome(n_loci = 50, n_chr = 2, chr_len_Mb = 50) |>
    define_founder_haplotypes(n_haplotypes = 10)
  close_pop(pop)

  pop2 <- restore_pop(tmp, pop_name = "replicate_01")
  expect_equal(pop2$pop_name, "replicate_01")

  close_pop(pop2)
  unlink(tmp)
})


test_that("restore_pop() errors when file does not exist", {
  expect_error(
    restore_pop("/tmp/definitely_does_not_exist_xyz123.duckdb"),
    "does not exist"
  )
})


test_that("restore_pop() errors for ':memory:' path", {
  expect_error(
    restore_pop(":memory:"),
    "cannot restore in-memory databases"
  )
})


test_that("restore_pop() errors when genome_meta table is missing", {
  tmp  <- tempfile(fileext = ".duckdb")
  drv  <- duckdb::duckdb()
  conn <- DBI::dbConnect(drv, dbdir = tmp)
  DBI::dbExecute(conn, "CREATE TABLE dummy (x INTEGER)")
  DBI::dbDisconnect(conn, shutdown = TRUE)

  expect_error(
    restore_pop(tmp),
    "ind_meta"
  )

  unlink(tmp)
})


# ── Gate 49: a pre-term-model file stops here, not inside add_tbv() ─────────

test_that("restore_pop() refuses the pre-v0.65.0 genome_effects shape", {
  tmp <- tempfile(fileext = ".duckdb")
  pop <- open_pop(pop_name = "old_shape", db_name = tmp) |>
    define_genome(n_loci = 10, n_chr = 1, chr_len_Mb = 10)
  close_pop(pop)

  # Rebuild the old table in place: one row per (locus, line, trait), with
  # locus_name and genome_effect_type, and no member/origin children.
  conn <- DBI::dbConnect(duckdb::duckdb(), dbdir = tmp)
  for (vw in c("genome_effect_loci", "genome_effect_terms"))
    DBI::dbExecute(conn, paste0("DROP VIEW ", vw))
  for (tb in c("genome_effect_member_origins", "genome_effect_members",
               "genome_effects"))
    DBI::dbExecute(conn, paste0("DROP TABLE ", tb))
  DBI::dbExecute(conn, "
    CREATE TABLE genome_effects (
      id_genome_effect   INTEGER PRIMARY KEY,
      locus_name         VARCHAR NOT NULL,
      line_name          VARCHAR,
      trait_name         VARCHAR NOT NULL,
      genome_effect_type VARCHAR NOT NULL,
      genome_value       DOUBLE  NOT NULL,
      base_allele_freq   DOUBLE
    )")
  DBI::dbDisconnect(conn, shutdown = TRUE)

  expect_error(restore_pop(tmp), "pre-v0\\.65\\.0")

  # The refusal must close its connection, or the file stays locked and the
  # message is followed by an unopenable database.
  conn2 <- DBI::dbConnect(duckdb::duckdb(), dbdir = tmp)
  expect_true("genome_effects" %in% DBI::dbListTables(conn2))
  DBI::dbDisconnect(conn2, shutdown = TRUE)

  unlink(tmp)
})

test_that("restore_pop() refuses a term table whose member children are gone", {
  tmp <- tempfile(fileext = ".duckdb")
  pop <- open_pop(pop_name = "half_shape", db_name = tmp) |>
    define_genome(n_loci = 10, n_chr = 1, chr_len_Mb = 10)
  close_pop(pop)

  conn <- DBI::dbConnect(duckdb::duckdb(), dbdir = tmp)
  for (vw in c("genome_effect_loci", "genome_effect_terms"))
    DBI::dbExecute(conn, paste0("DROP VIEW ", vw))
  DBI::dbExecute(conn, "DROP TABLE genome_effect_member_origins")
  DBI::dbExecute(conn, "DROP TABLE genome_effect_members")
  DBI::dbDisconnect(conn, shutdown = TRUE)

  expect_error(restore_pop(tmp), "genome_effect_members")

  unlink(tmp)
})

test_that("restore_pop() refuses a pre-0.68.0 phenotype_components shape", {
  # Same class as gate 49, one release later: a file written by 0.65.0-0.67.0
  # has the term/member genome effects but the old reserved column name, and
  # would fail inside define_phenotype(components = ) rather than here.
  tmp <- tempfile(fileext = ".duckdb")
  # A genome is defined so the genome-effects guard passes first: this test is
  # about the second check, not about reaching it by accident.
  pop <- open_pop(pop_name = "old_comp", db_name = tmp) |>
    define_genome(n_loci = 10, n_chr = 1, chr_len_Mb = 10)
  close_pop(pop)

  conn <- DBI::dbConnect(duckdb::duckdb(), dbdir = tmp)
  DBI::dbExecute(conn, paste0("ALTER TABLE phenotype_components ",
                              "RENAME COLUMN component_names TO ",
                              "genome_effect_types"))
  DBI::dbDisconnect(conn, shutdown = TRUE)

  expect_error(restore_pop(tmp), "pre-v0\\.68\\.0")

  conn2 <- DBI::dbConnect(duckdb::duckdb(), dbdir = tmp)
  expect_true("phenotype_components" %in% DBI::dbListTables(conn2))
  DBI::dbDisconnect(conn2, shutdown = TRUE)

  unlink(tmp)
})
