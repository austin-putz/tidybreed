# Sub-stage 1b checkpoint: the long-format tables and the chromosome-rule tables
# exist with the correct shape.

test_that("define_genome creates long ind_haplotype/ind_genotype and chr tables", {
  pop <- make_pop_base(n_loci = 10, n_chr = 2, chr_len_Mb = 100)
  on.exit(close_pop(pop), add = TRUE)

  for (tb in c("ind_haplotype", "ind_genotype",
               "chr_inheritance", "chr_recombination")) {
    expect_true(tb %in% pop$tables)
    expect_true(DBI::dbExistsTable(pop$db_conn, tb))
  }

  hap_cols <- DBI::dbListFields(pop$db_conn, "ind_haplotype")
  expect_setequal(
    hap_cols,
    c("id_ind", "parent_origin", "strand", "line_origin",
      "locus_id", "locus_name", "allele")
  )
  geno_cols <- DBI::dbListFields(pop$db_conn, "ind_genotype")
  expect_setequal(
    geno_cols,
    c("id_ind", "locus_id", "locus_name", "dosage_value")
  )

  # Long tables start empty (ind_genotype is on-demand; ind_haplotype not yet
  # populated).
  expect_equal(DBI::dbGetQuery(pop$db_conn,
    "SELECT COUNT(*) AS n FROM ind_haplotype")$n, 0)
  expect_equal(DBI::dbGetQuery(pop$db_conn,
    "SELECT COUNT(*) AS n FROM ind_genotype")$n, 0)
})

test_that("chr_inheritance/chr_recombination are seeded with default rows", {
  pop <- make_pop_base(n_loci = 12, n_chr = 3, chr_len_Mb = 100)
  on.exit(close_pop(pop), add = TRUE)

  inh <- DBI::dbGetQuery(pop$db_conn,
    "SELECT * FROM chr_inheritance ORDER BY chr_name")
  expect_equal(nrow(inh), 3)
  expect_true(all(is.na(inh$offspring_sex)))
  expect_true(all(is.na(inh$line_name)))
  expect_true(all(inh$from_parent_1 == 1))
  expect_true(all(inh$from_parent_2 == 1))

  rec <- DBI::dbGetQuery(pop$db_conn,
    "SELECT * FROM chr_recombination ORDER BY chr_name")
  expect_equal(nrow(rec), 3)
  expect_true(all(is.na(rec$parent_sex)))
  expect_true(all(is.na(rec$line_name)))
  expect_true(all(rec$recombines))
})

test_that("genome_meta stores pos_bp (BIGINT) and has no legacy columns", {
  pop <- make_pop_base(n_loci = 8, n_chr = 2, chr_len_Mb = 100)
  on.exit(close_pop(pop), add = TRUE)

  gm <- DBI::dbGetQuery(pop$db_conn, "SELECT * FROM genome_meta")
  expect_true("pos_bp" %in% names(gm))
  expect_false("pos_Mb"         %in% names(gm))
  expect_false("pos_cM"         %in% names(gm))  # genetic map is in genome_map
  expect_false("introduced_gen" %in% names(gm))
})

test_that("duplicate custom locus_names are rejected", {
  expect_error(
    open_pop(pop_name = "dup", db_name = ":memory:") |>
      define_genome(n_loci = 3, n_chr = 1, chr_len_Mb = 100,
                    locus_names = c("A", "B", "A")),
    "unique"
  )
})
