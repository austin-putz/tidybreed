# Sub-stage 1b checkpoint: the long-format tables and chr_meta exist with the
# correct shape alongside the (still-authoritative) wide tables.

test_that("define_genome creates long ind_haplotype/ind_genotype and chr_meta", {
  pop <- make_pop_base(n_loci = 10, n_chr = 2, chr_len_Mb = 100)
  on.exit(close_pop(pop), add = TRUE)

  for (tb in c("ind_haplotype", "ind_genotype", "chr_meta")) {
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

test_that("chr_meta is seeded with default diploid-autosome rows", {
  pop <- make_pop_base(n_loci = 12, n_chr = 3, chr_len_Mb = 100)
  on.exit(close_pop(pop), add = TRUE)

  cm <- DBI::dbGetQuery(pop$db_conn, "SELECT * FROM chr_meta ORDER BY chr_name")
  expect_equal(nrow(cm), 3)
  expect_true(all(cm$copy_mode_M == "full"))
  expect_true(all(cm$copy_mode_F == "full"))
  expect_true(all(cm$recombines))
  expect_true(all(is.na(cm$hemi_parent)))
})

test_that("genome_meta has introduced_gen (NULL for founding loci)", {
  pop <- make_pop_base(n_loci = 8, n_chr = 2, chr_len_Mb = 100)
  on.exit(close_pop(pop), add = TRUE)

  gm <- DBI::dbGetQuery(pop$db_conn, "SELECT * FROM genome_meta")
  expect_true("introduced_gen" %in% names(gm))
  expect_true(all(is.na(gm$introduced_gen)))
})

test_that("duplicate custom locus_names are rejected", {
  expect_error(
    open_pop(pop_name = "dup", db_name = ":memory:") |>
      define_genome(n_loci = 3, n_chr = 1, chr_len_Mb = 100,
                    locus_names = c("A", "B", "A")),
    "unique"
  )
})
