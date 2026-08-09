# Tests for add_dosage(): the on-demand ind_genotype dosage cache.

make_dosage_pop <- function(n_loci = 30) {
  open_pop(pop_name = "dose", db_name = ":memory:") |>
    define_genome(n_loci = n_loci, n_chr = 2, chr_len_Mb = 100) |>
    define_founder_haplotypes(n_haplotypes = 20) |>
    get_table("founder_haplotypes") |>
    add_founders(n_males = 4, n_females = 4, line_name = "A", gen = 0L)
}

# Dosage derived directly from ind_haplotype (the ground truth).
expected_dosage <- function(pop) {
  d <- DBI::dbGetQuery(pop$db_conn,
    "SELECT id_ind, locus_id, CAST(SUM(allele) AS INTEGER) AS dosage_value
     FROM ind_haplotype GROUP BY id_ind, locus_id
     ORDER BY id_ind, locus_id")
  d
}

test_that("add_dosage materializes correct dosages for all loci", {
  pop <- make_dosage_pop(n_loci = 30)
  on.exit(close_pop(pop), add = TRUE)

  pop <- pop |> get_table("ind_meta") |> add_dosage()

  got <- DBI::dbGetQuery(pop$db_conn,
    "SELECT id_ind, locus_id, CAST(dosage_value AS INTEGER) AS dosage_value
     FROM ind_genotype ORDER BY id_ind, locus_id")
  exp <- expected_dosage(pop)

  expect_equal(nrow(got), 8L * 30L)          # 8 founders x 30 loci
  expect_equal(got, exp)
  expect_true(all(got$dosage_value %in% c(0L, 1L, 2L)))
})

test_that("add_dosage is idempotent (re-running replaces, does not duplicate)", {
  pop <- make_dosage_pop(n_loci = 20)
  on.exit(close_pop(pop), add = TRUE)

  pop <- pop |> get_table("ind_meta") |> add_dosage()
  n1  <- DBI::dbGetQuery(pop$db_conn, "SELECT COUNT(*) AS n FROM ind_genotype")$n
  pop <- pop |> get_table("ind_meta") |> add_dosage()
  n2  <- DBI::dbGetQuery(pop$db_conn, "SELECT COUNT(*) AS n FROM ind_genotype")$n

  expect_equal(n1, n2)
  expect_equal(n2, 8L * 20L)
})

test_that("add_dosage restricts to a filtered individual set", {
  pop <- make_dosage_pop(n_loci = 15)
  on.exit(close_pop(pop), add = TRUE)

  pop <- pop |>
    get_table("ind_meta") |>
    dplyr::filter(sex == "M") |>
    add_dosage()

  ids <- DBI::dbGetQuery(pop$db_conn,
    "SELECT DISTINCT id_ind FROM ind_genotype")$id_ind
  males <- DBI::dbGetQuery(pop$db_conn,
    "SELECT id_ind FROM ind_meta WHERE sex = 'M'")$id_ind
  expect_setequal(ids, males)
})

test_that("add_dosage restricts to locus_names", {
  pop <- make_dosage_pop(n_loci = 40)
  on.exit(close_pop(pop), add = TRUE)

  pop <- pop |>
    get_table("ind_meta") |>
    add_dosage(locus_names = c("Locus_1", "Locus_5", "Locus_9"))

  ln <- DBI::dbGetQuery(pop$db_conn,
    "SELECT DISTINCT locus_name FROM ind_genotype")$locus_name
  expect_setequal(ln, c("Locus_1", "Locus_5", "Locus_9"))
})

test_that("add_dosage chip_name filters by is_<chip>; unknown chip errors", {
  pop <- make_dosage_pop(n_loci = 20)
  on.exit(close_pop(pop), add = TRUE)

  # Define a chip on the first 5 loci.
  chip_loci <- paste0("Locus_", 1:5)
  pop <- pop |>
    get_table("genome_meta") |>
    dplyr::filter(locus_name %in% chip_loci) |>
    define_chip("mychip")

  pop <- pop |> get_table("ind_meta") |> add_dosage(chip_name = "mychip")
  ln <- DBI::dbGetQuery(pop$db_conn,
    "SELECT DISTINCT locus_name FROM ind_genotype")$locus_name
  expect_setequal(ln, chip_loci)

  expect_error(
    pop |> get_table("ind_meta") |> add_dosage(chip_name = "nope"),
    "not found in genome_meta"
  )
})

test_that("overwrite_dosage clears prior rows for the candidate set", {
  pop <- make_dosage_pop(n_loci = 20)
  on.exit(close_pop(pop), add = TRUE)

  # Materialize all loci, then overwrite with just 3 loci.
  pop <- pop |> get_table("ind_meta") |> add_dosage()
  expect_equal(DBI::dbGetQuery(pop$db_conn,
    "SELECT COUNT(*) AS n FROM ind_genotype")$n, 8L * 20L)

  pop <- pop |>
    get_table("ind_meta") |>
    add_dosage(locus_names = c("Locus_1", "Locus_2"), overwrite_dosage = TRUE)

  # Prior rows cleared; only the 2 loci remain.
  ln <- DBI::dbGetQuery(pop$db_conn,
    "SELECT DISTINCT locus_name FROM ind_genotype")$locus_name
  expect_setequal(ln, c("Locus_1", "Locus_2"))
  expect_equal(DBI::dbGetQuery(pop$db_conn,
    "SELECT COUNT(*) AS n FROM ind_genotype")$n, 8L * 2L)
})

test_that("add_dosage errors before writing anything when an individual has ploidy != 2", {
  pop <- make_dosage_pop(n_loci = 10)
  on.exit(close_pop(pop), add = TRUE)

  ids <- DBI::dbGetQuery(pop$db_conn, "SELECT id_ind FROM ind_meta LIMIT 1")$id_ind
  DBI::dbExecute(pop$db_conn, paste0(
    "UPDATE ind_meta SET ploidy = 4 WHERE id_ind = '", ids, "'"))

  expect_error(pop |> get_table("ind_meta") |> add_dosage(), "ploidy")
  expect_equal(DBI::dbGetQuery(pop$db_conn,
    "SELECT COUNT(*) AS n FROM ind_genotype")$n, 0L)
})

test_that("add_dosage does not error when chr_inheritance has a sex-linked row", {
  pop <- make_dosage_pop(n_loci = 10)
  on.exit(close_pop(pop), add = TRUE)

  suppressMessages(
    pop <- define_chromosome(pop, "1", offspring_sex = "M",
                             from_parent_1 = 0, from_parent_2 = 1))

  expect_no_error(pop |> get_table("ind_meta") |> add_dosage())
})
