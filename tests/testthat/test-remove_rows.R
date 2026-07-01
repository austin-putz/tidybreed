# Helper: small population with founders + one trait + phenotypes
make_pop_for_remove <- function(n_males = 3, n_females = 3) {
  pop <- make_test_pop(n_males = n_males, n_females = n_females)

  pop <- pop |>
    define_trait(trait_name = "ADG", target_add_var = 100) |>
    define_phenotype(phenotype_name = "ADG", mean = 500, residual_var = 50)

  # Assign additive effects to all loci before calling add_phenotype
  pop <- pop |>
    get_table("genome_meta") |>
    define_additive_effects("ADG")

  pop |>
    get_table("ind_meta") |>
    add_phenotype("ADG")
}


# ── 1. Input validation ───────────────────────────────────────────────────────

test_that("remove_rows() errors when tbl is not a tidybreed_table", {
  expect_error(remove_rows(list()), "'tbl' must be a tidybreed_table")
})


# ── 2. No-filter guard ────────────────────────────────────────────────────────

test_that("remove_rows() stops when no filter applied and confirm_all = FALSE", {
  pop <- make_pop_for_remove()
  expect_error(
    get_table(pop, "ind_phenotype") |> remove_rows(),
    regexp = "no filter applied"
  )
  close_pop(pop)
})

test_that("remove_rows() error message includes row count and table name", {
  pop <- make_pop_for_remove()
  err <- tryCatch(
    get_table(pop, "ind_phenotype") |> remove_rows(),
    error = function(e) conditionMessage(e)
  )
  expect_match(err, "ind_phenotype")
  expect_match(err, "confirm_all = TRUE")
  close_pop(pop)
})

test_that("remove_rows(confirm_all = TRUE) deletes all rows without a filter", {
  pop <- make_pop_for_remove()
  n_before <- DBI::dbGetQuery(pop$db_conn,
    "SELECT COUNT(*) AS n FROM ind_phenotype")$n
  expect_gt(n_before, 0L)

  pop <- get_table(pop, "ind_phenotype") |> remove_rows(confirm_all = TRUE)

  n_after <- DBI::dbGetQuery(pop$db_conn,
    "SELECT COUNT(*) AS n FROM ind_phenotype")$n
  expect_equal(n_after, 0L)
  close_pop(pop)
})


# ── 3. 0-row filter warning ───────────────────────────────────────────────────

test_that("remove_rows() warns and returns pop unchanged when filter matches 0 rows", {
  pop <- make_pop_for_remove()
  n_before <- DBI::dbGetQuery(pop$db_conn,
    "SELECT COUNT(*) AS n FROM ind_phenotype")$n

  expect_warning(
    pop2 <- get_table(pop, "ind_phenotype") |>
      dplyr::filter(id_ind == "DOES_NOT_EXIST") |>
      remove_rows(),
    regexp = "matched 0 rows"
  )

  n_after <- DBI::dbGetQuery(pop2$db_conn,
    "SELECT COUNT(*) AS n FROM ind_phenotype")$n
  expect_equal(n_after, n_before)
  close_pop(pop)
})


# ── 4. Single-table mode ──────────────────────────────────────────────────────

test_that("remove_rows() deletes only filtered ind_phenotype rows (single-table)", {
  pop <- make_pop_for_remove()
  all_ids <- dplyr::collect(get_table(pop, "ind_meta"))$id_ind
  keep_id  <- all_ids[1]
  remove_ids <- all_ids[-1]

  pop <- get_table(pop, "ind_phenotype") |>
    dplyr::filter(id_ind %in% remove_ids) |>
    remove_rows(verbose = FALSE)

  remaining <- dplyr::collect(get_table(pop, "ind_phenotype"))
  expect_true(all(remaining$id_ind == keep_id))
  close_pop(pop)
})

test_that("remove_rows() in single-table mode does not touch other tables", {
  pop <- make_pop_for_remove()
  all_ids <- dplyr::collect(get_table(pop, "ind_meta"))$id_ind
  remove_id <- all_ids[1]

  n_meta_before <- DBI::dbGetQuery(pop$db_conn,
    "SELECT COUNT(*) AS n FROM ind_meta")$n

  pop <- get_table(pop, "ind_phenotype") |>
    dplyr::filter(id_ind == remove_id) |>
    remove_rows(verbose = FALSE)

  n_meta_after <- DBI::dbGetQuery(pop$db_conn,
    "SELECT COUNT(*) AS n FROM ind_meta")$n
  expect_equal(n_meta_after, n_meta_before)
  close_pop(pop)
})

test_that("remove_rows() deletes ind_tbv rows filtered by trait_name only", {
  pop <- make_pop_for_remove()

  # Add a second trait so we can check only ADG is deleted
  pop <- pop |>
    define_trait(trait_name = "BW", target_add_var = 50) |>
    define_phenotype(phenotype_name = "BW", mean = 100, residual_var = 10)
  pop <- pop |>
    get_table("genome_meta") |>
    define_additive_effects("BW")
  pop <- get_table(pop, "ind_meta") |> add_phenotype("BW")

  n_bw_before <- DBI::dbGetQuery(pop$db_conn,
    "SELECT COUNT(*) AS n FROM ind_tbv WHERE trait_name = 'BW'")$n
  expect_gt(n_bw_before, 0L)

  pop <- get_table(pop, "ind_tbv") |>
    dplyr::filter(trait_name == "ADG") |>
    remove_rows(verbose = FALSE)

  n_adg_after <- DBI::dbGetQuery(pop$db_conn,
    "SELECT COUNT(*) AS n FROM ind_tbv WHERE trait_name = 'ADG'")$n
  n_bw_after  <- DBI::dbGetQuery(pop$db_conn,
    "SELECT COUNT(*) AS n FROM ind_tbv WHERE trait_name = 'BW'")$n

  expect_equal(n_adg_after, 0L)
  expect_equal(n_bw_after, n_bw_before)
  close_pop(pop)
})


# ── 5. Cross-table mode ───────────────────────────────────────────────────────

test_that("remove_rows(tables = 'all') removes individual from all ind_* tables", {
  pop <- make_pop_for_remove()
  all_ids   <- dplyr::collect(get_table(pop, "ind_meta"))$id_ind
  remove_id <- all_ids[1]

  pop <- get_table(pop, "ind_meta") |>
    dplyr::filter(id_ind == remove_id) |>
    remove_rows(tables = "all", verbose = FALSE)

  for (tbl_name in c("ind_meta", "ind_phenotype", "ind_tbv",
                      "genome_haplotype", "genome_genotype")) {
    n <- DBI::dbGetQuery(
      pop$db_conn,
      paste0("SELECT COUNT(*) AS n FROM ", tbl_name,
             " WHERE id_ind = '", remove_id, "'")
    )$n
    expect_equal(n, 0L, label = paste(tbl_name, "row count for removed id"))
  }

  # Other individuals untouched
  n_remaining <- DBI::dbGetQuery(pop$db_conn,
    "SELECT COUNT(*) AS n FROM ind_meta")$n
  expect_equal(n_remaining, length(all_ids) - 1L)
  close_pop(pop)
})

test_that("remove_rows(tables = c(...)) removes from only listed tables", {
  pop <- make_pop_for_remove()
  all_ids   <- dplyr::collect(get_table(pop, "ind_meta"))$id_ind
  remove_id <- all_ids[1]

  n_meta_before <- DBI::dbGetQuery(pop$db_conn,
    "SELECT COUNT(*) AS n FROM ind_meta")$n

  pop <- get_table(pop, "ind_meta") |>
    dplyr::filter(id_ind == remove_id) |>
    remove_rows(tables = c("ind_phenotype", "ind_tbv"), verbose = FALSE)

  # ind_phenotype and ind_tbv cleaned
  n_pheno <- DBI::dbGetQuery(pop$db_conn,
    paste0("SELECT COUNT(*) AS n FROM ind_phenotype WHERE id_ind = '",
           remove_id, "'"))$n
  n_tbv   <- DBI::dbGetQuery(pop$db_conn,
    paste0("SELECT COUNT(*) AS n FROM ind_tbv WHERE id_ind = '",
           remove_id, "'"))$n
  expect_equal(n_pheno, 0L)
  expect_equal(n_tbv, 0L)

  # ind_meta untouched (not in tables list)
  n_meta_after <- DBI::dbGetQuery(pop$db_conn,
    "SELECT COUNT(*) AS n FROM ind_meta")$n
  expect_equal(n_meta_after, n_meta_before)
  close_pop(pop)
})

test_that("remove_rows() errors when tables arg names a table without id_ind", {
  pop <- make_pop_for_remove()
  target_id <- dplyr::collect(get_table(pop, "ind_meta"))$id_ind[[1]]
  expect_error(
    get_table(pop, "ind_meta") |>
      dplyr::filter(id_ind == target_id) |>
      remove_rows(tables = "trait_meta"),
    regexp = "do not have an 'id_ind' column"
  )
  close_pop(pop)
})

test_that("remove_rows() errors when tables arg names a non-existent table", {
  pop <- make_pop_for_remove()
  target_id <- dplyr::collect(get_table(pop, "ind_meta"))$id_ind[[1]]
  expect_error(
    get_table(pop, "ind_meta") |>
      dplyr::filter(id_ind == target_id) |>
      remove_rows(tables = "no_such_table"),
    regexp = "do not exist in this population"
  )
  close_pop(pop)
})

test_that("remove_rows() errors for cross-table mode when source has no id_ind", {
  pop <- make_pop_for_remove()
  expect_error(
    get_table(pop, "trait_meta") |>
      dplyr::filter(trait_name == "ADG") |>
      remove_rows(tables = "all"),
    regexp = "id_ind"
  )
  close_pop(pop)
})


# ── 6. dry_run ────────────────────────────────────────────────────────────────

test_that("remove_rows(dry_run = TRUE) does not modify the database", {
  pop <- make_pop_for_remove()
  all_ids   <- dplyr::collect(get_table(pop, "ind_meta"))$id_ind
  remove_id <- all_ids[1]

  n_before <- DBI::dbGetQuery(pop$db_conn,
    "SELECT COUNT(*) AS n FROM ind_phenotype")$n

  expect_message(
    get_table(pop, "ind_phenotype") |>
      dplyr::filter(id_ind == remove_id) |>
      remove_rows(dry_run = TRUE),
    regexp = "\\[dry_run\\]"
  )

  n_after <- DBI::dbGetQuery(pop$db_conn,
    "SELECT COUNT(*) AS n FROM ind_phenotype")$n
  expect_equal(n_after, n_before)
  close_pop(pop)
})


# ── 7. verbose ────────────────────────────────────────────────────────────────

test_that("remove_rows(verbose = TRUE) emits a deletion message", {
  pop <- make_pop_for_remove()
  target_id <- dplyr::collect(get_table(pop, "ind_meta"))$id_ind[[1]]

  expect_message(
    get_table(pop, "ind_phenotype") |>
      dplyr::filter(id_ind == target_id) |>
      remove_rows(verbose = TRUE),
    regexp = "Deleted"
  )
  close_pop(pop)
})

test_that("remove_rows(verbose = FALSE) emits no message", {
  pop <- make_pop_for_remove()
  target_id <- dplyr::collect(get_table(pop, "ind_meta"))$id_ind[[1]]

  expect_no_message(
    get_table(pop, "ind_phenotype") |>
      dplyr::filter(id_ind == target_id) |>
      remove_rows(verbose = FALSE)
  )
  close_pop(pop)
})


# ── 8. Return value ───────────────────────────────────────────────────────────

test_that("remove_rows() returns the pop invisibly", {
  pop <- make_pop_for_remove()
  target_id <- dplyr::collect(get_table(pop, "ind_meta"))$id_ind[[1]]

  result <- withVisible(
    get_table(pop, "ind_phenotype") |>
      dplyr::filter(id_ind == target_id) |>
      remove_rows(verbose = FALSE)
  )
  expect_false(result$visible)
  expect_s3_class(result$value, "tidybreed_pop")
  close_pop(pop)
})
