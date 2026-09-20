# resolve_subset_ids(): the single place a tidybreed_table becomes the set of
# individuals an action function acts on. The contract is "the individuals
# present in the (filtered) table", whatever the table and whether or not a
# filter is pending.

make_subset_pop <- function(pop_name) {
  pop <- open_pop(pop_name = pop_name, db_name = ":memory:") |>
    define_genome(n_loci = 100, n_chr = 2, chr_len_Mb = 100) |>
    define_founder_haplotypes(n_haplotypes = 50, method = "fixed")
  pop <- pop |>
    get_table("founder_haplotypes") |>
    add_founders(n_males = 10, n_females = 10, line_name = "A", gen = 0L)
  pop <- define_trait(pop, "ADG", target_add_var = 1)
  pop <- pop |> get_table("genome_meta") |> dplyr::filter(chr == 1) |>
    define_additive_effects("ADG")
  pop
}

# Hand-written ind_ebv rows for the first `n` ids, ebv_value = 1..n
seed_ebv <- function(pop, ids) {
  DBI::dbExecute(pop$db_conn, paste0(
    "INSERT INTO ind_ebv (id_ebv, id_ind, trait_name, model, ebv_value, ",
    "eval_number) VALUES ",
    paste(sprintf("(%d, '%s', 'ADG', 'm1', %d, 1)",
                  seq_along(ids), ids, seq_along(ids)), collapse = ", ")))
  invisible(pop)
}

all_ids <- function(pop) {
  sort(DBI::dbGetQuery(pop$db_conn, "SELECT id_ind FROM ind_meta")$id_ind)
}


test_that("unfiltered ind_meta is NULL by default and every id on request", {
  pop <- make_subset_pop("rs_meta")
  tbl <- get_table(pop, "ind_meta")

  expect_null(resolve_subset_ids(tbl))
  expect_identical(resolve_subset_ids(tbl, all_if_null = TRUE), all_ids(pop))

  close_pop(pop)
})


test_that("filtered ind_meta returns exactly the matching ids, sorted", {
  pop <- make_subset_pop("rs_meta_filter")
  ids <- pop |> get_table("ind_meta") |> dplyr::filter(sex == "M") |>
    resolve_subset_ids()

  expect_length(ids, 10)
  expect_identical(ids, sort(ids))
  meta <- dplyr::collect(get_table(pop, "ind_meta"))
  expect_true(all(meta$sex[match(ids, meta$id_ind)] == "M"))

  close_pop(pop)
})


test_that("unfiltered non-ind_meta table selects only the animals in it", {
  pop <- make_subset_pop("rs_ebv_nofilter")
  ebv_ids <- head(all_ids(pop), 5)
  seed_ebv(pop, ebv_ids)

  ids <- resolve_subset_ids(get_table(pop, "ind_ebv"))
  expect_identical(ids, sort(ebv_ids))

  close_pop(pop)
})


test_that("filtered ind_ebv respects the predicate", {
  pop <- make_subset_pop("rs_ebv_filter")
  ebv_ids <- head(all_ids(pop), 5)
  seed_ebv(pop, ebv_ids)                         # ebv_value 1..5

  ids <- pop |> get_table("ind_ebv") |> dplyr::filter(ebv_value > 2) |>
    resolve_subset_ids()
  expect_identical(ids, sort(ebv_ids[3:5]))

  close_pop(pop)
})


test_that("ind_genotype filtered by dosage gives distinct carriers only", {
  pop <- make_subset_pop("rs_geno")
  pop <- pop |> get_table("ind_meta") |> add_dosage()

  carriers <- DBI::dbGetQuery(pop$db_conn,
    "SELECT DISTINCT id_ind FROM ind_genotype
     WHERE locus_name = 'Locus_1' AND dosage_value = 2")$id_ind
  ids <- pop |> get_table("ind_genotype") |>
    dplyr::filter(locus_name == "Locus_1", dosage_value == 2L) |>
    resolve_subset_ids()

  expect_identical(ids, sort(carriers))
  expect_false(anyDuplicated(ids) > 0)

  close_pop(pop)
})


test_that("ind_haplotype (many rows per animal) reduces to distinct ids", {
  pop <- make_subset_pop("rs_hap")
  ids <- pop |> get_table("ind_haplotype") |>
    dplyr::filter(parent_origin == 1L) |>
    resolve_subset_ids()

  expect_identical(ids, all_ids(pop))            # every animal has a sire copy
  expect_false(anyDuplicated(ids) > 0)

  close_pop(pop)
})


test_that("a table without id_ind errors, filtered or not", {
  pop <- make_subset_pop("rs_no_id")

  expect_error(resolve_subset_ids(get_table(pop, "genome_meta"), "phenotyping"),
               "has no 'id_ind' column.*phenotyping")
  expect_error(
    pop |> get_table("genome_meta") |> dplyr::filter(chr == 1) |>
      resolve_subset_ids(),
    "has no 'id_ind' column")

  close_pop(pop)
})


test_that("select() before filter() cannot hide id_ind", {
  pop <- make_subset_pop("rs_select")
  ids <- pop |> get_table("ind_meta") |> dplyr::select(sex) |>
    dplyr::filter(sex == "F") |>
    resolve_subset_ids()

  expect_length(ids, 10)

  close_pop(pop)
})


test_that("a filter matching nobody returns character(0)", {
  pop <- make_subset_pop("rs_empty")
  ids <- pop |> get_table("ind_meta") |> dplyr::filter(sex == "X") |>
    resolve_subset_ids()

  expect_identical(ids, character(0))

  close_pop(pop)
})


test_that("ids absent from ind_meta are dropped by the semi-join", {
  pop <- make_subset_pop("rs_orphan")
  ebv_ids <- head(all_ids(pop), 2)
  seed_ebv(pop, c(ebv_ids, "ghost_1"))

  ids <- resolve_subset_ids(get_table(pop, "ind_ebv"))
  expect_identical(ids, sort(ebv_ids))

  close_pop(pop)
})
