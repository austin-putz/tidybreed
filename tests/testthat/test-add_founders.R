test_that("add_founders creates ind_meta table with correct structure", {
  pop <- open_pop(pop_name = "test", db_name = ":memory:") |>
    define_genome(n_loci = 100, n_chr = 5, chr_len_Mb = 100) |>
    define_founder_haplotypes(n_haplotypes = 50)

  pop <- pop |>
    get_table("founder_haplotypes") |>
    add_founders(n_males = 10, n_females = 20, line_name = "A")

  expect_true("ind_meta" %in% pop$tables)
  expect_true("ind_meta" %in% DBI::dbListTables(pop$db_conn))

  ind_meta <- get_table(pop, "ind_meta") %>% dplyr::collect()

  expect_equal(nrow(ind_meta), 30)
  expect_equal(ncol(ind_meta), 6)
  expect_true(all(c("id_ind", "id_parent_1", "id_parent_2", "line_name", "sex", "ploidy") %in% colnames(ind_meta)))

  expect_type(ind_meta$id_ind, "character")
  expect_type(ind_meta$id_parent_1, "character")
  expect_type(ind_meta$id_parent_2, "character")
  expect_type(ind_meta$line_name, "character")
  expect_type(ind_meta$sex, "character")
  expect_true(all(ind_meta$ploidy == 2L))

  close_pop(pop)
})


test_that("add_founders correctly assigns IDs and sex", {
  pop <- open_pop(pop_name = "test", db_name = ":memory:") |>
    define_genome(n_loci = 50, n_chr = 2, chr_len_Mb = 100) |>
    define_founder_haplotypes(n_haplotypes = 20)

  pop <- pop |>
    get_table("founder_haplotypes") |>
    add_founders(n_males = 5, n_females = 10, line_name = "A")

  ind_meta <- get_table(pop, "ind_meta") %>% dplyr::collect()

  expected_ids <- paste0("A_", 1:15)
  expect_equal(ind_meta$id_ind, expected_ids)

  expect_true(all(ind_meta$line_name == "A"))

  expect_equal(sum(ind_meta$sex == "M"), 5)
  expect_equal(sum(ind_meta$sex == "F"), 10)
  expect_equal(ind_meta$sex[1:5], rep("M", 5))
  expect_equal(ind_meta$sex[6:15], rep("F", 10))

  expect_true(all(is.na(ind_meta$id_parent_1)))
  expect_true(all(is.na(ind_meta$id_parent_2)))

  close_pop(pop)
})


test_that("add_founders creates correct long ind_haplotype table", {
  pop <- open_pop(pop_name = "test", db_name = ":memory:") |>
    define_genome(n_loci = 50, n_chr = 2, chr_len_Mb = 100) |>
    define_founder_haplotypes(n_haplotypes = 20)

  pop <- pop |>
    get_table("founder_haplotypes") |>
    add_founders(n_males = 5, n_females = 5, line_name = "A")

  haps <- get_table(pop, "ind_haplotype") %>% dplyr::collect()

  # 10 individuals x 2 haplotypes x 50 loci
  expect_equal(nrow(haps), 10L * 2L * 50L)
  expect_setequal(colnames(haps),
    c("id_ind", "parent_origin", "strand", "line_origin",
      "locus_id", "locus_name", "allele"))
  expect_true(all(haps$parent_origin %in% c(1L, 2L)))
  expect_true(all(haps$strand == 1L))
  expect_true(all(haps$line_origin == "A"))
  expect_true(all(haps$allele %in% c(0L, 1L)))

  # Exactly 2 haplotypes (parent_origin 1 and 2) x 50 loci per individual.
  for (i in 1:10) {
    ind_haps <- haps %>% dplyr::filter(id_ind == paste0("A_", i))
    expect_equal(nrow(ind_haps), 100L)
    expect_setequal(unique(ind_haps$parent_origin), c(1L, 2L))
  }

  close_pop(pop)
})


test_that("add_founders leaves ind_genotype empty (on-demand cache)", {
  pop <- open_pop(pop_name = "test", db_name = ":memory:") |>
    define_genome(n_loci = 50, n_chr = 2, chr_len_Mb = 100) |>
    define_founder_haplotypes(n_haplotypes = 20)

  pop <- pop |>
    get_table("founder_haplotypes") |>
    add_founders(n_males = 5, n_females = 5, line_name = "A")

  expect_equal(DBI::dbGetQuery(pop$db_conn,
    "SELECT COUNT(*) AS n FROM ind_genotype")$n, 0)

  # Dosage derived from haplotypes is 0/1/2.
  dose <- DBI::dbGetQuery(pop$db_conn,
    "SELECT SUM(allele) AS d FROM ind_haplotype GROUP BY id_ind, locus_id")
  expect_true(all(dose$d %in% c(0, 1, 2)))

  close_pop(pop)
})


test_that("dosage equals sum of the two haplotype alleles", {
  pop <- open_pop(pop_name = "test", db_name = ":memory:") |>
    define_genome(n_loci = 50, n_chr = 2, chr_len_Mb = 100) |>
    define_founder_haplotypes(n_haplotypes = 20)

  pop <- pop |>
    get_table("founder_haplotypes") |>
    add_founders(n_males = 5, n_females = 5, line_name = "A")

  haps <- get_table(pop, "ind_haplotype") %>% dplyr::collect()
  agg <- stats::aggregate(allele ~ id_ind + locus_id, data = haps, FUN = sum)
  # Two haplotype rows per (id_ind, locus): dosage in 0..2.
  expect_true(all(agg$allele %in% c(0L, 1L, 2L)))
  # Every individual x locus is present exactly once after aggregation.
  expect_equal(nrow(agg), 10L * 50L)

  close_pop(pop)
})


test_that("add_founders works with multiple lines", {
  pop <- open_pop(pop_name = "test", db_name = ":memory:") |>
    define_genome(n_loci = 50, n_chr = 2, chr_len_Mb = 100) |>
    define_founder_haplotypes(n_haplotypes = 20)

  pop <- pop |>
    get_table("founder_haplotypes") |>
    add_founders(n_males = 5, n_females = 5, line_name = "A")

  pop <- pop |>
    get_table("founder_haplotypes") |>
    add_founders(n_males = 3, n_females = 7, line_name = "B")

  ind_meta <- get_table(pop, "ind_meta") %>% dplyr::collect()

  expect_equal(nrow(ind_meta), 20)

  line_a <- ind_meta %>% dplyr::filter(line_name == "A")
  expect_equal(nrow(line_a), 10)
  expect_equal(line_a$id_ind, paste0("A_", 1:10))

  line_b <- ind_meta %>% dplyr::filter(line_name == "B")
  expect_equal(nrow(line_b), 10)
  expect_equal(line_b$id_ind, paste0("B_", 1:10))

  expect_equal(sum(line_a$sex == "M"), 5)
  expect_equal(sum(line_a$sex == "F"), 5)
  expect_equal(sum(line_b$sex == "M"), 3)
  expect_equal(sum(line_b$sex == "F"), 7)

  close_pop(pop)
})


test_that("sequential additions to same line continue numbering", {
  pop <- open_pop(pop_name = "test", db_name = ":memory:") |>
    define_genome(n_loci = 50, n_chr = 2, chr_len_Mb = 100) |>
    define_founder_haplotypes(n_haplotypes = 20)

  pop <- pop |>
    get_table("founder_haplotypes") |>
    add_founders(n_males = 5, n_females = 5, line_name = "A")

  pop <- pop |>
    get_table("founder_haplotypes") |>
    add_founders(n_males = 2, n_females = 3, line_name = "A")

  ind_meta <- get_table(pop, "ind_meta") %>% dplyr::collect()

  expect_equal(nrow(ind_meta), 15)
  expect_equal(ind_meta$id_ind, paste0("A_", 1:15))
  expect_true(all(ind_meta$line_name == "A"))

  close_pop(pop)
})


test_that("add_founders works with only males", {
  pop <- open_pop(pop_name = "test", db_name = ":memory:") |>
    define_genome(n_loci = 50, n_chr = 2, chr_len_Mb = 100) |>
    define_founder_haplotypes(n_haplotypes = 20)

  pop <- pop |>
    get_table("founder_haplotypes") |>
    add_founders(n_males = 10, n_females = 0, line_name = "A")

  ind_meta <- get_table(pop, "ind_meta") %>% dplyr::collect()

  expect_equal(nrow(ind_meta), 10)
  expect_true(all(ind_meta$sex == "M"))

  close_pop(pop)
})


test_that("add_founders works with only females", {
  pop <- open_pop(pop_name = "test", db_name = ":memory:") |>
    define_genome(n_loci = 50, n_chr = 2, chr_len_Mb = 100) |>
    define_founder_haplotypes(n_haplotypes = 20)

  pop <- pop |>
    get_table("founder_haplotypes") |>
    add_founders(n_males = 0, n_females = 10, line_name = "A")

  ind_meta <- get_table(pop, "ind_meta") %>% dplyr::collect()

  expect_equal(nrow(ind_meta), 10)
  expect_true(all(ind_meta$sex == "F"))

  close_pop(pop)
})


test_that("get_table errors if founder_haplotypes table does not exist", {
  pop <- open_pop(pop_name = "test", db_name = ":memory:") |>
    define_genome(n_loci = 50, n_chr = 2, chr_len_Mb = 100)
    # Note: no define_founder_haplotypes() call

  expect_error(
    get_table(pop, "founder_haplotypes"),
    "founder_haplotypes"
  )

  close_pop(pop)
})


test_that("add_founders errors if first arg is not a tidybreed_table", {
  expect_error(
    add_founders(list(), n_males = 5, n_females = 5, line_name = "A"),
    "tidybreed_table"
  )
})


test_that("add_founders errors if n_males + n_females = 0", {
  pop <- open_pop(pop_name = "test", db_name = ":memory:") |>
    define_genome(n_loci = 50, n_chr = 2, chr_len_Mb = 100) |>
    define_founder_haplotypes(n_haplotypes = 20)

  expect_error(
    pop |>
      get_table("founder_haplotypes") |>
      add_founders(n_males = 0, n_females = 0, line_name = "A"),
    "At least one founder must be specified"
  )

  close_pop(pop)
})


test_that("add_founders errors if line_name has invalid format", {
  pop <- open_pop(pop_name = "test", db_name = ":memory:") |>
    define_genome(n_loci = 50, n_chr = 2, chr_len_Mb = 100) |>
    define_founder_haplotypes(n_haplotypes = 20)

  expect_error(
    pop |>
      get_table("founder_haplotypes") |>
      add_founders(n_males = 5, n_females = 5, line_name = "1A"),
    "line_name must start with letter"
  )

  expect_error(
    pop |>
      get_table("founder_haplotypes") |>
      add_founders(n_males = 5, n_females = 5, line_name = "A B"),
    "line_name must start with letter"
  )

  close_pop(pop)
})


test_that("add_founders errors if n_males or n_females invalid", {
  pop <- open_pop(pop_name = "test", db_name = ":memory:") |>
    define_genome(n_loci = 50, n_chr = 2, chr_len_Mb = 100) |>
    define_founder_haplotypes(n_haplotypes = 20)

  expect_error(
    pop |>
      get_table("founder_haplotypes") |>
      add_founders(n_males = -1, n_females = 5, line_name = "A")
  )

  expect_error(
    pop |>
      get_table("founder_haplotypes") |>
      add_founders(n_males = 5, n_females = -1, line_name = "A")
  )

  expect_error(
    pop |>
      get_table("founder_haplotypes") |>
      add_founders(n_males = "five", n_females = 5, line_name = "A")
  )

  close_pop(pop)
})


test_that("integration: initialize_genome -> add_founders -> mutate_table", {
  pop <- open_pop(pop_name = "test", db_name = ":memory:") |>
    define_genome(n_loci = 100, n_chr = 5, chr_len_Mb = 100) |>
    define_founder_haplotypes(n_haplotypes = 50) |>
    get_table("founder_haplotypes") |>
    add_founders(n_males = 10, n_females = 100, line_name = "A") %>%
    get_table("ind_meta") %>%
    mutate_table(
      gen = 0,
      farm = "FarmA",
      date_birth = as.Date("2024-01-01")
    )

  ind_meta <- get_table(pop, "ind_meta") %>% dplyr::collect()

  expect_equal(nrow(ind_meta), 110)
  expect_true(all(c("id_ind", "id_parent_1", "id_parent_2", "line_name", "sex",
                    "gen", "farm", "date_birth") %in% colnames(ind_meta)))

  expect_true(all(ind_meta$gen == 0))
  expect_true(all(ind_meta$farm == "FarmA"))
  expect_true(all(ind_meta$date_birth == as.Date("2024-01-01")))

  close_pop(pop)
})


test_that("haplotypes are sampled with replacement", {
  pop <- open_pop(pop_name = "test", db_name = ":memory:") |>
    define_genome(n_loci = 10, n_chr = 1, chr_len_Mb = 100) |>
    define_founder_haplotypes(n_haplotypes = 5)  # Small pool

  pop <- pop |>
    get_table("founder_haplotypes") |>
    add_founders(n_males = 20, n_females = 20, line_name = "A")

  ind_meta <- get_table(pop, "ind_meta") %>% dplyr::collect()
  expect_equal(nrow(ind_meta), 40)

  close_pop(pop)
})


test_that("add_founders handles large lines efficiently", {
  skip_on_cran()

  pop <- open_pop(pop_name = "test", db_name = ":memory:") |>
    define_genome(n_loci = 1000, n_chr = 10, chr_len_Mb = 100) |>
    define_founder_haplotypes(n_haplotypes = 100)

  expect_no_error({
    pop <- pop |>
      get_table("founder_haplotypes") |>
      add_founders(n_males = 500, n_females = 500, line_name = "A")
  })

  ind_meta <- get_table(pop, "ind_meta") %>% dplyr::collect()
  expect_equal(nrow(ind_meta), 1000)

  # 1000 individuals x 2 haplotypes x 1000 loci in the long table.
  n_rows <- DBI::dbGetQuery(pop$db_conn,
    "SELECT COUNT(*) AS n FROM ind_haplotype")$n
  expect_equal(n_rows, 1000 * 2 * 1000)
  n_ind <- DBI::dbGetQuery(pop$db_conn,
    "SELECT COUNT(DISTINCT id_ind) AS n FROM ind_haplotype")$n
  expect_equal(n_ind, 1000)

  close_pop(pop)
})


# ─── ... forwarding tests ────────────────────────────────────────────────────

make_pop_for_extra <- function() {
  open_pop(pop_name = "t", db_name = ":memory:") |>
    define_genome(n_loci = 20, n_chr = 1, chr_len_Mb = 10) |>
    define_founder_haplotypes(n_haplotypes = 20)
}

test_that("add_founders() accepts scalar ... fields and writes them to ind_meta", {
  pop <- make_pop_for_extra()
  pop <- pop |>
    get_table("founder_haplotypes") |>
    add_founders(n_males = 5, n_females = 5, line_name = "A",
                 gen = 0L, farm = "Iowa")

  result <- get_table(pop, "ind_meta") |> dplyr::collect()
  expect_equal(nrow(result), 10L)
  expect_true("gen"  %in% names(result))
  expect_true("farm" %in% names(result))
  expect_true(all(result$gen  == 0L))
  expect_true(all(result$farm == "Iowa"))

  close_pop(pop)
})

test_that("add_founders() ... uses correct DuckDB type from R type", {
  pop <- make_pop_for_extra()
  pop <- pop |>
    get_table("founder_haplotypes") |>
    add_founders(n_males = 3, n_females = 3, line_name = "A",
                 gen_int = 1L, gen_dbl = 1.0, flag = TRUE)

  col_types <- DBI::dbGetQuery(pop$db_conn,
    "SELECT column_name, data_type FROM information_schema.columns
     WHERE table_name = 'ind_meta'
     AND column_name IN ('gen_int', 'gen_dbl', 'flag')")
  types <- setNames(col_types$data_type, col_types$column_name)

  expect_equal(types[["gen_int"]], "INTEGER")
  expect_equal(types[["gen_dbl"]], "DOUBLE")
  expect_equal(types[["flag"]],    "BOOLEAN")

  close_pop(pop)
})

test_that("add_founders() ... accepts a vector of correct length", {
  pop <- make_pop_for_extra()
  scores <- c(1L, 2L, 3L, 4L, 5L)
  pop <- pop |>
    get_table("founder_haplotypes") |>
    add_founders(n_males = 5, n_females = 0, line_name = "A",
                 rank = scores)

  result <- get_table(pop, "ind_meta") |> dplyr::collect()
  expect_equal(result$rank, scores)

  close_pop(pop)
})

test_that("add_founders() ... errors for vector of wrong length", {
  pop <- make_pop_for_extra()
  expect_error(
    pop |>
      get_table("founder_haplotypes") |>
      add_founders(n_males = 5, n_females = 5, line_name = "A",
                   rank = c(1L, 2L, 3L)),
    "must equal the number of rows"
  )
  close_pop(pop)
})

test_that("add_founders() ... blocks reserved column names", {
  pop <- make_pop_for_extra()
  expect_error(
    pop |>
      get_table("founder_haplotypes") |>
      add_founders(n_males = 3, n_females = 3, line_name = "A", sex = "X"),
    "reserved"
  )
  close_pop(pop)
})

test_that("add_founders() ... preserves values across two cohort calls", {
  pop <- make_pop_for_extra()
  pop <- pop |>
    get_table("founder_haplotypes") |>
    add_founders(n_males = 3, n_females = 3, line_name = "A", gen = 0L)
  pop <- pop |>
    get_table("founder_haplotypes") |>
    add_founders(n_males = 2, n_females = 2, line_name = "B", gen = 1L)

  result <- get_table(pop, "ind_meta") |> dplyr::collect()
  expect_equal(sum(result$gen == 0L), 6L)
  expect_equal(sum(result$gen == 1L), 4L)

  close_pop(pop)
})


# ─── haplotype pool filter tests ──────────────────────────────────────────────

test_that("add_founders uses haplotypes from an explicitly filtered pool", {
  pop <- open_pop(pop_name = "lnpool", db_name = ":memory:") |>
    define_genome(n_loci = 20, n_chr = 1, chr_len_Mb = 50) |>
    define_founder_haplotypes(n_haplotypes = 30, line_name = "A") |>
    define_founder_haplotypes(n_haplotypes = 40, line_name = "B")

  pop <- pop |>
    get_table("founder_haplotypes") |>
    dplyr::filter(line_name == "A") |>
    add_founders(n_males = 5, n_females = 5, line_name = "A")

  pop <- pop |>
    get_table("founder_haplotypes") |>
    dplyr::filter(line_name == "B") |>
    add_founders(n_males = 5, n_females = 5, line_name = "B")

  ind_meta <- get_table(pop, "ind_meta") |> dplyr::collect()
  expect_equal(nrow(ind_meta), 20L)

  # 20 individuals in ind_haplotype; line_origin split A/B by line.
  n_ind <- DBI::dbGetQuery(pop$db_conn,
    "SELECT COUNT(DISTINCT id_ind) AS n FROM ind_haplotype")$n
  expect_equal(n_ind, 20L)
  lo <- DBI::dbGetQuery(pop$db_conn,
    "SELECT DISTINCT id_ind, line_origin FROM ind_haplotype")
  expect_setequal(unique(lo$line_origin), c("A", "B"))
  expect_true(all(lo$line_origin == sub("_[0-9]+$", "", lo$id_ind)))

  close_pop(pop)
})


test_that("add_founders uses all haplotypes when no filter is applied", {
  pop <- open_pop(pop_name = "lnnofilter", db_name = ":memory:") |>
    define_genome(n_loci = 20, n_chr = 1, chr_len_Mb = 50) |>
    define_founder_haplotypes(n_haplotypes = 25)

  pop <- pop |>
    get_table("founder_haplotypes") |>
    add_founders(n_males = 5, n_females = 5, line_name = "A")

  ind_meta <- get_table(pop, "ind_meta") |> dplyr::collect()
  expect_equal(nrow(ind_meta), 10L)
  expect_true(all(ind_meta$line_name == "A"))

  close_pop(pop)
})


test_that("add_founders errors when filter produces zero rows", {
  pop <- open_pop(pop_name = "lnerr", db_name = ":memory:") |>
    define_genome(n_loci = 10, n_chr = 1, chr_len_Mb = 50) |>
    define_founder_haplotypes(n_haplotypes = 20, line_name = "A")

  expect_error(
    pop |>
      get_table("founder_haplotypes") |>
      dplyr::filter(line_name == "B") |>
      add_founders(n_males = 5, n_females = 5, line_name = "B"),
    "No haplotypes selected"
  )

  close_pop(pop)
})


test_that("add_founders with named pool uses correct hap_id prefix in pool", {
  set.seed(42)
  pop <- open_pop(pop_name = "lnpfx", db_name = ":memory:") |>
    define_genome(n_loci = 10, n_chr = 1, chr_len_Mb = 50) |>
    define_founder_haplotypes(n_haplotypes = 20, line_name = "Duroc")

  fh <- get_table(pop, "founder_haplotypes") |> dplyr::collect()
  expect_true(all(grepl("^Duroc_hap_", fh$hap_id)))
  expect_true(all(fh$line_name == "Duroc"))

  pop <- pop |>
    get_table("founder_haplotypes") |>
    dplyr::filter(line_name == "Duroc") |>
    add_founders(n_males = 3, n_females = 3, line_name = "Duroc")

  ind_meta <- get_table(pop, "ind_meta") |> dplyr::collect()
  expect_equal(nrow(ind_meta), 6L)

  close_pop(pop)
})
