# SQL-injection hardening for add_dosage() / extract_genotypes(): chip_name
# identifier validation (via the constructed is_<chip_name>/has_<chip_name>
# column) and sql_in_list() escaping for locus_name/id_ind IN-lists.

# ---------------------------------------------------------------------------
# sql_in_list() direct unit tests
# ---------------------------------------------------------------------------

test_that("sql_in_list() escapes embedded single quotes and joins values", {
  expect_equal(sql_in_list(c("A", "B's", "C")), "'A', 'B''s', 'C'")
})

test_that("sql_in_list() errors on empty input", {
  expect_error(sql_in_list(character(0)), "No.*provided")
})

test_that("sql_in_list() errors on NA", {
  expect_error(sql_in_list(c("A", NA_character_)), "NA")
})

# ---------------------------------------------------------------------------
# chip_name identifier validation
# ---------------------------------------------------------------------------

test_that("add_dosage() rejects an invalid chip_name before touching data", {
  pop <- make_pop_base(n_loci = 10, n_chr = 2, chr_len_Mb = 100) |>
    define_founder_haplotypes(n_haplotypes = 10) |>
    get_table("founder_haplotypes") |>
    add_founders(n_males = 2, n_females = 2, line_name = "A")
  on.exit(close_pop(pop), add = TRUE)

  expect_error(
    pop |> get_table("ind_meta") |> add_dosage(chip_name = "bad name"),
    "chip name"
  )
  expect_error(
    pop |> get_table("ind_meta") |> add_dosage(chip_name = "bad'name"),
    "chip name"
  )
  # "50k"-style chip names (digit-leading, industry-standard) must still work
  # because the constructed identifier ("is_50k") starts with a letter.
  pop <- pop |>
    get_table("genome_meta") |>
    dplyr::filter(locus_id <= 3) |>
    define_chip("50k")
  expect_no_error(pop |> get_table("ind_meta") |> add_dosage(chip_name = "50k"))
})

test_that("extract_genotypes() rejects an invalid chip_name before touching data", {
  pop <- make_pop_base(n_loci = 10, n_chr = 2, chr_len_Mb = 100) |>
    define_founder_haplotypes(n_haplotypes = 10) |>
    get_table("founder_haplotypes") |>
    add_founders(n_males = 2, n_females = 2, line_name = "A")
  on.exit(close_pop(pop), add = TRUE)

  # With the default col_name (derived from chip_name as has_<chip_name>),
  # the pre-existing col_name check already rejects an unsafe chip_name.
  expect_error(
    pop |> get_table("ind_meta") |> extract_genotypes(chip_name = "bad name"),
    "Invalid"
  )

  # The gap this stage closes: a caller can supply their OWN safe col_name,
  # which passes the col_name check regardless of chip_name — chip_name
  # itself is still used to build chip_col (is_<chip_name>) for the
  # genome_meta lookup, and previously reached that bare-identifier SQL
  # context completely unvalidated.
  expect_error(
    pop |>
      get_table("ind_meta") |>
      extract_genotypes(chip_name = "bad'name", col_name = "has_safe"),
    "chip name"
  )
})

# ---------------------------------------------------------------------------
# Realistic end-to-end escaping check: a real locus_name containing a quote
# (e.g. from a careless import) must not corrupt the IN-list query.
# ---------------------------------------------------------------------------

test_that("add_dosage() locus_names containing a quote do not corrupt the query", {
  pop <- make_pop_base(n_loci = 10, n_chr = 2, chr_len_Mb = 100) |>
    define_founder_haplotypes(n_haplotypes = 10) |>
    get_table("founder_haplotypes") |>
    add_founders(n_males = 2, n_females = 2, line_name = "A")
  on.exit(close_pop(pop), add = TRUE)

  expect_no_error(
    pop |>
      get_table("ind_meta") |>
      add_dosage(locus_names = c("Locus_1", "nonexistent'; DROP TABLE ind_genotype; --"))
  )

  expect_true(DBI::dbExistsTable(pop$db_conn, "ind_genotype"))
  ln <- DBI::dbGetQuery(pop$db_conn,
    "SELECT DISTINCT locus_name FROM ind_genotype")$locus_name
  expect_equal(ln, "Locus_1")
})

test_that("extract_genotypes() with a real quote-embedded locus_name still returns correct results", {
  pop <- open_pop(pop_name = "inj_test", db_name = ":memory:") |>
    define_genome(n_loci = 10, n_chr = 2, chr_len_Mb = 100) |>
    define_founder_haplotypes(n_haplotypes = 10) |>
    get_table("founder_haplotypes") |>
    add_founders(n_males = 2, n_females = 2, line_name = "A") |>
    define_trait("ADG", target_add_mean = 100, target_add_var = 400)
  on.exit(close_pop(pop), add = TRUE)

  pop <- pop |>
    get_table("genome_meta") |>
    dplyr::filter(locus_id <= 3) |>
    define_additive_effects("ADG")

  # Rename one QTL locus to embed a single quote, keeping genome_meta,
  # ind_haplotype, and genome_effects consistent (locus_name is denormalized
  # across all three).
  for (tbl in c("genome_meta", "ind_haplotype", "genome_effects")) {
    DBI::dbExecute(pop$db_conn, paste0(
      "UPDATE ", tbl, " SET locus_name = 'Locus_1''s' WHERE locus_name = 'Locus_1'"
    ))
  }

  geno <- pop |>
    get_table("ind_meta") |>
    extract_genotypes(
      effects_tbl = get_table(pop, "genome_effects") |>
        dplyr::filter(trait_name == "ADG")
    )

  expect_equal(nrow(geno), 4L)
  locus_cols <- grep("^locus_", names(geno), value = TRUE)
  expect_equal(length(locus_cols), 3L)  # still all 3 QTL loci present
})
