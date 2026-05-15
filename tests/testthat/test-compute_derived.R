# ── Helpers ───────────────────────────────────────────────────────────────────

# Minimal pop with 5 males + 10 females, birth_date in ind_meta,
# a count trait "AP" phenotyped for all females (pheno_number = 1).
make_ap_pop <- function() {
  pop <- initialize_genome(
    pop_name   = "cd_test",
    n_loci     = 50,
    n_chr      = 2,
    chr_len_Mb = 50,
    n_haplotypes = 50,
    db_path    = ":memory:"
  )
  pop <- add_founders(pop, n_males = 5, n_females = 10, line_name = "L")
  pop <- get_table(pop, "ind_meta") |>
    mutate_table(birth_date = as.Date("2024-01-01"))
  pop <- pop |>
    add_trait(
      "AP",
      trait_type     = "count",
      target_add_var = 1,
      residual_var   = 1,
      target_add_mean = 180,
      min_value      = 150,
      max_value      = 250,
      expressed_sex  = "F"
    )
  sel <- pop |> get_table("genome_meta") |> dplyr::collect() |>
    dplyr::slice_sample(n = 10) |> dplyr::pull(locus_name)
  pop <- pop |>
    get_table("genome_meta") |>
    dplyr::filter(locus_name %in% sel) |>
    define_qtl("AP") |>
    set_qtl_effects("AP")
  pop <- get_table(pop, "ind_meta") |>
    dplyr::filter(sex == "F") |>
    add_phenotype("AP")
  pop
}


# ── 1. Core functionality ─────────────────────────────────────────────────────

test_that("compute_derived() returns pop invisibly", {
  # run all code above as function
  pop <- make_ap_pop()
  
  # use `compute_derived()` function to add `puberty_date` to `ind_meta`
  result <- get_table(pop, "ind_phenotype") |>
    dplyr::filter(trait_name == "AP", pheno_number == 1L) |>
    compute_derived(
      compute    = \(df) df$birth_date + as.integer(df$pheno_value),
      join_table = "ind_meta",
      join_by    = "id_ind",
      write_to   = c(ind_meta = "puberty_date")
    )
  
  expect_s3_class(result, "tidybreed_pop")
  close_pop(pop)
})


test_that("compute_derived() writes to secondary table (ind_meta)", {
  pop <- make_ap_pop()
  pop <- get_table(pop, "ind_phenotype") |>
    dplyr::filter(trait_name == "AP", pheno_number == 1L) |>
    compute_derived(
      compute    = \(df) df$birth_date + as.integer(df$pheno_value),
      join_table = "ind_meta",
      join_by    = "id_ind",
      write_to   = c(ind_meta = "puberty_date")
    )

  meta_f <- get_table(pop, "ind_meta") |>
    dplyr::filter(sex == "F") |>
    dplyr::collect()
  ap_pheno <- get_table(pop, "ind_phenotype") |>
    dplyr::filter(trait_name == "AP", pheno_number == 1L) |>
    dplyr::collect()

  # All females should have puberty_date
  expect_false(any(is.na(meta_f$puberty_date)))
  expect_s3_class(meta_f$puberty_date, "Date")

  # Value should equal birth_date + AP value (in days)
  joined <- merge(
    meta_f[, c("id_ind", "birth_date", "puberty_date")],
    ap_pheno[, c("id_ind", "pheno_value")],
    by = "id_ind"
  )
  expected <- as.Date("2024-01-01") + as.integer(joined$pheno_value)
  expect_equal(joined$puberty_date, expected)

  close_pop(pop)
})


test_that("compute_derived() writes to primary table (ind_phenotype)", {
  pop <- make_ap_pop()
  pop <- get_table(pop, "ind_phenotype") |>
    dplyr::filter(trait_name == "AP", pheno_number == 1L) |>
    compute_derived(
      compute    = \(df) df$birth_date + as.integer(df$pheno_value),
      join_table = "ind_meta",
      join_by    = "id_ind",
      write_to   = c(ind_phenotype = "pheno_date")
    )

  ap_pheno <- get_table(pop, "ind_phenotype") |>
    dplyr::filter(trait_name == "AP", pheno_number == 1L) |>
    dplyr::collect()

  expect_true("pheno_date" %in% names(ap_pheno))
  expect_false(any(is.na(ap_pheno$pheno_date)))
  expect_s3_class(ap_pheno$pheno_date, "Date")

  # Verify value
  meta_f <- get_table(pop, "ind_meta") |>
    dplyr::filter(sex == "F") |>
    dplyr::collect()
  joined <- merge(ap_pheno[, c("id_ind", "pheno_value", "pheno_date")],
                  meta_f[, c("id_ind", "birth_date")], by = "id_ind")
  expect_equal(joined$pheno_date, as.Date("2024-01-01") + as.integer(joined$pheno_value))

  close_pop(pop)
})


test_that("compute_derived() writes to both tables in a single call", {
  pop <- make_ap_pop()
  pop <- get_table(pop, "ind_phenotype") |>
    dplyr::filter(trait_name == "AP", pheno_number == 1L) |>
    compute_derived(
      compute    = \(df) df$birth_date + as.integer(df$pheno_value),
      join_table = "ind_meta",
      join_by    = "id_ind",
      write_to   = c(ind_meta = "puberty_date", ind_phenotype = "pheno_date")
    )

  meta_f   <- get_table(pop, "ind_meta") |>
    dplyr::filter(sex == "F") |> dplyr::collect()
  ap_pheno <- get_table(pop, "ind_phenotype") |>
    dplyr::filter(trait_name == "AP", pheno_number == 1L) |>
    dplyr::collect()

  expect_true("puberty_date" %in% names(meta_f))
  expect_true("pheno_date" %in% names(ap_pheno))
  expect_false(any(is.na(meta_f$puberty_date)))
  expect_false(any(is.na(ap_pheno$pheno_date)))

  # puberty_date in ind_meta and pheno_date in ind_phenotype should match
  joined <- merge(
    meta_f[, c("id_ind", "puberty_date")],
    ap_pheno[, c("id_ind", "pheno_date")],
    by = "id_ind"
  )
  expect_equal(joined$puberty_date, joined$pheno_date)

  close_pop(pop)
})


test_that("compute_derived() updates an existing column (no duplicate column added)", {
  pop <- make_ap_pop()
  # First call creates puberty_date
  pop <- get_table(pop, "ind_phenotype") |>
    dplyr::filter(trait_name == "AP", pheno_number == 1L) |>
    compute_derived(
      compute    = \(df) df$birth_date + as.integer(df$pheno_value),
      join_table = "ind_meta",
      join_by    = "id_ind",
      write_to   = c(ind_meta = "puberty_date")
    )
  first_vals <- get_table(pop, "ind_meta") |>
    dplyr::filter(sex == "F") |> dplyr::collect() |>
    dplyr::pull(puberty_date)

  # Second call overwrites with a different constant
  pop <- get_table(pop, "ind_meta") |>
    compute_derived(
      compute  = \(df) as.Date("2025-06-01"),
      write_to = c(ind_meta = "puberty_date")
    )
  second_vals <- get_table(pop, "ind_meta") |>
    dplyr::collect() |> dplyr::pull(puberty_date)

  # Column still exists (not duplicated)
  expect_equal(sum(names(dplyr::collect(get_table(pop, "ind_meta"))) == "puberty_date"), 1L)
  # All rows updated to the new value
  expect_true(all(second_vals == as.Date("2025-06-01"), na.rm = TRUE))

  close_pop(pop)
})


test_that("NAs in result_vec are written as NULL (remain NA in R)", {
  pop <- make_ap_pop()
  # compute returns NA for everyone
  pop <- get_table(pop, "ind_meta") |>
    compute_derived(
      compute  = \(df) rep(NA_real_, nrow(df)),
      write_to = c(ind_meta = "some_score")
    )
  scores <- get_table(pop, "ind_meta") |> dplyr::collect() |>
    dplyr::pull(some_score)
  expect_true(all(is.na(scores)))
  close_pop(pop)
})


test_that("compute_derived() works without join_table (primary table only)", {
  pop <- make_ap_pop()
  pop <- get_table(pop, "ind_phenotype") |>
    dplyr::filter(trait_name == "AP") |>
    compute_derived(
      compute  = \(df) df$pheno_value * 0.1,   # scaled value
      write_to = c(ind_phenotype = "value_scaled")
    )
  ap <- get_table(pop, "ind_phenotype") |>
    dplyr::filter(trait_name == "AP") |>
    dplyr::collect()
  expect_true("value_scaled" %in% names(ap))
  expect_equal(ap$value_scaled, ap$pheno_value * 0.1)
  close_pop(pop)
})


# ── 2. Deduplication / conflict ───────────────────────────────────────────────

test_that("conflict error when two rows for same id_ind produce different values", {
  pop <- make_ap_pop()
  # Add a second repeatable AP record by temporarily making it repeatable;
  # we simulate this by directly inserting a duplicate ind_phenotype row
  ap_rows <- get_table(pop, "ind_phenotype") |>
    dplyr::filter(trait_name == "AP") |>
    dplyr::collect()
  # Create duplicate rows with a different value
  dup_rows <- ap_rows
  dup_rows$id_phenotype <- dup_rows$id_phenotype + max(ap_rows$id_phenotype)
  dup_rows$pheno_number <- 2L
  dup_rows$pheno_value  <- dup_rows$pheno_value + 10  # different value
  DBI::dbWriteTable(pop$db_conn, "ind_phenotype", dup_rows, append = TRUE)

  expect_error(
    get_table(pop, "ind_phenotype") |>
      dplyr::filter(trait_name == "AP") |>
      compute_derived(
        compute    = \(df) df$birth_date + as.integer(df$pheno_value),
        join_table = "ind_meta",
        join_by    = "id_ind",
        write_to   = c(ind_meta = "puberty_date")
      ),
    "maps to different computed results"
  )
  close_pop(pop)
})


test_that("no conflict when two rows produce the same value (silently deduped)", {
  pop <- make_ap_pop()
  # Add duplicate rows with the SAME value
  ap_rows <- get_table(pop, "ind_phenotype") |>
    dplyr::filter(trait_name == "AP") |>
    dplyr::collect()
  dup_rows <- ap_rows
  dup_rows$id_phenotype <- dup_rows$id_phenotype + max(ap_rows$id_phenotype)
  dup_rows$pheno_number <- 2L
  # value stays the same — constant compute = \(df) "same_label"
  DBI::dbWriteTable(pop$db_conn, "ind_phenotype", dup_rows, append = TRUE)

  expect_no_error(
    get_table(pop, "ind_phenotype") |>
      dplyr::filter(trait_name == "AP") |>
      compute_derived(
        compute  = \(df) "constant",   # same for every row
        write_to = c(ind_meta = "label_col")
      )
  )
  close_pop(pop)
})


# ── 3. Input validation ───────────────────────────────────────────────────────

test_that("error if tbl_obj is not a tidybreed_table", {
  pop <- make_ap_pop()
  expect_error(
    compute_derived(pop, compute = \(df) 1, write_to = c(ind_meta = "x")),
    "compute_derived\\(\\) must be called after get_table"
  )
  close_pop(pop)
})


test_that("error if compute is not a function", {
  pop <- make_ap_pop()
  expect_error(
    get_table(pop, "ind_meta") |>
      compute_derived(compute = "not_a_function", write_to = c(ind_meta = "x")),
    "'compute' must be a function"
  )
  close_pop(pop)
})


test_that("error if write_to is NULL", {
  pop <- make_ap_pop()
  expect_error(
    get_table(pop, "ind_meta") |>
      compute_derived(compute = \(df) 1, write_to = NULL),
    "'write_to' must be a named character vector"
  )
  close_pop(pop)
})


test_that("error if write_to has no names", {
  pop <- make_ap_pop()
  expect_error(
    get_table(pop, "ind_meta") |>
      compute_derived(compute = \(df) 1, write_to = "puberty_date"),
    "'write_to' must be a named character vector"
  )
  close_pop(pop)
})


test_that("error if compute() returns wrong-length vector (not 1 and not nrow)", {
  pop <- make_ap_pop()
  # ind_meta has 15 rows; returning length 3 is neither 1 nor 15
  expect_error(
    get_table(pop, "ind_meta") |>
      compute_derived(compute = \(df) c(1, 2, 3), write_to = c(ind_meta = "x")),
    "compute\\(\\) returned .* values but the joined data has"
  )
  close_pop(pop)
})


test_that("error if join_table not in pop$tables", {
  pop <- make_ap_pop()
  expect_error(
    get_table(pop, "ind_phenotype") |>
      compute_derived(
        compute    = \(df) 1,
        join_table = "no_such_table",
        write_to   = c(ind_phenotype = "x")
      ),
    "does not exist in this population"
  )
  close_pop(pop)
})


test_that("error if join_by not in primary table", {
  pop <- make_ap_pop()
  expect_error(
    get_table(pop, "ind_phenotype") |>
      compute_derived(
        compute    = \(df) 1,
        join_table = "ind_meta",
        join_by    = "no_such_col",
        write_to   = c(ind_phenotype = "x")
      ),
    "not found in primary table"
  )
  close_pop(pop)
})


test_that("error if join_by not in join_table", {
  pop <- make_ap_pop()
  expect_error(
    get_table(pop, "ind_meta") |>
      compute_derived(
        compute    = \(df) 1,
        join_table = "genome_meta",
        join_by    = "id_ind",
        write_to   = c(ind_meta = "x")
      ),
    "not found in join_table"
  )
  close_pop(pop)
})


test_that("error if destination table not in TABLE_PRIMARY_KEYS", {
  pop <- make_ap_pop()
  expect_error(
    get_table(pop, "ind_meta") |>
      compute_derived(
        compute  = \(df) 1,
        write_to = c(genome_haplotype = "x")
      ),
    "not registered in the primary-key map"
  )
  close_pop(pop)
})


test_that("error if dest_col is a reserved column", {
  pop <- make_ap_pop()
  expect_error(
    get_table(pop, "ind_phenotype") |>
      dplyr::filter(trait_name == "AP") |>
      compute_derived(
        compute  = \(df) df$pheno_value,
        write_to = c(ind_phenotype = "pheno_value")   # "pheno_value" is reserved
      ),
    "Cannot modify reserved column"
  )
  close_pop(pop)
})


test_that("early return with message when filter matches 0 rows", {
  pop <- make_ap_pop()
  expect_message(
    pop2 <- get_table(pop, "ind_phenotype") |>
      dplyr::filter(trait_name == "NONEXISTENT") |>
      compute_derived(
        compute  = \(df) df$pheno_value,
        write_to = c(ind_phenotype = "value_copy")
      ),
    "0 rows"
  )
  expect_s3_class(pop2, "tidybreed_pop")
  close_pop(pop)
})


# ── 4. Type handling ──────────────────────────────────────────────────────────

test_that("Date result creates DATE column", {
  pop <- make_ap_pop()
  pop <- get_table(pop, "ind_phenotype") |>
    dplyr::filter(trait_name == "AP", pheno_number == 1L) |>
    compute_derived(
      compute    = \(df) df$birth_date + as.integer(df$pheno_value),
      join_table = "ind_meta",
      join_by    = "id_ind",
      write_to   = c(ind_meta = "puberty_date")
    )
  result <- get_table(pop, "ind_meta") |>
    dplyr::filter(sex == "F") |> dplyr::collect()
  expect_s3_class(result$puberty_date, "Date")
  close_pop(pop)
})


test_that("integer result creates INTEGER column", {
  pop <- make_ap_pop()
  pop <- get_table(pop, "ind_meta") |>
    compute_derived(
      compute  = \(df) 1L,
      write_to = c(ind_meta = "gen")
    )
  result <- get_table(pop, "ind_meta") |> dplyr::collect()
  expect_type(result$gen, "integer")
  close_pop(pop)
})


test_that("double result creates DOUBLE column", {
  pop <- make_ap_pop()
  pop <- get_table(pop, "ind_meta") |>
    compute_derived(
      compute  = \(df) 1.5,
      write_to = c(ind_meta = "score_d")
    )
  result <- get_table(pop, "ind_meta") |> dplyr::collect()
  expect_type(result$score_d, "double")
  close_pop(pop)
})


test_that("character result creates VARCHAR column", {
  pop <- make_ap_pop()
  pop <- get_table(pop, "ind_meta") |>
    compute_derived(
      compute  = \(df) "group_A",
      write_to = c(ind_meta = "grp")
    )
  result <- get_table(pop, "ind_meta") |> dplyr::collect()
  expect_type(result$grp, "character")
  close_pop(pop)
})


# ── 5. Column name collision on join ─────────────────────────────────────────

test_that("primary table columns take precedence when join_table has same column", {
  pop <- make_ap_pop()
  # Both ind_phenotype and ind_meta have a 'line'-like column — add 'trait_name'
  # to ind_meta to create an artificial collision
  pop <- get_table(pop, "ind_meta") |>
    mutate_table(trait_name = "meta_side")

  pop <- get_table(pop, "ind_phenotype") |>
    dplyr::filter(pheno_number == 1L) |>
    compute_derived(
      compute    = \(df) df$trait_name,   # should be ind_phenotype's trait_name
      join_table = "ind_meta",
      join_by    = "id_ind",
      write_to   = c(ind_phenotype = "trait_name_check")
    )

  result <- get_table(pop, "ind_phenotype") |>
    dplyr::filter(pheno_number == 1L) |>
    dplyr::collect()
  # Should be the phenotype trait_name ("AP"), not "meta_side"
  expect_true(all(result$trait_name_check == "AP"))
  close_pop(pop)
})
