# Helper: create a population with an ind_meta table
create_pop_for_mutate <- function(n_males = 5, n_females = 5) {
  pop <- initialize_genome(
    pop_name = "test",
    n_loci   = 20,
    n_chr    = 2,
    chr_len_Mb = 50,
    db_path  = ":memory:"
  )
  n_ind  <- n_males + n_females
  ind_df <- tibble::tibble(
    id_ind      = paste0("A-", seq_len(n_ind)),
    id_parent_1 = NA_character_,
    id_parent_2 = NA_character_,
    line        = "A",
    sex         = c(rep("M", n_males), rep("F", n_females))
  )
  DBI::dbWriteTable(pop$db_conn, "ind_meta", ind_df, overwrite = TRUE)
  pop$tables <- c(pop$tables, "ind_meta")
  pop
}


# ── 1. tidybreed_table class ──────────────────────────────────────────────────

test_that("get_table() returns a tidybreed_table", {
  pop <- create_pop_for_mutate()
  tbl <- get_table(pop, "ind_meta")
  expect_s3_class(tbl, "tidybreed_table")
  close_pop(pop)
})


test_that("get_table() errors with clear message for non-existent table", {
  pop <- create_pop_for_mutate()
  expect_error(get_table(pop, "no_such_table"), "does not exist in this population")
  close_pop(pop)
})


test_that("collect(get_table()) returns a tibble (backward compat)", {
  pop <- create_pop_for_mutate()
  result <- dplyr::collect(get_table(pop, "ind_meta"))
  expect_s3_class(result, "data.frame")
  expect_equal(nrow(result), 10)
  close_pop(pop)
})


test_that("get_table() |> filter() |> collect() works (backward compat)", {
  pop <- create_pop_for_mutate()
  males <- get_table(pop, "ind_meta") |>
    dplyr::filter(sex == "M") |>
    dplyr::collect()
  expect_equal(nrow(males), 5)
  close_pop(pop)
})


test_that("get_table() |> pull() works", {
  pop <- create_pop_for_mutate()
  ids <- get_table(pop, "ind_meta") |> dplyr::pull(id_ind)
  expect_length(ids, 10)
  close_pop(pop)
})


test_that("get_table() |> filter() |> pull() works", {
  pop <- create_pop_for_mutate()
  male_ids <- get_table(pop, "ind_meta") |>
    dplyr::filter(sex == "M") |>
    dplyr::pull(id_ind)
  expect_length(male_ids, 5)
  close_pop(pop)
})


test_that("get_table() |> filter() |> count() |> pull(n) works", {
  pop <- create_pop_for_mutate()
  n <- get_table(pop, "ind_meta") |>
    dplyr::filter(sex == "F") |>
    dplyr::count() |>
    dplyr::pull(n)
  expect_equal(n, 5)
  close_pop(pop)
})


# ── 2. mutate_table() — no filter ────────────────────────────────────────────

test_that("mutate_table() adds a scalar INTEGER column to ind_meta", {
  pop <- create_pop_for_mutate()
  pop <- get_table(pop, "ind_meta") |> mutate_table(gen = 1L)
  result <- dplyr::collect(get_table(pop, "ind_meta"))
  expect_true("gen" %in% colnames(result))
  expect_true(all(result$gen == 1L))
  close_pop(pop)
})


test_that("mutate_table() adds a scalar column to genome_meta", {
  pop <- initialize_genome(pop_name = "t", n_loci = 10, n_chr = 1,
                           chr_len_Mb = 10, db_path = ":memory:")
  pop <- get_table(pop, "genome_meta") |> mutate_table(maf = 0.3)
  result <- dplyr::collect(get_table(pop, "genome_meta"))
  expect_true("maf" %in% colnames(result))
  expect_true(all(abs(result$maf - 0.3) < 1e-10))
  close_pop(pop)
})


test_that("mutate_table() adds multiple scalar columns at once", {
  pop <- create_pop_for_mutate()
  pop <- get_table(pop, "ind_meta") |> mutate_table(gen = 0L, farm = "A")
  result <- dplyr::collect(get_table(pop, "ind_meta"))
  expect_true(all(c("gen", "farm") %in% colnames(result)))
  expect_true(all(result$gen == 0L))
  expect_true(all(result$farm == "A"))
  close_pop(pop)
})


test_that("mutate_table() adds a vector column with correct per-row values", {
  pop <- create_pop_for_mutate(n_males = 5, n_females = 5)
  weights <- seq(50, 59, by = 1)
  pop <- get_table(pop, "ind_meta") |> mutate_table(weight = weights)
  result <- dplyr::collect(get_table(pop, "ind_meta"))
  expect_equal(sort(result$weight), sort(weights))
  close_pop(pop)
})


test_that("mutate_table() updates an existing column", {
  pop <- create_pop_for_mutate()
  pop <- get_table(pop, "ind_meta") |> mutate_table(gen = 0L)
  expect_warning(
    { pop <- get_table(pop, "ind_meta") |> mutate_table(gen = 1L) },
    "replaced"
  )
  result <- dplyr::collect(get_table(pop, "ind_meta"))
  expect_true(all(result$gen == 1L))
  close_pop(pop)
})


test_that("mutate_table() returns pop invisibly", {
  pop <- create_pop_for_mutate()
  result <- get_table(pop, "ind_meta") |> mutate_table(gen = 0L)
  expect_s3_class(result, "tidybreed_pop")
  close_pop(pop)
})


test_that("mutate_table() creates column schema on empty table without warning", {
  pop <- initialize_genome(pop_name = "t", n_loci = 5, n_chr = 1,
                           chr_len_Mb = 10, db_path = ":memory:")
  # ind_meta is empty immediately after initialize_genome()
  expect_message(
    get_table(pop, "ind_meta") |> mutate_table(gen = NA_integer_),
    "empty"
  )
  cols <- DBI::dbListFields(pop$db_conn, "ind_meta")
  expect_true("gen" %in% cols)
  # Column type should be INTEGER (DuckDB stores as INTEGER)
  col_types <- DBI::dbGetQuery(pop$db_conn,
    "SELECT data_type FROM information_schema.columns
     WHERE table_name = 'ind_meta' AND column_name = 'gen'")
  expect_equal(col_types$data_type, "INTEGER")
  close_pop(pop)
})


test_that("mutate_table() warns and returns early when no fields given", {
  pop <- create_pop_for_mutate()
  expect_warning(
    get_table(pop, "ind_meta") |> mutate_table(),
    "No fields specified"
  )
  close_pop(pop)
})


# ── 3. mutate_table() — with filter ──────────────────────────────────────────

test_that("new column + filter: matched rows get value, others get NULL", {
  pop <- create_pop_for_mutate(n_males = 5, n_females = 5)
  pop <- get_table(pop, "ind_meta") |>
    dplyr::filter(sex == "M") |>
    mutate_table(gen = 1L)
  result <- dplyr::collect(get_table(pop, "ind_meta"))
  expect_true("gen" %in% colnames(result))
  male_gen   <- result$gen[result$sex == "M"]
  female_gen <- result$gen[result$sex == "F"]
  expect_true(all(male_gen == 1L))
  expect_true(all(is.na(female_gen)))
  close_pop(pop)
})


test_that("existing column + filter: only matched rows updated", {
  pop <- create_pop_for_mutate(n_males = 5, n_females = 5)
  pop <- get_table(pop, "ind_meta") |> mutate_table(gen = 0L)
  expect_warning(
    {
      pop <- get_table(pop, "ind_meta") |>
        dplyr::filter(sex == "M") |>
        mutate_table(gen = 2L)
    },
    "replaced"
  )
  result <- dplyr::collect(get_table(pop, "ind_meta"))
  expect_true(all(result$gen[result$sex == "M"] == 2L))
  expect_true(all(result$gen[result$sex == "F"] == 0L))
  close_pop(pop)
})


test_that("vector + filter: vector of length n_filtered succeeds", {
  pop <- create_pop_for_mutate(n_males = 5, n_females = 5)
  scores <- c(10.1, 10.2, 10.3, 10.4, 10.5)
  pop <- get_table(pop, "ind_meta") |>
    dplyr::filter(sex == "M") |>
    mutate_table(score = scores)
  result <- dplyr::collect(get_table(pop, "ind_meta"))
  expect_equal(sort(result$score[result$sex == "M"]), sort(scores))
  expect_true(all(is.na(result$score[result$sex == "F"])))
  close_pop(pop)
})


test_that("vector + filter: wrong length errors with clear message", {
  pop <- create_pop_for_mutate(n_males = 5, n_females = 5)
  expect_error(
    get_table(pop, "ind_meta") |>
      dplyr::filter(sex == "M") |>
      mutate_table(score = 1:4),  # 4, not 5
    "5"
  )
  close_pop(pop)
})


# ── 4. Error paths ────────────────────────────────────────────────────────────

test_that("mutate_table() errors if not called after get_table()", {
  expect_error(mutate_table(list(a = 1), gen = 1L), "must be called after get_table")
})


test_that("reserved columns are blocked in ind_meta", {
  pop <- create_pop_for_mutate()
  for (col in c("id_ind", "id_parent_1", "id_parent_2", "line", "sex")) {
    args <- setNames(list(1L), col)
    expect_error(
      do.call(mutate_table, c(list(tbl_obj = get_table(pop, "ind_meta")), args)),
      "reserved"
    )
  }
  close_pop(pop)
})


test_that("reserved columns are blocked in genome_meta", {
  pop <- initialize_genome(pop_name = "t", n_loci = 5, n_chr = 1,
                           chr_len_Mb = 10, db_path = ":memory:")
  for (col in c("locus_id", "locus_name", "chr", "chr_name", "pos_Mb")) {
    args <- setNames(list(1L), col)
    expect_error(
      do.call(mutate_table, c(list(tbl_obj = get_table(pop, "genome_meta")), args)),
      "reserved"
    )
  }
  close_pop(pop)
})


test_that("invalid SQL identifier errors", {
  pop <- create_pop_for_mutate()
  expect_error(
    get_table(pop, "ind_meta") |> mutate_table(`1bad` = 1L),
    "Invalid"
  )
  close_pop(pop)
})


test_that("SQL reserved keywords are blocked", {
  pop <- create_pop_for_mutate()
  expect_error(
    get_table(pop, "ind_meta") |> mutate_table(SELECT = 1L),
    "reserved keyword"
  )
  close_pop(pop)
})


test_that("vector length mismatch errors with clear message", {
  pop <- create_pop_for_mutate(n_males = 5, n_females = 5)
  expect_error(
    get_table(pop, "ind_meta") |> mutate_table(gen = 1:9),   # 9, not 10
    "10"
  )
  close_pop(pop)
})


# ── 5. Type handling ──────────────────────────────────────────────────────────

test_that("logical value maps to BOOLEAN", {
  pop <- create_pop_for_mutate()
  pop <- get_table(pop, "ind_meta") |> mutate_table(is_sel = TRUE)
  result <- dplyr::collect(get_table(pop, "ind_meta"))
  expect_type(result$is_sel, "logical")
  close_pop(pop)
})


test_that("integer value maps to INTEGER", {
  pop <- create_pop_for_mutate()
  pop <- get_table(pop, "ind_meta") |> mutate_table(gen = 1L)
  result <- dplyr::collect(get_table(pop, "ind_meta"))
  expect_type(result$gen, "integer")
  close_pop(pop)
})


test_that("numeric value maps to DOUBLE", {
  pop <- create_pop_for_mutate()
  pop <- get_table(pop, "ind_meta") |> mutate_table(weight = 75.5)
  result <- dplyr::collect(get_table(pop, "ind_meta"))
  expect_type(result$weight, "double")
  close_pop(pop)
})


test_that("Date value maps to DATE", {
  pop <- create_pop_for_mutate()
  pop <- get_table(pop, "ind_meta") |> mutate_table(dob = as.Date("2024-01-15"))
  result <- dplyr::collect(get_table(pop, "ind_meta"))
  expect_s3_class(result$dob, "Date")
  close_pop(pop)
})


test_that("character value maps to VARCHAR", {
  pop <- create_pop_for_mutate()
  pop <- get_table(pop, "ind_meta") |> mutate_table(farm = "Iowa")
  result <- dplyr::collect(get_table(pop, "ind_meta"))
  expect_type(result$farm, "character")
  expect_true(all(result$farm == "Iowa"))
  close_pop(pop)
})
