# Helper: tiny pop with two traits and fake EBVs for 10 individuals
make_index_pop <- function(pop_name = "idx") {
  pop <- initialize_genome(
    pop_name          = pop_name,
    n_loci            = 100,
    n_chr             = 2,
    chr_len_Mb        = 100,
    n_haplotypes      = 50,
    db_path           = ":memory:",
    fixed_allele_freq = 0.5
  )
  pop <- add_founders(pop, n_males = 5, n_females = 5, line_name = "A")

  # Seed EBVs directly (bypasses BLUPF90)
  ind_ids <- DBI::dbGetQuery(pop$db_conn,
    "SELECT id_ind FROM ind_meta")$id_ind

  ebv_adg <- data.frame(
    id_ind      = ind_ids,
    trait_name  = "ADG",
    model       = "test_model",
    ebv         = seq(100, 190, length.out = length(ind_ids)),
    acc         = NA_real_,
    se          = NA_real_,
    eval_number = 1L,
    stringsAsFactors = FALSE
  )
  ebv_fcr <- data.frame(
    id_ind      = ind_ids,
    trait_name  = "FCR",
    model       = "test_model",
    ebv         = seq(2.5, 3.4, length.out = length(ind_ids)),
    acc         = NA_real_,
    se          = NA_real_,
    eval_number = 1L,
    stringsAsFactors = FALSE
  )
  DBI::dbWriteTable(pop$db_conn, "ind_ebv",
                    rbind(ebv_adg, ebv_fcr), append = TRUE)
  pop
}


# ============================================================
# define_index() tests
# ============================================================

test_that("define_index() creates index_meta rows", {
  pop <- make_index_pop("di_basic")
  on.exit(close_pop(pop))

  pop <- define_index(pop, "terminal",
                      trait_names = c("ADG", "FCR"),
                      index_wts   = c(1.2, -0.8))

  rows <- DBI::dbGetQuery(pop$db_conn,
    "SELECT * FROM index_meta ORDER BY trait_name")
  expect_equal(nrow(rows), 2L)
  expect_equal(rows$index_name, c("terminal", "terminal"))
  expect_equal(rows$trait_name, c("ADG", "FCR"))
  expect_equal(rows$index_wt,   c(1.2, -0.8), tolerance = 1e-9)
})


test_that("define_index() upserts — re-running updates weights", {
  pop <- make_index_pop("di_upsert")
  on.exit(close_pop(pop))

  pop <- define_index(pop, "terminal",
                      trait_names = c("ADG", "FCR"),
                      index_wts   = c(1.2, -0.8))
  pop <- define_index(pop, "terminal",
                      trait_names = c("ADG", "FCR"),
                      index_wts   = c(2.0, -1.0))

  rows <- DBI::dbGetQuery(pop$db_conn,
    "SELECT * FROM index_meta WHERE index_name = 'terminal' ORDER BY trait_name")
  expect_equal(nrow(rows), 2L)
  expect_equal(rows$index_wt, c(2.0, -1.0), tolerance = 1e-9)
})


test_that("define_index() supports multiple named indexes", {
  pop <- make_index_pop("di_multi")
  on.exit(close_pop(pop))

  pop <- define_index(pop, "terminal",  trait_names = "ADG",  index_wts = 1.0)
  pop <- define_index(pop, "maternal",  trait_names = "FCR",  index_wts = -0.5)

  n <- DBI::dbGetQuery(pop$db_conn,
    "SELECT COUNT(*) AS n FROM index_meta")$n
  expect_equal(n, 2L)

  idx_names <- DBI::dbGetQuery(pop$db_conn,
    "SELECT DISTINCT index_name FROM index_meta ORDER BY index_name")$index_name
  expect_equal(idx_names, c("maternal", "terminal"))
})


test_that("define_index() supports extra user columns (scalar broadcast)", {
  pop <- make_index_pop("di_extra")
  on.exit(close_pop(pop))

  pop <- define_index(pop, "terminal",
                      trait_names = c("ADG", "FCR"),
                      index_wts   = c(1.2, -0.8),
                      source      = "gen_team_v1")

  cols <- DBI::dbListFields(pop$db_conn, "index_meta")
  expect_true("source" %in% cols)

  vals <- DBI::dbGetQuery(pop$db_conn,
    "SELECT source FROM index_meta ORDER BY trait_name")$source
  expect_equal(vals, c("gen_team_v1", "gen_team_v1"))
})


test_that("define_index() errors on invalid index_name", {
  pop <- make_index_pop("di_err_name")
  on.exit(close_pop(pop))

  expect_error(
    define_index(pop, "bad name", trait_names = "ADG", index_wts = 1.0),
    regexp = "Must start with a letter"
  )
  expect_error(
    define_index(pop, "SELECT", trait_names = "ADG", index_wts = 1.0),
    regexp = "SQL reserved keyword"
  )
})


test_that("define_index() errors when index_wts length mismatches trait_names", {
  pop <- make_index_pop("di_err_len")
  on.exit(close_pop(pop))

  expect_error(
    define_index(pop, "terminal",
                 trait_names = c("ADG", "FCR"),
                 index_wts   = 1.0),
    regexp = "index_wts"
  )
})


# ============================================================
# add_index() tests
# ============================================================

test_that("add_index() computes correct index values", {
  pop <- make_index_pop("ai_basic")
  on.exit(close_pop(pop))

  pop <- define_index(pop, "terminal",
                      trait_names = c("ADG", "FCR"),
                      index_wts   = c(1.0, -1.0))

  expect_warning(
    pop <- pop |>
      get_table("ind_ebv") |>
      add_index("terminal"),
    regexp = "No filter applied"
  )

  result <- DBI::dbGetQuery(pop$db_conn,
    "SELECT * FROM ind_index ORDER BY id_ind")
  expect_equal(nrow(result), 10L)
  expect_equal(result$index_name, rep("terminal", 10L))
  expect_equal(result$index_number, rep(1L, 10L))

  # Verify one computed value manually
  adg_vals <- DBI::dbGetQuery(pop$db_conn,
    "SELECT id_ind, ebv FROM ind_ebv WHERE trait_name = 'ADG' ORDER BY id_ind")
  fcr_vals <- DBI::dbGetQuery(pop$db_conn,
    "SELECT id_ind, ebv FROM ind_ebv WHERE trait_name = 'FCR' ORDER BY id_ind")
  expected <- adg_vals$ebv * 1.0 + fcr_vals$ebv * (-1.0)
  expect_equal(result$index_value, expected, tolerance = 1e-9)
})


test_that("add_index() issues warning when no filter applied", {
  pop <- make_index_pop("ai_warn")
  on.exit(close_pop(pop))

  pop <- define_index(pop, "terminal",
                      trait_names = "ADG", index_wts = 1.0)

  expect_warning(
    pop |> get_table("ind_ebv") |> add_index("terminal"),
    regexp = "No filter applied"
  )
})


test_that("add_index() does NOT warn when filter is applied", {
  pop <- make_index_pop("ai_no_warn")
  on.exit(close_pop(pop))

  pop <- define_index(pop, "terminal",
                      trait_names = "ADG", index_wts = 1.0)

  expect_no_warning(
    pop |>
      get_table("ind_ebv") |>
      dplyr::filter(model == "test_model", eval_number == 1L) |>
      add_index("terminal")
  )
})


test_that("add_index() increments index_number on successive runs", {
  pop <- make_index_pop("ai_version")
  on.exit(close_pop(pop))

  pop <- define_index(pop, "terminal", trait_names = "ADG", index_wts = 1.0)

  suppressWarnings({
    pop <- pop |> get_table("ind_ebv") |> add_index("terminal")
    pop <- pop |> get_table("ind_ebv") |> add_index("terminal")
    pop <- pop |> get_table("ind_ebv") |> add_index("terminal")
  })

  nums <- DBI::dbGetQuery(pop$db_conn,
    "SELECT DISTINCT index_number FROM ind_index ORDER BY index_number")$index_number
  expect_equal(nums, 1:3)
})


test_that("add_index() overwrite_index = TRUE resets to index_number 1", {
  pop <- make_index_pop("ai_replace")
  on.exit(close_pop(pop))

  pop <- define_index(pop, "terminal", trait_names = "ADG", index_wts = 1.0)

  suppressWarnings({
    pop <- pop |> get_table("ind_ebv") |> add_index("terminal")
    pop <- pop |> get_table("ind_ebv") |> add_index("terminal")
    pop <- pop |> get_table("ind_ebv") |>
      add_index("terminal", overwrite_index = TRUE)
  })

  nums <- DBI::dbGetQuery(pop$db_conn,
    "SELECT DISTINCT index_number FROM ind_index ORDER BY index_number")$index_number
  expect_equal(nums, 1L)
})


test_that("add_index() delete_all = TRUE clears entire ind_index", {
  pop <- make_index_pop("ai_del_all")
  on.exit(close_pop(pop))

  pop <- define_index(pop, "t1", trait_names = "ADG", index_wts = 1.0)
  pop <- define_index(pop, "t2", trait_names = "FCR", index_wts = 1.0)

  suppressWarnings({
    pop <- pop |> get_table("ind_ebv") |> add_index("t1")
    pop <- pop |> get_table("ind_ebv") |> add_index("t2")
    pop <- pop |> get_table("ind_ebv") |>
      add_index("t1", delete_all = TRUE)
  })

  # t2 rows are gone; only t1 rows remain with index_number = 1
  remaining <- DBI::dbGetQuery(pop$db_conn,
    "SELECT DISTINCT index_name FROM ind_index")$index_name
  expect_equal(remaining, "t1")

  nums <- DBI::dbGetQuery(pop$db_conn,
    "SELECT DISTINCT index_number FROM ind_index")$index_number
  expect_equal(nums, 1L)
})


test_that("add_index() supports extra user columns", {
  pop <- make_index_pop("ai_extra")
  on.exit(close_pop(pop))

  pop <- define_index(pop, "terminal", trait_names = "ADG", index_wts = 1.0)

  suppressWarnings(
    pop <- pop |>
      get_table("ind_ebv") |>
      add_index("terminal", gen = 1L)
  )

  cols <- DBI::dbListFields(pop$db_conn, "ind_index")
  expect_true("gen" %in% cols)

  gen_vals <- DBI::dbGetQuery(pop$db_conn,
    "SELECT DISTINCT gen FROM ind_index")$gen
  expect_equal(gen_vals, 1L)
})


test_that("add_index() errors when tbl is not ind_ebv", {
  pop <- make_index_pop("ai_err_tbl")
  on.exit(close_pop(pop))

  pop <- define_index(pop, "terminal", trait_names = "ADG", index_wts = 1.0)

  expect_error(
    pop |> get_table("ind_meta") |> add_index("terminal"),
    regexp = "ind_ebv"
  )
})


test_that("add_index() errors when index_name not in index_meta", {
  pop <- make_index_pop("ai_err_no_idx")
  on.exit(close_pop(pop))

  expect_error(
    suppressWarnings(
      pop |> get_table("ind_ebv") |> add_index("nonexistent")
    ),
    regexp = "not found in index_meta"
  )
})


test_that("add_index() errors when a required trait has no EBVs for some individuals", {
  pop <- make_index_pop("ai_err_partial")
  on.exit(close_pop(pop))

  pop <- define_index(pop, "terminal",
                      trait_names = c("ADG", "FCR"),
                      index_wts   = c(1.0, -1.0))

  # Remove FCR EBVs for 3 individuals to simulate missing data
  ind_ids <- DBI::dbGetQuery(pop$db_conn,
    "SELECT id_ind FROM ind_meta LIMIT 3")$id_ind
  for (id in ind_ids) {
    DBI::dbExecute(pop$db_conn,
      paste0("DELETE FROM ind_ebv WHERE id_ind = '", id,
             "' AND trait_name = 'FCR'"))
  }

  expect_error(
    pop |>
      get_table("ind_ebv") |>
      dplyr::filter(model == "test_model") |>
      add_index("terminal"),
    regexp = "missing EBVs for required traits"
  )
})


test_that("add_index() errors when duplicate (id_ind, trait_name) rows remain", {
  pop <- make_index_pop("ai_err_dup")
  on.exit(close_pop(pop))

  # Add a second evaluation for all individuals
  ind_ids <- DBI::dbGetQuery(pop$db_conn,
    "SELECT id_ind FROM ind_meta")$id_ind
  ebv2 <- data.frame(
    id_ind      = ind_ids,
    trait_name  = "ADG",
    model       = "model_2",
    ebv         = rnorm(length(ind_ids)),
    acc         = NA_real_,
    se          = NA_real_,
    eval_number = 1L,
    stringsAsFactors = FALSE
  )
  DBI::dbWriteTable(pop$db_conn, "ind_ebv", ebv2, append = TRUE)

  pop <- define_index(pop, "terminal", trait_names = "ADG", index_wts = 1.0)

  # No filter — auto-select latest should still give 1 per (id_ind, trait_name)
  # since both have eval_number = 1. This should fail because 2 models tie.
  expect_error(
    suppressWarnings(
      pop |> get_table("ind_ebv") |> add_index("terminal")
    ),
    regexp = "more than one EBV row"
  )
})


test_that("index_meta and ind_index tables exist after initialize_genome()", {
  pop <- initialize_genome(
    pop_name = "idx_init_test",
    n_loci   = 50,
    n_chr    = 2,
    chr_len_Mb = 50,
    db_path  = ":memory:"
  )
  on.exit(close_pop(pop))

  tables <- DBI::dbListTables(pop$db_conn)
  expect_true("index_meta" %in% tables)
  expect_true("ind_index"  %in% tables)
  expect_true("index_meta" %in% pop$tables)
  expect_true("ind_index"  %in% pop$tables)
})
