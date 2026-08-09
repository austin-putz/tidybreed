# Stage 1: direct long write + RAM-aware batching + one-transaction-per-call.
# These tests lock in (a) the memory/batch-size helpers, (b) batch invariance
# (any batch_size -> byte-identical ind_haplotype for a fixed seed), (c) the
# ind_crossover schema, and (d) transaction atomicity on a mid-write failure.

# ---- memory + batch-size helpers ------------------------------------------

test_that("parse_mem_size parses byte counts and suffixed strings", {
  expect_equal(parse_mem_size(1024), 1024)
  expect_equal(parse_mem_size("1KB"), 1024)
  expect_equal(parse_mem_size("512MB"), 512 * 1024^2)
  expect_equal(parse_mem_size("2g"), 2 * 1024^3)
  expect_equal(parse_mem_size("1.5G"), 1.5 * 1024^3)
  expect_error(parse_mem_size("banana"))
  expect_error(parse_mem_size(-1))
})

test_that("detect_available_memory returns a positive number or NA, never errors", {
  m <- detect_available_memory()
  expect_true(is.na(m) || (is.numeric(m) && length(m) == 1L && is.finite(m) && m > 0))
})

test_that("resolve_batch_size honors explicit batch_size, memory budget, and clamps", {
  # explicit batch_size wins and is clamped to n_total
  expect_equal(resolve_batch_size(500, 5000, batch_size = 7), 7L)
  expect_equal(resolve_batch_size(5, 5000, batch_size = 100), 5L)
  expect_error(resolve_batch_size(500, 5000, batch_size = 0))
  expect_error(resolve_batch_size(500, 5000, batch_size = -3))
  # explicit memory budget derives B and never drops below 1
  b <- resolve_batch_size(10000, 50000, max_batch_mem = "128MB")
  expect_true(b >= 1L && b <= 10000L)
  expect_equal(resolve_batch_size(10, 5000, max_batch_mem = 1), 1L)  # tiny budget -> 1
  # n_total 0 -> 0
  expect_equal(resolve_batch_size(0, 5000), 0L)
  # default path (both NULL) returns something in [1, n_total]
  d <- resolve_batch_size(200, 1000)
  expect_true(d >= 1L && d <= 200L)
})

# ---- batch invariance ------------------------------------------------------

build_pop <- function(batch_size = NULL) {
  set.seed(4242)
  pop <- open_pop(pop_name = "bi", db_name = ":memory:") |>
    define_genome(n_loci = 300, n_chr = 3, chr_len_Mb = 100)
  pop <- pop |> define_founder_haplotypes(n_haplotypes = 40, line_name = "A",
                                          method = "uniform")
  pop <- pop |> get_table("founder_haplotypes") |>
    add_founders(n_males = 6, n_females = 6, line_name = "A",
                 batch_size = batch_size)
  matings <- tibble::tibble(
    id_parent_1 = rep(paste0("A_", 1:6), length.out = 24),
    id_parent_2 = rep(paste0("A_", 7:12), length.out = 24),
    sex = rep(c("M", "F"), length.out = 24), line_name = "A")
  pop |> add_offspring(matings, batch_size = batch_size)
}

canon_hap <- function(pop) {
  h <- get_table(pop, "ind_haplotype") |> dplyr::collect() |> as.data.frame()
  h <- h[do.call(order, h[order(names(h))]), order(names(h))]
  rownames(h) <- NULL
  h
}

test_that("ind_haplotype is byte-identical across batch sizes (founders + offspring)", {
  ref <- NULL
  for (bs in list(1L, 5L, 24L, NULL)) {
    pop <- build_pop(bs)
    h <- canon_hap(pop)
    if (is.null(ref)) ref <- h else expect_equal(h, ref)
    close_pop(pop)
  }
})

test_that("batch invariance holds with a sex chromosome configured", {
  make <- function(bs) {
    set.seed(99)
    pop <- open_pop(pop_name = "bx", db_name = ":memory:") |>
      define_genome(n_loci = 160, n_chr = 4, chr_len_Mb = 100) |>
      define_chromosome("4", offspring_sex = "M",
                        from_parent_1 = 0, from_parent_2 = 1) |>
      define_founder_haplotypes(n_haplotypes = 30, line_name = "A", method = "uniform")
    pop <- pop |> get_table("founder_haplotypes") |>
      add_founders(n_males = 5, n_females = 5, line_name = "A", batch_size = bs)
    matings <- tibble::tibble(
      id_parent_1 = rep(paste0("A_", 1:5), length.out = 16),
      id_parent_2 = rep(paste0("A_", 6:10), length.out = 16),
      sex = rep(c("M", "F"), length.out = 16), line_name = "A")
    pop |> add_offspring(matings, batch_size = bs)
  }
  p1 <- make(1L);   h1 <- canon_hap(p1); close_pop(p1)
  p2 <- make(16L);  h2 <- canon_hap(p2); close_pop(p2)
  expect_equal(h1, h2)
})

# ---- ind_crossover schema --------------------------------------------------

test_that("ind_crossover is created empty, registered, and documented", {
  pop <- open_pop(pop_name = "xc", db_name = ":memory:") |>
    define_genome(n_loci = 50, n_chr = 2, chr_len_Mb = 100)
  expect_true("ind_crossover" %in% pop$tables)
  expect_true("ind_crossover" %in% DBI::dbListTables(pop$db_conn))
  xc <- get_table(pop, "ind_crossover") |> dplyr::collect()
  expect_equal(nrow(xc), 0L)
  expect_setequal(names(xc),
    c("id_crossover", "id_ind", "parent_origin", "chr", "chr_name", "pos_cM"))
  n_meta <- DBI::dbGetQuery(pop$db_conn,
    "SELECT COUNT(*) n FROM _schema_meta WHERE table_name = 'ind_crossover'")$n
  expect_gt(n_meta, 0L)
  close_pop(pop)
})

# ---- transaction atomicity -------------------------------------------------

test_that("a mid-write failure in add_offspring rolls back the whole generation", {
  set.seed(7)
  pop <- open_pop(pop_name = "tx", db_name = ":memory:") |>
    define_genome(n_loci = 80, n_chr = 2, chr_len_Mb = 100) |>
    define_founder_haplotypes(n_haplotypes = 20, line_name = "A", method = "uniform")
  pop <- pop |> get_table("founder_haplotypes") |>
    add_founders(n_males = 5, n_females = 5, line_name = "A")
  meta0 <- nrow(get_table(pop, "ind_meta") |> dplyr::collect())
  hap0  <- nrow(get_table(pop, "ind_haplotype") |> dplyr::collect())

  orig  <- get(".write_long_haplotypes", envir = asNamespace("tidybreed"))
  calls <- 0L
  assignInNamespace(".write_long_haplotypes", function(conn, long_frame) {
    calls <<- calls + 1L
    if (calls == 2L) stop("injected mid-write failure")
    orig(conn, long_frame)
  }, "tidybreed")
  on.exit(assignInNamespace(".write_long_haplotypes", orig, "tidybreed"), add = TRUE)

  matings <- tibble::tibble(
    id_parent_1 = rep("A_1", 4), id_parent_2 = paste0("A_", 6:9),
    sex = rep(c("M", "F"), 2), line_name = "A")
  expect_error(add_offspring(pop, matings, batch_size = 1L), "injected")

  expect_equal(nrow(get_table(pop, "ind_meta") |> dplyr::collect()), meta0)
  expect_equal(nrow(get_table(pop, "ind_haplotype") |> dplyr::collect()), hap0)
  close_pop(pop)
})
