# Individual ids reach SQL through `.with_id_view()`, never as text
# (plans/individual_id_sql_followup.md), and the reads that became unordered
# or order-dependent in that change are pinned here: explicit row order, the
# duplicate source-row refusal, and build_data_file()'s fixed-order mean.
#
# Selections use non-id predicates: an explicit filter(id_ind %in% ...) is
# rendered into SQL by dbplyr, which is the documented exception to the rule.
# Each capture also asserts the intended statement ran; "no id in the SQL"
# proves nothing about a query that never executed.

IV_LINE <- "Qzline_F1"
iv_id <- function(n) paste0(IV_LINE, "_", n)
expect_no_ids <- function(sql) {
  expect_gt(length(sql), 0L)
  expect_false(any(grepl("Qzline_F1_[0-9]", sql)))
}
ran <- function(sql, ...) {
  hit <- rep(TRUE, length(sql))
  for (p in c(...)) hit <- hit & grepl(p, sql, fixed = TRUE)
  any(hit)
}
radix <- function(x) sort(x, method = "radix")

# Founders 1-4 are male, 5-8 female, in pens 1, 2, 1, 2, ...; `chip` covers
# loci 1-5 of 20. The ids pass _9, so radix order differs from numeric order.
iv_pop <- function(n_males = 4, n_females = 4) {
  set.seed(11)
  n <- n_males + n_females
  pop <- suppressMessages(
    open_pop(pop_name = "iv", db_name = ":memory:") |>
      define_genome(n_loci = 20, n_chr = 1, chr_len_Mb = 50) |>
      define_founder_haplotypes(n_haplotypes = 20) |>
      get_table("founder_haplotypes") |>
      add_founders(n_males = n_males, n_females = n_females,
                   line_name = IV_LINE, pen = rep_len(1:2, n)))
  suppressMessages(pop |> get_table("genome_meta") |>
    dplyr::filter(locus_id <= 5L) |> define_chip("chip"))
  pop
}

# Three generations on iv_pop(): 9 = 1x5, 10 = 2x6, 11 = 1x6, 12 = 9x10,
# 13 = 9x11, so 9 and 1 are shared ancestors of 12 and 13.
iv_pedigree <- function(pop) {
  mate <- function(pop, p1, p2, sex, seed) {
    suppressMessages(add_offspring(pop, data.frame(
      id_parent_1 = iv_id(p1), id_parent_2 = iv_id(p2), sex = sex,
      line_name = IV_LINE), seed = seed))
  }
  pop <- mate(pop, c(1, 2, 1), c(5, 6, 6), c("M", "F", "F"), 1)
  mate(pop, c(9, 9), c(10, 11), c("F", "M"), 2)
}

meta_ids <- function(pop, where) DBI::dbGetQuery(pop$db_conn, paste0(
  "SELECT id_ind FROM ind_meta WHERE ", where, " ORDER BY id_ind"))$id_ind

dosages <- function(pop, ids, loci) {
  d <- DBI::dbGetQuery(pop$db_conn, paste0(
    "SELECT id_ind, locus_id, CAST(SUM(allele) AS INTEGER) AS d ",
    "FROM ind_haplotype WHERE locus_id IN (", paste(loci, collapse = ", "), ") ",
    "GROUP BY id_ind, locus_id"))
  d <- d[d$id_ind %in% ids, ]
  m <- matrix(NA_integer_, length(ids), length(loci))
  m[cbind(match(d$id_ind, ids), match(d$locus_id, loci))] <- d$d
  m
}

# ---- .with_id_view() -------------------------------------------------------

iv_conn <- function() {
  conn <- DBI::dbConnect(duckdb::duckdb())
  DBI::dbExecute(conn, "CREATE TABLE t (id_ind VARCHAR, v INTEGER)")
  DBI::dbExecute(conn, "INSERT INTO t VALUES ('a', 1), ('a', 2), ('b', 4), ('c', 8)")
  conn
}
iv_close <- function(conn) DBI::dbDisconnect(conn, shutdown = TRUE)
sum_of <- function(conn, view) DBI::dbGetQuery(conn, paste0(
  "SELECT COALESCE(SUM(v), 0) AS s FROM t JOIN ", view, " USING (id_ind)"))$s
live_views <- function(conn) {
  v <- DBI::dbGetQuery(conn,
    "SELECT view_name FROM duckdb_views() WHERE NOT internal")$view_name
  v[startsWith(v, "__id_view_")]
}

test_that(".with_id_view() has set semantics, even through a join", {
  conn <- iv_conn()
  on.exit(iv_close(conn), add = TRUE)
  # 'a' repeated must not double a's rows: 1 + 2 + 4, not 1 + 2 + 1 + 2 + 4
  expect_equal(.with_id_view(conn, c("a", "a", "b"),
                             function(v) sum_of(conn, v)), 7)
})

test_that(".with_id_view() treats empty as nobody and refuses bad input", {
  conn <- iv_conn()
  on.exit(iv_close(conn), add = TRUE)

  expect_equal(.with_id_view(conn, character(0), function(v) sum_of(conn, v)), 0)
  type <- .with_id_view(conn, character(0), function(v) DBI::dbGetQuery(
    conn, paste0("DESCRIBE ", v))$column_type)
  expect_identical(type, "VARCHAR")

  never <- function(v) stop("fun must not run")
  expect_error(.with_id_view(conn, NULL, never), "character vector")
  expect_error(.with_id_view(conn, 1:2, never), "character vector")
  expect_error(.with_id_view(conn, factor("a"), never), "character vector")
  expect_error(.with_id_view(conn, c("a", NA), never), "NA not allowed")
  expect_length(live_views(conn), 0L)
})

test_that(".with_id_view() carries an id SQL text could not", {
  conn <- iv_conn()
  on.exit(iv_close(conn), add = TRUE)
  DBI::dbExecute(conn, "INSERT INTO t VALUES ('o''brien', 16)")
  expect_equal(.with_id_view(conn, "o'brien", function(v) sum_of(conn, v)), 16)
})

test_that("nested .with_id_view() calls each read their own ids", {
  conn <- iv_conn()
  on.exit(iv_close(conn), add = TRUE)
  duckdb::duckdb_register(conn, "unrelated", data.frame(x = 1L))
  on.exit(duckdb::duckdb_unregister(conn, "unrelated"), add = TRUE, after = FALSE)

  # success: the outer view outlives the inner one
  res <- .with_id_view(conn, "a", function(outer) {
    before <- sum_of(conn, outer)
    inner  <- .with_id_view(conn, "b", function(v) sum_of(conn, v))
    c(before, inner, sum_of(conn, outer), length(live_views(conn)))
  })
  expect_equal(res, c(3, 4, 3, 1))

  # an inner error caught by the outer callback: inner gone, outer usable
  res <- .with_id_view(conn, "a", function(outer) {
    tryCatch(.with_id_view(conn, "b", function(v) stop("inner")),
             error = function(e) NULL)
    c(length(live_views(conn)), sum_of(conn, outer))
  })
  expect_equal(res, c(1, 3))

  # an error that escapes both: neither view remains, the connection still works
  expect_error(.with_id_view(conn, "a", function(outer) {
    .with_id_view(conn, "b", function(inner) stop("escaped"))
  }), "escaped")
  expect_length(live_views(conn), 0L)
  expect_equal(.with_id_view(conn, "c", function(v) sum_of(conn, v)), 8)
  expect_equal(DBI::dbGetQuery(conn, "SELECT x FROM unrelated")$x, 1L)
})

test_that("a failed registration leaves a view it does not own alone", {
  conn <- iv_conn()
  on.exit(iv_close(conn), add = TRUE)
  taken <- paste0("__id_view_", .id_view_state$n + 1L)
  duckdb::duckdb_register(conn, taken, data.frame(id_ind = "c"))

  expect_error(.with_id_view(conn, "a", function(v) stop("fun must not run")),
               "already exists")
  expect_equal(sum_of(conn, taken), 8)
  duckdb::duckdb_unregister(conn, taken)
})

test_that(".with_id_view() draws no random numbers and keeps the caller's error", {
  conn <- iv_conn()
  on.exit(iv_close(conn), add = TRUE)

  set.seed(3)
  s <- .Random.seed
  .with_id_view(conn, "a", function(v) sum_of(conn, v))
  expect_identical(.Random.seed, s)
  try(.with_id_view(conn, "a", function(v) stop("x")), silent = TRUE)
  expect_identical(.Random.seed, s)

  custom <- structure(class = c("iv_error", "error", "condition"),
                      list(message = "custom", call = NULL))
  got <- tryCatch(.with_id_view(conn, "a", function(v) stop(custom)),
                  condition = identity)
  expect_identical(got, custom)

  # A cleanup failure cannot replace it: the view is already gone at exit
  got <- tryCatch(.with_id_view(conn, "a", function(v) {
    duckdb::duckdb_unregister(conn, v)
    stop(custom)
  }), condition = identity)
  expect_identical(got, custom)
})

# ---- add_dosage() ----------------------------------------------------------

test_that("add_dosage() overwrites only the selected animals' cache", {
  pop <- iv_pop()
  on.exit(close_pop(pop), add = TRUE)
  suppressMessages(pop |> get_table("ind_meta") |> add_dosage())
  males <- meta_ids(pop, "sex = 'M'")

  sql <- record_sql(suppressMessages(pop |> get_table("ind_meta") |>
    dplyr::filter(sex == "M") |>
    add_dosage(locus_names = c("Locus_1", "Locus_2"), overwrite_dosage = TRUE)))
  expect_no_ids(sql)
  expect_true(ran(sql, "DELETE FROM ind_genotype", "__id_view_"))
  expect_true(ran(sql, "INSERT OR REPLACE INTO ind_genotype", "__id_view_"))

  cached <- DBI::dbGetQuery(pop$db_conn, "SELECT id_ind, locus_id, dosage_value
    FROM ind_genotype")
  n <- table(cached$id_ind)
  expect_true(all(n[males] == 2L))
  expect_true(all(n[setdiff(names(n), males)] == 20L))
  expect_length(n, 8L)

  sel <- cached[cached$id_ind %in% males, ]
  want <- dosages(pop, males, 1:2)
  expect_equal(as.integer(sel$dosage_value),
               want[cbind(match(sel$id_ind, males), sel$locus_id)])
})

test_that("add_dosage() unrestricted overwrite replaces the whole cache", {
  pop <- iv_pop()
  on.exit(close_pop(pop), add = TRUE)
  suppressMessages(pop |> get_table("ind_meta") |> add_dosage())
  suppressMessages(pop |> get_table("ind_meta") |>
    add_dosage(locus_names = "Locus_3", overwrite_dosage = TRUE))
  got <- DBI::dbGetQuery(pop$db_conn, "SELECT DISTINCT locus_id FROM ind_genotype")
  expect_equal(got$locus_id, 3L)
  expect_equal(DBI::dbGetQuery(pop$db_conn,
    "SELECT COUNT(*) AS n FROM ind_genotype")$n, 8)
})

test_that("add_dosage() selecting from ind_genotype itself survives its own DELETE", {
  pop <- iv_pop()
  on.exit(close_pop(pop), add = TRUE)
  suppressMessages(pop |> get_table("ind_meta") |>
    dplyr::filter(sex == "F") |> add_dosage())
  females <- meta_ids(pop, "sex = 'F'")

  # A live subquery on the selection would empty it at the DELETE
  suppressMessages(pop |> get_table("ind_genotype") |>
    add_dosage(locus_names = "Locus_1", overwrite_dosage = TRUE))
  got <- DBI::dbGetQuery(pop$db_conn,
    "SELECT id_ind, locus_id FROM ind_genotype ORDER BY id_ind")
  expect_identical(got$id_ind, females)
  expect_true(all(got$locus_id == 1L))
})

# ---- add_genotypes() -------------------------------------------------------

test_that("add_genotypes() marks exactly the requested animals, additively", {
  pop <- iv_pop()
  on.exit(close_pop(pop), add = TRUE)
  flagged <- function() meta_ids(pop, "has_chip = TRUE")

  sql <- record_sql(suppressMessages(pop |> get_table("ind_meta") |>
    dplyr::filter(sex == "M", pen == 1L) |> add_genotypes("chip")))
  expect_no_ids(sql)
  expect_true(ran(sql, "UPDATE ind_meta SET has_chip", "__id_view_"))
  first <- meta_ids(pop, "sex = 'M' AND pen = 1")
  expect_identical(flagged(), first)

  suppressMessages(pop |> get_table("ind_meta") |>
    dplyr::filter(sex == "F", pen == 2L) |> add_genotypes("chip"))
  expect_identical(flagged(),
                   radix(c(first, meta_ids(pop, "sex = 'F' AND pen = 2"))))
})

# ---- extract_genotypes() ---------------------------------------------------

test_that("extract_genotypes() intersects chip and filter, in id order", {
  pop <- iv_pop(6, 6)
  on.exit(close_pop(pop), add = TRUE)
  suppressMessages(pop |> get_table("ind_meta") |>
    dplyr::filter(pen == 1L) |> add_genotypes("chip"))

  got <- NULL
  sql <- record_sql(got <- pop |> get_table("ind_meta") |>
    dplyr::filter(sex == "F") |> extract_genotypes(chip_name = "chip"))
  expect_no_ids(sql)
  expect_true(ran(sql, "FROM ind_meta WHERE has_chip", "__id_view_", "ORDER BY id_ind"))
  expect_true(ran(sql, "FROM ind_haplotype", "__id_view_"))

  want <- meta_ids(pop, "pen = 1 AND sex = 'F'")
  expect_gt(length(want), 1L)
  expect_identical(got$id_ind, want)
  expect_identical(names(got), c("id_ind", paste0("locus_", 1:5)))
  expect_identical(unname(as.matrix(got[, -1])), dosages(pop, want, 1:5))
})

test_that("extract_genotypes() orders the unfiltered chip and haplotype paths", {
  pop <- iv_pop(6, 6)
  on.exit(close_pop(pop), add = TRUE)
  suppressMessages(pop |> get_table("ind_meta") |> add_genotypes("chip"))

  all_ids <- meta_ids(pop, "TRUE")
  expect_identical(all_ids, radix(all_ids))
  expect_identical(
    (pop |> get_table("ind_meta") |> extract_genotypes(chip_name = "chip"))$id_ind,
    all_ids)

  sql <- record_sql(got <- pop |> get_table("ind_meta") |>
    extract_genotypes(loci_tbl = get_table(pop, "genome_meta") |>
                        dplyr::filter(locus_id %in% c(2L, 7L))))
  expect_true(ran(sql, "SELECT DISTINCT id_ind FROM ind_haplotype ORDER BY id_ind"))
  expect_identical(got$id_ind, all_ids)
  expect_identical(names(got), c("id_ind", "locus_2", "locus_7"))
  expect_identical(unname(as.matrix(got[, -1])), dosages(pop, all_ids, c(2L, 7L)))
})

test_that("extract_genotypes() reads the effects_tbl loci for a filtered subset", {
  pop <- iv_pop(6, 6)
  on.exit(close_pop(pop), add = TRUE)
  pop <- suppressMessages(with_additive_target(pop, "T1", 1))
  suppressMessages(pop |> get_table("genome_meta") |>
    dplyr::filter(locus_id %in% c(4L, 9L, 16L)) |> define_additive_effects("T1"))

  sql <- record_sql(got <- pop |> get_table("ind_meta") |>
    dplyr::filter(sex == "M") |>
    extract_genotypes(effects_tbl = get_table(pop, "genome_effect_loci")))
  expect_no_ids(sql)
  want <- meta_ids(pop, "sex = 'M'")
  expect_identical(got$id_ind, want)
  expect_identical(names(got), c("id_ind", "locus_4", "locus_9", "locus_16"))
  expect_identical(unname(as.matrix(got[, -1])), dosages(pop, want, c(4L, 9L, 16L)))
})

# ---- assert_ploidy_2() -----------------------------------------------------

test_that("assert_ploidy_2() reads only the given animals, through a view", {
  pop <- iv_pop()
  on.exit(close_pop(pop), add = TRUE)
  DBI::dbExecute(pop$db_conn, paste0(
    "UPDATE ind_meta SET ploidy = 4 WHERE id_ind = '", iv_id(1), "'"))

  sql <- record_sql(expect_null(assert_ploidy_2(pop, iv_id(2:3))))
  expect_no_ids(sql)
  expect_true(ran(sql, "ploidy != 2", "__id_view_"))

  expect_error(assert_ploidy_2(pop, iv_id(1:2)), "1 individual\\(s\\)")
  expect_error(assert_ploidy_2(pop), "1 individual\\(s\\)")
  expect_null(assert_ploidy_2(pop, character(0)))
})

# ---- ebv_parent_avg() ------------------------------------------------------

test_that("parent average reads each parent's latest EBV within the model", {
  pop <- iv_pedigree(iv_pop())
  on.exit(close_pop(pop), add = TRUE)
  pop <- suppressMessages(with_additive_target(pop, "ADG", 1))

  # Parents 9 and 10: two evaluations under m1, and a decoy under m2
  e <- data.frame(
    id_ebv      = 1:6,
    id_ind      = iv_id(c(9, 9, 10, 10, 9, 10)),
    trait_name  = "ADG",
    model       = c("m1", "m1", "m1", "m1", "m2", "m2"),
    ebv_value   = c(100, 2, 100, 4, -50, -50),
    eval_number = c(1L, 2L, 1L, 2L, 3L, 3L))
  DBI::dbAppendTable(pop$db_conn, "ind_ebv", e)

  sql <- record_sql(suppressMessages(suppressWarnings(
    pop |> get_table("ind_meta") |>
      dplyr::filter(sex == "F", !is.na(id_parent_1)) |>
      add_ebv("ADG", parent_avg = TRUE, model = "m1"))))
  expect_no_ids(sql)
  expect_true(ran(sql, "SELECT id_ind, id_parent_1, id_parent_2", "__id_view_"))
  expect_true(ran(sql, "FROM ind_ebv e", "__id_view_"))

  # The new evaluation is number 4: eval_number counts across models
  got <- DBI::dbGetQuery(pop$db_conn,
    "SELECT id_ind, ebv_value FROM ind_ebv WHERE eval_number = 4 AND model = 'm1'")
  # Animal 12 = 9 x 10: mean of the latest m1 EBVs. 10 and 11 have parents
  # with no EBV at all, so they are stored as NA, not as a wrong number.
  expect_equal(got$ebv_value[got$id_ind == iv_id(12)], 3)
  expect_setequal(got$id_ind, iv_id(c(10, 11, 12)))
  expect_true(all(is.na(got$ebv_value[got$id_ind != iv_id(12)])))
})

# ---- trace_pedigree() ------------------------------------------------------

test_that("trace_pedigree() follows each generation, once, in id order", {
  pop <- iv_pedigree(iv_pop())
  on.exit(close_pop(pop), add = TRUE)

  for (case in list(list(0L, c(12, 13)),
                    list(1L, c(12, 13, 9, 10, 11)),
                    list(2L, c(12, 13, 9, 10, 11, 1, 2, 5, 6)),
                    list(9L, c(12, 13, 9, 10, 11, 1, 2, 5, 6)))) {
    sql <- record_sql(ped <- trace_pedigree(pop, iv_id(c(12, 13)), case[[1]]))
    expect_no_ids(sql)
    expect_true(ran(sql, "SELECT id_ind, id_parent_1", "__id_view_", "ORDER BY id_ind"))
    expect_identical(ped$id_ind, radix(iv_id(case[[2]])))
  }

  ped <- trace_pedigree(pop, iv_id(12), 2L)
  expect_identical(ped$id_parent_1[ped$id_ind == iv_id(12)], iv_id(9))
  expect_identical(ped$id_parent_2[ped$id_ind == iv_id(12)], iv_id(10))
  expect_identical(ped$id_parent_1[ped$id_ind == iv_id(1)], "0")
  expect_identical(ped$id_parent_2[ped$id_ind == iv_id(1)], "0")
})

# ---- build_data_file() -----------------------------------------------------

put_pheno <- function(pop, id_ind, value, name = "ADG") {
  DBI::dbAppendTable(pop$db_conn, "ind_phenotype", data.frame(
    id_phenotype = seq_along(id_ind), id_ind = id_ind, phenotype_name = name,
    pheno_value = value, pheno_number = 1L))
}

# ADG with a class and a covariate from a non-ind_meta table. `herd_info` has
# one row per animal, with values that differ from animal to animal.
iv_effects_pop <- function() {
  pop <- iv_pop() |> with_additive_target("ADG", 1) |> suppressMessages()
  pop <- define_phenotype(pop, "ADG", type = "continuous", mean = 1,
                          residual_var = 0.5) |> suppressMessages()
  n <- 1:8
  DBI::dbWriteTable(pop$db_conn, "herd_info", data.frame(
    id_ind = iv_id(n), herd = paste("herd", n), wt = n * 10))
  pop <- define_effect_fixed_class(pop, "ADG", effect_name = "herd",
    source_column = "herd", source_table = "herd_info",
    levels = stats::setNames(rep(0, 8), paste("herd", n))) |> suppressMessages()
  define_effect_fixed_cov(pop, "ADG", effect_name = "wt", source_column = "wt",
    slope = 0, source_table = "herd_info") |> suppressMessages()
}

test_that("build_data_file() reads the selected animals and their fixed effects", {
  pop <- iv_effects_pop()
  on.exit(close_pop(pop), add = TRUE)
  eval_dir <- withr::local_tempdir()
  # Animals 1-3 selected; 4 phenotyped but not selected; record 3 is excluded
  put_pheno(pop, iv_id(c(1, 1, 1, 2, 3, 4)), c(1, 3, 100, 5, 7, 9))

  sql <- record_sql(res <- suppressMessages(
    build_data_file(pop, iv_id(1:3), "ADG", eval_dir,
                    pheno_ids = c(1L, 2L, 4L, 5L, 6L))))
  expect_no_ids(sql)
  expect_true(ran(sql, "FROM ind_phenotype", "__id_view_", "list_avg"))
  expect_true(ran(sql, "FROM herd_info", "__id_view_"))

  d <- res$data
  expect_identical(d$id_ind, iv_id(1:3))
  expect_identical(d$herd, c("herd_1", "herd_2", "herd_3"))
  expect_equal(d$wt, c(10, 20, 30))
  expect_equal(d$ADG, c(2, 5, 7))
  expect_identical(names(res$col_map), c("mu", "id_ind", "herd", "wt", "ADG"))
  expect_equal(unname(res$col_map), 1:5)
  file <- read.table(file.path(eval_dir, "data.txt"), header = TRUE)
  expect_equal(file$ADG, d$ADG)
  expect_identical(file$herd, d$herd)
})

test_that("build_data_file() refuses duplicate source rows it would have to read", {
  pop <- iv_effects_pop()
  on.exit(close_pop(pop), add = TRUE)
  eval_dir <- withr::local_tempdir()
  put_pheno(pop, iv_id(1:3), c(1, 2, 3))

  # A duplicate outside the animals read is not this call's concern
  DBI::dbExecute(pop$db_conn, paste0(
    "INSERT INTO herd_info VALUES ('", iv_id(8), "', 'herd 8b', 1)"))
  expect_no_error(suppressMessages(build_data_file(pop, iv_id(1:3), "ADG", eval_dir)))

  DBI::dbExecute(pop$db_conn, paste0(
    "INSERT INTO herd_info VALUES ('", iv_id(2), "', 'herd 2b', 1)"))
  expect_error(suppressMessages(build_data_file(pop, iv_id(1:3), "ADG", eval_dir)),
               "Effect 'herd', 'wt': source table 'herd_info' has duplicate")
})

test_that("build_data_file() averages records exactly as AVG() would", {
  pop <- iv_pop()
  on.exit(close_pop(pop), add = TRUE)
  eval_dir <- withr::local_tempdir()
  put_pheno(pop, iv_id(c(1, 2, 2, 2, 3, 3)),
            c(0.30000000000000004, NA, 2, 4, NA, NA))

  res <- suppressMessages(build_data_file(pop, iv_id(1:3), "ADG", eval_dir))
  d <- res$data
  # data.txt would round this to 0.3, so it is asserted on the frame
  expect_identical(d$ADG[d$id_ind == iv_id(1)], 0.30000000000000004)
  expect_identical(d$ADG[d$id_ind == iv_id(2)], 3)  # the NULL is skipped
  expect_identical(d$ADG[d$id_ind == iv_id(3)], 0)  # all NULL: missing value 0
})

test_that("build_data_file() treats empty pheno_ids as no records, never all", {
  pop <- iv_pop()
  on.exit(close_pop(pop), add = TRUE)
  eval_dir <- withr::local_tempdir()
  put_pheno(pop, iv_id(1:3), c(1, 2, 3))

  expect_error(build_data_file(pop, iv_id(1:3), "ADG", eval_dir,
                               pheno_ids = integer(0)),
               "No phenotypic records found")
  expect_false(file.exists(file.path(eval_dir, "data.txt")))
})

test_that("build_data_file() reports animals with records per trait", {
  pop <- iv_pop()
  on.exit(close_pop(pop), add = TRUE)
  eval_dir <- withr::local_tempdir()
  # Only ADG is recorded; LS (say, litter size) is not measured yet
  put_pheno(pop, iv_id(1:3), c(1, 2, 3))

  msgs <- testthat::capture_messages(
    res <- build_data_file(pop, iv_id(1:3), c("ADG", "LS"), eval_dir))
  expect_true(any(grepl("Records [ADG]: 3 animals", msgs, fixed = TRUE)))
  expect_true(any(grepl(
    "Records [LS]: 0 animals (EBVs come only from genetic correlations",
    msgs, fixed = TRUE)))
  expect_false(any(grepl("Records [ADG]: 3 animals (", msgs, fixed = TRUE)))
  # ADG's line comes first: trait order
  expect_lt(grep("[ADG]", msgs, fixed = TRUE), grep("[LS]", msgs, fixed = TRUE))
  expect_identical(res$data$LS, c(0, 0, 0))
  expect_equal(res$data$ADG, c(1, 2, 3))
})

test_that("build_data_file()'s count skips animals whose records are all missing", {
  pop <- iv_pop()
  on.exit(close_pop(pop), add = TRUE)
  eval_dir <- withr::local_tempdir()
  put_pheno(pop, iv_id(c(1, 2, 3, 3)), c(1, 2, NA, NA))

  msgs <- testthat::capture_messages(
    res <- build_data_file(pop, iv_id(1:3), "ADG", eval_dir))
  expect_true(any(grepl("Records [ADG]: 2 animals", msgs, fixed = TRUE)))
  expect_identical(nrow(res$data), 3L)  # animal 3 is still a (missing) row
})

test_that("build_data_file() refuses when no trait has a non-missing value", {
  pop <- iv_pop()
  on.exit(close_pop(pop), add = TRUE)
  eval_dir <- withr::local_tempdir()
  put_pheno(pop, iv_id(1:3), c(NA, NA, NA))

  expect_error(build_data_file(pop, iv_id(1:3), c("ADG", "LS"), eval_dir),
               "No phenotypic records found")
  expect_false(file.exists(file.path(eval_dir, "data.txt")))
})

test_that("add_ebv() refuses a phenotype filter that matches no records", {
  pop <- iv_pop() |> with_additive_target("ADG", 1) |> suppressMessages()
  on.exit(close_pop(pop), add = TRUE)
  put_pheno(pop, iv_id(1:3), c(1, 2, 3))

  expect_error(
    pop |> get_table("ind_meta") |>
      add_ebv("ADG", software = "blupf90",
              phenotype = pop |> get_table("ind_phenotype") |>
                dplyr::filter(pheno_value > 1e9)),
    "`phenotype` filter matched no records in ind_phenotype")
  expect_equal(
    DBI::dbGetQuery(pop$db_conn, "SELECT COUNT(*) AS n FROM ind_ebv")$n, 0)
})

test_that("build_data_file()'s mean is bit-identical across thread counts", {
  conn_pop <- suppressMessages(
    open_pop(pop_name = "iv_mean", db_name = ":memory:") |>
      define_genome(n_loci = 2, n_chr = 1, chr_len_Mb = 1))
  on.exit(close_pop(conn_pop), add = TRUE)
  conn <- conn_pop$db_conn

  # Interleaved by record, full-precision values: raw runif() cannot fail (a
  # few 32-bit values sum exactly) and records contiguous by animal barely can.
  set.seed(1)
  ids <- sprintf("X_%06d", seq_len(100000L))
  put_pheno(conn_pop, rep(ids, times = 5L), rnorm(500000L, 10, 2))

  threads <- DBI::dbGetQuery(conn, "SELECT current_setting('threads') AS t")$t
  on.exit(DBI::dbExecute(conn, paste0("SET threads = ", threads)),
          add = TRUE, after = FALSE)
  build <- function(n) {
    DBI::dbExecute(conn, paste0("SET threads = ", n))
    dir <- withr::local_tempdir()
    res <- suppressMessages(build_data_file(conn_pop, ids, "ADG", dir))
    list(data = res$data, lines = readLines(file.path(dir, "data.txt")))
  }
  one  <- build(1L)
  four <- build(4L)
  expect_identical(four$data, one$data)
  expect_identical(four$lines, one$lines)
})

# ---- write_geno_file() -----------------------------------------------------

test_that("write_geno_file() writes only genotyped animals, in id order", {
  pop <- iv_pop(6, 6)
  on.exit(close_pop(pop), add = TRUE)
  eval_dir <- withr::local_tempdir()
  suppressMessages(pop |> get_table("ind_meta") |>
    dplyr::filter(pen == 1L) |> add_genotypes("chip"))

  # A broader pedigree than the genotyped set, repeated ids included
  ped_ids <- iv_id(c(12, 11, 10, 9, 2, 1, 1))
  sql <- record_sql(res <- suppressMessages(
    write_geno_file(pop, ped_ids, "chip", eval_dir)))
  expect_no_ids(sql)
  expect_true(ran(sql, "FROM ind_meta", "__id_view_", "ORDER BY id_ind"))
  expect_true(ran(sql, "FROM ind_haplotype", "__id_view_"))

  # pen alternates, so the odd numbers are genotyped
  want <- radix(iv_id(c(11, 9, 1)))
  expect_identical(res$geno_ids, want)
  expect_identical(res$n_loci, 5L)

  lines <- readLines(res$geno_file)
  m <- dosages(pop, want, 1:5)
  expect_identical(lines, sprintf("%-*s%s", res$id_width, want,
                                  apply(m, 1, paste0, collapse = "")))
})

test_that("write_geno_file() returns empty when nobody in the set is genotyped", {
  pop <- iv_pop()
  on.exit(close_pop(pop), add = TRUE)
  eval_dir <- withr::local_tempdir()
  DBI::dbExecute(pop$db_conn,
    "ALTER TABLE ind_meta ADD COLUMN has_chip BOOLEAN DEFAULT FALSE")

  res <- suppressMessages(write_geno_file(pop, iv_id(1:2), "chip", eval_dir))
  expect_null(res$geno_file)
  expect_identical(res$geno_ids, character(0))
  expect_false(file.exists(file.path(eval_dir, "genotype.txt")))
})

# ---- the views do not outlive the calls ------------------------------------

test_that("no id view remains after the public paths, success or failure", {
  pop <- iv_pop()
  on.exit(close_pop(pop), add = TRUE)
  suppressMessages(pop |> get_table("ind_meta") |>
    dplyr::filter(sex == "M") |> add_dosage())
  expect_length(live_views(pop$db_conn), 0L)
  expect_error(pop |> get_table("ind_meta") |> dplyr::filter(sex == "Z") |>
                 extract_genotypes(loci_tbl = get_table(pop, "genome_meta")),
               "No individuals")
  expect_length(live_views(pop$db_conn), 0L)
})
