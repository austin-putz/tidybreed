# extract_allele_freq(): the single place a population selection becomes p.
# plans/update_genome_effects_base_tbl.md §2.2-2.3, §3.1-3.2; tests 1-9.

# Two named pools plus founders from each; deterministic pool contents so the
# hand computations below are exact.
make_freq_pop <- function(pop_name = "af", n_loci = 6) {
  pop <- open_pop(pop_name = pop_name, db_name = ":memory:") |>
    define_genome(n_loci = n_loci, n_chr = 1, chr_len_Mb = 50)
  pop <- define_founder_haplotypes(pop, n_haplotypes = 20, line_name = "Duroc",
                                   method = "fixed", allele_freq = 0.5)
  pop <- define_founder_haplotypes(pop, n_haplotypes = 20, line_name = "Landrace",
                                   method = "fixed", allele_freq = 0.5)
  pop <- pop |> get_table("founder_haplotypes") |>
    dplyr::filter(line_name == "Duroc") |>
    add_founders(n_males = 2, n_females = 2, line_name = "Duroc", gen = 0L)
  pop <- pop |> get_table("founder_haplotypes") |>
    dplyr::filter(line_name == "Landrace") |>
    add_founders(n_males = 2, n_females = 2, line_name = "Landrace", gen = 0L)
  pop
}

# Independent recomputation: mean allele over the selected rows, per locus,
# in plain R rather than SQL, so the two sides cannot share a bug.
hand_freq <- function(rows, n_loci, key = "locus_id") {
  out <- rep(NA_real_, n_loci)
  agg <- tapply(rows$allele, rows[[key]], mean)
  out[as.integer(names(agg))] <- as.numeric(agg)
  out
}

# ── shapes ─────────────────────────────────────────────────────────────────

test_that("founder_haplotypes: pooled and line-filtered frequencies are exact", {
  set.seed(11)
  pop <- make_freq_pop("af_fh")
  on.exit(close_pop(pop))
  fh <- DBI::dbGetQuery(pop$db_conn,
    "SELECT fh.line_name, gm.locus_id, fh.allele FROM founder_haplotypes fh
     JOIN genome_meta gm USING (locus_name)")

  got_all <- pop |> get_table("founder_haplotypes") |> extract_allele_freq()
  expect_equal(got_all$allele_freq, hand_freq(fh, 6))

  got_d <- pop |> get_table("founder_haplotypes") |>
    dplyr::filter(line_name == "Duroc") |> extract_allele_freq()
  expect_equal(got_d$allele_freq, hand_freq(fh[fh$line_name == "Duroc", ], 6))
  # method = "fixed" at 0.5 is exact by construction
  expect_true(all(got_d$allele_freq == 0.5))
})

test_that("ind_haplotype: filtered copies are averaged directly", {
  set.seed(12)
  pop <- make_freq_pop("af_ih")
  on.exit(close_pop(pop))
  ih <- DBI::dbGetQuery(pop$db_conn,
    "SELECT locus_id, allele, line_origin FROM ind_haplotype")

  got <- pop |> get_table("ind_haplotype") |>
    dplyr::filter(line_origin == "Duroc") |> extract_allele_freq()
  expect_equal(got$allele_freq, hand_freq(ih[ih$line_origin == "Duroc", ], 6))
})

test_that("id_ind tables: frequency depends on which individuals, not how many rows", {
  set.seed(13)
  pop <- make_freq_pop("af_ids")
  on.exit(close_pop(pop))
  ih <- DBI::dbGetQuery(pop$db_conn,
    "SELECT id_ind, locus_id, allele FROM ind_haplotype")
  duroc_ids <- DBI::dbGetQuery(pop$db_conn,
    "SELECT id_ind FROM ind_meta WHERE line_name = 'Duroc'")$id_ind
  expected <- hand_freq(ih[ih$id_ind %in% duroc_ids, ], 6)

  via_meta <- pop |> get_table("ind_meta") |>
    dplyr::filter(line_name == "Duroc") |> extract_allele_freq()
  expect_equal(via_meta$allele_freq, expected)

  # Repeated rows per animal: five phenotype records each. Must not weight.
  pop <- define_trait(pop, "ADG", target_add_var = 1)
  pop <- pop |> get_table("genome_meta") |>
    define_additive_effects("ADG", effects = rep(1, 6), line_name = "Duroc")
  pop <- define_phenotype(pop, "ADG", mean = 0, residual_var = 1,
                          repeatable = TRUE)
  for (i in 1:5) {
    pop <- pop |> get_table("ind_meta") |> dplyr::filter(line_name == "Duroc") |>
      add_phenotype("ADG")
  }
  n_rec <- DBI::dbGetQuery(pop$db_conn,
    "SELECT COUNT(*) AS n FROM ind_phenotype")$n
  expect_equal(n_rec, 5L * length(duroc_ids))
  via_pheno <- pop |> get_table("ind_phenotype") |>
    dplyr::filter(phenotype_name == "ADG") |> extract_allele_freq()
  expect_equal(via_pheno$allele_freq, expected)

  # And the direct-copy shape agrees when it selects the same copies.
  via_hap <- pop |> get_table("ind_haplotype") |>
    dplyr::filter(line_origin == "Duroc") |> extract_allele_freq()
  expect_equal(via_hap$allele_freq, expected)
})

test_that("founder pool and founders' own copies agree on founders-only data", {
  # Not a general identity (founders are a sample of the pool); here the pool
  # is method = "fixed" at exactly 0.5 per column, so every haplotype is a
  # 0/1 pattern with the same marginal -- the sampled founders need not match.
  # What must hold: the *pool* selection reads pool rows and the *copy*
  # selection reads ind_haplotype rows, and each equals its own hand value.
  set.seed(14)
  pop <- make_freq_pop("af_agree")
  on.exit(close_pop(pop))
  fh <- DBI::dbGetQuery(pop$db_conn,
    "SELECT gm.locus_id, fh.allele FROM founder_haplotypes fh
     JOIN genome_meta gm USING (locus_name) WHERE fh.line_name = 'Landrace'")
  ih <- DBI::dbGetQuery(pop$db_conn,
    "SELECT locus_id, allele FROM ind_haplotype WHERE line_origin = 'Landrace'")
  pool <- pop |> get_table("founder_haplotypes") |>
    dplyr::filter(line_name == "Landrace") |> extract_allele_freq()
  copies <- pop |> get_table("ind_haplotype") |>
    dplyr::filter(line_origin == "Landrace") |> extract_allele_freq()
  expect_equal(pool$allele_freq,   hand_freq(fh, 6))
  expect_equal(copies$allele_freq, hand_freq(ih, 6))
})

# ── contract ───────────────────────────────────────────────────────────────

test_that("partial coverage gives NA at absent loci only; 0 and 1 are kept", {
  set.seed(15)
  pop <- make_freq_pop("af_partial")
  on.exit(close_pop(pop))

  got <- pop |> get_table("ind_haplotype") |>
    dplyr::filter(locus_id %in% c(1L, 2L)) |> extract_allele_freq()
  expect_equal(nrow(got), 6L)
  expect_false(anyNA(got$allele_freq[1:2]))
  expect_true(all(is.na(got$allele_freq[3:6])))

  # Force a fixed and an absent allele at two loci; 0 / 1 must survive as
  # values and not be confused with "no copies".
  DBI::dbExecute(pop$db_conn, "UPDATE ind_haplotype SET allele = 0 WHERE locus_id = 1")
  DBI::dbExecute(pop$db_conn, "UPDATE ind_haplotype SET allele = 1 WHERE locus_id = 2")
  got <- pop |> get_table("ind_haplotype") |>
    dplyr::filter(locus_id %in% c(1L, 2L)) |> extract_allele_freq()
  expect_identical(got$allele_freq[1:2], c(0, 1))
  expect_true(all(is.na(got$allele_freq[3:6])))
})

test_that("select() removing a required column is a package error, not a DuckDB one", {
  set.seed(16)
  pop <- make_freq_pop("af_select")
  on.exit(close_pop(pop))

  expect_error(
    pop |> get_table("ind_haplotype") |> dplyr::select(locus_id) |>
      extract_allele_freq(),
    "ind_haplotype.*missing column.*allele")
  expect_error(
    pop |> get_table("founder_haplotypes") |> dplyr::select(allele) |>
      extract_allele_freq(),
    "founder_haplotypes.*missing column.*locus_name")
  expect_error(
    pop |> get_table("ind_meta") |> dplyr::select(sex) |>
      extract_allele_freq(),
    "ind_meta.*missing column.*id_ind")

  # line_name is NOT required for founder_haplotypes (plan v4 amendment 1)
  got <- pop |> get_table("founder_haplotypes") |>
    dplyr::select(locus_name, allele) |> extract_allele_freq()
  expect_equal(nrow(got), 6L)
})

test_that("a table without id_ind and a non-table both error naming the shapes", {
  set.seed(17)
  pop <- make_freq_pop("af_shape")
  on.exit(close_pop(pop))
  expect_error(pop |> get_table("genome_meta") |> extract_allele_freq(),
               "id_ind.*Accepted shapes")
  expect_error(extract_allele_freq(data.frame(id_ind = "x")),
               "must be a tidybreed_table.*Accepted shapes")
})

test_that("empty selections: founder diagnostic lists pools; others are generic", {
  set.seed(18)
  pop <- make_freq_pop("af_empty")
  on.exit(close_pop(pop))
  expect_error(
    pop |> get_table("founder_haplotypes") |> dplyr::filter(line_name == "Nope") |>
      extract_allele_freq(),
    "No founder_haplotypes rows for this selection. Available: 'Duroc', 'Landrace'")
  expect_error(
    pop |> get_table("ind_meta") |> dplyr::filter(line_name == "Nope") |>
      extract_allele_freq(),
    "filtered base \\(ind_meta\\) contains no allele copies")
  expect_error(
    pop |> get_table("ind_haplotype") |> dplyr::filter(locus_id > 999L) |>
      extract_allele_freq(),
    "filtered base \\(ind_haplotype\\) contains no allele copies")
})

test_that("the unnamed pool is named as such in the founder diagnostic", {
  pop <- open_pop(pop_name = "af_null", db_name = ":memory:") |>
    define_genome(n_loci = 4, n_chr = 1, chr_len_Mb = 10)
  on.exit(close_pop(pop))
  pop <- define_founder_haplotypes(pop, n_haplotypes = 10)          # line_name = NULL
  expect_error(
    pop |> get_table("founder_haplotypes") |> dplyr::filter(line_name == "A") |>
      extract_allele_freq(),
    "Available: an unnamed \\(line_name = NULL\\) pool")
})

test_that("extract_allele_freq() never warns", {
  set.seed(19)
  pop <- make_freq_pop("af_nowarn")
  on.exit(close_pop(pop))
  # A pooled two-line founder table is exactly what define_additive_effects()
  # warns about on its default path; the helper is a computation and stays silent.
  expect_no_warning(pop |> get_table("founder_haplotypes") |> extract_allele_freq())
  expect_no_warning(pop |> get_table("ind_meta") |> extract_allele_freq())
  expect_no_warning(pop |> get_table("ind_haplotype") |> extract_allele_freq())
})

test_that("output is locus_id-ordered, one row per locus, typed, and read-only", {
  set.seed(20)
  pop <- make_freq_pop("af_contract")
  on.exit(close_pop(pop))
  before <- vapply(pop$tables, function(t)
    DBI::dbGetQuery(pop$db_conn, paste0("SELECT COUNT(*) AS n FROM ", t))$n,
    numeric(1))

  got <- pop |> get_table("ind_meta") |> extract_allele_freq()
  expect_s3_class(got, "tbl_df")
  expect_identical(names(got), c("locus_id", "locus_name", "allele_freq"))
  expect_type(got$locus_id, "integer")
  expect_type(got$locus_name, "character")
  expect_type(got$allele_freq, "double")
  expect_identical(got$locus_id, 1:6)
  expect_identical(got$locus_name, paste0("Locus_", 1:6))

  after <- vapply(pop$tables, function(t)
    DBI::dbGetQuery(pop$db_conn, paste0("SELECT COUNT(*) AS n FROM ", t))$n,
    numeric(1))
  expect_identical(before, after)
})

test_that("with F1s present, line_origin selects copies and ind_meta.line_name selects animals", {
  # Plan test 17. Realized Duroc-origin copies are a finite transmitted sample
  # of the Duroc pool, so they are compared with a hand computation over
  # exactly those copies -- never with the pool frequency.
  set.seed(21)
  pop <- make_freq_pop("af_f1")
  on.exit(close_pop(pop))
  pop <- define_trait(pop, "ADG", target_add_var = 1)
  pop <- pop |> get_table("genome_meta") |>
    define_additive_effects("ADG", effects = rep(1, 6),
                            base_tbl = get_table(pop, "founder_haplotypes"))
  f1 <- tibble::tibble(id_parent_1 = c("Duroc_1", "Duroc_2"),
                       id_parent_2 = c("Landrace_3", "Landrace_4"),
                       sex = c("M", "F"), line_name = "F1", gen = 1L)
  pop <- add_offspring(pop, f1)
  ih <- DBI::dbGetQuery(pop$db_conn,
    "SELECT h.id_ind, h.locus_id, h.allele, h.line_origin, m.line_name
     FROM ind_haplotype h JOIN ind_meta m USING (id_ind)")

  by_origin <- pop |> get_table("ind_haplotype") |>
    dplyr::filter(line_origin == "Duroc") |> extract_allele_freq()
  expect_equal(by_origin$allele_freq, hand_freq(ih[ih$line_origin == "Duroc", ], 6))

  by_label <- pop |> get_table("ind_meta") |>
    dplyr::filter(line_name == "F1") |> extract_allele_freq()
  f1_rows <- ih[ih$line_name == "F1", ]
  expect_equal(by_label$allele_freq, hand_freq(f1_rows, 6))
  expect_setequal(unique(f1_rows$line_origin), c("Duroc", "Landrace"))

  # Make the two selections provably different at one locus: fix every
  # Landrace-origin copy at allele 1 and every Duroc-origin copy at 0 there.
  DBI::dbExecute(pop$db_conn,
    "UPDATE ind_haplotype SET allele = CASE WHEN line_origin = 'Duroc' THEN 0 ELSE 1 END
     WHERE locus_id = 1")
  by_origin <- pop |> get_table("ind_haplotype") |>
    dplyr::filter(line_origin == "Duroc") |> extract_allele_freq()
  by_label <- pop |> get_table("ind_meta") |>
    dplyr::filter(line_name == "F1") |> extract_allele_freq()
  expect_equal(by_origin$allele_freq[1], 0)
  expect_equal(by_label$allele_freq[1], 0.5)   # one Duroc + one Landrace copy each
})
