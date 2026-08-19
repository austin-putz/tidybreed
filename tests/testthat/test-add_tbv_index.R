# Coverage for add_tbv()'s index_names / true-index block (R/add_tbv.R) and its
# scattered input guards. The crossbreeding line-precedence join is covered
# separately in test-add_tbv.R -- keep the two concerns in separate files.

make_index_pop <- function(pop_name = "ix", n_ind = 10, traits = c("ADG", "BW")) {
  pop <- open_pop(pop_name = pop_name, db_name = ":memory:") |>
    define_genome(n_loci = 200, n_chr = 2, chr_len_Mb = 100) |>
    define_founder_haplotypes(n_haplotypes = 40, method = "beta")
  pop <- pop |>
    get_table("founder_haplotypes") |>
    add_founders(n_males = n_ind / 2, n_females = n_ind / 2, line_name = "A")
  for (t in traits) {
    pop <- define_trait(pop, t, target_add_var = 0.25)
    sel <- get_table(pop, "genome_meta") |>
      dplyr::collect() |> dplyr::slice_sample(n = 30) |> dplyr::pull(locus_name)
    pop <- pop |>
      get_table("genome_meta") |>
      dplyr::filter(locus_name %in% sel) |>
      define_additive_effects(t)
  }
  pop
}

# Wide TBV lookup keyed by id_ind -- reshape() mangles the value column name.
tbv_wide <- function(pop, ids, trait) {
  tb <- dplyr::collect(get_table(pop, "ind_tbv"))
  tb$tbv_value[match(paste(ids, trait), paste(tb$id_ind, tb$trait_name))]
}

true_index <- function(pop, ids, weight_type) {
  ti <- dplyr::collect(get_table(pop, "ind_true_index"))
  ti <- ti[ti$weight_type == weight_type, ]
  ti$true_index_value[match(ids, ti$id_ind)]
}

all_ids <- function(pop) {
  sort(dplyr::pull(dplyr::collect(get_table(pop, "ind_meta")), id_ind))
}


# ── 1. The index_names block ────────────────────────────────────────────────

test_that("add_tbv(type = 'index') writes the index-weighted sum of TBVs", {
  set.seed(11001)
  pop <- make_index_pop("ix_index")
  pop <- define_index(pop, "sel", c("ADG", "BW"), index_wts = c(2, 3))
  pop <- suppressMessages(
    pop |> get_table("ind_meta") |> add_tbv(index_names = "sel", type = "index"))

  ids <- all_ids(pop)
  ti  <- dplyr::collect(get_table(pop, "ind_true_index"))
  expect_equal(nrow(ti), length(ids))
  expect_true(all(ti$weight_type == "index"))
  expect_true(all(ti$index_name == "sel"))
  expect_equal(true_index(pop, ids, "index"),
               2 * tbv_wide(pop, ids, "ADG") + 3 * tbv_wide(pop, ids, "BW"))

  close_pop(pop)
})


test_that("add_tbv(type = 'economic') uses economic_weight, not index_weight", {
  set.seed(11002)
  pop <- make_index_pop("ix_econ")
  # Deliberately different index vs economic weights so a column mix-up fails.
  pop <- define_index(pop, "sel", c("ADG", "BW"),
                      index_wts = c(2, 3), economic_wts = c(10, -1))
  pop <- suppressMessages(
    pop |> get_table("ind_meta") |> add_tbv(index_names = "sel", type = "economic"))

  ids <- all_ids(pop)
  ti  <- dplyr::collect(get_table(pop, "ind_true_index"))
  expect_true(all(ti$weight_type == "economic"))
  expect_equal(true_index(pop, ids, "economic"),
               10 * tbv_wide(pop, ids, "ADG") - 1 * tbv_wide(pop, ids, "BW"))

  close_pop(pop)
})


test_that("add_tbv(type = 'both') writes one index and one economic row per individual", {
  set.seed(11003)
  pop <- make_index_pop("ix_both")
  pop <- define_index(pop, "sel", c("ADG", "BW"),
                      index_wts = c(2, 3), economic_wts = c(10, -1))
  pop <- suppressMessages(
    pop |> get_table("ind_meta") |> add_tbv(index_names = "sel", type = "both"))

  ids <- all_ids(pop)
  ti  <- dplyr::collect(get_table(pop, "ind_true_index"))
  expect_equal(nrow(ti), 2L * length(ids))
  expect_setequal(unique(ti$weight_type), c("index", "economic"))

  expect_equal(true_index(pop, ids, "index"),
               2 * tbv_wide(pop, ids, "ADG") + 3 * tbv_wide(pop, ids, "BW"))
  expect_equal(true_index(pop, ids, "economic"),
               10 * tbv_wide(pop, ids, "ADG") - 1 * tbv_wide(pop, ids, "BW"))

  close_pop(pop)
})


test_that("add_tbv() type defaults to 'index'", {
  set.seed(11004)
  pop <- make_index_pop("ix_default")
  pop <- define_index(pop, "sel", c("ADG", "BW"),
                      index_wts = c(2, 3), economic_wts = c(10, -1))
  # type omitted entirely -- exercises the match.arg() default.
  pop <- suppressMessages(
    pop |> get_table("ind_meta") |> add_tbv(index_names = "sel"))

  ids <- all_ids(pop)
  ti  <- dplyr::collect(get_table(pop, "ind_true_index"))
  expect_true(all(ti$weight_type == "index"))
  expect_equal(true_index(pop, ids, "index"),
               2 * tbv_wide(pop, ids, "ADG") + 3 * tbv_wide(pop, ids, "BW"))

  close_pop(pop)
})


test_that("add_tbv(overwrite_index = FALSE) skips individuals that already have a value", {
  set.seed(11005)
  pop <- make_index_pop("ix_noover")
  pop <- define_index(pop, "sel", c("ADG", "BW"), index_wts = c(2, 3))
  pop <- suppressMessages(
    pop |> get_table("ind_meta") |> add_tbv(index_names = "sel"))

  ids    <- all_ids(pop)
  before <- true_index(pop, ids, "index")

  # ind_true_index has no SQL UNIQUE constraint -- this skip is the only thing
  # standing between a re-run and silently duplicated rows.
  pop <- suppressMessages(
    pop |> get_table("ind_meta") |> add_tbv(index_names = "sel"))

  ti <- dplyr::collect(get_table(pop, "ind_true_index"))
  expect_equal(nrow(ti), length(ids))
  expect_equal(true_index(pop, ids, "index"), before)

  close_pop(pop)
})


test_that("add_tbv(overwrite_index = TRUE) recomputes in place after weights change", {
  set.seed(11006)
  pop <- make_index_pop("ix_over")
  pop <- define_index(pop, "sel", c("ADG", "BW"), index_wts = c(2, 3))
  pop <- suppressMessages(
    pop |> get_table("ind_meta") |> add_tbv(index_names = "sel"))

  ids <- all_ids(pop)
  expect_equal(true_index(pop, ids, "index"),
               2 * tbv_wide(pop, ids, "ADG") + 3 * tbv_wide(pop, ids, "BW"))

  pop <- suppressMessages(
    define_index(pop, "sel", c("ADG", "BW"), index_wts = c(-1, 5), overwrite = TRUE))
  pop <- suppressMessages(
    pop |> get_table("ind_meta") |> add_tbv(index_names = "sel", overwrite_index = TRUE))

  ti <- dplyr::collect(get_table(pop, "ind_true_index"))
  expect_equal(nrow(ti), length(ids))          # replaced, not appended
  expect_equal(true_index(pop, ids, "index"),
               -1 * tbv_wide(pop, ids, "ADG") + 5 * tbv_wide(pop, ids, "BW"))

  close_pop(pop)
})


test_that("add_tbv() handles multiple index_names in one call independently", {
  set.seed(11007)
  pop <- make_index_pop("ix_multi")
  pop <- define_index(pop, "idxa", c("ADG", "BW"), index_wts = c(1, 0))
  pop <- define_index(pop, "idxb", c("ADG", "BW"), index_wts = c(0, 4))
  pop <- suppressMessages(
    pop |> get_table("ind_meta") |> add_tbv(index_names = c("idxa", "idxb")))

  ids <- all_ids(pop)
  ti  <- dplyr::collect(get_table(pop, "ind_true_index"))
  expect_setequal(unique(ti$index_name), c("idxa", "idxb"))
  expect_equal(nrow(ti), 2L * length(ids))

  pick <- function(nm) {
    sub <- ti[ti$index_name == nm, ]
    sub$true_index_value[match(ids, sub$id_ind)]
  }
  expect_equal(pick("idxa"), 1 * tbv_wide(pop, ids, "ADG"))
  expect_equal(pick("idxb"), 4 * tbv_wide(pop, ids, "BW"))

  close_pop(pop)
})


test_that("add_tbv() writes true index rows only for the filtered subset", {
  set.seed(11008)
  pop <- make_index_pop("ix_subset")
  pop <- define_index(pop, "sel", c("ADG", "BW"), index_wts = c(2, 3))
  pop <- suppressMessages(
    pop |> get_table("ind_meta") |> dplyr::filter(sex == "M") |>
      add_tbv(index_names = "sel"))

  males <- sort(dplyr::pull(
    dplyr::collect(dplyr::filter(get_table(pop, "ind_meta"), sex == "M")), id_ind))
  ti <- dplyr::collect(get_table(pop, "ind_true_index"))

  expect_setequal(ti$id_ind, males)
  expect_equal(true_index(pop, males, "index"),
               2 * tbv_wide(pop, males, "ADG") + 3 * tbv_wide(pop, males, "BW"))

  close_pop(pop)
})


test_that("add_tbv() errors on an index name that is not in index_meta", {
  set.seed(11009)
  pop <- make_index_pop("ix_unknown")
  expect_error(
    suppressMessages(pop |> get_table("ind_meta") |> add_tbv(index_names = "nope")),
    "not found in index_meta"
  )
  close_pop(pop)
})


test_that("add_tbv(type = 'economic') errors when economic_weight was never supplied", {
  set.seed(11010)
  pop <- make_index_pop("ix_naecon")
  # No economic_wts -> define_index() omits the column from the INSERT and
  # index_meta.economic_weight has no SQL DEFAULT, so it reads back as NA.
  pop <- define_index(pop, "noecon", c("ADG", "BW"), index_wts = c(1, 2))

  em <- DBI::dbGetQuery(pop$db_conn,
    "SELECT economic_weight FROM index_meta WHERE index_name = 'noecon'")
  expect_true(all(is.na(em$economic_weight)))

  # define_trait() also writes a global index_name IS NULL row with
  # economic_weight = 0. Different row class -- must not rescue the named index.
  glob <- DBI::dbGetQuery(pop$db_conn,
    "SELECT economic_weight FROM index_meta WHERE index_name IS NULL")
  expect_true(all(glob$economic_weight == 0))

  expect_error(
    suppressMessages(pop |> get_table("ind_meta") |>
                       add_tbv(index_names = "noecon", type = "economic")),
    "Supply them via define_index"
  )
  close_pop(pop)
})


test_that("add_tbv() rejects an index name that is not a valid SQL identifier", {
  set.seed(11011)
  pop <- make_index_pop("ix_badname")
  expect_error(
    suppressMessages(pop |> get_table("ind_meta") |> add_tbv(index_names = "3bad")),
    "Invalid index name"
  )
  close_pop(pop)
})


test_that("add_tbv() errors instead of silently writing nothing when index traits have no TBVs", {
  set.seed(11012)
  pop <- make_index_pop("ix_notbv")
  # Index covers BW only, but the call computes ADG only -- so no individual has
  # a BW row in ind_tbv. Before the guard this wrote zero rows and returned
  # cleanly, violating the documented "all index traits must be in trait_name".
  pop <- define_index(pop, "sel", "BW", index_wts = 1)

  expect_error(
    suppressMessages(pop |> get_table("ind_meta") |>
                       add_tbv("ADG", index_names = "sel")),
    "No TBVs found for index 'sel'"
  )
  expect_equal(nrow(dplyr::collect(get_table(pop, "ind_true_index"))), 0L)

  close_pop(pop)
})


test_that("add_tbv() errors when only some index traits have TBVs", {
  set.seed(11013)
  pop <- make_index_pop("ix_partialtbv")
  pop <- define_index(pop, "sel", c("ADG", "BW"), index_wts = c(1, 1))
  # ADG TBVs exist; BW ones do not -> the per-trait NA check fires.
  pop <- suppressMessages(pop |> get_table("ind_meta") |> add_tbv("ADG"))

  expect_error(
    suppressMessages(pop |> get_table("ind_meta") |>
                       add_tbv("ADG", index_names = "sel")),
    "TBVs missing for index 'sel'"
  )
  close_pop(pop)
})


# ── 2. Scattered input guards ───────────────────────────────────────────────

test_that("add_tbv() defaults trait_name to every trait in trait_meta", {
  set.seed(12001)
  pop <- make_index_pop("ix_alltraits")
  pop <- suppressMessages(pop |> get_table("ind_meta") |> add_tbv())

  tb <- dplyr::collect(get_table(pop, "ind_tbv"))
  expect_setequal(unique(tb$trait_name), c("ADG", "BW"))
  expect_equal(nrow(tb), 2L * length(all_ids(pop)))

  close_pop(pop)
})


test_that("add_tbv() errors when trait_meta is empty and no trait_name is given", {
  pop <- make_test_pop("ix_notraits", n_loci = 20, n_chr = 1)
  expect_error(
    pop |> get_table("ind_meta") |> add_tbv(),
    "No traits found in trait_meta"
  )
  close_pop(pop)
})


test_that("add_tbv() rejects a non-scalar custom field", {
  set.seed(12002)
  pop <- make_index_pop("ix_vecfield")
  expect_error(
    pop |> get_table("ind_meta") |> add_tbv("ADG", note = c("a", "b")),
    "must be a scalar"
  )
  close_pop(pop)
})


test_that("add_tbv() errors when the filtered table has no id_ind column", {
  set.seed(12003)
  pop <- make_index_pop("ix_noid")
  expect_error(
    pop |> get_table("genome_meta") |> dplyr::filter(chr == 1L) |> add_tbv("ADG"),
    "must contain 'id_ind'"
  )
  close_pop(pop)
})


test_that("add_tbv() warns and no-ops when the filter matches no individuals", {
  set.seed(12004)
  pop <- make_index_pop("ix_nomatch")
  expect_warning(
    res <- pop |> get_table("ind_meta") |>
      dplyr::filter(id_ind == "does_not_exist") |> add_tbv("ADG"),
    "No individuals matched"
  )
  expect_s3_class(res, "tidybreed_pop")
  expect_equal(nrow(dplyr::collect(get_table(pop, "ind_tbv"))), 0L)
  close_pop(pop)
})


test_that("add_tbv() errors on a trait that is not in trait_meta", {
  set.seed(12005)
  pop <- make_index_pop("ix_notrait")
  expect_error(
    pop |> get_table("ind_meta") |> add_tbv("NOPE"),
    "Traits not found: NOPE"
  )
  close_pop(pop)
})


test_that("add_tbv() errors on a trait with no additive effects", {
  set.seed(12006)
  pop <- make_index_pop("ix_noqtl")
  pop <- define_trait(pop, "NOQTL", target_add_var = 0.25)
  expect_error(
    pop |> get_table("ind_meta") |> add_tbv("NOQTL"),
    "No additive effects found"
  )
  close_pop(pop)
})


test_that("add_tbv() errors, naming chromosome inheritance, when an individual carries no QTL-bearing chromosome", {
  set.seed(12007)
  pop <- open_pop(pop_name = "ix_ychr", db_name = ":memory:") |>
    define_genome(n_loci = 10, n_chr = 2, chr_names = c("1", "Y"), chr_len_Mb = 50) |>
    define_chromosome("Y", offspring_sex = "M", from_parent_1 = 1, from_parent_2 = 0) |>
    define_chromosome("Y", offspring_sex = "F", from_parent_1 = 0, from_parent_2 = 0) |>
    define_chromosome("Y", recombines = FALSE) |>
    define_founder_haplotypes(n_haplotypes = 20, method = "fixed", allele_freq = 0.5)

  pop <- pop |> get_table("founder_haplotypes") |>
    add_founders(n_males = 2, n_females = 2, line_name = "A")

  pop <- define_trait(pop, "YTRAIT")
  y_loci <- get_table(pop, "genome_meta") |> dplyr::collect() |>
    dplyr::filter(chr_name == "Y") |> dplyr::pull(locus_name)
  expect_gt(length(y_loci), 0)

  # QTL on Y only. scale_to_target = FALSE is required: assert_qtl_autosomal()
  # rejects non-autosomal QTL on the rescaling path.
  pop <- pop |> get_table("genome_meta") |>
    dplyr::filter(locus_name %in% y_loci) |>
    define_additive_effects("YTRAIT", effects = rep(1.0, length(y_loci)),
                            scale_to_target = FALSE)

  # Females carry no Y at all -> zero ind_haplotype rows at every QTL.
  n_f_rows <- DBI::dbGetQuery(pop$db_conn, paste0(
    "SELECT COUNT(*) AS n FROM ind_haplotype h JOIN ind_meta m USING (id_ind) ",
    "WHERE m.sex = 'F' AND h.locus_name IN ('",
    paste(y_loci, collapse = "', '"), "')"))$n
  expect_equal(n_f_rows, 0)

  expect_error(
    pop |> get_table("ind_meta") |> add_tbv("YTRAIT"),
    "chr_inheritance"
  )
  # Males do carry Y, so restricting to them succeeds.
  pop <- suppressMessages(
    pop |> get_table("ind_meta") |> dplyr::filter(sex == "M") |> add_tbv("YTRAIT"))
  expect_equal(nrow(dplyr::collect(get_table(pop, "ind_tbv"))), 2L)

  close_pop(pop)
})


test_that("add_tbv() writes scalar custom fields to ind_tbv", {
  set.seed(12008)
  pop <- make_index_pop("ix_extracol")
  pop <- suppressMessages(
    pop |> get_table("ind_meta") |> add_tbv("ADG", eval_note = "x"))

  tb <- dplyr::collect(get_table(pop, "ind_tbv"))
  expect_true("eval_note" %in% names(tb))
  expect_true(all(tb$eval_note == "x"))

  close_pop(pop)
})
