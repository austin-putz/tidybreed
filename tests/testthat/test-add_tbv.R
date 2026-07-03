# Stage 2 (plans/refactor_haplotype.md) -- add_tbv() using line_origin for
# crossbreeding additive TBV: line-specific effects with population-wide
# fallback, imprinting, and per-line base_allele_freq centering.

make_lines_pop <- function(pop_name, n_loci = 10, n_chr = 1) {
  pop <- open_pop(pop_name = pop_name, db_name = ":memory:") |>
    define_genome(n_loci = n_loci, n_chr = n_chr, chr_len_Mb = 50)
  pop <- define_founder_haplotypes(pop, n_haplotypes = 20, line_name = "Duroc",
                                   method = "fixed", allele_freq = 0.5)
  pop <- define_founder_haplotypes(pop, n_haplotypes = 20, line_name = "Landrace",
                                   method = "fixed", allele_freq = 0.5)
  pop
}

# Independently recompute expected TBV from stored ind_haplotype +
# genome_effects rows, mirroring (but not literally copying) the centered SQL
# under test in add_tbv(): line-specific effect first, population-wide
# fallback otherwise, per locus.
independent_tbv <- function(pop, trait, ids) {
  hap <- DBI::dbGetQuery(pop$db_conn, paste0(
    "SELECT id_ind, locus_name, allele, line_origin FROM ind_haplotype ",
    "WHERE id_ind IN (", paste0("'", ids, "'", collapse = ", "), ")"))
  eff <- DBI::dbGetQuery(pop$db_conn, paste0(
    "SELECT locus_name, line_name, genome_value, base_allele_freq FROM genome_effects ",
    "WHERE trait_name = '", trait, "' AND genome_effect_type = 'additive'"))

  stats::setNames(vapply(ids, function(id) {
    rows <- hap[hap$id_ind == id, , drop = FALSE]
    sum(vapply(seq_len(nrow(rows)), function(i) {
      r <- rows[i, ]
      e_line <- eff[eff$locus_name == r$locus_name &
                    !is.na(eff$line_name) & !is.na(r$line_origin) &
                    eff$line_name == r$line_origin, , drop = FALSE]
      e <- if (nrow(e_line) == 1L) {
        e_line
      } else {
        eff[eff$locus_name == r$locus_name & is.na(eff$line_name), , drop = FALSE]
      }
      if (nrow(e) != 1L) stop("test helper: ambiguous or missing effect match")
      p <- if (is.na(e$base_allele_freq)) 0 else e$base_allele_freq
      (r$allele - p) * e$genome_value
    }, numeric(1)))
  }, numeric(1)), ids)
}


test_that("add_tbv() computes correct F1 crossbred TBV with line-specific effects", {
  set.seed(1001)
  pop <- make_lines_pop("tbv_f1")
  pop <- pop |> get_table("founder_haplotypes") |> dplyr::filter(line_name == "Duroc") |>
    add_founders(n_males = 2, n_females = 2, line_name = "Duroc", gen = 0L)
  pop <- pop |> get_table("founder_haplotypes") |> dplyr::filter(line_name == "Landrace") |>
    add_founders(n_males = 2, n_females = 2, line_name = "Landrace", gen = 0L)

  pop <- define_trait(pop, "ADG")
  loci <- pop |> get_table("genome_meta") |> dplyr::collect() |>
    dplyr::arrange(.data$locus_id) |> dplyr::pull(locus_name)

  pop <- pop |> get_table("genome_meta") |>
    define_additive_effects("ADG", effects = seq(0.1, 1.0, length.out = length(loci)),
                            line_name = "Duroc")
  pop <- pop |> get_table("genome_meta") |>
    define_additive_effects("ADG", effects = seq(-0.5, 0.5, length.out = length(loci)),
                            line_name = "Landrace")

  matings <- tibble::tibble(id_parent_1 = "Duroc_1", id_parent_2 = "Landrace_3",
                            sex = "M", line_name = "F1", gen = 1L)
  set.seed(2001)
  pop <- add_offspring(pop, matings)

  # F1 of two pure-line parents: parent_origin=1 rows are ALL "Duroc",
  # parent_origin=2 rows are ALL "Landrace" (recombining within a pure line
  # cannot change the line label), so this exercises the plain line-specific
  # match branch (no fallback needed).
  hap_lines <- DBI::dbGetQuery(pop$db_conn,
    "SELECT DISTINCT parent_origin, line_origin FROM ind_haplotype WHERE id_ind = 'F1_1'")
  expect_equal(hap_lines$line_origin[hap_lines$parent_origin == 1], "Duroc")
  expect_equal(hap_lines$line_origin[hap_lines$parent_origin == 2], "Landrace")

  pop <- pop |> get_table("ind_meta") |> dplyr::filter(id_ind == "F1_1") |> add_tbv("ADG")

  expected <- independent_tbv(pop, "ADG", "F1_1")
  actual <- DBI::dbGetQuery(pop$db_conn,
    "SELECT tbv_value FROM ind_tbv WHERE id_ind = 'F1_1' AND trait_name = 'ADG'")$tbv_value
  expect_equal(actual, unname(expected), tolerance = 1e-8)

  close_pop(pop)
})


test_that("add_tbv() falls back to population-wide effect when line_origin is NULL", {
  pop <- make_lines_pop("tbv_null_fallback", n_loci = 6, n_chr = 1)
  pop <- pop |> get_table("founder_haplotypes") |> dplyr::filter(line_name == "Duroc") |>
    add_founders(n_males = 1, n_females = 1, line_name = "Duroc", gen = 0L)

  pop <- define_trait(pop, "ADG")
  loci <- pop |> get_table("genome_meta") |> dplyr::collect() |>
    dplyr::arrange(.data$locus_id) |> dplyr::pull(locus_name)

  # Population-wide effect only (line_name IS NULL) -- no line-specific rows.
  pop <- pop |> get_table("genome_meta") |>
    define_additive_effects("ADG", effects = rep(1.0, length(loci)))

  # Simulate untracked line origin for this individual's haplotypes.
  DBI::dbExecute(pop$db_conn,
    "UPDATE ind_haplotype SET line_origin = NULL WHERE id_ind = 'Duroc_1'")

  pop <- pop |> get_table("ind_meta") |> dplyr::filter(id_ind == "Duroc_1") |> add_tbv("ADG")

  expected <- independent_tbv(pop, "ADG", "Duroc_1")
  actual <- DBI::dbGetQuery(pop$db_conn,
    "SELECT tbv_value FROM ind_tbv WHERE id_ind = 'Duroc_1' AND trait_name = 'ADG'")$tbv_value
  expect_equal(actual, unname(expected), tolerance = 1e-8)

  close_pop(pop)
})


test_that("add_tbv() prefers line-specific effect over population-wide for the same locus", {
  pop <- make_lines_pop("tbv_dual_row", n_loci = 6, n_chr = 1)
  pop <- pop |> get_table("founder_haplotypes") |> dplyr::filter(line_name == "Duroc") |>
    add_founders(n_males = 1, n_females = 1, line_name = "Duroc", gen = 0L)

  pop <- define_trait(pop, "ADG")
  loci <- pop |> get_table("genome_meta") |> dplyr::collect() |>
    dplyr::arrange(.data$locus_id) |> dplyr::pull(locus_name)

  # Population-wide effect at ALL loci ...
  pop <- pop |> get_table("genome_meta") |>
    define_additive_effects("ADG", effects = rep(1.0, length(loci)))
  # ... and a DIFFERENT line-specific effect at the SAME loci for "Duroc".
  pop <- pop |> get_table("genome_meta") |>
    define_additive_effects("ADG", effects = rep(5.0, length(loci)), line_name = "Duroc")

  n_eff_rows <- DBI::dbGetQuery(pop$db_conn,
    "SELECT COUNT(*) AS n FROM genome_effects WHERE trait_name = 'ADG'")$n
  expect_equal(n_eff_rows, 2L * length(loci))  # both rows coexist, per locus

  pop <- pop |> get_table("ind_meta") |> dplyr::filter(id_ind == "Duroc_1") |> add_tbv("ADG")

  expected <- independent_tbv(pop, "ADG", "Duroc_1")
  actual <- DBI::dbGetQuery(pop$db_conn,
    "SELECT tbv_value FROM ind_tbv WHERE id_ind = 'Duroc_1' AND trait_name = 'ADG'")$tbv_value
  expect_equal(actual, unname(expected), tolerance = 1e-8)
  # Sanity: the line-specific value (5.0), not a doubled/summed value, drove the
  # result -- join per-locus since base_allele_freq varies by locus.
  hap <- DBI::dbGetQuery(pop$db_conn,
    "SELECT locus_name, allele FROM ind_haplotype WHERE id_ind = 'Duroc_1'")
  eff_duroc <- DBI::dbGetQuery(pop$db_conn,
    "SELECT locus_name, base_allele_freq FROM genome_effects WHERE trait_name = 'ADG' AND line_name = 'Duroc'")
  hap <- merge(hap, eff_duroc, by = "locus_name")
  expect_equal(actual, sum((hap$allele - hap$base_allele_freq) * 5.0), tolerance = 1e-8)

  close_pop(pop)
})


test_that("add_tbv() falls back per-locus when a line has effects at only some loci", {
  pop <- make_lines_pop("tbv_partial", n_loci = 10, n_chr = 1)
  pop <- pop |> get_table("founder_haplotypes") |> dplyr::filter(line_name == "Duroc") |>
    add_founders(n_males = 1, n_females = 1, line_name = "Duroc", gen = 0L)

  pop <- define_trait(pop, "ADG")
  all_loci <- pop |> get_table("genome_meta") |> dplyr::collect() |>
    dplyr::arrange(.data$locus_id) |> dplyr::pull(locus_name)
  half_loci <- all_loci[1:5]

  pop <- pop |> get_table("genome_meta") |>
    define_additive_effects("ADG", effects = rep(1.0, length(all_loci)))
  pop <- pop |> get_table("genome_meta") |> dplyr::filter(locus_name %in% half_loci) |>
    define_additive_effects("ADG", effects = rep(3.0, length(half_loci)), line_name = "Duroc")

  pop <- pop |> get_table("ind_meta") |> dplyr::filter(id_ind == "Duroc_1") |> add_tbv("ADG")

  expected <- independent_tbv(pop, "ADG", "Duroc_1")
  actual <- DBI::dbGetQuery(pop$db_conn,
    "SELECT tbv_value FROM ind_tbv WHERE id_ind = 'Duroc_1' AND trait_name = 'ADG'")$tbv_value
  expect_equal(actual, unname(expected), tolerance = 1e-8)

  close_pop(pop)
})


test_that("add_tbv() centers each allele with its own line's base_allele_freq", {
  set.seed(5001)
  pop <- make_lines_pop("tbv_base_by_line", n_loci = 6, n_chr = 1)
  pop <- pop |> get_table("founder_haplotypes") |> dplyr::filter(line_name == "Duroc") |>
    add_founders(n_males = 3, n_females = 3, line_name = "Duroc", gen = 0L)
  pop <- pop |> get_table("founder_haplotypes") |> dplyr::filter(line_name == "Landrace") |>
    add_founders(n_males = 3, n_females = 3, line_name = "Landrace", gen = 0L)

  pop <- define_trait(pop, "ADG")
  loci <- pop |> get_table("genome_meta") |> dplyr::collect() |>
    dplyr::arrange(.data$locus_id) |> dplyr::pull(locus_name)

  duroc_tbl    <- get_table(pop, "ind_meta") |> dplyr::filter(line_name == "Duroc")
  landrace_tbl <- get_table(pop, "ind_meta") |> dplyr::filter(line_name == "Landrace")

  pop <- pop |> get_table("genome_meta") |>
    define_additive_effects("ADG", effects = rep(2.0, length(loci)), line_name = "Duroc",
                            base = "current_pop", base_tbl = duroc_tbl)
  pop <- pop |> get_table("genome_meta") |>
    define_additive_effects("ADG", effects = rep(2.0, length(loci)), line_name = "Landrace",
                            base = "current_pop", base_tbl = landrace_tbl)

  base_duroc <- DBI::dbGetQuery(pop$db_conn,
    "SELECT base_allele_freq FROM genome_effects WHERE trait_name='ADG' AND line_name='Duroc'")$base_allele_freq
  base_landrace <- DBI::dbGetQuery(pop$db_conn,
    "SELECT base_allele_freq FROM genome_effects WHERE trait_name='ADG' AND line_name='Landrace'")$base_allele_freq
  # Not a degenerate test: the two lines' realized founder allele frequencies differ.
  expect_true(any(abs(base_duroc - base_landrace) > 1e-6))

  pop <- pop |> get_table("ind_meta") |> add_tbv("ADG")

  all_ids <- DBI::dbGetQuery(pop$db_conn, "SELECT id_ind FROM ind_meta")$id_ind
  expected <- independent_tbv(pop, "ADG", all_ids)
  actual_df <- DBI::dbGetQuery(pop$db_conn,
    "SELECT id_ind, tbv_value FROM ind_tbv WHERE trait_name = 'ADG'")
  actual <- stats::setNames(actual_df$tbv_value, actual_df$id_ind)[all_ids]
  expect_equal(unname(actual), unname(expected[all_ids]), tolerance = 1e-8)

  close_pop(pop)
})


test_that("add_tbv() respects parent_origin restriction with line-specific imprinted effects", {
  pop <- make_lines_pop("tbv_imprint", n_loci = 6, n_chr = 1)
  pop <- pop |> get_table("founder_haplotypes") |> dplyr::filter(line_name == "Duroc") |>
    add_founders(n_males = 2, n_females = 2, line_name = "Duroc", gen = 0L)
  pop <- pop |> get_table("founder_haplotypes") |> dplyr::filter(line_name == "Landrace") |>
    add_founders(n_males = 2, n_females = 2, line_name = "Landrace", gen = 0L)

  pop <- define_trait(pop, "IMP", expressed_parent = "parent_1")
  loci <- pop |> get_table("genome_meta") |> dplyr::collect() |>
    dplyr::arrange(.data$locus_id) |> dplyr::pull(locus_name)

  pop <- pop |> get_table("genome_meta") |>
    define_additive_effects("IMP", effects = rep(1.0, length(loci)), line_name = "Duroc")
  pop <- pop |> get_table("genome_meta") |>
    define_additive_effects("IMP", effects = rep(4.0, length(loci)), line_name = "Landrace")

  matings <- tibble::tibble(id_parent_1 = "Duroc_1", id_parent_2 = "Landrace_3",
                            sex = "M", line_name = "F1", gen = 1L)
  set.seed(3001)
  pop <- add_offspring(pop, matings)

  pop <- pop |> get_table("ind_meta") |> dplyr::filter(id_ind == "F1_1") |> add_tbv("IMP")

  hap_po1 <- DBI::dbGetQuery(pop$db_conn,
    "SELECT locus_name, allele, line_origin FROM ind_haplotype WHERE id_ind = 'F1_1' AND parent_origin = 1")
  expect_true(all(hap_po1$line_origin == "Duroc"))  # sire's line only

  # base_allele_freq varies by locus -- join per-locus rather than a single scalar.
  eff_duroc <- DBI::dbGetQuery(pop$db_conn,
    "SELECT locus_name, base_allele_freq FROM genome_effects WHERE trait_name='IMP' AND line_name='Duroc'")
  hap_po1 <- merge(hap_po1, eff_duroc, by = "locus_name")
  expected <- sum((hap_po1$allele - hap_po1$base_allele_freq) * 1.0)
  actual <- DBI::dbGetQuery(pop$db_conn,
    "SELECT tbv_value FROM ind_tbv WHERE id_ind='F1_1' AND trait_name='IMP'")$tbv_value
  expect_equal(actual, expected, tolerance = 1e-8)

  close_pop(pop)
})


test_that("add_tbv() correctly follows line_origin through F2 recombination", {
  set.seed(4001)
  pop <- make_lines_pop("tbv_f2", n_loci = 40, n_chr = 3)
  pop <- pop |> get_table("founder_haplotypes") |> dplyr::filter(line_name == "Duroc") |>
    add_founders(n_males = 4, n_females = 4, line_name = "Duroc", gen = 0L)
  pop <- pop |> get_table("founder_haplotypes") |> dplyr::filter(line_name == "Landrace") |>
    add_founders(n_males = 4, n_females = 4, line_name = "Landrace", gen = 0L)

  pop <- define_trait(pop, "ADG")
  loci <- pop |> get_table("genome_meta") |> dplyr::collect() |>
    dplyr::arrange(.data$locus_id) |> dplyr::pull(locus_name)
  set.seed(4002)
  pop <- pop |> get_table("genome_meta") |>
    define_additive_effects("ADG", effects = stats::rnorm(length(loci)), line_name = "Duroc")
  set.seed(4003)
  pop <- pop |> get_table("genome_meta") |>
    define_additive_effects("ADG", effects = stats::rnorm(length(loci)), line_name = "Landrace")

  matings_f1 <- tibble::tibble(
    id_parent_1 = c("Duroc_1", "Duroc_2"),
    id_parent_2 = c("Landrace_5", "Landrace_6"),
    sex         = c("M", "F"),
    line_name   = "F1",
    gen         = 1L
  )
  set.seed(4004)
  pop <- add_offspring(pop, matings_f1)

  matings_f2 <- tibble::tibble(
    id_parent_1 = "F1_1", id_parent_2 = "F1_2", sex = "M", line_name = "F2", gen = 2L
  )
  set.seed(4005)
  pop <- add_offspring(pop, matings_f2)

  # Not a degenerate test: the F2 individual carries a genuine mosaic of both
  # lines (proves line_origin survived two generations of recombination).
  f2_lines <- DBI::dbGetQuery(pop$db_conn,
    "SELECT DISTINCT line_origin FROM ind_haplotype WHERE id_ind = 'F2_1'")$line_origin
  expect_true(all(c("Duroc", "Landrace") %in% f2_lines))

  pop <- pop |> get_table("ind_meta") |> dplyr::filter(id_ind == "F2_1") |> add_tbv("ADG")

  expected <- independent_tbv(pop, "ADG", "F2_1")
  actual <- DBI::dbGetQuery(pop$db_conn,
    "SELECT tbv_value FROM ind_tbv WHERE id_ind = 'F2_1' AND trait_name = 'ADG'")$tbv_value
  expect_equal(actual, unname(expected), tolerance = 1e-8)

  close_pop(pop)
})


test_that("add_tbv() correctly computes a backcross (F1 x parental line) TBV", {
  set.seed(6001)
  pop <- make_lines_pop("tbv_backcross", n_loci = 40, n_chr = 3)
  pop <- pop |> get_table("founder_haplotypes") |> dplyr::filter(line_name == "Duroc") |>
    add_founders(n_males = 4, n_females = 4, line_name = "Duroc", gen = 0L)
  pop <- pop |> get_table("founder_haplotypes") |> dplyr::filter(line_name == "Landrace") |>
    add_founders(n_males = 4, n_females = 4, line_name = "Landrace", gen = 0L)

  pop <- define_trait(pop, "ADG")
  loci <- pop |> get_table("genome_meta") |> dplyr::collect() |>
    dplyr::arrange(.data$locus_id) |> dplyr::pull(locus_name)
  set.seed(6002)
  pop <- pop |> get_table("genome_meta") |>
    define_additive_effects("ADG", effects = stats::rnorm(length(loci)), line_name = "Duroc")
  set.seed(6003)
  pop <- pop |> get_table("genome_meta") |>
    define_additive_effects("ADG", effects = stats::rnorm(length(loci)), line_name = "Landrace")

  matings_f1 <- tibble::tibble(
    id_parent_1 = "Duroc_1", id_parent_2 = "Landrace_5",
    sex = "M", line_name = "F1", gen = 1L
  )
  set.seed(6004)
  pop <- add_offspring(pop, matings_f1)

  # Backcross: F1 sire x a pure Landrace dam.
  matings_bc <- tibble::tibble(
    id_parent_1 = "F1_1", id_parent_2 = "Landrace_6", sex = "F", line_name = "BC", gen = 2L
  )
  set.seed(6005)
  pop <- add_offspring(pop, matings_bc)

  pop <- pop |> get_table("ind_meta") |> dplyr::filter(id_ind == "BC_1") |> add_tbv("ADG")

  expected <- independent_tbv(pop, "ADG", "BC_1")
  actual <- DBI::dbGetQuery(pop$db_conn,
    "SELECT tbv_value FROM ind_tbv WHERE id_ind = 'BC_1' AND trait_name = 'ADG'")$tbv_value
  expect_equal(actual, unname(expected), tolerance = 1e-8)

  close_pop(pop)
})
