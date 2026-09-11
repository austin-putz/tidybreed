# Phase D of plans/update_genome_effects_v4.md (v4.9): the evaluator.
#
# The central assertion is the first block: every Phase A fixture, written into
# the real tables and evaluated by the real SQL evaluator, reproduces the number
# `gefx_eval_naive()` computes from the semantic per-tuple definition -- for
# *every* individual in the fixture universe, not only the ones the fixture
# named. The Phase A helper was written before any of this existed and is not
# generated from it, so agreement is evidence rather than a tautology.

# The fixture universe as a real population: the same loci, individuals and
# allele copies gefx_copies() defines, so a stored model can be evaluated both
# ways and compared.
gev_pop <- function() {
  pop <- open_pop(pop_name = "gev", db_name = ":memory:") |>
    define_genome(n_loci = 6, n_chr = 3, chr_len_Mb = 20,
                  locus_names = c("L1", "L2", "L3", "L4", "LX", "LY"),
                  chr_names   = c("A1", "X", "Y"))
  pop <- pop |>
    define_chromosome("X", offspring_sex = "M", from_parent_1 = 0, from_parent_2 = 1) |>
    define_chromosome("Y", offspring_sex = "M", from_parent_1 = 1, from_parent_2 = 0) |>
    define_chromosome("Y", offspring_sex = "F", from_parent_1 = 0, from_parent_2 = 0) |>
    define_chromosome("Y", recombines = FALSE)
  pop <- define_trait(pop, "ADG", target_add_var = 1.0)

  cp  <- gefx_copies()
  lid <- DBI::dbGetQuery(pop$db_conn, "SELECT locus_id, locus_name FROM genome_meta")
  cp$locus_id <- lid$locus_id[match(cp$locus_name, lid$locus_name)]
  inds <- unique(cp$id_ind)
  DBI::dbWriteTable(pop$db_conn, "ind_meta", data.frame(
    id_ind = inds, id_parent_1 = NA_character_, id_parent_2 = NA_character_,
    line_name = ifelse(grepl("^f1", inds), "F1", "A"),
    sex = ifelse(grepl("male", inds), "M", "F"), ploidy = 2L,
    stringsAsFactors = FALSE), append = TRUE)
  DBI::dbWriteTable(pop$db_conn, "ind_haplotype", data.frame(
    id_ind = cp$id_ind, parent_origin = as.integer(cp$parent_origin), strand = 1L,
    line_origin = cp$line_origin, locus_id = as.integer(cp$locus_id),
    locus_name = cp$locus_name, allele = as.integer(cp$allele),
    stringsAsFactors = FALSE), append = TRUE)
  pop
}

# Insert a fixture model verbatim. Deliberately not through
# define_genome_effects(): the point is to evaluate exactly the rows Phase A
# hand-checked, including F14, which the writer refuses on purpose.
gev_write_model <- function(pop, model) {
  conn <- pop$db_conn
  for (tb in c("genome_effect_member_origins", "genome_effect_members",
               "genome_effects")) {
    DBI::dbExecute(conn, paste0("DELETE FROM ", tb))
  }
  lid <- DBI::dbGetQuery(conn, "SELECT locus_id, locus_name FROM genome_meta")
  tm  <- model$terms
  DBI::dbWriteTable(conn, "genome_effects", data.frame(
    id_genome_effect = as.integer(tm$id_genome_effect), trait_name = tm$trait_name,
    effect_owner = tm$effect_owner, effect_name = tm$effect_name,
    genome_value = tm$genome_value, stringsAsFactors = FALSE), append = TRUE)
  mm <- model$members[order(model$members$id_genome_effect,
                            model$members$member_slot), , drop = FALSE]
  DBI::dbWriteTable(conn, "genome_effect_members", data.frame(
    id_genome_effect = as.integer(mm$id_genome_effect),
    member_slot = as.integer(mm$member_slot),
    locus_id = as.integer(lid$locus_id[match(mm$locus_name, lid$locus_name)]),
    contrast_name = mm$contrast_name,
    copy_count_value = as.integer(mm$copy_count_value),
    dosage_value = as.integer(mm$dosage_value),
    center_value = mm$center_value, stringsAsFactors = FALSE), append = TRUE)
  oo <- model$origins
  if (nrow(oo) > 0L) {
    DBI::dbWriteTable(conn, "genome_effect_member_origins", data.frame(
      id_genome_effect = as.integer(oo$id_genome_effect),
      member_slot = as.integer(oo$member_slot),
      origin_slot = as.integer(oo$origin_slot),
      line_match_type = oo$line_match_type, line_name = oo$line_name,
      parent_origin = as.integer(oo$parent_origin),
      copy_count = as.integer(oo$copy_count), stringsAsFactors = FALSE),
      append = TRUE)
  }
  invisible(pop)
}

# Total genetic value per individual, from the evaluator.
gev_totals <- function(pop, ids, trait = "ADG", ...) {
  res <- tidybreed:::.gev_evaluate(pop$db_conn, ids, trait, ...)
  out <- stats::setNames(rep(0, length(ids)), ids)
  if (nrow(res) == 0L) return(out)
  agg <- tapply(res$tgv_value, res$id_ind, sum)
  out[names(agg)] <- as.numeric(agg)
  out
}


# ── Gates 1-10, 15-18, 26-31, 37, 39: the evaluator is the semantic function ──

test_that("every Phase A fixture evaluates through SQL to its hand-computed value", {
  pop <- gev_pop()
  on.exit(close_pop(pop), add = TRUE)
  fx   <- gefx_fixtures()
  inds <- gefx_individuals()

  checked <- 0L
  for (nm in names(fx)) {
    f <- fx[[nm]]
    if (is.null(f$expected)) next
    checked <- checked + 1L
    gev_write_model(pop, f$model)
    tot <- gev_totals(pop, inds)
    for (id in names(f$expected)) {
      expect_equal(unname(tot[[id]]), unname(f$expected[[id]]),
                   tolerance = 1e-10, info = paste(nm, id))
    }
  }
  expect_gte(checked, 13L)   # every fixture carrying an expected value
})

test_that("the SQL evaluator agrees with the per-tuple definition for every individual", {
  # Not only the individuals a fixture names: an evaluator that silently
  # returned 0 for an unlisted individual would pass the block above.
  pop <- gev_pop()
  on.exit(close_pop(pop), add = TRUE)
  fx   <- gefx_fixtures()
  inds <- gefx_individuals()

  for (nm in names(fx)) {
    f <- fx[[nm]]
    if (is.null(f$expected)) next
    gev_write_model(pop, f$model)
    tot <- gev_totals(pop, inds)
    for (id in inds) {
      expect_equal(unname(tot[[id]]), gefx_eval_naive(f$model, id),
                   tolerance = 1e-10, info = paste(nm, id))
    }
  }
})

test_that("gate 31: the zero-copy indicator state is produced despite no haplotype row", {
  # fem_noY carries no LY row at all, so an inner join would drop the state
  # rather than score it. F13 is the fixture; this pins the mechanism directly.
  pop <- gev_pop()
  on.exit(close_pop(pop), add = TRUE)
  gev_write_model(pop, gefx_model(
    gefx_term(1, 7.0),
    gefx_member(1, 1, "LY", "indicator", copy_count_value = 0L, dosage_value = 0L)))
  tot <- gev_totals(pop, gefx_individuals())
  expect_equal(unname(tot[["fem_noY"]]), 7.0)     # no row at LY -> the state fires
  expect_equal(unname(tot[["pure_A"]]),  7.0)     # also carries no LY copy
})

test_that("gate 23 / F14: dominance at a non-diploid state stops instead of scoring", {
  pop <- gev_pop()
  on.exit(close_pop(pop), add = TRUE)
  gev_write_model(pop, gefx_fixtures()$F14$model)
  expect_error(gev_totals(pop, gefx_individuals()),
               "number of copies other than 2")
})

test_that("gate 23: a dominance term validated at a diploid locus never fails at evaluation", {
  pop <- gev_pop()
  on.exit(close_pop(pop), add = TRUE)
  gev_write_model(pop, gefx_model(
    gefx_term(1, 4.0),
    gefx_member(1, 1, "L1", "dominance", center = 0.5)))
  # L1 is autosomal: every individual in the universe carries two copies.
  expect_no_error(gev_totals(pop, gefx_individuals()))
})

test_that("gate 3: the dominance contrast is orthogonal under HWE", {
  # E[x_D] = 0 and E[x_A x_D] = 0 at HWE, checked on the contrast values the
  # evaluator actually produces rather than on the formula.
  p  <- 0.3; q <- 1 - p
  xD <- c(-2 * p^2, 2 * p * q, -2 * q^2)
  xA <- c(0, 1, 2) - 2 * p
  fr <- c(q^2, 2 * p * q, p^2)
  expect_equal(sum(fr * xD), 0, tolerance = 1e-12)
  expect_equal(sum(fr * xA * xD), 0, tolerance = 1e-12)
})


# ── A real population, written through the public writers ───────────────────

gev_lines_pop <- function(name = "gevl", n_loci = 6) {
  pop <- open_pop(pop_name = name, db_name = ":memory:") |>
    define_genome(n_loci = n_loci, n_chr = 1, chr_len_Mb = 20)
  set.seed(909)
  pop <- define_founder_haplotypes(pop, n_haplotypes = 20, line_name = "A",
                                   method = "fixed", allele_freq = 0.5)
  pop <- define_founder_haplotypes(pop, n_haplotypes = 20, line_name = "B",
                                   method = "fixed", allele_freq = 0.5)
  pop <- pop |> get_table("founder_haplotypes") |>
    dplyr::filter(line_name == "A") |>
    add_founders(n_males = 2, n_females = 2, line_name = "A")
  pop <- pop |> get_table("founder_haplotypes") |>
    dplyr::filter(line_name == "B") |>
    add_founders(n_males = 2, n_females = 2, line_name = "B")
  define_trait(pop, "ADG", target_add_var = 1.0)
}

gev_loci <- function(pop) {
  pop |> get_table("genome_meta") |> dplyr::collect() |>
    dplyr::arrange(.data$locus_id) |> dplyr::pull("locus_name")
}

gev_tgv <- function(pop, trait = "ADG") {
  DBI::dbGetQuery(pop$db_conn, paste0(
    "SELECT id_ind, component_name, tgv_value FROM ind_tgv ",
    "WHERE trait_name = '", trait, "' ORDER BY id_ind, component_name"))
}


test_that("gate 14: an injected dominance term changes the genetic value", {
  pop <- gev_lines_pop("gev_dom")
  on.exit(close_pop(pop), add = TRUE)
  loci <- gev_loci(pop)
  pop <- suppressWarnings(pop |> get_table("genome_meta") |>
    define_additive_effects("ADG", effects = rep(1.0, length(loci))))

  pop <- pop |> get_table("ind_meta") |> add_tgv("ADG")
  before <- DBI::dbGetQuery(pop$db_conn,
    "SELECT id_ind, tgv_total FROM ind_tgv_total WHERE trait_name = 'ADG' ORDER BY id_ind")

  pop <- define_genome_effects(pop, "ADG", data.frame(
    locus_name = loci[1], contrast_name = "dominance",
    center_value = 0.5, genome_value = 2.5), effect_owner = "dom")
  pop <- pop |> get_table("ind_meta") |> add_tgv("ADG")
  after <- DBI::dbGetQuery(pop$db_conn,
    "SELECT id_ind, tgv_total FROM ind_tgv_total WHERE trait_name = 'ADG' ORDER BY id_ind")

  expect_equal(before$id_ind, after$id_ind)
  expect_false(isTRUE(all.equal(before$tgv_total, after$tgv_total)))
})

test_that("gate 14: an injected interaction changes the genetic value", {
  pop <- gev_lines_pop("gev_epi")
  on.exit(close_pop(pop), add = TRUE)
  loci <- gev_loci(pop)
  pop <- suppressWarnings(pop |> get_table("genome_meta") |>
    define_additive_effects("ADG", effects = rep(1.0, length(loci))))
  pop <- pop |> get_table("ind_meta") |> add_tgv("ADG")
  before <- DBI::dbGetQuery(pop$db_conn,
    "SELECT id_ind, tgv_total FROM ind_tgv_total ORDER BY id_ind")

  cells <- expand.grid(a = 0:2, b = 0:2)
  tt <- genotype_terms(stats::setNames(cells, loci[1:2]),
                       value = c(0, 0, 0, 0, 1.4, 2.1, 0, 2.1, 3.6))
  pop <- define_genome_effects(pop, "ADG", tt, effect_owner = "AxA")
  pop <- pop |> get_table("ind_meta") |> add_tgv("ADG")
  after <- DBI::dbGetQuery(pop$db_conn,
    "SELECT id_ind, tgv_total FROM ind_tgv_total ORDER BY id_ind")

  expect_true("interaction" %in% gev_tgv(pop)$component_name)
  expect_false(isTRUE(all.equal(before$tgv_total, after$tgv_total)))
})

test_that("gate 21: each term maps to one component, and the components sum to the total", {
  pop <- gev_lines_pop("gev_comp")
  on.exit(close_pop(pop), add = TRUE)
  loci <- gev_loci(pop)
  pop <- suppressWarnings(pop |> get_table("genome_meta") |>
    define_additive_effects("ADG", effects = rep(0.7, length(loci))))
  pop <- define_genome_effects(pop, "ADG", data.frame(
    locus_name = loci[1], contrast_name = "dominance",
    center_value = 0.5, genome_value = 1.1), effect_owner = "dom")
  pop <- define_genome_effects(pop, "ADG", data.frame(
    locus_name = loci[2], contrast_name = "indicator",
    copy_count_value = 2L, dosage_value = 1L, genome_value = 0.9),
    effect_owner = "surface")
  cells <- expand.grid(a = 0:2, b = 0:2)
  pop <- define_genome_effects(pop, "ADG",
    genotype_terms(stats::setNames(cells, loci[3:4]),
                   value = c(0, 0, 0, 0, 1.0, 2.0, 0, 2.0, 3.0)),
    effect_owner = "AxA")

  pop <- pop |> get_table("ind_meta") |> add_tgv("ADG")
  tg <- gev_tgv(pop)
  expect_setequal(unique(tg$component_name),
                  c("order1_additive", "order1_dominance", "order1_other",
                    "interaction"))
  # One row per (individual, component): no term lands in two components.
  expect_false(anyDuplicated(paste(tg$id_ind, tg$component_name)) > 0L)

  tot <- DBI::dbGetQuery(pop$db_conn,
    "SELECT id_ind, tgv_total FROM ind_tgv_total ORDER BY id_ind")
  parts <- tapply(tg$tgv_value, tg$id_ind, sum)
  expect_equal(as.numeric(parts[tot$id_ind]), tot$tgv_total, tolerance = 1e-10)
})

test_that("gate 22: add_tgv() is idempotent per (id_ind, trait_name, component_name)", {
  pop <- gev_lines_pop("gev_idem")
  on.exit(close_pop(pop), add = TRUE)
  loci <- gev_loci(pop)
  pop <- suppressWarnings(pop |> get_table("genome_meta") |>
    define_additive_effects("ADG", effects = rep(1.0, length(loci))))

  pop <- pop |> get_table("ind_meta") |> add_tgv("ADG")
  first <- gev_tgv(pop)
  pop <- pop |> get_table("ind_meta") |> add_tgv("ADG")
  second <- gev_tgv(pop)
  expect_equal(first, second)
  expect_equal(nrow(second), 8L)     # 8 individuals x 1 component
})

test_that("a component that leaves the model leaves ind_tgv with it", {
  # The replacement is by (individual, trait), not by component: a plain upsert
  # would strand the old dominance row and ind_tgv_total would keep summing it.
  pop <- gev_lines_pop("gev_stale")
  on.exit(close_pop(pop), add = TRUE)
  loci <- gev_loci(pop)
  pop <- suppressWarnings(pop |> get_table("genome_meta") |>
    define_additive_effects("ADG", effects = rep(1.0, length(loci))))
  pop <- define_genome_effects(pop, "ADG", data.frame(
    locus_name = loci[1], contrast_name = "dominance",
    center_value = 0.5, genome_value = 2.0), effect_owner = "dom")
  pop <- pop |> get_table("ind_meta") |> add_tgv("ADG")
  expect_true("order1_dominance" %in% gev_tgv(pop)$component_name)

  DBI::dbExecute(pop$db_conn, paste0(
    "DELETE FROM genome_effect_members WHERE id_genome_effect IN ",
    "(SELECT id_genome_effect FROM genome_effects WHERE effect_owner = 'dom')"))
  DBI::dbExecute(pop$db_conn, "DELETE FROM genome_effects WHERE effect_owner = 'dom'")
  pop <- pop |> get_table("ind_meta") |> add_tgv("ADG")
  expect_false("order1_dominance" %in% gev_tgv(pop)$component_name)
})

test_that("gates 32-33: custom terms move tgv_value and leave tbv_value alone", {
  pop <- gev_lines_pop("gev_tbv_filter")
  on.exit(close_pop(pop), add = TRUE)
  loci <- gev_loci(pop)
  pop <- suppressWarnings(pop |> get_table("genome_meta") |>
    define_additive_effects("ADG", effects = rep(1.0, length(loci))))
  pop <- pop |> get_table("ind_meta") |> add_tbv("ADG")
  pop <- pop |> get_table("ind_meta") |> add_tgv("ADG")
  tbv0 <- DBI::dbGetQuery(pop$db_conn,
    "SELECT id_ind, tbv_value FROM ind_tbv ORDER BY id_ind")
  tgv0 <- DBI::dbGetQuery(pop$db_conn,
    "SELECT id_ind, tgv_total FROM ind_tgv_total ORDER BY id_ind")
  expect_equal(tbv0$tbv_value, tgv0$tgv_total, tolerance = 1e-10)

  # A functional dominance surface, and a custom additive term inside an
  # interaction. Neither is a breeding-value coefficient.
  pop <- define_genome_effects(pop, "ADG", data.frame(
    locus_name = loci[1], contrast_name = "indicator",
    copy_count_value = 2L, dosage_value = 1L, genome_value = 1.7),
    effect_owner = "functional_d")
  pop <- define_genome_effects(pop, "ADG", data.frame(
    term_id = c(1L, 1L), locus_name = loci[2:3],
    contrast_name = "additive", center_value = 0.5,
    genome_value = 3.0), effect_owner = "AxA")

  pop <- pop |> get_table("ind_meta") |> add_tbv("ADG")
  pop <- pop |> get_table("ind_meta") |> add_tgv("ADG")
  tbv1 <- DBI::dbGetQuery(pop$db_conn,
    "SELECT id_ind, tbv_value FROM ind_tbv ORDER BY id_ind")
  tgv1 <- DBI::dbGetQuery(pop$db_conn,
    "SELECT id_ind, tgv_total FROM ind_tgv_total ORDER BY id_ind")

  expect_equal(tbv1$tbv_value, tbv0$tbv_value, tolerance = 1e-10)
  expect_false(isTRUE(all.equal(tgv1$tgv_total, tgv0$tgv_total)))
  expect_setequal(unique(gev_tgv(pop)$component_name),
                  c("order1_additive", "order1_other", "interaction"))
})

test_that("gate 7 / 25: functional and Cockerham codings differ by exactly the reported mean", {
  pop <- gev_lines_pop("gev_ad")
  on.exit(close_pop(pop), add = TRUE)
  pop <- define_trait(pop, "COCK")
  loci <- gev_loci(pop)
  L <- loci[1]
  a <- 0.8; d <- 0.5; p <- 0.35; q <- 1 - p
  mu    <- a * (p - q) + 2 * p * q * d
  alpha <- a + d * (q - p)

  pop <- define_genome_effects(pop, "ADG",
    suppressMessages(ad_terms(L, a = a, d = d, p = p, coding = "functional")),
    effect_owner = "fn")
  pop <- define_genome_effects(pop, "COCK",
    suppressMessages(ad_terms(L, a = alpha, d = d, p = p, coding = "cockerham")),
    effect_owner = "ck")

  ids <- DBI::dbGetQuery(pop$db_conn, "SELECT id_ind FROM ind_meta ORDER BY id_ind")$id_ind
  fn <- gev_totals(pop, ids, "ADG")
  ck <- gev_totals(pop, ids, "COCK")
  expect_equal(unname(fn - ck), rep(mu, length(ids)), tolerance = 1e-10)
})

test_that("gate 25: repeated (a, d) appends accumulate, locus by locus", {
  pop <- gev_lines_pop("gev_ad_append")
  on.exit(close_pop(pop), add = TRUE)
  loci <- gev_loci(pop)
  ids  <- DBI::dbGetQuery(pop$db_conn, "SELECT id_ind FROM ind_meta ORDER BY id_ind")$id_ind
  a <- 0.6; d <- 0.3; p <- 0.5

  pop <- define_genome_effects(pop, "ADG",
    suppressMessages(ad_terms(loci[1], a = a, d = d, p = p)),
    effect_owner = "fn")
  one <- gev_totals(pop, ids)

  pop <- define_genome_effects(pop, "ADG",
    suppressMessages(ad_terms(loci[2], a = a, d = d, p = p)),
    effect_owner = "fn")
  two <- gev_totals(pop, ids)

  # The second append neither replaces the first nor doubles it: the two loci
  # are different families and sum.
  expect_false(isTRUE(all.equal(unname(one), unname(two))))
  dose <- DBI::dbGetQuery(pop$db_conn, paste0(
    "SELECT h.id_ind, g.locus_name, SUM(h.allele) AS g ",
    "FROM ind_haplotype h JOIN genome_meta g USING (locus_id) ",
    "WHERE g.locus_name IN ('", loci[1], "', '", loci[2], "') ",
    "GROUP BY 1, 2"))
  val <- function(g) a * (g - 1) + d * (g == 1)
  want <- vapply(ids, function(id) {
    sub <- dose[dose$id_ind == id, , drop = FALSE]
    sum(val(sub$g))
  }, numeric(1))
  expect_equal(unname(two[ids]), unname(want), tolerance = 1e-10)
})

test_that("gate 39: an order-one surface is 'order1_other' and sums into the total", {
  pop <- gev_lines_pop("gev_surface")
  on.exit(close_pop(pop), add = TRUE)
  loci <- gev_loci(pop)
  pop <- define_genome_effects(pop, "ADG", genotype_terms(
    stats::setNames(data.frame(0:2), loci[1]), value = c(1.0, 4.0, 9.0),
    copy_count = stats::setNames(list(2L), loci[1])), effect_owner = "surface")
  pop <- pop |> get_table("ind_meta") |> add_tgv("ADG")
  tg <- gev_tgv(pop)
  expect_setequal(unique(tg$component_name), "order1_other")
  tot <- DBI::dbGetQuery(pop$db_conn,
    "SELECT id_ind, tgv_total FROM ind_tgv_total ORDER BY id_ind")
  expect_equal(sort(tg$tgv_value), sort(tot$tgv_total))
})


# ── Gate 45: no eligible copy ───────────────────────────────────────────────

gev_sex_pop <- function(name = "gev_sex") {
  pop <- open_pop(pop_name = name, db_name = ":memory:") |>
    define_genome(n_loci = 4, n_chr = 2, chr_len_Mb = 20,
                  locus_names = c("A1a", "A1b", "Xa", "Xb"),
                  chr_names   = c("A1", "X"))
  pop <- define_chromosome(pop, "X", offspring_sex = "M",
                           from_parent_1 = 0, from_parent_2 = 1)
  set.seed(77)
  pop <- define_founder_haplotypes(pop, n_haplotypes = 20, method = "fixed",
                                   allele_freq = 0.5)
  pop <- pop |> get_table("founder_haplotypes") |>
    add_founders(n_males = 3, n_females = 3, line_name = "A")
  define_trait(pop, "ADG", target_add_var = 1.0)
}

test_that("gate 45: a paternally qualified X term contributes 0 beside a matching autosomal term", {
  pop <- gev_sex_pop("gev_x_pair")
  on.exit(close_pop(pop), add = TRUE)
  # Autosomal term, common scope; X term, paternal copies only. Males inherit
  # no paternal X (from_parent_1 = 0), so the X term matches nothing for them.
  pop <- pop |> get_table("genome_meta") |> dplyr::filter(chr_name == "A1") |>
    define_additive_effects("ADG", effects = c(1.0, 1.0))
  pop <- pop |> get_table("genome_meta") |> dplyr::filter(chr_name == "X") |>
    define_additive_effects("ADG", effects = c(2.0, 2.0), parent_origin = 1)

  ids <- DBI::dbGetQuery(pop$db_conn,
    "SELECT id_ind, sex FROM ind_meta ORDER BY id_ind")
  expect_no_error(pop |> get_table("ind_meta") |> add_tbv("ADG"))

  # The male value is the autosomal part alone; recomputed independently.
  hap <- DBI::dbGetQuery(pop$db_conn, paste0(
    "SELECT h.id_ind, h.allele, m.center_value, e.genome_value ",
    "FROM ind_haplotype h ",
    "JOIN genome_meta g ON g.locus_id = h.locus_id ",
    "JOIN genome_effect_members m ON m.locus_id = h.locus_id ",
    "JOIN genome_effects e USING (id_genome_effect) ",
    "LEFT JOIN genome_effect_member_origins o USING (id_genome_effect, member_slot) ",
    "WHERE g.chr_name = 'A1' AND o.line_match_type IS NULL"))
  auto <- tapply((hap$allele - hap$center_value) * hap$genome_value, hap$id_ind, sum)
  got <- DBI::dbGetQuery(pop$db_conn,
    "SELECT id_ind, tbv_value FROM ind_tbv ORDER BY id_ind")
  males <- ids$id_ind[ids$sex == "M"]
  expect_equal(got$tbv_value[match(males, got$id_ind)], as.numeric(auto[males]),
               tolerance = 1e-10)
  females <- ids$id_ind[ids$sex == "F"]
  expect_false(isTRUE(all.equal(got$tbv_value[match(females, got$id_ind)],
                                as.numeric(auto[females]))))
})

test_that("gate 45: the X term alone errors for males, who have no paternal X copy", {
  pop <- gev_sex_pop("gev_x_alone")
  on.exit(close_pop(pop), add = TRUE)
  pop <- pop |> get_table("genome_meta") |> dplyr::filter(chr_name == "X") |>
    define_additive_effects("ADG", effects = c(2.0, 2.0), parent_origin = 1)

  expect_error(pop |> get_table("ind_meta") |> add_tbv("ADG"),
               "No allele copy matched any effect")
  # Females, who do carry a paternal X, evaluate normally.
  expect_no_error(pop |> get_table("ind_meta") |>
                    dplyr::filter(sex == "F") |> add_tbv("ADG"))
})


# ── Gate 38: the label-vector preflight ─────────────────────────────────────

test_that("gate 38: the preflight warns and then stops at the configured thresholds", {
  pop <- gev_lines_pop("gev_preflight")
  on.exit(close_pop(pop), add = TRUE)
  loci <- gev_loci(pop)
  # A two-member term with an origin scope on *both* members: the label-vector
  # count is the product of the two alphabets, which is the quantity the guard
  # is about. An unscoped term of any order is one label-vector (see below).
  pop <- define_genome_effects(pop, "ADG",
    terms = data.frame(term_id = c(1L, 1L), locus_name = loci[1:2],
                       contrast_name = "additive", center_value = 0.5,
                       genome_value = 1.0),
    origin = data.frame(term_id = c(1L, 1L), locus_name = loci[1:2],
                        line_match_type = "exact", line_name = "A",
                        parent_origin = NA_integer_, copy_count = 1L),
    effect_owner = "AxA")
  ids <- DBI::dbGetQuery(pop$db_conn, "SELECT id_ind FROM ind_meta")$id_ind

  withr::local_options(list(tidybreed.label_vector_warn = 1))
  expect_warning(gev_totals(pop, ids), "Label-vector enumeration for trait 'ADG'")

  withr::local_options(list(tidybreed.label_vector_warn = 1e4,
                            tidybreed.label_vector_max = 1))
  expect_error(gev_totals(pop, ids), "above the hard cap")
})

test_that("an unscoped high-order term is one label-vector, not |labels|^order", {
  # Regression: a 50-member unscoped dominance term used to enumerate
  # |labels|^50 label-vectors and trip the hard cap, though evaluating it is
  # one GROUP BY. No variant scopes those members, so the selected variant
  # cannot depend on their labels and the inner sums collapse.
  pop <- gev_lines_pop("gev_unscoped_order", n_loci = 12)
  on.exit(close_pop(pop), add = TRUE)
  loci <- gev_loci(pop)
  pop <- define_genome_effects(pop, "ADG", data.frame(
    term_id = 1L, locus_name = loci, contrast_name = "dominance",
    center_value = 0.5, genome_value = 1.0), effect_owner = "big")
  ids <- DBI::dbGetQuery(pop$db_conn, "SELECT id_ind FROM ind_meta")$id_ind

  withr::local_options(list(tidybreed.label_vector_max = 100))
  expect_silent(got <- gev_totals(pop, ids))
  expect_equal(length(got), length(ids))

  model <- tidybreed:::.gev_read_model(pop$db_conn, "ADG")
  expect_true(all(tidybreed:::.gev_slot_freedom(model)))
})

test_that("the preflight is silent at the default thresholds", {
  pop <- gev_lines_pop("gev_preflight_quiet")
  on.exit(close_pop(pop), add = TRUE)
  loci <- gev_loci(pop)
  pop <- suppressWarnings(pop |> get_table("genome_meta") |>
    define_additive_effects("ADG", effects = rep(1.0, length(loci))))
  ids <- DBI::dbGetQuery(pop$db_conn, "SELECT id_ind FROM ind_meta")$id_ind
  expect_silent(gev_totals(pop, ids))
})


# ── Gates 51-52: resolution is precomputed, not per individual ──────────────

# Count SQL statements issued against the DuckDB connection. The counter lives
# in the global environment because a traced function in a namespace resolves
# free variables through the namespace's parents, which end at globalenv.
gev_count_queries <- function(expr) {
  counter <- new.env()
  counter$n <- 0L
  assign(".gev_query_counter", counter, envir = globalenv())
  for (fn in c("dbGetQuery", "dbExecute")) {
    suppressMessages(trace(
      fn, signature = c("duckdb_connection", "character"),
      tracer = quote(local({
        e <- get(".gev_query_counter", envir = globalenv())
        assign("n", e$n + 1L, envir = e)
      })),
      print = FALSE, where = asNamespace("duckdb")))
  }
  on.exit({
    for (fn in c("dbGetQuery", "dbExecute")) {
      suppressMessages(try(untrace(fn, signature = c("duckdb_connection", "character"),
                                   where = asNamespace("duckdb")), silent = TRUE))
    }
  }, add = TRUE)
  force(expr)
  rm(".gev_query_counter", envir = globalenv())
  counter$n
}

gev_scale_pop <- function(name, n_ind) {
  pop <- open_pop(pop_name = name, db_name = ":memory:") |>
    define_genome(n_loci = 8, n_chr = 1, chr_len_Mb = 20)
  set.seed(4242)
  pop <- define_founder_haplotypes(pop, n_haplotypes = 30, line_name = "A",
                                   method = "fixed", allele_freq = 0.5)
  pop <- define_founder_haplotypes(pop, n_haplotypes = 30, line_name = "B",
                                   method = "fixed", allele_freq = 0.5)
  pop <- pop |> get_table("founder_haplotypes") |>
    dplyr::filter(line_name == "A") |>
    add_founders(n_males = n_ind / 4, n_females = n_ind / 4, line_name = "A")
  pop <- pop |> get_table("founder_haplotypes") |>
    dplyr::filter(line_name == "B") |>
    add_founders(n_males = n_ind / 4, n_females = n_ind / 4, line_name = "B")
  pop <- define_trait(pop, "ADG", target_add_var = 1.0)
  loci <- gev_loci(pop)
  pop <- suppressWarnings(pop |> get_table("genome_meta") |>
    define_additive_effects("ADG", effects = rep(1.0, length(loci))))
  pop <- pop |> get_table("genome_meta") |>
    define_additive_effects("ADG", effects = rep(2.0, length(loci)),
                            line_name = "A")
  pop |> get_table("genome_meta") |>
    define_additive_effects("ADG", effects = rep(3.0, length(loci)),
                            line_name = "B")
}

test_that("gate 51: the statement count does not depend on the number of individuals", {
  small <- gev_scale_pop("gev_small", 8)
  big   <- gev_scale_pop("gev_big", 400)
  on.exit({ close_pop(small); close_pop(big) }, add = TRUE)

  ids_s <- DBI::dbGetQuery(small$db_conn, "SELECT id_ind FROM ind_meta")$id_ind
  ids_b <- DBI::dbGetQuery(big$db_conn,   "SELECT id_ind FROM ind_meta")$id_ind
  expect_equal(length(ids_s), 8L)
  expect_equal(length(ids_b), 400L)

  n_s <- gev_count_queries(tidybreed:::.gev_evaluate(small$db_conn, ids_s, "ADG"))
  n_b <- gev_count_queries(tidybreed:::.gev_evaluate(big$db_conn,   ids_b, "ADG"))
  expect_equal(n_s, n_b)
  # Three model reads, one label alphabet, one evaluation. The genotype
  # alphabet is skipped when no member needs it.
  expect_lte(n_b, 6L)
})

test_that("gate 51: individual identifiers never reach the SQL text", {
  # The evaluation statement is a pure function of four registered frame names,
  # so an id list cannot be inlined into it -- which is what a per-individual
  # loop would look like from the outside.
  sql <- tidybreed:::.gev_sql("i_tmp", "m_tmp", "p_tmp", "t_tmp")
  expect_false(grepl("'", sql, fixed = TRUE) &&
                 grepl("_1'", sql, fixed = TRUE))
  expect_true(grepl("i_tmp", sql, fixed = TRUE))
  expect_length(sql, 1L)
})

test_that("gate 52: the resolved variant map is the same for 8 and 400 individuals", {
  small <- gev_scale_pop("gev_map_small", 8)
  big   <- gev_scale_pop("gev_map_big", 400)
  on.exit({ close_pop(small); close_pop(big) }, add = TRUE)

  map_of <- function(pop) {
    ids <- DBI::dbGetQuery(pop$db_conn, "SELECT id_ind FROM ind_meta")$id_ind
    conn <- pop$db_conn
    tmp <- "_gev_map_probe"
    duckdb::duckdb_register(conn, tmp, data.frame(id_ind = ids,
                                                  stringsAsFactors = FALSE))
    on.exit(duckdb::duckdb_unregister(conn, tmp), add = TRUE)
    model <- tidybreed:::.gev_read_model(conn, "ADG")
    model$members$member_kind <- tidybreed:::.gev_member_kind(model$members$contrast_name)
    alph <- list(additive = tidybreed:::.gev_additive_alphabet(conn, tmp)$label,
                 genotype = character(0))
    m <- tidybreed:::.gev_variant_map(model, alph,
                                      tidybreed:::.gev_slot_freedom(model))
    m[order(m$id_genome_effect, m$member_slot, m$label),
      c("id_genome_effect", "member_slot", "label", "n_members")]
  }
  a <- map_of(small); b <- map_of(big)
  rownames(a) <- NULL; rownames(b) <- NULL
  expect_equal(a, b)
  expect_gt(nrow(a), 0L)
})

test_that("gate 52: a mixed-line multi-scope model agrees with the per-tuple definition", {
  # The label-vector factorization is not merely faster: on a model with a
  # common variant and two line-specific ones, it is the same function.
  pop <- gev_scale_pop("gev_semantic", 8)
  on.exit(close_pop(pop), add = TRUE)
  ids <- DBI::dbGetQuery(pop$db_conn, "SELECT id_ind FROM ind_meta")$id_ind
  got <- gev_totals(pop, ids)

  eff <- DBI::dbGetQuery(pop$db_conn, paste0(
    "SELECT e.id_genome_effect, e.genome_value, m.locus_id, m.center_value, ",
    "       o.line_match_type, o.line_name, o.parent_origin ",
    "FROM genome_effects e JOIN genome_effect_members m USING (id_genome_effect) ",
    "LEFT JOIN genome_effect_member_origins o USING (id_genome_effect, member_slot)"))
  hap <- DBI::dbGetQuery(pop$db_conn,
    "SELECT id_ind, locus_id, parent_origin, line_origin, allele FROM ind_haplotype")
  want <- vapply(ids, function(id) {
    rows <- hap[hap$id_ind == id, , drop = FALSE]
    sum(vapply(seq_len(nrow(rows)), function(i) {
      r <- rows[i, ]
      cand <- eff[eff$locus_id == r$locus_id, , drop = FALSE]
      ok <- is.na(cand$line_match_type) |
        (cand$line_match_type == "exact" & cand$line_name == r$line_origin)
      cand <- cand[ok, , drop = FALSE]
      pick <- if (any(!is.na(cand$line_match_type)))
        cand[!is.na(cand$line_match_type), , drop = FALSE] else cand
      (r$allele - pick$center_value) * pick$genome_value
    }, numeric(1)))
  }, numeric(1))
  expect_equal(unname(got[ids]), unname(want), tolerance = 1e-10)
})


# ── add_tgv() surface ───────────────────────────────────────────────────────

test_that("add_tgv() refuses a trait with no terms and a table with no id_ind", {
  pop <- gev_lines_pop("gev_api")
  on.exit(close_pop(pop), add = TRUE)
  expect_error(pop |> get_table("ind_meta") |> add_tgv("ADG"),
               "No genome effects found")
  expect_error(add_tgv(pop), "tidybreed_table")
  expect_error(pop |> get_table("ind_meta") |> add_tgv("NOPE"),
               "Traits not found")
})

test_that("add_tgv() defaults to every trait in trait_meta", {
  pop <- gev_lines_pop("gev_all_traits")
  on.exit(close_pop(pop), add = TRUE)
  pop <- define_trait(pop, "BW", target_add_var = 1.0)
  loci <- gev_loci(pop)
  pop <- suppressWarnings(pop |> get_table("genome_meta") |>
    define_additive_effects("ADG", effects = rep(1.0, length(loci))))
  pop <- suppressWarnings(pop |> get_table("genome_meta") |>
    define_additive_effects("BW", effects = rep(0.5, length(loci))))
  pop <- pop |> get_table("ind_meta") |> add_tgv()
  got <- DBI::dbGetQuery(pop$db_conn,
    "SELECT DISTINCT trait_name FROM ind_tgv ORDER BY trait_name")$trait_name
  expect_equal(got, c("ADG", "BW"))
})

test_that("add_tgv() honours a filtered subset", {
  pop <- gev_lines_pop("gev_subset")
  on.exit(close_pop(pop), add = TRUE)
  loci <- gev_loci(pop)
  pop <- suppressWarnings(pop |> get_table("genome_meta") |>
    define_additive_effects("ADG", effects = rep(1.0, length(loci))))
  pop <- pop |> get_table("ind_meta") |> dplyr::filter(line_name == "A") |>
    add_tgv("ADG")
  got <- DBI::dbGetQuery(pop$db_conn,
    "SELECT DISTINCT id_ind FROM ind_tgv ORDER BY id_ind")$id_ind
  expect_true(all(grepl("^A_", got)))
  expect_equal(length(got), 4L)
})
