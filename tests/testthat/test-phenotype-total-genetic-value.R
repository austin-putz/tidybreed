# Step 3a (plans/import_qtl_effect_methods.md §6A, precondition P2, gates
# PH1-PH7): every phenotype path reads the total genetic value
# (ind_tgv_total) by default, or the components phenotype_components lists.
# Before 0.74.0 they read the additive breeding value only, so a dominance or
# A x A term never reached a phenotype.

# 40 founders (20 M, 20 F), 12 loci; trait T with generated additive effects
# at loci 1-8 (target 1) and, with `mixed = TRUE`, user-owner dominance,
# A x A and indicator terms (the same fixture as test-tgv-consolidation.R).
ph_pop <- function(name, mixed = FALSE, n_males = 20L, n_females = 20L,
                   db_name = ":memory:") {
  set.seed(3401)
  pop <- open_pop(pop_name = name, db_name = db_name) |>
    define_genome(n_loci = 12, n_chr = 1, chr_len_Mb = 50) |>
    define_founder_haplotypes(n_haplotypes = 60, method = "beta")
  pop <- pop |> get_table("founder_haplotypes") |>
    add_founders(n_males = n_males, n_females = n_females, line_name = "A")
  pop <- with_additive_target(pop, "T", 1)
  pop <- suppressMessages(pop |> get_table("genome_meta") |>
    dplyr::filter(locus_id <= 8L) |>
    define_additive_effects("T", warn_bounds = NULL))
  if (!mixed) return(pop)
  loci <- pop |> get_table("genome_meta") |> dplyr::collect() |>
    dplyr::arrange(.data$locus_id) |> dplyr::pull("locus_name")
  p <- extract_allele_freq(get_table(pop, "founder_haplotypes"))
  dom <- loci[2:4]
  pop <- define_genome_effect_terms(pop, "T",
    ad_terms(dom, a = 0, d = c(1.6, -1.4, 1.3),
             p = p$allele_freq[match(dom, p$locus_name)],
             coding = "cockerham", report = FALSE),
    effect_owner = "dom")
  pop <- define_genome_effect_terms(pop, "T", data.frame(
    term_id = c(1L, 1L), locus_name = loci[9:10], contrast_name = "additive",
    center_value = 0.5, genome_value = 0.7), effect_owner = "epi")
  define_genome_effect_terms(pop, "T",
    genotype_terms(stats::setNames(data.frame(c(0L, 2L)), loci[11]),
                   value = c(-0.5, 0.8)),
    effect_owner = "surf")
}

ph_records <- function(pop, phenotype) {
  DBI::dbGetQuery(pop$db_conn, paste0(
    "SELECT id_ind, pheno_value, residual_value FROM ind_phenotype ",
    "WHERE phenotype_name = '", phenotype, "' ORDER BY id_ind"))
}

ph_total <- function(pop, trait) {
  tot <- dplyr::collect(get_table(pop, "ind_tgv_total"))
  tot <- tot[tot$trait_name == trait, ]
  stats::setNames(tot$tgv_total, tot$id_ind)
}


test_that("PH1: a simple phenotype is mean + total genetic value + residual", {
  # Additive-only: the total is the breeding value, bit for bit.
  pop <- ph_pop("ph1_add")
  on.exit(close_pop(pop), add = TRUE)
  pop <- define_phenotype(pop, "T", mean = 10, residual_var = 1)
  set.seed(1)
  pop <- suppressMessages(pop |> get_table("ind_meta") |> add_phenotype("T"))
  r <- ph_records(pop, "T")
  a <- tgv_additive(pop, "T")
  g <- a$tgv_value[match(r$id_ind, a$id_ind)]
  expect_identical(r$pheno_value, 10 + 0 + 0 + g + r$residual_value)

  # With dominance, A x A and an indicator: the total, which is not the
  # additive value.
  pop2 <- ph_pop("ph1_mixed", mixed = TRUE)
  on.exit(close_pop(pop2), add = TRUE)
  pop2 <- define_phenotype(pop2, "T", mean = 10, residual_var = 1)
  set.seed(1)
  pop2 <- suppressMessages(pop2 |> get_table("ind_meta") |> add_phenotype("T"))
  r2 <- ph_records(pop2, "T")
  tot <- ph_total(pop2, "T")[r2$id_ind]
  expect_identical(r2$pheno_value, 10 + 0 + 0 + unname(tot) + r2$residual_value)
  a2 <- tgv_additive(pop2, "T")
  expect_false(isTRUE(all.equal(unname(tot), a2$tgv_value[match(r2$id_ind, a2$id_ind)])))
})

test_that("PH1/PH3: composite contributors read the total by default, or the listed components", {
  pop <- ph_pop("ph3_comp", mixed = TRUE)
  on.exit(close_pop(pop), add = TRUE)
  pop <- define_trait(pop, "Z")    # a component trait with no terms of its own kind
  pop <- suppressMessages(pop |> get_table("genome_meta") |>
    dplyr::filter(locus_id == 12L) |>
    define_additive_effects("Z", G = 0.5, warn_bounds = NULL))

  comps <- function(names) tibble::tribble(
    ~source_trait_name, ~contributor_type, ~weight, ~component_names,
    "T",                "self",            1,       names[1],
    "Z",                "self",            2,       names[2])
  pop <- define_phenotype(pop, "P_tot", residual_var = 1,
                          components = comps(c(NA, NA)))
  pop <- define_phenotype(pop, "P_add", residual_var = 1,
                          components = comps(c("additive", "dominance")))
  pop <- define_phenotype(pop, "P_ad", residual_var = 1,
                          components = comps(c("additive,dominance", "additive")))

  stored <- DBI::dbGetQuery(pop$db_conn, paste0(
    "SELECT phenotype_name, source_trait_name, component_names ",
    "FROM phenotype_components ORDER BY phenotype_name, source_trait_name"))
  expect_equal(stored$component_names,
               c("additive,dominance", "additive", "additive", "dominance",
                 "total", "total"))

  set.seed(2)
  pop <- suppressMessages(pop |> get_table("ind_meta") |>
    add_phenotype(c("P_tot", "P_add", "P_ad")))

  tg <- dplyr::collect(get_table(pop, "ind_tgv"))
  comp_of <- function(trait, cn, ids) {
    v <- vapply(ids, function(i) sum(tg$tgv_value[tg$trait_name == trait &
                                                    tg$id_ind == i &
                                                    tg$component_name %in% cn]),
                numeric(1))
    unname(v)
  }
  for (ph in c("P_tot", "P_add", "P_ad")) {
    r <- ph_records(pop, ph)
    expect_equal(nrow(r), 40L)
    g <- switch(ph,
      P_tot = unname(ph_total(pop, "T")[r$id_ind]) + 2 * unname(ph_total(pop, "Z")[r$id_ind]),
      # Z has no dominance terms: the listed component contributes 0.
      P_add = comp_of("T", "additive", r$id_ind) + 2 * 0,
      P_ad  = comp_of("T", c("additive", "dominance"), r$id_ind) +
              2 * comp_of("Z", "additive", r$id_ind))
    expect_equal(r$pheno_value - r$residual_value, g, tolerance = 1e-12,
                 label = ph)
  }
})

test_that("PH3: define_phenotype() refuses a bad component_names value", {
  pop <- ph_pop("ph3_bad")
  on.exit(close_pop(pop), add = TRUE)
  bad <- function(x) tibble::tibble(source_trait_name = "T",
                                    contributor_type = "self",
                                    component_names = x)
  expect_error(define_phenotype(pop, "P", residual_var = 1, components = bad("bogus")),
               "unknown component\\(s\\) 'bogus'")
  expect_error(define_phenotype(pop, "P", residual_var = 1,
                                components = bad("total,additive")),
               "cannot be combined")
  expect_error(define_phenotype(pop, "P", residual_var = 1,
                                components = bad("additive,additive")),
               "lists a component twice")
  expect_error(define_phenotype(pop, "P", residual_var = 1,
                                components = bad("order1_additive")),
               "unknown component")
  # Refused before anything is written: no half-defined phenotype is left.
  expect_equal(nrow(dplyr::collect(get_table(pop, "phenotype_components"))), 0L)
  expect_equal(nrow(dplyr::collect(get_table(pop, "phenotype_meta"))), 0L)

  # With overwrite = TRUE a refused call keeps the existing definition.
  pop <- define_phenotype(pop, "P", residual_var = 1, components = bad("additive"))
  expect_error(define_phenotype(pop, "P", residual_var = 1, overwrite = TRUE,
                                components = bad("bogus")), "unknown component")
  expect_equal(dplyr::collect(get_table(pop, "phenotype_components"))$component_names,
               "additive")
})

test_that("PH3: a failing composite add_phenotype() leaves the record tables unchanged (D7)", {
  pop <- ph_pop("ph3_d7", mixed = TRUE)
  on.exit(close_pop(pop), add = TRUE)
  pop <- define_phenotype(pop, "WW", residual_var = 1,
    missing_component_action = "error",
    components = tibble::tribble(
      ~source_trait_name, ~contributor_type,
      "T",                "self",
      "T",                "dam"))
  # Founders have no dam: the missing component is an error.
  expect_error(suppressMessages(pop |> get_table("ind_meta") |>
                                  add_phenotype("WW")),
               "missing components")
  expect_equal(nrow(dplyr::collect(get_table(pop, "ind_phenotype"))), 0L)
  expect_equal(nrow(dplyr::collect(get_table(pop, "phenotype_random_effects"))), 0L)
})

test_that("PH4: phenotype_meta.mean is an intercept, not the realised base mean", {
  pop <- ph_pop("ph4_mean", mixed = TRUE)
  on.exit(close_pop(pop), add = TRUE)
  pop <- define_phenotype(pop, "T", mean = 100, residual_var = 1)
  set.seed(4)
  pop <- suppressMessages(pop |> get_table("ind_meta") |> add_phenotype("T"))
  r   <- ph_records(pop, "T")
  tot <- unname(ph_total(pop, "T")[r$id_ind])
  # The raw sum: nothing is added to centre the genetic value.
  expect_equal(mean(r$pheno_value) - 100,
               mean(tot) + mean(r$residual_value), tolerance = 1e-12)
  # This base's genetic mean is not zero (functional indicator surface, A x A
  # centred at 0.5, finite sample), so `mean` is not the realised mean.
  expect_gt(abs(mean(tot)), 1e-3)
})

test_that("PH5: the total and the phenotype built from it are bit-identical at 1 and 8 threads", {
  run <- function(name, threads) {
    pop <- ph_pop(name, mixed = TRUE, n_males = 200L, n_females = 200L)
    on.exit(close_pop(pop), add = TRUE)
    DBI::dbExecute(pop$db_conn, paste0("SET threads = ", threads))
    pop <- define_phenotype(pop, "T", mean = 10, residual_var = 1)
    set.seed(5)
    pop <- suppressMessages(pop |> get_table("ind_meta") |> add_phenotype("T"))
    list(total = DBI::dbGetQuery(pop$db_conn, paste0(
           "SELECT id_ind, tgv_total FROM ind_tgv_total ORDER BY id_ind")),
         pheno = ph_records(pop, "T"))
  }
  one <- run("ph5_t1", 1L)
  expect_equal(nrow(one$total), 400L)
  expect_identical(run("ph5_t8", 8L), one)
})

test_that("PH6: a trait with only user-owner terms records phenotypes", {
  pop <- ph_pop("ph6")
  on.exit(close_pop(pop), add = TRUE)
  pop <- define_trait(pop, "H")
  loci <- pop |> get_table("genome_meta") |> dplyr::collect() |>
    dplyr::arrange(.data$locus_id) |> dplyr::pull("locus_name")
  pop <- define_genome_effect_terms(pop, "H", data.frame(
    locus_name = loci[1:5], contrast_name = "additive", genome_value = 0.4),
    base_tbl = get_table(pop, "founder_haplotypes"), effect_owner = "hand")
  pop <- define_phenotype(pop, "H", residual_var = 1)
  pop <- suppressMessages(pop |> get_table("ind_meta") |> add_phenotype("H"))
  expect_equal(nrow(ph_records(pop, "H")), 40L)

  # A trait with no terms at all still errors.
  pop <- define_trait(pop, "EMPTY")
  pop <- define_phenotype(pop, "EMPTY", residual_var = 1)
  expect_error(pop |> get_table("ind_meta") |> add_phenotype("EMPTY"),
               "No genome effects found for phenotype 'EMPTY'")
})

# PH7, the parts that need no owner rule (that is step 3b). The rule: the
# threshold sums the stored population-wide diagonals of the kinds of terms
# the model has, and refuses kinds no target describes.
test_that("PH7: the prevalence threshold uses the active blocks only", {
  pop <- ph_pop("ph7_blocks")
  on.exit(close_pop(pop), add = TRUE)
  loci <- pop |> get_table("genome_meta") |> dplyr::collect() |>
    dplyr::arrange(.data$locus_id) |> dplyr::pull("locus_name")
  p <- extract_allele_freq(get_table(pop, "founder_haplotypes"))

  # Additive only: the additive target alone.
  expect_equal(.ap_prevalence_genetic_var(pop, "T"), 1)

  # Stored dominance and A x A targets, but no terms of either kind yet: the
  # threshold still uses the additive target only.
  pop <- suppressMessages(define_effect_cov_matrix(pop, "dominance", 0.3,
                                                   trait_name = "T"))
  pop <- suppressMessages(define_effect_cov_matrix(pop, "additive_by_additive",
                                                   0.2, trait_name = "T"))
  expect_equal(.ap_prevalence_genetic_var(pop, "T"), 1)

  # Dominance terms (planted under the reserved owner, test-only, as a step-5
  # generator would write them): additive + dominance.
  dom <- loci[2:4]
  pop <- .ge_write_terms(pop, "T",
    ad_terms(dom, a = 0, d = c(0.6, -0.4, 0.3),
             p = p$allele_freq[match(dom, p$locus_name)],
             coding = "cockerham", report = FALSE),
    effect_owner = GE_GENERATED_OWNER, allow_reserved_owner = TRUE)
  expect_equal(.ap_prevalence_genetic_var(pop, "T"), 1 + 0.3)

  # An A x A pair: all three.
  pop <- .ge_write_terms(pop, "T", data.frame(
    term_id = c(1L, 1L), locus_name = loci[9:10], contrast_name = "additive",
    center_value = 0.5, genome_value = 0.7),
    effect_owner = GE_GENERATED_OWNER, allow_reserved_owner = TRUE)
  expect_equal(.ap_prevalence_genetic_var(pop, "T"), 1 + 0.3 + 0.2)

  # End to end: the phenotype runs, and its threshold uses that sum.
  pop <- define_phenotype(pop, "T", type = "categorical", prevalence = 0.2,
                          residual_var = 1, store_liability = TRUE)
  set.seed(7)
  pop <- suppressMessages(pop |> get_table("ind_meta") |> add_phenotype("T"))
  ph <- DBI::dbGetQuery(pop$db_conn,
    "SELECT liability_value, pheno_value FROM ind_phenotype")
  thr <- stats::qnorm(0.8) * sqrt(1.5 + 1)
  expect_identical(ph$pheno_value, as.numeric(ifelse(ph$liability_value > thr, 2, 1)))
})

test_that("PH7: kinds no stored target describes are refused, naming thresholds", {
  pop <- ph_pop("ph7_refuse", mixed = TRUE)
  on.exit(close_pop(pop), add = TRUE)
  pop <- define_phenotype(pop, "T", type = "categorical", prevalence = 0.1,
                          residual_var = 1)
  # The mixed model has an indicator surface: no target can describe it.
  seed_before <- .Random.seed
  expect_error(pop |> get_table("ind_meta") |> add_phenotype("T"),
               "another kind \\(indicator\\).*thresholds")
  expect_identical(.Random.seed, seed_before)
  expect_equal(nrow(dplyr::collect(get_table(pop, "ind_phenotype"))), 0L)
  expect_equal(nrow(dplyr::collect(get_table(pop, "ind_tgv"))), 0L)

  # Explicit thresholds work on the same model.
  pop <- define_phenotype(pop, "T", type = "categorical", thresholds = 0.5,
                          residual_var = 1, overwrite = TRUE)
  pop <- suppressMessages(pop |> get_table("ind_meta") |> add_phenotype("T"))
  expect_equal(nrow(ph_records(pop, "T")), 40L)

  # Dominance terms with no stored dominance target are refused too.
  pop2 <- ph_pop("ph7_nodomtarget")
  on.exit(close_pop(pop2), add = TRUE)
  loci <- pop2 |> get_table("genome_meta") |> dplyr::collect() |>
    dplyr::arrange(.data$locus_id) |> dplyr::pull("locus_name")
  pop2 <- define_genome_effect_terms(pop2, "T",
    ad_terms(loci[2], a = 0, d = 0.5, p = 0.3, coding = "cockerham",
             report = FALSE), effect_owner = "dom")
  expect_error(.ap_prevalence_genetic_var(pop2, "T"),
               "no population-wide row for 'dominance'.*thresholds")
})
