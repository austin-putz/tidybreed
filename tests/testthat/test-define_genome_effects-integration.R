# define_genome_effects() end to end: the 5c gates of
# plans/import_qtl_effect_methods_phase_5_plan.md (C13, prevalence,
# remove_generated_effects(), C19, C20). The generator's own gates are in
# test-define_genome_effects.R. Expectations come from the individuals'
# records and dosages, or from exact identities of the fixture; the
# extractor is the measuring instrument, never the expected value's source
# for the quantity it measures.

dgi_pop <- function(name, n_loci = 40, n_ind = 200, n_hap = 300,
                    traits = c("T1", "T2", "T3")) {
  set.seed(8101)
  suppressMessages(suppressWarnings({
    p <- open_pop(pop_name = name, db_name = ":memory:") |>
      define_genome(n_loci = n_loci, n_chr = 2, chr_len_Mb = 100) |>
      define_founder_haplotypes(n_haplotypes = n_hap)
    p <- p |> get_table("founder_haplotypes") |>
      add_founders(n_males = n_ind / 2, n_females = n_ind / 2, line_name = "A")
    for (t in traits) p <- define_trait(p, t)
    p
  }))
}

quiet  <- function(expr) suppressWarnings(suppressMessages(expr))
all_gm <- function(pop) get_table(pop, "genome_meta")
gm_to  <- function(pop, n) get_table(pop, "genome_meta") |> dplyr::filter(locus_id <= n)
inds   <- function(pop) get_table(pop, "ind_meta")
dgi_ids <- function(pop) {
  DBI::dbGetQuery(pop$db_conn, "SELECT id_ind FROM ind_meta ORDER BY id_ind")$id_ind
}
nm <- function(M, t) { dimnames(M) <- list(t, t); M }

# Assign dosages: `G` is n x m for the individuals `ids` at loci `loci`.
dgi_set_genotypes <- function(pop, ids, loci, G) {
  df <- data.frame(id_ind = rep(ids, length(loci)),
                   locus_name = rep(loci, each = length(ids)),
                   g = as.vector(G), stringsAsFactors = FALSE)
  duckdb::duckdb_register(pop$db_conn, "dgi_set_g", df)
  on.exit(duckdb::duckdb_unregister(pop$db_conn, "dgi_set_g"))
  DBI::dbExecute(pop$db_conn, paste0(
    "UPDATE ind_haplotype h SET allele = CASE WHEN h.parent_origin = 1 ",
    "THEN CAST(s.g >= 1 AS INTEGER) ELSE CAST(s.g = 2 AS INTEGER) END ",
    "FROM dgi_set_g s JOIN genome_meta m USING (locus_name) ",
    "WHERE h.id_ind = s.id_ind AND h.locus_id = m.locus_id"))
  invisible(pop)
}

# An exact HWE + LE cohort: the full factorial of per-locus genotype counts
# (1, 2, 1) at p = 1/2 and (9, 6, 1) at p = 1/4. Every function of one locus
# is uncorrelated with every function of the others, and each locus is in
# exact HWE, so every realised block is n / (n - 1) times its genic value.
HWE_COUNTS <- list(half = c(1, 2, 1), quarter = c(9, 6, 1))
hwe_le_dosages <- function(kinds) {
  per <- lapply(kinds, function(k) rep(0:2, HWE_COUNTS[[k]]))
  as.matrix(expand.grid(per, KEEP.OUT.ATTRS = FALSE))
}
hwe_le_pop <- function(name, traits = "T1") {
  kinds <- c("half", "quarter", "half")
  X <- hwe_le_dosages(kinds)
  pop <- dgi_pop(name, n_loci = length(kinds), n_ind = nrow(X), n_hap = 60,
                 traits = traits)
  dgi_set_genotypes(pop, dgi_ids(pop), paste0("Locus_", seq_along(kinds)), X)
  colnames(X) <- paste0("Locus_", seq_along(kinds))
  list(pop = pop, X = X, n = nrow(X))
}

egv <- function(res, effect, t1, t2 = t1) {
  v <- res$cov_value[res$effect_name == effect & res$trait_name_1 == t1 &
                       res$trait_name_2 == t2]
  if (length(v) == 0L) NA_real_ else v
}

tgv_by <- function(pop, trait, ids, component = NULL) {
  q <- if (is.null(component)) {
    paste0("SELECT id_ind, tgv_total AS v FROM ind_tgv_total WHERE trait_name = '",
           trait, "'")
  } else {
    paste0("SELECT id_ind, tgv_value AS v FROM ind_tgv WHERE trait_name = '",
           trait, "' AND component_name = '", component, "'")
  }
  r <- DBI::dbGetQuery(pop$db_conn, q)
  out <- numeric(length(ids))
  out[match(r$id_ind, ids)] <- r$v
  out
}


# -- C13: phenotypes carry the dominance variance --------------------------

test_that("C13: phenotype minus residual is the total genetic value; its variance is V_A + V_D", {
  pop <- dgi_pop("c13", n_ind = 400)
  on.exit(close_pop(pop))
  set.seed(201)
  quiet(define_genome_effects(all_gm(pop), "T1", G_A = 1, G_D = 2,
                              anchor = "realised", base_tbl = inds(pop)))
  pop <- define_phenotype(pop, "T1", mean = 10, residual_var = 1)
  set.seed(202)
  pop <- quiet(inds(pop) |> add_phenotype("T1"))
  ph <- DBI::dbGetQuery(pop$db_conn,
    "SELECT id_ind, pheno_value, residual_value FROM ind_phenotype ORDER BY id_ind")
  ids <- dgi_ids(pop)
  expect_identical(ph$id_ind, ids)
  g <- tgv_by(pop, "T1", ids)

  # Deterministic: the record is mean + total genetic value + residual.
  expect_equal(ph$pheno_value - ph$residual_value - 10, g, tolerance = 1e-12)
  ev <- quiet(extract_genetic_variance(inds(pop), "T1", anchor = "realised"))
  expect_equal(stats::var(g), egv(ev, "total", "T1"), tolerance = 1e-10)
  # The realised calibration on these individuals: blocks are the targets.
  expect_equal(egv(ev, "additive", "T1"), 1, tolerance = 1e-10)
  expect_equal(egv(ev, "dominance", "T1"), 2, tolerance = 1e-10)

  # Coarse: the phenotypic variance less the residual's is near V_A + V_D = 3,
  # far from V_A = 1 (the genetic total also carries LD covariances).
  v_g <- stats::var(ph$pheno_value) - stats::var(ph$residual_value)
  expect_gt(v_g, 2)
  expect_lt(v_g, 4)
})


# -- Prevalence on a generated A + D + A x A model --------------------------

test_that("prevalence (a): at exact HWE + LE the genic blocks are the targets and the realised total is n/(n-1) times their sum", {
  f <- hwe_le_pop("prev_a")
  pop <- f$pop
  on.exit(close_pop(pop))
  set.seed(211)
  quiet(define_genome_effects(all_gm(pop), "T1", G_A = 1, G_D = 0.3, G_AA = 0.2,
                              pairs = data.frame(locus_name_1 = "Locus_1",
                                                 locus_name_2 = "Locus_2"),
                              base_tbl = inds(pop), warn_bounds = NULL))
  # Generation base = the cohort; centres are its exact frequencies.
  cen <- DBI::dbGetQuery(pop$db_conn, paste0(
    "SELECT DISTINCT m.locus_id, m.center_value FROM genome_effect_members m ",
    "ORDER BY 1"))
  expect_equal(cen$center_value, c(0.5, 0.25, 0.5))

  ge <- quiet(extract_genetic_variance(inds(pop), "T1", anchor = "genic",
                                       base_tbl = inds(pop)))
  expect_equal(egv(ge, "additive", "T1"), 1, tolerance = 1e-10)
  expect_equal(egv(ge, "dominance", "T1"), 0.3, tolerance = 1e-10)
  expect_equal(egv(ge, "additive_by_additive", "T1"), 0.2, tolerance = 1e-10)

  n  <- f$n
  re <- quiet(extract_genetic_variance(inds(pop), "T1", anchor = "realised"))
  expect_equal(egv(re, "between_components", "T1"), 0, tolerance = 1e-12)
  expect_equal(egv(re, "total", "T1"), n / (n - 1) * 1.5, tolerance = 1e-10)

  # The threshold sums the stored diagonals: the realised total as a
  # population variance.
  v_thr <- .ap_prevalence_genetic_var(pop, "T1")
  expect_equal(v_thr, 1.5)
  expect_equal(v_thr, (n - 1) / n * egv(re, "total", "T1"), tolerance = 1e-10)

  # And the phenotype is recorded with that cutpoint.
  pop <- define_phenotype(pop, "T1", type = "categorical", prevalence = 0.2,
                          residual_var = 1, store_liability = TRUE)
  set.seed(212)
  pop <- quiet(inds(pop) |> add_phenotype("T1"))
  ph <- DBI::dbGetQuery(pop$db_conn,
    "SELECT pheno_value, liability_value FROM ind_phenotype")
  cut <- stats::qnorm(0.8) * sqrt(1.5 + 1)
  expect_identical(ph$pheno_value == 2, ph$liability_value > cut)
})

test_that("prevalence (b): off equilibrium the realised total departs from the summed targets by between_components and drift", {
  pop <- dgi_pop("prev_b")
  on.exit(close_pop(pop))
  set.seed(221)
  quiet(define_genome_effects(gm_to(pop, 30), "T1", G_A = 1, G_D = 0.3, G_AA = 0.2))
  re <- quiet(extract_genetic_variance(inds(pop), "T1", anchor = "realised"))
  blocks <- c("additive", "dominance", "additive_by_additive")
  b <- vapply(blocks, function(e) egv(re, e, "T1"), numeric(1))
  between <- egv(re, "between_components", "T1")
  # The B3 identity holds on any cohort ...
  expect_equal(sum(b) + between, egv(re, "total", "T1"), tolerance = 1e-10)
  # ... and here the founders (300 sampled haplotypes) are off the genic
  # base: both parts of the departure are visible.
  expect_gt(abs(between), 1e-4)
  expect_gt(max(abs(b - c(1, 0.3, 0.2))), 1e-3)
  expect_equal(.ap_prevalence_genetic_var(pop, "T1"), 1.5)
})

test_that("prevalence (c): a generated D / A x A model cannot get a parent-scoped additive variant", {
  pop <- dgi_pop("prev_c")
  on.exit(close_pop(pop))
  set.seed(231)
  quiet(define_genome_effects(gm_to(pop, 20), "T1", G_A = 1, G_D = 0.3, G_AA = 0.2))
  before <- DBI::dbGetQuery(pop$db_conn,
    "SELECT COUNT(*) AS n FROM genome_effects")$n
  expect_error(quiet(gm_to(pop, 20) |> define_additive_effects("T1", parent_origin = 1)),
               "define_genome_effects")
  expect_identical(DBI::dbGetQuery(pop$db_conn,
    "SELECT COUNT(*) AS n FROM genome_effects")$n, before)
  # Common scope only: the one-parent-scope rule is never tripped.
  expect_equal(.ap_prevalence_genetic_var(pop, "T1"), 1.5)
})


# -- remove_generated_effects() on a real define_genome_effects() model -----

# The kinds of generated term: one member's contrast, or "interaction".
gen_kinds <- function(pop) {
  DBI::dbGetQuery(pop$db_conn, paste0(
    "SELECT DISTINCT CASE WHEN COUNT(*) > 1 THEN 'interaction' ",
    "ELSE MIN(m.contrast_name) END AS kind ",
    "FROM genome_effects e JOIN genome_effect_members m USING (id_genome_effect) ",
    "WHERE e.effect_owner = 'generated' GROUP BY e.id_genome_effect ORDER BY 1"))$kind
}

test_that("remove_generated_effects() removes the whole A + D + A x A model, keeps targets and custom terms", {
  pop <- dgi_pop("rge")
  on.exit(close_pop(pop))
  set.seed(241)
  quiet(define_genome_effects(gm_to(pop, 20), "T1", G_A = 1, G_D = 0.3, G_AA = 0.2))
  pop <- quiet(define_genome_effect_terms(pop, "T1",
    ad_terms("Locus_30", a = 0.5, d = 0.1, p = 0.5), effect_owner = "custom"))
  tvc_before <- DBI::dbGetQuery(pop$db_conn,
    "SELECT * FROM trait_var_comp ORDER BY effect_name, trait_name_1, trait_name_2")
  expect_identical(gen_kinds(pop), c("additive", "dominance", "interaction"))

  pop <- quiet(remove_generated_effects(pop, "T1"))
  left <- DBI::dbGetQuery(pop$db_conn,
    "SELECT effect_owner, COUNT(*) AS n FROM genome_effects GROUP BY 1")
  expect_identical(left$effect_owner, "custom")
  expect_identical(DBI::dbGetQuery(pop$db_conn,
    "SELECT * FROM trait_var_comp ORDER BY effect_name, trait_name_1, trait_name_2"),
    tvc_before)
  # Nothing generated is left, so the additive generator accepts the trait
  # without the non-additive-model refusal. The stored D and A x A targets
  # stay, so it is told to calibrate the additive block alone.
  set.seed(242)
  only_a <- get_table(pop, "trait_var_comp") |>
    dplyr::filter(effect_name == "additive", is.na(line_name))
  expect_no_error(quiet(gm_to(pop, 20) |>
                          define_additive_effects("T1", trait_var_comp_tbl = only_a)))
  expect_identical(gen_kinds(pop), "additive")
})


# -- C19: the genic anchor measured back ------------------------------------

test_that("C19: a genic model's genic blocks are the targets; at exact HWE + LE ind_tgv's additive is the realised block", {
  f <- hwe_le_pop("c19", traits = c("T1", "T2"))
  pop <- f$pop
  on.exit(close_pop(pop))
  GA <- nm(matrix(c(1, 0.3, 0.3, 2), 2), c("T1", "T2"))
  GD <- nm(matrix(c(0.3, 0.1, 0.1, 0.5), 2), c("T1", "T2"))
  set.seed(251)
  quiet(define_genome_effects(all_gm(pop), c("T1", "T2"), G_A = GA, G_D = GD,
                              base_tbl = inds(pop), warn_bounds = NULL))
  ge <- quiet(extract_genetic_variance(inds(pop), c("T1", "T2"), anchor = "genic",
                                       base_tbl = inds(pop)))
  for (e in c("additive", "dominance")) {
    G <- if (e == "additive") GA else GD
    expect_equal(egv(ge, e, "T1"), G[1, 1], tolerance = 1e-10)
    expect_equal(egv(ge, e, "T2"), G[2, 2], tolerance = 1e-10)
    expect_equal(egv(ge, e, "T1", "T2"), G[1, 2], tolerance = 1e-10)
  }

  # The stored HWE-centred additive component is the realised additive
  # contrast exactly when the cohort is at the stored centres in HWE + LE.
  ids <- dgi_ids(pop)
  pop <- quiet(inds(pop) |> add_tgv(c("T1", "T2")))
  A_tgv <- stats::cov(cbind(tgv_by(pop, "T1", ids, "additive"),
                            tgv_by(pop, "T2", ids, "additive")))
  re <- quiet(extract_genetic_variance(inds(pop), c("T1", "T2"), anchor = "realised"))
  expect_equal(A_tgv[1, 1], egv(re, "additive", "T1"), tolerance = 1e-10)
  expect_equal(A_tgv[2, 2], egv(re, "additive", "T2"), tolerance = 1e-10)
  expect_equal(A_tgv[1, 2], egv(re, "additive", "T1", "T2"), tolerance = 1e-10)
  expect_equal(egv(re, "additive", "T1"), f$n / (f$n - 1) * GA[1, 1],
               tolerance = 1e-10)
})

test_that("C19: away from HWE + LE the stored additive component is not the realised block", {
  pop <- dgi_pop("c19b")
  on.exit(close_pop(pop))
  # Half the individuals fully inbred: far from HWE at the stored centres.
  ids <- dgi_ids(pop)
  inb <- ids[seq(1, length(ids), by = 2)]
  duckdb::duckdb_register(pop$db_conn, "c19_ids", data.frame(id_ind = inb))
  DBI::dbExecute(pop$db_conn, paste0(
    "UPDATE ind_haplotype AS h SET allele = s.allele FROM ind_haplotype AS s ",
    "WHERE h.id_ind = s.id_ind AND h.locus_id = s.locus_id AND h.parent_origin = 2 ",
    "AND s.parent_origin = 1 AND h.id_ind IN (SELECT id_ind FROM c19_ids)"))
  duckdb::duckdb_unregister(pop$db_conn, "c19_ids")
  set.seed(261)
  quiet(define_genome_effects(gm_to(pop, 30), "T1", G_A = 1, G_D = 0.5,
                              base_tbl = inds(pop)))
  pop <- quiet(inds(pop) |> add_tgv("T1"))
  v_tgv <- stats::var(tgv_by(pop, "T1", ids, "additive"))
  re <- quiet(extract_genetic_variance(inds(pop), "T1", anchor = "realised"))
  # The realised contrast uses the cohort's own heterozygote regression b,
  # not HWE's q - p: the two additive measures differ.
  expect_gt(abs(v_tgv - egv(re, "additive", "T1")), 1e-3)
})


# -- C20: the realised round trip -------------------------------------------

test_that("C20: a realised model measured on its base returns G_A, G_D and G_AA", {
  pop <- dgi_pop("c20")
  on.exit(close_pop(pop))
  tr <- c("T1", "T2")
  GA  <- nm(matrix(c(1, 0.3, 0.3, 2), 2), tr)
  GD  <- nm(matrix(c(0.3, 0.1, 0.1, 0.5), 2), tr)
  GAA <- nm(matrix(c(0.2, 0.05, 0.05, 0.3), 2), tr)
  set.seed(271)
  quiet(define_genome_effects(gm_to(pop, 30), tr, G_A = GA, G_D = GD, G_AA = GAA,
                              anchor = "realised", base_tbl = inds(pop)))
  re <- quiet(extract_genetic_variance(inds(pop), tr, anchor = "realised"))
  for (e in c("additive", "dominance", "additive_by_additive")) {
    G <- switch(e, additive = GA, dominance = GD, additive_by_additive = GAA)
    expect_equal(egv(re, e, "T1"), G[1, 1], tolerance = 1e-10)
    expect_equal(egv(re, e, "T2"), G[2, 2], tolerance = 1e-10)
    expect_equal(egv(re, e, "T1", "T2"), G[1, 2], tolerance = 1e-10)
  }
  # The rows join the stored targets one to one (B8): three blocks, both
  # orders of every trait pair.
  tvc <- dplyr::collect(get_table(pop, "trait_var_comp"))
  j <- dplyr::inner_join(re, tvc, by = c("effect_name", "trait_name_1", "trait_name_2"),
                         suffix = c("_measured", "_target"))
  expect_equal(nrow(j), 12L)
  expect_equal(j$cov_value_measured, j$cov_value_target, tolerance = 1e-10)
  # warn_bounds on an inbred panel and NULL: test-define_genome_effects.R,
  # "16 / D6".
})
