# The `prevalence` threshold is mean + qnorm(1 - prevalence) * sqrt(V), a
# Gaussian approximation in which V must be the variance of the whole
# liability: the stored genetic target(s), every named random effect and the
# residual. These gates come from the step-3 Codex review
# (plans/import_qtl_effect_methods_phase_3_codex_review.md, findings 3, 4, 9);
# each checks the deterministic cutoff or refusal, not a sampled fraction.

pt_pop <- function(name, n = 30L) {
  set.seed(3401)
  open_pop(pop_name = name, db_name = ":memory:") |>
    define_genome(n_loci = 12, n_chr = 1, chr_len_Mb = 100) |>
    define_founder_haplotypes(n_haplotypes = 100) |>
    get_table("founder_haplotypes") |>
    add_founders(n_males = n, n_females = n, line_name = "A")
}

# sum p q a^2 over the trait's order-one additive terms (one eligible copy)
pt_pq_a2 <- function(pop, trait) {
  x <- DBI::dbGetQuery(pop$db_conn,
    "SELECT m.center_value AS p, e.genome_value AS a
     FROM genome_effect_members m JOIN genome_effects e USING (id_genome_effect)
     WHERE e.trait_name = ?", params = list(trait))
  sum(x$p * (1 - x$p) * x$a^2)
}


test_that("prevalence refuses generated variants for two parent scopes (finding 3)", {
  pop <- pt_pop("pt_parents")
  on.exit(close_pop(pop), add = TRUE)
  pop <- with_additive_target(pop, "T", 1)
  gm  <- get_table(pop, "genome_meta")

  # One parent-only variant alone is described by the target.
  pop <- suppressMessages(gm |> define_additive_effects("T", parent_origin = 1,
                                                        seed = 31, warn_bounds = NULL))
  expect_equal(pt_pq_a2(pop, "T"), 1, tolerance = 1e-8)
  expect_equal(.ap_prevalence_genetic_var(pop, "T"), 1)

  # Reciprocal: each scope calibrated to 1, their combination has variance 2.
  pop <- suppressMessages(gm |> define_additive_effects("T", parent_origin = 2,
                                                        seed = 32, warn_bounds = NULL))
  expect_equal(pt_pq_a2(pop, "T"), 2, tolerance = 1e-8)
  expect_error(.ap_prevalence_genetic_var(pop, "T"),
               "parent_origin 1 only \\+ parent_origin 2 only.*thresholds")

  # Refused in PLAN: nothing written, no draw.
  pop <- define_phenotype(pop, "T", type = "categorical", prevalence = 0.1,
                          residual_var = 1)
  seed_before <- .Random.seed
  expect_error(suppressMessages(get_table(pop, "ind_meta") |> add_phenotype("T")),
               "more than one parent-of-origin scope")
  expect_identical(.Random.seed, seed_before)
  expect_equal(DBI::dbGetQuery(pop$db_conn,
    "SELECT COUNT(*) AS n FROM ind_phenotype")$n, 0)
})

test_that("prevalence refuses a common variant with a parent-only fallback (finding 3)", {
  pop <- pt_pop("pt_fallback")
  on.exit(close_pop(pop), add = TRUE)
  pop <- with_additive_target(pop, "T", 1)
  gm  <- get_table(pop, "genome_meta")
  pop <- suppressMessages(gm |> define_additive_effects("T", seed = 1,
                                                        warn_bounds = NULL))
  expect_warning(pop <- suppressMessages(gm |>
    define_additive_effects("T", parent_origin = 1, seed = 2, warn_bounds = NULL)),
    "differ only in the parent dimension.*remove_generated_effects")
  expect_error(.ap_prevalence_genetic_var(pop, "T"),
               "both parents \\+ parent_origin 1 only")

  # The route the warning names: remove the parent-only variant, and the
  # common variant alone is described by the target again.
  pop <- suppressMessages(remove_generated_effects(pop, "T", parent_origin = 1))
  expect_equal(.ap_prevalence_genetic_var(pop, "T"), 1)
})

test_that("the threshold includes named random-effect variance (finding 4)", {
  # Genetic variance ~0, residual 1, a normal PE effect of 9: the liability is
  # Gaussian with variance 10, so the cutoff is qnorm(0.9) * sqrt(10).
  pop <- pt_pop("pt_random")
  on.exit(close_pop(pop), add = TRUE)
  pop <- define_trait(pop, "T")
  pop <- suppressMessages(get_table(pop, "genome_meta") |>
    define_additive_effects("T", G = 1e-12, seed = 1, warn_bounds = NULL))
  pop <- define_phenotype(pop, "T", type = "categorical", prevalence = 0.1,
                          residual_var = 1, store_liability = TRUE)
  pop <- suppressMessages(define_effect_random(pop, "T", "pe",
                                               source_column = "id_ind",
                                               variance = 9))
  expect_equal(.ap_prevalence_env_var(pop, "T"), 9)
  pop <- suppressMessages(get_table(pop, "ind_meta") |> add_phenotype("T", seed = 15))
  ph  <- DBI::dbGetQuery(pop$db_conn,
    "SELECT liability_value, pheno_value FROM ind_phenotype")
  thr <- stats::qnorm(0.9) * sqrt(1e-12 + 9 + 1)
  expect_identical(ph$pheno_value,
                   as.numeric(ifelse(ph$liability_value > thr, 2, 1)))

  # A uniform effect is centred with the declared variance and counts too.
  pop <- suppressMessages(define_effect_random(pop, "T", "herd",
    source_column = "id_ind", variance = 2, distribution = "uniform"))
  expect_equal(.ap_prevalence_env_var(pop, "T"), 11)
})

test_that("prevalence refuses a gamma effect and conditional residual strata (finding 4)", {
  pop <- pt_pop("pt_refuse_env")
  on.exit(close_pop(pop), add = TRUE)
  pop <- with_additive_target(pop, "T", 1)
  pop <- suppressMessages(get_table(pop, "genome_meta") |>
    define_additive_effects("T", seed = 1, warn_bounds = NULL))
  pop <- define_phenotype(pop, "T", type = "categorical", prevalence = 0.1,
                          residual_var = 1)
  pop <- suppressMessages(define_effect_random(pop, "T", "pen",
    source_column = "id_ind", variance = 1, distribution = "gamma"))
  seed_before <- .Random.seed
  expect_error(suppressMessages(get_table(pop, "ind_meta") |> add_phenotype("T")),
               "'pen' is a gamma effect.*thresholds")
  expect_identical(.Random.seed, seed_before)

  pop2 <- pt_pop("pt_strata")
  on.exit(close_pop(pop2), add = TRUE)
  pop2 <- with_additive_target(pop2, "T", 1)
  pop2 <- suppressMessages(get_table(pop2, "genome_meta") |>
    define_additive_effects("T", seed = 1, warn_bounds = NULL))
  pop2 <- define_phenotype(pop2, "T", type = "categorical", prevalence = 0.1,
                           residual_var = 1)
  for (s in c("M", "F")) {
    pop2 <- define_residual_cov(pop2, "T", matrix(if (s == "M") 1 else 4,
                                                  dimnames = list("T", "T")),
                                condition_column = "sex", condition_level = s)
  }
  expect_error(suppressMessages(get_table(pop2, "ind_meta") |> add_phenotype("T")),
               "conditional strata \\(sex\\).*thresholds")
})

test_that("a liability exactly on a cutpoint stays in the lower category (finding 9)", {
  # A 0/1/2 dosage surface with no residual puts many liabilities exactly on
  # the cutpoint 1; "above the threshold" means strictly above.
  expect_identical(liability_to_categorical(c(0, 1, 1 + 1e-12, 2), 1),
                   c(1L, 1L, 2L, 2L))
  expect_identical(liability_to_categorical(c(-1, 0, 0.5, 1, 3), c(0, 1)),
                   c(1L, 1L, 2L, 2L, 3L))

  pop <- pt_pop("pt_ties")
  on.exit(close_pop(pop), add = TRUE)
  pop <- define_trait(pop, "T")
  pop <- define_genome_effect_terms(pop, "T",
    genotype_terms(data.frame(Locus_1 = 0:2), c(0, 1, 2)))
  pop <- define_phenotype(pop, "T", type = "categorical", thresholds = 1,
                          residual_var = 0, store_liability = TRUE)
  pop <- suppressMessages(get_table(pop, "ind_meta") |> add_phenotype("T", seed = 1))
  ph <- DBI::dbGetQuery(pop$db_conn,
    "SELECT liability_value, pheno_value FROM ind_phenotype")
  expect_gt(sum(ph$liability_value == 1), 0)
  expect_identical(ph$pheno_value,
                   as.numeric(ifelse(ph$liability_value > 1, 2, 1)))
})

test_that("prevalence on a liability with zero variance is an error, not a fraction", {
  pop <- pt_pop("pt_zero")
  on.exit(close_pop(pop), add = TRUE)
  pop <- define_trait(pop, "T")
  pop <- suppressMessages(get_table(pop, "genome_meta") |>
    define_additive_effects("T", G = 0, seed = 1, warn_bounds = NULL))
  pop <- define_phenotype(pop, "T", type = "categorical", prevalence = 0.1,
                          residual_var = 0)
  expect_error(suppressMessages(get_table(pop, "ind_meta") |>
                                  add_phenotype("T", seed = 1)),
               "positive variance.*point mass")
  expect_equal(DBI::dbGetQuery(pop$db_conn,
    "SELECT COUNT(*) AS n FROM ind_phenotype")$n, 0)
})
