make_pheno_base_pop <- function(pop_name = "defph") {
  pop <- initialize_genome(
    pop_name          = pop_name,
    n_loci            = 200,
    n_chr             = 2,
    chr_len_Mb        = 100,
    n_haplotypes      = 50,
    db_path           = ":memory:",
    fixed_allele_freq = 0.5
  )
  pop <- add_founders(pop, n_males = 20, n_females = 20, line_name = "A")
  pop
}


test_that("define_phenotype() inserts a row into phenotype_meta", {
  pop <- make_pheno_base_pop("dp_basic")
  on.exit(close_pop(pop))

  pop <- define_trait(pop, "ADG", target_add_var = 1)
  pop <- define_phenotype(pop, "ADG",
                          trait_type   = "continuous",
                          mean         = 500,
                          residual_var = 120)

  row <- DBI::dbGetQuery(pop$db_conn,
    "SELECT * FROM phenotype_meta WHERE phenotype_name = 'ADG'")
  expect_equal(nrow(row), 1L)
  expect_equal(row$trait_type,    "continuous")
  expect_equal(row$mean,          500)
  expect_equal(row$expressed_sex, "both")
  expect_false(row$repeatable)
})


test_that("define_phenotype() residual_var writes to phenotype_residual_cov", {
  pop <- make_pheno_base_pop("dp_resid")
  on.exit(close_pop(pop))

  pop <- define_trait(pop, "ADG", target_add_var = 1)
  pop <- define_phenotype(pop, "ADG", residual_var = 150)

  rcov <- DBI::dbGetQuery(pop$db_conn,
    "SELECT * FROM phenotype_residual_cov WHERE phenotype_name_1 = 'ADG'")
  expect_equal(nrow(rcov), 1L)
  expect_equal(rcov$cov_value, 150)
  expect_true(is.na(rcov$condition_column))
})


test_that("define_phenotype() overwrite = FALSE errors on duplicate", {
  pop <- make_pheno_base_pop("dp_dup")
  on.exit(close_pop(pop))

  pop <- define_trait(pop, "ADG", target_add_var = 1)
  pop <- define_phenotype(pop, "ADG", residual_var = 100)
  expect_error(define_phenotype(pop, "ADG", residual_var = 200), "already exists")
})


test_that("define_phenotype() overwrite = TRUE replaces the row", {
  pop <- make_pheno_base_pop("dp_overwrite")
  on.exit(close_pop(pop))

  pop <- define_trait(pop, "ADG", target_add_var = 1)
  pop <- define_phenotype(pop, "ADG", mean = 100, residual_var = 50)
  pop <- define_phenotype(pop, "ADG", mean = 999, residual_var = 77,
                          overwrite = TRUE)

  row <- DBI::dbGetQuery(pop$db_conn,
    "SELECT mean FROM phenotype_meta WHERE phenotype_name = 'ADG'")
  expect_equal(row$mean, 999)

  rcov <- DBI::dbGetQuery(pop$db_conn,
    "SELECT cov_value FROM phenotype_residual_cov
     WHERE phenotype_name_1 = 'ADG' AND condition_column IS NULL")
  expect_equal(rcov$cov_value, 77)
})


test_that("define_phenotype() categorical with prevalence validates correctly", {
  pop <- make_pheno_base_pop("dp_cat_prev")
  on.exit(close_pop(pop))

  pop <- define_trait(pop, "mort", target_add_var = 1)

  # Must supply thresholds OR prevalence
  expect_error(
    define_phenotype(pop, "mort", trait_type = "categorical", residual_var = 1),
    "thresholds"
  )
  # Cannot supply both
  expect_error(
    define_phenotype(pop, "mort", trait_type = "categorical",
                     prevalence = 0.1, thresholds = c(0), residual_var = 1),
    "not both"
  )

  pop <- define_phenotype(pop, "mort",
                          trait_type      = "categorical",
                          prevalence      = 0.10,
                          cat_values      = c(0, 1),
                          cat_names       = c("Alive", "Dead"),
                          store_liability = TRUE,
                          residual_var    = 9)

  row <- DBI::dbGetQuery(pop$db_conn,
    "SELECT * FROM phenotype_meta WHERE phenotype_name = 'mort'")
  expect_equal(row$trait_type,      "categorical")
  expect_equal(row$prevalence,      0.10)
  expect_equal(row$cat_values,      "0,1")
  expect_equal(row$cat_names,       "Alive,Dead")
  expect_true(row$store_liability)
})


test_that("define_phenotype() categorical with thresholds stores them correctly", {
  pop <- make_pheno_base_pop("dp_cat_thresh")
  on.exit(close_pop(pop))

  pop <- define_trait(pop, "score", target_add_var = 1)
  pop <- define_phenotype(pop, "score",
                          trait_type  = "categorical",
                          thresholds  = c(-1, 0, 1),
                          cat_values  = c(1, 2, 3, 4),
                          cat_names   = c("thin", "fair", "good", "excellent"),
                          residual_var = 1)

  row <- DBI::dbGetQuery(pop$db_conn,
    "SELECT thresholds, cat_values, cat_names FROM phenotype_meta
     WHERE phenotype_name = 'score'")
  expect_equal(row$thresholds, "-1,0,1")
  expect_equal(row$cat_values, "1,2,3,4")
  expect_equal(row$cat_names,  "thin,fair,good,excellent")
})


test_that("define_phenotype() validates cat_values and cat_names lengths", {
  pop <- make_pheno_base_pop("dp_cat_len")
  on.exit(close_pop(pop))

  pop <- define_trait(pop, "x", target_add_var = 1)

  # 1 threshold → 2 categories; cat_values must be length 2
  expect_error(
    define_phenotype(pop, "x", trait_type = "categorical",
                     thresholds = c(0), cat_values = c(1, 2, 3), residual_var = 1),
    "length 2"
  )
  # cat_names must match cat_values count
  expect_error(
    define_phenotype(pop, "x", trait_type = "categorical",
                     thresholds = c(0), cat_values = c(0, 1),
                     cat_names = c("No"), residual_var = 1),
    "length 2"
  )
  # cat_names cannot contain commas
  expect_error(
    define_phenotype(pop, "x", trait_type = "categorical",
                     thresholds = c(0),
                     cat_names = c("Yes, really", "No"), residual_var = 1),
    "commas"
  )
})


test_that("define_phenotype() expressed_sex is stored", {
  pop <- make_pheno_base_pop("dp_sex")
  on.exit(close_pop(pop))

  pop <- define_trait(pop, "milk", target_add_var = 1)
  pop <- define_phenotype(pop, "milk",
                          expressed_sex = "F",
                          residual_var  = 5)

  row <- DBI::dbGetQuery(pop$db_conn,
    "SELECT expressed_sex FROM phenotype_meta WHERE phenotype_name = 'milk'")
  expect_equal(row$expressed_sex, "F")
})


test_that("define_phenotype() stores components in phenotype_components", {
  pop <- make_pheno_base_pop("dp_comp")
  on.exit(close_pop(pop))

  pop <- define_trait(pop, "WWD", target_add_var = 200)
  pop <- define_trait(pop, "WWM", target_add_var = 80)

  pop <- define_phenotype(pop, "WW",
                          trait_type   = "continuous",
                          mean         = 230,
                          residual_var = 180,
                          components   = tibble::tribble(
                            ~source_trait_name, ~contributor_type,
                            "WWD",              "self",
                            "WWM",              "dam"
                          ))

  comp <- DBI::dbGetQuery(pop$db_conn,
    "SELECT * FROM phenotype_components WHERE phenotype_name = 'WW'
     ORDER BY id_phenotype_comp")
  expect_equal(nrow(comp), 2L)
  expect_equal(comp$source_trait_name, c("WWD", "WWM"))
  expect_equal(comp$contributor_type,  c("self", "dam"))
  expect_equal(comp$weight,            c(1, 1))
  expect_equal(comp$weight_type,       c("fixed", "fixed"))

  # phenotype_meta row exists
  pm <- DBI::dbGetQuery(pop$db_conn,
    "SELECT * FROM phenotype_meta WHERE phenotype_name = 'WW'")
  expect_equal(nrow(pm), 1L)
  expect_equal(pm$mean, 230)
})


test_that("add_residual_cov() writes conditional residual rows", {
  pop <- make_pheno_base_pop("dp_cond_resid")
  on.exit(close_pop(pop))

  pop <- define_trait(pop, "BW", target_add_var = 100)
  pop <- define_phenotype(pop, "BW", residual_var = 600)

  pop <- add_residual_cov(pop,
                          phenotype_names  = "BW",
                          cov_matrix       = matrix(400, 1, 1,
                                                    dimnames = list("BW","BW")),
                          condition_column = "sex",
                          condition_level  = "M")
  pop <- add_residual_cov(pop,
                          phenotype_names  = "BW",
                          cov_matrix       = matrix(800, 1, 1,
                                                    dimnames = list("BW","BW")),
                          condition_column = "sex",
                          condition_level  = "F")

  rcov <- DBI::dbGetQuery(pop$db_conn,
    "SELECT * FROM phenotype_residual_cov WHERE phenotype_name_1 = 'BW'
     ORDER BY condition_level")
  expect_equal(nrow(rcov), 3L)  # unconditional + M + F
  conds <- rcov$condition_level[!is.na(rcov$condition_level)]
  expect_true("M" %in% conds && "F" %in% conds)
})


test_that("define_effect_cov_matrix() with effect_name='residual' routes to phenotype_residual_cov", {
  pop <- make_pheno_base_pop("dp_cov_mat_route")
  on.exit(close_pop(pop))

  pop <- define_trait(pop, "ADG", target_add_var = 1)
  pop <- define_trait(pop, "BW",  target_add_var = 2)
  pop <- define_phenotype(pop, "ADG", residual_var = NULL)
  pop <- define_phenotype(pop, "BW",  residual_var = NULL)

  R <- matrix(c(100, 30, 30, 200), 2, 2,
              dimnames = list(c("ADG","BW"), c("ADG","BW")))
  pop <- define_effect_cov_matrix(pop, "residual", R)

  rcov <- DBI::dbGetQuery(pop$db_conn,
    "SELECT * FROM phenotype_residual_cov WHERE condition_column IS NULL
     ORDER BY phenotype_name_1, phenotype_name_2")
  # Should have 4 rows: (ADG,ADG), (ADG,BW), (BW,ADG), (BW,BW)
  expect_equal(nrow(rcov), 4L)
  adg_var <- rcov$cov_value[rcov$phenotype_name_1 == "ADG" &
                             rcov$phenotype_name_2 == "ADG"]
  expect_equal(adg_var, 100)
})


test_that("missing_component_action stored in phenotype_meta, defaults to 'skip'", {
  pop <- make_pheno_base_pop("dp_mca")
  on.exit(close_pop(pop))

  pop <- define_trait(pop, "ADG", target_add_var = 1)

  # Default: skip
  pop <- define_phenotype(pop, "ADG", residual_var = 1)
  row <- DBI::dbGetQuery(pop$db_conn,
    "SELECT missing_component_action FROM phenotype_meta WHERE phenotype_name = 'ADG'")
  expect_equal(row$missing_component_action, "skip")

  # Explicit: error
  pop <- define_trait(pop, "BW", target_add_var = 1)
  pop <- define_phenotype(pop, "BW", residual_var = 1,
                          missing_component_action = "error")
  row2 <- DBI::dbGetQuery(pop$db_conn,
    "SELECT missing_component_action FROM phenotype_meta WHERE phenotype_name = 'BW'")
  expect_equal(row2$missing_component_action, "error")
})
