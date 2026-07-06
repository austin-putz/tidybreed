# Tests for formula-based phenotype specification:
#   Group A  — Schema / define_phenotype validation
#   Group B  — formula_tbv maternal model (WW = WWD + 0.5 * dam(WWM))
#   Group C  — formula_tbv SGE model (ADG_obs = ADG_direct + group_sum(...))
#   Group D  — derived_formula arithmetic (FCR = ADFI / ADG)
#   Group E  — Backward compatibility (existing components path unchanged)

# ── Shared helpers ─────────────────────────────────────────────────────────────

make_formula_pop <- function(pop_name = "fp", n_males = 20, n_females = 20,
                              n_loci = 300) {
  pop <- open_pop(pop_name = pop_name, db_name = ":memory:") |>
    define_genome(n_loci = n_loci, n_chr = 3, chr_len_Mb = 100) |>
    define_founder_haplotypes(n_haplotypes = 100, method = "fixed")
  pop |>
    get_table("founder_haplotypes") |>
    add_founders( n_males = n_males, n_females = n_females, line_name = "A")
}

add_offspring_gen <- function(pop, n_off = 40) {
  ind   <- dplyr::collect(get_table(pop, "ind_meta"))
  sires <- ind$id_ind[ind$sex == "M"]
  dams  <- ind$id_ind[ind$sex == "F"]
  matings <- data.frame(
    id_parent_1 = sample(sires, n_off, replace = TRUE),
    id_parent_2 = sample(dams,  n_off, replace = TRUE),
    sex         = sample(c("M", "F"), n_off, replace = TRUE),
    line_name   = "A",
    stringsAsFactors = FALSE
  )
  pop <- add_offspring(pop, matings = matings)
  pop <- get_table(pop, "ind_meta") |>
    dplyr::filter(!is.na(id_parent_1)) |>
    mutate_table(gen = 1L)
  pop
}

setup_ww_traits <- function(pop) {
  pop <- define_trait(pop, "WWD", target_add_var = 200)
  pop <- define_trait(pop, "WWM", target_add_var = 80)
  G_ww <- matrix(c(200, 40, 40, 80), 2, 2,
                 dimnames = list(c("WWD", "WWM"), c("WWD", "WWM")))
  sel <- pop |> get_table("genome_meta") |> dplyr::collect() |>
    dplyr::slice_sample(n = 100) |> dplyr::pull(locus_name)
  pop <- pop |>
    get_table("genome_meta") |>
    dplyr::filter(locus_name %in% sel) |>
    define_additive_effects(c("WWD", "WWM"), G = G_ww)
  pop
}

setup_sge_traits <- function(pop) {
  pop <- define_trait(pop, "ADG_direct", target_add_var = 0.4)
  pop <- define_trait(pop, "ADG_SGE",    target_add_var = 0.15)
  G_sge <- matrix(c(0.4, -0.1, -0.1, 0.15), 2, 2,
                  dimnames = list(c("ADG_direct", "ADG_SGE"),
                                  c("ADG_direct", "ADG_SGE")))
  sel <- pop |> get_table("genome_meta") |> dplyr::collect() |>
    dplyr::slice_sample(n = 100) |> dplyr::pull(locus_name)
  pop <- pop |>
    get_table("genome_meta") |>
    dplyr::filter(locus_name %in% sel) |>
    define_additive_effects(c("ADG_direct", "ADG_SGE"), G = G_sge)
  # Assign pens
  ids     <- dplyr::collect(get_table(pop, "ind_meta"))$id_ind
  pen_ids <- rep(paste0("pen", seq_len(ceiling(length(ids) / 10))), each = 10)[seq_along(ids)]
  pop <- get_table(pop, "ind_meta") |> mutate_table(
    pen_id = tibble::tibble(id_ind = ids, pen_id = pen_ids)
  )
  pop
}


# ── Group A: Schema / define validation ───────────────────────────────────────

test_that("formula_tbv is stored in phenotype_meta", {
  set.seed(1)
  pop <- make_formula_pop("fA1")
  on.exit(close_pop(pop))
  pop <- setup_ww_traits(pop)

  pop <- define_phenotype(pop, "WW",
    type        = "continuous",
    mean        = 230,
    residual_var = 180,
    formula_tbv = "WWD + 0.5 * dam(WWM)")

  pm <- dplyr::collect(get_table(pop, "phenotype_meta"))
  ww <- pm[pm$phenotype_name == "WW", ]
  expect_equal(nrow(ww), 1L)
  expect_equal(ww$formula_tbv, "WWD + 0.5 * dam(WWM)")
  expect_true(is.na(ww$formula))
})

test_that("formula is stored in phenotype_meta for derived_formula", {
  set.seed(2)
  pop <- make_formula_pop("fA2")
  on.exit(close_pop(pop))
  pop <- define_trait(pop, "ADFI", target_add_var = 0.1)
  pop <- define_trait(pop, "ADG",  target_add_var = 0.1)
  G_fcr <- matrix(c(0.1, 0.06, 0.06, 0.1), 2, 2,
                  dimnames = list(c("ADFI", "ADG"), c("ADFI", "ADG")))
  pop <- get_table(pop, "genome_meta") |>
    dplyr::slice_sample(n = 100) |>
    define_additive_effects(c("ADFI", "ADG"), G = G_fcr)
  pop <- define_phenotype(pop, "ADFI", type = "continuous",
                          mean = 2.2, residual_var = 0.1)
  pop <- define_phenotype(pop, "ADG",  type = "continuous",
                          mean = 0.9, residual_var = 0.08)
  pop <- define_phenotype(pop, "FCR",  type = "derived_formula",
                          formula = "ADFI / ADG")

  pm  <- dplyr::collect(get_table(pop, "phenotype_meta"))
  fcr <- pm[pm$phenotype_name == "FCR", ]
  expect_equal(fcr$formula,   "ADFI / ADG")
  expect_equal(fcr$type, "derived_formula")
  expect_true(is.na(fcr$formula_tbv))
})

test_that("formula_tbv and components together raises an error", {
  set.seed(3)
  pop <- make_formula_pop("fA3")
  on.exit(close_pop(pop))
  pop <- setup_ww_traits(pop)

  expect_error(
    define_phenotype(pop, "WW",
      type        = "continuous",
      mean        = 230,
      residual_var = 180,
      formula_tbv = "WWD + 0.5 * dam(WWM)",
      components  = tibble::tribble(
        ~source_trait_name, ~contributor_type,
        "WWD", "self", "WWM", "dam"
      )),
    "formula_tbv.*OR.*components"
  )
})

test_that("derived_formula with residual_var raises an error", {
  # Validation fires before any QTL or phenotype setup; only need a pop object
  set.seed(4)
  pop <- make_formula_pop("fA4")
  on.exit(close_pop(pop))

  expect_error(
    define_phenotype(pop, "FCR", type = "derived_formula",
                     formula = "ADFI / ADG", residual_var = 0.1),
    "derived_formula.*no independent residual"
  )
})

test_that("unknown trait in formula_tbv raises error with close-match suggestion", {
  set.seed(5)
  pop <- make_formula_pop("fA5")
  on.exit(close_pop(pop))
  # Only WWD is defined; WWM is not — should error when formula references dam(WWM)
  pop <- define_trait(pop, "WWD", target_add_var = 200)

  expect_error(
    define_phenotype(pop, "WW2",
      type        = "continuous",
      mean        = 230,
      residual_var = 180,
      formula_tbv = "WWD + 0.5 * dam(WWM)"),
    "Unknown trait name"
  )
})

test_that("formula_tbv with invalid R syntax raises a parse error", {
  set.seed(6)
  pop <- make_formula_pop("fA6")
  on.exit(close_pop(pop))
  pop <- define_trait(pop, "WWD", target_add_var = 200)

  expect_error(
    define_phenotype(pop, "WW3",
      type        = "continuous",
      mean        = 230,
      residual_var = 180,
      formula_tbv = "WWD + ("),
    "parse"
  )
})

test_that("formula_tbv with scalar constant warns", {
  set.seed(7)
  pop <- make_formula_pop("fA7")
  on.exit(close_pop(pop))
  pop <- setup_ww_traits(pop)

  expect_warning(
    define_phenotype(pop, "WW4",
      type        = "continuous",
      mean        = 230,
      residual_var = 180,
      formula_tbv = "WWD + 0.5 * dam(WWM) - 10"),
    "Scalar arithmetic constant"
  )
})

test_that("derived_formula missing formula raises error", {
  set.seed(8)
  pop <- make_formula_pop("fA8")
  on.exit(close_pop(pop))

  expect_error(
    define_phenotype(pop, "FCR2", type = "derived_formula"),
    "derived_formula.*requires.*formula"
  )
})


# ── Group B: formula_tbv maternal model ───────────────────────────────────────

test_that("formula_tbv maternal WW: correct record count and no WW in ind_tbv", {
  set.seed(101)
  pop <- make_formula_pop("fB1")
  on.exit(close_pop(pop))
  pop <- setup_ww_traits(pop)

  pop <- define_phenotype(pop, "WW",
    type        = "continuous",
    mean        = 230,
    residual_var = 180,
    formula_tbv = "WWD + 0.5 * dam(WWM)")

  pop <- add_offspring_gen(pop, n_off = 40)

  pop <- get_table(pop, "ind_meta") |>
    dplyr::filter(!is.na(id_parent_2)) |>
    add_phenotype("WW")

  ph  <- dplyr::collect(get_table(pop, "ind_phenotype"))
  expect_equal(nrow(ph), 40L)
  tbv_traits <- unique(dplyr::collect(get_table(pop, "ind_tbv"))$trait_name)
  expect_false("WW" %in% tbv_traits)
  expect_true("WWD" %in% tbv_traits)
  expect_true("WWM" %in% tbv_traits)
})

test_that("formula_tbv maternal WW: mean phenotype approximately correct", {
  set.seed(102)
  pop <- make_formula_pop("fB2")
  on.exit(close_pop(pop))
  pop <- setup_ww_traits(pop)

  pop <- define_phenotype(pop, "WW",
    type        = "continuous",
    mean        = 230,
    residual_var = 180,
    formula_tbv = "WWD + 0.5 * dam(WWM)")

  pop <- add_offspring_gen(pop, n_off = 80)
  pop <- get_table(pop, "ind_meta") |>
    dplyr::filter(!is.na(id_parent_2)) |>
    add_phenotype("WW")

  ph <- dplyr::collect(get_table(pop, "ind_phenotype"))
  expect_equal(mean(ph$pheno_value), 230, tolerance = 40)
})

test_that("formula_tbv maternal model produces no phenotype_components rows", {
  set.seed(103)
  pop <- make_formula_pop("fB3")
  on.exit(close_pop(pop))
  pop <- setup_ww_traits(pop)

  pop <- define_phenotype(pop, "WW",
    type        = "continuous",
    mean        = 230,
    residual_var = 180,
    formula_tbv = "WWD + 0.5 * dam(WWM)")

  pc <- dplyr::collect(get_table(pop, "phenotype_components"))
  expect_equal(nrow(pc[pc$phenotype_name == "WW", ]), 0L)
})

test_that("formula_tbv founders (no dam) are excluded with skip warning", {
  set.seed(104)
  pop <- make_formula_pop("fB4")
  on.exit(close_pop(pop))
  pop <- setup_ww_traits(pop)

  pop <- define_phenotype(pop, "WW",
    type        = "continuous",
    mean        = 230,
    residual_var = 180,
    formula_tbv = "WWD + 0.5 * dam(WWM)",
    missing_component_action = "skip")

  # Founders have no parents — dam TBV will be NA → all excluded with warning
  expect_warning(
    get_table(pop, "ind_meta") |> add_phenotype("WW")
  )
})


# ── Group C: formula_tbv SGE model ────────────────────────────────────────────

test_that("formula_tbv SGE basic end-to-end: 40 records, ADG_obs not in ind_tbv", {
  set.seed(201)
  pop <- make_formula_pop("fC1")
  on.exit(close_pop(pop))
  pop <- setup_sge_traits(pop)

  pop <- define_phenotype(pop, "ADG_obs",
    type        = "continuous",
    mean        = 850,
    residual_var = 300,
    formula_tbv = "ADG_direct + group_sum(ADG_SGE, pen_id)")

  pop <- get_table(pop, "ind_meta") |> add_phenotype("ADG_obs")
  ph  <- dplyr::collect(get_table(pop, "ind_phenotype"))
  expect_equal(nrow(ph), 40L)
  tbv_traits <- unique(dplyr::collect(get_table(pop, "ind_tbv"))$trait_name)
  expect_false("ADG_obs" %in% tbv_traits)
})

test_that("formula_tbv group_mean vs group_sum produce different phenotype values", {
  set.seed(202)
  pop_sum  <- make_formula_pop("fC2s")
  pop_mean <- make_formula_pop("fC2m")

  pop_sum  <- setup_sge_traits(pop_sum)
  pop_mean <- setup_sge_traits(pop_mean)

  pop_sum <- define_phenotype(pop_sum, "ADG_obs",
    type = "continuous", mean = 850, residual_var = 300,
    formula_tbv = "ADG_direct + group_sum(ADG_SGE, pen_id)")
  pop_mean <- define_phenotype(pop_mean, "ADG_obs",
    type = "continuous", mean = 850, residual_var = 300,
    formula_tbv = "ADG_direct + group_mean(ADG_SGE, pen_id)")

  set.seed(301); pop_sum  <- get_table(pop_sum,  "ind_meta") |> add_phenotype("ADG_obs", seed = 99L)
  set.seed(301); pop_mean <- get_table(pop_mean, "ind_meta") |> add_phenotype("ADG_obs", seed = 99L)

  ph_sum  <- dplyr::collect(get_table(pop_sum,  "ind_phenotype"))$pheno_value
  ph_mean <- dplyr::collect(get_table(pop_mean, "ind_phenotype"))$pheno_value

  close_pop(pop_sum); close_pop(pop_mean)
  expect_false(isTRUE(all.equal(ph_sum, ph_mean)))
})

test_that("formula_tbv multiple group terms in one formula work", {
  set.seed(203)
  pop <- make_formula_pop("fC3")
  on.exit(close_pop(pop))

  pop <- define_trait(pop, "ADG_direct",  target_add_var = 0.4)
  pop <- define_trait(pop, "SGE_pen",     target_add_var = 0.1)
  pop <- define_trait(pop, "SGE_barn",    target_add_var = 0.05)
  G3 <- diag(c(0.4, 0.1, 0.05))
  dimnames(G3) <- list(c("ADG_direct","SGE_pen","SGE_barn"),
                       c("ADG_direct","SGE_pen","SGE_barn"))
  pop <- get_table(pop, "genome_meta") |>
    dplyr::slice_sample(n = 100) |>
    define_additive_effects(c("ADG_direct", "SGE_pen", "SGE_barn"), G = G3)

  # Assign pen and barn
  ids      <- dplyr::collect(get_table(pop, "ind_meta"))$id_ind
  pen_ids  <- rep(paste0("pen",  seq_len(4)), each = 10)[seq_along(ids)]
  barn_ids <- rep(paste0("barn", seq_len(2)), each = 20)[seq_along(ids)]
  pop <- get_table(pop, "ind_meta") |>
    mutate_table(pen_id  = tibble::tibble(id_ind = ids, pen_id  = pen_ids),
                 barn_id = tibble::tibble(id_ind = ids, barn_id = barn_ids))

  pop <- define_phenotype(pop, "ADG_multi",
    type        = "continuous",
    mean        = 850,
    residual_var = 300,
    formula_tbv = "ADG_direct + group_sum(SGE_pen, pen_id) + group_sum(SGE_barn, barn_id)")

  pop <- get_table(pop, "ind_meta") |> add_phenotype("ADG_multi")
  ph  <- dplyr::collect(get_table(pop, "ind_phenotype"))
  expect_equal(nrow(ph), 40L)
})

test_that("formula_tbv explicit table= argument works", {
  set.seed(204)
  pop <- make_formula_pop("fC4")
  on.exit(close_pop(pop))
  pop <- setup_sge_traits(pop)

  pop <- define_phenotype(pop, "ADG_obs",
    type        = "continuous",
    mean        = 850,
    residual_var = 300,
    formula_tbv = "ADG_direct + group_sum(ADG_SGE, pen_id, table = \"ind_meta\")")

  pop <- get_table(pop, "ind_meta") |> add_phenotype("ADG_obs")
  ph  <- dplyr::collect(get_table(pop, "ind_phenotype"))
  expect_equal(nrow(ph), 40L)
})

test_that("formula_tbv missing group column errors at add_phenotype with clear message", {
  set.seed(205)
  pop <- make_formula_pop("fC5")
  on.exit(close_pop(pop))
  pop <- setup_sge_traits(pop)

  # Define formula with a non-existent column name
  pop <- define_phenotype(pop, "ADG_obs",
    type        = "continuous",
    mean        = 850,
    residual_var = 300,
    formula_tbv = "ADG_direct + group_sum(ADG_SGE, nonexistent_col)")

  expect_error(
    get_table(pop, "ind_meta") |> add_phenotype("ADG_obs"),
    "nonexistent_col.*not found"
  )
})


# ── Group D: derived_formula arithmetic ───────────────────────────────────────

setup_fcr_pop <- function(pop_name, n = 50, seed = 999) {
  set.seed(seed)
  pop <- make_formula_pop(pop_name, n_males = n %/% 2, n_females = n %/% 2)
  pop <- define_trait(pop, "ADFI", target_add_var = 0.1)
  pop <- define_trait(pop, "ADG",  target_add_var = 0.1)
  G_fcr <- matrix(c(0.1, 0.06, 0.06, 0.1), 2, 2,
                  dimnames = list(c("ADFI", "ADG"), c("ADFI", "ADG")))
  sel  <- dplyr::collect(get_table(pop, "genome_meta")) |>
    dplyr::slice_sample(n = 100) |>
    dplyr::pull(locus_name)
  pop <- get_table(pop, "genome_meta") |>
    dplyr::filter(locus_name %in% sel) |>
    define_additive_effects(c("ADFI", "ADG"), G = G_fcr)
  pop <- define_phenotype(pop, "ADFI", type = "continuous",
                          mean = 2.2, residual_var = 0.1)
  pop <- define_phenotype(pop, "ADG",  type = "continuous",
                          mean = 0.9, residual_var = 0.08)
  pop
}

test_that("FCR derived_formula: correct record count and FCR not in ind_tbv", {
  pop <- setup_fcr_pop("fD1")
  on.exit(close_pop(pop))
  pop <- define_phenotype(pop, "FCR", type = "derived_formula",
                          formula = "ADFI / ADG")
  pop <- get_table(pop, "ind_meta") |> add_phenotype()

  ph  <- dplyr::collect(get_table(pop, "ind_phenotype"))
  fcr <- ph[ph$phenotype_name == "FCR", ]
  expect_equal(nrow(fcr), 50L)
  tbv_traits <- unique(dplyr::collect(get_table(pop, "ind_tbv"))$trait_name)
  expect_false("FCR" %in% tbv_traits)
  expect_true("ADFI" %in% tbv_traits)
})

test_that("FCR central tendency is approximately ADFI_mean / ADG_mean", {
  pop <- setup_fcr_pop("fD2")
  on.exit(close_pop(pop))
  pop <- define_phenotype(pop, "FCR", type = "derived_formula",
                          formula = "ADFI / ADG")
  pop <- get_table(pop, "ind_meta") |> add_phenotype()

  ph  <- dplyr::collect(get_table(pop, "ind_phenotype"))
  fcr <- ph[ph$phenotype_name == "FCR", ]
  # Use the MEDIAN, not the mean: FCR = ADFI / ADG is a ratio whose denominator
  # (ADG, mean 0.9) can approach or cross zero for individual animals, making the
  # mean explode from a single outlier. The median is the stable central-tendency
  # statistic for this ratio and tracks ADFI_mean / ADG_mean.
  expect_equal(stats::median(fcr$pheno_value, na.rm = TRUE), 2.2 / 0.9,
               tolerance = 0.4)
})

test_that("topological ordering: FCR_pct evaluated after FCR when add_phenotype called with no name", {
  pop <- setup_fcr_pop("fD3")
  on.exit(close_pop(pop))
  pop <- define_phenotype(pop, "FCR",     type = "derived_formula",
                          formula = "ADFI / ADG")
  pop <- define_phenotype(pop, "FCR_pct", type = "derived_formula",
                          formula = "FCR * 100")

  pop <- get_table(pop, "ind_meta") |> add_phenotype()

  ph      <- dplyr::collect(get_table(pop, "ind_phenotype"))
  fcr     <- ph[ph$phenotype_name == "FCR",     ]
  fcr_pct <- ph[ph$phenotype_name == "FCR_pct", ]
  expect_equal(nrow(fcr_pct), 50L)
  # FCR_pct should be approximately 100 × FCR for matched individuals
  merged <- merge(
    fcr[,     c("id_ind", "pheno_value")],
    fcr_pct[, c("id_ind", "pheno_value")],
    by = "id_ind", suffixes = c("_fcr", "_pct")
  )
  expect_equal(merged$pheno_value_pct, merged$pheno_value_fcr * 100,
               tolerance = 1e-6)
})

test_that("division by zero produces NA with warning", {
  pop <- setup_fcr_pop("fD4", seed = 555)
  on.exit(close_pop(pop))

  # Force ADG = 0 for all individuals by using user_values
  pop <- define_phenotype(pop, "ADG_zero", type = "derived_formula",
                          formula = "ADG - ADG")

  pop <- get_table(pop, "ind_meta") |> add_phenotype()
  pop <- define_phenotype(pop, "BadFCR", type = "derived_formula",
                          formula = "ADFI / ADG_zero")

  expect_warning(
    get_table(pop, "ind_meta") |> add_phenotype("BadFCR"),
    "non-finite"
  )
  ph     <- dplyr::collect(get_table(pop, "ind_phenotype"))
  badfcr <- ph[ph$phenotype_name == "BadFCR", ]
  expect_true(all(is.na(badfcr$pheno_value)))
})

test_that("RFI with scalar coefficients: ADFI - 0.036*ADG - 0.0072*MBW", {
  set.seed(701)
  pop <- make_formula_pop("fD5", n_males = 25, n_females = 25)
  on.exit(close_pop(pop))

  pop <- define_trait(pop, "ADFI", target_add_var = 0.1)
  pop <- define_trait(pop, "ADG",  target_add_var = 0.1)
  pop <- define_trait(pop, "MBW",  target_add_var = 5)
  G3rfi <- diag(c(0.1, 0.1, 5))
  dimnames(G3rfi) <- list(c("ADFI","ADG","MBW"), c("ADFI","ADG","MBW"))
  sel3 <- dplyr::collect(get_table(pop, "genome_meta")) |>
    dplyr::slice_sample(n = 100) |>
    dplyr::pull(locus_name)
  pop <- get_table(pop, "genome_meta") |>
    dplyr::filter(locus_name %in% sel3) |>
    define_additive_effects(c("ADFI", "ADG", "MBW"), G = G3rfi)

  pop <- define_phenotype(pop, "ADFI", type = "continuous",
                          mean = 2.2, residual_var = 0.1)
  pop <- define_phenotype(pop, "ADG",  type = "continuous",
                          mean = 0.9, residual_var = 0.08)
  pop <- define_phenotype(pop, "MBW",  type = "continuous",
                          mean = 50,  residual_var = 5)
  pop <- define_phenotype(pop, "RFI",  type = "derived_formula",
                          formula = "ADFI - 0.036*ADG - 0.0072*MBW")

  pop <- get_table(pop, "ind_meta") |> add_phenotype()

  ph  <- dplyr::collect(get_table(pop, "ind_phenotype"))
  rfi <- ph[ph$phenotype_name == "RFI", ]
  expect_equal(nrow(rfi), 50L)
  expect_false(any(is.infinite(rfi$pheno_value)))
})

test_that("whitelisted math function sqrt() works in derived formula", {
  pop <- setup_fcr_pop("fD6")
  on.exit(close_pop(pop))
  pop <- define_phenotype(pop, "sqrtADG", type = "derived_formula",
                          formula = "sqrt(ADG)")
  pop <- get_table(pop, "ind_meta") |> add_phenotype()

  ph      <- dplyr::collect(get_table(pop, "ind_phenotype"))
  sqrtadg <- ph[ph$phenotype_name == "sqrtADG", ]
  expect_equal(nrow(sqrtadg), 50L)
  adg     <- ph[ph$phenotype_name == "ADG", ]
  merged  <- merge(sqrtadg[, c("id_ind", "pheno_value")],
                   adg[,     c("id_ind", "pheno_value")],
                   by = "id_ind", suffixes = c("_sqrt", "_adg"))
  expect_equal(merged$pheno_value_sqrt, sqrt(merged$pheno_value_adg),
               tolerance = 1e-9)
})


# ── Group E: Backward compatibility ───────────────────────────────────────────

test_that("existing components-based maternal WW still works", {
  set.seed(501)
  pop <- make_formula_pop("fE1")
  on.exit(close_pop(pop))
  pop <- setup_ww_traits(pop)

  pop <- define_phenotype(pop, "WW",
    type         = "continuous",
    mean         = 230,
    residual_var = 180,
    components   = tibble::tribble(
      ~source_trait_name, ~contributor_type,
      "WWD",              "self",
      "WWM",              "dam"
    ))

  pop <- add_offspring_gen(pop, n_off = 40)
  pop <- get_table(pop, "ind_meta") |>
    dplyr::filter(!is.na(id_parent_2)) |>
    add_phenotype("WW")

  ph <- dplyr::collect(get_table(pop, "ind_phenotype"))
  expect_equal(nrow(ph), 40L)
  expect_equal(mean(ph$pheno_value), 230, tolerance = 40)
})

test_that("existing components-based SGE model still works", {
  set.seed(502)
  pop <- make_formula_pop("fE2")
  on.exit(close_pop(pop))
  pop <- setup_sge_traits(pop)

  pop <- define_phenotype(pop, "ADG_obs",
    type        = "continuous",
    mean        = 850,
    residual_var = 300,
    components  = tibble::tribble(
      ~source_trait_name, ~contributor_type, ~group_column,
      "ADG_direct",       "self",            NA_character_,
      "ADG_SGE",          "group",           "pen_id"
    ))

  pop <- get_table(pop, "ind_meta") |> add_phenotype("ADG_obs")
  ph  <- dplyr::collect(get_table(pop, "ind_phenotype"))
  expect_equal(nrow(ph), 40L)
})

test_that("simple phenotype (no formula, no components) still works", {
  set.seed(503)
  pop <- make_formula_pop("fE3")
  on.exit(close_pop(pop))
  pop <- define_trait(pop, "ADG", target_add_var = 0.4)
  sel <- dplyr::collect(get_table(pop, "genome_meta")) |>
    dplyr::slice_sample(n = 100) |>
    dplyr::pull(locus_name)
  pop <- get_table(pop, "genome_meta") |>
    dplyr::filter(locus_name %in% sel) |>
    define_additive_effects("ADG")
  pop <- define_phenotype(pop, "ADG", type = "continuous",
                          mean = 850, residual_var = 300)
  pop <- get_table(pop, "ind_meta") |> add_phenotype("ADG")
  ph  <- dplyr::collect(get_table(pop, "ind_phenotype"))
  expect_equal(nrow(ph), 40L)
})

test_that("old DB without formula_tbv/formula columns gets migrated by ensure_trait_tables", {
  set.seed(504)
  pop <- make_formula_pop("fE4")
  on.exit(close_pop(pop))

  # Manually drop the new columns to simulate a pre-v0.34.0 database
  DBI::dbExecute(pop$db_conn,
    "ALTER TABLE phenotype_meta DROP COLUMN formula_tbv")
  DBI::dbExecute(pop$db_conn,
    "ALTER TABLE phenotype_meta DROP COLUMN formula")

  cols_before <- DBI::dbListFields(pop$db_conn, "phenotype_meta")
  expect_false("formula_tbv" %in% cols_before)
  expect_false("formula"     %in% cols_before)

  # Calling define_trait triggers ensure_trait_tables which should re-add them
  pop <- define_trait(pop, "ADG", target_add_var = 0.4)

  cols_after <- DBI::dbListFields(pop$db_conn, "phenotype_meta")
  expect_true("formula_tbv" %in% cols_after)
  expect_true("formula"     %in% cols_after)
})
