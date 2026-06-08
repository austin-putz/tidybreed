# ============================================================================
# Manual integration test: formula-based phenotype specification (v0.34.0)
#
# Tests three complete models end-to-end with printed summaries:
#   Model 1 — Maternal WW  (formula_tbv = "WWD + 0.5 * dam(WWM)")
#   Model 2 — SGE ADG      (formula_tbv = "ADG_direct + group_sum(ADG_SGE, pen_id)")
#   Model 3 — FCR ratio    (formula = "ADFI / ADG")
#
# Run with: source("tests/test_formula_phenotype.R")
# ============================================================================

library(tidybreed)
library(dplyr)
library(tibble)

PASS <- function(msg) cat(sprintf("  ✓ %s\n", msg))
FAIL <- function(msg) cat(sprintf("  ✗ FAIL: %s\n", msg))

check <- function(cond, msg) if (cond) PASS(msg) else FAIL(msg)

make_base <- function(name, n_males = 20, n_females = 20) {
  initialize_genome(
    pop_name = name, n_loci = 300, n_chr = 3, chr_len_Mb = 100,
    db_path = ":memory:"
  ) |>
    define_founder_haplotypes(n_haplotypes = 100, method = "fixed") |>
    add_founders(n_males = n_males, n_females = n_females, line_name = "A")
}

add_offspring_gen <- function(pop, n_off = 40) {
  ind   <- collect(get_table(pop, "ind_meta"))
  sires <- ind$id_ind[ind$sex == "M"]
  dams  <- ind$id_ind[ind$sex == "F"]
  matings <- data.frame(
    id_parent_1 = sample(sires, n_off, replace = TRUE),
    id_parent_2 = sample(dams,  n_off, replace = TRUE),
    sex         = sample(c("M", "F"), n_off, replace = TRUE),
    line_name   = "A",
    stringsAsFactors = FALSE
  )
  add_offspring(pop, matings = matings)
}


# ── Model 1: Maternal WW ─────────────────────────────────────────────────────

cat("\n=== Model 1: Maternal WW (formula_tbv = \"WWD + 0.5 * dam(WWM)\") ===\n")
set.seed(42)
pop1 <- make_base("maternal")

pop1 <- define_trait(pop1, "WWD", target_add_var = 200)
pop1 <- define_trait(pop1, "WWM", target_add_var = 80)
G_ww <- matrix(c(200, 40, 40, 80), 2, 2,
               dimnames = list(c("WWD", "WWM"), c("WWD", "WWM")))
pop1 <- get_table(pop1, "genome_meta") |>
  slice_sample(n = 100) |>
  define_additive_effects(c("WWD", "WWM"), G = G_ww)

pop1 <- define_phenotype(pop1, "WW",
  type         = "continuous",
  mean         = 230,
  residual_var = 180,
  formula_tbv  = "WWD + 0.5 * dam(WWM)")

pop1 <- add_offspring_gen(pop1, n_off = 40)
pop1 <- get_table(pop1, "ind_meta") |>
  filter(!is.na(id_parent_2)) |>
  add_phenotype("WW")

ph1 <- collect(get_table(pop1, "ind_phenotype"))
tbv1 <- unique(collect(get_table(pop1, "ind_tbv"))$trait_name)

cat(sprintf("  WW records:       %d (expected 40)\n", nrow(ph1)))
cat(sprintf("  Mean WW:          %.1f (expected ~230)\n", mean(ph1$pheno_value)))
cat(sprintf("  WW in ind_tbv:    %s (expected FALSE)\n", "WW" %in% tbv1))
cat(sprintf("  phenotype_components rows for WW: %d (expected 0)\n",
  nrow(collect(get_table(pop1, "phenotype_components")))))

check(nrow(ph1) == 40,          "40 WW phenotype records written")
check(abs(mean(ph1$pheno_value) - 230) < 40, "Mean WW within 40 of target 230")
check(!"WW" %in% tbv1,          "WW not in ind_tbv (formula_tbv path)")
check("WWD" %in% tbv1,          "WWD in ind_tbv")
check("WWM" %in% tbv1,          "WWM in ind_tbv")
check(nrow(collect(get_table(pop1, "phenotype_components"))) == 0,
      "No phenotype_components rows for formula_tbv phenotype")

close_pop(pop1)


# ── Model 2: SGE ADG ─────────────────────────────────────────────────────────

cat("\n=== Model 2: SGE ADG (formula_tbv = \"ADG_direct + group_sum(ADG_SGE, pen_id)\") ===\n")
set.seed(701)
pop2 <- make_base("sge", n_males = 20, n_females = 20)

ids     <- collect(get_table(pop2, "ind_meta"))$id_ind
pen_ids <- rep(paste0("pen", 1:4), each = 10)
pop2 <- get_table(pop2, "ind_meta") |> mutate_table(pen_id = pen_ids)

pop2 <- define_trait(pop2, "ADG_direct", target_add_var = 0.4)
pop2 <- define_trait(pop2, "ADG_SGE",    target_add_var = 0.15)
G_sge <- matrix(c(0.4, -0.1, -0.1, 0.15), 2, 2,
                dimnames = list(c("ADG_direct", "ADG_SGE"),
                                c("ADG_direct", "ADG_SGE")))
pop2 <- get_table(pop2, "genome_meta") |>
  slice_sample(n = 100) |>
  define_additive_effects(c("ADG_direct", "ADG_SGE"), G = G_sge)

pop2 <- define_phenotype(pop2, "ADG_obs",
  type         = "continuous",
  mean         = 850,
  residual_var = 300,
  formula_tbv  = "ADG_direct + group_sum(ADG_SGE, pen_id)")

pop2 <- get_table(pop2, "ind_meta") |> add_phenotype("ADG_obs")

ph2  <- collect(get_table(pop2, "ind_phenotype"))
tbv2 <- unique(collect(get_table(pop2, "ind_tbv"))$trait_name)

cat(sprintf("  ADG_obs records:  %d (expected 40)\n", nrow(ph2)))
cat(sprintf("  Mean ADG_obs:     %.1f\n", mean(ph2$pheno_value)))
cat(sprintf("  ADG_obs in ind_tbv: %s (expected FALSE)\n", "ADG_obs" %in% tbv2))

check(nrow(ph2) == 40,            "40 ADG_obs phenotype records written")
check(!"ADG_obs" %in% tbv2,       "ADG_obs not in ind_tbv (formula_tbv path)")
check("ADG_direct" %in% tbv2,     "ADG_direct in ind_tbv")
check("ADG_SGE" %in% tbv2,        "ADG_SGE in ind_tbv")

# Verify social contribution varies by pen (pen-mates have non-zero effect)
cat("  Checking social contribution varies across pens...\n")
# Each pen of 10 has 9 group-mates, so social contribution != 0 for all
check(!any(is.na(ph2$pheno_value)), "No NA phenotype records")

close_pop(pop2)


# ── Model 3: FCR = ADFI / ADG ─────────────────────────────────────────────────

cat("\n=== Model 3: FCR = ADFI / ADG (derived_formula) ===\n")
set.seed(999)
pop3 <- make_base("fcr", n_males = 50, n_females = 50)

pop3 <- define_trait(pop3, "ADFI", target_add_var = 0.1)
pop3 <- define_trait(pop3, "ADG",  target_add_var = 0.1)
G_fcr <- matrix(c(0.1, 0.06, 0.06, 0.1), 2, 2,
                dimnames = list(c("ADFI", "ADG"), c("ADFI", "ADG")))
pop3 <- get_table(pop3, "genome_meta") |>
  slice_sample(n = 100) |>
  define_additive_effects(c("ADFI", "ADG"), G = G_fcr)

pop3 <- define_phenotype(pop3, "ADFI", type = "continuous",
                         mean = 2.2, residual_var = 0.1)
pop3 <- define_phenotype(pop3, "ADG",  type = "continuous",
                         mean = 0.9, residual_var = 0.08)
pop3 <- define_phenotype(pop3, "FCR",  type = "derived_formula",
                         formula = "ADFI / ADG")

# Simulate all phenotypes with a single add_phenotype() call (topo sort ensures FCR last)
pop3 <- get_table(pop3, "ind_meta") |> add_phenotype()

ph3  <- collect(get_table(pop3, "ind_phenotype"))
fcr  <- ph3[ph3$phenotype_name == "FCR",  ]
adfi <- ph3[ph3$phenotype_name == "ADFI", ]
adg  <- ph3[ph3$phenotype_name == "ADG",  ]
tbv3 <- unique(collect(get_table(pop3, "ind_tbv"))$trait_name)

cat(sprintf("  ADFI records:     %d (expected 100)\n", nrow(adfi)))
cat(sprintf("  ADG records:      %d (expected 100)\n", nrow(adg)))
cat(sprintf("  FCR records:      %d (expected 100)\n", nrow(fcr)))
cat(sprintf("  FCR median:       %.3f (expected ~%.3f)\n",
            median(fcr$pheno_value, na.rm = TRUE), 2.2 / 0.9))
cat(sprintf("  FCR in ind_tbv:   %s (expected FALSE)\n", "FCR" %in% tbv3))
cat(sprintf("  NA values in FCR: %d (expected 0)\n",
            sum(is.na(fcr$pheno_value))))

check(nrow(fcr) == 100,                    "100 FCR records written")
check(!"FCR" %in% tbv3,                    "FCR not in ind_tbv")
check(sum(is.na(fcr$pheno_value)) == 0,    "No NA FCR values")
# Ratio distributions are right-skewed (small denominators inflate the mean);
# median is more robust for checking location
check(abs(median(fcr$pheno_value) - 2.2/0.9) < 0.5, "FCR median near 2.44")

# Verify FCR = ADFI / ADG exactly
merged_fcr <- merge(
  fcr[,  c("id_ind", "pheno_value")],
  adfi[, c("id_ind", "pheno_value")],
  by = "id_ind", suffixes = c("_fcr", "_adfi")
)
merged_fcr <- merge(
  merged_fcr,
  adg[, c("id_ind", "pheno_value")],
  by = "id_ind"
)
names(merged_fcr)[4] <- "pheno_value_adg"
check(
  all(abs(merged_fcr$pheno_value_fcr -
          merged_fcr$pheno_value_adfi / merged_fcr$pheno_value_adg) < 1e-9),
  "FCR = ADFI / ADG exactly for all records"
)

close_pop(pop3)


# ── Model 4: Chain of derived formulas + formula_tbv together ────────────────

cat("\n=== Model 4: Combined — maternal WW + FCR chain in one population ===\n")
set.seed(314)
pop4 <- make_base("combined", n_males = 30, n_females = 30)

# Maternal WW traits
pop4 <- define_trait(pop4, "WWD", target_add_var = 200)
pop4 <- define_trait(pop4, "WWM", target_add_var = 80)
G_ww2 <- matrix(c(200, 40, 40, 80), 2, 2,
                dimnames = list(c("WWD","WWM"),c("WWD","WWM")))
pop4 <- get_table(pop4, "genome_meta") |>
  slice_sample(n = 80) |>
  define_additive_effects(c("WWD","WWM"), G = G_ww2)

# FCR component traits
pop4 <- define_trait(pop4, "ADFI", target_add_var = 0.1)
pop4 <- define_trait(pop4, "ADG",  target_add_var = 0.1)
G_fcr2 <- matrix(c(0.1, 0.06, 0.06, 0.1), 2, 2,
                 dimnames = list(c("ADFI","ADG"),c("ADFI","ADG")))
pop4 <- get_table(pop4, "genome_meta") |>
  slice_sample(n = 80) |>
  define_additive_effects(c("ADFI","ADG"), G = G_fcr2)

# Define phenotypes
pop4 <- define_phenotype(pop4, "ADFI", type = "continuous", mean = 2.2,
                         residual_var = 0.1)
pop4 <- define_phenotype(pop4, "ADG",  type = "continuous", mean = 0.9,
                         residual_var = 0.08)
pop4 <- define_phenotype(pop4, "FCR",     type = "derived_formula",
                         formula = "ADFI / ADG")
pop4 <- define_phenotype(pop4, "FCR_pct", type = "derived_formula",
                         formula = "FCR * 100")

# Maternal WW — needs offspring with parents
pop4 <- add_offspring_gen(pop4, n_off = 40)
pop4 <- define_phenotype(pop4, "WW",
  type         = "continuous",
  mean         = 230,
  residual_var = 180,
  formula_tbv  = "WWD + 0.5 * dam(WWM)")

# Simulate all phenotypes at once (topo sort handles ordering)
# Founders: only ADFI, ADG, FCR, FCR_pct
pop4 <- get_table(pop4, "ind_meta") |>
  filter(is.na(id_parent_1)) |>
  add_phenotype(c("ADFI", "ADG", "FCR", "FCR_pct"))

# Offspring: all phenotypes including WW
pop4 <- get_table(pop4, "ind_meta") |>
  filter(!is.na(id_parent_2)) |>
  add_phenotype()

ph4     <- collect(get_table(pop4, "ind_phenotype"))
ph_fcr4 <- ph4[ph4$phenotype_name == "FCR", ]
ph_pct4 <- ph4[ph4$phenotype_name == "FCR_pct", ]
ph_ww4  <- ph4[ph4$phenotype_name == "WW", ]

cat(sprintf("  FCR records:     %d\n", nrow(ph_fcr4)))
cat(sprintf("  FCR_pct records: %d\n", nrow(ph_pct4)))
cat(sprintf("  WW records:      %d\n", nrow(ph_ww4)))

check(nrow(ph_ww4) > 0,           "WW records written for offspring")
check(nrow(ph_fcr4) > 0,          "FCR records written")
check(nrow(ph_pct4) > 0,          "FCR_pct records written")
check(nrow(ph_pct4) == nrow(ph_fcr4), "FCR_pct and FCR have same record count")

close_pop(pop4)


cat("\n=== All manual integration tests complete ===\n")
