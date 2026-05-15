### test_binary_mortality.R
### Demonstrates 2-category threshold trait simulation.
###
### Setup:
###   - Trait: neonatal mortality (cat_values = c(0,1); 1 = died, 0 = survived)
###   - Target prevalence: 10%  (used to compute the liability threshold)
###   - Additive variance (VA) = 1, Residual variance (VR) = 9
###   - Heritability h² = VA / (VA + VR) = 0.10
###   - Total liability variance = 10
###   - Fixed threshold = target_add_mean + qnorm(0.90) * sqrt(10) ≈ 4.05
###   - store_liability = TRUE: raw liability also written to ind_phenotype
###
### Expected result: observed mortality rate near 10% (sampling variation expected).

devtools::load_all(rprojroot::find_root(rprojroot::has_file("DESCRIPTION")), quiet = TRUE)
library(dplyr)

cat("=== Categorical Trait (Mortality) Simulation Test ===\n\n")

# ---- Parameters -------------------------------------------------------
N_MALES   <- 2000
N_FEMALES <- 2000
N_LOCI    <- 2000
N_QTL     <- 200
VA        <- 1   # additive genetic variance on liability scale
VR        <- 9   # residual variance on liability scale
PREV      <- 0.10  # target prevalence (mortality rate)
SEED      <- 2024

# ---- Expected threshold -----------------------------------------------
# total liability variance
total_var <- VA + VR
# expected threshold with qnorm()
expected_threshold <- stats::qnorm(1 - PREV) * sqrt(total_var)

# message
cat(sprintf("Liability model: mean=0, VA=%.1f, VR=%.1f, total_var=%.1f, h2=%.2f\n",
            VA, VR, total_var, VA / total_var))
cat(sprintf("Fixed threshold on liability scale: %.4f\n", expected_threshold))
cat(sprintf("  (qnorm(%.2f) * sqrt(%.1f) = %.4f * %.4f)\n\n",
            1 - PREV, total_var,
            stats::qnorm(1 - PREV), sqrt(total_var)))

# ---- Build population -------------------------------------------------
set.seed(SEED)

pop <- initialize_genome(
  pop_name     = "mort_test",
  n_loci       = N_LOCI,
  n_chr        = 5,
  chr_len_Mb   = 100,
  n_haplotypes = 400,
  db_path      = ":memory:"
)

pop <- add_founders(pop,
                    n_males   = N_MALES,
                    n_females = N_FEMALES,
                    line_name = "Duroc")

cat(sprintf("Founders: %d total (%d M + %d F)\n\n", N_MALES + N_FEMALES, N_MALES, N_FEMALES))

# ---- Define trait on liability scale ----------------------------------
# VA and VR are set in trait_effect_cov via target_add_var / residual_var.
# prevalence controls the threshold: threshold = qnorm(1 - prevalence) * sqrt(VA + VR)
# (when target_add_mean = 0, which is the default).
pop <- pop |>
  add_trait(
    "mortality",
    trait_type      = "categorical",
    prevalence      = PREV,
    cat_values      = c(0, 1),
    cat_names       = c("Survived", "Died"),
    store_liability = TRUE,
    description     = "Neonatal mortality (0 = survived, 1 = died)",
    target_add_var  = VA,
    residual_var    = VR
  )
sel_mort <- get_table(pop, "genome_meta") |> collect() |>
  slice_sample(n = N_QTL) |> pull(locus_name)
pop <- pop |>
  get_table("genome_meta") |>
  filter(locus_name %in% sel_mort) |>
  define_qtl("mortality") |>
  add_additive_effects("mortality")

# ---- Simulate phenotypes ----------------------------------------------
pop <- pop |>
  get_table("ind_meta") |>
  add_phenotype("mortality")

# ---- Results ----------------------------------------------------------
ph  <- pop |> get_table("ind_phenotype") |> collect()
print(ph)
tbv <- pop |> get_table("ind_tbv")       |> collect()

n_total   <- nrow(ph)
n_dead    <- sum(ph$pheno_value == 1)
obs_rate  <- mean(ph$pheno_value)

cat("=== Phenotype results ===\n")
cat(sprintf("Total animals phenotyped: %d\n", n_total))
cat(sprintf("Animals dead (value=1):   %d\n", n_dead))
cat(sprintf("Animals alive (value=0):  %d\n", n_total - n_dead))
cat(sprintf("Observed mortality rate:  %.4f  (%.2f%%)\n", obs_rate, obs_rate * 100))
cat(sprintf("Target  mortality rate:   %.4f  (%.2f%%)\n", PREV, PREV * 100))
cat(sprintf("Difference:               %+.4f (%+.2f%%)\n",
            obs_rate - PREV, (obs_rate - PREV) * 100))
cat(sprintf("All values are 0 or 1:    %s\n", all(ph$pheno_value %in% c(0, 1))))
cat(sprintf("liability_value written:  %s\n", "liability_value" %in% names(ph)))
cat(sprintf("cat_name written:         %s\n\n", "cat_name" %in% names(ph)))
if ("cat_name" %in% names(ph)) {
  print(table(ph$cat_name))
  cat("\n")
}

cat("=== TBV (liability scale) ===\n")
cat(sprintf("Mean TBV:     % .4f  (expected ~0)\n",  mean(tbv$tbv_value)))
cat(sprintf("Var(TBV):     %.4f  (target VA = %.1f)\n", var(tbv$tbv_value), VA))
cat(sprintf("SD(TBV):      %.4f  (expected ~%.4f)\n",
            sd(tbv$tbv_value), sqrt(VA)))

# ---- Verify threshold is approximately correct -------------------------
cat("\n=== Threshold verification ===\n")
# The threshold applied was: target_add_mean + qnorm(1-prev) * sqrt(VA + VR)
# We can back-check: P(N(0, total_var) > threshold) should equal PREV
p_above <- 1 - stats::pnorm(expected_threshold, mean = 0, sd = sqrt(total_var))
cat(sprintf("P(N(0, %.1f) > %.4f) = %.4f  (should be %.4f)\n",
            total_var, expected_threshold, p_above, PREV))

# ---- Simple assertion -------------------------------------------------
tol <- 0.015  # within 1.5 pp for n=4000
cat(sprintf("\nPass/Fail (tolerance ±%.1f pp): %s\n",
            tol * 100,
            ifelse(abs(obs_rate - PREV) <= tol, "PASS", "FAIL")))

close_pop(pop)
