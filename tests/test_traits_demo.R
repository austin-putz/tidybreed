# tests/test_traits_demo.R
#
# End-to-end exercise of the 0.2.0 trait + phenotype system.
# Run interactively to confirm the API works under both in-memory and
# file-based DuckDB.

suppressPackageStartupMessages({
  if (requireNamespace("devtools", quietly = TRUE) &&
      file.exists("DESCRIPTION")) {
    devtools::load_all(".", quiet = TRUE)
  } else {
    library(tidybreed)
  }
  library(dplyr)
})

set.seed(20260420)

# ---- 1. Initialise genome + founders + one generation of offspring --------

pop <- initialize_genome(
  pop_name        = "demo",
  n_loci          = 600,
  n_chr           = 6,
  chr_len_Mb      = 100,
  n_haplotypes    = 200,
  db_path         = ":memory:",
  fixed_allele_freq = 0.5
)
pop <- add_founders(pop, n_males = 50, n_females = 50, line_name = "A", gen = 0L)

# Produce one generation of offspring
sires <- get_table(pop, "ind_meta") |> filter(sex == "M", gen == 0L) |>
  collect() |> pull(id_ind)
dams  <- get_table(pop, "ind_meta") |> filter(sex == "F", gen == 0L) |>
  collect() |> pull(id_ind)

matings <- tibble::tibble(
  id_parent_1 = sample(sires, 200, replace = TRUE),
  id_parent_2 = sample(dams,  200, replace = TRUE),
  sex         = sample(c("M", "F"), 200, replace = TRUE),
  line_name   = "A",
  gen         = 1L
)
pop <- add_offspring(pop, matings)

# ---- 2. Two correlated continuous traits ----------------------------------

G <- matrix(c(0.25, 0.10, 0.10, 0.30), 2, 2)
R <- matrix(c(0.75, 0.20, 0.20, 0.70), 2, 2)

pop <- define_trait(pop, "ADG", target_add_var = 0.25, residual_var = 0.75,
                 mean = 850)
pop <- define_trait(pop, "BW",  target_add_var = 0.30, residual_var = 0.70,
                 mean = 450)

sel_qtl <- get_table(pop, "genome_meta") |> collect() |>
  slice_sample(n = 150) |> pull(locus_name)
pop <- pop |>
  get_table("genome_meta") |>
  filter(locus_name %in% sel_qtl) |>
  define_additive_effects(c("ADG", "BW"), G = G, seed = 42)
pop <- define_effect_cov_matrix(pop, "residual", R)

pop <- define_effect_fixed_class(pop, "ADG", "sex",
  source_column = "sex",
  levels = c(M = 30, F = 0))

pop <- pop |>
  get_table("ind_meta") |>
  filter(gen == 1L) |>
  add_phenotype(c("ADG", "BW"))

# ---- 3. Binary trait (mortality) ------------------------------------------

pop <- define_trait(pop, "mort", trait_type = "binary", prevalence = 0.08,
                 target_add_var = 1, residual_var = 1)
sel_mort <- get_table(pop, "genome_meta") |> collect() |>
  slice_sample(n = 60) |> pull(locus_name)
pop <- pop |>
  get_table("genome_meta") |>
  filter(locus_name %in% sel_mort) |>
  define_qtl("mort") |>
  define_additive_effects("mort")
pop <- pop |>
  get_table("ind_meta") |>
  filter(gen == 1L) |>
  add_phenotype("mort")

# ---- 4. Summary ----------------------------------------------------------

summary_tbl <- get_table(pop, "ind_phenotype") |>
  group_by(trait_name) |>
  summarise(
    n    = n(),
    mean = mean(value),
    var  = var(value),
    .groups = "drop"
  ) |>
  collect()

cat("\n=== Phenotype summary ===\n")
print(summary_tbl)

# Realised correlation between ADG and BW phenotypes
pheno <- get_table(pop, "ind_phenotype") |>
  filter(trait_name %in% c("ADG", "BW")) |>
  collect()
wide <- reshape(as.data.frame(pheno[, c("id_ind", "trait_name", "pheno_value")]),
                idvar = "id_ind", timevar = "trait_name",
                direction = "wide")
if (all(c("pheno_value.ADG", "pheno_value.BW") %in% names(wide))) {
  cat("\n=== Phenotype correlation diagnostic ===\n")
  cat(sprintf("cor(ADG, BW) = %.3f\n",
              cor(wide$pheno_value.ADG, wide$pheno_value.BW, use = "complete.obs")))
}

close_pop(pop)
cat("\nDemo complete.\n")
