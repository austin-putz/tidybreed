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
pop <- add_founders(pop, n_males = 50, n_females = 50, line_name = "A")
pop <- mutate_ind_meta(pop, gen = 0L)

# Produce one generation of offspring
sires <- get_table(pop, "ind_meta") |> filter(sex == "M", gen == 0L) |>
  collect() |> pull(id_ind)
dams  <- get_table(pop, "ind_meta") |> filter(sex == "F", gen == 0L) |>
  collect() |> pull(id_ind)

matings <- tibble::tibble(
  id_parent_1 = sample(sires, 200, replace = TRUE),
  id_parent_2 = sample(dams,  200, replace = TRUE),
  sex         = sample(c("M", "F"), 200, replace = TRUE),
  line        = "A",
  gen         = 1L
)
pop <- add_offspring(pop, matings)

# ---- 2. Two correlated continuous traits ----------------------------------

G <- matrix(c(0.25, 0.10, 0.10, 0.30), 2, 2)
R <- matrix(c(0.75, 0.20, 0.20, 0.70), 2, 2)

pop <- add_trait(pop, "ADG", target_add_var = 0.25, residual_var = 0.75,
                 mean = 850)
pop <- add_trait(pop, "BW",  target_add_var = 0.30, residual_var = 0.70,
                 mean = 450)

pop <- define_qtl(pop, "ADG", n = 150, method = "random")
qtl_tf <- get_table(pop, "genome_meta") |> collect() |> pull(is_QTL_ADG)
pop <- define_qtl(pop, "BW", locus_tf = qtl_tf)

pop <- set_qtl_effects_multi(pop, c("ADG", "BW"), G = G, seed = 42)
pop <- set_residual_cov(pop, c("ADG", "BW"), R = R)

pop <- add_trait_covariate(pop, "ADG", "sex",
  effect_class = "fixed", source_column = "sex",
  levels = c(M = 30, F = 0))

pop <- pop |>
  filter(gen == 1L) |>
  add_phenotype(c("ADG", "BW"))

# ---- 3. Binary trait (mortality) ------------------------------------------

pop <- add_trait(pop, "mort", trait_type = "binary", prevalence = 0.08,
                 target_add_var = 1, residual_var = 1)
pop <- define_qtl(pop, "mort", n = 60)
pop <- set_qtl_effects(pop, "mort")
pop <- pop |> filter(gen == 1L) |> add_phenotype("mort")

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
wide <- reshape(as.data.frame(pheno[, c("id_ind", "trait_name", "value")]),
                idvar = "id_ind", timevar = "trait_name",
                direction = "wide")
if (all(c("value.ADG", "value.BW") %in% names(wide))) {
  cat("\n=== Phenotype correlation diagnostic ===\n")
  cat(sprintf("cor(ADG, BW) = %.3f\n",
              cor(wide$value.ADG, wide$value.BW, use = "complete.obs")))
}

close_pop(pop)
cat("\nDemo complete.\n")
