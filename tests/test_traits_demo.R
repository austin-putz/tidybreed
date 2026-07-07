# tests/test_traits_demo.R
#
# End-to-end exercise of the trait + phenotype system on the current
# long-format schema. Run interactively to confirm the API works under
# in-memory DuckDB.

suppressPackageStartupMessages({
  library(dplyr)
})

set.seed(20260420)

# ---- 1. Initialise genome + founders + one generation of offspring --------

pop <- open_pop(pop_name = "demo", db_name = ":memory:") |>
  define_genome(n_loci = 600, n_chr = 6, chr_len_Mb = 100) |>
  define_founder_haplotypes(n_haplotypes = 200, method = "fixed")

pop <- pop |>
  get_table("founder_haplotypes") |>
  add_founders(n_males = 50, n_females = 50, line_name = "A", gen = 0L)

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
# Genetic layer (define_trait) carries additive variance only; the observation
# layer (define_phenotype) carries mean / residual / fixed effects.

G <- matrix(c(0.25, 0.10, 0.10, 0.30), 2, 2)
R <- matrix(c(0.75, 0.20, 0.20, 0.70), 2, 2,
            dimnames = list(c("ADG", "BW"), c("ADG", "BW")))

pop <- define_trait(pop, "ADG", target_add_var = 0.25)
pop <- define_trait(pop, "BW",  target_add_var = 0.30)

sel_qtl <- get_table(pop, "genome_meta") |> collect() |>
  slice_sample(n = 150) |> pull(locus_name)
pop <- pop |>
  get_table("genome_meta") |>
  filter(locus_name %in% sel_qtl) |>
  define_additive_effects(c("ADG", "BW"), G = G, seed = 42)

pop <- define_phenotype(pop, "ADG", type = "continuous", mean = 850)
pop <- define_phenotype(pop, "BW",  type = "continuous", mean = 450)

# Correlated residual (co)variance for the two phenotypes
pop <- define_effect_cov_matrix(pop, "residual", R)

# Fixed sex effect on ADG (observation layer -> phenotype_name)
pop <- define_effect_fixed_class(pop, "ADG", "sex",
  source_column = "sex",
  levels = c(M = 30, F = 0))

pop <- pop |>
  get_table("ind_meta") |>
  filter(gen == 1L) |>
  add_phenotype(c("ADG", "BW"))

# ---- 3. Binary trait (mortality) ------------------------------------------
# Genetic layer: additive variance on the liability scale.
# Observation layer: categorical 2-category via prevalence.

pop <- define_trait(pop, "mort", target_add_var = 1)
sel_mort <- get_table(pop, "genome_meta") |> collect() |>
  slice_sample(n = 60) |> pull(locus_name)
pop <- pop |>
  get_table("genome_meta") |>
  filter(locus_name %in% sel_mort) |>
  define_additive_effects("mort")

pop <- define_phenotype(pop, "mort",
  type         = "categorical",
  prevalence   = 0.08,
  cat_values   = c(0, 1),
  cat_names    = c("Alive", "Dead"),
  residual_var = 1)

pop <- pop |>
  get_table("ind_meta") |>
  filter(gen == 1L) |>
  add_phenotype("mort")

# ---- 4. Summary ----------------------------------------------------------

summary_tbl <- get_table(pop, "ind_phenotype") |>
  collect() |>
  group_by(phenotype_name) |>
  summarise(
    n    = n(),
    mean = mean(pheno_value),
    var  = var(pheno_value),
    .groups = "drop"
  )

cat("\n=== Phenotype summary ===\n")
print(summary_tbl)

stopifnot(all(c("ADG", "BW", "mort") %in% summary_tbl$phenotype_name))

# Realised correlation between ADG and BW phenotypes
pheno <- get_table(pop, "ind_phenotype") |>
  filter(phenotype_name %in% c("ADG", "BW")) |>
  collect()
wide <- reshape(as.data.frame(pheno[, c("id_ind", "phenotype_name", "pheno_value")]),
                idvar = "id_ind", timevar = "phenotype_name",
                direction = "wide")
if (all(c("pheno_value.ADG", "pheno_value.BW") %in% names(wide))) {
  cat("\n=== Phenotype correlation diagnostic ===\n")
  cat(sprintf("cor(ADG, BW) = %.3f\n",
              cor(wide$pheno_value.ADG, wide$pheno_value.BW, use = "complete.obs")))
}

close_pop(pop)
cat("\nDemo complete.\n")
