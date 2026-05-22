# test_sge_adg_selection.R
#
# 20-generation discrete-generation pig ADG selection experiment
# using the Social Genetic Effects (SGE / Bijma) model.
#
# Genetic model
#   ADG_phenotype = ADG_direct (own) + sum(pen-mates' ADG_SGE) + residual
#
# Genetic evaluation
#   Generation 0 founders: EBV = h2_eff * (y - generation_mean)
#   Generation 1+        : parent_avg via add_ebv(parent_avg = TRUE)
#
# Variance components used (single-trait animal model approximation)
#   VA_direct = 90    (h2_direct  ≈ 0.30)
#   VA_social = 15    (h2_social  ≈ 0.05)
#   r_ds      = -0.30 (Bijma competition)
#   VE        = 210   (residual)
#   pen_size  = 10
#   Vp_total  = VA_direct + (n-1)^2*VA_social + 2*(n-1)*r*sqrt(Vd*Vs) + VE
#
# Run from repo root:
#   Rscript tests/test_sge_adg_selection.R

suppressPackageStartupMessages({
  library(tidybreed)
  library(dplyr)
  library(tibble)
})

# ── Parameters ────────────────────────────────────────────────────────────────

set.seed(20250520)

N_MALES_SELECT   <- 10L    # sires selected per generation
N_FEMALES_SELECT <- 50L    # dams selected per generation
N_OFFSPRING      <- 300L   # piglets born per generation
PEN_SIZE         <- 10L    # pen size (constant)
N_GEN            <- 20L    # generations to run

MEAN_ADG         <- 850    # population mean (g/day)
VA_DIRECT        <- 90     # additive variance, direct effect
VA_SOCIAL        <- 15     # additive variance, social (SGE) effect
R_DS             <- -0.30  # genetic correlation direct-social
VE               <- 210    # residual variance

# Approx total phenotypic variance in a pen-mate SGE model (Bijma 2010)
n1 <- PEN_SIZE - 1L
Vp <- VA_DIRECT +
      n1^2 * VA_SOCIAL +
      2 * n1 * R_DS * sqrt(VA_DIRECT * VA_SOCIAL) +
      VE

# Effective h2 for phenotypic regression EBV (gen 0 only)
h2_eff <- VA_DIRECT / Vp

G_mat <- matrix(
  c(VA_DIRECT,
    R_DS * sqrt(VA_DIRECT * VA_SOCIAL),
    R_DS * sqrt(VA_DIRECT * VA_SOCIAL),
    VA_SOCIAL),
  nrow = 2, ncol = 2,
  dimnames = list(c("ADG_direct", "ADG_SGE"),
                  c("ADG_direct", "ADG_SGE"))
)

cat("=================================================================\n")
cat(" Pig ADG — Social Genetic Effects selection experiment\n")
cat("=================================================================\n")
cat(sprintf(" Pens: %d pigs/pen | Select %dM + %dF from %d offspring/gen\n",
            PEN_SIZE, N_MALES_SELECT, N_FEMALES_SELECT, N_OFFSPRING))
cat(sprintf(" VA_direct=%.0f  VA_social=%.0f  r_ds=%.2f  VE=%.0f\n",
            VA_DIRECT, VA_SOCIAL, R_DS, VE))
cat(sprintf(" Vp_approx=%.1f  h2_eff(dir)=%.3f\n\n", Vp, h2_eff))

# ── Helpers ───────────────────────────────────────────────────────────────────

# Globally unique pen IDs so no two generations share a pen label
make_pen_ids <- function(gen, n) {
  n_pens <- ceiling(n / PEN_SIZE)
  rep(paste0("g", gen, "_pen", seq_len(n_pens)), each = PEN_SIZE)[seq_len(n)]
}

# Write EBVs directly to ind_ebv.
# add_ebv() validates trait_name against trait_meta, but "ADG" is a composite
# phenotype that only lives in phenotype_meta — it has no trait_meta row.
# Until add_ebv() is extended to accept phenotype_meta names, we write the
# ind_ebv rows directly. The table format is identical to what add_ebv() produces.
write_ebvs <- function(pop, id_inds, ebv_values, model = "pa_v0", eval_num = 0L) {
  n        <- length(id_inds)
  first_id <- tidybreed:::next_int_id(pop$db_conn, "ind_ebv", "id_ebv")
  rows <- tibble::tibble(
    id_ebv      = seq(first_id, first_id + n - 1L),
    id_ind      = id_inds,
    trait_name  = "ADG",
    model       = model,
    ebv_value   = as.numeric(ebv_values),
    acc         = as.numeric(NA),
    se          = as.numeric(NA),
    eval_number = as.integer(eval_num)
  )
  DBI::dbWriteTable(pop$db_conn, "ind_ebv", rows, append = TRUE)
  pop
}

# Gen 0: phenotypic regression EBV = h2_eff * (y - generation_mean)
seed_pheno_ebvs <- function(pop, pheno_df, model = "pa_v0", eval_num = 0L) {
  mu   <- mean(pheno_df$pheno_value)
  ebvs <- h2_eff * (pheno_df$pheno_value - mu)
  write_ebvs(pop, pheno_df$id_ind, ebvs, model = model, eval_num = eval_num)
}

# Gen 1+: parent average EBV = 0.5 * (sire_EBV + dam_EBV)
compute_parent_avg_ebvs <- function(pop, new_ids, model = "pa_v0", eval_num = 1L) {
  meta <- collect(get_table(pop, "ind_meta") |> filter(id_ind %in% new_ids))
  ebv  <- collect(
    get_table(pop, "ind_ebv") |> filter(trait_name == "ADG", model == !!model)
  ) |>
    group_by(id_ind) |>
    slice_max(eval_number, n = 1L, with_ties = FALSE) |>
    ungroup()

  ebv_map <- setNames(ebv$ebv_value, ebv$id_ind)

  pa <- vapply(seq_len(nrow(meta)), function(i) {
    s <- ebv_map[meta$id_parent_1[i]]
    d <- ebv_map[meta$id_parent_2[i]]
    if (is.na(s) && is.na(d)) 0
    else if (is.na(s)) d / 2
    else if (is.na(d)) s / 2
    else (s + d) / 2
  }, numeric(1L))

  write_ebvs(pop, meta$id_ind, pa, model = model, eval_num = eval_num)
}

# Summarise one generation into one results row
summarise_gen <- function(pop, gen, candidate_ids) {
  ph <- collect(
    get_table(pop, "ind_phenotype") |>
      filter(phenotype_name == "ADG", pheno_number == 1L,
             id_ind %in% candidate_ids)
  )
  tbv <- collect(get_table(pop, "ind_tbv") |> filter(id_ind %in% candidate_ids))
  d   <- filter(tbv, trait_name == "ADG_direct")$tbv_value
  s   <- filter(tbv, trait_name == "ADG_SGE")$tbv_value

  tibble::tibble(
    gen              = gen,
    n_phenotyped     = nrow(ph),
    mean_pheno       = mean(ph$pheno_value),
    sd_pheno         = sd(ph$pheno_value),
    mean_tbv_direct  = mean(d),
    mean_tbv_social  = mean(s),
    # Total genetic value in the pen = direct + (n-1)*social
    mean_total_gv    = mean(d) + n1 * mean(s)
  )
}

# Select top animals by EBV from a generation, return id vectors
select_parents <- function(pop, gen, model = "pa_v0") {
  meta <- collect(get_table(pop, "ind_meta") |> filter(gen == !!gen))
  ebv  <- collect(
    get_table(pop, "ind_ebv") |>
      filter(id_ind %in% meta$id_ind, trait_name == "ADG",
             model == !!model)
  ) |>
    group_by(id_ind) |>
    slice_max(eval_number, n = 1L, with_ties = FALSE) |>
    ungroup()

  ebv_m <- filter(ebv, id_ind %in% filter(meta, sex == "M")$id_ind)
  ebv_f <- filter(ebv, id_ind %in% filter(meta, sex == "F")$id_ind)

  sires <- slice_max(ebv_m, ebv_value, n = N_MALES_SELECT)$id_ind
  dams  <- slice_max(ebv_f, ebv_value, n = N_FEMALES_SELECT)$id_ind
  list(sires = sires, dams = dams)
}

# ── Initialise population ─────────────────────────────────────────────────────

pop <- initialize_genome(
  pop_name          = "pig_sge_sel",
  n_loci            = 2000L,
  n_chr             = 18L,
  chr_len_Mb        = 100,
  n_haplotypes      = 400L,
  db_path           = ":memory:",
  min_allele_freq   = 0.05,
  max_allele_freq   = 0.95
)

# Pre-declare typed columns on empty tables
pop <- get_table(pop, "ind_meta") |> mutate_table(gen    = NA_integer_,
                                                   pen_id = NA_character_)

# ── Define genetic model ──────────────────────────────────────────────────────

pop <- define_trait(pop, "ADG_direct", target_add_var = VA_DIRECT)
pop <- define_trait(pop, "ADG_SGE",   target_add_var = VA_SOCIAL)

pop <- get_table(pop, "genome_meta") |>
  define_additive_effects(c("ADG_direct", "ADG_SGE"), G = G_mat)

pop <- define_phenotype(pop, "ADG",
  trait_type   = "continuous",
  mean         = MEAN_ADG,
  residual_var = VE,
  components   = tibble::tribble(
    ~source_trait_name, ~contributor_type, ~group_column,
    "ADG_direct",       "self",            NA_character_,
    "ADG_SGE",          "group",           "pen_id"
  )
)

# ── Generation 0: founders ────────────────────────────────────────────────────

# Extra founders so we have enough candidates for first selection round
n_founders_m <- N_MALES_SELECT   * 10L
n_founders_f <- N_FEMALES_SELECT * 3L

pop <- add_founders(pop,
  n_males   = n_founders_m,
  n_females = n_founders_f,
  line_name = "Pig",
  gen       = 0L)

g0_ids  <- collect(get_table(pop, "ind_meta") |> filter(gen == 0L))$id_ind
g0_pens <- make_pen_ids(0L, length(g0_ids))
pop <- get_table(pop, "ind_meta") |> filter(gen == 0L) |>
  mutate_table(pen_id = g0_pens)

# Phenotype founders
pop <- get_table(pop, "ind_meta") |> filter(gen == 0L) |>
  add_phenotype("ADG")

# Seed EBVs from phenotype regression (no parents for gen 0)
g0_ph  <- collect(
  get_table(pop, "ind_phenotype") |>
    filter(phenotype_name == "ADG", pheno_number == 1L,
           id_ind %in% g0_ids)
)
pop <- seed_pheno_ebvs(pop, g0_ph, model = "pa_v0", eval_num = 0L)

# Collect gen 0 results
results <- summarise_gen(pop, gen = 0L, candidate_ids = g0_ids)

cat(sprintf(
  "Gen %2d | n=%3d | mean_pheno=%6.1f | TBV_dir=%+6.2f | TBV_sge=%+6.3f | Total_GV=%+6.2f\n",
  results$gen, results$n_phenotyped, results$mean_pheno,
  results$mean_tbv_direct, results$mean_tbv_social, results$mean_total_gv
))

# ── Selection loop ────────────────────────────────────────────────────────────

for (gen in seq_len(N_GEN)) {

  prev_gen <- gen - 1L

  # 1. Select parents from previous generation by EBV
  sel <- select_parents(pop, gen = prev_gen, model = "pa_v0")

  # 2. Random matings among selected parents
  matings <- tibble::tibble(
    id_parent_1 = sample(sel$sires, N_OFFSPRING, replace = TRUE),
    id_parent_2 = sample(sel$dams,  N_OFFSPRING, replace = TRUE),
    sex         = sample(c("M", "F"), N_OFFSPRING, replace = TRUE),
    line_name   = "Pig"
  )
  pop <- add_offspring(pop, matings = matings)

  # 3. Tag new offspring with current gen + pen assignment
  new_ids <- collect(get_table(pop, "ind_meta") |> filter(is.na(gen)))$id_ind
  pop <- get_table(pop, "ind_meta") |>
    filter(id_ind %in% new_ids) |>
    mutate_table(gen = as.integer(gen))

  new_pens <- make_pen_ids(gen, length(new_ids))
  pop <- get_table(pop, "ind_meta") |>
    filter(id_ind %in% new_ids) |>
    mutate_table(pen_id = new_pens)

  # 4. Phenotype current generation
  pop <- get_table(pop, "ind_meta") |>
    filter(gen == !!as.integer(gen)) |>
    add_phenotype("ADG")

  # 5. Parent-average EBV (0.5 * sire_EBV + 0.5 * dam_EBV)
  #    Note: add_ebv(parent_avg=TRUE) validates against trait_meta; "ADG" is a
  #    composite phenotype (phenotype_meta only), so we compute PA manually.
  pop <- compute_parent_avg_ebvs(pop, new_ids, model = "pa_v0", eval_num = gen)

  # 6. Collect results
  row <- summarise_gen(pop, gen = as.integer(gen), candidate_ids = new_ids)
  results <- bind_rows(results, row)

  # Selected parents' mean EBV for reporting selection differential
  g0_sel_ebv <- collect(
    get_table(pop, "ind_ebv") |>
      filter(id_ind %in% c(sel$sires, sel$dams),
             trait_name == "ADG", model == "pa_v0")
  ) |>
    group_by(id_ind) |>
    slice_max(eval_number, n = 1L, with_ties = FALSE) |>
    ungroup()

  cat(sprintf(
    "Gen %2d | n=%3d | mean_pheno=%6.1f | TBV_dir=%+6.2f | TBV_sge=%+6.3f | Total_GV=%+6.2f | sel_EBV=%+6.2f\n",
    gen, row$n_phenotyped, row$mean_pheno,
    row$mean_tbv_direct, row$mean_tbv_social, row$mean_total_gv,
    mean(g0_sel_ebv$ebv_value)
  ))
}

# ── Summary table ─────────────────────────────────────────────────────────────

cat("\n=================================================================\n")
cat(" Results by generation\n")
cat("=================================================================\n")
cat(sprintf("%-4s  %-6s  %-10s  %-10s  %-10s  %-10s\n",
            "Gen", "n", "mean_pheno", "TBV_direct", "TBV_social", "Total_GV"))
cat(strrep("-", 62), "\n")
for (i in seq_len(nrow(results))) {
  r <- results[i, ]
  cat(sprintf("%-4d  %-6d  %-10.1f  %-+10.3f  %-+10.4f  %-+10.3f\n",
              r$gen, r$n_phenotyped, r$mean_pheno,
              r$mean_tbv_direct, r$mean_tbv_social, r$mean_total_gv))
}

# Genetic gain per generation (linear regression on TBV over gens 1-N_GEN)
fit_dir   <- lm(mean_tbv_direct ~ gen,  data = filter(results, gen >= 1L))
fit_soc   <- lm(mean_tbv_social ~ gen,  data = filter(results, gen >= 1L))
fit_total <- lm(mean_total_gv   ~ gen,  data = filter(results, gen >= 1L))
fit_pheno <- lm(mean_pheno      ~ gen,  data = filter(results, gen >= 1L))

cat("\nLinear genetic trend (gen 1 to", N_GEN, "):\n")
cat(sprintf("  Delta TBV_direct  per gen: %+.3f g/day\n", coef(fit_dir)[2]))
cat(sprintf("  Delta TBV_social  per gen: %+.4f g/day\n", coef(fit_soc)[2]))
cat(sprintf("  Delta Total_GV    per gen: %+.3f g/day\n", coef(fit_total)[2]))
cat(sprintf("  Delta mean_pheno  per gen: %+.3f g/day\n", coef(fit_pheno)[2]))

cat("\nNote: Total_GV = TBV_direct + (pen_size-1) * TBV_social\n")
cat("      Negative TBV_social trend = selection against social aggression/competition\n")
cat("      (expected with r_ds < 0 when selecting on direct performance)\n")

close_pop(pop)

# ── Plots ─────────────────────────────────────────────────────────────────────
#
# Saved to: tests/test_sge_adg_selection.pdf
#
# Panel A: Mean phenotype +/- 1 SD over generations
# Panel B: Genetic trends — TBV_direct, (n-1)*TBV_social, Total_GV (all zero-centred)
# Panel C: TBV_social raw values (shows antagonistic response to direct selection)
# Panel D: Phenotypic mean vs Total_GV — illustrates social drag masking true gain

if (!requireNamespace("ggplot2",  quietly = TRUE)) stop("Install ggplot2 to produce plots.")
if (!requireNamespace("patchwork", quietly = TRUE)) stop("Install patchwork to produce plots.")

library(ggplot2)
library(patchwork)

plot_file <- "tests/test_sge_adg_selection.pdf"

# Tidy data for plotting
res_plot <- results |>
  mutate(
    social_contribution = n1 * mean_tbv_social,   # (pen_size-1) * mean_TBV_social
    # zero-centre all genetic values relative to gen 0
    d_direct  = mean_tbv_direct  - mean_tbv_direct[gen == 0],
    d_social  = mean_tbv_social  - mean_tbv_social[gen == 0],
    d_soc_con = social_contribution - social_contribution[gen == 0],
    d_total   = mean_total_gv    - mean_total_gv[gen == 0],
    d_pheno   = mean_pheno       - mean_pheno[gen == 0]
  )

theme_sge <- theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank(),
        plot.title = element_text(size = 11, face = "bold"),
        legend.position = "bottom",
        legend.key.width = unit(1.2, "cm"))

pal <- c(
  "TBV direct"           = "#1b7837",
  "(n-1) × TBV social"   = "#d6604d",
  "Total GV"             = "#2166ac",
  "Mean phenotype"       = "#762a83"
)

# Panel A: Phenotype mean ± 1 SD
pA <- ggplot(res_plot, aes(gen, mean_pheno)) +
  geom_ribbon(aes(ymin = mean_pheno - sd_pheno,
                  ymax = mean_pheno + sd_pheno),
              fill = "#762a83", alpha = 0.12) +
  geom_line(colour = "#762a83", linewidth = 0.8) +
  geom_point(colour = "#762a83", size = 1.8) +
  geom_hline(yintercept = MEAN_ADG, linetype = "dashed", colour = "grey50") +
  labs(title = "A  Phenotypic mean ± 1 SD",
       x = "Generation", y = "Mean ADG (g/day)") +
  theme_sge

# Panel B: Zero-centred genetic trends on one axis
long_B <- tidyr::pivot_longer(
  res_plot,
  cols      = c(d_direct, d_soc_con, d_total),
  names_to  = "component",
  values_to = "delta"
) |>
  mutate(component = recode(component,
    d_direct  = "TBV direct",
    d_soc_con = "(n-1) × TBV social",
    d_total   = "Total GV"
  ))

pB <- ggplot(long_B, aes(gen, delta, colour = component, linetype = component)) +
  geom_hline(yintercept = 0, colour = "grey70") +
  geom_line(linewidth = 0.8) +
  geom_point(size = 1.6) +
  scale_colour_manual(values = pal, name = NULL) +
  scale_linetype_manual(
    values = c("TBV direct" = "solid",
               "(n-1) × TBV social" = "dashed",
               "Total GV"  = "solid"),
    name = NULL
  ) +
  labs(title = "B  Genetic trends (zero-centred at gen 0)",
       x = "Generation", y = "Delta genetic value (g/day)") +
  theme_sge

# Panel C: Raw TBV_social — shows antagonistic correlated response
pC <- ggplot(res_plot, aes(gen, mean_tbv_social)) +
  geom_hline(yintercept = 0, colour = "grey70") +
  geom_line(colour = "#d6604d", linewidth = 0.8) +
  geom_point(colour = "#d6604d", size = 1.8) +
  labs(title = "C  Mean TBV_social (raw) -- correlated response",
       x = "Generation",
       y = "Mean TBV_social (g/day)",
       caption = sprintf("r(direct, social) = %.2f; pen size n = %d", R_DS, PEN_SIZE)) +
  theme_sge +
  theme(plot.caption = element_text(size = 8, colour = "grey40"))

# Panel D: Phenotypic mean vs Total GV (both zero-centred) — shows social drag
long_D <- tidyr::pivot_longer(
  res_plot,
  cols      = c(d_pheno, d_total),
  names_to  = "series",
  values_to = "delta"
) |>
  mutate(series = recode(series,
    d_pheno = "Mean phenotype",
    d_total = "Total GV"
  ))

pD <- ggplot(long_D, aes(gen, delta, colour = series)) +
  geom_hline(yintercept = 0, colour = "grey70") +
  geom_line(linewidth = 0.8) +
  geom_point(size = 1.6) +
  scale_colour_manual(values = pal, name = NULL) +
  labs(title = "D  Phenotype vs Total GV -- social drag",
       x = "Generation", y = "Delta value (g/day)",
       caption = paste0(
         sprintf("VA_direct=%.0f  VA_social=%.0f  VE=%.0f  Vp~%.0f\n",
                 VA_DIRECT, VA_SOCIAL, VE, Vp),
         "Total GV = TBV_direct + (n-1) × TBV_social"
       )) +
  theme_sge +
  theme(plot.caption = element_text(size = 8, colour = "grey40"))

# Assemble and save
fig <- (pA | pB) / (pC | pD) +
  plot_annotation(
    title    = "Pig ADG: 20 generations of phenotypic selection under SGE (Bijma model)",
    subtitle = sprintf(
      "10 sires × 50 dams × 300 offspring/gen | %d pigs/pen | r(d,s)=%.2f | seed=20250520",
      PEN_SIZE, R_DS
    ),
    theme = theme(
      plot.title    = element_text(size = 13, face = "bold"),
      plot.subtitle = element_text(size = 9,  colour = "grey30")
    )
  )

ggsave(plot_file, fig, width = 11, height = 8.5, device = "pdf")
cat("\nPlot saved to:", plot_file, "\n")
