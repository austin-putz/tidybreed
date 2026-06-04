#------------------------------------------------------------------------------#
# Load Packages
#------------------------------------------------------------------------------#

library(pak)
pak::pak("austin-putz/tidybreed")

library(tidybreed)
library(tidyverse)
library(yaml)

#------------------------------------------------------------------------------#
# Set Input Options
#------------------------------------------------------------------------------#

n_reps <- 3
n_gens <- 10

#------------------------------------------------------------------------------#
# yaml inputs
#------------------------------------------------------------------------------#

config <- yaml::read_yaml("~/Claude/tidybreed-test/scenarios/sires_n_3.yaml")

if (dir.exists(config$output$save_dir)) {
  warning("Output directory already exists")
} else {
  warning("Output directory will be created for user")
  dir.create(config$output$save_dir, recursive = TRUE, showWarnings = TRUE)
}

#------------------------------------------------------------------------------#
# Start Genome + Population Object with database
#------------------------------------------------------------------------------#

pop <- initialize_genome(
  pop_name     = "Swine",
  n_loci       = 10000,
  n_chr        = 18,
  chr_len_Mb   = 50,
  n_haplotypes = 1000,
  overwrite    = TRUE
)

#------------------------------------------------------------------------------#
# Pre-declare custom column schemas
#------------------------------------------------------------------------------#

pop <- get_table(pop, "ind_meta") |>
  mutate_table(
    rep      = NA_integer_,
    gen_born = NA_integer_
  )

pop <- get_table(pop, "ind_meta") |>
  mutate_table(
    active       = TRUE,
    .set_default = TRUE    # new animals are active by default
  )

pop <- get_table(pop, "ind_phenotype") |>
  mutate_table(
    rep       = NA_integer_,
    gen_pheno = NA_integer_,
    .set_default = TRUE
  )

pop <- get_table(pop, "ind_tbv") |>
  mutate_table(rep = NA_integer_)

pop <- get_table(pop, "ind_ebv") |>
  mutate_table(
    gen_eval     = NA_integer_,
    rep          = NA_integer_,
    .set_default = TRUE
  )

pop <- get_table(pop, "ind_index") |>
  mutate_table(
    gen_eval = NA_integer_,
    rep      = NA_integer_
  )

#==============================================================================
# REPLICATE LOOP
#==============================================================================

for (repl in seq_len(n_reps)) {

  warning("\n ----------    REPLICATE: ", repl, "    --------------------\n")

  #----------------------------------------------------------------------------#
  # Add Founders
  #----------------------------------------------------------------------------#

  pop <- pop |>
    add_founders(
      n_males   = 50,
      n_females = 100,
      line_name = "Libra",
      rep       = repl,
      gen_born  = 0L,
      active    = TRUE
    )

  #----------------------------------------------------------------------------#
  # Add 9k SNP Chip
  #----------------------------------------------------------------------------#

  warning("Add 9k chip")

  chip_loci <- get_table(pop, "genome_meta") |>
    collect() |>
    slice_sample(n = 9000) |>
    pull(locus_name)

  pop <- get_table(pop, "genome_meta") |>
    filter(locus_name %in% chip_loci) |>
    define_chip("9k")

  #----------------------------------------------------------------------------#
  # Define Genetic (Co)Variance Matrices
  #----------------------------------------------------------------------------#

  warning("Define covariance matrices")

  # ── Additive genetic (co)variances: NW, ADG, BF ──────────────────────────
  mat.add.gen.vars <- matrix(
    c(0.90,    3.07,  0.21,
      3.07, 1050.00, 17.80,
      0.21,   17.80,  1.20),
    nrow = 3, byrow = TRUE,
    dimnames = list(c("NW", "ADG", "BF"), c("NW", "ADG", "BF"))
  )

  # ── Residual (co)variances: NW, ADG, BF ──────────────────────────────────
  mat.res.vars <- matrix(
    c(8.10,   10.00,  0.65,
     10.00, 2100.00, 10.00,
      0.65,   10.00,  1.30),
    nrow = 3, byrow = TRUE,
    dimnames = list(c("NW", "ADG", "BF"), c("NW", "ADG", "BF"))
  )

  pop <- pop |>
    define_effect_cov_matrix(effect_name = "gen_add",  cov_matrix = mat.add.gen.vars)

  pop <- pop |>
    define_effect_cov_matrix(effect_name = "residual", cov_matrix = mat.res.vars)

  # ── Additive genetic (co)variances: WWD + WWM (weaning weight) ───────────
  # Pig WW at 21 days: direct h2 ~10%, maternal h2 ~15%, VP_total ~3 kg^2
  # Slight negative covariance between direct and maternal effects (common in pigs)
  mat.add.ww <- matrix(
    c( 0.30, -0.05,
      -0.05,  0.45),
    nrow = 2, byrow = TRUE,
    dimnames = list(c("WWD", "WWM"), c("WWD", "WWM"))
  )

  pop <- pop |>
    define_effect_cov_matrix(effect_name = "gen_add", cov_matrix = mat.add.ww)

  #----------------------------------------------------------------------------#
  # Define Selection Index
  # Maternal pig index: NW + ADG + BF + WW
  #----------------------------------------------------------------------------#

  warning("Define maternal index")

  pop <- pop |>
    define_index(
      index_name  = "maternal",
      trait_names = c("NW", "ADG", "BF", "WW"),
      index_wts   = c(93.0, 1.5, -30.0, 25.0),
      overwrite   = TRUE
    )

  get_table(pop, "index_meta")

  #----------------------------------------------------------------------------#
  # Trait: ADG
  #   Genetic layer  — define_trait()
  #   Phenotype layer — define_phenotype()
  #   QTL + effects  — get_table("genome_meta") |> filter() |> define_additive_effects()
  #----------------------------------------------------------------------------#

  warning("Define ADG")

  pop <- pop |>
    define_trait(
      trait_name      = "ADG",
      description     = "Average Daily Gain",
      units           = "g/d",
      target_add_mean = 0,
      overwrite       = TRUE
    )

  pop <- pop |>
    define_phenotype(
      phenotype_name = "ADG",
      type           = "continuous",
      mean           = 1000,
      residual_var   = 2100,
      overwrite      = TRUE
    )

  # 1000 random QTL; effects scaled to target_add_var stored in trait_effect_cov
  qtl_adg <- get_table(pop, "genome_meta") |> collect() |>
    slice_sample(n = 1000) |> pull(locus_name)

  pop <- get_table(pop, "genome_meta") |>
    filter(locus_name %in% qtl_adg) |>
    define_additive_effects("ADG")

  # Compute TBVs for founders
  pop <- get_table(pop, "ind_meta") |>
    filter(rep == repl) |>
    add_tbv("ADG", rep = repl)

  get_table(pop, "ind_tbv")

  #----------------------------------------------------------------------------#
  # Trait: BF
  #----------------------------------------------------------------------------#

  warning("Define BF")

  pop <- pop |>
    define_trait(
      trait_name      = "BF",
      description     = "Ultrasound Backfat",
      units           = "mm",
      target_add_mean = 0,
      overwrite       = TRUE
    )

  pop <- pop |>
    define_phenotype(
      phenotype_name = "BF",
      type           = "continuous",
      mean           = 10,
      residual_var   = 1.30,
      min_value      = 0,
      overwrite      = TRUE
    )

  qtl_bf <- get_table(pop, "genome_meta") |> collect() |>
    slice_sample(n = 1000) |> pull(locus_name)

  pop <- get_table(pop, "genome_meta") |>
    filter(locus_name %in% qtl_bf) |>
    define_additive_effects("BF")

  pop <- get_table(pop, "ind_meta") |>
    filter(rep == repl) |>
    add_tbv("BF", rep = repl)

  get_table(pop, "ind_tbv") |> filter(trait_name == "BF")

  #----------------------------------------------------------------------------#
  # Trait: NW (Number Weaned)
  #   expressed_sex = "F": only females can express this phenotype
  #----------------------------------------------------------------------------#

  warning("Define NW")

  pop <- pop |>
    define_trait(
      trait_name      = "NW",
      description     = "Number Weaned",
      units           = "count",
      target_add_mean = 0,
      overwrite       = TRUE
    )

  pop <- pop |>
    define_phenotype(
      phenotype_name = "NW",
      type           = "count",
      mean           = 10,
      residual_var   = 8.10,
      min_value      = 1,
      expressed_sex  = "F",
      overwrite      = TRUE
    )

  qtl_nw <- get_table(pop, "genome_meta") |> collect() |>
    slice_sample(n = 1000) |> pull(locus_name)

  pop <- get_table(pop, "genome_meta") |>
    filter(locus_name %in% qtl_nw) |>
    define_additive_effects("NW")

  pop <- get_table(pop, "ind_meta") |>
    filter(rep == repl) |>
    add_tbv("NW", rep = repl)

  get_table(pop, "ind_tbv") |> filter(trait_name == "NW")

  #----------------------------------------------------------------------------#
  # Trait: WW (Weaning Weight) — maternal model
  #
  #   Two genetic components in trait_meta:
  #     WWD — direct additive effect (piglet's own genes)
  #     WWM — maternal additive effect (dam's genes that affect piglet environment)
  #
  #   Observed phenotype model:
  #     WW = mu + WWD_self + WWM_dam + residual
  #
  #   formula_tbv assembles the composite TBV at add_phenotype() time.
  #   Founders (no dam) are skipped automatically (missing_component_action = "skip").
  #----------------------------------------------------------------------------#

  warning("Define WW — direct (WWD) and maternal (WWM) genetic components")

  pop <- pop |>
    define_trait(
      trait_name      = "WWD",
      description     = "Weaning Weight — direct genetic component",
      units           = "kg",
      target_add_mean = 0,
      overwrite       = TRUE
    )

  pop <- pop |>
    define_trait(
      trait_name      = "WWM",
      description     = "Weaning Weight — maternal genetic component",
      units           = "kg",
      target_add_mean = 0,
      overwrite       = TRUE
    )

  # Shared QTL for WWD + WWM with correlated effects (reads G from trait_effect_cov)
  qtl_ww <- get_table(pop, "genome_meta") |> collect() |>
    slice_sample(n = 1000) |> pull(locus_name)

  pop <- get_table(pop, "genome_meta") |>
    filter(locus_name %in% qtl_ww) |>
    define_additive_effects(c("WWD", "WWM"))

  # Compute founder TBVs for both components (used when offspring add_phenotype runs)
  pop <- get_table(pop, "ind_meta") |>
    filter(rep == repl) |>
    add_tbv(c("WWD", "WWM"), rep = repl)

  # Observed WW phenotype — composite TBV from formula_tbv DSL
  pop <- pop |>
    define_phenotype(
      phenotype_name           = "WW",
      type                     = "continuous",
      mean                     = 6.5,     # mean pig WW at 21 days (kg)
      residual_var             = 2.30,
      min_value                = 0.5,
      formula_tbv              = "WWD + dam(WWM)",
      missing_component_action = "skip",  # founders with no dam are excluded
      overwrite                = TRUE
    )

  #----------------------------------------------------------------------------#
  # Founder Phenotypes (generation 0)
  #----------------------------------------------------------------------------#

  warning("Add founder phenotypes: ADG + BF")

  pop <- get_table(pop, "ind_meta") |>
    filter(rep == repl) |>
    add_phenotype(c("ADG", "BF"), rep = repl, gen_pheno = 0L)

  get_table(pop, "ind_phenotype")

  warning("Add NW for founder females")

  pop <- get_table(pop, "ind_meta") |>
    filter(rep == repl, sex == "F") |>
    add_phenotype("NW", rep = repl, gen_pheno = 0L)

  get_table(pop, "ind_phenotype") |> filter(phenotype_name == "NW")

  # Note: WW not added for founders — they have no dam so all would be skipped

  #----------------------------------------------------------------------------#
  # Founder EBVs
  #----------------------------------------------------------------------------#

  warning("Run founder EBVs")

  pop <- get_table(pop, "ind_meta") |>
    filter(rep == repl) |>
    add_ebv("ADG", software = "blupf90", model = "blup", gen_eval = 0L, rep = repl)

  pop <- get_table(pop, "ind_meta") |>
    filter(rep == repl) |>
    add_ebv("BF",  software = "blupf90", model = "blup", gen_eval = 0L, rep = repl)

  pop <- get_table(pop, "ind_meta") |>
    filter(rep == repl) |>
    add_ebv("NW",  software = "blupf90", model = "blup", gen_eval = 0L, rep = repl)

  # WW EBV at gen 0: no phenotype records yet, so BLUP returns near-zero estimates
  # from pedigree only — these improve rapidly once offspring have WW phenotypes
  pop <- get_table(pop, "ind_meta") |>
    filter(rep == repl) |>
    add_ebv("WW",  software = "blupf90", model = "blup", gen_eval = 0L, rep = repl)

  #----------------------------------------------------------------------------#
  # Founder Index
  #----------------------------------------------------------------------------#

  warning("Calculate founder indexes")

  pop <- get_table(pop, "ind_ebv") |>
    filter(rep == repl, gen_eval == 0L) |>
    add_index("maternal", gen_eval = 0L, rep = repl)

  get_table(pop, "ind_index") |> filter(index_name == "maternal", rep == repl)

  #----------------------------------------------------------------------------#
  # Founder Summary
  #----------------------------------------------------------------------------#

  warning("Founder generation means")

  get_table(pop, "ind_tbv") |>
    collect() |>
    group_by(rep, trait_name) |>
    summarise(MeanTBV = mean(tbv_value), .groups = "drop_last") |>
    print(n = 20)

  get_table(pop, "ind_ebv") |>
    collect() |>
    group_by(rep, trait_name) |>
    summarise(MeanEBV = mean(ebv_value), .groups = "drop_last") |>
    print(n = 20)

  get_table(pop, "ind_phenotype") |>
    collect() |>
    group_by(rep, phenotype_name) |>
    summarise(MeanP = mean(pheno_value), .groups = "drop_last") |>
    print(n = 20)

  #============================================================================
  # GENERATION LOOP
  #============================================================================

  warning("Begin generation loop")

  for (genl in seq_len(n_gens)) {

    warning("\n ----------    GEN: ", genl, "    --------------------\n")

    #---------- SUBSET CANDIDATES ----------------------------------------#

    warning("Select candidates")

    male_candidates <- get_table(pop, "ind_meta") |>
      filter(rep == repl, sex == "M", active == TRUE) |>
      pull(id_ind)

    female_candidates <- get_table(pop, "ind_meta") |>
      filter(rep == repl, sex == "F", active == TRUE) |>
      pull(id_ind)

    #---------- SELECT PARENTS -------------------------------------------#

    warning("Select parents on maternal index")

    males_selected <- get_table(pop, "ind_index") |>
      filter(
        index_name == "maternal",
        id_ind     %in% male_candidates,
        gen_eval   == genl - 1L,
        rep        == repl
      ) |>
      slice_max(index_value, n = 3) |>
      pull(id_ind)

    females_selected <- get_table(pop, "ind_index") |>
      filter(
        index_name == "maternal",
        id_ind     %in% female_candidates,
        gen_eval   == genl - 1L,
        rep        == repl
      ) |>
      slice_max(index_value, n = 15) |>
      pull(id_ind)

    #---------- RESET ACTIVES --------------------------------------------#

    warning("Reset actives")

    pop <- get_table(pop, "ind_meta") |>
      filter(rep == repl) |>
      mutate_table(active = FALSE)

    pop <- get_table(pop, "ind_meta") |>
      filter(rep == repl, id_ind %in% c(males_selected, females_selected)) |>
      mutate_table(active = TRUE)

    #---------- NW PHENOTYPE ON SELECTED DAMS ----------------------------#

    warning("Sample NW phenotype on selected dams")

    pop <- get_table(pop, "ind_meta") |>
      filter(rep == repl, id_ind %in% females_selected) |>
      add_phenotype("NW", rep = repl, gen_pheno = genl)

    warning("Extract NW phenotype")

    data.nw <- get_table(pop, "ind_phenotype") |>
      filter(
        rep            == repl,
        phenotype_name == "NW",
        id_ind         %in% females_selected,
        gen_pheno      == genl
      ) |>
      collect() |>
      select(id_ind, pheno_value)

    #---------- PRODUCE PROGENY ------------------------------------------#

    warning("Build progeny mating plan")

    data.new.matings <- tibble(
      id_parent_1 = sample(males_selected,
                           size = sum(data.nw$pheno_value), replace = TRUE),
      id_parent_2 = rep(data.nw$id_ind,
                        times = as.integer(data.nw$pheno_value)),
      line_name   = "Libra",
      sex         = sample(c("M", "F"),
                           size = sum(data.nw$pheno_value), replace = TRUE),
      rep         = repl,
      gen_born    = genl
    )

    warning("Add offspring")

    pop <- pop |> add_offspring(data.new.matings)

    #---------- ADD GENOTYPES --------------------------------------------#

    warning("Add genotypes (males only)")

    pop <- get_table(pop, "ind_meta") |>
      filter(rep == repl, sex == "M") |>
      add_genotypes(chip_name = "9k")

    #---------- ADD PHENOTYPES -------------------------------------------#

    warning("Add ADG + BF phenotypes for current gen")

    pop <- get_table(pop, "ind_meta") |>
      filter(rep == repl, gen_born == genl) |>
      add_phenotype(c("ADG", "BF"), rep = repl, gen_pheno = genl)

    warning("Add WW phenotypes for current gen offspring")

    # add_phenotype() automatically computes WWD + WWM TBVs (including dam TBVs)
    # via the formula_tbv DSL; individuals without a dam are skipped with a warning
    pop <- get_table(pop, "ind_meta") |>
      filter(rep == repl, gen_born == genl) |>
      add_phenotype("WW", rep = repl, gen_pheno = genl)

    #---------- RUN BLUPF90 ----------------------------------------------#

    warning("Run single-trait BLUPF90")

    pop <- get_table(pop, "ind_meta") |>
      filter(rep == repl) |>
      add_ebv("ADG", software = "blupf90", model = "blup", gen_eval = genl, rep = repl)

    pop <- get_table(pop, "ind_meta") |>
      filter(rep == repl) |>
      add_ebv("BF",  software = "blupf90", model = "blup", gen_eval = genl, rep = repl)

    pop <- get_table(pop, "ind_meta") |>
      filter(rep == repl) |>
      add_ebv("NW",  software = "blupf90", model = "blup", gen_eval = genl, rep = repl)

    # WW EBV: single-trait BLUP ignores maternal structure but is sufficient
    # for index ranking; upgrade to maternal model for unbiased genetic trend
    pop <- get_table(pop, "ind_meta") |>
      filter(rep == repl) |>
      add_ebv("WW",  software = "blupf90", model = "blup", gen_eval = genl, rep = repl)

    #---------- CALCULATE INDEX ------------------------------------------#

    warning("Calculate maternal index")

    pop <- get_table(pop, "ind_ebv") |>
      filter(gen_eval == genl, rep == repl) |>
      add_index("maternal", gen_eval = genl, rep = repl)

  } # END GENERATION LOOP

  warning("End gen loop — ensure all animals have TBVs for summary")

  pop <- get_table(pop, "ind_meta") |>
    filter(rep == repl) |>
    add_tbv(c("ADG", "BF", "NW", "WWD", "WWM"), rep = repl)

} # END REPLICATE LOOP


#==============================================================================
# Summarize Population
#==============================================================================

#------------------------------------------------------------#
# Combine ind_meta + ind_ebv + ind_tbv into one tibble
#------------------------------------------------------------#

data.ind <- get_table(pop, "ind_meta") |> collect()

# Use last generation EBVs for trends
data.ebv <- get_table(pop, "ind_ebv") |>
  filter(gen_eval == n_gens) |>
  collect()

data.tbv <- get_table(pop, "ind_tbv") |> collect()

# Join: EBV + ind_meta (gen_born, sex, rep), then attach TBV
data.tbv.ebv <- left_join(data.ebv,     data.ind, by = "id_ind") |>
                left_join(data.tbv,               by = c("id_ind", "trait_name", "rep"))

#------------------------------------------------------------#
# EBV Trends
#------------------------------------------------------------#

data.tbv.ebv |>
  group_by(rep, gen_born, trait_name) |>
  summarise(MeanEBV = mean(ebv_value), .groups = "drop_last") |>
  ggplot(aes(x = gen_born, y = MeanEBV,
             color = trait_name, shape = as.factor(rep), group = rep)) +
  geom_line() +
  geom_point() +
  facet_wrap(~trait_name, scales = "free_y") +
  labs(
    title   = "Estimated Breeding Values (EBVs)",
    x       = "Generation Born",
    y       = "Mean EBV",
    caption = "Swine Breeding Program — 10 generations"
  )

#------------------------------------------------------------#
# TBV Trends
#------------------------------------------------------------#

data.tbv.ebv |>
  group_by(rep, gen_born, trait_name) |>
  summarise(MeanTBV = mean(tbv_value), .groups = "drop_last") |>
  ggplot(aes(x = gen_born, y = MeanTBV,
             color = trait_name, shape = as.factor(rep), group = rep)) +
  geom_line() +
  geom_point() +
  facet_wrap(~trait_name, scales = "free_y") +
  labs(
    title   = "True Breeding Values (TBVs)",
    x       = "Generation Born",
    y       = "Mean TBV",
    caption = "Swine Breeding Program — 10 generations"
  )

#------------------------------------------------------------#
# EBV + TBV ribbon overlay
#------------------------------------------------------------#

data.tbv.ebv |>
  pivot_longer(cols = c("ebv_value", "tbv_value"),
               names_to = "type", values_to = "value") |>
  ggplot(aes(x = gen_born, y = value, fill = type, color = type)) +
  stat_summary(fun.min = min, fun.max = max,
               geom = "ribbon", alpha = 0.5) +
  stat_summary(fun.data = mean_sdl, fun.args = list(mult = 1),
               geom = "ribbon", alpha = 0.5) +
  stat_summary(fun = mean, geom = "line", linewidth = 1, alpha = 0.7) +
  facet_wrap(~trait_name, scales = "free_y", ncol = 3) +
  scale_fill_manual("Value",  values = c("steelblue3", "grey50"),
                    labels = c("EBV", "TBV")) +
  scale_color_manual("Value", values = c("steelblue3", "grey50"),
                    labels = c("EBV", "TBV")) +
  scale_x_discrete(limits = factor(0:n_gens)) +
  labs(
    title   = "EBVs and TBVs over Time",
    x       = "Generation Born",
    y       = "Value (EBV or TBV)",
    caption = "Swine Breeding Program — 10 generations"
  )

#------------------------------------------------------------#
# Animal counts by rep + generation
#------------------------------------------------------------#

get_table(pop, "ind_meta") |>
  collect() |>
  ggplot(aes(x = as.factor(gen_born), fill = as.factor(rep))) +
  geom_bar(position = "dodge") +
  scale_x_discrete(limits = factor(0:n_gens)) +
  scale_fill_manual("Replicate",
    values = c("dodgerblue3", "cadetblue3", "mediumorchid2")) +
  labs(
    title   = "Animal Records by Rep + Generation",
    x       = "Generation Born",
    y       = "Count",
    caption = "Swine Breeding Program — 10 generations"
  )

#------------------------------------------------------------#
# Phenotype records by rep + generation
#------------------------------------------------------------#

get_table(pop, "ind_phenotype") |>
  collect() |>
  ggplot(aes(x = as.factor(gen_pheno), fill = as.factor(rep))) +
  geom_bar(position = "dodge") +
  facet_wrap(~phenotype_name) +
  scale_x_discrete(limits = factor(0:n_gens)) +
  scale_fill_manual("Replicate",
    values = c("dodgerblue3", "cadetblue3", "mediumorchid2")) +
  labs(
    title   = "Phenotype Records by Rep + Generation",
    x       = "Generation Phenotyped",
    y       = "Count",
    caption = "Swine Breeding Program — 10 generations"
  )

#------------------------------------------------------------#
# TBV ribbon by trait + rep
#------------------------------------------------------------#

get_table(pop, "ind_tbv") |>
  collect() |>
  left_join(data.ind, by = "id_ind") |>
  group_by(rep, gen_born, trait_name) |>
  summarise(
    MinTBV = min(tbv_value),
    Q1TBV  = quantile(tbv_value, 0.25),
    Q2TBV  = quantile(tbv_value, 0.50),
    Q3TBV  = quantile(tbv_value, 0.75),
    MaxTBV = max(tbv_value),
    .groups = "drop_last"
  ) |>
  mutate(gen_born = as.factor(gen_born), rep = as.factor(rep)) |>
  ggplot(aes(x = gen_born, fill = rep, color = rep, group = rep)) +
  geom_hline(yintercept = 0, color = "red3", linewidth = 0.75, linetype = 3) +
  geom_ribbon(aes(ymin = MinTBV, ymax = MaxTBV), alpha = 0.05) +
  geom_ribbon(aes(ymin = Q1TBV,  ymax = Q3TBV),  alpha = 0.40) +
  geom_line(aes(y = Q2TBV)) +
  facet_wrap(~trait_name, scales = "free_y") +
  scale_color_manual("Replicate",
    values = c("dodgerblue3", "cadetblue3", "mediumorchid2")) +
  scale_fill_manual("Replicate",
    values = c("dodgerblue3", "cadetblue3", "mediumorchid2")) +
  labs(
    title    = "TBV Trends by Trait and Rep",
    subtitle = "Median + middle 50% + min/max",
    x        = "Generation Born",
    y        = "TBV",
    caption  = "Swine Breeding Program — 10 generations"
  )

close_pop(pop)
