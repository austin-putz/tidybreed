
library(yaml)
library(tidybreed)
library(tidyverse)

n_reps <- 5
n_gens <- 10

#------------------------------------------------------------------------------#
# Start Genome + Population Object with database
#------------------------------------------------------------------------------#

# start population by building a genome
pop <- initialize_genome(
         pop_name     = "Swine",  # different than "line_name" below, this is the entire genome + DB name
		     n_loci       = 10000,    # number of loci (nothing assigned as SNP or QTL yet...)
		     n_chr        = 18,       # number of chromosomes
		     chr_len_Mb   = 50,       # length in Mb (1,000,000 bp) (e.g. 1.20 and not 1_200_000)
		     n_haplotypes = 1000,     # number of random haplotypes generated (no LD generated)
		     overwrite    = TRUE      # overwrite the database if it exists from a previous run
		  )

#------------------------------------------------------------------------------#
# Add custom fields!
#------------------------------------------------------------------------------#

# create THREE (3) new CUSTOM Columns/Fields
pop %>% 
  get_table("ind_meta") %>%
  mutate_table(
    rep = NA_integer_,    # replication number as 'integer'
    gen = NA_integer_     # generation number as 'integer'
  )

# add 'active' field, default to FALSE
pop %>%
  get_table("ind_meta") %>%
  mutate_table(
    active = TRUE,
    .set_default = TRUE
  )

#------------------------------------------------------------------------------#
# Add Founders
#------------------------------------------------------------------------------#

#for (repl in 1:n_reps){

# add founders
pop <- pop %>%
  add_founders(          # add founders
    n_males   = 50,      # sample male founders
    n_females = 100,     # sample female founders
    line_name = "Libra", # name this population
    rep = 1L,          # USER DEFINED - replicate
    gen = 0L             # USER DEFINED - generation
  )

#------------------------------------------------------------------------------#
# Add SNP Chip
#------------------------------------------------------------------------------#

# add 9k SNP Chip
pop %>%
  define_chip(
    chip_name = "9k",
    n = 9000,
    method = "random"
  )

#------------------------------------------------------------------------------#
# Add Traits
#------------------------------------------------------------------------------#

# additive genetic (CO)VARIANCES
mat.add.gen.vars <- matrix(c(0.90,    3.07,  0.21,
                             3.07, 1050.00, 17.80,
                             0.21,   17.80,  1.20), 
                      nrow = 3, byrow=TRUE, 
                      dimnames = list(c("NW", "ADG", "BF"), c("NW", "ADG", "BF")))

# residual (CO)VARIANCES
mat.res.vars <- matrix(c(8.10,   10.00,  0.65,
                        10.00, 2100.00, 10.00,
                         0.65,   10.00,  1.30), 
                      nrow = 3, byrow=TRUE, 
                      dimnames = list(c("NW", "ADG", "BF"), c("NW", "ADG", "BF")))

# add this additive genetic covariance matrix to a table with function
pop <- pop %>%
  add_effect_cov_matrix(
    effect_name = "gen_add",
    cov_matrix  = mat.add.gen.vars
  )

# add this residual covariance matrix to a table with function
pop <- pop %>%
  add_effect_cov_matrix(
    effect_name = "residual",
    cov_matrix  = mat.res.vars
  )

#------------------------------------------------------------#
# Trait: ADG
#------------------------------------------------------------#

# add ADG as a trait
pop <- pop %>%
  add_trait(
    trait_name  = "ADG",
    description = "Average Daily Gain", 
    units       = "g/d",                  # grams per day during testing period
    trait_type  = "continuous",           # not categorical or binary
    repeatable  = FALSE,                  # only 1 phenotype per individual
    recorded_on = "self",                 # recorded on itself only
    target_add_mean = 0,      # mean TBV in 'base'
    #target_add_var  = 1082,  # additive variance target in 'base'
    #residual_var    = 2500,  # residual variance (target h2 ~0.30) - 'fixed'
    index_weight    = 1.49,   # index weight (not implemented yet)
    economic_value  = 1.49,   # economic value (not implemented yet)
    overwrite       = TRUE    # wipe this row if it exists and replace with this new data
  )

# add which loci are QTL and their effects
pop %>%
  # set loci as QTL for this trait
  define_qtl(
    trait_name = "ADG",         # name of trait
    n = 1000,                   # set number of QTL for this trait (will set locations)
    method = "random"           # allocate to loci randomly
  ) %>%
  # set QTL effects
  set_qtl_effects(
    trait_name      = "ADG",        # trait name
    distribution    = "normal",     # distribution of QTL effects
    scale_to_target = TRUE,         # scale to meet additive variance target
    base            = "current_pop" # use all animals in pop to standardized (or if filtered)
  )

# add all TBV for ADG
pop <- pop %>%
  get_table("ind_meta") %>% # here we specify the 'ind_meta' table so all animals will have their TBV calculated
  add_tbv("ADG")

# add overall mean for ADG
pop %>%
  add_effect_int(
    trait_name = "ADG",              # trait (need to change to "trait_name")
    mean = 1000
  )

# test `add_phenotype()` function
pop %>%
  get_table("ind_meta") %>%     # will phenotype all individuals in this table with no filter
  add_phenotype(                # add rows to 'ind_phenotye' table
    trait = "ADG"              # trait name
  )  

#------------------------------------------------------------#
# Trait: BF
#------------------------------------------------------------#

# add ADG as a trait
pop <- pop %>%
  add_trait(
    trait_name  = "BF",
    description = "Ultrasound Backfat", 
    units       = "mm",                  # grams per day during testing period
    trait_type  = "continuous",           # not categorical or binary
    repeatable  = FALSE,                  # only 1 phenotype per individual
    recorded_on = "self",                 # recorded on itself only
    target_add_mean = 0,      # mean TBV in 'base'
    min_value       = 0,      # cannot be negative
    #target_add_var  = 1082,  # additive variance target in 'base'
    #residual_var    = 2500,  # residual variance (target h2 ~0.30) - 'fixed'
    index_weight    = -28.61, # index weight (not implemented yet)
    economic_value  = -28.61, # economic value (not implemented yet)
    overwrite       = TRUE    # wipe this row if it exists and replace with this new data
  )

# add which loci are QTL and their effects
pop %>%
  # set loci as QTL for this trait
  define_qtl(
    trait_name = "BF",         # name of trait
    n = 1000,                   # set number of QTL for this trait (will set locations)
    method = "random"           # allocate to loci randomly
  ) %>%
  # set QTL effects
  set_qtl_effects(
    trait_name      = "BF",        # trait name
    distribution    = "normal",     # distribution of QTL effects
    scale_to_target = TRUE,         # scale to meet additive variance target
    base            = "current_pop" # use all animals in pop to standardized (or if filtered)
  )

# add all TBV for ADG
pop <- pop %>%
  get_table("ind_meta") %>% # here we specify the 'ind_meta' table so all animals will have their TBV calculated
  add_tbv("BF")

# add overall mean for ADG
pop %>%
  add_effect_int(
    trait_name = "BF",              # trait (need to change to "trait_name")
    mean = 10
  )

# test `add_phenotype()` function
pop %>%
  get_table("ind_meta") %>%     # will phenotype all individuals in this table with no filter
  add_phenotype(                # add rows to 'ind_phenotye' table
    trait = "BF"              # trait name
  )

#------------------------------------------------------------#
# Trait: NW
#------------------------------------------------------------#

# add ADG as a trait
pop <- pop %>%
  add_trait(
    trait_name  = "NW",
    description = "Number Weaned", 
    units       = "count",                  # grams per day during testing period
    trait_type  = "count",           # not categorical or binary
    repeatable  = FALSE,                  # only 1 phenotype per individual
    recorded_on = "dam",                 # recorded on itself only
    target_add_mean = 0,      # mean TBV in 'base'
    min_value       = 1,      # cannot be negative
    #target_add_var  = 1082,  # additive variance target in 'base'
    #residual_var    = 2500,  # residual variance (target h2 ~0.30) - 'fixed'
    index_weight    = 91.93, # index weight (not implemented yet)
    economic_value  = 92.93, # economic value (not implemented yet)
    overwrite       = TRUE    # wipe this row if it exists and replace with this new data
  )

# add which loci are QTL and their effects
pop %>%
  # set loci as QTL for this trait
  define_qtl(
    trait_name = "NW",         # name of trait
    n = 1000,                   # set number of QTL for this trait (will set locations)
    method = "random"           # allocate to loci randomly
  ) %>%
  # set QTL effects
  set_qtl_effects(
    trait_name      = "NW",        # trait name
    distribution    = "normal",     # distribution of QTL effects
    scale_to_target = TRUE,         # scale to meet additive variance target
    base            = "current_pop" # use all animals in pop to standardized (or if filtered)
  )

# add all TBV for ADG
pop <- pop %>%
  get_table("ind_meta") %>% # here we specify the 'ind_meta' table so all animals will have their TBV calculated
  add_tbv("NW")

# add overall mean for ADG
pop %>%
  add_effect_int(
    trait_name = "NW",              # trait (need to change to "trait_name")
    mean = 10
  )

# test `add_phenotype()` function
pop %>%
  get_table("ind_meta") %>%     # will phenotype all individuals in this table with no filter
  filter(sex == "F") %>%        # only females get NW phenotype
  add_phenotype(                # add rows to 'ind_phenotye' table
    trait = "NW"                # trait name
  )

#------------------------------------------------------------------------------#
# Run EBVs
#------------------------------------------------------------------------------#

# first add gen as field in 'ind_ebv'
pop <- pop %>%
  get_table("ind_ebv") %>%
  mutate_table(
    gen = NA_integer_
  )

# run EBV for ADG
pop <- pop %>%
  get_table("ind_meta") %>%
  add_ebv("ADG", software="blupf90", model="blup", gen=0L)

# run EBV for BF
pop <- pop %>%
  get_table("ind_meta") %>%
  add_ebv("BF", software="blupf90", model="blup", gen=0L)

# run EBV for NW
pop <- pop %>%
  get_table("ind_meta") %>%
  add_ebv("NW", software="blupf90", model="blup", gen=0L)

#------------------------------------------------------------------------------#
# Checks
#------------------------------------------------------------------------------#

#------------------------------------------------------------#
# True Breeding Values
#------------------------------------------------------------#

# print mean of TBV by trait
pop %>%
  get_table("ind_tbv") %>%
  collect() %>%
  group_by(trait_name) %>%
  summarise(
    MeanTBV = mean(tbv)
  ) %>%
  print(n=10)

#------------------------------------------------------------#
# EBVs
#------------------------------------------------------------#

# print mean of TBV by trait
pop %>%
  get_table("ind_ebv") %>%
  collect() %>%
  group_by(trait_name) %>%
  summarise(
    MeanEBV = mean(ebv)
  ) %>%
  print(n=10)

#------------------------------------------------------------#
# Phenotypes
#------------------------------------------------------------#

# print mean of TBV by trait
pop %>%
  get_table("ind_phenotype") %>%
  collect() %>%
  group_by(trait_name) %>%
  summarise(
    MeanP = mean(value)
  ) %>%
  print(n=10)


#------------------------------------------------------------------------------#
# Run by generation
#------------------------------------------------------------------------------#

for (genl in 1L:n_gens){
  
  #---------------- SUBSET CANDIDATES ----------------#
  
  # subset male candidates
  male_candidates <- pop %>%
    get_table("ind_meta") %>%
    filter(sex == "M", active==TRUE) %>%
    pull(id_ind)
  
  # subset female candidates
  female_candidates <- pop %>%
    get_table("ind_meta") %>%
    filter(sex == "F", active==TRUE) %>%
    pull(id_ind)
  
  #---------------- SELECT CANDIDATES ----------------#
  
  # select males
  males_selected <- pop %>%
    get_table("ind_ebv") %>%
    filter(
      trait_name=="ADG", 
      gen == genl - 1L,
      id_ind %in% male_candidates
    ) %>%
    slice_max(ebv, n=5) %>%
    pull(id_ind)
  
  # select females
  females_selected <- pop %>%
    get_table("ind_ebv") %>%
    filter(
      trait_name=="ADG", 
      gen == genl - 1L,
      id_ind %in% female_candidates
    ) %>%
    slice_max(ebv, n=30) %>%
    pull(id_ind)
  
  #---------------- RESET ACTIVES ----------------#
  
  # RESET ALL to ACTIVE=FALSE
  pop <- pop %>%
    get_table("ind_meta") %>%
    mutate_table(
      active = FALSE
    )
  
  # set selected parents as ACTIVE=TRUE
  pop <- pop %>%
    get_table("ind_meta") %>%
    filter(
      id_ind %in% c(males_selected, females_selected)
    ) %>%
    mutate_table(
      active = TRUE
    )
  
  #---------------- SAMPLE NW PHENOTYPE ON DAMS ----------------#
  
  # phenotype selected dams first
  pop %>%
    get_table("ind_meta") %>%
    filter(id_ind %in% females_selected) %>%
    add_phenotype(
      "NW",
      gen = genl
    )
  
  # extract NW phenotype to produce offspring numbers correctly
  data.nw <- pop %>%
    get_table("ind_phenotype") %>%
    filter(
      trait_name == "NW",
      id_ind %in% females_selected,
      gen == genl
    ) %>%
    collect() %>%
    select(id_ind, value)
  
  #---------------- PRODUCE PROGENY ----------------#
  
  # use new phenotype to build mating plan
  data.new.matings <- tibble(
    id_parent_1 = sample(males_selected, size=sum(data.nw$value), replace=TRUE),
    id_parent_2 = rep(c(data.nw$id_ind), time=data.nw$value),
    line = "Libra",
    sex = sample(c("M", "F"), size=sum(data.nw$value), replace=TRUE, prob=c(0.5, 0.5)),
    gen = genl
  )
  
  # add new offspring based on tibble mating plan (1 row per offspring)
  pop %>%
    add_offspring(
      data.new.matings
    )
  
  #---------------- ADD GENOTYPES ----------------#
  
  # all Males get a genotype for 9k
  pop <- pop %>%
    get_table("ind_meta") %>%
    filter(sex=="M") %>%            # MALES ONLY
    add_genotypes(chip_name="9k")
  
  #---------------- ADD PHENOTYPES ----------------#
  
  # all offspring get a genotype for 9k
  pop <- pop %>%
    get_table("ind_meta") %>%
    filter(gen == genl) %>%             # Current Gen phenotyped only!
    add_phenotype(c("ADG", "BF"))
  
  #---------------- RUN BLUPF90 ----------------#
  
  # run EBV for ADG
  pop <- pop %>%
    get_table("ind_meta") %>%
    add_ebv("ADG", software="blupf90", model="blup", gen=genl)
  
  # run EBV for BF
  pop <- pop %>%
    get_table("ind_meta") %>%
    add_ebv("BF", software="blupf90", model="blup", gen=genl)
  
  # run EBV for NW
  pop <- pop %>%
    get_table("ind_meta") %>%
    add_ebv("NW", software="blupf90", model="blup", gen=genl)
  
  #---------------- RUN OVERLAP GENERATION ----------------#
  
  # Offspring already ACTIVE=TRUE by default now, no need to set
  
} # END GENERATION LOOP
#} # END REPLICATE LOOP






pop %>%
  get_table("ind_meta") %>%
  add_tbv(c("ADG", "BF", "NW"))








#------------------------------------------------------------------------------#
# Summarize Pop
#------------------------------------------------------------------------------#

#------------------------------------------------------------#
# EBVs
#------------------------------------------------------------#

data.ind <- pop %>%
  get_table("ind_meta") %>%
  collect() %>%
  rename(gen_born = gen)

data.ebv <- pop %>%
  get_table("ind_ebv") %>%
  filter(gen == genl) %>%
  collect()

# join ind data to EBV data
data.tbv.ebv <- left_join(data.ebv, data.ind)

# plot EBV Trends
data.tbv.ebv %>%
  group_by(gen_born, trait_name) %>%
  summarise(
    MeanEBV = mean(ebv),
    .groups = "drop_last"
  ) %>%
ggplot(., aes(x=gen_born, y=MeanEBV)) +
  geom_point() +
  facet_wrap(~trait_name) +
  labs(
    title = "Estimated Breeding Values (EBVs)",
    x = "Generation Born",
    y = "Mean Estimated Breeding Value (EBV)",
    caption = "Swine Breeding Program - 10 generations"
  )

#------------------------------------------------------------#
# TBVs
#------------------------------------------------------------#

# pull TBVs
data.tbv <- pop %>%
  get_table("ind_tbv") %>%
  collect()

# join TBV to EBV tibble
data.tbv.ebv <- left_join(data.tbv.ebv, data.tbv)

# plot TBV trends
data.tbv.ebv %>%
  group_by(gen_born, trait_name) %>%
  summarise(
    MeanTBV = mean(tbv),
    .groups = "drop_last"
  ) %>%
ggplot(., aes(x=gen_born, y=MeanTBV)) +
  geom_point() +
  facet_wrap(~trait_name) +
  labs(
    title = "True Breeding Values (TBVs)",
    x = "Generation Born",
    y = "Mean True Breeding Value (TBV)",
    caption = "Swine Breeding Program - 10 generations"
  )

# data.tbv.ebv %>%
# ggplot(., aes(x = gen_born, y = ebv)) +
#   geom_ribbon(aes(ymin = min_val, ymax = max_val), alpha = 0.3) +
#   geom_line(aes(y = mean_val)) +
#   facet_wrap(~trait_name)

# plot EBVs
data.tbv.ebv %>%
ggplot(., aes(x = gen_born, y = ebv)) +
  stat_summary(fun.min = min, fun.max = max,
               geom = "ribbon", alpha = 0.2, fill = "steelblue") +
  stat_summary(fun.data = mean_sdl, fun.args = list(mult = 1),
               geom = "ribbon", alpha = 0.3, fill = "steelblue") +
  stat_summary(fun = mean, geom = "line", color = "steelblue", linewidth = 1) +
  facet_wrap(~ trait_name, scales = "free_y") +
  labs(
    title = "EBVs",
    x = "Generation Born",
    y = "Estimated Breeding Values (EBVs)",
    caption = "Swine Breeding Program - 10 Gens"
  )

# plot TBVs
data.tbv.ebv %>%
ggplot(., aes(x = gen_born, y = tbv)) +
  stat_summary(fun.min = min, fun.max = max,
               geom = "ribbon", alpha = 0.2, fill = "steelblue") +
  stat_summary(fun.data = mean_sdl, fun.args = list(mult = 1),
               geom = "ribbon", alpha = 0.3, fill = "steelblue") +
  stat_summary(fun = mean, geom = "line", color = "steelblue", linewidth = 1) +
  facet_wrap(~ trait_name, scales = "free_y") +
  labs(
    title = "TBVs",
    x = "Generation Born",
    y = "True Breeding Values (TBVs)",
    caption = "Swine Breeding Program - 10 Gens"
  )

# plot TBVs
data.tbv.ebv %>%
  pivot_longer(
    cols = c("tbv", "ebv")
  ) %>%
ggplot(., aes(x = gen_born, y = value, fill=name, color=name)) +
  stat_summary(fun.min = min, fun.max = max,
               geom = "ribbon", alpha = 0.5) +
  stat_summary(fun.data = mean_sdl, fun.args = list(mult = 1),
               geom = "ribbon", alpha = 0.5) +
  stat_summary(fun = mean, geom = "line", linewidth = 1, alpha=0.7) +
  facet_wrap(~ trait_name, scales = "free_y", ncol=3) +
  scale_fill_manual("Value", values = c("steelblue3", "grey50")) +
  scale_color_manual("Value", values = c("steelblue3", "grey50")) +
  scale_x_discrete(limits = factor(0:10)) +
  labs(
    title = "EBVs and TBVs over Time",
    x = "Generation Born",
    y = "Values (EBVs or TBVs)",
    caption = "Swine Breeding Program - 10 Gens"
  )

# pop %>%
#   get_table("ind_phenotype") %>%
#   collect() %>%
#   count(trait_name, gen) %>%
# ggplot(., aes(x=gen, y=n, fill=trait_name)) +
#   geom_col() +
#   facet_wrap(~trait_name)

# ADG EBV
data.tbv.ebv %>%
  filter(trait_name=="ADG") %>%
ggplot(., aes(x = gen_born, y = ebv)) +
  stat_summary(fun.min = min, fun.max = max,      # shaded = full range
               geom = "ribbon", alpha = 0.2, fill = "steelblue") +
  stat_summary(fun.data = mean_sdl,               # darker band = ±1 SD
               fun.args = list(mult = 1),
               geom = "ribbon", alpha = 0.3, fill = "steelblue") +
  stat_summary(fun = mean, geom = "line",          # line = mean
               color = "steelblue", linewidth = 1)

# ADG TBV
data.tbv.ebv %>%
  filter(trait_name=="ADG") %>%
ggplot(., aes(x = gen_born, y = tbv)) +
  stat_summary(fun.min = min, fun.max = max,      # shaded = full range
               geom = "ribbon", alpha = 0.2, fill = "steelblue") +
  stat_summary(fun.data = mean_sdl,               # darker band = ±1 SD
               fun.args = list(mult = 1),
               geom = "ribbon", alpha = 0.3, fill = "steelblue") +
  stat_summary(fun = mean, geom = "line",          # line = mean
               color = "steelblue", linewidth = 1)


# ADG EBV
data.tbv.ebv %>%
  filter(trait_name=="NW") %>%
ggplot(., aes(x = gen_born, y = ebv)) +
  stat_summary(fun.min = min, fun.max = max,      # shaded = full range
               geom = "ribbon", alpha = 0.2, fill = "steelblue") +
  stat_summary(fun.data = mean_sdl,               # darker band = ±1 SD
               fun.args = list(mult = 1),
               geom = "ribbon", alpha = 0.3, fill = "steelblue") +
  stat_summary(fun = mean, geom = "line",          # line = mean
               color = "steelblue", linewidth = 1)

# ADG TBV
data.tbv.ebv %>%
  filter(trait_name=="NW") %>%
ggplot(., aes(x = gen_born, y = tbv)) +
  stat_summary(fun.min = min, fun.max = max,      # shaded = full range
               geom = "ribbon", alpha = 0.2, fill = "steelblue") +
  stat_summary(fun.data = mean_sdl,               # darker band = ±1 SD
               fun.args = list(mult = 1),
               geom = "ribbon", alpha = 0.3, fill = "steelblue") +
  stat_summary(fun = mean, geom = "line",          # line = mean
               color = "steelblue", linewidth = 1)

# close pop object for database
#close_pop(pop)

