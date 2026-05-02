#------------------------------------------------------------------------------#
# Load Packages
#------------------------------------------------------------------------------#

library(yaml)
library(tidybreed)
library(tidyverse)

#------------------------------------------------------------------------------#
# Set Input Options
#------------------------------------------------------------------------------#

n_reps <- 3
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

# create 2 custom fields for simulation pipeline
pop %>% 
  get_table("ind_meta") %>%
  mutate_table(
    rep      = NA_integer_,    # rep (replication) number as 'integer'
    gen_born = NA_integer_     # gen (generation) number as 'integer'
  )

# add 'active' field, default to FALSE
pop %>%
  get_table("ind_meta") %>%
  mutate_table(
    active = TRUE,          # this is the default value (no rows yet)
    .set_default = TRUE     # if TRUE, sets default of given value
    # here so when animals are born they are "active" by default, later culled
  )

# add 'rep' + 'gen_pheno' field to 'ind_phenotype'
pop %>%
  get_table("ind_phenotype") %>%
  mutate_table(
    rep = NA_integer_,          # add rep to phenotypes
    gen_pheno = NA_integer_,    # this is the default value (no rows yet)
    .set_default = TRUE         # if TRUE, sets default of given value
    # here we need to identify which generation they were phenotyped
  )

pop %>% get_table("ind_phenotype")

# add 'rep' to 'ind_tbv'
pop %>%
  get_table("ind_tbv") %>%
  mutate_table(
    rep = NA_integer_          # add rep
    #.set_default = TRUE         # if TRUE, sets default of given value
  )

pop %>% get_table("ind_tbv")

# add 'rep' to 'ind_ebv'
pop %>%
  get_table("ind_ebv") %>%
  mutate_table(
    gen_eval = NA_integer_,     # add gen of evaluation
    rep = NA_integer_,          # add rep
    .set_default = TRUE         # if TRUE, sets default of given value
  )

pop %>% get_table("ind_ebv")

pop %>% get_table("ind_index")
pop %>% get_table("index_meta")

#------------------------------------------------------------------------------#
# Add Founders
#------------------------------------------------------------------------------#

#for (repl in 1:n_reps){
repl = 1

# add founders
pop <- pop %>%
  add_founders(          # add founders
    n_males   = 50,      # sample male founders
    n_females = 100,     # sample female founders
    line_name = "Libra", # name this population
    rep = repl,          # USER DEFINED - replicate
    gen_born = 0L        # USER DEFINED - generation
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
    effect_name = "gen_add",          # fixed term for additive genetic (co)variance matrix
    cov_matrix  = mat.add.gen.vars    # name of matrix with row/col names
  )

# add this residual covariance matrix to a table with function
pop <- pop %>%
  add_effect_cov_matrix(
    effect_name = "residual",      # fixed term for residual (co)variance matrix
    cov_matrix  = mat.res.vars     # name for matrix with row/col names
  )

#------------------------------------------------------------#
# Index
#------------------------------------------------------------#

pop %>%
  define_index(
    index_name = "maternal",
    trait_names = c("NW", "ADG", "BF"),
    index_wts = c(93, 1.5, -30)
  )

pop %>% get_table("index_meta")

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
  filter(rep == repl) %>%
  add_tbv("ADG", rep = repl)

pop %>% get_table("ind_tbv")

# add overall mean for ADG
pop %>%
  add_effect_int(
    trait_name = "ADG",              # trait (need to change to "trait_name")
    mean = 1000
  )

# test `add_phenotype()` function
pop %>%
  get_table("ind_meta") %>%     # will phenotype all individuals in this table with no filter
  filter(rep == repl) %>%
  add_phenotype(                # add rows to 'ind_phenotye' table
    trait = "ADG",              # trait name
    rep = repl,                 # set rep number
    gen_pheno = 0L              # founder gen
  )  

pop %>% get_table("ind_phenotype")

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
  filter(rep == repl) %>%
  add_tbv("BF", rep = repl)

pop %>% get_table("ind_tbv") %>% filter(trait_name == "ADG")
pop %>% get_table("ind_tbv") %>% filter(trait_name == "BF")

# add overall mean for ADG
pop %>%
  add_effect_int(
    trait_name = "BF",              # trait (need to change to "trait_name")
    mean = 10
  )

# test `add_phenotype()` function
pop %>%
  get_table("ind_meta") %>%     # will phenotype all individuals in this table with no filter
  filter(rep == repl) %>%
  add_phenotype(                # add rows to 'ind_phenotye' table
    trait = "BF",              # trait name
    rep = repl,                 # add rep
    gen_pheno = 0L              # founder generation
  )

pop %>% get_table("ind_phenotype") %>% filter(trait_name == "ADG")
pop %>% get_table("ind_phenotype") %>% filter(trait_name == "BF")

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
  filter(rep == repl) %>%
  add_tbv("NW", rep = repl)

pop %>% get_table("ind_tbv") %>% filter(trait_name == "NW")

# add overall mean for ADG
pop %>%
  add_effect_int(
    trait_name = "NW",              # trait (need to change to "trait_name")
    mean = 10
  )

# test `add_phenotype()` function
pop %>%
  get_table("ind_meta") %>%     # will phenotype all individuals in this table with no filter
  filter(rep == repl, sex == "F") %>%        # only females get NW phenotype
  add_phenotype(                # add rows to 'ind_phenotye' table
    trait = "NW",                # trait name
    rep = repl,                 # add rep number
    gen_pheno = 0L              # founder gen
  )

pop %>% get_table("ind_phenotype") %>% filter(trait_name == "NW")

#------------------------------------------------------------------------------#
# Run EBVs
#------------------------------------------------------------------------------#

# run EBV for ADG
pop <- pop %>%
  get_table("ind_meta") %>%
  filter(rep == repl) %>%
  add_ebv("ADG", software="blupf90", model="blup", gen_eval=0L, rep = repl)

pop %>% get_table("ind_ebv") %>% filter(trait_name == "ADG")

# run EBV for BF
pop <- pop %>%
  get_table("ind_meta") %>%
  filter(rep == repl) %>%
  add_ebv("BF", software="blupf90", model="blup", gen_eval=0L, rep = repl)

pop %>% get_table("ind_ebv") %>% filter(trait_name == "BF")

# run EBV for NW
pop <- pop %>%
  get_table("ind_meta") %>%
  filter(rep == repl) %>%
  add_ebv("NW", software="blupf90", model="blup", gen_eval=0L, rep = repl)

#------------------------------------------------------------------------------#
# Run Index
#------------------------------------------------------------------------------#

# run index calculation
pop %>%
  get_table("ind_ebv") %>%    # must pass 'ind_ebv' because it contains the EBVs needed
  add_index("maternal")       # just give the index name and it will grab weights

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
  group_by(rep, trait_name) %>%
  summarise(
    MeanTBV = mean(tbv),
    .groups = "drop_last"
  ) %>%
  print(n=10)

#------------------------------------------------------------#
# EBVs
#------------------------------------------------------------#

# print mean of TBV by trait
pop %>%
  get_table("ind_ebv") %>%
  collect() %>%
  group_by(rep, trait_name) %>%
  summarise(
    MeanEBV = mean(ebv),
    .groups = "drop_last"
  ) %>%
  print(n=10)

#------------------------------------------------------------#
# Phenotypes
#------------------------------------------------------------#

# print mean of TBV by trait
pop %>%
  get_table("ind_phenotype") %>%
  collect() %>%
  group_by(rep, trait_name) %>%
  summarise(
    MeanP = mean(value),
    .groups = "drop_last"
  ) %>%
  print(n=10)


#------------------------------------------------------------------------------#
# Run by generation
#------------------------------------------------------------------------------#

for (genl in 1:n_gens){
  
  #---------------- SUBSET CANDIDATES ----------------#
  
  # subset male candidates
  male_candidates <- pop %>%
    get_table("ind_meta") %>%
    filter(rep == repl, sex == "M", active==TRUE) %>%
    pull(id_ind)
  
  # subset female candidates
  female_candidates <- pop %>%
    get_table("ind_meta") %>%
    filter(rep == repl, sex == "F", active==TRUE) %>%
    pull(id_ind)
  
  #---------------- SELECT CANDIDATES ----------------#
  
  # select males
  males_selected <- pop %>%
    get_table("ind_ebv") %>%
    filter(
      trait_name=="ADG", 
      gen_eval == genl - 1L,
      id_ind %in% male_candidates,
      rep == repl
    ) %>%
    slice_max(ebv, n=3) %>%
    pull(id_ind)
  
  # select females
  females_selected <- pop %>%
    get_table("ind_ebv") %>%
    filter(
      trait_name=="ADG", 
      gen_eval == genl - 1L,
      id_ind %in% female_candidates,
      rep == repl
    ) %>%
    slice_max(ebv, n=15) %>%
    pull(id_ind)
  
  #---------------- RESET ACTIVES ----------------#
  
  # RESET ALL to ACTIVE=FALSE
  pop <- pop %>%
    get_table("ind_meta") %>%
    filter(rep == repl) %>%
    mutate_table(
      active = FALSE
    )
  
  # set selected parents as ACTIVE=TRUE
  pop <- pop %>%
    get_table("ind_meta") %>%
    filter(
      rep == repl, 
      id_ind %in% c(males_selected, females_selected)
    ) %>%
    mutate_table(
      active = TRUE
    )
  
  #---------------- SAMPLE NW PHENOTYPE ON DAMS ----------------#
  
  # phenotype selected dams first
  pop %>%
    get_table("ind_meta") %>%
    filter(rep == repl, id_ind %in% females_selected) %>%
    add_phenotype(
      "NW",
      rep = repl,          # add rep number
      gen_pheno = genl     # add generation phenotyped as current gen
    )
  
  # extract NW phenotype to produce offspring numbers correctly
  data.nw <- pop %>%
    get_table("ind_phenotype") %>%
    filter(
      rep == repl, 
      trait_name == "NW",
      id_ind %in% females_selected,
      gen_pheno == genl        # filter to current gen as that's what we set above
    ) %>%
    collect() %>%
    select(id_ind, value)
  
  #---------------- PRODUCE PROGENY ----------------#
  
  # HOW:
  # boars: select males randomly like pooled semen I guess (possible but wouldn't know pedigree)
  # gilts/sows: fill number of rows based on NW phenotype just sampled (rows will vary by chance)
  # sex: randomly assign 50/50
  
  # use new phenotype to build mating plan
  data.new.matings <- tibble(
    id_parent_1 = sample(males_selected, size=sum(data.nw$value), replace=TRUE),
    id_parent_2 = rep(c(data.nw$id_ind), time=data.nw$value),
    line     = "Libra",
    sex      = sample(c("M", "F"), size=sum(data.nw$value), replace=TRUE, prob=c(0.5, 0.5)),
    rep      = repl,
    gen_born = genl
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
    filter(rep == repl, sex=="M") %>%            # MALES ONLY
    add_genotypes(chip_name="9k")
  
  #---------------- ADD PHENOTYPES ----------------#
  
  # all offspring get a genotype for 9k
  pop <- pop %>%
    get_table("ind_meta") %>%
    filter(rep == repl, gen_born == genl) %>%             # Current Gen phenotyped only!
    add_phenotype(
      c("ADG", "BF"),
      rep = repl, 
      gen_pheno = genl
    )
  
  #---------------- RUN BLUPF90 ----------------#
  
  # run EBV for ADG
  pop <- pop %>%
    get_table("ind_meta") %>%
    filter(rep == repl) %>%
    add_ebv("ADG", software="blupf90", model="blup", gen_eval=genl, rep = repl)
  
  # run EBV for BF
  pop <- pop %>%
    get_table("ind_meta") %>%
    filter(rep == repl) %>%
    add_ebv("BF", software="blupf90", model="blup", gen_eval=genl, rep = repl)
  
  # run EBV for NW
  pop <- pop %>%
    get_table("ind_meta") %>%
    filter(rep == repl) %>%
    add_ebv("NW", software="blupf90", model="blup", gen_eval=genl, rep = repl)
  
  #---------------- RUN OVERLAP GENERATION ----------------#
  
  # Offspring already ACTIVE=TRUE by default now, no need to set
  
} # END GENERATION LOOP

# make sure all animals have TBVs for later
pop %>%
  get_table("ind_meta") %>%
  filter(rep == repl) %>%
  add_tbv(c("ADG", "BF", "NW"), rep = repl)

} # END REPLICATE LOOP













#------------------------------------------------------------------------------#
# Summarize Pop
#------------------------------------------------------------------------------#

#------------------------------------------------------------#
# Combine Individual + EBV + TBV data into 1 tibble
#------------------------------------------------------------#

# extract individual data
data.ind <- pop %>%
  get_table("ind_meta") %>%
  collect()

# extract EBV data (only use last generation EBVs)
data.ebv <- pop %>%
  get_table("ind_ebv") %>%
  filter(gen_eval == genl) %>%
  collect()

# pull TBVs
data.tbv <- pop %>%
  get_table("ind_tbv") %>%
  collect()

# join ind data to EBV data
data.tbv.ebv <- left_join(data.ebv, data.ind)

# join TBV to EBV tibble
data.tbv.ebv <- left_join(data.tbv.ebv, data.tbv)

#------------------------------------------------------------#
# EBVs
#------------------------------------------------------------#

# plot EBV Trends
data.tbv.ebv %>%
  group_by(rep, gen_born, trait_name) %>%
  summarise(
    MeanEBV = mean(ebv),
    .groups = "drop_last"
  ) %>%
ggplot(., aes(x=gen_born, y=MeanEBV, color=trait_name, shape=as.factor(rep), group=rep)) +
  geom_line() +
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

# plot EBV Trends
data.tbv.ebv %>%
  group_by(rep, gen_born, trait_name) %>%
  summarise(
    MeanTBV = mean(tbv),
    .groups = "drop_last"
  ) %>%
ggplot(., aes(x=gen_born, y=MeanTBV, color=trait_name, shape=as.factor(rep), group=rep)) +
  geom_line() +
  geom_point() +
  facet_wrap(~trait_name) +
  labs(
    title = "True Breeding Values (TBVs)",
    x = "Generation Born",
    y = "Mean True Breeding Value (TBV)",
    caption = "Swine Breeding Program - 10 generations"
  )



#------------------------------------------------------------------------------#
# Ribbins
#------------------------------------------------------------------------------#

# data.tbv.ebv %>%
# ggplot(., aes(x = gen_born, y = ebv)) +
#   geom_ribbon(aes(ymin = min_val, ymax = max_val), alpha = 0.3) +
#   geom_line(aes(y = mean_val)) +
#   facet_wrap(~trait_name)

#------------------------------------------------------------#
# EBVs
#------------------------------------------------------------#

# plot EBVs
data.tbv.ebv %>%
ggplot(., aes(x = gen_born, y = ebv)) +
  stat_summary(fun.min = min, fun.max = max,
               geom = "ribbon", alpha = 0.2, fill = "steelblue") +
  stat_summary(fun.data = mean_sdl, fun.args = list(mult = 1),
               geom = "ribbon", alpha = 0.5, fill = "steelblue") +
  stat_summary(fun = mean, geom = "line", color = "steelblue", linewidth = 1) +
  facet_wrap(~ trait_name, scales = "free_y") +
  labs(
    title = "EBVs",
    x = "Generation Born",
    y = "Estimated Breeding Values (EBVs)",
    caption = "Swine Breeding Program - 10 Gens"
  )

#------------------------------------------------------------#
# TBVs
#------------------------------------------------------------#

# plot TBVs
data.tbv.ebv %>%
ggplot(., aes(x = gen_born, y = tbv)) +
  stat_summary(fun.min = min, fun.max = max,
               geom = "ribbon", alpha = 0.2, fill = "steelblue") +
  stat_summary(fun.data = mean_sdl, fun.args = list(mult = 1),
               geom = "ribbon", alpha = 0.5, fill = "steelblue") +
  stat_summary(fun = mean, geom = "line", color = "steelblue", linewidth = 1) +
  facet_wrap(~ trait_name, scales = "free_y") +
  labs(
    title = "TBVs",
    x = "Generation Born",
    y = "True Breeding Values (TBVs)",
    caption = "Swine Breeding Program - 10 Gens"
  )

#------------------------------------------------------------#
# EBVs + TBVs
#------------------------------------------------------------#

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
  scale_fill_manual("Value", values = c("steelblue3", "grey50"), labels = c("EBV", "TBV")) +
  scale_color_manual("Value", values = c("steelblue3", "grey50"), labels = c("EBV", "TBV")) +
  scale_x_discrete(limits = factor(0:10)) +
  labs(
    title = "EBVs and TBVs over Time",
    x = "Generation Born",
    y = "Values (EBVs or TBVs)",
    caption = "Swine Breeding Program - 10 Gens"
  )





#------------------------------------------------------------------------------#
# Summary of Data
#------------------------------------------------------------------------------#

#------------------------------------------------------------#
# Animals - Rep + Gen
#------------------------------------------------------------#

# count by rep and gen
pop %>%
  get_table("ind_meta") %>%
  collect() %>%
ggplot(., aes(x=as.factor(gen_born), fill=as.factor(rep))) +
  geom_bar(position="dodge") +
  scale_x_discrete(limits = factor(0:10)) +
  scale_fill_manual("Replicate", 
      values=c("dodgerblue3", "cadetblue3", "mediumorchid2")) +
  labs(
    title = "Count Animal Records by Rep + Generation",
    x = "Generation Born",
    y = "Count",
    caption = "Swine Breeding Program - 10 Gens"
  )

#------------------------------------------------------------#
# Phenotypes - Rep + Gen
#------------------------------------------------------------#

# count by rep and gen
pop %>%
  get_table("ind_phenotype") %>%
  collect() %>%
ggplot(., aes(x=as.factor(gen_pheno), fill=as.factor(rep))) +
  geom_bar(position="dodge") +
  facet_wrap(~trait_name) +
  scale_x_discrete(limits = factor(0:10)) +
  scale_fill_manual("Replicate", 
      values=c("dodgerblue3", "cadetblue3", "mediumorchid2")) +
  labs(
    title = "Count Animal Records by Rep + Generation",
    x = "Generation Born",
    y = "Count",
    caption = "Swine Breeding Program - 10 Gens"
  )

#------------------------------------------------------------#
# TBV - Rep + Gen
#------------------------------------------------------------#

# count by rep and gen
pop %>%
  get_table("ind_tbv") %>%
  collect() %>%
  left_join(., data.ind) %>%
  group_by(rep, gen_born, trait_name) %>%
  summarise(
    MinTBV = min(tbv),
    Q1TBV = quantile(tbv, prob=0.25),
    Q2TBV = quantile(tbv, prob=0.50),
    Q3TBV = quantile(tbv, prob=0.75),
    MaxTBV = max(tbv),
    .groups = "drop_last"
  ) %>%
  mutate(
    gen_born = as.factor(gen_born),
    rep = as.factor(rep)
  ) %>% 
ggplot(aes(x=gen_born, fill=rep, color=rep, group=rep)) +
  geom_ribbon(aes(ymin=MinTBV, ymax=MaxTBV), alpha=0.15) +
  geom_ribbon(aes(ymin=Q1TBV, ymax=Q3TBV), alpha=0.40) +
  geom_line(aes(y=Q2TBV)) +
  facet_wrap(~ trait_name, scale="free_y") +
  scale_color_manual("Replicate", 
      values=c("dodgerblue3", "cadetblue3", "mediumorchid2")) +
  scale_fill_manual("Replicate", 
      values=c("dodgerblue3", "cadetblue3", "mediumorchid2")) +
  labs(
    title = "TBV Trends By Trait and Rep"
  )







# close pop object for database
close_pop(pop)

