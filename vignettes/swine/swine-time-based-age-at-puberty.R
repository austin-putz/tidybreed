#------------------------------------------------------------------------------#
# Load Packages
#------------------------------------------------------------------------------#

#install.packages("pak")
#library(pak)
#pak::pak("austin-putz/tidybreed")

if ("tidybreed" %in% installed.packages()){
  message("tidybreed is installed")
} else {
  warning("tidybreed not installed, will install for you")
  if ("pak" %in% installed.packages()){
    message("pak installed")
    library(pak)
  } else {
    install.packages("pak")
    library(pak)
  }
  pak::pak("austin-putz/tidybreed")
}

library(yaml)
library(tidyverse)
library(tidybreed)

#------------------------------------------------------------------------------#
# Set other input options outside yaml file
#------------------------------------------------------------------------------#

hg_green      <- rgb(150, 200,  40, max=255)   # "yellowgreen"
hg_blue       <- rgb(  0, 150, 180, max=255)   # "deepskyblue3"
hypor_blue    <- rgb(  0,  65, 120, max=255)   # "navy"

#------------------------------------------------------------------------------#
# yaml inputs
#------------------------------------------------------------------------------#

#args <- commandArgs(trailingOnly = TRUE)
#config_path <- args[1]
config_path <- "~/Claude/tidybreed-test/scenarios/age_at_puberty.yaml"

# check for yaml file name
if (is.na(config_path)) {
  stop("No yaml file provided")
} else {
  message("yaml file name provided")
}

# check if yaml file exists
if (!file.exists(config_path)) {
  stop("yaml file does not exist")
} else {
  message("yaml file exists")
}

# read config file
message("reading yaml file now...")
config <- yaml::read_yaml(config_path)

summary(config)

# create output directory
if (dir.exists(config$output$save_dir)){
  warning("Output directory already exists")
} else {
  warning("Output directory will be created for user")
  # create directory if it doesn't exist
  dir.create(config$output$save_dir, recursive=TRUE, showWarnings = TRUE)
}

start_date <- as.Date(config$general$start_date)

#------------------------------------------------------------------------------#
# Start Genome + Population Object with database
#------------------------------------------------------------------------------#

# start population by building a genome
pop <- initialize_genome(
         pop_name     = config$scenario_name,  # different than "line_name" below, this is the entire genome + DB name
		     n_loci       = config$genome$n_loci,       # number of loci (nothing assigned as SNP or QTL yet...)
		     n_chr        = config$genome$n_chr,        # number of chromosomes
		     chr_len_Mb   = config$genome$chr_len,      # length in Mb (1,000,000 bp) (e.g. 1.20 and not 1_200_000)
		     n_haplotypes = config$genome$n_haplotypes, # number of random haplotypes generated (no LD generated)
		     overwrite    = TRUE      # overwrite the database if it exists from a previous run
		  )

#------------------------------------------------------------------------------#
# Add custom fields!
#------------------------------------------------------------------------------#

# create 2 custom fields for simulation pipeline
pop %>% 
  get_table("ind_meta") %>%
  mutate_table(
    rep          = NA_integer_,    # rep (replication) number as 'integer'
    status       = NA_character_,  # status [e.g. 'juvenile', 'off-test-gilt', 'gest', 'lact', etc]
    conc_date    = as.Date(NA),
    birth_date   = as.Date(NA),
    on_test_date = as.Date(NA),
    off_test_date = as.Date(NA),
    puberty_date = as.Date(NA),
    mate_date   = as.Date(NA),
    farrow_date = as.Date(NA),
    wean_date   = as.Date(NA),
    cull_date   = as.Date(NA),
    death_date  = as.Date(NA)
  )

# add 'active' field, default to FALSE
pop %>%
  get_table("ind_meta") %>%
  mutate_table(
    alive = TRUE,
    active = FALSE,          # this is the default value (no rows yet)
    .set_default = TRUE     # if TRUE, sets default of given value
    # here so when animals are born they are "active" by default, later culled
  )

# add 'rep' + 'gen_pheno' field to 'ind_phenotype'
pop %>%
  get_table("ind_phenotype") %>%
  mutate_table(
    rep = NA_integer_,          # add rep to phenotypes
    pheno_date = as.Date(NA),
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
    rep = NA_integer_,          # add rep
    eval_date = as.Date(NA),
    .set_default = TRUE         # if TRUE, sets default of given value
  )

pop %>% get_table("ind_ebv")

# add 'rep' to 'ind_ebv'
pop %>%
  get_table("ind_index") %>%
  mutate_table(
    index_date = as.Date(NA),     # add gen of evaluation
    rep = NA_integer_           # add rep
  )

pop %>% get_table("ind_index")
pop %>% get_table("index_meta")

#------------------------------------------------------------------------------#
# Add Founders
#------------------------------------------------------------------------------#

#for (repl in 1:config$general$n_reps){
repl = 1
  
warning("\n ----------    REPLICATE: ", repl, "    --------------------\n")

# sample birth dates of founders
min_birth_date <- as.Date(config$general$start_date) - config$testing$off_test_age
max_birth_date <- as.Date(config$general$start_date) + 7
sampled_birth_dates <- min_birth_date + 
  round(runif(n = (config$population$n_founder_male + config$population$n_founder_female),
         min = 0, max = config$testing$off_test_age + 7))

hist(sampled_birth_dates, breaks = 50)

# add founders
pop <- pop %>%
  add_founders(          # add founders
    n_males   = config$population$n_founder_male,      # sample male founders
    n_females = config$population$n_founder_female,     # sample female founders
    line_name = "Libra",         # name this line
    rep       = repl,            # USER DEFINED - replicate
    conc_date     = sampled_birth_dates - config$general$gest_len, 
    birth_date    = sampled_birth_dates, 
    on_test_date  = sampled_birth_dates + config$testing$on_test_age,
    off_test_date = sampled_birth_dates + config$testing$off_test_age,
    alive     = TRUE,
    active    = FALSE
  )

#------------------------------------------------------------------------------#
# Add SNP Chip
#------------------------------------------------------------------------------#

warning("Add 9k Chip")

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

warning("Add cov matrices")

# additive genetic (CO)VARIANCES
mat.add.gen.vars <- matrix(c(200.00,  0.00,    0.00,  0.00,
                               0.00,  0.90,    3.07,  0.21,
                               0.00,  3.07, 1050.00, 17.80,
                               0.00,  0.21,   17.80,  1.20), 
                      nrow = 4, byrow=TRUE, 
                      dimnames = list(c("AP", "NW", "ADG", "BF"), 
                                      c("AP", "NW", "ADG", "BF")))

# residual (CO)VARIANCES
mat.res.vars <- matrix(c(400,  0.00,    0.00,  0.00,
                         0.0,  8.10,   10.00,  0.65,
                         0.0, 10.00, 2100.00, 10.00,
                         0.0,  0.65,   10.00,  1.30), 
                      nrow = 4, byrow=TRUE, 
                      dimnames = list(c("AP", "NW", "ADG", "BF"), 
                                      c("AP", "NW", "ADG", "BF")))

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

warning("Define index")

# add generic terminal index
pop %>%
  define_index(
    index_name  = "maternal",
    trait_names = c("AP", "NW", "ADG", "BF"),
    index_wts   = c(-14, 93, 1.5, -30)
  )

# print table with index values
pop %>% get_table("index_meta")

#------------------------------------------------------------#
# Trait: AP - Age at Puberty
#------------------------------------------------------------#

warning("Define AP - Age at Puberty")

# add ADG as a trait
pop <- pop %>%
  add_trait(
    trait_name  = "AP",
    description = "Age at Puberty", 
    units       = "days",                  # grams per day during testing period
    trait_type  = "count",           # not categorical or binary
    repeatable  = FALSE,                  # only 1 phenotype per individual
    recorded_on = "self",                 # recorded on itself only
    target_add_mean = 0,      # mean TBV in 'base'
    #target_add_var  = 1082,  # additive variance target in 'base'
    #residual_var    = 2500,  # residual variance (target h2 ~0.30) - 'fixed'
    index_weight    = -14,   # index weight (not implemented yet)
    economic_value  = -14,   # economic value (not implemented yet)
    overwrite       = TRUE    # wipe this row if it exists and replace with this new data
  )

# add which loci are QTL and their effects
pop %>%
  # set loci as QTL for this trait
  define_qtl(
    trait_name = "AP",         # name of trait
    n = 1000,                   # set number of QTL for this trait (will set locations)
    method = "random"           # allocate to loci randomly
  ) %>%
  # set QTL effects
  add_additive_effects(
    trait_name      = "AP",        # trait name
    distribution    = "normal",     # distribution of QTL effects
    scale_to_target = TRUE,         # scale to meet additive variance target
    base            = "current_pop" # use all animals in pop to standardized (or if filtered)
  )

# add all TBV for AP
pop <- pop %>%
  get_table("ind_meta") %>% # here we specify the 'ind_meta' table so all animals will have their TBV calculated
  filter(rep == repl) %>%
  add_tbv("AP", rep = repl)

pop %>% get_table("ind_tbv")

# add overall mean for AP
pop %>%
  add_effect_int(
    trait_name = "AP",              # trait (need to change to "trait_name")
    mean = config$general$mean_puberty_age
  )

# test `add_phenotype()` function
pop %>%
  get_table("ind_meta") %>%     # will phenotype all individuals in this table with no filter
  filter(
    rep == repl,
    sex == "F",
  ) %>%
  add_phenotype(                # add rows to 'ind_phenotye' table
    trait = "AP",              # trait name
    rep = repl                 # set rep number
  )  

data.age.puberty <- pop %>%
  get_table("ind_phenotype") %>%
  filter(
    rep == repl,
    trait_name == "AP"
  ) %>%
  select(id_ind, value) %>%
  collect()

data.birth.date <- pop %>%
  get_table("ind_meta") %>%
  filter(
    rep == repl,
    sex == "F"
  ) %>%
  select(id_ind, birth_date) %>%
  collect() %>%
  left_join(., data.age.puberty) %>%
  mutate(
    pheno_date = birth_date + value
  )

# add phenotype date to phenotype table
pop <- pop %>%
  get_table("ind_meta") %>%
  filter(
    rep == repl,
    sex == "F"
  ) %>%
  mutate_table(
    puberty_date = data.birth.date$pheno_date
  )

# add phenotype date to phenotype table
pop <- pop %>%
  get_table("ind_phenotype") %>%
  filter(
    rep == repl,
    trait_name == "AP"
  ) %>%
  mutate_table(
    pheno_date = data.birth.date$pheno_date
  )

pop %>% get_table("ind_phenotype") %>% filter(trait_name == "AP")

#dbExecute(pop$db_conn, "DELETE FROM ind_phenotype WHERE trait_name = 'AP'")

#------------------------------------------------------------#
# Trait: ADG
#------------------------------------------------------------#

warning("Define ADG")

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
  add_additive_effects(
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
  filter(
    rep == repl,
    off_test_date < start_date
    ) %>%
  add_phenotype(                # add rows to 'ind_phenotye' table
    trait = "ADG",              # trait name
    rep = repl                  # set rep number
  )  

pop %>% get_table("ind_phenotype") %>% filter(trait_name == "ADG")

#------------------------------------------------------------#
# Trait: BF
#------------------------------------------------------------#

warning("Define BF")

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
  add_additive_effects(
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
  filter(
    rep == repl,
    off_test_date < start_date
  ) %>%
  add_phenotype(                # add rows to 'ind_phenotye' table
    trait = "BF",               # trait name
    rep = repl                  # add rep
  )

pop %>% get_table("ind_phenotype") %>% filter(trait_name == "BF")

#------------------------------------------------------------#
# Trait: NW
#------------------------------------------------------------#

warning("Define NW")

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
  add_additive_effects(
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
# pop %>%
#   get_table("ind_meta") %>%     # will phenotype all individuals in this table with no filter
#   filter(
#     rep == repl, 
#     sex == "F"
#   ) %>%        # only females get NW phenotype
#   add_phenotype(                 # add rows to 'ind_phenotye' table
#     trait = "NW",                # trait name
#     rep = repl                   # add rep number
#   )

pop %>% get_table("ind_phenotype") %>% filter(trait_name == "NW")

#------------------------------------------------------------------------------#
# Run EBVs
#------------------------------------------------------------------------------#

# warning("Run founder EBVs")
# 
# # run EBV for AP
# pop <- pop %>%
#   get_table("ind_meta") %>%
#   filter(rep == repl) %>%
#   add_ebv("AP", software="blupf90", model="blup", 
#           eval_date = as.Date(config$general$start_date), rep = repl)
# 
# pop %>% get_table("ind_ebv") %>%
#   filter(trait_name == "AP", rep == repl)
# 
# # run EBV for ADG
# pop <- pop %>%
#   get_table("ind_meta") %>%
#   filter(rep == repl) %>%
#   add_ebv("ADG", software="blupf90", model="blup",
#           eval_date = as.Date(config$general$start_date), rep = repl)
# 
# pop %>% get_table("ind_ebv") %>%
#   filter(trait_name == "ADG", rep == repl)
# 
# # run EBV for BF
# pop <- pop %>%
#   get_table("ind_meta") %>%
#   filter(rep == repl) %>%
#   add_ebv("BF", software="blupf90", model="blup", 
#           eval_date = as.Date(config$general$start_date), rep = repl)
# 
# pop %>% get_table("ind_ebv") %>%
#   filter(trait_name == "BF", rep == repl)
# 
# # run EBV for NW
# pop <- pop %>%
#   get_table("ind_meta") %>%
#   filter(rep == repl) %>%
#   add_ebv("NW", software="blupf90", model="blup", 
#           eval_date = as.Date(config$general$start_date), rep = repl)
# 
# pop %>% get_table("ind_ebv") %>%
#   filter(trait_name == "NW", rep == repl)

#------------------------------------------------------------------------------#
# Run Index
#------------------------------------------------------------------------------#

# warning("Calculate founder Indexes")
# 
# # run index calculation
# pop %>%
#   get_table("ind_ebv") %>%    # must pass 'ind_ebv' because it contains the EBVs needed
#   filter(
#     rep == repl, 
#     eval_date == as.Date(config$general$start_date)) %>%
#   add_index(
#     "maternal",          # just give the index name and it will grab weights
#     index_date = as.Date(config$general$start_date),
#     rep = repl
#   )
# 
# pop %>% get_table("ind_index") %>%
#   filter(index_name == "maternal", rep == repl)
# 
# pop %>% get_table("ind_index") %>%
#   filter(index_name == "maternal", rep == repl) %>%
#   collect() %>%
#   count()

#------------------------------------------------------------------------------#
# Checks
#------------------------------------------------------------------------------#

warning("Calculate Means in Founder Generation")

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

# # print mean EBVs
# pop %>%
#   get_table("ind_ebv") %>%
#   collect() %>%
#   group_by(rep, trait_name) %>%
#   summarise(
#     MeanEBV = mean(ebv),
#     .groups = "drop_last"
#   ) %>%
#   print(n=10)

#------------------------------------------------------------#
# Phenotypes
#------------------------------------------------------------#

# print mean phenotype by trait
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
# Run by Time
#------------------------------------------------------------------------------#

if (1 > 2){

  # delete rows in tables if needed!  
  dbExecute(pop$db_conn, "DELETE FROM ind_phenotype WHERE trait_name = 'AP'")
  dbExecute(pop$db_conn, "DELETE FROM ind_phenotype WHERE trait_name = 'ADG'")
  dbExecute(pop$db_conn, "DELETE FROM ind_phenotype WHERE trait_name = 'BF'")
  dbExecute(pop$db_conn, "DELETE FROM ind_phenotype WHERE trait_name = 'NW'")
  # delete rows in tables if needed!  
  dbExecute(pop$db_conn, "DELETE FROM ind_ebv WHERE trait_name = 'AP'")
  dbExecute(pop$db_conn, "DELETE FROM ind_ebv WHERE trait_name = 'ADG'")
  dbExecute(pop$db_conn, "DELETE FROM ind_ebv WHERE trait_name = 'BF'")
  dbExecute(pop$db_conn, "DELETE FROM ind_ebv WHERE trait_name = 'NW'")
  # delete index rows
  dbExecute(pop$db_conn, "DELETE FROM ind_index WHERE index_name = 'maternal'")
  # delete new offspring
  dbExecute(pop$db_conn, "DELETE FROM ind_meta WHERE birth_date > '2026-08-01'")
}

# set starting selection date
selection_date <- as.Date("2026-01-29") + 28 

#------------------------------------------------------------------------------#
# Run by Time
#------------------------------------------------------------------------------#

warning("Begin Date Loop")

# big date "loop" but continuous
for (cur_date in seq(as.Date(start_date), as.Date("2026-08-01"))){

# convert to date, not integer
cur_date = as.Date(cur_date)

# get day of week (run evals on certain days for instance)
cur_day_of_week = weekdays(cur_date)

# print "current date" (simulated advance by 1 day)
message("Current Date: ", cur_date)

# ------------------------------ PHENOTYPE ----------------------------------- #

# ---------- AP ---------- #

pop %>%
  get_table("ind_meta") %>%
  filter(
    rep == repl,
    puberty_date  == cur_date
  ) %>%
  add_phenotype("AP", rep = repl)

# ---------- ADG ---------- #

pop %>%
  get_table("ind_meta") %>%
  filter(
    rep == repl,
    off_test_date  == cur_date
  ) %>%
  add_phenotype("ADG", rep = repl)

# ---------- BF ---------- #

pop %>%
  get_table("ind_meta") %>%
  filter(
    rep == repl,
    off_test_date  == cur_date
  ) %>%
  add_phenotype("BF", rep = repl)

# ---------- NW ---------- #

pop %>%
  get_table("ind_meta") %>%
  filter(
    rep == repl,
    puberty_date  == cur_date
  ) %>%
  add_phenotype("NW", rep = repl)


# ------------------------------ GENOTYPE ------------------------------------ #


# ------------------------------ RUN EVALUATIONS ----------------------------- #

if (cur_date > as.Date("2026-04-01") & cur_day_of_week == "Friday"){
  
  message("Calculate EBVs by Trait")

  # run EBV for AP
  pop <- pop %>%
    get_table("ind_meta") %>%
    filter(rep == repl) %>%
    add_ebv("AP", software="blupf90", model="blup",
            eval_date = cur_date, rep = repl)
  
  pop %>% get_table("ind_ebv") %>%
    filter(trait_name == "AP", rep == repl)
  
  # run EBV for ADG
  pop <- pop %>%
    get_table("ind_meta") %>%
    filter(rep == repl) %>%
    add_ebv("ADG", software="blupf90", model="blup",
            eval_date = cur_date, rep = repl)
  
  pop %>% get_table("ind_ebv") %>%
    filter(trait_name == "ADG", rep == repl)
  
  # run EBV for BF
  pop <- pop %>%
    get_table("ind_meta") %>%
    filter(rep == repl) %>%
    add_ebv("BF", software="blupf90", model="blup",
            eval_date = cur_date, rep = repl)
  
  pop %>% get_table("ind_ebv") %>%
    filter(trait_name == "BF", rep == repl)
  
  # run EBV for NW
  pop <- pop %>%
    get_table("ind_meta") %>%
    filter(rep == repl) %>%
    add_ebv("NW", software="blupf90", model="blup",
            eval_date = cur_date, rep = repl)

  pop %>% get_table("ind_ebv") %>%
    filter(trait_name == "NW", rep == repl)

} # END CALCULATE EBVs

# ------------------------------ CALC INDEX ---------------------------------- #

if (cur_date > as.Date("2026-04-01") & cur_day_of_week == "Friday"){
  
  message("Calculate Indexes")
  
  # run index calculation
  pop %>%
    get_table("ind_ebv") %>%    # must pass 'ind_ebv' because it contains the EBVs needed
    filter(
      rep == repl,
      eval_date == cur_date
    ) %>%
    add_index(
      "maternal",          # just give the index name and it will grab weights
      index_date = cur_date,
      rep = repl
    )
} # END INDEX CALCULATION


# ------------------------------ SELECT -------------------------------------- #

if (cur_date == selection_date){

  message("Selection Date!")
  
  boar_cull_date <- cur_date - 50
  gilt_cull_date <- cur_date - 75
  
  # | status == "breeding-boar"
  #  & off_test_date > (cur_date - 50)
  male_candidates <- pop %>%
    get_table("ind_meta") %>%
    filter(
      rep == repl,
      sex == "M",
      cur_date > off_test_date & off_test_date > boar_cull_date | status == "breeding-boar"
    ) %>%
    pull(id_ind)
  
  message("male candidates:")
  print(male_candidates)
  
  female_candidates <- pop %>%
    get_table("ind_meta") %>%
    filter(
      rep == repl,
      sex == "F",
      cur_date > puberty_date,
      #off_test_date > gilt_cull_date,
      #status != "gestation"
    ) %>%
    pull(id_ind)
  
  message("female candidates:")
  print(female_candidates)
  
  # random selection early on before EBVs running
  if (cur_date < as.Date("2026-04-01")){
    
    message("Using RANDOM")
    
    selected_males <- pop %>%
      get_table("ind_meta") %>%
      filter(
        id_ind %in% male_candidates
      ) %>%
      slice_sample(n=3) %>%
      pull(id_ind)
    
    message("selected males:")
    print(selected_males)
  
    selected_females <- pop %>%
      get_table("ind_meta") %>%
      filter(
        id_ind %in% female_candidates
      ) %>%
      slice_sample(n=10) %>%
      pull(id_ind)
    
    message("selected females:")
    print(selected_females)
    
  # index selection after EBVs are running
  } else {
    
    message("Using INDEX")
    
    latest_index_date = pop %>%
      get_table("ind_index") %>%
      collect() %>%
      pull(index_date) %>%
      max()
    
    selected_males <- pop %>%
      get_table("ind_index") %>%
      filter(
        id_ind %in% male_candidates,
        index_date == latest_index_date
      ) %>%
      slice_sample(n=3) %>%
      pull(id_ind)
    
    message("selected males:")
    print(selected_males)
  
    selected_females <- pop %>%
      get_table("ind_index") %>%
      filter(
        id_ind %in% female_candidates,
        index_date == latest_index_date
      ) %>%
      slice_sample(n=10) %>%
      pull(id_ind)
    
    message("selected females:")
    print(selected_females)
  
  } # end index selection
  
# set males to boar
pop %>% 
  get_table("ind_meta") %>%
  filter(
    id_ind %in% selected_males
  ) %>%
  mutate_table(
    status = "breeding-boar"
  )

# set females to gest
pop %>% 
  get_table("ind_meta") %>%
  filter(
    id_ind %in% selected_females
  ) %>%
  mutate_table(
    status = "gestation"
  )


} # END SELECTION STEP

# ------------------------------ MATE ---------------------------------------- #

if (cur_date == selection_date){
  
  #---------------- SAMPLE NW PHENOTYPE ON DAMS ----------------#
  
  message("Sample NW phenotype on selected dams")
  
  # phenotype selected dams first
  pop %>%
    get_table("ind_meta") %>%
    filter(
      rep == repl, 
      id_ind %in% selected_females
    ) %>%
    add_phenotype(
      "NW",
      rep = repl,          # add rep number
      pheno_date = cur_date + config$general$gest_len
    )
  
  #---------------- EXTRACT NW PHENOTYPE ----------------#
  
  message("Extract NW phenotype")
  
  # extract NW phenotype to produce offspring numbers correctly
  data.nw <- pop %>%
    get_table("ind_phenotype") %>%
    filter(
      rep == repl, 
      trait_name == "NW",
      id_ind %in% selected_females,
      pheno_date == cur_date + config$general$gest_len
    ) %>%
    collect() %>%
    select(id_ind, value)
    
  #---------------- BUILD PROGENY MATRIX ----------------#
  
  message("Build progeny matrix")
  
  # HOW:
  # boars sampled randomly, 1 per dam (1 mating per sire/dam pair)
  # offspring repeated based on the number of "NW" sampled above
  
  # sample males randomly so 1 sire per dam (1 mating but multiple offspring)
  list_selected_males <- sample(selected_males, size=nrow(data.nw), replace=TRUE)
  
  # use new phenotype to build mating plan
  data.new.matings <- tibble(
    # rep sires by "NW" phenotype from dam so they match the same rows (1 litter)
    id_parent_1 = rep(list_selected_males, time=data.nw$value),
    # rep dams by "NW" phenotype to get a full litter (1 row / offspring)
    id_parent_2 = rep(c(data.nw$id_ind), time=data.nw$value),
    line     = "Libra",
    sex      = sample(c("M", "F"), size=sum(data.nw$value), replace=TRUE, prob=c(0.5, 0.5)),
    rep      = repl,
    conc_date     = cur_date,
    birth_date    = cur_date + config$general$gest_len,
    on_test_date  = cur_date + config$general$gest_len + config$testing$on_test_age,
    off_test_date = cur_date + config$general$gest_len + config$testing$off_test_age
  )

  #---------------- ADD OFFSPRING ----------------#
  
  message("Add offspring")
  
  # add new offspring based on tibble mating plan (1 row per offspring)
  pop %>%
    add_offspring(
      data.new.matings
    )

  #---------------- UPDATE DATES AND STATUS ----------------#
  
  message("Add latest mating date")
  
  # convert to dates just in case
  cur_farrow_date = as.Date(cur_date + config$general$gest_len)
  cur_wean_date   = as.Date(cur_date + config$general$gest_len + config$general$lact_len)
  
  # set dates after mating
  pop %>%
    get_table("ind_meta") %>%
    filter(
      rep == repl,
      id_ind %in% selected_females
    ) %>%
    mutate_table(
      mate_date   = format(cur_date, "%Y-%m-%d"),        # add mating date
      farrow_date = format(cur_farrow_date, "%Y-%m-%d"), # add farrowing date
      wean_date   = format(cur_wean_date, "%Y-%m-%d")    # add weaning date
    )
  
} # # END MATING STEP


# increase selection_date (every 28 days)
if (cur_date == selection_date){
  selection_date = selection_date + 28
}

} # END SIMULATION LOOP FOR DATE








#------------------------------------------------------------------------------#
# Data Summary
#------------------------------------------------------------------------------#

# count status by sex
pop %>%
  get_table("ind_meta") %>%
  collect() %>%
  group_by(sex, status) %>%
  summarise(
    Status = n()
  )

# calculate mean EBVs by trait
pop %>%
  get_table("ind_ebv") %>%
  collect() %>%
  group_by(trait_name) %>%
  summarise(
    MeanEBV = mean(ebv)
  )

# calculate mean phenotypes by trait
pop %>%
  get_table("ind_phenotype") %>%
  collect() %>%
  group_by(trait_name) %>%
  summarise(
    MeanPhenotype = mean(value)
  )

# histogram of EBVs by Trait
pop %>%
  get_table("ind_ebv") %>%
  collect() %>%
  filter(eval_date == max(eval_date)) %>%
ggplot(aes(x=ebv)) +
  geom_histogram(color="white", fill=hg_blue) +
  facet_wrap(~trait_name, scale="free") +
  labs(
    title = "EBVs (latest)"
  )

# histogram of phenotypes by Trait
pop %>%
  get_table("ind_phenotype") %>%
  collect() %>%
ggplot(aes(x=value)) +
  geom_histogram(color="white", fill=hg_green, bins=17) +
  facet_wrap(~trait_name, scale="free_x") +
  labs(
    title = "Phenotypes"
  )

# histogram of indexes
pop %>% get_table("ind_index") %>%
  filter(
    index_name == "maternal", 
    rep == repl, 
    index_date == as.Date("2026-06-12")
  ) %>%
  collect() %>%
ggplot(aes(x=index_value)) +
  geom_histogram(color="white", fill="aquamarine3") +
  labs(
    title = "Index Values (latest)"
  )

# count rows in index table
pop %>% get_table("ind_index") %>%
  filter(index_name == "maternal", rep == repl) %>%
  collect() %>%
  group_by(index_name, rep, index_date) %>%
  count()







#------------------------------------------------------------------------------#
# Summarize Pop - Tables and Figures
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
# plot options
#------------------------------------------------------------#

trait_colors = c("cadetblue2", "darkolivegreen2", "darkorange", 
                 "darkorchid", "dodgerblue2")

#------------------------------------------------------------#
# Genetic Variances
#------------------------------------------------------------#

data.gen.vars <- pop %>%
  get_table("trait_effect_cov") %>%
  filter(
    effect_name == "gen_add",
    trait_1 == trait_2
  ) %>%
  collect() %>%
  select(trait_1, cov) %>%
  rename(trait_name = trait_1, var = cov) %>%
  mutate(
    sd = sqrt(var)
  )

#------------------------------------------------------------#
# EBVs
#------------------------------------------------------------#

# plot EBV Trends
p <- data.tbv.ebv %>%
  group_by(rep, gen_born, trait_name) %>%
  summarise(
    MeanEBV = mean(ebv),
    .groups = "drop_last"
  ) %>%
  mutate(
    rep = as.factor(rep),
    gen_born = as.factor(gen_born)
  ) %>%
ggplot(., aes(x=gen_born, y=MeanEBV, color=trait_name, shape=as.factor(rep), group=rep)) +
  geom_hline(yintercept = 0, color="grey50", linetype=2) +
  geom_line() +
  geom_point() +
  facet_wrap(~trait_name, scales="free_y") +
  scale_shape_discrete("Rep Number") +
  #scale_color_manual("Trait", values = trait_colors) +
  scale_color_discrete("Trait") +
  scale_x_discrete(limits=factor(0:10)) +
  labs(
    title = "Estimated Breeding Values (EBVs)",
    x = "Generation Born",
    y = "Mean Estimated Breeding Value (EBV)",
    caption = config$scenario_name
  )

print(p)

ggsave(
  filename = file.path(config$output$save_dir, "ebv_mean_on_gen_born_facet_trait_color_trait.png"),
  plot = p, 
  width = 8,
  height = 5,
  units = "in",
  dpi = 100,
  bg = "white"
)

#------------------------------------------------------------#
# TBVs
#------------------------------------------------------------#

# plot EBV Trends
p <- data.tbv.ebv %>%
  group_by(rep, gen_born, trait_name) %>%
  summarise(
    MeanTBV = mean(tbv),
    .groups = "drop_last"
  ) %>%
  mutate(
    rep = as.factor(rep),
    gen_born = as.factor(gen_born)
  ) %>%
ggplot(., aes(x=gen_born, y=MeanTBV, color=trait_name, shape=as.factor(rep), group=rep)) +
  geom_hline(yintercept = 0, color="grey50", linetype=2) +
  geom_line() +
  geom_point() +
  facet_wrap(~trait_name, scales="free_y") +
  scale_shape_discrete("Rep Number") +
  #scale_color_manual("Trait", values = trait_colors) +
  scale_color_discrete("Trait") +
  scale_x_discrete(limits=factor(0:10)) +
  labs(
    title = "True Breeding Values (TBVs)",
    x = "Generation Born",
    y = "Mean True Breeding Value (TBV)",
    caption = config$scenario_name
  )

print(p)

ggsave(
  filename = file.path(config$output$save_dir, "tbv_mean_on_gen_born_facet_trait_color_trait.png"),
  plot = p, 
  width = 8,
  height = 5,
  units = "in",
  dpi = 100,
  bg = "white"
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
p <- data.tbv.ebv %>%
  # mutate(
  #   rep = as.factor(rep),
  #   gen_born = as.factor(gen_born)
  # ) %>%
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
    caption = config$scenario_name
  )

print(p)

ggsave(
  filename = file.path(config$output$save_dir, "ribbon_ebv_on_gen_born_facet_trait.png"),
  plot = p, 
  width = 8,
  height = 5,
  units = "in",
  dpi = 100,
  bg = "white"
)

#------------------------------------------------------------#
# TBVs
#------------------------------------------------------------#

# plot TBVs
p <- data.tbv.ebv %>%
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
    caption = config$scenario_name
  )

print(p)

ggsave(
  filename = file.path(config$output$save_dir, "ribbon_tbv_on_gen_born_facet_trait.png"),
  plot = p, 
  width = 8,
  height = 5,
  units = "in",
  dpi = 100,
  bg = "white"
)

#------------------------------------------------------------#
# EBVs + TBVs
#------------------------------------------------------------#

# plot TBVs
p <- data.tbv.ebv %>%
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
    caption = config$scenario_name
  )

print(p)

ggsave(
  filename = file.path(config$output$save_dir, "ribbon_ebv_and_tbv_on_gen_born_facet_trait.png"),
  plot = p, 
  width = 8,
  height = 5,
  units = "in",
  dpi = 100,
  bg = "white"
)




#------------------------------------------------------------------------------#
# Summary of Data
#------------------------------------------------------------------------------#

#------------------------------------------------------------#
# Animals - Rep + Gen
#------------------------------------------------------------#

# count by rep and gen
p <- pop %>%
  get_table("ind_meta") %>%
  collect() %>%
ggplot(., aes(x=as.factor(gen_born), fill=as.factor(rep))) +
  geom_bar(position="dodge") +
  scale_x_discrete(limits = factor(0:10)) +
  #scale_fill_manual("Replicate", values=trait_colors)) +
  scale_fill_discrete("Replicate") +
  labs(
    title = "Count Animal Records by Rep + Generation",
    x = "Generation Born",
    y = "Count",
    caption = config$scenario_name
  )

print(p)

ggsave(
  filename = file.path(config$output$save_dir, "bar_animal_count_on_gen_born_fill_rep.png"),
  plot = p, 
  width = 8,
  height = 5,
  units = "in",
  dpi = 100,
  bg = "white"
)

#------------------------------------------------------------#
# Phenotypes - Rep + Gen
#------------------------------------------------------------#

# count by rep and gen
p <- pop %>%
  get_table("ind_phenotype") %>%
  collect() %>%
ggplot(., aes(x=as.factor(gen_pheno), fill=as.factor(rep))) +
  geom_bar(position="dodge") +
  facet_wrap(~trait_name) +
  scale_x_discrete(limits = factor(0:10)) +
  #scale_fill_manual("Replicate", values=trait_colors) +
  scale_fill_discrete("Replicate") +
  labs(
    title = "Phenotype Count by Rep + Generation",
    x = "Generation Born",
    y = "Phenotype Count",
    caption = config$scenario_name
  )

print(p)

ggsave(
  filename = file.path(config$output$save_dir, "bar_phenotype_count_on_gen_born_fill_rep.png"),
  plot = p, 
  width = 8,
  height = 5,
  units = "in",
  dpi = 100,
  bg = "white"
)

#------------------------------------------------------------#
# TBV - Rep + Gen
#------------------------------------------------------------#

# count by rep and gen
p <- pop %>%
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
  geom_hline(aes(yintercept = 0), color="red3", linewidth=0.75, linetype=3) +
  geom_ribbon(aes(ymin=MinTBV, ymax=MaxTBV), alpha=0.05) +
  geom_ribbon(aes(ymin=Q1TBV, ymax=Q3TBV), alpha=0.40) +
  geom_line(aes(y=Q2TBV)) +
  facet_wrap(~ trait_name, scale="free_y") +
  #scale_color_manual("Replicate", values=trait_colors) +
  scale_color_discrete("Replicate") +
  #scale_fill_manual("Replicate", values=trait_colors) +
  scale_fill_discrete("Replicate") +
  labs(
    title = "TBV Trends By Trait and Rep",
    subtitle = "Median Trend + middle 50 percentile + min/max",
    x = "Generation Born",
    y = "TBV (by rep)",
    caption = config$scenario_name
  )

print(p)

ggsave(
  filename = file.path(config$output$save_dir, "ribbon_tbv_on_gen_born_facet_trait_fill_trait.png"),
  plot = p, 
  width = 8,
  height = 5,
  units = "in",
  dpi = 100,
  bg = "white"
)







# close pop object for database
close_pop(pop)

