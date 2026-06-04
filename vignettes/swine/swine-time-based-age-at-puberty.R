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
  pak::pak("austin-putz/tidybreed", upgrade=TRUE)
}

library(DBI)
library(glue)
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

data.config <- summary(config) %>% as.data.frame()
summary(config)

# create output directory
if (dir.exists(config$output$save_dir)){
  warning("Output directory already exists")
} else {
  warning("Output directory will be created for user")
  # create directory if it doesn't exist
  dir.create(config$output$save_dir, recursive=TRUE, showWarnings = TRUE)
}

#------------------------------------------------------------------------------#
# general inputs
#------------------------------------------------------------------------------#

# format for printing time
format_elapsed <- function(pt_diff) {
  secs <- pt_diff["elapsed"]
  if (secs < 60) {
    sprintf("%.1f sec", secs)
  } else {
    sprintf("%.1f min", secs / 60)
  }
}

# start time
time_start_total <- proc.time()

# start date of simluation
start_date <- as.Date(config$general$start_date)
end_date   <- as.Date(config$general$end_date)

# SEPARATE BOAR AND GILT/SOW SELECTION STEPS

# set starting selection date
female_selection_date <- as.Date(config$general$start_date_selection)
male_selection_date   <- as.Date(config$general$start_date_selection)

# add to data frame
data.timing <- tibble(
  sim_date       = as.Date(character()),
  real_date_time = as.POSIXct(character()),
  type           = character(),
  elapsed_sec    = numeric(),
  cumulative_sec = numeric()
)

# add to data frame
data.timing <- add_row(data.timing,
  sim_date       = start_date,
  real_date_time = Sys.time(),
  type           = "begin-simulation",
  elapsed_sec    = NA,
  cumulative_sec = NA
)

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
    rep = NA_integer_,           # add rep to phenotypes
    current_date = as.Date(NA),  # add current date (to subset ind later)
    pheno_date = as.Date(NA),    # add phenotype date
    .set_default = TRUE          # if TRUE, sets default of given value
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

# set min and max birth dates
min_birth_date <- as.Date(config$general$start_date) - config$general$mean_puberty_age
max_birth_date <- as.Date(config$general$start_date_selection) + 
                    config$general$gest_len +
                    config$general$lact_len + 
                    config$general$w2e_int + 
                    config$testing$off_test_age 

days_between_founders <- as.numeric(max_birth_date - min_birth_date)
  
# sample birth dates
sampled_birth_dates <- min_birth_date + 
  round(runif(n = (config$population$n_founder_male + config$population$n_founder_female),
         min = 0, max = days_between_founders))

# hist of birth dates
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
  get_table("genome_meta") %>%
  slice_sample(n=9000) %>%
  define_chip(chip_name = "9k")

pop %>% get_table("genome_meta")

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
  define_effect_cov_matrix(
    effect_name = "gen_add",          # fixed term for additive genetic (co)variance matrix
    cov_matrix  = mat.add.gen.vars    # name of matrix with row/col names
  )

# add this residual covariance matrix to a table with function
pop <- pop %>%
  define_effect_cov_matrix(
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
    index_name   = "maternal",
    trait_names  = c("AP", "NW", "ADG", "BF"),
    index_wts    = c(-3, 93, 1.5, -30),
    economic_wts = c(0, 80, 2.0, -20)
  )

# print table with index values
pop %>% get_table("index_meta")

#------------------------------------------------------------#
# Trait: AP - Age at Puberty
#------------------------------------------------------------#

warning("Define AP - Age at Puberty")

# add ADG as a trait
pop <- pop %>%
  define_trait(
    trait_name  = "AP",
    description = "Age at Puberty", 
    units       = "days",                  # grams per day during testing period
    type        = "count",           # not categorical or binary
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
  get_table("genome_meta") %>%
  filter(is_9k != TRUE) %>%
  define_additive_effects(
    trait_name = "AP", 
    distribution = "normal", 
    scale_to_target = TRUE, 
    base = "current_pop"
  )

# print 'genome_effects'
pop %>% get_table("genome_effects")

# add all TBV for AP
pop <- pop %>%
  get_table("ind_meta") %>% # here we specify the 'ind_meta' table so all animals will have their TBV calculated
  filter(rep == repl) %>%
  add_tbv("AP", rep = repl)

pop %>% get_table("ind_tbv")

pop %>%
  get_table("ind_meta") %>%
  add_genotypes(chip_name = "9k")

pop %>%
  get_table("ind_meta") %>%
  extract_genotypes(chip_name = "9k")

pop %>%
  get_table("ind_meta") %>%
  filter(
    id_ind %in% c("Libra_1", "Libra_2")
  ) %>%
  extract_genotypes(
    effects_tbl = pop %>% get_table("genome_effects") %>% filter(trait_name=="AP")
  ) %>%
  collect()


# ----- ADD OVERALL MEAN ----- #

# add overall mean for AP
pop %>%
  define_effect_intercept(
    trait_name = "AP",              # trait (need to change to "trait_name")
    mean = config$general$mean_puberty_age
  )

# ----- SAMPLE AP PHENOTYPE ON FEMALES ----- #

# sample phenotype for all founder females
pop %>%
  get_table("ind_meta") %>%     # will phenotype all individuals in this table with no filter
  filter(
    rep == repl,
    sex == "F",
  ) %>%
  add_phenotype(                # add rows to 'ind_phenotye' table
    trait_name = "AP",              # trait name
    rep = repl                 # set rep number
  )

pop %>% get_table("ind_phenotype")

# ----- EXTRACT PHENOTYPE ----- #

# pull phenotype sampled
data.age.puberty <- pop %>%
  get_table("ind_phenotype") %>%
  filter(
    rep == repl,
    trait_name == "AP"
  ) %>%
  select(id_ind, pheno_value) %>%
  collect()

list_founder_AP_ids <- as.character(data.age.puberty$id_ind)

# ----- ADD PUBERTY DATE ----- #

data.birth.date <- pop %>%
  get_table("ind_meta") %>%
  filter(
    rep == repl,
    id_ind %in% list_founder_AP_ids
  ) %>%
  select(id_ind, birth_date) %>%
  collect() %>%
  left_join(., data.age.puberty) %>%
  mutate(
    pheno_date = birth_date + pheno_value
  )

# ----- VERY IMPORTANT! DON'T MESS UP ORDER! ----- #

# check if order is the same
if (all(list_founder_AP_ids == data.birth.date$id_ind) == FALSE){
  stop("IDs don't match up!!!")
} else {
  message("IDs match up")
}

# ----- UPDATE 'puberty_date' in 'ind_meta' ----- #

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

# ----- UPDATE 'pheno_date' in 'ind_phenotype' ----- #

# add phenotype date to phenotype table
pop <- pop %>%
  get_table("ind_phenotype") %>%
  filter(
    rep == repl,
    trait_name == "AP",
    id_ind %in% list_founder_AP_ids
  ) %>%
  mutate_table(
    pheno_date = data.birth.date$pheno_date
  )

# ----- CHECK 'ind_phenotype' ----- #

# print table
pop %>% get_table("ind_phenotype") %>% filter(trait_name == "AP")

# count how many rows are before the start date
pop %>% 
  get_table("ind_phenotype") %>% 
  filter(trait_name == "AP", pheno_date < start_date) %>% 
  collect() %>% 
  count()

# ------------------------------ UPDATE STATUS: FEMALES ---------------------- #

message("Change Status - 'after-test-gilt'")

# change status based on reaching off-test
pop %>%
  get_table("ind_meta") %>%
  filter(
    rep == repl,
    sex == "F",
    off_test_date < start_date
  ) %>%
  mutate_table(
    status = "after-test-gilt"
  )

message("Change Status - 'puberty-gilt'")

# change status based on reaching puberty
pop %>%
  get_table("ind_meta") %>%
  filter(
    rep == repl,
    sex == "F",
    puberty_date < start_date
  ) %>%
  mutate_table(
    status = "puberty-gilt"
  )

message("Change Status - 'cull-gilt'")

# create new gilt cull date
cur_gilt_cull_date <- start_date - config$culling$gilt_cull_days_after_off_test

# change status to cull if gilt is not selected
pop %>%
  get_table("ind_meta") %>%
  filter(
    rep == repl,
    sex == "F",
    status == "after-test-gilt" | status == "puberty-gilt",
    off_test_date < cur_gilt_cull_date       # paste X days, cull them if not in puberty
  ) %>%
  mutate_table(
    status = "cull-gilt"
  )

#------------------------------------------------------------#
# Trait: ADG
#------------------------------------------------------------#

warning("Define ADG")

# add ADG as a trait
pop <- pop %>%
  define_trait(
    trait_name  = "ADG",
    description = "Average Daily Gain", 
    units       = "g/d",                  # grams per day during testing period
    type        = "continuous",           # not categorical or binary
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
  get_table("genome_meta") %>%
  filter(is_9k != TRUE) %>%
  # set loci as QTL for this trait
  define_additive_effects(
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
  define_effect_intercept(
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

# check ADG phenotype count
pop %>% get_table("ind_phenotype") %>% filter(trait_name == "ADG")
pop %>% get_table("ind_phenotype") %>% filter(trait_name == "ADG") %>% collect() %>% count()

#------------------------------------------------------------#
# Trait: BF
#------------------------------------------------------------#

warning("Define BF")

# add ADG as a trait
pop <- pop %>%
  define_trait(
    trait_name  = "BF",
    description = "Ultrasound Backfat", 
    units       = "mm",                  # grams per day during testing period
    type        = "continuous",           # not categorical or binary
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
  get_table("genome_meta") %>%
  filter(is_9k != TRUE) %>%
  # set loci as QTL for this trait
  define_additive_effects(
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

# look at TBV table
pop %>% get_table("ind_tbv") %>% filter(trait_name == "BF")

# add overall mean for ADG
pop %>%
  define_effect_intercept(
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
    trait_name = "BF",               # trait name
    rep = repl                  # add rep
  )

# print table
pop %>% get_table("ind_phenotype") %>% filter(trait_name == "BF")

# count new phenotypes
pop %>% get_table("ind_phenotype") %>% filter(trait_name == "BF") %>% collect() %>% count()

#------------------------------------------------------------#
# Trait: NW
#------------------------------------------------------------#

warning("Define NW")

# add ADG as a trait
pop <- pop %>%
  define_trait(
    trait_name  = "NW",
    description = "Number Weaned", 
    units       = "count",                  # grams per day during testing period
    type        = "count",           # not categorical or binary
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
  get_table("genome_meta") %>%
  filter(is_9k != TRUE) %>%
  # set loci as QTL for this trait
  define_additive_effects(
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
  define_effect_intercept(
    trait_name = "NW",              # trait (need to change to "trait_name")
    mean = 10
  )

# check phenotype table for NW
pop %>% get_table("ind_phenotype") %>% filter(trait_name == "NW")

# count NW phenotypes 
pop %>% get_table("ind_phenotype") %>% filter(trait_name == "NW") %>% collect() %>% count()

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
    MeanTBV = mean(tbv_value),
    .groups = "drop_last"
  ) %>%
  print(n=10)

#------------------------------------------------------------#
# Phenotypes
#------------------------------------------------------------#

# print mean phenotype by trait
pop %>%
  get_table("ind_phenotype") %>%
  collect() %>%
  group_by(rep, trait_name) %>%
  summarise(
    MeanP = mean(pheno_value),
    .groups = "drop_last"
  ) %>%
  print(n=10)



pop %>% get_table("ind_meta") %>% add_tbv(index_names = "maternal")

pop %>% get_table("ind_true_index")






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

#------------------------------------------------------------------------------#
# Run DATE Loop
#------------------------------------------------------------------------------#

# just set it for now
repl = 1

# start loop time
time_start_loop <- proc.time()

# elapsed since start
startup_elapsed <- (time_start_loop - time_start_total)["elapsed"]
message(sprintf("✔ Startup complete | %s", format_elapsed(startup_elapsed)))

# add to timing data frame
data.timing <- add_row(data.timing,
  sim_date       = start_date,
  real_date_time = Sys.time(),
  type           = "end-founder-setup",
  elapsed_sec    = NA,
  cumulative_sec = startup_elapsed
)

warning("Begin Date Loop")

# big date "loop" but continuous
for (cur_date in seq(as.Date(start_date), as.Date(end_date))){
  
# loop start time
loop_start <- proc.time()

# convert to date, not integer
cur_date = as.Date(cur_date)

# get day of week (run evals on certain days for instance)
cur_day_of_week = weekdays(cur_date)

# print "current date" (simulated advance by 1 day)
warning("Current Date: ", cur_date)

# ------------------------------ PHENOTYPE ----------------------------------- #

# ---------- AP ---------- #

message("Add AP Phenotype")

warning("Removed adding AP phenotype at this point because we added this phenotype in founders or after adding offspring")

warning("Phenotypes will be added to the evaluation by filtering by 'pheno_date'")

# pop %>%
#   get_table("ind_meta") %>%
#   filter(
#     rep == repl,
#     sex == "F", 
#     puberty_date == cur_date       # pulled and added puberty date to founders only (UPDATE)
#   ) %>%
#   add_phenotype("AP", rep = repl)

# ERROR: this is wrong, puberty date was sampled above already, this would add 
# another separate (unobserved phenotype) to the current animals. 
# 
# STEPS:
#   1. add phenotype on day of creation (founders or mating/birth)
#   2. add birth_date to sampled phenotype (age at puberty in days...) and add to `ind_meta` update `puberty_age`
#   3. add pheno_date from sum of birth_date + age at puberty to `ind_phenotype` table

# CHECKS: 
#   - do not mess up the order... 

# ---------- ADG ---------- #

message("Add ADG Phenotype")

pop %>%
  get_table("ind_meta") %>%
  filter(
    rep == repl,
    off_test_date  == cur_date
  ) %>%
  add_phenotype("ADG", rep = repl, pheno_date = cur_date)

# ---------- BF ---------- #

message("Add BF Phenotype")

pop %>%
  get_table("ind_meta") %>%
  filter(
    rep == repl,
    off_test_date  == cur_date
  ) %>%
  add_phenotype("BF", rep = repl, pheno_date = cur_date)






# ------------------------------ UPDATE STATUS: MALES ------------------------ #

message("Change Status - Males - 'after-test-boar'")

# change status based on reaching off-test
pop %>%
  get_table("ind_meta") %>%
  filter(
    rep == repl,
    sex == "M",
    off_test_date == cur_date
  ) %>%
  mutate_table(
    status = "after-test-boar"
  )

message("Change Status - Males - 'cull-juvenile-boar'")

# change status based on after-test-boar and culling age of those not selected already
pop %>%
  get_table("ind_meta") %>%
  filter(
    rep == repl,
    sex == "M",
    status == "after-test-boar",
    off_test_date < as.Date(cur_date - config$culling$boar_cull_days_after_off_test)
  ) %>%
  mutate_table(
    status = "cull-juvenile-boar"
  )

# ------------------------------ UPDATE STATUS: FEMALES ---------------------- #

message("Change Status - 'after-test-gilt'")

# change status based on reaching off-test
pop %>%
  get_table("ind_meta") %>%
  filter(
    rep == repl,
    sex == "F",
    off_test_date == cur_date,
    puberty_date > off_test_date         # skip if puberty age < off-test age
  ) %>%
  mutate_table(
    status = "after-test-gilt"
  )

message("Change Status - 'puberty-gilt'")

# change status based on reaching puberty
pop %>%
  get_table("ind_meta") %>%
  filter(
    rep == repl,
    sex == "F",
    puberty_date == cur_date
  ) %>%
  mutate_table(
    status = "puberty-gilt"
  )

message("Change Status - 'cull-gilt'")

# create new gilt cull date
cur_gilt_cull_date <- cur_date - config$culling$gilt_cull_days_after_off_test

# change status to cull if gilt is not selected
pop %>%
  get_table("ind_meta") %>%
  filter(
    rep == repl,
    sex == "F",
    status == "after-test-gilt" | status == "puberty-gilt",
    off_test_date < cur_gilt_cull_date       # past X days, cull them if not in puberty
  ) %>%
  mutate_table(
    status = "cull-gilt"       # change status to 'cull-gilt'
  )

# ------------------------------ GENOTYPE ------------------------------------ #


# ------------------------------ RUN EVALUATIONS ----------------------------- #

if (cur_date >= as.Date(config$general$start_date_evaluations) & cur_day_of_week == "Friday"){
  
  message("Calculate EBVs by Trait")

  # run EBV for AP
  pop <- pop %>%
    get_table("ind_meta") %>%
    filter(rep == repl) %>%
    add_ebv("AP", software="blupf90", model="blup",
            phenotype = pop %>% get_table("ind_phenotype") %>% 
              filter(pheno_date < cur_date | is.na(pheno_date)),
            eval_date = cur_date, rep = repl)
  
  pop %>% get_table("ind_ebv") %>%
    filter(trait_name == "AP", rep == repl)
  
  # run EBV for ADG
  pop <- pop %>%
    get_table("ind_meta") %>%
    filter(rep == repl) %>%
    add_ebv("ADG", software="blupf90", model="blup",
            phenotype = pop %>% get_table("ind_phenotype") %>% 
              filter(pheno_date < cur_date | is.na(pheno_date)),
            eval_date = cur_date, rep = repl)
  
  pop %>% get_table("ind_ebv") %>%
    filter(trait_name == "ADG", rep == repl)
  
  # run EBV for BF
  pop <- pop %>%
    get_table("ind_meta") %>%
    filter(rep == repl) %>%
    add_ebv("BF", software="blupf90", model="blup",
            phenotype = pop %>% get_table("ind_phenotype") %>% 
              filter(pheno_date < cur_date | is.na(pheno_date)),
            eval_date = cur_date, rep = repl)
  
  pop %>% get_table("ind_ebv") %>%
    filter(trait_name == "BF", rep == repl)
  
  # run EBV for NW
  pop <- pop %>%
    get_table("ind_meta") %>%
    filter(rep == repl) %>%
    add_ebv("NW", software="blupf90", model="blup",
            phenotype = pop %>% get_table("ind_phenotype") %>% 
              filter(pheno_date < cur_date | is.na(pheno_date)),
            eval_date = cur_date, rep = repl)

  pop %>% get_table("ind_ebv") %>%
    filter(trait_name == "NW", rep == repl)

} else {  # END CALCULATE EBVs
  warning("EBVs NOT RUN TODAY")
}

# ------------------------------ CALC INDEX ---------------------------------- #

if (cur_date >= as.Date(config$general$start_date_evaluations) & cur_day_of_week == "Friday"){
  
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
} else {  # END INDEX CALCULATION
  warning("No new EBVs, no need to calculate INDEXES")
}


# ------------------------------ SELECT MALES -------------------------------- #

if (cur_date == male_selection_date){
  
  message("MALE Selection Date!")
  
  message("male candidates:")
  
  # pull male candidates to select from
  male_candidates <- pop %>%
    get_table("ind_meta") %>%
    filter(
      rep == repl,
      sex == "M",
      status %in% c("after-test-boar", "breeding-boar")
    ) %>%
    pull(id_ind)
  
  print(male_candidates)
  
  # random selection early on before EBVs running with more data!
  if (cur_date < as.Date(config$general$start_date_evaluations)){
    
    message("Using RANDOM")
    
    selected_males <- pop %>%
      get_table("ind_meta") %>%
      filter(
        id_ind %in% male_candidates
      ) %>%
      slice_sample(n=config$selection$n_sires) %>%
      pull(id_ind)

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
      slice_max(index_value, n=config$selection$n_sires) %>%
      pull(id_ind)
    
    message("selected males:")
    print(selected_males)

  } # end index selection
  
  # ---------- UPDATE BOAR STATUS ---------- #
  
  # set current 'breeding-boar' to 'cull-boar' if was already a 'breeding-boar' and not in list
  list_culled_boars <- pop %>%
    get_table("ind_meta") %>%
    filter(
      status == "breeding-boar",
      !id_ind %in% selected_males
    ) %>%
    mutate_table(
      status = "cull-boar",
      cull_date = cur_date
    )
  
  # set new selected males to 'breeding-boar'
  pop %>% 
    get_table("ind_meta") %>%
    filter(
      id_ind %in% selected_males
    ) %>%
    mutate_table(
      status = "breeding-boar"
    )

} # END MALE SELECTION STEP

# ------------------------------ SELECT FEMALES ------------------------------ #

if (cur_date == female_selection_date){

  message("FEMALE Selection Date!")
  
  female_candidates <- pop %>%
    get_table("ind_meta") %>%
    filter(
      rep == repl,
      sex == "F",
      status %in% c("puberty-gilt", "open-sow")
    ) %>%
    pull(id_ind)
  
  message("female candidates:")
  print(female_candidates)
  
  # random selection early on before EBVs running
  if (cur_date < as.Date(config$general$start_date_evaluations)){
    
    selected_females = pop %>%
      get_table("ind_meta") %>%
      filter(
        id_ind %in% female_candidates
      ) %>%
      slice_sample(n=config$selection$n_dams_per_breeding) %>%
      pull(id_ind)
    
    message("randomly selected females:")
    print(selected_females)
    
  # index selection after EBVs are running
  } else {
    
    message("Using INDEX")
    
    latest_index_date = pop %>%
      get_table("ind_index") %>%
      collect() %>%
      pull(index_date) %>%
      max()
    
    selected_females <- pop %>%
      get_table("ind_index") %>%
      filter(
        id_ind %in% female_candidates,
        index_date == latest_index_date
      ) %>%
      slice_max(index_value, n=config$selection$n_dams_per_breeding) %>%
      pull(id_ind)

    message("index selected females:")
    print(selected_females)

  } # end index selection

  #-------------------- SET STATUS: FEMALES -----------------------------------#
  
  message("Change selected females status = 'selected-female'")
  
  # set females to selected to pull later
  pop <- pop %>%
    get_table("ind_meta") %>%
    filter(
      id_ind %in% selected_females
    ) %>%
    mutate_table(
      status = "selected-female"
    )

} # END FEMALE SELECTION STEP





# ------------------------------ MATE ---------------------------------------- #

if (cur_date == female_selection_date){
  
  #---------------- SAMPLE NW PHENOTYPE ON DAMS ----------------#
  
  message("Sample NW phenotype on selected dams")
  
  # phenotype date will be in the future
  cur_NW_pheno_date <- as.Date(cur_date + config$general$gest_len)
  
  # phenotype selected dams first
  pop <- pop %>%
    get_table("ind_meta") %>%
    filter(
      #rep == repl, 
      status == "selected-female"
    ) %>%
    add_phenotype(
      "NW",
      rep = repl,                      # add rep number
      pheno_date = cur_NW_pheno_date   # add future phenotype date for 'NW'
    )
  
  #---------------- EXTRACT NW PHENOTYPE ----------------#
  
  message("Extract NW phenotype")
  
  # list the selected females
  list_cur_selected_females <- pop %>%
    get_table("ind_meta") %>%
    filter(
      rep == repl, 
      status %in% c("selected-female")
    ) %>%
    pull(id_ind)
  
  # extract NW phenotype to produce offspring numbers correctly
  data.nw <- pop %>%
    get_table("ind_phenotype") %>%
    filter(
      rep == repl, 
      trait_name == "NW",
      id_ind %in% list_cur_selected_females,
      pheno_date == cur_NW_pheno_date
    ) %>%
    collect() %>%
    select(id_ind, pheno_value)
  
  message("Number of litters: ", nrow(data.nw))
    
  #---------------- BUILD PROGENY MATRIX ----------------#
  
  message("Build progeny matrix")
  
  # HOW:
  # boars sampled randomly, 1 per dam (1 mating per sire/dam pair)
  # offspring repeated based on the number of "NW" sampled above
  
  # pull list of selected males
  cur_active_boars <- pop %>%
    get_table("ind_meta") %>%
    filter(
      rep == repl,
      sex == "M",
      status == "breeding-boar"
    ) %>%
    pull(id_ind)
  
  # sample males randomly so 1 sire per dam (1 mating but multiple offspring)
  list_sampled_boar_matings <- sample(cur_active_boars, size=nrow(data.nw), replace=TRUE)
  
  # use new phenotype to build mating plan
  data.new.matings <- tibble(
    # rep sires by "NW" phenotype from dam so they match the same rows (1 litter)
    id_parent_1 = rep(list_sampled_boar_matings, time=data.nw$pheno_value),
    # rep dams by "NW" phenotype to get a full litter (1 row / offspring)
    id_parent_2 = rep(c(data.nw$id_ind), time=data.nw$pheno_value),
    line_name     = "Libra",
    sex           = sample(c("M", "F"), size=sum(data.nw$pheno_value), replace=TRUE, prob=c(0.5, 0.5)),
    rep           = repl,
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
  
  #---------------- UPDATE SOW DATES ----------------#
  
  message("Add latest mating date")
  
  # convert to dates just in case
  cur_farrow_date = as.Date(cur_date + config$general$gest_len)
  cur_wean_date   = as.Date(cur_date + config$general$gest_len + config$general$lact_len)
  
  # set dates after mating
  pop %>%
    get_table("ind_meta") %>%
    filter(
      rep == repl,
      status == "selected-female"
    ) %>%
    mutate_table(
      mate_date   = cur_date,        # add mating date
      farrow_date = cur_farrow_date, # add farrowing date
      wean_date   = cur_wean_date    # add weaning date
      #mate_date   = format(cur_date, "%Y-%m-%d"),        # add mating date
      #farrow_date = format(cur_farrow_date, "%Y-%m-%d"), # add farrowing date
      #wean_date   = format(cur_wean_date, "%Y-%m-%d")    # add weaning date
    )
  
  #-------------------- SET STATUS: FEMALES -----------------------------------#
  
  message("Update mated sows to status = 'gestation'")
  
  # set females to gest
  pop %>%
    get_table("ind_meta") %>%
    filter(
      status == "selected-female"
    ) %>%
    mutate_table(
      status = "gestation"          # update to gestation since just bred
    )
  
  #-------------------- SAMPLE PUBERTY AGE/DATE ON NEW GILTS ------------------#
  
  message("Sample 'AP' phenotype on new gilts")
  
  # ----- SAMPLE AP PHENOTYPE ON FEMALES ----- #
  
  # sample phenotype for new gilts
  pop %>%
    get_table("ind_meta") %>%     # will phenotype all individuals in this table with no filter
    filter(
      #rep == repl,
      sex == "F",
      conc_date == cur_date
    ) %>%
    add_phenotype(                # add rows to 'ind_phenotye' table
      trait = "AP",               # trait name
      rep = repl,                 # set rep number
      current_date = cur_date     # set current date (ONLY to filter right below, new phenotypes!)
    )
  
  #pop %>% get_table("ind_phenotype")
  
  # ----- EXTRACT PHENOTYPE ----- #
  
  message("Extract 'AP' phenotype on new gilts")
  
  data.age.puberty <- pop %>%
    get_table("ind_phenotype") %>%
    filter(
      current_date == cur_date
    ) %>%
    select(id_ind, pheno_value) %>%
    collect()
  
  message("list new gilt ids:")
  
  # extract list for 'ind_meta' table
  list_cur_AP_ids <- as.character(data.age.puberty$id_ind)
  
  print(list_cur_AP_ids)
  
  # ----- ADD PUBERTY DATE ----- #
  
  message("Calculate 'pheno_date' on new gilts")
  
  data.birth.date <- pop %>%
    get_table("ind_meta") %>%
    filter(
      #rep == repl,
      id_ind %in% list_cur_AP_ids
    ) %>%
    select(id_ind, birth_date) %>%
    collect() %>%
    left_join(., data.age.puberty) %>%
    mutate(
      pheno_date = birth_date + pheno_value
    )
  
  # ----- VERY IMPORTANT! DON'T MESS UP ORDER! ----- #
  
  warning("CHECK: Do your ID lists match??")
  
  # check if order is the same
  if (all(list_cur_AP_ids == data.birth.date$id_ind) == FALSE){
    stop("IDs don't match up!!!")
  } else {
    message("IDs match up")
  }
  
  # ----- UPDATE 'puberty_date' in 'ind_meta' ----- #
  
  message("Add 'puberty_date' on new gilts")
  
  # add puberty date to 'ind_meta' table
  pop <- pop %>%
    get_table("ind_meta") %>%
    filter(
      #rep == repl,
      #sex == "F",
      #conc_date == cur_date
      id_ind %in% list_cur_AP_ids
    ) %>%
    mutate_table(
      puberty_date = data.birth.date$pheno_date
    )
  
  # ----- UPDATE 'pheno_date' in 'ind_phenotype' ----- #
  
  message("Update 'pheno_date' on new gilts")
  
  # add phenotype date to phenotype table
  pop <- pop %>%
    get_table("ind_phenotype") %>%
    filter(
      #rep == repl,
      trait_name == "AP",
      id_ind %in% list_cur_AP_ids,
      current_date == cur_date
    ) %>%
    mutate_table(
      pheno_date = data.birth.date$pheno_date
    )
  
  # ----- CHECK 'ind_phenotype' ----- #
  
  # print table
  #pop %>% get_table("ind_phenotype") %>% filter(trait_name == "AP")
  
  message("Count new 'AP' phenotypes")
  
  # count how many rows added to AP phenotype data
  pop %>% 
    get_table("ind_phenotype") %>% 
    filter(trait_name == "AP", current_date == cur_date, id_ind %in% list_cur_AP_ids) %>% 
    collect() %>% 
    count()
  
} # # END MATING STEP

# ------------------------------ UPDATE STATUS ------------------------------- #

message("Update farrowed sows to status = 'lactation'")

# convert sow status to open if just weaned a litter
pop %>%
  get_table("ind_meta") %>%
  filter(
    rep == repl,
    farrow_date == cur_date
  ) %>%
  mutate_table(
    status = "lactation"  # if today is the farrow date, now considered in lactation phase
  )

message("Update weaned sows to status = 'post-wean-sow'")

# convert sow status to open if just weaned a litter
pop %>%
  get_table("ind_meta") %>%
  filter(
    rep == repl,
    wean_date == cur_date
  ) %>%
  mutate_table(
    status = "post-wean-sow"  # if today is the wean date, now considered in post-wean phase before mating
  )

message("Update weaned sows to status = 'open-sow'")

# convert sow status to open if just weaned a litter
pop %>%
  get_table("ind_meta") %>%
  filter(
    rep == repl,
    wean_date == (cur_date - config$general$w2e_int)
  ) %>%
  mutate_table(
    status = "open-sow" # now open if past weaning + wean-2-estrus interval (usually 4-5 days in sows)
  )

message("Increase selection dates for males and females")

# increase male selection_date (every 14 days)
if (cur_date == male_selection_date){
  male_selection_date = male_selection_date +
                            config$selection$male_selection_interval
}
# increase female selection_date (every 28 days)
if (cur_date == female_selection_date){
  female_selection_date = female_selection_date +
                            config$selection$female_selection_interval
}

# ------------------------------ TIMING -------------------------------------- #

# calculate loop and total time elapsed
loop_elapsed   <- (proc.time() - loop_start)["elapsed"]
total_elapsed  <- (proc.time() - time_start_total)["elapsed"]

# print message
message(sprintf("  Date %s | loop: %s | total: %s", 
                cur_date, 
                format_elapsed(loop_elapsed), 
                format_elapsed(total_elapsed)))

# add to timing data frame
data.timing <- add_row(data.timing,
  sim_date       = cur_date,
  real_date_time = Sys.time(),
  type           = "date-loop",
  elapsed_sec    = loop_elapsed,
  cumulative_sec = total_elapsed
)

} # END SIMULATION LOOP FOR DATE






# close pop
#close_pop(pop)

#stop("stopping after finishing")







#------------------------------------------------------------------------------#
# RESTORE DUCKDB DATABASE
#------------------------------------------------------------------------------#

# libraries needed
library(DBI)
library(duckdb)
library(tidyverse)

# duckdb file path
cur_duck_db_file_name = "~/Claude/tidybreed-test/age_at_puberty_tidybreed.duckdb"

# restore population object for tidybreed
pop <- restore_pop(
  db_path = cur_duck_db_file_name
)

# add TBVs you forgot...
pop <- pop %>%
  get_table("ind_meta") %>%
  add_tbv(c("ADG", "BF", "AP", "NW"))

# close pop
close_pop(pop)

# how to use DBI directly
if (1 > 2) {
  
  # reconnect to DB
  pop <- dbConnect(duckdb(), 
            dbdir = cur_duck_db_file_name)
  
  # extract latest EBVs
  data.ebv.latest <- tbl(pop, "ind_ebv") %>%
    filter(eval_date == max(eval_date)) %>%
    collect()
  
  # extract TBVs
  data.tbv <- tbl(pop, "ind_tbv") %>%
    collect()
  
  # extract TBVs
  data.ind.meta <- tbl(pop, "ind_meta") %>%
    collect()

  # close DB connection
  dbDisconnect(pop, shutdown = TRUE)

}







#------------------------------------------------------------------------------#
# UPDATE TIME TIBBLE
#------------------------------------------------------------------------------#

# ---------- MALES ---------- #

male_selection_start_date <- as.Date(config$general$start_date_selection)

#male_selection_dates <- seq(male_selection_start_date, by = "56 days", length.out = ceiling(365 * 2 / 56))

male_selection_dates <- seq.Date(
  from   = male_selection_start_date,
  to     = male_selection_start_date + (365 * 2),
  by     = as.numeric(config$selection$male_selection_interval)
)

# ---------- FEMALES ---------- #

female_selection_start_date <- as.Date(config$general$start_date_selection)

#female_selection_dates <- seq(female_selection_start_date, by = "56 days", length.out = ceiling(365 * 2 / 56))

female_selection_dates <- seq.Date(
  from   = female_selection_start_date,
  to     = female_selection_start_date + (365 * 2),
  by     = as.numeric(config$selection$female_selection_interval)
)

# ---------- EVALUATIONS ---------- #

first_friday <- as.Date(config$general$start_date_evaluations) + 
  (5 - wday(as.Date(config$general$start_date_evaluations), week_start = 1)) %% 7

# Generate every 7 days
selection_dates <- seq.Date(
  from = first_friday,
  to   = first_friday + (365 * 2),
  by   = 7
)

# ---------- UPDATE TIBBLE ---------- #

data.timing <- data.timing %>%
  mutate(
    male_selection_date = sim_date %in% male_selection_dates,
    female_selection_date = sim_date %in% female_selection_dates,
    eval_date = sim_date %in% selection_dates
  ) %>%
  mutate(
    day_type = case_when(
      male_selection_date == TRUE & female_selection_date == FALSE & eval_date == FALSE ~ "male-sel",
      male_selection_date == FALSE & female_selection_date == TRUE & eval_date == FALSE ~ "female-sel",
      male_selection_date == TRUE & female_selection_date == TRUE & eval_date == FALSE ~ "male-sel-female-sel",
      male_selection_date == TRUE & female_selection_date == FALSE & eval_date == TRUE ~ "male-sel-eval-date",
      male_selection_date == FALSE & female_selection_date == TRUE & eval_date == TRUE ~ "female-sel-eval-date",
      male_selection_date == TRUE & female_selection_date == TRUE & eval_date == TRUE ~ "male-sel-female-sel-eval-date",
      male_selection_date == FALSE & female_selection_date == FALSE & eval_date == TRUE ~ "eval-date",
      .default = "regular-day"
    )
  )
#------------------------------------------------------------------------------#
# TIMING PER LOOP
#------------------------------------------------------------------------------#

# ---------- CUMULATIVE TIME ON DATE ---------- #

data.timing %>%
  mutate(
    cumulative_min  = cumulative_sec / 60,
    cumulative_hour = cumulative_min / 60
  ) %>%
ggplot(aes(x=sim_date, y=cumulative_min, group=1)) +
  geom_line(color=hg_blue) +
  theme_bw() +
  labs(
    title = "Cumulative Minutes Total (1-day)",
    subtitle = "Daily Swine Breeding Program",
    x = "Simulation Date",
    y = "Total Minutes",
    caption = "tidybreed timing"
  )

# ---------- LOOP TIME ON DATE ---------- #

data.timing %>%
  mutate(
    cumulative_min  = cumulative_sec / 60,
    cumulative_hour = cumulative_min / 60
  ) %>%
  #filter(!is.na(elapsed_sec)) %>%
ggplot(aes(x=sim_date, y=elapsed_sec, fill=day_type, color=day_type)) +
  geom_col() +
  theme_bw() +
  labs(
    title = "Elapsed Seconds Per Loop (1-day)",
    subtitle = "Daily Swine Breeding Program",
    x = "Simulation Date",
    y = "Elapsed Seconds",
    caption = "tidybreed timing"
  )

# ---------- LOOP TIME HISTOGRAM ---------- #

data.timing %>%
  mutate(
    cumulative_min  = cumulative_sec / 60,
    cumulative_hour = cumulative_min / 60
  ) %>%
  #filter(!is.na(elapsed_sec)) %>%
ggplot(aes(x=elapsed_sec, fill=day_type)) +
  geom_histogram(color="grey30") +
  theme_bw() +
  labs(
    title = "Elapsed Seconds Per Loop (1-day)",
    subtitle = "Daily Swine Breeding Program",
    x = "Simulation Date",
    y = "Elapsed Seconds",
    caption = "tidybreed timing"
  )




#------------------------------------------------------------------------------#
# TEST HOW TO REMOVE ROWS WITH DBI
#------------------------------------------------------------------------------#

# WITHOUT GLUE PACKAGE
# id_list <- c(1, 2, 3)
# id_string <- paste(id_list, collapse = ", ")
# sql <- glue("DELETE FROM ind_phenotype WHERE id_ind IN ({id_string})")
# dbExecute(pop$db_conn, sql)

# GLUE DIDN'T WORK - should work now with removal of "ID-01" with dashed (clashes with SQL)
# sql <- glue("DELETE FROM ind_phenotype WHERE id_ind IN ({glue_collapse(id_list, sep = ', ')})")
# dbExecute(pop$db_conn, sql)

if (1 > 2){
  
  # check how many records first
  pop %>% 
    get_table("ind_phenotype") %>% 
    collect() %>% 
    filter(pheno_date > cur_date) 
  
  # list rows to delete in phenotype dataset
  list_pheno_ids <- pop %>% 
    get_table("ind_phenotype") %>% 
    collect() %>% 
    filter(pheno_date > cur_date) %>% 
    pull(id_record)
  
  # put R list into a correct string for SQL
  id_string <- paste0("'", list_pheno_ids, "'", collapse = ", ")
  
  # sql code
  sql <- glue("DELETE FROM ind_phenotype WHERE id_record IN ({id_string})")
  
  # run SQL code with `dbExecute`
  dbExecute(pop$db_conn, sql)
  
  # VERIFY the rows were deleted! 
  pop %>% 
    get_table("ind_phenotype") %>% 
    collect() %>% 
    filter(pheno_date > cur_date)
  
  # YUP - No rows left that are into the future! 

}




#------------------------------------------------------------------------------#
# Summary - Tables and Plots
#------------------------------------------------------------------------------#

# add TBVs
pop <- pop %>%
  get_table("ind_meta") %>%
  add_tbv(c("ADG", "BF", "AP", "NW"))

# calculate latest EBV evaluation date
latest_eval_date <- pop %>%
  get_table("ind_ebv") %>%
  collect() %>%
  pull(eval_date) %>%
  max()

# ---------- BAR - STATUS color SEX ---------- #

# count status by sex
pop %>%
  get_table("ind_meta") %>%
  collect() %>%
  group_by(sex, status) %>%
  summarise(
    n = n()
  ) %>%
  filter(!is.na(status)) %>%
ggplot(aes(x=status, y=n, fill=sex)) +
  geom_col() +
  geom_label(aes(label=n), fill="white") +
  scale_fill_manual("Sex", values = c("magenta3", "dodgerblue3")) +
  labs(
    title = "Status Counts",
    subtitle = "Weekly Evaluations",
    x = "(Latest) Simulation Status",
    y = "Count",
    caption = "tidybreed - 'Daily' Swine Breeding Program"
  )

# ---------- COMPARE EBVS AND TBVS ---------- #

list_old_animals <- pop %>%
  get_table("ind_meta") %>%
  filter(birth_date < as.Date("2026-01-01")) %>%
  pull(id_ind)

list_young_animals <- pop %>%
  get_table("ind_meta") %>%
  filter(birth_date > as.Date("2027-06-01")) %>%
  pull(id_ind)

# calculate mean EBVs by trait - OLD ANIMALS
pop %>%
  get_table("ind_ebv") %>%
  filter(
    eval_date == latest_eval_date,
    id_ind %in% list_old_animals) %>%
  collect() %>%
  group_by(trait_name) %>%
  summarise(
    MeanEBV = mean(ebv)
  )

# calculate mean EBVs by trait - YOUNG ANIMALS
pop %>%
  get_table("ind_ebv") %>%
  filter(
    eval_date == latest_eval_date,
    id_ind %in% list_young_animals) %>%
  collect() %>%
  group_by(trait_name) %>%
  summarise(
    MeanEBV = mean(ebv)
  )

# calculate mean TBVs by trait - OLD ANIMALS
pop %>%
  get_table("ind_tbv") %>%
  filter(
    id_ind %in% list_old_animals
  ) %>%
  collect() %>%
  group_by(trait_name) %>%
  summarise(
    MeanTBV = mean(tbv)
  )

# calculate mean TBVs by trait - YOUNG ANIMALS
pop %>%
  get_table("ind_tbv") %>%
  filter(
    id_ind %in% list_young_animals
  ) %>%
  collect() %>%
  group_by(trait_name) %>%
  summarise(
    MeanTBV = mean(tbv)
  )

# ---------- TABLE - MEAN PHENOTYPE ---------- #

# calculate mean phenotypes by trait
pop %>%
  get_table("ind_phenotype") %>%
  collect() %>%
  group_by(trait_name) %>%
  summarise(
    MeanPhenotype = mean(value)
  )

# ---------- HISTOGRAM EBVs ---------- #

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

# time-series of EBVs by Trait
# pop %>%
#   get_table("ind_ebv") %>%
#   collect() %>%
# ggplot(aes(x=eval_date, y=ebv, color=trait_name)) +
#   geom_line(aes(group=id_ind), alpha=0.05) +
#   geom_point(color=hg_blue, alpha=0.1) +
#   geom_smooth(aes(group = NULL), method = "loess", se = TRUE) +
#   facet_wrap(~trait_name, scale="free") +
#   labs(
#     title = "Timeseries EBVs (latest)",
#     subtitle = "Weekly Evaluations",
#     x = "Evaluation Date (BLUPF90)",
#     y = "Estimated Breeding Value",
#     caption = "tidybreed - 'Daily' Swine Breeding Program"
#   )

pop %>%
  get_table("ind_meta") %>%
  collect() %>%
  count()

# 8381 animals/rows

pop %>%
  get_table("ind_tbv") %>%
  collect() %>%
  count()

# 33,524 rows
# 33,524 / 4 = 8381

# TBVs tibble
data.tbvs <- pop %>%
  get_table("ind_tbv") %>%
  collect() %>%
  select(id_ind, trait_name, tbv)

# EBVs dataset (latest eval)
data.ebvs.latest <- pop %>%
  get_table("ind_ebv") %>%
  filter(eval_date == latest_eval_date) %>%
  collect()

# pull individual data
data.ind.meta <- pop %>%
  get_table("ind_meta") %>%
  collect()

# join meta data to EBVs data
data.ebvs.latest <- left_join(data.ebvs.latest, data.ind.meta)

# join meta data to EBVs data
data.tbvs <- left_join(data.tbvs, data.ind.meta)

# ---------- EBVs on Birth Date ---------- #

# EBVs on Birth Date
data.ebvs.latest %>%
ggplot(aes(x=birth_date, y=ebv, color=trait_name)) +
  #geom_point(color=hg_blue, alpha=0.1) +
  geom_hex() +
  geom_smooth(aes(group = NULL), method = "loess", se = TRUE, linewidth=2) +
  facet_wrap(~trait_name, scales="free") +
  scale_color_discrete("Trait Name") +
  scale_fill_gradient(
    low  = "grey80",
    high = "grey20"
  ) +
  labs(
    title    = "EBV (latest) trend by birth date",
    subtitle = "Weekly Evaluations",
    x        = "Birth Date",
    y        = "Estimated Breeding Value (EBV)",
    caption  = "tidybreed - 'Daily' Swine Breeding Program"
  )

# ---------- TBVs on Birth Date ---------- #

# TBV on birth date by Trait
data.tbvs %>%
ggplot(aes(x=birth_date, y=tbv, color=trait_name)) +
  #geom_point(color="grey80", alpha=0.5) +
  geom_hex() +
  geom_smooth(method = "loess", se = TRUE, linewidth=2) +
  facet_wrap(~trait_name, scales="free") +
  scale_color_discrete("Trait Name") +
  scale_fill_gradient(
    low  = "grey80",
    high = "grey20"
  ) +
  labs(
    title    = "TBV trend by birth date",
    subtitle = "Weekly Evaluations",
    x        = "Birth Date",
    y        = "True Breeding Value (TBV)",
    caption  = "tidybreed - 'Daily' Swine Breeding Program"
  )


# ---------- MEAN TBVs on Birth Date ---------- #

# TBV on birth date by Trait
data.tbvs %>%
  group_by(trait_name, birth_date) %>%
  summarise(
    MeanTBV = mean(tbv)
  ) %>%
ggplot(aes(x=birth_date, y=MeanTBV, color=trait_name)) +
  #geom_line() +
  geom_point() +
  geom_smooth(method = "loess", se = TRUE) +
  facet_wrap(~trait_name, scale="free") +
  scale_color_discrete("Trait Name") +
  labs(
    title    = "TBV trend by birth date",
    subtitle = "Weekly Evaluations",
    x        = "Birth Date",
    y        = "True Breeding Value (TBV)",
    caption  = "tidybreed - 'Daily' Swine Breeding Program"
  )

# ----------  ---------- #

# count matings by date
pop %>%
  get_table("ind_meta") %>%
  collect() %>%
  mutate(mate_date = as.factor(mate_date)) %>%
  filter(!is.na(mate_date)) %>%
ggplot(aes(x=mate_date, group=id_ind)) +
  geom_bar(fill=hg_blue) +
  labs(
    title = "Timeseries Mating Counts",
    subtitle = "28 day batch",
    x = "Mating Date",
    y = "Count Mated",
    caption = "tidybreed - 'Daily' Swine Breeding Program"
  )

# count age at puberty by date
pop %>%
  get_table("ind_meta") %>%
  collect() %>%
  #mutate(puberty_date = as.factor(puberty_date)) %>%
  filter(!is.na(puberty_date)) %>%
ggplot(aes(x=puberty_date, group=id_ind)) +
  geom_bar(fill=hg_blue) +
  labs(
    title = "Timeseries Puberty Date Counts",
    subtitle = "28 day batch",
    x = "Puberty Date (daily)",
    y = "Count",
    caption = "tidybreed - 'Daily' Swine Breeding Program"
  )

# data frame - count age at puberty by year-week
df <- pop %>%
  get_table("ind_meta") %>%
  collect() %>%
  #mutate(puberty_date = as.factor(puberty_date)) %>%
  filter(!is.na(puberty_date)) %>%
  mutate(
    puberty_year = year(puberty_date),
    puberty_week = str_pad(week(puberty_date), 2, side="left",pad = "0"),
    puberty_yw = paste(puberty_year, puberty_week, sep="_")
  )

# plot - count age at puberty by year-week
df %>%
  count(puberty_yw) %>%
  mutate(puberty_yw = as.factor(puberty_yw)) %>%
ggplot(aes(x=puberty_yw, y=n)) +
  geom_col(fill=hg_blue) +
  #scale_x_date(breaks=df$puberty_date, date_labels = "%d") +
  labs(
    title = "Timeseries Puberty Date Counts",
    subtitle = "28 day batch",
    x = "Puberty Year-Week",
    y = "Count",
    caption = "tidybreed - 'Daily' Swine Breeding Program"
  ) + 
  theme(axis.text.x = element_text(angle=45))

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

