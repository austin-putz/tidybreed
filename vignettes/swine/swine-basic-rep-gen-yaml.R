#------------------------------------------------------------------------------#
# Load Packages
#------------------------------------------------------------------------------#

#install.packages("pak")
#library(pak)
#pak::pak("austin-putz/tidybreed")

library(yaml)
library(tidyverse)
library(tidybreed)

#------------------------------------------------------------------------------#
# Set Input Options
#------------------------------------------------------------------------------#

#n_reps <- 3
#n_gens <- 10

#------------------------------------------------------------------------------#
# yaml inputs
#------------------------------------------------------------------------------#

#args <- commandArgs(trailingOnly = TRUE)
#config_path <- args[1]
config_path <- "~/Claude/tidybreed-test/scenarios/age_at_puberty.yaml"

# check for yaml file name
if (is.na(config_path)) stop("No yaml file provided")
# check if yaml file exists
if (!file.exists(config_path)) stop("yaml file does not exist")

# read config file
config <- yaml::read_yaml(config_path)

# create output directory
if (dir.exists(config$output$save_dir)){
  warning("Output directory already exists")
} else {
  warning("Output directory will be created for user")
  # create directory if it doesn't exist
  dir.create(config$output$save_dir, recursive=TRUE, showWarnings = TRUE)
}

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
    rep       = NA_integer_,    # rep (replication) number as 'integer'
    status    = NA_character_,   # status [juv, gest, lact, ]
    conc_date = NA_Date_,       # conception date
    birth_date = NA_Date_,      # birth date
    on_test_date = NA_Date_,    # on-test date
    off_test_date = NA_Date_,   # off-test date
    puberty_date = NA_Date_,    # puberty date
    mate_date = NA_Date_,       # mate date
    farr_date = NA_Date_,       # farrowing date
    wean_date = NA_Date_,       # wean date
    cull_date = NA_Date_,       # cull date
    death_date = NA_Date_       # death date
  )

# add 'active' field, default to FALSE
pop %>%
  get_table("ind_meta") %>%
  mutate_table(
    alive  = TRUE,          # this is the default value only
    active = FALSE,          # default value only
    .set_default = TRUE     # if TRUE, sets default of given value
    # here so when animals are born they are "active" by default, later culled
  )

# add 'rep' + 'gen_pheno' field to 'ind_phenotype'
pop %>%
  get_table("ind_phenotype") %>%
  mutate_table(
    rep = NA_integer_,          # add rep to phenotypes
    pheno_date = NA_integer_,    # this is the default value (no rows yet)
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
    eval_date = NA_Date_,       # add gen of evaluation
    .set_default = TRUE         # if TRUE, sets default of given value
  )

pop %>% get_table("ind_ebv")

# add 'rep' to 'ind_ebv'
pop %>%
  get_table("ind_index") %>%
  mutate_table(
    rep = NA_integer_,           # add rep
    index_date = NA_Date_       # add gen of evaluation
  )

pop %>% get_table("ind_index")
pop %>% get_table("index_meta")

#------------------------------------------------------------------------------#
# Add Founders
#------------------------------------------------------------------------------#

#for (repl in 1:config$general$n_reps){
repl = 1
  
warning("\n ----------    REPLICATE: ", repl, "    --------------------\n")

min_birth_date <- as.Date(config$general$start_date) - config$testing$off_test_age
max_birth_date <- as.Date(config$general$start_date) + 7
samp_birth_dates <- as.Date(min_birth_date) + 
  round(runif(n = config$population$n_founder_male + config$population$n_founder_female, 
        min = 0, max = config$testing$off_test_age + 7)

# add founders
pop <- pop %>%
  add_founders(          # add founders
    n_males   = config$population$n_founder_male,      # sample male founders
    n_females = config$population$n_founder_female,     # sample female founders
    line_name = "Libra",         # name this line
    rep       = repl,            # USER DEFINED - replicate
    
    gen_born  = 0L,              # USER DEFINED - generation
    active    = TRUE
  )

#------------------------------------------------------------------------------#
# Add SNP Chip
#------------------------------------------------------------------------------#

warning("Add 9k Chip")

# add 9k SNP Chip
pop %>%
  define_chip(
    chip_name = "9k",
    n = 900,
    method = "random"
  )

#------------------------------------------------------------------------------#
# Add Traits
#------------------------------------------------------------------------------#

warning("Add cov matrices")

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

warning("Define index")

# add generic terminal index
pop %>%
  define_index(
    index_name  = "maternal",
    trait_names = c("NW", "ADG", "BF"),
    index_wts   = c(93, 1.5, -30)
  )

# print table with index values
pop %>% get_table("index_meta")

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

warning("Run founder EBVs")

# run EBV for ADG
pop <- pop %>%
  get_table("ind_meta") %>%
  filter(rep == repl) %>%
  add_ebv("ADG", software="blupf90", model="blup", gen_eval=0L, rep = repl)

pop %>% get_table("ind_ebv") %>%
  filter(trait_name == "ADG", rep == repl)

# run EBV for BF
pop <- pop %>%
  get_table("ind_meta") %>%
  filter(rep == repl) %>%
  add_ebv("BF", software="blupf90", model="blup", gen_eval=0L, rep = repl)

pop %>% get_table("ind_ebv") %>%
  filter(trait_name == "BF", rep == repl)

# run EBV for NW
pop <- pop %>%
  get_table("ind_meta") %>%
  filter(rep == repl) %>%
  add_ebv("NW", software="blupf90", model="blup", gen_eval=0L, rep = repl)

pop %>% get_table("ind_ebv") %>%
  filter(trait_name == "NW", rep == repl)

#------------------------------------------------------------------------------#
# Run Index
#------------------------------------------------------------------------------#

warning("Calculate founder Indexes")

# run index calculation
pop %>%
  get_table("ind_ebv") %>%    # must pass 'ind_ebv' because it contains the EBVs needed
  filter(rep == repl, gen_eval == 0L) %>%
  add_index(
    "maternal",          # just give the index name and it will grab weights
    gen_eval = 0L,
    rep = repl
  )

pop %>% get_table("ind_index") %>%
  filter(index_name == "maternal", rep == repl)

pop %>% get_table("ind_index") %>%
  filter(index_name == "maternal", rep == repl) %>%
  collect() %>%
  count()


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

warning("Begin Generation Loop")

for (genl in 1:config$general$n_gens){
    
  warning("\n ----------    GEN: ", genl, "    --------------------\n")

  #---------------- SUBSET CANDIDATES ----------------#
    
  warning("Select Candidates")

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
  
  #---------------- SELECT PARENTS ----------------#
    
  warning("Select Parents")

  #---------- SELECT ON INDEX
  
  # select males
  males_selected <- pop %>%
    get_table("ind_index") %>%
    filter(
      index_name=="maternal",
      id_ind %in% male_candidates,
      gen_eval == genl - 1L,
      rep == repl
    ) %>%
    slice_max(index_value, n=config$selection$n_sires_per_gen) %>%
    pull(id_ind)

  # select females
  females_selected <- pop %>%
    get_table("ind_index") %>%
    filter(
      index_name=="maternal",
      id_ind %in% female_candidates,
      gen_eval == genl - 1L,
      rep == repl
    ) %>%
    slice_max(index_value, n=config$selection$n_dams_per_gen) %>%
    pull(id_ind)
  
  #---------- SELECT ON ADG
  
  # # select males
  # males_selected <- pop %>%
  #   get_table("ind_ebv") %>%
  #   filter(
  #     trait_name=="ADG", 
  #     gen_eval == genl - 1L,
  #     id_ind %in% male_candidates,
  #     rep == repl
  #   ) %>%
  #   slice_max(ebv, n=3) %>%
  #   pull(id_ind)
  # 
  # # select females
  # females_selected <- pop %>%
  #   get_table("ind_ebv") %>%
  #   filter(
  #     trait_name=="ADG", 
  #     gen_eval == genl - 1L,
  #     id_ind %in% female_candidates,
  #     rep == repl
  #   ) %>%
  #   slice_max(ebv, n=15) %>%
  #   pull(id_ind)
  
  #---------------- RESET ACTIVES ----------------#
    
  warning("Reset Actives")

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
    
  warning("Sample NW phenotype on selected dams")

  # phenotype selected dams first
  pop %>%
    get_table("ind_meta") %>%
    filter(rep == repl, id_ind %in% females_selected) %>%
    add_phenotype(
      "NW",
      rep = repl,          # add rep number
      gen_pheno = genl     # add generation phenotyped as current gen
    )
    
  warning("Extract NW phenotype")
  
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
    
  warning("Build progeny matrix")
  
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
    
  warning("add offspring")
  
  # add new offspring based on tibble mating plan (1 row per offspring)
  pop %>%
    add_offspring(
      data.new.matings
    )
  
  #---------------- ADD GENOTYPES ----------------#
    
  warning("add genotypes")
  
  # all Males get a genotype for 9k
  pop <- pop %>%
    get_table("ind_meta") %>%
    filter(rep == repl, sex=="M") %>%            # MALES ONLY
    add_genotypes(chip_name="9k")
  
  #---------------- ADD PHENOTYPES ----------------#
    
  warning("add phenotypes for ADG + BF")
  
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
    
  warning("run single-trait BLUPF90")
  
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

  #---------------- CALC INDEX ----------------#
    
  warning("calculate index")
  
  pop <- pop %>%
    get_table("ind_ebv") %>%
    filter(
      gen_eval == genl,
      rep == repl
    ) %>%
    add_index(
      "maternal",
      gen_eval = genl, 
      rep = repl
    )
  
  #---------------- RUN OVERLAP GENERATION ----------------#
  
  # Offspring already ACTIVE=TRUE by default now, no need to set
  
} # END GENERATION LOOP

warning("end gen loop, add all TBVs if needed")

# make sure all animals have TBVs for later
pop %>%
  get_table("ind_meta") %>%
  filter(rep == repl) %>%
  add_tbv(c("ADG", "BF", "NW"), rep = repl)

} # END REPLICATE LOOP















#------------------------------------------------------------------------------#
# Summarize Pop
#------------------------------------------------------------------------------#

#as.data.frame(config)

#tibble(data = list(config)) %>% unnest_wider(data)

#fromJSON(toJSON(config, auto_unbox = TRUE), flatten=TRUE)

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

trait_colors = c("cadetblue2", "darkolivegreen2", "darkorange", "darkorchid", "dodgerblue2")

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

