#------------------------------------------------------------------------------#
# Age at Puberty + Sexed Semen
#------------------------------------------------------------------------------#

# Description:
#    - provided myself a new "challenge" or feature, sexed semen for top 
#      10 / 25 females mated will produce ALL males, rest will produce females

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

# load libraries
library(hrbrthemes)    # ggplot2 themes to try
library(DBI)           # DBI -> main database connection/execution package
library(glue)          # glue package for commands into DBI
library(yaml)          # load/read yaml files for input options
library(tidyverse)     # tidyverse
library(tidybreed)     # tidybreed

#------------------------------------------------------------------------------#
# Options
#------------------------------------------------------------------------------#

# REMEMBER: Set your own paths here, go into a directory to setup this project

# set based on yaml in the future or command line argument
cur_scenario_name = "age_at_puberty"

# set options
options(
  tidybreed.pop_name = "swine",
  tidybreed.base_dir = "~/Claude/tidybreed/vignettes/swine/",  # default is 'getwd()'
  tidybreed.output   = "tidybreed_output",
  tidybreed.scenario = cur_scenario_name,
  tidybreed.tools    = c("blupf90", "plink"),
  tidybreed.db_name  = "sim.duckdb",
  tidybreed.replicate = 1L,                             # set with input later
  tidybreed.archive_path = "~/Claude/tidybreed/vignettes/swine/results/",
  tidybreed.db_name_archive = paste0(cur_scenario_name, "_all_reps.duckdb") # added with archive_rep() later
)

#------------------------------------------------------------------------------#
# Set other input options 
#------------------------------------------------------------------------------#

tb_theme <- function(base_size = 12) {
  theme_minimal(base_size = base_size) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.major.y = element_line(linewidth = 0.3),
      axis.title = element_text(face = "bold"),
      plot.title = element_text(face = "bold", size = rel(1.3)),
      plot.subtitle = element_text(color = "grey40"),
      strip.text = element_text(face = "bold"),
      legend.position = "bottom",
      legend.title = element_text(face = "bold"),
      plot.title.position = "plot"
    )
}

tb_colors <- c("#E67E22","#4E79A7","#59A14F","#E15759","#B07AA1","#EDC948",
               "#2C3E50","#9D9D9D")
tb_science <- c("#D55E00","#0072B2","#009E73","#CC79A7","#F0E442","#56B4E9",
                "#999999")

# scale_fill_gradientn(
#   colours = c("#2C3E50","#4E79A7","#F5F5F5","#E67E22","#A04000")
# )

#------------------------------------------------------------------------------#
# yaml inputs
#------------------------------------------------------------------------------#

# yaml file name
config_path <- "~/Claude/tidybreed/vignettes/swine/age_at_puberty_sexed_semen.yaml"

# check for yaml file name
if (is.na(config_path)) {
  stop("No yaml file provided")
} else {
  message("yaml file name provided")
}

# check if yaml file exists
if (!file.exists(config_path)) {
  stop("yaml file does not exist, change the path above to find your local path")
} else {
  message("yaml file exists")
}

# read config file
message("reading yaml file now...")
config <- yaml::read_yaml(config_path)

# summary of config file
data.config <- summary(config) %>% as.data.frame()

# summarize yaml config file
summary(config)

# create output directory
if (dir.exists(config$output$save_dir)){
  warning("Output directory already exists")
} else {
  warning("Output directory will be created for user")
  # create directory if it doesn't exist
  dir.create(purrr::chuck(config, "output", "save_dir"), recursive=TRUE, showWarnings = TRUE)
}

# why use purrr::chuck() function? 
#  - answer: purrr::chuck() will throw an error instead of 
#    'config$level_1$level_2' will silently break and return NULL and not alert users

# print numbers
message("Number of matings per cycle: ", 
        purrr::chuck(config, "selection", "n_dams_per_breeding"))
message("Number of top females to produce males only (rest female): ", 
        purrr::chuck(config, "sexed_semen", "n_females_produce_males_per_breeding"))

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

message("Start/End Date: ", start_date, " / ", end_date)

# SEPARATE BOAR AND GILT/SOW SELECTION STEPS

# set starting selection date
female_selection_date <- as.Date(config$general$start_date_selection)
male_selection_date   <- as.Date(config$general$start_date_selection)

message("Female/Male Start Selection Date: ", female_selection_date, " / ",
                                              male_selection_date)

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

# print
print(data.timing)

#------------------------------------------------------------------------------#
# Open Population Object with database
#------------------------------------------------------------------------------#

# start/open the population object and initialize the .duckdb database
pop <- open_pop(
  clean = TRUE
)

# print pop object
print(pop)

# summary of tables
schema(pop)

#------------------------------------------------------------------------------#
# Start Genome + Population Object with database
#------------------------------------------------------------------------------#

# start population by building a genome
#
# define_genome() is a ONE-SHOT call: it writes seven tables (genome_meta,
# genome_map, ind_haplotype, ind_genotype, ind_crossover, chr_inheritance,
# chr_recombination) inside a SINGLE transaction. If anything fails, everything
# rolls back and this same `pop` can be reused for a corrected call. Calling it
# a second time on a population that already has any of those tables is an
# error -- there is no partial re-definition.
#
# Arguments are validated BEFORE anything is written: n_loci / n_chr must be
# whole numbers with n_loci >= n_chr, and chr_len_Mb / cM_per_Mb must be finite
# and strictly positive (each may be a scalar or a vector of length n_chr).
pop <- pop %>%
  define_genome(
		n_loci       = config$genome$n_loci,  # number of loci (all* -> SNP/QTL/etc)
		n_chr        = config$genome$n_chr,   # number of chromosomes
		chr_len_Mb   = config$genome$chr_len, # length in Mb (1,000,000 bp) (e.g. 1.20 and not 1_200_000)
		cM_per_Mb    = 1.0                     # genetic-map rate: pos_cM = pos_bp/1e6 * cM_per_Mb.
		                                       #   Writes the DEFAULT map to the new `genome_map`
		                                       #   table. 50 Mb chr * 1.0 = 50 cM -> ~0.5 crossovers
		                                       #   per chromosome per meiosis (lambda = cM/100).
		                                       #   Sex-/line-specific maps can be layered in later as
		                                       #   extra `genome_map` rows without a schema change.

		# Optional arguments not used here:
		#   locus_names  = ...   # length n_loci; must be unique, non-NA, non-empty
		#                        #   (default "Locus_1" ... "Locus_n")
		#   chr_names    = ...   # length n_chr, same rules (default "1" ... "18")
		#   recombines_M = TRUE  # genome-wide per-PARENT-sex recombination defaults;
		#   recombines_F = TRUE  #   set one FALSE for a whole-genome achiasmatic sex.
		#                        #   Seeded into `chr_recombination`. Per-chromosome
		#                        #   exceptions (X/Y, MT) go through define_chromosome().
	)

# print genome info table (1 row per locus)
pop %>% get_table("genome_meta")
pop %>% get_table("genome_map")

#-------------------------------#
# Verify genome tables
#-------------------------------#

# verify total counts
pop |> 
  get_table("genome_meta") |>
  count()

# verify chr counts
pop |> 
  get_table("genome_meta") |>
  count(chr_name)

# verify pos_bp locations
pop |> 
  get_table("genome_meta") |>
  collect() |>
  group_by(chr_name) |>
  summarize(
    min_bp = min(pos_bp),
    max_bp = max(pos_bp)
  )

#------------------------------------------------------------------------------#
# Add founder haplotypes
#------------------------------------------------------------------------------#

# define_founder_haplotypes() fills the `founder_haplotypes` pool that
# add_founders() samples from. Six methods, in two families:
#
#   NO LD -- draw a per-locus allele frequency p, then sample alleles at each
#   locus independently:
#     "uniform"          min_allele_freq (0.01), max_allele_freq (0.99)
#     "fixed"            allele_freq (0.5)
#     "beta"             beta_shape1 (0.5), beta_shape2 (0.5)
#     "balding_nichols"  fst (0.1), mean_allele_freq (0.5)
#
#   LD -- build correlation along the GENETIC MAP (genome_map), resolved for
#   this pool's own line_name, so founder LD is built on the same map that will
#   later drive recombination for that line:
#     "mosaic"           n_templates, template_switch_rate (1.0 per cM)
#     "gaussian_copula"  ld_decay_rate (1.0; rho = exp(-lambda * d_cM))
#
# Each argument belongs to exactly one method (except exact_freq, below) --
# passing one to the wrong method is a hard error naming the method it belongs
# to. Scalars are validated strictly: no NA/Inf, no fractional counts,
# n_templates must be a whole number in [2, n_haplotypes].
#
# Calling with line_name = NULL stores ONE shared pool that add_founders() falls
# back to when no named pool exists for the line it is building. Re-using a
# line_name that already has a pool (or NULL twice) is an error.
#
# Every call also (re)writes genome_meta.founder_allele_freq with the per-locus
# frequency of the pool written MOST RECENTLY. It is informational only -- no
# other tidybreed function reads it -- so with six lines below it describes only
# line F. For per-line Falconer centering use
#   define_additive_effects(..., base = "current_pop", base_tbl = <line-filtered>)
# rather than base = "founder_haplotypes", which pools ALL lines together.

# line A
pop <- pop %>%
  define_founder_haplotypes(
    line_name        = "A",
    n_haplotypes     = config$genome$n_haplotypes,  # number of haplotypes generated
    method           = "uniform",                   # Uniform(min, max), no LD
    min_allele_freq  = 0.01,                        # min allele freq
    max_allele_freq  = 0.99                         # max allele freq
  )

# print genome info table (1 row per locus)
pop %>% get_table("founder_haplotypes")
pop %>% get_table("founder_haplotypes") |> collect() |> count(locus_name)

# per-locus frequency of the pool written most recently (informational)
pop %>% get_table("genome_meta") |> select(locus_name, founder_allele_freq)

# line B
pop <- pop %>%
  define_founder_haplotypes(
    line_name        = "B",
    n_haplotypes     = config$genome$n_haplotypes,
    method           = "fixed",     # one frequency at every locus (default = 0.5)
    allele_freq      = 0.5          # all p = 0.50, realized EXACTLY: exactly
                                    #   round(0.5 * n_haplotypes) haplotypes carry
                                    #   the 1-allele at each locus, on an
                                    #   independently drawn subset per locus (so no
                                    #   LD is induced). This is exact_freq = TRUE,
                                    #   the default for method = "fixed". Frequencies
                                    #   live on a 1/n_haplotypes grid; an off-grid
                                    #   request warns and names what was used.
  )

# line C
pop <- pop %>%
  define_founder_haplotypes(
    line_name        = "C",
    n_haplotypes     = config$genome$n_haplotypes,
    method           = "beta",       # Beta(shape1, shape2) per-locus freq, no LD
    beta_shape1      = 0.5,          # beta shape 1 parameter value
    beta_shape2      = 0.5           # beta shape 2 parameter value
                                     #   0.5/0.5 = Jeffreys prior -> U-shaped MAF
                                     #   spectrum (many rare + many common alleles)
  )

# line D
pop <- pop %>%
  define_founder_haplotypes(
    line_name        = "D",
    n_haplotypes     = config$genome$n_haplotypes,
    method           = "balding_nichols",     # Balding-Nichols Method (no LD)
    fst              = 0.1,                   # fst value (larger -> more extreme freqs)
    mean_allele_freq = 0.5                    # ancestral mean allele freq
  )

# exact_freq: shared by all four frequency-based methods above.
#   TRUE  -> each locus gets exactly round(p * n_haplotypes) 1-alleles on an
#            independently drawn random subset, so the REALIZED pool frequency
#            equals the target p (no binomial scatter, still no LD).
#   FALSE -> alleles are independent Bernoulli(p) draws, so realized frequencies
#            scatter around p with sd = sqrt(p(1-p)/n_haplotypes).
# Defaults: TRUE for "fixed" (a fixed frequency that drifts is not fixed),
# FALSE for "uniform"/"beta"/"balding_nichols" (drawing p and then sampling
# binomially is their correct generative model). Set TRUE on those when you want
# a drift-free base frequency.
#
# pop <- pop %>%
#   define_founder_haplotypes(
#     line_name    = "G",
#     n_haplotypes = config$genome$n_haplotypes,
#     method       = "beta",
#     beta_shape1  = 0.5,
#     beta_shape2  = 0.5,
#     exact_freq   = TRUE                     # realize each drawn p exactly
#   )

# line E
pop <- pop %>%
  define_founder_haplotypes(
    line_name        = "E",
    n_haplotypes     = config$genome$n_haplotypes,
    method           = "mosaic",              # Li-Stephens block copying -> LD blocks
    n_templates      = ceiling(sqrt(config$genome$n_haplotypes)),
                                              #   n_templates also controls the MAF
                                              #   SPECTRUM, not just block length:
                                              #   templates are the only source of
                                              #   variation, so ~2/(n_templates + 1) of
                                              #   loci are monomorphic REGARDLESS of
                                              #   n_haplotypes, and MAF is quantized to
                                              #   multiples of 1/n_templates. A warning
                                              #   fires above 10% monomorphic -- QTL
                                              #   placed there add nothing to sum(2pq a^2).
                                              #   Raise it, or use "gaussian_copula".
    template_switch_rate = 1.0                # template re-draws per cM. The re-draw is
                                              #   uniform over ALL templates including the
                                              #   current one (the standard Li-Stephens
                                              #   kernel -- it makes realized LD invariant
                                              #   to marker density), so OBSERVABLE changes
                                              #   occur at rate
                                              #   template_switch_rate * (n_templates-1)/n_templates.
                                              #   0 = never switch (complete LD within a chr).
  )

# line F
pop <- pop %>%
  define_founder_haplotypes(
    line_name        = "F",
    n_haplotypes     = config$genome$n_haplotypes,
    method           = "gaussian_copula",     # AR(1) latent-normal LD; fully vectorized
                                              #   and, unlike "mosaic", gives an
                                              #   UNQUANTIZED MAF spectrum
    ld_decay_rate    = 0.25                   # LD decay rate per cM: rho = exp(-0.25 * d_cM)
                                              #   -> rho ~ 0.78 at 1 cM (slow decay / long
                                              #   LD blocks). 0 = no decay (complete LD
                                              #   within a chromosome).
  )

#------------------------------------------------------------------------------#
# Add custom fields to each table
#------------------------------------------------------------------------------#

#-------------------- ind_meta --------------------#

# ind_meta stores 1 row per individual and would track the pedigree, sex, etc

# add custom fields for simulation pipeline
pop %>% 
  get_table("ind_meta") %>% 
  mutate_table(
    #rep           = NA_integer_,    # rep (replication) number as 'integer'
    status        = NA_character_,  # status [e.g. 'juvenile', 'off-test-gilt', 'gest', 'lact', etc]
    conc_date     = as.Date(NA),    # conception date
    birth_date    = as.Date(NA),    # birth date
    on_test_date  = as.Date(NA),    # on-test date (~ 70 days old)
    off_test_date = as.Date(NA),    # off-test date (~ 160 days old)
    puberty_date  = as.Date(NA),    # puberty date (females only)
    mate_date     = as.Date(NA),    # mating date
    farrow_date   = as.Date(NA),    # farrowing date
    wean_date     = as.Date(NA),    # weaning date
    cull_date     = as.Date(NA),    # culling date
    death_date    = as.Date(NA)     # death date
  )

# add 'active' & 'alive' fields, set DEFAULT value (not missing)
pop %>%
  get_table("ind_meta") %>%
  mutate_table(
    alive = TRUE,
    active = FALSE,          # this is the default value (no rows yet)
    .set_default = TRUE      # if TRUE, sets default of given value
  )

# define descriptions of user define fields for schema() and describe_table()
pop %>%
  get_table("ind_meta") %>%
  #define_schema_description("rep",          "Replicate number (1 to n)") %>%
  define_schema_description("status",       "Production status of the animals (e.g. 'gestation')") %>%
  define_schema_description("conc_date",    "Conception Date") %>%
  define_schema_description("birth_date",   "Birth date") %>%
  define_schema_description("on_test_date", "On-Test Date") %>%
  define_schema_description("off_test_date", "Off-Test Date (often slaughter weight)") %>%
  define_schema_description("puberty_date", "Puberty/estrus date (females only)") %>%
  define_schema_description("mate_date",    "Mating date") %>%
  define_schema_description("farrow_date",  "Farrow date (e.g. gave birth)") %>%
  define_schema_description("wean_date",    "Weaning date") %>%
  define_schema_description("cull_date",    "Cull date") %>%
  define_schema_description("death_date",   "Death date") %>%
  define_schema_description("alive",        "Is alive? (logical)") %>%
  define_schema_description("active",       "Is active (reproductively)? (logical)")

pop %>% describe_table("ind_meta")

# print schema
schema(pop)

#-------------------- ind_phenotype --------------------#

# ind_phenotype stores individual phenotypes

# add custom fields to ind_phenotype table
pop %>%
  get_table("ind_phenotype") %>%
  mutate_table(
    #rep          = NA_integer_,    # add rep to phenotypes
    current_date = as.Date(NA),    # add current date (to subset ind later)
    pheno_date   = as.Date(NA),    # add phenotype date
    .set_default = TRUE            # if TRUE, sets default of given value
  )

# set new schema descriptions
pop %>%
  get_table("ind_phenotype") %>%
  #define_schema_description("rep", "Replicate number (1 to n)") %>%
  define_schema_description("current_date", "Current date") %>%
  define_schema_description("pheno_date", "Phenotype date (could be in the future...)")

# print field descriptions
pop %>% describe_table("ind_phenotype")

# print table (no rows)
pop %>% get_table("ind_phenotype")

#-------------------- ind_tbv --------------------#

# ind_tbv stores the calculated TBV for each individual (only calculated once!)

# print field descriptions
pop %>% describe_table("ind_tbv")

# print table (no rows)
pop %>% get_table("ind_tbv")

#-------------------- ind_ebv --------------------#

# ind_ebv stores your individual EBVs for each trait, long format

# add 'rep' to 'ind_ebv'
pop %>%
  get_table("ind_ebv") %>%
  mutate_table(
    #rep          = NA_integer_,     # add replicate
    eval_date    = as.Date(NA),     # add evaluation date
    .set_default = TRUE             # if TRUE, sets default of given value
  )

# set new schema descriptions
pop %>%
  get_table("ind_ebv") %>%
  #define_schema_description("rep", "Replicate number (1 to n)") %>%
  define_schema_description("eval_date", "Evaluation date")

# print field descriptions
pop %>% describe_table("ind_ebv")

# print table (no rows)
pop %>% get_table("ind_ebv")

#-------------------- ind_index --------------------#

# add 'rep' to 'ind_index' table
pop %>%
  get_table("ind_index") %>%
  mutate_table(
    index_date = as.Date(NA)     # index calculation date
    #rep        = NA_integer_      # replicate
  )

# set new schema descriptions
pop %>%
  get_table("ind_index") %>%
  #define_schema_description("rep", "Replicate number (1 to n)") %>%
  define_schema_description("index_date", "Index date of calculation")

# print field descriptions
pop %>% describe_table("ind_index")

# print table (no rows)
pop %>% get_table("ind_index")

#------------------------------------------------------------------------------#
# Add Founders
#------------------------------------------------------------------------------#

#for (repl in 1:config$general$n_reps){
repl = 1

# ----- REPRODUCIBILITY (recombination refactor, v0.53.0) ----- #
# One set.seed() here pins the entire replicate: founder sampling, birth-date
# jitter, sire sampling AND the new per-gamete dqrng recombination streams.
# Because add_offspring(seed = NULL) draws its base seed from this base-R stream,
# every mating event downstream is reproducible from this single seed. (Below we
# ALSO pass an explicit per-date seed to add_offspring() so recombination is
# pinned independently of unrelated RNG churn in the loop.)
set.seed(config$general$seed %||% (1000L + repl))

warning("\n ----------    REPLICATE: ", repl, "    --------------------\n")

# set min birth dates
min_birth_date <- as.Date(config$general$start_date) - config$general$mean_puberty_age - 60

# set max birth date
max_birth_date <- as.Date(config$general$start_date_selection) + 
                    config$general$gest_len +
                    config$general$lact_len + 
                    config$general$w2e_int + 
                    config$testing$off_test_age + 100

# set days between founders
days_between_founders <- as.numeric(max_birth_date - min_birth_date)

# print message
message("min/max birth date: ", min_birth_date, " / ", max_birth_date, 
        " (", days_between_founders, " days)")

# sample birth dates for founders to simulate starting a farm with a ladder 
# of different ages
sampled_birth_dates <- min_birth_date + 
  round(runif(n = (config$population$n_founder_male + config$population$n_founder_female),
         min = 0, max = days_between_founders))

# now create a tibble with sampled birth dates
data.founder.birth.dates <- tibble(
  birth_date = sampled_birth_dates
)

# plot birth dates
data.founder.birth.dates %>%
ggplot(aes(x=birth_date)) +
  geom_histogram(fill="darkorange1", color="white") + 
  labs(
    title = "Founder Birth Dates"
  )

# add founders
pop <- pop %>%
  get_table("founder_haplotypes") %>%      # filter specific haplotypes
  filter(
    line_name == "A"
  ) %>%
  add_founders(                                            # add founders
    n_males       = config$population$n_founder_male,      # sample male founders
    n_females     = config$population$n_founder_female,    # sample female founders
    line_name     = "A",                                   # name this line
    #rep           = repl,                                  # USER DEFINED - replicate
    conc_date     = sampled_birth_dates - config$general$gest_len, 
    birth_date    = sampled_birth_dates, 
    on_test_date  = sampled_birth_dates + config$testing$on_test_age,
    off_test_date = sampled_birth_dates + config$testing$off_test_age,
    alive         = TRUE,
    active        = FALSE
  )

pop |> get_table("ind_meta")

#------------------------------------------------------------------------------#
# Add SNP Chip
#------------------------------------------------------------------------------#

warning("Add 9k Chip")

# add 9k SNP Chip
pop %>%
  get_table("genome_meta") %>%      # pass table with loci info (all loci)
    slice_sample(n=9000) %>%        # sample 9000 SNP randomly
  define_chip(chip_name = "9k")     # define SNP Chip (give name -> assign loci)

# print genome_meta table
pop %>% get_table("genome_meta")

# define_chip() added field name "is_9k" to genome_meta

#------------------------------------------------------------------------------#
# Add Traits
#------------------------------------------------------------------------------#

#------------------------------------------------------------#
# Additive Genetic Covariance
#------------------------------------------------------------#

warning("Add genetic covariance matrix")

# additive genetic (CO)VARIANCES
vars.mat.add <- matrix(c(200.00,  0.00,    0.00,  0.00,    0.00, 0.00, 0.00,
                           0.00,  0.90,    3.07,  0.00,    0.21, 0.00, 0.00,
                           0.00,  3.07,  0.0045,  0.0058,  0.01, 0.00, 0.00,
                           0.00,  0.00,  0.0058,  0.03,    0.00, 0.00, 0.00,
                           0.00,  0.21,    0.01,  0.00,    1.20, 0.00, 0.00,
                           0.00,  0.00,    0.00,  0.00,    0.00, 0.04, 0.00,
                           0.00,  0.00,    0.00,  0.00,    0.00, 0.00, 0.13), 
                      nrow = 7, byrow=TRUE, 
                      dimnames = list(c("AP", "NW", "ADG", "ADFI", "BF", "WWD", "WWM"), 
                                      c("AP", "NW", "ADG", "ADFI", "BF", "WWD", "WWM")))

if (isSymmetric(vars.mat.add)){
  message("Check: Additive genetic (co)variance matrix is symmetric")
} else {
  stop("Check: Additive genetic (co)varianc matrix NOT symmetric")
}

# add this additive genetic covariance matrix to a table with function
pop <- pop %>%
  define_effect_cov_matrix(
    effect_name = "gen_add",    # fixed term for additive genetic (co)variance matrix
    cov_matrix  = vars.mat.add  # name of matrix with row/col names
  )

# print additive variance components
pop %>% get_table("trait_var_comp")

#------------------------------------------------------------#
# Residual Covariance
#------------------------------------------------------------#

warning("Add residual covariance matrix")

# NOTE: Only 6 phenotypes, while 7 "traits" above

# residual (CO)VARIANCES
vars.mat.res <- matrix(c(400,  0.00,  0.00,   0.0000,  0.00, 0.00,
                         0.0,  8.10,  0.00,   0.0000,  0.65, 0.00,
                         0.0,  0.00,  0.0067, 0.0077,  0.05, 0.00,
                         0.0,  0.00,  0.0077, 0.0560,  0.00, 0.00,
                         0.0,  0.65,  0.05,   0.0000,  1.30, 0.00,
                         0.0,  0.00,  0.00,   0.0000,  0.00, 0.45), 
                      nrow = 6, byrow=TRUE, 
                      dimnames = list(c("AP", "NW", "ADG", "ADFI", "BF", "WW"), 
                                      c("AP", "NW", "ADG", "ADFI", "BF", "WW")))

if (isSymmetric(vars.mat.res)){
  message("Check: Residual (co)variance matrix is symmetric")
} else {
  stop("Check: Residual (co)varianc matrix NOT symmetric")
}

# add this residual covariance matrix to a table with function
pop <- pop %>%
  define_effect_cov_matrix(
    effect_name = "residual",      # fixed term for residual (co)variance matrix
    cov_matrix  = vars.mat.res     # name for matrix with row/col names
  )

# print residual variance components
pop %>% get_table("phenotype_var_comp")

# phenotype_var_comp will also store other random effect cov matrices

#------------------------------------------------------------#
# Pen Covariance
#------------------------------------------------------------#

warning("Add pen covariance matrix")

# NOTE: Only 2 phenotypes have pen effects

# residual (CO)VARIANCES
vars.mat.pen <- matrix(c(0.0005, 0.00,   0.00,
                         0.00,   0.0072, 0.00,
                         0.00,   0.00,   0.3125), 
                      nrow = 3, byrow=TRUE, 
                      dimnames = list(c("ADG", "ADFI", "BF"), 
                                      c("ADG", "ADFI", "BF")))

if (isSymmetric(vars.mat.pen)){
  message("Check: Pen (co)variance matrix is symmetric")
} else {
  stop("Check: Pen (co)varianc matrix NOT symmetric")
}

# add this pen covariance matrix to a table with function
pop <- pop %>%
  define_effect_cov_matrix(
    effect_name = "pen",      # fixed term for residual (co)variance matrix
    cov_matrix  = vars.mat.pen     # name for matrix with row/col names
  )

# print residual variance components
pop %>% get_table("phenotype_var_comp") %>% 
  filter(effect_name=="pen")


#------------------------------------------------------------------------------#
# Index
#------------------------------------------------------------------------------#

warning("Define index")

# genetic SD per trait (WW = sqrt(Var_WWD + Var_WWM) = sqrt(0.04 + 0.13))
gen_sd <- sqrt(diag(vars.mat.add))

# percent inclusion
pct <- c(
  AP   = 0.10,  # lower is better
  NW   = 0.15,  # higher is better
  ADG  = 0.20,  # higher is better
  ADFI = 0.10,  # lower is better
  BF   = 0.10,  # lower is better
  WWD  = 0.15,  # higher is better (weaning weight direct)
  WWM  = 0.20   # higher is better (weaning weight maternal)
)
barplot(pct, col="steelblue", main="Target Percent of Index")

# direction required
direction <- c(
  AP   = -1,
  NW   =  1,
  ADG  =  1,
  ADFI = -1,
  BF   = -1,
  WWD  =  1,
  WWM  =  1
)
barplot(direction, col="steelblue", main="Direction (+1 or -1)")

acc <- c(
  AP   = 0.50,
  NW   = 0.30,
  ADG  = 0.75,
  ADFI = 0.65,
  BF   = 0.75,
  WWD  = 0.50,
  WWM  = 0.45
)
barplot(acc, col="steelblue", main="Accuracy (Guess)")

# G for EBV (not TBV)
G_ebv <- diag(acc) %*% vars.mat.add %*% diag(acc)
G_ebv %>% round(4)

# first-pass raw weights
b <- pct * direction / gen_sd
barplot(b, col="steelblue", main="b raw weights")

# expected response in original trait units
response_units <- G_ebv %*% b
rownames(response_units) <- rownames(vars.mat.add)
barplot(response_units[,1], col="steelblue", main="response (phenotype units)")

# convert response to genetic SD units
response_sd_units <- response_units / gen_sd
barplot(response_sd_units[,1], col="steelblue", main="response (sd units)")

# realized standardized emphasis
realized_pct <- abs(response_sd_units) / sum(abs(response_sd_units))
barplot(realized_pct[,1], col="steelblue", main="Realized Percent")

# diagonal matrix with inverse of genetic SD on diagonal
D_inv <- diag(1 / gen_sd)

# maps raw weights to response in SD units
M <- D_inv %*% G_ebv

# target %
target <- pct * direction
barplot(target, col="steelblue", main="Target Percent (+/-)")

b_adjusted <- solve(M, target)
names(b_adjusted) <- rownames(vars.mat.add)
barplot(b_adjusted, col="steelblue", main="b weights (adjusted)")

# check if it gives the same relative weights
check <- M %*% b_adjusted
abs(check) / sum(abs(check))

# add maternal index (WW = composite weaning weight phenotype EBV)
pop %>%
  define_index(
    index_name   = "maternal",
    trait_names  = names(b_adjusted),
    index_wts    = b_adjusted,
    economic_wts = b_adjusted
  )

# print table with index values
pop %>% get_table("index_meta")

#------------------------------------------------------------------------------#
# Add Traits + Phenotypes
#------------------------------------------------------------------------------#

#------------------------------------------------------------#
# Trait + Phenotype: Age at Puberty (AP)
#------------------------------------------------------------#

warning("Define AP - Age at Puberty")

# add ADG as a trait
pop <- pop %>%
  define_trait(
    trait_name      = "AP",
    target_add_mean = 0,
    description     = "Age at Puberty",
    units           = "days",
    overwrite       = TRUE
  ) %>%
  define_phenotype(
    phenotype_name = "AP",
    type           = "count",
    mean           = config$general$mean_puberty_age,
    expressed_sex  = "F", 
    repeatable     = FALSE,
    min_value      = 20,
    overwrite      = TRUE
  )

pop %>% get_table("trait_meta") %>% collect() %>% print(width=Inf)
pop %>% get_table("phenotype_meta") %>% collect() %>% print(width=Inf)
pop %>% get_table("phenotype_components") %>% collect() %>% print(width=Inf)

# add which loci are QTL and their effects
pop %>%
  get_table("genome_meta") %>%
    filter(is_9k != TRUE) %>%            # QTL will not be on 9k SNP chip
  define_additive_effects(
    trait_name      = "AP", 
    distribution    = "normal", 
    scale_to_target = TRUE, 
    base            = "current_pop"
  )

# print 'genome_effects'
pop %>% get_table("genome_effects") |>
  count(locus_name)

# calculate all TBV for AP
pop <- pop %>%
  get_table("ind_meta") %>% # here we specify the 'ind_meta' table so all animals will have their TBV calculated
    #filter(
    #  rep == repl
    #) %>%
  add_tbv(
    trait_name = "AP"
  )

# print TBV table
pop %>% get_table("ind_tbv")

# add 9k genotypes to all animals
pop %>%
  get_table("ind_meta") %>%
  add_genotypes(chip_name = "9k")

# print ind_meta table with new field for genotyped or not
pop %>% get_table("ind_meta") %>% slice_sample(n=2) %>% collect() %>% print(width=Inf)

# extract genotypes on 9k (all with 9k genotype)
pop %>%
  get_table("ind_meta") %>%
  extract_genotypes(chip_name = "9k")

# extract QTL for this trait (2 animals)
pop %>%
  get_table("ind_meta") %>%
    filter(
      id_ind %in% c("A_1", "A_2")
    ) %>%
  extract_genotypes(
    effects_tbl = pop %>% get_table("genome_effects") %>% filter(trait_name=="AP")
  ) %>%
  collect()


# ----- ADD OVERALL MEAN ----- #

# UPDATE: Version 0.31.0 we add 'mean' argument inside `define_phenotype()` 
# to replace this function. 

# # add overall mean for AP
# pop %>%
#   add_effect_int(
#     trait_name = "AP",              # trait (need to change to "trait_name")
#     mean = config$general$mean_puberty_age
#   )

# ----- SAMPLE AP PHENOTYPE ON FEMALES ----- #

# sample phenotype for all founder females
pop %>%
  get_table("ind_meta") %>%     # will phenotype all individuals in this table with no filter
    filter(
      #rep == repl,
      sex == "F"
    ) %>%
  add_phenotype(              # add rows to 'ind_phenotye' table
    phenotype_name = "AP"     # phenotype name
  )

# print phenotype table
pop %>% get_table("ind_phenotype")

# ----- EXTRACT PHENOTYPE ----- #

# pull phenotype sampled
data.age.puberty <- pop %>%
  get_table("ind_phenotype") %>%
  filter(
    #rep == repl,
    phenotype_name == "AP"
  ) %>%
  select(id_ind, pheno_value) %>%
  collect()

# extract IDs
list_founder_AP_ids <- as.character(data.age.puberty$id_ind)

# ----- ADD PUBERTY DATE ----- #

# extract birth dates and join to puberty age and sum to get phenotype date
data.birth.date <- pop %>%
  get_table("ind_meta") %>%
    filter(
      #rep == repl,
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

# pull IDs of animals
puberty_ids <- pop %>%
  get_table("ind_meta") %>%
  filter(
    #rep == repl,
    sex == "F"
  ) %>%
  pull(id_ind)

# add back phenotype date as "puberty_date" in ind_meta
pop <- pop %>%
  get_table("ind_meta") %>%
  mutate_table(
    puberty_date = tibble::tibble(
      id_ind = puberty_ids,
      puberty_date = data.birth.date$pheno_date
    )
  )

# ----- UPDATE 'pheno_date' in 'ind_phenotype' ----- #

# add phenotype date to phenotype table
pheno_rec_ids <- pop %>%
  get_table("ind_phenotype") %>%
  filter(
    #rep == repl,
    phenotype_name == "AP",
    id_ind %in% list_founder_AP_ids
  ) %>%
  pull(id_phenotype)

# add phenotype date now based on birth date + age at puberty
pop <- pop %>%
  get_table("ind_phenotype") %>%
  mutate_table(
    pheno_date = tibble::tibble(
      id_phenotype = pheno_rec_ids,
      pheno_date = data.birth.date$pheno_date
    )
  )

# ----- CHECK 'ind_phenotype' ----- #

# print table
pop %>% get_table("ind_phenotype") %>% filter(phenotype_name == "AP")

# count how many rows are before the start date
pop %>% 
  get_table("ind_phenotype") %>% 
    filter(phenotype_name == "AP", pheno_date < start_date) %>% 
  collect() %>% 
  count()

# ------------------------------ UPDATE STATUS: FEMALES ---------------------- #

message("Change Status - 'after-test-gilt'")

# change status based on reaching off-test
pop %>%
  get_table("ind_meta") %>%
  filter(
    #rep == repl,
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
    #rep == repl,
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
    #rep == repl,
    sex == "F",
    status == "after-test-gilt" | status == "puberty-gilt",
    off_test_date < cur_gilt_cull_date       # paste X days, cull them if not in puberty
  ) %>%
  mutate_table(
    status = "cull-gilt"
  )

#------------------------------------------------------------#
# Trait + Phenotype: Average Daily Gain (ADG)
#------------------------------------------------------------#

warning("Define ADG")

# add ADG as a trait
pop <- pop %>%
  define_trait(
    trait_name      = "ADG",
    description     = "Average Daily Gain",   # 
    units           = "kg/d",                  # grams per day during testing period
    target_add_mean = 0,                      # mean TBV in 'base'
    overwrite       = TRUE                    # wipe this row if it exists and replace with this new data
  ) %>%
  define_phenotype(
    phenotype_name = "ADG", 
    type           = "continuous",
    mean           = 0.92,
    expressed_sex  = "both",
    min_value      = 0, 
    overwrite      = TRUE
  ) 

pop %>% get_table("trait_meta")
pop %>% get_table("phenotype_meta") %>% collect() %>% print(width=Inf)

# add which loci are QTL and their effects
pop %>%
  get_table("genome_meta") %>%
    filter(is_9k != TRUE) %>%
  define_additive_effects(
    trait_name      = "ADG",        # trait name
    distribution    = "normal",     # distribution of QTL effects
    scale_to_target = TRUE,         # scale to meet additive variance target
    base            = "current_pop" # use all animals in pop to standardized (or if filtered)
  )

# add all TBV for ADG
pop <- pop %>%
  get_table("ind_meta") %>% # here we specify the 'ind_meta' table so all animals will have their TBV calculated
    #filter(rep == repl) %>%
  add_tbv(
    trait_name = "ADG"
  )

# print TBV table
pop %>% get_table("ind_tbv")

# # add overall mean for ADG
# pop %>%
#   define_effect_int(
#     trait_name = "ADG",              # trait (need to change to "trait_name")
#     mean = 1000
#   )

# add sex effect for ADG
pop %>%
  define_effect_fixed_class(
    "ADG",
    effect_name = "sex",
    source_column = "sex",
    levels = c(M = 0.08, F = 0),
    source_table = "ind_meta",
    overwrite = TRUE
  )

# trait effects
pop |> get_table("trait_effects")            # fixed effects
pop |> get_table("trait_random_effects")     # random effects sampled
pop |> get_table("phenotype_components")     # for special cases but I forget now why we left it

# test `add_phenotype()` function
pop %>%
  get_table("ind_meta") %>%     # will phenotype all individuals in this table with no filter
    filter(
      #rep == repl,
      off_test_date < start_date
      ) %>%
  add_phenotype(                # add rows to 'ind_phenotye' table
    phenotype_name = "ADG"      # trait name
  )

# check ADG phenotype count
pop %>% get_table("ind_phenotype") %>% filter(phenotype_name == "ADG")
pop %>% get_table("ind_phenotype") %>% filter(phenotype_name == "ADG") %>% 
  collect() %>% count()

#------------------------------------------------------------#
# Trait + Phenotype: Backfat (BF)
#------------------------------------------------------------#

warning("Define BF")

# add BF as a trait + phenotype
pop <- pop %>%
  define_trait(
    trait_name      = "BF",
    description     = "Ultrasound Backfat", 
    units           = "mm", 
    target_add_mean = 0,                 # mean TBV in 'base'
    overwrite       = TRUE
  ) %>%
  define_phenotype(
    phenotype_name = "BF", 
    type           = "continuous",
    mean           = 10,
    expressed_sex  = "both",
    min_value      = 0, 
    overwrite      = TRUE
  ) 

pop %>% get_table("trait_meta") %>% collect() %>% print(width=Inf)
pop %>% get_table("phenotype_meta") %>% collect() %>% print(width=Inf)

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

# add all TBV for BF
pop <- pop %>%
  get_table("ind_meta") %>% # here we specify the 'ind_meta' table so all animals will have their TBV calculated
    #filter(rep == repl) %>%
  add_tbv("BF", rep = repl)

# look at TBV table
pop %>% get_table("ind_tbv") %>% filter(trait_name == "BF")

# add overall mean for BF
pop %>%
  define_effect_intercept(
    phenotype_name = "BF",              # trait (need to change to "trait_name")
    mean = 10
  )

# test `add_phenotype()` function
pop %>%
  get_table("ind_meta") %>%     # will phenotype all individuals in this table with no filter
    filter(
      #rep == repl,
      off_test_date < start_date
    ) %>%
  add_phenotype(                # add rows to 'ind_phenotype' table
    phenotype_name = "BF"       # phenotype name
  )

# print table
pop %>% get_table("ind_phenotype") %>% filter(phenotype_name == "BF")

# count new phenotypes
pop %>% get_table("ind_phenotype") %>% filter(phenotype_name == "BF") %>% 
  collect() %>% count()

#------------------------------------------------------------#
# Trait + Phenotype: Average Daily Feed Intake (ADFI)
#------------------------------------------------------------#

warning("Define ADFI")

# add ADFI as a trait + phenotype
pop <- pop %>%
  define_trait(
    trait_name      = "ADFI",
    description     = "Average Daily Feed Intake", 
    units           = "kg/d", 
    target_add_mean = 0,                 # mean TBV in 'base'
    overwrite       = TRUE
  ) %>%
  define_phenotype(
    phenotype_name = "ADFI", 
    type           = "continuous",
    mean           = 2.52,
    expressed_sex  = "both",
    min_value      = 0, 
    overwrite      = TRUE
  ) 

pop %>% get_table("trait_meta") %>% collect() %>% print(width=Inf)
pop %>% get_table("phenotype_meta") %>% collect() %>% print(width=Inf)

# add which loci are QTL and their effects
pop %>%
  get_table("genome_meta") %>%
    filter(is_9k != TRUE) %>%
  # set loci as QTL for this trait
  define_additive_effects(
    trait_name      = "ADFI",        # trait name
    distribution    = "normal",     # distribution of QTL effects
    scale_to_target = TRUE,         # scale to meet additive variance target
    base            = "current_pop" # use all animals in pop to standardized (or if filtered)
  )

# add all TBV for ADFI
pop <- pop %>%
  get_table("ind_meta") %>% # here we specify the 'ind_meta' table so all animals will have their TBV calculated
    #filter(rep == repl) %>%
  add_tbv("ADFI")

# look at TBV table
pop %>% get_table("ind_tbv") %>% filter(trait_name == "ADFI")

# test `add_phenotype()` function
pop %>%
  get_table("ind_meta") %>%     # will phenotype all individuals in this table with no filter
    filter(
      #rep == repl,
      off_test_date < start_date
    ) %>%
  add_phenotype(                # add rows to 'ind_phenotype' table
    phenotype_name = "ADFI"     # phenotype name
  )

# print table
pop %>% get_table("ind_phenotype") %>% filter(phenotype_name == "ADFI")

# count new phenotypes
pop %>% get_table("ind_phenotype") %>% filter(phenotype_name == "ADFI") %>% 
  collect() %>% count()

#------------------------------------------------------------#
# Phenotype: Feed Conversion Ratio
#------------------------------------------------------------#

warning("Define FCR")

# add ADFI as a trait + phenotype
pop <- pop %>%
  define_phenotype(
    phenotype_name = "FCR", 
    formula        = "ADFI / ADG", 
    type           = "derived_formula",
    #mean           = 2.52,
    expressed_sex  = "both",
    min_value      = 0, 
    overwrite      = TRUE
  ) 

pop %>% get_table("trait_meta") %>% collect() %>% print(width=Inf)
pop %>% get_table("phenotype_meta") %>% collect() %>% print(width=Inf)

# test `add_phenotype()` function
pop %>%
  get_table("ind_meta") %>%     # will phenotype all individuals in this table with no filter
    filter(
      #rep == repl,
      off_test_date < start_date
    ) %>%
  add_phenotype(                # add rows to 'ind_phenotype' table
    phenotype_name = "FCR"      # phenotype name
  )

# print table
pop %>% get_table("ind_phenotype") %>% filter(phenotype_name == "FCR")

# count new phenotypes
pop %>% get_table("ind_phenotype") %>% filter(phenotype_name == "FCR") %>% 
  collect() %>% count()

#------------------------------------------------------------------------------#
# Phenotype Weaning Weight (WW)
#------------------------------------------------------------------------------#

# Description:
#   - weaning weight is the combination of weaning weight direct and maternal

#------------------------------------------------------------#
# Trait: WWD - Weaning Weight Direct
#------------------------------------------------------------#

warning("Define WWD - Weaning Weight Direct")

pop <- pop %>%
  define_trait(
    trait_name       = "WWD",
    description      = "Weaning Weight - Direct Genetic Effect",
    units            = "kg",
    expressed_parent = "both",
    target_add_mean  = 0,
    overwrite        = TRUE
  )

pop %>%
  get_table("genome_meta") %>%
    filter(is_9k != TRUE) %>%
  define_additive_effects(
    trait_name       = "WWD",
    distribution     = "normal",
    scale_to_target  = TRUE,
    base             = "current_pop"
  )

pop <- pop %>%
  get_table("ind_meta") %>%
    #filter(rep == repl) %>%
  add_tbv("WWD")

pop %>% get_table("ind_tbv") %>% filter(trait_name == "WWD")

#------------------------------------------------------------#
# Trait: WWM - Weaning Weight Maternal
#------------------------------------------------------------#

warning("Define WWM - Weaning Weight Maternal")

pop <- pop %>%
  define_trait(
    trait_name       = "WWM",
    description      = "Weaning Weight - Maternal Genetic Effect",
    units            = "kg",
    expressed_parent = "both",
    target_add_mean  = 0,
    overwrite        = TRUE
  )

pop %>%
  get_table("genome_meta") %>%
    filter(is_9k != TRUE) %>%
  define_additive_effects(
    trait_name       = "WWM",
    distribution     = "normal",
    scale_to_target  = TRUE,
    base             = "current_pop"
  )

pop <- pop %>%
  get_table("ind_meta") %>%
    #filter(rep == repl) %>%
  add_tbv("WWM")

pop %>% get_table("ind_tbv") %>% filter(trait_name == "WWM")

pop %>% get_table("trait_meta") %>% collect() %>% print(n=Inf, width=Inf)

#------------------------------------------------------------#
# Add additive effects for both WWD and WWM
#------------------------------------------------------------#

# add which loci are QTL and their effects
pop %>%
  get_table("genome_meta") %>%
    filter(is_9k != TRUE) %>%
  # set loci as QTL for this trait
  define_additive_effects(
    trait_name      = c("WWD", "WWM"), # trait names
    distribution    = "normal",        # distribution of QTL effects
    scale_to_target = TRUE,            # scale to meet additive variance target
    base            = "current_pop"    # use all animals in pop to standardized (or if filtered)
  )

#------------------------------------------------------------#
# Phenotype: WW - Weaning Weight (composite: WWD + dam(WWM))
#------------------------------------------------------------#

warning("Define WW phenotype (formula_tbv DSL)")

pop <- pop %>%
  define_phenotype(
    phenotype_name           = "WW",
    type                     = "continuous",
    formula_tbv              = "WWD + dam(WWM)",
    mean                     = config$general$wean_weight_mean,
    expressed_sex            = "both",
    repeatable               = FALSE,
    min_value                = 0,
    missing_component_action = "skip",    # founders have no dam — will be skipped
    overwrite                = TRUE
  )

pop %>% get_table("phenotype_meta") %>% collect() %>% print(n=Inf, width=Inf)
pop %>% get_table("phenotype_components") %>% collect() %>% print(width=Inf)

#------------------------------------------------------------#
# Trait: NW
#------------------------------------------------------------#

warning("Define NW")

# add NW as a trait + phenotype
pop <- pop %>%
  define_trait(
    trait_name      = "NW",
    description     = "Number Weaned", 
    units           = "count",
    target_add_mean = 0,
    overwrite       = TRUE
  ) %>%
  define_phenotype(
    phenotype_name = "NW", 
    type           = "count",
    mean           = 10,
    expressed_sex  = "F",
    min_value      = 0, 
    repeatable     = TRUE,
    overwrite      = TRUE
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
    #filter(rep == repl) %>%
  add_tbv("NW")

pop %>% get_table("ind_tbv") %>% filter(trait_name == "NW")

# add overall mean for ADG
pop %>%
  define_effect_intercept(
    phenotype_name = "NW",
    mean = 10
  )

# check phenotype table for NW
pop %>% get_table("ind_phenotype") %>% filter(phenotype_name == "NW")

# count NW phenotypes 
pop %>% get_table("ind_phenotype") %>% filter(phenotype_name == "NW") %>% 
  collect() %>% count()

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
    MeanTBV = round(mean(tbv_value), 3),
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
  group_by(phenotype_name) %>%
  summarise(
    MeanP = mean(pheno_value),
    .groups = "drop_last"
  ) %>%
  print(n=10)


#------------------------------------------------------------#
# Calculate true index value
#------------------------------------------------------------#

# add true index value given index weights
pop %>% get_table("ind_meta") %>% 
  add_tbv(index_names = "maternal")

pop %>% get_table("ind_true_index") 
pop %>% get_table("ind_true_index") %>% collect() %>% glimpse()






#------------------------------------------------------------------------------#
# Run by Time
#------------------------------------------------------------------------------#

if (1 > 2){

  # delete rows in tables if needed!  
  dbExecute(pop$db_conn, "DELETE FROM ind_phenotype WHERE phenotype_name = 'AP'")
  dbExecute(pop$db_conn, "DELETE FROM ind_phenotype WHERE phenotype_name = 'ADG'")
  dbExecute(pop$db_conn, "DELETE FROM ind_phenotype WHERE phenotype_name = 'BF'")
  dbExecute(pop$db_conn, "DELETE FROM ind_phenotype WHERE phenotype_name = 'NW'")
  dbExecute(pop$db_conn, "DELETE FROM ind_phenotype WHERE phenotype_name = 'WW'")
  # delete rows in tables if needed!
  dbExecute(pop$db_conn, "DELETE FROM ind_ebv WHERE trait_name = 'AP'")
  dbExecute(pop$db_conn, "DELETE FROM ind_ebv WHERE trait_name = 'ADG'")
  dbExecute(pop$db_conn, "DELETE FROM ind_ebv WHERE trait_name = 'BF'")
  dbExecute(pop$db_conn, "DELETE FROM ind_ebv WHERE trait_name = 'NW'")
  dbExecute(pop$db_conn, "DELETE FROM ind_ebv WHERE trait_name = 'WW'")
  # delete index rows
  dbExecute(pop$db_conn, "DELETE FROM ind_index WHERE index_name = 'maternal'")
  # delete new offspring
  dbExecute(pop$db_conn, "DELETE FROM ind_meta WHERE birth_date > '2026-08-01'")
}

#------------------------------------------------------------------------------#
# Run DATE Loop
#------------------------------------------------------------------------------#

# just set it for now
#repl = 1

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

#----------------- START LOOP -------------------------------------------------#

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

# ---------- ADG + BF ---------- #

message("Add ADG Phenotype")

pop %>%
  get_table("ind_meta") %>%
  filter(
    #rep == repl,
    off_test_date  == cur_date
  ) %>%
  add_phenotype(
    phenotype_name = c("ADG", "BF"),
    #rep = repl, 
    pheno_date = cur_date
  )

# ---------- WW ---------- #

message("Add WW Phenotype (weaning weight on piglets weaned today)")

pop <- pop %>%
  get_table("ind_meta") %>%
  filter(
    #rep == repl,
    birth_date == (cur_date - config$general$lact_len)  # piglets being weaned today
  ) %>%
  add_phenotype("WW", pheno_date = cur_date)


# ------------------------------ UPDATE STATUS: MALES ------------------------ #

message("Change Status - Males - 'after-test-boar'")

# change status based on reaching off-test
pop %>%
  get_table("ind_meta") %>%
  filter(
    #rep == repl,
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
    #rep == repl,
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
    #rep == repl,
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
    #rep == repl,
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
    #rep == repl,
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
      #filter(rep == repl) %>%
    add_ebv("AP", software="blupf90", model="blup",
      phenotype = pop %>% 
        get_table("ind_phenotype") %>% 
        filter(
          pheno_date < cur_date |    # make sure to remove future observations
          is.na(pheno_date)          # or if phenotype date is NULL/NA
        ),
      eval_date = cur_date           # add eval date
    )
  
  pop %>% get_table("ind_ebv") %>%
    filter(trait_name == "AP")
  
  # run EBV for ADG
  pop <- pop %>%
    get_table("ind_meta") %>%
      #filter(rep == repl) %>%
    add_ebv("ADG", software="blupf90", model="blup",
      phenotype = pop %>% 
        get_table("ind_phenotype") %>% 
        filter(
          pheno_date < cur_date |    # make sure to remove future observations
          is.na(pheno_date)          # or if phenotype date is NULL/NA
        ),
      eval_date = cur_date           # add data to output table
    )
  
  pop %>% get_table("ind_ebv") %>%
    filter(trait_name == "ADG")
  
  # run EBV for BF
  pop <- pop %>%
    get_table("ind_meta") %>%
      #filter(rep == repl) %>%
    add_ebv("BF", software="blupf90", model="blup",
      phenotype = pop %>% 
        get_table("ind_phenotype") %>% 
        filter(
          pheno_date < cur_date |    # make sure to remove future observations
          is.na(pheno_date)          # or if phenotype date is NULL/NA
        ),
      eval_date = cur_date           # add data to output table
    )
  
  pop %>% get_table("ind_ebv") %>%
    filter(trait_name == "BF")
  
  # run EBV for NW
  pop <- pop %>%
    get_table("ind_meta") %>%
      #filter(rep == repl) %>%
    add_ebv("NW", software="blupf90", model="blup",
      phenotype = pop %>% 
        get_table("ind_phenotype") %>% 
        filter(
          pheno_date < cur_date |    # make sure to remove future observations
          is.na(pheno_date)          # or if phenotype date is NULL/NA
        ),
      eval_date = cur_date           # add data to output table
    )

  pop %>% get_table("ind_ebv") %>%
    filter(trait_name == "NW")
  
  
  
  
  #---------------- Weaning Weight (WW) --------------------#
  
  # create run directory
  run_dir <- tidybreed:::.create_run_dir(pop, tool = "blupf90")
  
  # extract pedigree
  ped_df <- pop %>%
    get_table("ind_meta") %>%
    collect() %>%
    select(id_ind, id_parent_1, id_parent_2)
  
  # extract dataset
  data_ind_meta <- pop %>%
    get_table("ind_meta")
  
  # extract phenotypes
  data_phenotype <- pop %>%
    get_table("ind_phenotype") %>%
    filter(
      pheno_date < cur_date,
      phenotype_name == "WW"
    ) %>%
    collect()
  
  # function to convet matrix to string for glue
  matrix_to_string <- function(matrix) {
    rows <- apply(matrix, 1, function(row) paste(row, collapse = " "))
    paste(rows, collapse = "\n")
  }
  
  # Then in your glue block:
  G_WW <- vars.mat.add[c("WWD", "WWM"), c("WWD", "WWM")]
  R_WW <- as.matrix(vars.mat.res[c("WW"), c("WW")])
  
  matrix_to_string(G_WW)
  matrix_to_string(R_WW)
  
  # write paramter file directly within R
  par_lines <- glue("
DATAFILE
data.txt
SKIP_HEADER
1
TRAITS
4
FIELDS_PASSED TO OUTPUT
2
WEIGHT(S)

RESIDUAL_VARIANCE
{matrix_to_string(R_WW)}
EFFECT
5 cov
EFFECT
2 cross alpha
RANDOM
animal
OPTIONAL
mat
FILE
pedigree.txt
SKIP_HEADER
1
FILE_POS
1 2 3 0 0
PED_DEPTH
0
(CO)VARIANCES
{matrix_to_string(G_WW)}
OPTION origID
OPTION missing 0
OPTION method BLUP
")
  
  # WRITE OUT FILES
  
  # write pedigree
  write_delim(ped_df, file.path(run_dir, "pedigree.txt"), delim = " ", 
              na = "0")
  # write phenotype file
  write_delim(data_phenotype, file.path(run_dir, "data.txt"), delim = " ")
  # write parameter file
  writeLines(as.character(par_lines), file.path(run_dir, "renum.par"))
  
  # set current working directory
  old_wd <- getwd()
  
  # change to new folder we created to run the evaluation (temp directory)
  setwd(run_dir)
  
  # RUN BLUPF90
  
  # run renumf90
  system2("renumf90", args = "renum.par", stdout = "renumf90.out", stderr = "renumf90.err")
  # run blupf90+ with no VC, just BLUP
  system2("blupf90+", args = "renf90.par", stdout = "blupf90.out", stderr = "blupf90.err")
  
  # read in table of solutions
  solutions.ww <- read.table(file.path(run_dir, "solutions.orig"),
                             skip = 1) %>%
    tibble() %>%
    select(trait_renum = V1, effect_renum = V2, 
           level_renum = V3, id_ind = V4, 
           ebv_value = V5) %>%
    filter(effect_renum == 2 | effect_renum == 3) %>%
    mutate(
      trait_name = case_when(
        effect_renum == 2 ~ "WWD",
        effect_renum == 3 ~ "WWM",
        .default = "missing"
      )
    ) %>% 
    mutate(
      model = "blupf90", acc = NA, se = NA,
      eval_number = 1, eval_date = cur_date
    ) %>%
    select(-trait_renum, -level_renum, -effect_renum)
  
  # read in solutions
  solutions.ww <- solutions.ww %>%
    mutate(
      id_ebv = seq.int(
        tidybreed:::next_int_id(pop$db_conn, "ind_ebv", "id_ebv"),
        length.out = nrow(solutions.ww)
      )
    ) %>%
    relocate(id_ebv, .before = id_ind)
    
  # Insert rows into 'ind_ebv' table via DBI::dbAppendTable()
  DBI::dbAppendTable(pop$db_conn, "ind_ebv", solutions.ww)
  
  # print ebv table
  pop %>% get_table("ind_ebv") %>% print(n=5)
  
  # reset to original wd
  setwd(old_wd)
  
  # # run EBV for WW
  # pop <- pop %>%
  #   get_table("ind_meta") %>%
  #     filter(rep == repl) %>%
  #   add_ebv("WW", software="blupf90", model="blup",
  #     phenotype = pop %>% 
  #       get_table("ind_phenotype") %>% 
  #       filter(
  #         pheno_date < cur_date |    # make sure to remove future observations
  #         is.na(pheno_date)          # or if phenotype date is NULL/NA
  #       ),
  #     eval_date = cur_date           # add data to output table
  #   )

  pop %>% get_table("ind_ebv") %>%
    filter(trait_name == "WW")

} else {  # END CALCULATE EBVs
  warning("EVALUATIONS NOT RUN TODAY")
}

# ------------------------------ CALC INDEX ---------------------------------- #

if (cur_date >= as.Date(config$general$start_date_evaluations) & cur_day_of_week == "Friday"){
  
  message("Calculate Indexes")
  
  # run index calculation
  pop %>%
    get_table("ind_ebv") %>%    # must pass 'ind_ebv' because it contains the EBVs needed
    filter(
      #rep == repl,
      eval_date == cur_date
    ) %>%
    add_index(
      "maternal",          # just give the index name and it will grab weights
      index_date = cur_date
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
      #rep == repl,
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
      #rep == repl,
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
    
    message("  n Females Selected (Index): ", length(selected_females))
    
    non_selected_females = pop %>%
      get_table("ind_meta") %>%
      filter(
        id_ind %in% female_candidates,
        !(id_ind %in% selected_females)
      ) %>%
      pull(id_ind)
    
    message("  n Non-Selected Females: ", length(non_selected_females))
    
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
    
    message("  n Females Selected (Index): ", length(selected_females))
    
    non_selected_females <- pop %>%
      get_table("ind_index") %>%
      filter(
        id_ind %in% selected_females_candidates,
        !(id_ind %in% selected_females)
      ) %>%
      pull(id_ind)
    
    message("  n Non-Selected Females: ", length(non_selected_females))
    
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
  
  message("Change non-selected females status = 'cull-sow'")
  
  # set females to selected to pull later
  pop <- pop %>%
    get_table("ind_meta") %>%
    filter(
      status %in% "open-sow",
      id_ind %in% non_selected_females
    ) %>%
    mutate_table(
      status = "cull-sow"
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
      #rep = repl,                      # add rep number
      pheno_date = cur_NW_pheno_date   # add future phenotype date for 'NW'
    )
  
  #---------------- EXTRACT NW PHENOTYPE ----------------#
  
  message("Extract NW phenotype")
  
  # # random selection early on before EBVs running
  # if (cur_date < as.Date(config$general$start_date_evaluations)){
  #   
  #   # list the selected females
  #   list_cur_selected_females <- pop %>%
  #     get_table("ind_meta") %>%
  #     filter(
  #       #rep == repl,
  #       status %in% c("selected-female")
  #     ) %>%
  #     pull(id_ind)
  #   
  # } else {
  #   
  # }
  
  # 
  if (cur_date < as.Date(config$general$start_date_evaluations)){
      
    # list all selected females
    list_cur_selected_females <- pop %>%
      get_table("ind_meta") %>%
      filter(
        #rep == repl, 
        status %in% c("selected-female")
      ) %>%
      collect() %>%
      pull(id_ind)
    
    # top 10 selected females produce FEMALES ONLY
    list_cur_selected_females_sexed_males <- pop %>%
      get_table("ind_meta") %>%
      filter(
        #rep == repl, 
        status %in% c("selected-female")
      ) %>%
      collect() %>% 
      slice_sample(n=config$sexed_semen$n_females_produce_males_per_breeding) %>%
      pull(id_ind)
    
    message(" Number of MALE ONLY litters: ", length(list_cur_selected_females_sexed_males))
    
    # select females to produce FEMALES ONLY
    list_cur_selected_females_sexed_females <- list_cur_selected_females[!list_cur_selected_females %in% list_cur_selected_females_sexed_males]
    
    message(" Number of FEMALE ONLY litters: ", length(list_cur_selected_females_sexed_females))
    
  } else {
      
    # list all selected females
    list_cur_selected_females <- pop %>%
      get_table("ind_meta") %>%
      filter(
        #rep == repl, 
        status %in% c("selected-female")
      ) %>%
      collect() %>%
      pull(id_ind)
    
    # select females to produce MALES ONLY based on index
    list_cur_selected_females_sexed_males <- pop %>%
      get_table("ind_index") %>%
      filter(
        id_ind %in% list_cur_selected_females,
        index_date == latest_index_date
      ) %>%
      slice_max(index_value, n=config$sexed_semen$n_females_produce_males_per_breeding) %>%
      pull(id_ind)
    
    message(" Number of MALE ONLY litters: ", length(list_cur_selected_females_sexed_males))
    
    # select females to produce FEMALES ONLY
    list_cur_selected_females_sexed_females <- list_cur_selected_females[!list_cur_selected_females %in% list_cur_selected_females_sexed_males]
    
    message(" Number of FEMALE ONLY litters: ", length(list_cur_selected_females_sexed_females))
    
  }

  # extract NW phenotype to produce offspring numbers correctly
  data.nw <- pop %>%
    get_table("ind_phenotype") %>%
    filter(
      #rep == repl, 
      phenotype_name == "NW",
      id_ind %in% list_cur_selected_females,
      pheno_date == cur_NW_pheno_date
    ) %>%
    collect() %>%
    select(id_ind, pheno_value)
  
  data.nw.males <- data.nw %>%
    filter(id_ind %in% list_cur_selected_females_sexed_males)
  
  data.nw.females <- data.nw %>%
    filter(id_ind %in% list_cur_selected_females_sexed_females)
  
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
      #rep == repl,
      sex == "M",
      status == "breeding-boar"
    ) %>%
    pull(id_ind)
  
  # sample males randomly so 1 sire per dam (1 mating but multiple offspring)
  list_sampled_boar_matings_sexed_males <- sample(cur_active_boars, 
          size=nrow(data.nw.males), replace=TRUE)
  
  list_sampled_boar_matings_sexed_females <- sample(cur_active_boars, 
          size=nrow(data.nw.females), replace=TRUE)
  
    # use new phenotype to build mating plan
  data.new.matings.males <- tibble(
    # rep sires by "NW" phenotype from dam so they match the same rows (1 litter)
    id_parent_1 = rep(list_sampled_boar_matings_sexed_males, time=data.nw.males$pheno_value),
    # rep dams by "NW" phenotype to get a full litter (1 row / offspring)
    id_parent_2 = rep(c(data.nw.males$id_ind), time=data.nw.males$pheno_value),
    line_name     = "A",
    sex           = "M",       # SEXED SEMEN -> males only
    #rep           = repl,
    conc_date     = cur_date,
    birth_date    = cur_date + config$general$gest_len,
    on_test_date  = cur_date + config$general$gest_len + config$testing$on_test_age,
    off_test_date = cur_date + config$general$gest_len + config$testing$off_test_age
  )
  
  # use new phenotype to build mating plan
  data.new.matings.females <- tibble(
    # rep sires by "NW" phenotype from dam so they match the same rows (1 litter)
    id_parent_1 = rep(list_sampled_boar_matings_sexed_females, time=data.nw.females$pheno_value),
    # rep dams by "NW" phenotype to get a full litter (1 row / offspring)
    id_parent_2 = rep(c(data.nw.females$id_ind), time=data.nw.females$pheno_value),
    line_name     = "A",
    sex           = "F",         # SEXED SEMEN -> females only
    #rep           = repl,
    conc_date     = cur_date,
    birth_date    = cur_date + config$general$gest_len,
    on_test_date  = cur_date + config$general$gest_len + config$testing$on_test_age,
    off_test_date = cur_date + config$general$gest_len + config$testing$off_test_age
  )
  
  # stack both new matings
  data.new.matings <- bind_rows(data.new.matings.males, data.new.matings.females)
  
  #---------------- ADD OFFSPRING ----------------#
  
  message("Add offspring")
  
  # add new offspring based on tibble mating plan (1 row per offspring)
  #
  # RECOMBINATION (v0.53.0 refactor): each offspring gamete is drawn on its own
  # deterministic dqrng sub-stream keyed on (seed, offspring index, parent role),
  # using the genetic map written by define_genome() (genome_map.pos_cM). The
  # kernel runs in compiled C++ by default; force the pure-R reference with
  # Sys.setenv(TIDYBREED_KERNEL = "r") if you ever need to cross-check.
  #
  #   seed             : explicit per-DATE base seed -> this day's matings are
  #                      byte-reproducible and independent of other loop RNG.
  #   store_crossovers : TRUE also logs every crossover to the `ind_crossover`
  #                      table (id_ind, parent_origin, chr, chr_name, pos_cM).
  #   batch_size       : bounds peak memory to ~batch_size x n_loci long rows;
  #                      output is byte-identical for any batch size / same seed.
  pop %>%
    add_offspring(
      data.new.matings,
      seed             = as.integer(cur_date),   # days-since-epoch: unique per date, < int32
      store_crossovers = FALSE,                  # set TRUE to populate ind_crossover
      batch_size       = NULL                    # NULL = one batch; e.g. 2000L to cap memory
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
      #rep == repl,
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
      phenotype_name = "AP",               # trait name
      #rep = repl,                 # set rep number
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
  gilt_ids <- pop %>%
    get_table("ind_meta") %>%
    filter(
      #rep == repl,
      #sex == "F",
      #conc_date == cur_date
      id_ind %in% list_cur_AP_ids
    ) %>%
    pull(id_ind)

  pop <- pop %>%
    get_table("ind_meta") %>%
    mutate_table(
      puberty_date = tibble::tibble(
        id_ind = gilt_ids,
        puberty_date = data.birth.date$pheno_date
      )
    )
  
  # ----- UPDATE 'pheno_date' in 'ind_phenotype' ----- #
  
  message("Update 'pheno_date' on new gilts")
  
  # add phenotype date to phenotype table
  gilt_pheno_ids <- pop %>%
    get_table("ind_phenotype") %>%
    filter(
      #rep == repl,
      phenotype_name == "AP",
      id_ind %in% list_cur_AP_ids,
      current_date == cur_date
    ) %>%
    pull(id_phenotype)

  pop <- pop %>%
    get_table("ind_phenotype") %>%
    mutate_table(
      pheno_date = tibble::tibble(
        id_phenotype = gilt_pheno_ids,
        pheno_date = data.birth.date$pheno_date
      )
    )
  
  # ----- CHECK 'ind_phenotype' ----- #
  
  # print table
  #pop %>% get_table("ind_phenotype") %>% filter(trait_name == "AP")
  
  message("Count new 'AP' phenotypes")
  
  # count how many rows added to AP phenotype data
  pop %>% 
    get_table("ind_phenotype") %>% 
    filter(
      phenotype_name == "AP", 
      current_date == cur_date, 
      id_ind %in% list_cur_AP_ids
    ) %>% 
    collect() %>% 
    count()
  
} # # END MATING STEP

# ------------------------------ UPDATE STATUS ------------------------------- #

message("Update farrowed sows to status = 'lactation'")

# convert sow status to open if just weaned a litter
pop %>%
  get_table("ind_meta") %>%
  filter(
    #rep == repl,
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
    #rep == repl,
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
    #rep == repl,
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

# how to use DBI directly
if (1 > 2) {
  
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

}

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
    y = "Count",
    x = "Elapsed Seconds / Round",
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

# ---------- BAR - SEX ---------- #

pop %>%
  get_table("ind_meta") %>%
  collect() %>%
  count(sex)

# count status by sex
pop %>%
  get_table("ind_meta") %>%
  collect() %>%
  mutate(
    birth_year    = year(birth_date),
    birth_quarter = quarter(birth_date),
    birth_yq      = paste(birth_year, birth_quarter, sep="_")
  ) %>%
  group_by(sex, birth_yq) %>%
  summarise(
    n = n()
  ) %>%
ggplot(aes(x=birth_yq, y=n, fill=sex)) +
  geom_col(position="dodge") +
  #geom_label(aes(label=n), fill="white") +
  scale_fill_manual("Sex", values = c("magenta3", "dodgerblue3")) +
  labs(
    title = "Sex Count by Birth Year-Quarter",
    subtitle = "Daily Loop / Weekly Evaluations",
    x = "Birth Year-Quarter",
    y = "Count",
    caption = "tidybreed - 'Daily' Swine Breeding Program"
  ) + 
  theme(axis.text.x = element_text(angle=70))

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
    MeanEBV = mean(ebv_value)
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
    MeanEBV = mean(ebv_value)
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
    MeanTBV = mean(tbv_value)
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
    MeanTBV = mean(tbv_value)
  )

# ---------- TABLE - MEAN PHENOTYPE ---------- #

# calculate mean phenotypes by trait
pop %>%
  get_table("ind_phenotype") %>%
  collect() %>%
  group_by(trait_name) %>%
  summarise(
    MeanPhenotype = mean(pheno_value)
  )

# ---------- HISTOGRAM EBVs ---------- #

# histogram of EBVs by Trait
pop %>%
  get_table("ind_ebv") %>%
  collect() %>%
  filter(eval_date == max(eval_date)) %>%
ggplot(aes(x=ebv_value)) +
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
  select(id_ind, trait_name, tbv_value)

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
ggplot(aes(x=birth_date, y=ebv_value, color=trait_name)) +
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
ggplot(aes(x=birth_date, y=tbv_value, color=trait_name)) +
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
    MeanTBV = mean(tbv_value)
  ) %>%
ggplot(aes(x=birth_date, y=MeanTBV, color=trait_name)) +
  #geom_line() +
  geom_point() +
  geom_smooth(method = "loess", se = TRUE) +
  facet_wrap(~trait_name, scale="free") +
  scale_color_discrete("Trait Name") +
  labs(
    title    = "TBV Mean Trend by Birth Date",
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
ggplot(aes(x=pheno_value)) +
  geom_histogram(color="white", fill=hg_green, bins=17) +
  facet_wrap(~trait_name, scale="free_x") +
  labs(
    title = "Histogram - Phenotypes",
    subtitle = "Daily Loop / Sexed Semen",
    x = "Phenotype Value",
    y = "Count",
    caption = "tidybreed - 'Daily' Swine Breeding Program"
  ) 

# histogram of indexes
pop %>% get_table("ind_index") %>%
  filter(
    index_name == "maternal"
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

