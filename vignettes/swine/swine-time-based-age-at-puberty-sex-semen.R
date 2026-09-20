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
# define_genome() is a ONE-SHOT call: it writes ten tables (genome_meta,
# genome_map, ind_haplotype, ind_genotype, ind_crossover, chr_inheritance,
# chr_recombination, plus the empty genome_effects / genome_effect_members /
# genome_effect_member_origins effect tables) and two views (genome_effect_terms,
# genome_effect_loci) inside a SINGLE transaction. If anything fails, everything
# rolls back and this same `pop` can be reused for a corrected call. Calling it
# a second time on a population that already has any of those tables is an
# error -- there is no partial re-definition.
#
# The effect tables live here rather than in open_pop() because
# genome_effect_members carries a foreign key to genome_meta.locus_id, which
# has to exist first.
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
# line F. define_additive_effects() never reads it: with line_name set it
# centers on that line's own founder pool by default, and any other base is a
# filtered table passed as base_tbl (see extract_allele_freq()).

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

pop |>
  get_table("genome_meta") |>
  count(chr) |> 
  arrange(chr)

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
    base_tbl        = get_table(pop, "ind_meta") # all animals currently in pop define p
  )

# print the causal loci (one row per term x locus; a locus-level view over
# genome_effects / genome_effect_members)
pop %>% get_table("genome_effect_loci") |>
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

# extract the causal loci for this trait (2 animals)
pop %>%
  get_table("ind_meta") %>%
    filter(
      id_ind %in% c("A_1", "A_2")
    ) %>%
  extract_genotypes(
    effects_tbl = pop %>% get_table("genome_effect_loci") %>% filter(trait_name=="AP")
  ) %>%
  collect()


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
    base_tbl        = get_table(pop, "ind_meta") # all animals currently in pop define p
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
pop |> get_table("phenotype_effects")            # fixed effects
pop |> get_table("phenotype_random_effects")     # random effects sampled
pop |> get_table("phenotype_components")     # composite phenotypes only (e.g. WW = WWD + dam(WWM), defined below)

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
    base_tbl        = get_table(pop, "ind_meta") # all animals currently in pop define p
  )

# add all TBV for BF
pop <- pop %>%
  get_table("ind_meta") %>% # here we specify the 'ind_meta' table so all animals will have their TBV calculated
    #filter(rep == repl) %>%
  add_tbv("BF")

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
    base_tbl        = get_table(pop, "ind_meta") # all animals currently in pop define p
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
    target_add_mean  = 0,
    overwrite        = TRUE
  )

# QTL effects + TBVs for WWD and WWM are written together, below, from one
# multi-trait define_additive_effects() call.

#------------------------------------------------------------#
# Trait: WWM - Weaning Weight Maternal
#------------------------------------------------------------#

warning("Define WWM - Weaning Weight Maternal")

pop <- pop %>%
  define_trait(
    trait_name       = "WWM",
    description      = "Weaning Weight - Maternal Genetic Effect",
    units            = "kg",
    target_add_mean  = 0,
    overwrite        = TRUE
  )

pop %>% get_table("trait_meta") %>% collect() %>% print(n=Inf, width=Inf)

#------------------------------------------------------------#
# Add additive effects for both WWD and WWM
#------------------------------------------------------------#

# One multi-trait call draws the two traits' QTL effects jointly from
# MVN(0, G), with G read from trait_var_comp (the WWD/WWM block of
# vars.mat.add). method = "shared" (the default) puts both traits on the same
# QTL set. Calling define_additive_effects() again for the same trait and scope
# REPLACES that trait's effects, so this must be the only call for WWD/WWM --
# otherwise any TBVs written before it go stale.
pop %>%
  get_table("genome_meta") %>%
    filter(is_9k != TRUE) %>%
  define_additive_effects(
    trait_name      = c("WWD", "WWM"), # trait names
    distribution    = "normal",        # distribution of QTL effects
    scale_to_target = TRUE,            # scale to meet additive variance target
    base_tbl        = get_table(pop, "ind_meta") # all animals currently in pop define p
  )

# causal loci per trait (locus grain; one row per term x locus)
pop %>% get_table("genome_effect_loci") %>% count(trait_name)

# TBVs for both traits in one call
pop <- pop %>%
  get_table("ind_meta") %>%
  add_tbv(c("WWD", "WWM"))

pop %>% get_table("ind_tbv") %>% filter(trait_name %in% c("WWD", "WWM"))

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
    base_tbl        = get_table(pop, "ind_meta") # all animals currently in pop define p
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
  group_by(trait_name) %>%
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

# Re-running the loop on a population that already went through it? Reset the
# loop's output first with remove_rows() (filter -> delete; refuses to wipe a
# whole table unless confirm_all = TRUE). Not run by default.
if (FALSE) {

  # drop every EBV, index value and phenotype record
  pop %>% get_table("ind_ebv")       %>% remove_rows(confirm_all = TRUE)
  pop %>% get_table("ind_index")     %>% remove_rows(confirm_all = TRUE)
  pop %>% get_table("ind_phenotype") %>% remove_rows(confirm_all = TRUE)

  # drop every animal born after the founders, from EVERY ind_* table
  # (ind_meta, ind_haplotype, ind_tbv, ind_true_index, ...)
  pop %>%
    get_table("ind_meta") %>%
    filter(!is.na(id_parent_1)) %>%
    remove_rows(tables = "all")
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

# ---------- OFF-TEST: ADG + BF + ADFI + FCR ---------- #

message("Add off-test phenotypes (ADG, BF, ADFI, FCR)")

# ADG / BF / ADFI share the subset, so their residuals are drawn jointly from
# the residual R matrix; FCR is a derived_formula phenotype (ADFI / ADG) and
# add_phenotype() sorts it after the two records it reads.
pop %>%
  get_table("ind_meta") %>%
  filter(
    #rep == repl,
    off_test_date  == cur_date
  ) %>%
  add_phenotype(
    phenotype_name = c("ADG", "BF", "ADFI", "FCR"),
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
  
  # Only records observed BEFORE today enter an evaluation: pheno_date is set
  # when the record is made, and for AP / NW it lies in the future (birth date +
  # age at puberty; mating date + gestation). The same filter is passed to
  # add_ebv() below.
  observed_phenotypes <- pop %>%
    get_table("ind_phenotype") %>%
    filter(
      pheno_date < cur_date |    # make sure to remove future observations
      is.na(pheno_date)          # or if phenotype date is NULL/NA
    )
  
  # simple (single-trait, animal-model) phenotypes with at least one observed
  # record today; the rest are skipped until records exist
  traits_to_evaluate <- observed_phenotypes %>%
    collect() %>%
    distinct(phenotype_name) %>%
    pull(phenotype_name) %>%
    intersect(c("AP", "ADG", "BF", "ADFI", "NW"))
  
  message("  Traits with observed records: ", paste(traits_to_evaluate, collapse = ", "))
  
  # one BLUPF90 animal-model run per trait (pedigree traced by add_ebv())
  for (cur_trait in traits_to_evaluate) {
    
    pop <- pop %>%
      get_table("ind_meta") %>%
        #filter(rep == repl) %>%
      add_ebv(cur_trait, software="blupf90", model="blup",
        phenotype = observed_phenotypes,   # records observed before today only
        eval_date = cur_date               # custom column on ind_ebv
      )
    
  }
  
  # today's EBVs
  pop %>% get_table("ind_ebv") %>%
    filter(eval_date == cur_date) %>%
    count(trait_name)
  
  #---------------- Weaning Weight (WW) --------------------#
  
  # add_ebv() fits single-trait animal models only, so the maternal model for
  # WW (direct WWD + maternal WWM, from the WW phenotype) is written and run
  # against BLUPF90 by hand here. Solutions are stored in ind_ebv as separate
  # WWD and WWM rows, which is what the 'maternal' index reads.
  
  # extract phenotypes (weaned before today)
  data_phenotype <- pop %>%
    get_table("ind_phenotype") %>%
    filter(
      pheno_date < cur_date,
      phenotype_name == "WW"
    ) %>%
    collect()
  
  # the first litters are weaned ~ gest_len + lact_len days after the first
  # mating; until then there is nothing to evaluate
  if (nrow(data_phenotype) == 0) {
    
    warning("No WW phenotypes yet; skipping WW maternal evaluation")
    
  } else {
  
  # create run directory
  run_dir <- tidybreed:::.create_run_dir(pop, tool = "blupf90")
  
  # extract pedigree
  ped_df <- pop %>%
    get_table("ind_meta") %>%
    collect() %>%
    select(id_ind, id_parent_1, id_parent_2)
  
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
      model = "blupf90", acc = NA_real_, se = NA_real_,
      eval_date = cur_date
    ) %>%
    select(-trait_renum, -level_renum, -effect_renum)
  
  # eval_number the way add_ebv() assigns it: per trait, MAX(eval_number) + 1
  next_eval_ww <- DBI::dbGetQuery(pop$db_conn,
    "SELECT trait_name, COALESCE(MAX(eval_number), 0) + 1 AS eval_number
       FROM ind_ebv WHERE trait_name IN ('WWD', 'WWM') GROUP BY trait_name")
  
  # assign ids + eval_number
  solutions.ww <- solutions.ww %>%
    left_join(next_eval_ww, by = "trait_name") %>%
    mutate(
      eval_number = coalesce(as.integer(eval_number), 1L),
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
  
  pop %>% get_table("ind_ebv") %>%
    filter(trait_name %in% c("WWD", "WWM"), eval_date == cur_date)
  
  } # END WW maternal evaluation

} else {  # END CALCULATE EBVs
  warning("EVALUATIONS NOT RUN TODAY")
}

# ------------------------------ CALC INDEX ---------------------------------- #

if (cur_date >= as.Date(config$general$start_date_evaluations) & cur_day_of_week == "Friday"){
  
  message("Calculate Indexes")
  
  # add_index() requires every animal to have exactly one EBV for EVERY trait
  # in the index. Early on some index traits have no EBVs yet (NW records are
  # dated at farrowing, WWD/WWM need weaned litters), so only compute the index
  # once today's evaluation covered all of them.
  index_traits <- pop %>%
    get_table("index_meta") %>%
    filter(index_name == "maternal") %>%
    pull(trait_name)
  
  evaluated_today <- pop %>%
    get_table("ind_ebv") %>%
    filter(eval_date == cur_date) %>%
    collect() %>%
    distinct(trait_name) %>%
    pull(trait_name)
  
  if (all(index_traits %in% evaluated_today)) {
    
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
    
  } else {
    warning("Index NOT calculated: no EBVs yet for ",
            paste(setdiff(index_traits, evaluated_today), collapse = ", "))
  }
  
} else {  # END INDEX CALCULATION
  warning("No new EBVs, no need to calculate INDEXES")
}

# ------------------------------ SELECTION MODE ------------------------------ #

# Selection is random until the first index has been computed, then by index.
# (The index, not the evaluation start date, is the switch: see CALC INDEX.)
index_available <- (pop %>% get_table("ind_index") %>% collect() %>% nrow()) > 0


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
  
  # random selection early on before an index exists
  if (!index_available){
    
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
  pop %>%
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
  
  # random selection early on before an index exists
  if (!index_available){
    
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
        id_ind %in% female_candidates,
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
  
  # random split early on before an index exists
  if (!index_available){
      
    # list all selected females
    list_cur_selected_females <- pop %>%
      get_table("ind_meta") %>%
      filter(
        #rep == repl, 
        status %in% c("selected-female")
      ) %>%
      collect() %>%
      pull(id_ind)
    
    # randomly chosen selected females produce MALES ONLY (sexed semen)
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
    id_parent_1 = rep(list_sampled_boar_matings_sexed_males, times = data.nw.males$pheno_value),
    # rep dams by "NW" phenotype to get a full litter (1 row / offspring)
    id_parent_2 = rep(c(data.nw.males$id_ind), times = data.nw.males$pheno_value),
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
    id_parent_1 = rep(list_sampled_boar_matings_sexed_females, times = data.nw.females$pheno_value),
    # rep dams by "NW" phenotype to get a full litter (1 row / offspring)
    id_parent_2 = rep(c(data.nw.females$id_ind), times = data.nw.females$pheno_value),
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

#------------------------------------------------------------------------------#
# After the loop: fill in TBV / TGV / true index for every animal
#------------------------------------------------------------------------------#

warning("Fill in TBVs, TGVs and true index for every animal")

# TBVs are written by add_phenotype() (for the animals it phenotypes) and by
# add_tbv() -- NOT by add_offspring(). Animals born inside the loop that never
# received a phenotype (young boars, piglets still on test) have none yet. One
# call with no trait_name evaluates every trait in trait_meta; rows that already
# exist are upserted, so this is safe to run on everyone.
pop <- pop %>%
  get_table("ind_meta") %>%
  add_tbv()

# every animal x every trait
pop %>% get_table("ind_tbv") %>% count(trait_name)

# true index for everyone (overwrite_index = TRUE recomputes the founders too)
pop <- pop %>%
  get_table("ind_meta") %>%
  add_tbv(index_names = "maternal", overwrite_index = TRUE)

pop %>% get_table("ind_true_index")

# True GENETIC values: add_tgv() evaluates EVERY stored term of a trait (additive,
# dominance, epistatic, ...) and writes ind_tgv, one row per (animal x trait x
# component_name). This model is purely additive, so the only component is
# 'order1_additive' and the total (view ind_tgv_total) equals the TBV -- a
# cheap consistency check on the evaluator.
pop <- pop %>%
  get_table("ind_meta") %>%
  add_tgv()

pop %>% get_table("ind_tgv") %>% count(trait_name, component_name)

# TBV vs TGV total: identical under an additive-only model
pop %>%
  get_table("ind_tgv_total") %>%
  collect() %>%
  inner_join(
    pop %>% get_table("ind_tbv") %>% collect(),
    by = c("id_ind", "trait_name")
  ) %>%
  summarise(max_abs_diff = max(abs(tgv_total - tbv_value)))

# summary of the population (row counts per table, on-disk size)
print(pop)
schema(pop)

#------------------------------------------------------------------------------#
# Restoring the database later
#------------------------------------------------------------------------------#

# Everything lives in the .duckdb file, so a finished (or interrupted) run can be
# reopened in a fresh R session with restore_pop(). The path is stored on the
# pop object; with the options set at the top of this script it resolves to
#   <base_dir>/<output>/<scenario>/<db_name>
message("Database file: ", pop$db_path)

# Not run here -- example of reopening in a new session
if (FALSE) {

  library(tidyverse)
  library(tidybreed)

  # same options as the top of this script (tools -> run directories)
  options(tidybreed.tools = c("blupf90", "plink"))

  # restore population object (refuses a database written by an older schema)
  pop <- restore_pop(
    db_path = "~/Claude/tidybreed/vignettes/swine/tidybreed_output/age_at_puberty/sim.duckdb"
  )

  # continue as usual, e.g. fill in TBVs
  pop <- pop %>%
    get_table("ind_meta") %>%
    add_tbv()

  # close pop
  close_pop(pop)

}

#------------------------------------------------------------------------------#
# Timing: label each simulated day
#------------------------------------------------------------------------------#

# ---------- MALES ---------- #

male_selection_start_date <- as.Date(config$general$start_date_selection)

male_selection_dates <- seq.Date(
  from   = male_selection_start_date,
  to     = end_date,
  by     = as.numeric(config$selection$male_selection_interval)
)

# ---------- FEMALES ---------- #

female_selection_start_date <- as.Date(config$general$start_date_selection)

female_selection_dates <- seq.Date(
  from   = female_selection_start_date,
  to     = end_date,
  by     = as.numeric(config$selection$female_selection_interval)
)

# ---------- EVALUATIONS ---------- #

# first Friday on/after the evaluation start date; evaluations run every Friday
first_friday <- as.Date(config$general$start_date_evaluations) + 
  (5 - wday(as.Date(config$general$start_date_evaluations), week_start = 1)) %% 7

eval_dates <- seq.Date(
  from = first_friday,
  to   = end_date,
  by   = 7
)

# ---------- UPDATE TIBBLE ---------- #

data.timing <- data.timing %>%
  mutate(
    male_selection_date   = sim_date %in% male_selection_dates,
    female_selection_date = sim_date %in% female_selection_dates,
    eval_date             = sim_date %in% eval_dates
  ) %>%
  mutate(
    day_type = case_when(
      male_selection_date == TRUE  & female_selection_date == FALSE & eval_date == FALSE ~ "male-sel",
      male_selection_date == FALSE & female_selection_date == TRUE  & eval_date == FALSE ~ "female-sel",
      male_selection_date == TRUE  & female_selection_date == TRUE  & eval_date == FALSE ~ "male-sel-female-sel",
      male_selection_date == TRUE  & female_selection_date == FALSE & eval_date == TRUE  ~ "male-sel-eval-date",
      male_selection_date == FALSE & female_selection_date == TRUE  & eval_date == TRUE  ~ "female-sel-eval-date",
      male_selection_date == TRUE  & female_selection_date == TRUE  & eval_date == TRUE  ~ "male-sel-female-sel-eval-date",
      male_selection_date == FALSE & female_selection_date == FALSE & eval_date == TRUE  ~ "eval-date",
      .default = "regular-day"
    )
  )

# seconds per day type
data.timing %>%
  filter(type == "date-loop") %>%
  group_by(day_type) %>%
  summarise(
    n_days       = n(),
    mean_sec     = mean(elapsed_sec),
    max_sec      = max(elapsed_sec),
    total_min    = sum(elapsed_sec) / 60,
    .groups      = "drop"
  ) %>%
  arrange(desc(total_min))

#------------------------------------------------------------------------------#
# Timing: plots
#------------------------------------------------------------------------------#

# ---------- CUMULATIVE TIME ON DATE ---------- #

data.timing %>%
  mutate(
    cumulative_min  = cumulative_sec / 60,
    cumulative_hour = cumulative_min / 60
  ) %>%
ggplot(aes(x=sim_date, y=cumulative_min, group=1)) +
  geom_line(color=tb_colors[2]) +
  tb_theme() +
  labs(
    title = "Cumulative Minutes Total (1-day)",
    subtitle = "Daily Swine Breeding Program",
    x = "Simulation Date",
    y = "Total Minutes",
    caption = "tidybreed timing"
  )

# ---------- LOOP TIME ON DATE ---------- #

data.timing %>%
  filter(type == "date-loop") %>%
ggplot(aes(x=sim_date, y=elapsed_sec, fill=day_type, color=day_type)) +
  geom_col() +
  scale_fill_manual("Day type", values = tb_colors, aesthetics = c("fill", "colour")) +
  tb_theme() +
  labs(
    title = "Elapsed Seconds Per Loop (1-day)",
    subtitle = "Daily Swine Breeding Program",
    x = "Simulation Date",
    y = "Elapsed Seconds",
    caption = "tidybreed timing"
  )

# ---------- LOOP TIME HISTOGRAM ---------- #

data.timing %>%
  filter(type == "date-loop") %>%
ggplot(aes(x=elapsed_sec, fill=day_type)) +
  geom_histogram(color="grey30") +
  scale_fill_manual("Day type", values = tb_colors) +
  tb_theme() +
  labs(
    title = "Elapsed Seconds Per Loop (1-day)",
    subtitle = "Daily Swine Breeding Program",
    y = "Count",
    x = "Elapsed Seconds / Round",
    caption = "tidybreed timing"
  )

#------------------------------------------------------------------------------#
# Summary - Tables and Plots
#------------------------------------------------------------------------------#

warning("Summary tables and plots")

# latest EBV evaluation date
latest_eval_date <- pop %>%
  get_table("ind_ebv") %>%
  collect() %>%
  pull(eval_date) %>%
  max()

message("Latest evaluation date: ", latest_eval_date)

# ---------- ONE TIBBLE: individuals + latest EBVs + TBVs ---------- #

# individuals, with birth year-quarter as the time axis for trend plots
data.ind.meta <- pop %>%
  get_table("ind_meta") %>%
  collect() %>%
  mutate(
    birth_year    = year(birth_date),
    birth_quarter = quarter(birth_date),
    birth_yq      = paste(birth_year, birth_quarter, sep = "_")
  )

# latest EBVs (one row per animal x trait)
data.ebvs.latest <- pop %>%
  get_table("ind_ebv") %>%
  filter(eval_date == latest_eval_date) %>%
  collect() %>%
  select(id_ind, trait_name, ebv_value)

# TBVs (one row per animal x trait, every animal after add_tbv() above)
data.tbvs <- pop %>%
  get_table("ind_tbv") %>%
  collect() %>%
  select(id_ind, trait_name, tbv_value)

# join: every animal x trait, EBV is NA for animals not in the latest evaluation
data.tbv.ebv <- data.tbvs %>%
  left_join(data.ebvs.latest, by = c("id_ind", "trait_name")) %>%
  left_join(data.ind.meta,    by = "id_ind")

glimpse(data.tbv.ebv)

# ---------- COUNTS ---------- #

# animals by sex
data.ind.meta %>% count(sex)

# animals by sex and birth year-quarter
data.ind.meta %>%
  count(sex, birth_yq) %>%
ggplot(aes(x=birth_yq, y=n, fill=sex)) +
  geom_col(position="dodge") +
  scale_fill_manual("Sex", values = c(F = "magenta3", M = "dodgerblue3")) +
  tb_theme() +
  labs(
    title = "Sex Count by Birth Year-Quarter",
    subtitle = "Daily Loop / Weekly Evaluations",
    x = "Birth Year-Quarter",
    y = "Count",
    caption = "tidybreed - 'Daily' Swine Breeding Program"
  ) + 
  theme(axis.text.x = element_text(angle=70, hjust=1))

# status by sex (latest status of every animal)
data.ind.meta %>%
  filter(!is.na(status)) %>%
  count(sex, status) %>%
ggplot(aes(x=status, y=n, fill=sex)) +
  geom_col() +
  geom_label(aes(label=n), fill="white") +
  scale_fill_manual("Sex", values = c(F = "magenta3", M = "dodgerblue3")) +
  tb_theme() +
  labs(
    title = "Status Counts",
    subtitle = "Weekly Evaluations",
    x = "(Latest) Simulation Status",
    y = "Count",
    caption = "tidybreed - 'Daily' Swine Breeding Program"
  ) + 
  theme(axis.text.x = element_text(angle=45, hjust=1))

# phenotype records by phenotype and year-quarter phenotyped
pop %>%
  get_table("ind_phenotype") %>%
  collect() %>%
  mutate(
    pheno_yq = paste(year(pheno_date), quarter(pheno_date), sep = "_")
  ) %>%
  count(phenotype_name, pheno_yq) %>%
  pivot_wider(names_from = phenotype_name, values_from = n, values_fill = 0) %>%
  arrange(pheno_yq) %>%
  print(n = Inf)

# ---------- MEANS: OLD vs YOUNG ANIMALS ---------- #

# "old" = born before the simulation started (founders); "young" = born in the
# final year of the simulation
list_old_animals   <- data.ind.meta %>% filter(birth_date <  start_date)       %>% pull(id_ind)
list_young_animals <- data.ind.meta %>% filter(birth_date >= end_date - 365)   %>% pull(id_ind)

message("n old / young animals: ", length(list_old_animals), " / ", length(list_young_animals))

# mean EBV (latest) and TBV by trait, old vs young
data.tbv.ebv %>%
  mutate(
    age_group = case_when(
      id_ind %in% list_old_animals   ~ "old (founders)",
      id_ind %in% list_young_animals ~ "young (last year)",
      .default = "middle"
    )
  ) %>%
  group_by(age_group, trait_name) %>%
  summarise(
    n        = n(),
    mean_ebv = mean(ebv_value, na.rm = TRUE),
    mean_tbv = mean(tbv_value),
    .groups  = "drop"
  ) %>%
  arrange(trait_name, age_group) %>%
  print(n = Inf)

# mean phenotype by phenotype name
pop %>%
  get_table("ind_phenotype") %>%
  collect() %>%
  group_by(phenotype_name) %>%
  summarise(
    n              = n(),
    mean_phenotype = mean(pheno_value),
    sd_phenotype   = sd(pheno_value),
    .groups        = "drop"
  )

# ---------- HISTOGRAMS ---------- #

# histogram of EBVs (latest) by trait
data.ebvs.latest %>%
ggplot(aes(x=ebv_value)) +
  geom_histogram(color="white", fill=tb_colors[2]) +
  facet_wrap(~trait_name, scales="free") +
  tb_theme() +
  labs(
    title = "EBVs (latest evaluation)",
    subtitle = paste("Evaluation date:", latest_eval_date),
    x = "Estimated Breeding Value (EBV)",
    y = "Count"
  )

# histogram of phenotypes by phenotype
pop %>%
  get_table("ind_phenotype") %>%
  collect() %>%
ggplot(aes(x=pheno_value)) +
  geom_histogram(color="white", fill=tb_colors[3], bins=17) +
  facet_wrap(~phenotype_name, scales="free_x") +
  tb_theme() +
  labs(
    title = "Histogram - Phenotypes",
    subtitle = "Daily Loop / Sexed Semen",
    x = "Phenotype Value",
    y = "Count",
    caption = "tidybreed - 'Daily' Swine Breeding Program"
  ) 

# histogram of index values (latest)
pop %>% get_table("ind_index") %>%
  filter(index_name == "maternal") %>%
  collect() %>%
  filter(index_date == max(index_date)) %>%
ggplot(aes(x=index_value)) +
  geom_histogram(color="white", fill="aquamarine3") +
  tb_theme() +
  labs(
    title = "Index Values (latest)"
  )

# index rows per calculation date
pop %>% get_table("ind_index") %>%
  filter(index_name == "maternal") %>%
  collect() %>%
  count(index_name, index_date) %>%
  print(n = Inf)

# ---------- TRENDS ON BIRTH DATE ---------- #

# EBVs (latest) on birth date
data.tbv.ebv %>%
  filter(!is.na(ebv_value)) %>%
ggplot(aes(x=birth_date, y=ebv_value)) +
  geom_hex() +
  geom_smooth(method = "loess", se = TRUE, linewidth=2, color=tb_colors[1]) +
  facet_wrap(~trait_name, scales="free") +
  scale_fill_gradient(low = "grey80", high = "grey20") +
  tb_theme() +
  labs(
    title    = "EBV (latest) trend by birth date",
    subtitle = "Weekly Evaluations",
    x        = "Birth Date",
    y        = "Estimated Breeding Value (EBV)",
    caption  = "tidybreed - 'Daily' Swine Breeding Program"
  )

# TBVs on birth date
data.tbv.ebv %>%
ggplot(aes(x=birth_date, y=tbv_value)) +
  geom_hex() +
  geom_smooth(method = "loess", se = TRUE, linewidth=2, color=tb_colors[1]) +
  facet_wrap(~trait_name, scales="free") +
  scale_fill_gradient(low = "grey80", high = "grey20") +
  tb_theme() +
  labs(
    title    = "TBV trend by birth date",
    subtitle = "Weekly Evaluations",
    x        = "Birth Date",
    y        = "True Breeding Value (TBV)",
    caption  = "tidybreed - 'Daily' Swine Breeding Program"
  )

# mean TBV by birth date
data.tbv.ebv %>%
  group_by(trait_name, birth_date) %>%
  summarise(MeanTBV = mean(tbv_value), .groups = "drop") %>%
ggplot(aes(x=birth_date, y=MeanTBV, color=trait_name)) +
  geom_point(alpha=0.5) +
  geom_smooth(method = "loess", se = TRUE) +
  facet_wrap(~trait_name, scales="free") +
  scale_color_manual("Trait Name", values = tb_colors) +
  tb_theme() +
  labs(
    title    = "TBV Mean Trend by Birth Date",
    subtitle = "Weekly Evaluations",
    x        = "Birth Date",
    y        = "True Breeding Value (TBV)",
    caption  = "tidybreed - 'Daily' Swine Breeding Program"
  )

# ---------- MATINGS AND PUBERTY OVER TIME ---------- #

# count matings by date
data.ind.meta %>%
  filter(!is.na(mate_date)) %>%
  count(mate_date) %>%
ggplot(aes(x=mate_date, y=n)) +
  geom_col(fill=tb_colors[2]) +
  tb_theme() +
  labs(
    title = "Timeseries Mating Counts",
    subtitle = paste(config$selection$female_selection_interval, "day batch"),
    x = "Mating Date",
    y = "Count Mated",
    caption = "tidybreed - 'Daily' Swine Breeding Program"
  )

# count puberty dates by year-week
data.ind.meta %>%
  filter(!is.na(puberty_date)) %>%
  mutate(
    puberty_year = year(puberty_date),
    puberty_week = str_pad(week(puberty_date), 2, side="left", pad = "0"),
    puberty_yw   = paste(puberty_year, puberty_week, sep="_")
  ) %>%
  count(puberty_yw) %>%
ggplot(aes(x=puberty_yw, y=n)) +
  geom_col(fill=tb_colors[2]) +
  tb_theme() +
  labs(
    title = "Timeseries Puberty Date Counts",
    subtitle = "Gilts reaching puberty per year-week",
    x = "Puberty Year-Week",
    y = "Count",
    caption = "tidybreed - 'Daily' Swine Breeding Program"
  ) + 
  theme(axis.text.x = element_text(angle=45, hjust=1))

#------------------------------------------------------------------------------#
# Saved figures: genetic trend by birth year-quarter
#------------------------------------------------------------------------------#

# A time-based simulation has no "generation"; birth year-quarter is the
# grouping used for trend figures written to config$output$save_dir.

warning("Save trend figures to: ", config$output$save_dir)

# helper: save a plot with consistent settings
save_fig <- function(p, name) {
  ggsave(
    filename = file.path(config$output$save_dir, name),
    plot   = p,
    width  = 8,
    height = 5,
    units  = "in",
    dpi    = 100,
    bg     = "white"
  )
}

# ---------- MEAN EBV + TBV BY BIRTH YEAR-QUARTER ---------- #

p <- data.tbv.ebv %>%
  pivot_longer(cols = c(tbv_value, ebv_value), names_to = "value_type", values_to = "value") %>%
  filter(!is.na(value)) %>%
  mutate(value_type = recode(value_type, tbv_value = "TBV", ebv_value = "EBV")) %>%
  group_by(value_type, trait_name, birth_yq) %>%
  summarise(mean_value = mean(value), n = n(), .groups = "drop") %>%
ggplot(aes(x=birth_yq, y=mean_value, color=value_type, group=value_type)) +
  geom_hline(yintercept = 0, color="grey50", linetype=2) +
  geom_line() +
  geom_point() +
  facet_wrap(~trait_name, scales="free_y") +
  scale_color_manual("Value", values = c(EBV = tb_colors[2], TBV = tb_colors[1])) +
  tb_theme() +
  labs(
    title = "Mean EBV (latest) and TBV by Birth Year-Quarter",
    x = "Birth Year-Quarter",
    y = "Mean Breeding Value",
    caption = config$scenario_name
  ) +
  theme(axis.text.x = element_text(angle=70, hjust=1))

print(p)
save_fig(p, "mean_ebv_tbv_on_birth_yq_facet_trait.png")

# ---------- TBV DISTRIBUTION BY BIRTH YEAR-QUARTER ---------- #

p <- data.tbv.ebv %>%
  group_by(trait_name, birth_yq) %>%
  summarise(
    MinTBV = min(tbv_value),
    Q1TBV  = quantile(tbv_value, prob=0.25),
    Q2TBV  = quantile(tbv_value, prob=0.50),
    Q3TBV  = quantile(tbv_value, prob=0.75),
    MaxTBV = max(tbv_value),
    .groups = "drop"
  ) %>% 
ggplot(aes(x=birth_yq, group=1)) +
  geom_hline(aes(yintercept = 0), color="red3", linewidth=0.75, linetype=3) +
  geom_ribbon(aes(ymin=MinTBV, ymax=MaxTBV), fill=tb_colors[2], alpha=0.10) +
  geom_ribbon(aes(ymin=Q1TBV,  ymax=Q3TBV),  fill=tb_colors[2], alpha=0.40) +
  geom_line(aes(y=Q2TBV), color=tb_colors[2]) +
  facet_wrap(~ trait_name, scales="free_y") +
  tb_theme() +
  labs(
    title = "TBV Trends By Trait",
    subtitle = "Median + middle 50 percent + min/max, by birth year-quarter",
    x = "Birth Year-Quarter",
    y = "TBV",
    caption = config$scenario_name
  ) +
  theme(axis.text.x = element_text(angle=70, hjust=1))

print(p)
save_fig(p, "ribbon_tbv_on_birth_yq_facet_trait.png")

# ---------- ANIMAL COUNT BY BIRTH YEAR-QUARTER ---------- #

p <- data.ind.meta %>%
  count(birth_yq, sex) %>%
ggplot(aes(x=birth_yq, y=n, fill=sex)) +
  geom_col(position="dodge") +
  scale_fill_manual("Sex", values = c(F = "magenta3", M = "dodgerblue3")) +
  tb_theme() +
  labs(
    title = "Count Animals by Birth Year-Quarter and Sex",
    x = "Birth Year-Quarter",
    y = "Count",
    caption = config$scenario_name
  ) +
  theme(axis.text.x = element_text(angle=70, hjust=1))

print(p)
save_fig(p, "bar_animal_count_on_birth_yq_fill_sex.png")

# ---------- PHENOTYPE COUNT BY PHENOTYPE YEAR-QUARTER ---------- #

p <- pop %>%
  get_table("ind_phenotype") %>%
  collect() %>%
  mutate(pheno_yq = paste(year(pheno_date), quarter(pheno_date), sep = "_")) %>%
  count(phenotype_name, pheno_yq) %>%
ggplot(aes(x=pheno_yq, y=n, fill=phenotype_name)) +
  geom_col() +
  facet_wrap(~phenotype_name, scales="free_y") +
  scale_fill_manual("Phenotype", values = tb_colors) +
  tb_theme() +
  labs(
    title = "Phenotype Count by Year-Quarter Phenotyped",
    x = "Phenotype Year-Quarter",
    y = "Phenotype Count",
    caption = config$scenario_name
  ) +
  theme(axis.text.x = element_text(angle=70, hjust=1))

print(p)
save_fig(p, "bar_phenotype_count_on_pheno_yq_facet_phenotype.png")

#------------------------------------------------------------------------------#
# Total time + close
#------------------------------------------------------------------------------#

total_elapsed <- (proc.time() - time_start_total)["elapsed"]
message(sprintf("Simulation finished | total: %s", format_elapsed(total_elapsed)))

# close pop object for database
close_pop(pop)
