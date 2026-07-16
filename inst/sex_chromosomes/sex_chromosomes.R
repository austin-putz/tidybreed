
library(DBI)
library(tidyverse)
library(tidybreed)

pop <- open_pop(pop_name = "xy", db_name = ":memory:")
    
pop <- pop |> 
  define_genome(
    n_loci = 12, 
    n_chr = 3, 
    chr_names = c("1", "X", "Y"), 
    chr_len_Mb = 50
) 

pop |> get_table("genome_meta")
pop |> get_table("chr_meta")

# list chromsome 1 names
locus_names_1 <- pop |> get_table("genome_meta") |>
  filter(chr_name == "1") |>
  pull(locus_name)

# list X and Y chromosomes loci names
locus_names_X <- pop |> get_table("genome_meta") |>
  filter(chr_name == "X") |>
  pull(locus_name)
locus_names_Y <- pop |> get_table("genome_meta") |>
  filter(chr_name == "Y") |>
  pull(locus_name)

pop <- pop |>
  define_chr("X", copy_mode_M = "half", copy_mode_F = "full",
               hemi_parent = "parent_2", recombines = TRUE) |>
  define_chr("Y", copy_mode_M = "half", copy_mode_F = "none",
               hemi_parent = "parent_1", recombines = FALSE)

# updated the chr_meta table to fix the Y chromosome
pop |> get_table("chr_meta")

# define haplotypes and sample a few M/F individuals to test this
pop |>
  define_founder_haplotypes(n_haplotypes = 20, method = "fixed") |>
    get_table("founder_haplotypes") |>
  add_founders(n_males = 4, n_females = 4, line_name = "A")

pop |> get_table("ind_meta")

pop |> 
  get_table("ind_haplotype") |> 
  count(id_ind)

# count of chromosome 1 per individual
pop |> 
  get_table("ind_haplotype") |>
  filter(
    locus_name %in% locus_names_1
  ) |>
  count(id_ind)

# count of X chromosomes per individual
pop |> 
  get_table("ind_haplotype") |>
  filter(
    locus_name %in% locus_names_X
  ) |>
  count(id_ind)

# count of Y chromosomes per individual
pop |> 
  get_table("ind_haplotype") |>
  filter(
    locus_name %in% locus_names_Y
  ) |>
  count(id_ind)

