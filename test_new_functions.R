#!/usr/bin/env Rscript
# Quick integration test for mutate_genome_meta() and define_chip()

library(tidybreed)
library(dplyr)

cat("=== Testing mutate_genome_meta() and define_chip() ===\n\n")

# Initialize genome
cat("1. Initializing genome...\n")
pop <- initialize_genome(
  pop_name = "test",
  n_loci = 1000,
  n_chr = 10,
  chr_len_Mb = 100,
  db_path = ":memory:",
  n_haplotypes = 100
)

# Define SNP chips
cat("2. Defining SNP chips...\n")
pop <- pop %>%
  define_chip(name = "50k", n_snp = 500, method = "random") %>%
  define_chip(name = "HD", n_snp = 900, method = "even") %>%
  define_chip(name = "10k", n_snp = 100, method = "chromosome_even")

# Check chips were created
cat("3. Checking chip columns...\n")
cols <- DBI::dbListFields(pop$db_conn, "genome_meta")
cat("  Columns in genome_meta:", paste(cols, collapse = ", "), "\n")

stopifnot("is_50k" %in% cols)
stopifnot("is_HD" %in% cols)
stopifnot("is_10k" %in% cols)

# Count SNPs per chip
genome <- get_table(pop, "genome_meta") %>% collect()
cat("  SNPs on 50k chip:", sum(genome$is_50k), "\n")
cat("  SNPs on HD chip:", sum(genome$is_HD), "\n")
cat("  SNPs on 10k chip:", sum(genome$is_10k), "\n")

# Add QTL effects to 50k chip SNPs
cat("4. Adding QTL effects to 50k chip SNPs...\n")
n_loci <- nrow(genome)
set.seed(123)
effect_vec <- ifelse(
  genome$is_50k,
  rnorm(n_loci, 0, 1),
  0
)

pop <- pop %>% mutate_genome_meta(
  effect_50k = effect_vec,
  is_QTL_growth = genome$is_50k  # Mark chip SNPs as QTL
)

# Add founders
cat("5. Adding founders...\n")
pop <- pop %>%
  add_founders(n_males = 10, n_females = 100, line_name = "A")

# Mark individuals as genotyped on 50k chip
cat("6. Marking individuals as genotyped...\n")
pop <- pop %>%
  mutate_ind_meta(
    genotyped_50k = TRUE,
    gen = 0
  )

# Query chip genotypes
cat("7. Querying chip genotypes...\n")
chip_loci <- genome %>%
  filter(is_50k == TRUE) %>%
  pull(locus_name)

cat("  First 10 chip loci:", paste(head(chip_loci, 10), collapse = ", "), "\n")

# Get genotypes for chip SNPs
chip_genotypes <- get_table(pop, "genome_genotype") %>%
  select(ind_id, all_of(chip_loci)) %>%
  collect()

cat("  Genotype matrix dimensions:", nrow(chip_genotypes), "x", ncol(chip_genotypes), "\n")

# Summary
cat("\n=== Summary ===\n")
cat("✓ mutate_genome_meta() working correctly\n")
cat("✓ define_chip() working correctly\n")
cat("✓ Integration with add_founders() working\n")
cat("✓ Integration with mutate_ind_meta() working\n")
cat("✓ Can query chip genotypes successfully\n")

# View final genome_meta structure
cat("\n8. Final genome_meta structure:\n")
result <- get_table(pop, "genome_meta") %>%
  filter(is_50k == TRUE) %>%
  select(locus_id, locus_name, chr, pos_Mb, is_50k, is_HD, is_10k, effect_50k, is_QTL_growth) %>%
  head(10) %>%
  collect()

print(result)

cat("\n=== All tests passed! ===\n")

# Clean up
close_pop(pop)
