# Test script for founder haplotype generation
library(tidybreed)
library(dplyr)

# Test 1: Create population with uniform allele frequencies
cat("Test 1: Uniform allele frequencies\n")
cat("===================================\n")

pop1 <- initialize_genome(
  pop_name = "test_uniform",
  n_loci = 100,
  n_chr = 5,
  chr_len_Mb = 100,
  n_haplotypes = 50,
  min_allele_freq = 0.05,
  max_allele_freq = 0.95
)

# Inspect the population
print(pop1)

# Check genome_meta table
genome_meta <- get_table(pop1, "genome_meta") %>% collect()
cat("\nGenome metadata columns:", paste(colnames(genome_meta), collapse = ", "), "\n")
cat("Allele frequency range:", min(genome_meta$founder_allele_freq), "-",
    max(genome_meta$founder_allele_freq), "\n")

# Check founder_haplotypes table
founder_haps <- get_table(pop1, "founder_haplotypes") %>% collect()
cat("Founder haplotypes dimensions:", nrow(founder_haps), "x", ncol(founder_haps), "\n")
cat("First few haplotype IDs:", head(founder_haps$hap_id), "\n")

# Check observed allele frequencies match expected
observed_freqs <- founder_haps %>%
  select(starts_with("locus_")) %>%
  summarise(across(everything(), mean)) %>%
  as.numeric()

cat("\nExpected vs Observed allele frequencies:\n")
comparison <- data.frame(
  locus = 1:10,  # Just first 10 loci
  expected = genome_meta$founder_allele_freq[1:10],
  observed = observed_freqs[1:10]
)
print(comparison, digits = 3)

close_pop(pop1)

# Test 2: Create population with fixed 50% allele frequency
cat("\n\nTest 2: Fixed 50% allele frequency\n")
cat("====================================\n")

pop2 <- initialize_genome(
  pop_name = "test_fixed",
  n_loci = 100,
  n_chr = 5,
  chr_len_Mb = 100,
  n_haplotypes = 100,
  fixed_allele_freq = 0.5
)

print(pop2)

# Check all allele frequencies are 0.5
genome_meta2 <- get_table(pop2, "genome_meta") %>% collect()
cat("\nAll allele frequencies equal 0.5:",
    all(genome_meta2$founder_allele_freq == 0.5), "\n")

# Check observed frequencies
founder_haps2 <- get_table(pop2, "founder_haplotypes") %>% collect()
observed_freqs2 <- founder_haps2 %>%
  select(starts_with("locus_")) %>%
  summarise(across(everything(), mean)) %>%
  as.numeric()

cat("Observed allele frequency range:",
    round(min(observed_freqs2), 3), "-",
    round(max(observed_freqs2), 3), "\n")
cat("Mean observed allele frequency:",
    round(mean(observed_freqs2), 3), "\n")

close_pop(pop2)

# Test 3: Create population without founder haplotypes
cat("\n\nTest 3: No founder haplotypes\n")
cat("==============================\n")

pop3 <- initialize_genome(
  pop_name = "test_no_founders",
  n_loci = 100,
  n_chr = 5,
  chr_len_Mb = 100
)

print(pop3)

cat("\nTables:", paste(pop3$tables, collapse = ", "), "\n")
cat("Has founder_haplotypes table:",
    "founder_haplotypes" %in% pop3$tables, "\n")

# Check genome_meta does not have founder_allele_freq
genome_meta3 <- get_table(pop3, "genome_meta") %>% collect()
cat("Has founder_allele_freq column:",
    "founder_allele_freq" %in% colnames(genome_meta3), "\n")

close_pop(pop3)

cat("\nAll tests completed!\n")
