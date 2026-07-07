# Test script for founder haplotype generation
library(dplyr)

# Test 1: Create population with uniform allele frequencies
cat("Test 1: Uniform allele frequencies\n")
cat("===================================\n")

pop1 <- open_pop(pop_name = "test_uniform", db_name = ":memory:") |>
  define_genome(
    n_loci = 100,
    n_chr = 5,
    chr_len_Mb = 100
  ) |>
  define_founder_haplotypes(
    n_haplotypes    = 50,
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

# Check founder_haplotypes table (LONG: line_name, haplotype_id, locus_name, allele)
founder_haps <- get_table(pop1, "founder_haplotypes") %>% collect()
cat("Founder haplotypes rows x cols:", nrow(founder_haps), "x", ncol(founder_haps), "\n")
cat("Columns:", paste(colnames(founder_haps), collapse = ", "), "\n")
cat("Unique haplotype IDs:", length(unique(founder_haps$haplotype_id)), "\n")

# Observed per-locus allele frequency = mean allele across haplotypes
observed_freqs <- founder_haps %>%
  group_by(locus_name) %>%
  summarise(observed = mean(allele), .groups = "drop")

cat("\nExpected vs Observed allele frequencies (first 10 loci):\n")
comparison <- genome_meta %>%
  arrange(locus_id) %>%
  slice_head(n = 10) %>%
  select(locus = locus_id, locus_name, expected = founder_allele_freq) %>%
  left_join(observed_freqs, by = "locus_name")
print(comparison, digits = 3)

close_pop(pop1)

# Test 2: Create population with fixed 50% allele frequency
cat("\n\nTest 2: Fixed 50% allele frequency\n")
cat("====================================\n")

pop2 <- open_pop(pop_name = "test_fixed", db_name = ":memory:") |>
  define_genome(
    n_loci = 100,
    n_chr = 5,
    chr_len_Mb = 100
  ) |>
  define_founder_haplotypes(
    n_haplotypes      = 100,
    method = "fixed"
  )

print(pop2)

# Check all allele frequencies are 0.5
genome_meta2 <- get_table(pop2, "genome_meta") %>% collect()
cat("\nAll allele frequencies equal 0.5:",
    all(genome_meta2$founder_allele_freq == 0.5), "\n")

# Check observed frequencies (LONG format)
founder_haps2 <- get_table(pop2, "founder_haplotypes") %>% collect()
observed_freqs2 <- founder_haps2 %>%
  group_by(locus_name) %>%
  summarise(freq = mean(allele), .groups = "drop") %>%
  pull(freq)

cat("Observed allele frequency range:",
    round(min(observed_freqs2), 3), "-",
    round(max(observed_freqs2), 3), "\n")
cat("Mean observed allele frequency:",
    round(mean(observed_freqs2), 3), "\n")

close_pop(pop2)

# Test 3: Create population without founder haplotypes
cat("\n\nTest 3: No founder haplotypes\n")
cat("==============================\n")

pop3 <- open_pop(pop_name = "test_no_founders", db_name = ":memory:") |>
  define_genome(
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
