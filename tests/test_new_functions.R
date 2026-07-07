#!/usr/bin/env Rscript
# Quick integration test for mutate_table() on genome_meta and define_chip()

library(dplyr)

cat("=== Testing mutate_table() on genome_meta and define_chip() ===\n\n")

# Initialize genome (current long-format API)
cat("1. Initializing genome...\n")
pop <- open_pop(pop_name = "test", db_name = ":memory:") |>
  define_genome(n_loci = 1000, n_chr = 10, chr_len_Mb = 100) |>
  define_founder_haplotypes(n_haplotypes = 100)

# Define SNP chips
cat("2. Defining SNP chips...\n")
set.seed(42)
sel_50k <- get_table(pop, "genome_meta") |> collect() |>
  slice_sample(n = 500) |> pull(locus_name)
pop <- pop |> get_table("genome_meta") |>
  filter(locus_name %in% sel_50k) |> define_chip("50k")

sel_HD <- get_table(pop, "genome_meta") |> collect() |>
  slice_sample(n = 900) |> pull(locus_name)
pop <- pop |> get_table("genome_meta") |>
  filter(locus_name %in% sel_HD) |> define_chip("HD")

sel_10k <- get_table(pop, "genome_meta") |> collect() |>
  slice_sample(n = 100) |> pull(locus_name)
pop <- pop |> get_table("genome_meta") |>
  filter(locus_name %in% sel_10k) |> define_chip("10k")

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

stopifnot(sum(genome$is_50k) == 500)
stopifnot(sum(genome$is_HD) == 900)
stopifnot(sum(genome$is_10k) == 100)

# Add custom columns to genome_meta via mutate_table()
# Per-row values must be passed as a tibble keyed on the primary key (locus_id);
# plain vectors are rejected because DB row order is not guaranteed.
cat("4. Adding custom columns to genome_meta via mutate_table()...\n")
n_loci <- nrow(genome)
set.seed(123)
effect_vec <- ifelse(genome$is_50k, rnorm(n_loci, 0, 1), 0)

pop <- pop %>%
  get_table("genome_meta") %>%
  mutate_table(
    effect_50k = tibble::tibble(locus_id = genome$locus_id,
                                effect_50k = effect_vec)
  )
pop <- pop %>%
  get_table("genome_meta") %>%
  mutate_table(
    is_QTL_growth = tibble::tibble(locus_id = genome$locus_id,
                                   is_QTL_growth = genome$is_50k)  # chip SNPs as QTL
  )

# Add founders (current API: pipe founder_haplotypes into add_founders)
cat("5. Adding founders...\n")
pop <- pop %>%
  get_table("founder_haplotypes") %>%
  add_founders(n_males = 10, n_females = 100, line_name = "A")

# Mark individuals via mutate_table() on ind_meta
cat("6. Marking individuals with custom columns via mutate_table()...\n")
pop <- pop %>%
  get_table("ind_meta") %>%
  mutate_table(
    genotyped_50k = TRUE,
    gen           = 0L
  )

# Query chip genotypes
cat("7. Querying chip genotypes...\n")
chip_loci <- genome %>%
  filter(is_50k == TRUE) %>%
  pull(locus_name)

cat("  First 10 chip loci:", paste(head(chip_loci, 10), collapse = ", "), "\n")

# Materialize dosages for chip SNPs (on-demand cache: ind_genotype is empty
# until add_dosage() populates it), then query them in long format.
pop <- pop %>%
  get_table("ind_meta") %>%
  add_dosage(chip_name = "50k")

chip_genotypes <- get_table(pop, "ind_genotype") %>%
  filter(locus_name %in% chip_loci) %>%
  collect()

n_ind  <- length(unique(chip_genotypes$id_ind))
n_loc  <- length(unique(chip_genotypes$locus_name))
cat("  Genotype rows:", nrow(chip_genotypes),
    "| individuals:", n_ind, "| loci:", n_loc, "\n")

stopifnot(n_ind == 110)
stopifnot(n_loc == 500)
stopifnot(all(chip_genotypes$dosage_value %in% c(0L, 1L, 2L)))

# Summary
cat("\n=== Summary ===\n")
cat("get_table('genome_meta') |> mutate_table() working correctly\n")
cat("define_chip() working correctly\n")
cat("Integration with add_founders() working\n")
cat("get_table('ind_meta') |> mutate_table() working\n")
cat("Can query chip genotypes successfully (long-format ind_genotype)\n")

# View final genome_meta structure
cat("\n8. Final genome_meta structure:\n")
result <- get_table(pop, "genome_meta") %>%
  filter(is_50k == TRUE) %>%
  select(locus_id, locus_name, chr, pos_bp, is_50k, is_HD, is_10k, effect_50k, is_QTL_growth) %>%
  collect()

print(head(result, 10))

cat("\n=== All tests passed! ===\n")

# Clean up
close_pop(pop)
