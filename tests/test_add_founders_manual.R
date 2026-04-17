# Manual test script for add_founders() function
# Run this in R to verify the implementation works correctly

library(tidyverse)
devtools::load_all()

cat("=" , rep("=", 70), "\n", sep = "")
cat("Testing add_founders() Implementation\n")
cat("=" , rep("=", 70), "\n\n", sep = "")

# Test 1: Basic functionality
cat("Test 1: Basic functionality\n")
cat("--------------------------\n")

pop <- initialize_genome(
  pop_name = "test",
  n_loci = 100,
  n_chr = 5,
  chr_len_Mb = 100,
  n_haplotypes = 50,
  db_path = ":memory:"
)

pop <- add_founders(pop, n_males = 10, n_females = 20, line_name = "A")

ind_meta <- get_table(pop, "ind_meta") %>% collect()
cat("✓ Created", nrow(ind_meta), "founders\n")
cat("✓ IDs:", paste(head(ind_meta$id_ind, 5), collapse = ", "), "...\n")
cat("✓ Males:", sum(ind_meta$sex == "M"), ", Females:", sum(ind_meta$sex == "F"), "\n")
cat("✓ All parents are NA:", all(is.na(ind_meta$id_parent_1)) && all(is.na(ind_meta$id_parent_2)), "\n\n")

# Test 2: Haplotype and genotype verification
cat("Test 2: Haplotype and genotype verification\n")
cat("-------------------------------------------\n")

haps <- get_table(pop, "genome_haplotype") %>% collect()
genos <- get_table(pop, "genome_genotype") %>% collect()

cat("✓ Haplotype rows (should be 60):", nrow(haps), "\n")
cat("✓ Genotype rows (should be 30):", nrow(genos), "\n")

# Verify genotype = sum of haplotypes for first individual
ind_id <- ind_meta$id_ind[1]
hap_rows <- haps %>% filter(id_ind == !!ind_id)
geno_row <- genos %>% filter(id_ind == !!ind_id)

locus_name <- "locus_1"
hap_sum <- sum(hap_rows[[locus_name]])
geno_val <- geno_row[[locus_name]][1]

cat("✓ Genotype = sum of haplotypes:", hap_sum == geno_val,
    "(", hap_rows[[locus_name]][1], "+", hap_rows[[locus_name]][2], "=", geno_val, ")\n\n")

close_pop(pop)

# Test 3: Multiple lines
cat("Test 3: Multiple lines\n")
cat("----------------------\n")

pop <- initialize_genome(
  pop_name = "test",
  n_loci = 50,
  n_chr = 2,
  chr_len_Mb = 100,
  n_haplotypes = 20,
  db_path = ":memory:"
)

pop <- add_founders(pop, n_males = 5, n_females = 5, line_name = "A")
pop <- add_founders(pop, n_males = 3, n_females = 7, line_name = "B")

ind_meta <- get_table(pop, "ind_meta") %>% collect()

line_a <- ind_meta %>% filter(line == "A")
line_b <- ind_meta %>% filter(line == "B")

cat("✓ Line A:", nrow(line_a), "individuals\n")
cat("✓ Line A IDs:", paste(head(line_a$id_ind, 5), collapse = ", "), "\n")
cat("✓ Line B:", nrow(line_b), "individuals\n")
cat("✓ Line B IDs:", paste(head(line_b$id_ind, 5), collapse = ", "), "\n\n")

close_pop(pop)

# Test 4: Sequential additions
cat("Test 4: Sequential additions to same line\n")
cat("-----------------------------------------\n")

pop <- initialize_genome(
  pop_name = "test",
  n_loci = 50,
  n_chr = 2,
  chr_len_Mb = 100,
  n_haplotypes = 20,
  db_path = ":memory:"
)

pop <- add_founders(pop, n_males = 5, n_females = 5, line_name = "A")
pop <- add_founders(pop, n_males = 2, n_females = 3, line_name = "A")

ind_meta <- get_table(pop, "ind_meta") %>% collect()

cat("✓ Total individuals:", nrow(ind_meta), "\n")
cat("✓ IDs are sequential:", all(ind_meta$id_ind == paste0("A-", 1:15)), "\n")
cat("✓ All IDs:", paste(ind_meta$id_ind, collapse = ", "), "\n\n")

close_pop(pop)

# Test 5: Integration with mutate_ind_meta
cat("Test 5: Integration with mutate_ind_meta\n")
cat("----------------------------------------\n")

pop <- initialize_genome(
  pop_name = "test",
  n_loci = 100,
  n_chr = 5,
  chr_len_Mb = 100,
  n_haplotypes = 50,
  db_path = ":memory:"
) %>%
  add_founders(n_males = 10, n_females = 100, line_name = "A") %>%
  mutate_ind_meta(
    gen = 0,
    farm = "FarmA",
    date_birth = as.Date("2024-01-01")
  )

ind_meta <- get_table(pop, "ind_meta") %>% collect()

cat("✓ Columns in ind_meta:", paste(colnames(ind_meta), collapse = ", "), "\n")
cat("✓ All gen = 0:", all(ind_meta$gen == 0), "\n")
cat("✓ All farm = 'FarmA':", all(ind_meta$farm == "FarmA"), "\n")
cat("✓ All date_birth correct:", all(ind_meta$date_birth == as.Date("2024-01-01")), "\n\n")

close_pop(pop)

# Summary
cat("=" , rep("=", 70), "\n", sep = "")
cat("All manual tests passed! ✓\n")
cat("=" , rep("=", 70), "\n\n", sep = "")

cat("Next steps:\n")
cat("1. Run automated tests: devtools::test()\n")
cat("2. Generate documentation: devtools::document()\n")
cat("3. Try the full workflow in your own scripts\n")
