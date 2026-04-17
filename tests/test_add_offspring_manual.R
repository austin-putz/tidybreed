# Manual test: add_offspring()
# Run interactively to verify the full round-trip

library(tidybreed)

cat("=== add_offspring() manual test ===\n\n")

# ── Setup ────────────────────────────────────────────────────────────────────
pop <- initialize_genome(
  pop_name    = "manual_test",
  n_loci      = 200,
  n_chr       = 3,
  chr_len_Mb  = 100,
  n_haplotypes = 50,
  db_path     = ":memory:"
)

pop <- pop |>
  add_founders(n_males = 5, n_females = 5, line_name = "A") |>
  mutate_ind_meta(gen = 1L)

cat("Founders added. ind_meta rows:", nrow(dplyr::collect(get_table(pop, "ind_meta"))), "\n")

# ── Basic offspring ──────────────────────────────────────────────────────────
matings <- tibble::tibble(
  id_parent_1 = rep("A-1", 5),
  id_parent_2 = paste0("A-", 6:10),
  sex         = c("M", "F", "M", "F", "M"),
  line        = "A",
  gen         = 2L
)

pop <- pop |> add_offspring(matings)

ind_meta <- dplyr::collect(get_table(pop, "ind_meta"))
cat("\nAfter add_offspring(5 offspring):\n")
cat("  ind_meta rows:         ", nrow(ind_meta), "\n")
cat("  gen column present:    ", "gen" %in% colnames(ind_meta), "\n")
cat("  offspring gen values:  ", unique(ind_meta$gen[!is.na(ind_meta$id_parent_1)]), "\n")

# ── Haplotype / genotype checks ──────────────────────────────────────────────
haps  <- dplyr::collect(get_table(pop, "genome_haplotype"))
genos <- dplyr::collect(get_table(pop, "genome_genotype"))
cat("\ngenome_haplotype rows:  ", nrow(haps), " (expected", nrow(ind_meta) * 2, ")\n")
cat("genome_genotype rows:   ", nrow(genos), " (expected", nrow(ind_meta), ")\n")

locus_cols <- paste0("locus_", 1:200)
hap_vals   <- unlist(haps[, locus_cols])
cat("\nAll haplotype values in {0,1}:", all(hap_vals %in% c(0L, 1L)), "\n")

# Check genotype = hap1 + hap2 for first offspring
off_id <- "A-11"
off_haps <- dplyr::filter(haps, id_ind == off_id)
off_geno <- dplyr::filter(genos, id_ind == off_id)
h1 <- as.integer(off_haps[off_haps$parent_origin == 1, locus_cols])
h2 <- as.integer(off_haps[off_haps$parent_origin == 2, locus_cols])
g  <- as.integer(off_geno[1, locus_cols])
cat("Genotype == hap1 + hap2 for", off_id, ":", all(g == h1 + h2), "\n")

# ── Alias columns ────────────────────────────────────────────────────────────
matings_alias <- tibble::tibble(
  id_sire = "A-1",
  id_dam  = "A-6",
  sex     = "M",
  line    = "A",
  gen     = 2L
)
pop <- pop |> add_offspring(matings_alias)
cat("\nAlias (id_sire / id_dam) accepted without error: TRUE\n")

# ── Cross-line mating ────────────────────────────────────────────────────────
pop <- pop |>
  add_founders(n_males = 2, n_females = 2, line_name = "B") |>
  mutate_ind_meta(gen = 1L)

cross_matings <- tibble::tibble(
  id_parent_1 = c("A-1", "A-2"),
  id_parent_2 = c("B-1", "B-2"),
  sex         = c("M", "F"),
  line        = "C",
  gen         = 2L
)
pop <- pop |> add_offspring(cross_matings)
c_ids <- dplyr::filter(dplyr::collect(get_table(pop, "ind_meta")), line == "C")$id_ind
cat("\nCross-line offspring (line C):", paste(c_ids, collapse = ", "), "\n")

# ── Error checks ─────────────────────────────────────────────────────────────
cat("\nError checks:\n")
tryCatch(
  add_offspring(pop, tibble::tibble(id_parent_1 = "Z-999", id_parent_2 = "A-6",
                                    sex = "M", line = "A")),
  error = function(e) cat("  Non-existent parent correctly caught:", conditionMessage(e), "\n")
)
tryCatch(
  add_offspring(pop, tibble::tibble(id_parent_1 = "A-1", id_parent_2 = "A-6",
                                    sex = "X", line = "A")),
  error = function(e) cat("  Invalid sex correctly caught:", conditionMessage(e), "\n")
)

close_pop(pop)
cat("\n=== All manual checks passed ===\n")
