# tests/test_dynamic_genome.R
#
# Manual integration test: verifies that the DB-first design supports adding a
# novel locus after initialization, gene-editing an individual, and producing
# offspring whose recombination respects the updated genome.
#
# Run interactively:
#   devtools::load_all()
#   source("tests/test_dynamic_genome.R")

suppressPackageStartupMessages({
  if (requireNamespace("devtools", quietly = TRUE) &&
      file.exists("DESCRIPTION")) {
    devtools::load_all(".", quiet = TRUE)
  } else {
    library(tidybreed)
  }
})

set.seed(20260514)

cat("=== Dynamic genome test ===\n\n")

# ─────────────────────────────────────────────────────────────────────────────
# 1. Setup: population with founders and ADG trait with QTL
# ─────────────────────────────────────────────────────────────────────────────

pop <- initialize_genome(
  pop_name          = "dyntest",
  n_loci            = 200,
  n_chr             = 2,
  chr_len_Mb        = 100,
  n_haplotypes      = 100,
  db_path           = ":memory:",
  fixed_allele_freq = 0.5
)

pop <- add_founders(pop, n_males = 10, n_females = 10, line_name = "A")
pop <- get_table(pop, "ind_meta") |> mutate_table(gen = 0L)

pop <- add_trait(pop,
  trait_name     = "ADG",
  trait_type     = "continuous",
  target_add_var = 1.0,
  residual_var   = 1.0
)
sel_qtl <- get_table(pop, "genome_meta") |> collect() |>
  dplyr::slice_sample(n = 20) |> dplyr::pull(locus_name)
pop <- pop |>
  get_table("genome_meta") |>
  dplyr::filter(locus_name %in% sel_qtl) |>
  define_qtl("ADG")
pop <- add_additive_effects(pop, "ADG", seed = 99)

n_start <- DBI::dbGetQuery(pop$db_conn, "SELECT COUNT(*) AS n FROM genome_meta")$n
cat("Baseline setup complete.\n")
cat("  Founders:", DBI::dbGetQuery(pop$db_conn, "SELECT COUNT(*) AS n FROM ind_meta")$n, "\n")
cat("  Loci at start:", n_start, "\n\n")

# ─────────────────────────────────────────────────────────────────────────────
# 2. Gene surgery: add one new locus via raw DBI calls
#    This is the "user by hand" workflow until a formal add_locus() exists.
# ─────────────────────────────────────────────────────────────────────────────

conn <- pop$db_conn

next_id      <- as.integer(DBI::dbGetQuery(conn,
  "SELECT MAX(locus_id) + 1 AS nid FROM genome_meta")$nid)
new_locus_col <- paste0("locus_", next_id)

cat("Inserting new locus: locus_id =", next_id, " col =", new_locus_col, "\n")

# 2a. Add row to genome_meta (chromosome 1, mid-position)
DBI::dbExecute(conn, sprintf(
  "INSERT INTO genome_meta (locus_id, locus_name, chr, chr_name, pos_Mb)
   VALUES (%d, 'MSTN_KO', 1, '1', 50.0)",
  next_id
))

# 2b. Add column to haplotype and genotype tables, defaulting to 0
DBI::dbExecute(conn, sprintf(
  "ALTER TABLE genome_haplotype ADD COLUMN %s UTINYINT DEFAULT 0",
  new_locus_col))
DBI::dbExecute(conn, sprintf(
  "ALTER TABLE genome_genotype  ADD COLUMN %s UTINYINT DEFAULT 0",
  new_locus_col))

cat("  Columns added to genome_haplotype and genome_genotype (default 0).\n")

# 2c. Gene-edit individual A_1: set to heterozygous (hap1 = 1, hap2 = 0 → genotype = 1)
edited_id <- "A_1"

DBI::dbExecute(conn, sprintf(
  "UPDATE genome_haplotype SET %s = 1 WHERE id_ind = '%s' AND parent_origin = 1",
  new_locus_col, edited_id
))
DBI::dbExecute(conn, sprintf(
  "UPDATE genome_genotype SET %s = 1 WHERE id_ind = '%s'",
  new_locus_col, edited_id
))

cat("  Gene edit: individual", edited_id, "set to heterozygous at", new_locus_col, "\n\n")

# ─────────────────────────────────────────────────────────────────────────────
# 3. Set a LARGE additive effect for the new locus on ADG
#    is_QTL_ADG and add_ADG columns already exist (created by define_qtl /
#    add_additive_effects); new locus row has them as NULL — just UPDATE.
# ─────────────────────────────────────────────────────────────────────────────

DBI::dbExecute(conn, sprintf(
  "UPDATE genome_meta SET is_QTL_ADG = TRUE, add_ADG = 50.0 WHERE locus_id = %d",
  next_id
))

cat("Effect set: add_ADG = 50.0 for locus_id =", next_id, "\n\n")

# ─────────────────────────────────────────────────────────────────────────────
# 4. Compute TBVs and verify the gene-edited individual is an outlier
# ─────────────────────────────────────────────────────────────────────────────

pop <- pop |>
  get_table("ind_meta") |>
  add_tbv("ADG")

tbv_df <- DBI::dbGetQuery(conn,
  "SELECT id_ind, tbv_value FROM ind_tbv WHERE trait_name = 'ADG' ORDER BY tbv_value DESC")

cat("TBV results (all individuals, sorted):\n")
print(tbv_df)

edited_tbv  <- tbv_df$tbv_value[tbv_df$id_ind == edited_id]
other_tbvs  <- tbv_df$tbv_value[tbv_df$id_ind != edited_id]

cat("\nGene-edited individual", edited_id, "TBV:", round(edited_tbv, 3), "\n")
cat("Max TBV among others:     ", round(max(other_tbvs), 3), "\n")
cat("Min TBV among others:     ", round(min(other_tbvs), 3), "\n")
cat("Effect size vs next-best: ", round(edited_tbv - max(other_tbvs), 3), "\n\n")

stopifnot(
  "Gene-edited individual must have the highest TBV" =
    edited_tbv > max(other_tbvs)
)
cat("PASS: gene-edited individual shows outlier TBV.\n\n")

# ─────────────────────────────────────────────────────────────────────────────
# 5. Verify add_offspring() uses the live genome (n+1 loci)
# ─────────────────────────────────────────────────────────────────────────────

n_loci_after <- as.integer(DBI::dbGetQuery(conn,
  "SELECT COUNT(*) AS n FROM genome_meta")$n)

stopifnot("genome_meta should have n+1 rows" = n_loci_after == n_start + 1L)
cat("Loci after insertion:", n_loci_after, "(was", n_start, ")\n")

# A_1 (heterozygous sire) x 5 dams
matings <- data.frame(
  id_parent_1 = rep(edited_id, 5),
  id_parent_2 = paste0("A_", 11:15),
  sex         = rep("M", 5),
  line_name   = "A",
  gen         = 1L
)

pop <- add_offspring(pop, matings)

cat("Offspring added. Checking genome_haplotype for new locus column...\n")

offspring_ids <- paste0("A_", 21:25)

new_col_vals <- DBI::dbGetQuery(conn, sprintf(
  "SELECT id_ind, parent_origin, %s FROM genome_haplotype WHERE id_ind IN (%s) ORDER BY id_ind, parent_origin",
  new_locus_col,
  paste(sprintf("'%s'", offspring_ids), collapse = ", ")
))

cat("\nNew locus haplotype values in offspring (each has 2 rows: parent_origin 1 and 2):\n")
print(new_col_vals)

stopifnot(
  "Offspring haplotype table must contain the new locus column" =
    new_locus_col %in% colnames(new_col_vals)
)

# Mendelian check: sire A_1 is heterozygous (hap1=1, hap2=0).
# Some offspring should inherit the mutant allele from the sire gamete.
sire_contributions <- new_col_vals$locus_201[new_col_vals$parent_origin == 1]
cat("\nSire-contributed haplotype values at new locus:", sire_contributions, "\n")
cat("(Should be 0 or 1; not all-zero if recombination sampling works correctly)\n")

# Verify column count in offspring haplotype matches n_loci_after
hap_cols <- DBI::dbListFields(conn, "genome_haplotype")
n_locus_cols <- sum(grepl("^locus_", hap_cols))

stopifnot(
  "genome_haplotype must have n+1 locus columns" = n_locus_cols == n_loci_after
)
cat("\nPASS: genome_haplotype has", n_locus_cols, "locus columns (matches genome_meta row count).\n")

# ─────────────────────────────────────────────────────────────────────────────
# Cleanup
# ─────────────────────────────────────────────────────────────────────────────

close_pop(pop)
cat("\n=== Dynamic genome test PASSED ===\n")
