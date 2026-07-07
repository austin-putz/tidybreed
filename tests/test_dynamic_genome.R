# tests/test_dynamic_genome.R
#
# Manual integration test: verifies that the DB-first design supports adding a
# novel locus after initialization, gene-editing an individual, and producing
# offspring whose recombination respects the updated genome — all expressed in
# the CURRENT long-format schema (ind_haplotype / genome_map / genome_effects).
#
# The old wide-format "ALTER TABLE ADD COLUMN locus_N" premise is gone. Adding a
# locus is now: one row in genome_meta, one row in genome_map (its genetic
# position), one ind_haplotype row per allele copy for every existing individual,
# and one genome_effects row to give it a QTL effect.
#
# Run interactively:
#   Rscript -e 'suppressMessages(pkgload::load_all(".", quiet=TRUE)); source("tests/test_dynamic_genome.R")'

suppressPackageStartupMessages(library(dplyr))

set.seed(20260514)

cat("=== Dynamic genome test ===\n\n")

# ─────────────────────────────────────────────────────────────────────────────
# 1. Setup: population with founders and an ADG trait with QTL
# ─────────────────────────────────────────────────────────────────────────────

pop <- open_pop(pop_name = "dyntest", db_name = ":memory:") |>
  define_genome(n_loci = 200, n_chr = 2, chr_len_Mb = 100) |>
  define_founder_haplotypes(n_haplotypes = 100, method = "fixed")

pop <- pop |>
  get_table("founder_haplotypes") |>
  add_founders(n_males = 10, n_females = 10, line_name = "A")
pop <- get_table(pop, "ind_meta") |> mutate_table(gen = 0L)

# Genetic layer: register the trait, then place ~20 random QTL and draw effects.
pop <- define_trait(pop, trait_name = "ADG", target_add_var = 1.0)

sel_qtl <- get_table(pop, "genome_meta") |> collect() |>
  slice_sample(n = 20) |> pull(locus_name)
pop <- pop |>
  get_table("genome_meta") |>
  filter(locus_name %in% sel_qtl) |>
  define_additive_effects("ADG", seed = 99)

conn    <- pop$db_conn
n_start <- DBI::dbGetQuery(conn, "SELECT COUNT(*) AS n FROM genome_meta")$n
cat("Baseline setup complete.\n")
cat("  Founders:", DBI::dbGetQuery(conn, "SELECT COUNT(*) AS n FROM ind_meta")$n, "\n")
cat("  Loci at start:", n_start, "\n\n")

# ─────────────────────────────────────────────────────────────────────────────
# 2. Gene surgery: add one new locus in long form via raw DBI calls.
#    This is the "user by hand" workflow until a formal add_locus() exists.
# ─────────────────────────────────────────────────────────────────────────────

next_id     <- as.integer(DBI::dbGetQuery(conn,
  "SELECT MAX(locus_id) + 1 AS nid FROM genome_meta")$nid)
new_locus   <- "MSTN_KO"
edited_id   <- "A_1"

cat("Inserting new locus:", new_locus, " (locus_id =", next_id, ")\n")

# 2a. genome_meta row. The new locus sits at the terminal end of chromosome 1.
#     resolve_genome_map() orders loci by locus_id and requires pos_cM to be
#     monotonic non-decreasing within a chromosome; because the new locus_id is
#     the largest, it sorts LAST within chr 1, so its position must be >= the
#     current chr-1 maximum. We place it at that maximum.
chr1_max <- DBI::dbGetQuery(conn,
  "SELECT MAX(pos_bp) AS mbp FROM genome_meta WHERE chr = 1")
new_pos_bp <- as.numeric(chr1_max$mbp)

DBI::dbExecute(conn, sprintf(
  "INSERT INTO genome_meta (locus_id, locus_name, chr, chr_name, pos_bp)
   VALUES (%d, '%s', 1, '1', %d)",
  next_id, new_locus, as.integer(new_pos_bp)))

# 2b. genome_map row (its genetic position). pos_cM matches the chr-1 maximum so
#     the resolved map stays monotonic.
chr1_max_cM <- DBI::dbGetQuery(conn, "
  SELECT MAX(gm.pos_cM) AS mcm
  FROM genome_map gm JOIN genome_meta m ON gm.locus_id = m.locus_id
  WHERE m.chr = 1 AND gm.map_name = 'default'")$mcm
new_map_id <- next_int_id(conn, "genome_map", "id_genome_map")
DBI::dbWriteTable(conn, "genome_map",
  data.frame(id_genome_map = new_map_id, locus_id = next_id,
             locus_name = new_locus, sex = NA_character_,
             line_name = NA_character_, map_name = "default",
             pos_cM = chr1_max_cM, stringsAsFactors = FALSE),
  append = TRUE)

# 2c. ind_haplotype rows for the new locus: every existing individual gets one
#     row per (parent_origin, strand) copy, allele = 0, line_origin = its line.
inds <- DBI::dbGetQuery(conn, "SELECT id_ind, line_name FROM ind_meta")
new_haps <- do.call(rbind, lapply(c(1L, 2L), function(po) {
  data.frame(id_ind = inds$id_ind, parent_origin = po, strand = 1L,
             line_origin = inds$line_name, locus_id = next_id,
             locus_name = new_locus, allele = 0L, stringsAsFactors = FALSE)
}))
DBI::dbWriteTable(conn, "ind_haplotype", new_haps, append = TRUE)

cat("  Added genome_meta + genome_map + ind_haplotype rows for the new locus.\n")

# 2d. Gene-edit A_1: set its parent_origin = 1 copy to the mutant allele (1),
#     leaving parent_origin = 2 at 0 → heterozygous carrier.
DBI::dbExecute(conn, sprintf(
  "UPDATE ind_haplotype SET allele = 1
   WHERE id_ind = '%s' AND locus_id = %d AND parent_origin = 1",
  edited_id, next_id))

cat("  Gene edit:", edited_id, "set to heterozygous (mutant on parent_origin = 1).\n\n")

# ─────────────────────────────────────────────────────────────────────────────
# 3. Give the new locus a LARGE additive effect on ADG (population-wide row).
# ─────────────────────────────────────────────────────────────────────────────

new_eff_id <- next_int_id(conn, "genome_effects", "id_genome_effect")
DBI::dbWriteTable(conn, "genome_effects",
  data.frame(id_genome_effect = new_eff_id, locus_name = new_locus,
             line_name = NA_character_, trait_name = "ADG",
             genome_effect_type = "additive", genome_value = 50.0,
             base_allele_freq = 0.5, stringsAsFactors = FALSE),
  append = TRUE)

cat("Effect set: genome_value = 50.0 for", new_locus, "on ADG.\n\n")

# ─────────────────────────────────────────────────────────────────────────────
# 4. Compute TBVs and verify the gene-edited individual is an outlier.
# ─────────────────────────────────────────────────────────────────────────────

pop <- pop |>
  get_table("ind_meta") |>
  add_tbv("ADG")

tbv_df <- DBI::dbGetQuery(conn,
  "SELECT id_ind, tbv_value FROM ind_tbv WHERE trait_name = 'ADG' ORDER BY tbv_value DESC")

cat("TBV results (all individuals, sorted):\n")
print(tbv_df)

edited_tbv <- tbv_df$tbv_value[tbv_df$id_ind == edited_id]
other_tbvs <- tbv_df$tbv_value[tbv_df$id_ind != edited_id]

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
# 5. Verify add_offspring() uses the live genome (n+1 loci).
# ─────────────────────────────────────────────────────────────────────────────

n_loci_after <- as.integer(DBI::dbGetQuery(conn,
  "SELECT COUNT(*) AS n FROM genome_meta")$n)

stopifnot("genome_meta should have n+1 rows" = n_loci_after == n_start + 1L)
cat("Loci after insertion:", n_loci_after, "(was", n_start, ")\n")

# A_1 (heterozygous sire, male) x the 10 founder dams (A_11..A_20, female).
matings <- data.frame(
  id_parent_1 = rep(edited_id, 10),
  id_parent_2 = paste0("A_", 11:20),
  sex         = rep("M", 10),
  line_name   = "A",
  gen         = 1L,
  stringsAsFactors = FALSE
)

pop <- add_offspring(pop, matings)

cat("Offspring added. Checking ind_haplotype for the new locus...\n")

offspring_ids <- paste0("A_", 21:30)

new_col_vals <- DBI::dbGetQuery(conn, sprintf(
  "SELECT id_ind, parent_origin, allele FROM ind_haplotype
   WHERE locus_name = '%s' AND id_ind IN (%s)
   ORDER BY id_ind, parent_origin",
  new_locus,
  paste(sprintf("'%s'", offspring_ids), collapse = ", ")))

cat("\nNew-locus haplotype values in offspring (each has parent_origin 1 and 2):\n")
print(new_col_vals)

stopifnot(
  "Every offspring has haplotype rows at the new locus (2 per individual)" =
    nrow(new_col_vals) == length(offspring_ids) * 2L,
  "All offspring have a row for the new locus" =
    setequal(unique(new_col_vals$id_ind), offspring_ids)
)

# Mendelian check: the sire A_1 is heterozygous (mutant on parent_origin = 1).
# Offspring inherit their parent_origin = 1 haplotype from the sire's gamete, so
# some of those copies must carry the mutant allele (1).
sire_contributions <- new_col_vals$allele[new_col_vals$parent_origin == 1]
cat("\nSire-contributed alleles at new locus:", sire_contributions, "\n")
n_mutant <- sum(sire_contributions == 1L)
cat("Offspring inheriting the mutant allele from the sire:", n_mutant, "of",
    length(sire_contributions), "\n")

stopifnot(
  "Some offspring must inherit the mutant allele from the heterozygous sire" =
    n_mutant >= 1L,
  "Sire-contributed alleles are all 0/1" =
    all(sire_contributions %in% c(0L, 1L))
)

# The offspring genome carries the full n+1 loci in ind_haplotype.
n_distinct_loci <- as.integer(DBI::dbGetQuery(conn, sprintf(
  "SELECT COUNT(DISTINCT locus_id) AS n FROM ind_haplotype WHERE id_ind = '%s'",
  offspring_ids[1]))$n)
stopifnot(
  "offspring haplotype must span all n+1 loci" =
    n_distinct_loci == n_loci_after)
cat("\nPASS: offspring haplotypes span", n_distinct_loci,
    "loci (matches genome_meta row count).\n")

# ─────────────────────────────────────────────────────────────────────────────
# Cleanup
# ─────────────────────────────────────────────────────────────────────────────

close_pop(pop)
cat("\n=== Dynamic genome test PASSED ===\n")
