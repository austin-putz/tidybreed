# Manual test: add_offspring()
# Run interactively to verify the full round-trip against the current
# long-format schema (ind_haplotype / ind_genotype).
#
#   Rscript -e 'suppressMessages(pkgload::load_all(".", quiet=TRUE)); source("tests/test_add_offspring_manual.R")'

library(dplyr)

cat("=== add_offspring() manual test ===\n\n")

# ── Setup ────────────────────────────────────────────────────────────────────
pop <- open_pop(pop_name = "manual_test", db_name = ":memory:") |>
  define_genome(n_loci = 200, n_chr = 3, chr_len_Mb = 100) |>
  define_founder_haplotypes(n_haplotypes = 50)

pop <- pop |>
  get_table("founder_haplotypes") |>
  add_founders(n_males = 5, n_females = 5, line_name = "A", gen = 1L)

cat("Founders added. ind_meta rows:", nrow(collect(get_table(pop, "ind_meta"))), "\n")

# ── Basic offspring ──────────────────────────────────────────────────────────
matings <- tibble::tibble(
  id_parent_1 = rep("A_1", 5),
  id_parent_2 = paste0("A_", 6:10),
  sex         = c("M", "F", "M", "F", "M"),
  line_name   = "A",
  gen         = 2L
)

pop <- pop |> add_offspring(matings)

ind_meta <- collect(get_table(pop, "ind_meta"))
cat("\nAfter add_offspring(5 offspring):\n")
cat("  ind_meta rows:         ", nrow(ind_meta), "\n")
cat("  gen column present:    ", "gen" %in% colnames(ind_meta), "\n")
cat("  offspring gen values:  ", unique(ind_meta$gen[!is.na(ind_meta$id_parent_1)]), "\n")

stopifnot(
  "10 founders + 5 offspring = 15 individuals" = nrow(ind_meta) == 15L,
  "offspring gen recorded"                     = all(ind_meta$gen[!is.na(ind_meta$id_parent_1)] == 2L)
)

# ── Haplotype / genotype checks (long format) ────────────────────────────────
# ind_haplotype: 2 rows per (individual x locus) for diploid autosomes.
# ind_genotype: on-demand cache, empty until add_dosage() is called.
haps <- collect(get_table(pop, "ind_haplotype"))
cat("\nind_haplotype rows:     ", nrow(haps),
    " (expected", nrow(ind_meta) * 200 * 2, ")\n")

stopifnot(
  "ind_haplotype has 2 rows per (ind x locus)" =
    nrow(haps) == nrow(ind_meta) * 200 * 2,
  "all haplotype alleles in {0,1}" = all(haps$allele %in% c(0L, 1L))
)

# ind_genotype starts empty; populate dosages on demand.
genos_before <- collect(get_table(pop, "ind_genotype"))
cat("ind_genotype rows before add_dosage:", nrow(genos_before), "(expected 0)\n")
stopifnot("ind_genotype is empty before add_dosage()" = nrow(genos_before) == 0L)

pop   <- get_table(pop, "ind_meta") |> add_dosage()
genos <- collect(get_table(pop, "ind_genotype"))
cat("ind_genotype rows after add_dosage:  ", nrow(genos),
    " (expected", nrow(ind_meta) * 200, ")\n")
stopifnot("ind_genotype has 1 row per (ind x locus)" =
            nrow(genos) == nrow(ind_meta) * 200)

# Check dosage_value == sum(allele over the ind's 2 strands) for one offspring.
off_id   <- "A_11"
off_haps <- filter(haps,  id_ind == off_id)
off_geno <- filter(genos, id_ind == off_id)

hap_sum <- off_haps |>
  group_by(locus_id) |>
  summarise(dsum = sum(allele), .groups = "drop") |>
  arrange(locus_id)
geno_ord <- off_geno |> arrange(locus_id)

cat("dosage_value == sum(allele) for", off_id, ":",
    all(geno_ord$dosage_value == hap_sum$dsum), "\n")
stopifnot("dosage equals summed alleles" =
            all(geno_ord$dosage_value == hap_sum$dsum))

# ── Alias columns (id_sire / id_dam) ─────────────────────────────────────────
matings_alias <- tibble::tibble(
  id_sire   = "A_1",
  id_dam    = "A_6",
  sex       = "M",
  line_name = "A",
  gen       = 2L
)
pop <- pop |> add_offspring(matings_alias)
cat("\nAlias (id_sire / id_dam) accepted without error: TRUE\n")

# ── Cross-line mating ────────────────────────────────────────────────────────
pop <- pop |>
  define_founder_haplotypes(n_haplotypes = 50, line_name = "B")
pop <- pop |>
  get_table("founder_haplotypes") |>
  filter(line_name == "B") |>
  add_founders(n_males = 2, n_females = 2, line_name = "B", gen = 1L)

cross_matings <- tibble::tibble(
  id_parent_1 = c("A_1", "A_2"),
  id_parent_2 = c("B_1", "B_2"),
  sex         = c("M", "F"),
  line_name   = "C",
  gen         = 2L
)
pop <- pop |> add_offspring(cross_matings)
c_ids <- filter(collect(get_table(pop, "ind_meta")), line_name == "C")$id_ind
cat("\nCross-line offspring (line C):", paste(c_ids, collapse = ", "), "\n")
stopifnot("two line-C offspring created" = length(c_ids) == 2L)

# Cross-line offspring carry line_origin from BOTH parent lines in ind_haplotype.
first_c <- c_ids[1]
c_lo <- collect(get_table(pop, "ind_haplotype") |>
                  filter(id_ind == first_c)) |>
  pull(line_origin) |> unique() |> sort()
cat("line_origin values in first line-C offspring:", paste(c_lo, collapse = ", "), "\n")
stopifnot("line-C offspring inherits both parent lines" =
            identical(c_lo, c("A", "B")))

# ── Error checks ─────────────────────────────────────────────────────────────
cat("\nError checks:\n")
tryCatch(
  add_offspring(pop, tibble::tibble(id_parent_1 = "Z_999", id_parent_2 = "A_6",
                                    sex = "M", line_name = "A")),
  error = function(e) cat("  Non-existent parent correctly caught:", conditionMessage(e), "\n")
)
tryCatch(
  add_offspring(pop, tibble::tibble(id_parent_1 = "A_1", id_parent_2 = "A_6",
                                    sex = "X", line_name = "A")),
  error = function(e) cat("  Invalid sex correctly caught:", conditionMessage(e), "\n")
)

close_pop(pop)
cat("\n=== All manual checks passed ===\n")
