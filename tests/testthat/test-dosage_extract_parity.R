# Stage 3 exit-criteria tests spanning add_dosage() / extract_genotypes() /
# add_tbv() / add_phenotype(): the ind_genotype cache must never diverge from
# a direct extraction, and its (partial) population must never leak into
# TBV/phenotype results.

# Build two identically-seeded pops so a divergence between them can only be
# caused by whichever step differs between the two call sites (mirrors
# helper-parity.R's run_parity_sim() per-step re-seeding discipline).
make_identical_pop <- function(pop_name) {
  pop <- open_pop(pop_name = pop_name, db_name = ":memory:") |>
    define_genome(n_loci = 40, n_chr = 2, chr_len_Mb = 100)

  set.seed(501)
  pop <- define_founder_haplotypes(pop, n_haplotypes = 20, method = "uniform")

  set.seed(502)
  pop <- pop |>
    get_table("founder_haplotypes") |>
    add_founders(n_males = 6, n_females = 6, line_name = "A")

  pop <- define_trait(pop, "ADG", target_add_var = 1.0)

  set.seed(503)
  pop <- pop |>
    get_table("genome_meta") |>
    dplyr::filter(chr == 1L) |>
    define_additive_effects("ADG", distribution = "normal")

  pop <- define_phenotype(pop, "ADG", type = "continuous", mean = 100,
                          residual_var = 0.5)

  pop
}

test_that("add_dosage() cache matches extract_genotypes() direct computation", {
  pop <- make_identical_pop("dose_extract_parity")
  on.exit(close_pop(pop), add = TRUE)

  pop <- pop |>
    get_table("genome_meta") |>
    dplyr::filter(chr == 1L) |>
    define_chip("chip1")
  pop <- pop |> get_table("ind_meta") |> add_genotypes("chip1")

  # Cache path
  pop <- pop |> get_table("ind_meta") |> add_dosage(chip_name = "chip1")
  cache_vals <- DBI::dbGetQuery(pop$db_conn,
    "SELECT id_ind, locus_name, CAST(dosage_value AS INTEGER) AS dosage_value
     FROM ind_genotype ORDER BY id_ind, locus_name")

  # Direct path
  export_df <- pop |> get_table("ind_meta") |> extract_genotypes("chip1")
  direct_vals <- canonicalize_export(pop, export_df)$values
  names(direct_vals)[names(direct_vals) == "dosage_value"] <- "dosage_value"
  direct_vals <- as.data.frame(direct_vals)
  direct_vals <- direct_vals[order(direct_vals$id_ind, direct_vals$locus_name), ]
  rownames(direct_vals) <- NULL

  expect_equal(cache_vals, as.data.frame(direct_vals))
})

test_that("populating ind_genotype for a subset does not change add_tbv()/add_phenotype() results", {
  pop_a <- make_identical_pop("cache_a")
  pop_b <- make_identical_pop("cache_b")
  on.exit({close_pop(pop_a); close_pop(pop_b)}, add = TRUE)

  # pop_a: partially populate ind_genotype (males only, one locus subset).
  pop_a <- pop_a |>
    get_table("ind_meta") |>
    dplyr::filter(sex == "M") |>
    add_dosage(locus_names = c("Locus_1", "Locus_2"))
  expect_true(DBI::dbGetQuery(pop_a$db_conn,
    "SELECT COUNT(*) AS n FROM ind_genotype")$n > 0)

  # pop_b: ind_genotype stays empty.
  expect_equal(DBI::dbGetQuery(pop_b$db_conn,
    "SELECT COUNT(*) AS n FROM ind_genotype")$n, 0L)

  set.seed(601)
  pop_a <- pop_a |> get_table("ind_meta") |> add_tbv("ADG")
  set.seed(601)
  pop_b <- pop_b |> get_table("ind_meta") |> add_tbv("ADG")

  tbv_a <- DBI::dbGetQuery(pop_a$db_conn,
    "SELECT id_ind, trait_name, tbv_value FROM ind_tbv ORDER BY id_ind, trait_name")
  tbv_b <- DBI::dbGetQuery(pop_b$db_conn,
    "SELECT id_ind, trait_name, tbv_value FROM ind_tbv ORDER BY id_ind, trait_name")
  expect_equal(tbv_a, tbv_b)

  set.seed(602)
  pop_a <- pop_a |> get_table("ind_meta") |> add_phenotype("ADG")
  set.seed(602)
  pop_b <- pop_b |> get_table("ind_meta") |> add_phenotype("ADG")

  pheno_a <- DBI::dbGetQuery(pop_a$db_conn,
    "SELECT id_ind, phenotype_name, pheno_value FROM ind_phenotype ORDER BY id_ind, phenotype_name")
  pheno_b <- DBI::dbGetQuery(pop_b$db_conn,
    "SELECT id_ind, phenotype_name, pheno_value FROM ind_phenotype ORDER BY id_ind, phenotype_name")
  expect_equal(pheno_a, pheno_b)
})
