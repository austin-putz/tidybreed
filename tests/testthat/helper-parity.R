# ---------------------------------------------------------------------------
# Parity harness for the long-format haplotype/genotype schema.
#
# The goal: prove the simulation is reproducible for a fixed seed. A
# deterministic, seeded 2-line -> F1 -> F2 simulation is run and its semantic
# content (haplotypes, dosages, TBVs, exported genotype matrix) is canonicalized
# into long tibbles. These are captured as golden files on first run; every later
# run is compared against them. See plans/refactor_haplotype.md.
# ---------------------------------------------------------------------------

parity_golden_dir <- function() {
  testthat::test_path("parity_golden")
}

# Canonical haplotypes: (id_ind, parent_origin, locus_name, allele), sorted,
# read from the long ind_haplotype table.
canonicalize_haplotypes <- function(pop) {
  conn <- pop$db_conn

  df <- DBI::dbGetQuery(
    conn,
    "SELECT id_ind, parent_origin, locus_name, allele FROM ind_haplotype"
  )

  df$allele        <- as.integer(df$allele)
  df$parent_origin <- as.integer(df$parent_origin)
  df <- df[order(df$id_ind, df$parent_origin, df$locus_name), , drop = FALSE]
  rownames(df) <- NULL
  tibble::as_tibble(df)
}

# Canonical dosage: (id_ind, locus_name, dosage_value), sorted. Derived from
# haplotypes (dosage = sum of alleles across parent_origin), which matches
# add_dosage() output, so it stays available even though ind_genotype is
# on-demand/empty until add_dosage() is called.
canonicalize_dosage <- function(pop) {
  hap <- canonicalize_haplotypes(pop)
  agg <- stats::aggregate(
    allele ~ id_ind + locus_name, data = hap, FUN = sum
  )
  names(agg)[names(agg) == "allele"] <- "dosage_value"
  agg$dosage_value <- as.integer(agg$dosage_value)
  agg <- agg[order(agg$id_ind, agg$locus_name), , drop = FALSE]
  rownames(agg) <- NULL
  tibble::as_tibble(agg)
}

# Canonicalize an extract_genotypes() result. extract_genotypes() returns a wide
# matrix with columns named locus_<locus_id> (ordered by locus_id). Map those
# headers to locus_name and return BOTH a sorted long value tibble and the ordered
# locus_name vector (to verify column ordering is preserved).
canonicalize_export <- function(pop, export_df) {
  conn <- pop$db_conn
  gm   <- DBI::dbGetQuery(conn, "SELECT locus_id, locus_name FROM genome_meta")

  geno_cols <- setdiff(names(export_df), "id_ind")
  # Map each genotype column header to a locus_name.
  col_locus_name <- vapply(geno_cols, function(cc) {
    if (grepl("^locus_[0-9]+$", cc)) {
      lid <- as.integer(sub("^locus_", "", cc))
      gm$locus_name[match(lid, gm$locus_id)]
    } else {
      cc  # already a locus_name (long-format output)
    }
  }, character(1))

  mat <- as.matrix(export_df[, geno_cols, drop = FALSE])
  n   <- nrow(export_df)
  p   <- length(geno_cols)
  long <- data.frame(
    id_ind       = rep(export_df$id_ind, times = p),
    locus_name   = rep(col_locus_name, each = n),
    dosage_value = as.integer(mat),
    stringsAsFactors = FALSE
  )
  long <- long[order(long$id_ind, long$locus_name), , drop = FALSE]
  rownames(long) <- NULL

  list(
    values     = tibble::as_tibble(long),
    col_order  = unname(col_locus_name)   # ordering semantics (by locus_id)
  )
}

# Run the deterministic parity simulation and return a list of canonical
# artifacts. Every stochastic step is seeded independently so a failure
# localizes to one artifact and each artifact is reproducible on its own.
run_parity_sim <- function() {
  pop <- open_pop(pop_name = "parity", db_name = ":memory:") |>
    define_genome(n_loci = 60, n_chr = 3, chr_len_Mb = 100)

  set.seed(101)
  pop <- define_founder_haplotypes(pop, n_haplotypes = 30, line_name = "A",
                                   method = "uniform")
  set.seed(102)
  pop <- define_founder_haplotypes(pop, n_haplotypes = 30, line_name = "B",
                                   method = "uniform")

  # Founders: males first (A_1..A_4 / B_1..B_4), then females (A_5..A_8 / B_5..B_8)
  set.seed(201)
  pop <- pop |>
    get_table("founder_haplotypes") |>
    dplyr::filter(line_name == "A") |>
    add_founders(n_males = 4, n_females = 4, line_name = "A", gen = 0L)
  set.seed(202)
  pop <- pop |>
    get_table("founder_haplotypes") |>
    dplyr::filter(line_name == "B") |>
    add_founders(n_males = 4, n_females = 4, line_name = "B", gen = 0L)

  # One non-imprinted trait and one imprinted (parent_1) trait.
  pop <- define_trait(pop, "ADG", target_add_var = 1.0)
  pop <- define_trait(pop, "IMP", target_add_var = 1.0,
                      expressed_parent = "parent_1")

  set.seed(301)
  pop <- pop |>
    get_table("genome_meta") |>
    define_additive_effects("ADG", distribution = "normal")
  set.seed(302)
  pop <- pop |>
    get_table("genome_meta") |>
    define_additive_effects("IMP", distribution = "normal")

  # F1: A sires x B dams (crossbred).
  matings_f1 <- tibble::tibble(
    id_parent_1 = c("A_1", "A_2", "A_3", "A_4"),
    id_parent_2 = c("B_5", "B_6", "B_7", "B_8"),
    sex         = c("M", "F", "M", "F"),
    line_name   = "F1",
    gen         = 1L
  )
  set.seed(401)
  pop <- add_offspring(pop, matings_f1)

  # F2: F1 x F1.
  matings_f2 <- tibble::tibble(
    id_parent_1 = c("F1_1", "F1_3"),
    id_parent_2 = c("F1_2", "F1_4"),
    sex         = c("M", "F"),
    line_name   = "F2",
    gen         = 2L
  )
  set.seed(402)
  pop <- add_offspring(pop, matings_f2)

  # TBVs for everyone.
  pop <- pop |>
    get_table("ind_meta") |>
    add_tbv(c("ADG", "IMP"))

  # Exported genotype matrix over the ADG QTL set.
  export_df <- pop |>
    get_table("ind_meta") |>
    extract_genotypes(
      effects_tbl = get_table(pop, "genome_effects") |>
        dplyr::filter(trait_name == "ADG")
    )

  tbv <- DBI::dbGetQuery(
    pop$db_conn,
    "SELECT id_ind, trait_name, tbv_value FROM ind_tbv"
  )
  tbv <- tbv[order(tbv$id_ind, tbv$trait_name), , drop = FALSE]
  rownames(tbv) <- NULL

  export_canon <- canonicalize_export(pop, export_df)

  artifacts <- list(
    haplotypes   = canonicalize_haplotypes(pop),
    dosage       = canonicalize_dosage(pop),
    tbv          = tibble::as_tibble(tbv),
    export_vals  = export_canon$values,
    export_order = export_canon$col_order
  )

  close_pop(pop)
  artifacts
}
