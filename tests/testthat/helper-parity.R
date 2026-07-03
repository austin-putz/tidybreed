# ---------------------------------------------------------------------------
# Parity harness for the wide -> long haplotype/genotype refactor (Stage 1).
#
# The goal: prove the long-format rewrite preserves current simulation behavior
# exactly. A deterministic, seeded 2-line -> F1 -> F2 simulation is run and its
# semantic content (haplotypes, dosages, TBVs, exported genotype matrix) is
# canonicalized into format-agnostic long tibbles. On the first run (against the
# current WIDE code) these are captured as golden files; every later run (against
# the LONG code) is compared against them. See plans/refactor_haplotype.md.
#
# The canonicalize_*() helpers accept EITHER the wide tables
# (genome_haplotype / genome_genotype with locus_1..locus_n columns) OR the long
# tables (ind_haplotype / ind_genotype), so the same harness works before and
# after the refactor.
# ---------------------------------------------------------------------------

parity_golden_dir <- function() {
  testthat::test_path("parity_golden")
}

# Melt a wide haplotype/genotype data frame (id cols + locus_<id> columns) into
# a long data frame keyed by locus_id. `id_cols` are the non-locus columns.
.parity_melt_wide <- function(wide, id_cols) {
  locus_cols <- grep("^locus_[0-9]+$", names(wide), value = TRUE)
  locus_ids  <- as.integer(sub("^locus_", "", locus_cols))
  mat <- as.matrix(wide[, locus_cols, drop = FALSE])
  n   <- nrow(wide)
  p   <- length(locus_cols)
  out <- data.frame(
    lapply(id_cols, function(cc) rep(wide[[cc]], times = p)),
    locus_id = rep(locus_ids, each = n),
    allele   = as.integer(mat),   # column-major flatten matches rep(each = n)
    stringsAsFactors = FALSE
  )
  names(out)[seq_along(id_cols)] <- id_cols
  out
}

# Canonical haplotypes: (id_ind, parent_origin, locus_name, allele), sorted.
# Reads ind_haplotype if present (long), else genome_haplotype (wide).
# During the transition both wide (genome_haplotype) and long (ind_haplotype)
# tables can exist. Prefer the long table only when it is FULLY populated (a row
# for every individual in ind_meta) so that partial dual-write states (e.g.
# founders migrated but offspring not yet) correctly fall back to the still-
# authoritative wide table. Once the wide table is dropped, use long
# unconditionally.
.parity_use_long <- function(conn) {
  if (!DBI::dbExistsTable(conn, "ind_haplotype")) return(FALSE)
  if (!DBI::dbExistsTable(conn, "genome_haplotype")) return(TRUE)
  n_ind  <- DBI::dbGetQuery(conn, "SELECT COUNT(*) AS n FROM ind_meta")$n
  n_long <- DBI::dbGetQuery(
    conn, "SELECT COUNT(DISTINCT id_ind) AS n FROM ind_haplotype")$n
  n_long > 0 && n_long == n_ind
}

canonicalize_haplotypes <- function(pop) {
  conn <- pop$db_conn
  gm   <- DBI::dbGetQuery(conn, "SELECT locus_id, locus_name FROM genome_meta")

  if (.parity_use_long(conn)) {
    df <- DBI::dbGetQuery(
      conn,
      "SELECT id_ind, parent_origin, locus_name, allele FROM ind_haplotype"
    )
  } else {
    wide <- DBI::dbGetQuery(conn, "SELECT * FROM genome_haplotype")
    long <- .parity_melt_wide(wide, id_cols = c("id_ind", "parent_origin"))
    long <- merge(long, gm, by = "locus_id", all.x = TRUE)
    df   <- long[, c("id_ind", "parent_origin", "locus_name", "allele")]
  }

  df$allele        <- as.integer(df$allele)
  df$parent_origin <- as.integer(df$parent_origin)
  df <- df[order(df$id_ind, df$parent_origin, df$locus_name), , drop = FALSE]
  rownames(df) <- NULL
  tibble::as_tibble(df)
}

# Canonical dosage: (id_ind, locus_name, dosage_value), sorted. ALWAYS derived
# from haplotypes (dosage = sum of alleles across parent_origin). This equals the
# wide genome_genotype (which is that same sum) and the future add_dosage()
# output, so it stays available even after ind_genotype becomes on-demand/empty.
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

# Canonicalize an extract_genotypes() result. The wide code returns columns named
# locus_<locus_id>; the long code returns columns named by locus_name. Map both
# to locus_name and return BOTH a sorted long value tibble and the ordered
# locus_name vector (to verify column ordering survives the PIVOT).
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
