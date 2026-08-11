make_effects_pop <- function(pop_name = "eff", n_ind = 500, n_loci = 500) {
  pop <- open_pop(pop_name = pop_name, db_name = ":memory:") |>
    define_genome(n_loci = n_loci, n_chr = 5, chr_len_Mb = 100) |>
    define_founder_haplotypes(n_haplotypes = 200, method = "fixed")
  pop <- pop |>
    get_table("founder_haplotypes") |>
    add_founders(n_males = n_ind / 2, n_females = n_ind / 2,
                 line_name = "A")
  pop
}

# Dosage matrix (individuals x loci, locus_id order) from the long ind_haplotype.
.dosage_matrix <- function(pop) {
  agg <- DBI::dbGetQuery(pop$db_conn,
    "SELECT id_ind, locus_id, CAST(SUM(allele) AS INTEGER) AS d FROM ind_haplotype GROUP BY id_ind, locus_id")
  n_loci <- DBI::dbGetQuery(pop$db_conn, "SELECT COUNT(*) AS n FROM genome_meta")$n
  ids <- unique(agg$id_ind)
  m <- matrix(0L, nrow = length(ids), ncol = n_loci, dimnames = list(ids, NULL))
  m[cbind(match(agg$id_ind, ids), agg$locus_id)] <- as.integer(agg$d)
  m
}


test_that("define_additive_effects() rescales to target_add_var within tolerance", {
  set.seed(42)
  pop <- make_effects_pop("eff_scale", n_ind = 600, n_loci = 600)

  pop <- define_trait(pop, "ADG", target_add_var = 0.5)
  sel <- pop |> get_table("genome_meta") |> dplyr::collect() |>
    dplyr::slice_sample(n = 100) |> dplyr::pull(locus_name)
  pop <- pop |>
    get_table("genome_meta") |>
    dplyr::filter(locus_name %in% sel) |>
    define_additive_effects("ADG", distribution = "normal", seed = 1)

  # Effects are now in genome_effects, not genome_meta columns
  eff <- DBI::dbGetQuery(pop$db_conn,
    "SELECT locus_name, genome_value FROM genome_effects WHERE trait_name = 'ADG'")
  locus_order <- DBI::dbGetQuery(pop$db_conn,
    "SELECT locus_id, locus_name FROM genome_meta ORDER BY locus_id")
  a <- rep(0, nrow(locus_order))
  idx <- match(eff$locus_name, locus_order$locus_name)
  a[idx] <- eff$genome_value

  # The rescaler's actual contract: the Falconer expected additive variance
  # under the base allele frequencies, sum(2 p q a^2), equals target_add_var.
  # This is deterministic given the effects, so it is asserted tightly.
  p <- DBI::dbGetQuery(pop$db_conn,
    "SELECT base_allele_freq AS p, genome_value AS a
       FROM genome_effects WHERE trait_name = 'ADG'")
  expect_equal(sum(2 * p$p * (1 - p$p) * p$a^2), 0.5, tolerance = 1e-8)

  # The variance *realised* in the sampled founders is a noisy estimate of that
  # expectation: 600 individuals drawn from a 200-haplotype pool carry drift and
  # LD, so per-locus realised 2pq deviates from the base p. Measured spread
  # across seeds is roughly +/-25% of target, so this bound is a sanity check on
  # the order of magnitude, not a precision test -- do not tighten it to chase a
  # lucky seed.
  X <- .dosage_matrix(pop)
  realised <- var(as.numeric(X %*% a))

  expect_equal(realised, 0.5, tolerance = 0.35)
  close_pop(pop)
})


test_that("TBV mean is approximately 0 for founder population", {
  set.seed(7)
  pop <- make_effects_pop("eff_mean", n_ind = 500, n_loci = 500)

  pop <- define_trait(pop, "ADG", target_add_var = 100, target_add_mean = 0)
  sel <- pop |> get_table("genome_meta") |> dplyr::collect() |>
    dplyr::slice_sample(n = 200) |> dplyr::pull(locus_name)
  pop <- pop |>
    get_table("genome_meta") |>
    dplyr::filter(locus_name %in% sel) |>
    define_additive_effects("ADG", distribution = "normal",
                          base = "founder_haplotypes", seed = 3)

  pop <- pop |> get_table("ind_meta") |> add_tbv("ADG")

  tbv_df <- dplyr::collect(get_table(pop, "ind_tbv"))
  expect_equal(nrow(tbv_df), 500)
  # mean TBV should be close to 0 (within ~2 SE = 2*sqrt(100/500) ≈ 0.9)
  expect_equal(mean(tbv_df$tbv_value), 0, tolerance = 2.0)
  # var TBV should be close to target_add_var = 100
  expect_equal(var(tbv_df$tbv_value), 100, tolerance = 15)

  close_pop(pop)
})


test_that("base_allele_freq written to genome_effects, not genome_meta", {
  pop <- make_effects_pop("eff_base_col")

  pop <- define_trait(pop, "ADG", target_add_var = 1)
  sel <- pop |> get_table("genome_meta") |> dplyr::collect() |>
    dplyr::slice_sample(n = 50) |> dplyr::pull(locus_name)
  pop <- pop |>
    get_table("genome_meta") |>
    dplyr::filter(locus_name %in% sel) |>
    define_additive_effects("ADG", distribution = "normal")

  # base_allele_freq is in genome_effects, not genome_meta
  genome_cols <- DBI::dbListFields(pop$db_conn, "genome_meta")
  expect_false("base_allele_freq_ADG" %in% genome_cols)
  expect_false("add_ADG"              %in% genome_cols)
  expect_false("is_QTL_ADG"          %in% genome_cols)

  eff <- DBI::dbGetQuery(pop$db_conn,
    "SELECT base_allele_freq FROM genome_effects WHERE trait_name = 'ADG'")
  expect_equal(nrow(eff), 50)
  expect_true(all(eff$base_allele_freq >= 0 & eff$base_allele_freq <= 1))

  close_pop(pop)
})


test_that("base = 'current_pop' via base_tbl argument works", {
  set.seed(17)
  pop <- make_effects_pop("eff_currpop", n_ind = 200, n_loci = 300)

  pop <- get_table(pop, "ind_meta") |> mutate_table(gen = 0L)
  pop <- define_trait(pop, "ADG", target_add_var = 50)
  sel <- pop |> get_table("genome_meta") |> dplyr::collect() |>
    dplyr::slice_sample(n = 100) |> dplyr::pull(locus_name)

  gen0_tbl <- get_table(pop, "ind_meta") |> dplyr::filter(gen == 0L)
  pop <- pop |>
    get_table("genome_meta") |>
    dplyr::filter(locus_name %in% sel) |>
    define_additive_effects("ADG", base = "current_pop", base_tbl = gen0_tbl,
                          distribution = "normal", seed = 5)

  eff <- DBI::dbGetQuery(pop$db_conn,
    "SELECT locus_name, genome_value, base_allele_freq FROM genome_effects WHERE trait_name = 'ADG'")
  expect_equal(nrow(eff), 100)

  # TBV mean should be ≈ 0
  pop <- pop |> get_table("ind_meta") |> add_tbv("ADG")
  tbv_df <- dplyr::collect(get_table(pop, "ind_tbv"))
  expect_equal(mean(tbv_df$tbv_value), 0, tolerance = 3.0)

  close_pop(pop)
})


test_that("define_additive_effects() accepts manual effects", {
  pop <- make_effects_pop("eff_manual")
  pop <- define_trait(pop, "ADG", target_add_var = 1)
  sel <- pop |> get_table("genome_meta") |> dplyr::collect() |>
    dplyr::slice_sample(n = 10) |> dplyr::pull(locus_name)
  pop <- pop |>
    get_table("genome_meta") |>
    dplyr::filter(locus_name %in% sel) |>
    define_additive_effects("ADG", effects = rep(2.0, 10))

  eff <- DBI::dbGetQuery(pop$db_conn,
    "SELECT genome_value FROM genome_effects WHERE trait_name = 'ADG'")
  expect_equal(nrow(eff), 10)
  expect_true(all(eff$genome_value == 2.0))

  close_pop(pop)
})


test_that("re-calling define_additive_effects() replaces existing rows", {
  pop <- make_effects_pop("eff_replace")
  pop <- define_trait(pop, "ADG", target_add_var = 1)
  sel <- pop |> get_table("genome_meta") |> dplyr::collect() |>
    dplyr::slice_sample(n = 20) |> dplyr::pull(locus_name)

  pop <- pop |>
    get_table("genome_meta") |>
    dplyr::filter(locus_name %in% sel) |>
    define_additive_effects("ADG", effects = rep(1.0, 20))

  n_before <- DBI::dbGetQuery(pop$db_conn,
    "SELECT COUNT(*) AS n FROM genome_effects WHERE trait_name = 'ADG'")$n
  expect_equal(n_before, 20L)

  # Call again with different loci set
  sel2 <- pop |> get_table("genome_meta") |> dplyr::collect() |>
    dplyr::slice_sample(n = 30) |> dplyr::pull(locus_name)
  pop <- pop |>
    get_table("genome_meta") |>
    dplyr::filter(locus_name %in% sel2) |>
    define_additive_effects("ADG", effects = rep(3.0, 30))

  n_after <- DBI::dbGetQuery(pop$db_conn,
    "SELECT COUNT(*) AS n FROM genome_effects WHERE trait_name = 'ADG'")$n
  expect_equal(n_after, 30L)
  eff_vals <- DBI::dbGetQuery(pop$db_conn,
    "SELECT genome_value FROM genome_effects WHERE trait_name = 'ADG'")$genome_value
  expect_true(all(eff_vals == 3.0))

  close_pop(pop)
})


test_that("define_additive_effects() hits target variances per trait (multi-trait)", {
  set.seed(123)
  pop <- make_effects_pop("eff_multi", n_ind = 800, n_loci = 600)

  pop <- define_trait(pop, "ADG", target_add_var = 0.25)
  pop <- define_trait(pop, "BW",  target_add_var = 0.50)

  # Same QTL for both traits (full pleiotropy via method = "shared")
  sel <- pop |> get_table("genome_meta") |> dplyr::collect() |>
    dplyr::slice_sample(n = 150) |> dplyr::pull(locus_name)

  G <- matrix(c(0.25, 0.10, 0.10, 0.50), 2, 2)
  pop <- pop |>
    get_table("genome_meta") |>
    dplyr::filter(locus_name %in% sel) |>
    define_additive_effects(trait_name = c("ADG", "BW"), G = G,
                             method = "shared", seed = 7)

  locus_order <- DBI::dbGetQuery(pop$db_conn,
    "SELECT locus_id, locus_name FROM genome_meta ORDER BY locus_id")
  n_loci <- nrow(locus_order)

  load_eff <- function(t) {
    e <- DBI::dbGetQuery(pop$db_conn, paste0(
      "SELECT locus_name, genome_value FROM genome_effects WHERE trait_name = '", t, "'"))
    a <- rep(0, n_loci)
    idx <- match(e$locus_name, locus_order$locus_name)
    a[idx] <- e$genome_value
    a
  }
  aA <- load_eff("ADG")
  aB <- load_eff("BW")

  X <- .dosage_matrix(pop)

  bv_A <- as.numeric(X %*% aA)
  bv_B <- as.numeric(X %*% aB)

  expect_equal(var(bv_A), 0.25, tolerance = 0.15)
  expect_equal(var(bv_B), 0.50, tolerance = 0.16)
  expect_gt(cor(bv_A, bv_B), 0.05)

  close_pop(pop)
})


test_that("define_additive_effects() errors on bare tidybreed_pop", {
  pop <- make_effects_pop("eff_err_pop")
  pop <- define_trait(pop, "ADG", target_add_var = 1)

  expect_error(
    define_additive_effects(pop, "ADG"),
    "tidybreed_table"
  )
  close_pop(pop)
})


test_that("define_additive_effects() errors when filter returns zero rows", {
  pop <- make_effects_pop("eff_err_empty")
  pop <- define_trait(pop, "ADG", target_add_var = 1)

  expect_error(
    pop |> get_table("genome_meta") |> dplyr::filter(locus_name == "NONEXISTENT") |>
      define_additive_effects("ADG"),
    "zero rows"
  )
  close_pop(pop)
})


# ---------------------------------------------------------------------------
# scale_to_target guard for sex-linked/organelle QTL (Stage 4)
# ---------------------------------------------------------------------------

make_effects_pop_with_x <- function(pop_name = "eff_x", n_ind = 20, n_loci = 20) {
  pop <- open_pop(pop_name = pop_name, db_name = ":memory:") |>
    define_genome(n_loci = n_loci, n_chr = 2, chr_names = c("1", "X"), chr_len_Mb = 100) |>
    define_chromosome("X", offspring_sex = "M", from_parent_1 = 0, from_parent_2 = 1) |>
    define_founder_haplotypes(n_haplotypes = 20, method = "fixed")
  pop |>
    get_table("founder_haplotypes") |>
    add_founders(n_males = n_ind / 2, n_females = n_ind / 2, line_name = "A")
}

test_that("scale_to_target = TRUE errors when QTL set includes a sex-linked locus", {
  pop <- make_effects_pop_with_x()
  on.exit(close_pop(pop), add = TRUE)
  pop <- define_trait(pop, "ADG", target_add_var = 1)

  expect_error(
    pop |> get_table("genome_meta") |> dplyr::filter(chr_name == "X") |>
      define_additive_effects("ADG", scale_to_target = TRUE),
    "Falconer variance scaling"
  )
})

test_that("scale_to_target = FALSE with manual effects works fine for sex-linked QTL", {
  pop <- make_effects_pop_with_x()
  on.exit(close_pop(pop), add = TRUE)
  pop <- define_trait(pop, "ADG", target_add_var = 1)

  n_x <- DBI::dbGetQuery(pop$db_conn,
    "SELECT COUNT(*) AS n FROM genome_meta WHERE chr_name = 'X'")$n

  expect_no_error(
    pop |> get_table("genome_meta") |> dplyr::filter(chr_name == "X") |>
      define_additive_effects("ADG", effects = rep(1, n_x))
  )
})

test_that("scale_to_target = TRUE still works for purely autosomal QTL on a genome that also has a sex chromosome", {
  pop <- make_effects_pop_with_x()
  on.exit(close_pop(pop), add = TRUE)
  pop <- define_trait(pop, "ADG", target_add_var = 1)

  expect_no_error(
    pop |> get_table("genome_meta") |> dplyr::filter(chr_name == "1") |>
      define_additive_effects("ADG", scale_to_target = TRUE)
  )
})


# ============================================================
# base_line_name / per-line Falconer centering
# ============================================================

# Two lines fixed for OPPOSITE alleles at every locus. Within-line 2pq is 0 at
# every locus; pooling them gives p = 0.5 and an apparent 2pq = 0.5. This is the
# sharpest possible statement of the Wahlund effect.
make_two_line_pop <- function(pop_name, n_loci = 40, n_hap = 20) {
  pop <- open_pop(pop_name = pop_name, db_name = ":memory:") |>
    define_genome(n_loci = n_loci, n_chr = 2, chr_len_Mb = 50)
  gm <- DBI::dbGetQuery(pop$db_conn,
    "SELECT locus_name FROM genome_meta ORDER BY locus_id")
  for (ln in c("A", "B")) {
    fh <- data.frame(
      line_name    = ln,
      haplotype_id = rep(seq_len(n_hap), times = nrow(gm)),
      locus_name   = rep(gm$locus_name, each = n_hap),
      allele       = if (ln == "A") 0L else 1L,
      stringsAsFactors = FALSE
    )
    DBI::dbWriteTable(pop$db_conn, "founder_haplotypes", fh, append = TRUE)
  }
  pop$tables <- unique(c(pop$tables, "founder_haplotypes"))
  pop
}

test_that("base_line_name inherits line_name so each line centers on its own pool", {
  pop <- make_two_line_pop("bln_inherit")
  pop <- define_trait(pop, "ADG", target_add_var = 1)

  stored <- function(ln) {
    DBI::dbGetQuery(pop$db_conn, paste0(
      "SELECT DISTINCT base_allele_freq FROM genome_effects ",
      "WHERE trait_name = 'ADG' AND line_name ",
      if (is.null(ln)) "IS NULL" else paste0("= '", ln, "'")))$base_allele_freq
  }

  pop |> get_table("genome_meta") |>
    define_additive_effects("ADG", effects = rep(1, 40), line_name = "A")
  pop |> get_table("genome_meta") |>
    define_additive_effects("ADG", effects = rep(1, 40), line_name = "B")

  # Line A is fixed at allele 0, line B at allele 1 -- each sees its own.
  expect_equal(stored("A"), 0)
  expect_equal(stored("B"), 1)

  # A population-wide effect is still centered on the pooled base, which is the
  # right answer for an effect that applies to the whole founder base.
  expect_warning(
    pop |> get_table("genome_meta") |>
      define_additive_effects("ADG", effects = rep(1, 40)),
    "pooled across"
  )
  expect_equal(stored(NULL), 0.5)

  close_pop(pop)
})


test_that("base_line_name = NULL forces pooling even for a line-specific effect", {
  pop <- make_two_line_pop("bln_forced")
  pop <- define_trait(pop, "ADG", target_add_var = 1)

  suppressWarnings(
    pop |> get_table("genome_meta") |>
      define_additive_effects("ADG", effects = rep(1, 40), line_name = "A",
                              base_line_name = NULL))

  expect_equal(
    DBI::dbGetQuery(pop$db_conn,
      "SELECT DISTINCT base_allele_freq FROM genome_effects
         WHERE trait_name = 'ADG' AND line_name = 'A'")$base_allele_freq,
    0.5
  )

  close_pop(pop)
})


test_that("per-line centering recovers target_add_var that pooling misses", {
  # Line A: allele 0 fixed at half the loci, polymorphic at the rest, so the
  # within-line and pooled frequencies genuinely differ.
  set.seed(404)
  pop <- open_pop(pop_name = "bln_var", db_name = ":memory:") |>
    define_genome(n_loci = 60, n_chr = 2, chr_len_Mb = 50)
  gm <- DBI::dbGetQuery(pop$db_conn,
    "SELECT locus_name FROM genome_meta ORDER BY locus_id")
  n_hap <- 40
  mk_line <- function(ln, p) {
    alle <- as.integer(stats::runif(nrow(gm) * n_hap) < p)
    DBI::dbWriteTable(pop$db_conn, "founder_haplotypes", data.frame(
      line_name    = ln,
      haplotype_id = rep(seq_len(n_hap), times = nrow(gm)),
      locus_name   = rep(gm$locus_name, each = n_hap),
      allele       = alle,
      stringsAsFactors = FALSE
    ), append = TRUE)
  }
  mk_line("A", 0.1)   # rare allele in A
  mk_line("B", 0.9)   # common allele in B -- pooled sits near 0.5
  pop$tables <- unique(c(pop$tables, "founder_haplotypes"))
  pop <- define_trait(pop, "ADG", target_add_var = 2)

  falconer <- function(ln) {
    e <- DBI::dbGetQuery(pop$db_conn, paste0(
      "SELECT base_allele_freq p, genome_value a FROM genome_effects ",
      "WHERE trait_name = 'ADG' AND line_name = '", ln, "'"))
    sum(2 * e$p * (1 - e$p) * e$a^2)
  }
  # Realized within-line variance implied by the stored effects and the line's
  # OWN allele frequencies -- what the simulation actually delivers.
  realised <- function(ln) {
    e <- DBI::dbGetQuery(pop$db_conn, paste0(
      "SELECT e.locus_name, e.genome_value a, f.p FROM genome_effects e ",
      "JOIN (SELECT locus_name, AVG(CAST(allele AS DOUBLE)) p ",
      "        FROM founder_haplotypes WHERE line_name = '", ln, "' ",
      "        GROUP BY locus_name) f ON f.locus_name = e.locus_name ",
      "WHERE e.trait_name = 'ADG' AND e.line_name = '", ln, "'"))
    sum(2 * e$p * (1 - e$p) * e$a^2)
  }

  set.seed(1)
  pop |> get_table("genome_meta") |>
    define_additive_effects("ADG", line_name = "A", seed = 1)
  expect_equal(falconer("A"), 2, tolerance = 1e-8)
  expect_equal(realised("A"), 2, tolerance = 1e-8)

  # Now force the old pooled behaviour: the Falconer bookkeeping still "hits"
  # the target, but the variance actually realised within line A does not.
  set.seed(1)
  suppressWarnings(
    pop |> get_table("genome_meta") |>
      define_additive_effects("ADG", line_name = "A", base_line_name = NULL,
                              seed = 1))
  expect_equal(falconer("A"), 2, tolerance = 1e-8)
  expect_lt(realised("A"), 1.0)   # pooling under-scales: well short of 2

  close_pop(pop)
})


test_that("base_line_name validates its input", {
  pop <- make_two_line_pop("bln_valid")
  pop <- define_trait(pop, "ADG", target_add_var = 1)
  on.exit(close_pop(pop), add = TRUE)

  # A typo must be loud: base frequencies are zero-initialised, so a silent
  # miss would center every allele at 0 and contribute nothing to V_A.
  expect_error(
    pop |> get_table("genome_meta") |>
      define_additive_effects("ADG", effects = rep(1, 40), base_line_name = "NOPE"),
    "No founder_haplotypes rows for line"
  )
  expect_error(
    pop |> get_table("genome_meta") |>
      define_additive_effects("ADG", effects = rep(1, 40),
                              base = "current_pop", base_line_name = "A"),
    "applies only to base"
  )
  expect_error(
    pop |> get_table("genome_meta") |>
      define_additive_effects("ADG", effects = rep(1, 40),
                              base_line_name = "A'; DROP TABLE genome_meta; --"),
    "Invalid line name"
  )
  expect_true(DBI::dbExistsTable(pop$db_conn, "genome_meta"))
})


test_that("defining effects for one line does not clobber another line's rows", {
  pop <- make_two_line_pop("bln_clobber")
  pop <- define_trait(pop, "ADG", target_add_var = 1)

  pop |> get_table("genome_meta") |>
    define_additive_effects("ADG", effects = rep(1, 40), line_name = "A")
  pop |> get_table("genome_meta") |>
    define_additive_effects("ADG", effects = rep(2, 40), line_name = "B")
  expect_warning(
    pop |> get_table("genome_meta") |>
      define_additive_effects("ADG", effects = rep(3, 40)),
    "pooled across"
  )

  counts <- DBI::dbGetQuery(pop$db_conn,
    "SELECT line_name, COUNT(*) AS n, MIN(genome_value) AS v
       FROM genome_effects WHERE trait_name = 'ADG'
       GROUP BY line_name ORDER BY line_name NULLS LAST")

  expect_equal(nrow(counts), 3L)
  expect_true(all(counts$n == 40L))
  expect_equal(counts$v[counts$line_name == "A" & !is.na(counts$line_name)], 1)
  expect_equal(counts$v[counts$line_name == "B" & !is.na(counts$line_name)], 2)
  expect_equal(counts$v[is.na(counts$line_name)], 3)

  close_pop(pop)
})
