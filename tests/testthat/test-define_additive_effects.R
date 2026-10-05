# These tests were written against the flat genome_effects table: one row per
# (locus, line), carrying the coefficient and its centring frequency together.
# Effects now live as terms over members with an origin scope, so the flat shape
# is reconstructed as the test-only `additive_flat` view (ge_flat_view(), in
# helper-genome-effects-db.R). The package deliberately ships no such view --
# the point of the new schema is that a term is not a locus -- but every
# assertion below is about *generated additive* effects, which are exactly the
# order-one, single-origin case the flat shape described correctly.
make_effects_pop <- function(pop_name = "eff", n_ind = 500, n_loci = 500) {
  pop <- open_pop(pop_name = pop_name, db_name = ":memory:") |>
    define_genome(n_loci = n_loci, n_chr = 5, chr_len_Mb = 100) |>
    define_founder_haplotypes(n_haplotypes = 200, method = "fixed")
  pop <- pop |>
    get_table("founder_haplotypes") |>
    add_founders(n_males = n_ind / 2, n_females = n_ind / 2,
                 line_name = "A")
  ge_flat_view(pop)
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


test_that("define_additive_effects() rescales to the stored additive target", {
  set.seed(42)
  pop <- make_effects_pop("eff_scale", n_ind = 600, n_loci = 600)

  pop <- with_additive_target(pop, "ADG", 0.5)
  sel <- pop |> get_table("genome_meta") |> dplyr::collect() |>
    dplyr::slice_sample(n = 100) |> dplyr::pull(locus_name)
  pop <- pop |>
    get_table("genome_meta") |>
    dplyr::filter(locus_name %in% sel) |>
    define_additive_effects("ADG", distribution = "normal", seed = 1)

  # Effects are now in genome_effects, not genome_meta columns
  eff <- DBI::dbGetQuery(pop$db_conn,
    "SELECT locus_name, genome_value FROM additive_flat WHERE trait_name = 'ADG'")
  locus_order <- DBI::dbGetQuery(pop$db_conn,
    "SELECT locus_id, locus_name FROM genome_meta ORDER BY locus_id")
  a <- rep(0, nrow(locus_order))
  idx <- match(eff$locus_name, locus_order$locus_name)
  a[idx] <- eff$genome_value

  # The rescaler's actual contract: the Falconer expected additive variance
  # under the base allele frequencies, sum(2 p q a^2), equals the target.
  # This is deterministic given the effects, so it is asserted tightly.
  p <- DBI::dbGetQuery(pop$db_conn,
    "SELECT center_value AS p, genome_value AS a
       FROM additive_flat WHERE trait_name = 'ADG'")
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

  pop <- with_additive_target(pop, "ADG", 100)
  sel <- pop |> get_table("genome_meta") |> dplyr::collect() |>
    dplyr::slice_sample(n = 200) |> dplyr::pull(locus_name)
  pop <- pop |>
    get_table("genome_meta") |>
    dplyr::filter(locus_name %in% sel) |>
    define_additive_effects("ADG", distribution = "normal", seed = 3)

  pop <- pop |> get_table("ind_meta") |> add_tgv("ADG")

  tbv_df <- tgv_additive(pop)
  expect_equal(nrow(tbv_df), 500)
  # mean TBV should be close to 0 (within ~2 SE = 2*sqrt(100/500) ≈ 0.9)
  expect_equal(mean(tbv_df$tgv_value), 0, tolerance = 2.0)
  # var TBV should be close to the target, 100
  expect_equal(var(tbv_df$tgv_value), 100, tolerance = 15)

  close_pop(pop)
})


test_that("center_value written to genome_effect_members, not genome_meta", {
  pop <- make_effects_pop("eff_base_col")

  pop <- with_additive_target(pop, "ADG", 1)
  sel <- pop |> get_table("genome_meta") |> dplyr::collect() |>
    dplyr::slice_sample(n = 50) |> dplyr::pull(locus_name)
  pop <- pop |>
    get_table("genome_meta") |>
    dplyr::filter(locus_name %in% sel) |>
    define_additive_effects("ADG", distribution = "normal")

  # center_value is in genome_effect_members, not a genome_meta column
  genome_cols <- DBI::dbListFields(pop$db_conn, "genome_meta")
  expect_false("base_allele_freq_ADG" %in% genome_cols)
  expect_false("add_ADG"              %in% genome_cols)
  expect_false("is_QTL_ADG"          %in% genome_cols)

  eff <- DBI::dbGetQuery(pop$db_conn,
    "SELECT center_value FROM additive_flat WHERE trait_name = 'ADG'")
  expect_equal(nrow(eff), 50)
  expect_true(all(eff$center_value >= 0 & eff$center_value <= 1))

  close_pop(pop)
})


test_that("base = 'current_pop' via base_tbl argument works", {
  set.seed(17)
  pop <- make_effects_pop("eff_currpop", n_ind = 200, n_loci = 300)

  pop <- get_table(pop, "ind_meta") |> mutate_table(gen = 0L)
  pop <- with_additive_target(pop, "ADG", 50)
  sel <- pop |> get_table("genome_meta") |> dplyr::collect() |>
    dplyr::slice_sample(n = 100) |> dplyr::pull(locus_name)

  gen0_tbl <- get_table(pop, "ind_meta") |> dplyr::filter(gen == 0L)
  pop <- pop |>
    get_table("genome_meta") |>
    dplyr::filter(locus_name %in% sel) |>
    define_additive_effects("ADG", base_tbl = gen0_tbl,
                          distribution = "normal", seed = 5)

  eff <- DBI::dbGetQuery(pop$db_conn,
    "SELECT locus_name, genome_value, center_value FROM additive_flat WHERE trait_name = 'ADG'")
  expect_equal(nrow(eff), 100)

  # TBV mean should be ≈ 0
  pop <- pop |> get_table("ind_meta") |> add_tgv("ADG")
  tbv_df <- tgv_additive(pop)
  expect_equal(mean(tbv_df$tgv_value), 0, tolerance = 3.0)

  close_pop(pop)
})


test_that("the generator has no manual or unscaled mode (Q21)", {
  # 0.74.1: generated means calibrated. Known coefficients go through
  # define_genome_effect_terms(); there is no `effects =` or
  # `scale_to_target =` to bypass the calibration.
  pop <- make_effects_pop("eff_manual")
  on.exit(close_pop(pop), add = TRUE)
  pop <- with_additive_target(pop, "ADG", 1)
  gm  <- pop |> get_table("genome_meta")
  expect_error(gm |> define_additive_effects("ADG", effects = rep(2.0, 10)),
               "unused argument \\(effects")
  expect_error(gm |> define_additive_effects("ADG", scale_to_target = FALSE),
               "unused argument \\(scale_to_target")
  expect_equal(DBI::dbGetQuery(pop$db_conn,
    "SELECT COUNT(*) AS n FROM genome_effects")$n, 0)
})


test_that("re-calling define_additive_effects() replaces existing rows", {
  pop <- make_effects_pop("eff_replace")
  pop <- with_additive_target(pop, "ADG", 1)
  sel <- pop |> get_table("genome_meta") |> dplyr::collect() |>
    dplyr::slice_sample(n = 20) |> dplyr::pull(locus_name)

  pop <- pop |>
    get_table("genome_meta") |>
    dplyr::filter(locus_name %in% sel) |>
    define_additive_effects("ADG", seed = 1)

  n_before <- DBI::dbGetQuery(pop$db_conn,
    "SELECT COUNT(*) AS n FROM additive_flat WHERE trait_name = 'ADG'")$n
  expect_equal(n_before, 20L)

  # Call again with different loci set
  sel2 <- pop |> get_table("genome_meta") |> dplyr::collect() |>
    dplyr::slice_sample(n = 30) |> dplyr::pull(locus_name)
  pop <- pop |>
    get_table("genome_meta") |>
    dplyr::filter(locus_name %in% sel2) |>
    define_additive_effects("ADG", seed = 2)

  n_after <- DBI::dbGetQuery(pop$db_conn,
    "SELECT COUNT(*) AS n FROM additive_flat WHERE trait_name = 'ADG'")$n
  expect_equal(n_after, 30L)
  # Exactly the second call's loci remain.
  got <- DBI::dbGetQuery(pop$db_conn,
    "SELECT locus_name FROM additive_flat WHERE trait_name = 'ADG'")$locus_name
  expect_setequal(got, sel2)

  close_pop(pop)
})


test_that("define_additive_effects() hits target variances per trait (multi-trait)", {
  set.seed(123)
  pop <- make_effects_pop("eff_multi", n_ind = 800, n_loci = 600)

  pop <- define_trait(pop, "ADG")
  pop <- define_trait(pop, "BW")

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
      "SELECT locus_name, genome_value FROM additive_flat WHERE trait_name = '", t, "'"))
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
  pop <- with_additive_target(pop, "ADG", 1)

  expect_error(
    define_additive_effects(pop, "ADG"),
    "tidybreed_table"
  )
  close_pop(pop)
})


test_that("define_additive_effects() errors when filter returns zero rows", {
  pop <- make_effects_pop("eff_err_empty")
  pop <- with_additive_target(pop, "ADG", 1)

  expect_error(
    pop |> get_table("genome_meta") |> dplyr::filter(locus_name == "NONEXISTENT") |>
      define_additive_effects("ADG"),
    "zero rows"
  )
  close_pop(pop)
})


# ---------------------------------------------------------------------------
# Calibration guard for sex-linked/organelle QTL (Stage 4)
# ---------------------------------------------------------------------------

make_effects_pop_with_x <- function(pop_name = "eff_x", n_ind = 20, n_loci = 20) {
  pop <- open_pop(pop_name = pop_name, db_name = ":memory:") |>
    define_genome(n_loci = n_loci, n_chr = 2, chr_names = c("1", "X"), chr_len_Mb = 100) |>
    define_chromosome("X", offspring_sex = "M", from_parent_1 = 0, from_parent_2 = 1) |>
    define_founder_haplotypes(n_haplotypes = 20, method = "fixed")
  pop <- pop |>
    get_table("founder_haplotypes") |>
    add_founders(n_males = n_ind / 2, n_females = n_ind / 2, line_name = "A")
  ge_flat_view(pop)
}

test_that("the generator errors when the QTL set includes a sex-linked locus", {
  pop <- make_effects_pop_with_x()
  on.exit(close_pop(pop), add = TRUE)
  pop <- with_additive_target(pop, "ADG", 1)

  expect_error(
    pop |> get_table("genome_meta") |> dplyr::filter(chr_name == "X") |>
      define_additive_effects("ADG"),
    "assume diploid/autosomal QTL.*define_genome_effect_terms"
  )
})

test_that("sex-linked QTL with known values go through the writer", {
  pop <- make_effects_pop_with_x()
  on.exit(close_pop(pop), add = TRUE)
  pop <- with_additive_target(pop, "ADG", 1)

  n_x <- DBI::dbGetQuery(pop$db_conn,
    "SELECT COUNT(*) AS n FROM genome_meta WHERE chr_name = 'X'")$n

  expect_no_error(
    pop |> get_table("genome_meta") |> dplyr::filter(chr_name == "X") |>
      with_additive_terms("ADG", effects = rep(1, n_x))
  )
})

test_that("the generator still works for purely autosomal QTL on a genome that also has a sex chromosome", {
  pop <- make_effects_pop_with_x()
  on.exit(close_pop(pop), add = TRUE)
  pop <- with_additive_target(pop, "ADG", 1)

  expect_no_error(
    pop |> get_table("genome_meta") |> dplyr::filter(chr_name == "1") |>
      define_additive_effects("ADG")
  )
})


# ============================================================
# base_tbl / per-line Falconer centering
# (plans/update_genome_effects_base_tbl.md §3.3, §3.8 tests 10-15)
# ============================================================

# Two lines fixed for OPPOSITE alleles at every locus. Within-line 2pq is 0 at
# every locus; pooling them gives p = 0.5 and an apparent 2pq = 0.5. This is the
# sharpest possible statement of the Wahlund effect.
make_two_line_pop <- function(pop_name, n_loci = 40, n_hap = 20,
                              lines = c("A", "B")) {
  pop <- open_pop(pop_name = pop_name, db_name = ":memory:") |>
    define_genome(n_loci = n_loci, n_chr = 2, chr_len_Mb = 50)
  gm <- DBI::dbGetQuery(pop$db_conn,
    "SELECT locus_name FROM genome_meta ORDER BY locus_id")
  for (ln in lines) {
    fh <- data.frame(
      line_name    = if (is.na(ln)) NA_character_ else ln,
      haplotype_id = rep(seq_len(n_hap), times = nrow(gm)),
      locus_name   = rep(gm$locus_name, each = n_hap),
      allele       = if (identical(ln, "A")) 0L else 1L,
      stringsAsFactors = FALSE
    )
    DBI::dbWriteTable(pop$db_conn, "founder_haplotypes", fh, append = TRUE)
  }
  pop$tables <- unique(c(pop$tables, "founder_haplotypes"))
  ge_flat_view(pop)
}

stored_center <- function(pop, ln) {
  DBI::dbGetQuery(pop$db_conn, paste0(
    "SELECT DISTINCT center_value FROM additive_flat ",
    "WHERE trait_name = 'ADG' AND line_name ",
    if (is.null(ln)) "IS NULL" else paste0("= '", ln, "'")))$center_value
}

# Each line of make_two_line_pop() is fixed (p = 0 or 1), where no variance
# can be calibrated, so the line-scoped centring checks write known values with
# with_additive_terms(); it resolves the base through the generator's own
# .dae_resolve_base(). The pooled (p = 0.5) calls use the generator itself.
test_that("default base is the effect's own line pool; population-wide pools and warns", {
  pop <- make_two_line_pop("bt_inherit")
  on.exit(close_pop(pop), add = TRUE)
  pop <- with_additive_target(pop, "ADG", 1)

  pop |> get_table("genome_meta") |>
    with_additive_terms("ADG", effects = rep(1, 40), line_name = "A")
  pop |> get_table("genome_meta") |>
    with_additive_terms("ADG", effects = rep(1, 40), line_name = "B")

  # Line A is fixed at allele 0, line B at allele 1 -- each sees its own.
  expect_equal(stored_center(pop, "A"), 0)
  expect_equal(stored_center(pop, "B"), 1)

  # A population-wide effect on the default path is centered on the pooled
  # base -- the right answer for an effect that applies to the whole founder
  # base -- and warns, because pooling was not asked for.
  expect_warning(
    pop |> get_table("genome_meta") |>
      define_additive_effects("ADG", seed = 1),
    "pooled across")
  expect_equal(stored_center(pop, NULL), 0.5)
})

test_that("an explicit whole founder table pools on purpose and never warns", {
  pop <- make_two_line_pop("bt_forced")
  on.exit(close_pop(pop), add = TRUE)
  pop <- with_additive_target(pop, "ADG", 1)

  # Line-specific effect deliberately centered on the pooled base.
  expect_no_warning(
    pop |> get_table("genome_meta") |>
      define_additive_effects("ADG", line_name = "A", seed = 1,
                              base_tbl = get_table(pop, "founder_haplotypes")))
  expect_equal(stored_center(pop, "A"), 0.5)

  # Population-wide effect, explicit pooled base: same numbers, no warning.
  expect_no_warning(
    pop |> get_table("genome_meta") |>
      define_additive_effects("ADG", seed = 1,
                              base_tbl = get_table(pop, "founder_haplotypes")))
  expect_equal(stored_center(pop, NULL), 0.5)
})

test_that("default resolution: line pool -> shared pool -> error", {
  # (b) only a shared pool: a line-scoped effect falls back to it, silently.
  pop <- make_two_line_pop("bt_shared", lines = NA_character_)
  pop <- with_additive_target(pop, "ADG", 1)
  expect_no_warning(
    pop |> get_table("genome_meta") |>
      with_additive_terms("ADG", effects = rep(1, 40), line_name = "A"))
  expect_equal(stored_center(pop, "A"), 1)   # the NA-line pool is allele 1

  # (c) both exist: the named pool wins.
  DBI::dbWriteTable(pop$db_conn, "founder_haplotypes", data.frame(
    line_name = "A", haplotype_id = rep(1:20, times = 40),
    locus_name = rep(DBI::dbGetQuery(pop$db_conn,
      "SELECT locus_name FROM genome_meta ORDER BY locus_id")$locus_name,
      each = 20),
    allele = 0L, stringsAsFactors = FALSE), append = TRUE)
  pop |> get_table("genome_meta") |>
    with_additive_terms("ADG", effects = rep(1, 40), line_name = "A")
  expect_equal(stored_center(pop, "A"), 0)

  # Population-wide on the default path now sees two pools (NULL counts).
  expect_warning(
    pop |> get_table("genome_meta") |>
      define_additive_effects("ADG", seed = 1),
    "holds 2 pools")
  close_pop(pop)

  # A line with no named pool falls back to the shared pool when one exists --
  # pool-level fallback, the same rule as (b).
  pop <- make_two_line_pop("bt_fallback", lines = c("A", NA_character_))
  on.exit(close_pop(pop), add = TRUE)
  pop <- with_additive_target(pop, "ADG", 1)
  expect_no_warning(
    pop |> get_table("genome_meta") |>
      with_additive_terms("ADG", effects = rep(1, 40), line_name = "B"))
  expect_equal(stored_center(pop, "B"), 1)         # the NA-line pool is allele 1

  # (d) neither: loud, listing what exists.
  pop2 <- make_two_line_pop("bt_none")            # A and B, no NULL pool
  on.exit(close_pop(pop2), add = TRUE)
  pop2 <- with_additive_target(pop2, "ADG", 1)
  expect_error(
    pop2 |> get_table("genome_meta") |>
      define_additive_effects("ADG", line_name = "NOPE"),
    "No founder_haplotypes rows for line 'NOPE'. Available: 'A', 'B'")
})

test_that("per-line centering recovers the target that pooling misses", {
  # Line A: allele 0 fixed at half the loci, polymorphic at the rest, so the
  # within-line and pooled frequencies genuinely differ.
  set.seed(404)
  pop <- open_pop(pop_name = "bt_var", db_name = ":memory:") |>
    define_genome(n_loci = 60, n_chr = 2, chr_len_Mb = 50)
  on.exit(close_pop(pop), add = TRUE)
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
  pop <- ge_flat_view(pop)
  pop <- with_additive_target(pop, "ADG", 2)

  falconer <- function(ln) {
    e <- DBI::dbGetQuery(pop$db_conn, paste0(
      "SELECT center_value p, genome_value a FROM additive_flat ",
      "WHERE trait_name = 'ADG' AND line_name = '", ln, "'"))
    sum(2 * e$p * (1 - e$p) * e$a^2)
  }
  # Realized within-line variance implied by the stored effects and the line's
  # OWN allele frequencies -- what the simulation actually delivers.
  realised <- function(ln) {
    e <- DBI::dbGetQuery(pop$db_conn, paste0(
      "SELECT e.locus_name, e.genome_value a, f.p FROM additive_flat e ",
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

  # Now force pooling explicitly: the Falconer bookkeeping still "hits" the
  # target, but the variance actually realised within line A does not.
  set.seed(1)
  pop |> get_table("genome_meta") |>
    define_additive_effects("ADG", line_name = "A", seed = 1,
                            base_tbl = get_table(pop, "founder_haplotypes"))
  expect_equal(falconer("A"), 2, tolerance = 1e-8)
  expect_lt(realised("A"), 1.0)   # pooling under-scales: well short of 2
})

test_that("base_tbl is validated: class, same pop, and column projection", {
  pop <- make_two_line_pop("bt_valid")
  on.exit(close_pop(pop), add = TRUE)
  pop <- with_additive_target(pop, "ADG", 1)
  gm  <- pop |> get_table("genome_meta")

  expect_error(gm |> define_additive_effects("ADG", base_tbl = "A"),
               "must be a tidybreed_table")
  other <- make_two_line_pop("bt_other")
  on.exit(close_pop(other), add = TRUE)
  expect_error(gm |> define_additive_effects("ADG",
                 base_tbl = get_table(other, "founder_haplotypes")),
               "same pop as 'tbl'")
  expect_error(gm |> define_additive_effects("ADG",
                 base_tbl = get_table(pop, "founder_haplotypes") |>
                   dplyr::select(line_name)),
               "missing column\\(s\\) locus_name, allele")
  # A line name that is not a valid identifier is still refused up front.
  expect_error(gm |> define_additive_effects("ADG",
                 line_name = "A'; DROP TABLE genome_meta; --"),
               "Invalid line name")
  expect_true(DBI::dbExistsTable(pop$db_conn, "genome_meta"))
})

test_that("a selected QTL with no base copies errors, per trait under union", {
  set.seed(77)
  pop <- open_pop(pop_name = "bt_gap", db_name = ":memory:") |>
    define_genome(n_loci = 8, n_chr = 1, chr_len_Mb = 20)
  on.exit(close_pop(pop), add = TRUE)
  pop <- define_founder_haplotypes(pop, n_haplotypes = 20, line_name = "A",
                                   method = "fixed", allele_freq = 0.5)
  pop <- pop |> get_table("founder_haplotypes") |>
    add_founders(n_males = 2, n_females = 2, line_name = "A")
  pop <- define_trait(pop, "ADG")
  pop <- define_trait(pop, "BW")

  # Copies exist at loci 1-4 only.
  half <- get_table(pop, "ind_haplotype") |> dplyr::filter(locus_id <= 4L)

  # Single trait selecting loci 3..6: 5 and 6 have no base copies -> error,
  # naming them, before any effect is written.
  expect_error(
    pop |> get_table("genome_meta") |> dplyr::filter(locus_id %in% 3:6) |>
      define_additive_effects("ADG", G = 1, base_tbl = half),
    "no allele copies at 2 selected QTL \\(Locus_5, Locus_6\\)")
  expect_equal(DBI::dbGetQuery(pop$db_conn,
    "SELECT COUNT(*) n FROM genome_effects")$n, 0)

  # Fully covered selection is fine. Planted (uncalibrated, test-only) so the
  # union call below meets ADG's existing generated QTL and no stored target.
  pop |> get_table("genome_meta") |> dplyr::filter(locus_id %in% 1:4) |>
    plant_generated_additive("ADG", effects = rep(1, 4), base_tbl = half)

  # Union: ADG's existing QTL are 1-4 (covered); BW has none at this scope, so
  # its positive target variance cannot be delivered and nothing is stored.
  G <- diag(2); dimnames(G) <- list(c("ADG", "BW"), c("ADG", "BW"))
  expect_error(
    pop |> get_table("genome_meta") |>
      define_additive_effects(c("ADG", "BW"), G = G, method = "union",
                              base_tbl = half),
    "no existing generated additive effects.*cannot deliver its target variance 1")
  expect_equal(DBI::dbGetQuery(pop$db_conn,
    "SELECT COUNT(*) n FROM trait_var_comp")$n, 0)
  # Shared: every trait writes every candidate, so the gap now bites, per trait.
  expect_error(
    pop |> get_table("genome_meta") |>
      define_additive_effects(c("ADG", "BW"), G = G, method = "shared",
                              base_tbl = half),
    "no allele copies at 4 selected QTL for trait 'ADG'")
})

test_that("same seed reproduces itself under the new surface", {
  run <- function() {
    pop <- open_pop(pop_name = "bt_seed", db_name = ":memory:") |>
      define_genome(n_loci = 20, n_chr = 1, chr_len_Mb = 20)
    on.exit(close_pop(pop), add = TRUE)
    set.seed(3)
    pop <- define_founder_haplotypes(pop, n_haplotypes = 30, line_name = "A")
    pop <- ge_flat_view(pop)
    pop <- with_additive_target(pop, "ADG", 1)
    pop |> get_table("genome_meta") |>
      define_additive_effects("ADG", line_name = "A", seed = 42)
    DBI::dbGetQuery(pop$db_conn,
      "SELECT locus_name, genome_value FROM additive_flat ORDER BY locus_name")
  }
  expect_identical(run(), run())
})


test_that("defining effects for one line does not clobber another line's rows", {
  pop <- make_two_line_pop("bln_clobber")
  pop <- with_additive_target(pop, "ADG", 1)

  # Pooled base (p = 0.5): the lines alone are fixed and carry no variance.
  pooled <- get_table(pop, "founder_haplotypes")
  vals <- function(ln) DBI::dbGetQuery(pop$db_conn, paste0(
    "SELECT locus_name, genome_value FROM additive_flat WHERE trait_name = 'ADG' ",
    "AND line_name ", if (is.null(ln)) "IS NULL" else paste0("= '", ln, "'"),
    " ORDER BY locus_name"))
  pop |> get_table("genome_meta") |>
    define_additive_effects("ADG", line_name = "A", base_tbl = pooled, seed = 1)
  a1 <- vals("A")
  pop |> get_table("genome_meta") |>
    define_additive_effects("ADG", line_name = "B", base_tbl = pooled, seed = 2)
  pop |> get_table("genome_meta") |>
    define_additive_effects("ADG", base_tbl = pooled, seed = 3)

  counts <- DBI::dbGetQuery(pop$db_conn,
    "SELECT line_name, COUNT(*) AS n
       FROM additive_flat WHERE trait_name = 'ADG'
       GROUP BY line_name ORDER BY line_name NULLS LAST")

  expect_equal(nrow(counts), 3L)
  expect_true(all(counts$n == 40L))
  # Line A's variant is untouched by the line-B and common calls.
  expect_identical(vals("A"), a1)
  expect_false(isTRUE(all.equal(vals("B")$genome_value, a1$genome_value)))

  close_pop(pop)
})


test_that("centres and effects follow locus_id whatever order or projection tbl has", {
  # The written members are in locus_id order and p_base is indexed the same
  # way. A QTL table arranged descending and stripped of locus_id must not
  # shift a locus's centre onto its neighbour.
  pop <- open_pop(pop_name = "bt_order", db_name = ":memory:") |>
    define_genome(n_loci = 12, n_chr = 1, chr_len_Mb = 12)
  on.exit(close_pop(pop), add = TRUE)
  set.seed(7)
  pop <- define_founder_haplotypes(pop, n_haplotypes = 40, method = "uniform",
                                   line_name = "A")
  pop <- with_additive_target(pop, "ADG", 1)

  p <- extract_allele_freq(get_table(pop, "founder_haplotypes"))
  expect_gt(length(unique(round(p$allele_freq, 6))), 3L)   # p really varies

  qtl <- pop |> get_table("genome_meta") |>
    dplyr::filter(locus_id %% 2L == 0L) |>
    dplyr::arrange(dplyr::desc(pos_bp)) |>
    dplyr::select(locus_name, chr)
  expect_false("locus_id" %in% colnames(qtl$tbl))

  want <- p[p$locus_id %% 2L == 0L, ]           # locus_id order
  qtl |> define_additive_effects("ADG", line_name = "A", seed = 1)

  got <- DBI::dbGetQuery(pop$db_conn, "
    SELECT gm.locus_id, m.center_value
      FROM genome_effect_members m
      JOIN genome_effects e USING (id_genome_effect)
      JOIN genome_meta gm USING (locus_id)
     WHERE e.trait_name = 'ADG' ORDER BY gm.locus_id")
  expect_equal(got$locus_id, want$locus_id)
  expect_equal(got$center_value, want$allele_freq)
})


# ── Step-3 review (findings 1, 2): the stored target describes every retained
# generated term at its scope ─────────────────────────────────────────────────

dae_review_pop <- function(name) {
  set.seed(3401)
  open_pop(pop_name = name, db_name = ":memory:") |>
    define_genome(n_loci = 12, n_chr = 1, chr_len_Mb = 100) |>
    define_founder_haplotypes(n_haplotypes = 100) |>
    get_table("founder_haplotypes") |>
    add_founders(n_males = 30, n_females = 30, line_name = "A")
}

# sum n_eligible p q a^2 over a trait's generated additive terms
dae_genic <- function(pop, trait) {
  x <- DBI::dbGetQuery(pop$db_conn,
    "SELECT m.center_value AS p, e.genome_value AS a,
            COUNT(o.origin_slot) AS n_origin
     FROM genome_effect_members m
     JOIN genome_effects e USING (id_genome_effect)
     LEFT JOIN genome_effect_member_origins o
       ON o.id_genome_effect = m.id_genome_effect
      AND o.member_slot = m.member_slot AND o.parent_origin IS NOT NULL
     WHERE e.trait_name = ?
     GROUP BY m.id_genome_effect, m.member_slot, m.center_value, e.genome_value",
    params = list(trait))
  sum(ifelse(x$n_origin > 0, 1, 2) * x$p * (1 - x$p) * x$a^2)
}

test_that("trait_var_comp_tbl must select the block the call's scope reads (finding 1)", {
  pop <- dae_review_pop("dae_tvc_scope")
  on.exit(close_pop(pop), add = TRUE)
  pop <- define_trait(pop, "T")
  pop <- suppressMessages(define_effect_cov_matrix(pop, "additive", 1, trait_name = "T"))
  pop <- suppressMessages(define_effect_cov_matrix(pop, "additive", 100,
                                                   trait_name = "T", line_name = "A"))
  pop <- suppressMessages(define_effect_cov_matrix(pop, "additive", 7,
                                                   trait_name = "T", line_name = "B"))
  tvc <- function(...) get_table(pop, "trait_var_comp") |>
    dplyr::filter(effect_name == "additive", ...)
  gm <- get_table(pop, "genome_meta")
  n_terms <- function() DBI::dbGetQuery(pop$db_conn,
    "SELECT COUNT(*) AS n FROM genome_effects")$n

  seed_before <- .Random.seed
  # Common call, line-A rows.
  expect_error(gm |> define_additive_effects("T", seed = 11, warn_bounds = NULL,
                 trait_var_comp_tbl = tvc(line_name == "A")),
               "selects the line 'A'.*described by the population-wide block")
  # Line-A call, line-B rows.
  expect_error(gm |> define_additive_effects("T", line_name = "A", seed = 11,
                 warn_bounds = NULL, trait_var_comp_tbl = tvc(line_name == "B")),
               "selects the line 'B'.*line 'A' block")
  # Line-A call, population-wide rows although line A has its own block.
  expect_error(gm |> define_additive_effects("T", line_name = "A", seed = 11,
                 warn_bounds = NULL, trait_var_comp_tbl = tvc(is.na(line_name))),
               "selects the population-wide.*has a block of its own")
  expect_identical(.Random.seed, seed_before)
  expect_equal(n_terms(), 0)

  # Same scope: accepted and calibrated to that block.
  pop <- suppressMessages(gm |> define_additive_effects("T", seed = 11,
    warn_bounds = NULL, trait_var_comp_tbl = tvc(is.na(line_name))))
  expect_equal(dae_genic(pop, "T"), 1, tolerance = 1e-8)
})

test_that("trait_var_comp_tbl accepts the line -> population-wide fallback (finding 1)", {
  pop <- dae_review_pop("dae_tvc_fallback")
  on.exit(close_pop(pop), add = TRUE)
  pop <- with_additive_target(pop, "T", 1)
  pop <- suppressMessages(get_table(pop, "genome_meta") |>
    define_additive_effects("T", line_name = "A", seed = 1, warn_bounds = NULL,
      trait_var_comp_tbl = get_table(pop, "trait_var_comp") |>
        dplyr::filter(effect_name == "additive", is.na(line_name))))
  expect_equal(DBI::dbGetQuery(pop$db_conn,
    "SELECT COUNT(*) AS n FROM genome_effect_member_origins
     WHERE line_name = 'A'")$n, 12)
})

test_that("a zero-target union trait loses its old terms at the scope (finding 2)", {
  pop <- dae_review_pop("dae_union_zero")
  on.exit(close_pop(pop), add = TRUE)
  pop <- define_trait(pop, "T")
  pop <- define_trait(pop, "U")
  gm <- get_table(pop, "genome_meta")
  pop <- suppressMessages(gm |> dplyr::filter(locus_id <= 6L) |>
    define_additive_effects("T", G = 1, seed = 1, warn_bounds = NULL))
  pop <- suppressMessages(gm |> dplyr::filter(locus_id > 6L) |>
    define_additive_effects("U", G = 1, seed = 2, warn_bounds = NULL))
  expect_equal(dae_genic(pop, "U"), 1, tolerance = 1e-8)
  pop <- suppressMessages(get_table(pop, "trait_var_comp") |>
    remove_rows(confirm_all = TRUE))

  G <- diag(c(1, 0)); dimnames(G) <- list(c("T", "U"), c("T", "U"))
  expect_message(pop <- suppressWarnings(gm |> dplyr::filter(locus_id <= 6L) |>
    define_additive_effects(c("T", "U"), G = G, method = "union", seed = 13,
                            warn_bounds = NULL)),
    "Trait 'U' has no QTL.*generated effects at this scope are removed")
  # The stored diagonal (0) describes what U's retained model delivers (0).
  expect_equal(get_trait_var(pop, "additive", "U"), 0)
  expect_equal(DBI::dbGetQuery(pop$db_conn,
    "SELECT COUNT(*) AS n FROM genome_effects WHERE trait_name = 'U'")$n, 0)
  expect_equal(dae_genic(pop, "U"), 0)
  expect_equal(dae_genic(pop, "T"), 1, tolerance = 1e-8)
})
