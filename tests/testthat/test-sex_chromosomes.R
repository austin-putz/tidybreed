# Stage 4 (plans/refactor_haplotype.md) -- sex chromosomes and organelles via
# chr_meta/define_chr(). Hand-verified row counts and inheritance patterns for
# X/Y, Z/W, MT, and X0, plus mixed-genome dosage/export and RNG-neutrality of
# non-recombining inheritance.

# Row counts of ind_haplotype for one individual on one chromosome, by
# parent_origin.
.chr_rows <- function(pop, id, chr) {
  DBI::dbGetQuery(pop$db_conn, sprintf(
    "SELECT h.parent_origin, COUNT(*) AS n FROM ind_haplotype h
     JOIN genome_meta gm USING(locus_id)
     WHERE h.id_ind = '%s' AND gm.chr_name = '%s'
     GROUP BY h.parent_origin ORDER BY h.parent_origin", id, chr))
}

.chr_alleles <- function(pop, id, chr) {
  DBI::dbGetQuery(pop$db_conn, sprintf(
    "SELECT gm.locus_id, h.allele, h.line_origin FROM ind_haplotype h
     JOIN genome_meta gm USING(locus_id)
     WHERE h.id_ind = '%s' AND gm.chr_name = '%s'
     ORDER BY gm.locus_id", id, chr))
}


# ---------------------------------------------------------------------------
# X/Y (mammalian)
# ---------------------------------------------------------------------------

test_that("X/Y: founders get correct copy counts, offspring inherit patrilineal Y", {
  set.seed(101)
  pop <- open_pop(pop_name = "xy", db_name = ":memory:") |>
    define_genome(n_loci = 12, n_chr = 3, chr_names = c("1", "X", "Y"), chr_len_Mb = 50) |>
    define_chr("X", copy_mode_M = "half", copy_mode_F = "full",
               hemi_parent = "parent_2", recombines = TRUE) |>
    define_chr("Y", copy_mode_M = "half", copy_mode_F = "none",
               hemi_parent = "parent_1", recombines = FALSE) |>
    define_founder_haplotypes(n_haplotypes = 20, method = "fixed") |>
    get_table("founder_haplotypes") |>
    add_founders(n_males = 4, n_females = 4, line_name = "A")
  on.exit(close_pop(pop), add = TRUE)

  ids <- DBI::dbGetQuery(pop$db_conn, "SELECT id_ind, sex FROM ind_meta")
  males   <- ids$id_ind[ids$sex == "M"]
  females <- ids$id_ind[ids$sex == "F"]

  for (m in males) {
    x <- .chr_rows(pop, m, "X")
    expect_equal(x$parent_origin, 2L)
    expect_equal(x$n, 4L)
    y <- .chr_rows(pop, m, "Y")
    expect_equal(y$parent_origin, 1L)
    expect_equal(y$n, 4L)
  }
  for (f in females) {
    x <- .chr_rows(pop, f, "X")
    expect_equal(x$parent_origin, c(1L, 2L))
    expect_equal(x$n, c(4L, 4L))
    y <- .chr_rows(pop, f, "Y")
    expect_equal(nrow(y), 0L)
  }

  sire <- males[1]
  dam  <- females[1]
  matings_f1 <- tibble::tibble(
    id_parent_1 = sire, id_parent_2 = dam, sex = c("M", "F"), line_name = "F1"
  )
  set.seed(202)
  pop <- add_offspring(pop, matings_f1)

  f1 <- DBI::dbGetQuery(pop$db_conn, "SELECT id_ind, sex FROM ind_meta WHERE line_name = 'F1'")
  son      <- f1$id_ind[f1$sex == "M"]
  daughter <- f1$id_ind[f1$sex == "F"]

  sire_y <- .chr_alleles(pop, sire, "Y")
  son_y  <- .chr_alleles(pop, son, "Y")
  expect_equal(son_y$allele, sire_y$allele)
  expect_equal(son_y$line_origin, sire_y$line_origin)
  expect_equal(nrow(.chr_rows(pop, daughter, "Y")), 0L)

  # Grandson: Y still identical two generations later (no recombination, ever)
  granddam <- females[2]
  matings_f2 <- tibble::tibble(id_parent_1 = son, id_parent_2 = granddam,
                               sex = "M", line_name = "F2")
  set.seed(303)
  pop <- add_offspring(pop, matings_f2)
  grandson <- DBI::dbGetQuery(pop$db_conn, "SELECT id_ind FROM ind_meta WHERE line_name = 'F2'")$id_ind
  grandson_y <- .chr_alleles(pop, grandson, "Y")
  expect_equal(grandson_y$allele, sire_y$allele)
  expect_true(all(grandson_y$line_origin == "A"))
})


# ---------------------------------------------------------------------------
# Z/W (avian) -- mirror image: females heterogametic (matrilineal W)
# ---------------------------------------------------------------------------

test_that("Z/W: founders get correct copy counts, offspring inherit matrilineal W", {
  set.seed(111)
  pop <- open_pop(pop_name = "zw", db_name = ":memory:") |>
    define_genome(n_loci = 12, n_chr = 3, chr_names = c("1", "Z", "W"), chr_len_Mb = 50) |>
    define_chr("Z", copy_mode_M = "full", copy_mode_F = "half",
               hemi_parent = "parent_1", recombines = TRUE) |>
    define_chr("W", copy_mode_M = "none", copy_mode_F = "half",
               hemi_parent = "parent_2", recombines = FALSE) |>
    define_founder_haplotypes(n_haplotypes = 20, method = "fixed") |>
    get_table("founder_haplotypes") |>
    add_founders(n_males = 4, n_females = 4, line_name = "A")
  on.exit(close_pop(pop), add = TRUE)

  ids <- DBI::dbGetQuery(pop$db_conn, "SELECT id_ind, sex FROM ind_meta")
  males   <- ids$id_ind[ids$sex == "M"]
  females <- ids$id_ind[ids$sex == "F"]

  for (m in males) {
    z <- .chr_rows(pop, m, "Z")
    expect_equal(z$parent_origin, c(1L, 2L))
    expect_equal(z$n, c(4L, 4L))
    expect_equal(nrow(.chr_rows(pop, m, "W")), 0L)
  }
  for (f in females) {
    z <- .chr_rows(pop, f, "Z")
    expect_equal(z$parent_origin, 1L)
    expect_equal(z$n, 4L)
    w <- .chr_rows(pop, f, "W")
    expect_equal(w$parent_origin, 2L)
    expect_equal(w$n, 4L)
  }

  dam  <- females[1]
  sire <- males[1]
  matings_f1 <- tibble::tibble(
    id_parent_1 = sire, id_parent_2 = dam, sex = c("M", "F"), line_name = "F1"
  )
  set.seed(212)
  pop <- add_offspring(pop, matings_f1)

  f1 <- DBI::dbGetQuery(pop$db_conn, "SELECT id_ind, sex FROM ind_meta WHERE line_name = 'F1'")
  son      <- f1$id_ind[f1$sex == "M"]
  daughter <- f1$id_ind[f1$sex == "F"]

  dam_w <- .chr_alleles(pop, dam, "W")
  daughter_w <- .chr_alleles(pop, daughter, "W")
  expect_equal(daughter_w$allele, dam_w$allele)
  expect_equal(nrow(.chr_rows(pop, son, "W")), 0L)

  granddam_line <- males[2]  # a second sire, to mate with the F1 daughter
  matings_f2 <- tibble::tibble(id_parent_1 = granddam_line, id_parent_2 = daughter,
                               sex = "F", line_name = "F2")
  set.seed(313)
  pop <- add_offspring(pop, matings_f2)
  granddaughter <- DBI::dbGetQuery(pop$db_conn, "SELECT id_ind FROM ind_meta WHERE line_name = 'F2'")$id_ind
  granddaughter_w <- .chr_alleles(pop, granddaughter, "W")
  expect_equal(granddaughter_w$allele, dam_w$allele)
})


# ---------------------------------------------------------------------------
# Mitochondria: both sexes get exactly 1 copy, always maternal, never
# recombined, and inheritance consumes zero RNG draws.
# ---------------------------------------------------------------------------

test_that("MT: both sexes get 1 copy, strictly maternal, never recombined", {
  set.seed(121)
  pop <- open_pop(pop_name = "mt", db_name = ":memory:") |>
    define_genome(n_loci = 12, n_chr = 2, chr_names = c("1", "MT"), chr_len_Mb = 50) |>
    define_chr("MT", copy_mode_M = "half", copy_mode_F = "half",
               hemi_parent = "parent_2", recombines = FALSE) |>
    define_founder_haplotypes(n_haplotypes = 20, method = "fixed") |>
    get_table("founder_haplotypes") |>
    add_founders(n_males = 3, n_females = 3, line_name = "A")
  on.exit(close_pop(pop), add = TRUE)

  ids <- DBI::dbGetQuery(pop$db_conn, "SELECT id_ind, sex FROM ind_meta")
  for (id in ids$id_ind) {
    mt <- .chr_rows(pop, id, "MT")
    expect_equal(mt$parent_origin, 2L)
    expect_equal(mt$n, 6L)  # 12 loci / 2 chr
  }

  sire <- ids$id_ind[ids$sex == "M"][1]
  dam  <- ids$id_ind[ids$sex == "F"][1]
  matings <- tibble::tibble(
    id_parent_1 = sire, id_parent_2 = dam, sex = c("M", "F"), line_name = "F1"
  )
  set.seed(222)
  pop <- add_offspring(pop, matings)

  f1 <- DBI::dbGetQuery(pop$db_conn, "SELECT id_ind FROM ind_meta WHERE line_name = 'F1'")$id_ind
  dam_mt <- .chr_alleles(pop, dam, "MT")
  sire_mt <- .chr_alleles(pop, sire, "MT")
  for (id in f1) {
    off_mt <- .chr_alleles(pop, id, "MT")
    expect_equal(off_mt$allele, dam_mt$allele)
    # Never matches the father's MT (unless coincidentally identical alleles —
    # guard with a check that dam/sire actually differ somewhere first)
    if (!identical(dam_mt$allele, sire_mt$allele)) {
      expect_false(identical(off_mt$allele, sire_mt$allele))
    }
  }
})

test_that("pass_through_gamete() makes zero RNG draws when only one copy exists (k == 1)", {
  # This is the exact mechanism behind non-recombining, strictly hemizygous
  # inheritance (patrilineal Y, matrilineal MT): the contributing parent has
  # only one stored copy, so there is nothing to choose between.
  hap <- matrix(c(1L, 0L, 1L), nrow = 1)
  lo  <- matrix(c("A", "A", "A"), nrow = 1)

  set.seed(4242)
  before <- runif(1)

  set.seed(4242)
  g <- pass_through_gamete(hap, lo)
  after <- runif(1)

  expect_equal(before, after)
  expect_equal(g$allele, c(1L, 0L, 1L))
  expect_equal(g$line_origin, c("A", "A", "A"))
})

test_that("pass_through_gamete() draws exactly one value when two copies exist (k == 2, non-recombining)", {
  hap <- matrix(c(1L, 0L, 1L, 0L, 1L, 0L), nrow = 2, byrow = TRUE)
  lo  <- matrix(c("A", "A", "A", "B", "B", "B"), nrow = 2, byrow = TRUE)

  set.seed(99)
  chosen_row <- sample.int(2, size = 1)

  set.seed(99)
  g <- pass_through_gamete(hap, lo)

  expect_equal(g$allele, hap[chosen_row, ])
  expect_equal(g$line_origin, lo[chosen_row, ])
})


# ---------------------------------------------------------------------------
# X0 (no partner chromosome defined at all)
# ---------------------------------------------------------------------------

test_that("X0: a hemizygous chromosome works with no partner chromosome defined", {
  set.seed(131)
  pop <- open_pop(pop_name = "x0", db_name = ":memory:") |>
    define_genome(n_loci = 12, n_chr = 3, chr_names = c("1", "2", "X"), chr_len_Mb = 50) |>
    define_chr("X", copy_mode_M = "half", copy_mode_F = "full",
               hemi_parent = "parent_2", recombines = TRUE) |>
    define_founder_haplotypes(n_haplotypes = 20, method = "fixed") |>
    get_table("founder_haplotypes") |>
    add_founders(n_males = 3, n_females = 3, line_name = "A")
  on.exit(close_pop(pop), add = TRUE)

  ids <- DBI::dbGetQuery(pop$db_conn, "SELECT id_ind, sex FROM ind_meta")
  for (m in ids$id_ind[ids$sex == "M"]) {
    x <- .chr_rows(pop, m, "X")
    expect_equal(x$parent_origin, 2L)
    expect_equal(x$n, 4L)
  }
  for (f in ids$id_ind[ids$sex == "F"]) {
    x <- .chr_rows(pop, f, "X")
    expect_equal(x$parent_origin, c(1L, 2L))
    expect_equal(x$n, c(4L, 4L))
  }

  matings <- tibble::tibble(
    id_parent_1 = ids$id_ind[ids$sex == "M"][1],
    id_parent_2 = ids$id_ind[ids$sex == "F"][1],
    sex = c("M", "F"), line_name = "F1"
  )
  expect_no_error(pop <- add_offspring(pop, matings))
})


# ---------------------------------------------------------------------------
# Mixed autosome + sex-chromosome dosage / export
# ---------------------------------------------------------------------------

test_that("add_dosage()/extract_genotypes() handle a mixed autosome+X/Y genome correctly", {
  set.seed(141)
  pop <- open_pop(pop_name = "mixed_export", db_name = ":memory:") |>
    define_genome(n_loci = 12, n_chr = 3, chr_names = c("1", "X", "Y"), chr_len_Mb = 50) |>
    define_chr("X", copy_mode_M = "half", copy_mode_F = "full",
               hemi_parent = "parent_2", recombines = TRUE) |>
    define_chr("Y", copy_mode_M = "half", copy_mode_F = "none",
               hemi_parent = "parent_1", recombines = FALSE) |>
    define_founder_haplotypes(n_haplotypes = 20, method = "fixed") |>
    get_table("founder_haplotypes") |>
    add_founders(n_males = 3, n_females = 3, line_name = "A")
  on.exit(close_pop(pop), add = TRUE)

  ids <- DBI::dbGetQuery(pop$db_conn, "SELECT id_ind, sex FROM ind_meta")
  males   <- ids$id_ind[ids$sex == "M"]
  females <- ids$id_ind[ids$sex == "F"]

  pop <- pop |> get_table("ind_meta") |> add_dosage()

  auto_locus_ids <- DBI::dbGetQuery(pop$db_conn,
    "SELECT locus_id FROM genome_meta WHERE chr_name = '1'")$locus_id
  x_locus_ids <- DBI::dbGetQuery(pop$db_conn,
    "SELECT locus_id FROM genome_meta WHERE chr_name = 'X'")$locus_id
  y_locus_ids <- DBI::dbGetQuery(pop$db_conn,
    "SELECT locus_id FROM genome_meta WHERE chr_name = 'Y'")$locus_id

  geno <- DBI::dbGetQuery(pop$db_conn, "SELECT * FROM ind_genotype")

  # Autosome dosage: 0/1/2, present for everyone (regression, unaffected)
  auto_geno <- geno[geno$locus_id %in% auto_locus_ids, ]
  expect_true(all(auto_geno$dosage_value %in% c(0L, 1L, 2L)))
  expect_equal(nrow(auto_geno), length(auto_locus_ids) * 6L)

  # X dosage: males 0/1, females 0/1/2
  male_x <- geno[geno$locus_id %in% x_locus_ids & geno$id_ind %in% males, ]
  female_x <- geno[geno$locus_id %in% x_locus_ids & geno$id_ind %in% females, ]
  expect_true(all(male_x$dosage_value %in% c(0L, 1L)))
  expect_true(all(female_x$dosage_value %in% c(0L, 1L, 2L)))

  # Y: no ind_genotype row at all for females (partial-cache semantics)
  female_y <- geno[geno$locus_id %in% y_locus_ids & geno$id_ind %in% females, ]
  expect_equal(nrow(female_y), 0L)
  male_y <- geno[geno$locus_id %in% y_locus_ids & geno$id_ind %in% males, ]
  expect_equal(nrow(male_y), length(y_locus_ids) * 3L)

  # extract_genotypes(): dense matrix, 0 (not NA) for females at Y loci
  pop2 <- pop |> get_table("genome_meta") |> define_chip("all")
  pop2 <- pop2 |> get_table("ind_meta") |> add_genotypes("all")
  ext <- pop2 |> get_table("ind_meta") |> extract_genotypes("all")

  y_cols <- paste0("locus_", y_locus_ids)
  female_rows <- ext[ext$id_ind %in% females, y_cols, drop = FALSE]
  expect_true(all(!is.na(unlist(female_rows))))
  expect_true(all(unlist(female_rows) == 0L))
})


# ---------------------------------------------------------------------------
# Ploidy guard interaction (Stage 4 semantic change)
# ---------------------------------------------------------------------------

test_that("a sex-linked chr_meta configuration does not trip the ploidy-2 guard", {
  set.seed(151)
  pop <- open_pop(pop_name = "guard_interaction", db_name = ":memory:") |>
    define_genome(n_loci = 8, n_chr = 2, chr_names = c("1", "X"), chr_len_Mb = 50) |>
    define_chr("X", copy_mode_M = "half", copy_mode_F = "full", hemi_parent = "parent_2") |>
    define_founder_haplotypes(n_haplotypes = 20, method = "fixed") |>
    get_table("founder_haplotypes") |>
    add_founders(n_males = 3, n_females = 3, line_name = "A")
  on.exit(close_pop(pop), add = TRUE)

  expect_no_error(pop <- pop |> get_table("ind_meta") |> add_dosage())
  expect_no_error(
    pop2 <- pop |> get_table("genome_meta") |> define_chip("all")
  )
})
