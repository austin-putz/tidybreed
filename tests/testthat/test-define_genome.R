test_that("define_genome() adds genome tables to pop$tables", {
  pop <- open_pop(db_name = ":memory:") |>
    define_genome(n_loci = 100, n_chr = 2, chr_len_Mb = 100)
  on.exit(close_pop(pop))

  expect_true("genome_meta"   %in% pop$tables)
  expect_true("genome_map"    %in% pop$tables)
  expect_true("ind_haplotype" %in% pop$tables)
  expect_true("ind_genotype"  %in% pop$tables)
  expect_true("chr_meta"      %in% pop$tables)
})

test_that("define_genome() genome_meta has correct row count and columns", {
  n <- 200
  pop <- open_pop(db_name = ":memory:") |>
    define_genome(n_loci = n, n_chr = 4, chr_len_Mb = 50)
  on.exit(close_pop(pop))

  gm <- get_table(pop, "genome_meta") |> dplyr::collect()
  expect_equal(nrow(gm), n)
  # Physical coordinate is pos_bp (BIGINT); genetic pos_cM lives in genome_map.
  expect_true(all(c("locus_id", "locus_name", "chr", "chr_name", "pos_bp") %in%
                    names(gm)))
  # Old columns are gone.
  expect_false("pos_Mb"         %in% names(gm))
  expect_false("pos_cM"         %in% names(gm))
  expect_false("introduced_gen" %in% names(gm))

  # pos_bp is a genuine BIGINT column (not DOUBLE/INTEGER).
  pt <- DBI::dbGetQuery(pop$db_conn,
    "SELECT data_type FROM information_schema.columns
     WHERE table_name = 'genome_meta' AND column_name = 'pos_bp'")$data_type
  expect_equal(pt, "BIGINT")
})

test_that("define_genome() distributes loci evenly across chromosomes", {
  pop <- open_pop(db_name = ":memory:") |>
    define_genome(n_loci = 100, n_chr = 4, chr_len_Mb = 100)
  on.exit(close_pop(pop))

  gm <- get_table(pop, "genome_meta") |> dplyr::collect()
  chr_counts <- table(gm$chr)
  expect_equal(length(chr_counts), 4)
  # Each chromosome should have approximately 25 loci
  expect_true(all(chr_counts >= 24 & chr_counts <= 26))
})

test_that("define_genome() scalar chr_len_Mb applies to all chromosomes", {
  pop <- open_pop(db_name = ":memory:") |>
    define_genome(n_loci = 100, n_chr = 3, chr_len_Mb = 150)
  on.exit(close_pop(pop))

  gm <- get_table(pop, "genome_meta") |> dplyr::collect()
  for (c in 1:3) {
    chr_pos <- gm$pos_bp[gm$chr == c]        # base pairs now
    expect_true(max(chr_pos) <= 150 * 1e6)
    expect_true(min(chr_pos) > 0)            # 1-based, strictly positive
    expect_false(anyDuplicated(chr_pos) > 0) # strictly increasing, no dup
  }
})

test_that("define_genome() vector chr_len_Mb respected per chromosome", {
  lengths <- c(50, 100, 200)
  pop <- open_pop(db_name = ":memory:") |>
    define_genome(n_loci = 90, n_chr = 3, chr_len_Mb = lengths)
  on.exit(close_pop(pop))

  gm <- get_table(pop, "genome_meta") |> dplyr::collect()
  for (i in seq_along(lengths)) {
    chr_pos <- gm$pos_bp[gm$chr == i]        # base pairs
    expect_true(max(chr_pos) <= lengths[i] * 1e6)
  }
})

test_that("define_genome() custom locus_names are stored correctly", {
  names_vec <- paste0("SNP_", 1:50)
  pop <- open_pop(db_name = ":memory:") |>
    define_genome(n_loci = 50, n_chr = 2, chr_len_Mb = 100,
                  locus_names = names_vec)
  on.exit(close_pop(pop))

  gm <- get_table(pop, "genome_meta") |> dplyr::collect()
  expect_equal(sort(gm$locus_name), sort(names_vec))
})

test_that("define_genome() custom chr_names are stored correctly", {
  chr_names_vec <- c("chr1", "chr2", "chrX")
  pop <- open_pop(db_name = ":memory:") |>
    define_genome(n_loci = 60, n_chr = 3, chr_len_Mb = 100,
                  chr_names = chr_names_vec)
  on.exit(close_pop(pop))

  gm <- get_table(pop, "genome_meta") |> dplyr::collect()
  expect_equal(sort(unique(gm$chr_name)), sort(chr_names_vec))
})

test_that("define_genome() ind_haplotype and ind_genotype are empty tables", {
  pop <- open_pop(db_name = ":memory:") |>
    define_genome(n_loci = 50, n_chr = 2, chr_len_Mb = 100)
  on.exit(close_pop(pop))

  n_hap  <- DBI::dbGetQuery(pop$db_conn,
                             "SELECT COUNT(*) AS n FROM ind_haplotype")$n
  n_geno <- DBI::dbGetQuery(pop$db_conn,
                             "SELECT COUNT(*) AS n FROM ind_genotype")$n
  expect_equal(n_hap,  0L)
  expect_equal(n_geno, 0L)
})

test_that("define_genome() registers genome table descriptions in _schema_meta", {
  pop <- open_pop(db_name = ":memory:") |>
    define_genome(n_loci = 50, n_chr = 2, chr_len_Mb = 100)
  on.exit(close_pop(pop))

  sm <- DBI::dbGetQuery(
    pop$db_conn,
    "SELECT DISTINCT table_name FROM _schema_meta WHERE object_type = 'table'"
  )$table_name
  expect_true("genome_meta"   %in% sm)
  expect_true("ind_haplotype" %in% sm)
  expect_true("ind_genotype"  %in% sm)
  expect_true("chr_meta"      %in% sm)
})

test_that("define_genome() is pipe-friendly and returns pop", {
  pop <- open_pop(db_name = ":memory:") |>
    define_genome(n_loci = 50, n_chr = 2, chr_len_Mb = 100)
  on.exit(close_pop(pop))

  expect_s3_class(pop, "tidybreed_pop")
})
