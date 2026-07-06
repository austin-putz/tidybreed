# Genome position coordinates & genetic map (pos_bp + genome_map migration).
# Forward-only reproducibility (pre-1.0.0): we assert same-seed determinism of the
# CURRENT implementation, never a match to any prior version.

# ---------------------------------------------------------------------------
# Position generation & pos_cM derivation
# ---------------------------------------------------------------------------

test_that("define_genome() writes a default genome_map with pos_cM = pos_bp/1e6 * rate", {
  pop <- open_pop(pop_name = "gm", db_name = ":memory:") |>
    define_genome(n_loci = 30, n_chr = 3, chr_len_Mb = c(20, 15, 10))
  on.exit(close_pop(pop), add = TRUE)

  m <- resolve_genome_map(pop$db_conn)
  expect_equal(nrow(m), 30L)                       # one row per locus
  expect_false(anyNA(m$pos_cM))                    # every locus resolved
  expect_equal(m$locus_id, sort(m$locus_id))       # locus_id-ordered
  expect_equal(m$pos_cM, m$pos_bp / 1e6, tolerance = 1e-9)   # default rate 1

  # genome_map default rows: sex/line NULL, map_name 'default'
  raw <- DBI::dbGetQuery(pop$db_conn, "SELECT * FROM genome_map")
  expect_equal(nrow(raw), 30L)
  expect_true(all(is.na(raw$sex)))
  expect_true(all(is.na(raw$line_name)))
  expect_true(all(raw$map_name == "default"))
})

test_that("cM_per_Mb scales pos_cM; per-chromosome rate honored", {
  pop <- open_pop(pop_name = "gm2", db_name = ":memory:") |>
    define_genome(n_loci = 20, n_chr = 2, chr_len_Mb = 10, cM_per_Mb = c(2, 0.5))
  on.exit(close_pop(pop), add = TRUE)

  m <- resolve_genome_map(pop$db_conn)
  chr1 <- m[m$chr == 1, ]; chr2 <- m[m$chr == 2, ]
  expect_equal(chr1$pos_cM, chr1$pos_bp / 1e6 * 2,   tolerance = 1e-9)
  expect_equal(chr2$pos_cM, chr2$pos_bp / 1e6 * 0.5, tolerance = 1e-9)
})

test_that("cM_per_Mb must be finite and strictly positive", {
  pop <- open_pop(pop_name = "gm3", db_name = ":memory:")
  on.exit(close_pop(pop), add = TRUE)
  expect_error(define_genome(pop, n_loci = 10, n_chr = 1, chr_len_Mb = 10, cM_per_Mb = 0),
               "positive")
  expect_error(define_genome(pop, n_loci = 10, n_chr = 1, chr_len_Mb = 10, cM_per_Mb = -1),
               "positive")
  expect_error(define_genome(pop, n_loci = 10, n_chr = 1, chr_len_Mb = 10, cM_per_Mb = NA_real_),
               "positive")
})

test_that("pos_bp validation errors on a too-short chromosome and writes nothing", {
  pop <- open_pop(pop_name = "gm4", db_name = ":memory:")
  on.exit(close_pop(pop), add = TRUE)
  # 1000 bp span cannot hold 100000 strictly-increasing integer positions
  expect_error(
    define_genome(pop, n_loci = 100000, n_chr = 1, chr_len_Mb = 0.001),
    "too short"
  )
  expect_false("genome_meta" %in% DBI::dbListTables(pop$db_conn))
})

# ---------------------------------------------------------------------------
# Map resolution precedence & validator strictness
# ---------------------------------------------------------------------------

test_that("resolve_genome_map() applies sex precedence, default fallback", {
  pop <- open_pop(pop_name = "gm5", db_name = ":memory:") |>
    define_genome(n_loci = 3, n_chr = 1, chr_len_Mb = 10)
  on.exit(close_pop(pop), add = TRUE)
  con <- pop$db_conn

  def <- resolve_genome_map(con)
  # Add a MONOTONIC male-specific override for locus 1 only.
  DBI::dbExecute(con, sprintf(
    "INSERT INTO genome_map VALUES (100, 1, 'Locus_1', 'M', NULL, 'default', %f)",
    def$pos_cM[1] - 0.1))
  validate_genome_map(con)

  mM <- resolve_genome_map(con, sex = "M")
  mF <- resolve_genome_map(con, sex = "F")
  expect_equal(mM$pos_cM[1], def$pos_cM[1] - 0.1)  # M gets the override
  expect_equal(mF$pos_cM[1], def$pos_cM[1])         # F falls back to default
  expect_equal(mM$pos_cM[2:3], def$pos_cM[2:3])     # other loci fall back
})

test_that("resolve_genome_map() errors when a sex override makes the map non-monotonic", {
  pop <- open_pop(pop_name = "gm6", db_name = ":memory:") |>
    define_genome(n_loci = 3, n_chr = 1, chr_len_Mb = 10)
  on.exit(close_pop(pop), add = TRUE)
  con <- pop$db_conn
  # Override locus 1 with a huge cM that breaks monotonicity after fallback.
  DBI::dbExecute(con,
    "INSERT INTO genome_map VALUES (100, 1, 'Locus_1', 'M', NULL, 'default', 999.0)")
  expect_error(resolve_genome_map(con, sex = "M"), "monotonic")
})

test_that("resolve_genome_map() applies line precedence over sex-agnostic default", {
  pop <- open_pop(pop_name = "gm5b", db_name = ":memory:") |>
    define_genome(n_loci = 3, n_chr = 1, chr_len_Mb = 10)
  on.exit(close_pop(pop), add = TRUE)
  con <- pop$db_conn
  def <- resolve_genome_map(con)
  # line-specific override for locus 1 (monotonic)
  DBI::dbExecute(con, sprintf(
    "INSERT INTO genome_map VALUES (100, 1, 'Locus_1', NULL, 'LineX', 'default', %f)",
    def$pos_cM[1] - 0.1))
  validate_genome_map(con)
  mX <- resolve_genome_map(con, line_name = "LineX")
  mDef <- resolve_genome_map(con, line_name = "OtherLine")
  expect_equal(mX$pos_cM[1], def$pos_cM[1] - 0.1)   # LineX gets its row
  expect_equal(mDef$pos_cM[1], def$pos_cM[1])        # unknown line falls back to default
})

test_that("resolve_genome_map() errors when a locus is missing from the resolved map", {
  pop <- open_pop(pop_name = "gm5c", db_name = ":memory:") |>
    define_genome(n_loci = 4, n_chr = 1, chr_len_Mb = 10)
  on.exit(close_pop(pop), add = TRUE)
  con <- pop$db_conn
  DBI::dbExecute(con, "DELETE FROM genome_map WHERE locus_id = 2")
  expect_error(resolve_genome_map(con), "no genome_map row")
})

test_that("validate_genome_map() catches invalid sex, NULL map_name, non-finite pos_cM, dup id", {
  mk <- function(nm) open_pop(pop_name = nm, db_name = ":memory:") |>
    define_genome(n_loci = 4, n_chr = 1, chr_len_Mb = 10)

  p1 <- mk("gmv1"); on.exit(close_pop(p1), add = TRUE)
  DBI::dbExecute(p1$db_conn,
    "INSERT INTO genome_map VALUES (100, 1, 'Locus_1', 'X', NULL, 'default', 1.0)")  # sex 'X'
  expect_error(validate_genome_map(p1$db_conn), "invalid sex")

  p2 <- mk("gmv2"); on.exit(close_pop(p2), add = TRUE)
  DBI::dbExecute(p2$db_conn,
    "INSERT INTO genome_map VALUES (100, 1, 'Locus_1', 'M', NULL, NULL, 1.0)")        # NULL map_name
  expect_error(validate_genome_map(p2$db_conn), "invalid sex|map_name")

  p3 <- mk("gmv3"); on.exit(close_pop(p3), add = TRUE)
  DBI::dbExecute(p3$db_conn,
    "INSERT INTO genome_map VALUES (100, 1, 'Locus_1', 'M', NULL, 'default', 'NaN'::DOUBLE)")
  expect_error(validate_genome_map(p3$db_conn), "non-finite")

  p4 <- mk("gmv4"); on.exit(close_pop(p4), add = TRUE)
  # duplicate id_genome_map (reuse an existing id with a distinct logical key)
  DBI::dbExecute(p4$db_conn,
    "INSERT INTO genome_map VALUES (1, 1, 'Locus_1', 'M', NULL, 'default', 1.0)")
  expect_error(validate_genome_map(p4$db_conn), "id_genome_map")
})

test_that("genome_map is registered in schema metadata", {
  pop <- open_pop(pop_name = "gmmeta", db_name = ":memory:") |>
    define_genome(n_loci = 5, n_chr = 1, chr_len_Mb = 10)
  on.exit(close_pop(pop), add = TRUE)
  meta <- DBI::dbGetQuery(pop$db_conn,
    "SELECT column_name FROM _schema_meta WHERE table_name = 'genome_map' AND column_name IS NOT NULL")
  expect_true(all(c("id_genome_map", "locus_id", "sex", "line_name", "map_name", "pos_cM")
                  %in% meta$column_name))
})

test_that("add_offspring() uses the gamete-producing parent's resolved map (sex-specific)", {
  build <- function(seed, add_male_override) {
    set.seed(seed)
    pop <- open_pop(pop_name = "rm", db_name = ":memory:") |>
      define_genome(n_loci = 60, n_chr = 2, chr_len_Mb = c(30, 30)) |>
      define_founder_haplotypes(n_haplotypes = 30, method = "uniform") |>
      get_table("founder_haplotypes") |>
      add_founders(n_males = 4, n_females = 4, line_name = "A")
    if (add_male_override) {
      base <- resolve_genome_map(pop$db_conn)
      nid  <- next_int_id(pop$db_conn, "genome_map", "id_genome_map")
      df <- data.frame(id_genome_map = seq.int(nid, length.out = nrow(base)),
                       locus_id = base$locus_id, locus_name = base$locus_name,
                       sex = "M", line_name = NA_character_, map_name = "default",
                       pos_cM = base$pos_cM * 8)   # sires recombine ~8x more
      duckdb::duckdb_register(pop$db_conn, "__mo", df)
      DBI::dbExecute(pop$db_conn, "INSERT INTO genome_map SELECT * FROM __mo")
      duckdb::duckdb_unregister(pop$db_conn, "__mo")
      validate_genome_map(pop$db_conn)
    }
    ids <- DBI::dbGetQuery(pop$db_conn, "SELECT id_ind, sex FROM ind_meta")
    matings <- data.frame(
      id_parent_1 = ids$id_ind[ids$sex == "M"][1:4],
      id_parent_2 = ids$id_ind[ids$sex == "F"][1:4],
      sex = "M", line_name = "B", stringsAsFactors = FALSE)
    set.seed(seed + 1)
    pop <- add_offspring(pop, matings)
    h <- DBI::dbGetQuery(pop$db_conn,
      "SELECT allele FROM ind_haplotype WHERE id_ind LIKE 'B_%' AND parent_origin = 1
       ORDER BY id_ind, locus_id")$allele
    close_pop(pop); h
  }
  default_run  <- build(7, FALSE)
  override_run <- build(7, TRUE)
  # A male-specific genetic map must change sire-derived recombination.
  expect_false(identical(default_run, override_run))
})

test_that("validate_genome_map() catches duplicate key, orphan, and name mismatch", {
  pop <- open_pop(pop_name = "gm7", db_name = ":memory:") |>
    define_genome(n_loci = 4, n_chr = 1, chr_len_Mb = 10)
  on.exit(close_pop(pop), add = TRUE)
  con <- pop$db_conn

  # duplicate (locus_id, NULL, NULL, 'default')
  DBI::dbExecute(con,
    "INSERT INTO genome_map VALUES (100, 1, 'Locus_1', NULL, NULL, 'default', 1.0)")
  expect_error(validate_genome_map(con), "duplicate")

  # orphan locus_id (not in genome_meta)
  pop2 <- open_pop(pop_name = "gm7b", db_name = ":memory:") |>
    define_genome(n_loci = 4, n_chr = 1, chr_len_Mb = 10)
  on.exit(close_pop(pop2), add = TRUE)
  DBI::dbExecute(pop2$db_conn,
    "INSERT INTO genome_map VALUES (100, 999, 'ghost', NULL, NULL, 'altmap', 1.0)")
  expect_error(validate_genome_map(pop2$db_conn), "not present in genome_meta")

  # locus_id <-> locus_name mismatch
  pop3 <- open_pop(pop_name = "gm7c", db_name = ":memory:") |>
    define_genome(n_loci = 4, n_chr = 1, chr_len_Mb = 10)
  on.exit(close_pop(pop3), add = TRUE)
  DBI::dbExecute(pop3$db_conn,
    "INSERT INTO genome_map VALUES (100, 1, 'WRONG_NAME', 'M', NULL, 'default', 1.0)")
  expect_error(validate_genome_map(pop3$db_conn), "disagrees with genome_meta")
})

# ---------------------------------------------------------------------------
# Forward reproducibility + type preservation through the pipeline
# ---------------------------------------------------------------------------

.run_pipeline <- function(seed) {
  set.seed(seed)
  pop <- open_pop(pop_name = "rp", db_name = ":memory:") |>
    define_genome(n_loci = 40, n_chr = 2, chr_len_Mb = c(15, 10)) |>
    define_founder_haplotypes(n_haplotypes = 20, method = "mosaic") |>
    get_table("founder_haplotypes") |>
    add_founders(n_males = 4, n_females = 4, line_name = "A")
  ids <- DBI::dbGetQuery(pop$db_conn, "SELECT id_ind, sex FROM ind_meta")
  matings <- data.frame(
    id_parent_1 = ids$id_ind[ids$sex == "M"][c(1, 2, 1, 2)],
    id_parent_2 = ids$id_ind[ids$sex == "F"][c(1, 2, 3, 4)],
    sex = c("M", "F", "M", "F"), line_name = "B", stringsAsFactors = FALSE)
  pop <- add_offspring(pop, matings)
  hap <- DBI::dbGetQuery(pop$db_conn,
    "SELECT * FROM ind_haplotype ORDER BY id_ind, parent_origin, locus_id")
  pt <- DBI::dbGetQuery(pop$db_conn,
    "SELECT data_type FROM information_schema.columns
     WHERE table_name = 'genome_meta' AND column_name = 'pos_bp'")$data_type
  close_pop(pop)
  list(hap = hap, pos_bp_type = pt)
}

test_that("same seed reproduces founder+offspring haplotypes on the current code", {
  a <- .run_pipeline(2024)
  b <- .run_pipeline(2024)
  expect_identical(a$hap, b$hap)
})

test_that("pos_bp stays BIGINT through the full pipeline (no table-rewrite demotion)", {
  a <- .run_pipeline(7)
  expect_equal(a$pos_bp_type, "BIGINT")
})
