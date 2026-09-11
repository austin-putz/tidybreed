# Phase B of plans/update_genome_effects_v4.md (v4.9): the term/member/origin
# tables, their declared constraints, the views, and the cross-row validator.
#
# Gates 46-48 and 54. The evaluator is Phase D; nothing here computes a value.

ge_pop <- function(n_loci = 20, n_chr = 2) {
  open_pop(pop_name = "ge", db_name = ":memory:") |>
    define_genome(n_loci = n_loci, n_chr = n_chr, chr_len_Mb = 50)
}

# Insert one minimal valid term. Returns id_genome_effect.
ge_seed_term <- function(pop, id = 1L, trait = "ADG", owner = "custom",
                         value = 2.0, locus_id = 1L, contrast = "additive",
                         center = 0.5) {
  DBI::dbExecute(pop$db_conn, sprintf(
    "INSERT INTO genome_effects VALUES (%d, '%s', '%s', NULL, %f)",
    id, trait, owner, value))
  DBI::dbExecute(pop$db_conn, sprintf(
    "INSERT INTO genome_effect_members VALUES (%d, 1, %d, '%s', NULL, NULL, %f)",
    id, locus_id, contrast, center))
  invisible(id)
}


# ── Gate 47: genome_meta primary key and the FK it makes possible ───────────

test_that("genome_meta has a primary key and refuses a duplicate locus_id", {
  pop <- ge_pop()
  on.exit(close_pop(pop), add = TRUE)

  expect_error(
    DBI::dbExecute(pop$db_conn,
                   "INSERT INTO genome_meta VALUES (1, 'dup', 1, '1', 999)"),
    "[Pp]rimary key|Duplicate key"
  )
})

test_that("a member naming a nonexistent locus_id is refused by the FK", {
  pop <- ge_pop(n_loci = 20)
  on.exit(close_pop(pop), add = TRUE)
  ge_seed_term(pop)

  expect_error(
    DBI::dbExecute(pop$db_conn, paste0(
      "INSERT INTO genome_effect_members VALUES (1, 2, 9999, 'additive', ",
      "NULL, NULL, 0.5)")),
    "foreign key"
  )
})

test_that("genome_meta stays writable by ALTER + UPDATE under the new PK", {
  # define_chip(), define_founder_haplotypes() and mutate_table() all reach
  # genome_meta this way and never rewrite the table.
  pop <- ge_pop(n_loci = 30, n_chr = 3)
  on.exit(close_pop(pop), add = TRUE)

  pop <- pop |>
    get_table("genome_meta") |>
    dplyr::filter(chr %in% 1:2) |>
    define_chip("panel")

  got <- pop |> get_table("genome_meta") |> dplyr::collect()
  expect_true("is_panel" %in% names(got))
  expect_equal(sum(got$is_panel), sum(got$chr %in% 1:2))
})


# ── Gate 46 + row-local constraints: every one is declared in SQL ───────────

test_that("member row-local invariants are enforced by the database", {
  pop <- ge_pop()
  on.exit(close_pop(pop), add = TRUE)
  ge_seed_term(pop)
  ins <- function(sql) DBI::dbExecute(pop$db_conn, sql)

  # Gate 46: a bare CHECK (center BETWEEN 0 AND 1) accepts NULL, because SQL
  # accepts UNKNOWN. The branch must also say IS NOT NULL.
  expect_error(ins("INSERT INTO genome_effect_members VALUES (1,2,2,'additive',NULL,NULL,NULL)"),
               "CHECK")
  expect_error(ins("INSERT INTO genome_effect_members VALUES (1,2,2,'dominance',NULL,NULL,NULL)"),
               "CHECK")
  expect_error(ins("INSERT INTO genome_effect_members VALUES (1,2,2,'additive',NULL,NULL,1.5)"),
               "CHECK")
  expect_error(ins("INSERT INTO genome_effect_members VALUES (1,2,2,'nonsense',NULL,NULL,0.5)"),
               "CHECK")
  # An indicator state is the (copy_count, dosage) pair, and dosage cannot
  # exceed the copy count.
  expect_error(ins("INSERT INTO genome_effect_members VALUES (1,2,2,'indicator',NULL,0,NULL)"),
               "CHECK")
  expect_error(ins("INSERT INTO genome_effect_members VALUES (1,2,2,'indicator',1,2,NULL)"),
               "CHECK")
  expect_error(ins("INSERT INTO genome_effect_members VALUES (1,2,2,'indicator',2,1,0.5)"),
               "CHECK")
  # Valid ones.
  expect_silent(ins("INSERT INTO genome_effect_members VALUES (1,2,2,'indicator',2,1,NULL)"))
})

test_that("origin row-local invariants are enforced by the database", {
  pop <- ge_pop()
  on.exit(close_pop(pop), add = TRUE)
  ge_seed_term(pop)
  ins <- function(sql) DBI::dbExecute(pop$db_conn, sql)

  expect_error(ins("INSERT INTO genome_effect_member_origins VALUES (1,1,1,'nope',NULL,1,1)"),
               "CHECK")
  expect_error(ins("INSERT INTO genome_effect_member_origins VALUES (1,1,1,'exact',NULL,1,1)"),
               "CHECK")
  expect_error(ins("INSERT INTO genome_effect_member_origins VALUES (1,1,1,'unknown','A',1,1)"),
               "CHECK")
  expect_error(ins("INSERT INTO genome_effect_member_origins VALUES (1,1,1,'exact','A',1,0)"),
               "CHECK")
  expect_error(ins("INSERT INTO genome_effect_member_origins VALUES (1,1,1,'exact','A',3,1)"),
               "CHECK")
  # Gate 41 (row-local half): 'any' without a parent is the common scope
  # written the long way, and is refused by the CHECK, not only the validator.
  expect_error(ins("INSERT INTO genome_effect_member_origins VALUES (1,1,1,'any',NULL,NULL,1)"),
               "CHECK")
  expect_silent(ins("INSERT INTO genome_effect_member_origins VALUES (1,1,1,'any',NULL,1,1)"))
})

test_that("orphans inside the effect set are caught by the R validator", {
  # There is deliberately no FOREIGN KEY from members to effects or from
  # origins to members. DuckDB 1.5.5 refuses to delete a parent row inside an
  # explicit transaction whose children were deleted earlier in the same
  # transaction, which would make every replace mode of
  # define_genome_effects() unwritable atomically. The rule is enforced by
  # validate_genome_effects() instead, which runs before every COMMIT.
  pop <- ge_pop()
  on.exit(close_pop(pop), add = TRUE)
  ge_seed_term(pop)
  DBI::dbExecute(pop$db_conn,
                 "INSERT INTO genome_effect_member_origins VALUES (1,1,1,'exact','A',NULL,1)")
  expect_silent(tidybreed:::validate_genome_effects(pop$db_conn))

  # A term whose members are gone.
  DBI::dbExecute(pop$db_conn, "DELETE FROM genome_effect_member_origins")
  DBI::dbExecute(pop$db_conn, "DELETE FROM genome_effect_members")
  expect_error(tidybreed:::validate_genome_effects(pop$db_conn), "has no members")

  # A member with no term, and an origin row with no member.
  DBI::dbExecute(pop$db_conn, "DELETE FROM genome_effects")
  DBI::dbExecute(pop$db_conn, paste0(
    "INSERT INTO genome_effect_members VALUES (1, 1, 1, 'additive', NULL, NULL, 0.5)"))
  expect_error(tidybreed:::validate_genome_effects(pop$db_conn),
               "member rows with no term")
  DBI::dbExecute(pop$db_conn,
                 "INSERT INTO genome_effect_member_origins VALUES (1,9,1,'exact','A',NULL,1)")
  expect_error(tidybreed:::validate_genome_effects(pop$db_conn),
               "origin rows with no member")
})

test_that("a parent delete inside a transaction is what forced the FKs out", {
  # Pin the DuckDB behaviour the schema comment cites, so a future version that
  # fixes it is noticed rather than silently leaving the workaround in place.
  conn <- DBI::dbConnect(duckdb::duckdb())
  on.exit(DBI::dbDisconnect(conn, shutdown = TRUE), add = TRUE)
  DBI::dbExecute(conn, "CREATE TABLE p (id INTEGER PRIMARY KEY)")
  DBI::dbExecute(conn, paste0(
    "CREATE TABLE c (id INTEGER PRIMARY KEY, ",
    "FOREIGN KEY (id) REFERENCES p(id))"))
  DBI::dbExecute(conn, "INSERT INTO p VALUES (1)")
  DBI::dbExecute(conn, "INSERT INTO c VALUES (1)")

  DBI::dbExecute(conn, "BEGIN TRANSACTION")
  DBI::dbExecute(conn, "DELETE FROM c WHERE id = 1")
  expect_error(DBI::dbExecute(conn, "DELETE FROM p WHERE id = 1"),
               "foreign key")
  DBI::dbExecute(conn, "ROLLBACK")

  # The identical sequence in autocommit is fine, which is why this reads as a
  # transaction-visibility limitation and not as a modelling error.
  DBI::dbExecute(conn, "DELETE FROM c WHERE id = 1")
  expect_equal(DBI::dbExecute(conn, "DELETE FROM p WHERE id = 1"), 1L)
})

test_that("remove_rows() refuses the effect tables and says why", {
  pop <- ge_pop()
  on.exit(close_pop(pop), add = TRUE)
  ge_seed_term(pop)

  expect_error(
    pop |> get_table("genome_effects") |>
      dplyr::filter(id_genome_effect == 1L) |> remove_rows(),
    "define_genome_effects"
  )
  expect_error(
    pop |> get_table("genome_effect_members") |>
      dplyr::filter(id_genome_effect == 1L) |> remove_rows(),
    "define_genome_effects"
  )
})


# ── The validator: cross-row rules SQL cannot express ──────────────────────

test_that("a well-formed model validates clean", {
  pop <- ge_pop()
  on.exit(close_pop(pop), add = TRUE)
  ge_seed_term(pop, id = 1L, value = 2.0, center = 0.5)                 # common
  ge_seed_term(pop, id = 2L, value = 3.0, center = 0.4)                 # A-specific
  DBI::dbExecute(pop$db_conn,
                 "INSERT INTO genome_effect_member_origins VALUES (2,1,1,'exact','A',NULL,1)")

  expect_silent(validate_genome_effects(pop$db_conn))
})

test_that("duplicate family + scope identity is rejected", {
  pop <- ge_pop()
  on.exit(close_pop(pop), add = TRUE)
  ge_seed_term(pop, id = 1L, value = 2.0, center = 0.5)
  ge_seed_term(pop, id = 2L, value = 1.0, center = 0.5)

  expect_error(validate_genome_effects(pop$db_conn), "duplicate family")
})

test_that("mixed centring at one locus collides as a duplicate (Q2)", {
  # center_value is outside the family signature on purpose: a line-specific
  # variant legitimately carries a different frequency from its fallback.
  pop <- ge_pop()
  on.exit(close_pop(pop), add = TRUE)
  ge_seed_term(pop, id = 1L, value = 2.0, center = 0.5)
  ge_seed_term(pop, id = 2L, value = 1.0, center = 0.2)

  expect_error(validate_genome_effects(pop$db_conn), "duplicate family")
})

test_that("overlapping but incomparable scopes in one family are rejected", {
  pop <- ge_pop()
  on.exit(close_pop(pop), add = TRUE)
  ge_seed_term(pop, id = 1L)
  ge_seed_term(pop, id = 2L)
  DBI::dbExecute(pop$db_conn,
                 "INSERT INTO genome_effect_member_origins VALUES (1,1,1,'exact','A',NULL,1)")
  DBI::dbExecute(pop$db_conn,
                 "INSERT INTO genome_effect_member_origins VALUES (2,1,1,'any',NULL,1,1)")

  expect_error(validate_genome_effects(pop$db_conn), "incomparable")
})

test_that("'any' is refused on a genotype member (cross-table, so R)", {
  pop <- ge_pop()
  on.exit(close_pop(pop), add = TRUE)
  ge_seed_term(pop, id = 1L, contrast = "dominance")
  DBI::dbExecute(pop$db_conn,
                 "INSERT INTO genome_effect_member_origins VALUES (1,1,1,'any',NULL,1,1)")

  expect_error(validate_genome_effects(pop$db_conn), "only on additive members")
})

test_that("additive members take at most one origin row, with copy_count 1", {
  pop <- ge_pop()
  on.exit(close_pop(pop), add = TRUE)
  ge_seed_term(pop, id = 1L)
  DBI::dbExecute(pop$db_conn,
                 "INSERT INTO genome_effect_member_origins VALUES (1,1,1,'exact','A',NULL,1)")
  DBI::dbExecute(pop$db_conn,
                 "INSERT INTO genome_effect_member_origins VALUES (1,1,2,'exact','B',NULL,1)")

  expect_error(validate_genome_effects(pop$db_conn), "at most one origin row")

  pop2 <- ge_pop()
  on.exit(close_pop(pop2), add = TRUE)
  ge_seed_term(pop2, id = 1L)
  DBI::dbExecute(pop2$db_conn,
                 "INSERT INTO genome_effect_member_origins VALUES (1,1,1,'exact','A',NULL,2)")
  expect_error(validate_genome_effects(pop2$db_conn), "copy_count = 1")
})

test_that("members must be canonicalized by ascending locus_id", {
  pop <- ge_pop()
  on.exit(close_pop(pop), add = TRUE)
  DBI::dbExecute(pop$db_conn, "INSERT INTO genome_effects VALUES (1,'ADG','custom',NULL,1.0)")
  DBI::dbExecute(pop$db_conn,
                 "INSERT INTO genome_effect_members VALUES (1,1,5,'additive',NULL,NULL,0.5)")
  DBI::dbExecute(pop$db_conn,
                 "INSERT INTO genome_effect_members VALUES (1,2,3,'additive',NULL,NULL,0.5)")

  expect_error(validate_genome_effects(pop$db_conn), "ascending locus_id")
})


# ── Containment: the worked cases from the plan's lattices ─────────────────

test_that("additive containment matches the plan's worked cases", {
  mk <- function(...) {
    o <- list(...)
    origins <- if (!length(o)) {
      data.frame(member_slot = integer(0), origin_slot = integer(0),
                 line_match_type = character(0), line_name = character(0),
                 parent_origin = integer(0), copy_count = integer(0))
    } else {
      data.frame(member_slot = 1L, origin_slot = 1L,
                 line_match_type = o[[1]], line_name = o[[2]],
                 parent_origin = o[[3]], copy_count = 1L,
                 stringsAsFactors = FALSE)
    }
    .ge_predicate(
      data.frame(member_slot = 1L, locus_id = 1L, contrast_name = "additive",
                 stringsAsFactors = FALSE),
      origins)
  }
  common  <- mk()
  exactA  <- mk("exact", "A", NA_integer_)
  exactA1 <- mk("exact", "A", 1L)
  exactB2 <- mk("exact", "B", 2L)
  anyP1   <- mk("any", NA_character_, 1L)
  unknown <- mk("unknown", NA_character_, NA_integer_)

  expect_true(.ge_pred_leq(exactA, common))     # (ANY,ANY) contains (exact A, ANY)
  expect_false(.ge_pred_leq(common, exactA))
  expect_true(.ge_pred_leq(exactA1, exactA))    # the tie v4.1 got wrong
  expect_false(.ge_pred_leq(exactA, exactA1))
  expect_false(.ge_pred_leq(exactA1, exactB2))  # disjoint: both apply
  expect_false(.ge_pred_overlap(exactA1, exactB2))
  expect_false(.ge_pred_leq(exactA, anyP1))     # overlap, neither contains
  expect_false(.ge_pred_leq(anyP1, exactA))
  expect_true(.ge_pred_overlap(exactA, anyP1))
  expect_true(.ge_pred_leq(anyP1, common))
  expect_true(.ge_pred_leq(unknown, common))
  expect_false(.ge_pred_overlap(unknown, exactA))
})

test_that("genotype containment handles multisets and reciprocals", {
  mk <- function(rows) {
    origins <- if (is.null(rows)) {
      data.frame(member_slot = integer(0), origin_slot = integer(0),
                 line_match_type = character(0), line_name = character(0),
                 parent_origin = integer(0), copy_count = integer(0))
    } else rows
    .ge_predicate(
      data.frame(member_slot = 1L, locus_id = 1L, contrast_name = "dominance",
                 stringsAsFactors = FALSE),
      origins)
  }
  row <- function(slot, line, parent, count = 1L) {
    data.frame(member_slot = 1L, origin_slot = slot, line_match_type = "exact",
               line_name = line, parent_origin = parent, copy_count = count,
               stringsAsFactors = FALSE)
  }
  common <- mk(NULL)
  ab     <- mk(rbind(row(1L, "A", NA_integer_), row(2L, "B", NA_integer_)))
  recip1 <- mk(rbind(row(1L, "A", 1L),          row(2L, "B", 2L)))
  recip2 <- mk(rbind(row(1L, "A", 2L),          row(2L, "B", 1L)))
  aa     <- mk(row(1L, "A", NA_integer_, 2L))

  expect_true(.ge_pred_leq(ab, common))
  expect_true(.ge_pred_leq(recip1, ab))         # reciprocal override works
  expect_false(.ge_pred_leq(ab, recip1))
  expect_false(.ge_pred_leq(recip1, recip2))    # disjoint, both valid
  expect_false(.ge_pred_overlap(recip1, recip2))
  expect_false(.ge_pred_leq(aa, ab))            # {A:2} vs {A:1,B:1}
  expect_false(.ge_pred_overlap(aa, ab))
})

test_that("multiset matching is exhaustive, not greedy", {
  # {A@p1, A@ANY} is satisfiable by {A@p1, A@p2}: a greedy pass that lets the
  # ANY row eat the p1 copy first would wrongly call it unsatisfiable. Found by
  # the Phase A fixtures.
  demands <- data.frame(line = c("A", "A"), parent = c(1L, NA_integer_),
                        stringsAsFactors = FALSE)
  copies  <- data.frame(line = c("A", "A"), parent = c(2L, 1L),
                        stringsAsFactors = FALSE)
  expect_true(.ge_bijection(demands, copies, mode = "compatible"))
  expect_true(.ge_bijection(copies, demands, mode = "refines"))
})


# ── Views ──────────────────────────────────────────────────────────────────

test_that("genome_effect_terms derives order, signature and family_key", {
  pop <- ge_pop()
  on.exit(close_pop(pop), add = TRUE)
  ge_seed_term(pop, id = 1L, value = 2.0, center = 0.5)   # common additive @1
  ge_seed_term(pop, id = 2L, value = 3.0, center = 0.4)   # A-specific additive @1
  DBI::dbExecute(pop$db_conn,
                 "INSERT INTO genome_effect_member_origins VALUES (2,1,1,'exact','A',NULL,1)")
  ge_seed_term(pop, id = 3L, value = 4.0, contrast = "dominance", center = 0.5)
  DBI::dbExecute(pop$db_conn, "INSERT INTO genome_effects VALUES (4,'ADG','custom','axa',1.0)")
  DBI::dbExecute(pop$db_conn,
                 "INSERT INTO genome_effect_members VALUES (4,1,1,'additive',NULL,NULL,0.5)")
  DBI::dbExecute(pop$db_conn,
                 "INSERT INTO genome_effect_members VALUES (4,2,2,'additive',NULL,NULL,0.5)")

  v <- pop |> get_table("genome_effect_terms") |> dplyr::collect()
  v <- v[order(v$id_genome_effect), ]

  expect_equal(v$effect_order, c(1L, 1L, 1L, 2L))
  expect_equal(v$contrast_signature,
               c("additive", "additive", "dominance", "additive+additive"))
  expect_equal(v$scope_description,
               c("common", "1:exact(A)x1", "common", "common"))

  # Gate 54: scope variants share a family and compete; different contrasts and
  # different member sets are different families and sum.
  expect_equal(v$family_key[1], v$family_key[2])
  expect_false(v$family_key[1] == v$family_key[3])
  expect_false(v$family_key[1] == v$family_key[4])
})

test_that("the view's family_key agrees with the R family signature", {
  pop <- ge_pop()
  on.exit(close_pop(pop), add = TRUE)
  ge_seed_term(pop, id = 1L)
  DBI::dbExecute(pop$db_conn, "INSERT INTO genome_effects VALUES (2,'ADG','custom',NULL,1.0)")
  DBI::dbExecute(pop$db_conn,
                 "INSERT INTO genome_effect_members VALUES (2,1,1,'indicator',2,1,NULL)")

  v <- pop |> get_table("genome_effect_terms") |> dplyr::collect()
  members <- DBI::dbGetQuery(pop$db_conn, "SELECT * FROM genome_effect_members")
  expect_equal(v$family_key, .ge_family_keys(v, members))
})

test_that("genome_effect_loci joins locus_name back in", {
  pop <- ge_pop()
  on.exit(close_pop(pop), add = TRUE)
  ge_seed_term(pop, id = 1L, locus_id = 7L)

  got <- pop |> get_table("genome_effect_loci") |> dplyr::collect()
  expect_equal(got$locus_id, 7L)
  expect_equal(got$locus_name, "Locus_7")
  expect_equal(got$contrast_name, "additive")
  expect_false("locus_name" %in% DBI::dbListFields(pop$db_conn, "genome_effect_members"))
})

test_that("ind_tgv_total sums components and is never stored", {
  pop <- ge_pop()
  on.exit(close_pop(pop), add = TRUE)
  DBI::dbExecute(pop$db_conn,
                 "INSERT INTO ind_tgv VALUES (1,'A_1','ADG','order1_additive',2.0)")
  DBI::dbExecute(pop$db_conn,
                 "INSERT INTO ind_tgv VALUES (2,'A_1','ADG','order1_dominance',-0.5)")

  tot <- pop |> get_table("ind_tgv_total") |> dplyr::collect()
  expect_equal(nrow(tot), 1L)
  expect_equal(tot$tgv_total, 1.5)
  expect_false("replicate" %in% DBI::dbListFields(pop$db_conn, "ind_tgv"))
})


# ── Gate 48: fresh vs restored parity ──────────────────────────────────────

test_that("a fresh population and a restored one list the same objects", {
  dir <- withr::local_tempdir()
  f   <- file.path(dir, "ge.duckdb")

  pop <- open_pop(pop_name = "GE", db_name = f) |>
    define_genome(n_loci = 20, n_chr = 2, chr_len_Mb = 50)
  built <- schema(pop, show_empty = TRUE, include_system = TRUE)$table_name
  close_pop(pop)

  pop2 <- restore_pop(f)
  on.exit(close_pop(pop2), add = TRUE)
  restored <- schema(pop2, show_empty = TRUE, include_system = TRUE)$table_name

  expect_identical(built, restored)
  expect_true(all(c("genome_effect_terms", "genome_effect_loci", "ind_tgv_total")
                  %in% built))
  # Every listed object carries a description; a view that fell through the
  # registries would print "(no description)".
  s <- schema(pop2, show_empty = TRUE, include_system = TRUE)
  expect_false(any(is.na(s$description) | s$description == ""))
})


# ── define_genome() preflight ──────────────────────────────────────────────

test_that("define_genome() still refuses to run twice", {
  pop <- ge_pop()
  on.exit(close_pop(pop), add = TRUE)
  expect_error(define_genome(pop, n_loci = 10, n_chr = 1, chr_len_Mb = 50),
               "Genome already defined")
})

test_that("open_pop() alone creates no genome-effect table", {
  # The two edits are one commit: while open_pop() still created
  # genome_effects, listing it in GENOME_TABLES fired the preflight on every
  # fresh population.
  pop <- open_pop(pop_name = "bare", db_name = ":memory:")
  on.exit(close_pop(pop), add = TRUE)
  live <- DBI::dbListTables(pop$db_conn)

  expect_false(any(c("genome_effects", "genome_effect_members",
                     "genome_effect_member_origins") %in% live))
  expect_true("ind_tgv" %in% live)   # a result table, created with ind_tbv
  expect_message(define_genome(pop, n_loci = 10, n_chr = 1, chr_len_Mb = 50),
                 "Defined genome")
})


# ── Dominance is refused where the locus is not reliably diploid ───────────

test_that("dominance is refused at a non-diploid locus, and provable on the member", {
  pop <- ge_pop(n_loci = 20, n_chr = 2)
  on.exit(close_pop(pop), add = TRUE)

  # Chromosome 2 becomes a male-hemizygous X.
  pop <- define_chromosome(pop, "2", offspring_sex = "M",
                           from_parent_1 = 0, from_parent_2 = 1)
  x_locus <- DBI::dbGetQuery(
    pop$db_conn, "SELECT locus_id FROM genome_meta WHERE chr_name = '2' LIMIT 1")$locus_id

  ge_seed_term(pop, id = 1L, contrast = "dominance", locus_id = x_locus, center = 0.5)
  expect_error(validate_genome_effects(pop$db_conn), "needs a diploid locus")

  # An exact multiset demanding two copies on that same member proves diploidy.
  DBI::dbExecute(pop$db_conn,
                 "INSERT INTO genome_effect_member_origins VALUES (1,1,1,'exact','A',NULL,2)")
  expect_silent(validate_genome_effects(pop$db_conn))
})

test_that("a diploid-proving scope on a different member proves nothing", {
  pop <- ge_pop(n_loci = 20, n_chr = 2)
  on.exit(close_pop(pop), add = TRUE)
  pop <- define_chromosome(pop, "2", offspring_sex = "M",
                           from_parent_1 = 0, from_parent_2 = 1)
  auto <- DBI::dbGetQuery(
    pop$db_conn, "SELECT locus_id FROM genome_meta WHERE chr_name = '1' LIMIT 1")$locus_id
  x_locus <- DBI::dbGetQuery(
    pop$db_conn, "SELECT locus_id FROM genome_meta WHERE chr_name = '2' LIMIT 1")$locus_id

  DBI::dbExecute(pop$db_conn, "INSERT INTO genome_effects VALUES (1,'ADG','custom',NULL,1.0)")
  DBI::dbExecute(pop$db_conn, sprintf(
    "INSERT INTO genome_effect_members VALUES (1,1,%d,'dominance',NULL,NULL,0.5)", auto))
  DBI::dbExecute(pop$db_conn, sprintf(
    "INSERT INTO genome_effect_members VALUES (1,2,%d,'dominance',NULL,NULL,0.5)", x_locus))
  # Two-copy scope on member 1 (the autosome) only.
  DBI::dbExecute(pop$db_conn,
                 "INSERT INTO genome_effect_member_origins VALUES (1,1,1,'exact','A',NULL,2)")

  expect_error(validate_genome_effects(pop$db_conn), "member 2")
})


# ── Regressions found reviewing Phase B against the plan ──────────────────

test_that("'unknown' inside a genotype multiset compares like any other line", {
  # The line dimension is a token, not a raw line_name: comparing raw names made
  # every predicate containing an 'unknown' demand incomparable with itself, so
  # a duplicate scope went undetected. 'unknown' is legal on a genotype member
  # (it constrains: it matches copies whose founding line is NULL); only 'any'
  # is not.
  pop <- ge_pop()
  on.exit(close_pop(pop), add = TRUE)
  seed <- function(id, value) {
    DBI::dbExecute(pop$db_conn, sprintf(
      "INSERT INTO genome_effects VALUES (%d,'ADG','custom',NULL,%f)", id, value))
    DBI::dbExecute(pop$db_conn, sprintf(
      "INSERT INTO genome_effect_members VALUES (%d,1,1,'dominance',NULL,NULL,0.5)", id))
    DBI::dbExecute(pop$db_conn, sprintf(
      "INSERT INTO genome_effect_member_origins VALUES (%d,1,1,'exact','A',NULL,1)", id))
    DBI::dbExecute(pop$db_conn, sprintf(
      "INSERT INTO genome_effect_member_origins VALUES (%d,1,2,'unknown',NULL,NULL,1)", id))
  }
  seed(1L, 2.0)
  expect_silent(validate_genome_effects(pop$db_conn))
  seed(2L, 5.0)
  expect_error(validate_genome_effects(pop$db_conn), "duplicate family")
})

test_that("a line literally named 'unknown' is not the unknown-line scope", {
  pop <- ge_pop()
  on.exit(close_pop(pop), add = TRUE)
  ge_seed_term(pop, id = 1L)
  ge_seed_term(pop, id = 2L)
  DBI::dbExecute(pop$db_conn,
                 "INSERT INTO genome_effect_member_origins VALUES (1,1,1,'exact','unknown',NULL,1)")
  DBI::dbExecute(pop$db_conn,
                 "INSERT INTO genome_effect_member_origins VALUES (2,1,1,'unknown',NULL,NULL,1)")

  # Disjoint, not duplicate: one matches a line called "unknown", the other
  # matches copies with no line at all.
  expect_silent(validate_genome_effects(pop$db_conn))
})

test_that("a genotype origin multiset must account for every copy of the state", {
  pop <- ge_pop()
  on.exit(close_pop(pop), add = TRUE)

  ge_seed_term(pop, id = 1L, contrast = "dominance", center = 0.5)
  DBI::dbExecute(pop$db_conn,
                 "INSERT INTO genome_effect_member_origins VALUES (1,1,1,'exact','A',NULL,1)")
  expect_error(validate_genome_effects(pop$db_conn), "defined over 2")

  DBI::dbExecute(pop$db_conn,
                 "INSERT INTO genome_effect_member_origins VALUES (1,1,2,'exact','B',NULL,1)")
  expect_silent(validate_genome_effects(pop$db_conn))
})

test_that("an indicator multiset must sum to its declared copy_count_value", {
  pop <- ge_pop()
  on.exit(close_pop(pop), add = TRUE)
  DBI::dbExecute(pop$db_conn, "INSERT INTO genome_effects VALUES (1,'ADG','custom',NULL,3.0)")
  DBI::dbExecute(pop$db_conn,
                 "INSERT INTO genome_effect_members VALUES (1,1,1,'indicator',1,0,NULL)")
  DBI::dbExecute(pop$db_conn,
                 "INSERT INTO genome_effect_member_origins VALUES (1,1,1,'exact','A',NULL,2)")

  expect_error(validate_genome_effects(pop$db_conn), "defined over 1")
})

test_that("member rows with no term are reported, not skipped", {
  # The empty-model early return used to swallow this: with no term rows there
  # was nothing to iterate, so orphans passed silently.
  expect_true(any(grepl(
    "no term",
    .ge_validate_frames(
      terms = data.frame(id_genome_effect = integer(0), trait_name = character(0),
                         effect_owner = character(0), effect_name = character(0),
                         genome_value = numeric(0)),
      members = data.frame(id_genome_effect = 1L, member_slot = 1L, locus_id = 1L,
                           contrast_name = "additive",
                           copy_count_value = NA_integer_, dosage_value = NA_integer_,
                           center_value = 0.5),
      origins = data.frame(id_genome_effect = integer(0), member_slot = integer(0),
                           origin_slot = integer(0), line_match_type = character(0),
                           line_name = character(0), parent_origin = integer(0),
                           copy_count = integer(0))))))
})


# ── Views are registered like any other system object ──────────────────────

test_that("the views are registered in the row-key and reserved-column registries", {
  # Regression: adding the views to SYSTEM_TABLES without a TABLE_ROW_KEYS or
  # TABLE_NO_ROW_DELETE entry broke the registry-completeness guard in
  # test-schema-registries.R, which exists precisely to force that decision.
  for (vw in c("genome_effect_terms", "genome_effect_loci", "ind_tgv_total")) {
    expect_true(vw %in% names(TABLE_NO_ROW_DELETE), info = vw)
    expect_true(vw %in% names(TABLE_RESERVED_COLS), info = vw)
  }
})

test_that("a view lists exactly the columns it reserves, and refuses writes", {
  pop <- ge_pop()
  on.exit(close_pop(pop), add = TRUE)

  for (vw in c("genome_effect_terms", "genome_effect_loci", "ind_tgv_total")) {
    expect_setequal(DBI::dbListFields(pop$db_conn, vw), TABLE_RESERVED_COLS[[vw]])
  }
  expect_error(
    pop |> get_table("genome_effect_terms") |> mutate_table(family_key = "x"),
    "reserved"
  )
  # remove_rows() checks for an empty/unfiltered selection before it consults
  # TABLE_NO_ROW_DELETE, so the view needs a row for the refusal to be reached.
  DBI::dbExecute(pop$db_conn,
                 "INSERT INTO ind_tgv VALUES (1,'A_1','ADG','order1_additive',2.0)")
  expect_error(
    pop |> get_table("ind_tgv_total") |>
      dplyr::filter(trait_name == "ADG") |> remove_rows(),
    "derived view"
  )
})
