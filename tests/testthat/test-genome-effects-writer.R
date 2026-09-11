# Phase C of plans/update_genome_effects_v4.md (v4.9): define_genome_effects(),
# the terms/origin input format, the (a, d) and genotype-table builders, and
# define_additive_effects() rebuilt on top of them.
#
# Gates 34-35, 41-44, 50 and 53. Nothing here computes a genetic value; the
# evaluator is Phase D. What is asserted is what reaches storage.

gew_pop <- function(n_loci = 20, n_chr = 2, traits = "ADG", var = 0.25) {
  pop <- open_pop(pop_name = "gew", db_name = ":memory:") |>
    define_genome(n_loci = n_loci, n_chr = n_chr, chr_len_Mb = 50)
  for (t in traits) pop <- define_trait(pop, t, target_add_var = var)
  pop
}

# A genome shaped like the Phase A fixture context: L1-L4 autosomal, LX on a
# mammalian X, LY on a Y. The filler loci exist only because define_genome()
# spreads loci evenly across chromosomes.
gew_fixture_pop <- function() {
  pop <- open_pop(pop_name = "gewfix", db_name = ":memory:") |>
    define_genome(n_loci = 8, n_chr = 4, chr_len_Mb = 20,
                  locus_names = c("L1", "L2", "L3", "L4",
                                  "LX", "LXb", "LY", "LYb"),
                  chr_names   = c("A1", "A2", "X", "Y"))
  pop <- pop |>
    define_chromosome("X", offspring_sex = "M", from_parent_1 = 0, from_parent_2 = 1) |>
    define_chromosome("Y", offspring_sex = "M", from_parent_1 = 1, from_parent_2 = 0) |>
    define_chromosome("Y", offspring_sex = "F", from_parent_1 = 0, from_parent_2 = 0) |>
    define_chromosome("Y", recombines = FALSE)
  define_trait(pop, "ADG", target_add_var = 1.0)
}

gew_lines_pop <- function(n_loci = 6) {
  pop <- open_pop(pop_name = "gewl", db_name = ":memory:") |>
    define_genome(n_loci = n_loci, n_chr = 1, chr_len_Mb = 20)
  set.seed(11)
  pop <- define_founder_haplotypes(pop, n_haplotypes = 20, line_name = "A")
  pop <- define_founder_haplotypes(pop, n_haplotypes = 20, line_name = "B")
  define_trait(pop, "ADG", target_add_var = 1.0)
}

gew_scopes <- function(pop, trait = "ADG") {
  DBI::dbGetQuery(pop$db_conn, paste0(
    "SELECT t.id_genome_effect, t.effect_owner, t.genome_value, ",
    "t.scope_description FROM genome_effect_terms t ",
    "WHERE t.trait_name = '", trait, "' ORDER BY t.id_genome_effect"))
}


# ── Gate 53: the terms format round-trips ──────────────────────────────────

test_that("worked example 1 (one dominance term) round-trips", {
  pop <- gew_pop()
  on.exit(close_pop(pop), add = TRUE)

  pop <- define_genome_effects(pop, "ADG", terms = data.frame(
    locus_name = "Locus_10", contrast_name = "dominance",
    center_value = 0.3, genome_value = 0.8))

  got <- DBI::dbGetQuery(pop$db_conn, paste0(
    "SELECT l.locus_name, m.contrast_name, m.center_value, m.copy_count_value, ",
    "m.dosage_value, t.genome_value, t.effect_owner, t.effect_order, ",
    "t.scope_description FROM genome_effect_terms t ",
    "JOIN genome_effect_loci l USING (id_genome_effect) ",
    "JOIN genome_effect_members m USING (id_genome_effect, member_slot)"))
  expect_equal(nrow(got), 1L)
  expect_equal(got$locus_name, "Locus_10")
  expect_equal(got$contrast_name, "dominance")
  expect_equal(got$center_value, 0.3)
  expect_true(is.na(got$copy_count_value) && is.na(got$dosage_value))
  expect_equal(got$genome_value, 0.8)
  expect_equal(got$effect_owner, "custom")
  expect_equal(got$effect_order, 1L)
  expect_equal(got$scope_description, "common")
  # No origin rows at all is the common scope, not an empty scope.
  expect_equal(DBI::dbGetQuery(pop$db_conn,
    "SELECT COUNT(*) n FROM genome_effect_member_origins")$n, 0)
})

test_that("worked example 2 (a 3x3 A x A surface) round-trips, sparse", {
  pop <- gew_pop()
  on.exit(close_pop(pop), add = TRUE)

  cells <- expand.grid(g1 = 0:2, g2 = 0:2)
  cells$value <- c(0, 0, 0,  0, 1.4, 2.1,  0, 2.1, 3.6)
  surface <- rbind(
    data.frame(term_id = seq_len(nrow(cells)), locus_name = "Locus_10",
               contrast_name = "indicator", dosage_value = cells$g1,
               genome_value  = cells$value),
    data.frame(term_id = seq_len(nrow(cells)), locus_name = "Locus_14",
               contrast_name = "indicator", dosage_value = cells$g2,
               genome_value  = cells$value))
  pop <- define_genome_effects(pop, "ADG", surface[surface$genome_value != 0, ],
                               effect_owner = "epistasis_AxA")

  got <- gew_scopes(pop)
  expect_equal(nrow(got), 4L)                   # the four non-zero cells
  expect_equal(sort(got$genome_value), c(1.4, 2.1, 2.1, 3.6))
  expect_true(all(got$effect_owner == "epistasis_AxA"))

  # copy_count_value is inferred at an ordinary autosome: the state is the pair,
  # but there is only one copy count it can be, so the user never types it.
  cc <- DBI::dbGetQuery(pop$db_conn,
    "SELECT DISTINCT copy_count_value FROM genome_effect_members")
  expect_equal(cc$copy_count_value, 2L)

  # Every term is two members, so each is an interaction, and the four cells are
  # four families -- they sum, they do not compete.
  ord <- DBI::dbGetQuery(pop$db_conn,
    "SELECT effect_order, family_key FROM genome_effect_terms")
  expect_true(all(ord$effect_order == 2L))
  expect_equal(length(unique(ord$family_key)), 4L)
})

test_that("worked example 3 (reciprocal dominance) round-trips", {
  pop <- gew_pop()
  on.exit(close_pop(pop), add = TRUE)

  pop <- define_genome_effects(
    pop, "ADG",
    terms = data.frame(term_id = 1L, locus_name = "Locus_10",
                       contrast_name = "dominance", center_value = 0.3,
                       genome_value = 1.2),
    origin = data.frame(term_id = 1L, locus_name = "Locus_10",
                        line_match_type = "exact",
                        line_name     = c("Duroc", "Landrace"),
                        parent_origin = c(1L, 2L),
                        copy_count    = c(1L, 1L)),
    effect_owner = "reciprocal")

  o <- DBI::dbGetQuery(pop$db_conn, paste0(
    "SELECT * FROM genome_effect_member_origins ORDER BY origin_slot"))
  expect_equal(nrow(o), 2L)
  expect_equal(o$origin_slot, c(1L, 2L))
  expect_equal(o$line_name, c("Duroc", "Landrace"))
  expect_equal(o$parent_origin, c(1L, 2L))
  expect_equal(sum(o$copy_count), 2L)           # sums to the dominance state

  # The mirror image is a second, disjoint variant: both apply, to different
  # individuals, so they never compete.
  pop <- define_genome_effects(
    pop, "ADG",
    terms = data.frame(term_id = "mirror", locus_name = "Locus_10",
                       contrast_name = "dominance", center_value = 0.3,
                       genome_value = 0.4),
    origin = data.frame(term_id = "mirror", locus_name = "Locus_10",
                        line_match_type = "exact",
                        line_name     = c("Duroc", "Landrace"),
                        parent_origin = c(2L, 1L),
                        copy_count    = c(1L, 1L)),
    effect_owner = "reciprocal")
  fam <- DBI::dbGetQuery(pop$db_conn,
    "SELECT family_key FROM genome_effect_terms WHERE effect_owner = 'reciprocal'")
  expect_equal(length(unique(fam$family_key)), 1L)   # same family
  expect_equal(nrow(fam), 2L)                        # both stored
})

test_that("origin rows are canonicalized, so input order does not matter", {
  mk <- function(rev) {
    pop <- gew_pop()
    o <- data.frame(term_id = 1L, locus_name = "Locus_10",
                    line_match_type = "exact",
                    line_name = c("Duroc", "Landrace"),
                    parent_origin = c(1L, 2L), copy_count = c(1L, 1L))
    if (rev) o <- o[rev(seq_len(nrow(o))), ]
    pop <- define_genome_effects(
      pop, "ADG",
      terms = data.frame(term_id = 1L, locus_name = "Locus_10",
                         contrast_name = "dominance", center_value = 0.3,
                         genome_value = 1.2),
      origin = o)
    out <- DBI::dbGetQuery(pop$db_conn,
      "SELECT * FROM genome_effect_member_origins ORDER BY origin_slot")
    close_pop(pop)
    out
  }
  expect_equal(mk(FALSE), mk(TRUE))
})

test_that("members are canonicalized by ascending locus_id whatever the input order", {
  pop <- gew_pop()
  on.exit(close_pop(pop), add = TRUE)
  pop <- define_genome_effects(pop, "ADG", data.frame(
    term_id = 1L, locus_name = c("Locus_14", "Locus_3"),
    contrast_name = "additive", center_value = 0.5, genome_value = 2.5))

  m <- DBI::dbGetQuery(pop$db_conn, paste0(
    "SELECT m.member_slot, l.locus_name FROM genome_effect_members m ",
    "JOIN genome_effect_loci l USING (id_genome_effect, member_slot) ",
    "ORDER BY m.member_slot"))
  expect_equal(m$member_slot, c(1L, 2L))
  expect_equal(m$locus_name, c("Locus_3", "Locus_14"))
})

test_that("every valid Phase A fixture round-trips through the writer", {
  # The fixtures were hand-computed before any DDL existed and are the semantic
  # reference for the whole feature. Writing each one back through the public
  # entry point is what turns "representable" into "writable".
  pop <- gew_fixture_pop()
  on.exit(close_pop(pop), add = TRUE)

  # expect_invalid fixtures are rejected by the validator; expect_eval_error
  # ones are rejected by the *writer* now that a chr_inheritance table exists --
  # Phase A had none, so F14 could only be pinned at the evaluator.
  fx <- gefx_fixtures()
  valid <- Filter(function(f) is.null(f$expect_invalid) &&
                    is.null(f$expect_eval_error), fx)
  expect_gt(length(valid), 10L)

  for (nm in names(valid)) {
    model <- valid[[nm]]$model
    tt <- merge(model$members, model$terms[, c("id_genome_effect", "genome_value")],
                by = "id_genome_effect")
    terms <- data.frame(
      term_id          = tt$id_genome_effect,
      locus_name       = tt$locus_name,
      contrast_name    = tt$contrast_name,
      center_value     = tt$center_value,
      copy_count_value = tt$copy_count_value,
      dosage_value     = tt$dosage_value,
      genome_value     = tt$genome_value,
      stringsAsFactors = FALSE)
    org <- model$origins
    origin <- if (nrow(org) == 0L) NULL else {
      slots <- merge(org, model$members[, c("id_genome_effect", "member_slot",
                                            "locus_name")],
                     by = c("id_genome_effect", "member_slot"))
      data.frame(term_id = slots$id_genome_effect, locus_name = slots$locus_name,
                 line_match_type = slots$line_match_type,
                 line_name = slots$line_name,
                 parent_origin = slots$parent_origin,
                 copy_count = slots$copy_count, stringsAsFactors = FALSE)
    }

    expect_message(
      pop <- define_genome_effects(pop, "ADG", terms, origin = origin,
                                   effect_owner = "custom",
                                   mode = "replace_owner"),
      "Wrote", label = nm)

    # Read back and compare to the fixture's own canonical form, keyed by the
    # family the term lands in and its scope -- not by a surrogate id.
    back <- DBI::dbGetQuery(pop$db_conn, paste0(
      "SELECT t.family_key, t.scope_description, t.genome_value ",
      "FROM genome_effect_terms t ORDER BY t.family_key, t.scope_description"))
    want <- do.call(rbind, lapply(model$terms$id_genome_effect, function(id) {
      mm <- model$members[model$members$id_genome_effect == id, , drop = FALSE]
      data.frame(n_members = nrow(mm),
                 value = model$terms$genome_value[
                   model$terms$id_genome_effect == id])
    }))
    expect_equal(nrow(back), nrow(model$terms), label = nm)
    expect_equal(sort(back$genome_value), sort(want$value), label = nm)
    expect_equal(length(unique(back$family_key)),
                 length(gefx_families(model, "ADG")), label = nm)
  }

  # F14 -- dominance at a hemizygous locus -- was pinned at the evaluator in
  # Phase A for want of a chr_inheritance table. It is now refused at the
  # writer, which is where the fixture's own note said it belongs.
  expect_error(
    define_genome_effects(pop, "ADG", data.frame(
      locus_name = "LX", contrast_name = "dominance", center_value = 0.5,
      genome_value = 3.0), mode = "replace_owner"),
    "needs a diploid locus")
})

test_that("a malformed terms call names the user's term_id, not a surrogate id", {
  pop <- gew_pop()
  on.exit(close_pop(pop), add = TRUE)

  # Varying genome_value within one term.
  expect_error(define_genome_effects(pop, "ADG", data.frame(
    term_id = "AxB", locus_name = c("Locus_1", "Locus_2"),
    contrast_name = "additive", center_value = 0.5,
    genome_value = c(1, 2))), "term_id 'AxB'.*genome_value")

  # The same locus twice in one term.
  expect_error(define_genome_effects(pop, "ADG", data.frame(
    term_id = "dup", locus_name = c("Locus_1", "Locus_1"),
    contrast_name = "additive", center_value = 0.5,
    genome_value = 1)), "term_id 'dup'.*more than once")

  # An unknown column is rejected, never silently dropped.
  expect_error(define_genome_effects(pop, "ADG", data.frame(
    locus_name = "Locus_1", contrast_name = "additive", center_value = 0.5,
    genome_value = 1, base_allele_freq = 0.3)),
    "Unknown column in 'terms': 'base_allele_freq'")

  # effect_name varying within a term.
  expect_error(define_genome_effects(pop, "ADG", data.frame(
    term_id = "lbl", locus_name = c("Locus_1", "Locus_2"),
    contrast_name = "additive", center_value = 0.5, genome_value = 1,
    effect_name = c("x", "y"))), "term_id 'lbl'.*effect_name")

  # A contrast/state mismatch is caught before any SQL, naming the locus.
  expect_error(define_genome_effects(pop, "ADG", data.frame(
    term_id = "nc", locus_name = "Locus_1", contrast_name = "additive",
    genome_value = 1)), "term_id 'nc' locus 'Locus_1'.*center_value")
  expect_error(define_genome_effects(pop, "ADG", data.frame(
    term_id = "ic", locus_name = "Locus_1", contrast_name = "indicator",
    dosage_value = 1, center_value = 0.5, genome_value = 1)),
    "must not carry 'center_value'")

  # And a locus that is not in the genome.
  expect_error(define_genome_effects(pop, "ADG", data.frame(
    locus_name = "nope", contrast_name = "additive", center_value = 0.5,
    genome_value = 1)), "not in genome_meta")

  # Nothing was written by any of them.
  expect_equal(DBI::dbGetQuery(pop$db_conn,
    "SELECT COUNT(*) n FROM genome_effects")$n, 0)
})

test_that("a duplicate family + scope identity is rejected and rolls back", {
  pop <- gew_pop()
  on.exit(close_pop(pop), add = TRUE)
  pop <- define_genome_effects(pop, "ADG", data.frame(
    locus_name = "Locus_1", contrast_name = "additive", center_value = 0.5,
    genome_value = 1))

  expect_error(define_genome_effects(pop, "ADG", data.frame(
    term_id = "again", locus_name = "Locus_1", contrast_name = "additive",
    center_value = 0.5, genome_value = 9)), "same logical term at the same scope")
  # The failed call left the first term untouched and added nothing.
  got <- gew_scopes(pop)
  expect_equal(nrow(got), 1L)
  expect_equal(got$genome_value, 1)
})


# ── require_complete, and copy-count inference ─────────────────────────────

test_that("require_complete rejects a sparse surface and names what is missing", {
  pop <- gew_pop()
  on.exit(close_pop(pop), add = TRUE)
  cells <- expand.grid(g1 = 0:2, g2 = 0:2)
  cells$value <- c(0, 0, 0,  0, 1.4, 2.1,  0, 2.1, 3.6)
  surface <- rbind(
    data.frame(term_id = seq_len(nrow(cells)), locus_name = "Locus_10",
               contrast_name = "indicator", dosage_value = cells$g1,
               genome_value  = cells$value),
    data.frame(term_id = seq_len(nrow(cells)), locus_name = "Locus_14",
               contrast_name = "indicator", dosage_value = cells$g2,
               genome_value  = cells$value))

  expect_error(
    define_genome_effects(pop, "ADG", surface[surface$genome_value != 0, ],
                          require_complete = TRUE),
    "missing 5 of 9 reachable")
  # The complete surface is accepted, zeros and all.
  expect_message(
    define_genome_effects(pop, "ADG", surface, require_complete = TRUE),
    "Wrote 9")
})

test_that("require_complete counts the absent state where a chromosome can be absent", {
  pop <- gew_fixture_pop()
  on.exit(close_pop(pop), add = TRUE)

  # Y is (1,0) in males and (0,0) in females, so the reachable states are
  # copy_count 0 and 1: three states, not two, and one of them has no
  # haplotype row to join at evaluation time.
  states <- data.frame(term_id = 1:2, locus_name = "LY",
                       contrast_name = "indicator",
                       copy_count_value = c(1L, 1L), dosage_value = c(0L, 1L),
                       genome_value = c(1.0, 2.0))
  expect_error(define_genome_effects(pop, "ADG", states, require_complete = TRUE),
               "missing 1 of 3 reachable")

  full <- rbind(states, data.frame(term_id = 3L, locus_name = "LY",
                                   contrast_name = "indicator",
                                   copy_count_value = 0L, dosage_value = 0L,
                                   genome_value = 3.0))
  expect_message(define_genome_effects(pop, "ADG", full, require_complete = TRUE),
                 "Wrote 3")
})

test_that("copy_count_value cannot be inferred at a variable-copy locus", {
  pop <- gew_fixture_pop()
  on.exit(close_pop(pop), add = TRUE)

  # The X carries 1 copy in males and 2 in females, so dosage 0 is two different
  # biological situations and the writer must refuse to guess.
  expect_error(define_genome_effects(pop, "ADG", data.frame(
    locus_name = "LX", contrast_name = "indicator", dosage_value = 0L,
    genome_value = 1)), "cannot be inferred.*1 or 2 copies")

  expect_message(define_genome_effects(pop, "ADG", data.frame(
    locus_name = "LX", contrast_name = "indicator", copy_count_value = 1L,
    dosage_value = 0L, genome_value = 1)), "Wrote 1")
})


# ── Gate 41: 'any' constraints ─────────────────────────────────────────────

test_that("'any' without a parent_origin is refused by the writer and by SQL", {
  pop <- gew_pop()
  on.exit(close_pop(pop), add = TRUE)

  # Through the writer: a scope that constrains nothing is not the same thing
  # as the common scope, and asking for it is a mistake, not a synonym. The
  # message arrives before any SQL, naming the user's own term.
  expect_error(define_genome_effects(pop, "ADG", data.frame(
    term_id = "t1", locus_name = "Locus_1", contrast_name = "additive",
    center_value = 0.5, genome_value = 1),
    origin = list(line_match_type = "any", line_name = "A")),
    "term_id 't1' locus 'Locus_1'.*takes no line_name")
  expect_error(define_genome_effects(pop, "ADG", data.frame(
    term_id = "t2", locus_name = "Locus_1", contrast_name = "additive",
    center_value = 0.5, genome_value = 1),
    origin = data.frame(term_id = "t2", locus_name = "Locus_1",
                        line_match_type = "any", copy_count = 1L)),
    "constrains nothing and is not a second spelling")
  expect_equal(DBI::dbGetQuery(pop$db_conn,
    "SELECT COUNT(*) n FROM genome_effects")$n, 0)

  # And by direct insert, because it is also a row-local CHECK -- the R check
  # is for the message, the constraint is what makes it true of the table.
  DBI::dbExecute(pop$db_conn,
    "INSERT INTO genome_effects VALUES (1, 'ADG', 'custom', NULL, 2.0)")
  DBI::dbExecute(pop$db_conn, paste0(
    "INSERT INTO genome_effect_members VALUES (1, 1, 1, 'additive', ",
    "NULL, NULL, 0.5)"))
  expect_error(DBI::dbExecute(pop$db_conn,
    "INSERT INTO genome_effect_member_origins VALUES (1,1,1,'any',NULL,NULL,1)"),
    "CHECK")
})

test_that("'any' on a genotype member is refused (cross-table, so R-only)", {
  pop <- gew_pop()
  on.exit(close_pop(pop), add = TRUE)

  expect_error(define_genome_effects(pop, "ADG", data.frame(
    locus_name = "Locus_1", contrast_name = "dominance", center_value = 0.4,
    genome_value = 1), origin = list(line_match_type = "any",
                                     parent_origin = 1, copy_count = 2)),
    "'any' is permitted only on additive members")
})

test_that("a scoped genotype member must sum to its state's copy count", {
  pop <- gew_pop()
  on.exit(close_pop(pop), add = TRUE)

  # The scalar-list form gives one origin row, which is all an additive member
  # may have and is not enough for a dominance member.
  expect_error(define_genome_effects(pop, "ADG", data.frame(
    locus_name = "Locus_1", contrast_name = "dominance", center_value = 0.4,
    genome_value = 1), origin = list(line_name = "A")),
    "demands 1 copies but the member's state is defined over 2")
})


# ── Gates 34, 42: owners and the wrapper scope matrix ──────────────────────

test_that("the reserved owner is refused to the general writer", {
  pop <- gew_pop()
  on.exit(close_pop(pop), add = TRUE)
  expect_error(define_genome_effects(pop, "ADG", data.frame(
    locus_name = "Locus_1", contrast_name = "additive", center_value = 0.5,
    genome_value = 1), effect_owner = "generated_additive_tbv"),
    "reserved effect owner")
})

test_that("gate 34: rerunning define_additive_effects() cannot delete custom terms", {
  pop <- gew_lines_pop()
  on.exit(close_pop(pop), add = TRUE)

  pop <- define_genome_effects(pop, "ADG", data.frame(
    locus_name = "Locus_1", contrast_name = "dominance", center_value = 0.4,
    genome_value = 5), effect_owner = "custom")
  pop <- define_genome_effects(pop, "ADG", data.frame(
    locus_name = "Locus_2", contrast_name = "additive", center_value = 0.5,
    genome_value = 7), effect_owner = "my_model")

  for (i in 1:3) {
    pop <- pop |> get_table("genome_meta") |>
      define_additive_effects("ADG", effects = rep(i, 6),
                              base_line_name = "A")
  }

  owners <- DBI::dbGetQuery(pop$db_conn, paste0(
    "SELECT effect_owner, COUNT(*) n FROM genome_effect_terms ",
    "GROUP BY 1 ORDER BY 1"))
  expect_equal(owners$effect_owner,
               c("custom", "generated_additive_tbv", "my_model"))
  expect_equal(owners$n, c(1, 6, 1))
  # And the custom values are untouched.
  expect_equal(DBI::dbGetQuery(pop$db_conn, paste0(
    "SELECT genome_value FROM genome_effect_terms ",
    "WHERE effect_owner = 'my_model'"))$genome_value, 7)
})

test_that("replace_trait refuses to take the reserved owner with it", {
  pop <- gew_lines_pop()
  on.exit(close_pop(pop), add = TRUE)
  pop <- pop |> get_table("genome_meta") |>
    define_additive_effects("ADG", effects = rep(1, 6), base_line_name = "A")

  expect_error(define_genome_effects(pop, "ADG", data.frame(
    locus_name = "Locus_1", contrast_name = "dominance", center_value = 0.4,
    genome_value = 5), mode = "replace_trait"),
    "reserved owner 'generated_additive_tbv'")
  expect_equal(DBI::dbGetQuery(pop$db_conn,
    "SELECT COUNT(*) n FROM genome_effects")$n, 6)
})

test_that("gate 42: all four wrapper scopes round-trip, and replace_scope isolates them", {
  pop <- gew_lines_pop()
  on.exit(close_pop(pop), add = TRUE)

  # Each scope gets its own locus, so each lands in its own fallback family.
  # They cannot share one: ('exact' A, parent ANY) and ('any', parent 1)
  # overlap without either containing the other -- a line-A paternal copy
  # matches both -- and the validator refuses that pair by design. Family
  # separation is what lets all four coexist here; replace_scope's isolation is
  # per scope within (trait, owner), not per locus, so the test is still real.
  eff <- function(v, locus, ...) {
    pop <<- pop |> get_table("genome_meta") |>
      dplyr::filter(locus_name == !!locus) |>
      define_additive_effects("ADG", effects = v, base_line_name = "A", ...)
  }
  eff(1, "Locus_1")                                      # NULL / NULL
  eff(2, "Locus_2", line_name = "A")                     # "A"  / NULL
  eff(3, "Locus_3", parent_origin = 1)                   # NULL / 1
  eff(4, "Locus_4", line_name = "A", parent_origin = 2)  # "A"  / 2

  scopes <- DBI::dbGetQuery(pop$db_conn, paste0(
    "SELECT genome_value, scope_description FROM genome_effect_terms ",
    "ORDER BY genome_value"))
  expect_equal(scopes$genome_value, c(1, 2, 3, 4))
  expect_equal(scopes$scope_description,
               c("common", "1:exact(A)x1", "1:any@p1x1", "1:exact(A)@p2x1"))

  # Every scope carries exactly the origin rows the mapping table promises: one
  # row in every scoped case, none at all for the common one.
  n_org <- DBI::dbGetQuery(pop$db_conn, paste0(
    "SELECT e.genome_value, COUNT(o.origin_slot) n FROM genome_effects e ",
    "LEFT JOIN genome_effect_member_origins o USING (id_genome_effect) ",
    "GROUP BY 1 ORDER BY 1"))
  expect_equal(n_org$n, c(0, 1, 1, 1))

  # Replacing the ('exact' A, parent ANY) scope leaves the other three exactly
  # as they were, including the one that differs from it only by a parent.
  eff(99, "Locus_5", line_name = "A")
  after <- DBI::dbGetQuery(pop$db_conn, paste0(
    "SELECT genome_value, scope_description FROM genome_effect_terms ",
    "ORDER BY genome_value"))
  expect_equal(after$genome_value, c(1, 3, 4, 99))
  expect_equal(after$scope_description,
               c("common", "1:any@p1x1", "1:exact(A)@p2x1", "1:exact(A)x1"))
})

test_that("two scopes that overlap without nesting are refused, whichever order", {
  # ('exact' A, parent ANY) vs ('any', parent 1): a line-A paternal copy matches
  # both and neither is more specific, so no variant could be selected for it.
  for (order in list(c(NA, 1L), c(1L, NA))) {
    pop <- gew_lines_pop()
    first  <- if (is.na(order[1])) list(line_name = "A") else list(parent_origin = order[1])
    second <- if (is.na(order[2])) list(line_name = "A") else list(parent_origin = order[2])
    pop <- do.call(define_additive_effects, c(
      list(pop |> get_table("genome_meta"), "ADG", effects = rep(1, 6),
           base_line_name = "A"), first))
    expect_error(do.call(define_additive_effects, c(
      list(pop |> get_table("genome_meta"), "ADG", effects = rep(2, 6),
           base_line_name = "A"), second)),
      "overlapping but incomparable")
    expect_equal(DBI::dbGetQuery(pop$db_conn,
      "SELECT COUNT(*) n FROM genome_effects")$n, 6)
    close_pop(pop)
  }
})

test_that("gate 35: successive common, A and B calls all stand, in one family", {
  pop <- gew_lines_pop()
  on.exit(close_pop(pop), add = TRUE)
  for (a in list(list(1, NULL), list(2, "A"), list(3, "B"))) {
    pop <- if (is.null(a[[2]])) {
      pop |> get_table("genome_meta") |>
        define_additive_effects("ADG", effects = rep(a[[1]], 6),
                                base_line_name = "A")
    } else {
      pop |> get_table("genome_meta") |>
        define_additive_effects("ADG", effects = rep(a[[1]], 6),
                                line_name = a[[2]])
    }
  }
  # The two views share trait_name, effect_owner and genome_value, so a join
  # between them has to qualify every shared column.
  got <- DBI::dbGetQuery(pop$db_conn, paste0(
    "SELECT l.locus_name, t.scope_description, t.genome_value ",
    "FROM genome_effect_terms t ",
    "JOIN genome_effect_loci l USING (id_genome_effect) ",
    "WHERE l.locus_name = 'Locus_1' ORDER BY t.genome_value"))
  expect_equal(got$genome_value, c(1, 2, 3))
  expect_equal(got$scope_description,
               c("common", "1:exact(A)x1", "1:exact(B)x1"))

  # All three compete rather than sum: same family, different scopes. That is
  # the whole point -- separate owners would sum and give no common fallback.
  fam <- DBI::dbGetQuery(pop$db_conn, paste0(
    "SELECT DISTINCT t.family_key FROM genome_effect_terms t ",
    "JOIN genome_effect_loci l USING (id_genome_effect) ",
    "WHERE l.locus_name = 'Locus_1'"))
  expect_equal(nrow(fam), 1L)

  # Re-running the common call drops loci absent from the new set -- the whole
  # scope is replaced, not merged into.
  pop <- pop |> get_table("genome_meta") |> dplyr::filter(locus_id <= 3) |>
    define_additive_effects("ADG", effects = rep(8, 3), base_line_name = "A")
  by_scope <- DBI::dbGetQuery(pop$db_conn, paste0(
    "SELECT scope_description, COUNT(*) n FROM genome_effect_terms ",
    "GROUP BY 1 ORDER BY 1"))
  expect_equal(by_scope$scope_description,
               c("1:exact(A)x1", "1:exact(B)x1", "common"))
  expect_equal(by_scope$n, c(6, 6, 3))
})


# ── Gate 50: the parent-only re-run warning ────────────────────────────────

test_that("gate 50: a parent-only re-run warns, and the fallback pair still stands", {
  pop <- gew_lines_pop()
  on.exit(close_pop(pop), add = TRUE)
  pop <- pop |> get_table("genome_meta") |>
    define_additive_effects("ADG", effects = rep(2, 6), line_name = "A")

  expect_warning(
    pop <- pop |> get_table("genome_meta") |>
      define_additive_effects("ADG", effects = rep(5, 6), line_name = "A",
                              parent_origin = 1),
    "differ only in the parent dimension")

  got <- DBI::dbGetQuery(pop$db_conn, paste0(
    "SELECT DISTINCT genome_value, scope_description FROM genome_effect_terms ",
    "ORDER BY genome_value"))
  expect_equal(got$genome_value, c(2, 5))
  expect_equal(got$scope_description, c("1:exact(A)x1", "1:exact(A)@p1x1"))
})

test_that("gate 50: the common/A/B sequence and a reciprocal pair stay silent", {
  pop <- gew_lines_pop()
  on.exit(close_pop(pop), add = TRUE)

  # Members differ in *line*, not in parent: never the confusable case.
  pop <- pop |> get_table("genome_meta") |>
    define_additive_effects("ADG", effects = rep(1, 6), base_line_name = "A")
  expect_warning(
    pop <- pop |> get_table("genome_meta") |>
      define_additive_effects("ADG", effects = rep(2, 6), line_name = "A"),
    NA)
  expect_warning(
    pop <- pop |> get_table("genome_meta") |>
      define_additive_effects("ADG", effects = rep(3, 6), line_name = "B"),
    NA)

  # A reciprocal pair has disjoint parents rather than nested ones, so neither
  # contains the other and no fallback is implied.
  pop2 <- gew_lines_pop()
  on.exit(close_pop(pop2), add = TRUE)
  pop2 <- pop2 |> get_table("genome_meta") |>
    define_additive_effects("ADG", effects = rep(1, 6), line_name = "A",
                            parent_origin = 1)
  expect_warning(
    pop2 <- pop2 |> get_table("genome_meta") |>
      define_additive_effects("ADG", effects = rep(2, 6), line_name = "A",
                              parent_origin = 2),
    NA)
  expect_equal(DBI::dbGetQuery(pop2$db_conn,
    "SELECT COUNT(*) n FROM genome_effects")$n, 12)
})


# ── Gates 43, 44: variance target and the multi-trait origin contract ──────

test_that("gate 43: scale_to_target lands on V for parent-qualified effects too", {
  realized <- function(po) {
    pop <- gew_lines_pop(n_loci = 30)
    on.exit(close_pop(pop), add = TRUE)
    set.seed(909)
    pop <- suppressWarnings(
      pop |> get_table("genome_meta") |>
        define_additive_effects("ADG", distribution = "normal",
                                parent_origin = po, base_line_name = "A"))
    p <- tidybreed:::compute_base_allele_freq(pop, "founder_haplotypes", NULL, "A")
    a <- DBI::dbGetQuery(pop$db_conn, paste0(
      "SELECT m.locus_id, e.genome_value FROM genome_effects e ",
      "JOIN genome_effect_members m USING (id_genome_effect) ORDER BY m.locus_id"))
    n_elig <- if (is.null(po)) 2 else 1
    pq <- p[a$locus_id] * (1 - p[a$locus_id])
    sum(n_elig * pq * a$genome_value^2)
  }
  # The unparented model and the parent-qualified one both land at the target.
  # Before this, a parent-qualified model asked for V landed at V/2, because the
  # 2 in 2pq counts copies the term does not read.
  expect_equal(realized(NULL), 1.0, tolerance = 1e-8)
  expect_equal(realized(1L),   1.0, tolerance = 1e-8)
})

test_that("gate 44: parent_origin resolves per trait, in all three input forms", {
  mk <- function(po) {
    pop <- gew_lines_pop()
    pop <- define_trait(pop, "BW", target_add_var = 1.0)
    pop <- suppressWarnings(
      pop |> get_table("genome_meta") |>
        define_additive_effects(c("ADG", "BW"), effects = NULL,
                                G = diag(2), parent_origin = po,
                                base_line_name = "A"))
    out <- DBI::dbGetQuery(pop$db_conn, paste0(
      "SELECT DISTINCT e.trait_name, t.scope_description ",
      "FROM genome_effect_terms t JOIN genome_effects e USING (id_genome_effect) ",
      "ORDER BY e.trait_name"))
    close_pop(pop)
    out
  }
  # Scalar, recycled.
  expect_equal(mk(1L)$scope_description, c("1:any@p1x1", "1:any@p1x1"))
  # Positional vector, and named -- both uniform, since a mixed call with G is
  # rejected below.
  expect_equal(mk(c(2L, 2L))$scope_description, c("1:any@p2x1", "1:any@p2x1"))
  expect_equal(mk(c(ADG = 2L, BW = 2L))$scope_description,
               c("1:any@p2x1", "1:any@p2x1"))
})

test_that("gate 44: a mixed-origin correlated call is rejected, with the reason", {
  pop <- gew_lines_pop()
  on.exit(close_pop(pop), add = TRUE)
  pop <- define_trait(pop, "BW", target_add_var = 1.0)

  expect_error(
    pop |> get_table("genome_meta") |>
      define_additive_effects(c("ADG", "BW"), G = diag(2),
                              parent_origin = c(ADG = 1L, BW = 2L)),
    "covariance between a paternal-only and a maternal-only trait is zero")
  # Nothing was written.
  expect_equal(DBI::dbGetQuery(pop$db_conn,
    "SELECT COUNT(*) n FROM genome_effects")$n, 0)
})

test_that("parent_origin rejects values that are not 1 or 2", {
  pop <- gew_lines_pop()
  on.exit(close_pop(pop), add = TRUE)
  expect_error(pop |> get_table("genome_meta") |>
    define_additive_effects("ADG", effects = rep(1, 6), parent_origin = 3),
    "must be 1 .*or 2")
  expect_error(pop |> get_table("genome_meta") |>
    define_additive_effects("ADG", effects = rep(1, 6), parent_origin = "sire"),
    "must be 1 .*2 .*or NULL")
})


# ── The terms builders ─────────────────────────────────────────────────────

test_that("ad_terms() expands functional (a, d) and reports the mean without writing it", {
  expect_message(
    tt <- ad_terms("Locus_10", a = 0.4, d = 0.2, p = 0.3),
    "written to no table")
  expect_equal(nrow(tt), 2L)
  expect_equal(tt$contrast_name, c("additive", "indicator"))
  expect_equal(tt$center_value[1], 0.5)          # functional additive
  expect_equal(tt$copy_count_value[2], 2L)       # functional dominance
  expect_equal(tt$dosage_value[2], 1L)           # is 1[g = 1]
  expect_equal(tt$genome_value, c(0.4, 0.2))

  # mu = a(p - q) + 2pq d
  mu <- 0.4 * (0.3 - 0.7) + 2 * 0.3 * 0.7 * 0.2
  expect_message(ad_terms("Locus_10", a = 0.4, d = 0.2, p = 0.3),
                 format(round(mu, 6), nsmall = 0), fixed = TRUE)

  pop <- gew_pop()
  on.exit(close_pop(pop), add = TRUE)
  pop <- define_genome_effects(pop, "ADG", suppressMessages(
    ad_terms("Locus_10", a = 0.4, d = 0.2, p = 0.3)), effect_owner = "func")
  # Two terms, two families: an additive main effect and a dominance main
  # effect at one locus sum, they do not compete.
  fam <- DBI::dbGetQuery(pop$db_conn,
    "SELECT DISTINCT family_key FROM genome_effect_terms")
  expect_equal(nrow(fam), 2L)

  # Nothing reached phenotype_meta.
  expect_equal(DBI::dbGetQuery(pop$db_conn,
    "SELECT COUNT(*) n FROM phenotype_meta")$n, 0)
})

test_that("ad_terms() Cockerham coding centres on p and uses the dominance contrast", {
  tt <- suppressMessages(ad_terms("Locus_10", a = 0.4, d = 0.2, p = 0.3,
                                  coding = "cockerham"))
  expect_equal(tt$contrast_name, c("additive", "dominance"))
  expect_equal(tt$center_value, c(0.3, 0.3))
  expect_true(all(is.na(tt$copy_count_value)))
  # A centred model has mean zero by construction, and says so.
  expect_message(ad_terms("Locus_10", a = 0.4, d = 0.2, p = 0.3,
                          coding = "cockerham"), "mu = 0")
})

test_that("ad_terms() drops a zero coefficient rather than writing an inert term", {
  tt <- suppressMessages(ad_terms("Locus_10", a = 0.4, d = 0, p = 0.3))
  expect_equal(nrow(tt), 1L)
  expect_equal(tt$contrast_name, "additive")
  expect_error(suppressMessages(ad_terms("Locus_10", a = 0, d = 0, p = 0.3)),
               "no term to write")
  expect_error(ad_terms("Locus_10", a = 0.4, d = 0.2), "'p' .*is required")
})

test_that("genotype_terms() turns a genotype table into indicator terms", {
  cells <- expand.grid(Locus_10 = 0:2, Locus_14 = 0:2)
  vals  <- c(0, 0, 0,  0, 1.4, 2.1,  0, 2.1, 3.6)
  tt <- genotype_terms(cells, vals)
  expect_equal(nrow(tt), 8L)                     # 4 non-zero cells x 2 members
  expect_true(all(tt$contrast_name == "indicator"))
  expect_equal(length(unique(tt$term_id)), 4L)

  expect_equal(nrow(genotype_terms(cells, vals, drop_zero = FALSE)), 18L)
  expect_error(genotype_terms(cells, vals[1:3]), "one non-missing entry per row")
  expect_error(genotype_terms(cells, rep(0, 9)), "no term to write")

  pop <- gew_pop()
  on.exit(close_pop(pop), add = TRUE)
  expect_message(define_genome_effects(pop, "ADG", tt,
                                       effect_owner = "surface"), "Wrote 4")
})

test_that("genotype_terms() carries copy_count through for a variable-copy locus", {
  cells <- data.frame(LX = c(0L, 1L))
  tt <- genotype_terms(cells, c(1.5, 2.5), copy_count = list(LX = 1L))
  expect_equal(tt$copy_count_value, c(1L, 1L))

  pop <- gew_fixture_pop()
  on.exit(close_pop(pop), add = TRUE)
  expect_message(define_genome_effects(pop, "ADG", tt), "Wrote 2")
  expect_error(genotype_terms(cells, c(1, 2), copy_count = list(LZ = 1L)),
               "not a column of 'genotypes'")
})


# ── Replacement modes on the general writer ────────────────────────────────

test_that("replace_owner and replace_trait clear what they say they clear", {
  pop <- gew_pop()
  on.exit(close_pop(pop), add = TRUE)
  add <- function(owner, locus, v) {
    pop <<- define_genome_effects(pop, "ADG", data.frame(
      locus_name = locus, contrast_name = "additive", center_value = 0.5,
      genome_value = v), effect_owner = owner)
  }
  add("one", "Locus_1", 1); add("one", "Locus_2", 2); add("two", "Locus_3", 3)

  pop <- define_genome_effects(pop, "ADG", data.frame(
    locus_name = "Locus_5", contrast_name = "additive", center_value = 0.5,
    genome_value = 9), effect_owner = "one", mode = "replace_owner")
  got <- gew_scopes(pop)
  expect_equal(sort(got$genome_value), c(3, 9))

  # Children go with their parents: no orphan member or origin rows are left.
  expect_equal(DBI::dbGetQuery(pop$db_conn, paste0(
    "SELECT COUNT(*) n FROM genome_effect_members m WHERE NOT EXISTS (",
    "SELECT 1 FROM genome_effects e ",
    "WHERE e.id_genome_effect = m.id_genome_effect)"))$n, 0)

  pop <- define_genome_effects(pop, "ADG", data.frame(
    locus_name = "Locus_7", contrast_name = "additive", center_value = 0.5,
    genome_value = 4), effect_owner = "three", mode = "replace_trait")
  expect_equal(gew_scopes(pop)$genome_value, 4)
})

test_that("replace_scope refuses a per-member origin data frame", {
  pop <- gew_pop()
  on.exit(close_pop(pop), add = TRUE)
  expect_error(define_genome_effects(
    pop, "ADG",
    terms = data.frame(term_id = 1L, locus_name = "Locus_1",
                       contrast_name = "dominance", center_value = 0.4,
                       genome_value = 1),
    origin = data.frame(term_id = 1L, locus_name = "Locus_1",
                        line_match_type = "exact", line_name = c("A", "B"),
                        copy_count = c(1L, 1L)),
    mode = "replace_scope"),
    "no single scope")
})

test_that("a rejected write leaves the database exactly as it was", {
  pop <- gew_pop()
  on.exit(close_pop(pop), add = TRUE)
  pop <- define_genome_effects(pop, "ADG", data.frame(
    locus_name = "Locus_1", contrast_name = "additive", center_value = 0.5,
    genome_value = 1), effect_owner = "keep")
  before <- lapply(c("genome_effects", "genome_effect_members",
                     "genome_effect_member_origins"),
                   function(t) DBI::dbGetQuery(pop$db_conn,
                                               paste0("SELECT * FROM ", t)))

  # An overlapping-but-incomparable pair: neither scope is more specific, so no
  # variant could ever be selected. Rejected at write time, inside the
  # transaction that also deleted the old rows.
  expect_error(define_genome_effects(
    pop, "ADG",
    terms = data.frame(term_id = c("a", "b"), locus_name = "Locus_1",
                       contrast_name = "additive", center_value = 0.5,
                       genome_value = c(2, 3)),
    origin = data.frame(term_id = c("a", "b"), locus_name = "Locus_1",
                        line_match_type = c("exact", "any"),
                        line_name = c("A", NA), parent_origin = c(NA, 1L),
                        copy_count = 1L),
    effect_owner = "keep", mode = "replace_owner"),
    "overlapping but incomparable")

  after <- lapply(c("genome_effects", "genome_effect_members",
                    "genome_effect_member_origins"),
                  function(t) DBI::dbGetQuery(pop$db_conn,
                                              paste0("SELECT * FROM ", t)))
  expect_equal(after, before)
})

test_that("the writer needs define_genome() and an existing trait", {
  pop <- open_pop(pop_name = "bare", db_name = ":memory:")
  on.exit(close_pop(pop), add = TRUE)
  expect_error(define_genome_effects(pop, "ADG", data.frame(
    locus_name = "Locus_1", contrast_name = "additive", center_value = 0.5,
    genome_value = 1)), "genome-effect tables do not exist")

  pop2 <- gew_pop()
  on.exit(close_pop(pop2), add = TRUE)
  expect_error(define_genome_effects(pop2, "NOPE", data.frame(
    locus_name = "Locus_1", contrast_name = "additive", center_value = 0.5,
    genome_value = 1)), "not found in trait_meta")
})

test_that("a dominance member at a non-diploid locus is refused", {
  pop <- gew_fixture_pop()
  on.exit(close_pop(pop), add = TRUE)
  expect_error(define_genome_effects(pop, "ADG", data.frame(
    locus_name = "LX", contrast_name = "dominance", center_value = 0.4,
    genome_value = 1)), "needs a diploid locus")

  # The escape is an exact multiset on that same member demanding two copies.
  expect_message(define_genome_effects(
    pop, "ADG",
    terms = data.frame(term_id = 1L, locus_name = "LX",
                       contrast_name = "dominance", center_value = 0.4,
                       genome_value = 1),
    origin = data.frame(term_id = 1L, locus_name = "LX",
                        line_match_type = "exact", line_name = c("A", "B"),
                        copy_count = c(1L, 1L))), "Wrote 1")
})

test_that("a duplicate caused only by a different centring names the combined term", {
  # center_value is outside the family signature by design, so functional
  # additive@0.5 and Cockerham additive@p at one locus collide. They are
  # genuinely combinable, so the rejection says what single term to write.
  pop <- gew_pop()
  on.exit(close_pop(pop), add = TRUE)
  pop <- define_genome_effects(pop, "ADG", data.frame(
    locus_name = "Locus_1", contrast_name = "additive", center_value = 0.5,
    genome_value = 2))

  expect_error(define_genome_effects(pop, "ADG", data.frame(
    term_id = "cockerham", locus_name = "Locus_1", contrast_name = "additive",
    center_value = 0.3, genome_value = 6)),
    "genome_value = 8, center_value = 0.35")

  # Cancelling coefficients leave a constant, which this model has no place for.
  expect_error(define_genome_effects(pop, "ADG", data.frame(
    term_id = "cancel", locus_name = "Locus_1", contrast_name = "additive",
    center_value = 0.3, genome_value = -2)),
    "coefficients cancel")

  # A different centre at a *different scope* is not a duplicate at all: that
  # is exactly the line-specific case the exclusion exists for.
  pop2 <- gew_pop()
  on.exit(close_pop(pop2), add = TRUE)
  pop2 <- define_genome_effects(pop2, "ADG", data.frame(
    locus_name = "Locus_1", contrast_name = "additive", center_value = 0.5,
    genome_value = 2))
  expect_message(define_genome_effects(pop2, "ADG", data.frame(
    locus_name = "Locus_1", contrast_name = "additive", center_value = 0.3,
    genome_value = 6), origin = list(line_name = "A")), "Wrote 1")
})

test_that("silent coercions are refused rather than stored as plausible values", {
  pop <- gew_pop()
  on.exit(close_pop(pop), add = TRUE)
  # as.numeric() would turn this into a silent NA.
  expect_error(define_genome_effects(pop, "ADG", data.frame(
    locus_name = "Locus_1", contrast_name = "additive", center_value = 0.5,
    genome_value = "0.4")), "'genome_value' must be numeric")
  # as.integer() would truncate this to a different genotype state.
  expect_error(define_genome_effects(pop, "ADG", data.frame(
    locus_name = "Locus_1", contrast_name = "indicator", dosage_value = 1.5,
    genome_value = 1)), "must be a non-negative whole number")
  expect_equal(DBI::dbGetQuery(pop$db_conn,
    "SELECT COUNT(*) n FROM genome_effects")$n, 0)
})
