# The error-order contract of the genome-effect writer (decision D7 of
# plans/import_qtl_effect_methods_phase_5_plan.md).
#
# Written against the 0.75.2 code, before the step-5a vectorisation, and kept
# unchanged through it. When several things are wrong, the writer reports the
# lexicographically first violation by (term position in user order, rule
# position, member position) -- not the first term that fails whichever rule a
# vectorised version happens to check first. Validators that return a vector
# of violations keep the whole vector, in order.
#
# The expectations are hand-derived from the rules, not captured output.

gwo_pop <- function() {
  pop <- open_pop(pop_name = "gwo", db_name = ":memory:") |>
    define_genome(n_loci = 20, n_chr = 2, chr_len_Mb = 50)
  with_additive_target(pop, "ADG", 1)
}

gwo_error <- function(expr) {
  tryCatch({ force(expr); NA_character_ },
           error = function(e) conditionMessage(e))
}


# ── .ge_build(): term, then rule, then member ───────────────────────────────

test_that("two terms broken by the same rule: the first in user order is named", {
  pop <- gwo_pop()
  on.exit(close_pop(pop), add = TRUE)
  # term_id 'zeta' comes first in the input although it sorts last.
  terms <- data.frame(
    term_id       = c("zeta", "zeta", "alpha", "alpha"),
    locus_name    = c("Locus_1", "Locus_2", "Locus_3", "Locus_4"),
    contrast_name = "additive", center_value = 0.5,
    genome_value  = c(1, 2, 3, 4))
  msg <- gwo_error(define_genome_effect_terms(pop, "ADG", terms = terms))
  expect_identical(msg, paste0(
    "term_id 'zeta': 'genome_value' must be one non-missing value for the ",
    "whole term (got 1, 2). A term is one coefficient over its members."))
})

test_that("a later rule failing on an earlier term beats an earlier rule on a later term", {
  pop <- gwo_pop()
  on.exit(close_pop(pop), add = TRUE)
  # Term 'first' repeats a locus (rule 3); term 'second' has two coefficients
  # (rule 1). Term order decides, so 'first' is reported.
  terms <- data.frame(
    term_id       = c("first", "first", "second", "second"),
    locus_name    = c("Locus_5", "Locus_5", "Locus_6", "Locus_7"),
    contrast_name = "additive", center_value = 0.5,
    genome_value  = c(1, 1, 2, 3))
  msg <- gwo_error(define_genome_effect_terms(pop, "ADG", terms = terms))
  expect_identical(msg, paste0(
    "term_id 'first': locus 'Locus_5' appears more than once in one term. ",
    "Each locus contributes at most one member."))

  # Rule 2 (effect_name) on the first term beats rule 1 on the second.
  terms <- data.frame(
    term_id       = c("first", "first", "second", "second"),
    locus_name    = c("Locus_5", "Locus_6", "Locus_7", "Locus_8"),
    contrast_name = "additive", center_value = 0.5,
    effect_name   = c("additive", "dominance", "additive", "additive"),
    genome_value  = c(1, 1, 2, 3))
  msg <- gwo_error(define_genome_effect_terms(pop, "ADG", terms = terms))
  expect_identical(msg, paste0(
    "term_id 'first': 'effect_name' must be constant within a term (got ",
    "'additive', 'dominance')."))

  # Within one term, rules run in order: two coefficients beat a repeated locus.
  terms <- data.frame(
    term_id       = "only",
    locus_name    = c("Locus_9", "Locus_9"),
    contrast_name = "additive", center_value = 0.5,
    genome_value  = c(1, 2))
  msg <- gwo_error(define_genome_effect_terms(pop, "ADG", terms = terms))
  expect_match(msg, "^term_id 'only': 'genome_value' must be one")
})

test_that("member-field violations are listed whole, by term then locus_id", {
  pop <- gwo_pop()
  on.exit(close_pop(pop), add = TRUE)
  # Members are canonicalised by locus_id inside a term, so Locus_3 is listed
  # before Locus_12 in term 'b' although it was typed second. Term 'b' comes
  # before term 'a' because it was typed first.
  terms <- data.frame(
    term_id       = c("b", "b", "a"),
    locus_name    = c("Locus_12", "Locus_3", "Locus_2"),
    contrast_name = c("additive", "dominance", "indicator"),
    center_value  = c(NA, 1.5, 0.5),
    dosage_value  = c(NA, NA, NA),
    copy_count_value = c(NA, NA, 2),
    genome_value  = c(1, 1, 2))
  msg <- gwo_error(define_genome_effect_terms(pop, "ADG", terms = terms))
  expect_identical(msg, paste0(
    "Invalid 'terms':\n",
    "  - term_id 'b' locus 'Locus_3': center_value must be between 0 and 1 (got 1.5)\n",
    "  - term_id 'b' locus 'Locus_12' is 'additive' and needs 'center_value' ",
    "(p under Cockerham coding, 0.5 under functional)\n",
    "  - term_id 'a' locus 'Locus_2' is an indicator and needs 'dosage_value'\n",
    "  - term_id 'a' locus 'Locus_2' is an indicator and must not carry ",
    "'center_value' (an indicator is a state, not a centred contrast)"))
})

test_that("origin-field violations are listed whole, in origin-row order", {
  pop <- gwo_pop()
  on.exit(close_pop(pop), add = TRUE)
  terms <- data.frame(
    term_id       = c("t2", "t1"),
    locus_name    = c("Locus_4", "Locus_1"),
    contrast_name = "additive", center_value = 0.5,
    genome_value  = c(1, 2))
  origin <- data.frame(
    term_id         = c("t1", "t2"),
    locus_name      = c("Locus_1", "Locus_4"),
    line_match_type = c("unknown", "exact"),
    line_name       = c("A", NA),
    parent_origin   = c(3L, NA),
    copy_count      = c(0L, 1L))
  msg <- gwo_error(define_genome_effect_terms(pop, "ADG", terms = terms,
                                              origin = origin))
  # Origins are canonicalised by local term index: 't2' was typed first.
  expect_identical(msg, paste0(
    "Invalid 'origin':\n",
    "  - term_id 't2' locus 'Locus_4' is line_match_type 'exact' and needs a line_name\n",
    "  - term_id 't1' locus 'Locus_1': line_match_type 'unknown' takes no line_name (got 'A')\n",
    "  - term_id 't1' locus 'Locus_1': parent_origin must be 1 (sire) or 2 (dam), or NA for either (got 3)\n",
    "  - term_id 't1' locus 'Locus_1' needs copy_count >= 1"))
})


# ── .ge_validate_frames(): the whole ordered violation vector ───────────────

test_that(".ge_validate_frames() returns every violation, in its fixed order", {
  terms <- data.frame(
    # Row order differs from id order on purpose: the per-term phase walks rows.
    id_genome_effect = c(2L, 1L, 3L, 4L, 5L, 6L, 7L, 8L, 10L),
    trait_name = "ADG", effect_owner = "custom", effect_name = NA_character_,
    genome_value = 1, stringsAsFactors = FALSE)
  m <- function(id, slot, locus, kind = "additive", cc = NA_integer_,
                dv = NA_integer_, centre = 0.5) {
    data.frame(id_genome_effect = id, member_slot = slot, locus_id = locus,
               contrast_name = kind, copy_count_value = cc, dosage_value = dv,
               center_value = centre, stringsAsFactors = FALSE)
  }
  members <- rbind(
    # term 1: no members.
    m(2L, 1L, 5L), m(2L, 1L, 5L),                 # same locus, bad slots
    m(3L, 1L, 9L), m(3L, 2L, 3L),                 # not canonicalised
    m(4L, 1L, 7L, "dominance", centre = 0.4),     # genotype member
    m(5L, 1L, 8L),                                # additive member
    m(6L, 1L, 11L), m(7L, 1L, 11L),               # duplicate common terms
    m(8L, 1L, 12L), m(10L, 1L, 12L),              # incomparable scopes
    m(9L, 1L, 13L))                               # orphan member
  o <- function(id, slot, os, lmt, ln = NA_character_, po = NA_integer_,
                cc = 1L) {
    data.frame(id_genome_effect = id, member_slot = slot, origin_slot = os,
               line_match_type = lmt, line_name = ln, parent_origin = po,
               copy_count = cc, stringsAsFactors = FALSE)
  }
  origins <- rbind(
    o(4L, 1L, 1L, "exact", "A"), o(4L, 1L, 3L, "any"),
    o(5L, 1L, 1L, "exact", "A", cc = 2L), o(5L, 1L, 2L, "exact", "B", cc = 2L),
    o(5L, 2L, 1L, "exact", "A"),                  # no such member
    o(8L, 1L, 1L, "exact", "A"),
    o(10L, 1L, 1L, "any", po = 1L))
  labels <- c("1" = "one", "2" = "two", "3" = "three", "4" = "four",
              "5" = "five", "6" = "six", "7" = "seven", "8" = "eight")

  v <- tidybreed:::.ge_validate_frames(terms, members, origins, labels = labels)
  expect_identical(v, c(
    "member rows with no term: 9",
    "origin rows with no member: 5 2",
    "term_id 'two' names the same locus more than once",
    "term_id 'two' member_slot must be 1..n in ascending locus_id order",
    "term_id 'one' has no members",
    "term_id 'three' members are not canonicalized by ascending locus_id",
    "term_id 'four' member 1 origin_slot must be 1..n",
    paste("term_id 'four' member 1 'any' is permitted only on additive members:",
          "a genotype scope must name lines and sum to the",
          "realized copy count, so an 'any' entry constrains",
          "nothing"),
    paste("term_id 'five' member 1 additive members take at most one origin row",
          "(expand alternatives into separate variants)"),
    paste("term_id 'five' member 1 additive origin requires copy_count = 1",
          "(matching is per allele copy)"),
    paste0("term_id 'six' and term_id 'seven' are the same logical term at the ",
           "same scope (duplicate family + scope identity)"),
    paste0("term_id 'eight' and term 10 have overlapping but incomparable ",
           "scopes in one family: neither is more specific, so no variant can ",
           "be selected")))
})

test_that("a whole-table conflict with stored rows names the user's term_id", {
  pop <- gwo_pop()
  on.exit(close_pop(pop), add = TRUE)
  pop <- suppressMessages(define_genome_effect_terms(pop, "ADG", terms = data.frame(
    term_id = "old", locus_name = "Locus_3", contrast_name = "additive",
    center_value = 0.5, genome_value = 1)))
  msg <- gwo_error(suppressMessages(define_genome_effect_terms(
    pop, "ADG", terms = data.frame(
      term_id = c("new_b", "new_a"), locus_name = c("Locus_4", "Locus_3"),
      contrast_name = "additive", center_value = 0.5, genome_value = c(2, 3)))))
  expect_identical(msg, paste0(
    "Invalid genome effects:\n",
    "  - term 1 and term_id 'new_a' are the same logical term at the same scope ",
    "(duplicate family + scope identity)"))
  # Rolled back: only the stored term remains.
  expect_identical(DBI::dbGetQuery(pop$db_conn,
    "SELECT COUNT(*) AS n FROM genome_effects")$n, 1)
})


# ── .ge_pair_key(): distinct pairs, distinct keys ──────────────────────────

test_that(".ge_pair_key() gives distinct pairs distinct keys, at the 2^53 boundary too", {
  key <- tidybreed:::.ge_pair_key
  # Small keys: the exact numeric form x * (max(y) + 1) + y.
  expect_identical(key(c(1L, 1L, 2L), c(1L, 2L, 1L)), c(4, 5, 7))
  expect_false(anyDuplicated(key(c(1L, 1L, 2L), c(1L, 2L, 1L))) > 0L)

  # Radix 2^31 with x = 2^22 puts the keys at 2^53, where (2^22, 3) and
  # (2^22, 4) would round to the same double (Codex review of 5a).
  x <- c(4194304L, 4194304L, 1L)
  y <- c(3L, 4L, 2147483647L)
  k <- key(x, y)
  expect_identical(anyDuplicated(k), 0L)
  expect_identical(duplicated(k), duplicated(paste(x, y)))

  # Equal pairs stay equal on either path; NA is a value, as in duplicated().
  expect_identical(duplicated(key(c(x, 4194304L), c(y, 3L))),
                   c(FALSE, FALSE, FALSE, TRUE))
  expect_identical(duplicated(key(c(1L, 1L, 1L), c(NA, NA, 2L))),
                   c(FALSE, TRUE, FALSE))
})
