# Phase A of plans/update_genome_effects_v4.md (v4.9): the origin truth table,
# every case hand-computed before any table exists. Fixtures and the two
# evaluators live in helper-genome-effects.R; the hand computations are worked
# out line by line in plans/update_genome_effects_phase_A.md.
#
# No database, no package writer, no DDL. What these tests establish is that
# every case the plan promises is (a) representable in the three proposed
# tables, (b) accepted or rejected by the stated rules, and (c) has a value
# derivable from the stored rows alone that equals a number computed by hand.

fx <- gefx_fixtures()

test_that("every fixture is either evaluable or a declared rejection", {
  for (nm in names(fx)) {
    f <- fx[[nm]]
    kinds <- c(!is.null(f$expected), !is.null(f$expect_invalid),
               !is.null(f$expect_eval_error))
    expect_equal(sum(kinds), 1L, info = nm)
  }
})

test_that("valid fixtures pass the row-local and cross-row rules", {
  for (nm in names(fx)) {
    f <- fx[[nm]]
    if (is.null(f$expected) && is.null(f$expect_eval_error)) next
    expect_equal(gefx_validate(f$model), character(0), info = nm)
  }
})

test_that("fixture values match the hand computation", {
  for (nm in names(fx)) {
    f <- fx[[nm]]
    if (is.null(f$expected)) next
    for (id in names(f$expected)) {
      got <- gefx_eval_naive(f$model, id)
      expect_equal(got, unname(f$expected[[id]]),
                   tolerance = 1e-12, info = paste(nm, id))
    }
  }
})

test_that("the label-vector factorization equals the per-tuple definition", {
  # The central claim of the plan's Evaluation strategy section, checked before
  # any SQL is written against it.
  for (nm in names(fx)) {
    f <- fx[[nm]]
    if (is.null(f$expected)) next
    for (id in names(f$expected)) {
      expect_equal(gefx_eval_grouped(f$model, id), gefx_eval_naive(f$model, id),
                   tolerance = 1e-12, info = paste(nm, id))
    }
  }
})

test_that("every fixture is evaluable for every individual, not just the named ones", {
  # A fixture that silently errored on an unlisted individual would hide a
  # hole in the predicate system.
  for (nm in names(fx)) {
    f <- fx[[nm]]
    if (is.null(f$expected)) next
    for (id in gefx_individuals()) {
      expect_no_error(gefx_eval_naive(f$model, id))
    }
  }
})

test_that("rejection fixtures are rejected, with the stated reason", {
  for (nm in names(fx)) {
    f <- fx[[nm]]
    if (is.null(f$expect_invalid)) next
    v <- gefx_validate(f$model)
    expect_true(length(v) > 0L, info = nm)
    expect_true(any(grepl(f$expect_invalid, v)), info = paste(nm, ":", paste(v, collapse = "; ")))
  }
})

test_that("dominance at a non-diploid state stops instead of scoring", {
  f <- fx$F14
  expect_error(gefx_eval_naive(f$model, f$eval_ind), f$expect_eval_error)
})

test_that("the combined equivalent for mixed centring is well defined (Q2)", {
  # The writer must reject F16 *with this term in the message*, not with a bare
  # duplicate error: a1(g - c1) + a2(g - c2) = (a1 + a2)(g - c').
  a <- c(2.0, 1.0); cc <- c(0.5, 0.2)
  expect_equal(sum(a), unname(fx$F16$combined[["genome_value"]]))
  expect_equal(sum(a * cc) / sum(a), unname(fx$F16$combined[["center_value"]]))
})

test_that("families partition terms by signature, not by scope", {
  # Scope variants of one term share a family and compete.
  expect_length(gefx_families(fx$F01$model), 1L)
  expect_length(gefx_families(fx$F12$model), 1L)
  # An additive main effect, a dominance main effect and an interaction at the
  # same locus are three families and sum.
  expect_length(gefx_families(fx$F19$model), 3L)
  # Indicator states differing only in copy count are different basis functions.
  expect_length(gefx_families(fx$F09$model), 3L)
})

test_that("containment order is the one the plan's worked cases state", {
  m <- function(origins = NULL) {
    gefx_model(gefx_term(1, 1.0),
               gefx_member(1, 1, "L1", "additive", center = 0.5),
               origins)
  }
  common <- gefx_variant_pred(m(), 1L)
  exactA <- gefx_variant_pred(m(gefx_origin(1, 1, 1, "exact", "A")), 1L)
  exactA1 <- gefx_variant_pred(m(gefx_origin(1, 1, 1, "exact", "A", parent_origin = 1L)), 1L)
  exactB2 <- gefx_variant_pred(m(gefx_origin(1, 1, 1, "exact", "B", parent_origin = 2L)), 1L)
  anyP1 <- gefx_variant_pred(m(gefx_origin(1, 1, 1, "any", parent_origin = 1L)), 1L)

  expect_true(gefx_pred_leq(exactA, common))    # (ANY,ANY) contains (exact A, ANY)
  expect_false(gefx_pred_leq(common, exactA))
  expect_true(gefx_pred_leq(exactA1, exactA))   # the tie v4.1 got wrong
  expect_false(gefx_pred_leq(exactA, exactA1))
  expect_false(gefx_pred_leq(exactA1, exactB2)) # disjoint, both apply
  expect_false(gefx_pred_leq(exactB2, exactA1))
  expect_false(gefx_pred_leq(exactA, anyP1))    # overlap, neither contains
  expect_false(gefx_pred_leq(anyP1, exactA))
  expect_true(gefx_pred_leq(anyP1, common))
})

test_that("genotype containment handles multisets and reciprocals", {
  m <- function(origins = NULL) {
    gefx_model(gefx_term(1, 1.0),
               gefx_member(1, 1, "L1", "dominance", center = 0.5),
               origins)
  }
  common <- gefx_variant_pred(m(), 1L)
  ab <- gefx_variant_pred(m(rbind(gefx_origin(1, 1, 1, "exact", "A"),
                                  gefx_origin(1, 1, 2, "exact", "B"))), 1L)
  recip1 <- gefx_variant_pred(m(rbind(
    gefx_origin(1, 1, 1, "exact", "A", parent_origin = 1L),
    gefx_origin(1, 1, 2, "exact", "B", parent_origin = 2L))), 1L)
  recip2 <- gefx_variant_pred(m(rbind(
    gefx_origin(1, 1, 1, "exact", "A", parent_origin = 2L),
    gefx_origin(1, 1, 2, "exact", "B", parent_origin = 1L))), 1L)
  aa <- gefx_variant_pred(m(gefx_origin(1, 1, 1, "exact", "A", copy_count = 2L)), 1L)

  expect_true(gefx_pred_leq(ab, common))
  expect_true(gefx_pred_leq(recip1, ab))
  expect_false(gefx_pred_leq(ab, recip1))
  expect_false(gefx_pred_leq(recip1, recip2))
  expect_false(gefx_pred_leq(recip2, recip1))
  expect_false(gefx_pred_leq(aa, ab))
  expect_false(gefx_pred_leq(ab, aa))
})
