# remove_generated_effects(): the one route that deletes generated terms, one
# scope at a time (step-3 review; decided 2026-10-05).

rge_pop <- function(name) {
  set.seed(3401)
  pop <- open_pop(pop_name = name, db_name = ":memory:") |>
    define_genome(n_loci = 12, n_chr = 1, chr_len_Mb = 100) |>
    define_founder_haplotypes(n_haplotypes = 100) |>
    get_table("founder_haplotypes") |>
    add_founders(n_males = 10, n_females = 10, line_name = "A")
  pop <- with_additive_target(pop, "T", 1)
  pop <- with_additive_target(pop, "U", 1)
  gm <- get_table(pop, "genome_meta")
  pop <- suppressMessages(gm |> define_additive_effects("T", seed = 1,
                                                        warn_bounds = NULL))
  pop <- suppressWarnings(suppressMessages(gm |> define_additive_effects("T",
    parent_origin = 1, seed = 2, warn_bounds = NULL)))
  pop <- suppressWarnings(suppressMessages(gm |> define_additive_effects("T",
    line_name = "A", parent_origin = 1, seed = 3, warn_bounds = NULL)))
  suppressMessages(gm |> define_additive_effects("U", seed = 4, warn_bounds = NULL))
}

rge_scopes <- function(pop) {
  DBI::dbGetQuery(pop$db_conn,
    "SELECT e.trait_name, o.line_name, o.parent_origin, COUNT(*) AS n
     FROM genome_effects e
     LEFT JOIN genome_effect_member_origins o USING (id_genome_effect)
     GROUP BY ALL ORDER BY ALL")
}

test_that("remove_generated_effects() deletes exactly one scope", {
  pop <- rge_pop("rge_scope")
  on.exit(close_pop(pop), add = TRUE)
  before <- rge_scopes(pop)
  expect_message(pop <- remove_generated_effects(pop, "T", parent_origin = 1),
                 "Removed 12 generated term\\(s\\).*parent_origin 1 only")
  after <- rge_scopes(pop)
  gone <- before$trait_name == "T" & !is.na(before$parent_origin) &
    is.na(before$line_name)
  expect_identical(after, before[!gone, , drop = FALSE], ignore_attr = TRUE)
  expect_equal(nrow(dplyr::collect(get_table(pop, "trait_var_comp"))), 2L)

  # Line scope, then the common scope.
  expect_error(remove_generated_effects(pop, "T", line_name = "A"),
               "no generated effects at scope line A, both parents")
  pop <- suppressMessages(remove_generated_effects(pop, "T", line_name = "A",
                                                   parent_origin = 1))
  pop <- suppressMessages(remove_generated_effects(pop, "T"))
  expect_equal(DBI::dbGetQuery(pop$db_conn,
    "SELECT COUNT(*) AS n FROM genome_effects WHERE trait_name = 'T'")$n, 0)
  expect_equal(DBI::dbGetQuery(pop$db_conn,
    "SELECT COUNT(*) AS n FROM genome_effects WHERE trait_name = 'U'")$n, 12)
})

test_that("remove_generated_effects() refuses an empty scope or bad input, deleting nothing", {
  pop <- rge_pop("rge_refuse")
  on.exit(close_pop(pop), add = TRUE)
  before <- rge_scopes(pop)
  expect_error(remove_generated_effects(pop, "T", parent_origin = 2),
               "no generated effects at scope all lines, parent_origin 2")
  # All-or-nothing across traits.
  expect_error(remove_generated_effects(pop, c("U", "T"), line_name = "A"),
               "Trait 'U' has no generated")
  expect_error(remove_generated_effects(pop, "T", parent_origin = 3), "1 \\(sire\\)")
  expect_error(remove_generated_effects(pop, "T", line_name = c("A", "B")),
               "one line name")
  expect_identical(rge_scopes(pop), before)

  # User-owner terms are never touched.
  pop <- define_genome_effect_terms(pop, "U", data.frame(
    locus_name = "Locus_1", contrast_name = "additive", center_value = 0.5,
    genome_value = 0.3))
  pop <- suppressMessages(remove_generated_effects(pop, "U"))
  expect_identical(DBI::dbGetQuery(pop$db_conn,
    "SELECT effect_owner FROM genome_effects WHERE trait_name = 'U'")$effect_owner,
    "custom")
})

test_that("remove_generated_effects() removes every kind at the scope (decided 2026-10-05)", {
  # A jointly calibrated model is removed as a whole: generated dominance and
  # A x A terms at the common scope go with the common additive terms, while
  # a line variant at another scope survives. The non-additive terms are
  # planted under the reserved owner (test-only), as step 5's generator would
  # write them.
  pop <- rge_pop("rge_kinds")
  on.exit(close_pop(pop), add = TRUE)
  pop <- suppressMessages(remove_generated_effects(pop, "T", parent_origin = 1))
  p <- extract_allele_freq(get_table(pop, "founder_haplotypes"))
  dom <- c("Locus_2", "Locus_3")
  pop <- .ge_write_terms(pop, "T",
    ad_terms(dom, a = 0, d = c(0.4, -0.3),
             p = p$allele_freq[match(dom, p$locus_name)],
             coding = "cockerham", report = FALSE),
    effect_owner = GE_GENERATED_OWNER, allow_reserved_owner = TRUE)
  pop <- .ge_write_terms(pop, "T", data.frame(
    term_id = c(1L, 1L), locus_name = c("Locus_4", "Locus_5"),
    contrast_name = "additive",
    center_value = p$allele_freq[match(c("Locus_4", "Locus_5"), p$locus_name)],
    genome_value = 0.2),
    effect_owner = GE_GENERATED_OWNER, allow_reserved_owner = TRUE)
  kinds <- function() {
    m <- .gev_read_model(pop$db_conn, "T")
    table(.gev_target_kind(m), .gev_term_line(m), useNA = "ifany")
  }
  expect_setequal(rownames(kinds()), c("additive", "dominance",
                                        "additive_by_additive"))

  pop <- suppressMessages(remove_generated_effects(pop, "T"))
  m <- .gev_read_model(pop$db_conn, "T")
  # Only the line-A, parent-1 additive variant remains.
  expect_equal(unique(.gev_target_kind(m)), "additive")
  expect_equal(unique(.gev_term_line(m)), "A")
  expect_equal(unique(.gev_term_parent(m)), "1")
  expect_equal(nrow(m$terms), 12L)
})
