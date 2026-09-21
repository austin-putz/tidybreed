# The residual adapter of add_phenotype() (Phase 5 of
# plans/sample_correlated_effects.md): sequential conditional residuals across
# phenotypes, per-record strata, user_residual as fixed coordinates, D2/D6 at
# sampling, ordinal pairing on pheno_number. See §5.3–§5.5 and §7 "Residual
# blocks", "Defect 4", "Fixed coordinates", "Repeated records".
#
# Exact tests replay the resolver's contract: `n * m` standard normals in
# entity (sorted id_ind, pheno_number) order and coordinate order within an
# entity, applied through chol(C) (draws = z %*% U). Distributional tests use
# residual_value directly, so no TBV or model noise enters them.

sym <- function(names, values) {
  n <- length(names)
  matrix(values, n, n, byrow = TRUE, dimnames = list(names, names))
}

# Founders A_1..A_n with traits (QTL placed) and one continuous repeatable
# phenotype per trait; residual variance is declared by the test.
make_resid_pop <- function(pop_name, traits = c("A", "B"), n_males = 6,
                           n_females = 6, seed = 11, ...) {
  set.seed(seed)
  pop <- suppressMessages(make_test_pop(pop_name, n_loci = 60, n_chr = 1,
                                        n_males = n_males, n_females = n_females))
  for (t in traits) {
    pop <- suppressMessages(define_trait(pop, t, target_add_var = 1))
    pop <- suppressMessages(
      pop |> get_table("genome_meta") |> define_additive_effects(t))
    pop <- suppressMessages(define_phenotype(
      pop, t, type = "continuous", mean = 10, repeatable = TRUE, ...))
  }
  pop
}

set_resid <- function(pop, phenos, R, ...) {
  suppressMessages(define_residual_cov(pop, phenos, R, ...))
}

# residual_value / residual_condition_level of one phenotype at one
# pheno_number, named by id_ind, sorted
resid_rows_of <- function(pop, t, pheno_number = 1L) {
  DBI::dbGetQuery(pop$db_conn, paste0(
    "SELECT id_ind, residual_value, residual_condition_level FROM ind_phenotype ",
    "WHERE phenotype_name = '", t, "' AND pheno_number = ", pheno_number,
    " ORDER BY id_ind"))
}
resid_of <- function(pop, t, pheno_number = 1L) {
  r <- resid_rows_of(pop, t, pheno_number)
  stats::setNames(r$residual_value, r$id_ind)
}
level_of <- function(pop, t, pheno_number = 1L) {
  r <- resid_rows_of(pop, t, pheno_number)
  stats::setNames(r$residual_condition_level, r$id_ind)
}

phen_ids <- function(pop, t) sort(DBI::dbGetQuery(pop$db_conn, paste0(
  "SELECT DISTINCT id_ind FROM ind_phenotype WHERE phenotype_name = '", t, "'"))$id_ind)

ids_of <- function(pop, where = "") sort(DBI::dbGetQuery(pop$db_conn,
  paste0("SELECT id_ind FROM ind_meta ", where))$id_ind)

sex_of <- function(pop) {
  m <- DBI::dbGetQuery(pop$db_conn, "SELECT id_ind, sex FROM ind_meta")
  stats::setNames(m$sex, m$id_ind)
}

state_after <- function(seed, expr) {
  set.seed(seed)
  force(expr)
  .Random.seed
}

add_ph <- function(pop, phenos, tbl = get_table(pop, "ind_meta"), ...) {
  suppressMessages(add_phenotype(tbl, phenos, ...))
}

# Standard normals in resolver order for n entities x m coordinates.
z_mat <- function(n, m) matrix(stats::rnorm(n * m), n, m, byrow = TRUE)

R_AB <- sym(c("A", "B"), c(1, .8, .8, 1))
U_AB <- chol(R_AB)   # A = z1; B = .8 z1 + .6 z2


# ── Residual blocks: sequential conditioning ────────────────────────────────

test_that("sequential A -> B on the same individuals: B is drawn conditional on the stored A (exact and in distribution)", {
  pop <- make_resid_pop("rs_seq", n_males = 200, n_females = 200)
  on.exit(close_pop(pop))
  pop <- set_resid(pop, c("A", "B"), R_AB)
  n <- 400L

  pop <- add_ph(pop, "A", seed = 1)
  a <- resid_of(pop, "A")
  set.seed(1)
  expect_equal(unname(a), stats::rnorm(n))              # marginal: sd 1, entity order
  expect_true(all(is.na(level_of(pop, "A"))))

  pop <- add_ph(pop, "B", seed = 2)
  b <- resid_of(pop, "B")
  set.seed(2)
  expect_equal(unname(b), .8 * unname(a) + .6 * stats::rnorm(n))   # E[B|A] + sd .6
  expect_equal(names(b), names(a))

  # Distribution on 400 animals: cor ~ .8 (SE ~ .018), var(B) ~ 1
  expect_lt(abs(stats::cor(a, b) - .8), .07)
  expect_lt(abs(stats::var(b) - 1), .25)
  # Nothing but residual_value / residual_condition_level was affected
  ph <- dplyr::collect(get_table(pop, "ind_phenotype"))
  expect_false(anyNA(ph$residual_value))
  expect_true(all(is.na(ph$residual_condition_level)))
})

test_that("A -> B on a culled subset: survivors condition on their A, the culled get no B", {
  pop <- make_resid_pop("rs_cull", n_males = 150, n_females = 150)
  on.exit(close_pop(pop))
  pop <- set_resid(pop, c("A", "B"), R_AB)

  pop <- add_ph(pop, "A", seed = 3)
  a <- resid_of(pop, "A")
  males <- ids_of(pop, "WHERE sex = 'M'")
  pop <- add_ph(pop, "B", tbl = get_table(pop, "ind_meta") |> dplyr::filter(sex == "M"),
                seed = 4)
  b <- resid_of(pop, "B")
  expect_identical(names(b), males)
  set.seed(4)
  expect_equal(unname(b), .8 * unname(a[males]) + .6 * stats::rnorm(length(males)))
  expect_lt(abs(stats::cor(a[males], b) - .8), .09)
  # Untouched: the females' A rows are exactly as written
  expect_identical(resid_of(pop, "A"), a)
})

test_that("A and B in one call on partially overlapping planned sets: joint on the overlap, marginal elsewhere", {
  pop <- make_resid_pop("rs_partial", n_males = 6, n_females = 6)
  on.exit(close_pop(pop))
  pop <- suppressMessages(define_phenotype(pop, "B", expressed_sex = "M",
                                           repeatable = TRUE, overwrite = TRUE))
  pop <- set_resid(pop, c("A", "B"), R_AB)

  ids   <- ids_of(pop)
  males <- ids_of(pop, "WHERE sex = 'M'")
  pop <- add_ph(pop, c("A", "B"), seed = 5)
  a <- resid_of(pop, "A"); b <- resid_of(pop, "B")
  expect_identical(names(a), ids)
  expect_identical(names(b), males)

  # Groups in sorted key order: (unconditional, "A") before (unconditional, "A,B")
  set.seed(5)
  z_f <- stats::rnorm(length(ids) - length(males))          # A-only entities
  z_m <- z_mat(length(males), 2L) %*% U_AB                  # joint entities
  expect_equal(unname(a[setdiff(ids, males)]), z_f)
  expect_equal(unname(a[males]), z_m[, 1L])
  expect_equal(unname(b),        z_m[, 2L])
  # exactly n_A + n_B normals consumed
  expect_identical(state_after(5, add_ph(pop, c("A", "B"))),
                   state_after(5, stats::rnorm(length(ids) + length(males))))
})

test_that("B for a mixed group: individuals with a stored A condition on it, the others get the marginal draw", {
  pop <- make_resid_pop("rs_mixed", n_males = 150, n_females = 150)
  on.exit(close_pop(pop))
  pop <- set_resid(pop, c("A", "B"), R_AB)
  males   <- ids_of(pop, "WHERE sex = 'M'")
  females <- ids_of(pop, "WHERE sex = 'F'")

  pop <- add_ph(pop, "A", tbl = get_table(pop, "ind_meta") |> dplyr::filter(sex == "M"),
                seed = 6)
  a <- resid_of(pop, "A")
  pop <- add_ph(pop, "B", seed = 7)
  b <- resid_of(pop, "B")
  expect_length(b, 300L)

  # One resolver call (one stratum, sample set {B}); the observed pattern
  # grouping happens inside it, so the stream is entity order regardless.
  set.seed(7)
  z <- stats::rnorm(300L)
  names(z) <- names(b)
  expect_equal(unname(b[males]),   .8 * unname(a[males]) + .6 * unname(z[males]))
  expect_equal(unname(b[females]), unname(z[females]))
  expect_lt(abs(stats::cor(a[males], b[males]) - .8), .1)
  expect_lt(abs(stats::var(b[females]) - 1), .3)
})

test_that("B first, A later is symmetric", {
  pop <- make_resid_pop("rs_sym", n_males = 150, n_females = 150)
  on.exit(close_pop(pop))
  pop <- set_resid(pop, c("A", "B"), R_AB)
  pop <- add_ph(pop, "B", seed = 8)
  pop <- add_ph(pop, "A", seed = 9)
  a <- resid_of(pop, "A"); b <- resid_of(pop, "B")
  set.seed(9)
  expect_equal(unname(a), .8 * unname(b) + .6 * stats::rnorm(300L))
  expect_lt(abs(stats::cor(a, b) - .8), .08)
})

test_that("a three-phenotype block with several observation patterns: C conditions on whatever each entity has", {
  pop <- make_resid_pop("rs_abc", traits = c("A", "B", "C"),
                        n_males = 100, n_females = 100)
  on.exit(close_pop(pop))
  R <- sym(c("A", "B", "C"), c(1, .5, .3,
                               .5, 1, .6,
                               .3, .6, 1))
  pop <- set_resid(pop, c("A", "B", "C"), R)
  males   <- ids_of(pop, "WHERE sex = 'M'")
  females <- ids_of(pop, "WHERE sex = 'F'")

  pop <- add_ph(pop, "A", seed = 10)
  pop <- add_ph(pop, "B", tbl = get_table(pop, "ind_meta") |> dplyr::filter(sex == "M"),
                seed = 11)
  pop <- add_ph(pop, "C", seed = 12)
  a <- resid_of(pop, "A"); b <- resid_of(pop, "B"); c_ <- resid_of(pop, "C")

  set.seed(12)
  z <- stats::setNames(stats::rnorm(200L), names(c_))
  # females: C | A
  w_f  <- R["C", "A"] / R["A", "A"]
  sd_f <- sqrt(R["C", "C"] - w_f * R["A", "C"])
  expect_equal(unname(c_[females]), w_f * unname(a[females]) + sd_f * unname(z[females]))
  # males: C | A, B
  W_m  <- R["C", c("A", "B")] %*% solve(R[c("A", "B"), c("A", "B")])
  sd_m <- sqrt(drop(R["C", "C"] - W_m %*% R[c("A", "B"), "C"]))
  mu_m <- drop(cbind(a[males], b[males]) %*% t(W_m))
  expect_equal(unname(c_[males]), unname(mu_m) + sd_m * unname(z[males]))
})

test_that("disconnected components draw independently, in block order", {
  pop <- make_resid_pop("rs_disc", traits = c("A", "B", "C"))
  on.exit(close_pop(pop))
  pop <- set_resid(pop, c("A", "B"), R_AB)
  pop <- set_resid(pop, "C", matrix(4, 1, 1, dimnames = list("C", "C")))
  n <- 12L

  pop <- add_ph(pop, c("C", "A", "B"), seed = 13)   # call order != block order
  a <- resid_of(pop, "A"); b <- resid_of(pop, "B"); c_ <- resid_of(pop, "C")
  set.seed(13)
  z_ab <- z_mat(n, 2L) %*% U_AB                     # block {A, B} first
  z_c  <- stats::rnorm(n)                           # then block {C}
  expect_equal(unname(a),  z_ab[, 1L])
  expect_equal(unname(b),  z_ab[, 2L])
  expect_equal(unname(c_), 2 * z_c)
  expect_identical(state_after(13, add_ph(pop, c("A", "B", "C"))),
                   state_after(13, stats::rnorm(3L * n)))
})

test_that("an explicit zero covariance keeps two phenotypes in one block and they come out uncorrelated", {
  pop <- make_resid_pop("rs_zero", n_males = 150, n_females = 150)
  on.exit(close_pop(pop))
  pop <- set_resid(pop, c("A", "B"), sym(c("A", "B"), c(1, 0, 0, 1)))
  expect_length(find_covariance_blocks(pop$db_conn, "residual", "A")[[1L]]$phenotypes, 2L)

  pop <- add_ph(pop, "A", seed = 14)
  pop <- add_ph(pop, "B", seed = 15)
  a <- resid_of(pop, "A"); b <- resid_of(pop, "B")
  set.seed(15)
  expect_equal(unname(b), stats::rnorm(300L))        # conditioning on A changes nothing
  expect_lt(abs(stats::cor(a, b)), .12)
})

test_that("a phenotype in no residual block errors as before", {
  pop <- make_resid_pop("rs_none", traits = "A")
  on.exit(close_pop(pop))
  expect_error(add_ph(pop, "A"), "No residual variance found for phenotype 'A'")
  expect_equal(nrow(dplyr::collect(get_table(pop, "ind_phenotype"))), 0L)
})


# ── Defect 4: heterogeneous residuals on single-phenotype calls ─────────────

test_that("a single-phenotype call draws each record from its own stratum and stores the stratum", {
  pop <- make_resid_pop("rs_het", traits = "A", n_males = 6, n_females = 6)
  on.exit(close_pop(pop))
  pop <- set_resid(pop, "A", matrix(4, 1, 1, dimnames = list("A", "A")),
                   condition_column = "sex", condition_level = "M")
  pop <- set_resid(pop, "A", matrix(9, 1, 1, dimnames = list("A", "A")),
                   condition_column = "sex", condition_level = "F")
  # No unconditional stratum at all: the call still succeeds
  sex <- sex_of(pop)

  pop <- add_ph(pop, "A", seed = 16)
  a <- resid_of(pop, "A")
  expect_identical(unname(level_of(pop, "A")), unname(sex[names(a)]))
  # groups in sorted stratum order: F entities first, then M
  set.seed(16)
  f <- names(a)[sex[names(a)] == "F"]; m <- names(a)[sex[names(a)] == "M"]
  expect_equal(unname(a[f]), 3 * stats::rnorm(length(f)))
  expect_equal(unname(a[m]), 2 * stats::rnorm(length(m)))
})

test_that("a level matching no stratum falls back to the unconditional R (stored NULL) with a warning; NULL falls back silently; no unconditional R is an error", {
  pop <- make_resid_pop("rs_fallback", traits = "A", n_males = 6, n_females = 6)
  on.exit(close_pop(pop))
  conn <- pop$db_conn
  ids <- ids_of(pop)
  pop <- suppressWarnings(suppressMessages(
    pop |> get_table("ind_meta") |> mutate_table(farm = "F1")))
  pop <- suppressWarnings(suppressMessages(
    pop |> get_table("ind_meta") |> dplyr::filter(id_ind %in% c("A_2", "A_3")) |>
      mutate_table(farm = "F9")))
  DBI::dbExecute(conn, "UPDATE ind_meta SET farm = NULL WHERE id_ind = 'A_4'")

  pop <- set_resid(pop, "A", matrix(4, 1, 1, dimnames = list("A", "A")),
                   condition_column = "farm", condition_level = "F1")
  expect_error(add_ph(pop, "A"),
               "3 planned record\\(s\\).*farm = 'F9', NULL.*no unconditional stratum.*A_2, A_3, A_4")
  expect_equal(nrow(dplyr::collect(get_table(pop, "ind_phenotype"))), 0L)

  pop <- suppressMessages(define_phenotype(pop, "A", residual_var = 1, overwrite = TRUE))
  expect_warning(add_ph(pop, "A", seed = 17),
                 "2 planned record\\(s\\).*farm value matching no residual stratum \\('F9'\\).*A_2, A_3")
  a  <- resid_of(pop, "A")
  lv <- level_of(pop, "A")
  expect_identical(unname(lv[names(a) %in% c("A_2", "A_3", "A_4")]), rep(NA_character_, 3L))
  expect_true(all(lv[!names(a) %in% c("A_2", "A_3", "A_4")] == "F1"))
  set.seed(17)
  fb <- c("A_2", "A_3", "A_4")                      # unconditional group first
  expect_equal(unname(a[fb]), stats::rnorm(3L))
  expect_equal(unname(a[setdiff(names(a), fb)]), 2 * stats::rnorm(9L))
})


# ── Fixed coordinates: user_residual ────────────────────────────────────────

test_that("a supplied A conditions a generated B in the same call, and both are stored", {
  pop <- make_resid_pop("rs_fixed", n_males = 100, n_females = 100)
  on.exit(close_pop(pop))
  pop <- set_resid(pop, c("A", "B"), R_AB)
  n <- 200L
  set.seed(18); a_user <- stats::rnorm(n)

  pop <- add_ph(pop, c("A", "B"), user_residual = list(A = a_user), seed = 19)
  a <- resid_of(pop, "A"); b <- resid_of(pop, "B")
  expect_equal(unname(a), a_user)
  set.seed(19)
  expect_equal(unname(b), .8 * a_user + .6 * stats::rnorm(n))
  expect_lt(abs(stats::cor(a, b) - .8), .1)
  expect_identical(state_after(19, add_ph(pop, c("A", "B"), user_residual = list(A = a_user))),
                   state_after(19, stats::rnorm(n)))   # only B consumed RNG
})

test_that("both supplied: RNG-neutral, both stored, and a later C conditions on them", {
  pop <- make_resid_pop("rs_fixed2", traits = c("A", "B", "C"))
  on.exit(close_pop(pop))
  R <- sym(c("A", "B", "C"), c(1, .5, .3,
                               .5, 1, .6,
                               .3, .6, 1))
  pop <- set_resid(pop, c("A", "B", "C"), R)
  n <- 12L
  ur <- list(A = seq_len(n) / 10, B = -seq_len(n) / 20)
  expect_identical(state_after(20, add_ph(pop, c("A", "B"), user_residual = ur)),
                   state_after(20, NULL))
  expect_equal(unname(resid_of(pop, "A")), ur$A)
  expect_equal(unname(resid_of(pop, "B")), ur$B)

  pop <- add_ph(pop, "C", seed = 21)
  W  <- R["C", c("A", "B")] %*% solve(R[c("A", "B"), c("A", "B")])
  sd <- sqrt(drop(R["C", "C"] - W %*% R[c("A", "B"), "C"]))
  set.seed(21)
  expect_equal(unname(resid_of(pop, "C")),
               drop(cbind(ur$A, ur$B) %*% t(W)) + sd * stats::rnorm(n))
})

test_that("user_residual is validated against the plan", {
  pop <- make_resid_pop("rs_ur_val", traits = c("A", "B", "D"))
  on.exit(close_pop(pop))
  pop <- set_resid(pop, c("A", "B"), R_AB)
  pop <- suppressMessages(define_phenotype(pop, "D", type = "derived_formula",
                                           formula = "2 * A", overwrite = TRUE))

  expect_error(add_ph(pop, c("A", "B"), user_residual = seq_len(12)),
               "named list.*more than one phenotype.*\\{A, B\\}")
  expect_error(add_ph(pop, c("A", "B"), user_residual = list(A = 1:12, Z = 1:12)),
               "not generated from the model.*\\{Z\\}")
  expect_error(add_ph(pop, c("A", "D"), user_residual = list(D = 1:12)),
               "not generated from the model.*\\{D\\}")
  expect_error(add_ph(pop, "A", user_residual = 1:3),
               "must equal 12 \\(the planned records")
  expect_error(add_ph(pop, "A", user_residual = c(1:11, NA)), "must be finite")
  expect_error(add_ph(pop, "A", user_residual = list(1:12)), "named list")
  expect_error(add_ph(pop, "A", user_residual = stats::setNames(1:12, ids_of(pop))),
               "per-id_ind names are not supported")
  expect_error(add_ph(pop, "A", user_residual = 1:12, user_values = 1:12),
               "cannot be combined")
  expect_equal(nrow(dplyr::collect(get_table(pop, "ind_phenotype"))), 0L)

  # A plain vector is fine when only one phenotype is on the model path
  pop <- add_ph(pop, c("A", "D"), user_residual = seq_len(12) / 10)
  expect_equal(unname(resid_of(pop, "A")), seq_len(12) / 10)
  expect_true(all(is.na(DBI::dbGetQuery(pop$db_conn,
    "SELECT residual_value FROM ind_phenotype WHERE phenotype_name = 'D'")$residual_value)))
})

test_that("a supplied value off the support of a singular R errors in the resolver, before any draw", {
  pop <- make_resid_pop("rs_support")
  on.exit(close_pop(pop))
  pop <- set_resid(pop, c("A", "B"), sym(c("A", "B"), c(1, 0, 0, 0)))   # B has zero variance
  expect_error(add_ph(pop, c("A", "B"), user_residual = list(B = rep(.5, 12))),
               "outside the support")
  expect_equal(nrow(dplyr::collect(get_table(pop, "ind_phenotype"))), 0L)
  pop <- add_ph(pop, c("A", "B"), user_residual = list(B = rep(0, 12)))
  expect_equal(unname(resid_of(pop, "B")), rep(0, 12))
})


# ── Repeated records: ordinal pairing on pheno_number ───────────────────────

test_that("distinct phenotypes pair on equal pheno_number; repeats of one phenotype are independent", {
  pop <- make_resid_pop("rs_repeat")
  on.exit(close_pop(pop))
  pop <- set_resid(pop, c("A", "B"), R_AB)
  n <- 12L

  pop <- add_ph(pop, "A", seed = 22)
  pop <- add_ph(pop, "A", seed = 23)
  a1 <- resid_of(pop, "A", 1L); a2 <- resid_of(pop, "A", 2L)
  set.seed(23)
  expect_equal(unname(a2), stats::rnorm(n))          # no other member at pheno_number 2

  pop <- add_ph(pop, "B", seed = 24)                 # B's first record pairs with A(1)
  b1 <- resid_of(pop, "B", 1L)
  set.seed(24)
  expect_equal(unname(b1), .8 * unname(a1) + .6 * stats::rnorm(n))
})

test_that("unequal record counts never condition on the wrong record", {
  pop <- make_resid_pop("rs_unequal")
  on.exit(close_pop(pop))
  pop <- set_resid(pop, c("A", "B"), R_AB)
  males <- ids_of(pop, "WHERE sex = 'M'")

  pop <- add_ph(pop, "A", seed = 25)
  pop <- add_ph(pop, "A", tbl = get_table(pop, "ind_meta") |> dplyr::filter(sex == "M"),
                seed = 26)
  a1 <- resid_of(pop, "A", 1L)
  pop <- add_ph(pop, "B", seed = 27)
  b1 <- resid_of(pop, "B", 1L)
  set.seed(27)
  expect_equal(unname(b1), .8 * unname(a1) + .6 * stats::rnorm(12L))

  # A's third record for the males and B's second: entities (id, 2) for B see A(2)
  pop <- add_ph(pop, "B", tbl = get_table(pop, "ind_meta") |> dplyr::filter(sex == "M"),
                seed = 28)
  a2 <- resid_of(pop, "A", 2L); b2 <- resid_of(pop, "B", 2L)
  expect_identical(names(b2), males)
  set.seed(28)
  expect_equal(unname(b2), .8 * unname(a2[males]) + .6 * stats::rnorm(length(males)))
})


# ── D2 / D6 at sampling ─────────────────────────────────────────────────────

make_stratum_pop <- function(name, ...) {
  pop <- make_resid_pop(name, ...)
  pop <- suppressWarnings(suppressMessages(
    pop |> get_table("ind_meta") |> mutate_table(farm = "F1")))
  R1 <- R_AB; R2 <- 4 * R_AB
  pop <- set_resid(pop, c("A", "B"), R1, condition_column = "farm", condition_level = "F1")
  pop <- set_resid(pop, c("A", "B"), R2, condition_column = "farm", condition_level = "F2")
  pop
}

move_farm <- function(pop, ids, farm) {
  suppressWarnings(suppressMessages(
    pop |> get_table("ind_meta") |> dplyr::filter(id_ind %in% !!ids) |>
      mutate_table(farm = farm)))
}

test_that("D2: a stored residual under a different stratum is an error by default, naming the records", {
  pop <- make_stratum_pop("rs_d2_error")
  on.exit(close_pop(pop))
  pop <- add_ph(pop, "A", seed = 29)
  expect_true(all(level_of(pop, "A") == "F1"))

  pop <- move_farm(pop, c("A_3", "A_5"), "F2")
  expect_error(add_ph(pop, "B"),
               "2 record\\(s\\) in block \\{A, B\\}.*A_3 \\(A: stored under 'F1', now 'F2'\\).*condition_change_action = 'independent'")
  expect_equal(nrow(dplyr::collect(get_table(pop, "ind_phenotype"))), 12L)
})

test_that("D2: 'independent' drops only the mismatched stored coordinates, warns once, and draws from the current stratum", {
  # (Set at definition: the D6 check refuses to flip one member of a defined
  # block, so the value cannot be changed afterwards.)
  pop <- make_stratum_pop("rs_d2_indep", condition_change_action = "independent")
  on.exit(close_pop(pop))
  pop <- add_ph(pop, "A", seed = 30)
  a <- resid_of(pop, "A")
  moved <- c("A_3", "A_5")
  pop <- move_farm(pop, moved, "F2")

  expect_warning(pop <- add_ph(pop, "B", seed = 31),
                 "2 record\\(s\\) in block \\{A, B\\}.*dropped \\{A\\}.*A_3 \\(A: stored under 'F1', now 'F2'\\)")
  b  <- resid_of(pop, "B")
  lv <- level_of(pop, "B")
  expect_identical(unname(lv[names(b) %in% moved]), c("F2", "F2"))
  expect_true(all(lv[!names(b) %in% moved] == "F1"))
  # Groups in sorted order: F1 (conditioned on A) then F2 (marginal under 4R)
  set.seed(31)
  stay <- setdiff(names(b), moved)
  expect_equal(unname(b[stay]),  .8 * unname(a[stay]) + .6 * stats::rnorm(length(stay)))
  expect_equal(unname(b[moved]), 2 * stats::rnorm(2L))
})

test_that("D6 at sampling: disagreeing condition_change_action across the block is an error", {
  pop <- make_stratum_pop("rs_d6")
  on.exit(close_pop(pop))
  pop <- add_ph(pop, "A", seed = 32)
  pop <- move_farm(pop, "A_3", "F2")
  DBI::dbExecute(pop$db_conn,
    "UPDATE phenotype_meta SET condition_change_action = 'independent' WHERE phenotype_name = 'B'")
  expect_error(add_ph(pop, "B"),
               "add_phenotype\\(\\): `condition_change_action` must agree across residual covariance block \\{A, B\\}: A = 'error', B = 'independent'")
})

test_that("sex as a condition column never triggers the change path", {
  pop <- make_resid_pop("rs_sex", n_males = 6, n_females = 6)
  on.exit(close_pop(pop))
  pop <- set_resid(pop, c("A", "B"), R_AB, condition_column = "sex", condition_level = "M")
  pop <- set_resid(pop, c("A", "B"), 4 * R_AB, condition_column = "sex", condition_level = "F")
  sex <- sex_of(pop)
  pop <- add_ph(pop, "A", seed = 33)
  pop <- add_ph(pop, "B", seed = 34)
  a <- resid_of(pop, "A"); b <- resid_of(pop, "B")
  expect_identical(unname(level_of(pop, "B")), unname(sex[names(b)]))
  set.seed(34)
  f <- names(b)[sex[names(b)] == "F"]; m <- names(b)[sex[names(b)] == "M"]
  expect_equal(unname(b[f]), .8 * unname(a[f]) + 2 * .6 * stats::rnorm(length(f)))
  expect_equal(unname(b[m]), .8 * unname(a[m]) + .6 * stats::rnorm(length(m)))
})


# ── The realization lock is now real ────────────────────────────────────────

test_that("a residual block is locked by the records add_phenotype() writes", {
  pop <- make_resid_pop("rs_lock")
  on.exit(close_pop(pop))
  pop <- set_resid(pop, c("A", "B"), R_AB)
  pop <- add_ph(pop, "A")
  err <- tryCatch(define_residual_cov(pop, c("A", "B"), 2 * R_AB),
                  error = function(e) conditionMessage(e))
  expect_match(err, "has 12 realized draws in ind_phenotype")
  pop <- suppressMessages(
    pop |> get_table("ind_phenotype") |>
      dplyr::filter(phenotype_name %in% c("A", "B"), !is.na(residual_value)) |>
      remove_rows())
  pop <- set_resid(pop, c("A", "B"), 2 * R_AB)
  expect_equal(find_covariance_blocks(pop$db_conn, "residual", "A")[[1L]]$unconditional,
               2 * R_AB)
})
