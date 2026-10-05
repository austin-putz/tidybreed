# Step 3c (plans/import_qtl_effect_methods.md Q18, gate PH2): the
# formula_tgv DSL. Every reference reads the contributor's total genetic value
# by default, or the one component a named `component =` gives; the group
# calls read the group column from a named `table =`. define_phenotype()
# refuses anything else before writing.

# 20 founders (10 M, 10 F) and 30 offspring, 12 loci. Trait T has generated
# additive effects at loci 1-8 and user-owner dominance terms at loci 2-4, so
# its total is not its additive value. Offspring are in pens of 5 in
# ind_meta, and in pens of 10 in a separate table `pens`.
dsl_pop <- function(name) {
  set.seed(3501)
  pop <- open_pop(pop_name = name, db_name = ":memory:") |>
    define_genome(n_loci = 12, n_chr = 1, chr_len_Mb = 50) |>
    define_founder_haplotypes(n_haplotypes = 60, method = "beta")
  pop <- pop |> get_table("founder_haplotypes") |>
    add_founders(n_males = 10, n_females = 10, line_name = "A")
  pop <- with_additive_target(pop, "T", 1)
  pop <- suppressMessages(pop |> get_table("genome_meta") |>
    dplyr::filter(locus_id <= 8L) |>
    define_additive_effects("T", warn_bounds = NULL))
  loci <- pop |> get_table("genome_meta") |> dplyr::collect() |>
    dplyr::arrange(.data$locus_id) |> dplyr::pull("locus_name")
  p <- extract_allele_freq(get_table(pop, "founder_haplotypes"))
  dom <- loci[2:4]
  pop <- define_genome_effect_terms(pop, "T",
    ad_terms(dom, a = 0, d = c(1.6, -1.4, 1.3),
             p = p$allele_freq[match(dom, p$locus_name)],
             coding = "cockerham", report = FALSE),
    effect_owner = "dom")

  ind <- dplyr::collect(get_table(pop, "ind_meta"))
  matings <- data.frame(
    id_parent_1 = rep(ind$id_ind[ind$sex == "M"], length.out = 30),
    id_parent_2 = rep(ind$id_ind[ind$sex == "F"], length.out = 30),
    sex = rep(c("M", "F"), 15), line_name = "A", stringsAsFactors = FALSE)
  pop <- suppressMessages(add_offspring(pop, matings = matings))

  ids <- dplyr::collect(get_table(pop, "ind_meta"))$id_ind
  pen <- paste0("p", (seq_along(ids) - 1L) %/% 5L)
  pop <- get_table(pop, "ind_meta") |>
    mutate_table(pen = tibble::tibble(id_ind = ids, pen = pen))
  DBI::dbWriteTable(pop$db_conn, "pens", data.frame(
    id_ind = ids, pen = paste0("q", (seq_along(ids) - 1L) %/% 10L),
    stringsAsFactors = FALSE))
  pop
}

dsl_offspring <- function(pop) {
  pop |> get_table("ind_meta") |> dplyr::filter(!is.na(id_parent_1))
}

dsl_records <- function(pop, phenotype) {
  DBI::dbGetQuery(pop$db_conn, paste0(
    "SELECT p.id_ind, p.pheno_value, p.residual_value, m.id_parent_2 ",
    "FROM ind_phenotype p JOIN ind_meta m USING (id_ind) ",
    "WHERE p.phenotype_name = '", phenotype, "' ORDER BY p.id_ind"))
}

# A component's value per id (sum over the listed components; "total" = all).
dsl_value <- function(pop, ids, component = "total") {
  tg <- dplyr::collect(get_table(pop, "ind_tgv"))
  tg <- tg[tg$trait_name == "T", ]
  if (!identical(component, "total")) tg <- tg[tg$component_name %in% component, ]
  v <- tapply(tg$tgv_value, tg$id_ind, sum)
  out <- as.vector(v[ids])
  out[is.na(out)] <- 0
  out
}

# Hand-computed group-mate sum: the other members of the focal's group.
dsl_mate_sum <- function(pop, ids, table, component = "total") {
  g <- DBI::dbGetQuery(pop$db_conn, paste0("SELECT id_ind, pen FROM ", table))
  grp <- stats::setNames(g$pen, g$id_ind)
  vapply(ids, function(i) {
    mates <- setdiff(g$id_ind[g$pen == grp[[i]]], i)
    sum(dsl_value(pop, mates, component))
  }, numeric(1), USE.NAMES = FALSE)
}

phenotype_rows <- function(pop) {
  DBI::dbGetQuery(pop$db_conn, "SELECT COUNT(*) AS n FROM phenotype_meta")$n
}


test_that("PH2: references read the total by default and a named component when given", {
  pop <- dsl_pop("ph2_comp")
  on.exit(close_pop(pop), add = TRUE)
  pop <- define_phenotype(pop, "W_tot", residual_var = 1,
                          formula_tgv = "T + dam(T)")
  pop <- define_phenotype(pop, "W_add", residual_var = 1,
    formula_tgv = "self(T, component = \"total\") + dam(T, component = \"additive\")")
  pop <- define_phenotype(pop, "W_dom", residual_var = 1,
    formula_tgv = "sire(T, component = 'dominance')")
  stored <- DBI::dbGetQuery(pop$db_conn, paste0(
    "SELECT formula_tgv FROM phenotype_meta WHERE phenotype_name = 'W_add'"))
  expect_match(stored$formula_tgv, "component = \"additive\"", fixed = TRUE)

  set.seed(7)
  pop <- suppressMessages(dsl_offspring(pop) |>
    add_phenotype(c("W_tot", "W_add", "W_dom")))

  r <- dsl_records(pop, "W_tot")
  expect_equal(nrow(r), 30L)
  expect_equal(r$pheno_value - r$residual_value,
               dsl_value(pop, r$id_ind) + dsl_value(pop, r$id_parent_2),
               tolerance = 1e-12)

  r <- dsl_records(pop, "W_add")
  g <- dsl_value(pop, r$id_ind) + dsl_value(pop, r$id_parent_2, "additive")
  expect_equal(r$pheno_value - r$residual_value, g, tolerance = 1e-12)
  # The dam's additive value is not her total: dominance is present.
  expect_false(isTRUE(all.equal(dsl_value(pop, r$id_parent_2, "additive"),
                                dsl_value(pop, r$id_parent_2))))

  r <- dsl_records(pop, "W_dom")
  sires <- DBI::dbGetQuery(pop$db_conn, paste0(
    "SELECT id_ind, id_parent_1 FROM ind_meta WHERE id_parent_1 IS NOT NULL ",
    "ORDER BY id_ind"))
  expect_equal(r$pheno_value - r$residual_value,
               dsl_value(pop, sires$id_parent_1[match(r$id_ind, sires$id_ind)],
                         "dominance"), tolerance = 1e-12)
})

test_that("PH2: table = reads the named table; group terms differing only in table stay distinct", {
  pop <- dsl_pop("ph2_table")
  on.exit(close_pop(pop), add = TRUE)
  expect_warning(pop <- define_phenotype(pop, "G2", residual_var = 1,
    formula_tgv = "group_sum(T, pen) + 2 * group_sum(T, pen, table = \"pens\")"),
    "Scalar arithmetic constant")
  pop <- define_phenotype(pop, "Gc", residual_var = 1,
    formula_tgv = "group_mean(T, pen, table = pens, component = \"additive\")")
  set.seed(8)
  pop <- suppressMessages(suppressWarnings(dsl_offspring(pop) |>
    add_phenotype(c("G2", "Gc"))))

  r <- dsl_records(pop, "G2")
  expect_equal(r$pheno_value - r$residual_value,
               dsl_mate_sum(pop, r$id_ind, "ind_meta") +
                 2 * dsl_mate_sum(pop, r$id_ind, "pens"),
               tolerance = 1e-10)
  # The two groupings differ, so a mix-up of the two terms would show.
  expect_false(isTRUE(all.equal(dsl_mate_sum(pop, r$id_ind, "ind_meta"),
                                dsl_mate_sum(pop, r$id_ind, "pens"))))

  r <- dsl_records(pop, "Gc")
  expect_equal(r$pheno_value - r$residual_value,
               dsl_mate_sum(pop, r$id_ind, "pens", "additive") / 9,
               tolerance = 1e-10)
})

test_that("PH2: define_phenotype() refuses anything outside the DSL, before writing", {
  pop <- dsl_pop("ph2_bad")
  on.exit(close_pop(pop), add = TRUE)
  bad <- function(f, pattern) {
    expect_error(define_phenotype(pop, "W", residual_var = 1, formula_tgv = f),
                 pattern, label = f)
  }
  bad("dam(T, component = \"bogus\")", "`component = \"bogus\"` must be one of")
  bad("dam(T, component = 1)", "`component` must be a name or a single string")
  bad("dam(T, \"additive\")", "dam\\(\\) takes 1 positional argument")
  bad("self(T, table = \"pens\")", "self\\(\\) has no argument `table`")
  bad("dam(T, weight = 2)", "has no argument `weight`")
  bad("dam(dam(T))", "`trait` must be a name or a single string")
  bad("dam(T, component = \"additive\", component = \"total\")", "given twice")
  bad("group_sum(T, pen, \"ind_meta\")",
      "group_sum\\(\\) takes 2 positional arguments.*`table =` must be named")
  bad("group_sum(T)", "takes 2 positional arguments")
  bad("group_sum(T, \"pen-id\")", "Invalid group column 'pen-id'")
  bad("group_sum(T, pen, table = \"pens; DROP TABLE x\")", "Invalid group table")
  bad("group_sum(T, pen, table = \"nope\")", "table 'nope' does not exist")
  bad("group_sum(T, pen, table = \"genome_meta\")", "has no 'id_ind' column")
  bad("group_sum(T, litter)", "column 'litter' not found in table 'ind_meta'")
  bad("system(\"echo hi\")", "`system\\(\"echo hi\"\\)` is not allowed")
  bad("T + \"a\"", "the constant `\"a\"` is not allowed")
  bad("T + dam", "`dam` is a function name, not a trait")
  bad("T; T", "must be a single expression")
  bad("Q + dam(T)", "Unknown trait name\\(s\\)")
  expect_equal(phenotype_rows(pop), 0L)

  # The math whitelist and numbers still work (with the usual constant warning).
  expect_warning(
    pop <- define_phenotype(pop, "W", residual_var = 1,
                            formula_tgv = "sqrt(abs(T)) + 0.5 * dam(T)"),
    "Scalar arithmetic constant")
  expect_equal(phenotype_rows(pop), 1L)
})

test_that("3c.2b: Stage 1 evaluates contributor ids without SQL text, as an ind_meta filter would", {
  # Ids are never rendered into SQL, so an id with a quote and a large set
  # behave like any other: unknown and NA ids are dropped, repeats collapse.
  set.seed(3502)
  pop <- open_pop(pop_name = "tgv_ids", db_name = ":memory:") |>
    define_genome(n_loci = 12, n_chr = 1, chr_len_Mb = 50) |>
    define_founder_haplotypes(n_haplotypes = 60, method = "beta")
  pop <- pop |> get_table("founder_haplotypes") |>
    add_founders(n_males = 700, n_females = 700, line_name = "A")
  on.exit(close_pop(pop), add = TRUE)
  pop <- with_additive_target(pop, "T", 1)
  pop <- suppressMessages(pop |> get_table("genome_meta") |>
    define_additive_effects("T", warn_bounds = NULL))

  ids  <- dplyr::collect(get_table(pop, "ind_meta"))$id_ind
  pick <- ids[seq(1L, length(ids), by = 2L)]
  pop <- suppressMessages(.tgv_compute_ids(
    pop, c(rev(pick), pick[1:5], NA, "O'Brien_1", "x\"y"), "T"))
  via_ids <- DBI::dbGetQuery(pop$db_conn,
    "SELECT * FROM ind_tgv ORDER BY id_ind, component_name")
  expect_equal(sort(unique(via_ids$id_ind)), sort(pick))

  DBI::dbExecute(pop$db_conn, "DELETE FROM ind_tgv")
  pop <- suppressMessages(get_table(pop, "ind_meta") |>
    dplyr::filter(id_ind %in% !!pick) |> add_tgv("T"))
  via_filter <- DBI::dbGetQuery(pop$db_conn,
    "SELECT * FROM ind_tgv ORDER BY id_ind, component_name")
  expect_identical(via_ids[, c("id_ind", "trait_name", "component_name", "tgv_value")],
                   via_filter[, c("id_ind", "trait_name", "component_name", "tgv_value")])
})

# Every SQL statement DuckDB receives while `code` runs.
record_sql <- function(code) {
  rec <- new.env()
  rec$sql <- character()
  # do.call() passes the built expression: trace() quotes its `tracer`.
  suppressMessages(do.call(trace, list("dbSendQuery",
    signature = c("duckdb_connection", "character"),
    tracer = bquote(assign("sql", c(get("sql", envir = .(rec)), statement),
                           envir = .(rec))),
    where = asNamespace("DBI"), print = FALSE), quote = TRUE))
  on.exit(suppressMessages(untrace("dbSendQuery",
    signature = c("duckdb_connection", "character"),
    where = asNamespace("DBI"))), add = TRUE)
  force(code)
  rec$sql
}

test_that("3c.2b: no individual id appears in any SQL add_phenotype() sends", {
  # A filtered subset (the plan's ind_meta read), dams and group-mates (the
  # Stage 1 contributor sets), and every lookup after them.
  pop <- dsl_pop("ids_sql")
  on.exit(close_pop(pop), add = TRUE)
  pop <- define_phenotype(pop, "T", residual_var = 1)
  pop <- define_phenotype(pop, "W", residual_var = 1,
    formula_tgv = "T + dam(T, component = \"additive\") + group_sum(T, pen)")
  pop <- define_phenotype(pop, "C", residual_var = 1,
    components = tibble::tribble(~source_trait_name, ~contributor_type,
                                 "T",                "self",
                                 "T",                "sire"))
  ids <- dplyr::collect(get_table(pop, "ind_meta"))$id_ind
  set.seed(9)
  sql <- record_sql(pop <- suppressMessages(dsl_offspring(pop) |>
    add_phenotype(c("T", "W", "C"))))
  expect_gt(length(sql), 20L)
  expect_equal(nrow(dsl_records(pop, "W")), 30L)
  quoted <- paste0("'", ids, "'")
  leaked <- quoted[vapply(quoted, function(q) any(grepl(q, sql, fixed = TRUE)),
                          logical(1))]
  expect_identical(leaked, character(0))
})

test_that("derived formula: only phenotype names, numbers, operators and math functions", {
  pop <- dsl_pop("derived_grammar")
  on.exit(close_pop(pop), add = TRUE)
  pop <- define_phenotype(pop, "T", residual_var = 1)
  bad <- function(f, pattern) {
    expect_error(define_phenotype(pop, "D", type = "derived_formula", formula = f),
                 pattern, label = f)
  }
  bad("T + nchar(system(\"echo hi\", intern = TRUE))",
      "`nchar\\(system.*` is not allowed")
  bad("get(\"T\")", "`get\\(\"T\"\\)` is not allowed")
  bad("base::sqrt(T)", "is not allowed")
  bad("T + \"a\"", "the constant `\"a\"` is not allowed")
  bad("T + TRUE", "the constant `TRUE` is not allowed")
  bad("sqrt + T", "`sqrt` is not a phenotype name")
  bad("T; T", "must be a single expression")
  expect_equal(phenotype_rows(pop), 1L)
  pop <- define_phenotype(pop, "D", type = "derived_formula",
                          formula = "round(sqrt(abs(T)) * 2, 1)")
  expect_equal(phenotype_rows(pop), 2L)

  # A stored formula that skipped define_phenotype() is refused before eval:
  # the command never runs and no record is written.
  DBI::dbExecute(pop$db_conn, paste0(
    "UPDATE phenotype_meta SET formula = 'T + nchar(Sys.setenv(TB_PWNED = \"1\"))' ",
    "WHERE phenotype_name = 'D'"))
  Sys.unsetenv("TB_PWNED")
  set.seed(10)
  pop <- suppressMessages(dsl_offspring(pop) |> add_phenotype("T"))
  expect_error(suppressMessages(dsl_offspring(pop) |> add_phenotype("D")),
               "is not allowed")
  expect_identical(Sys.getenv("TB_PWNED"), "")
  expect_equal(nrow(dsl_records(pop, "D")), 0L)
})
