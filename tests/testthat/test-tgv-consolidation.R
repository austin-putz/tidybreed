# Step 3a (plans/import_qtl_effect_methods.md §6, gates T3-T9): ind_tgv is the
# one table of true genetic values. The breeding value is its 'additive'
# component; there is no second table and no second evaluator.

# 40 founders, 12 loci; generated additive effects at loci 1-8 (target 1).
cons_pop <- function(name, db_name = ":memory:") {
  set.seed(3301)
  pop <- open_pop(pop_name = name, db_name = db_name) |>
    define_genome(n_loci = 12, n_chr = 1, chr_len_Mb = 50) |>
    define_founder_haplotypes(n_haplotypes = 60, method = "beta")
  pop <- pop |> get_table("founder_haplotypes") |>
    add_founders(n_males = 20, n_females = 20, line_name = "A")
  pop <- with_additive_target(pop, "T", 1)
  suppressMessages(pop |> get_table("genome_meta") |>
    dplyr::filter(locus_id <= 8L) |>
    define_additive_effects("T", warn_bounds = NULL))
}

cons_loci <- function(pop) {
  pop |> get_table("genome_meta") |> dplyr::collect() |>
    dplyr::arrange(.data$locus_id) |> dplyr::pull("locus_name")
}

# Add one term of every other kind, through the writer, under user owners:
# Cockerham dominance at loci 2-4, an A x A pair at loci 9-10, and a one-locus
# indicator surface at locus 11.
cons_add_mixed <- function(pop) {
  loci <- cons_loci(pop)
  p <- extract_allele_freq(get_table(pop, "founder_haplotypes"))
  dom <- loci[2:4]
  pop <- define_genome_effect_terms(pop, "T",
    ad_terms(dom, a = 0, d = c(0.6, -0.4, 0.3),
             p = p$allele_freq[match(dom, p$locus_name)],
             coding = "cockerham", report = FALSE),
    effect_owner = "dom")
  pop <- define_genome_effect_terms(pop, "T", data.frame(
    term_id = c(1L, 1L), locus_name = loci[9:10], contrast_name = "additive",
    center_value = 0.5, genome_value = 0.7), effect_owner = "epi")
  define_genome_effect_terms(pop, "T",
    genotype_terms(stats::setNames(data.frame(c(0L, 2L)), loci[11]),
                   value = c(-0.5, 0.8)),
    effect_owner = "surf")
}

# Independent oracle: every term of a trait evaluated from ind_haplotype and the
# stored rows, in R, one individual at a time. Common-scope terms only (no
# origin rows), which is all these fixtures write.
cons_oracle <- function(pop, trait) {
  conn <- pop$db_conn
  hap <- DBI::dbGetQuery(conn,
    "SELECT id_ind, locus_id, allele FROM ind_haplotype")
  mem <- DBI::dbGetQuery(conn, paste0(
    "SELECT e.id_genome_effect, e.genome_value, m.member_slot, m.locus_id, ",
    "m.contrast_name, m.center_value, m.copy_count_value, m.dosage_value ",
    "FROM genome_effects e JOIN genome_effect_members m USING (id_genome_effect) ",
    "WHERE e.trait_name = '", trait, "'"))
  ids <- sort(unique(hap$id_ind))
  out <- list()
  for (id in ids) {
    h <- hap[hap$id_ind == id, ]
    for (g in split(mem, mem$id_genome_effect)) {
      x <- vapply(seq_len(nrow(g)), function(k) {
        a <- h$allele[h$locus_id == g$locus_id[k]]
        c0 <- g$center_value[k]
        switch(g$contrast_name[k],
          additive  = sum(a - c0),
          dominance = c(-2 * c0^2, 2 * c0 * (1 - c0),
                        -2 * (1 - c0)^2)[sum(a) + 1],
          indicator = as.numeric(length(a) == g$copy_count_value[k] &&
                                   sum(a) == g$dosage_value[k]))
      }, numeric(1))
      comp <- if (nrow(g) > 1L) "interaction" else g$contrast_name[1]
      out[[length(out) + 1L]] <- data.frame(
        id_ind = id, component_name = comp, value = g$genome_value[1] * prod(x))
    }
  }
  out <- do.call(rbind, out)
  stats::aggregate(value ~ id_ind + component_name, data = out, FUN = sum)
}


test_that("T3: an additive-only model writes one 'additive' row, equal to the total", {
  pop <- cons_pop("cons_t3")
  on.exit(close_pop(pop), add = TRUE)
  pop <- suppressMessages(pop |> get_table("ind_meta") |> add_tgv("T"))

  tg <- dplyr::collect(get_table(pop, "ind_tgv"))
  expect_setequal(unique(tg$component_name), "additive")
  expect_equal(nrow(tg), 40L)
  tot <- dplyr::collect(get_table(pop, "ind_tgv_total"))
  # The exact (DECIMAL) total of one row is that row, bit for bit.
  expect_identical(tot$tgv_total[match(tg$id_ind, tot$id_ind)], tg$tgv_value)
})

test_that("T4: every component of a mixed model equals the test's own computation", {
  pop <- cons_add_mixed(cons_pop("cons_t4"))
  on.exit(close_pop(pop), add = TRUE)
  pop <- suppressMessages(pop |> get_table("ind_meta") |> add_tgv("T"))

  tg <- dplyr::collect(get_table(pop, "ind_tgv"))
  expect_setequal(unique(tg$component_name), TGV_COMPONENT_NAMES)
  want <- cons_oracle(pop, "T")
  key  <- function(d) paste(d$id_ind, d$component_name)
  expect_equal(tg$tgv_value[match(key(want), key(tg))], want$value,
               tolerance = 1e-12)

  tot <- dplyr::collect(get_table(pop, "ind_tgv_total"))
  sums <- tapply(want$value, want$id_ind, sum)
  expect_equal(tot$tgv_total[match(names(sums), tot$id_ind)], as.numeric(sums),
               tolerance = 1e-12)

  # The consumer reader aggregates the selected individuals only, and must
  # agree with the view bit for bit.
  rd <- .tgv_read(pop$db_conn, tot$id_ind, "T", "total")
  expect_identical(rd$value[match(tot$id_ind, rd$id_ind)], tot$tgv_total)
})

test_that("re-evaluation keeps custom columns on surviving rows", {
  pop <- cons_add_mixed(cons_pop("cons_custom"))
  on.exit(close_pop(pop), add = TRUE)
  pop <- suppressMessages(
    pop |> get_table("ind_meta") |> add_tgv("T", run_label = "first"))
  pop <- suppressMessages(pop |> get_table("ind_meta") |> add_tgv("T"))
  tg <- dplyr::collect(get_table(pop, "ind_tgv"))
  expect_equal(nrow(tg), 4L * 40L)
  expect_true(all(tg$run_label == "first"))
  # A new value for the same column replaces it.
  pop <- suppressMessages(
    pop |> get_table("ind_meta") |> add_tgv("T", run_label = "second"))
  expect_true(all(dplyr::collect(get_table(pop, "ind_tgv"))$run_label == "second"))
})

test_that("T7: true indices on the breeding value and on the total coexist", {
  pop <- cons_add_mixed(cons_pop("cons_t7"))
  on.exit(close_pop(pop), add = TRUE)
  pop <- with_additive_target(pop, "U", 1)
  pop <- suppressMessages(pop |> get_table("genome_meta") |>
    dplyr::filter(locus_id >= 5L) |> define_additive_effects("U", warn_bounds = NULL))
  pop <- define_index(pop, "sel", c("T", "U"), index_wts = c(2, -1))

  # T's hand-written indicator and interaction terms: the additive index
  # warns that its 'additive' component is not the breeding value.
  expect_warning(pop <- suppressMessages(pop |> get_table("ind_meta") |>
    add_tgv(c("T", "U"), index_names = "sel")),
    "'additive' component is not the breeding value")
  pop <- suppressMessages(pop |> get_table("ind_meta") |>
    add_tgv(c("T", "U"), index_names = "sel", component_name = "total"))

  ti <- dplyr::collect(get_table(pop, "ind_true_index"))
  expect_equal(nrow(ti), 80L)
  expect_setequal(unique(ti$component_name), c("additive", "total"))

  a   <- tgv_additive(pop)
  tot <- dplyr::collect(get_table(pop, "ind_tgv_total"))
  ids <- sort(unique(a$id_ind))
  val <- function(d, col, tr) {
    d <- d[d$trait_name == tr, ]
    d[[col]][match(ids, d$id_ind)]
  }
  want_a <- 2 * val(a, "tgv_value", "T") - val(a, "tgv_value", "U")
  want_t <- 2 * val(tot, "tgv_total", "T") - val(tot, "tgv_total", "U")
  ti_a <- ti[ti$component_name == "additive", ]
  ti_t <- ti[ti$component_name == "total", ]
  expect_equal(ti_a$true_index_value[match(ids, ti_a$id_ind)], want_a,
               tolerance = 1e-12)
  expect_equal(ti_t$true_index_value[match(ids, ti_t$id_ind)], want_t,
               tolerance = 1e-12)
  expect_false(isTRUE(all.equal(want_a, want_t)))

  # A second additive run skips (overwrite_index = FALSE) and leaves the total.
  pop <- suppressMessages(pop |> get_table("ind_meta") |>
    add_tgv(c("T", "U"), index_names = "sel"))
  expect_equal(nrow(dplyr::collect(get_table(pop, "ind_true_index"))), 80L)

  expect_error(pop |> get_table("ind_meta") |>
                 add_tgv("T", index_names = "sel", component_name = "bogus"),
               "`component_name` must be one of")
})

test_that("T8: archive, remove_rows and restore_pop know ind_tgv, not the old table", {
  tmp <- tempfile(fileext = ".duckdb")
  arc <- tempfile(fileext = ".duckdb")
  on.exit(unlink(c(tmp, arc)), add = TRUE)
  pop <- cons_pop("cons_t8", db_name = tmp)
  pop <- suppressMessages(pop |> get_table("ind_meta") |> add_tgv("T"))

  # remove_rows() on ind_tgv deletes exactly the filtered rows.
  one <- sort(dplyr::collect(get_table(pop, "ind_meta"))$id_ind)[1]
  pop <- get_table(pop, "ind_tgv") |> dplyr::filter(id_ind == !!one) |>
    remove_rows(verbose = FALSE)
  expect_equal(nrow(dplyr::collect(get_table(pop, "ind_tgv"))), 39L)

  pop <- suppressMessages(archive_replicate(pop, replicate = 1L,
                                            archive_path = arc))
  arc_conn <- DBI::dbConnect(duckdb::duckdb(), dbdir = arc, read_only = TRUE)
  expect_true("replicate" %in% DBI::dbListFields(arc_conn, "ind_tgv"))
  expect_false("ind_tbv" %in% DBI::dbListTables(arc_conn))
  DBI::dbDisconnect(arc_conn, shutdown = TRUE)
  close_pop(pop)

  # A current file restores.
  pop <- restore_pop(tmp)
  close_pop(pop)

  run_sql <- function(sql) {
    conn <- DBI::dbConnect(duckdb::duckdb(), dbdir = tmp)
    DBI::dbExecute(conn, sql)
    DBI::dbDisconnect(conn, shutdown = TRUE)
  }
  # 0.74.3 (Q18): the formula_tbv column became formula_tgv.
  run_sql("ALTER TABLE phenotype_meta RENAME COLUMN formula_tgv TO formula_tbv")
  expect_error(restore_pop(tmp), "pre-v0\\.74\\.3 'phenotype_meta' shape")
  run_sql("ALTER TABLE phenotype_meta RENAME COLUMN formula_tbv TO formula_tgv")
  run_sql("INSERT INTO ind_tgv VALUES (999, 'X_1', 'T', 'order1_additive', 1.0)")
  expect_error(restore_pop(tmp), "pre-v0\\.74\\.0 component names in ind_tgv")
  run_sql("DELETE FROM ind_tgv WHERE id_tgv = 999")
  run_sql("ALTER TABLE ind_true_index DROP COLUMN component_name")
  expect_error(restore_pop(tmp), "pre-v0\\.74\\.0 'ind_true_index' shape")
  run_sql("CREATE TABLE ind_tbv (id_ind VARCHAR)")
  expect_error(restore_pop(tmp), "pre-v0\\.74\\.0 'ind_tbv'")
})

test_that("T9: schema() lists ind_tgv and ind_true_index.component_name, and no TBV table", {
  pop <- cons_pop("cons_t9")
  on.exit(close_pop(pop), add = TRUE)
  s <- suppressMessages(schema(pop, show_empty = TRUE))
  expect_true(all(c("ind_tgv", "ind_tgv_total", "ind_true_index") %in% s$table_name))
  expect_false(any(grepl("tbv", s$table_name)))
  d <- describe_table(pop, "ind_true_index")
  expect_true("component_name" %in% d$column_name)
  expect_true(all(nzchar(d$description[d$column_name == "component_name"])))
})

test_that("an unknown index is refused before anything is written", {
  pop <- cons_pop("cons_noindex")
  on.exit(close_pop(pop), add = TRUE)
  expect_error(pop |> get_table("ind_meta") |> add_tgv("T", index_names = "nope"),
               "Index 'nope' not found")
  expect_equal(nrow(dplyr::collect(get_table(pop, "ind_tgv"))), 0L)
})
