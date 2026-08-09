# v5 review #2/#3: the resolved copy total must equal the number of
# ind_haplotype rows actually emitted, and a storage-expressible-but-unimplemented
# transmission (uniparental disomy, 2,0) must fail explicitly at the kernel
# boundary rather than silently duplicating or dropping strands.

.rows_for <- function(pop, id, chr) {
  DBI::dbGetQuery(pop$db_conn, sprintf(
    "SELECT COUNT(*) AS n FROM ind_haplotype h JOIN genome_meta gm USING(locus_id)
     WHERE h.id_ind = '%s' AND gm.chr_name = '%s'", id, chr))$n
}

test_that("emitted ind_haplotype rows equal the resolved copy total (XY genome)", {
  set.seed(4001)
  suppressMessages({
    pop <- open_pop(pop_name = "kb_xy", db_name = ":memory:") |>
      define_genome(n_loci = 15, n_chr = 3, chr_names = c("1", "X", "Y"),
                    chr_len_Mb = 50) |>
      define_chromosome("X", offspring_sex = "M", from_parent_1 = 0, from_parent_2 = 1) |>
      define_chromosome("Y", offspring_sex = "M", from_parent_1 = 1, from_parent_2 = 0) |>
      define_chromosome("Y", offspring_sex = "F", from_parent_1 = 0, from_parent_2 = 0) |>
      define_chromosome("Y", recombines = FALSE) |>
      define_founder_haplotypes(n_haplotypes = 20, method = "fixed") |>
      get_table("founder_haplotypes") |>
      add_founders(n_males = 3, n_females = 3, line_name = "A")
  })
  on.exit(close_pop(pop), add = TRUE)
  conn <- pop$db_conn
  n_per_chr <- 5L   # 15 loci / 3 chr

  inh_M <- resolve_chr_inheritance(conn, "M")
  inh_F <- resolve_chr_inheritance(conn, "F")
  tot <- function(inh, chr) {
    r <- inh[inh$chr_name == chr, ]
    (r$from_parent_1 + r$from_parent_2) * n_per_chr
  }

  ids <- DBI::dbGetQuery(conn, "SELECT id_ind, sex FROM ind_meta")
  for (i in seq_len(nrow(ids))) {
    id <- ids$id_ind[i]; sx <- ids$sex[i]
    inh <- if (sx == "M") inh_M else inh_F
    for (chr in c("1", "X", "Y")) {
      expect_equal(.rows_for(pop, id, chr), tot(inh, chr),
                   info = sprintf("id=%s chr=%s sex=%s", id, chr, sx))
    }
  }
})

test_that("a stored 2,0 (uniparental disomy) fails explicitly in add_founders", {
  suppressMessages({
    pop <- open_pop(pop_name = "kb_20", db_name = ":memory:") |>
      define_genome(n_loci = 10, n_chr = 2, chr_len_Mb = 50)
    # 2,0 is storage-expressible (sum <= 2) ...
    pop <- define_chromosome(pop, "1", offspring_sex = "M",
                             from_parent_1 = 2, from_parent_2 = 0)
    pop <- define_founder_haplotypes(pop, n_haplotypes = 10, method = "fixed")
  })
  on.exit(close_pop(pop), add = TRUE)

  # ... but the founder kernel has no mechanism for two child copies from one
  # parent slot, so it errors.
  expect_error(
    pop |> get_table("founder_haplotypes") |>
      add_founders(n_males = 2, n_females = 2, line_name = "A"),
    "not implemented")
})

test_that("a stored 2,0 fails explicitly in add_offspring", {
  set.seed(4002)
  suppressMessages({
    pop <- open_pop(pop_name = "kb_20o", db_name = ":memory:") |>
      define_genome(n_loci = 10, n_chr = 2, chr_len_Mb = 50) |>
      define_founder_haplotypes(n_haplotypes = 10, method = "fixed") |>
      get_table("founder_haplotypes") |>
      add_founders(n_males = 2, n_females = 2, line_name = "A")
    # Introduce the unimplemented rule AFTER plain-autosome founders exist.
    pop <- define_chromosome(pop, "1", offspring_sex = "M",
                             from_parent_1 = 2, from_parent_2 = 0)
  })
  on.exit(close_pop(pop), add = TRUE)

  ids <- DBI::dbGetQuery(pop$db_conn, "SELECT id_ind, sex FROM ind_meta")
  matings <- tibble::tibble(
    id_parent_1 = ids$id_ind[ids$sex == "M"][1],
    id_parent_2 = ids$id_ind[ids$sex == "F"][1],
    sex = "M", line_name = "F1")
  expect_error(add_offspring(pop, matings), "not implemented")
})
