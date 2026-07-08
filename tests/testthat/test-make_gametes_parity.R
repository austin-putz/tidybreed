# Stage 3b: R <-> C++ recombination-kernel parity.
# Locks (a) the compiled kernel reproduces the R reference bit-for-bit for a
# fixed base seed across configs and store_crossovers, (b) wrapper-level parity
# through add_offspring() including a configured sex chromosome (C++ autosomes +
# R special path), and (c) the runtime kernel selector.

skip_if_not(exists("make_gametes_batch_cpp", mode = "function"),
            "compiled kernel unavailable")

# ---- direct kernel parity --------------------------------------------------

# Build packed kernel inputs for n_par parents, n_loci autosome loci laid out in
# `n_chr` contiguous chromosomes, and G gametes with random parent/o/origin.
make_kernel_inputs <- function(n_par, n_chr, loci_per_chr, G, seed) {
  set.seed(seed)
  n_loci <- n_chr * loci_per_chr
  parent_allele  <- matrix(sample(0:1, 2 * n_par * n_loci, TRUE), nrow = 2 * n_par)
  parent_lo_code <- matrix(sample(1:4, 2 * n_par * n_loci, TRUE), nrow = 2 * n_par)
  chr_start <- as.integer(seq(1L, n_loci, by = loci_per_chr))
  chr_end   <- as.integer(chr_start + loci_per_chr - 1L)
  chr_len_cM <- runif(n_chr, 40, 300)          # exercises varied lambda incl. > 1 Morgan
  chr_pos_cM <- numeric(n_loci)
  for (c in seq_len(n_chr)) {
    idx <- chr_start[c]:chr_end[c]
    chr_pos_cM[idx] <- sort(runif(loci_per_chr, 0, chr_len_cM[c]))
  }
  list(parent_allele = parent_allele, parent_lo_code = parent_lo_code,
       gamete_parent_idx = sample(seq_len(n_par), G, TRUE),
       gamete_o = sample(1:100000, G), gamete_origin = sample(1:2, G, TRUE),
       chr_start = chr_start, chr_end = chr_end,
       chr_pos_cM = chr_pos_cM, chr_len_cM = chr_len_cM)
}

call_kernel <- function(fn, k, base_seed, store) {
  fn(k$parent_allele, k$parent_lo_code, k$gamete_parent_idx, k$gamete_o,
     k$gamete_origin, k$chr_start, k$chr_end, k$chr_pos_cM, k$chr_len_cM,
     base_seed, store)
}

test_that("C++ kernel reproduces the R reference bit-for-bit", {
  configs <- list(
    list(n_par = 6, n_chr = 3, lpc = 20, G = 12),
    list(n_par = 4, n_chr = 1, lpc = 50, G = 8),    # single chromosome
    list(n_par = 8, n_chr = 5, lpc = 8,  G = 20)    # many short chromosomes
  )
  for (cf in configs) {
    k <- make_kernel_inputs(cf$n_par, cf$n_chr, cf$lpc, cf$G, seed = cf$G)
    for (base_seed in c(1L, 42L, 7L, 2147483647L)) {
      for (store in c(FALSE, TRUE)) {
        r_out   <- call_kernel(make_gametes_batch_r,   k, base_seed, store)
        cpp_out <- call_kernel(make_gametes_batch_cpp, k, base_seed, store)
        for (nm in c("parent_origin", "locus_idx", "allele", "line_origin_code",
                     "xover_gamete_o", "xover_parent_origin", "xover_chr")) {
          expect_identical(as.integer(r_out[[nm]]), as.integer(cpp_out[[nm]]),
                           info = paste(nm, "seed", base_seed, "store", store))
        }
        expect_identical(as.numeric(r_out$xover_pos_cM),
                         as.numeric(cpp_out$xover_pos_cM),
                         info = paste("xover_pos_cM seed", base_seed, "store", store))
      }
    }
  }
})

test_that("lambda above the ceiling errors identically in both kernels", {
  k <- make_kernel_inputs(2, 1, 10, 2, seed = 1)
  k$chr_len_cM <- 3100                          # 31 Morgans > ceiling
  expect_error(call_kernel(make_gametes_batch_r,   k, 1L, FALSE), "ceiling")
  expect_error(call_kernel(make_gametes_batch_cpp, k, 1L, FALSE), "ceiling")
})

# ---- wrapper-level parity through add_offspring() --------------------------

test_that("add_offspring() is identical under C++ and forced-R, incl. a sex chromosome", {
  build <- function() {
    set.seed(2024)
    pop <- open_pop(pop_name = "kp", db_name = ":memory:") |>
      define_genome(n_loci = 180, n_chr = 3, chr_len_Mb = 100)
    # Make chromosome 3 an X: males hemizygous (dam-supplied), females diploid.
    x_name <- DBI::dbGetQuery(pop$db_conn,
      "SELECT DISTINCT chr_name FROM genome_meta WHERE chr = 3")$chr_name
    pop <- define_chr(pop, x_name, copy_mode_M = "half", copy_mode_F = "full",
                      hemi_parent = "parent_2", recombines = TRUE)
    pop |>
      define_founder_haplotypes(n_haplotypes = 30, line_name = "A", method = "uniform") |>
      get_table("founder_haplotypes") |>
      add_founders(n_males = 5, n_females = 5, line_name = "A")
  }
  mat <- tibble::tibble(
    id_parent_1 = rep(paste0("A_", 1:5), 4),
    id_parent_2 = rep(paste0("A_", 6:10), 4),
    sex = rep(c("M", "F"), 10), line_name = "A")

  run <- function(kernel) {
    withr::local_options(tidybreed.kernel = kernel)
    p <- build()
    p <- suppressMessages(add_offspring(p, mat, seed = 314L, store_crossovers = TRUE))
    h  <- DBI::dbGetQuery(p$db_conn,
      "SELECT * FROM ind_haplotype ORDER BY id_ind, parent_origin, locus_id")
    xc <- DBI::dbGetQuery(p$db_conn,
      "SELECT id_ind, parent_origin, chr, chr_name, pos_cM FROM ind_crossover
       ORDER BY id_ind, parent_origin, chr, pos_cM")
    close_pop(p); list(h = h, xc = xc)
  }
  cpp_run <- run("auto"); r_run <- run("r")

  expect_equal(cpp_run$h, r_run$h)
  expect_equal(cpp_run$xc, r_run$xc)
  expect_gt(nrow(cpp_run$xc), 0L)               # crossovers actually stored
})

# ---- runtime selector ------------------------------------------------------

test_that("the kernel selector honours option / env and defaults to C++", {
  args <- make_kernel_inputs(3, 2, 10, 6, seed = 5)
  ref  <- call_kernel(make_gametes_batch_r, args, 99L, TRUE)

  withr::with_options(list(tidybreed.kernel = "r"), {
    out <- call_kernel(make_gametes_batch, args, 99L, TRUE)
    expect_identical(as.integer(out$allele), as.integer(ref$allele))
  })
  withr::with_options(list(tidybreed.kernel = "auto"), {
    out <- call_kernel(make_gametes_batch, args, 99L, TRUE)   # C++ path
    expect_identical(as.integer(out$allele), as.integer(ref$allele))
  })
})
