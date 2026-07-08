#' Pre-compute per-chromosome locus metadata for gamete simulation
#'
#' Called once per `add_offspring()` invocation. The returned list is passed to
#' every `make_gamete()` call, avoiding O(n_loci × n_chr) masking work inside
#' the per-offspring loop.
#'
#' @param locus_map Data frame with columns `chr` and `pos_cM` — the resolved
#'   genetic map (from `resolve_genome_map()`), one row per locus, ordered by
#'   `locus_id`.
#' @return Named list, one element per chromosome. Each element is a list with
#'   `locus_idx` (integer indices into the locus dimension), `pos_cM`
#'   (genetic positions in centiMorgans), and `chr_len` (derived as
#'   `MAX(pos_cM)` for the chromosome; crossovers beyond the last locus are
#'   invisible to inheritance).
#' @keywords internal
build_chr_info <- function(locus_map) {
  chrs <- sort(unique(locus_map$chr))
  stats::setNames(
    lapply(chrs, function(chr_id) {
      mask <- locus_map$chr == chr_id
      list(
        locus_idx = which(mask),
        pos_cM    = locus_map$pos_cM[mask],
        chr_len   = max(locus_map$pos_cM[mask])
      )
    }),
    as.character(chrs)
  )
}


# ---------------------------------------------------------------------------
# dqrng per-gamete stream model (Stage 2)
#
# The unit of randomness is a gamete identified by (offspring o, parent_origin
# r, stream_kind k): autosomes use k = 1, special chromosomes use k = 2 (an
# independent sub-stream, so the R special path never has to know where the
# autosome kernel left the counter). Each gamete's stream is seeded ONCE from
# the base seed and the folded stream id `sid`, then every draw for that gamete
# is consumed sequentially from it. This makes output independent of batch size,
# iteration order, and (later) thread count — gamete (o, r, k) always consumes
# stream (o, r, k) from its start. Scheme A from dev/spikes/dqrng_parity/: a
# 2-int seed `c(base_seed, sid)` folded by dqrng's own convert_seed, so the R
# reference here and the future C++ kernel seed bit-identically.
# ---------------------------------------------------------------------------

#' Fold (offspring, parent_origin, stream_kind) into a single dqrng stream id
#'
#' @param o Integer offspring index (the original `matings` row index; global,
#'   never batch-local — so any `batch_size` yields identical output).
#' @param parent_origin Integer 1 (parent_1/sire) or 2 (parent_2/dam).
#' @param kind Integer 1 (autosome) or 2 (special chromosome).
#' @return Integer stream id `sid`. Fits signed 32-bit for `o < ~5e8`.
#' @keywords internal
.gamete_stream_id <- function(o, parent_origin, kind) {
  as.integer((o * 2L + (parent_origin - 1L)) * 2L + (kind - 1L))
}

#' Seed the global dqrng stream for one gamete (Scheme A)
#'
#' @param base_seed Integer base seed for the whole `add_offspring()` call.
#' @param sid Integer stream id from [.gamete_stream_id()].
#' @return Invisibly `NULL`. Sets the dqrng generator kind + seed as a side
#'   effect; subsequent `dqrng::dqrunif()` draws consume this gamete's stream.
#' @keywords internal
.seed_gamete_stream <- function(base_seed, sid) {
  dqrng::dqRNGkind("Xoroshiro128++")
  dqrng::dqset.seed(c(as.integer(base_seed), as.integer(sid)))
  invisible(NULL)
}

#' Draw a Poisson count from the current dqrng stream (log-accumulation inversion)
#'
#' Knuth's multiplicative inversion in log space: accumulate `log(u)` and count
#' steps until the sum drops to `-lambda`. Consumes ONLY dqrng uniforms (no base
#' `rpois`, no `std::poisson_distribution`), so the R reference and the future
#' C++ kernel match by construction from the shared uniform stream. Log form
#' (not the raw product) avoids underflow when `exp(-lambda)` is denormal.
#'
#' @param lambda Non-negative Poisson rate. `lambda <= 0` returns `0L` and draws
#'   nothing. Errors above a documented ceiling (30 Morgans / 3000 cM) rather
#'   than spinning `O(lambda)` uniforms.
#' @return Integer crossover count.
#' @keywords internal
.draw_poisson_dqrng <- function(lambda) {
  if (lambda <= 0) return(0L)
  if (lambda > 30) {
    stop(
      "Recombination rate lambda = ", signif(lambda, 4),
      " exceeds the supported ceiling (30 Morgans / 3000 cM per chromosome). ",
      "Check cM_per_Mb / the genetic map — a per-chromosome genetic length above ",
      "3000 cM is almost certainly a units error.",
      call. = FALSE
    )
  }
  target <- -lambda
  s <- 0
  k <- 0L
  repeat {
    s <- s + log(dqrng::dqrunif(1L))
    if (s <= target) break
    k <- k + 1L
  }
  k
}

#' Draw one chromosome's recombination pattern from the current dqrng stream
#'
#' The per-chromosome core shared by autosomes ([make_gamete()]) and recombining
#' special chromosomes (the `add_offspring()` special path). Assumes the gamete's
#' stream is ALREADY seeded — draws sequentially in the fixed canonical order
#' (count, then start homolog, then sorted positions) so R↔C++ parity holds.
#' Purely a homolog-index computation: it never touches allele values, so the
#' caller gathers alleles/`line_origin` from its own matrices via `hap_idx_vec`.
#'
#' @param chr_pos Numeric locus genetic positions (cM) for this chromosome.
#' @param chr_len Numeric chromosome genetic length (cM) = `max(chr_pos)`.
#' @param store_crossovers Logical; when `TRUE` also return the drawn crossover
#'   positions (cM) for `ind_crossover`.
#' @return List with `hap_idx_vec` (length `length(chr_pos)`, values in `{1,2}`
#'   — which homolog donates each locus) and `xover_pos_cM` (sorted crossover
#'   positions in cM, or `numeric(0)`).
#' @keywords internal
.draw_chr_recombination <- function(chr_pos, chr_len, store_crossovers = FALSE) {
  n_cross <- .draw_poisson_dqrng(chr_len / 100)
  start   <- if (dqrng::dqrunif(1L) < 0.5) 1L else 2L
  if (n_cross == 0L) {
    hap_idx_vec <- rep(start, length(chr_pos))
    cross_pos   <- numeric(0)
  } else {
    cross_pos   <- sort(dqrng::dqrunif(n_cross, min = 0, max = chr_len))
    # findInterval counts crossovers strictly before each locus position; parity
    # of that count determines how many times the haplotype toggled.
    n_toggles   <- findInterval(chr_pos, cross_pos)
    hap_idx_vec <- (start - 1L + n_toggles %% 2L) %% 2L + 1L
  }
  list(hap_idx_vec  = hap_idx_vec,
       xover_pos_cM = if (store_crossovers) cross_pos else numeric(0))
}

#' Produce a gamete from a parent's two haplotypes via chromosomal crossovers
#'
#' Simulates chromosomal recombination using the Haldane map function, drawing
#' from this gamete's own dqrng stream (seeded here from `base_seed` + `sid`).
#' Crossover count per chromosome ~ Poisson(chr_len_cM / 100), i.e.
#' Poisson(genetic length in Morgans); positions uniform in genetic (cM)
#' distance. All chromosomes of the gamete are drawn sequentially from the one
#' seeded stream in ascending-chromosome order.
#'
#' @param hap_matrix 2 x n_loci integer matrix. Row 1 = haplotype from
#'   parent_origin 1, row 2 = parent_origin 2.
#' @param chr_info List from `build_chr_info()` — pre-computed locus indices,
#'   positions, and lengths per chromosome.
#' @param base_seed Integer base seed for the `add_offspring()` call.
#' @param sid Integer stream id (see [.gamete_stream_id()], kind = 1 autosome).
#' @param store_crossovers Logical; emit the crossover buffer when `TRUE`.
#' @return List: `allele` and `homolog` (length-`n_loci` integer vectors, the
#'   gamete alleles and the donating homolog per locus), plus `xover_chr_idx`
#'   (the `chr_info` position of each crossover event) and `xover_pos_cM` (its
#'   cM position) — both `length(0)` unless `store_crossovers`.
#' @keywords internal
make_gamete <- function(hap_matrix, chr_info, base_seed, sid,
                        store_crossovers = FALSE) {
  .seed_gamete_stream(base_seed, sid)

  n_loci  <- ncol(hap_matrix)
  gamete  <- integer(n_loci)
  homolog <- integer(n_loci)
  xchr    <- integer(0)
  xpos    <- numeric(0)

  for (ci_i in seq_along(chr_info)) {
    ci  <- chr_info[[ci_i]]
    idx <- ci$locus_idx
    d   <- .draw_chr_recombination(ci$pos_cM, ci$chr_len, store_crossovers)

    gamete[idx]  <- hap_matrix[cbind(d$hap_idx_vec, idx)]
    homolog[idx] <- d$hap_idx_vec

    if (store_crossovers && length(d$xover_pos_cM)) {
      xchr <- c(xchr, rep.int(ci_i, length(d$xover_pos_cM)))
      xpos <- c(xpos, d$xover_pos_cM)
    }
  }

  list(allele = gamete, homolog = homolog,
       xover_chr_idx = xchr, xover_pos_cM = xpos)
}

#' Generate autosome gametes for a batch of offspring (the swappable seam)
#'
#' The single seam the C++ kernel (Stage 3) drops into. Dependency-free — no
#' `tidybreed_pop`, no DBI, no `data.frame` in its numeric inputs — so it is
#' drop-in swappable and unit-testable in isolation.
#'
#' **RNG order:** offspring in the given order, **parent_1 (sire) gamete then
#' parent_2 (dam) gamete**; each gamete draws from its own dqrng stream keyed on
#' its global offspring index `o_index[i]` (autosome kind = 1), so output is
#' independent of batch size and iteration order.
#'
#' @param sire_ids,dam_ids Character vectors of parent IDs, length `n`; index
#'   into `parent_haps`/`parent_lo`.
#' @param parent_haps Named list; `parent_haps[[pid]]` is the parent's
#'   `2 x n_autosome_loci` integer allele matrix (homolog x locus).
#' @param parent_lo Named list; parallel `2 x n_autosome_loci` `line_origin`
#'   character matrix.
#' @param autosome_info_by_parent Named list; `[[pid]]` is the parent's resolved
#'   autosome `chr_info` (LOCAL indices) for its `(sex, line)` map.
#' @param n_autosome_loci Integer autosome locus count (matrix column count).
#' @param base_seed Integer base seed for the `add_offspring()` call.
#' @param o_index Integer vector length `n` — the global offspring index `o`
#'   (original `matings` row) per batch element; drives the per-gamete stream.
#' @param store_crossovers Logical; collect the crossover buffer when `TRUE`.
#' @return List of four `n x n_autosome_loci` matrices (`sire_allele`,
#'   `dam_allele`, `sire_lo`, `dam_lo`) plus `crossovers`: `NULL`, or a
#'   data.frame `(local_j, parent_origin, chr_idx, pos_cM)` (`chr_idx` indexes
#'   into the autosome `chr_info` order; the caller maps it to `chr`/`chr_name`).
#' @keywords internal
make_gametes_batch <- function(sire_ids, dam_ids, parent_haps, parent_lo,
                               autosome_info_by_parent, n_autosome_loci,
                               base_seed, o_index, store_crossovers = FALSE) {
  n       <- length(sire_ids)
  col_idx <- seq_len(n_autosome_loci)

  sire_allele <- matrix(0L, nrow = n, ncol = n_autosome_loci)
  dam_allele  <- matrix(0L, nrow = n, ncol = n_autosome_loci)
  sire_lo     <- matrix(NA_character_, nrow = n, ncol = n_autosome_loci)
  dam_lo      <- matrix(NA_character_, nrow = n, ncol = n_autosome_loci)

  # Preallocated to the upper bound (one crossover frame per gamete = 2n).
  xo <- if (store_crossovers) vector("list", 2L * n) else NULL
  xk <- 0L

  for (i in seq_len(n)) {
    o    <- o_index[i]
    sire <- sire_ids[i]
    dam  <- dam_ids[i]

    gs <- make_gamete(parent_haps[[sire]], autosome_info_by_parent[[sire]],
                      base_seed, .gamete_stream_id(o, 1L, 1L), store_crossovers)
    sire_allele[i, ] <- gs$allele
    sire_lo[i, ]     <- parent_lo[[sire]][cbind(gs$homolog, col_idx)]

    gd <- make_gamete(parent_haps[[dam]], autosome_info_by_parent[[dam]],
                      base_seed, .gamete_stream_id(o, 2L, 1L), store_crossovers)
    dam_allele[i, ] <- gd$allele
    dam_lo[i, ]     <- parent_lo[[dam]][cbind(gd$homolog, col_idx)]

    if (store_crossovers) {
      if (length(gs$xover_pos_cM)) {
        xk <- xk + 1L
        xo[[xk]] <- data.frame(local_j = i, parent_origin = 1L,
                               chr_idx = gs$xover_chr_idx, pos_cM = gs$xover_pos_cM)
      }
      if (length(gd$xover_pos_cM)) {
        xk <- xk + 1L
        xo[[xk]] <- data.frame(local_j = i, parent_origin = 2L,
                               chr_idx = gd$xover_chr_idx, pos_cM = gd$xover_pos_cM)
      }
    }
  }

  crossovers <- if (store_crossovers && xk > 0L) do.call(rbind, xo[seq_len(xk)]) else NULL

  list(sire_allele = sire_allele, dam_allele = dam_allele,
       sire_lo = sire_lo, dam_lo = dam_lo, crossovers = crossovers)
}

#' Pass a chromosome copy through to a gamete without recombination
#'
#' Used for non-recombining / single-copy special chromosomes (Y, W, MT, most
#' organelles). Draws from the **current, already-seeded** dqrng stream (the
#' caller seeds the gamete's `"special"` stream once, then processes all of that
#' gamete's special chromosomes in sequence). When `k == 2`, one dqrng uniform
#' chooses the homolog (`u < 0.5`); when `k == 1` there is nothing to choose, so
#' strictly hemizygous inheritance consumes **no** RNG.
#'
#' @param hap_matrix `k x n_chr_loci` integer matrix, `k` in `{1, 2}`.
#' @param lo_matrix `k x n_chr_loci` character `line_origin` matrix, parallel.
#' @return List with `allele` and `line_origin`, each length `n_chr_loci`.
#' @keywords internal
pass_through_gamete <- function(hap_matrix, lo_matrix) {
  k <- nrow(hap_matrix)
  chosen <- if (k == 1L) 1L else if (dqrng::dqrunif(1L) < 0.5) 1L else 2L
  list(
    allele      = hap_matrix[chosen, ],
    line_origin = lo_matrix[chosen, ]
  )
}
