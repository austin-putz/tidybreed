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

#' Reconstruct the per-chromosome `chr_info` list from flat kernel arrays
#'
#' The C++-ready seam ([make_gametes_batch_r()]) receives the autosome map as
#' flat arrays (`chr_start`/`chr_end`/`chr_pos_cM`/`chr_len_cM`) rather than the
#' nested `chr_info` list [make_gamete()] consumes. This rebuilds that list so
#' the unchanged per-gamete core can be reused verbatim. `chr_start`/`chr_end`
#' are 1-based inclusive locus-index ranges into `1..n_autosome_loci`.
#'
#' @param chr_start,chr_end Integer per-chromosome locus-index ranges.
#' @param chr_pos_cM Numeric per-locus genetic positions (cM), length
#'   `n_autosome_loci`.
#' @param chr_len_cM Numeric per-chromosome genetic length (cM).
#' @return A `chr_info`-shaped list (one element per chromosome).
#' @keywords internal
.chr_info_from_arrays <- function(chr_start, chr_end, chr_pos_cM, chr_len_cM) {
  lapply(seq_along(chr_start), function(c_i) {
    idx <- chr_start[c_i]:chr_end[c_i]
    list(locus_idx = idx, pos_cM = chr_pos_cM[idx], chr_len = chr_len_cM[c_i])
  })
}

#' Generate autosome gametes for one resolved-map group (the swappable seam)
#'
#' The single seam the C++ kernel (Stage 3) drops into: a dependency-free numeric
#' function over **packed integer arrays** — no `tidybreed_pop`, no DBI, no
#' `data.frame`, no character `line_origin`, no string-keyed lists — so the R
#' reference and the C++ kernel take byte-identical inputs and their outputs can
#' be compared directly for parity. This is the R reference / parity oracle;
#' [make_gametes_batch()] dispatches to it (or to the compiled kernel).
#'
#' **Gamete-flat / one-map-per-call (group-by-map seam).** Each gamete may use a
#' different genetic map (sex-specific, achiasmy, line-specific), so the caller
#' groups gametes by resolved-map key and calls this once per group with that
#' group's `chr_*` arrays. The unit is a **gamete**, identified by
#' `(gamete_o, gamete_origin)`; each draws from its own dqrng stream keyed on its
#' **global** offspring index `gamete_o` (autosome kind = 1), so output is
#' independent of batch size, iteration order, and how gametes are partitioned
#' into map-groups or calls.
#'
#' @param parent_allele Integer matrix `(2 * n_parents) x n_autosome_loci`,
#'   homolog-major: rows `2p-1` and `2p` are parent `p`'s homolog 1 and 2.
#' @param parent_lo_code Integer matrix, same shape: `line_origin` dictionary
#'   codes (the caller decodes back to `VARCHAR`).
#' @param gamete_parent_idx Integer length `G`: 1-based packed parent row that
#'   produces each gamete.
#' @param gamete_o Integer length `G`: global offspring index `o` per gamete
#'   (original `matings` row); drives the per-gamete stream.
#' @param gamete_origin Integer length `G`: `parent_origin` (1 = parent_1/sire,
#'   2 = parent_2/dam) per gamete.
#' @param chr_start,chr_end Integer per-chromosome locus-index ranges (1-based
#'   inclusive, contiguous, into `1..n_autosome_loci`) for THIS group's map.
#' @param chr_pos_cM Numeric per-locus genetic positions (cM), length
#'   `n_autosome_loci`, for this group's map.
#' @param chr_len_cM Numeric per-chromosome genetic length (cM) for this group.
#' @param base_seed Integer base seed for the `add_offspring()` call.
#' @param store_crossovers Logical; emit the ragged crossover buffer when `TRUE`.
#' @return List of long parallel vectors over the `G` gametes in the given order
#'   (gamete -> locus ascending): `parent_origin`, `locus_idx`
#'   (1..n_autosome_loci), `allele`, `line_origin_code`; plus the ragged
#'   crossover buffer `xover_gamete_o` (global `o`), `xover_parent_origin`,
#'   `xover_chr` (index into the autosome chromosome order), `xover_pos_cM`
#'   (length 0 unless `store_crossovers`).
#' @keywords internal
make_gametes_batch_r <- function(parent_allele, parent_lo_code,
                                 gamete_parent_idx, gamete_o, gamete_origin,
                                 chr_start, chr_end, chr_pos_cM, chr_len_cM,
                                 base_seed, store_crossovers = FALSE) {
  n_chr <- length(chr_start)
  # Contiguity contract: chromosomes tile 1..n_autosome_loci with no gaps, in
  # ascending order (holds by construction today; guards future out-of-order
  # locus_id additions — see plans/optimize_add_offspring.md Stage 0 caveat).
  if (n_chr > 1L &&
      !identical(as.integer(chr_start[-1L]), as.integer(chr_end[-n_chr]) + 1L)) {
    stop("make_gametes_batch_r(): chromosomes must be contiguous locus-index blocks.",
         call. = FALSE)
  }

  G       <- length(gamete_parent_idx)
  n_loci  <- ncol(parent_allele)
  col_idx <- seq_len(n_loci)
  chr_info <- .chr_info_from_arrays(chr_start, chr_end, chr_pos_cM, chr_len_cM)

  total         <- G * n_loci
  parent_origin <- integer(total)
  locus_idx     <- integer(total)
  allele        <- integer(total)
  lo_code       <- integer(total)

  # Ragged crossover accumulator, preallocated to the upper bound (G gametes).
  xo <- if (store_crossovers) vector("list", G) else NULL
  xk <- 0L

  pos <- 0L
  for (i in seq_len(G)) {
    o    <- gamete_o[i]
    r    <- gamete_origin[i]
    pidx <- gamete_parent_idx[i]
    rows <- c(2L * pidx - 1L, 2L * pidx)
    g    <- make_gamete(parent_allele[rows, , drop = FALSE], chr_info,
                        base_seed, .gamete_stream_id(o, r, 1L), store_crossovers)

    rng                <- (pos + 1L):(pos + n_loci)
    parent_origin[rng] <- r
    locus_idx[rng]     <- col_idx
    allele[rng]        <- g$allele
    lo_code[rng]       <- parent_lo_code[rows, , drop = FALSE][cbind(g$homolog, col_idx)]
    pos                <- pos + n_loci

    if (store_crossovers && length(g$xover_pos_cM)) {
      xk <- xk + 1L
      xo[[xk]] <- data.frame(o = o, parent_origin = r,
                             chr = g$xover_chr_idx, pos_cM = g$xover_pos_cM)
    }
  }

  x <- if (store_crossovers && xk > 0L) do.call(rbind, xo[seq_len(xk)]) else
    data.frame(o = integer(0), parent_origin = integer(0),
               chr = integer(0), pos_cM = numeric(0))

  list(parent_origin = parent_origin, locus_idx = locus_idx,
       allele = allele, line_origin_code = lo_code,
       xover_gamete_o = x$o, xover_parent_origin = x$parent_origin,
       xover_chr = x$chr, xover_pos_cM = x$pos_cM)
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
