#' Haplotype long-write + memory-batching helpers
#'
#' Internal helpers shared by the three haplotype write paths
#' (`add_offspring()`, `add_founders()`, and `define_founder_haplotypes()`'s
#' pool write). They convert the write path to direct **long** inserts and bound
#' peak memory with RAM-aware batching. No values change — only the write shape.
#'
#' @keywords internal
#' @name haplotype_write_helpers
NULL


#' Detect available system memory in bytes (cross-OS, best-effort)
#'
#' Returns available (not total) physical RAM in bytes on Linux, macOS
#' (Intel + Apple Silicon), and Windows, or `NA_real_` when detection fails or
#' the OS is unsupported. **Never errors** — callers fall back to a conservative
#' fixed budget when this returns `NA`.
#'
#' - Linux: `/proc/meminfo` `MemAvailable` (falls back to `MemFree`), in kB.
#' - macOS: `vm_stat` free + inactive + speculative pages x `hw.pagesize`.
#' - Windows: `wmic OS get FreePhysicalMemory` (kB), then a PowerShell fallback.
#'
#' @return Numeric scalar bytes, or `NA_real_`.
#' @keywords internal
detect_available_memory <- function() {
  sysname <- Sys.info()[["sysname"]]

  bytes <- tryCatch(
    suppressWarnings({
      if (identical(sysname, "Linux")) {
        mi   <- readLines("/proc/meminfo", warn = FALSE)
        line <- grep("^MemAvailable:", mi, value = TRUE)
        if (length(line) == 0L) line <- grep("^MemFree:", mi, value = TRUE)
        if (length(line) == 0L) return(NA_real_)
        as.numeric(gsub("[^0-9]", "", line[1])) * 1024

      } else if (identical(sysname, "Darwin")) {
        page <- as.numeric(system2("sysctl", c("-n", "hw.pagesize"),
                                   stdout = TRUE, stderr = FALSE))
        vs   <- system2("vm_stat", stdout = TRUE, stderr = FALSE)
        if (length(vs) == 0L || is.na(page)) return(NA_real_)
        getpages <- function(key) {
          l <- grep(key, vs, value = TRUE, fixed = TRUE)
          if (length(l) == 0L) return(0)
          as.numeric(gsub("[^0-9]", "", sub(".*:", "", l[1])))
        }
        (getpages("Pages free") + getpages("Pages inactive") +
           getpages("Pages speculative")) * page

      } else if (identical(sysname, "Windows")) {
        out  <- system2("wmic", c("OS", "get", "FreePhysicalMemory", "/Value"),
                        stdout = TRUE, stderr = FALSE)
        line <- grep("FreePhysicalMemory=", out, value = TRUE)
        if (length(line) > 0L) {
          kb <- as.numeric(gsub("[^0-9]", "", line[1]))
          if (!is.na(kb) && kb > 0) return(kb * 1024)
        }
        # PowerShell fallback (wmic is deprecated on newer Windows).
        ps <- system2("powershell",
                      c("-NoProfile", "-Command",
                        "(Get-CimInstance Win32_OperatingSystem).FreePhysicalMemory"),
                      stdout = TRUE, stderr = FALSE)
        kb <- as.numeric(gsub("[^0-9]", "", ps[1]))
        if (is.na(kb)) return(NA_real_)
        kb * 1024

      } else {
        NA_real_
      }
    }),
    error = function(e) NA_real_
  )

  if (length(bytes) != 1L || is.na(bytes) || !is.finite(bytes) || bytes <= 0)
    return(NA_real_)
  bytes
}


#' Parse a memory size (bytes number or a string like "512MB", "2GB")
#'
#' @param x Numeric bytes, or a string with an optional k/m/g/t suffix
#'   (case-insensitive; a trailing "B" is optional). Powers of 1024.
#' @return Numeric bytes.
#' @keywords internal
parse_mem_size <- function(x) {
  if (is.numeric(x)) {
    if (length(x) != 1L || is.na(x) || x <= 0)
      stop("max_batch_mem must be a positive size.", call. = FALSE)
    return(as.numeric(x))
  }
  s <- trimws(as.character(x))
  m <- regmatches(s, regexec("^([0-9.]+)\\s*([kKmMgGtT]?)[bB]?$", s))[[1]]
  if (length(m) != 3L)
    stop("Cannot parse memory size '", x, "'. Use e.g. 512e6 or \"512MB\".",
         call. = FALSE)
  num    <- as.numeric(m[2])
  suffix <- toupper(m[3])
  mult   <- if (suffix == "") 1 else
    switch(suffix, "K" = 1024, "M" = 1024^2, "G" = 1024^3, "T" = 1024^4)
  num * mult
}


#' Resolve the offspring/founders-per-batch size `B`
#'
#' Precedence: explicit `batch_size` > explicit `max_batch_mem` > RAM-aware
#' auto-pick (target `target_frac` of detected available memory) > conservative
#' fixed fallback (`fallback_bytes`) when detection fails. Always clamped to
#' `[1, n_total]`. Any `B` yields byte-identical output — this only bounds peak
#' memory, never changes values.
#'
#' @param n_total Integer. Total rows to process (offspring or founders).
#' @param n_loci Integer. Loci count (drives per-row memory).
#' @param batch_size Optional explicit batch size (wins when set).
#' @param max_batch_mem Optional explicit per-batch memory budget (bytes or
#'   `"512MB"`-style string).
#' @param bytes_per_offspring_row Heuristic peak bytes per (offspring x locus)
#'   during the in-R long build + register. Calibrated (~200) from Phase-0/Stage-1
#'   measurements: each offspring emits `2 x n_loci` long rows (two haplotypes),
#'   each ~90 bytes at peak including the `rbind`/register transient copy — so the
#'   per-(offspring x locus) constant is ~2 x 90 ~ 190, rounded up for headroom.
#' @param target_frac Fraction of available memory to target for one batch.
#' @param fallback_bytes Conservative per-batch budget when detection fails.
#' @return Integer `B` in `[1, n_total]` (or `0L` when `n_total == 0`).
#' @keywords internal
resolve_batch_size <- function(n_total, n_loci,
                               batch_size = NULL, max_batch_mem = NULL,
                               bytes_per_offspring_row = 200,
                               target_frac = 0.25,
                               fallback_bytes = 512 * 1024^2) {
  n_total <- as.integer(n_total)
  if (is.na(n_total) || n_total <= 0L) return(0L)

  if (!is.null(batch_size)) {
    b <- suppressWarnings(as.integer(batch_size))
    if (length(batch_size) != 1L || is.na(b) || b < 1L)
      stop("batch_size must be a single positive integer.", call. = FALSE)
    return(min(b, n_total))
  }

  target_bytes <- if (!is.null(max_batch_mem)) {
    parse_mem_size(max_batch_mem)
  } else {
    avail <- detect_available_memory()
    if (is.na(avail)) fallback_bytes else avail * target_frac
  }

  per_offspring <- as.numeric(n_loci) * bytes_per_offspring_row
  b <- floor(target_bytes / per_offspring)
  if (!is.finite(b) || b < 1) b <- 1
  as.integer(min(max(1, b), n_total))
}


#' Insert a long haplotype frame into `ind_haplotype` (no UNPIVOT)
#'
#' The single shared write core for `add_offspring()` and `add_founders()`.
#' Registers `long_frame` and does one `INSERT ... SELECT` that joins
#' `genome_meta` on `locus_id` to attach `locus_name`. `strand` is always `1L`
#' in this ploidy-2 version. RNG-neutral (register + INSERT, no `dbWriteTable`).
#'
#' Write mechanism (Stage-1 Appender vs register+INSERT decision): we use
#' `duckdb_register()` + `INSERT ... SELECT`, not the DuckDB Appender API. It is
#' RNG-neutral (unlike `dbWriteTable`, it advances no R RNG — the founder/offspring
#' RNG discipline relies on this), zero-copy over the R frame, and already the
#' fast path the pre-Stage-1 UNPIVOT write used. The Appender API would save only
#' the per-batch view registration (marginal, and batching already bounds that),
#' at the cost of a C++ boundary; not worth it here. Revisit if profiling ever
#' shows registration overhead dominating.
#'
#' @param conn DuckDB connection.
#' @param long_frame data.frame with columns `id_ind`, `parent_origin`,
#'   `strand`, `line_origin`, `locus_id`, `allele`.
#' @return Invisibly `NULL`.
#' @keywords internal
.write_long_haplotypes <- function(conn, long_frame) {
  if (is.null(long_frame) || nrow(long_frame) == 0L) return(invisible(NULL))
  duckdb::duckdb_register(conn, "__tmp_long_hap", as.data.frame(long_frame))
  on.exit(duckdb::duckdb_unregister(conn, "__tmp_long_hap"), add = TRUE)
  DBI::dbExecute(conn, paste0(
    "INSERT INTO ind_haplotype ",
    "(id_ind, parent_origin, strand, line_origin, locus_id, locus_name, allele) ",
    "SELECT f.id_ind, f.parent_origin, f.strand, f.line_origin, ",
    "f.locus_id, gm.locus_name, f.allele ",
    "FROM __tmp_long_hap f JOIN genome_meta gm ON f.locus_id = gm.locus_id"
  ))
  invisible(NULL)
}
