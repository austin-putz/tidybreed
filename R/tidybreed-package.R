#' @keywords internal
"_PACKAGE"

## usethis namespace: start
#' @importFrom tibble tibble
#' @importFrom dplyr filter mutate select tbl
#' @importFrom rlang .data
#' @importFrom DBI dbConnect dbDisconnect dbIsValid dbListTables
#' @importFrom DBI dbGetQuery dbExecute dbWriteTable dbListFields
#' @import duckdb
#' @useDynLib tidybreed, .registration = TRUE
#' @importFrom Rcpp sourceCpp
## usethis namespace: end
NULL

.onLoad <- function(libname, pkgname) {
  op      <- options()
  defaults <- list(
    tidybreed.base_dir        = getwd(),
    tidybreed.output          = "tidybreed_output",
    tidybreed.scenario        = NULL,
    tidybreed.tools           = NULL,
    tidybreed.db_name         = "sim.duckdb",
    tidybreed.pop_name        = "sim",
    tidybreed.quiet           = FALSE,
    tidybreed.replicate       = 1L,
    tidybreed.archive_path    = NULL,
    tidybreed.db_name_archive = NULL
  )
  toset <- !(names(defaults) %in% names(op))
  if (any(toset)) options(defaults[toset])
  invisible()
}

.onAttach <- function(libname, pkgname) {
  if (isTRUE(getOption("tidybreed.quiet", FALSE))) return(invisible())

  ver      <- as.character(utils::packageVersion("tidybreed"))
  duck_ver <- as.character(utils::packageVersion("duckdb"))

  fmt_opt <- function(name) {
    val <- getOption(name)
    if (is.null(val)) return("NULL")
    if (is.character(val)) return(paste0("\"", paste(val, collapse = "\", \""), "\""))
    paste(val, collapse = ", ")
  }

  opt_names <- c("tidybreed.pop_name", "tidybreed.base_dir", "tidybreed.output",
                 "tidybreed.scenario", "tidybreed.tools", "tidybreed.db_name",
                 "tidybreed.replicate", "tidybreed.archive_path",
                 "tidybreed.db_name_archive", "tidybreed.quiet")
  opt_width <- max(nchar(opt_names))
  opt_lines <- vapply(opt_names, function(name) {
    paste0("    ", formatC(name, width = -opt_width), " = ", fmt_opt(name))
  }, character(1))

  msg <- paste(c(
    cli::rule(
      left  = paste0(cli::style_bold(cli::col_cyan("tidybreed")), " ", ver),
      width = 72
    ),
    paste0("  ", cli::col_grey("Tidy Breeding Program Simulation")),
    paste0("  ", cli::col_grey("License: MIT · No warranty; use at your own risk")),
    "",
    paste0("  ", cli::style_bold("Citation:"), "  citation(\"tidybreed\")"),
    paste0("  ", cli::style_bold("Help:    "), "  ?tidybreed  ·  browseVignettes(\"tidybreed\")"),
    "",
    paste0("  ", cli::style_bold("Backend: "), "  DuckDB ", duck_ver),
    paste0("  ", cli::style_bold("Storage: "), "  All data lives on disk in a .duckdb file, not in R memory"),
    "",
    paste0("  ", cli::style_bold("Options"), " ", cli::col_grey("(current values, set with options()):")),
    opt_lines,
    cli::rule(width = 72)
  ), collapse = "\n")

  packageStartupMessage(msg)
}
