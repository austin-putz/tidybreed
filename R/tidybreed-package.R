#' @keywords internal
"_PACKAGE"

## usethis namespace: start
#' @importFrom tibble tibble
#' @importFrom dplyr filter mutate select tbl
#' @importFrom rlang .data
#' @importFrom DBI dbConnect dbDisconnect dbIsValid dbListTables
#' @importFrom DBI dbGetQuery dbExecute dbWriteTable dbListFields
#' @import duckdb
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
  base_dir <- getOption("tidybreed.base_dir", getwd())

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
    paste0("  ", cli::style_bold("Options"), " ", cli::col_grey("(set with options()):")),
    paste0("    tidybreed.base_dir  = \"", base_dir, "\""),
    paste0("    tidybreed.db_name   = \"sim.duckdb\""),
    paste0("    tidybreed.pop_name  = \"sim\""),
    paste0("    tidybreed.output    = \"tidybreed_output\""),
    paste0("    tidybreed.quiet     = FALSE"),
    cli::rule(width = 72)
  ), collapse = "\n")

  packageStartupMessage(msg)
}
