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
    tidybreed.base_dir  = getwd(),
    tidybreed.output    = "tidybreed_output",
    tidybreed.scenario  = NULL,
    tidybreed.tools     = NULL,
    tidybreed.db_name   = "sim.duckdb",
    tidybreed.pop_name  = "sim"
  )
  toset <- !(names(defaults) %in% names(op))
  if (any(toset)) options(defaults[toset])
  invisible()
}
