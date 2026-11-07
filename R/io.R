#' Resolve a path under the data directory
#'
#' All inputs and outputs of the simulation (census frame, covariate tables,
#' samples, results) live under one directory named by the `SAECOL_DATA_DIR`
#' environment variable. Set it in `.Renviron` (see `.Renviron.example`) before
#' using this package.
#' @section The covariate tables live in "aux_tables":
#' `AUX` is a reserved device name on Windows: no file or directory called
#' `aux` can be created there, on any file system, and a zip file holding an
#' entry of that name cannot be extracted. The covariate directory is therefore
#' called `aux_tables` on every platform, and a path component spelled exactly
#' `"aux"` is mapped onto `"aux_tables"`, so a path built through this function
#' never names the reserved directory.
#'
#' @param ... path components appended to the data root
#' @param create logical; create the resulting directory if missing
#' @return the resolved path (character)
#' @export
output_path <- function(..., create = FALSE) {
  base <- Sys.getenv("SAECOL_DATA_DIR", unset = NA)
  if (is.na(base) || !nzchar(base))
    stop("Set SAECOL_DATA_DIR to the data folder (see .Renviron.example).")
  parts <- list(...)
  # Applied on every platform, so a path lands in the same directory wherever
  # it is built (see the roxygen section above).
  parts <- lapply(parts, function(x) if (identical(x, "aux")) "aux_tables" else x)
  p <- do.call(file.path, c(list(base), parts))
  if (create) dir.create(p, recursive = TRUE, showWarnings = FALSE)
  p
}

#' @export
census_path <- function() output_path("Census2018")

#' @export
ecv_path <- function() output_path("ECV2024")

#' Write a parquet file, ensuring the parent directory exists
#' @export
write_parquet_safe <- function(df, path) {
  dir <- dirname(path)
  if (!dir.exists(dir)) dir.create(dir, recursive = TRUE)
  arrow::write_parquet(df, path)
  invisible(path)
}

#' @export
read_parquet_safe <- function(path) {
  if (!file.exists(path)) stop("Parquet file not found: ", path)
  arrow::read_parquet(path)
}
