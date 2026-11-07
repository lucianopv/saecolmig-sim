#' Parse the ECV 2024 design YAML into a design object
#'
#' The YAML's `strata:` block can be either:
#'   - a list of entries with `stratum` and `n_psu` (used in tests with synthetic data), or
#'   - a single mapping with `source_csv`, `stratum_col`, `n_psu_col` that
#'     points to a CSV of per-stratum allocations (the shipped YAML points to
#'     `ecv2024_strata_realised.csv`).
#'
#' @param path optional path to the YAML; defaults to the bundled
#'   `inst/extdata/ecv2024_design.yml`.
#' @return a list with class "ecv_design_spec".
#' @export
parse_ecv_design <- function(path = NULL) {
  if (is.null(path)) {
    path <- system.file("extdata", "ecv2024_design.yml", package = "saecol2026")
    if (!nzchar(path)) {
      # devtools::load_all() context: fall back to the in-repo path
      path <- file.path("inst", "extdata", "ecv2024_design.yml")
    }
  }
  if (!file.exists(path)) stop("Design YAML not found: ", path)
  raw <- yaml::read_yaml(path)
  raw$strata <- coerce_strata(raw$strata, yaml_path = path)
  structure(raw, class = c("ecv_design_spec", "list"))
}

# Internal: turn the YAML's strata block into a tibble(stratum, n_psu).
coerce_strata <- function(strata_block, yaml_path) {
  if (is.list(strata_block) && !is.null(strata_block$source_csv)) {
    # Resolve CSV path relative to the YAML's directory.
    csv_path <- strata_block$source_csv
    if (!file.exists(csv_path)) {
      candidate <- file.path(dirname(yaml_path), basename(csv_path))
      if (file.exists(candidate)) csv_path <- candidate
    }
    if (!file.exists(csv_path)) {
      stop("Strata source CSV not found: ", strata_block$source_csv)
    }
    df <- utils::read.csv(csv_path, stringsAsFactors = FALSE)
    sc <- .null_or(strata_block$stratum_col, "stratum")
    nc <- .null_or(strata_block$n_psu_col,   "n_psu")
    out <- data.frame(stratum = as.character(df[[sc]]),
                      n_psu   = as.integer(df[[nc]]),
                      stringsAsFactors = FALSE)
    # Optional n_seg column lets draw_replicate use the stratum's realised
    # segment count per PSU instead of a global stage2$n_per_psu fallback.
    if ("n_seg" %in% names(df)) {
      out$n_seg <- as.integer(df$n_seg)
      out$n_seg_per_psu <- pmax(1L, as.integer(round(out$n_seg / out$n_psu)))
    }
    return(out)
  }
  # Inline list of {stratum, n_psu}
  out <- do.call(rbind, lapply(strata_block, as.data.frame))
  out$stratum <- as.character(out$stratum)
  out$n_psu   <- as.integer(out$n_psu)
  rownames(out) <- NULL
  out
}

#' @export
print.ecv_design_spec <- function(x, ...) {
  cat("<ecv_design_spec>\n")
  cat("  stages:", x$stage1$unit, "->", x$stage2$unit, "->", x$stage3$unit, "\n")
  cat("  total PSUs across strata:", sum(x$strata$n_psu), "\n")
  cat("  strata:", nrow(x$strata), "\n")
  invisible(x)
}
