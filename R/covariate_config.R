#' Covariate selections, XGB hyperparameters and the sub-area universe
#'
#' The area-level configuration `geo_bld` is the census covariates of the
#' municipal matrix, a geo block and a buildings block; the sub-area
#' configuration `geo_blk_bld` is the census shares of each census block or
#' rural section, an urban flag, the same geo layers extracted at sub-area
#' resolution and the same buildings columns. Which geo and buildings columns
#' enter is fixed by two column lists shipped in `inst/extdata`, the selections
#' the published fits used. They came from a VIF prune (threshold 10) against
#' the other covariates of the parent matrix; a prune against a different set of
#' covariates can keep different columns, so the builders read the lists rather
#' than re-running the prune. `bld_prune_vif()` runs the prune itself.
#'
#' @name covariate_config
NULL

# Geo columns whose name carries the census year of the layer. The sub-area
# matrix names them after the layer alone; the values are extractions of the
# same layers.
.geo_year_labels <- c(
  pop_density_2018       = "pop_density",
  pop_change_2015_2018   = "pop_change",
  precip_annual_2018     = "precip_annual",
  precip_annual_2018_pw  = "precip_annual_pw",
  precip_anomaly_2018    = "precip_anomaly",
  precip_anomaly_2018_pw = "precip_anomaly_pw")

#' Sub-area names of geo columns
#'
#' Maps the column names `build_geo_aux()` returns to the names the sub-area
#' matrix uses: the census-year label is dropped (`pop_density_2018` ->
#' `pop_density`, `precip_annual_2018_pw` -> `precip_annual_pw`, ...). Names
#' without a year label are returned unchanged.
#'
#' @param x character vector of column names
#' @return character vector of the same length
#' @export
geo_subarea_names <- function(x) {
  x <- as.character(x)
  hit <- x %in% names(.geo_year_labels)
  x[hit] <- unname(.geo_year_labels[x[hit]])
  x
}

#' Geo and buildings columns of the published configurations
#'
#' @param which `"geo"` for the geo block, `"bld"` for the buildings block
#' @param level `"area"` (default) gives the municipal matrix's names;
#'   `"subarea"` the sub-area matrix's (see [geo_subarea_names()]; the
#'   buildings names are the same at both levels)
#' @return character vector of column names, in matrix order
#' @export
kept_cols <- function(which = c("geo", "bld"), level = c("area", "subarea")) {
  which <- match.arg(which)
  level <- match.arg(level)
  f <- system.file("extdata", sprintf("kept_%s_cols.txt", which),
                   package = "saecol2026", mustWork = TRUE)
  x <- trimws(readLines(f, warn = FALSE))
  x <- x[nzchar(x)]
  if (identical(level, "subarea")) x <- geo_subarea_names(x)
  x
}

# Hyperparameters held as integers by xgboost's interface; the rest are numeric.
.xgb_int_params <- c("nrounds", "max_depth")

#' Tuned hyperparameters of the sub-area XGB method
#'
#' One set per outcome, selected by the procedure of `data-raw/tune_xgb_sub.R`:
#' a five-fold cross-validation grid over learning rate, depth and minimum
#' child weight on each of the first five replicates, each parameter then set
#' to its most frequent winner. The shipped table was selected on the sub-area
#' matrix of the published fits; tuning on other covariates can select other
#' values. The values are passed to [fit_xgb_sub()] as `params`.
#'
#' @param config sub-area configuration whose shipped table is read
#' @param file path of a table in the same layout, read instead of the shipped
#'   one (for example one written by `data-raw/tune_xgb_sub.R` on other
#'   covariates)
#' @return named list, outcome -> list of hyperparameters (`nrounds`, `eta`,
#'   `max_depth`, ...), in `sae_outcomes` order, with the path of the table it
#'   was read from in the attribute `source`
#' @export
xgb_sub_params <- function(config = "geo_blk_bld", file = NULL) {
  if (is.null(file)) {
    file <- system.file("extdata", sprintf("xgb_sub_params_%s.tsv", config),
                        package = "saecol2026")
    if (!nzchar(file)) {
      stop("xgb_sub_params: no tuned XGB hyperparameters ship for configuration '",
           config, "'; tune with data-raw/tune_xgb_sub.R and pass its table as `file`")
    }
  }
  hdr <- strsplit(readLines(file, n = 1L, warn = FALSE), "\t", fixed = TRUE)[[1]]
  if (!"outcome" %in% hdr) stop("xgb_sub_params: ", file, " has no `outcome` column")
  need <- c("nrounds", "eta", "max_depth")
  if (!all(need %in% hdr)) {
    stop("xgb_sub_params: ", file, " lacks ", paste(setdiff(need, hdr), collapse = ", "))
  }
  d <- utils::read.delim(file, colClasses = ifelse(hdr == "outcome", "character", "numeric"),
                         stringsAsFactors = FALSE)
  unknown <- setdiff(d$outcome, sae_outcomes$outcome)
  missing <- setdiff(sae_outcomes$outcome, d$outcome)
  if (length(unknown) || length(missing) || anyDuplicated(d$outcome)) {
    stop("xgb_sub_params: ", file, " must hold one row per outcome;",
         if (length(missing)) paste0(" missing: ", paste(missing, collapse = ", ")),
         if (length(unknown)) paste0(" unknown: ", paste(unknown, collapse = ", ")),
         if (anyDuplicated(d$outcome)) " duplicated outcomes")
  }
  pn <- setdiff(hdr, "outcome")
  na_cols <- pn[vapply(pn, function(k) anyNA(d[[k]]), logical(1))]
  if (length(na_cols)) {
    stop("xgb_sub_params: ", file, " has missing values in ",
         paste(na_cols, collapse = ", "))
  }
  out <- lapply(seq_len(nrow(d)), function(i) {
    x <- as.list(d[i, pn, drop = FALSE])
    for (k in intersect(.xgb_int_params, pn)) x[[k]] <- as.integer(x[[k]])
    x
  })
  names(out) <- d$outcome
  out <- out[sae_outcomes$outcome]
  attr(out, "source") <- file
  out
}

#' Write sub-area XGB hyperparameters as a table
#'
#' The layout [xgb_sub_params()] reads: one tab-separated row per outcome, one
#' column per hyperparameter.
#'
#' @param params named list, outcome -> list of hyperparameters, every outcome
#'   holding the same names
#' @param file output path
#' @return `file`, invisibly
#' @export
write_xgb_sub_params <- function(params, file) {
  if (!is.list(params) || !length(params) || is.null(names(params)) ||
      any(!nzchar(names(params)))) {
    stop("write_xgb_sub_params: params must be a non-empty list named by outcome")
  }
  pn <- names(params[[1]])
  same <- vapply(params, function(x) identical(names(x), pn), logical(1))
  if (!all(same)) {
    stop("write_xgb_sub_params: every outcome must hold the same hyperparameters")
  }
  d <- data.frame(outcome = names(params), stringsAsFactors = FALSE)
  for (k in pn) d[[k]] <- vapply(params, function(x) as.numeric(x[[k]]), numeric(1))
  utils::write.table(d, file, sep = "\t", quote = FALSE, row.names = FALSE)
  invisible(file)
}

#' Covariate columns of a sub-area matrix
#'
#' Every column except the structural ones: the key `code`, `municipality`, the
#' class `stratum`, the census counts `n_households_census` and
#' `n_persons_census`, and the WorldPop roll-up weight `pop_worldpop` with the
#' polygon area `area_km2` it is derived from.
#'
#' @param pop_sub sub-area covariate matrix (data.frame or data.table)
#' @return character vector of covariate column names, in matrix order
#' @export
subarea_covariate_names <- function(pop_sub) {
  setdiff(names(pop_sub), .subarea_id_cols)
}

#' Sub-area prediction universe
#'
#' A sub-area enters the prediction universe when every covariate of the
#' sub-area matrix is observed: one that is not cannot be predicted, and the
#' municipal estimate is a weighted mean over the sub-areas the method
#' predicts. Structural columns (keys, census counts, area and WorldPop weight)
#' may be missing without excluding a row.
#'
#' @param pop_sub sub-area covariate matrix (data.frame or data.table)
#' @return `pop_sub` restricted to the rows complete on its covariate columns,
#'   of the same class
#' @export
subarea_universe <- function(pop_sub) {
  x_cols <- subarea_covariate_names(pop_sub)
  if (!length(x_cols)) {
    stop("subarea_universe: the sub-area matrix has no covariate columns besides ",
         "its structural ones (", paste(.subarea_id_cols, collapse = ", "), ")")
  }
  keep <- stats::complete.cases(as.data.frame(pop_sub)[, x_cols, drop = FALSE])
  if (data.table::is.data.table(pop_sub)) pop_sub[keep] else pop_sub[keep, , drop = FALSE]
}
