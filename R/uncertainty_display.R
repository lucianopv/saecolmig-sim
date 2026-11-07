#' Display summaries of the uncertainty measures
#'
#' Each function reads a per-municipality summary `a` (the output of
#' `ur_area_summary()`, or several of them stacked) for one `outcome` and
#' returns the aggregates a comparison figure prints. They filter by outcome
#' (and, where stated, to the in-sample municipalities) and by nothing else:
#' a table that holds several covariate configurations must be reduced to one
#' before the call, because the summaries would otherwise pool them.
#' @name uncertainty_display
NULL

# Rows of one outcome (optionally the in-sample municipalities only), reduced to
# `cols`. The selection is built from vectors taken with `$`, so the column
# `outcome` is never confused with the argument of the same name. The result is a
# new table: the caller may add columns to it.
.ur_rows_of <- function(a, outcome, cols, in_sample = FALSE) {
  a <- data.table::as.data.table(a)
  keep <- a$outcome == outcome
  if (in_sample) keep <- keep & a$area_set == "in_sample"
  keep[is.na(keep)] <- FALSE
  a[keep, cols, with = FALSE]
}

#' Ratio of estimated to empirical RMSE, summarised over municipalities
#'
#' Per municipality the ratio is `rep_rmse / sqrt(emp_mse)`: the mean over the
#' replicates of the root of the method's reported MSE, over the root of the
#' realised MSE. Municipalities whose realised MSE is 0 have no ratio and are
#' left out. A ratio below 1 means the method reports less uncertainty than it
#' realises.
#' @param a per-municipality summary with `method`, `area_set`, `outcome`,
#'   `rep_rmse` and `emp_mse`
#' @param outcome the outcome to summarise
#' @return data.table by method and area set: `p25`, `med`, `p75` (quartiles
#'   of the ratio, type 7) and `below` (share of municipalities with a ratio
#'   under 1)
#' @export
ur_ratio_summary <- function(a, outcome) {
  rr <- .ur_rows_of(a, outcome, c("method", "area_set", "rep_rmse", "emp_mse"))
  rr <- rr[emp_mse > 0]
  rr[, ratio := rep_rmse / sqrt(emp_mse)]
  rr[, .(p25 = stats::quantile(ratio, 0.25), med = stats::median(ratio),
         p75 = stats::quantile(ratio, 0.75), below = mean(ratio < 1)),
     by = .(method, area_set)]
}

#' How closely the estimated RMSE tracks the empirical RMSE across municipalities
#'
#' Every municipality of the outcome enters, including those with a realised
#' MSE of 0 (their empirical RMSE is 0).
#' @inheritParams ur_ratio_summary
#' @return data.table by method and area set: `spearman` (rank correlation of
#'   `rep_rmse` with `sqrt(emp_mse)`) and `below` (share of municipalities whose
#'   `rep_rmse` is under `sqrt(emp_mse)`)
#' @export
ur_tracking_summary <- function(a, outcome) {
  a <- .ur_rows_of(a, outcome, c("method", "area_set", "rep_rmse", "emp_mse"))
  a[, `:=`(emp = sqrt(emp_mse), est = rep_rmse)]
  a[, .(spearman = stats::cor(est, emp, method = "spearman"),
        below = mean(est < emp)),
    by = .(method, area_set)]
}

#' Sample-size groups: sparse sizes pooled upward
#'
#' Each distinct value of `n_dwell_med` is its own group while it holds at least
#' `min_n` municipalities. Values that hold fewer are pooled with the next larger
#' values until the group reaches `min_n`, and a short tail joins the last group.
#' Equal values always share a group.
#' @param n_dwell_med median dwellings drawn, one value per municipality
#' @param min_n smallest number of municipalities in a group
#' @return integer group number for each element of `n_dwell_med`, increasing
#'   with the size; an error if the municipalities in total are fewer than
#'   `min_n`
#' @export
ur_size_groups <- function(n_dwell_med, min_n = 25L) {
  sizes <- data.table::data.table(v = n_dwell_med)[, .N, keyby = v]
  grp <- integer(nrow(sizes)); g <- 1L; acc <- 0L
  for (i in seq_len(nrow(sizes))) {
    grp[i] <- g; acc <- acc + sizes$N[i]
    if (acc >= min_n) { g <- g + 1L; acc <- 0L }
  }
  if (acc > 0L && g > 1L) grp[grp == g] <- g - 1L
  if (any(tapply(sizes$N, grp, sum) < min_n)) {
    stop("ur_size_groups: the municipalities are fewer than min_n = ", min_n,
         "; a size group would hold fewer than that")
  }
  grp[match(n_dwell_med, sizes$v)]
}

#' Median empirical and estimated RMSE by sample-size group
#'
#' In-sample municipalities only. A municipality's size is its median dwellings
#' drawn over the replicates in which it was drawn. The groups come from
#' `ur_size_groups()` on the sizes of the first method's municipalities (the
#' first level of `method` for a factor, else the first one in the table); every
#' method is summarised on those groups, so all methods must carry the same
#' sizes, and rows of a method with a size the first method lacks are left out
#' with a warning. Per method, kind and group the size is the median size and the
#' RMSE the median RMSE.
#' @inheritParams ur_ratio_summary
#' @param a per-municipality summary with `method`, `area_set`, `outcome`,
#'   `emp_mse`, `rep_rmse` and `n_dwell_med`
#' @param min_n smallest number of municipalities in a size group
#' @return data.table: `method`, `kind` ("Empirical RMSE" is `sqrt(emp_mse)`,
#'   "Estimated RMSE" is `rep_rmse`), `size_group`, `n_dwell` (median size) and
#'   `rmse` (median RMSE)
#' @export
ur_size_medians <- function(a, outcome, min_n = 25L) {
  sz <- .ur_rows_of(a, outcome, c("method", "emp_mse", "rep_rmse", "n_dwell_med"),
                    in_sample = TRUE)
  sz <- rbind(sz[, .(method, n_dwell_med, kind = "Empirical RMSE", rmse = sqrt(emp_mse))],
              sz[, .(method, n_dwell_med, kind = "Estimated RMSE", rmse = rep_rmse)])
  ms <- if (is.factor(sz$method)) levels(droplevels(sz$method)) else unique(sz$method)
  first_method <- ms[1L]
  e <- sz[kind == "Empirical RMSE" & method == first_method]
  sizes <- unique(data.table::data.table(n_dwell_med = e$n_dwell_med,
                                         size_group = ur_size_groups(e$n_dwell_med, min_n)))
  m <- merge(sz, sizes, by = "n_dwell_med")
  if (nrow(m) < nrow(sz)) {
    warning("ur_size_medians: ", nrow(sz) - nrow(m), " rows have a size not carried by ",
            "the first method's municipalities and are left out", call. = FALSE)
  }
  m[, .(n_dwell = stats::median(n_dwell_med), rmse = stats::median(rmse)),
    by = .(method, kind, size_group)]
}

#' Coverage of the 95\% interval, summarised over in-sample municipalities
#'
#' `overall` is the coverage of every interval: the municipalities' coverages
#' weighted by their number of replicates. The two shares count municipalities.
#' @param a per-municipality summary with `method`, `area_set`, `outcome`,
#'   `n_reps` and `coverage`
#' @inheritParams ur_ratio_summary
#' @return data.table by method: `overall`, `share_ge_095` (share of
#'   municipalities with coverage of at least 0.95) and `share_lt_080` (share
#'   with coverage under 0.80)
#' @export
ur_coverage_summary <- function(a, outcome) {
  cvm <- .ur_rows_of(a, outcome, c("method", "n_reps", "coverage"), in_sample = TRUE)
  cvm[, .(overall = stats::weighted.mean(coverage, n_reps),
          share_ge_095 = mean(coverage >= 0.95),
          share_lt_080 = mean(coverage < 0.8)),
      by = method]
}
