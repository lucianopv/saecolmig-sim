# R/variance_reduction.R
#
# Efficiency of a model against the direct estimator: the model's reported MSE
# over the direct sampling variance, and the realised RMSE of the model over
# that of the direct estimator, area by area.

#' Build a "direct" pseudo-method panel for evaluate_estimates()
#'
#' Turns a panel of direct estimates into the layout of model estimates with
#' method = "direct", so the direct estimator can be evaluated by
#' evaluate_estimates() alongside the SAE methods (same bias/rmse/
#' coverage_95/calibration_ratio columns, directly comparable).
#'
#' @param direct_panel tibble/data.frame with columns area_id, outcome,
#'   replicate, theta_direct, var_direct_smoothed (the direct-estimate files of
#'   several replicates row-bound, with a replicate column added).
#' @return tibble with columns method, replicate, area_id, outcome,
#'   theta_hat, mse, ci_lo, ci_hi (method is always "direct").
#' @export
build_direct_pseudo_method <- function(direct_panel) {
  stopifnot(all(c("area_id", "outcome", "replicate", "theta_direct",
                  "var_direct_smoothed") %in% colnames(direct_panel)))
  se <- sqrt(direct_panel$var_direct_smoothed)
  tibble::tibble(
    method    = "direct",
    replicate = direct_panel$replicate,
    area_id   = direct_panel$area_id,
    outcome   = direct_panel$outcome,
    theta_hat = direct_panel$theta_direct,
    mse       = direct_panel$var_direct_smoothed,
    ci_lo     = direct_panel$theta_direct - 1.96 * se,
    ci_hi     = direct_panel$theta_direct + 1.96 * se
  )
}

#' Paired per-area-replicate ratio of model mse to direct variance
#'
#' @param results_panel model estimates (tibble/data.frame): method,
#'   replicate, area_id, outcome, mse (pass in-sample rows only — direct
#'   variance only exists for sampled areas).
#' @param direct_panel direct estimates row-bound across replicates
#'   (tibble/data.frame): area_id, outcome, replicate, var_direct,
#'   var_direct_smoothed, n_eff.
#' @return tibble with method, outcome, area_id, replicate, n_eff, mse,
#'   var_direct, var_direct_smoothed, ratio_raw, ratio_smoothed.
#'   Note: `ratio_raw` is set to `NA_real_` when `var_direct == 0` (undefined
#'   direct-variance improvement ratio; this is common for rare outcomes).
#'   `ratio_smoothed` remains valid in these cells since the GVF-smoothed
#'   variance is never exactly zero.
#' @export
compute_direct_mse_ratio <- function(results_panel, direct_panel) {
  results_panel <- as.data.frame(results_panel)
  direct_panel  <- as.data.frame(direct_panel)
  keys <- c("replicate", "outcome", "area_id")
  stopifnot(all(keys %in% colnames(results_panel)),
            all(keys %in% colnames(direct_panel)))
  d <- merge(
    results_panel[, c(keys, "method", "mse")],
    direct_panel[, c(keys, "var_direct", "var_direct_smoothed", "n_eff")],
    by = keys
  )
  d <- data.table::as.data.table(d)
  d[, ratio_raw      := ifelse(var_direct == 0, NA_real_, mse / var_direct)]
  d[, ratio_smoothed := mse / var_direct_smoothed]
  tibble::as_tibble(d)
}

#' Per-area model-vs-direct RMSE reduction, summarised by method x outcome
#'
#' For each area, RMSE_i = sqrt(mean_r (theta_hat - theta_true)^2) over the reps
#' it is in-sample, with the model estimate and the direct estimate PAIRED on
#' the same (replicate, area) draw -- so only municipalities that actually carry
#' a direct estimate that replicate contribute (the pairing IS the sample
#' restriction; no separate draw-count filter is needed or applied).
#' ratio_i = RMSE_model_i / RMSE_direct_i (< 1 means the model reduces error at
#' area i). Summarised per (method, outcome) by the MEDIAN (and mean, for
#' reference) of ratio_i over all areas with RMSE_direct_i > 0 -- the only areas
#' for which the ratio is defined -- plus pct_win = share of those areas with
#' ratio_i < 1.
#'
#' @param results_panel model estimates (pass in-sample, non-failed rows
#'   only): method, replicate, area_id, outcome, theta_hat.
#' @param direct_panel direct estimates row-bound across replicates:
#'   replicate, area_id, outcome, theta_direct.
#' @param truth outcome-long census truth: area_id, outcome, theta_true.
#' @return tibble(method, outcome, median_ratio, mean_ratio, pct_win, n_area,
#'   n_zero_direct). median_ratio/mean_ratio/pct_win are computed over the
#'   n_area - n_zero_direct areas with a positive direct-RMSE denominator.
#' @export
summarise_error_reduction <- function(results_panel, direct_panel, truth) {
  keys <- c("replicate", "area_id", "outcome")
  stopifnot(all(c(keys, "method", "theta_hat") %in% colnames(results_panel)),
            all(c(keys, "theta_direct") %in% colnames(direct_panel)),
            all(c("area_id", "outcome", "theta_true") %in% colnames(truth)))
  d <- merge(as.data.frame(results_panel)[, c(keys, "method", "theta_hat")],
             as.data.frame(direct_panel)[, c(keys, "theta_direct")], by = keys)
  d <- merge(d, as.data.frame(truth)[, c("area_id", "outcome", "theta_true")],
             by = c("area_id", "outcome"))
  d <- data.table::as.data.table(d)
  d[, sq_model  := (theta_hat    - theta_true)^2]
  d[, sq_direct := (theta_direct - theta_true)^2]
  area <- d[, .(rmse_model  = sqrt(mean(sq_model)),
                rmse_direct = sqrt(mean(sq_direct))),
            by = c("method", "outcome", "area_id")]
  area[, ratio_i := rmse_model / rmse_direct]
  out <- area[, {
    ok <- rmse_direct > 0 & is.finite(ratio_i)
    list(median_ratio  = stats::median(ratio_i[ok]),
         mean_ratio    = mean(ratio_i[ok]),
         pct_win       = mean(ratio_i[ok] < 1),
         n_area        = .N,
         n_zero_direct = sum(!ok))
  }, by = c("method", "outcome")]
  tibble::as_tibble(out)
}

#' Method x outcome x n_eff-quartile summary of the mse / direct-variance ratio
#'
#' @param ratio_panel output of compute_direct_mse_ratio()
#' @return tibble(method, outcome, n_eff_quartile, median_ratio_raw,
#'   median_ratio_smoothed, pct_below_1_raw, n)
#' @export
summarise_mse_ratio <- function(ratio_panel) {
  d <- data.table::as.data.table(ratio_panel)
  d[, n_eff_quartile := as.integer(
    cut(n_eff, stats::quantile(n_eff, probs = seq(0, 1, 0.25)),
        include.lowest = TRUE, labels = FALSE)),
    by = .(method, outcome)]
  out <- d[, .(
    median_ratio_raw      = stats::median(ratio_raw, na.rm = TRUE),
    median_ratio_smoothed = stats::median(ratio_smoothed, na.rm = TRUE),
    pct_below_1_raw        = mean(ratio_raw < 1, na.rm = TRUE),
    n                      = .N
  ), by = .(method, outcome, n_eff_quartile)]
  tibble::as_tibble(out)
}
