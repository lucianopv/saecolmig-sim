#' Uncertainty measures of small area estimates
#'
#' Every function takes estimate rows -- one per method x config x outcome x
#' replicate x municipality, with `theta_hat`, `mse`, `ci_lo`, `ci_hi` and
#' `theta_true` -- and returns aggregates only. `mse` is each method's own
#' reported uncertainty: emdi's MSE for FH-arcsin, the posterior variance for
#' Bayes-Beta, povmap's bootstrap variance for MEGB and XGB.
#' @name uncertainty_eval
NULL

#' Grouping keys shared by the uncertainty summaries
#' @export
UR_KEYS <- c("method", "config", "outcome", "area_set")

.ur_z <- stats::qnorm(0.975)

#' Interval score (Gneiting and Raftery 2007): width plus 2/alpha per unit of miss
#' @param lo,hi interval bounds
#' @param truth the true value
#' @param alpha nominal miss rate (0.05 for a 95\% interval)
#' @export
ur_interval_score <- function(lo, hi, truth, alpha = 0.05) {
  (hi - lo) + (2 / alpha) * pmax(lo - truth, 0) + (2 / alpha) * pmax(truth - hi, 0)
}

#' Drop failed and non-converged fits, and say what went
#'
#' A fit is (method, config, outcome, replicate). It is dropped whole
#' when any of its rows is `fit_failed`, or when its `max_rhat` exceeds
#' `rhat_max`; rows with a non-finite estimate, MSE or bound are dropped
#' individually and counted.
#' @param d estimate rows (data.table)
#' @param rhat_max convergence threshold for the Bayesian fits
#' @return list(d = kept rows, dropped = one row per dropped fit with its
#'   reason, rows_nonfinite = count)
#' @export
ur_drop_bad_fits <- function(d, rhat_max = 1.05) {
  fk <- intersect(c("method", "config", "outcome", "replicate"), names(d))
  ff <- if ("fit_failed" %in% names(d)) d$fit_failed %in% TRUE else rep(FALSE, nrow(d))
  rh <- if ("max_rhat" %in% names(d)) !is.na(d$max_rhat) & d$max_rhat > rhat_max else
    rep(FALSE, nrow(d))
  flag <- data.table::copy(d[, fk, with = FALSE])[, `:=`(f = ff, r = rh)][
    , .(failed = any(f), rhat = any(r)), by = fk]
  bad <- flag[failed | rhat][, reason := ifelse(failed, "failed", "rhat")][
    , c(fk, "reason"), with = FALSE]
  kept <- if (nrow(bad)) d[!bad, on = fk] else data.table::copy(d)
  ok <- with(kept, is.finite(theta_hat) & is.finite(mse) & is.finite(ci_lo) & is.finite(ci_hi))
  list(d = kept[ok], dropped = bad, rows_nonfinite = sum(!ok))
}

#' Per-row quantities every summary uses (added by reference)
#' @param d estimate rows (data.table) carrying `theta_true`
#' @param alpha nominal miss rate of the intervals
#' @export
ur_add_row_metrics <- function(d, alpha = 0.05) {
  stopifnot(data.table::is.data.table(d),
            all(c("theta_hat", "mse", "ci_lo", "ci_hi", "theta_true") %in% names(d)))
  d[, `:=`(sq      = (theta_hat - theta_true)^2,
           cov     = as.integer(ci_lo <= theta_true & theta_true <= ci_hi),
           miss_lo = as.integer(theta_true < ci_lo),
           miss_hi = as.integer(theta_true > ci_hi),
           width   = ci_hi - ci_lo,
           iscore  = ur_interval_score(ci_lo, ci_hi, theta_true, alpha))]
  invisible(d)
}

#' Per-municipality summary across replicates
#'
#' `ratio` = mean reported MSE / realised MSE; `rrmse_mse` = RMSE of the
#' reported MSE around the realised MSE, relative to it; `width_rel` = mean
#' width / (2 z realised RMSE), z = `qnorm(0.975)`: the mean width over that of
#' a normal 95\% interval at the realised RMSE, so 1 means intervals as wide as
#' the realised error calls for. All three are NA where the realised MSE is 0
#' (a municipality every replicate estimated exactly), never Inf.
#'
#' The realised MSE also splits as `emp_mse = bias^2 + emp_var`: `bias` is the
#' mean error over replicates and `emp_var` the variance of the estimates
#' around their own mean. `ratio_var` = mean reported uncertainty / `emp_var`
#' (NA where the estimates never move) and `bias_share` = bias^2 / `emp_mse`.
#' They matter because povmap's default MEGB and XGB uncertainty, like a
#' posterior variance, leaves each municipality's own bias out: a method can
#' price its variance well and still look overconfident against the MSE.
#'
#' `rb_rmse` is the relative bias of the estimated RMSE as the MERF/MEGB
#' literature defines it (Krennmair & Schmid 2022; Messer & Schmid 2024):
#' (mean over replicates of the root reported MSE - realised RMSE) / realised
#' RMSE, NA where the realised MSE is 0.
#'
#' When the rows carry `n_dwell` (dwellings drawn), `n_dwell_med` is the
#' municipality's median over its replicates in the set (0 out of sample).
#' @param d rows after `ur_add_row_metrics()`
#' @param by grouping keys
#' @export
ur_area_summary <- function(d, by = UR_KEYS) {
  a <- d[, {
    emp <- mean(sq)
    .(n_reps = .N, theta_true = theta_true[1], emp_mse = emp, rep_mse = mean(mse),
      rep_rmse = mean(sqrt(mse)),
      rrmse_mse = sqrt(mean((mse - emp)^2)) / emp,
      bias = mean(theta_hat - theta_true),
      emp_var = mean((theta_hat - mean(theta_hat))^2),
      coverage = mean(cov), width_mean = mean(width), iscore = mean(iscore))
  }, by = c(by, "area_id")]
  a[, ratio_var := ifelse(emp_var > 0, rep_mse / emp_var, NA_real_)]
  a[, bias_share := ifelse(emp_mse > 0, bias^2 / emp_mse, NA_real_)]
  a[, ratio := ifelse(emp_mse > 0, rep_mse / emp_mse, NA_real_)]
  a[, rrmse_mse := ifelse(emp_mse > 0, rrmse_mse, NA_real_)]
  a[, width_rel := ifelse(emp_mse > 0, width_mean / (2 * .ur_z * sqrt(emp_mse)), NA_real_)]
  a[, rb_rmse := ifelse(emp_mse > 0, rep_rmse / sqrt(emp_mse) - 1, NA_real_)]
  if ("n_dwell" %in% names(d)) {
    nd <- d[, .(n_dwell_med = as.numeric(stats::median(n_dwell))), by = c(by, "area_id")]
    a[nd, n_dwell_med := i.n_dwell_med, on = c(by, "area_id")]
  }
  a[]
}

#' Per-set summary: pooled row metrics plus the per-municipality distribution
#'
#' Per-municipality statistics use only municipalities with at least
#' `min_reps` replicates in the set, because a ratio of two means over a handful
#' of replicates is mostly noise; `n_areas` reports how many qualified.
#' @param d rows after `ur_add_row_metrics()`
#' @param area output of `ur_area_summary()` on the same rows
#' @param by grouping keys
#' @param min_reps minimum replicates for a municipality to enter the
#'   per-municipality statistics
#' @export
ur_set_summary <- function(d, area, by = UR_KEYS, min_reps = 20L) {
  rows <- d[, .(n_rows = .N, coverage = mean(cov), miss_lo = mean(miss_lo),
                miss_hi = mean(miss_hi), width_med = stats::median(width),
                iscore = mean(iscore), ratio_pooled = mean(mse) / mean(sq),
                rmse = sqrt(mean(sq)), rmse_reported = sqrt(mean(mse))), by = by]
  q <- function(x, p) if (all(is.na(x))) NA_real_ else
    stats::quantile(x, p, na.rm = TRUE, names = FALSE)
  sp <- function(x, y) {
    ok <- is.finite(x) & is.finite(y)
    if (sum(ok) < 3L) NA_real_ else stats::cor(x[ok], y[ok], method = "spearman")
  }
  ar <- area[n_reps >= min_reps, .(
    n_areas = .N, n_ratio_undefined = sum(is.na(ratio)),
    ratio_med = q(ratio, 0.5), ratio_q25 = q(ratio, 0.25), ratio_q75 = q(ratio, 0.75),
    share_ratio_lt05 = mean(ratio < 0.5, na.rm = TRUE),
    share_ratio_gt2 = mean(ratio > 2, na.rm = TRUE),
    rrmse_mse_med = q(rrmse_mse, 0.5),
    ratio_var_med = q(ratio_var, 0.5), bias_share_med = q(bias_share, 0.5),
    spearman_rmse = sp(sqrt(rep_mse), sqrt(emp_mse)),
    share_cov_lt090 = mean(coverage < 0.90), share_cov_lt080 = mean(coverage < 0.80),
    width_rel_med = q(width_rel, 0.5),
    rb_rmse_med = q(rb_rmse, 0.5),
    rb_rmse_mean = if (all(is.na(rb_rmse))) NA_real_ else mean(rb_rmse, na.rm = TRUE)), by = by]
  out <- merge(rows, ar, by = by, all.x = TRUE)
  out[is.na(n_areas), n_areas := 0L]
  out[]
}

#' Reliability curve: reported vs realised RMSE by reported-RMSE bin
#' @param d rows after `ur_add_row_metrics()`
#' @param by grouping keys
#' @param n_bins number of quantile bins of the reported RMSE
#' @export
ur_reliability <- function(d, by = UR_KEYS, n_bins = 10L) {
  x <- d[, c(by, "mse", "sq"), with = FALSE]
  x[, rbin := {
    r <- sqrt(mse)
    br <- unique(stats::quantile(r, seq(0, 1, length.out = n_bins + 1L), names = FALSE))
    if (length(br) < 2L) 1L else as.integer(cut(r, br, include.lowest = TRUE, labels = FALSE))
  }, by = by]
  x[, .(n = .N, rmse_reported = sqrt(mean(mse)), rmse_realised = sqrt(mean(sq))),
    by = c(by, "rbin")][order(rbin)]
}

#' Quantile breaks of dwellings drawn, over in-sample rows only
#' @param n_dwell dwellings drawn per municipality and replicate
#' @param n_bins number of in-sample bins
#' @export
ur_size_breaks <- function(n_dwell, n_bins = 5L) {
  x <- n_dwell[is.finite(n_dwell) & n_dwell > 0]
  unique(stats::quantile(x, seq(0, 1, length.out = n_bins + 1L), names = FALSE))
}

#' Sample-size bin: "0 (out of sample)" first, then labelled in-sample bins
#' @param n_dwell dwellings drawn (0 or NA out of sample)
#' @param area_set "in_sample" / "out_of_sample"
#' @param breaks from `ur_size_breaks()`
#' @export
ur_size_bin <- function(n_dwell, area_set, breaks) {
  k <- length(breaks) - 1L
  labs <- sprintf("%d-%d", as.integer(c(breaks[1], breaks[-c(1, k + 1L)] + 1)),
                  as.integer(breaks[-1]))
  b <- as.character(cut(n_dwell, breaks, labels = labs, include.lowest = TRUE))
  b[area_set == "out_of_sample"] <- "0 (out of sample)"
  factor(b, levels = c("0 (out of sample)", labs))
}

#' Coverage, width, pooled MSE ratio and interval score by a bin column
#'
#' Also the bin's empirical RMSE (root mean squared error against the truth),
#' estimated RMSE (root mean of the method's own MSE) and mean true value.
#' @param d rows after `ur_add_row_metrics()`
#' @param bin name of the bin column
#' @param by grouping keys
#' @param x optional numeric column whose median per bin is returned as
#'   `x_med` (the plotting position of a sample-size bin)
#' @export
ur_by_bin <- function(d, bin, by, x = NULL) {
  keys_ <- c(by, bin)   # combined before the data.table call: a column named like an argument would shadow it
  r <- d[, .(n = .N, coverage = mean(cov), width_med = stats::median(width),
             ratio_pooled = mean(mse) / mean(sq), iscore = mean(iscore),
             miss_lo = mean(miss_lo), miss_hi = mean(miss_hi),
             rmse_empirical = sqrt(mean(sq)), rmse_estimated = sqrt(mean(mse)),
             truth_mean = mean(theta_true)), by = keys_]
  if (!is.null(x)) {
    xm <- d[, .(x_med = as.numeric(stats::median(.SD[[1L]]))), by = keys_, .SDcols = x]
    r <- r[xm, on = keys_]
  }
  r
}

#' Publication flags: reported vs realised CV against a threshold
#'
#' Reported CV = sqrt(mse) / theta_hat per row; realised CV = realised RMSE /
#' true value per municipality. A non-positive estimate or truth is never
#' publishable. `false_publish` = share of rows flagged publishable whose
#' realised CV is at or above the threshold.
#' @param d rows after `ur_add_row_metrics()`
#' @param area output of `ur_area_summary()` on the same rows
#' @param threshold CV threshold
#' @param by grouping keys
#' @export
ur_cv_flags <- function(d, area, threshold = 0.20, by = UR_KEYS) {
  tru <- area[, c(by, "area_id", "emp_mse", "theta_true"), with = FALSE]
  tru[, pub_true := theta_true > 0 & sqrt(emp_mse) / theta_true < threshold]
  x <- merge(d[, c(by, "area_id", "theta_hat", "mse"), with = FALSE],
             tru[, c(by, "area_id", "pub_true"), with = FALSE], by = c(by, "area_id"))
  x[, pub_rep := theta_hat > 0 & sqrt(mse) / theta_hat < threshold]
  x[, .(n = .N, share_pub_reported = mean(pub_rep), share_pub_true = mean(pub_true),
        false_publish = if (any(pub_rep)) mean(!pub_true[pub_rep]) else NA_real_,
        missed_publish = if (any(pub_true)) mean(!pub_rep[pub_true]) else NA_real_),
    by = by]
}

#' Per-replicate sums for the Monte Carlo standard errors
#' @param d rows after `ur_add_row_metrics()`
#' @param by grouping keys
#' @export
ur_rep_sums <- function(d, by = UR_KEYS) {
  d[, .(n = .N, s_cov = sum(cov), s_mse = sum(mse), s_sq = sum(sq), s_is = sum(iscore)),
    by = c(by, "replicate")]
}

#' Monte Carlo SEs by replicate bootstrap (whole replicates, common draws)
#' @param sums output of `ur_rep_sums()`
#' @param by grouping keys
#' @param B bootstrap draws
#' @param seed RNG seed
#' @export
ur_mcse <- function(sums, by = UR_KEYS, B = 200L, seed = 20261009L) {
  reps <- sort(unique(sums$replicate))
  set.seed(seed)
  draws <- data.table::rbindlist(lapply(seq_len(B), function(b) {
    w <- data.table::as.data.table(table(replicate = sample(reps, replace = TRUE)))
    w[, replicate := as.integer(as.character(replicate))]
    data.table::setnames(w, "N", "w")
    x <- sums[w, on = "replicate", nomatch = 0L]
    x[, .(coverage = sum(w * s_cov) / sum(w * n),
          ratio_pooled = sum(w * s_mse) / sum(w * s_sq),
          iscore = sum(w * s_is) / sum(w * n)), by = by]
  }))
  draws[, .(mcse_coverage = stats::sd(coverage), mcse_ratio_pooled = stats::sd(ratio_pooled),
            mcse_iscore = stats::sd(iscore)), by = by]
}

#' Every method must see the same drawn municipalities in each replicate
#' @param d rows of several methods for one covariate configuration
#' @export
ur_assert_same_sets <- function(d) {
  fp <- d[area_set == "in_sample", .(fp = paste(sort(unique(area_id)), collapse = ",")),
          by = .(method, outcome, replicate)]
  n <- fp[, .(n_sets = data.table::uniqueN(fp)), by = .(outcome, replicate)]
  if (any(n$n_sets > 1L)) {
    stop("ur_assert_same_sets: in-sample area sets differ across methods in ",
         sum(n$n_sets > 1L), " (outcome, replicate) cells")
  }
  invisible(TRUE)
}
