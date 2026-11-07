#' Compute the census-derived ground truth per area
#'
#' For each area in the frame, returns the unweighted mean of the outcome
#' (`migrated_last_12m` by default). NA outcome rows are excluded from both
#' numerator and denominator.
#'
#' @param frame household-level frame (must contain `migrated_last_12m`)
#' @param area_var character; column name of the area identifier
#' @param outcome_var character; outcome column name
#' @return tibble(area_id, theta_true, n_total)
#' @export
compute_truth <- function(frame,
                          area_var    = "municipality",
                          outcome_var = "migrated_last_12m") {
  dt <- data.table::as.data.table(frame)
  out <- dt[, .(
    theta_true = mean(get(outcome_var), na.rm = TRUE),
    n_total    = sum(!is.na(get(outcome_var)))
  ), by = c(area_var)]
  data.table::setnames(out, area_var, "area_id")
  tibble::as_tibble(out)
}

#' Evaluate SAE estimates against census ground truth
#'
#' @param estimates_long tibble with columns method, replicate, area_id,
#'   theta_hat, mse, ci_lo, ci_hi (and optionally area_set, outcome).
#'   When both estimates_long and truth carry an `outcome` column, the merge
#'   is keyed on (area_id, outcome); grouping summaries include `outcome`.
#' @param truth tibble with columns area_id, theta_true (and optionally
#'   outcome, as produced by compute_truth_all)
#' @return list(by_method, by_method_set, by_method_area), each gaining
#'   rmse_reported (sqrt(mean(mse))) and calibration_ratio
#'   (mean(mse) / empirical mse, pooled ratio-of-means) alongside the existing
#'   bias/rmse/coverage_95/ci_length columns. `by_method_set` additionally
#'   gains calibration_ratio_med: the per-area median of
#'   mean_r(mse) / mean_r(squared error) (each area's estimated vs. realised
#'   MSE across reps), taken over all areas with a positive empirical-MSE
#'   denominator -- a robust, area-democratic companion to the tail-weighted
#'   pooled ratio.
#' @export
evaluate_estimates <- function(estimates_long, truth) {
  estimates_long <- as.data.frame(estimates_long)
  truth <- as.data.frame(truth)
  keys <- intersect(intersect(c("area_id", "outcome"),
                              colnames(estimates_long)),
                    colnames(truth))
  if ("outcome" %in% colnames(truth) && !"outcome" %in% keys) {
    stop("truth is outcome-long but estimates carry no `outcome` column; ",
         "filter truth to one outcome or tag the estimates.")
  }
  if ("outcome" %in% colnames(estimates_long) && !"outcome" %in% keys) {
    stop("estimates are outcome-long but truth carries no `outcome` column; ",
         "pass outcome-long truth (compute_truth_all) or drop `outcome` ",
         "from estimates.")
  }
  d <- merge(estimates_long, truth[, c(keys, "theta_true")], by = keys)
  d <- data.table::as.data.table(d)
  if (!"area_set" %in% colnames(d)) d[, area_set := "in_sample"]
  d[, err    := theta_hat - theta_true]
  d[, sq     := err^2]
  d[, cov    := as.integer(ci_lo <= theta_true & ci_hi >= theta_true)]
  d[, ci_len := ci_hi - ci_lo]

  g_area <- intersect(c("method", "outcome", "area_id", "area_set"),
                      colnames(d))
  g_set  <- intersect(c("method", "outcome", "area_set"), colnames(d))
  g_meth <- intersect(c("method", "outcome"), colnames(d))

  by_method_area <- d[, .(
    bias              = mean(err,  na.rm = TRUE),
    rmse              = sqrt(mean(sq, na.rm = TRUE)),
    rmse_reported     = sqrt(mean(mse, na.rm = TRUE)),
    calibration_ratio = mean(mse, na.rm = TRUE) / mean(sq, na.rm = TRUE),
    coverage_95       = mean(cov,  na.rm = TRUE),
    ci_length         = mean(ci_len, na.rm = TRUE),
    mean_theta_hat    = mean(theta_hat, na.rm = TRUE),
    n_reps            = .N
  ), by = g_area]

  by_method_set <- d[, .(
    bias              = mean(err,  na.rm = TRUE),
    rmse              = sqrt(mean(sq, na.rm = TRUE)),
    rmse_reported     = sqrt(mean(mse, na.rm = TRUE)),
    calibration_ratio = mean(mse, na.rm = TRUE) / mean(sq, na.rm = TRUE),
    mae               = mean(abs(err), na.rm = TRUE),
    coverage_95       = mean(cov, na.rm = TRUE),
    ci_length         = mean(ci_len, na.rm = TRUE),
    n_areas           = data.table::uniqueN(area_id),
    n_reps            = data.table::uniqueN(replicate)
  ), by = g_set]

  # Per-area median calibration (calibration_ratio_med): a robust, area-
  # democratic companion to the pooled `calibration_ratio` above. For each
  # area, c_i = mean_r(mse) / mean_r(sq) -- its mean estimated MSE over its
  # realised (empirical) MSE across the reps it appears in -- summarised by the
  # MEDIAN over all areas with a positive empirical-MSE denominator (the only
  # ones for which the ratio is defined). No draw-count filter is applied, so
  # every sampled municipality is represented, not just the frequently-drawn
  # (large) ones -- the median is already robust to the handful of areas drawn
  # in few reps. Unlike the pooled ratio-of-means (dominated by the high-
  # variance tail areas), and unlike the mean-of-ratios (which blows up on the
  # rare outcomes' small denominators), the median tracks the *typical* area.
  area_cal <- d[is.finite(sq) & is.finite(mse),
                .(mse_bar = mean(mse), emp = mean(sq)),
                by = c(g_set, "area_id")]
  area_cal[, c_i := mse_bar / emp]
  cal_med_set <- area_cal[emp > 0 & is.finite(c_i),
                          .(calibration_ratio_med = stats::median(c_i)),
                          by = g_set]
  by_method_set <- merge(by_method_set, cal_med_set, by = g_set, all.x = TRUE)

  rank_corr <- d[, .(rho = stats::cor(theta_hat, theta_true,
                                       method = "spearman")),
                 by = c(g_meth, "replicate")]
  rank_summary <- rank_corr[, .(rank_corr_mean = mean(rho, na.rm = TRUE)),
                             by = g_meth]

  by_method <- d[, .(
    bias              = mean(err,  na.rm = TRUE),
    rmse              = sqrt(mean(sq, na.rm = TRUE)),
    rmse_reported     = sqrt(mean(mse, na.rm = TRUE)),
    calibration_ratio = mean(mse, na.rm = TRUE) / mean(sq, na.rm = TRUE),
    mae               = mean(abs(err), na.rm = TRUE),
    coverage_95       = mean(cov, na.rm = TRUE),
    ci_length         = mean(ci_len, na.rm = TRUE)
  ), by = g_meth]
  by_method <- merge(by_method, rank_summary, by = g_meth)

  list(by_method      = tibble::as_tibble(by_method),
       by_method_set  = tibble::as_tibble(by_method_set),
       by_method_area = tibble::as_tibble(by_method_area))
}

#' Census ratio truth per area (person-level outcomes)
#'
#' theta_true = sum(numerator) / sum(denominator) over the frame.
#' @param frame household-level frame with the count columns
#' @param area_var character; area identifier column
#' @param numerator_var,denominator_var character; count columns
#' @return tibble(area_id, theta_true, n_total); n_total = sum(denominator)
#' @export
compute_truth_ratio <- function(frame,
                                area_var = "municipality",
                                numerator_var,
                                denominator_var) {
  dt <- data.table::as.data.table(frame)
  out <- dt[, .(
    theta_true = sum(get(numerator_var), na.rm = TRUE) /
                 sum(get(denominator_var), na.rm = TRUE),
    n_total    = sum(get(denominator_var), na.rm = TRUE)
  ), by = c(area_var)]
  data.table::setnames(out, area_var, "area_id")
  tibble::as_tibble(out)
}

#' Census truth for all registered outcomes, stacked long
#' @export
compute_truth_all <- function(frame,
                              area_var = "municipality",
                              outcomes = sae_outcomes) {
  parts <- lapply(seq_len(nrow(outcomes)), function(i) {
    o <- as.list(outcomes[i, ])
    t <- if (o$level == "household") {
      compute_truth(frame, area_var = area_var, outcome_var = o$indicator)
    } else {
      compute_truth_ratio(frame, area_var = area_var,
                          numerator_var   = o$numerator,
                          denominator_var = o$denominator)
    }
    t$outcome <- o$outcome
    t
  })
  tibble::as_tibble(do.call(rbind, parts))
}
