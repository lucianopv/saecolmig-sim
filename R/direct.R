#' Design-based direct estimator (Hajek + Taylor variance) per area
#'
#' Builds a survey design from the sample (strata `stratum`, clusters
#' `ssu_id`, weights `weight_cal`), then computes a per-area Hajek proportion
#' and Taylor variance
#' for the outcome. Effective sample size n_eff = (sum w)^2 / sum(w^2);
#' n_dwell = number of sample rows in the area.
#'
#' @param sample tibble from draw_replicate + apply_calibration
#' @param area_var character; area identifier column (default "municipality")
#' @param outcome_var character (default "migrated_last_12m")
#' @return tibble(area_id, theta_direct, var_direct, n_eff, n_dwell)
#' @export
direct_estimate <- function(sample,
                            area_var    = "municipality",
                            outcome_var = "migrated_last_12m") {
  s <- as.data.frame(sample)
  s$.area    <- s[[area_var]]
  s$.outcome <- s[[outcome_var]]
  # The design clusters on the segment (ids = ~ssu_id), the unit selected
  # within a municipality, so a lonely PSU here is a stratum in which a single
  # segment is drawn. Its Taylor variance contribution is not defined;
  # survey.lonely.psu = "average" gives such a stratum the average variance
  # contribution of the other strata.
  prev <- getOption("survey.lonely.psu", "fail")
  options(survey.lonely.psu = "average")
  on.exit(options(survey.lonely.psu = prev), add = TRUE)
  des <- survey::svydesign(ids = ~ssu_id, strata = ~stratum,
                           weights = ~weight_cal, nest = TRUE, data = s)
  by  <- survey::svyby(~.outcome, ~.area, design = des,
                       FUN = survey::svymean, na.rm = TRUE)
  out <- data.frame(
    area_id      = by$.area,
    theta_direct = by$.outcome,
    var_direct   = by$se^2,
    stringsAsFactors = FALSE
  )

  # n_eff and n_dwell aggregates per area
  s_dt <- data.table::as.data.table(s)
  agg  <- s_dt[, .(
    n_eff   = (sum(weight_cal[!is.na(.outcome)]))^2 /
              sum(weight_cal[!is.na(.outcome)]^2),
    n_dwell = .N
  ), by = .area]
  data.table::setnames(agg, ".area", "area_id")
  out <- merge(out, agg, by = "area_id", all.x = TRUE)
  tibble::as_tibble(out)
}

#' Smooth direct sampling variances across areas
#'
#' Two backends:
#'   "ghs" (default): a Generalised Heteroskedastic Smoother. Fits
#'     log(var_direct) ~ log(n_eff) + theta_direct + theta_direct^2 via OLS;
#'     predicted log-variances are exponentiated to give strictly positive
#'     smoothed variances. Follows the spirit of You & Hidiroglou (2023)
#'     while staying lightweight (no extra dependency beyond stats).
#'   "tipsae": delegates to tipsae::smoothing() when the package is
#'     installed; falls back to "ghs" with a warning otherwise.
#'
#' @param direct output of direct_estimate()
#' @param method one of "ghs", "tipsae"
#' @return direct with new column var_direct_smoothed (strictly positive for
#'   areas with a finite direct estimate; NA otherwise)
#' @export
smooth_variance <- function(direct, method = c("ghs", "tipsae")) {
  method <- match.arg(method)
  d <- direct
  ok <- with(d, var_direct > 0 & is.finite(var_direct) &
                n_eff > 0 & is.finite(n_eff))
  if (sum(ok) < 5L) {
    floor_var <- max(min(d$var_direct[ok], na.rm = TRUE), 1e-8)
    d$var_direct_smoothed <- ifelse(ok, d$var_direct, floor_var)
    d$var_direct_smoothed[!is.finite(d$theta_direct)] <- NA_real_
    return(tibble::as_tibble(d))
  }

  if (method == "tipsae") {
    if (!requireNamespace("tipsae", quietly = TRUE)) {
      warning("tipsae not installed; falling back to GHS smoother.")
    } else {
      smoothed <- try(tipsae::smoothing(formula  = var_direct ~ n_eff,
                                        data     = d[ok, ],
                                        method_var = "GVF"),
                      silent = TRUE)
      if (!inherits(smoothed, "try-error")) {
        d$var_direct_smoothed <- NA_real_
        d$var_direct_smoothed[ok]  <- smoothed$smoothed_var
        d$var_direct_smoothed[!ok] <- min(d$var_direct_smoothed[ok])
        return(tibble::as_tibble(d))
      }
      warning("tipsae::smoothing failed; falling back to GHS smoother.")
    }
  }

  fit <- stats::lm(log(var_direct) ~ log(n_eff) + theta_direct +
                                       I(theta_direct^2),
                   data = d[ok, ])
  d$var_direct_smoothed <- exp(stats::predict(fit, newdata = d))
  d$var_direct_smoothed[!is.finite(d$theta_direct)] <- NA_real_
  tibble::as_tibble(d)
}

#' Design-based ratio estimator (Hajek ratio + Taylor variance) per area
#'
#' For person-level outcomes: theta = sum(w * numerator) / sum(w *
#' denominator) via survey::svyratio under the same SSU-clustered design as
#' direct_estimate(). Areas whose weighted denominator is zero return NA
#' (count reported via message()). n_eff is computed on denominator
#' contributors (rows with denominator > 0); n_dwell counts all rows.
#'
#' @param sample tibble from draw_replicate + apply_calibration
#' @param area_var character; area identifier column (default "municipality")
#' @param numerator_var,denominator_var character; count columns
#' @return tibble(area_id, theta_direct, var_direct, n_eff, n_dwell)
#' @export
direct_estimate_ratio <- function(sample,
                                  area_var = "municipality",
                                  numerator_var,
                                  denominator_var) {
  s <- as.data.frame(sample)
  s$.area <- s[[area_var]]
  s$.num  <- s[[numerator_var]]
  s$.den  <- s[[denominator_var]]
  prev <- getOption("survey.lonely.psu", "fail")
  options(survey.lonely.psu = "average")
  on.exit(options(survey.lonely.psu = prev), add = TRUE)
  des <- survey::svydesign(ids = ~ssu_id, strata = ~stratum,
                           weights = ~weight_cal, nest = TRUE, data = s)
  by <- survey::svyby(~.num, by = ~.area, design = des,
                      FUN = survey::svyratio, denominator = ~.den,
                      na.rm = TRUE)
  out <- data.frame(
    area_id      = by$.area,
    theta_direct = as.numeric(stats::coef(by)),
    var_direct   = as.numeric(survey::SE(by))^2,
    stringsAsFactors = FALSE
  )
  bad <- !is.finite(out$theta_direct) | !is.finite(out$var_direct)
  if (any(bad)) {
    out$theta_direct[bad] <- NA_real_
    out$var_direct[bad]   <- NA_real_
    message(sum(bad), " area(s) with zero/NA weighted denominator -> NA (",
            numerator_var, "/", denominator_var, ")")
  }
  s_dt <- data.table::as.data.table(s)
  agg  <- s_dt[, .(
    n_eff   = (sum(weight_cal[!is.na(.den) & .den > 0]))^2 /
              sum(weight_cal[!is.na(.den) & .den > 0]^2),
    n_dwell = .N
  ), by = .area]
  data.table::setnames(agg, ".area", "area_id")
  out <- merge(out, agg, by = "area_id", all.x = TRUE)
  tibble::as_tibble(out)
}

#' Direct estimates for all registered outcomes, stacked long
#'
#' Dispatches each registry row to direct_estimate() (household indicators)
#' or direct_estimate_ratio() (person ratios) and stacks with an `outcome`
#' column.
#'
#' @param sample tibble from draw_replicate + apply_calibration
#' @param area_var character (default "municipality")
#' @param outcomes tibble; defaults to the full sae_outcomes registry
#' @return tibble(area_id, theta_direct, var_direct, n_eff, n_dwell, outcome)
#' @export
direct_estimate_all <- function(sample,
                                area_var = "municipality",
                                outcomes = sae_outcomes) {
  parts <- lapply(seq_len(nrow(outcomes)), function(i) {
    o <- as.list(outcomes[i, ])
    d <- if (o$level == "household") {
      direct_estimate(sample, area_var = area_var,
                      outcome_var = o$indicator)
    } else {
      direct_estimate_ratio(sample, area_var = area_var,
                            numerator_var   = o$numerator,
                            denominator_var = o$denominator)
    }
    d$outcome <- o$outcome
    d
  })
  tibble::as_tibble(do.call(rbind, parts))
}

#' Apply smooth_variance() independently within each outcome
#'
#' One GVF per outcome — pooling outcomes with different prevalence scales
#' would bias the log-linear fit. Areas with NA direct estimates keep NA
#' smoothed variances.
#'
#' @param direct_long output of direct_estimate_all()
#' @param method passed to smooth_variance()
#' @return direct_long with var_direct_smoothed
#' @export
smooth_variance_all <- function(direct_long, method = "ghs") {
  parts <- lapply(split(direct_long, direct_long$outcome),
                  smooth_variance, method = method)
  out <- do.call(rbind, parts)
  rownames(out) <- NULL
  tibble::as_tibble(out)
}
