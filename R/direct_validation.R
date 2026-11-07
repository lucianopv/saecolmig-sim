# R/direct_validation.R
#
# Inferential validation of the design-based direct estimator across
# Monte Carlo replicates of the ECV 2024 design draw: bias, RMSE, standard-error
# ratio and Wald coverage at the national, department and department x area
# levels, where the design gives the direct estimator enough sample.

#' Direct estimator at national and dept x area levels
#'
#' Builds a single survey design from the calibrated sample and returns the
#' Hajek mean (or weighted ratio) and Taylor-linearised SE of the outcome at
#' the national level and at each dept x area cell. Matches the survey design
#' conventions of R/direct.R (single-stage SSU clustering with
#' lonely.psu = "average").
#'
#' @param sample tibble or data.frame from draw_replicate + apply_calibration,
#'   must contain columns dept, area, stratum, ssu_id, weight_cal, and either
#'   the outcome column (indicator path) or the numerator/denominator columns
#'   (ratio path).
#' @param outcome character; name of the binary outcome column. Ignored when
#'   \code{numerator} is supplied.
#' @param numerator character or NULL; column name of the ratio numerator.
#'   When non-NULL, activates the ratio path (svyratio) and \code{outcome} is
#'   ignored.
#' @param denominator character or NULL; column name of the ratio denominator.
#'   Must be non-NULL when \code{numerator} is non-NULL.
#' @return tibble with columns
#'   level   (character: "national" or "dept_x_area"),
#'   label   (character: "national" or "<dept>_<area>"),
#'   estimate (numeric),
#'   se      (numeric),
#'   n_obs   (integer),
#'   df      (numeric; degrees of freedom from svymean / svyratio / svyby).
#' @export
direct_two_levels <- function(sample, outcome = "migrated_last_12m",
                               numerator = NULL, denominator = NULL) {
  is_ratio <- !is.null(numerator)

  s <- as.data.frame(sample)

  # Lonely-PSU treatment — survey default "fail" is too strict for small frames
  prev <- getOption("survey.lonely.psu", "fail")
  options(survey.lonely.psu = "average")
  on.exit(options(survey.lonely.psu = prev), add = TRUE)

  if (!is_ratio) {
    # ---- indicator (Hajek mean) path ----------------------------------------
    stopifnot(outcome %in% colnames(s))
    s$.outcome <- s[[outcome]]

    des <- survey::svydesign(ids = ~ssu_id, strata = ~stratum,
                             weights = ~weight_cal, nest = TRUE, data = s)

    nat <- survey::svymean(~.outcome, des, na.rm = TRUE)
    nat_row <- data.frame(
      level    = "national",
      label    = "national",
      estimate = as.numeric(nat),
      se       = as.numeric(survey::SE(nat)),
      n_obs    = sum(!is.na(s$.outcome)),
      df       = survey::degf(des),
      stringsAsFactors = FALSE
    )

    by <- survey::svyby(~.outcome, ~dept + area, design = des,
                        FUN = survey::svymean, na.rm = TRUE)
    da_n <- as.data.frame(
      data.table::as.data.table(s)[!is.na(.outcome),
                                    .(n_obs = .N), by = .(dept, area)]
    )
    da <- merge(by, da_n, by = c("dept", "area"), all.x = TRUE)
    da_row <- data.frame(
      level    = "dept_x_area",
      label    = paste(da$dept, da$area, sep = "_"),
      estimate = da$.outcome,
      se       = da$se,
      n_obs    = da$n_obs,
      df       = survey::degf(des),
      stringsAsFactors = FALSE
    )

    dept_by <- survey::svyby(~.outcome, ~dept, design = des,
                             FUN = survey::svymean, na.rm = TRUE)
    dept_n <- as.data.frame(
      data.table::as.data.table(s)[!is.na(.outcome),
                                    .(n_obs = .N), by = .(dept)]
    )
    dpt <- merge(dept_by, dept_n, by = "dept", all.x = TRUE)
    dept_row <- data.frame(
      level    = "dept",
      label    = as.character(dpt$dept),
      estimate = dpt$.outcome,
      se       = dpt$se,
      n_obs    = dpt$n_obs,
      df       = survey::degf(des),
      stringsAsFactors = FALSE
    )

  } else {
    # ---- ratio (svyratio) path -----------------------------------------------
    stopifnot(
      !is.null(denominator),
      numerator   %in% colnames(s),
      denominator %in% colnames(s)
    )
    s$.num <- s[[numerator]]
    s$.den <- s[[denominator]]

    des <- survey::svydesign(ids = ~ssu_id, strata = ~stratum,
                             weights = ~weight_cal, nest = TRUE, data = s)

    nat <- survey::svyratio(~.num, ~.den, des, na.rm = TRUE)
    nat_row <- data.frame(
      level    = "national",
      label    = "national",
      estimate = as.numeric(stats::coef(nat)),
      se       = as.numeric(survey::SE(nat)),
      n_obs    = sum(!is.na(s$.den) & s$.den > 0),
      df       = survey::degf(des),
      stringsAsFactors = FALSE
    )

    by <- survey::svyby(~.num, ~dept + area, design = des,
                        FUN = survey::svyratio, denominator = ~.den,
                        na.rm = TRUE)
    da_n <- as.data.frame(
      data.table::as.data.table(s)[!is.na(.den) & .den > 0,
                                   .(n_obs = .N), by = .(dept, area)]
    )
    da <- data.frame(dept = by$dept, area = by$area,
                     est  = as.numeric(stats::coef(by)),
                     se   = as.numeric(survey::SE(by)),
                     stringsAsFactors = FALSE)
    da <- merge(da, da_n, by = c("dept", "area"), all.x = TRUE)
    da_row <- data.frame(
      level    = "dept_x_area",
      label    = paste(da$dept, da$area, sep = "_"),
      estimate = da$est,
      se       = da$se,
      n_obs    = da$n_obs,
      df       = survey::degf(des),
      stringsAsFactors = FALSE
    )

    dept_by <- survey::svyby(~.num, ~dept, design = des,
                             FUN = survey::svyratio, denominator = ~.den,
                             na.rm = TRUE)
    dept_n <- as.data.frame(
      data.table::as.data.table(s)[!is.na(.den) & .den > 0,
                                   .(n_obs = .N), by = .(dept)]
    )
    dpt <- data.frame(dept = dept_by$dept,
                      est  = as.numeric(stats::coef(dept_by)),
                      se   = as.numeric(survey::SE(dept_by)),
                      stringsAsFactors = FALSE)
    dpt <- merge(dpt, dept_n, by = "dept", all.x = TRUE)
    dept_row <- data.frame(
      level    = "dept",
      label    = as.character(dpt$dept),
      estimate = dpt$est,
      se       = dpt$se,
      n_obs    = dpt$n_obs,
      df       = survey::degf(des),
      stringsAsFactors = FALSE
    )
  }

  tibble::as_tibble(rbind(nat_row, da_row, dept_row))
}

#' Census-derived truth at national and dept x area levels
#'
#' Computes the unweighted proportion of `outcome == 1` on the household-level
#' frame restricted to `in_sample_munis`, at the national and dept x area
#' levels used by direct_two_levels(). No row-level data is returned.
#'
#' @param frame household-level data.frame/tibble with dept, area,
#'   municipality, and the outcome column.
#' @param in_sample_munis character or integer vector of municipality codes
#'   to retain.
#' @param outcome character; column name of the binary outcome.
#' @return tibble(level, label, theta_true).
#' @export
compute_truth_two_levels <- function(frame, in_sample_munis,
                                     outcome = "migrated_last_12m") {
  stopifnot(outcome %in% colnames(frame))
  fr <- data.table::as.data.table(frame)
  fr <- fr[municipality %in% in_sample_munis]
  y  <- fr[[outcome]]

  nat <- data.frame(level = "national", label = "national",
                    theta_true = mean(y, na.rm = TRUE),
                    stringsAsFactors = FALSE)

  da <- fr[, .(theta_true = mean(get(outcome), na.rm = TRUE)),
           by = .(dept, area)]
  da_row <- data.frame(
    level      = "dept_x_area",
    label      = paste(da$dept, da$area, sep = "_"),
    theta_true = da$theta_true,
    stringsAsFactors = FALSE
  )

  dept_tb <- fr[, .(theta_true = mean(get(outcome), na.rm = TRUE)),
                by = .(dept)]
  dept_row <- data.frame(
    level      = "dept",
    label      = as.character(dept_tb$dept),
    theta_true = dept_tb$theta_true,
    stringsAsFactors = FALSE
  )
  tibble::as_tibble(rbind(nat, da_row, dept_row))
}

#' direct_two_levels for every registered outcome, stacked long
#'
#' Iterates over all rows of \code{outcomes} (default: \code{sae_outcomes}),
#' calls \code{direct_two_levels()} with the appropriate path (indicator or
#' ratio), and row-binds results with an \code{outcome} column prepended.
#'
#' @param sample tibble or data.frame from \code{draw_replicate} +
#'   \code{apply_calibration}.
#' @param outcomes data.frame/tibble with the column contract of
#'   \code{sae_outcomes}: outcome, level, indicator, numerator, denominator.
#' @return tibble with all columns of \code{direct_two_levels()} plus
#'   \code{outcome}.
#' @export
direct_two_levels_all <- function(sample, outcomes = sae_outcomes) {
  parts <- lapply(seq_len(nrow(outcomes)), function(i) {
    o <- as.list(outcomes[i, ])
    r <- if (o$level == "household") {
      direct_two_levels(sample, outcome = o$indicator)
    } else {
      direct_two_levels(sample, numerator = o$numerator,
                        denominator = o$denominator)
    }
    r$outcome <- o$outcome
    r
  })
  tibble::as_tibble(do.call(rbind, parts))
}

#' compute_truth_two_levels for every registered outcome, stacked long
#'
#' Iterates over all rows of \code{outcomes} (default: \code{sae_outcomes}),
#' computes census truth for each, and row-binds results with an \code{outcome}
#' column prepended.
#'
#' @param frame household-level data.frame/tibble with dept, area, municipality,
#'   and all indicator/numerator/denominator columns.
#' @param in_sample_munis character or integer vector of municipality codes to
#'   restrict the frame to before computing truth.
#' @param outcomes data.frame/tibble with the column contract of
#'   \code{sae_outcomes}.
#' @return tibble with columns level, label, theta_true, outcome.
#' @export
compute_truth_two_levels_all <- function(frame, in_sample_munis,
                                          outcomes = sae_outcomes) {
  fr <- data.table::as.data.table(frame)
  fr <- fr[municipality %in% in_sample_munis]
  parts <- lapply(seq_len(nrow(outcomes)), function(i) {
    o <- as.list(outcomes[i, ])
    if (o$level == "household") {
      t <- compute_truth_two_levels(fr, in_sample_munis = in_sample_munis,
                                    outcome = o$indicator)
    } else {
      nat <- data.frame(
        level      = "national",
        label      = "national",
        theta_true = sum(fr[[o$numerator]], na.rm = TRUE) /
                     sum(fr[[o$denominator]], na.rm = TRUE),
        stringsAsFactors = FALSE
      )
      da <- fr[, .(theta_true = sum(get(o$numerator),   na.rm = TRUE) /
                                sum(get(o$denominator), na.rm = TRUE)),
               by = .(dept, area)]
      da_row <- data.frame(
        level      = "dept_x_area",
        label      = paste(da$dept, da$area, sep = "_"),
        theta_true = da$theta_true,
        stringsAsFactors = FALSE
      )
      dept_tb <- fr[, .(theta_true = sum(get(o$numerator),   na.rm = TRUE) /
                                     sum(get(o$denominator), na.rm = TRUE)),
                    by = .(dept)]
      dept_row <- data.frame(
        level      = "dept",
        label      = as.character(dept_tb$dept),
        theta_true = dept_tb$theta_true,
        stringsAsFactors = FALSE
      )
      t <- tibble::as_tibble(rbind(nat, da_row, dept_row))
    }
    t$outcome <- o$outcome
    t
  })
  tibble::as_tibble(do.call(rbind, parts))
}

#' Monte Carlo summary metrics per cell
#'
#' Given a long tibble of per-replicate direct estimates (one row per
#' (rep_id, level, label)) and a truth tibble, computes the cell-level
#' frequentist summary: empirical bias, RMSE, empirical SE, mean estimated SE,
#' the SE ratio, Wald-CI coverage, MCSE of bias, and the participation rate
#' (fraction of reps in which the cell was estimable, i.e. non-NA estimate
#' and SE > 0).
#'
#' If \code{direct_long} contains an \code{outcome} column the function
#' groups by (outcome, level, label); otherwise it groups by (level, label).
#'
#' @param direct_long tibble with rep_id, level, label, estimate, se.
#'   An optional \code{outcome} column triggers multi-outcome grouping.
#' @param truth tibble with level, label, theta_true. Should also contain
#'   \code{outcome} when \code{direct_long} does.
#' @param alpha numeric; significance level for the Wald CI (default 0.05).
#' @return tibble keyed by (outcome?, level, label) with summary metrics.
#' @export
summarise_mc_metrics <- function(direct_long, truth, alpha = 0.05) {
  stopifnot(all(c("rep_id", "level", "label", "estimate", "se") %in%
                colnames(direct_long)))
  z <- stats::qnorm(1 - alpha / 2)
  d <- data.table::as.data.table(direct_long)
  t <- data.table::as.data.table(truth)

  # Build merge/group keys — include "outcome" only when present in d
  keys      <- intersect(c("outcome", "level", "label"), colnames(d))
  merge_by  <- intersect(keys, colnames(t))

  d <- merge(d, t, by = merge_by, all.x = TRUE)

  d[, estimable := !is.na(estimate) & !is.na(se) & se > 0]
  out <- d[, {
    R_total <- .N
    R_est   <- sum(estimable)
    if (R_est == 0) {
      list(R_total = R_total, R_estimable = 0L,
           participation = 0, bias = NA_real_, rmse = NA_real_,
           emp_se = NA_real_, mean_est_se = NA_real_, se_ratio = NA_real_,
           coverage = NA_real_, mcse_bias = NA_real_)
    } else {
      est  <- estimate[estimable]
      se_v <- se[estimable]
      tv   <- theta_true[estimable][1]
      bias <- mean(est) - tv
      rmse <- sqrt(mean((est - tv)^2))
      emp  <- stats::sd(est)
      mest <- mean(se_v)
      cov  <- mean(abs(est - tv) <= z * se_v)
      list(R_total = R_total, R_estimable = R_est,
           participation = R_est / R_total,
           bias = bias, rmse = rmse,
           emp_se = emp, mean_est_se = mest,
           se_ratio = if (!is.na(emp) && emp > 0) mest / emp else NA_real_,
           coverage = cov,
           mcse_bias = emp / sqrt(R_est))
    }
  }, by = keys]
  tibble::as_tibble(out)
}

#' Monte Carlo running diagnostic for the national-level direct estimator
#'
#' For each R in R_grid, computes the running bias, empirical SE, Wald-CI
#' coverage, and MCSE of bias using the first R replicates. It shows whether
#' the number of replicates resolves the bias and the coverage to the
#' precision a comparison needs.
#'
#' If \code{direct_long} contains an \code{outcome} column the function
#' iterates over outcomes and prepends an \code{outcome} column to the result.
#'
#' @param direct_long tibble; assumed to contain \code{level == "national"}
#'   rows ordered by rep_id.
#' @param truth tibble with level, label, theta_true (plus outcome when
#'   direct_long carries it).
#' @param R_grid integer vector of cumulative R values to evaluate.
#' @param alpha significance level (default 0.05).
#' @return tibble with R, running_bias, running_emp_se, running_coverage,
#'   running_mcse (and outcome when present in the input).
#' @export
mc_running_trace <- function(direct_long, truth,
                             R_grid = c(50, 100, 200, 500),
                             alpha = 0.05) {
  if ("outcome" %in% colnames(direct_long)) {
    outcomes <- unique(direct_long$outcome)
    parts <- lapply(outcomes, function(oc) {
      sub_d <- direct_long[direct_long$outcome == oc, ]
      sub_t <- truth[truth$outcome == oc, ]
      r <- .mc_running_trace_one(sub_d, sub_t, R_grid = R_grid, alpha = alpha)
      r$outcome <- oc
      r
    })
    return(tibble::as_tibble(do.call(rbind, parts)))
  }
  .mc_running_trace_one(direct_long, truth, R_grid = R_grid, alpha = alpha)
}

# Internal single-outcome worker for mc_running_trace
.mc_running_trace_one <- function(direct_long, truth,
                                   R_grid = c(50, 100, 200, 500),
                                   alpha = 0.05) {
  d <- data.table::as.data.table(direct_long)[level == "national"]
  data.table::setorder(d, rep_id)
  tv <- truth$theta_true[truth$level == "national" & truth$label == "national"]
  z  <- stats::qnorm(1 - alpha / 2)

  rows <- lapply(R_grid, function(R) {
    sub <- d[seq_len(min(R, .N))]
    est <- sub$estimate
    se  <- sub$se
    keep <- !is.na(est) & !is.na(se) & se > 0
    est <- est[keep]; se <- se[keep]
    if (length(est) == 0) {
      return(data.frame(R = R, running_bias = NA_real_,
                        running_emp_se = NA_real_,
                        running_coverage = NA_real_,
                        running_mcse = NA_real_))
    }
    emp <- if (length(est) > 1) stats::sd(est) else NA_real_
    data.frame(
      R                = R,
      running_bias     = mean(est) - tv,
      running_emp_se   = emp,
      running_coverage = mean(abs(est - tv) <= z * se),
      running_mcse     = if (!is.na(emp)) emp / sqrt(length(est)) else NA_real_,
      stringsAsFactors = FALSE
    )
  })
  tibble::as_tibble(do.call(rbind, rows))
}

#' Distributional summary of dept x area cell metrics
#'
#' Reduces the 66 dept x area cells to a compact median/IQR/min/max table
#' for the four headline metrics (bias, rmse, se_ratio, coverage),
#' optionally split by area (`"all"`, `"urban"` or `"rural"`). National rows
#' are dropped.
#'
#' If \code{eval_summary} contains an \code{outcome} column the function
#' iterates over outcomes and prepends an \code{outcome} column to the result.
#'
#' @param eval_summary output of summarise_mc_metrics() (or any tibble with
#'   columns level, label, bias, rmse, se_ratio, coverage).
#' @return tibble(outcome?, metric, area_split, median, q25, q75, min, max).
#' @export
cell_dist_summary <- function(eval_summary) {
  metrics <- c("bias", "rmse", "se_ratio", "coverage")
  stopifnot(all(metrics %in% colnames(eval_summary)))

  if ("outcome" %in% colnames(eval_summary)) {
    outcomes <- unique(eval_summary$outcome)
    parts <- lapply(outcomes, function(oc) {
      sub <- eval_summary[eval_summary$outcome == oc, ]
      r <- .cell_dist_summary_one(sub)
      r$outcome <- oc
      r
    })
    return(tibble::as_tibble(do.call(rbind, parts)))
  }
  .cell_dist_summary_one(eval_summary)
}

# Internal single-outcome worker for cell_dist_summary
.cell_dist_summary_one <- function(eval_summary) {
  metrics <- c("bias", "rmse", "se_ratio", "coverage")
  d <- data.table::as.data.table(eval_summary)[level == "dept_x_area"]
  d[, area := sub("^.*_", "", label)]

  splits <- list(all = d, urban = d[area == "urban"], rural = d[area == "rural"])

  rows <- list()
  for (split_name in names(splits)) {
    sub <- splits[[split_name]]
    for (m in metrics) {
      v <- sub[[m]]
      v <- v[!is.na(v)]
      rows[[length(rows) + 1L]] <- data.frame(
        metric     = m,
        area_split = split_name,
        median     = if (length(v)) stats::median(v) else NA_real_,
        q25        = if (length(v)) stats::quantile(v, 0.25, names = FALSE) else NA_real_,
        q75        = if (length(v)) stats::quantile(v, 0.75, names = FALSE) else NA_real_,
        min        = if (length(v)) min(v) else NA_real_,
        max        = if (length(v)) max(v) else NA_real_,
        stringsAsFactors = FALSE
      )
    }
  }
  tibble::as_tibble(do.call(rbind, rows))
}
