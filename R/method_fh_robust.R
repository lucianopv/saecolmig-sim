#' Fit a Sinha-Rao robust Fay-Herriot model via saeRobust::rfh
#'
#' Backend wraps `saeRobust::rfh()` with Huber's tuning constant k = 1.345
#' (95\% efficiency at the normal). The Huber influence function bounds the
#' weight of areas whose direct estimate lies far from the regression fit.
#'
#' Out-of-sample handling: `saeRobust::rfh()` requires complete cases and
#' `predict.fitrfh()` has no `newdata` argument. OOS rows are predicted
#' synthetically as `X_oos \%*\% beta_hat`, with the leading-term
#' synthetic-predictor MSE `sigma2_v + x_i' Var(beta_hat) x_i` (Rao & Molina
#' 2015, sec. 6.2.2). saeRobust exposes no `vcov()`, so `Var(beta_hat)` is
#' approximated by the GLS coefficient covariance `(X' V^-1 X)^-1`,
#' `V = diag(vardir + sigma2_v)`, at the robust variance estimate. The
#' leading-term form omits the (typically minor) g3 variance-of-variance
#' component a parametric bootstrap would capture.
#'
#' @param direct tibble (area_id, theta_direct, var_direct_smoothed, n_eff)
#' @param X tibble (municipality, x1, ..., xK) for in-sample areas
#' @param geo ignored
#' @param X_oos optional tibble (municipality, x1, ..., xK) for held-out areas
#' @param k Huber tuning constant; default 1.345 (saeRobust default)
#' @param B bootstrap reps for MSE; default 100 (saeRobust default)
#' @param standardize logical; centre and scale the covariates before calling
#'   saeRobust::rfh(), then back-transform the reported coefficients. Default
#'   TRUE. FALSE fits on the covariates as given.
#'
#'   Why this is needed. A covariate matrix that mixes columns whose norms span
#'   orders of magnitude (population counts against shares near one) has a large
#'   condition number. rfh() forms and inverts a design-derived cross-product,
#'   which squares it, so solve() can refuse a full-rank matrix as
#'   "computationally singular": a problem of conditioning, not of rank.
#'   Standardising brings the condition number down to that of the shares.
#'
#'   The transformation is affine, so a linear model is invariant to it: fitted
#'   values, MSEs and predictions are unchanged and only the factorisation
#'   succeeds or fails. Coefficients are back-transformed below so the reported
#'   beta_* diagnostics stay on the original covariate scale.
#' @param ... ignored
#' @return list(estimates, diagnostics, meta) with class c("sae_fit", "sae_fit_fh_robust")
#' @export
fit_fh_robust <- function(direct, X, geo = NULL, X_oos = NULL,
                          k = 1.345, B = 100L, standardize = TRUE, ...) {
  if (!requireNamespace("saeRobust", quietly = TRUE)) {
    stop("Package 'saeRobust' is required for fit_fh_robust().")
  }
  t0 <- Sys.time()

  x_cols <- setdiff(colnames(X), "municipality")
  rhs    <- paste(x_cols, collapse = " + ")
  fixed  <- stats::as.formula(paste0("theta_direct ~ ", rhs))

  in_df <- data.frame(
    area_id      = direct$area_id,
    theta_direct = direct$theta_direct,
    vardir       = pmax(direct$var_direct_smoothed, .Machine$double.eps),
    X[match(direct$area_id, X$municipality), x_cols, drop = FALSE]
  )
  in_df <- in_df[stats::complete.cases(in_df), , drop = FALSE]

  # Centre and scale on the IN-SAMPLE covariates only, and reuse exactly those
  # constants for the out-of-sample block, so both matrices sit in one frame.
  #
  # Constant columns get scale 1 rather than 0, which avoids a division by
  # zero. They are still centred, so a constant column becomes all-zero; it is
  # collinear with the intercept with or without standardising.
  ctr <- stats::setNames(rep(0, length(x_cols)), x_cols)
  scl <- stats::setNames(rep(1, length(x_cols)), x_cols)
  if (isTRUE(standardize)) {
    for (cc in x_cols) {
      m <- mean(in_df[[cc]]); sdv <- stats::sd(in_df[[cc]])
      if (!is.finite(sdv) || sdv == 0) sdv <- 1
      ctr[cc] <- m; scl[cc] <- sdv
      in_df[[cc]] <- (in_df[[cc]] - m) / sdv
    }
  }

  fit <- tryCatch(
    saeRobust::rfh(formula = fixed, data = in_df,
                   samplingVar = "vardir", k = k),
    error = function(e) e
  )
  failed <- inherits(fit, "error")

  # oos_df is needed by both the failure and success branches.
  oos_df <- NULL
  if (!is.null(X_oos) && nrow(X_oos) > 0L) {
    oos_df <- data.frame(area_id = X_oos$municipality,
                         X_oos[, x_cols, drop = FALSE])
    if (isTRUE(standardize))
      for (cc in x_cols) oos_df[[cc]] <- (oos_df[[cc]] - ctr[cc]) / scl[cc]
  }

  if (failed) {
    in_ids  <- in_df$area_id
    oos_ids <- if (!is.null(oos_df)) oos_df$area_id else integer(0)
    all_ids   <- c(in_ids, oos_ids)
    area_set  <- c(rep("in_sample", length(in_ids)),
                   rep("out_of_sample", length(oos_ids)))
    est <- tibble::tibble(
      area_id   = all_ids,
      area_set  = area_set,
      theta_hat = NA_real_,
      mse       = NA_real_,
      ci_lo     = NA_real_,
      ci_hi     = NA_real_
    )
    diag <- tibble::tibble(
      metric = c("rfh_failed", "rfh_error_message"),
      value  = c(1, NA_real_),
      character_value = c(NA_character_, conditionMessage(fit))
    )
    meta <- list(method = "fh_robust",
                 runtime_sec = as.numeric(difftime(Sys.time(), t0, units = "secs")),
                 seed = NA_integer_,
                 n_areas = nrow(est),
                 n_in_sample = length(in_ids),
                 n_out_of_sample = length(oos_ids))
    return(structure(list(estimates = est, diagnostics = diag, meta = meta),
                     class = c("sae_fit", "sae_fit_fh_robust")))
  }

  # In-sample EBLUP + bootstrap MSE
  # Use the S3 generic predict() — predict.fitrfh is registered but not
  # directly exported from saeRobust >= 0.5.0.
  pred_in   <- predict(fit)
  # saeRobust 0.5.0: predict.fitrfh() returns the fields "re", "direct" and
  # "reblup"; there is no bias-corrected "reblupbc". The point estimate and the
  # MSE must target the same predictor, so the uncorrected REBLUP is used with
  # its matching pseudo-MSE (predType = "reblup", column "pseudo").
  theta_in  <- as.numeric(pred_in$reblup)
  # In saeRobust >= 0.5.0, mse() returns a data.frame; with predType = "reblup"
  # the pseudo-linearisation MSE is in the "pseudo" column.
  mse_obj   <- saeRobust::mse(fit, type = "pseudo",
                              predType = "reblup", B = B)
  mse_in    <- as.numeric(mse_obj$pseudo)
  # In saeRobust >= 0.5.0, fit$coefficients is a named numeric vector.
  beta_hat  <- as.numeric(fit$coefficients)
  sigma2_v  <- as.numeric(fit$variance)

  in_est <- tibble::tibble(
    area_id   = in_df$area_id,
    area_set  = "in_sample",
    theta_hat = theta_in,
    mse       = pmax(mse_in, .Machine$double.eps)
  )

  # OOS via synthetic prediction. OOS areas get no random-effect prediction,
  # so the leading-term synthetic-predictor variance (Rao & Molina 2015,
  # sec. 6.2.2) is
  #   MSE_i = sigma2_v + x_i' Var(beta_hat) x_i.
  # saeRobust exposes no vcov(), so Var(beta_hat) is approximated by the GLS
  # coefficient covariance (X' V^-1 X)^-1, V = diag(vardir + sigma2_v),
  # evaluated at the robust variance estimate.
  if (!is.null(oos_df)) {
    X_mat_oos <- cbind(1, as.matrix(oos_df[, x_cols, drop = FALSE]))
    theta_oos <- as.numeric(X_mat_oos %*% beta_hat)
    X_mat_in  <- cbind(1, as.matrix(in_df[, x_cols, drop = FALSE]))
    v_in      <- in_df$vardir + sigma2_v
    # NO back-transform here, deliberately. Under the reparametrisation the
    # scaling matrix cancels in the quadratic form x' Var(beta) x: beta and x
    # are transformed inversely, so the OOS variance is invariant and is
    # already on the right scale. Back-transforming beta_vcov to "match" the
    # reported betas would double-apply the scaling and inflate every OOS MSE.
    beta_vcov <- tryCatch(solve(crossprod(X_mat_in, X_mat_in / v_in)),
                          error = function(e) NULL)
    mse_oos <- if (!is.null(beta_vcov)) {
      sigma2_v + rowSums((X_mat_oos %*% beta_vcov) * X_mat_oos)
    } else {
      stats::median(in_est$mse)  # defensive fallback
    }
    oos_est <- tibble::tibble(
      area_id   = oos_df$area_id,
      area_set  = "out_of_sample",
      theta_hat = theta_oos,
      mse       = pmax(mse_oos, .Machine$double.eps)
    )
    est <- rbind(in_est, oos_est)
  } else {
    est <- in_est
  }
  se_v <- sqrt(est$mse)
  est$ci_lo <- est$theta_hat - 1.96 * se_v
  est$ci_hi <- est$theta_hat + 1.96 * se_v

  # Report coefficients on the ORIGINAL covariate scale, so the beta_*
  # diagnostics mean the same thing whether or not standardize was used:
  #   beta_j        = beta_std_j / scale_j
  #   beta_intercept = beta_std_0 - sum_j beta_std_j * centre_j / scale_j
  beta_report <- beta_hat
  if (isTRUE(standardize)) {
    # The back-transform needs one coefficient per covariate plus the
    # intercept; a mismatch stops here, where the message names it.
    stopifnot(length(beta_hat) == length(x_cols) + 1L)
    bs <- beta_hat[-1]
    # unname(): ctr/scl are named vectors, so the arithmetic would otherwise
    # attach covariate names to the diagnostics `value` column, which the
    # unstandardised path leaves unnamed.
    beta_report <- unname(c(beta_hat[1] - sum(bs * ctr[x_cols] / scl[x_cols]),
                            bs / scl[x_cols]))
  }
  beta_metric_names <- paste0("beta_", c("intercept", x_cols))
  diag <- tibble::tibble(
    metric          = c("sigma2_v", "k_huber", "B_bootstrap",
                        "rfh_failed", beta_metric_names),
    value           = c(sigma2_v, k, B, 0, beta_report),
    character_value = NA_character_
  )

  meta <- list(method = "fh_robust",
               runtime_sec = as.numeric(difftime(Sys.time(), t0, units = "secs")),
               seed = NA_integer_,
               n_areas = nrow(est),
               n_in_sample = sum(est$area_set == "in_sample"),
               n_out_of_sample = sum(est$area_set == "out_of_sample"))
  structure(list(estimates = est, diagnostics = diag, meta = meta),
            class = c("sae_fit", "sae_fit_fh_robust"))
}
