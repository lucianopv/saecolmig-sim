#' Fit a Fay-Herriot model with arcsin transformation via emdi2::fh
#'
#' Wraps emdi2::fh(transformation = "arcsin", backtransformation = "bc"),
#' the bias-corrected back-transformation of Hadam et al. (2020), as emdi2
#' implements it for the arcsin transformation. The variance-stabilising
#' transformation forces the model to respect the [0, 1] support of the
#' outcome.
#'
#' Pre-clips theta_direct to the open (0, 1) interval at
#' 1/(2 * n_eff) and 1 - 1/(2 * n_eff) — the standard arcsin preprocessor
#' for empirical proportions. On the arcsine scale emdi2 takes the sampling
#' variance from the effective sample size, 1 / (4 n_eff), so
#' `var_direct_smoothed` has no effect on the fit.
#'
#' Out-of-sample handling: emdi2::fh natively supports OOS via
#' combined_data with NA vardir on held-out rows. Out of sample it returns the
#' back-transformed synthetic prediction sin^2(x'beta), without the bias
#' correction of the in-sample areas, with the bootstrap MSE of its own.
#'
#' @param direct tibble (area_id, theta_direct, var_direct_smoothed, n_eff)
#' @param X tibble (municipality, x1, ..., xK) for in-sample areas
#' @param geo ignored
#' @param X_oos optional tibble (municipality, x1, ..., xK) for held-out areas
#' @param backtransformation "bc" (Hadam et al. 2020, default) or "naive" (sin^2)
#' @param B bootstrap reps for back-transformed MSE; default 50
#' @param seed integer seed of the parametric bootstrap; `NA` leaves emdi2's
#'   own default
#' @param ... ignored
#' @return list(estimates, diagnostics, meta) with class c("sae_fit", "sae_fit_fh_arcsin")
#' @export
fit_fh_arcsin <- function(direct, X, geo = NULL, X_oos = NULL,
                          backtransformation = "bc", B = 50L,
                          seed = NA_integer_, ...) {
  if (!requireNamespace("emdi2", quietly = TRUE)) {
    stop("Package 'emdi2' is required for fit_fh_arcsin().")
  }
  t0 <- Sys.time()

  x_cols <- setdiff(colnames(X), "municipality")
  rhs    <- paste(x_cols, collapse = " + ")
  fixed  <- stats::as.formula(paste0("theta_direct ~ ", rhs))

  # Effective sample size: prefer n_eff, fall back to n with a one-shot message
  if ("n_eff" %in% colnames(direct)) {
    eff_n <- pmax(direct$n_eff, 1)
  } else if ("n" %in% colnames(direct)) {
    message("fit_fh_arcsin: 'n_eff' not in direct; falling back to 'n'.")
    eff_n <- pmax(direct$n, 1)
  } else {
    stop("fit_fh_arcsin: direct must contain 'n_eff' or 'n'.")
  }

  # Pre-clip theta_direct to the open (0, 1) interval.
  # Note: n_eff == 1 forces lo = hi = 0.5, clipping any direct
  # estimate from such areas to exactly 0.5; at ECV sample sizes n_eff is far
  # above 1, so this is rare.
  lo <- 1 / (2 * eff_n)
  hi <- 1 - 1 / (2 * eff_n)
  theta_clipped <- pmin(pmax(direct$theta_direct, lo), hi)

  in_df <- data.frame(
    area_id      = direct$area_id,
    theta_direct = theta_clipped,
    vardir       = pmax(direct$var_direct_smoothed, .Machine$double.eps),
    n_eff        = eff_n,
    X[match(direct$area_id, X$municipality), x_cols, drop = FALSE]
  )
  combined <- in_df
  if (!is.null(X_oos) && nrow(X_oos) > 0L) {
    oos_df <- data.frame(
      area_id      = X_oos$municipality,
      theta_direct = NA_real_,
      vardir       = NA_real_,
      n_eff        = NA_real_,
      X_oos[, x_cols, drop = FALSE]
    )
    combined <- rbind(in_df, oos_df)
  }

  fh_args <- list(
    fixed              = fixed,
    vardir             = "vardir",
    combined_data      = combined,
    domains            = "area_id",
    method             = "reml",
    transformation     = "arcsin",
    backtransformation = backtransformation,
    eff_smpsize        = "n_eff",
    MSE                = TRUE,
    mse_type           = "boot",
    B                  = c(B, 0)
  )
  # emdi2::fh rejects a non-integer seed; pass it only when supplied,
  # otherwise emdi2 uses its own default (reproducible by construction).
  if (!is.na(seed)) fh_args$seed <- as.integer(seed)
  fit <- tryCatch(do.call(emdi2::fh, fh_args), error = function(e) e)
  failed <- inherits(fit, "error")

  if (failed) {
    in_ids  <- in_df$area_id
    oos_ids <- if (!is.null(X_oos) && nrow(X_oos) > 0L) X_oos$municipality else integer(0)
    all_ids  <- c(in_ids, oos_ids)
    area_set <- c(rep("in_sample", length(in_ids)),
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
      metric = c("fh_arcsin_failed", "fh_arcsin_error_message"),
      value  = c(1, NA_real_),
      character_value = c(NA_character_, conditionMessage(fit))
    )
    meta <- list(method = "fh_arcsin",
                 runtime_sec = as.numeric(difftime(Sys.time(), t0, units = "secs")),
                 seed = seed,
                 n_areas = nrow(est),
                 n_in_sample = length(in_ids),
                 n_out_of_sample = length(oos_ids))
    return(structure(list(estimates = est, diagnostics = diag, meta = meta),
                     class = c("sae_fit", "sae_fit_fh_arcsin")))
  }

  ind <- fit$ind
  mse <- fit$MSE
  ind <- ind[match(combined$area_id, ind$Domain), , drop = FALSE]
  mse <- mse[match(combined$area_id, mse$Domain), , drop = FALSE]

  theta_hat_v <- as.numeric(ind$FH)
  theta_hat_v <- pmin(pmax(theta_hat_v, 0), 1)  # belt-and-braces [0,1] guard
  mse_v       <- as.numeric(mse$FH)
  area_set    <- ifelse(is.na(combined$theta_direct),
                        "out_of_sample", "in_sample")
  # Guard: emdi2 returns a finite per-area bootstrap MSE for every domain,
  # out-of-sample ones included. Should a degenerate bootstrap draw leave an MSE
  # missing, it is filled with the in-sample median MSE so that no interval is
  # broken.
  if (any(is.na(mse_v))) {
    in_med <- stats::median(mse_v[area_set == "in_sample"], na.rm = TRUE)
    if (is.na(in_med)) in_med <- .Machine$double.eps
    mse_v[is.na(mse_v)] <- in_med
  }
  mse_v <- pmax(mse_v, .Machine$double.eps)
  se_v  <- sqrt(mse_v)
  ci_lo <- pmax(0, theta_hat_v - 1.96 * se_v)
  ci_hi <- pmin(1, theta_hat_v + 1.96 * se_v)
  # Prevent ci_lo == ci_hi when theta_hat = 1 and mse is near machine precision.
  ci_hi <- pmin(1 - .Machine$double.eps, ci_hi)

  est <- tibble::tibble(
    area_id   = combined$area_id,
    area_set  = area_set,
    theta_hat = theta_hat_v,
    mse       = mse_v,
    ci_lo     = ci_lo,
    ci_hi     = ci_hi
  )

  v <- fit$model$variance
  sigma2_v_val <- if (length(v) >= 1L) as.numeric(v[1]) else NA_real_
  diag <- tibble::tibble(
    metric          = c("sigma2_v", "transformation", "backtransformation",
                        "B_bootstrap", "fh_arcsin_failed"),
    value           = c(sigma2_v_val, NA_real_, NA_real_, B, 0),
    character_value = c(NA_character_, "arcsin", backtransformation,
                        NA_character_, NA_character_)
  )
  meta <- list(method = "fh_arcsin",
               runtime_sec = as.numeric(difftime(Sys.time(), t0, units = "secs")),
               seed = seed,
               n_areas = nrow(est),
               n_in_sample = sum(area_set == "in_sample"),
               n_out_of_sample = sum(area_set == "out_of_sample"))
  structure(list(estimates = est, diagnostics = diag, meta = meta),
            class = c("sae_fit", "sae_fit_fh_arcsin"))
}
