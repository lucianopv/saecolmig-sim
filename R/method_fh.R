# Extract the random-effect variance sigma^2_v from an emdi2::fh model's
# `$variance`. Non-spatial FH returns a scalar; spatial FH returns a named
# list (correlation, variance, convergence), so position-based access would
# silently pick up rho. Always resolve by name when available.
.fh_sigma2_v <- function(v) {
  if (!is.null(names(v)) && "variance" %in% names(v)) {
    return(as.numeric(v[["variance"]]))
  }
  as.numeric(v[1])
}

#' Fit a Fay-Herriot model via emdi2::fh
#'
#' emdi2::fh handles out-of-sample areas natively through `combined_data`:
#' held-out rows carry an NA direct estimate and NA vardir and are predicted
#' from their covariates. The arcsine-transformed variant is [fit_fh_arcsin()].
#'
#' @param direct tibble (area_id, theta_direct, var_direct_smoothed, ...)
#' @param X tibble (municipality, x1, ..., xK) for in-sample areas
#' @param geo optional list; geo$W passed as corMatrix when correlation="spatial"
#' @param X_oos optional tibble (municipality, x1, ..., xK) for held-out areas
#' @param transformation "no" (default) or "arcsin"
#' @param correlation "no" (default) or "spatial"
#' @param ... ignored
#' @return list(estimates, diagnostics, meta) with class
#'   c("sae_fit", "sae_fit_fh") (or "sae_fit_fh_spatial" via the wrapper)
#' @export
fit_fh <- function(direct, X, geo = NULL, X_oos = NULL,
                   transformation = "no",
                   correlation = "no", ...) {
  if (!requireNamespace("emdi2", quietly = TRUE)) {
    stop("Package 'emdi2' is required for fit_fh().")
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
  combined <- in_df
  if (!is.null(X_oos) && nrow(X_oos) > 0L) {
    oos_df <- data.frame(
      area_id      = X_oos$municipality,
      theta_direct = NA_real_,
      vardir       = NA_real_,
      X_oos[, x_cols, drop = FALSE]
    )
    combined <- rbind(in_df, oos_df)
  }

  fit_args <- list(
    fixed          = fixed,
    vardir         = "vardir",
    combined_data  = combined,
    domains        = "area_id",
    method         = "reml",
    transformation = transformation,
    MSE            = TRUE
  )
  if (correlation == "spatial") {
    if (is.null(geo) || is.null(geo$W)) {
      stop("fit_fh(correlation='spatial') requires geo$W.")
    }
    # emdi requires corMatrix to span ONLY the in-sample areas; OOS rows
    # are predicted via the synthetic X*beta component (no spatial
    # borrowing for them).
    in_ids <- as.character(in_df$area_id)
    W_full <- as.matrix(geo$W)
    W_in   <- W_full[in_ids, in_ids, drop = FALSE]
    fit_args$correlation <- "spatial"
    fit_args$corMatrix   <- W_in
    # 'analytical' MSE for spatial FH - the parametric bootstrap variants
    # ('spatialparboot', 'spatialnonparboot') periodically hit a singular
    # (I - rho*W')*(I - rho*W) on small samples or degenerate W structures.
    # Analytical is faster and stable on samples of several hundred areas.
    fit_args$mse_type    <- "analytical"
  }
  fit <- do.call(emdi2::fh, fit_args)

  ind <- fit$ind                         # cols: Domain, Direct, FH, Out
  mse <- fit$MSE                         # cols: Domain, Direct, FH, Out
  ind <- ind[match(combined$area_id, ind$Domain), ]
  mse <- mse[match(combined$area_id, mse$Domain), ]

  theta_hat_v <- as.numeric(ind$FH)
  mse_v       <- as.numeric(mse$FH)
  area_set    <- ifelse(is.na(combined$theta_direct),
                        "out_of_sample", "in_sample")
  # emdi's analytical MSE covers in-sample domains only; under spatial FH the
  # out-of-sample rows come back NA. OOS areas receive no random-effect
  # prediction (synthetic X*beta only), so their MSE is the leading-term
  # synthetic-predictor variance (Rao & Molina 2015, sec. 6.2.2):
  #   MSE_i = sigma2_v + x_i' Var(beta_hat) x_i,
  # with Var(beta_hat) taken from emdi's `beta_vcov`.
  na_mse <- is.na(mse_v)
  if (any(na_mse)) {
    bvc <- fit$model$beta_vcov
    s2v <- .fh_sigma2_v(fit$model$variance)
    if (!is.null(bvc) && is.finite(s2v)) {
      Xna  <- cbind(1, as.matrix(combined[na_mse, x_cols, drop = FALSE]))
      quad <- rowSums((Xna %*% bvc) * Xna)
      mse_v[na_mse] <- pmax(s2v + quad, .Machine$double.eps)
    } else {
      # Defensive fallback if emdi did not expose beta_vcov / variance.
      in_med <- stats::median(mse_v[!na_mse], na.rm = TRUE)
      mse_v[na_mse] <- if (is.na(in_med)) .Machine$double.eps else in_med
    }
  }
  mse_v <- pmax(mse_v, .Machine$double.eps)
  se_v  <- sqrt(mse_v)

  est <- tibble::tibble(
    area_id   = combined$area_id,
    area_set  = area_set,
    theta_hat = theta_hat_v,
    mse       = mse_v,
    ci_lo     = theta_hat_v - 1.96 * se_v,
    ci_hi     = theta_hat_v + 1.96 * se_v
  )

  # Spatial FH returns variance as a named list (correlation, variance,
  # convergence); plain FH returns a scalar. Resolve by name: position-based
  # access would return rho in place of sigma2_v for the spatial case.
  v <- fit$model$variance
  sigma2_v_val <- .fh_sigma2_v(v)
  rho_val      <- if (!is.null(names(v)) && "correlation" %in% names(v)) {
    as.numeric(v[["correlation"]])
  } else NA_real_
  diag <- tibble::tibble(
    metric          = c("sigma2_v", "rho", "transformation", "correlation"),
    value           = c(sigma2_v_val, rho_val, NA_real_, NA_real_),
    character_value = c(NA_character_, NA_character_, transformation, correlation)
  )

  meta <- list(method = if (correlation == "spatial") "fh_spatial" else "fh",
               runtime_sec = as.numeric(difftime(Sys.time(), t0, units = "secs")),
               seed = NA_integer_,
               n_areas = nrow(est),
               n_in_sample = sum(area_set == "in_sample"),
               n_out_of_sample = sum(area_set == "out_of_sample"))
  cls <- if (correlation == "spatial") "sae_fit_fh_spatial" else "sae_fit_fh"
  structure(list(estimates = est, diagnostics = diag, meta = meta),
            class = c("sae_fit", cls))
}
