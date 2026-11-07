#' Compile (or re-use) the bayes_fh_heavy Stan model
#'
#' cmdstanr caches the compiled binary alongside the .stan source file
#' (model hash -> executable). Within one R process this helper also
#' memoises the compiled model object, saving the cmdstan_model() lookup.
#'
#' @keywords internal
.bayes_fh_heavy_model <- local({
  cached <- NULL
  function() {
    if (!is.null(cached)) return(cached)
    if (!requireNamespace("cmdstanr", quietly = TRUE)) {
      stop("Package 'cmdstanr' is required for fit_bayes_fh_heavy().")
    }
    stan_path <- system.file("stan", "bayes_fh_heavy.stan",
                             package = "saecol2026")
    if (!nzchar(stan_path)) {
      stop("Stan source 'bayes_fh_heavy.stan' not found in package.")
    }
    cached <<- cmdstanr::cmdstan_model(stan_file = stan_path,
                                        compile = TRUE, quiet = TRUE)
    cached
  }
})

#' Fit a Bayesian Fay-Herriot model with a half-Normal prior on the random-effect scale
#'
#' The method registered as `bayes_fh_heavy`. A Gaussian Fay-Herriot model in
#' Stan (`inst/stan/bayes_fh_heavy.stan`): y_d ~ N(theta_d, vardir_d),
#' theta_d = x_d' beta + v_d, v_d ~ N(0, sigma_v^2), with beta ~ N(0, 1) on
#' covariates standardised with the in-sample moments and a half-Normal(0,
#' `sigma_v_scale`) prior on sigma_v, sampled by cmdstanr.
#'
#' Out-of-sample areas are predicted in the generated-quantities block:
#' theta_oos = x_oos' beta + v_oos, with v_oos drawn from its prior
#' N(0, sigma_v^2) (an out-of-sample area has no direct estimate to inform its
#' random effect), so the interval carries both the regression uncertainty and
#' the random-effect variance.
#'
#' @param direct tibble (area_id, theta_direct, var_direct_smoothed, n_eff)
#' @param X tibble (municipality, x1, ..., xK) for in-sample areas
#' @param geo ignored
#' @param X_oos optional tibble (municipality, x1, ..., xK) for held-out areas
#' @param chains number of MCMC chains; default 4
#' @param iter_sampling post-warmup samples per chain; default 1000
#' @param iter_warmup warmup samples per chain; default 1000
#' @param sigma_v_scale scale of the half-Normal prior on sigma_v; default
#'   0.05
#' @param seed RNG seed; default NULL (cmdstanr generates one)
#' @param parallel_chains how many chains to run in parallel; default = chains
#' @param ... ignored
#' @return list(estimates, diagnostics, meta) with class
#'   c("sae_fit", "sae_fit_bayes_fh_heavy")
#' @export
fit_bayes_fh_heavy <- function(direct, X, geo = NULL, X_oos = NULL,
                               chains = 4L,
                               iter_sampling = 1000L,
                               iter_warmup = 1000L,
                               sigma_v_scale = 0.05,
                               seed = NULL,
                               parallel_chains = NULL, ...) {
  t0 <- Sys.time()
  fit_seed <- if (is.null(seed)) {
    sample.int(.Machine$integer.max, 1L)
  } else {
    as.integer(seed)
  }

  x_cols <- setdiff(colnames(X), "municipality")
  X_in_raw <- as.matrix(
    X[match(direct$area_id, X$municipality), x_cols, drop = FALSE])

  # Standardise covariates (centre + scale to unit sd) using in-sample
  # moments only, so OOS standardisation does not leak from held-out
  # data. Unstandardised count-scale columns (scale ~ 10^4) would blow up
  # the linear predictor.
  X_means <- colMeans(X_in_raw, na.rm = TRUE)
  X_sds   <- apply(X_in_raw, 2, stats::sd, na.rm = TRUE)
  X_sds[X_sds == 0 | !is.finite(X_sds)] <- 1
  X_in_std <- sweep(sweep(X_in_raw, 2, X_means, "-"), 2, X_sds, "/")

  X_in_mat <- cbind(1, X_in_std)
  colnames(X_in_mat) <- c("(Intercept)", x_cols)
  K <- ncol(X_in_mat)

  has_oos <- !is.null(X_oos) && nrow(X_oos) > 0L
  if (has_oos) {
    X_oos_raw <- as.matrix(X_oos[, x_cols, drop = FALSE])
    X_oos_std <- sweep(sweep(X_oos_raw, 2, X_means, "-"), 2, X_sds, "/")
    X_oos_mat <- cbind(1, X_oos_std)
    colnames(X_oos_mat) <- c("(Intercept)", x_cols)
    oos_ids <- X_oos$municipality
  } else {
    X_oos_mat <- matrix(0, nrow = 0, ncol = K)
    oos_ids <- integer(0)
  }

  in_ids <- direct$area_id
  data_list <- list(
    N_in     = length(in_ids),
    N_oos    = length(oos_ids),
    K        = K,
    X_in     = X_in_mat,
    X_oos    = X_oos_mat,
    y_direct = as.numeric(direct$theta_direct),
    vardir   = pmax(direct$var_direct_smoothed, .Machine$double.eps),
    sigma_v_scale = sigma_v_scale
  )

  model <- .bayes_fh_heavy_model()

  fit <- tryCatch(
    model$sample(
      data            = data_list,
      chains          = chains,
      iter_warmup     = iter_warmup,
      iter_sampling   = iter_sampling,
      parallel_chains = .null_or(parallel_chains, chains),
      seed            = fit_seed,
      init            = 0.1,           # tighter init around constrained-scale 0.1
      refresh         = 0,
      show_messages   = FALSE,
      show_exceptions = FALSE
    ),
    error = function(e) e
  )

  failed <- inherits(fit, "error")

  if (failed) {
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
      metric = c("bayes_fh_heavy_failed", "bayes_fh_heavy_error_message"),
      value  = c(1, NA_real_),
      character_value = c(NA_character_, conditionMessage(fit))
    )
    meta <- list(method = "bayes_fh_heavy",
                 runtime_sec = as.numeric(difftime(Sys.time(), t0, units = "secs")),
                 seed = fit_seed,
                 n_areas = nrow(est),
                 n_in_sample = length(in_ids),
                 n_out_of_sample = length(oos_ids))
    return(structure(list(estimates = est, diagnostics = diag, meta = meta),
                     class = c("sae_fit", "sae_fit_bayes_fh_heavy")))
  }

  # Posterior extraction for in-sample theta
  theta_in_draws  <- posterior::as_draws_matrix(fit$draws("theta_in"))
  theta_hat_in    <- as.numeric(colMeans(theta_in_draws))
  mse_in          <- as.numeric(apply(theta_in_draws, 2, stats::var))
  ci_in           <- apply(theta_in_draws, 2,
                            stats::quantile, probs = c(0.025, 0.975))

  in_est <- tibble::tibble(
    area_id   = in_ids,
    area_set  = "in_sample",
    theta_hat = theta_hat_in,
    mse       = pmax(mse_in, .Machine$double.eps),
    ci_lo     = as.numeric(ci_in[1, ]),
    ci_hi     = as.numeric(ci_in[2, ])
  )

  # OOS posterior (synthetic; uses generated quantities block)
  if (has_oos) {
    theta_oos_draws <- posterior::as_draws_matrix(fit$draws("theta_oos"))
    theta_hat_oos   <- as.numeric(colMeans(theta_oos_draws))
    mse_oos         <- as.numeric(apply(theta_oos_draws, 2, stats::var))
    ci_oos          <- apply(theta_oos_draws, 2,
                              stats::quantile, probs = c(0.025, 0.975))
    oos_est <- tibble::tibble(
      area_id   = oos_ids,
      area_set  = "out_of_sample",
      theta_hat = theta_hat_oos,
      mse       = pmax(mse_oos, .Machine$double.eps),
      ci_lo     = as.numeric(ci_oos[1, ]),
      ci_hi     = as.numeric(ci_oos[2, ])
    )
    est <- rbind(in_est, oos_est)
  } else {
    est <- in_est
  }

  # Diagnostics
  s   <- fit$summary(c("sigma_v", "beta"))
  ds  <- fit$diagnostic_summary(quiet = TRUE)
  sigma_v_mean <- s$mean[s$variable == "sigma_v"]
  max_rhat     <- max(s$rhat, na.rm = TRUE)
  min_ess_bulk <- min(s$ess_bulk, na.rm = TRUE)
  n_div <- sum(ds$num_divergent)
  total_iter <- chains * iter_sampling
  beta_rows  <- s[grepl("^beta\\[", s$variable), ]
  beta_metric_names <- paste0("beta_", c("intercept", x_cols))
  diag <- tibble::tibble(
    metric = c("sigma_v_mean", "max_rhat", "min_ess_bulk",
               "n_divergences", "frac_divergences",
               "bayes_fh_heavy_failed",
               beta_metric_names),
    value  = c(sigma_v_mean, max_rhat, min_ess_bulk,
               n_div, n_div / total_iter,
               0,
               beta_rows$mean),
    character_value = NA_character_
  )

  meta <- list(method = "bayes_fh_heavy",
               runtime_sec = as.numeric(difftime(Sys.time(), t0, units = "secs")),
               seed = fit_seed,
               n_areas = nrow(est),
               n_in_sample = sum(est$area_set == "in_sample"),
               n_out_of_sample = sum(est$area_set == "out_of_sample"))
  structure(list(estimates = est, diagnostics = diag, meta = meta),
            class = c("sae_fit", "sae_fit_bayes_fh_heavy"))
}
