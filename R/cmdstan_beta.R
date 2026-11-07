#' Capture tipsae's exact standata for the ExtBeta model
#'
#' tipsae assembles the Stan data (create_data + dummy_standata + X scaling +
#' island/spatial/temporal blocks) inside fit_sae, then calls
#' rstan::sampling(stanmodels$global, data = standata). Reimplementing that
#' assembly would be fragile, so a cheap 2-iteration fit_sae runs (the
#' data-prep is fast; only sampling is slow) and intercept the `data` argument
#' to rstan::sampling via a temporary namespace wrapper. Returns the exact
#' standata plus the in-/out-of-sample area_id ordering so posterior
#' theta[]/theta_oos[] map back.
#'
#' NOTE on interception strategy: rstan::sampling is an S4 generic whose body
#' delegates immediately to .local(object, data=..., ...). trace() on an S4
#' generic injects code before dispatch, where `...` does not exist, causing
#' "... used in a situation where it does not exist". Instead this function
#' temporarily replaces rstan::sampling in the rstan namespace with a plain wrapper that
#' captures the `data` argument and then delegates to the original S4 generic
#' (saved as orig_sampling in this function's environment). On exit the
#' namespace is restored unconditionally.
#'
#' Mirrors fit_bayes_fh_beta's ExtBeta inputs: smoothed phi via
#' tipsae::smoothing, type_disp = "neff", household size = n_dwell.
#'
#' The direct estimates reach tipsae unmodified: the ExtBeta likelihood has a
#' point mass at 0 and fits an exact zero directly.
#'
#' @param direct tibble (area_id, theta_direct, var_direct_smoothed, n_eff,
#'   n_dwell)
#' @param X tibble (municipality, x1..xK)
#' @param X_oos optional tibble (municipality, x1..xK)
#' @return list(standata, in_sample_ids, oos_ids); the ids are in tipsae's
#'   order (it sorts its data by area id), the order of `y`, `theta[]` and
#'   `theta_oos[]`
#' @keywords internal
capture_tipsae_standata <- function(direct, X, X_oos = NULL) {
  stopifnot(requireNamespace("tipsae", quietly = TRUE),
            requireNamespace("rstan", quietly = TRUE))
  x_cols <- setdiff(colnames(X), "municipality")
  rhs    <- paste(x_cols, collapse = " + ")
  fixed  <- stats::as.formula(paste0("theta_direct ~ ", rhs))

  in_df <- data.frame(
    area_id      = direct$area_id,
    theta_direct = direct$theta_direct,
    vardir       = pmax(direct$var_direct_smoothed, .Machine$double.eps),
    n_size       = as.numeric(direct$n_dwell),
    X[match(direct$area_id, X$municipality), x_cols, drop = FALSE])
  combined <- in_df
  if (!is.null(X_oos) && nrow(X_oos) > 0L) {
    oos_df <- data.frame(
      area_id = X_oos$municipality, theta_direct = NA_real_,
      vardir = NA_real_, n_size = NA_real_,
      X_oos[, x_cols, drop = FALSE])
    combined <- rbind(in_df, oos_df)
  }
  # smoothed dispersion (ExtBeta path), in-sample rows only
  in_only <- combined[!is.na(combined$theta_direct), ]
  sm <- tipsae::smoothing(direct_estimates = "theta_direct", data = in_only,
                          area_id = "area_id", raw_variance = "vardir",
                          areas_sample_sizes = "n_size", var_function = NULL,
                          method = "ols")
  combined$phi_smoothed <- NA_real_
  combined$phi_smoothed[!is.na(combined$theta_direct)] <- sm$phi

  # --- capture standata via namespace wrapper ----------------------------
  # tipsae calls rstan::sampling(stanmodels$global, data = standata, ...).
  # Temporarily replace rstan::sampling with a wrapper that saves `data`
  # into `captured` in this environment, then delegates to the real generic.
  captured     <- NULL
  orig_sampling <- rstan::sampling
  # Register the restore handler BEFORE patching so an interrupt between the
  # two cannot leave rstan::sampling monkey-patched.
  on.exit(
    suppressWarnings(utils::assignInNamespace("sampling", orig_sampling,
                                              ns = "rstan")),
    add = TRUE)
  env          <- environment()
  wrapper <- function(object, data = list(), ...) {
    env$captured <- data
    orig_sampling(object, data = data, ...)
  }
  suppressWarnings(utils::assignInNamespace("sampling", wrapper, ns = "rstan"))
  fit <- NULL
  invisible(utils::capture.output(fit <- suppressWarnings(suppressMessages(
    tipsae::fit_sae(formula_fixed = fixed, data = combined,
                    domains = "area_id", likelihood = "ExtBeta",
                    disp_direct = "phi_smoothed", type_disp = "neff",
                    household_size = "n_size",
                    chains = 1L, iter = 2L, warmup = 1L)))))
  if (is.null(captured)) stop("failed to capture tipsae standata")
  # tipsae sorts its data by the domain column, so y, theta[] and theta_oos[]
  # follow its order, not the caller's: the labels come from tipsae.
  dobj <- fit$data_obj
  list(standata      = captured,
       in_sample_ids = dobj$domains_names[!dobj$is_oos],
       oos_ids       = dobj$domains_names[dobj$is_oos])
}

#' Compile (or re-use) the vendored tipsae global model under cmdstanr
#'
#' Memoised: compiles global.stan once per R process and returns the cached
#' CmdStanModel thereafter. The model is compiled with stan_threads = TRUE;
#' cmdstanr >= 0.9 then requires threads_per_chain at sample time, so
#' fit_bayes_fh_beta_cmdstan() passes threads_per_chain = 1. Within-chain
#' threading is not used: the serial transformed-parameter loops dominate the
#' cost of this model.
#' @keywords internal
.tipsae_cmdstan_model <- local({
  cached <- NULL
  function() {
    if (!is.null(cached)) return(cached)
    if (!requireNamespace("cmdstanr", quietly = TRUE))
      stop("Package 'cmdstanr' is required for the cmdstan beta backend.")
    stan_path <- system.file("stan", "tipsae", "global.stan",
                             package = "saecol2026")
    inc <- system.file("stan", "tipsae", package = "saecol2026")
    if (stan_path == "") stop("vendored Stan model not found: global.stan")
    cached <<- cmdstanr::cmdstan_model(
      stan_file = stan_path, include_paths = inc,
      cpp_options = list(stan_threads = TRUE), compile = TRUE, quiet = TRUE)
    cached
  }
})

#' Fit the ExtBeta model via cmdstanr on tipsae's exact standata
#'
#' The `backend = "cmdstan"` path of [fit_bayes_fh_beta()]: tipsae's own Stan
#' model (vendored in `inst/stan/tipsae/`) on tipsae's own data
#' ([capture_tipsae_standata()]), sampled by CmdStan with tipsae's sampler
#' controls (`adapt_delta = 0.95`, `max_treedepth = 10`) and the start of
#' [extbeta_init_intercept()]. The default backend of [fit_bayes_fh_beta()] is
#' rstan through tipsae. Out of sample the estimate is tipsae's `theta_oos`.
#'
#' @param direct,X,X_oos as [capture_tipsae_standata()]
#' @param chains,iter chains and total iterations per chain (half warm-up)
#' @param adapt_delta,max_treedepth sampler controls
#' @param seed integer; `NA` leaves the sampler unseeded
#' @keywords internal
fit_bayes_fh_beta_cmdstan <- function(direct, X, X_oos = NULL,
                                      chains = 4L, iter = 2000L,
                                      adapt_delta = 0.95,
                                      max_treedepth = 10L,
                                      seed = NA_integer_) {
  t0  <- Sys.time()
  cap <- capture_tipsae_standata(direct, X, X_oos = X_oos)
  init <- rep(list(extbeta_init_intercept(cap$standata)), chains)
  # cmdstanr lists the model's zero-size parameters as "missing" from a
  # partial init; they have no values to set, so the message is suppressed.
  old_opt <- options(cmdstanr_warn_inits = FALSE)
  on.exit(options(old_opt), add = TRUE)
  mod <- .tipsae_cmdstan_model()
  # tipsae's Stan sampler controls: adapt_delta = 0.95 and max_treedepth = 10,
  # and the start of extbeta_init_intercept(). cmdstanr's defaults (random
  # inits on [-2, 2], adapt_delta = 0.8) do not mix on the near-boundary
  # ExtBeta posterior (the sampler hits the maximum tree depth).
  sm  <- mod$sample(
    data = cap$standata, chains = chains, parallel_chains = chains,
    threads_per_chain = 1L,
    iter_warmup = iter %/% 2L, iter_sampling = iter %/% 2L,
    adapt_delta = adapt_delta, max_treedepth = max_treedepth, init = init,
    seed = if (is.na(seed)) NULL else as.integer(seed),
    refresh = 0, show_messages = FALSE)

  pull <- function(var, ids) {
    if (!length(ids)) return(NULL)
    s <- sm$summary(var, "mean", "sd",
                    q = ~quantile(., probs = c(0.025, 0.975)))
    data.frame(area_id = ids, theta_hat = s$mean,
               mse = pmax(s$sd^2, .Machine$double.eps),
               ci_lo = s$`2.5%`, ci_hi = s$`97.5%`,
               stringsAsFactors = FALSE)
  }
  est_in  <- pull("theta", cap$in_sample_ids)
  est_in$area_set <- "in_sample"
  est <- est_in
  if (length(cap$oos_ids)) {
    est_oos <- pull("theta_oos", cap$oos_ids)
    est_oos$area_set <- "out_of_sample"
    est <- rbind(est_in, est_oos)
  }
  div <- tryCatch(sum(sm$diagnostic_summary()$num_divergent),
                  error = function(e) NA_real_)
  rh <- tryCatch(max(sm$summary(NULL, "rhat")$rhat, na.rm = TRUE),
                 error = function(e) NA_real_)
  diag <- tibble::tibble(metric = c("n_divergent", "max_rhat"),
                         value  = c(as.numeric(div), rh))
  meta <- list(method = "bayes_fh_beta", likelihood = "ExtBeta",
               backend = "cmdstan",
               runtime_sec = as.numeric(difftime(Sys.time(), t0, units="secs")),
               n_areas = nrow(est), n_in_sample = nrow(est_in),
               n_out_of_sample = nrow(est) - nrow(est_in))
  structure(list(estimates = tibble::as_tibble(est), diagnostics = diag,
                 meta = meta),
            class = c("sae_fit", "sae_fit_bayes_fh_beta"))
}

#' Initial values for ExtBeta
#'
#' tipsae's `init = 0` starts every mu at 0.5, so the zero mass p0 of an area
#' with an exact zero starts at 0.5^m. For an area with many sampled households
#' (m of about 1,075 or more) that underflows to 0 in double precision, giving
#' log(p0) = -Inf and a rejected start in every chain; sampling then does not
#' start and tipsae raises no R error. This returns the constrained-scale
#' equivalent of `init = 0` for every active parameter (sigma_v = exp(0) = 1,
#' lambda_star_EB = inv_logit(0) = 0.5) except the intercept, which starts at
#' logit(mean y). Every active parameter is listed so none falls back to a
#' random start, which does not mix on this posterior. Both backends of
#' [fit_bayes_fh_beta()] start ExtBeta here.
#'
#' @param standata tipsae's standata (from [capture_tipsae_standata()])
#' @return a named list for one chain
#' @keywords internal
extbeta_init_intercept <- function(standata) {
  sd <- standata
  ok <- identical(as.integer(sd$likelihood), 3L) &&
    identical(as.integer(sd$intercept), 1L) &&
    identical(as.integer(sd$prior_coeff), 0L) &&
    identical(as.integer(sd$prior_reff), 0L) &&
    identical(as.integer(sd$spatial_err), 0L) &&
    identical(as.integer(sd$temporal_err), 0L)
  if (!ok) {
    stop("extbeta_init_intercept: only the ExtBeta structure of ",
         "fit_bayes_fh_beta() is supported (intercept, normal coefficients and ",
         "random effect, no spatial or temporal error)")
  }
  ybar <- mean(sd$y)
  if (!(ybar > 0 && ybar < 1)) stop("extbeta_init_intercept: mean(y) not in (0, 1)")
  list(beta0 = array(stats::qlogis(ybar), dim = 1L),
       z_beta = array(0, dim = sd$P),
       v_raw = rep(0, sd$M_is),       # length_v = M_is without a temporal error
       sigma_v = array(1, dim = 1L),
       lambda_star_EB = array(0.5, dim = 1L))
}
