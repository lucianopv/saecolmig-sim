#' Fit a Bayesian Fay-Herriot model (Beta or Extended Beta; tipsae backend)
#'
#' Wraps tipsae::fit_sae. Two paths:
#'
#' - `likelihood = "beta"` (default, faster): variance-based dispersion via
#'   `type_disp = "var"` and `disp_direct = vardir`. Simpler model.
#' - `likelihood = "ExtBeta"`: De Nicolo et al. (2024) Extended Beta with
#'   smoothed phi via `tipsae::smoothing()` + `type_disp = "neff"` + the area
#'   sample size m_d (tipsae's `household_size`, required by ExtBeta), here the
#'   number of sampled households `n_dwell`.
#'
#' The direct estimates reach tipsae unmodified, for every likelihood and on
#' both backends. ExtBeta has a point mass at 0 and fits an exact zero
#' directly; the plain beta likelihoods take none, and tipsae stops with its own
#' message. Under ExtBeta both backends start the chains at the values of
#' [extbeta_init_intercept()] (every parameter at tipsae's `init = 0` except the
#' intercept, at the logit of the mean direct estimate): from tipsae's start the
#' zero mass of an area with an exact zero is 0.5^m, which is 0 in double
#' precision for about 1,075 sampled households or more, and the chains cannot
#' start. On rstan the start goes to the sampler call tipsae makes (its `init`
#' argument takes one list per chain only when `chains` is 1). The other
#' likelihoods, and ExtBeta with another prior for the coefficients on rstan,
#' use tipsae's default start.
#'
#' tipsae natively detects out-of-sample areas (rows with an NA direct
#' estimate) and saves theta_oos in the Stan output; that is the out-of-sample
#' estimate (for ExtBeta it is mu, not the in-sample estimand
#' theta = mu (1 - p0 - p1) + p1).
#' Posterior summaries map onto the method contract ([sae_fit_contract()]):
#' theta_hat = posterior mean, mse = posterior variance, ci_lo / ci_hi = 2.5 /
#' 97.5 posterior quantiles.
#'
#' @param direct tibble (area_id, theta_direct, var_direct_smoothed, n_eff,
#'   n_dwell, ...). ExtBeta takes n_dwell as the area sample size m_d.
#' @param X tibble (municipality, x1, ..., xK)
#' @param geo unused
#' @param X_oos optional tibble (municipality, x1, ..., xK) for held-out areas
#' @param chains integer; default 4
#' @param iter integer total iterations per chain; default 2000
#' @param warmup integer warmup iterations per chain; default half of `iter`
#' @param likelihood "beta" (default), "ExtBeta", "Infbeta01", "Infbeta0",
#'   "Infbeta1" or "flexbeta". "Infbeta01" does not start from tipsae's default
#'   start (tipsae warns of this); "flexbeta" does not save `theta_oos` and
#'   stops when `X_oos` is given.
#' @param prior_coeff "normal" (default) or "HorseShoe"
#' @param p0_HorseShoe integer; expected number of non-zero coefficients
#'   under HorseShoe (only used when prior_coeff = "HorseShoe")
#' @param adapt_delta numeric in (0,1); forwarded to tipsae::fit_sae (which
#'   passes it to rstan::sampling via HMC control). NULL uses the tipsae/rstan
#'   default (0.95 in tipsae; 0.8 in raw rstan).
#' @param backend "rstan" (default) runs the tipsae path; "cmdstan" runs the
#'   vendored tipsae global.stan via cmdstanr on tipsae's exact standata
#'   (always the ExtBeta model). See
#'   [fit_bayes_fh_beta_cmdstan()].
#' @param seed integer
#' @param ... ignored
#' @section Arguments by backend:
#' * `backend = "rstan"` honours `direct`, `X`, `X_oos`, `chains`, `iter`,
#'   `warmup`, `likelihood`, `prior_coeff`, `p0_HorseShoe`, `adapt_delta` and
#'   `seed`. The ExtBeta start is defined for `prior_coeff = "normal"`; with
#'   another prior the chains start at tipsae's default.
#' * `backend = "cmdstan"` honours `direct`, `X`, `X_oos`, `chains`, `iter` and
#'   `seed`. It always fits the ExtBeta model with a normal prior on the
#'   coefficients, warms up for half of `iter` iterations and samples at
#'   `adapt_delta = 0.95`. A call that asks for something else stops, naming
#'   the argument: a `likelihood` other than `"ExtBeta"` (it must be given),
#'   `prior_coeff` other than `"normal"`, a `p0_HorseShoe`, an `adapt_delta`
#'   other than 0.95, or a `warmup` other than half of `iter`.
#'   `adapt_delta = NULL` is not a request.
#'
#' The diagnostics `max_rhat` and `n_divergent` are those of the fit. `max_rhat`
#' is the largest rank-normalised split R-hat of the `posterior` package over
#' the variables the backend keeps. Under ExtBeta that is, on CmdStan, every
#' variable of the model and its generated quantities, and on rstan the
#' variables tipsae keeps (`beta0`, `beta`, `v`, `sigma_v`, `mu`, `p0`, `p1`,
#' `lambda_EB`, `theta`, `theta_oos`, `y_rep` and `log_lik`, with `lp__`);
#' other likelihoods keep other variables on rstan.
#'
#' tipsae's horseshoe prior (`prior_coeff = "HorseShoe"`) does not take exact
#' zeros in the direct estimates: rstan then returns no draws.
#'
#' `geo` and `...` are ignored by both backends.
#' @return list(estimates, diagnostics, meta) with class
#'   c("sae_fit", "sae_fit_bayes_fh_beta")
#' @export
fit_bayes_fh_beta <- function(direct, X, geo = NULL, X_oos = NULL,
                              chains = 4L, iter = 2000L, warmup = iter %/% 2L,
                              likelihood = "beta",
                              prior_coeff = "normal",
                              p0_HorseShoe = NULL,
                              adapt_delta = NULL,
                              backend = c("rstan", "cmdstan"),
                              seed = NA_integer_, ...) {
  if (!requireNamespace("tipsae", quietly = TRUE)) {
    stop("Package 'tipsae' is required for fit_bayes_fh_beta().")
  }
  if (!requireNamespace("rstan", quietly = TRUE)) {
    stop("Package 'rstan' is required for fit_bayes_fh_beta().")
  }
  backend <- match.arg(backend)
  if (backend == "cmdstan") {
    # This backend fits ExtBeta with a normal prior, warms up for iter %/% 2
    # iterations and samples at adapt_delta = 0.95. An argument that asks for
    # something else would be dropped without a trace, so the call stops.
    # `adapt_delta` at the value applied, and a `warmup` equal to the number
    # applied, are not requests.
    unused <- c(
      if (!identical(likelihood, "ExtBeta")) "likelihood",
      if (!identical(prior_coeff, "normal")) "prior_coeff",
      if (!is.null(p0_HorseShoe)) "p0_HorseShoe",
      if (!is.null(adapt_delta) && !isTRUE(all.equal(adapt_delta, 0.95))) "adapt_delta",
      if (!isTRUE(all.equal(warmup, iter %/% 2L))) "warmup")
    if (length(unused)) {
      stop("backend = 'cmdstan' does not use ", paste0("`", unused, "`", collapse = ", "),
           ": it fits ExtBeta with a normal prior, warms up for iter %/% 2 ",
           "iterations and samples at adapt_delta = 0.95. Use backend = 'rstan' ",
           "for these arguments",
           if ("likelihood" %in% unused) {
             ", or pass likelihood = \"ExtBeta\" to fit ExtBeta on CmdStan"
           },
           ".", call. = FALSE)
    }
    return(fit_bayes_fh_beta_cmdstan(direct, X, X_oos = X_oos,
                                     chains = chains, iter = iter,
                                     seed = seed))
  }
  if (!requireNamespace("posterior", quietly = TRUE)) {
    stop("Package 'posterior' is required for fit_bayes_fh_beta(backend = 'rstan').")
  }
  t0 <- Sys.time()

  x_cols <- setdiff(colnames(X), "municipality")
  rhs    <- paste(x_cols, collapse = " + ")
  fixed  <- stats::as.formula(paste0("theta_direct ~ ", rhs))

  needs_n <- likelihood %in% c("ExtBeta", "Infbeta0", "Infbeta1", "Infbeta01")
  if (needs_n && !"n_dwell" %in% colnames(direct)) {
    stop(sprintf(
      "likelihood = '%s' requires direct$n_dwell (sampled households per area).",
      likelihood
    ))
  }

  # The direct estimates go to tipsae as they are.
  in_df <- data.frame(
    area_id      = direct$area_id,
    theta_direct = direct$theta_direct,
    vardir       = pmax(direct$var_direct_smoothed, .Machine$double.eps),
    n_size       = if (needs_n) as.numeric(direct$n_dwell) else NA_real_,
    X[match(direct$area_id, X$municipality), x_cols, drop = FALSE]
  )
  combined <- in_df
  if (!is.null(X_oos) && nrow(X_oos) > 0L) {
    oos_df <- data.frame(
      area_id      = X_oos$municipality,
      theta_direct = NA_real_,
      vardir       = NA_real_,
      n_size       = NA_real_,
      X_oos[, x_cols, drop = FALSE]
    )
    combined <- rbind(in_df, oos_df)
  }

  fit_seed <- if (is.na(seed)) sample.int(.Machine$integer.max, 1L) else seed

  # tipsae::fit_sae forwards `...` to rstan::sampling, which accepts `cores`.
  # Running 4 chains in parallel cuts wall by ~2.5-3x vs the sequential default.
  fit_args <- list(
    formula_fixed = fixed,
    data          = combined,
    domains       = "area_id",
    likelihood    = likelihood,
    prior_coeff   = prior_coeff,
    chains        = chains,
    iter          = iter,
    warmup        = warmup,
    cores         = chains,
    seed          = fit_seed,
    refresh       = 0
  )
  if (!is.null(p0_HorseShoe)) fit_args$p0_HorseShoe <- p0_HorseShoe
  if (!is.null(adapt_delta)) fit_args$adapt_delta <- adapt_delta

  if (likelihood == "ExtBeta") {
    # Smooth phi via tipsae::smoothing on in-sample rows; OOS rows have an NA
    # direct estimate, which tipsae's create_data() flags as is_oos.
    in_only <- combined[!is.na(combined$theta_direct), ]
    sm <- tipsae::smoothing(
      direct_estimates   = "theta_direct",
      data               = in_only,
      area_id            = "area_id",
      raw_variance       = "vardir",
      areas_sample_sizes = "n_size",
      var_function       = NULL,
      method             = "ols"
    )
    combined$phi_smoothed <- NA_real_
    combined$phi_smoothed[!is.na(combined$theta_direct)] <- sm$phi
    fit_args$data           <- combined
    fit_args$disp_direct    <- "phi_smoothed"
    fit_args$type_disp      <- "neff"
    fit_args$household_size <- "n_size"
  } else if (likelihood %in% c("Infbeta0", "Infbeta1", "Infbeta01")) {
    # Inflated-beta likelihoods: tipsae requires type_disp = "neff".
    # n_size is already populated (checked above via needs_n).
    fit_args$disp_direct <- "n_size"
    fit_args$type_disp   <- "neff"
  } else {
    # 'beta' / 'flexbeta': variance-based dispersion.
    fit_args$disp_direct <- "vardir"
    fit_args$type_disp   <- "var"
  }

  fit <- if (likelihood == "ExtBeta" && identical(prior_coeff, "normal")) {
    .with_extbeta_start(do.call(tipsae::fit_sae, fit_args))
  } else {
    do.call(tipsae::fit_sae, fit_args)
  }
  # tipsae raises no error when rstan cannot sample (the chains cannot start,
  # or an argument is refused); rstan then returns an empty fit, and the
  # summaries below would fail obscurely.
  if (fit$stanfit@mode != 0L) {
    stop("rstan returned no draws; see its messages above")
  }

  post_in <- rstan::summary(fit$stanfit, pars = "theta")$summary

  # tipsae::fit_sae() unconditionally re-sorts its `data` argument by the
  # domains column before fitting (whenever domains is non-NULL and
  # temporal_error is off, which is always true here) -- so theta[]'s row
  # order follows that internal sort, not `combined`'s original row order.
  # Must pull labels from fit$data_obj$domains_names (post-sort), not
  # combined$area_id (pre-sort), or in-sample/OOS labels silently
  # misalign whenever the two blocks interleave under a global area_id
  # sort (i.e. whenever real OOS rows are present).
  est_in <- tibble::tibble(
    area_id   = fit$data_obj$domains_names[!fit$data_obj$is_oos],
    area_set  = "in_sample",
    theta_hat = unname(post_in[, "mean"]),
    mse       = pmax(unname(post_in[, "sd"])^2, .Machine$double.eps),
    ci_lo     = unname(post_in[, "2.5%"]),
    ci_hi     = unname(post_in[, "97.5%"])
  )

  est <- est_in
  if (any(fit$data_obj$is_oos)) {
    post_oos <- rstan::summary(fit$stanfit, pars = "theta_oos")$summary
    est_oos <- tibble::tibble(
      area_id   = fit$data_obj$domains_names[fit$data_obj$is_oos],
      area_set  = "out_of_sample",
      theta_hat = unname(post_oos[, "mean"]),
      mse       = pmax(unname(post_oos[, "sd"])^2, .Machine$double.eps),
      ci_lo     = unname(post_oos[, "2.5%"]),
      ci_hi     = unname(post_oos[, "97.5%"])
    )
    est <- rbind(est_in, est_oos)
  }

  # R-hat and bulk ESS: the rank-normalised statistics of the posterior
  # package, as on the CmdStan path, over the variables the fit keeps (R-hat)
  # and over theta (ESS).
  # (posterior warns about every capped ESS of a very short chain.)
  draw_stats <- withCallingHandlers(
    posterior::summarise_draws(posterior::as_draws_array(fit$stanfit), "rhat", "ess_bulk"),
    warning = function(w) {
      if (grepl("ESS has been capped", conditionMessage(w), fixed = TRUE)) {
        invokeRestart("muffleWarning")
      }
    })
  rhat_vec <- draw_stats$rhat
  ess_vec  <- draw_stats$ess_bulk[grepl("^theta\\[", draw_stats$variable)]
  n_div <- tryCatch(rstan::get_num_divergent(fit$stanfit),
                    error = function(e) NA_real_)
  diag <- tibble::tibble(
    metric = c("frac_rhat_above_1.05", "min_ess_bulk",
               "max_rhat", "median_ess_bulk", "n_divergent"),
    value  = c(mean(rhat_vec > 1.05, na.rm = TRUE),
               suppressWarnings(min(ess_vec, na.rm = TRUE)),
               suppressWarnings(max(rhat_vec, na.rm = TRUE)),
               stats::median(ess_vec, na.rm = TRUE),
               as.numeric(n_div))
  )
  meta <- list(method = "bayes_fh_beta",
               likelihood = likelihood,
               runtime_sec = as.numeric(difftime(Sys.time(), t0, units = "secs")),
               seed = fit_seed,
               adapt_delta = adapt_delta,
               n_areas = nrow(est),
               n_in_sample = sum(est$area_set == "in_sample"),
               n_out_of_sample = sum(est$area_set == "out_of_sample"))
  structure(list(estimates = est, diagnostics = diag, meta = meta),
            class = c("sae_fit", "sae_fit_bayes_fh_beta"))
}

#' Run a tipsae ExtBeta fit with the start of [extbeta_init_intercept()]
#'
#' `tipsae::fit_sae()` passes its `init` argument to `rstan::sampling()`, but
#' tests it first with `init == "0"` inside `if`, so a list with one element
#' per chain (several chains) stops with "the condition has length > 1". The
#' start is therefore set on the call `fit_sae()` makes to the sampler: for
#' the duration of `expr`, `rstan::sampling` is replaced by a wrapper that,
#' when it is called with tipsae's default `init = "0"`, computes the start
#' from the Stan data it receives and passes it on, one list per chain; a call
#' with any other `init` (rstan runs the chains of a parallel fit as one call
#' each, with their own start) passes through. The original is restored on
#' exit.
#'
#' @param expr expression that calls `tipsae::fit_sae()` for ExtBeta
#' @keywords internal
.with_extbeta_start <- function(expr) {
  orig_sampling <- rstan::sampling
  # Register the restore handler BEFORE patching so an interrupt between the
  # two cannot leave rstan::sampling replaced.
  on.exit(
    suppressWarnings(utils::assignInNamespace("sampling", orig_sampling,
                                              ns = "rstan")),
    add = TRUE)
  wrapper <- function(object, data = list(), ..., init = "random") {
    if (identical(init, "0")) {
      # rstan's default chain count when the call names none.
      init <- rep(list(extbeta_init_intercept(data)), .null_or(list(...)$chains, 4L))
    }
    orig_sampling(object, data = data, ..., init = init)
  }
  suppressWarnings(utils::assignInNamespace("sampling", wrapper, ns = "rstan"))
  force(expr)
}
