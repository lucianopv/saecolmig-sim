#' Sub-area mixed-effects gradient boosting via povmap::megb
#'
#' Companion to `fit_xgb_sub()` on identical inputs, so the two differ in
#' estimator only.
#'
#' `povmap::megb()` has no `sub_domains` argument, so the two-level structure
#' is expressed differently: the sub-area rows are the population cells,
#' `domains = "municipality"` names the target level, and `pop_weights` carries
#' the census counts that weight cells into a municipality estimate. The model
#' sees the same data as `fit_xgb_sub()`; only the interface for declaring the
#' hierarchy differs.
#'
#' Uncertainty is megb's parametric bootstrap (`mse = TRUE`, `mse_type =
#' "var"`): the variance of the centred bootstrap prediction errors and a basic
#' bootstrap interval. The estimator is cross-fitted: folds hold out whole
#' municipalities (`cv_folds = "domain"`, `cv_nfold = 10`), sampled
#' municipalities are predicted by the fold model that held them out and
#' unsampled ones by the fold-model average, and the bootstrap re-estimates the
#' random effects on out-of-fold residuals. These settings are povmap's
#' defaults but povmap reads them from `getOption()`, so they are passed
#' explicitly: a global option must not change the estimator. The cross-fitting
#' needs more drawn municipalities than folds.
#'
#' @param direct tibble of direct estimates for this outcome; used only to
#'   identify drawn municipalities, never fed to the model
#' @param smp_sub collapsed sub-area sample for this outcome
#' @param pop_sub sub-area covariate matrix over the whole universe
#' @param outcome outcome name from the registry
#' @param B bootstrap replications; default 100
#' @param seed integer seed passed to povmap; `NA` passes none
#' @param mse_type povmap's `mse_type`: "var" (default; the variance of the
#'   centred bootstrap prediction errors and a basic interval, which leave each
#'   domain's bootstrap bias out) or "mse" (the prediction MSE including that
#'   bias, with a normal-approximation interval)
#' @param bootstrap_cores cores for megb's bootstrap; default 0 (serial)
#' @param params optional named list forwarded as `gradient_params`
#' @param level `"subarea"` (default) or `"municipality"`; the latter
#'   aggregates the predictors to municipality and fits the same estimator
#'   there, which separates resolution from method.
#' @param weight `"census"` (default) or `"worldpop"`; the population weight,
#'   see `.subarea_inputs()`.
#' @param ... forwarded to `povmap::megb()`
#' @return list(estimates, diagnostics, meta), class c("sae_fit", "sae_fit_megb_sub")
#' @export
fit_megb_sub <- function(direct, smp_sub, pop_sub, outcome,
                         B = 100L, seed = NA_integer_, mse_type = c("var", "mse"),
                         bootstrap_cores = 0L,
                         params = NULL,
                         level = c("subarea", "municipality"),
                         weight = c("census", "worldpop"), ...) {
  if (!requireNamespace("povmap", quietly = TRUE)) {
    stop("Package 'povmap' is required for fit_megb_sub().")
  }
  t0 <- Sys.time()
  mse_type <- match.arg(mse_type)
  level <- match.arg(level)
  weight <- match.arg(weight)
  inp <- .subarea_inputs(direct, smp_sub, pop_sub, outcome, level = level,
                         weight = weight)

  args <- list(
    fixed = inp$fixed,
    smp_data = inp$sm, smp_weights = "wt",
    pop_data = inp$pop, pop_weights = inp$wpop,
    domains = "municipality",
    transformation = "no",
    mse = TRUE, mse_type = mse_type, B = as.integer(B),
    bootstrap_cores = as.integer(bootstrap_cores), na.rm = FALSE,
    seed = if (is.na(seed)) NULL else as.integer(seed),
    cv_folds = "domain", cv_nfold = 10L,
    predict_sampled = "crossfit", predict_unsampled = "foldmean",
    early_stopping = TRUE)
  if (!is.null(params)) args$gradient_params <- params

  fit_raw <- tryCatch(
    suppressMessages(do.call(povmap::megb, c(args, list(...)))),
    error = function(e) e)

  if (inherits(fit_raw, "error")) {
    return(.povmap_failure("megb_sub", inp$in_ids, inp$oos_ids, t0,
                           err_msg = conditionMessage(fit_raw)))
  }

  est <- .povmap_estimates_from_fit(fit_raw, inp$in_ids, inp$oos_ids)

  # Variance components of the mixed model behind megb, recorded per fit. megb
  # stores standard deviations; the variances are what the FH family reports as
  # sigma2_v, so both are emitted and the comparison is direct.
  #
  # A random-effect standard deviation near zero is a result, not a failure: it
  # says the sub-area covariates already explain the between-municipality
  # variation, leaving no area effect. It is the megb counterpart of fh's
  # sigma2_v at its interval limit, so it is recorded per fit.
  num <- function(x) if (is.null(x) || !length(x)) NA_real_ else as.numeric(x)[1]
  ran_sd <- num(fit_raw$megb_model$ran_eff_sd)
  err_sd <- num(fit_raw$megb_model$error_sd)
  # Fraction of bootstrap replicates whose area effect collapsed: present only
  # when megb ran its bootstrap and reported diagnostics.
  bcol <- if (!is.null(fit_raw$boot_diag$boot_ran_eff_sd)) {
    mean(fit_raw$boot_diag$boot_ran_eff_sd < 1e-4, na.rm = TRUE)
  } else NA_real_
  best_iter <- num(fit_raw$megb_model$best_iter)
  chr <- function(x) if (is.null(x) || !length(x)) NA_character_ else as.character(x)[1]

  diag <- tibble::tibble(
    metric = c("B", "n_subarea_train", "n_subarea_pop", "n_covariates",
               "degenerate_y_share", "level_municipality",
               "ran_eff_sd", "ran_eff_var", "error_sd", "error_var",
               "boot_ran_eff_collapsed_share", "cv_nfold", "best_iter",
               "cv_folds", "predict_sampled", "predict_unsampled", "mse_type"),
    value  = c(as.numeric(B), nrow(inp$sm), nrow(inp$pop), length(inp$x_cols),
               mean(inp$sm$y %in% c(0, 1)),
               as.numeric(identical(level, "municipality")),
               ran_sd, ran_sd^2, err_sd, err_sd^2, bcol,
               num(fit_raw$cv_nfold), best_iter, NA_real_, NA_real_, NA_real_,
               NA_real_),
    character_value = c(rep(NA_character_, 13),
                        args$cv_folds, chr(fit_raw$predict_sampled),
                        chr(fit_raw$predict_unsampled), mse_type))
  meta <- .povmap_meta("megb_sub", est, inp$in_ids, inp$oos_ids, t0)
  meta$seed <- if (is.na(seed)) NA_integer_ else as.integer(seed)

  fit <- structure(list(estimates = est, diagnostics = diag, meta = meta),
                   class = c("sae_fit", "sae_fit_megb_sub"))
  sae_fit_contract(fit)
  fit
}
