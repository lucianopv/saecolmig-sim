#' Sub-area XGBoost EBP via povmap::xgb
#'
#' Fits `povmap::xgb()` on census-block and rural-section aggregates and lets
#' povmap aggregate the sub-area predictions to municipality, where the
#' estimates are scored against the same census truth as the area-level
#' methods.
#'
#' `povmap::xgb()` is the povmap entry point that takes `sub_domains`, so the
#' two-level structure is declared explicitly: `domains = "municipality"`,
#' `sub_domains = "code"`. Aggregation from sub-area to municipality is
#' population-weighted internally by `pop_weights`.
#'
#' `transformation = "no"` is fixed: a realised sub-area is about one
#' ten-dwelling segment, so most have `y` exactly 0 or 1, and no logit or
#' arcsin transformation survives that.
#'
#' Uncertainty is povmap's residual bootstrap. `bootstrap = TRUE` is required:
#' povmap builds `$CI` only under bootstrap, so without it the method would
#' return point estimates with no interval. Estimates are unbenchmarked,
#' because povmap's bootstrap interval is built around the unbenchmarked
#' estimator and pairing a benchmarked point with it would not be coherent.
#'
#' povmap validates against the xgboost version pinned in DESCRIPTION and
#' warns, rather than stops, on a mismatch; check the version if results move
#' unexpectedly.
#'
#' @param direct tibble of direct estimates for this outcome; used only to
#'   identify drawn municipalities, never fed to the model
#' @param smp_sub collapsed sub-area sample for this outcome
#' @param pop_sub sub-area covariate matrix over the whole universe
#' @param outcome outcome name from the registry
#' @param B bootstrap replications; default 200
#' @param seed integer seed passed to povmap; `NA` passes none
#' @param cpus cores for povmap's bootstrap; default 1
#' @param params optional named list of xgboost hyperparameters (nrounds, eta,
#'   max_depth, min_child_weight, subsample, colsample_bytree, lambda, alpha);
#'   anything left out takes the default (200, 0.1, 4, 1, 0.8, 0.8, 1, 0)
#' @param level `"subarea"` (default) or `"municipality"`; the latter
#'   aggregates the predictors to municipality and fits the same estimator
#'   there, which separates resolution from method.
#' @param weight `"census"` (default) or `"worldpop"`; the population weight,
#'   see `.subarea_inputs()`.
#' @param ... forwarded to `povmap::xgb()`
#' @return list(estimates, diagnostics, meta), class c("sae_fit", "sae_fit_xgb_sub")
#' @export
fit_xgb_sub <- function(direct, smp_sub, pop_sub, outcome,
                        B = 200L, seed = NA_integer_, cpus = 1L,
                        params = list(),
                        level = c("subarea", "municipality"),
                        weight = c("census", "worldpop"), ...) {
  if (!requireNamespace("povmap", quietly = TRUE)) {
    stop("Package 'povmap' is required for fit_xgb_sub().")
  }
  t0 <- Sys.time()
  level <- match.arg(level)
  weight <- match.arg(weight)
  inp <- .subarea_inputs(direct, smp_sub, pop_sub, outcome, level = level,
                         weight = weight)

  dflt <- list(nrounds = 200L, eta = 0.1, max_depth = 4L,
               min_child_weight = 1, subsample = 0.8,
               colsample_bytree = 0.8, lambda = 1, alpha = 0)
  p <- utils::modifyList(dflt, params)

  fit_raw <- tryCatch(
    povmap::xgb(
      fixed = inp$fixed,
      smp_data = inp$sm, smp_weights = "wt",
      pop_data = inp$pop, pop_weights = inp$wpop,
      domains = "municipality",
      # povmap::xgb requires `sub_domains`. At municipality level
      # .subarea_inputs() sets `code` to the municipality itself, so each domain
      # is its own single cell and the aggregation step is the identity.
      sub_domains = "code",
      transformation = "no",
      bootstrap = TRUE, bootstrap_type = "residual", B = as.integer(B),
      cpus = as.integer(cpus), na.rm = FALSE,
      seed = if (is.na(seed)) NULL else as.integer(seed),
      nrounds = as.integer(p$nrounds), eta = p$eta,
      max_depth = as.integer(p$max_depth),
      min_child_weight = p$min_child_weight, subsample = p$subsample,
      colsample_bytree = p$colsample_bytree,
      lambda = p$lambda, alpha = p$alpha,
      verbose = FALSE, ...),
    error = function(e) e)

  if (inherits(fit_raw, "error")) {
    return(.povmap_failure("xgb_sub", inp$in_ids, inp$oos_ids, t0,
                           err_msg = conditionMessage(fit_raw)))
  }

  est  <- .povmap_estimates_from_fit(fit_raw, inp$in_ids, inp$oos_ids)
  diag <- tibble::tibble(
    metric = c("B", "n_subarea_train", "n_subarea_pop", "n_covariates",
               "degenerate_y_share", "level_municipality",
               "nrounds", "eta", "max_depth"),
    value  = c(as.numeric(B), nrow(inp$sm), nrow(inp$pop), length(inp$x_cols),
               mean(inp$sm$y %in% c(0, 1)),
               as.numeric(identical(level, "municipality")),
               as.numeric(p$nrounds),
               p$eta, as.numeric(p$max_depth)),
    character_value = NA_character_)
  meta <- .povmap_meta("xgb_sub", est, inp$in_ids, inp$oos_ids, t0)
  meta$seed <- if (is.na(seed)) NA_integer_ else as.integer(seed)
  meta$xgboost_version <- as.character(utils::packageVersion("xgboost"))

  fit <- structure(list(estimates = est, diagnostics = diag, meta = meta),
                   class = c("sae_fit", "sae_fit_xgb_sub"))
  sae_fit_contract(fit)
  fit
}
