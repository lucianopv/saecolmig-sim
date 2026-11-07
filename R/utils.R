# Internal utilities shared across the package.

# Default-fallback helper. Returns `a` unless it is NULL, in which case `b`.
# Dot-prefixed so exportPattern("^[^.]") in NAMESPACE does NOT export it
# (avoids conflict with rlang::`%||%`).
.null_or <- function(a, b) if (is.null(a)) b else a

#' Build the argument list for an area-level method's fit function
#'
#' Centralises per-method argument dispatch for the area-level methods, which
#' take area-level `direct` + muni-level X. The sub-area methods
#' (`sae_methods_subarea`) take a sample and a population at sub-area level
#' instead and are called directly, see [fit_xgb_sub()] and [fit_megb_sub()];
#' asking for one here is an error.
#'
#' @param method_name character; an area-level method of
#'   `sae_methods_implemented`, i.e. not one of `sae_methods_subarea`
#' @param direct one rep's direct estimates (tibble). When it contains an
#'   `outcome` column (outcome-long format), only rows matching `outcome` are
#'   kept before dispatch.
#' @param aux_in,aux_oos muni-level auxiliary X for in-sample / OOS areas
#' @param W_matrix_full queens-contiguity W (1,122 munis); used only by
#'   fh_spatial.
#' @param outcome character; registered outcome name (see `sae_outcomes`).
#'   `direct` is pre-filtered to rows matching this outcome when `direct`
#'   carries an `outcome` column.
#' @param seed integer; forwarded to `fh_arcsin` (bootstrap) and to the two
#'   Bayesian methods (sampler). NA leaves the Bayesian fits unseeded. To
#'   re-run a fit, use [sae_fit_seed()].
#' @return list ready to pass to `do.call(sae_methods[[method_name]], args)`
#' @keywords internal
method_args <- function(method_name, direct,
                        aux_in, aux_oos,
                        W_matrix_full = NULL,
                        outcome = "migrated_last_12m",
                        seed = NA_integer_) {
  if (method_name %in% sae_methods_subarea) {
    stop("method_args() is for the area-level methods; ", method_name,
         " is a sub-area method and takes its arguments directly ",
         "(direct, smp_sub, pop_sub, outcome), see fit_xgb_sub() and ",
         "fit_megb_sub().", call. = FALSE)
  }
  # Filter outcome-long direct to the requested outcome.
  if ("outcome" %in% colnames(direct)) {
    direct <- direct[direct$outcome == outcome, , drop = FALSE]
    if (nrow(direct) == 0) stop("No direct estimates for outcome ", outcome)
  }
  args <- list(direct = direct, X = aux_in, X_oos = aux_oos)
  if (method_name == "fh_spatial") {
    args$geo <- list(W = W_matrix_full)
  }
  if (method_name == "fh_arcsin") {
    # The emdi2 arcsin bootstrap is cheap (~0.1 s per fit), so B = 500 is
    # affordable and keeps the Monte Carlo noise of the per-area MSE small; the
    # per-replicate `seed` makes the bootstrap reproducible.
    args$B    <- 500L
    args$seed <- seed
  }
  if (method_name == "bayes_fh_beta") {
    # ExtBeta likelihood with a smoothed dispersion, which keeps the coverage of
    # rare outcomes close to nominal without divergent transitions.
    args$likelihood <- "ExtBeta"
    args$iter       <- 2000L
    # Three chains rather than the function default of four; chains run in
    # parallel, so the count does not change per-fit wall-clock time.
    args$chains     <- 3L
    # The sampler seed.
    args$seed       <- seed
  }
  if (method_name == "bayes_fh_heavy") {
    args$iter_sampling <- 1000L  # cmdstanr post-warmup samples per chain.
    # fit_bayes_fh_heavy's default is NULL, not NA.
    if (!is.na(seed)) args$seed <- seed
  }
  args
}
