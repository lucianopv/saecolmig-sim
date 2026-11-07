#' Registered SAE methods
#'
#' Each entry maps a method name to its fit function. The method names appear
#' in run outputs and evaluation summaries. Add a method by writing
#' R/method_<name>.R (conforming to the contract of [sae_fit_contract()]) and
#' adding one line here.
#'
#' The area-level methods (`fh`, `fh_spatial`, `fh_robust`, `fh_arcsin`,
#' `bayes_fh_beta`, `bayes_fh_heavy`) fit municipal direct estimates and
#' covariates. The sub-area methods (`xgb_sub`, `megb_sub`) take
#' `(direct, smp_sub, pop_sub, outcome)`: a sample and a population at the level
#' of census blocks and rural sections, aggregated to municipality; see
#' [sae_methods_subarea].
#'
#' @export
sae_methods <- list(
  fh             = fit_fh,
  fh_spatial     = fit_fh_spatial,
  bayes_fh_beta  = fit_bayes_fh_beta,
  fh_robust      = fit_fh_robust,
  fh_arcsin      = fit_fh_arcsin,
  bayes_fh_heavy = fit_bayes_fh_heavy,
  xgb_sub        = fit_xgb_sub,
  megb_sub       = fit_megb_sub
)

#' Methods that are implemented
#' @export
sae_methods_implemented <- c("fh", "fh_spatial",
                             "bayes_fh_beta", "fh_robust", "fh_arcsin",
                             "bayes_fh_heavy",
                             "xgb_sub", "megb_sub")

#' Sub-area methods
#'
#' These model census blocks and rural sections and aggregate to municipality,
#' where they are scored against the same truth as the area-level methods. They
#' take the sub-area sample and population instead of municipal covariates.
#' @export
sae_methods_subarea <- c("xgb_sub", "megb_sub")

#' Verify a fit object satisfies the method contract
#'
#' Used by tests and by the pipeline to fail fast if a method returns the
#' wrong shape.
#'
#' @param fit return value of any fit_<name>()
#' @return invisibly TRUE; throws on contract violation
#' @export
sae_fit_contract <- function(fit) {
  if (!inherits(fit, "sae_fit")) stop("not an sae_fit object")
  required_top <- c("estimates", "diagnostics", "meta")
  miss <- setdiff(required_top, names(fit))
  if (length(miss)) stop("missing top-level: ", paste(miss, collapse = ", "))
  required_est <- c("area_id", "theta_hat", "mse", "ci_lo", "ci_hi")
  miss <- setdiff(required_est, colnames(fit$estimates))
  if (length(miss)) stop("estimates missing cols: ", paste(miss, collapse = ", "))
  if (anyNA(fit$estimates$mse)) stop("NA mse in estimates")
  if (!all(fit$estimates$mse > 0)) stop("non-positive mse")
  if (!all(fit$estimates$ci_lo < fit$estimates$ci_hi)) stop("ci_lo >= ci_hi")
  required_meta <- c("method", "runtime_sec", "n_areas")
  miss <- setdiff(required_meta, names(fit$meta))
  if (length(miss)) stop("meta missing: ", paste(miss, collapse = ", "))
  invisible(TRUE)
}
