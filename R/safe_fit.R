#' Run a registered SAE method, converting errors into recorded failures
#'
#' A fit can fail in two ways: the fit function raises an error (a singular
#' fit, for example), or it returns an object that breaks the method contract
#' of [sae_fit_contract()] (a wrapper that catches its own fitting error and
#' returns NA estimates). Both the fit call and the contract check run inside
#' one error handler, so either failure returns one NA-estimate row per
#' requested area with `fit_failed = TRUE` and the first 200 characters of the
#' error message; on success the contract-checked estimates are returned with
#' `fit_failed = FALSE`. The failure rate by method and outcome is itself a
#' result, so a failure is recorded rather than raised.
#'
#' @param method_name name in `sae_methods`
#' @param args list passed to the fit via do.call (from method_args())
#' @param area_ids area ids used to shape the NA rows on failure
#' @param diagnostics logical; if TRUE attach the fit's diagnostics tibble as
#'   `attr(est, "diagnostics")`. On a failure it is the diagnostics of the
#'   object the fit function returned when that object failed the contract
#'   (where a wrapper reports its own error message), and NULL when the fit
#'   function raised.
#' @return tibble of estimates with fit_failed + error_msg columns
#' @export
safe_method_fit <- function(method_name, args, area_ids, diagnostics = FALSE) {
  fit <- NULL
  est <- tryCatch({
    fit <- do.call(sae_methods[[method_name]], args)
    sae_fit_contract(fit)
    e <- tibble::as_tibble(fit$estimates)
    if (!"area_set" %in% colnames(e)) e$area_set <- "in_sample"
    e$fit_failed <- FALSE
    e$error_msg  <- NA_character_
    if (isTRUE(diagnostics)) attr(e, "diagnostics") <- fit$diagnostics
    e
  }, error = function(err) {
    e <- tibble::tibble(
      area_id    = area_ids,
      area_set   = "in_sample",
      theta_hat  = NA_real_,
      mse        = NA_real_,
      ci_lo      = NA_real_,
      ci_hi      = NA_real_,
      fit_failed = TRUE,
      error_msg  = substr(conditionMessage(err), 1, 200))
    if (isTRUE(diagnostics) && inherits(fit, "sae_fit")) {
      attr(e, "diagnostics") <- fit$diagnostics
    }
    e
  })
  est
}
