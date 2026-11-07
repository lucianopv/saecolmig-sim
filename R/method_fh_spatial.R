#' Fit a spatial Fay-Herriot model (SAR random effects) via emdi2::fh
#'
#' Thin wrapper around fit_fh() with correlation = "spatial". emdi fits the
#' SAR random effects; the MSE is emdi's analytical MSE in sample and the
#' synthetic-predictor variance out of sample (see [fit_fh()]). The W matrix in
#' geo$W must span all areas in direct + X_oos.
#'
#' @param direct tibble (area_id, theta_direct, var_direct_smoothed, ...)
#' @param X tibble (municipality, x1, ..., xK)
#' @param geo list with element W (sparse Matrix, row-normalised)
#' @param X_oos optional tibble for held-out areas
#' @param ... passed through to fit_fh
#' @return list(estimates, diagnostics, meta) with class
#'   c("sae_fit", "sae_fit_fh_spatial")
#' @export
fit_fh_spatial <- function(direct, X, geo, X_oos = NULL, ...) {
  fit_fh(direct = direct, X = X, geo = geo, X_oos = X_oos,
         correlation = "spatial", ...)
}
