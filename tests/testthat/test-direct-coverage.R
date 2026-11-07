test_that("5-rep mini-simulation: FH coverage > 0.5 on synthetic frame", {
  skip_if_not_installed("sae")
  fr <- make_synthetic_frame()

  design <- structure(list(
    stage1 = list(unit = "municipality", selection = "pps_wor",
                  size_variable = "n_dwellings"),
    stage2 = list(unit = "segment", selection = "systematic_fractional",
                  target_segment_size = 10L, n_per_psu = 3L),
    stage3 = list(unit = "household", selection = "census"),
    strata = data.frame(stratum = c("99_urban", "99_rural"),
                        n_psu   = c(2L, 2L),
                        n_seg_per_psu = c(3L, 3L)),
    calibration = list(margins = "stratum")
  ), class = c("ecv_design_spec", "list"))

  pop_totals <- as.data.frame(table(fr$stratum), stringsAsFactors = FALSE)
  names(pop_totals) <- c("stratum", "N")
  truth <- compute_truth(fr, area_var = "municipality")
  X <- truth[, "area_id", drop = FALSE]
  X$municipality <- X$area_id
  set.seed(123)
  X$x1 <- rnorm(nrow(X))
  X <- X[, c("municipality", "x1")]

  set.seed(2026)
  results <- list()
  for (r in 1:5) {
    s <- draw_replicate(fr, design, seed = sample.int(1e6, 1))
    s <- apply_calibration(s, pop_totals)
    d <- suppressWarnings(direct_estimate(s, area_var = "municipality"))
    d <- smooth_variance(d, method = "ghs")
    # On some replicates of this small frame sigma2_v sits at the limit of
    # emdi2's search interval, which emdi2 reports with a warning.
    fit <- muffle_warnings(fit_fh(d, X), "interval limit")
    fit$estimates$replicate <- r
    fit$estimates$method <- "fh"
    results[[r]] <- fit$estimates
  }
  est <- do.call(rbind, results)
  res <- evaluate_estimates(est, truth)

  cov <- res$by_method$coverage_95[res$by_method$method == "fh"]
  expect_gt(cov, 0.5,
            label = sprintf("FH coverage_95 = %.2f", cov))
})
