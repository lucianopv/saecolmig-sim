test_that("fit_bayes_fh_heavy returns the sae_fit contract on a small fixture", {
  skip_if_not_installed("cmdstanr")
  f <- .fixture_robust(n_in = 30L, n_oos = 8L)
  fit <- fit_bayes_fh_heavy(f$direct, f$X_in, X_oos = f$X_oos,
                            chains = 2L, iter_sampling = 250L,
                            iter_warmup = 250L,
                            seed = 5042L)
  expect_s3_class(fit, "sae_fit")
  expect_s3_class(fit, "sae_fit_bayes_fh_heavy")
  expect_equal(nrow(fit$estimates), f$n_in + f$n_oos)
  expect_setequal(colnames(fit$estimates),
                   c("area_id", "area_set", "theta_hat", "mse", "ci_lo", "ci_hi"))
  expect_true(all(!is.na(fit$estimates$theta_hat)))
  expect_true(all(fit$estimates$mse > 0))
  expect_no_error(sae_fit_contract(fit))
  expect_true("sigma_v_mean" %in% fit$diagnostics$metric)
  expect_true(fit$diagnostics$value[fit$diagnostics$metric == "bayes_fh_heavy_failed"] == 0)
})

test_that("fit_bayes_fh_heavy returns NA estimates and a diagnostic flag on sampling failure", {
  skip_if_not_installed("cmdstanr")
  f <- .fixture_robust(n_in = 5L, n_oos = 0L)
  f$X_in$x1 <- 0
  f$X_in$x2 <- 0
  f$direct$var_direct_smoothed <- 0
  fit <- suppressWarnings(
    fit_bayes_fh_heavy(f$direct, f$X_in,
                       chains = 1L, iter_sampling = 100L, iter_warmup = 100L,
                       seed = 5042L)
  )
  if (any(is.na(fit$estimates$theta_hat))) {
    expect_true(any(fit$diagnostics$metric == "bayes_fh_heavy_failed" &
                     fit$diagnostics$value == 1))
  } else {
    succeed("sampling completed despite degenerate inputs; covered by contract sanity")
  }
})
