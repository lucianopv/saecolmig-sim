test_that("fit_fh_arcsin returns the sae_fit contract on a 60-area fixture", {
  skip_if_not_installed("emdi2")
  f <- .fixture_robust()
  fit <- fit_fh_arcsin(f$direct, f$X_in, X_oos = f$X_oos)
  expect_s3_class(fit, "sae_fit")
  expect_s3_class(fit, "sae_fit_fh_arcsin")
  expect_true(all(c("estimates", "diagnostics", "meta") %in% names(fit)))
  expect_equal(nrow(fit$estimates), f$n_in + f$n_oos)
  expect_setequal(colnames(fit$estimates),
                   c("area_id", "area_set", "theta_hat", "mse", "ci_lo", "ci_hi"))
  expect_equal(sum(fit$estimates$area_set == "in_sample"), f$n_in)
  expect_equal(sum(fit$estimates$area_set == "out_of_sample"), f$n_oos)
  expect_true(all(!is.na(fit$estimates$theta_hat)))
  expect_true(all(fit$estimates$mse > 0))
  expect_no_error(sae_fit_contract(fit))
})

test_that("fit_fh_arcsin theta_hat respects the [0, 1] support after back-transform", {
  skip_if_not_installed("emdi2")
  f <- .fixture_robust()
  fit <- fit_fh_arcsin(f$direct, f$X_in, X_oos = f$X_oos)
  expect_true(all(fit$estimates$theta_hat >= 0 - 1e-8))
  expect_true(all(fit$estimates$theta_hat <= 1 + 1e-8))
})

test_that("fit_fh_arcsin pre-clips theta_direct outside the open (0, 1) interval", {
  skip_if_not_installed("emdi2")
  f <- .fixture_robust()
  # Inject a zero direct estimate at row 1
  f$direct$theta_direct[1] <- 0
  expect_no_error(
    fit <- fit_fh_arcsin(f$direct, f$X_in, X_oos = f$X_oos)
  )
  # The pre-clip should leave theta_hat well-defined for all areas
  expect_true(all(!is.na(fit$estimates$theta_hat)))
})

test_that("fit_fh_arcsin falls back to n when n_eff is absent and emits a message", {
  skip_if_not_installed("emdi2")
  f <- .fixture_robust()
  # Drop n_eff; provide n instead
  direct_no_neff <- f$direct
  direct_no_neff$n <- direct_no_neff$n_eff
  direct_no_neff$n_eff <- NULL
  expect_message(
    fit <- fit_fh_arcsin(direct_no_neff, f$X_in, X_oos = f$X_oos),
    "n_eff.+falling back to.+n"
  )
  expect_s3_class(fit, "sae_fit_fh_arcsin")
  expect_true(all(!is.na(fit$estimates$theta_hat)))
})

test_that("fit_fh_arcsin has these arguments, diagnostics and the seed in its meta", {
  expect_setequal(names(formals(fit_fh_arcsin)),
                  c("direct", "X", "geo", "X_oos", "backtransformation", "B", "seed", "..."))
  skip_if_not_installed("emdi2")
  f <- .fixture_robust()
  fit <- fit_fh_arcsin(f$direct, f$X_in, X_oos = f$X_oos, seed = 5L)
  expect_setequal(fit$diagnostics$metric,
                  c("sigma2_v", "transformation", "backtransformation", "B_bootstrap",
                    "fh_arcsin_failed"))
  expect_identical(fit$meta$seed, 5L)
})
