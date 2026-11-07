# Helpers shared by the povmap-backed sub-area methods (R/method_subarea_common.R).

.povmap_raw <- function(domains = c("1", "2", "3"), mean = c(0.1, 0.2, 0.3),
                        var = c(1e-4, 4e-4, 9e-4), lower = mean - 0.02,
                        upper = mean + 0.02) {
  list(ind = data.frame(Domain = domains, Mean = mean),
       var = data.frame(Domain = domains, Mean = var),
       CI  = data.frame(Domain = domains, Lower = lower, Upper = upper))
}

test_that(".povmap_sha reports the installed commit", {
  skip_if_not_installed("povmap")
  expect_match(saecol2026:::.povmap_sha(), "^[0-9a-f]{40}$")
})

test_that(".povmap_estimates_from_fit maps povmap's tables to the contract shape", {
  raw <- .povmap_raw()
  est <- saecol2026:::.povmap_estimates_from_fit(raw, in_ids = c(1L, 3L), oos_ids = 2L)
  expect_s3_class(est, "tbl_df")
  expect_true(all(c("area_id", "theta_hat", "mse", "ci_lo", "ci_hi", "area_set")
                  %in% names(est)))
  # in-sample first, then out-of-sample, each ordered by area
  expect_identical(est$area_id, c(1L, 3L, 2L))
  expect_identical(est$area_set, c("in_sample", "in_sample", "out_of_sample"))
  expect_equal(est$theta_hat, c(0.1, 0.3, 0.2))
  expect_equal(est$mse, c(1e-4, 9e-4, 4e-4))
  expect_equal(est$ci_lo, c(0.08, 0.28, 0.18))
  expect_equal(est$ci_hi, c(0.12, 0.32, 0.22))
})

test_that(".povmap_estimates_from_fit replaces a missing or collapsed interval by a Wald one", {
  raw <- .povmap_raw(lower = c(0.08, 0.2, NA), upper = c(0.12, 0.2, NA))
  est <- saecol2026:::.povmap_estimates_from_fit(raw, in_ids = 1:3, oos_ids = integer())
  est <- est[order(est$area_id), ]
  expect_equal(est$ci_lo[1], 0.08)
  # area 2: interval collapsed to a point; area 3: missing; both use 1.96 * sqrt(mse)
  expect_equal(est$ci_lo[2:3], est$theta_hat[2:3] - 1.96 * sqrt(est$mse[2:3]))
  expect_equal(est$ci_hi[2:3], est$theta_hat[2:3] + 1.96 * sqrt(est$mse[2:3]))
  expect_true(all(est$ci_lo < est$ci_hi))
})

test_that(".povmap_estimates_from_fit keeps a positive variance floor for absent variances", {
  raw <- .povmap_raw()
  raw$var <- NULL
  raw$CI <- NULL
  est <- saecol2026:::.povmap_estimates_from_fit(raw, in_ids = 1:3, oos_ids = integer())
  expect_true(all(est$mse > 0))
  expect_true(all(est$ci_lo < est$ci_hi))
})

test_that(".povmap_failure returns NA estimates for every area and names the error", {
  fit <- saecol2026:::.povmap_failure("xgb_sub", in_ids = c(1L, 2L), oos_ids = 3L,
                                      t0 = Sys.time(), err_msg = "boom")
  expect_s3_class(fit, c("sae_fit", "sae_fit_xgb_sub"))
  expect_identical(fit$estimates$area_id, 1:3)
  expect_identical(fit$estimates$area_set,
                   c("in_sample", "in_sample", "out_of_sample"))
  expect_true(all(is.na(fit$estimates[c("theta_hat", "mse", "ci_lo", "ci_hi")])))
  expect_identical(fit$diagnostics$metric, c("xgb_sub_failed", "xgb_sub_error_message"))
  expect_identical(fit$diagnostics$character_value[2], "boom")
  expect_identical(fit$meta$n_in_sample, 2L)
  expect_identical(fit$meta$n_out_of_sample, 1L)
})
