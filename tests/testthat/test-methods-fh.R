.tiny_method_inputs <- function(seed = 1, n_areas = 30) {
  set.seed(seed)
  area_ids <- as.integer(seq_len(n_areas))
  X <- tibble::tibble(
    municipality = area_ids,
    x1 = rnorm(n_areas),
    x2 = rnorm(n_areas)
  )
  beta <- c(intercept = 0.05, x1 = 0.02, x2 = -0.01)
  theta_true <- pmin(pmax(beta["intercept"] + beta["x1"] * X$x1 +
                            beta["x2"] * X$x2 + rnorm(n_areas, 0, 0.01),
                          0), 0.5)
  n_eff <- runif(n_areas, 20, 200)
  var_d <- theta_true * (1 - theta_true) / n_eff
  direct <- tibble::tibble(
    area_id              = area_ids,
    theta_direct         = theta_true + rnorm(n_areas, 0, sqrt(var_d)),
    var_direct           = var_d,
    var_direct_smoothed  = var_d,
    n_eff                = n_eff,
    n_dwell              = round(n_eff)
  )
  list(direct = direct, X = X, theta_true = theta_true)
}

test_that("fit_fh returns the contract schema and shrinks toward the model fit", {
  skip_if_not_installed("emdi2")
  inp <- .tiny_method_inputs()
  # On this fixture the REML estimate of sigma2_v sits at the limit of emdi2's
  # search interval, which emdi2 reports with a warning.
  expect_warning(res <- fit_fh(inp$direct, inp$X), "interval limit")

  expect_s3_class(res, "sae_fit")
  expect_s3_class(res, "sae_fit_fh")
  for (col in c("area_id", "theta_hat", "mse", "ci_lo", "ci_hi")) {
    expect_true(col %in% colnames(res$estimates), info = col)
  }
  expect_equal(nrow(res$estimates), nrow(inp$direct))
  expect_true(all(res$estimates$mse > 0))
  expect_true(all(res$estimates$ci_lo < res$estimates$ci_hi))
  expect_true(res$meta$method == "fh")
  expect_true(res$meta$runtime_sec > 0)

  mae_fh     <- mean(abs(res$estimates$theta_hat - inp$theta_true))
  mae_direct <- mean(abs(inp$direct$theta_direct - inp$theta_true))
  expect_lt(mae_fh, mae_direct)
})

test_that("fit_fh handles X_oos via emdi combined_data", {
  skip_if_not_installed("emdi2")
  inp <- .tiny_method_inputs(n_areas = 30)
  in_ids  <- inp$direct$area_id[1:20]
  oos_ids <- inp$direct$area_id[21:30]
  direct_in <- inp$direct[inp$direct$area_id %in% in_ids, ]
  X_in      <- inp$X[inp$X$municipality %in% in_ids, ]
  X_oos     <- inp$X[inp$X$municipality %in% oos_ids, ]

  expect_warning(res <- fit_fh(direct_in, X_in, X_oos = X_oos), "interval limit")

  expect_true("area_set" %in% colnames(res$estimates))
  expect_setequal(unique(res$estimates$area_set),
                  c("in_sample", "out_of_sample"))
  expect_equal(sum(res$estimates$area_set == "in_sample"), 20L)
  expect_equal(sum(res$estimates$area_set == "out_of_sample"), 10L)
  expect_true(all(res$estimates$mse > 0))
  expect_true(all(res$estimates$ci_lo < res$estimates$ci_hi))
})
