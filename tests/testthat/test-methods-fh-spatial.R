test_that("fit_fh_spatial returns the contract schema with a W matrix", {
  skip_if_not_installed("emdi2")
  skip_if_not_installed("Matrix")

  set.seed(1)
  n_areas  <- 30L
  area_ids <- as.integer(c(rep(5L, 10), rep(11L, 10), rep(13L, 10)) * 1000L +
                          rep(seq_len(10), 3))
  X <- tibble::tibble(municipality = area_ids,
                      x1 = rnorm(n_areas),
                      x2 = rnorm(n_areas))
  theta_true <- 0.05 + 0.02 * X$x1 - 0.01 * X$x2 + rnorm(n_areas, 0, 0.01)
  theta_true <- pmin(pmax(theta_true, 0.01), 0.5)
  n_eff <- runif(n_areas, 50, 200)
  var_d <- theta_true * (1 - theta_true) / n_eff
  direct <- tibble::tibble(
    area_id = area_ids,
    theta_direct = theta_true + rnorm(n_areas, 0, sqrt(var_d)),
    var_direct = var_d, var_direct_smoothed = var_d,
    n_eff = n_eff, n_dwell = round(n_eff)
  )

  W <- build_dept_adjacency(area_ids)
  res <- fit_fh_spatial(direct, X, geo = list(W = W))

  expect_s3_class(res, "sae_fit_fh_spatial")
  for (col in c("area_id", "theta_hat", "mse", "ci_lo", "ci_hi")) {
    expect_true(col %in% colnames(res$estimates), info = col)
  }
  expect_equal(nrow(res$estimates), n_areas)
  expect_true(all(res$estimates$mse > 0))
  expect_true(all(res$estimates$ci_lo < res$estimates$ci_hi))
  expect_equal(res$meta$method, "fh_spatial")
})

test_that("fit_fh_spatial passes through to fit_fh with correlation='spatial' and OOS", {
  skip_if_not_installed("emdi2")
  skip_if_not_installed("Matrix")

  set.seed(7)
  n <- 30L
  area_ids <- as.integer(c(rep(5L, 10), rep(11L, 10), rep(13L, 10)) * 1000L +
                          rep(seq_len(10), 3))
  X <- tibble::tibble(municipality = area_ids,
                      x1 = rnorm(n), x2 = rnorm(n))
  theta_true <- pmin(pmax(0.05 + 0.02 * X$x1 - 0.01 * X$x2 +
                            rnorm(n, 0, 0.01), 0.005), 0.5)
  n_eff <- runif(n, 50, 200)
  var_d <- theta_true * (1 - theta_true) / n_eff

  in_ids <- area_ids[1:20]
  oos_ids <- area_ids[21:30]
  direct_in <- tibble::tibble(
    area_id              = in_ids,
    theta_direct         = theta_true[1:20] + rnorm(20, 0, sqrt(var_d[1:20])),
    var_direct           = var_d[1:20],
    var_direct_smoothed  = var_d[1:20]
  )
  W <- build_dept_adjacency(area_ids)
  X_in  <- X[X$municipality %in% in_ids, ]
  X_oos <- X[X$municipality %in% oos_ids, ]

  res <- suppressWarnings(
    fit_fh_spatial(direct_in, X_in, geo = list(W = W), X_oos = X_oos)
  )

  expect_s3_class(res, "sae_fit_fh_spatial")
  expect_equal(res$meta$method, "fh_spatial")
  expect_equal(sum(res$estimates$area_set == "out_of_sample"), 10L)
  expect_true(all(res$estimates$mse > 0))

  # OOS MSE is the per-area synthetic-predictor variance (sigma2_v + x'Vx),
  # so it varies across the OOS areas.
  oos_mse <- res$estimates$mse[res$estimates$area_set == "out_of_sample"]
  expect_gt(length(unique(round(oos_mse, 12))), 1L)
})
