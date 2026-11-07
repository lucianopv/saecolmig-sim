test_that("fit_bayes_fh_beta returns the contract schema and converges on small synthetic input", {
  skip_if_not_installed("tipsae")
  skip_if_not_installed("rstan")
  skip_on_ci()

  set.seed(2026)
  n_areas <- 30L
  area_ids <- as.integer(seq_len(n_areas))
  X <- tibble::tibble(municipality = area_ids,
                      x1 = rnorm(n_areas),
                      x2 = rnorm(n_areas))
  theta_true <- pmin(pmax(0.05 + 0.02 * X$x1 - 0.01 * X$x2 +
                            rnorm(n_areas, 0, 0.01), 0.005), 0.5)
  n_eff <- runif(n_areas, 50, 200)
  var_d <- theta_true * (1 - theta_true) / n_eff
  direct <- tibble::tibble(
    area_id              = area_ids,
    theta_direct         = pmin(pmax(theta_true + rnorm(n_areas, 0, sqrt(var_d)),
                                     0.001), 0.999),
    var_direct           = var_d,
    var_direct_smoothed  = var_d,
    n_eff                = n_eff,
    n_dwell              = round(n_eff)
  )

  res <- suppressWarnings(
    fit_bayes_fh_beta(direct, X,
                      chains = 2L, iter = 500L, warmup = 250L,
                      seed = 1L)
  )

  expect_s3_class(res, "sae_fit")
  expect_s3_class(res, "sae_fit_bayes_fh_beta")
  for (col in c("area_id", "theta_hat", "mse", "ci_lo", "ci_hi")) {
    expect_true(col %in% colnames(res$estimates), info = col)
  }
  expect_equal(nrow(res$estimates), n_areas)
  expect_true(all(res$estimates$mse > 0))
  expect_true(all(res$estimates$ci_lo < res$estimates$ci_hi))
  expect_equal(res$meta$method, "bayes_fh_beta")
  # Allow generous time on the tiny fixture (Stan compile is amortised
  # across the suite but a fresh run can hit ~120 s)
  expect_true(res$meta$runtime_sec > 0,
              info = sprintf("runtime = %.1fs", res$meta$runtime_sec))
})

test_that("fit_bayes_fh_beta returns native OOS predictions via likelihood='beta'", {
  skip_if_not_installed("tipsae")
  skip_if_not_installed("rstan")
  skip_on_ci()

  set.seed(2026)
  n_areas <- 30L
  area_ids <- as.integer(seq_len(n_areas))
  X <- tibble::tibble(municipality = area_ids,
                      x1 = rnorm(n_areas), x2 = rnorm(n_areas))
  theta_true <- pmin(pmax(0.05 + 0.02 * X$x1 - 0.01 * X$x2 +
                            rnorm(n_areas, 0, 0.01), 0.005), 0.5)
  n_eff <- runif(n_areas, 50, 200)
  var_d <- theta_true * (1 - theta_true) / n_eff

  in_ids  <- area_ids[1:20]
  oos_ids <- area_ids[21:30]
  direct_in <- tibble::tibble(
    area_id              = in_ids,
    theta_direct         = pmin(pmax(theta_true[1:20] +
                                      rnorm(20, 0, sqrt(var_d[1:20])),
                                    0.001), 0.999),
    var_direct           = var_d[1:20],
    var_direct_smoothed  = var_d[1:20]
  )
  X_in  <- X[X$municipality %in% in_ids, ]
  X_oos <- X[X$municipality %in% oos_ids, ]

  res <- suppressWarnings(
    fit_bayes_fh_beta(direct_in, X_in, X_oos = X_oos,
                      chains = 2L, iter = 500L, warmup = 250L, seed = 1L)
  )

  expect_true("area_set" %in% colnames(res$estimates))
  expect_equal(sum(res$estimates$area_set == "in_sample"), 20L)
  expect_equal(sum(res$estimates$area_set == "out_of_sample"), 10L)
  expect_true(all(res$estimates$mse > 0))
})

test_that("fit_bayes_fh_beta forwards adapt_delta and reports divergences", {
  skip_if_not_installed("tipsae")
  skip_if_not_installed("rstan")
  skip_on_ci()

  # .fixture_robust() holds out at least one area; X_oos = NULL leaves the fit
  # without out-of-sample areas.
  fx <- .fixture_robust(n_in = 30L, n_oos = 1L, seed = 7L)
  d  <- fx$direct
  d$n_dwell <- 40L
  fit <- suppressWarnings(
    fit_bayes_fh_beta(d, fx$X_in, X_oos = NULL,
                      chains = 1L, iter = 600L,
                      warmup = 300L, seed = 42L,
                      adapt_delta = 0.95)
  )
  expect_s3_class(fit, "sae_fit")
  expect_true("n_divergent" %in% fit$diagnostics$metric)
  expect_true(is.finite(
    fit$diagnostics$value[fit$diagnostics$metric == "n_divergent"]))
})

test_that("Infbeta0 fits exact zeros", {
  skip_if_not_installed("tipsae")
  skip_if_not_installed("rstan")
  skip_on_ci()

  # n_oos = 0L would break .fixture_robust indexing; use n_oos = 1L, X_oos = NULL.
  fx <- .fixture_robust(n_in = 30L, n_oos = 1L, seed = 8L)
  d  <- fx$direct
  d$n_dwell <- 40L
  d$theta_direct[1:12] <- 0           # rare-outcome pattern
  fit <- suppressWarnings(
    fit_bayes_fh_beta(d, fx$X_in, X_oos = NULL,
                      chains = 1L, iter = 600L,
                      warmup = 300L, seed = 43L,
                      likelihood = "Infbeta0")
  )
  expect_s3_class(fit, "sae_fit")
  expect_true(all(fit$estimates$theta_hat >= 0))
  expect_true(all(fit$estimates$mse > 0))
})
