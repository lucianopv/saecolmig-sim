test_that("compute_truth returns one row per area with theta_true in [0,1]", {
  fr <- make_synthetic_frame()
  truth <- compute_truth(fr, area_var = "municipality")
  expect_equal(nrow(truth), length(unique(fr$municipality)))
  expect_setequal(unique(truth$area_id), unique(fr$municipality))
  expect_true(all(truth$theta_true >= 0 & truth$theta_true <= 1))
  # theta_true equals the unweighted mean of migrated_last_12m per muni
  for (m in unique(fr$municipality)[1:3]) {
    expected <- mean(fr$migrated_last_12m[fr$municipality == m], na.rm = TRUE)
    expect_equal(truth$theta_true[truth$area_id == m], expected, tolerance = 1e-9)
  }
})

.tiny_estimates <- function() {
  truth <- tibble::tibble(area_id    = c("A", "B", "C"),
                          theta_true = c(0.10, 0.20, 0.30))
  set.seed(42)
  rows <- list()
  for (m in c("A", "B")) for (r in 1:5) for (a in c("A", "B", "C")) {
    bias <- if (m == "B") 0.10 else 0
    th_true <- truth$theta_true[truth$area_id == a]
    rows[[length(rows) + 1L]] <- tibble::tibble(
      method     = m,
      replicate  = r,
      area_id    = a,
      theta_hat  = th_true + bias + rnorm(1, 0, 0.01),
      mse        = 0.0001,
      ci_lo      = th_true + bias - 0.05,
      ci_hi      = th_true + bias + 0.05
    )
  }
  list(estimates = do.call(rbind, rows), truth = truth)
}

test_that("evaluate_estimates returns by_method and by_method_area summaries", {
  fx <- .tiny_estimates()
  res <- evaluate_estimates(fx$estimates, fx$truth)

  expect_setequal(res$by_method$method, c("A", "B"))
  bias_A <- res$by_method$bias[res$by_method$method == "A"]
  bias_B <- res$by_method$bias[res$by_method$method == "B"]
  expect_lt(abs(bias_A), 0.02)
  expect_lt(abs(bias_B - 0.10), 0.02)

  expect_equal(nrow(res$by_method_area), 6)
  for (col in c("bias", "rmse", "coverage_95", "ci_length")) {
    expect_true(col %in% colnames(res$by_method_area), info = col)
  }

  cov_A <- mean(res$by_method_area$coverage_95[res$by_method_area$method == "A"])
  cov_B <- mean(res$by_method_area$coverage_95[res$by_method_area$method == "B"])
  expect_gt(cov_A, 0.9)
  expect_lt(cov_B, 0.1)

  # Rank correlation: theta_true is monotone, so both methods should hit ~1.
  expect_gt(res$by_method$rank_corr_mean[res$by_method$method == "A"], 0.9)
  expect_gt(res$by_method$rank_corr_mean[res$by_method$method == "B"], 0.9)
})

test_that("evaluate_estimates returns by_method_set when area_set is present", {
  set.seed(7)
  truth <- tibble::tibble(area_id    = c("A", "B", "C", "D"),
                          theta_true = c(0.10, 0.20, 0.30, 0.40))
  rows <- list()
  for (m in c("X", "Y")) for (r in 1:5) for (a in c("A", "B", "C", "D")) {
    set <- if (a %in% c("A", "B")) "in_sample" else "out_of_sample"
    th_true <- truth$theta_true[truth$area_id == a]
    bias <- if (set == "out_of_sample") 0.05 else 0
    rows[[length(rows) + 1L]] <- tibble::tibble(
      method     = m,
      replicate  = r,
      area_id    = a,
      area_set   = set,
      theta_hat  = th_true + bias + rnorm(1, 0, 0.01),
      mse        = 0.0001,
      ci_lo      = th_true + bias - 0.05,
      ci_hi      = th_true + bias + 0.05
    )
  }
  est <- do.call(rbind, rows)
  res <- evaluate_estimates(est, truth)

  expect_true("by_method_set" %in% names(res))
  expect_equal(nrow(res$by_method_set), 4L)
  expect_setequal(unique(res$by_method_set$area_set),
                  c("in_sample", "out_of_sample"))
  bias_oos <- res$by_method_set$bias[res$by_method_set$area_set == "out_of_sample"]
  expect_true(all(abs(bias_oos - 0.05) < 0.02))
})

test_that("evaluate_estimates adds rmse_reported and calibration_ratio at all three grains", {
  fx <- .tiny_estimates()
  res <- evaluate_estimates(fx$estimates, fx$truth)

  for (col in c("rmse_reported", "calibration_ratio")) {
    expect_true(col %in% colnames(res$by_method),      info = col)
    expect_true(col %in% colnames(res$by_method_set),  info = col)
    expect_true(col %in% colnames(res$by_method_area), info = col)
  }

  # mse is a constant 0.0001 in .tiny_estimates(), so rmse_reported is exact.
  expect_equal(res$by_method$rmse_reported, rep(sqrt(0.0001), 2))

  # calibration_ratio = mean(mse) / empirical_mse; empirical_mse = rmse^2.
  expect_equal(res$by_method$calibration_ratio,
               0.0001 / res$by_method$rmse^2, tolerance = 1e-12)
})

test_that("compute_truth_all returns 6 outcomes with correct ratios", {
  fr <- make_synthetic_frame()
  tr <- compute_truth_all(fr, area_var = "municipality")
  expect_setequal(unique(tr$outcome), sae_outcomes$outcome)
  m  <- unique(fr$municipality)[1]
  sub <- fr[fr$municipality == m, ]
  expect_equal(
    tr$theta_true[tr$outcome == "pers_migrated_12m" & tr$area_id == m],
    sum(sub$n_mig_total) / sum(sub$n_pers_resp))
  expect_equal(
    tr$theta_true[tr$outcome == "migrated_last_12m" & tr$area_id == m],
    mean(sub$migrated_last_12m))
})

test_that("evaluate_estimates merges on outcome when both sides carry it", {
  est <- tibble::tibble(
    method    = "fh", replicate = 1L,
    area_id   = c(1L, 1L), outcome = c("a", "b"),
    theta_hat = c(0.10, 0.50),
    mse = 0.01, ci_lo = c(0.05, 0.45), ci_hi = c(0.15, 0.55))
  tr <- tibble::tibble(area_id = c(1L, 1L), outcome = c("a", "b"),
                       theta_true = c(0.10, 0.40))
  ev <- evaluate_estimates(est, tr)
  bm <- ev$by_method_set
  expect_equal(nrow(bm), 2L)               # one row per outcome
  expect_true("outcome" %in% colnames(bm))
  expect_equal(bm$bias[bm$outcome == "a"], 0)
  expect_equal(bm$bias[bm$outcome == "b"], 0.1)
})

test_that("evaluate_estimates errors when truth is outcome-long but estimates have no outcome", {
  est <- tibble::tibble(
    method = "fh", replicate = 1L,
    area_id = c(1L, 2L),
    theta_hat = c(0.1, 0.2),
    mse = 0.01, ci_lo = c(0.05, 0.15), ci_hi = c(0.15, 0.25))
  tr <- tibble::tibble(
    area_id   = c(1L, 1L, 2L, 2L),
    outcome   = c("a", "b", "a", "b"),
    theta_true = c(0.1, 0.4, 0.2, 0.5))
  expect_error(evaluate_estimates(est, tr), "outcome")
})

test_that("evaluate_estimates stops when estimates are outcome-long but truth is not", {
  est <- tibble::tibble(
    method = "fh", replicate = 1L,
    area_id = c(1L, 1L), outcome = c("a", "b"),
    theta_hat = c(0.10, 0.50),
    mse = 0.01, ci_lo = c(0.05, 0.45), ci_hi = c(0.15, 0.55))
  tr <- tibble::tibble(area_id = 1L, theta_true = 0.10)
  expect_error(evaluate_estimates(est, tr), "truth carries no")
})
