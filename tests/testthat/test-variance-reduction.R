test_that("build_direct_pseudo_method constructs a method=direct panel with Wald CIs", {
  dp <- tibble::tibble(
    area_id = c(1L, 2L), outcome = "a", replicate = c(1L, 1L),
    theta_direct = c(0.10, 0.20),
    var_direct_smoothed = c(0.0001, 0.0004)
  )
  out <- build_direct_pseudo_method(dp)
  expect_equal(unique(out$method), "direct")
  expect_equal(out$theta_hat, dp$theta_direct)
  expect_equal(out$mse, dp$var_direct_smoothed)
  expect_equal(out$ci_lo, dp$theta_direct - 1.96 * sqrt(dp$var_direct_smoothed),
               tolerance = 1e-12)
  expect_equal(out$ci_hi, dp$theta_direct + 1.96 * sqrt(dp$var_direct_smoothed),
               tolerance = 1e-12)
})

test_that("compute_direct_mse_ratio computes paired ratios against direct variance", {
  results <- tibble::tibble(
    method = "fh", replicate = c(1L, 1L), outcome = "a",
    area_id = c(1L, 2L), mse = c(0.00005, 0.0002)
  )
  direct <- tibble::tibble(
    replicate = c(1L, 1L), outcome = "a", area_id = c(1L, 2L),
    var_direct = c(0.0001, 0.0002), var_direct_smoothed = c(0.00008, 0.00025),
    n_eff = c(40L, 60L)
  )
  out <- compute_direct_mse_ratio(results, direct)
  expect_equal(nrow(out), 2L)
  expect_equal(out$ratio_raw, c(0.00005 / 0.0001, 0.0002 / 0.0002))
  expect_equal(out$ratio_smoothed, c(0.00005 / 0.00008, 0.0002 / 0.00025))
})

test_that("summarise_mse_ratio reports median ratio and pct below 1 per group", {
  ratio_panel <- tibble::tibble(
    method = "fh", outcome = "a",
    area_id = 1:8, replicate = 1L,
    n_eff = c(10, 20, 30, 40, 50, 60, 70, 80),
    ratio_raw      = c(0.5, 0.6, 0.7, 0.8, 1.2, 1.3, 1.4, 1.5),
    ratio_smoothed = c(0.5, 0.6, 0.7, 0.8, 1.2, 1.3, 1.4, 1.5) * 0.9
  )
  out <- summarise_mse_ratio(ratio_panel)
  expect_true(all(c("method", "outcome", "n_eff_quartile", "median_ratio_raw",
                    "median_ratio_smoothed", "pct_below_1_raw", "n") %in%
                  colnames(out)))
  expect_equal(sum(out$n), 8L)
  low_q <- out[out$n_eff_quartile %in% c(1L, 2L), ]
  expect_true(all(low_q$pct_below_1_raw == 1))
})

test_that("summarise_error_reduction: per-area median ratio, positive-denom only, no draw filter", {
  # area 1: model closer than direct in both reps -> ratio 0.1 (win)
  # area 2: model worse than direct              -> ratio 5   (loss)
  # area 3: direct exact in every rep -> rmse_direct = 0 -> dropped (undefined)
  results <- tibble::tibble(
    method = "fh", outcome = "a",
    replicate = c(1L, 2L, 1L, 2L, 1L, 2L),
    area_id   = c(1L, 1L, 2L, 2L, 3L, 3L),
    theta_hat = c(0.11, 0.09, 0.30, 0.30, 0.50, 0.50)
  )
  direct <- tibble::tibble(
    outcome = "a",
    replicate = c(1L, 2L, 1L, 2L, 1L, 2L),
    area_id   = c(1L, 1L, 2L, 2L, 3L, 3L),
    theta_direct = c(0.20, 0.00, 0.22, 0.18, 0.40, 0.40)
  )
  truth <- tibble::tibble(area_id = c(1L, 2L, 3L), outcome = "a",
                          theta_true = c(0.10, 0.20, 0.40))
  out <- summarise_error_reduction(results, direct, truth)
  expect_equal(out$n_area, 3L)
  expect_equal(out$n_zero_direct, 1L)                 # area 3 dropped
  expect_equal(out$median_ratio, stats::median(c(0.1, 5)))
  expect_equal(out$mean_ratio, mean(c(0.1, 5)))
  expect_equal(out$pct_win, 0.5)                       # area 1 wins, area 2 loses
})

test_that("compute_direct_mse_ratio sets ratio_raw to NA when var_direct is zero", {
  results <- tibble::tibble(
    method = "fh", replicate = c(1L, 1L), outcome = "a",
    area_id = c(1L, 2L), mse = c(0.00005, 0.0002)
  )
  # Row 1: normal case with non-zero var_direct
  # Row 2: zero var_direct case (common for rare outcomes)
  direct <- tibble::tibble(
    replicate = c(1L, 1L), outcome = "a", area_id = c(1L, 2L),
    var_direct = c(0.0001, 0), var_direct_smoothed = c(0.00008, 0.00025),
    n_eff = c(40L, 60L)
  )
  out <- compute_direct_mse_ratio(results, direct)
  expect_equal(nrow(out), 2L)
  # Row 1: normal case, ratio_raw should be computed
  expect_equal(out$ratio_raw[1], 0.00005 / 0.0001)
  expect_true(is.na(out$ratio_raw[2]))
  # ratio_smoothed should always be computed (never NA even when var_direct == 0)
  expect_equal(out$ratio_smoothed[1], 0.00005 / 0.00008)
  expect_equal(out$ratio_smoothed[2], 0.0002 / 0.00025)
  expect_false(is.na(out$ratio_smoothed[2]))
})
