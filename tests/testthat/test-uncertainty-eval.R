library(data.table)

.ur_rows <- function() {
  # two areas x four replicates, one method; truth 0.10 and 0.00
  d <- data.table(method = "m", config = "c", outcome = "o",
                  area_set = "in_sample", area_id = rep(1:2, each = 4L),
                  replicate = rep(1:4, 2L),
                  theta_hat = c(0.12, 0.08, 0.10, 0.14, 0, 0, 0, 0),
                  mse = c(4e-4, 4e-4, 4e-4, 4e-4, 1e-6, 1e-6, 1e-6, 1e-6),
                  ci_lo = c(0.09, 0.05, 0.07, 0.11, 0, 0, 0, 0),
                  ci_hi = c(0.15, 0.11, 0.13, 0.17, 0.01, 0.01, 0.01, 0.01),
                  theta_true = rep(c(0.10, 0.00), each = 4L))
  d
}

test_that("grouping keys carry no version column and the bootstrap seed is fixed", {
  expect_identical(UR_KEYS, c("method", "config", "outcome", "area_set"))
  expect_identical(formals(ur_mcse)$seed, 20261009L)
})

test_that("interval score is width inside and adds 2/alpha per unit of miss", {
  expect_equal(ur_interval_score(0.1, 0.3, 0.2), 0.2)
  expect_equal(ur_interval_score(0.1, 0.3, 0.05), 0.2 + 40 * 0.05)
  expect_equal(ur_interval_score(0.1, 0.3, 0.4), 0.2 + 40 * 0.1)
})

test_that("bad fits are dropped and counted, by reason", {
  d <- rbind(.ur_rows()[, `:=`(fit_failed = FALSE, max_rhat = 1.01)],
             .ur_rows()[, `:=`(method = "f", fit_failed = TRUE, max_rhat = NA_real_)],
             .ur_rows()[, `:=`(method = "r", fit_failed = FALSE, max_rhat = 1.30)])
  out <- ur_drop_bad_fits(d)
  expect_setequal(unique(out$d$method), "m")
  expect_identical(nrow(out$dropped), 8L)           # 4 replicates x 2 methods
  expect_setequal(out$dropped$reason, c("failed", "rhat"))
  d2 <- copy(d)[method == "m"][1, theta_hat := NA_real_]
  expect_identical(ur_drop_bad_fits(d2)$rows_nonfinite, 1L)
})

test_that("row and area metrics are right, and a zero realised MSE gives NA ratios", {
  d <- ur_add_row_metrics(.ur_rows())
  expect_equal(d$cov, c(1, 1, 1, 0, 1, 1, 1, 1))
  a <- ur_area_summary(d)
  a1 <- a[area_id == 1L]
  expect_equal(a1$emp_mse, mean(c(0.02, -0.02, 0, 0.04)^2))
  expect_equal(a1$ratio, 4e-4 / a1$emp_mse)
  expect_equal(a1$coverage, 0.75)
  a2 <- a[area_id == 2L]
  expect_identical(a2$emp_mse, 0)
  expect_true(is.na(a2$ratio) && is.na(a2$rrmse_mse) && is.na(a2$width_rel))
})

test_that("set summary ignores undefined ratios and areas below min_reps", {
  d <- ur_add_row_metrics(.ur_rows())
  a <- ur_area_summary(d)
  s <- ur_set_summary(d, a, min_reps = 4L)
  expect_identical(s$n_ratio_undefined, 1L)
  expect_equal(s$ratio_med, a[area_id == 1L, ratio])
  expect_equal(s$coverage, 7 / 8)
  expect_true(is.finite(s$ratio_pooled))
  s5 <- ur_set_summary(d, a, min_reps = 5L)
  expect_identical(s5$n_areas, 0L)
  expect_true(is.na(s5$ratio_med))
})

test_that("reliability bins: reported equals realised when mse equals the squared error", {
  set.seed(1)
  e <- rnorm(2000, 0, rep(c(0.01, 0.05), each = 1000))
  d <- data.table(method = "m", config = "c", outcome = "o",
                  area_set = "in_sample", area_id = 1L, replicate = seq_along(e),
                  theta_hat = 0.1 + e, theta_true = 0.1, mse = e^2,
                  ci_lo = 0, ci_hi = 1)
  r <- ur_reliability(ur_add_row_metrics(d), n_bins = 10L)
  expect_identical(nrow(r), 10L)
  expect_equal(r$rmse_reported, r$rmse_realised)
})

test_that("size bins put out-of-sample first and cover every in-sample row", {
  br <- ur_size_breaks(c(1, 5, 10, 20, 40, 80, 160, 320, 640, 1280), n_bins = 5L)
  b <- ur_size_bin(c(0, 1, 30, 1280), c("out_of_sample", rep("in_sample", 3)), br)
  expect_identical(levels(b)[1], "0 (out of sample)")
  expect_identical(as.character(b[1]), "0 (out of sample)")
  expect_false(anyNA(b))
})

test_that("CV flags: zero truth is never truly publishable, zero estimate never reported publishable", {
  d <- ur_add_row_metrics(.ur_rows())
  a <- ur_area_summary(d)
  f <- ur_cv_flags(d, a, threshold = 0.20)
  # area 1: reported CV 0.02/0.08..0.14 -> publishable only where theta_hat > 0.1;
  # realised CV sqrt(emp)/0.10 = 0.245 -> never; area 2: theta 0 -> neither
  expect_equal(f$share_pub_true, 0)
  expect_equal(f$false_publish, 1)
  expect_true(is.na(f$missed_publish))
  expect_equal(f$share_pub_reported, mean(c(0.02 / c(0.12, 0.08, 0.10, 0.14) < 0.2, rep(FALSE, 4))))
})

test_that("MCSE is zero when every replicate is identical", {
  d <- ur_add_row_metrics(.ur_rows()[area_id == 1L][, `:=`(theta_hat = 0.12, ci_lo = 0.09, ci_hi = 0.15)])
  m <- ur_mcse(ur_rep_sums(d), B = 50L)
  expect_equal(m$mcse_coverage, 0)
  expect_equal(m$mcse_ratio_pooled, 0)
})

test_that("the set gate catches methods that drew different municipalities", {
  a <- .ur_rows(); b <- copy(a)[, method := "n"]
  expect_invisible(ur_assert_same_sets(rbind(a, b)))
  b[area_id == 2L, area_set := "out_of_sample"]
  expect_error(ur_assert_same_sets(rbind(a, b)), "in-sample area sets differ")
})

test_that("realised MSE splits into squared bias plus variance", {
  d <- ur_add_row_metrics(.ur_rows())
  a <- ur_area_summary(d)
  a1 <- a[area_id == 1L]          # errors 0.02, -0.02, 0, 0.04
  expect_equal(a1$bias, 0.01)
  expect_equal(a1$emp_var, 5e-4)
  expect_equal(a1$emp_mse, a1$bias^2 + a1$emp_var)
  expect_equal(a1$ratio_var, 4e-4 / 5e-4)
  expect_equal(a1$bias_share, 1e-4 / 6e-4)
  a2 <- a[area_id == 2L]          # estimated exactly every time
  expect_true(is.na(a2$ratio_var) && is.na(a2$bias_share))
  # constant but wrong: no variance to compare against, all of the MSE is bias
  d3 <- ur_add_row_metrics(.ur_rows()[area_id == 1L][, theta_hat := 0.13])
  a3 <- ur_area_summary(d3)
  expect_identical(a3$emp_var, 0)
  expect_true(is.na(a3$ratio_var))
  expect_equal(a3$bias_share, 1)
  s <- ur_set_summary(d, a, min_reps = 4L)
  expect_equal(s$ratio_var_med, 0.8)
  expect_equal(s$bias_share_med, 1 / 6)
})

test_that("per-bin empirical and estimated RMSE, mean truth and the bin's median size", {
  d <- ur_add_row_metrics(.ur_rows())
  d[, `:=`(n_dwell = rep(c(40L, 100L), each = 4L), bin = rep(c("a", "b"), each = 4L))]
  r <- ur_by_bin(d, "bin", by = "method", x = "n_dwell")[order(bin)]
  # area 1 errors 0.02, -0.02, 0, 0.04 -> mean square 6e-4; its reported MSE is 4e-4
  expect_equal(r$rmse_empirical, c(sqrt(6e-4), 0))
  expect_equal(r$rmse_estimated, c(sqrt(4e-4), sqrt(1e-6)))
  expect_equal(r$truth_mean, c(0.10, 0))
  expect_equal(r$x_med, c(40, 100))
  # without `x` (the prevalence bins) there is no size column
  expect_false("x_med" %in% names(ur_by_bin(d, "bin", by = "method")))
})

test_that("relative bias of the estimated RMSE uses the mean of the root MSE, as the literature does", {
  d <- ur_add_row_metrics(.ur_rows())
  d[area_id == 1L, mse := c(1e-4, 4e-4, 9e-4, 4e-4)]   # mean root 0.02; root of the mean 0.0212
  a <- ur_area_summary(d)
  expect_equal(a[area_id == 1L, rep_rmse], 0.02)
  expect_equal(a[area_id == 1L, rb_rmse], 0.02 / sqrt(6e-4) - 1)
  expect_true(is.na(a[area_id == 2L, rb_rmse]))          # empirical MSE 0: undefined, not Inf
  s <- ur_set_summary(d, a, min_reps = 1L)
  expect_equal(s$rb_rmse_med, 0.02 / sqrt(6e-4) - 1)
  expect_equal(s$rb_rmse_mean, 0.02 / sqrt(6e-4) - 1)
})

test_that("area summary carries each municipality's median dwellings drawn, when the rows have them", {
  d <- ur_add_row_metrics(.ur_rows())
  # without n_dwell there is no size column, and no error
  expect_false("n_dwell_med" %in% names(ur_area_summary(d)))
  d[, n_dwell := c(40L, 50L, 60L, 70L, 100L, 100L, 120L, 100L)]
  a <- ur_area_summary(d)
  expect_equal(a[area_id == 1L, n_dwell_med], 55)    # even count: a double, not an integer
  expect_equal(a[area_id == 2L, n_dwell_med], 100)
  expect_type(a$n_dwell_med, "double")
})
