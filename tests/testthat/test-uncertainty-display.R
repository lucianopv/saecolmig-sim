library(data.table)

# One method per row block: area-level summary rows as ur_area_summary() returns them
.disp_area <- function() {
  data.table(
    method = rep(c("a", "b"), each = 6L), config = "c", outcome = "o",
    area_set = rep(rep(c("in_sample", "out_of_sample"), each = 3L), 2L),
    area_id = rep(1:3, 4L),
    n_reps = rep(c(10L, 20L, 30L), 4L),
    rep_rmse = c(0.02, 0.04, 0.03, 0.05, 0.01, 0.02,
                 0.03, 0.01, 0.06, 0.02, 0.04, 0.05),
    emp_mse = c(4e-4, 9e-4, 0, 1e-4, 4e-4, 9e-4,
                9e-4, 1e-4, 4e-4, 4e-4, 4e-4, 9e-4),
    coverage = c(0.9, 0.96, 0.5, 1, 0.8, 0.95, 0.99, 0.7, 0.95, 1, 1, 0.2),
    n_dwell_med = rep(c(30, 40, 50), 4L))
}

test_that("size groups pool sparse sizes upward and the tail joins the last group", {
  n <- c(rep(10, 30), rep(11, 5), rep(12, 5), rep(13, 20), rep(14, 3))
  g <- ur_size_groups(n, min_n = 25L)
  expect_equal(as.vector(tapply(g, n, unique)), c(1L, 2L, 2L, 2L, 2L))
  expect_true(all(table(g) >= 25L))
})

test_that("coverage summary weights municipalities by replicates", {
  a <- data.table::data.table(method = "m", config = "c", outcome = "o",
    area_set = "in_sample", area_id = 1:2, n_reps = c(10L, 30L), coverage = c(1, 0.5))
  s <- ur_coverage_summary(a, "o")
  expect_equal(s$overall, (10 * 1 + 30 * 0.5) / 40)
  expect_equal(s$share_ge_095, 0.5); expect_equal(s$share_lt_080, 0.5)
})

test_that("ratio summary ignores municipalities with zero empirical MSE", {
  a <- data.table::data.table(method = "m", config = "c", outcome = "o",
    area_set = "in_sample", area_id = 1:3, rep_rmse = c(1, 2, 5), emp_mse = c(1, 1, 0))
  s <- ur_ratio_summary(a, "o")
  expect_equal(s$med, 1.5); expect_equal(s$below, 0)
})

test_that("ratio summary has the quartiles of the ratio and the share under one, per method and set", {
  a <- .disp_area()
  s <- ur_ratio_summary(a, "o")
  expect_named(s, c("method", "area_set", "p25", "med", "p75", "below"))
  expect_equal(nrow(s), 4L)
  # method a, in sample: area 3 has zero empirical MSE and drops out
  r <- c(0.02 / 0.02, 0.04 / 0.03)
  x <- s[method == "a" & area_set == "in_sample"]
  expect_equal(x$p25, unname(quantile(r, 0.25)))
  expect_equal(x$med, median(r))
  expect_equal(x$p75, unname(quantile(r, 0.75)))
  expect_equal(x$below, mean(r < 1))
  # method b, out of sample: ratios 0.02/0.02, 0.04/0.02, 0.05/0.03
  r <- c(1, 2, 0.05 / 0.03)
  x <- s[method == "b" & area_set == "out_of_sample"]
  expect_equal(x$med, median(r)); expect_equal(x$below, 0)
})

test_that("the summaries read one outcome and ignore the others", {
  a <- .disp_area()
  b <- copy(a)[, outcome := "p"][, `:=`(rep_rmse = 9, emp_mse = 1, coverage = 0, n_reps = 1L)]
  ab <- rbind(a, b)
  expect_equal(ur_ratio_summary(ab, "o"), ur_ratio_summary(a, "o"))
  expect_equal(ur_tracking_summary(ab, "o"), ur_tracking_summary(a, "o"))
  expect_equal(ur_size_medians(ab, "o", min_n = 3L), ur_size_medians(a, "o", min_n = 3L))
  expect_equal(ur_coverage_summary(ab, "o"), ur_coverage_summary(a, "o"))
  # an outcome called like the argument name must not select every row
  a2 <- copy(a)[, outcome := rep(c("o", "outcome"), length.out = .N)]
  expect_equal(unique(ur_ratio_summary(a2, "outcome")$method), c("a", "b"))
  expect_equal(ur_coverage_summary(a2, "outcome")$overall,
               a2[outcome == "outcome" & area_set == "in_sample", weighted.mean(coverage, n_reps), by = method]$V1)
})

test_that("the summaries need no version column and tolerate a column column", {
  a <- .disp_area()[, column := "geo"]
  expect_false("version" %in% names(a))
  expect_s3_class(ur_ratio_summary(a, "o"), "data.table")
  expect_s3_class(ur_tracking_summary(a, "o"), "data.table")
  expect_s3_class(ur_size_medians(a, "o", min_n = 2L), "data.table")
  expect_s3_class(ur_coverage_summary(a, "o"), "data.table")
  expect_false("column" %in% names(ur_coverage_summary(a, "o")))
})

test_that("tracking summary keeps every municipality: rank correlation and share under the line", {
  a <- .disp_area()
  s <- ur_tracking_summary(a, "o")
  expect_named(s, c("method", "area_set", "spearman", "below"))
  x <- a[method == "a" & area_set == "in_sample"]       # includes the zero-MSE area
  expect_equal(s[method == "a" & area_set == "in_sample", spearman],
               cor(x$rep_rmse, sqrt(x$emp_mse), method = "spearman"))
  expect_equal(s[method == "a" & area_set == "in_sample", below],
               mean(x$rep_rmse < sqrt(x$emp_mse)))
  # a perfectly tracking method: rank correlation 1, nothing strictly below
  tr <- data.table(method = "m", area_set = "in_sample", outcome = "o",
                   rep_rmse = c(1, 2, 3), emp_mse = c(1, 4, 9))
  st <- ur_tracking_summary(tr, "o")
  expect_equal(st$spearman, 1); expect_equal(st$below, 0)
})

test_that("size groups: every distinct size stays whole, equal sizes share a group", {
  n <- c(rep(5, 30), rep(6, 25), rep(7, 40))
  expect_identical(ur_size_groups(n, min_n = 25L), rep(1:3, c(30L, 25L, 40L)))
  expect_type(ur_size_groups(n), "integer")
  # sizes arrive unsorted: the groups follow the sorted distinct values
  g <- ur_size_groups(rev(n), min_n = 25L)
  expect_identical(g, rev(rep(1:3, c(30L, 25L, 40L))))
  # a group that cannot reach the minimum is an error, not a silent short group
  expect_error(ur_size_groups(1:5, min_n = 25L), "fewer than")
})

test_that("size medians pool sizes and give the median size and RMSE per group, empirical and estimated", {
  n <- c(rep(10, 3), rep(11, 1), rep(12, 1), rep(13, 2))
  a <- data.table(method = "m", config = "c", outcome = "o", area_set = "in_sample",
                  area_id = seq_along(n), n_dwell_med = n,
                  emp_mse = (seq_along(n) / 100)^2, rep_rmse = seq_along(n) / 50)
  s <- ur_size_medians(a, "o", min_n = 3L)
  expect_named(s, c("method", "kind", "size_group", "n_dwell", "rmse"))
  # groups: {10} (3 areas), {11, 12, 13} (4 areas)
  e <- s[kind == "Empirical RMSE"][order(size_group)]
  expect_equal(e$size_group, 1:2)
  expect_equal(e$n_dwell, c(10, median(c(11, 12, 13, 13))))
  expect_equal(e$rmse, c(median((1:3) / 100), median((4:7) / 100)))
  m <- s[kind == "Estimated RMSE"][order(size_group)]
  expect_equal(m$rmse, c(median((1:3) / 50), median((4:7) / 50)))
  expect_setequal(s$kind, c("Empirical RMSE", "Estimated RMSE"))
})

test_that("size medians use in-sample rows only", {
  a <- .disp_area()
  a[area_set == "out_of_sample", n_dwell_med := 0]
  s <- ur_size_medians(a, "o", min_n = 3L)
  expect_setequal(s$method, c("a", "b"))
  expect_true(all(s$n_dwell > 0))
  expect_equal(s[method == "b" & kind == "Estimated RMSE", rmse], median(c(0.03, 0.01, 0.06)))
})

test_that("sizes the first method does not carry are left out of the medians, with a warning", {
  a <- .disp_area()[area_set == "in_sample"]
  a[method == "b" & area_id == 3L, n_dwell_med := 60]
  expect_warning(s <- ur_size_medians(a, "o", min_n = 3L), "not carried by")
  expect_equal(s[method == "b" & kind == "Estimated RMSE", rmse], median(c(0.03, 0.01)))
  expect_equal(s[method == "a" & kind == "Estimated RMSE", rmse], median(c(0.02, 0.04, 0.03)))
})

test_that("coverage summary is in sample only, per method, with thresholds 0.95 and 0.80", {
  a <- .disp_area()
  s <- ur_coverage_summary(a, "o")
  expect_named(s, c("method", "overall", "share_ge_095", "share_lt_080"))
  x <- a[method == "a" & area_set == "in_sample"]
  expect_equal(s[method == "a", overall], weighted.mean(x$coverage, x$n_reps))
  expect_equal(s[method == "a", share_ge_095], mean(x$coverage >= 0.95))
  expect_equal(s[method == "a", share_lt_080], mean(x$coverage < 0.8))
  # a coverage of exactly 0.95 counts as at least 0.95, exactly 0.80 does not count as below
  tc <- data.table(method = "m", outcome = "o", area_set = "in_sample", n_reps = 1L,
                   coverage = c(0.95, 0.8))
  st <- ur_coverage_summary(tc, "o")
  expect_equal(st$share_ge_095, 0.5); expect_equal(st$share_lt_080, 0)
})
