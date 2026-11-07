test_that("aggregate_theta_to_dept computes the population-weighted mean per dept", {
  area_id   <- c(11001L, 11002L, 25001L, 25002L, 25003L)
  theta_hat <- c(0.10,   0.20,   0.05,   0.15,   0.25)
  n_persons <- c(1000,   3000,   500,    500,    1000)
  out <- aggregate_theta_to_dept(area_id, theta_hat, n_persons)

  expect_setequal(out$dept, c(11L, 25L))
  expected_11 <- (0.10 * 1000 + 0.20 * 3000) / (1000 + 3000)
  expected_25 <- (0.05 * 500 + 0.15 * 500 + 0.25 * 1000) / (500 + 500 + 1000)
  expect_equal(out$theta_dept[out$dept == 11L], expected_11, tolerance = 1e-12)
  expect_equal(out$theta_dept[out$dept == 25L], expected_25, tolerance = 1e-12)
  expect_equal(out$n_munis[out$dept == 11L], 2L)
  expect_equal(out$n_munis[out$dept == 25L], 3L)
})

test_that("combine_in_sample_oos unions and tags area_set with no overlap", {
  in_s <- tibble::tibble(area_id = c(1L, 2L), outcome = "a", replicate = 1L,
                         theta_hat = c(0.1, 0.2), mse = c(0.001, 0.002))
  oos  <- tibble::tibble(area_id = c(3L, 4L), outcome = "a", replicate = 1L,
                         theta_hat = c(0.3, 0.4), mse = c(0.003, 0.004))
  out <- combine_in_sample_oos(in_s, oos)
  expect_equal(nrow(out), 4L)
  expect_setequal(out$area_set, c("in_sample", "out_of_sample"))
  expect_equal(sort(out$area_id), 1:4)
})

test_that("combine_in_sample_oos errors on overlapping area_id/outcome/replicate", {
  in_s <- tibble::tibble(area_id = 1L, outcome = "a", replicate = 1L,
                         theta_hat = 0.1, mse = 0.001)
  oos  <- tibble::tibble(area_id = 1L, outcome = "a", replicate = 1L,
                         theta_hat = 0.2, mse = 0.002)
  expect_error(combine_in_sample_oos(in_s, oos), "overlapping")
})

test_that("summarise_dept_comparison recovers zero bias when model equals truth", {
  dept_truth <- tibble::tibble(dept = c(11L, 25L), outcome = "a",
                               theta_true = c(0.10, 0.20))
  dept_direct <- tibble::tibble(
    dept = rep(c(11L, 25L), each = 3), outcome = "a",
    replicate = rep(1:3, 2),
    estimate = c(0.11, 0.09, 0.10,  0.22, 0.18, 0.20)
  )
  dept_model <- tibble::tibble(
    dept = rep(c(11L, 25L), each = 3), outcome = "a",
    replicate = rep(1:3, 2),
    theta_dept = c(0.10, 0.10, 0.10,  0.20, 0.20, 0.20)
  )
  out <- summarise_dept_comparison(dept_direct, dept_model, dept_truth)
  expect_equal(out$bias_model, c(0, 0), tolerance = 1e-12)
  expect_true(all(abs(out$bias_direct) < 0.02))
  expect_equal(out$median_gap[out$dept == 11L],
               stats::median(c(0.10, 0.10, 0.10) - c(0.11, 0.09, 0.10)),
               tolerance = 1e-12)
})
