test_that("safe_method_fit returns contract estimates + flags on success", {
  fx <- .fixture_robust(n_in = 40L, n_oos = 5L, seed = 11L)
  # sigma2_v at the limit of emdi2's search interval on this fixture: a warning,
  # not a failure.
  expect_warning(
    est <- safe_method_fit("fh", list(direct = fx$direct, X = fx$X_in,
                                      X_oos = fx$X_oos),
                           area_ids = fx$direct$area_id),
    "interval limit")
  expect_true(all(c("area_id", "theta_hat", "mse", "ci_lo", "ci_hi",
                    "fit_failed", "error_msg") %in% colnames(est)))
  expect_true(all(!est$fit_failed))
  expect_true(all(is.na(est$error_msg)))
})

test_that("safe_method_fit converts errors into NA rows with fit_failed", {
  bad_args <- list(direct = data.frame(nonsense = 1), X = data.frame())
  est <- safe_method_fit("fh", bad_args, area_ids = c(101L, 102L, 103L))
  expect_equal(nrow(est), 3L)
  expect_true(all(est$fit_failed))
  expect_true(all(is.na(est$theta_hat)))
  expect_true(all(nzchar(est$error_msg)))
  expect_equal(est$area_id, c(101L, 102L, 103L))
})

test_that("safe_method_fit optionally attaches the fit's diagnostics", {
  skip_if_not_installed("emdi2")
  set.seed(7)
  n <- 30L
  X <- tibble::tibble(municipality = 1:n, x1 = rnorm(n), x2 = rnorm(n))
  th <- pmin(pmax(0.1 + 0.02 * X$x1 + rnorm(n, 0, 0.01), 0.01), 0.5)
  v  <- th * (1 - th) / 80
  direct <- tibble::tibble(area_id = 1:n, theta_direct = th,
                           var_direct = v, var_direct_smoothed = v,
                           n_eff = rep(80, n), n_dwell = rep(80, n),
                           outcome = "migrated_last_12m")
  args <- method_args("fh", direct = direct, aux_in = X, aux_oos = NULL,
                      outcome = "migrated_last_12m")
  # sigma2_v at the limit of emdi2's search interval on this fixture.
  expect_warning(est_plain <- safe_method_fit("fh", args, area_ids = 1:n), "interval limit")
  expect_null(attr(est_plain, "diagnostics"))
  expect_warning(est_diag <- safe_method_fit("fh", args, area_ids = 1:n, diagnostics = TRUE),
                 "interval limit")
  dg <- attr(est_diag, "diagnostics")
  expect_s3_class(dg, "tbl_df")
  expect_true(all(c("metric", "value") %in% colnames(dg)))
  # failure path: unknown covariate column forces an error -> NA rows, NULL attr
  bad <- args; bad$X <- args$X[, "municipality", drop = FALSE]
  est_fail <- safe_method_fit("fh", bad, area_ids = 1:n, diagnostics = TRUE)
  expect_true(all(est_fail$fit_failed))
  expect_null(attr(est_fail, "diagnostics"))
})

test_that("sae_fit_contract rejects NA mse with an NA-specific message", {
  # An NA mse must fail with the contract's own message: `if (!all(mse > 0))`
  # on an NA vector would instead evaluate `if (NA)` and error with "missing
  # value where TRUE/FALSE needed" (the check is in R/zz_method_registry.R).
  fit <- structure(
    list(
      estimates = tibble::tibble(area_id = 1:3, theta_hat = c(0.1, 0.2, NA),
                                 mse = c(0.01, NA, 0.02),
                                 ci_lo = c(0.05, 0.1, 0.1),
                                 ci_hi = c(0.15, 0.3, 0.3)),
      diagnostics = tibble::tibble(metric = "x", value = 1),
      meta = list(method = "fake", runtime_sec = 0.1, n_areas = 3L)
    ),
    class = "sae_fit"
  )
  expect_error(sae_fit_contract(fit), "NA mse")
})

test_that("safe_method_fit records a contract violation as a fit failure, not a crash", {
  # fit_fh_arcsin catches an emdi2::fh error and returns an all-NA sae_fit
  # instead of raising (R/method_fh_arcsin.R). The contract check of that
  # object runs inside safe_method_fit's error handler, so the failure is
  # recorded rather than raised. An NA covariate on one in-sample area makes
  # emdi2::fh fail deterministically.
  skip_if_not_installed("emdi2")
  fx <- .fixture_robust(n_in = 40L, n_oos = 5L, seed = 11L)
  X_in_bad <- fx$X_in
  X_in_bad$x1[1] <- NA_real_
  args <- method_args("fh_arcsin", direct = fx$direct, aux_in = X_in_bad,
                      aux_oos = fx$X_oos, outcome = "migrated_last_12m",
                      seed = 123L)

  # Sanity check: the underlying wrapper really does return a contract-
  # violating (all-NA mse) sae_fit rather than raising an R error.
  raw_fit <- do.call(saecol2026::sae_methods[["fh_arcsin"]], args)
  expect_true(all(is.na(raw_fit$estimates$mse)))

  # safe_method_fit must return recorded failure rows, not raise.
  est <- safe_method_fit("fh_arcsin", args, area_ids = fx$direct$area_id)
  expect_equal(nrow(est), length(fx$direct$area_id))
  expect_true(all(est$fit_failed))
  expect_true(all(is.na(est$theta_hat)))
  expect_true(all(is.na(est$mse)))
  expect_true(all(grepl("NA mse", est$error_msg)))
})

test_that("safe_method_fit keeps the diagnostics of a fit that fails the contract", {
  failed_fit <- function(...) structure(
    list(estimates = tibble::tibble(area_id = 1:2, theta_hat = NA_real_, mse = NA_real_,
                                    ci_lo = NA_real_, ci_hi = NA_real_),
         diagnostics = tibble::tibble(metric = c("m_failed", "m_error_message"),
                                      value = c(1, NA_real_),
                                      character_value = c(NA_character_, "no convergence")),
         meta = list(method = "m", runtime_sec = 0, n_areas = 2L)),
    class = "sae_fit")
  reg <- sae_methods
  reg$fh <- failed_fit
  local_mocked_bindings(sae_methods = reg)
  est <- safe_method_fit("fh", list(), area_ids = 1:2, diagnostics = TRUE)
  expect_true(all(est$fit_failed))
  expect_match(est$error_msg[1], "NA mse")
  expect_identical(attr(est, "diagnostics")$character_value[2], "no convergence")
  expect_null(attr(safe_method_fit("fh", list(), area_ids = 1:2), "diagnostics"))
})
