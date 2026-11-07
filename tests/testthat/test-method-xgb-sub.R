# fit_xgb_sub (R/method_xgb_sub.R)

# Small enough to fit in seconds: few boosting rounds and bootstrap replicates.
.xgb_sub_fast <- list(nrounds = 15L)

test_that("fit_xgb_sub satisfies the contract on the fixture", {
  skip_if_not_installed("povmap")
  f <- .fixture_sub()
  fit <- .quiet_fit(fit_xgb_sub(f$direct, f$smp, f$pop, "migrated_last_12m",
                                B = 5L, seed = 1L, params = .xgb_sub_fast))
  expect_s3_class(fit, c("sae_fit", "sae_fit_xgb_sub"))
  expect_false(any(fit$diagnostics$metric == "xgb_sub_failed"))
  expect_silent(sae_fit_contract(fit))
  expect_identical(fit$meta$method, "xgb_sub")
  expect_equal(nrow(fit$estimates), f$n_muni)
  expect_identical(sort(fit$estimates$area_id), seq_len(f$n_muni))
  expect_identical(sum(fit$estimates$area_set == "in_sample"), length(f$drawn))
  expect_identical(sum(fit$estimates$area_set == "out_of_sample"),
                   f$n_muni - length(f$drawn))
  expect_true(all(fit$estimates$area_id[fit$estimates$area_set == "in_sample"]
                  %in% f$drawn))
})

test_that("fit_xgb_sub records what it trained on and how", {
  skip_if_not_installed("povmap")
  f <- .fixture_sub()
  fit <- .quiet_fit(fit_xgb_sub(f$direct, f$smp, f$pop, "migrated_last_12m",
                                B = 5L, seed = 7L,
                                params = list(nrounds = 15L, eta = 0.2)))
  d <- stats::setNames(fit$diagnostics$value, fit$diagnostics$metric)
  expect_identical(unname(d["B"]), 5)
  expect_identical(unname(d["n_subarea_train"]), as.numeric(nrow(f$smp)))
  expect_identical(unname(d["n_subarea_pop"]), as.numeric(nrow(f$pop)))
  expect_identical(unname(d["n_covariates"]), 3)
  expect_identical(unname(d["level_municipality"]), 0)
  expect_identical(unname(d["nrounds"]), 15)
  expect_identical(unname(d["eta"]), 0.2)
  expect_identical(unname(d["max_depth"]), 4)
  expect_identical(fit$meta$seed, 7L)
  expect_identical(fit$meta$xgboost_version,
                   as.character(utils::packageVersion("xgboost")))
  expect_identical(fit$meta$povmap_version,
                   as.character(utils::packageVersion("povmap")))
})

test_that("fit_xgb_sub is reproducible for a given seed", {
  skip_if_not_installed("povmap")
  f <- .fixture_sub()
  run <- function(seed) .quiet_fit(
    fit_xgb_sub(f$direct, f$smp, f$pop, "migrated_last_12m",
                B = 5L, seed = seed, params = .xgb_sub_fast))
  a <- run(3L); b <- run(3L)
  expect_identical(a$estimates$theta_hat, b$estimates$theta_hat)
  expect_identical(a$estimates$mse, b$estimates$mse)
})

test_that("fit_xgb_sub declares the two-level structure and the bootstrap to povmap", {
  skip_if_not_installed("povmap")
  f <- .fixture_sub()
  seen <- NULL
  local_mocked_bindings(
    xgb = function(...) { seen <<- list(...); stop("captured") },
    .package = "povmap")
  fit_xgb_sub(f$direct, f$smp, f$pop, "migrated_last_12m", B = 7L, seed = 5L,
              cpus = 2L, params = list(max_depth = 3L, lambda = 2))
  expect_identical(seen$domains, "municipality")
  expect_identical(seen$sub_domains, "code")
  expect_identical(seen$transformation, "no")
  expect_true(seen$bootstrap)
  expect_identical(seen$bootstrap_type, "residual")
  expect_identical(seen$B, 7L)
  expect_identical(seen$seed, 5L)
  expect_identical(seen$cpus, 2L)
  expect_identical(seen$smp_weights, "wt")
  expect_identical(seen$pop_weights, "n_households_census")
  expect_identical(deparse(seen$fixed), "y ~ z1 + z2 + z3")
  # defaults fill what the caller leaves out
  expect_identical(seen$nrounds, 200L)
  expect_identical(seen$eta, 0.1)
  expect_identical(seen$max_depth, 3L)
  expect_identical(seen$min_child_weight, 1)
  expect_identical(seen$subsample, 0.8)
  expect_identical(seen$colsample_bytree, 0.8)
  expect_identical(seen$lambda, 2)
  expect_identical(seen$alpha, 0)
})

test_that("fit_xgb_sub passes no seed to povmap when none is given", {
  skip_if_not_installed("povmap")
  f <- .fixture_sub()
  seen <- NULL
  local_mocked_bindings(
    xgb = function(...) { seen <<- list(...); stop("captured") },
    .package = "povmap")
  fit_xgb_sub(f$direct, f$smp, f$pop, "migrated_last_12m", B = 2L)
  expect_true("seed" %in% names(seen))
  expect_null(seen$seed)
})

test_that("a povmap error is recorded as NA rows for every municipality, not thrown", {
  skip_if_not_installed("povmap")
  f <- .fixture_sub()
  local_mocked_bindings(
    xgb = function(...) stop("no usable boosting round"),
    .package = "povmap")
  fit <- fit_xgb_sub(f$direct, f$smp, f$pop, "migrated_last_12m", B = 2L, seed = 1L)
  expect_s3_class(fit, "sae_fit")
  expect_equal(nrow(fit$estimates), f$n_muni)
  expect_true(all(is.na(fit$estimates$theta_hat)))
  expect_true(all(is.na(fit$estimates$mse)))
  d <- fit$diagnostics
  expect_identical(d$value[d$metric == "xgb_sub_failed"], 1)
  expect_identical(d$character_value[d$metric == "xgb_sub_error_message"],
                   "no usable boosting round")
  expect_identical(fit$meta$method, "xgb_sub")
  expect_error(sae_fit_contract(fit), "NA mse")
})

test_that("fit_xgb_sub at municipality level fits the same estimator on rolled-up inputs", {
  skip_if_not_installed("povmap")
  f <- .fixture_sub(n_muni = 12L, per_muni = 4L)
  fit <- .quiet_fit(
    fit_xgb_sub(f$direct, f$smp, f$pop, "migrated_last_12m", B = 5L, seed = 1L,
                params = .xgb_sub_fast, level = "municipality"))
  expect_false(any(fit$diagnostics$metric == "xgb_sub_failed"))
  expect_silent(sae_fit_contract(fit))
  d <- stats::setNames(fit$diagnostics$value, fit$diagnostics$metric)
  expect_identical(unname(d["level_municipality"]), 1)
  expect_identical(unname(d["n_subarea_pop"]), 12)
  expect_identical(unname(d["n_subarea_train"]), 6)
})

test_that("fit_xgb_sub refuses an unknown level or weight", {
  skip_if_not_installed("povmap")
  f <- .fixture_sub()
  expect_error(fit_xgb_sub(f$direct, f$smp, f$pop, "migrated_last_12m",
                           level = "block"), "'arg'")
  expect_error(fit_xgb_sub(f$direct, f$smp, f$pop, "migrated_last_12m",
                           weight = "other"), "'arg'")
})
