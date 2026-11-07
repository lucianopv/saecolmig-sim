# fit_megb_sub (R/method_megb_sub.R)

test_that("fit_megb_sub satisfies the contract and records its settings", {
  skip_if_not_installed("povmap")
  # 15 drawn municipalities: megb's domain-wise cross-fitting needs more drawn
  # municipalities than its ten folds.
  f <- .fixture_sub(n_muni = 30L, per_muni = 8L)
  fit <- .quiet_fit(fit_megb_sub(f$direct, f$smp, f$pop, "migrated_last_12m",
                                 B = 5L, seed = 1L))
  expect_s3_class(fit, c("sae_fit", "sae_fit_megb_sub"))
  expect_false(any(fit$diagnostics$metric == "megb_sub_failed"))
  expect_silent(sae_fit_contract(fit))
  expect_identical(fit$meta$method, "megb_sub")
  expect_setequal(unique(fit$estimates$area_set), c("in_sample", "out_of_sample"))
  expect_equal(nrow(fit$estimates), f$n_muni)
  expect_identical(sum(fit$estimates$area_set == "in_sample"), length(f$drawn))
  expect_match(fit$meta$povmap_sha, "^[0-9a-f]{40}$")
  expect_identical(fit$meta$seed, 1L)
  dg <- fit$diagnostics
  num <- function(m) dg$value[dg$metric == m]
  chr <- function(m) dg$character_value[dg$metric == m]
  expect_identical(num("cv_nfold"), 10)
  expect_gt(num("best_iter"), 0)
  expect_identical(num("B"), 5)
  expect_identical(num("n_subarea_train"), as.numeric(nrow(f$smp)))
  expect_identical(num("level_municipality"), 0)
  expect_identical(chr("cv_folds"), "domain")
  expect_identical(chr("predict_sampled"), "crossfit")
  expect_identical(chr("predict_unsampled"), "foldmean")
  expect_identical(chr("mse_type"), "var")
  # the variance components are the squares of the standard deviations
  expect_equal(num("ran_eff_var"), num("ran_eff_sd")^2)
  expect_equal(num("error_var"), num("error_sd")^2)
})

test_that("fit_megb_sub passes povmap's fold options explicitly", {
  skip_if_not_installed("povmap")
  f <- .fixture_sub()
  seen <- NULL
  local_mocked_bindings(
    megb = function(...) { seen <<- list(...); stop("captured") },
    .package = "povmap")
  # A global option must not change what the wrapper runs.
  withr::local_options(megb.cv_folds = "rows", megb.cv_nfold = 5L,
                       megb.early_stopping = FALSE)
  fit_megb_sub(f$direct, f$smp, f$pop, "migrated_last_12m", B = 2L, seed = 1L)
  expect_identical(seen$cv_folds, "domain")
  expect_identical(seen$cv_nfold, 10L)
  expect_identical(seen$predict_sampled, "crossfit")
  expect_identical(seen$predict_unsampled, "foldmean")
  expect_identical(seen$early_stopping, TRUE)
  expect_identical(seen$mse_type, "var")
  expect_identical(seen$mse, TRUE)
  expect_identical(seen$transformation, "no")
  expect_identical(seen$domains, "municipality")
  expect_identical(seen$smp_weights, "wt")
  expect_identical(seen$pop_weights, "n_households_census")
  expect_false("sub_domains" %in% names(seen))
  expect_false("gradient_params" %in% names(seen))
})

test_that("fit_megb_sub forwards mse_type, defaulting to \"var\"", {
  skip_if_not_installed("povmap")
  f <- .fixture_sub()
  seen <- NULL
  local_mocked_bindings(
    megb = function(...) { seen <<- list(...); stop("captured") },
    .package = "povmap")
  fit_megb_sub(f$direct, f$smp, f$pop, "migrated_last_12m", B = 2L, seed = 1L,
               mse_type = "mse")
  expect_identical(seen$mse_type, "mse")
  fit_megb_sub(f$direct, f$smp, f$pop, "migrated_last_12m", B = 2L, seed = 1L)
  expect_identical(seen$mse_type, "var")
  expect_error(fit_megb_sub(f$direct, f$smp, f$pop, "migrated_last_12m",
                            mse_type = "bias"), "'arg'")
})

test_that("fit_megb_sub forwards params as gradient_params and the bootstrap cores", {
  skip_if_not_installed("povmap")
  f <- .fixture_sub()
  seen <- NULL
  local_mocked_bindings(
    megb = function(...) { seen <<- list(...); stop("captured") },
    .package = "povmap")
  fit_megb_sub(f$direct, f$smp, f$pop, "pers_migrated_12m", B = 3L,
               bootstrap_cores = 2L, params = list(max_depth = 3L, eta = 0.2))
  expect_identical(seen$gradient_params, list(max_depth = 3L, eta = 0.2))
  expect_identical(seen$bootstrap_cores, 2L)
  expect_identical(seen$B, 3L)
  expect_identical(seen$pop_weights, "n_persons_census")
  # no seed given: nothing is passed on
  expect_null(seen$seed)
})

test_that("a povmap error in fit_megb_sub is recorded as NA rows, not thrown", {
  skip_if_not_installed("povmap")
  f <- .fixture_sub()
  local_mocked_bindings(
    megb = function(...) stop("fewer domains than folds"),
    .package = "povmap")
  fit <- fit_megb_sub(f$direct, f$smp, f$pop, "migrated_last_12m", B = 2L, seed = 1L)
  expect_s3_class(fit, "sae_fit")
  expect_equal(nrow(fit$estimates), f$n_muni)
  expect_true(all(is.na(fit$estimates$mse)))
  d <- fit$diagnostics
  expect_identical(d$value[d$metric == "megb_sub_failed"], 1)
  expect_identical(d$character_value[d$metric == "megb_sub_error_message"],
                   "fewer domains than folds")
  expect_identical(fit$meta$method, "megb_sub")
})
