test_that("method_args filters outcome-long direct to the requested outcome", {
  direct <- tibble::tibble(
    area_id = rep(1:3, 2),
    theta_direct = 0.1, var_direct_smoothed = 0.01, n_eff = 50,
    outcome = rep(c("migrated_last_12m", "migrated_intl_12m"), each = 3))
  a <- method_args("fh", direct = direct,
                   aux_in = data.frame(municipality = 1:3),
                   aux_oos = data.frame(municipality = 4:5),
                   outcome = "migrated_intl_12m")
  expect_equal(unique(a$direct$outcome), "migrated_intl_12m")
  expect_equal(nrow(a$direct), 3L)
  expect_error(
    method_args("fh", direct = direct,
                aux_in = data.frame(), aux_oos = data.frame(),
                outcome = "not_an_outcome"),
    "No direct estimates")
})

test_that("every implemented sae_method satisfies the contract on synthetic input", {
  skip_if_not_installed("emdi2")
  skip_if_not_installed("xgboost")
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

  in_ids  <- area_ids[1:20]
  oos_ids <- area_ids[21:30]
  direct_in <- tibble::tibble(
    area_id              = in_ids,
    theta_direct         = pmin(pmax(theta_true[1:20] +
                                      rnorm(20, 0, sqrt(var_d[1:20])),
                                    0.001), 0.999),
    var_direct           = var_d[1:20],
    var_direct_smoothed  = var_d[1:20],
    n_eff                = n_eff[1:20],
    n_dwell              = round(n_eff[1:20]),
    stratum              = paste0("S", (in_ids %/% 1000L))
  )
  X_in  <- X[X$municipality %in% in_ids, ]
  X_oos <- X[X$municipality %in% oos_ids, ]
  W <- build_dept_adjacency(area_ids)
  geo <- list(W = W)

  for (m in sae_methods_implemented) {
    if (m == "bayes_fh_beta") next  # tested in its own file with skip_on_ci()
    if (m %in% sae_methods_subarea) next  # (direct, smp_sub, pop_sub, outcome) signature; tested below
    fn   <- sae_methods[[m]]
    args <- list(direct = direct_in, X = X_in, X_oos = X_oos, geo = geo)
    fit <- suppressWarnings(do.call(fn, args))
    expect_silent(sae_fit_contract(fit))
    expect_equal(fit$meta$method, m)
    expect_equal(nrow(fit$estimates), n)
    expect_true("area_set" %in% colnames(fit$estimates))
    expect_equal(sum(fit$estimates$area_set == "out_of_sample"), 10L)
  }
})

test_that("method_args carries the ExtBeta settings for bayes_fh_beta", {
  direct <- tibble::tibble(
    area_id = 1:3, theta_direct = 0.1,
    var_direct_smoothed = 0.01, n_eff = 50, n_dwell = 40L,
    outcome = "migrated_last_12m")
  a <- method_args("bayes_fh_beta", direct = direct,
                   aux_in = data.frame(municipality = 1:3),
                   aux_oos = data.frame(municipality = 4:5))
  expect_equal(a$likelihood, "ExtBeta")
  expect_equal(a$iter, 2000L)
})

test_that("the sub-area methods are registered as a pair", {
  expect_identical(sae_methods_subarea, c("xgb_sub", "megb_sub"))
  expect_true(all(sae_methods_subarea %in% names(sae_methods)))
  expect_true(all(sae_methods_subarea %in% sae_methods_implemented))
  expect_identical(sae_methods$xgb_sub, fit_xgb_sub)
  expect_identical(sae_methods$megb_sub, fit_megb_sub)
})

test_that("every sub-area method satisfies the contract on synthetic input", {
  skip_if_not_installed("povmap")
  skip_if_not_installed("xgboost")
  f <- .fixture_sub(n_muni = 30L, per_muni = 8L)
  for (m in sae_methods_subarea) {
    fn <- sae_methods[[m]]
    fit <- .quiet_fit(
      fn(f$direct, f$smp, f$pop, "migrated_last_12m", B = 5L, seed = 1L,
         params = if (m == "xgb_sub") list(nrounds = 15L)))
    expect_false(any(fit$diagnostics$metric == paste0(m, "_failed")), info = m)
    expect_silent(sae_fit_contract(fit))
    expect_equal(fit$meta$method, m)
    expect_equal(nrow(fit$estimates), f$n_muni, info = m)
    expect_setequal(unique(fit$estimates$area_set),
                    c("in_sample", "out_of_sample"))
    expect_identical(sum(fit$estimates$area_set == "out_of_sample"),
                     f$n_muni - length(f$drawn), info = m)
  }
})

test_that("method_args refuses the sub-area methods, which take their arguments directly", {
  direct <- tibble::tibble(area_id = 1:3, theta_direct = 0.1,
                           var_direct_smoothed = 0.01, n_eff = 50,
                           outcome = "migrated_last_12m")
  for (m in sae_methods_subarea) {
    expect_error(
      method_args(m, direct = direct, aux_in = data.frame(municipality = 1:3),
                  aux_oos = data.frame(municipality = 4:5)),
      paste0(m, " is a sub-area method and takes its arguments directly"),
      fixed = TRUE)
  }
  expect_error(
    method_args("xgb_sub", direct = direct, aux_in = NULL, aux_oos = NULL),
    "fit_xgb_sub()", fixed = TRUE)
  # an area-level method is still dispatched
  expect_named(method_args("fh", direct = direct,
                           aux_in = data.frame(municipality = 1:3),
                           aux_oos = data.frame(municipality = 4:5)),
               c("direct", "X", "X_oos"))
})
