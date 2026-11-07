test_that("fit_fh_robust returns the sae_fit contract on a 60-area fixture", {
  skip_if_not_installed("saeRobust")
  f <- .fixture_robust()
  fit <- fit_fh_robust(f$direct, f$X_in, X_oos = f$X_oos)
  expect_s3_class(fit, "sae_fit")
  expect_s3_class(fit, "sae_fit_fh_robust")
  expect_true(all(c("estimates", "diagnostics", "meta") %in% names(fit)))
  expect_equal(nrow(fit$estimates), f$n_in + f$n_oos)
  expect_setequal(colnames(fit$estimates),
                   c("area_id", "area_set", "theta_hat", "mse", "ci_lo", "ci_hi"))
  expect_equal(sum(fit$estimates$area_set == "in_sample"), f$n_in)
  expect_equal(sum(fit$estimates$area_set == "out_of_sample"), f$n_oos)
  expect_true(all(!is.na(fit$estimates$theta_hat)))
  expect_true(all(fit$estimates$mse > 0))
  expect_true(all(fit$estimates$ci_lo < fit$estimates$ci_hi))
  expect_no_error(sae_fit_contract(fit))
})

test_that("fit_fh_robust OOS theta_hat equals X_oos %*% beta_hat (synthetic prediction)", {
  skip_if_not_installed("saeRobust")
  f <- .fixture_robust()
  fit <- fit_fh_robust(f$direct, f$X_in, X_oos = f$X_oos)
  beta_hat <- fit$diagnostics$value[fit$diagnostics$metric == "beta_intercept"]
  beta_x1  <- fit$diagnostics$value[fit$diagnostics$metric == "beta_x1"]
  beta_x2  <- fit$diagnostics$value[fit$diagnostics$metric == "beta_x2"]
  oos_rows <- fit$estimates[fit$estimates$area_set == "out_of_sample", ]
  expected <- beta_hat + beta_x1 * f$X_oos$x1 + beta_x2 * f$X_oos$x2
  ord <- match(oos_rows$area_id, f$X_oos$municipality)
  expect_equal(oos_rows$theta_hat, expected[ord], tolerance = 1e-8)
})

test_that("fit_fh_robust OOS MSE is the per-area synthetic-predictor variance", {
  skip_if_not_installed("saeRobust")
  f <- .fixture_robust()
  fit <- fit_fh_robust(f$direct, f$X_in, X_oos = f$X_oos)
  sigma2_v <- fit$diagnostics$value[fit$diagnostics$metric == "sigma2_v"]
  oos_mse  <- fit$estimates$mse[fit$estimates$area_set == "out_of_sample"]
  # synthetic-predictor variance = sigma2_v + x'Var(beta)x: varies per area
  # and is bounded below by sigma2_v (the quadratic form is non-negative).
  expect_gt(length(unique(round(oos_mse, 12))), 1L)
  expect_true(all(oos_mse >= sigma2_v - 1e-12))
})

test_that("fit_fh_robust returns NA estimates and a diagnostic flag on rfh failure", {
  skip_if_not_installed("saeRobust")
  f <- .fixture_robust(n_in = 8L, n_oos = 0L)  # too small for stable rfh
  f$X_in$x1 <- 0  # zero-variance covariate triggers degenerate fit
  f$X_in$x2 <- 0
  fit <- suppressWarnings(fit_fh_robust(f$direct, f$X_in))
  expect_s3_class(fit, "sae_fit")
  expect_true(all(is.na(fit$estimates$theta_hat)))
  expect_true(any(fit$diagnostics$metric == "rfh_failed"))
})

test_that("fit_fh_robust uses the REBLUP variant consistent with the MSE (not reblupbc)", {
  skip_if_not_installed("saeRobust")
  f <- .fixture_robust()
  # standardize = FALSE so the hand re-fit below is an exact reference: this
  # test is about WHICH predictor variant is used, not about covariate scaling.
  fit <- fit_fh_robust(f$direct, f$X_in, X_oos = f$X_oos, standardize = FALSE)

  # Re-fit by hand to read both the corrected and uncorrected variants.
  x_cols <- setdiff(colnames(f$X_in), "municipality")
  in_df <- data.frame(
    area_id      = f$direct$area_id,
    theta_direct = f$direct$theta_direct,
    vardir       = pmax(f$direct$var_direct_smoothed, .Machine$double.eps),
    f$X_in[match(f$direct$area_id, f$X_in$municipality), x_cols, drop = FALSE]
  )
  rhs   <- paste(x_cols, collapse = " + ")
  fixed <- stats::as.formula(paste0("theta_direct ~ ", rhs))
  hand_fit <- saeRobust::rfh(formula = fixed, data = in_df,
                             samplingVar = "vardir", k = 1.345)
  pred <- predict(hand_fit)
  in_rows <- fit$estimates[fit$estimates$area_set == "in_sample", ]
  ord <- match(in_rows$area_id, in_df$area_id)

  # saeRobust 0.5.0 does not expose reblupbc in predict(); theta_hat must
  # match the plain reblup (the consistent choice given predType = "reblup"
  # MSE). If a future version adds reblupbc, update both and flip this test.
  expect_equal(in_rows$theta_hat, as.numeric(pred$reblup[ord]),
               tolerance = 1e-8)
  # Confirm reblupbc is absent in this version (guards against silent API
  # changes that would reintroduce the inconsistency).
  expect_false("reblupbc" %in% names(pred),
               label = "reblupbc absent from predict() in installed saeRobust")
})

test_that("the two paths agree to solver tolerance on a WELL-CONDITIONED design", {
  skip_if_not_installed("saeRobust")
  f <- .fixture_robust()
  # SCOPE. This establishes agreement where both paths can run, which is a
  # weaker property than "standardize is inert" and is the one that is true.
  # The transformation is exactly affine, so in exact arithmetic the paths are
  # identical; what separates them is the numerical error of the
  # unstandardised solve, which grows with cond(X). On this well-conditioned
  # fixture it is negligible; on a badly scaled matrix it is not, so the bound
  # used here is a property of a well-conditioned design.
  #
  # The rescue property -- that standardising makes an ill-conditioned design
  # solvable at all, which is why the flag exists -- is the test below.
  raw <- fit_fh_robust(f$direct, f$X_in, X_oos = f$X_oos, standardize = FALSE)
  std <- fit_fh_robust(f$direct, f$X_in, X_oos = f$X_oos, standardize = TRUE)
  a <- raw$estimates[order(raw$estimates$area_id), ]
  b <- std$estimates[order(std$estimates$area_id), ]
  expect_equal(a$area_id, b$area_id)
  # Judge on absolute difference scaled by the SPREAD of each quantity. A
  # relative test is the wrong measure here because ci_lo legitimately crosses
  # zero on rare outcomes, so dividing by it inflates the ratio without the
  # numbers moving.
  for (cc in c("theta_hat", "mse", "ci_lo", "ci_hi")) {
    sd_cc <- stats::sd(a[[cc]], na.rm = TRUE)
    if (!is.finite(sd_cc) || sd_cc == 0) next
    # 1e-5 rather than tighter because rfh is iterative: on this small fixture
    # the two paths converge to points differing by ~1.2e-6 of the spread.
    expect_lt(max(abs(a[[cc]] - b[[cc]]), na.rm = TRUE) / sd_cc, 1e-5)
  }
  # Coefficients are back-transformed, so the reported betas match too.
  ga <- raw$diagnostics; gb <- std$diagnostics
  bi <- grep("^beta_", ga$metric)
  expect_equal(ga$value[bi], gb$value[bi], tolerance = 1e-6)
})

# An ill-conditioned twin of the fixture: one covariate on a scale many orders
# of magnitude above the others, the shape of a matrix that mixes population
# counts with shares near 1.
.fixture_illcond <- function(scale = 1e9, ...) {
  f <- .fixture_robust(...)
  f$X_in$x2  <- f$X_in$x2  * scale
  f$X_oos$x2 <- f$X_oos$x2 * scale
  f
}

.rcond_xtx <- function(X) {
  M <- cbind(1, as.matrix(X[, setdiff(colnames(X), "municipality"), drop = FALSE]))
  rcond(crossprod(M))
}

test_that("standardising RESCUES an ill-conditioned design", {
  skip_if_not_installed("saeRobust")
  f <- .fixture_illcond()

  # rfh inverts a design-derived cross-product, so the quantity that decides
  # whether it factorises is rcond(X'X), not rcond(X). Assert that first,
  # because it is deterministic arithmetic and travels across machines.
  x_cols <- setdiff(colnames(f$X_in), "municipality")
  std_X  <- f$X_in
  for (cc in x_cols) std_X[[cc]] <- scale(std_X[[cc]])[, 1]
  expect_lt(.rcond_xtx(f$X_in), 1e-16)
  expect_gt(.rcond_xtx(std_X), 1e-6)
  expect_gt(log10(.rcond_xtx(std_X)) - log10(.rcond_xtx(f$X_in)), 10)

  # The standardised path must return a usable fit on a design this badly
  # scaled. Unconditional: this is the property the flag exists for.
  std <- fit_fh_robust(f$direct, f$X_in, X_oos = f$X_oos, standardize = TRUE)
  expect_equal(
    std$diagnostics$value[std$diagnostics$metric == "rfh_failed"], 0)
  expect_true(all(is.finite(std$estimates$theta_hat)))
  expect_true(all(std$estimates$mse > 0))
  expect_true(all(std$estimates$ci_lo < std$estimates$ci_hi))
  expect_no_error(sae_fit_contract(std))

  # The other half of "rescue": the unstandardised path must NOT survive, or
  # the fixture is not demonstrating anything. Guarded by skip rather than
  # asserted, because where rfh gives up is a property of the LAPACK build.
  # rcond(X'X) of this fixture by `scale`, and whether the unstandardised fit
  # fails, on a reference LAPACK build:
  #
  #   scale   rcond(X'X)   raw fails?
  #   1e6       7.7e-13        no
  #   1e7       7.7e-15        no
  #   1e8       7.7e-17        YES
  #   1e9       7.7e-19        YES
  #
  # so the boundary sits between 1e-15 and 1e-17, near .Machine$double.eps as
  # expected, and 1e9 clears it by two orders. A build that copes anyway skips
  # here rather than turning the suite red.
  raw <- fit_fh_robust(f$direct, f$X_in, X_oos = f$X_oos, standardize = FALSE)
  raw_failed <- raw$diagnostics$value[raw$diagnostics$metric == "rfh_failed"]
  skip_if(raw_failed == 0,
          "this LAPACK factorised rcond 7.7e-19; rescue not demonstrable here")
  expect_equal(raw_failed, 1)
  expect_true(all(is.na(raw$estimates$theta_hat[
    raw$estimates$area_set == "in_sample"])))
})
