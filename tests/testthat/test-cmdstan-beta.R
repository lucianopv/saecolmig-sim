test_that("capture_tipsae_standata returns the 41-field ExtBeta contract", {
  skip_if_not_installed("tipsae")
  # .fixture_robust() holds out at least one area; X_oos = NULL exercises the
  # path without out-of-sample areas.
  fx <- .fixture_robust(n_in = 40L, n_oos = 1L, seed = 21L)
  d  <- fx$direct; d$n_dwell <- 40L
  cap <- capture_tipsae_standata(d, fx$X_in, X_oos = NULL)
  sd <- cap$standata
  expect_true(all(c("y","X","X_oos","M_is","M_oos","P","disp","likelihood",
                    "deff","m_d") %in% names(sd)))
  expect_equal(sd$likelihood, 3L)   # ExtBeta
  expect_equal(sd$deff, 1L)         # neff parameterisation
  expect_equal(sd$M_is, 40L)
  expect_equal(length(cap$in_sample_ids), 40L)  # area_id order for theta[]
})

test_that("vendored tipsae global.stan compiles under cmdstanr with threads", {
  skip_if_not_installed("cmdstanr")
  stan_path <- system.file("stan", "tipsae", "global.stan",
                           package = "saecol2026")
  skip_if(stan_path == "", "vendored global.stan not installed")
  # Include directives are #include /include/x.stan, so include_paths must be
  # the parent of include/ (i.e., .../stan/tipsae), not include/ itself.
  inc <- system.file("stan", "tipsae", package = "saecol2026")
  m <- cmdstanr::cmdstan_model(stan_file = stan_path,
                               include_paths = inc,
                               cpp_options = list(stan_threads = TRUE),
                               compile = TRUE, quiet = TRUE)
  expect_s3_class(m, "CmdStanModel")
})

test_that("fit_bayes_fh_beta(backend='cmdstan') matches rstan within MC tol", {
  skip_if_not_installed("cmdstanr"); skip_if_not_installed("tipsae")
  # .fixture_robust() holds out at least one area; X_oos = NULL in both fits
  # leaves neither with an out-of-sample area.
  fx <- .fixture_robust(n_in = 50L, n_oos = 1L, seed = 31L)
  d  <- fx$direct; d$n_dwell <- 40L
  # Run BOTH backends on the SAME (ExtBeta / vendored global) model for an
  # apples-to-apples comparison: rstan via tipsae, cmdstan via vendored model.
  # Both warm up for half of `iter` (400 warmup + 400 sampling per chain).
  # rstan reports a low effective sample size for these short chains.
  r  <- muffle_warnings(
    fit_bayes_fh_beta(d, fx$X_in, X_oos = NULL, backend = "rstan",
                      likelihood = "ExtBeta",
                      chains = 2L, iter = 800L, seed = 7L),
    "Effective Samples Size|The largest R-hat")
  c2 <- fit_bayes_fh_beta(d, fx$X_in, X_oos = NULL, backend = "cmdstan",
                          likelihood = "ExtBeta",
                          chains = 2L, iter = 800L, seed = 7L)
  expect_s3_class(c2, "sae_fit")
  sae_fit_contract(c2)
  j <- merge(r$estimates[, c("area_id","theta_hat")],
             c2$estimates[, c("area_id","theta_hat")],
             by = "area_id", suffixes = c("_r","_c"))
  # Monte Carlo noise at the test's reduced chains and iterations: at the
  # runner's settings the two backends agree to about 0.002 in max |theta_hat
  # difference|; 0.03 bounds the noise of these short chains, not a model
  # difference.
  expect_lt(max(abs(j$theta_hat_r - j$theta_hat_c)), 0.03)
})

test_that("cmdstan backend tags in/OOS area sets disjointly and completely", {
  skip_if_not_installed("cmdstanr")
  skip_if_not_installed("tipsae")
  skip_if(is.null(tryCatch(cmdstanr::cmdstan_path(),
                           error = function(e) NULL)),
          "CmdStan toolchain not installed")
  inp <- .fixture_robust(n_in = 20L, n_oos = 5L)
  inp$direct$n_dwell <- as.numeric(inp$direct$n_eff)
  res <- fit_bayes_fh_beta(inp$direct, inp$X_in, X_oos = inp$X_oos,
                           chains = 1L, iter = 400L,
                           likelihood = "ExtBeta",
                           backend = "cmdstan", seed = 42L)
  est <- res$estimates
  expect_identical(nrow(est), 25L)
  in_ids  <- est$area_id[est$area_set == "in_sample"]
  oos_ids <- est$area_id[est$area_set == "out_of_sample"]
  expect_setequal(in_ids,  inp$direct$area_id)
  expect_setequal(oos_ids, inp$X_oos$municipality)
  expect_length(intersect(in_ids, oos_ids), 0L)
  # n_divergent/max_rhat must flow to diagnostics for downstream run metadata
  expect_true("n_divergent" %in% res$diagnostics$metric)
  expect_true("max_rhat" %in% res$diagnostics$metric)
  # Out of sample the estimate is tipsae's theta_oos; nothing else is reported.
  expect_setequal(names(est), c("area_id", "theta_hat", "mse", "ci_lo", "ci_hi", "area_set"))
  expect_setequal(res$diagnostics$metric, c("n_divergent", "max_rhat"))
})

# Inputs of the next two tests; the sampler itself is replaced by a stub, so
# neither test needs CmdStan.
.cmdstan_guard_inputs <- function() {
  list(d = data.frame(area_id = 1:3, theta_direct = c(0.05, 0.10, 0.20),
                      var_direct_smoothed = 0.01, n_dwell = 50L),
       X = data.frame(municipality = 1:3, x1 = c(0.1, 0.5, 0.9)))
}

test_that("backend = 'cmdstan' stops on a sampler argument it does not use", {
  skip_if_not_installed("tipsae"); skip_if_not_installed("rstan")
  inp <- .cmdstan_guard_inputs()
  # A stub that fails differently from the guard: an error that does not name
  # the argument means the guard did not stop the call.
  local_mocked_bindings(fit_bayes_fh_beta_cmdstan = function(...) stop("sampler reached"))
  fit <- function(...) fit_bayes_fh_beta(inp$d, inp$X, likelihood = "ExtBeta",
                                         backend = "cmdstan", ...)
  # The message puts the offending arguments, and only those, in backticks; the
  # explanation that follows does not.
  offenders <- function(...) {
    msg <- tryCatch(fit(...), error = conditionMessage)
    regmatches(msg, gregexpr("`[A-Za-z0-9_]+`", msg))[[1]]
  }
  expect_identical(offenders(adapt_delta = 0.99), "`adapt_delta`")
  expect_identical(offenders(warmup = 500L), "`warmup`")
  # the backend warms up for 2000 iterations here, whatever the default of warmup is
  expect_identical(offenders(iter = 4000L, warmup = 1000L), "`warmup`")
  expect_identical(offenders(prior_coeff = "HorseShoe"), "`prior_coeff`")
  expect_identical(offenders(p0_HorseShoe = 5L), "`p0_HorseShoe`")
  expect_setequal(offenders(adapt_delta = 0.99, warmup = 500L),
                  c("`adapt_delta`", "`warmup`"))
  expect_error(fit_bayes_fh_beta(inp$d, inp$X, likelihood = "Infbeta0",
                                 backend = "cmdstan"), "`likelihood`", fixed = TRUE)
  # the backend alone does not choose the model: ExtBeta has to be asked for
  expect_error(fit_bayes_fh_beta(inp$d, inp$X, backend = "cmdstan"),
               "`likelihood`", fixed = TRUE)
})

test_that("backend = 'cmdstan' accepts the arguments it uses and the values it applies", {
  skip_if_not_installed("tipsae"); skip_if_not_installed("rstan")
  inp <- .cmdstan_guard_inputs()
  seen <- NULL
  local_mocked_bindings(fit_bayes_fh_beta_cmdstan = function(direct, X, ...) {
    seen <<- list(...)
    "sampler reached"
  })
  # the runner's call: method_args() overlaid with the run settings
  args <- utils::modifyList(
    method_args("bayes_fh_beta", direct = inp$d, aux_in = inp$X, aux_oos = NULL, seed = 7L),
    sim_settings()$bayes_fh_beta)
  expect_identical(do.call(fit_bayes_fh_beta, args), "sampler reached")
  expect_setequal(names(seen), c("X_oos", "chains", "iter", "seed"))
  expect_identical(seen$chains, 3L)
  expect_identical(seen$seed, 7L)
  # the defaults, and values equal to what the backend applies, are not a request
  # for something else
  expect_identical(fit_bayes_fh_beta(inp$d, inp$X, likelihood = "ExtBeta", backend = "cmdstan"),
                   "sampler reached")
  expect_identical(fit_bayes_fh_beta(inp$d, inp$X, likelihood = "ExtBeta", backend = "cmdstan",
                                     iter = 800L, warmup = 400L, adapt_delta = 0.95,
                                     prior_coeff = "normal"), "sampler reached")
  # an explicit warmup equal to half of iter, with the default iter
  expect_identical(fit_bayes_fh_beta(inp$d, inp$X, likelihood = "ExtBeta", backend = "cmdstan",
                                     warmup = 1000L), "sampler reached")
})

# ---- the direct estimates, the start and the output ---------------------------

# 40 areas with a rare outcome, 6 of them exact zeros, each with 1,200 sampled
# households; 4 areas out of sample.
.zero_fixture <- function(m = 1200L) {
  withr::local_seed(11)
  n <- 40L
  X <- data.frame(municipality = 1:(n + 4L), x1 = rnorm(n + 4L), x2 = rnorm(n + 4L))
  y <- rbinom(n, m, plogis(-4 + 0.5 * X$x1[1:n])) / m
  y[1:6] <- 0
  list(d = data.frame(area_id = 1:n, theta_direct = y,
                      var_direct_smoothed = pmax(y * (1 - y) / m, 1e-8),
                      n_dwell = m, n_eff = m),
       X_in = X[1:n, ], X_oos = X[(n + 1L):(n + 4L), ])
}

test_that("the direct estimates reach tipsae unmodified, exact zeros included", {
  skip_if_not_installed("tipsae"); skip_if_not_installed("rstan")
  fx <- .zero_fixture(m = 40L)
  cap <- capture_tipsae_standata(fx$d, fx$X_in, X_oos = fx$X_oos)
  expect_identical(sort(unname(cap$standata$y)), sort(fx$d$theta_direct))
  expect_identical(sum(cap$standata$y == 0), sum(fx$d$theta_direct == 0))
  expect_gte(sum(cap$standata$y == 0), 6L)
  expect_identical(cap$oos_ids, fx$X_oos$municipality)
})

test_that("a plain beta fit on direct estimates with zeros stops with tipsae's message", {
  skip_if_not_installed("tipsae"); skip_if_not_installed("rstan")
  fx <- .zero_fixture(m = 40L)
  for (lik in c("beta", "flexbeta")) {
    expect_error(fit_bayes_fh_beta(fx$d, fx$X_in, X_oos = fx$X_oos, likelihood = lik,
                                   backend = "rstan", chains = 1L, iter = 20L, warmup = 10L),
                 "direct estimates equal to 0 or 1")
  }
})

test_that("the Bayes-Beta functions, settings and per-fit columns are exactly these", {
  expect_setequal(names(formals(fit_bayes_fh_beta)),
                  c("direct", "X", "geo", "X_oos", "chains", "iter", "warmup", "likelihood",
                    "prior_coeff", "p0_HorseShoe", "adapt_delta", "backend", "seed", "..."))
  expect_setequal(names(formals(fit_bayes_fh_beta_cmdstan)),
                  c("direct", "X", "X_oos", "chains", "iter", "adapt_delta", "max_treedepth", "seed"))
  expect_setequal(names(formals(capture_tipsae_standata)), c("direct", "X", "X_oos"))
  expect_setequal(names(sim_settings()$bayes_fh_beta), c("likelihood", "chains", "iter", "backend"))
  expect_identical(sim_fit_columns,
                   c("area_id", "theta_hat", "mse", "ci_lo", "ci_hi", "area_set", "fit_failed",
                     "error_msg", "config", "method", "replicate", "outcome", "seed", "wall_secs",
                     "weight", "n_divergent", "max_rhat", "sigma2_v", "povmap_sha", "mse_type",
                     "pkg_versions"))
})

test_that("extbeta_init_intercept starts every active parameter at 0 except beta0", {
  sd <- list(likelihood = 3L, intercept = 1L, prior_coeff = 0L, prior_reff = 0L,
             spatial_err = 0L, temporal_err = 0L, P = 4L, M_is = 10L,
             y = c(rep(0, 6), 0.02, 0.05, 0.01, 0.03))
  i <- extbeta_init_intercept(sd)
  expect_named(i, c("beta0", "z_beta", "v_raw", "sigma_v", "lambda_star_EB"))
  expect_equal(as.numeric(i$beta0), qlogis(mean(sd$y)))
  expect_equal(as.numeric(i$z_beta), rep(0, 4)); expect_equal(i$v_raw, rep(0, 10))
  # one covariate: z_beta is still an array, which rstan requires
  sd1 <- modifyList(sd, list(P = 1L, prior_coeff = 0L))
  expect_identical(dim(extbeta_init_intercept(sd1)$z_beta), 1L)
  expect_equal(as.numeric(i$sigma_v), 1); expect_equal(as.numeric(i$lambda_star_EB), 0.5)
  sd$prior_coeff <- 1L                       # horseshoe: unsupported, refuse
  expect_error(extbeta_init_intercept(sd), "only the ExtBeta structure of fit_bayes_fh_beta")
})

test_that("rstan::sampling is restored after an ExtBeta start, also when the fit fails", {
  skip_if_not_installed("rstan")
  before <- rstan::sampling
  expect_error(.with_extbeta_start(stop("stop inside")), "stop inside")
  expect_identical(rstan::sampling, before)
})

test_that("ExtBeta fits an exact zero at a large sample size on both backends", {
  skip_if_not_installed("tipsae"); skip_if_not_installed("rstan")
  fx <- .zero_fixture(m = 1200L)
  # From tipsae's start every mu is 0.5, so the zero mass of a zero area is
  # 0.5^1200, which is 0 in double precision: the chains cannot start there.
  expect_identical(0.5^1200, 0)
  before <- rstan::sampling
  rs <- muffle_warnings(
    fit_bayes_fh_beta(fx$d, fx$X_in, X_oos = fx$X_oos, likelihood = "ExtBeta",
                      backend = "rstan", chains = 2L, iter = 300L, seed = 3L),
    "Effective Samples Size|The largest R-hat")
  expect_identical(rstan::sampling, before)
  expect_identical(nrow(rs$estimates), 44L)
  expect_identical(sum(rs$estimates$area_set == "out_of_sample"), 4L)
  expect_true(all(is.finite(rs$estimates$theta_hat)) && all(rs$estimates$mse > 0))
  expect_setequal(names(rs$estimates), c("area_id", "area_set", "theta_hat", "mse", "ci_lo", "ci_hi"))
  expect_setequal(rs$diagnostics$metric,
                  c("frac_rhat_above_1.05", "min_ess_bulk", "max_rhat", "median_ess_bulk", "n_divergent"))
  expect_true(is.finite(rs$diagnostics$value[rs$diagnostics$metric == "max_rhat"]))
  skip_if_not_installed("cmdstanr")
  skip_if(is.null(tryCatch(cmdstanr::cmdstan_path(), error = function(e) NULL)),
          "CmdStan toolchain not installed")
  cs <- fit_bayes_fh_beta(fx$d, fx$X_in, X_oos = fx$X_oos, likelihood = "ExtBeta",
                          backend = "cmdstan", chains = 2L, iter = 300L, seed = 3L)
  expect_identical(nrow(cs$estimates), 44L)
  j <- merge(rs$estimates[, c("area_id", "theta_hat")], cs$estimates[, c("area_id", "theta_hat")],
             by = "area_id", suffixes = c("_r", "_c"))
  # the two backends sample the same posterior; short chains bound the difference
  expect_lt(max(abs(j$theta_hat_r - j$theta_hat_c)), 0.005)
})

test_that("area ids follow tipsae's sort when direct and X_oos are not sorted by area id", {
  skip_if_not_installed("tipsae"); skip_if_not_installed("rstan")
  fx <- .zero_fixture(m = 40L)
  rev_d   <- fx$d[rev(seq_len(nrow(fx$d))), ]
  rev_oos <- fx$X_oos[rev(seq_len(nrow(fx$X_oos))), ]
  cap <- capture_tipsae_standata(rev_d, fx$X_in, X_oos = rev_oos)
  expect_identical(unname(cap$standata$y), rev_d$theta_direct[match(cap$in_sample_ids, rev_d$area_id)])
  expect_identical(cap$in_sample_ids, sort(fx$d$area_id))
  expect_identical(cap$oos_ids, sort(fx$X_oos$municipality))
  skip_if_not_installed("cmdstanr")
  skip_if(is.null(tryCatch(cmdstanr::cmdstan_path(), error = function(e) NULL)),
          "CmdStan toolchain not installed")
  # tipsae sorts its data, so the data reach Stan in the same row order (the dispersion
  # smoothing can differ in the last bits); the fits agree by id within Monte Carlo noise
  fit <- function(d, xo) fit_bayes_fh_beta(d, fx$X_in, X_oos = xo, likelihood = "ExtBeta",
                                           backend = "cmdstan", chains = 2L, iter = 200L, seed = 5L)$estimates
  j <- merge(fit(fx$d, fx$X_oos), fit(rev_d, rev_oos), by = c("area_id", "area_set"),
             suffixes = c("_s", "_r"))
  expect_identical(nrow(j), 44L)
  expect_lt(max(abs(j$theta_hat_s - j$theta_hat_r)), 0.03)
})

test_that("ExtBeta runs with one covariate on both backends", {
  skip_if_not_installed("tipsae"); skip_if_not_installed("rstan")
  fx <- .zero_fixture(m = 40L)
  X1 <- fx$X_in[, c("municipality", "x1")]; Xo1 <- fx$X_oos[, c("municipality", "x1")]
  rs <- muffle_warnings(
    fit_bayes_fh_beta(fx$d, X1, X_oos = Xo1, likelihood = "ExtBeta", backend = "rstan",
                      chains = 2L, iter = 200L, seed = 1L),
    "Effective Samples Size|The largest R-hat")
  expect_identical(nrow(rs$estimates), 44L)
  skip_if_not_installed("cmdstanr")
  skip_if(is.null(tryCatch(cmdstanr::cmdstan_path(), error = function(e) NULL)),
          "CmdStan toolchain not installed")
  cs <- fit_bayes_fh_beta(fx$d, X1, X_oos = Xo1, likelihood = "ExtBeta", backend = "cmdstan",
                          chains = 2L, iter = 200L, seed = 1L)
  expect_identical(nrow(cs$estimates), 44L)
})

test_that("a fit for which rstan returns no draws stops with a message that points to rstan", {
  skip_if_not_installed("tipsae"); skip_if_not_installed("rstan")
  withr::local_seed(3)
  n <- 30L
  X <- data.frame(municipality = 1:n, x1 = rnorm(n))
  y <- rbinom(n, 4000L, 0.5) / 4000L; y[1] <- 0
  d <- data.frame(area_id = 1:n, theta_direct = y, var_direct_smoothed = pmax(y * (1 - y) / 4000, 1e-8),
                  n_dwell = 4000L, n_eff = 4000)
  # y about 0.5 at 4,000 households: even the intercept start gives the zero area an
  # underflowing zero mass, so no chain starts.
  utils::capture.output(
    expect_error(fit_bayes_fh_beta(d, X, likelihood = "ExtBeta", backend = "rstan",
                                   chains = 1L, iter = 40L, seed = 1L),
                 "rstan returned no draws; see its messages above"),
    type = "message")
})

# 40 areas without exact zeros, unequal sizes, a variance that is not exactly
# proportional to mu (1 - mu) / n (tipsae::smoothing() cannot fit an exactly
# proportional one), 4 areas out of sample, 3 covariates.
.free_fixture <- function() {
  withr::local_seed(12)
  n <- 40L
  X <- data.frame(municipality = 1:(n + 4L), x1 = rnorm(n + 4L), x2 = rnorm(n + 4L), x3 = rnorm(n + 4L))
  m <- sample(50:100, n, replace = TRUE)
  y <- pmax(rbinom(n, m, plogis(-2 + 0.5 * X$x1[1:n])) / m, 0.01)
  list(d = data.frame(area_id = 1:n, theta_direct = y,
                      var_direct_smoothed = y * (1 - y) / m * runif(n, 0.8, 1.2),
                      n_dwell = m, n_eff = m),
       X_in = X[1:n, ], X_oos = X[(n + 1L):(n + 4L), ])
}

test_that("the ExtBeta start is used with the normal prior and not with another prior", {
  skip_if_not_installed("tipsae"); skip_if_not_installed("rstan"); skip_if_not_installed("posterior")
  fx <- .free_fixture()
  calls <- 0L; f <- NULL
  # the mock evaluates the fit, keeps tipsae's return value and passes it on
  local_mocked_bindings(.with_extbeta_start = function(expr) { calls <<- calls + 1L; f <<- expr; f })
  rs <- muffle_warnings(
    fit_bayes_fh_beta(fx$d, fx$X_in, X_oos = fx$X_oos, likelihood = "ExtBeta", backend = "rstan",
                      chains = 2L, iter = 200L, seed = 1L),
    "Effective Samples Size|The largest R-hat")
  expect_identical(calls, 1L)
  # max_rhat is the rank-normalised split R-hat of the posterior package over the
  # variables the fit keeps
  rhat <- max(posterior::summarise_draws(posterior::as_draws_array(f$stanfit), "rhat")$rhat,
              na.rm = TRUE)
  expect_equal(rs$diagnostics$value[rs$diagnostics$metric == "max_rhat"], rhat)
  hs <- muffle_warnings(
    fit_bayes_fh_beta(fx$d, fx$X_in, X_oos = fx$X_oos, likelihood = "ExtBeta", backend = "rstan",
                      prior_coeff = "HorseShoe", p0_HorseShoe = 1L,
                      chains = 1L, iter = 100L, seed = 1L),
    "Effective Samples Size|The largest R-hat|divergent|Bulk|Tail|ESS")
  expect_identical(calls, 1L)
  expect_identical(nrow(hs$estimates), 44L)
})
