fit_path <- function(out, method, config, rep) {
  file.path(out, "fits", paste0(method, "__", config), sprintf("rep_%04d.parquet", rep))
}

# A fit in the method contract: fixed estimates for `in_ids` (in sample) and
# `oos_ids` (out of sample), with the given diagnostics.
fake_fit <- function(in_ids, oos_ids, diagnostics = tibble::tibble(metric = "x", value = 1)) {
  ids <- c(in_ids, oos_ids)
  est <- tibble::tibble(
    area_id = ids, theta_hat = 0.1, mse = 0.01, ci_lo = 0.05, ci_hi = 0.15,
    area_set = rep(c("in_sample", "out_of_sample"), c(length(in_ids), length(oos_ids))))
  structure(list(estimates = est, diagnostics = diagnostics,
                 meta = list(method = "fake", runtime_sec = 0, n_areas = length(ids))),
            class = "sae_fit")
}

# Replace the registry entry of `method` by `fun` for the calling test.
local_method <- function(method, fun, env = parent.frame()) {
  reg <- sae_methods
  reg[[method]] <- fun
  local_mocked_bindings(sae_methods = reg, .env = env)
}

test_that("a replicate writes direct and per-fit files with the documented columns", {
  out <- withr::local_tempdir(); inp <- fixture_sim_inputs()
  suppressMessages(
    sim_run_replicate(1L, inp, methods = "fh", outcomes = "migrated_last_12m", out_dir = out))
  f <- arrow::read_parquet(file.path(out, "fits", "fh__geo_bld", "rep_0001.parquet"))
  expect_true(all(c("area_id", "theta_hat", "mse", "area_set", "seed", "pkg_versions") %in% names(f)))
  expect_equal(nrow(f), length(inp$universe))
  expect_equal(unique(f$seed), sae_fit_seed("geo_bld", "fh", "migrated_last_12m", 1L))
  expect_identical(names(f), sim_fit_columns)
  d <- arrow::read_parquet(file.path(out, "direct", "rep_0001.parquet"))
  expect_identical(names(d), c("area_id", "theta_direct", "var_direct", "n_eff", "n_dwell",
                               "outcome", "var_direct_smoothed"))
  expect_setequal(unique(d$outcome), sae_outcomes$outcome)
  expect_setequal(f$area_id[f$area_set == "in_sample"], unique(d$area_id))
  expect_false(any(grepl("\\.tmp", list.files(out, recursive = TRUE))))
})

test_that("the sample of replicate r uses seed 100000 + r", {
  inp <- fixture_sim_inputs()
  a <- sim_draw(3L, inp)
  b <- draw_replicate(inp$frame, inp$design, seed = 100003L, segments = inp$segments)
  expect_identical(a$sample, apply_calibration(b, inp$pop_totals))
  # Segments built inside the draw select the same households.
  b0 <- draw_replicate(inp$frame, inp$design, seed = 100003L)
  expect_identical(sort(a$sample$hh_id), sort(b0$hh_id))
})

test_that("the draw holds smoothed direct estimates and the sub-area sample of every outcome", {
  inp <- fixture_sim_inputs()
  a <- sim_draw(2L, inp)
  expect_identical(a$direct, smooth_variance_all(
    direct_estimate_all(a$sample, area_var = "municipality"), method = "ghs"))
  for (o in sae_outcomes$outcome) {
    want <- collapse_sample_to_subarea(a$sample, inp$pop_sub, o)
    got <- a$smp_sub[a$smp_sub$outcome == o, setdiff(names(a$smp_sub), "outcome")]
    rownames(got) <- NULL
    expect_equal(got, as.data.frame(want), ignore_attr = TRUE)
  }
})

test_that("a failing fit is recorded and the replicate continues", {
  out <- withr::local_tempdir(); inp <- fixture_sim_inputs()
  local_method("fh", function(...) stop("singular"))
  suppressMessages(sim_run_replicate(1L, inp, methods = c("fh", "fh_arcsin"),
                                     outcomes = "migrated_last_12m", out_dir = out))
  f <- arrow::read_parquet(file.path(out, "fits", "fh__geo_bld", "rep_0001.parquet"))
  expect_true(all(f$fit_failed)); expect_match(f$error_msg[1], "singular")
  expect_true(file.exists(file.path(out, "fits", "fh_arcsin__geo_bld", "rep_0001.parquet")))
  g <- arrow::read_parquet(fit_path(out, "fh_arcsin", "geo_bld", 1L))
  expect_false(any(g$fit_failed))
  expect_equal(nrow(g), length(inp$universe))
})

test_that("a rerun skips complete replicates and rewrites partial ones", {
  out <- withr::local_tempdir(); inp <- fixture_sim_inputs()
  suppressMessages(sim_run(1:2, inp, methods = "fh", outcomes = "migrated_last_12m", out_dir = out))
  p <- file.path(out, "fits", "fh__geo_bld", "rep_0002.parquet")
  writeBin(raw(10), p)                                 # simulate a broken file
  t1 <- file.mtime(file.path(out, "fits", "fh__geo_bld", "rep_0001.parquet"))
  suppressMessages(sim_run(1:2, inp, methods = "fh", outcomes = "migrated_last_12m", out_dir = out))
  expect_equal(file.mtime(file.path(out, "fits", "fh__geo_bld", "rep_0001.parquet")), t1)
  expect_equal(nrow(arrow::read_parquet(p)), length(inp$universe))
})

test_that("a per-fit file holding another replicate is refitted", {
  out <- withr::local_tempdir(); inp <- fixture_sim_inputs()
  suppressMessages(sim_run(1:2, inp, methods = "fh", outcomes = "migrated_last_12m", out_dir = out))
  file.copy(fit_path(out, "fh", "geo_bld", 1L), fit_path(out, "fh", "geo_bld", 2L), overwrite = TRUE)
  suppressMessages(sim_run(2L, inp, methods = "fh", outcomes = "migrated_last_12m", out_dir = out))
  f <- arrow::read_parquet(fit_path(out, "fh", "geo_bld", 2L))
  expect_identical(unique(f$replicate), 2L)
  expect_identical(unique(f$seed), sae_fit_seed("geo_bld", "fh", "migrated_last_12m", 2L))
})

test_that("a rerun with more outcomes keeps the fitted rows and adds the others", {
  out <- withr::local_tempdir(); inp <- fixture_sim_inputs()
  suppressMessages(sim_run_replicate(1L, inp, methods = "fh", outcomes = "migrated_last_12m",
                                     out_dir = out))
  p <- fit_path(out, "fh", "geo_bld", 1L)
  first <- arrow::read_parquet(p)
  suppressMessages(sim_run_replicate(1L, inp, methods = "fh",
                                     outcomes = c("migrated_intl_12m", "migrated_last_12m"),
                                     out_dir = out))
  both <- arrow::read_parquet(p)
  expect_identical(unique(both$outcome), c("migrated_last_12m", "migrated_intl_12m"))
  kept <- both[both$outcome == "migrated_last_12m", ]
  expect_identical(kept$wall_secs, first$wall_secs)
  expect_identical(kept$theta_hat, first$theta_hat)
  expect_equal(sum(both$outcome == "migrated_intl_12m"), length(inp$universe))
})

test_that("a broken direct file is rewritten and complete fits are kept", {
  out <- withr::local_tempdir(); inp <- fixture_sim_inputs()
  suppressMessages(sim_run_replicate(1L, inp, methods = "fh", outcomes = "migrated_last_12m",
                                     out_dir = out))
  dp <- file.path(out, "direct", "rep_0001.parquet")
  d1 <- arrow::read_parquet(dp)
  t1 <- file.mtime(fit_path(out, "fh", "geo_bld", 1L))
  writeBin(raw(10), dp)
  suppressMessages(sim_run_replicate(1L, inp, methods = "fh", outcomes = "migrated_last_12m",
                                     out_dir = out))
  expect_identical(arrow::read_parquet(dp), d1)
  expect_equal(file.mtime(fit_path(out, "fh", "geo_bld", 1L)), t1)
})

test_that("replicates run in forked workers write the same files as one process", {
  skip_on_os("windows")
  inp <- fixture_sim_inputs()
  o1 <- withr::local_tempdir(); o2 <- withr::local_tempdir()
  suppressMessages(sim_run(1:3, inp, methods = "fh", outcomes = "migrated_last_12m", out_dir = o1))
  st <- sim_run(1:3, inp, methods = "fh", outcomes = "migrated_last_12m", out_dir = o2,
                workers = 2L, verbose = FALSE)
  expect_identical(st$status, rep("done", 3L))
  for (r in 1:3) {
    a <- arrow::read_parquet(fit_path(o1, "fh", "geo_bld", r))
    b <- arrow::read_parquet(fit_path(o2, "fh", "geo_bld", r))
    a$wall_secs <- b$wall_secs <- NULL
    expect_identical(a, b)
  }
})

test_that("a replicate that errors is reported and the others are written", {
  inp <- fixture_sim_inputs(); out <- withr::local_tempdir()
  draw <- sim_draw
  local_mocked_bindings(sim_draw = function(rep, inputs) {
    if (rep == 2L) stop("frame unreadable")
    draw(rep, inputs)
  })
  expect_error(suppressMessages(sim_run(1:3, inp, methods = "fh", outcomes = "migrated_last_12m",
                                        out_dir = out)),
               "replicate 2: frame unreadable")
  expect_true(file.exists(fit_path(out, "fh", "geo_bld", 1L)))
  expect_true(file.exists(fit_path(out, "fh", "geo_bld", 3L)))
  expect_false(file.exists(fit_path(out, "fh", "geo_bld", 2L)))
})

test_that("the published settings of each method are fixed in one place", {
  s <- sim_settings(xgb_sub_params("geo_blk_bld"))
  expect_identical(s$fh_arcsin, list(B = 500L))
  expect_identical(s$bayes_fh_beta,
                   list(likelihood = "ExtBeta", chains = 3L, iter = 2000L, backend = "cmdstan"))
  expect_identical(s$megb_sub, list(B = 100L, weight = "census", mse_type = "var"))
  expect_identical(s$xgb_sub$B, 200L)
  expect_identical(s$xgb_sub$weight, "census")
  expect_identical(s$xgb_sub$params, xgb_sub_params("geo_blk_bld"))
  expect_identical(sim_v3_methods, c("fh_arcsin", "bayes_fh_beta", "xgb_sub"))
  # The defaults of the fit functions and of method_args() agree with the settings.
  expect_identical(formals(fit_xgb_sub)$B, s$xgb_sub$B)
  expect_identical(formals(fit_megb_sub)$B, s$megb_sub$B)
  d <- data.frame(area_id = 1L, theta_direct = 0.1, outcome = "migrated_last_12m")
  for (m in c("fh_arcsin", "bayes_fh_beta")) {
    a <- method_args(m, direct = d, aux_in = NULL, aux_oos = NULL, seed = 1L)
    shared <- intersect(names(a), names(s[[m]]))
    expect_identical(a[shared], s[[m]][shared])
  }
})

test_that("Bayes-Beta receives the published arguments and its run columns", {
  inp <- fixture_sim_inputs(); dr <- sim_draw(1L, inp)
  o <- "migrated_intl_12m"
  seen <- new.env()
  local_method("bayes_fh_beta", function(...) {
    seen$args <- list(...)
    fake_fit(seen$args$X$municipality, seen$args$X_oos$municipality,
             tibble::tibble(metric = c("n_divergent", "max_rhat"), value = c(2, 1.004)))
  })
  f <- sim_fit_one("bayes_fh_beta", o, 1L, dr$direct, dr$smp_sub, inp)
  a <- seen$args
  drawn <- sort(unique(dr$direct$area_id[dr$direct$outcome == o]))
  X <- inp$X_area
  expect_identical(a$X, X[X$municipality %in% drawn, , drop = FALSE])
  expect_identical(a$X_oos, X[X$municipality %in% setdiff(inp$universe, drawn), , drop = FALSE])
  expect_identical(a$direct, as.data.frame(dr$direct)[dr$direct$outcome == o, , drop = FALSE])
  expect_identical(a[c("likelihood", "iter", "chains", "seed", "backend")],
                   list(likelihood = "ExtBeta", iter = 2000L, chains = 3L,
                        seed = sae_fit_seed("geo_bld", "bayes_fh_beta", o, 1L),
                        backend = "cmdstan"))
  expect_setequal(names(a), c("direct", "X", "X_oos", "likelihood", "iter", "chains", "seed",
                              "backend"))
  expect_identical(names(f), sim_fit_columns)
  expect_equal(nrow(f), length(inp$universe))
  expect_identical(unique(f$n_divergent), 2)
  expect_identical(unique(f$max_rhat), 1.004)
  expect_true(all(is.na(f$weight)) && all(is.na(f$povmap_sha)) && all(is.na(f$mse_type)))
  expect_identical(unique(f$config), "geo_bld")
  expect_identical(unique(f$seed), sae_fit_seed("geo_bld", "bayes_fh_beta", o, 1L))
  # settings that override the seed: the row records the seed the fit used
  st <- sim_settings(); st$bayes_fh_beta$seed <- 99L
  f2 <- sim_fit_one("bayes_fh_beta", o, 1L, dr$direct, dr$smp_sub, inp, settings = st)
  expect_identical(seen$args$seed, 99L)
  expect_identical(unique(f2$seed), 99L)
})

test_that("FH-arcsin receives the bootstrap size and its seed", {
  inp <- fixture_sim_inputs(); dr <- sim_draw(1L, inp)
  seen <- new.env()
  local_method("fh_arcsin", function(...) {
    seen$args <- list(...)
    fake_fit(seen$args$X$municipality, seen$args$X_oos$municipality,
             tibble::tibble(metric = "sigma2_v", value = 0.003))
  })
  f <- sim_fit_one("fh_arcsin", "pers_migrated_12m", 4L, dr$direct, NULL, inp)
  expect_setequal(names(seen$args), c("direct", "X", "X_oos", "B", "seed"))
  expect_identical(seen$args$B, 500L)
  expect_identical(seen$args$seed, sae_fit_seed("geo_bld", "fh_arcsin", "pers_migrated_12m", 4L))
  expect_identical(unique(f$sigma2_v), 0.003)
})

test_that("the sub-area methods receive the sub-area sample, the matrix and their settings", {
  inp <- fixture_sim_inputs(); dr <- sim_draw(1L, inp)
  o <- "pers_migrated_internal_12m"
  seen <- new.env()
  recorder <- function(tag) function(...) {
    a <- list(...)
    seen[[tag]] <- a
    in_ids <- sort(unique(a$direct$area_id))
    fake_fit(in_ids, sort(setdiff(unique(a$pop_sub$municipality), in_ids)),
             tibble::tibble(metric = "ran_eff_var", value = 0.002))
  }
  local_method("xgb_sub", recorder("xgb_sub"))
  local_method("megb_sub", recorder("megb_sub"))
  fx <- sim_fit_one("xgb_sub", o, 2L, dr$direct, dr$smp_sub, inp)
  fm <- sim_fit_one("megb_sub", o, 2L, dr$direct, dr$smp_sub, inp)
  ax <- seen$xgb_sub; am <- seen$megb_sub
  want_direct <- as.data.frame(dr$direct)[dr$direct$outcome == o, , drop = FALSE]
  want_smp <- dr$smp_sub[dr$smp_sub$outcome == o, , drop = FALSE]
  for (a in list(ax, am)) {
    expect_identical(a$direct, want_direct)
    expect_identical(a$smp_sub, want_smp)
    expect_identical(a$pop_sub, inp$pop_sub)
    expect_identical(a$outcome, o)
    expect_identical(a$weight, "census")
  }
  expect_identical(ax$seed, sae_fit_seed("geo_blk_bld", "xgb_sub", o, 2L))
  expect_identical(am$seed, sae_fit_seed("geo_blk_bld", "megb_sub", o, 2L))
  expect_identical(ax$params, inp$tune[[o]])
  expect_identical(ax$B, 200L)
  expect_setequal(names(ax), c("direct", "smp_sub", "pop_sub", "outcome", "seed", "B", "weight",
                               "params"))
  expect_identical(am$mse_type, "var")
  expect_identical(am$B, 100L)
  expect_setequal(names(am), c("direct", "smp_sub", "pop_sub", "outcome", "seed", "B", "weight",
                               "mse_type"))
  expect_identical(unique(fx$config), "geo_blk_bld")
  expect_identical(unique(fx$weight), "census")
  expect_true(all(is.na(fx$mse_type)))
  expect_identical(unique(fm$mse_type), "var")
  expect_identical(unique(fm$sigma2_v), 0.002)
  expect_equal(nrow(fx), length(inp$universe))
})

test_that("XGB is refused without tuned hyperparameters", {
  inp <- fixture_sim_inputs(); inp$tune <- NULL
  out <- withr::local_tempdir()
  expect_error(sim_run_replicate(1L, inp, methods = "xgb_sub", outcomes = "migrated_last_12m",
                                 out_dir = out), "xgb-params")
  expect_error(sim_run(1L, inp, methods = c("fh", "xgb_sub"), outcomes = "migrated_last_12m",
                       out_dir = out), "xgb_sub_params\\(file")
  expect_false(dir.exists(file.path(out, "direct")))
  dr <- sim_draw(1L, fixture_sim_inputs())
  expect_error(sim_fit_one("xgb_sub", "migrated_last_12m", 1L, dr$direct, dr$smp_sub, inp),
               "hyperparameters")
})

test_that("unknown methods and outcomes are refused before any draw", {
  inp <- fixture_sim_inputs(); out <- withr::local_tempdir()
  expect_error(sim_run(1L, inp, methods = "fh_typo", out_dir = out), "fh_typo")
  expect_error(sim_run(1L, inp, methods = "fh", outcomes = "migrated", out_dir = out), "migrated")
  expect_error(sim_run(0L, inp, methods = "fh", out_dir = out), "replicate")
  expect_false(dir.exists(file.path(out, "direct")))
})

test_that("a fit that reports its own failure records that message", {
  inp <- fixture_sim_inputs(); dr <- sim_draw(1L, inp)
  local_method("xgb_sub", function(direct, smp_sub, pop_sub, outcome, ...) {
    in_ids <- sort(unique(direct$area_id))
    .povmap_failure("xgb_sub", in_ids, sort(setdiff(unique(pop_sub$municipality), in_ids)),
                    Sys.time(), "too few sub-areas for the bootstrap")
  })
  f <- sim_fit_one("xgb_sub", "migrated_last_12m", 1L, dr$direct, dr$smp_sub, inp)
  expect_true(all(f$fit_failed))
  expect_identical(unique(f$error_msg), "too few sub-areas for the bootstrap")
  expect_true(all(is.na(f$theta_hat)))
  expect_true(all(is.na(f$sigma2_v)))
})

test_that("package versions name the toolchain of the fits", {
  pv <- sim_pkg_versions()
  expect_match(pv, "^emdi2=[^;]+;povmap=[^;]+;xgboost=[^;]+;tipsae=[^;]+;cmdstan=[^;]+$")
  if (requireNamespace("povmap", quietly = TRUE)) {
    sha <- utils::packageDescription("povmap")$RemoteSha
    if (!is.null(sha)) expect_match(pv, paste0("povmap=[^;@]+@", substr(sha, 1, 7), ";"))
  }
  cmd <- withr::local_tempdir()
  writeLines(c("# makefile", "CMDSTAN_VERSION := 2.99.1", "other := 1"), file.path(cmd, "makefile"))
  withr::local_envvar(CMDSTAN = cmd)
  expect_match(sim_pkg_versions(), "cmdstan=2.99.1$")
})

test_that("CmdStan is found in the default install folder when CMDSTAN is unset", {
  home <- withr::local_tempdir()
  dir.create(file.path(home, ".cmdstan", "cmdstan-2.9.0"), recursive = TRUE)
  dir.create(file.path(home, ".cmdstan", "cmdstan-2.38.0"), recursive = TRUE)
  withr::local_envvar(CMDSTAN = NA)
  want <- file.path(home, ".cmdstan", "cmdstan-2.38.0")
  expect_identical(sim_cmdstan_env(home = home), want)
  expect_identical(Sys.getenv("CMDSTAN"), want)
  withr::local_envvar(CMDSTAN = "/opt/cmdstan-x")
  expect_identical(sim_cmdstan_env(home = home), "/opt/cmdstan-x")
})

test_that("a missing CmdStan stops a run that needs it and names CMDSTAN", {
  home <- withr::local_tempdir()
  withr::local_envvar(CMDSTAN = NA)
  expect_error(sim_cmdstan_env(home = home, required = TRUE), "CMDSTAN")
  expect_identical(sim_cmdstan_env(home = home), "")
  expect_identical(Sys.getenv("CMDSTAN"), "")
})

test_that("sub-area XGB runs end to end on the fixture", {
  skip_if_not_installed("povmap")
  skip_if_not_installed("xgboost")
  inp <- fixture_sim_inputs(); dr <- sim_draw(1L, inp)
  s <- sim_settings(inp$tune); s$xgb_sub$B <- 5L
  utils::capture.output(f <- suppressWarnings(suppressMessages(
    sim_fit_one("xgb_sub", "migrated_last_12m", 1L, dr$direct, dr$smp_sub, inp, settings = s))))
  expect_false(any(f$fit_failed))
  expect_equal(nrow(f), length(inp$universe))
  expect_true(all(is.finite(f$theta_hat)) && all(f$mse > 0))
  expect_identical(unique(f$povmap_sha), .povmap_sha())
})

test_that("a method whose packages do not load is refused before any draw, naming them", {
  inp <- fixture_sim_inputs(); out <- withr::local_tempdir()
  local_mocked_bindings(.sim_has_pkg = function(pkg) !pkg %in% c("povmap", "emdi2"))
  msg <- "emdi2, povmap; needed by fh_arcsin \\(emdi2\\), xgb_sub \\(povmap\\)"
  expect_error(sim_run(1L, inp, methods = c("fh_arcsin", "xgb_sub"),
                       outcomes = "migrated_last_12m", out_dir = out), msg)
  expect_error(sim_run_replicate(1L, inp, methods = c("fh_arcsin", "xgb_sub"),
                                 outcomes = "migrated_last_12m", out_dir = out), msg)
  expect_identical(list.files(out, recursive = TRUE), character(0))
  s <- sim_settings(inp$tune)
  expect_identical(.sim_method_packages("bayes_fh_beta", s), c("tipsae", "rstan", "cmdstanr"))
  s$bayes_fh_beta$backend <- "rstan"
  expect_identical(.sim_method_packages("bayes_fh_beta", s), c("tipsae", "rstan", "posterior"))
  expect_identical(.sim_method_packages("megb_sub", s), c("povmap", "xgboost"))
  expect_identical(.sim_method_packages("fh", s), "emdi2")
})

test_that("a Stan model that does not compile stops the run before anything is written", {
  inp <- fixture_sim_inputs(); out <- withr::local_tempdir()
  withr::local_envvar(CMDSTAN = withr::local_tempdir())
  local_mocked_bindings(.sim_has_pkg = function(pkg) TRUE,
                        .tipsae_cmdstan_model = function() stop("g++: command not found"))
  expect_error(sim_run(1:2, inp, methods = c("fh", "bayes_fh_beta"),
                       outcomes = "migrated_last_12m", out_dir = out),
               "bayes_fh_beta could not be compiled: g\\+\\+: command not found")
  expect_identical(list.files(out, recursive = TRUE), character(0))
})

test_that("the Stan model is compiled once for a range of replicates", {
  inp <- fixture_sim_inputs(); out <- withr::local_tempdir()
  withr::local_envvar(CMDSTAN = withr::local_tempdir())
  n <- new.env(); n$calls <- 0L
  local_mocked_bindings(.sim_has_pkg = function(pkg) TRUE,
                        .tipsae_cmdstan_model = function() n$calls <- n$calls + 1L)
  local_method("bayes_fh_beta", function(...) {
    a <- list(...)
    fake_fit(a$X$municipality, a$X_oos$municipality)
  })
  st <- suppressMessages(sim_run(1:3, inp, methods = "bayes_fh_beta",
                                 outcomes = "migrated_last_12m", out_dir = out))
  expect_identical(n$calls, 1L)
  expect_identical(st$status, rep("done", 3L))
})

test_that("failed fits are counted per method in the messages and in the run's result", {
  inp <- fixture_sim_inputs(); out <- withr::local_tempdir()
  local_method("fh", function(...) stop("singular"))
  oc <- c("migrated_last_12m", "migrated_intl_12m")
  msgs <- capture_messages(sim_run_replicate(1L, inp, methods = c("fh", "fh_arcsin"),
                                             outcomes = oc, out_dir = out))
  expect_true(any(grepl("^replicate 1: 4 fit\\(s\\) in .*; failed fits: fh 2, fh_arcsin 0", msgs)))
  st <- suppressMessages(sim_run(1:2, inp, methods = c("fh_arcsin", "fh"), outcomes = oc,
                                 out_dir = out))
  expect_identical(names(st), c("replicate", "status", "message", "failed_fh", "failed_fh_arcsin"))
  expect_identical(st$status, c("skipped", "done"))
  expect_identical(st$failed_fh, c(2L, 2L))
  expect_identical(st$failed_fh_arcsin, c(0L, 0L))
  msgs <- capture_messages(sim_run(1L, inp, methods = c("fh", "fh_arcsin"), outcomes = oc,
                                   out_dir = out))
  expect_true(any(grepl("^replicate 1: complete, skipped; failed fits: fh 2, fh_arcsin 0", msgs)))
})

test_that("failed fits are kept on a rerun and fitted again with retry_failed", {
  inp <- fixture_sim_inputs(); out <- withr::local_tempdir()
  oc <- c("migrated_last_12m", "migrated_intl_12m")
  real <- sae_methods$fh
  reg <- sae_methods
  reg$fh <- function(...) {
    a <- list(...)
    if (all(a$direct$outcome == "migrated_intl_12m")) stop("singular") else real(...)
  }
  with_mocked_bindings(
    suppressMessages(sim_run_replicate(1L, inp, methods = "fh", outcomes = oc, out_dir = out)),
    sae_methods = reg)
  p <- fit_path(out, "fh", "geo_bld", 1L)
  f1 <- arrow::read_parquet(p); t1 <- file.mtime(p)
  expect_true(all(f1$fit_failed[f1$outcome == "migrated_intl_12m"]))
  expect_false(any(f1$fit_failed[f1$outcome == "migrated_last_12m"]))
  st <- suppressMessages(sim_run(1L, inp, methods = "fh", outcomes = oc, out_dir = out))
  expect_identical(st$status, "skipped")
  expect_identical(st$failed_fh, 1L)
  expect_equal(file.mtime(p), t1)
  st <- suppressMessages(sim_run(1L, inp, methods = "fh", outcomes = oc, out_dir = out,
                                 retry_failed = TRUE))
  f2 <- arrow::read_parquet(p)
  expect_identical(st$failed_fh, 0L)
  expect_false(any(f2$fit_failed))
  expect_equal(nrow(f2), 2L * length(inp$universe))
  expect_identical(f2$wall_secs[f2$outcome == "migrated_last_12m"],
                   f1$wall_secs[f1$outcome == "migrated_last_12m"])
})

test_that("retry_failed refits failed outcomes that are requested and keeps the others", {
  inp <- fixture_sim_inputs(); out <- withr::local_tempdir()
  oc <- c("migrated_last_12m", "migrated_intl_12m")
  reg <- sae_methods
  reg$fh <- function(...) stop("singular")
  with_mocked_bindings(
    suppressMessages(sim_run_replicate(1L, inp, methods = "fh", outcomes = oc, out_dir = out)),
    sae_methods = reg)
  p <- fit_path(out, "fh", "geo_bld", 1L)
  expect_true(all(arrow::read_parquet(p)$fit_failed))
  suppressMessages(sim_run(1L, inp, methods = "fh", outcomes = "migrated_last_12m",
                           out_dir = out, retry_failed = TRUE))
  f <- arrow::read_parquet(p)
  expect_setequal(unique(f$outcome), oc)
  expect_false(any(f$fit_failed[f$outcome == "migrated_last_12m"]))
  expect_true(all(f$fit_failed[f$outcome == "migrated_intl_12m"]))
})

test_that("hyperparameters given as a list hash alike whatever file they were read from", {
  p <- xgb_sub_params()
  expect_identical(.sim_md5_tune(structure(p, source = "a.tsv")),
                   .sim_md5_tune(structure(p, source = "b.tsv")))
  q <- p; q$migrated_last_12m$eta <- 0.5
  expect_false(identical(.sim_md5_tune(p), .sim_md5_tune(q)))
})

test_that("a folder refuses results computed from other inputs or settings", {
  inp <- fixture_sim_inputs(); out <- withr::local_tempdir()
  o <- "migrated_last_12m"
  suppressMessages(sim_run_replicate(1L, inp, methods = "fh", outcomes = o, out_dir = out))
  m <- jsonlite::fromJSON(file.path(out, "run_manifest.json"), simplifyVector = FALSE)
  expect_identical(m$area_config, "geo_bld")
  expect_identical(m$subarea_config, "geo_blk_bld")
  expect_identical(m$settings$fh_arcsin$B, 500L)
  expect_true(all(c("frame", "stratum_map", "area_matrix", "subarea_matrix", "extra_area",
                    "extra_subarea", "spatial_weights", "xgb_params") %in% names(m$sources)))
  expect_identical(m$sources$area_matrix$md5, inp$sources$area_matrix$md5)
  expect_match(m$sources$xgb_params$path, "xgb_sub_params_geo_blk_bld\\.tsv$")
  expect_identical(m$pkg_versions$tipsae, as.character(utils::packageVersion("tipsae")))
  # The same inputs read again (other temporary paths, same contents) are accepted.
  suppressMessages(sim_run_replicate(2L, fixture_sim_inputs(), methods = "fh", outcomes = o,
                                     out_dir = out))
  inp2 <- inp; inp2$X_area$share_urban[1] <- 0.5
  expect_error(sim_run_replicate(3L, inp2, methods = "fh", outcomes = o, out_dir = out),
               "differs in: data.area_matrix_md5")
  s <- sim_settings(inp$tune); s$fh_arcsin$B <- 50L
  expect_error(sim_run(3L, inp, methods = "fh", outcomes = o, out_dir = out, settings = s),
               "differs in: settings")
  inp3 <- inp; inp3$sources$extra_area <- list(path = "x.csv", size = 10, md5 = "abc")
  expect_error(sim_run(3L, inp3, methods = "fh", outcomes = o, out_dir = out),
               "differs in: sources.extra_area")
  expect_false(file.exists(fit_path(out, "fh", "geo_bld", 3L)))
  # A package version that is not known is not compared.
  local_mocked_bindings(sim_pkg_versions = function() {
    "emdi2=NA;povmap=NA;xgboost=NA;tipsae=NA;cmdstan=NA"
  })
  suppressMessages(sim_run(3L, inp, methods = "fh", outcomes = o, out_dir = out))
  expect_true(file.exists(fit_path(out, "fh", "geo_bld", 3L)))
})

test_that("old temporary files of a target are removed when it is written, recent ones kept", {
  inp <- fixture_sim_inputs(); out <- withr::local_tempdir()
  dir.create(file.path(out, "fits", "fh__geo_bld"), recursive = TRUE)
  dir.create(file.path(out, "direct"))
  stale <- c(paste0(fit_path(out, "fh", "geo_bld", 1L), ".tmp.12345"),
             file.path(out, "direct", "rep_0001.parquet.tmp.999"))
  # A recent temporary file of the same target may belong to a process that is
  # still writing it.
  live <- paste0(fit_path(out, "fh", "geo_bld", 1L), ".tmp.777")
  other <- paste0(fit_path(out, "fh", "geo_bld", 2L), ".tmp.12345")
  file.create(c(stale, live, other))
  Sys.setFileTime(c(stale, other), Sys.time() - 2 * 3600)
  suppressMessages(sim_run_replicate(1L, inp, methods = "fh", outcomes = "migrated_last_12m",
                                     out_dir = out))
  expect_false(any(file.exists(stale)))
  expect_true(file.exists(live))
  expect_true(file.exists(other))
})

test_that("the script's arguments are read into the runner's", {
  a <- sim_parse_args(c("--reps=1-3,7,9-10", "--out=res"))
  expect_identical(a$reps, c(1L, 2L, 3L, 7L, 9L, 10L))
  expect_identical(a$out, "res")
  expect_identical(a$methods, sim_v3_methods)
  expect_identical(a$outcomes, sae_outcomes$outcome)
  expect_identical(a$workers, 1L)
  expect_false(a$retry_failed)
  expect_null(a$xgb_params)
  expect_null(a$area_matrix)
  f <- withr::local_tempfile(fileext = ".tsv"); writeLines("outcome", f)
  b <- sim_parse_args(c("--reps=5", "--out=res", "--methods=fh,xgb_sub",
                        "--outcomes=migrated_intl_12m", "--workers=4",
                        paste0("--xgb-params=", f), paste0("--extra-area=", f), "--retry-failed"))
  expect_identical(b$reps, 5L)
  expect_identical(b$methods, c("fh", "xgb_sub"))
  expect_identical(b$outcomes, "migrated_intl_12m")
  expect_identical(b$workers, 4L)
  expect_identical(b$xgb_params, f)
  expect_identical(b$extra_area, f)
  expect_true(b$retry_failed)
})

test_that("malformed script arguments are refused", {
  base <- c("--reps=1", "--out=res")
  expect_error(sim_parse_args(character(0)), "no arguments")
  expect_error(sim_parse_args("--reps=1-3"), "--out is required")
  expect_error(sim_parse_args(c(base, "--colour=red")), "unknown option\\(s\\): --colour")
  expect_error(sim_parse_args(c(base, "res")), "not of the form")
  expect_error(sim_parse_args(c(base, "--reps=2")), "more than once")
  expect_error(sim_parse_args(c("--reps=3-1", "--out=res")), "decreasing range 3-1")
  expect_error(sim_parse_args(c("--reps=1-a", "--out=res")), "cannot read 1-a")
  expect_error(sim_parse_args(c("--reps=0-2", "--out=res")), "1 to 9999")
  expect_error(sim_parse_args(c("--reps=1-999999999", "--out=res")), "1 to 9999")
  expect_error(sim_parse_args(c(base, "--workers=0")), "--workers")
  expect_error(sim_parse_args(c(base, "--methods=fh,nope")), "unknown nope")
  expect_error(sim_parse_args(c(base, "--outcomes=migrated")), "unknown migrated")
  expect_error(sim_parse_args(c(base, "--xgb-params=/no/such/file.tsv")),
               "--xgb-params: file not found: /no/such/file.tsv")
  expect_error(sim_parse_args(c(base, "--retry-failed=yes")), "takes no value")
})
