eval_fit_path <- function(run, method, config, rep) {
  file.path(run, "fits", paste0(method, "__", config), sprintf("rep_%04d.parquet", rep))
}

# Rewrite one per-fit file of a fixture run through `fun(rows)`.
edit_fit <- function(run, method, config, rep, fun) {
  p <- eval_fit_path(run, method, config, rep)
  arrow::write_parquet(fun(as.data.frame(arrow::read_parquet(p))), p)
}

# Record the fit of outcome `o` in a per-fit file as failed, as the runner does:
# NA rows of the drawn municipalities with `fit_failed = TRUE`.
fail_outcome <- function(o) {
  function(d) {
    d <- d[d$outcome != o | d$area_set == "in_sample", ]
    f <- d$outcome == o
    d[f, c("theta_hat", "mse", "ci_lo", "ci_hi")] <- NA_real_
    d$fit_failed[f] <- TRUE
    d$error_msg[f] <- "singular fit"
    d
  }
}

# Evaluate `expr` and return its value, expecting exactly the warnings that
# `patterns` describe: each pattern matches at least one warning, and every
# warning matches a pattern.
expect_warnings_matching <- function(expr, patterns) {
  msgs <- character(0)
  val <- withCallingHandlers(expr, warning = function(w) {
    msgs <<- c(msgs, conditionMessage(w))
    invokeRestart("muffleWarning")
  })
  for (p in patterns) expect_true(any(grepl(p, msgs)), label = paste0("a warning matching '", p, "'"))
  hit <- Reduce(`|`, lapply(patterns, grepl, x = msgs), rep(FALSE, length(msgs)))
  expect_identical(msgs[!hit], character(0))
  val
}

UNEVEN <- "methods keep different numbers of fits of the same outcome"
NOT_CARRIED <- "size not carried by the first method"

test_that("external duplicates are refused", {
  e <- fixture_external(); e <- rbind(e, e[1, ])
  expect_error(sim_read_external(e, "MEGB"), "duplicate")
})

test_that("external replicates must match the evaluated ones", {
  run <- fixture_run_dir(reps = 1:2)
  expect_error(sim_evaluate(run, fixture_truth(), reps = 1:2,
                            external = fixture_external(reps = 1:3)), "replicates")
})

test_that("a fit without any drawn municipality stops the evaluation", {
  run <- fixture_run_dir(reps = 1:2); e <- fixture_external(reps = 1:2)
  e$sampled[e$replicate == 2 & e$outcome == "migrated_intl_12m"] <- FALSE
  expect_error(sim_evaluate(run, fixture_truth(), reps = 1:2, external = e),
               "1 fit\\(s\\) have no drawn municipality \\(first: MEGB, migrated_intl_12m, replicate 2\\)")
})

test_that("different drawn municipalities across methods stop the evaluation", {
  run <- fixture_run_dir(reps = 1:2); e <- fixture_external(reps = 1:2)
  e$sampled[1] <- !e$sampled[1]
  expect_error(sim_evaluate(run, fixture_truth(), reps = 1:2, external = e), "drawn")
})

test_that("outputs carry no version column and the documented keys", {
  res <- sim_evaluate(fixture_run_dir(reps = 1:2), fixture_truth(), reps = 1:2)
  expect_false("version" %in% names(res$area_summary))
  expect_true(all(c("method", "config", "outcome", "area_set", "column") %in% names(res$set_summary)))
})

test_that("the external file maps to the per-fit layout", {
  e <- fixture_external(reps = 1)
  x <- sim_read_external(e, "MEGB")
  expect_identical(names(x), c(sim_fit_columns, "truth"))
  expect_identical(x$area_id, as.integer(e$municipality))
  expect_identical(x$area_set, ifelse(e$sampled, "in_sample", "out_of_sample"))
  expect_identical(x$theta_hat, e$estimate)
  expect_identical(x$mse, e$estimated_variance)
  expect_identical(c(x$ci_lo, x$ci_hi), c(e$lower, e$upper))
  expect_identical(x$truth, e$truth)
  expect_true(all(!x$fit_failed))
  expect_identical(unique(x$method), "MEGB")
  expect_identical(unique(x$config), "external")
  expect_identical(unique(sim_read_external(e, "M2", config = "external_cfg")$config),
                   "external_cfg")
})

test_that("an external CSV file reads as the same table", {
  e <- fixture_external(reps = 1)
  f <- withr::local_tempfile(fileext = ".csv")
  data.table::fwrite(e, f)
  expect_equal(sim_read_external(f, "MEGB"), sim_read_external(e, "MEGB"))
})

test_that("malformed external files are refused, naming the problem", {
  e <- fixture_external(reps = 1)
  expect_error(sim_read_external(e[, names(e) != "upper"], "MEGB"), "missing column\\(s\\): upper")
  bad <- e; bad$estimate[2] <- NA
  expect_error(sim_read_external(bad, "MEGB"), "1 row has a missing value")
  bad <- e; bad$estimated_variance[3] <- -1e-6
  expect_error(sim_read_external(bad, "MEGB"), "negative variance")
  bad <- e; bad$lower[4] <- bad$upper[4] + 1
  expect_error(sim_read_external(bad, "MEGB"), "lower bound above")
  bad <- e; bad$outcome[5] <- "migrated"
  expect_error(sim_read_external(bad, "MEGB"), "unknown outcome")
  bad <- e; bad$sampled <- "yes"
  expect_error(sim_read_external(bad, "MEGB"), "`sampled`")
  bad <- e; bad$municipality[6] <- "99001"
  expect_error(sim_read_external(bad, "MEGB", universe = fixture_eval_universe()),
               "1 municipality code is not in the universe")
  expect_no_error(sim_read_external(bad, "MEGB"))
})

test_that("the external truth must equal the truth file", {
  run <- fixture_run_dir(reps = 1:2); e <- fixture_external(reps = 1:2)
  e$truth[10] <- e$truth[10] + 1e-6
  expect_error(sim_evaluate(run, fixture_truth(), reps = 1:2, external = e),
               "truth .* differs from the truth file by up to 1e-06")
})

test_that("every table is returned, keyed and labelled with the covariate set", {
  res <- sim_evaluate(fixture_run_dir(reps = 1:2), fixture_truth(), reps = 1:2,
                      external = fixture_external(reps = 1:2), column = "geo")
  expect_identical(names(res), c("area_summary", "set_summary", "import_check", "fit_check", "by_bin",
                                 "cv_flags", "mcse", "reliability", "ratio_summary",
                                 "tracking_summary", "size_medians", "coverage_summary"))
  for (k in names(res)) {
    expect_true(data.table::is.data.table(res[[k]]), label = k)
    expect_identical(unique(res[[k]]$column), "geo", label = k)
  }
  for (k in c("area_summary", "set_summary", "by_bin", "cv_flags", "mcse", "reliability")) {
    expect_true(all(UR_KEYS %in% names(res[[k]])), label = k)
  }
  s <- res$set_summary
  expect_setequal(s$method, c("FH-arcsin", "Bayes-Beta", "XGB", "MEGB"))
  expect_identical(unique(s$config[s$method == "FH-arcsin"]), "geo_bld")
  expect_identical(unique(s$config[s$method == "XGB"]), "geo_blk_bld")
  expect_identical(unique(s$config[s$method == "MEGB"]), "external")
  expect_setequal(s$outcome, sae_outcomes$outcome)
  expect_identical(nrow(s), 4L * 6L * 2L)
})

test_that("the measures are those of the per-fit rows against the truth", {
  run <- fixture_run_dir(reps = 1:2)
  tr <- fixture_truth()
  res <- sim_evaluate(run, tr, reps = 1:2)
  rows <- do.call(rbind, lapply(1:2, function(r) {
    as.data.frame(arrow::read_parquet(eval_fit_path(run, "fh_arcsin", "geo_bld", r)))
  }))
  rows <- merge(rows[rows$outcome == "migrated_intl_12m", ], tr, by = c("area_id", "outcome"))
  for (s in c("in_sample", "out_of_sample")) {
    x <- rows[rows$area_set == s, ]
    got <- res$set_summary[method == "FH-arcsin" & outcome == "migrated_intl_12m" & area_set == s]
    expect_equal(got$n_rows, nrow(x))
    expect_equal(got$rmse, sqrt(mean((x$theta_hat - x$theta_true)^2)))
    expect_equal(got$coverage, mean(x$ci_lo <= x$theta_true & x$theta_true <= x$ci_hi))
    expect_equal(got$ratio_pooled, mean(x$mse) / mean((x$theta_hat - x$theta_true)^2))
  }
  # Dwellings drawn come from the direct files; out of sample they are 0.
  dr <- do.call(rbind, lapply(1:2, function(r) {
    d <- as.data.frame(arrow::read_parquet(file.path(run, "direct", sprintf("rep_%04d.parquet", r))))
    d$replicate <- r
    d
  }))
  dr <- dr[dr$outcome == "migrated_intl_12m", ]
  a <- res$area_summary[method == "XGB" & outcome == "migrated_intl_12m"]
  want <- tapply(dr$n_dwell, dr$area_id, stats::median)
  ins <- a[area_set == "in_sample"]
  expect_equal(ins$n_dwell_med, unname(as.numeric(want[as.character(ins$area_id)])))
  expect_true(all(a[area_set == "out_of_sample"]$n_dwell_med == 0))
  expect_identical(res$import_check$rows_kept, 3L * 2L * 6L * 60L)
  expect_true(is.na(res$import_check$source))
})

test_that("the import check counts the external rows and the truth difference", {
  run <- fixture_run_dir(reps = 1:2)
  f <- withr::local_tempfile(fileext = ".parquet")
  arrow::write_parquet(fixture_external(reps = 1:2), f)
  ic <- sim_evaluate(run, fixture_truth(), reps = 1:2, external = f)$import_check
  expect_identical(ic$source, basename(f))
  expect_identical(ic$rows, 2L * 6L * 60L)
  expect_identical(ic$truth_max_abs_diff, 0)
  expect_identical(ic$rows_kept, 4L * 2L * 6L * 60L)
  expect_identical(c(ic$fits_dropped, ic$fits_dropped_failed, ic$fits_dropped_rhat,
                     ic$rows_nonfinite), c(0L, 0L, 0L, 0L))
})

test_that("failed and non-converged fits are dropped whole and counted", {
  run <- fixture_run_dir(reps = 1:2)
  edit_fit(run, "xgb_sub", "geo_blk_bld", 2L, fail_outcome("migrated_last_12m"))
  edit_fit(run, "bayes_fh_beta", "geo_bld", 1L, function(d) {
    d$max_rhat[d$outcome == "pers_migrated_12m"] <- 1.2
    d
  })
  edit_fit(run, "fh_arcsin", "geo_bld", 1L, function(d) {
    d$mse[d$outcome == "migrated_intl_12m"][1:2] <- Inf
    d
  })
  res <- expect_warnings_matching(sim_evaluate(run, fixture_truth(), reps = 1:2),
                                  c(UNEVEN, NOT_CARRIED))
  ic <- res$import_check
  expect_identical(c(ic$fits_dropped, ic$fits_dropped_failed, ic$fits_dropped_rhat), c(2L, 1L, 1L))
  expect_identical(ic$rows_nonfinite, 2L)
  s <- res$set_summary
  expect_identical(sum(s[method == "XGB" & outcome == "migrated_last_12m"]$n_rows), 60L)
  expect_identical(sum(s[method == "Bayes-Beta" & outcome == "pers_migrated_12m"]$n_rows), 60L)
  expect_identical(sum(s[method == "FH-arcsin" & outcome == "migrated_intl_12m"]$n_rows), 118L)
  fc <- res$fit_check
  expect_identical(nrow(fc), 3L * 6L)
  x <- fc[method == "XGB" & outcome == "migrated_last_12m"]
  expect_identical(c(x$fits_read, x$failed, x$dropped_rhat, x$dropped_nonfinite, x$kept),
                   c(2L, 1L, 0L, 0L, 1L))
  x <- fc[method == "Bayes-Beta" & outcome == "pers_migrated_12m"]
  expect_identical(c(x$failed, x$dropped_rhat, x$kept), c(0L, 1L, 1L))
  x <- fc[method == "FH-arcsin" & outcome == "migrated_intl_12m"]
  expect_identical(c(x$kept, x$rows_nonfinite), c(2L, 2L))
  expect_identical(sum(fc$kept), 3L * 6L * 2L - 2L)
})

test_that("uneven numbers of fits are warned about once, naming methods and counts", {
  run <- fixture_run_dir(reps = 1:3)
  edit_fit(run, "xgb_sub", "geo_blk_bld", 3L, fail_outcome("migrated_intl_12m"))
  edit_fit(run, "fh_arcsin", "geo_bld", 1L, fail_outcome("migrated_intl_12m"))
  edit_fit(run, "fh_arcsin", "geo_bld", 2L, fail_outcome("pers_migrated_12m"))
  msgs <- character(0)
  withCallingHandlers(sim_evaluate(run, fixture_truth(), reps = 1:3), warning = function(w) {
    msgs <<- c(msgs, conditionMessage(w))
    invokeRestart("muffleWarning")
  })
  uneven <- msgs[grepl(UNEVEN, msgs)]
  expect_length(uneven, 1L)
  expect_match(uneven, "migrated_intl_12m \\(FH-arcsin 2, Bayes-Beta 3, XGB 2\\)")
  expect_match(uneven, "pers_migrated_12m \\(FH-arcsin 2, Bayes-Beta 3, XGB 3\\)")
  expect_true(all(grepl(paste(UNEVEN, NOT_CARRIED, sep = "|"), msgs)))
})

test_that("a method that keeps no fit of an outcome stops the evaluation, saying why", {
  run <- fixture_run_dir(reps = 1:2)
  for (r in 1:2) edit_fit(run, "xgb_sub", "geo_blk_bld", r, fail_outcome("migrated_last_12m"))
  expect_error(sim_evaluate(run, fixture_truth(), reps = 1:2),
               "XGB keeps no fit of migrated_last_12m \\(of 2 fits: 2 failed, 0 with a split R-hat")
  run <- fixture_run_dir(reps = 1:2)
  for (r in 1:2) {
    edit_fit(run, "bayes_fh_beta", "geo_bld", r, function(d) {
      d$max_rhat[d$outcome == "migrated_intl_12m"] <- 1.3
      d
    })
  }
  expect_error(sim_evaluate(run, fixture_truth(), reps = 1:2),
               "Bayes-Beta keeps no fit of migrated_intl_12m .*0 failed, 2 with a split R-hat above 1.05")
  run <- fixture_run_dir(reps = 1:2)
  edit_fit(run, "fh_arcsin", "geo_bld", 1L, fail_outcome("pers_migrated_12m"))
  edit_fit(run, "fh_arcsin", "geo_bld", 2L, function(d) {
    d$theta_hat[d$outcome == "pers_migrated_12m"] <- Inf
    d
  })
  expect_error(sim_evaluate(run, fixture_truth(), reps = 1:2),
               "FH-arcsin keeps no fit of pers_migrated_12m \\(of 2 fits: 1 failed, 0 .*, 1 without a finite row\\)")
  e <- fixture_external(reps = 1:2)
  e$estimate[e$outcome == "migrated_internal_12m"] <- Inf
  expect_error(sim_evaluate(fixture_run_dir(reps = 1:2), fixture_truth(), reps = 1:2, external = e),
               "MEGB keeps no fit of migrated_internal_12m .*2 without a finite row")
})

test_that("a fit without one row per municipality stops the evaluation", {
  run <- fixture_run_dir(reps = 1:2)
  edit_fit(run, "fh_arcsin", "geo_bld", 2L, function(d) d[-1, ])
  expect_error(sim_evaluate(run, fixture_truth(), reps = 1:2),
               "fh_arcsin.*without one row per municipality.*replicate 2, 59 rows for 60")
  run <- fixture_run_dir(reps = 1:2)
  edit_fit(run, "fh_arcsin", "geo_bld", 2L, function(d) { d$area_id[1] <- d$area_id[2]; d })
  expect_error(sim_evaluate(run, fixture_truth(), reps = 1:2), "one row per municipality")
})

test_that("missing replicates, outcomes and direct files stop the evaluation", {
  run <- fixture_run_dir(reps = 1:2)
  expect_error(sim_evaluate(run, fixture_truth(), reps = 1:3),
               "xgb_sub__geo_blk_bld: no file for replicate 3")
  run <- fixture_run_dir(reps = 1:2)
  edit_fit(run, "xgb_sub", "geo_blk_bld", 1L, function(d) d[d$outcome != "migrated_intl_12m", ])
  expect_error(sim_evaluate(run, fixture_truth(), reps = 1:2),
               "xgb_sub.*no fit for 1 of the \\(outcome, replicate\\).*migrated_intl_12m, replicate 1")
  run <- fixture_run_dir(reps = 1:2)
  unlink(file.path(run, "direct", "rep_0002.parquet"))
  expect_error(sim_evaluate(run, fixture_truth(), reps = 1:2), "direct.*replicate 2")
  run <- fixture_run_dir(reps = 1:2)
  edit_fit(run, "fh_arcsin", "geo_bld", 1L, function(d) { d$replicate <- 2L; d })
  expect_error(sim_evaluate(run, fixture_truth(), reps = 1:2), "rows of another")
})

test_that("methods selects the run's methods, and a label taken twice is refused", {
  run <- fixture_run_dir(reps = 1:2, methods = c(sim_v3_methods, "megb_sub"))
  res <- sim_evaluate(run, fixture_truth(), reps = 1:2, methods = c("xgb_sub", "fh_arcsin"))
  expect_setequal(res$set_summary$method, c("XGB", "FH-arcsin"))
  expect_error(sim_evaluate(run, fixture_truth(), reps = 1:2, external = fixture_external()),
               "label MEGB .*megb_sub")
  expect_error(sim_evaluate(run, fixture_truth(), reps = 1:2, methods = "nope"),
               "no per-fit results for method nope")
  res <- sim_evaluate(run, fixture_truth(), reps = 1:2, methods = sim_v3_methods,
                      external = fixture_external())
  expect_setequal(res$set_summary$method, c("FH-arcsin", "Bayes-Beta", "XGB", "MEGB"))
  dir.create(file.path(run, "fits", "fh_arcsin__geo"))
  expect_error(sim_evaluate(run, fixture_truth(), reps = 1:2), "more than one configuration")
})

test_that("without methods, the run's methods come in the order of the display labels", {
  res <- sim_evaluate(fixture_run_dir(reps = 1:2), fixture_truth(), reps = 1:2,
                      external = fixture_external())
  expect_identical(unique(res$area_summary$method), c("FH-arcsin", "Bayes-Beta", "XGB", "MEGB"))
})

test_that("the display summaries are computed for every outcome", {
  res <- sim_evaluate(fixture_run_dir(reps = 1:2), fixture_truth(), reps = 1:2,
                      external = fixture_external(reps = 1:2))
  a <- res$area_summary
  for (k in c("ratio_summary", "tracking_summary", "size_medians", "coverage_summary")) {
    expect_setequal(res[[k]]$outcome, sae_outcomes$outcome)
    expect_type(res[[k]]$method, "character")
  }
  o <- "migrated_internal_12m"
  got <- res$ratio_summary[outcome == o, setdiff(names(res$ratio_summary), c("outcome", "column")),
                           with = FALSE]
  expect_equal(got[order(method, area_set)], ur_ratio_summary(a, o)[order(method, area_set)])
  got <- res$coverage_summary[outcome == o]
  expect_equal(got[order(method)]$overall, ur_coverage_summary(a, o)[order(method)]$overall)
  expect_setequal(res$size_medians$kind, c("Empirical RMSE", "Estimated RMSE"))
})

test_that("sample-size groups do not depend on the order in which methods are named", {
  run <- fixture_run_dir(reps = 1:2)
  # A failed XGB fit gives XGB other median sizes than the other methods.
  edit_fit(run, "xgb_sub", "geo_blk_bld", 1L, fail_outcome("migrated_internal_12m"))
  a <- expect_warnings_matching(
    sim_evaluate(run, fixture_truth(), reps = 1:2, methods = c("xgb_sub", "fh_arcsin", "bayes_fh_beta")),
    c(UNEVEN, NOT_CARRIED))
  b <- expect_warnings_matching(
    sim_evaluate(run, fixture_truth(), reps = 1:2, methods = c("fh_arcsin", "bayes_fh_beta", "xgb_sub")),
    c(UNEVEN, NOT_CARRIED))
  expect_equal(a$size_medians[order(outcome, method, kind, size_group)],
               b$size_medians[order(outcome, method, kind, size_group)])
})

test_that("bins hold the sample sizes, out of sample first, and the truth quintiles", {
  res <- sim_evaluate(fixture_run_dir(reps = 1:2), fixture_truth(), reps = 1:2)
  b <- res$by_bin
  expect_setequal(b$bin_var, c("size", "prevalence"))
  sz <- b[bin_var == "size"]
  expect_true(all(sz[bin == "0 (out of sample)"]$area_set == "out_of_sample"))
  expect_true(all(sz[bin == "0 (out of sample)"]$bin_order == 1L))
  expect_true(all(sz[bin != "0 (out of sample)"]$area_set == "in_sample"))
  expect_true(all(sz[bin != "0 (out of sample)"]$bin_order > 1L))
  expect_lte(max(sz$bin_order), 9L)
  pv <- b[bin_var == "prevalence"]
  expect_setequal(pv$bin_order, 1:5)
  expect_true(all(is.na(pv$x_med)))
  # Each bin variable partitions the rows of a key.
  n_rows <- res$set_summary[, .(n = sum(n_rows)), by = .(method, outcome)]
  for (v in c("size", "prevalence")) {
    n_bin <- b[bin_var == v, .(n = sum(n)), by = .(method, outcome)]
    expect_equal(n_bin[order(method, outcome)]$n, n_rows[order(method, outcome)]$n)
  }
})

test_that("the evaluation script's arguments are read", {
  run <- withr::local_tempdir()
  tr <- withr::local_tempfile(fileext = ".parquet"); writeLines("x", tr)
  a <- sim_parse_eval_args(c(paste0("--run=", run), paste0("--truth=", tr),
                             "--reps=1-3,5", "--out=res"))
  expect_identical(a$run, run)
  expect_identical(a$truth, tr)
  expect_identical(a$reps, c(1L, 2L, 3L, 5L))
  expect_identical(a$out, "res")
  expect_null(a$external)
  expect_identical(a$external_label, "MEGB")
  expect_null(a$methods)
  expect_identical(a$column, "geo")
  b <- sim_parse_eval_args(c(paste0("--run=", run), paste0("--truth=", tr), "--reps=1-100",
                             "--out=res", paste0("--external=", tr), "--external-label=M2",
                             "--methods=fh_arcsin,xgb_sub", "--column=set_a"))
  expect_identical(b$external, tr)
  expect_identical(b$external_label, "M2")
  expect_identical(b$methods, c("fh_arcsin", "xgb_sub"))
  expect_identical(b$column, "set_a")
  base <- c(paste0("--run=", run), paste0("--truth=", tr), "--reps=1-2", "--out=res")
  expect_error(sim_parse_eval_args(character(0)), "no arguments")
  expect_error(sim_parse_eval_args(base[-4]), "--out is required")
  expect_error(sim_parse_eval_args(base[-3]), "--reps is required")
  expect_error(sim_parse_eval_args(c(base, "--colour=red")), "unknown option\\(s\\): --colour")
  expect_error(sim_parse_eval_args(c(base, "--reps=3")), "more than once")
  expect_error(sim_parse_eval_args(c(base[-3], "--reps=3-1")), "decreasing range 3-1")
  expect_error(sim_parse_eval_args(c(base[-1], "--run=/no/such/dir")),
               "--run: folder not found: /no/such/dir")
  expect_error(sim_parse_eval_args(c(base, "--external=/no/such/file.parquet")),
               "--external: file not found")
  # `--external-label` alone is refused, not read as `--external`.
  expect_error(sim_parse_eval_args(c(base, "--external-label=M2")),
               "--external-label names the method of --external, which is not given")
  expect_null(sim_parse_eval_args(c(base, "--column=set_a"))$external)
})
