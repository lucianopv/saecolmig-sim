# tests/testthat/test-direct-validation.R

tiny_design <- function() {
  structure(list(
    stage1 = list(unit = "municipality", selection = "pps_wor",
                  size_variable = "n_dwellings"),
    stage2 = list(unit = "segment", selection = "systematic_fractional",
                  target_segment_size = 10L, n_per_psu = 2L),
    stage3 = list(unit = "household", selection = "census"),
    strata = data.frame(stratum = c("99_urban", "99_rural"),
                        n_psu   = c(2L, 2L)),
    calibration = list(margins = "stratum")
  ), class = c("ecv_design_spec", "list"))
}

.tiny_sample_two_dept <- function() {
  make_stratum <- function(dept, stratum_id, weight, y) {
    data.frame(
      dept = dept, area = "all", stratum = stratum_id,
      ssu_id = rep(paste0(stratum_id, "_psu", 1:2), each = 3),
      weight_cal = weight,
      migrated_last_12m = y,
      stringsAsFactors = FALSE
    )
  }
  rbind(
    make_stratum(11L, "11_all", 2.0, c(1L, 0L, 0L, 1L, 0L, 0L)),
    make_stratum(25L, "25_all", 3.0, c(1L, 1L, 0L, 0L, 1L, 0L))
  )
}

.tiny_frame_two_dept <- function() {
  data.frame(
    dept = c(11L, 11L, 11L, 25L, 25L, 25L, 25L),
    area = "all",
    municipality = c(11001L, 11001L, 11002L, 25001L, 25001L, 25002L, 25002L),
    migrated_last_12m = c(1L, 0L, 1L, 0L, 0L, 1L, 1L)
  )
}

test_that("direct_two_levels adds a dept-only level (population-weighted across area)", {
  s <- .tiny_sample_two_dept()
  out <- direct_two_levels(s, outcome = "migrated_last_12m")
  expect_equal(sum(out$level == "dept"), 2L)
  expect_setequal(out$label[out$level == "dept"], c("11", "25"))

  for (d in c(11L, 25L)) {
    sub <- s[s$dept == d, ]
    hajek <- sum(sub$weight_cal * sub$migrated_last_12m) / sum(sub$weight_cal)
    expect_equal(
      out$estimate[out$level == "dept" & out$label == as.character(d)],
      hajek, tolerance = 1e-9)
  }
})

test_that("compute_truth_two_levels adds a dept-only level", {
  fr <- .tiny_frame_two_dept()
  out <- compute_truth_two_levels(fr, in_sample_munis = unique(fr$municipality),
                                  outcome = "migrated_last_12m")
  expect_equal(sum(out$level == "dept"), 2L)
  expect_equal(
    out$theta_true[out$level == "dept" & out$label == "11"],
    mean(fr$migrated_last_12m[fr$dept == 11L]), tolerance = 1e-12)
  expect_equal(
    out$theta_true[out$level == "dept" & out$label == "25"],
    mean(fr$migrated_last_12m[fr$dept == 25L]), tolerance = 1e-12)
})

test_that("direct_two_levels_all / compute_truth_two_levels_all carry the dept level for ratio outcomes too", {
  fr   <- make_synthetic_frame()
  design <- .tiny_design_for_direct()
  s    <- draw_replicate(fr, design, seed = 100001)
  pop_totals <- as.data.frame(table(fr$stratum), stringsAsFactors = FALSE)
  names(pop_totals) <- c("stratum", "N")
  s <- apply_calibration(s, pop_totals)

  r  <- direct_two_levels_all(s)
  tr <- compute_truth_two_levels_all(fr, in_sample_munis = unique(fr$municipality))
  expect_true(all(r$outcome[r$level == "dept"] %in% sae_outcomes$outcome))
  expect_true(all(tr$outcome[tr$level == "dept"] %in% sae_outcomes$outcome))
  expect_true(sum(tr$level == "dept" & tr$outcome == "pers_migrated_12m") >= 1L)
})

test_that("direct_two_levels returns one national row plus dept_x_area rows", {
  fr   <- make_synthetic_frame()
  design <- tiny_design()
  s    <- draw_replicate(fr, design, seed = 100001)
  pop_totals <- as.data.frame(table(fr$stratum), stringsAsFactors = FALSE)
  names(pop_totals) <- c("stratum", "N")
  s <- apply_calibration(s, pop_totals)

  out <- direct_two_levels(s, outcome = "migrated_last_12m")

  expect_true(all(c("level", "label", "estimate", "se",
                    "n_obs", "df") %in% colnames(out)))
  expect_equal(sum(out$level == "national"), 1L)
  expect_true(sum(out$level == "dept_x_area") >= 1L)
  expect_true(all(out$estimate >= 0 & out$estimate <= 1))
  expect_true(all(out$se >= 0))
  hajek <- sum(s$weight_cal * s$migrated_last_12m, na.rm = TRUE) /
           sum(s$weight_cal[!is.na(s$migrated_last_12m)])
  expect_equal(out$estimate[out$level == "national"], hajek, tolerance = 1e-9)
})

test_that("compute_truth_two_levels reproduces planted cell proportions", {
  fr <- make_synthetic_frame(seed = 99)
  munis <- unique(fr$municipality)
  out <- compute_truth_two_levels(fr,
                                  in_sample_munis = munis,
                                  outcome = "migrated_last_12m")
  expect_true(all(c("level", "label", "theta_true") %in% colnames(out)))
  expect_equal(sum(out$level == "national"), 1L)

  nat <- mean(fr$migrated_last_12m)
  expect_equal(out$theta_true[out$level == "national"], nat, tolerance = 1e-12)

  expected_cells <- aggregate(migrated_last_12m ~ dept + area, data = fr,
                              FUN = mean)
  for (i in seq_len(nrow(expected_cells))) {
    lab <- paste(expected_cells$dept[i], expected_cells$area[i], sep = "_")
    expect_equal(out$theta_true[out$label == lab],
                 expected_cells$migrated_last_12m[i],
                 tolerance = 1e-12)
  }
})

test_that("compute_truth_two_levels restricts to in_sample_munis", {
  fr <- make_synthetic_frame()
  keep <- unique(fr$municipality)[1:4]
  out <- compute_truth_two_levels(fr, in_sample_munis = keep,
                                  outcome = "migrated_last_12m")
  sub_nat <- mean(fr$migrated_last_12m[fr$municipality %in% keep])
  expect_equal(out$theta_true[out$level == "national"], sub_nat,
               tolerance = 1e-12)
})

test_that("summarise_mc_metrics returns zero bias / full coverage on perfect data", {
  R <- 200L
  direct_long <- tibble::tibble(
    rep_id   = rep(seq_len(R), each = 1),
    level    = "national",
    label    = "national",
    estimate = rep(0.1, R),
    se       = rep(0.01, R),
    n_obs    = 1000L,
    df       = 50
  )
  truth <- tibble::tibble(level = "national", label = "national",
                          theta_true = 0.1)
  out <- summarise_mc_metrics(direct_long, truth)
  expect_equal(out$bias, 0)
  expect_equal(out$rmse, 0)
  expect_equal(out$coverage, 1)
  expect_equal(out$participation, 1)
})

test_that("summarise_mc_metrics recovers the correct sign of bias", {
  set.seed(101)
  R <- 500L
  truth_val <- 0.10
  bias_true <- 0.02
  direct_long <- tibble::tibble(
    rep_id   = seq_len(R),
    level    = "national",
    label    = "national",
    estimate = truth_val + bias_true + rnorm(R, sd = 0.01),
    se       = 0.01,
    n_obs    = 1000L,
    df       = 50
  )
  truth <- tibble::tibble(level = "national", label = "national",
                          theta_true = truth_val)
  out <- summarise_mc_metrics(direct_long, truth)
  expect_true(out$bias > 0)
  expect_lt(abs(out$bias - bias_true), 5 * out$mcse_bias)
})

test_that("summarise_mc_metrics drops NA rows and lowers participation", {
  R <- 100L
  direct_long <- tibble::tibble(
    rep_id   = seq_len(R),
    level    = "national",
    label    = "national",
    estimate = c(rep(NA_real_, 10), rep(0.1, R - 10)),
    se       = c(rep(NA_real_, 10), rep(0.01, R - 10)),
    n_obs    = 1000L,
    df       = 50
  )
  truth <- tibble::tibble(level = "national", label = "national",
                          theta_true = 0.1)
  out <- summarise_mc_metrics(direct_long, truth)
  expect_equal(out$participation, (R - 10) / R)
  expect_equal(out$bias, 0)
})

test_that("mc_running_trace returns one row per R in the grid", {
  set.seed(13)
  R <- 500L
  direct_long <- tibble::tibble(
    rep_id   = seq_len(R),
    level    = "national",
    label    = "national",
    estimate = 0.1 + rnorm(R, sd = 0.01),
    se       = 0.01
  )
  truth <- tibble::tibble(level = "national", label = "national",
                          theta_true = 0.1)
  grid <- c(50, 100, 200, 500)
  out <- mc_running_trace(direct_long, truth, R_grid = grid)
  expect_equal(nrow(out), length(grid))
  expect_equal(out$R, grid)
  expect_true(all(c("running_bias", "running_emp_se",
                    "running_coverage", "running_mcse") %in% colnames(out)))
  full <- summarise_mc_metrics(direct_long, truth)
  expect_equal(out$running_bias[out$R == 500], full$bias, tolerance = 1e-12)
})

test_that("cell_dist_summary collapses dept_x_area cells to quartile bands", {
  set.seed(7)
  K <- 66L
  cells <- tibble::tibble(
    level    = "dept_x_area",
    label    = paste(rep(11:43, each = 2), c("urban", "rural"), sep = "_")[1:K],
    bias     = rnorm(K, 0, 0.005),
    rmse     = runif(K, 0.005, 0.02),
    se_ratio = runif(K, 0.85, 1.15),
    coverage = pmin(pmax(0.93 + rnorm(K, 0, 0.02), 0), 1)
  )
  nat <- tibble::tibble(level = "national", label = "national",
                        bias = 0.001, rmse = 0.003, se_ratio = 1.02,
                        coverage = 0.95)
  out <- cell_dist_summary(rbind(nat, cells))
  expect_true(all(c("metric", "area_split", "median",
                    "q25", "q75", "min", "max") %in% colnames(out)))
  expect_setequal(out$metric, c("bias", "rmse", "se_ratio", "coverage"))
  expect_setequal(out$area_split, c("all", "urban", "rural"))
  expect_equal(
    out$median[out$metric == "rmse" & out$area_split == "all"],
    median(cells$rmse)
  )
})

test_that("direct_two_levels supports ratio outcomes", {
  fr   <- make_synthetic_frame()
  design <- .tiny_design_for_direct()
  s    <- draw_replicate(fr, design, seed = 100001)
  pop_totals <- as.data.frame(table(fr$stratum), stringsAsFactors = FALSE)
  names(pop_totals) <- c("stratum", "N")
  s <- apply_calibration(s, pop_totals)

  r <- direct_two_levels(s, numerator = "n_mig_total",
                         denominator = "n_pers_resp")
  nat <- r[r$level == "national", ]
  expect_equal(nat$estimate,
               sum(s$weight_cal * s$n_mig_total) /
                 sum(s$weight_cal * s$n_pers_resp),
               tolerance = 1e-9)
})

test_that("direct_two_levels_all stacks 6 outcomes; truth wrapper matches", {
  fr   <- make_synthetic_frame()
  design <- .tiny_design_for_direct()
  s    <- draw_replicate(fr, design, seed = 100001)
  pop_totals <- as.data.frame(table(fr$stratum), stringsAsFactors = FALSE)
  names(pop_totals) <- c("stratum", "N")
  s <- apply_calibration(s, pop_totals)

  r <- direct_two_levels_all(s)
  expect_setequal(unique(r$outcome), sae_outcomes$outcome)

  tr <- compute_truth_two_levels_all(fr,
          in_sample_munis = unique(fr$municipality))
  expect_setequal(unique(tr$outcome), sae_outcomes$outcome)
  expect_equal(
    tr$theta_true[tr$outcome == "pers_migrated_12m" &
                  tr$level == "national"],
    sum(fr$n_mig_total) / sum(fr$n_pers_resp))
})

test_that("summarise_mc_metrics groups by outcome when present", {
  dl <- tibble::tibble(
    rep_id = rep(1:3, 2),
    outcome = rep(c("a", "b"), each = 3),
    level = "national", label = "national",
    estimate = c(0.1, 0.11, 0.09, 0.5, 0.52, 0.48),
    se = 0.01)
  tr <- tibble::tibble(outcome = c("a", "b"), level = "national",
                       label = "national", theta_true = c(0.1, 0.5))
  out <- summarise_mc_metrics(dl, tr)
  expect_equal(nrow(out), 2L)
  expect_true("outcome" %in% colnames(out))
})
