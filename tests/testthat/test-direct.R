# .tiny_design_for_direct() lives in helper-fixtures.R (shared with
# test-direct-validation.R).

test_that("direct_estimate returns a tidy per-area tibble", {
  fr   <- make_synthetic_frame()
  design <- .tiny_design_for_direct()
  s    <- draw_replicate(fr, design, seed = 100001)
  pop_totals <- as.data.frame(table(fr$stratum), stringsAsFactors = FALSE)
  names(pop_totals) <- c("stratum", "N")
  s <- apply_calibration(s, pop_totals)

  d <- direct_estimate(s, area_var = "municipality")

  expect_true(all(c("area_id", "theta_direct", "var_direct",
                    "n_eff", "n_dwell") %in% colnames(d)))
  expect_equal(nrow(d), length(unique(s$municipality)))
  # Hajek proportion = sum(w * y) / sum(w * !is.na(y))
  for (m in d$area_id[1:2]) {
    sm    <- s[s$municipality == m, ]
    wsum  <- sum(sm$weight_cal[!is.na(sm$migrated_last_12m)])
    psum  <- sum(sm$weight_cal * sm$migrated_last_12m, na.rm = TRUE)
    expect_equal(d$theta_direct[d$area_id == m], psum / wsum, tolerance = 1e-9)
  }
  expect_true(all(d$var_direct >= 0))
  for (m in d$area_id[1:2]) {
    expect_equal(d$n_dwell[d$area_id == m], sum(s$municipality == m))
  }
})

test_that("direct_estimate gives non-zero variance for single-PSU strata", {
  fr <- make_synthetic_frame()
  design <- .tiny_design_for_direct()
  design$strata$n_psu <- rep(1L, nrow(design$strata))
  design$strata$n_seg_per_psu <- rep(3L, nrow(design$strata))

  s <- draw_replicate(fr, design, seed = 100001)
  pop_totals <- as.data.frame(table(fr$stratum), stringsAsFactors = FALSE)
  names(pop_totals) <- c("stratum", "N")
  s <- apply_calibration(s, pop_totals)

  d <- direct_estimate(s, area_var = "municipality")

  # Taylor with survey.lonely.psu='average' gives non-zero variance from
  # multiple-segments-per-stratum even when stage 1 has only 1 PSU.
  expect_true(median(d$var_direct) > 1e-8,
              info = sprintf("median var_direct = %.3e", median(d$var_direct)))
})

test_that("smooth_variance produces strictly positive smoothed variances", {
  fr   <- make_synthetic_frame()
  design <- .tiny_design_for_direct()
  s    <- draw_replicate(fr, design, seed = 100001)
  pop_totals <- as.data.frame(table(fr$stratum), stringsAsFactors = FALSE)
  names(pop_totals) <- c("stratum", "N")
  s <- apply_calibration(s, pop_totals)

  d  <- direct_estimate(s)
  sm <- smooth_variance(d, method = "ghs")

  expect_true("var_direct_smoothed" %in% colnames(sm))
  expect_equal(nrow(sm), nrow(d))
  expect_true(all(sm$var_direct_smoothed > 0))
  expect_true(all(is.finite(sm$var_direct_smoothed)))
})

test_that("direct_estimate_ratio matches the weighted count ratio", {
  fr   <- make_synthetic_frame()
  design <- .tiny_design_for_direct()
  s    <- draw_replicate(fr, design, seed = 100001)
  pop_totals <- as.data.frame(table(fr$stratum), stringsAsFactors = FALSE)
  names(pop_totals) <- c("stratum", "N")
  s <- apply_calibration(s, pop_totals)

  d <- direct_estimate_ratio(s, area_var = "municipality",
                             numerator_var   = "n_mig_total",
                             denominator_var = "n_pers_resp")
  expect_true(all(c("area_id", "theta_direct", "var_direct",
                    "n_eff", "n_dwell") %in% colnames(d)))
  for (m in d$area_id[1:2]) {
    sm <- s[s$municipality == m, ]
    expect_equal(d$theta_direct[d$area_id == m],
                 sum(sm$weight_cal * sm$n_mig_total) /
                   sum(sm$weight_cal * sm$n_pers_resp),
                 tolerance = 1e-9)
  }
  expect_true(all(d$var_direct >= 0, na.rm = TRUE))
})

test_that("direct_estimate_all stacks the 6 registered outcomes", {
  fr   <- make_synthetic_frame()
  design <- .tiny_design_for_direct()
  s    <- draw_replicate(fr, design, seed = 100001)
  pop_totals <- as.data.frame(table(fr$stratum), stringsAsFactors = FALSE)
  names(pop_totals) <- c("stratum", "N")
  s <- apply_calibration(s, pop_totals)

  d <- direct_estimate_all(s, area_var = "municipality")
  expect_setequal(unique(d$outcome), sae_outcomes$outcome)
  n_areas <- length(unique(s$municipality))
  expect_equal(nrow(d), 6L * n_areas)
  # headline outcome identical to the single-outcome path
  d1 <- direct_estimate(s, area_var = "municipality")
  dh <- d[d$outcome == "migrated_last_12m", ]
  expect_equal(dh$theta_direct[order(dh$area_id)],
               d1$theta_direct[order(d1$area_id)], tolerance = 1e-12)
})

test_that("smooth_variance_all smooths within outcome", {
  fr   <- make_synthetic_frame()
  design <- .tiny_design_for_direct()
  s    <- draw_replicate(fr, design, seed = 100001)
  pop_totals <- as.data.frame(table(fr$stratum), stringsAsFactors = FALSE)
  names(pop_totals) <- c("stratum", "N")
  s  <- apply_calibration(s, pop_totals)
  d  <- direct_estimate_all(s)
  sm <- smooth_variance_all(d, method = "ghs")
  expect_equal(nrow(sm), nrow(d))
  expect_true(all(sm$var_direct_smoothed > 0, na.rm = TRUE))
  # per-outcome smoothing == smoothing one outcome's slice alone
  one  <- smooth_variance(d[d$outcome == "migrated_intl_12m", ],
                          method = "ghs")
  both <- sm[sm$outcome == "migrated_intl_12m", ]
  expect_equal(both$var_direct_smoothed[order(both$area_id)],
               one$var_direct_smoothed[order(one$area_id)],
               tolerance = 1e-12)
})

test_that("direct_estimate_ratio returns NA (with message) for zero-denominator areas", {
  fr   <- make_synthetic_frame()
  # force one municipality to have no responding persons
  zero_muni <- unique(fr$municipality)[1]
  fr$n_pers_resp[fr$municipality == zero_muni] <- 0L
  fr$n_mig_total[fr$municipality == zero_muni] <- 0L
  design <- .tiny_design_for_direct()
  s    <- draw_replicate(fr, design, seed = 100001)
  pop_totals <- as.data.frame(table(fr$stratum), stringsAsFactors = FALSE)
  names(pop_totals) <- c("stratum", "N")
  s <- apply_calibration(s, pop_totals)

  if (zero_muni %in% s$municipality) {
    expect_message(
      d <- direct_estimate_ratio(s, area_var = "municipality",
                                 numerator_var   = "n_mig_total",
                                 denominator_var = "n_pers_resp"),
      "zero/NA weighted denominator")
    expect_true(is.na(d$theta_direct[d$area_id == zero_muni]))
    expect_true(is.na(d$var_direct[d$area_id == zero_muni]))
    # the area row is retained with its dwelling count
    expect_equal(d$n_dwell[d$area_id == zero_muni],
                 sum(s$municipality == zero_muni))
    # and smoothing keeps it NA, all others positive
    d$outcome <- "x"
    sm <- smooth_variance_all(d)
    expect_true(is.na(sm$var_direct_smoothed[sm$area_id == zero_muni]))
    expect_true(all(sm$var_direct_smoothed[sm$area_id != zero_muni] > 0))
  } else {
    skip("zero-denominator municipality not drawn into this seed's sample")
  }
})

test_that("apply_calibration hits household AND person margins when given N_persons", {
  fr   <- make_synthetic_frame()
  design <- .tiny_design_for_direct()
  s    <- draw_replicate(fr, design, seed = 100001)

  dt <- data.table::as.data.table(fr)
  pt <- dt[, .(N = .N, N_persons = sum(n_persons)), by = stratum]
  pop_totals <- as.data.frame(pt[order(stratum)])

  s2 <- apply_calibration(s, pop_totals)
  chk <- data.table::as.data.table(s2)[, .(
    hh   = sum(weight_cal),
    pers = sum(weight_cal * n_persons)
  ), by = stratum]
  chk <- merge(chk, pt, by = "stratum")
  expect_equal(as.numeric(chk$hh),   as.numeric(chk$N),         tolerance = 1e-6)
  expect_equal(as.numeric(chk$pers), as.numeric(chk$N_persons), tolerance = 1e-6)
})

test_that("apply_calibration tolerates NA n_persons (NA contributes 0 to person margin)", {
  fr   <- make_synthetic_frame()
  design <- .tiny_design_for_direct()
  s    <- draw_replicate(fr, design, seed = 100001)
  # poison two households with NA size, as in the real census frame
  s$n_persons[c(3L, 17L)] <- NA_integer_

  dt <- data.table::as.data.table(fr)
  pt <- dt[, .(N = .N, N_persons = sum(n_persons)), by = stratum]
  pop_totals <- as.data.frame(pt[order(stratum)])

  s2 <- apply_calibration(s, pop_totals)
  expect_false(".n_persons_cal" %in% colnames(s2))
  expect_true(all(is.finite(s2$weight_cal)))
  chk <- data.table::as.data.table(s2)[, .(
    hh   = sum(weight_cal),
    pers = sum(weight_cal * ifelse(is.na(n_persons), 0, n_persons))
  ), by = stratum]
  chk <- merge(chk, pt, by = "stratum")
  expect_equal(as.numeric(chk$hh), as.numeric(chk$N), tolerance = 1e-6)
  # person margin: the frame totals here INCLUDE the persons of the households
  # given a missing size (only the sample was altered), so allow the small gap
  # those two households create rather than exact equality
  expect_true(all(abs(chk$pers - chk$N_persons) / chk$N_persons < 0.15))
})

test_that("apply_calibration without N_persons calibrates to the household counts only", {
  fr   <- make_synthetic_frame()
  design <- .tiny_design_for_direct()
  s    <- draw_replicate(fr, design, seed = 100001)
  pop_totals <- as.data.frame(table(fr$stratum), stringsAsFactors = FALSE)
  names(pop_totals) <- c("stratum", "N")
  s2 <- apply_calibration(s, pop_totals)
  chk <- data.table::as.data.table(s2)[, .(hh = sum(weight_cal)),
                                       by = stratum]
  chk <- merge(chk, pop_totals, by = "stratum")
  expect_equal(as.numeric(chk$hh), as.numeric(chk$N), tolerance = 1e-6)
})
