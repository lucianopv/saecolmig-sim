test_that("compute_psu_sizes counts dwellings per (stratum, municipality)", {
  fr <- make_synthetic_frame()
  sizes <- compute_psu_sizes(fr)
  expect_true(all(c("stratum", "psu_id", "n_dwellings") %in% colnames(sizes)))
  # Each municipality has 6 manzanas * 5 households = 30 dwellings
  expect_true(all(sizes$n_dwellings == 30L))
  # 8 municipalities total
  expect_equal(nrow(sizes), 8L)
})

test_that("compute_psu_sizes errors when the PSU column is missing", {
  fr <- make_synthetic_frame()
  fr$municipality <- NULL
  expect_error(compute_psu_sizes(fr), "PSU column 'municipality' not found")
})

test_that("draw_psus selects exactly n_psu PSUs from the stratum with PPS", {
  fr <- make_synthetic_frame()
  sizes <- compute_psu_sizes(fr)
  urban <- sizes[sizes$stratum == "99_urban", ]
  drawn <- draw_psus(urban, n_psu = 2, seed = 1)
  expect_equal(nrow(drawn), 2)
  expect_true(all(drawn$psu_id %in% urban$psu_id))
  expect_true(all(drawn$pi_psu > 0 & drawn$pi_psu <= 1))
})

test_that("draw_psus short-circuits forced-inclusion (n_psu == N) strata", {
  fr <- make_synthetic_frame()
  sizes <- compute_psu_sizes(fr)
  urban <- sizes[sizes$stratum == "99_urban", ]
  drawn <- draw_psus(urban, n_psu = nrow(urban), seed = 1)
  expect_equal(nrow(drawn), nrow(urban))
  expect_true(all(drawn$pi_psu == 1))
})

test_that("draw_psus short-circuits n_psu == 1 strata with PPS multinomial", {
  fr <- make_synthetic_frame()
  sizes <- compute_psu_sizes(fr)
  urban <- sizes[sizes$stratum == "99_urban", ]
  # All PSUs in the synthetic stratum have the same size (30), so the inclusion
  # prob for the single drawn PSU should equal 1 / nrow(urban).
  drawn <- draw_psus(urban, n_psu = 1, seed = 7)
  expect_equal(nrow(drawn), 1)
  expect_equal(drawn$pi_psu, 1 / nrow(urban))
})

test_that("build_segments chunks dwellings into target_size groups", {
  fr <- make_synthetic_frame()
  muni_id <- unique(fr$municipality)[1]
  one_muni <- fr[fr$municipality == muni_id, ]
  segmented <- build_segments(one_muni, target_size = 10L)
  expect_true("segment_id" %in% colnames(segmented))
  expect_equal(nrow(segmented), nrow(one_muni))
  # 30 dwellings / 10 -> 3 segments of exactly 10 each (no tail)
  expect_equal(length(unique(segmented$segment_id)), 3L)
  expect_true(all(table(segmented$segment_id) == 10L))
})

test_that("draw_segments returns dwellings from exactly n_seg distinct segments", {
  fr <- make_synthetic_frame()
  muni_id <- unique(fr$municipality)[1]
  segmented <- build_segments(fr[fr$municipality == muni_id, ], target_size = 10L)
  drawn <- draw_segments(segmented, n_seg = 2, seed = 7)
  expect_equal(length(unique(drawn$segment_id)), 2L)
  # pi_seg = n_seg / N
  expect_true(all(drawn$pi_seg == 2 / 3))
})

# ---------------------------------------------------------------------------
# Tiny ECV-style design for orchestrator tests
# ---------------------------------------------------------------------------
.tiny_design <- function() {
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

test_that("draw_replicate returns a sample with valid weights and stable seed", {
  fr   <- make_synthetic_frame()
  design <- .tiny_design()

  s1 <- draw_replicate(fr, design, seed = 100001)
  s2 <- draw_replicate(fr, design, seed = 100001)

  for (col in c("hh_id", "stratum", "municipality", "psu_id", "ssu_id",
                "segment_id", "pi_psu", "pi_seg", "weight",
                "migrated_last_12m", "rep_seed")) {
    expect_true(col %in% colnames(s1), info = paste("missing column:", col))
  }

  # Same seed -> same sample
  expect_equal(sort(s1$hh_id), sort(s2$hh_id))

  # 2 PSUs per stratum * 2 segments per PSU * 10 dwellings per segment = 40
  # (per stratum). Two strata -> 80 total.
  expect_equal(nrow(s1), 80L)

  expect_true(all(s1$weight > 0 & is.finite(s1$weight)))
  expect_equal(s1$weight, 1 / (s1$pi_psu * s1$pi_seg))
})

test_that("apply_calibration matches stratum population totals", {
  fr   <- make_synthetic_frame()
  design <- .tiny_design()
  s    <- draw_replicate(fr, design, seed = 100001)

  pop_totals <- as.data.frame(table(fr$stratum), stringsAsFactors = FALSE)
  names(pop_totals) <- c("stratum", "N")
  sc <- apply_calibration(s, pop_totals)

  for (st in pop_totals$stratum) {
    expect_equal(sum(sc$weight_cal[sc$stratum == st]),
                 pop_totals$N[pop_totals$stratum == st],
                 tolerance = 1e-6)
  }
})

test_that("build_all_segments returns one entry per municipality, each pre-segmented", {
  fr <- make_synthetic_frame()
  cache <- build_all_segments(fr, target_size = 10L)

  expect_type(cache, "list")
  expect_setequal(names(cache), as.character(unique(fr$municipality)))
  for (m in names(cache)[1:3]) {
    expect_true("segment_id" %in% colnames(cache[[m]]))
    expect_equal(nrow(cache[[m]]), sum(fr$municipality == as.integer(m)))
    expect_equal(length(unique(cache[[m]]$segment_id)),
                 ceiling(sum(fr$municipality == as.integer(m)) / 10))
  }
})

test_that("draw_replicate(segments = cache) is identical to without-cache for same seed", {
  fr <- make_synthetic_frame()
  design <- .tiny_design()

  s_nocache <- draw_replicate(fr, design, seed = 100001)
  cache     <- build_all_segments(fr, target_size = design$stage2$target_segment_size)
  s_cached  <- draw_replicate(fr, design, seed = 100001, segments = cache)

  expect_equal(sort(s_nocache$hh_id), sort(s_cached$hh_id))
  ord_no  <- order(s_nocache$hh_id)
  ord_yes <- order(s_cached$hh_id)
  expect_equal(s_nocache$pi_psu[ord_no], s_cached$pi_psu[ord_yes])
  expect_equal(s_nocache$pi_seg[ord_no], s_cached$pi_seg[ord_yes])
  expect_equal(s_nocache$weight[ord_no], s_cached$weight[ord_yes])
})

test_that("draw_psus handles a non-degenerate PPS-WOR draw without crashing", {
  set.seed(11L)
  sizes <- tibble::tibble(
    stratum = rep(1L, 5L),
    psu_id = paste0("P", 1:5),
    n_dwellings = c(100, 110, 90, 95, 105)
  )
  out <- draw_psus(sizes, n_psu = 2L, seed = 11L)
  expect_equal(nrow(out), 2L)
  expect_true("pi_psu" %in% colnames(out))
  expect_true(all(out$pi_psu > 0 & out$pi_psu < 1))
})

test_that("draw_psus survives a stratum with a certainty unit (pik == 1)", {
  # One PSU has 1e6 dwellings; the rest have 100. With n_psu = 2 the
  # first PSU's pik = 1.0 and the remaining slot is shared by 4 PSUs.
  # Pre-fix, this crashed `sampling::UPmaxentropy` with
  # "subscript out of bounds".
  set.seed(22L)
  sizes <- tibble::tibble(
    stratum = rep(2L, 5L),
    psu_id = paste0("Q", 1:5),
    n_dwellings = c(1e6, 100, 100, 100, 100)
  )
  out <- draw_psus(sizes, n_psu = 2L, seed = 22L)
  expect_equal(nrow(out), 2L)
  expect_true("pi_psu" %in% colnames(out))
  expect_true("Q1" %in% out$psu_id)  # certainty unit always drawn
  expect_equal(out$pi_psu[out$psu_id == "Q1"], 1)
  # The non-certainty pi values sum to (n_psu - n_certain) = 1.
  non_cert <- out$pi_psu[out$psu_id != "Q1"]
  expect_equal(length(non_cert), 1L)
  expect_true(non_cert > 0 & non_cert < 1)
})


test_that("end-to-end: average direct estimate is close to truth on synthetic frame", {
  fr   <- make_synthetic_frame(seed = 1)
  design <- structure(list(
    stage1 = list(unit = "municipality", selection = "pps_wor",
                  size_variable = "n_dwellings"),
    stage2 = list(unit = "segment", selection = "systematic_fractional",
                  target_segment_size = 10L, n_per_psu = 3L),
    stage3 = list(unit = "household", selection = "census"),
    strata = data.frame(stratum = c("99_urban", "99_rural"),
                        n_psu   = c(3L, 3L)),
    calibration = list(margins = "stratum")
  ), class = c("ecv_design_spec", "list"))

  truth <- mean(fr$migrated_last_12m, na.rm = TRUE)
  pop_totals <- as.data.frame(table(fr$stratum), stringsAsFactors = FALSE)
  names(pop_totals) <- c("stratum", "N")

  set.seed(2026)
  ests <- replicate(30, {
    s <- draw_replicate(fr, design, seed = sample.int(1e6, 1))
    s <- apply_calibration(s, pop_totals)
    sum(s$weight_cal * s$migrated_last_12m, na.rm = TRUE) /
      sum(s$weight_cal[!is.na(s$migrated_last_12m)])
  })

  # Coarse sanity check: mean of 30 replicate direct estimates within 7 pp
  # of truth on the small fixture (lots of MC noise expected).
  expect_lt(abs(mean(ests) - truth), 0.07)
})
