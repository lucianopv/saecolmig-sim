test_that("zscore_features standardises numeric columns and preserves identifier", {
  df <- data.frame(
    municipality = c(1L, 2L, 3L, 4L),
    a = c(1, 2, 3, 4),
    b = c(10, 20, 30, 40)
  )
  out <- zscore_features(df, id_col = "municipality",
                         feature_cols = c("a", "b"))
  expect_equal(colnames(out), c("municipality", "a", "b"))
  expect_equal(out$municipality, df$municipality)
  expect_equal(mean(out$a), 0, tolerance = 1e-12)
  expect_equal(sd(out$a), 1, tolerance = 1e-12)
  expect_equal(mean(out$b), 0, tolerance = 1e-12)
  expect_equal(sd(out$b), 1, tolerance = 1e-12)
})

test_that("zscore_features errors when feature has zero variance", {
  df <- data.frame(municipality = 1:3, a = c(1, 1, 1))
  expect_error(
    zscore_features(df, id_col = "municipality", feature_cols = "a"),
    "zero variance"
  )
})

test_that("predict_loo_nn classifies a clean within-cell case", {
  # Two cells (cell 1 and cell 2), each with 3 munis. Within cell 1,
  # munis A and B share stratum 100, muni C is stratum 200.
  # Distances in the contrived feature x make A the NN of B and vice versa.
  features <- data.frame(
    municipality = c("A", "B", "C", "D", "E", "F"),
    cell         = c(1L, 1L, 1L, 2L, 2L, 2L),
    stratum      = c(100L, 100L, 200L, 300L, 300L, 400L),
    x            = c(0.0, 0.1, 5.0, 0.0, 0.1, 5.0)
  )
  res <- predict_loo_nn(features,
                        id_col = "municipality",
                        cell_col = "cell",
                        label_col = "stratum",
                        feature_cols = "x")
  expect_equal(nrow(res), 6)
  # A's NN within cell 1 (after hold-out) is B; B's stratum is 100 — match.
  expect_equal(res$predicted_stratum[res$municipality == "A"], 100L)
  # B's NN is A — match.
  expect_equal(res$predicted_stratum[res$municipality == "B"], 100L)
  # C's NN is A (distance 5 vs 4.9 to B); A's stratum is 100, C's true is
  # 200 — mismatch.
  expect_equal(res$predicted_stratum[res$municipality == "C"], 100L)
  expect_true(all(c("predicted_stratum", "candidate_count",
                    "distance_to_nn", "tie_flag") %in% colnames(res)))
})

test_that("predict_loo_nn returns NA when candidate set is empty", {
  # Muni A is the only member of cell 1; held-out, candidate set empty.
  features <- data.frame(
    municipality = c("A", "B", "C"),
    cell         = c(1L, 2L, 2L),
    stratum      = c(100L, 200L, 200L),
    x            = c(0.0, 1.0, 2.0)
  )
  res <- predict_loo_nn(features,
                        id_col = "municipality",
                        cell_col = "cell",
                        label_col = "stratum",
                        feature_cols = "x")
  a_row <- res[res$municipality == "A", ]
  expect_true(is.na(a_row$predicted_stratum))
  expect_equal(a_row$candidate_count, 0L)
  expect_true(is.na(a_row$distance_to_nn))
})

test_that("predict_loo_nn records ties", {
  # Equidistant candidates; which.min picks the first. Tie flag = TRUE.
  features <- data.frame(
    municipality = c("A", "B", "C"),
    cell         = c(1L, 1L, 1L),
    stratum      = c(100L, 200L, 300L),
    x            = c(0.0, 1.0, -1.0)
  )
  res <- predict_loo_nn(features,
                        id_col = "municipality",
                        cell_col = "cell",
                        label_col = "stratum",
                        feature_cols = "x")
  a_row <- res[res$municipality == "A", ]
  expect_true(a_row$tie_flag)
  expect_equal(a_row$candidate_count, 2L)
})

test_that("loo_nn_metrics computes exact-match, non-trivial match, and tau", {
  # 6 munis. Two cells. Cell 1 (4 munis): 3 strata. Cell 2 (2 munis):
  # both share stratum 400 (so all predictions in cell 2 are trivial-correct).
  res <- data.frame(
    municipality = c("A", "B", "C", "D", "E", "F"),
    cell = c(1L, 1L, 1L, 1L, 2L, 2L),
    true_stratum = c(100L, 100L, 200L, 300L, 400L, 400L),
    predicted_stratum = c(100L, 200L, 100L, 300L, 400L, 400L),
    candidate_count = c(3L, 3L, 3L, 3L, 1L, 1L),
    distance_to_nn = c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6),
    tie_flag = c(FALSE, FALSE, FALSE, FALSE, FALSE, FALSE),
    stringsAsFactors = FALSE
  )
  m <- loo_nn_metrics(res,
                      label_col = "stratum",
                      cell_col = "cell",
                      candidate_count_col = "candidate_count")
  # 4/6 exact matches (A, D, E, F).
  expect_equal(m$exact_match_rate, 4 / 6)
  # Non-trivial cells = cell 1 (3 strata). A & D match; B & C don't.
  expect_equal(m$non_trivial_exact_match_rate, 2 / 4)
  # Kendall's tau between predicted and true on predictable munis.
  expect_equal(m$kendall_tau,
               stats::cor(res$true_stratum, res$predicted_stratum,
                          method = "kendall"),
               tolerance = 1e-12)
  expect_equal(m$n_predictable, 6L)
  expect_equal(m$n_singleton_cell, 0L)
})

test_that("loo_nn_metrics excludes singleton-cell rows from headline metrics", {
  res <- data.frame(
    municipality = c("A", "B", "C"),
    cell = c(1L, 2L, 2L),
    true_stratum = c(100L, 200L, 200L),
    predicted_stratum = c(NA_integer_, 200L, 200L),
    candidate_count = c(0L, 1L, 1L),
    distance_to_nn = c(NA_real_, 0.1, 0.2),
    tie_flag = c(FALSE, FALSE, FALSE),
    stringsAsFactors = FALSE
  )
  # Both predictable rows carry the same label, so Kendall's tau is undefined
  # and cor() warns (its message is translated, so it is not matched).
  expect_warning(m <- loo_nn_metrics(res,
                                     label_col = "stratum",
                                     cell_col = "cell",
                                     candidate_count_col = "candidate_count"))
  expect_true(is.na(m$kendall_tau))
  expect_equal(m$exact_match_rate, 2 / 2)
  expect_equal(m$n_predictable, 2L)
  expect_equal(m$n_singleton_cell, 1L)
})

test_that("extend_strata_nn assigns OOS munis via within-dept NN", {
  # 4 munis. A, B known (stratum 100, 200). C, D OOS in same dept as A, B.
  # Feature x: A=0, B=10, C=0.1 (near A), D=9.9 (near B).
  aux <- data.frame(
    municipality = 1:4,
    dept = c(1L, 1L, 1L, 1L),
    n_persons = c(100, 100, 100, 100),
    n_households = c(20, 20, 20, 20),
    share_urban = c(0.5, 0.5, 0.5, 0.5),
    share_working_age = c(0.6, 0.6, 0.6, 0.6),
    share_elderly = c(0.1, 0.1, 0.1, 0.1),
    share_secondary_completed = c(0.4, 0.4, 0.4, 0.4),
    share_illiterate = c(0.1, 0.1, 0.1, 0.1),
    share_economically_active = c(0.6, 0.6, 0.6, 0.6),
    share_afro = c(0.05, 0.05, 0.05, 0.05),
    share_indigenous = c(0.05, 0.05, 0.05, 0.05),
    share_piped_water = c(0.9, 0.9, 0.9, 0.9),
    share_electricity = c(0.95, 0.95, 0.95, 0.95),
    discriminator = c(0.0, 10.0, 0.1, 9.9)
  )
  known <- data.frame(
    municipality = c(1L, 2L),
    stratum = c(100L, 200L)
  )
  out <- extend_strata_nn(aux, known,
                          id_col = "municipality",
                          dept_col = "dept",
                          stratum_col = "stratum",
                          feature_cols = "discriminator")
  expect_equal(nrow(out), 4L)
  expect_equal(sort(colnames(out)),
               sort(c("municipality", "dept", "stratum", "source")))
  # 1 and 2 keep their true strata.
  expect_equal(out$stratum[out$municipality == 1L], 100L)
  expect_equal(out$source[out$municipality == 1L], "known")
  # 3's nearest in dept 1 is A (stratum 100). 4's nearest is B (stratum 200).
  expect_equal(out$stratum[out$municipality == 3L], 100L)
  expect_equal(out$source[out$municipality == 3L], "nn")
  expect_equal(out$stratum[out$municipality == 4L], 200L)
})

test_that("extend_strata_nn errors when an OOS muni has no in-dept candidate", {
  # OOS muni in dept 2 has no in-frame (known) muni in dept 2.
  aux <- data.frame(
    municipality = 1:3,
    dept = c(1L, 1L, 2L),
    n_persons = c(100, 100, 100),
    n_households = c(20, 20, 20),
    share_urban = c(0.5, 0.5, 0.5),
    share_working_age = c(0.6, 0.6, 0.6),
    share_elderly = c(0.1, 0.1, 0.1),
    share_secondary_completed = c(0.4, 0.4, 0.4),
    share_illiterate = c(0.1, 0.1, 0.1),
    share_economically_active = c(0.6, 0.6, 0.6),
    share_afro = c(0.05, 0.05, 0.05),
    share_indigenous = c(0.05, 0.05, 0.05),
    share_piped_water = c(0.9, 0.9, 0.9),
    share_electricity = c(0.95, 0.95, 0.95),
    discriminator = c(0.0, 1.0, 5.0)
  )
  known <- data.frame(
    municipality = c(1L, 2L),
    stratum = c(100L, 200L)
  )
  expect_error(
    extend_strata_nn(aux, known,
                     id_col = "municipality",
                     dept_col = "dept",
                     stratum_col = "stratum",
                     feature_cols = "discriminator"),
    "no in-frame candidate"
  )
})

test_that("simulate_stage1_pps draws n_psu PSUs per stratum each rep", {
  # Stratum 100 has 4 munis with equal size; n_psu = 2. Stratum 200 has
  # 1 muni with size = 50; n_psu = 1 (take-all).
  universe <- data.frame(
    municipality = 1:5,
    stratum = c(100L, 100L, 100L, 100L, 200L),
    size = c(10, 10, 10, 10, 50)
  )
  alloc <- data.frame(
    stratum = c(100L, 200L),
    n_psu = c(2L, 1L)
  )
  set.seed(42L)
  sim <- simulate_stage1_pps(universe, alloc,
                             id_col = "municipality",
                             stratum_col = "stratum",
                             size_col = "size",
                             n_psu_col = "n_psu",
                             S = 20L,
                             base_seed = 100000L)
  expect_equal(nrow(sim), 5L)
  expect_true(all(c("municipality", "stratum", "emp_pi", "n_drawn",
                    "theoretical_pi", "n_reps") %in% colnames(sim)))
  expect_equal(sim$n_reps, rep(20L, 5L))
  # Muni 5 is take-all in its stratum -> emp_pi = 1.
  expect_equal(sim$emp_pi[sim$municipality == 5], 1)
  # Stratum 100: each muni's emp_pi ~ 0.5 (theoretical), but realised
  # over 20 reps may differ. Check theoretical_pi values.
  expect_equal(sim$theoretical_pi[sim$municipality == 1], 0.5)
  expect_equal(sim$theoretical_pi[sim$municipality == 5], 1)
})

test_that("simulate_stage1_pps respects per-stratum n_psu count exactly", {
  # In each rep, exactly n_psu PSUs per stratum should be drawn.
  universe <- data.frame(
    municipality = 1:6,
    stratum = c(rep(100L, 3), rep(200L, 3)),
    size = c(10, 10, 10, 5, 5, 5)
  )
  alloc <- data.frame(
    stratum = c(100L, 200L),
    n_psu = c(2L, 1L)
  )
  set.seed(7L)
  sim <- simulate_stage1_pps(universe, alloc,
                             id_col = "municipality",
                             stratum_col = "stratum",
                             size_col = "size",
                             n_psu_col = "n_psu",
                             S = 10L,
                             base_seed = 100000L)
  # Each rep draws exactly 2 PSUs in stratum 100 and 1 in stratum 200.
  # Aggregate: over S = 10 reps, sum of n_drawn within each stratum
  # = n_psu * S.
  draws_100 <- sum(sim$n_drawn[sim$stratum == 100L])
  draws_200 <- sum(sim$n_drawn[sim$stratum == 200L])
  expect_equal(draws_100, 2L * 10L)
  expect_equal(draws_200, 1L * 10L)
})
