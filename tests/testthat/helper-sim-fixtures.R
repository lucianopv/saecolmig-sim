# Synthetic inputs of the replicate runner: 12 municipalities in two departments,
# a household frame, the stratum map, the census truth, an area-level covariate
# matrix and a sub-area matrix. Nothing here is drawn from real data.

# DIVIPOLA-style codes: department 5 and department 8, six municipalities each.
fixture_universe <- function() c(5001L:5006L, 8001L:8006L)

# Area-level covariate matrix: the key first, then covariates only.
fixture_area_matrix <- function(universe = fixture_universe(), seed = 61L) {
  withr::local_seed(seed)
  n <- length(universe)
  data.frame(municipality = universe,
             share_urban = round(stats::runif(n, 0.2, 0.9), 4),
             share_elderly = round(stats::runif(n, 0.05, 0.2), 4),
             ntl_mean_pw = round(stats::rlnorm(n, 0.5, 0.6), 4),
             stringsAsFactors = FALSE)
}

# Sub-area matrix of the fixture frame: the census block and rural-section
# populations of the frame (build_subarea_pop()), then two covariates. Each
# municipality has two urban blocks and one rural section; sample households match
# these codes, so a sample collapses to sub-areas.
fixture_subarea_matrix <- function(universe = fixture_universe(), seed = 62L) {
  withr::local_seed(seed)
  pop <- as.data.frame(build_subarea_pop(fixture_frame(universe)))
  pop <- pop[order(pop$code), ]
  rownames(pop) <- NULL
  pop$z1 <- round(stats::rnorm(nrow(pop)), 4)
  pop$z2 <- round(stats::runif(nrow(pop)), 4)
  pop
}

# Household frame: 20 households per municipality in four blocks of five (two
# urban, two rural). One household has no recorded size.
fixture_frame <- function(universe = fixture_universe(), seed = 63L) {
  withr::local_seed(seed)
  hh <- 20L
  n <- length(universe) * hh
  muni <- rep(universe, each = hh)
  blk <- rep(rep(1:4, each = 5L), times = length(universe))
  area <- ifelse(blk <= 2L, "urban", "rural")
  size <- sample(1L:6L, n, replace = TRUE)
  mig_int <- stats::rbinom(n, size, 0.08)
  mig_intl <- stats::rbinom(n, size, 0.02)
  fr <- data.frame(
    dept = as.integer(muni %/% 1000L),
    muni = as.integer(muni %% 1000L),
    municipality = as.integer(muni),
    area = area,
    seccion = sprintf("%05d%s%04d%04d", muni, ifelse(area == "urban", "1", "3"), 1L, blk),
    manzana = sprintf("%05d%s%04d%04d%04d%04d", muni, ifelse(area == "urban", "1", "3"),
                      1L, blk, 1L, 1L),
    hh_id = sprintf("H%06d", seq_len(n)),
    n_persons = size,
    migrated_last_12m = as.integer(mig_int + mig_intl > 0L),
    migrated_internal_12m = as.integer(mig_int > 0L),
    migrated_intl_12m = as.integer(mig_intl > 0L),
    n_pers_resp = size,
    n_mig_total = mig_int + mig_intl,
    n_mig_internal = mig_int,
    n_mig_intl = mig_intl,
    stringsAsFactors = FALSE)
  fr$n_persons[3L] <- NA_integer_
  fr
}

# Municipality -> stratum map, as the strata builder writes it. The stratum codes
# are three of the design's realised strata, four municipalities each.
fixture_strata_map <- function(universe = fixture_universe()) {
  data.frame(municipality = universe,
             dept = as.integer(universe %/% 1000L),
             stratum = rep(c(500L, 501L, 502L), each = length(universe) %/% 3L),
             source = "known",
             stringsAsFactors = FALSE)
}

# Census truth in the layout of the truth builder: one row per area and outcome.
fixture_truth_full <- function(universe = fixture_universe(), seed = 64L) {
  withr::local_seed(seed)
  out <- sae_outcomes$outcome
  d <- expand.grid(area_id = universe, outcome = out, stringsAsFactors = FALSE)
  d$theta_true <- round(stats::runif(nrow(d), 0.01, 0.2), 5)
  d$n_total <- sample(100L:900L, nrow(d), replace = TRUE)
  d[, c("area_id", "theta_true", "n_total", "outcome")]
}

# A data root in the layout the runner reads (see DATA_DICTIONARY.md); the files
# are written into a temporary directory that is removed when the caller exits.
fixture_data_dir <- function(universe = fixture_universe(),
                             strata_map = fixture_strata_map(universe),
                             .local_envir = parent.frame()) {
  root <- withr::local_tempdir(.local_envir = .local_envir)
  w <- function(d, ...) write_parquet_safe(d, file.path(root, ...))
  w(fixture_frame(universe), "frame", "national_household_frame.parquet")
  w(strata_map, "aux_tables", "muni_stratum_nn_extended.parquet")
  w(fixture_truth_full(universe), "aux_tables", "truth_full.parquet")
  w(fixture_area_matrix(universe), "aux_tables", "municipal_X_geo_bld.parquet")
  w(fixture_subarea_matrix(universe), "subarea", "X_geo_blk_bld.parquet")
  root
}

# Runner inputs on which an area-level model is identifiable: 24 municipalities in
# four design strata that select 8, 4, 3 and 1 municipalities, so a replicate
# draws 16 of them (three covariates and an intercept are fitted on 16 areas) and
# leaves 8 out of sample. Every municipality has two segments of ten households,
# fewer than any of these strata selects, so a drawn municipality is enumerated
# whole.
fixture_sim_universe <- function() c(5001L:5012L, 8001L:8012L)

fixture_sim_strata_map <- function(universe = fixture_sim_universe()) {
  data.frame(municipality = universe,
             dept = as.integer(universe %/% 1000L),
             stratum = c(rep(599L, 10), rep(524L, 6), rep(504L, 4), rep(500L, 4)),
             source = "known",
             stringsAsFactors = FALSE)
}

fixture_sim_inputs <- function() {
  u <- fixture_sim_universe()
  root <- fixture_data_dir(u, strata_map = fixture_sim_strata_map(u))
  suppressMessages(sim_inputs(root))
}

# Evaluation fixtures: the files of a replicate run, its truth and an external
# per-fit file, all synthetic. Every method and the external file see the same
# drawn municipalities in a replicate, and the direct file of a replicate holds
# exactly those municipalities.

# Municipalities of the evaluation fixtures: 60 in two departments, enough for a
# sample-size group of 25 in-sample municipalities.
fixture_eval_universe <- function() c(5001L:5030L, 8001L:8030L)

# Census truth of the evaluation universe (layout of `truth_full.parquet`).
fixture_truth <- function(universe = fixture_eval_universe()) {
  fixture_truth_full(universe, seed = 65L)
}

# The 30 municipalities drawn in replicate `rep`.
fixture_drawn <- function(rep, universe = fixture_eval_universe()) {
  withr::with_seed(7000L + rep, sort(sample(universe, 30L)))
}

# Dwellings drawn per municipality of replicate `rep`, the same for every outcome.
fixture_n_dwell <- function(rep, ids) {
  withr::with_seed(8000L + rep, sample(12L:240L, length(ids), replace = TRUE))
}

# Per-fit rows of one method and replicate, every outcome: one row per
# municipality of the universe, estimates scattered around the truth with a
# method-specific bias, intervals of +- 1.96 reported standard errors.
fixture_fit_rows <- function(method, config, rep, truth = fixture_truth(),
                             universe = fixture_eval_universe()) {
  k <- match(method, sae_methods_implemented)
  withr::local_seed(9000L + 100L * k + rep)
  drawn <- fixture_drawn(rep, universe)
  n <- length(universe)
  ins <- universe %in% drawn
  rows <- lapply(sae_outcomes$outcome, function(o) {
    tr <- truth$theta_true[truth$outcome == o][match(universe, truth$area_id[truth$outcome == o])]
    sd <- ifelse(ins, 0.02, 0.04)
    theta <- tr + 0.002 * k + stats::rnorm(n, 0, sd)
    mse <- (sd * stats::runif(n, 0.5, 1.5))^2
    data.frame(
      area_id = universe, theta_hat = theta, mse = mse,
      ci_lo = theta - 1.96 * sqrt(mse), ci_hi = theta + 1.96 * sqrt(mse),
      area_set = ifelse(ins, "in_sample", "out_of_sample"), fit_failed = FALSE,
      error_msg = NA_character_, config = config, method = method, replicate = as.integer(rep),
      outcome = o, seed = 1L, wall_secs = 0.1, weight = NA_character_,
      n_divergent = if (method == "bayes_fh_beta") 0 else NA_real_,
      max_rhat = if (method == "bayes_fh_beta") 1.002 else NA_real_,
      sigma2_v = NA_real_, povmap_sha = NA_character_, mse_type = NA_character_,
      pkg_versions = "fixture",
      stringsAsFactors = FALSE)
  })
  do.call(rbind, rows)[, sim_fit_columns]
}

# A run folder in the runner's layout: `direct/rep_NNNN.parquet` and
# `fits/<method>__<config>/rep_NNNN.parquet` for each replicate of `reps`.
# Written into a temporary directory removed when the caller exits.
fixture_run_dir <- function(reps, methods = sim_v3_methods, truth = fixture_truth(),
                            universe = fixture_eval_universe(),
                            .local_envir = parent.frame()) {
  root <- withr::local_tempdir(.local_envir = .local_envir)
  for (r in reps) {
    drawn <- fixture_drawn(r, universe)
    nd <- fixture_n_dwell(r, drawn)
    dr <- do.call(rbind, lapply(sae_outcomes$outcome, function(o) {
      data.frame(area_id = drawn, theta_direct = 0.1, var_direct = 1e-3, n_eff = nd / 2,
                 n_dwell = nd, outcome = o, var_direct_smoothed = 1e-3,
                 stringsAsFactors = FALSE)
    }))
    write_parquet_safe(dr, file.path(root, "direct", sprintf("rep_%04d.parquet", r)))
    for (m in methods) {
      cfg <- if (m %in% sae_methods_subarea) "geo_blk_bld" else "geo_bld"
      write_parquet_safe(fixture_fit_rows(m, cfg, r, truth, universe),
                         file.path(root, "fits", paste0(m, "__", cfg),
                                   sprintf("rep_%04d.parquet", r)))
    }
  }
  root
}

# An external per-fit file of the same replicates and drawn municipalities, with
# zero-padded text municipality codes.
fixture_external <- function(reps = 1:2, truth = fixture_truth(),
                             universe = fixture_eval_universe()) {
  withr::local_seed(9900L)
  rows <- lapply(reps, function(r) {
    drawn <- fixture_drawn(r, universe)
    do.call(rbind, lapply(sae_outcomes$outcome, function(o) {
      tr <- truth$theta_true[truth$outcome == o][match(universe, truth$area_id[truth$outcome == o])]
      ins <- universe %in% drawn
      v <- (ifelse(ins, 0.02, 0.04) * stats::runif(length(universe), 0.5, 1.5))^2
      est <- tr + stats::rnorm(length(universe), 0, sqrt(v))
      data.frame(replicate = as.integer(r), outcome = o, municipality = sprintf("%05d", universe),
                 sampled = ins, estimate = est, estimated_variance = v,
                 lower = est - 1.96 * sqrt(v), upper = est + 1.96 * sqrt(v), truth = tr,
                 stringsAsFactors = FALSE)
    }))
  })
  do.call(rbind, rows)
}
