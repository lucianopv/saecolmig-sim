#!/usr/bin/env Rscript
# data-raw/tune_xgb_sub.R
# Per-outcome hyperparameters of the sub-area XGB method (fit_xgb_sub()).
#
#   input   subarea/X_geo_blk_bld.parquet (build_subarea_geo_blk_bld.R), or the
#           matrix named by XFILE; the national frame, the NN strata and the
#           ECV design, from which the tuning replicates are drawn
#   output  aux_tables/xgb_sub_params_<CONFIG>.tsv   (read with
#           xgb_sub_params(file = ...))
#           aux_tables/tune_xgb_sub_<CONFIG>.rds     (per-replicate winners)
#
# The shipped table inst/extdata/xgb_sub_params_geo_blk_bld.tsv holds the values
# the replicate runner passes to fit_xgb_sub() for configuration geo_blk_bld
# (sim_inputs()$tune); fit_xgb_sub() itself has generic defaults. This script
# repeats the selection on another matrix; pass its table to the runner with
# --xgb-params.
#
# Method. The cross-validated optimum moves with the outcome, so each outcome
# is tuned separately. For each outcome and each of the first TUNE_REPS
# replicates (samples drawn with seed 100000 + r, calibrated as in the
# simulation), povmap::xgb_tune() scores the whole grid by five-fold
# cross-validation on the design-weighted sub-area sample and reports its
# winning configuration. Each hyperparameter is then set to its most frequent
# value among the replicates' winners (the marginal mode), with a fixed 200
# rounds. Selecting on several replicates keeps one draw's noise from deciding
# the configuration. MEGB runs at its documented defaults (povmap has no MEGB
# tuner).
#
# Env: CONFIG (default geo_blk_bld; names the outputs), XFILE (path of another
#      sub-area matrix), TUNE_REPS (default 5).
# Runs 6 x TUNE_REPS grid searches of 32 configurations each. Prints selected
# values and cross-validated R2 only.

suppressPackageStartupMessages({
  library(saecol2026); library(arrow); library(data.table); library(povmap)
})
Sys.setenv(OMP_NUM_THREADS = 1, OPENBLAS_NUM_THREADS = 1)

CONFIG    <- Sys.getenv("CONFIG", "geo_blk_bld")
XFILE     <- Sys.getenv("XFILE", file.path(output_path("subarea"), "X_geo_blk_bld.parquet"))
TUNE_REPS <- as.integer(Sys.getenv("TUNE_REPS", "5"))
FOLDS     <- 5L
NROUND    <- 200L
GRID <- list(eta = c(0.02, 0.05, 0.10, 0.20),
             max_depth = c(2L, 3L, 4L, 6L),
             min_child_weight = c(1, 5))

## ---- sub-area matrix: the universe and the covariates the fit uses ----------
pop <- subarea_universe(as.data.table(read_parquet(XFILE)))
xc  <- subarea_covariate_names(pop)
fx  <- stats::as.formula(paste("y ~", paste(xc, collapse = " + ")))
cat("covariate matrix:", basename(XFILE), "|", format(nrow(pop), big.mark = ","),
    "sub-areas |", length(xc), "covariates\n")

## ---- the sampling inputs of the simulation -----------------------------------
# The frame with its stratum, the design, the segments and the calibration
# totals are built by the function the replicate runner uses.
samp <- sim_sampling_inputs()
draw <- function(r) {
  s <- draw_replicate(samp$frame, samp$design, seed = 100000L + r, segments = samp$segments)
  as.data.table(apply_calibration(s, samp$pop_totals))
}

## ---- one xgb_tune call per (outcome, replicate) over the whole grid ----------
t_all <- Sys.time()
smps  <- lapply(seq_len(TUNE_REPS), draw)
rows  <- list()
for (o in sae_outcomes$outcome) {
  for (i in seq_len(TUNE_REPS)) {
    sa <- collapse_sample_to_subarea(smps[[i]], pop, o)
    sm <- merge(sa[, .(code, y, wt)], pop, by = "code")
    sm[, municipality := as.character(as.integer(municipality))]
    t0 <- Sys.time()
    # Unspecified dials are fixed scalars, so the grid is the one stated above.
    r <- tryCatch(povmap::xgb_tune(
           fixed = fx, smp_data = as.data.frame(sm), smp_weights = "wt",
           domains = "municipality", transformation = "no", folds = FOLDS,
           eta = GRID$eta, max_depth = GRID$max_depth,
           min_child_weight = GRID$min_child_weight,
           nround = NROUND, subsample = 0.8, colsample_bytree = 0.8,
           lambda = 1, alpha = 0, gamma = 0,
           seed = 1000L + i, cpus = 1L, verbose = FALSE),
         error = function(e) e)
    el <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
    if (inherits(r, "error")) {
      cat(sprintf("  %-28s rep %d  ERROR %s\n", o, i, substr(conditionMessage(r), 1, 60)))
      next
    }
    rows[[length(rows) + 1L]] <- data.table(
      outcome = o, rep = i, eta = r$eta, max_depth = r$max_depth,
      min_child_weight = r$min_child_weight, nround = r$nround,
      mse_oos = r$mse_oos, r2_oos = r$r2_oos, secs = el)
    cat(sprintf("  %-28s rep %d  eta %.2f depth %d mcw %g | R2 %.4f | %.0f s\n",
                o, i, r$eta, r$max_depth, r$min_child_weight, r$r2_oos, el))
    utils::flush.console()
  }
}
R <- rbindlist(rows)

## ---- each hyperparameter on its modal value across replicates ----------------
marg <- function(col) {
  m <- R[, .N, by = c("outcome", col)][order(outcome, -N)][, .SD[1], by = outcome]
  setnames(m, c(col, "N"), c(col, paste0("n_", col)))[]
}
pick <- Reduce(function(a, b) merge(a, b, by = "outcome"),
               list(marg("eta"), marg("max_depth"), marg("min_child_weight")))
pick[, nround := NROUND]
pick <- merge(pick, R[, .(r2 = mean(r2_oos)), by = outcome], by = "outcome")
cat("\n--- selected configuration per outcome (agreement across replicates) ---\n")
print(pick[, .(outcome, eta, agree_eta = paste0(n_eta, "/", TUNE_REPS),
               max_depth, agree_depth = paste0(n_max_depth, "/", TUNE_REPS),
               min_child_weight, agree_mcw = paste0(n_min_child_weight, "/", TUNE_REPS),
               r2 = round(r2, 4))])

params <- stats::setNames(lapply(seq_len(nrow(pick)), function(i) list(
  nrounds = as.integer(pick$nround[i]), eta = pick$eta[i],
  max_depth = as.integer(pick$max_depth[i]),
  min_child_weight = pick$min_child_weight[i],
  subsample = 0.8, colsample_bytree = 0.8, lambda = 1, alpha = 0)), pick$outcome)
params <- params[sae_outcomes$outcome]

aux <- output_path("aux_tables", create = TRUE)
f_tsv <- write_xgb_sub_params(params, file.path(aux, sprintf("xgb_sub_params_%s.tsv", CONFIG)))
saveRDS(list(params = params, trace = R, selected = pick, grid = GRID, folds = FOLDS,
             tune_reps = TUNE_REPS,
             povmap = as.character(utils::packageVersion("povmap")),
             povmap_sha = utils::packageDescription("povmap")$RemoteSha,
             xgboost = as.character(utils::packageVersion("xgboost"))),
        file.path(aux, sprintf("tune_xgb_sub_%s.rds", CONFIG)))
cat("\nwrote:", f_tsv, "\n")
cat("total:", round(as.numeric(difftime(Sys.time(), t_all, units = "mins")), 1), "min\n")
