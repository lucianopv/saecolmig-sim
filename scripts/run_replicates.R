#!/usr/bin/env Rscript
# Run replicates of the simulation: for each replicate r, draw the sample (seed
# 100000 + r), compute and smooth the direct estimates, fit each method to each
# outcome with the fit's own seed, and write the direct estimates and the
# per-fit rows (layouts in DATA_DICTIONARY.md, section 4). A rerun skips what is
# complete, so the same command can be repeated after an interruption, and a
# range of replicates can be split across jobs of any scheduler.
#
# Usage:
#   Rscript scripts/run_replicates.R --reps=1-100 --out=DIR
#       [--methods=fh_arcsin,bayes_fh_beta,xgb_sub] [--outcomes=NAME,...]
#       [--workers=N] [--area-matrix=PATH] [--subarea-matrix=PATH]
#       [--extra-area=PATH] [--extra-subarea=PATH] [--xgb-params=PATH]
#       [--retry-failed]
#
#   --reps            replicate numbers: a range, a list or both ("1-100", "1,4,7-9")
#   --out             output directory
#   --methods         methods to fit (default: the three refitted methods)
#   --outcomes        outcomes to fit (default: all six)
#   --workers         replicates run at a time in forked processes (default 1)
#   --area-matrix     area-level covariate matrix used in place of the built one
#   --subarea-matrix  sub-area covariate matrix used in place of the built one
#   --extra-area      extra area-level covariates, joined on `municipality`
#   --extra-subarea   extra sub-area covariates, joined on `code`
#   --xgb-params      table of tuned XGB hyperparameters (data-raw/tune_xgb_sub.R),
#                     read in place of the shipped table
#   --retry-failed    fit again the outcomes whose fit failed in an earlier run
#
# Cores: each worker runs one fit at a time on one thread, except a Bayes-Beta
# fit, which runs its three chains in parallel; with bayes_fh_beta requested,
# --workers=N can occupy about 3N cores.
#
# The data root is the folder named by SAECOL_DATA_DIR. Bayes-Beta samples with
# CmdStan through cmdstanr; CMDSTAN names the installation and, when it is
# unset, is set to the newest ~/.cmdstan/cmdstan-<version>. The output directory
# keeps a run manifest (inputs, settings, package versions); a run whose
# manifest differs from it is refused.
#
# The run ends with one line: replicates run and skipped, and per method the
# number of replicate-outcome fits that failed (recorded with fit_failed = TRUE).

# One BLAS and OpenMP thread per process, set before any package loads: these
# libraries read the variables once, and several threads per fit oversubscribe
# the cores when replicates run in parallel.
Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1", MKL_NUM_THREADS = "1")

suppressPackageStartupMessages(library(saecol2026))
opt <- tryCatch(sim_parse_args(commandArgs(trailingOnly = TRUE)), error = function(e) {
  cat("run_replicates.R: ", conditionMessage(e), "\n", sep = "", file = stderr())
  quit(status = 2L)
})

# CmdStan location, before anything loads cmdstanr (the package itself does not).
needs_stan <- any(opt$methods %in% c("bayes_fh_beta", "bayes_fh_heavy"))
cmdstan <- sim_cmdstan_env(required = needs_stan)
if (needs_stan) cat("CmdStan:", cmdstan, "\n")

inputs <- sim_inputs(area_matrix = opt$area_matrix,
                     subarea_matrix = opt$subarea_matrix,
                     extra_area = opt$extra_area,
                     extra_subarea = opt$extra_subarea)
if (!is.null(opt$xgb_params)) inputs$tune <- xgb_sub_params(file = opt$xgb_params)

cat(sprintf("replicates %d-%d (%d) | methods %s | outcomes %d | workers %d | output %s\n",
            min(opt$reps), max(opt$reps), length(opt$reps), paste(opt$methods, collapse = ","),
            length(opt$outcomes), opt$workers, opt$out))
cat(sim_pkg_versions(), "\n")

st <- sim_run(opt$reps, inputs, methods = opt$methods, outcomes = opt$outcomes,
              out_dir = opt$out, workers = opt$workers, retry_failed = opt$retry_failed)
fc <- st[grep("^failed_", names(st))]
cat(sprintf("done: %d replicate(s), %d run and %d skipped as complete; failed fits: %s\n",
            nrow(st), sum(st$status == "done"), sum(st$status == "skipped"),
            paste0(sub("^failed_", "", names(fc)), " ", colSums(fc), collapse = ", ")))
