#!/usr/bin/env Rscript
# Evaluate the per-fit results of a replicate run against the census truth: read
# the direct and per-fit files of the run folder (DATA_DICTIONARY.md, sections
# 4.1-4.2) and, optionally, an external per-fit file of a method fitted
# elsewhere (section 4.3); score every estimate in sample and out of sample over
# the replicates; write the aggregate tables of section 4.4 as
# <out>/<table>.parquet.
#
# Usage:
#   Rscript scripts/evaluate.R --run=DIR --truth=PATH --reps=1-100 --out=DIR
#       [--external=PATH] [--external-label=MEGB] [--methods=NAME,...] [--column=geo]
#
#   --run             output folder of scripts/run_replicates.R
#   --truth           census truth (aux_tables/truth_full.parquet under the data root)
#   --reps            replicates to evaluate: a range, a list or both ("1-100", "1,4,7-9")
#   --out             folder for the tables
#   --external        external per-fit file (parquet or CSV)
#   --external-label  method label of the external rows (default MEGB)
#   --methods         methods of the run to evaluate (default: all in --run)
#   --column          label of the covariate set (default geo)
#
# The evaluation stops, saying why, when the external truth differs from the
# truth file, the external replicates differ from --reps, a method lacks a fit
# or a fit lacks municipalities, the methods saw different drawn municipalities,
# or a method keeps no fit of an outcome; it warns when the methods keep
# different numbers of fits of an outcome. The tables and the printed summary
# hold aggregates only.

suppressPackageStartupMessages(library(saecol2026))
opt <- tryCatch(sim_parse_eval_args(commandArgs(trailingOnly = TRUE)), error = function(e) {
  cat("evaluate.R: ", conditionMessage(e), "\n", sep = "", file = stderr())
  quit(status = 2L)
})

cat(sprintf("run %s | replicates %d-%d (%d) | external %s | output %s\n", opt$run,
            min(opt$reps), max(opt$reps), length(opt$reps),
            if (is.null(opt$external)) "none" else paste0(opt$external, " as ", opt$external_label),
            opt$out))
res <- sim_evaluate(opt$run, opt$truth, reps = opt$reps, methods = opt$methods,
                    external = opt$external, external_label = opt$external_label,
                    column = opt$column)
dir.create(opt$out, recursive = TRUE, showWarnings = FALSE)
for (k in names(res)) {
  write_parquet_safe(tibble::as_tibble(res[[k]]), file.path(opt$out, paste0(k, ".parquet")))
  cat(sprintf("wrote %s (%s rows)\n", k, format(nrow(res[[k]]), big.mark = ",")))
}

## ---- printed summary (aggregates only) ---------------------------------------
cat("\nImport check:\n")
print(as.data.frame(res$import_check), row.names = FALSE)
fc <- res$fit_check[kept < fits_read]
if (nrow(fc)) {
  cat("\nFits dropped, by method and outcome:\n")
  print(as.data.frame(fc), row.names = FALSE)
} else {
  cat("\nEvery fit read enters the tables.\n")
}
s <- res$set_summary
cat("\nMean over the outcomes, by method and set:\n")
print(as.data.frame(s[, .(coverage = mean(coverage), rmse_x1000 = 1e3 * mean(rmse),
                          ratio_pooled = mean(ratio_pooled), ratio_med = mean(ratio_med)),
                      by = .(method, area_set)][order(area_set, method)]),
      row.names = FALSE, digits = 3)
