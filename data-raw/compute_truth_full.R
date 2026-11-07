# Census "ground truth" per municipality and outcome over the FULL national
# frame: every municipality of the universe (1,122), not only the 636 that the
# ECV 2024 sampled (data-raw/compute_truth.R restricts to those).
#
# The runner scores every municipality of the universe: the drawn set changes
# from replicate to replicate, and the municipalities not drawn are scored out
# of sample. The truth therefore covers the whole universe, and sim_inputs()
# refuses a truth file that lacks a municipality or an outcome.
#
# Writes aux_tables/truth_full.parquet. Reads/writes under SAECOL_DATA_DIR only.
suppressPackageStartupMessages({
  library(saecol2026)
  library(arrow)
  library(data.table)
})

frame_path <- file.path(output_path("frame"), "national_household_frame.parquet")
if (!file.exists(frame_path)) {
  stop("Run data-raw/build_frame.R first to materialise the frame.")
}
frame <- as.data.table(read_parquet_safe(frame_path))

truth_full <- compute_truth_all(frame, area_var = "municipality")
out <- output_path("aux_tables", create = TRUE)
write_parquet_safe(truth_full, file.path(out, "truth_full.parquet"))

cat("Rows in truth_full (areas x outcomes):", nrow(truth_full), "\n")
cat("Municipalities covered:", length(unique(truth_full$area_id)), "\n")
for (oc in unique(truth_full$outcome)) {
  v <- truth_full$theta_true[truth_full$outcome == oc]
  cat(sprintf("%-28s mean %.4f  q5/q50/q95 %.4f/%.4f/%.4f  zero-areas %d\n",
              oc, mean(v, na.rm = TRUE),
              quantile(v, .05, na.rm = TRUE), quantile(v, .5, na.rm = TRUE),
              quantile(v, .95, na.rm = TRUE), sum(v == 0, na.rm = TRUE)))
}

# National rates: the population-level value over the whole frame (household
# means; person-level count ratios).
cat("\nNational rates over the full frame:\n")
for (i in seq_len(nrow(sae_outcomes))) {
  o <- as.list(sae_outcomes[i, ])
  v <- if (o$level == "household") {
    mean(frame[[o$indicator]], na.rm = TRUE)
  } else {
    sum(frame[[o$numerator]], na.rm = TRUE) / sum(frame[[o$denominator]], na.rm = TRUE)
  }
  cat(sprintf("%-28s %.5f\n", o$outcome, v))
}
