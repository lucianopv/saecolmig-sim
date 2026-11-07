# Compute census "ground truth" per municipality/outcome from the frame
# built by data-raw/build_frame.R, restricted to the ECV 2024 muni universe
# (bundled DANE strata map). Reads/writes under SAECOL_DATA_DIR only; no
# raw census access here (that happened when the frame was built).
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

muni_map_path <- system.file("extdata", "ecv2024_muni_to_stratum.csv",
                             package = "saecol2026")
if (!nzchar(muni_map_path)) {
  muni_map_path <- file.path("inst", "extdata", "ecv2024_muni_to_stratum.csv")
}
muni_map <- fread(file = muni_map_path)
frame <- frame[municipality %in% muni_map$municipality]

truth <- compute_truth_all(frame, area_var = "municipality")
out <- output_path("aux_tables", create = TRUE)
write_parquet_safe(truth, file.path(out, "truth.parquet"))

cat("Rows in truth (areas x outcomes):", nrow(truth), "\n")
for (oc in unique(truth$outcome)) {
  v <- truth$theta_true[truth$outcome == oc]
  cat(sprintf("%-28s mean %.4f  q5/q50/q95 %.4f/%.4f/%.4f  zero-areas %d\n",
              oc, mean(v, na.rm = TRUE),
              quantile(v, .05, na.rm = TRUE), quantile(v, .5, na.rm = TRUE),
              quantile(v, .95, na.rm = TRUE), sum(v == 0, na.rm = TRUE)))
}
