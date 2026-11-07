# Build the full national CNPV 2018 household frame and write it under the
# data directory (SAECOL_DATA_DIR). Run once after all department ZIPs are
# unzipped under SAECOL_DATA_DIR/Census2018 (see data-raw/README.md).
suppressPackageStartupMessages({
  library(saecol2026)
  library(arrow)
})

frame <- build_national_frame()
out <- output_path("frame", create = TRUE)
write_parquet_safe(frame, file.path(out, "national_household_frame.parquet"))

# Privacy: only print aggregates.
cat("Households:    ", nrow(frame), "\n")
cat("Departments:   ", length(unique(frame$dept)), "\n")
cat("Municipalities:", nrow(unique(frame[, c("dept", "muni")])), "\n")
cat("Manzanas:      ", length(unique(frame$manzana)), "\n")
cat("Urban share:   ", round(mean(frame$area == "urban"), 4), "\n")
cat("Migrant share (national, household-level): ",
    round(mean(frame$migrated_last_12m, na.rm = TRUE), 4), "\n")
cat("NA migrant flag share: ",
    round(mean(is.na(frame$migrated_last_12m)), 6), "\n")
cat("HH internal-migrant share: ",
    round(mean(frame$migrated_internal_12m, na.rm = TRUE), 4), "\n")
cat("HH intl-migrant share:     ",
    round(mean(frame$migrated_intl_12m, na.rm = TRUE), 4), "\n")
cat("Person migration shares (total/internal/intl): ",
    paste(round(c(sum(frame$n_mig_total,    na.rm = TRUE),
                  sum(frame$n_mig_internal, na.rm = TRUE),
                  sum(frame$n_mig_intl,     na.rm = TRUE)) /
                sum(frame$n_pers_resp, na.rm = TRUE), 4), collapse = " / "),
    "\n")
