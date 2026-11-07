# Build the census-derived municipal covariate matrix: one row per
# municipality, aggregated from the CNPV 2018 department files under
# SAECOL_DATA_DIR/Census2018 (via census_path()) — see data-raw/README.md.
suppressPackageStartupMessages({
  library(saecol2026)
  library(arrow)
})

aux <- build_aux_municipal()

out_dir <- output_path("aux_tables", create = TRUE)
write_parquet_safe(aux, file.path(out_dir, "municipal_X.parquet"))

cat("Auxiliary rows (municipalities): ", nrow(aux), "\n")
cat("Columns: ", ncol(aux), "\n")
cat("Mean share_urban:    ", round(mean(aux$share_urban, na.rm = TRUE), 4), "\n")
cat("Mean share_indigenous:", round(mean(aux$share_indigenous, na.rm = TRUE), 4), "\n")
