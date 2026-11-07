# data-raw/build_nn_strata.R
# Build the NN-extended 1,122-municipality -> ECV 2024 stratum map
# (muni_stratum_nn_extended.parquet). DANE's public ECV 2024 design only
# assigns 636 of the 1,122 CNPV 2018 municipalities to a stratum; the
# remaining ("out-of-sample", OOS) municipalities are assigned the stratum
# of their nearest in-frame neighbour (by department + a 12-feature census
# profile) via extend_strata_nn(). The sampling design (draw_psus()) reads
# this map, so every municipality in the frame has a stratum, not just the 636
# DANE names explicitly.
#
# Reads: SAECOL_DATA_DIR/aux_tables/municipal_X.parquet (12 census features for
#   all 1,122 munis, built by data-raw/build_auxiliary.R) and the DANE-known
#   636-muni map shipped in inst/extdata/ecv2024_muni_to_stratum.csv.
# Writes: SAECOL_DATA_DIR/aux_tables/muni_stratum_nn_extended.parquet (1,122
#   rows).

suppressPackageStartupMessages({
  library(saecol2026)
  library(arrow)
  library(data.table)
})

# 1. Load aux features and known map -------------------------------------
aux <- as.data.table(arrow::read_parquet(
  file.path(output_path("aux_tables"), "municipal_X.parquet")))
stopifnot(nrow(aux) == 1122L, "dept" %in% colnames(aux))
aux[, log_n_persons := log(pmax(n_persons, 1))]

known <- as.data.table(read.csv(
  system.file("extdata", "ecv2024_muni_to_stratum.csv",
              package = "saecol2026")))
stopifnot(nrow(known) == 636L)

# 2. Build NN-extended 1,122-muni -> stratum map --------------------------
feature_cols <- c(
  "share_urban", "share_working_age", "share_elderly",
  "share_secondary_completed", "share_illiterate",
  "share_economically_active", "share_afro", "share_indigenous",
  "share_piped_water", "share_electricity",
  "n_households", "log_n_persons"
)
ext_map <- extend_strata_nn(
  as.data.frame(aux),
  as.data.frame(known),
  id_col = "municipality",
  dept_col = "dept",
  stratum_col = "stratum",
  feature_cols = feature_cols)
stopifnot(nrow(ext_map) == 1122L, !anyNA(ext_map$stratum))

ext_path <- file.path(output_path("aux_tables", create = TRUE),
                      "muni_stratum_nn_extended.parquet")
arrow::write_parquet(ext_map, ext_path)
cat(sprintf("Wrote %s (1122 rows)\n", ext_path))

# Aggregate diagnostics only (no row-level data).
cat("Source breakdown (known DANE stratum vs NN-imputed):\n")
print(table(ext_map$source))
