#!/usr/bin/env Rscript
# data-raw/assemble_geo_X.R
# The area-level geo matrix: the census covariates of municipal_X.parquet plus
# the population-weighted geo covariates of the published selection,
# kept_cols("geo"). It is the parent of the geo_bld matrix
# (build_area_geo_bld.R).
#
#   inputs  aux_tables/municipal_X.parquet    (build_auxiliary.R)
#           aux_tables/municipal_geo.parquet  (build_municipal_geo.R)
#   output  aux_tables/municipal_X_geo.parquet
#
# The selection is a VIF prune (threshold 10, bld_prune_vif()) of ten candidate
# layers against the other covariates of the matrix. The prune is re-run here
# against the census covariates and printed beside the published selection,
# which is the one written: a prune against other covariates can keep other
# columns.
#
# Prints counts and VIFs only.

suppressPackageStartupMessages({ library(saecol2026); library(arrow) })

aux_dir <- output_path("aux_tables")
X   <- as.data.frame(read_parquet(file.path(aux_dir, "municipal_X.parquet")))
geo <- as.data.frame(read_parquet(file.path(aux_dir, "municipal_geo.parquet")))

candidates <- c("ntl_mean_pw", "travel_time_city_pw", "terrain_ruggedness_pw",
                "dist_venezuela_border_km", "built_up_share_pw", "pop_density_2018",
                "cropland_share_pw", "forest_share_pw", "precip_annual_2018_pw",
                "precip_anomaly_2018_pw")
kept_geo <- kept_cols("geo")
stopifnot(all(kept_geo %in% candidates), all(candidates %in% names(geo)))

Xg <- merge(X, geo[, c("municipality", candidates)], by = "municipality", all.x = TRUE)

# A municipality left without a geo value (a raster with no data over it, e.g.
# accessibility over the islands) takes the median over municipalities, so the
# matrix is complete: the area-level models drop rows with a missing covariate.
for (g in candidates) {
  if (anyNA(Xg[[g]])) {
    n_na <- sum(is.na(Xg[[g]]))
    Xg[[g]][is.na(Xg[[g]])] <- stats::median(Xg[[g]], na.rm = TRUE)
    cat(sprintf("median-imputed %d NA in %s\n", n_na, g))
  }
}

# The prune against the census covariates, for comparison only.
num  <- names(Xg)[vapply(Xg, is.numeric, logical(1))]
base <- setdiff(num, c("municipality", "dept", "muni", "n_persons", "n_households",
                       candidates))
r <- bld_prune_vif(Xg, base = base, prunable = candidates)
cat(sprintf("\nVIF prune against %d census covariates keeps %d of %d candidates\n",
            length(base), length(r$kept), length(candidates)))
if (nrow(r$dropped)) {
  cat("  dropped:", paste(sprintf("%s (%.1f)", r$dropped$column, r$dropped$vif),
                          collapse = ", "), "\n")
}
cat("  same as the published selection:", identical(r$kept, kept_geo), "\n")
if (!identical(r$kept, kept_geo)) {
  cat("  only in the prune     :", paste(setdiff(r$kept, kept_geo), collapse = ", "), "\n")
  cat("  only in the selection :", paste(setdiff(kept_geo, r$kept), collapse = ", "), "\n")
}

# Census covariates plus the published geo selection, in its order.
Xg <- Xg[, c(setdiff(names(Xg), candidates), kept_geo)]
stopifnot(nrow(Xg) == nrow(X), !anyNA(Xg[, kept_geo]))
write_parquet_safe(Xg, file.path(aux_dir, "municipal_X_geo.parquet"))
cat(sprintf("\nmunicipal_X_geo.parquet written: %d rows x %d cols (geo: %d)\n",
            nrow(Xg), ncol(Xg), length(kept_geo)))
