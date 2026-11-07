#!/usr/bin/env Rscript
# data-raw/build_area_geo_bld.R
# The area-level matrix of configuration `geo_bld`: the geo matrix
# (assemble_geo_X.R) plus the buildings columns of the published selection,
# kept_cols("bld").
#
#   municipal  population-weighted building covariates per municipality, on two
#              WorldPop weightings
#              -> aux_tables/municipal_bld_<year>.parquet     (unconstrained 2018)
#              -> aux_tables/municipal_bld_<year>con.parquet  (constrained R2024B 2018)
#   assemble   -> aux_tables/municipal_X_geo_bld.parquet
#   vif        prints the VIF prune of the four building columns against the
#              geo matrix beside the published selection; writes nothing
#
#   Rscript data-raw/build_area_geo_bld.R <municipal|assemble|vif>
#
# The recipe is the other `_pw` layers': a WorldPop-weighted zonal mean of the
# four building bands (extract_weighted_mean_multi()), then
# bld_derive_covariates() for the presence-weighted height and SD. The matrix
# takes the unconstrained weights, as its geo block does; the constrained table
# is the municipal value for sub-areas without one
# (build_subarea_geo_blk_bld.R).
#
# Inputs under SAECOL_DATA_DIR: geo_cache/buildings_<year>_3as.tif
# (build_buildings_raster.R), geo_cache/worldpop_2018.tif,
# geo_cache/worldpop_r2024b_2018.tif, aux_tables/municipal_X_geo.parquet.
# Env: BLD_YEAR (2018), BLD_TAG (suffix of the raster and of every output, as in
#      build_buildings_raster.R).
#
# Prints counts, quantiles and VIFs only.

suppressPackageStartupMessages({
  library(saecol2026); library(sf); library(terra); library(arrow); library(data.table)
})
args <- commandArgs(trailingOnly = TRUE)
mode <- if (length(args)) args[[1]] else stop("usage: build_area_geo_bld.R <municipal|assemble|vif>")
YEAR <- as.integer(Sys.getenv("BLD_YEAR", "2018"))
TAG  <- Sys.getenv("BLD_TAG", "")
GEO  <- output_path("geo_cache")
AUX  <- output_path("aux_tables", create = TRUE)
say  <- function(...) { cat(format(Sys.time(), "%H:%M:%S"), " ", ..., "\n", sep = ""); flush.console() }
out_p    <- function(stem) file.path(AUX, paste0(stem, TAG, ".parquet"))
stem_unc <- sprintf("municipal_bld_%d", YEAR)
stem_con <- sprintf("municipal_bld_%dcon", YEAR)
kept_bld <- kept_cols("bld")

describe <- function(d, label) {
  q <- rbindlist(lapply(bld_cov_names, function(k) data.table(
    covariate = k, n_na = sum(is.na(d[[k]])),
    p05 = stats::quantile(d[[k]], .05, na.rm = TRUE),
    median = stats::median(d[[k]], na.rm = TRUE),
    p95 = stats::quantile(d[[k]], .95, na.rm = TRUE))))
  say(label); print(as.data.frame(q), digits = 4)
  stopifnot(all(d$building_presence_pw >= 0 & d$building_presence_pw <= 1, na.rm = TRUE),
            all(d$building_density_pw >= 0, na.rm = TRUE),
            all(d$building_height_pw >= 0, na.rm = TRUE))
}
parent <- function() {
  p <- file.path(AUX, "municipal_X_geo.parquet")
  if (!file.exists(p)) stop("run data-raw/assemble_geo_X.R first: ", p)
  as.data.table(read_parquet(p))
}

if (mode == "municipal") {
  # Planar geometry for the municipal polygons (the zonal statistics run on the
  # raster grid). Only this mode reads polygons.
  sf::sf_use_s2(FALSE)
  bp <- file.path(GEO, sprintf("buildings_%d_3as%s.tif", YEAR, TAG))
  if (!file.exists(bp)) stop("missing building raster: ", bp)
  b <- terra::rast(bp); names(b) <- bld_band_names
  shp <- sf::st_as_sf(colmaps::municipios)
  shp$municipality <- as.integer(as.character(shp$id))
  shp <- shp[order(shp$municipality), c("municipality", "geometry")]
  for (v in list(list(stem = stem_unc, pop = "worldpop_2018.tif"),
                 list(stem = stem_con, pop = "worldpop_r2024b_2018.tif"))) {
    t0 <- Sys.time()
    d  <- as.data.table(bld_pw_covariates(shp, b, terra::rast(file.path(GEO, v$pop))))
    # A municipality left without a value (e.g. islands outside the raster)
    # takes the median over municipalities, as in the area-level geo block.
    if (!nzchar(TAG)) for (k in bld_cov_names) {
      n_na <- sum(is.na(d[[k]]))
      if (n_na) {
        set(d, which(is.na(d[[k]])), k, stats::median(d[[k]], na.rm = TRUE))
        say(sprintf("median-imputed %d NA in %s", n_na, k))
      }
    }
    describe(d, sprintf("%s (%d municipalities, %.1f min)", v$stem, nrow(d),
                        as.numeric(difftime(Sys.time(), t0, units = "mins"))))
    stopifnot(nrow(d) == 1122L, !anyDuplicated(d$municipality))
    write_parquet_safe(d, out_p(v$stem))
    say("wrote ", out_p(v$stem))
  }

} else if (mode == "assemble") {
  P <- parent()
  b <- as.data.table(read_parquet(out_p(stem_unc)))
  stopifnot(all(kept_bld %in% bld_cov_names), all(kept_bld %in% names(b)))
  X <- merge(P, b[, c("municipality", kept_bld), with = FALSE], by = "municipality",
             all.x = TRUE, sort = FALSE)
  setcolorder(X, c(names(P), kept_bld))
  # The parent's columns and rows are carried unchanged, and no added column
  # may be missing.
  stopifnot(identical(X$municipality, P$municipality), nrow(X) == nrow(P))
  for (cc in setdiff(names(P), "municipality")) stopifnot(identical(X[[cc]], P[[cc]]))
  na <- vapply(X[, ..kept_bld], function(v) sum(is.na(v)), 0L)
  say(sprintf("geo_bld: rows %s | buildings: %s | NA in buildings: %d",
              format(nrow(X), big.mark = ","), paste(kept_bld, collapse = ", "), sum(na)))
  stopifnot(sum(na) == 0L)
  f <- file.path(AUX, "municipal_X_geo_bld.parquet")
  write_parquet_safe(tibble::as_tibble(X), f)
  say("wrote ", f)

} else if (mode == "vif") {
  P  <- parent()
  bm <- as.data.table(read_parquet(out_p(stem_unc)))
  base <- setdiff(names(P)[vapply(P, is.numeric, logical(1))],
                  c("municipality", "dept", "muni", "n_persons", "n_households"))
  stopifnot(!any(bld_cov_names %in% base))
  D <- merge(P, bm, by = "municipality")
  stopifnot(nrow(D) == nrow(P))
  r <- bld_prune_vif(as.data.frame(D), base = base, prunable = bld_cov_names)
  say(sprintf("prune against %d covariates keeps: %s | dropped: %s", length(base),
              paste(r$kept, collapse = ", "),
              if (nrow(r$dropped)) paste(sprintf("%s (%.1f)", r$dropped$column,
                                                 r$dropped$vif), collapse = ", ")
              else "none"))
  print(round(r$final_vif, 2))
  say("same as the published selection: ", identical(r$kept, kept_bld))

} else stop("unknown mode: ", mode)
