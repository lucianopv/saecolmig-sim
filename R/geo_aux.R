# R/geo_aux.R
# Municipality-level geospatial covariates (2018) for the saecol2026 aux X
# matrix. Pure, area-weighted zonal-statistic helpers + a build_geo_aux()
# orchestrator. All public functions return tibbles keyed by `municipality`
# (integer DIVIPOLA = dept * 1000 + muni, matching build_dept_aux()).
#
# Output is always aggregated to municipality; no row-level / pixel-level data
# is ever returned or printed.

NULL

#' Area-weighted zonal statistic of a raster over municipality polygons
#'
#' @param shp sf of municipalities; must contain an integer `municipality` column.
#' @param rast a terra SpatRaster (single layer).
#' @param name character; output column name.
#' @param fun "mean" or "sum" (area-weighted via exactextractr coverage).
#' @return tibble(municipality, <name>), one row per polygon, input order.
#' @export
extract_raster_stat <- function(shp, rast, name, fun = c("mean", "sum")) {
  fun <- match.arg(fun)
  shp <- sf::st_transform(shp, terra::crs(rast))
  if (fun == "sum") {
    vals <- exactextractr::exact_extract(rast, shp, function(v, cov) {
      if (all(is.na(v))) NA_real_ else sum(v * cov, na.rm = TRUE)
    }, progress = FALSE)
  } else {
    vals <- exactextractr::exact_extract(rast, shp, fun, progress = FALSE)
  }
  out <- tibble::tibble(municipality = shp$municipality)
  out[[name]] <- as.numeric(vals)
  out
}

#' Population-weighted mean of a raster over municipality polygons
#'
#' Weights each value cell by its area-coverage fraction times the population
#' from `weights` (e.g. WorldPop) — i.e. the covariate value the average
#' resident experiences, rather than the average hectare. The weights raster may
#' be at a different resolution than `rast` (exactextractr aligns them); it is
#' reprojected to `rast`'s CRS if they differ.
#'
#' @param shp sf of municipalities with an integer `municipality` column.
#' @param rast a terra SpatRaster (single layer) of covariate values.
#' @param weights a terra SpatRaster of population weights.
#' @param name character; output column name.
#' @return tibble(municipality, <name>), one row per polygon, input order.
#' @export
extract_weighted_mean <- function(shp, rast, weights, name) {
  shp <- sf::st_transform(shp, terra::crs(rast))
  if (!terra::same.crs(rast, weights)) weights <- terra::project(weights, rast)
  # exactextractr weighted ops require the value and weight rasters to share a
  # grid (extent AND alignment). Real layers differ (global CHIRPS vs Colombia
  # WorldPop; offset origins) -> "Incompatible extents". Restrict the value to
  # the weights' extent, then resample the population onto the value grid.
  # method = "sum" aggregates the finer WorldPop counts into each (coarser) value
  # cell = the population that cell represents (terra's recommended pop method;
  # falls back to bilinear if "sum" is unavailable for the given grids).
  rast <- terra::crop(rast, terra::ext(weights))
  w <- tryCatch(terra::resample(weights, rast, method = "sum"),
                error = function(e) terra::resample(weights, rast, method = "bilinear"))
  # WorldPop no-data (ocean / outside Colombia) is NA; exact_extract's
  # weighted_mean returns NA for ANY polygon touching an NA-weight cell, which
  # NA's out every coastal/border muni (worse on finer grids). No-data = zero
  # population, so set it to 0 and weight by the populated cells only.
  w <- terra::subst(w, NA, 0)
  vals <- exactextractr::exact_extract(rast, shp, "weighted_mean",
                                       weights = w, progress = FALSE)
  out <- tibble::tibble(municipality = shp$municipality)
  out[[name]] <- as.numeric(vals)
  out
}

#' Mean Terrain Ruggedness Index per municipality, from an elevation raster
#'
#' @param shp sf of municipalities with `municipality`.
#' @param dem terra SpatRaster of elevation.
#' @return tibble(municipality, terrain_ruggedness)
#' @export
terrain_ruggedness <- function(shp, dem) {
  tri <- terra::terrain(dem, v = "TRI")
  extract_raster_stat(shp, tri, "terrain_ruggedness", fun = "mean")
}

#' 2018 precipitation-anomaly raster (2018 minus long-term normal)
#'
#' The normal is reprojected (CRS mismatch) or resampled (grid mismatch) onto the
#' 2018 grid before subtraction. Exposed so both the area-weighted and
#' population-weighted aggregations consume one shared anomaly raster.
#'
#' @param precip_2018 terra SpatRaster of 2018 annual precipitation.
#' @param precip_normal terra SpatRaster of the climatological annual normal.
#' @return terra SpatRaster of the anomaly (2018 - normal).
#' @export
precip_anomaly_raster <- function(precip_2018, precip_normal) {
  if (!terra::same.crs(precip_2018, precip_normal)) {
    precip_normal <- terra::project(precip_normal, precip_2018)
  } else if (!terra::compareGeom(precip_2018, precip_normal, stopOnError = FALSE)) {
    precip_normal <- terra::resample(precip_normal, precip_2018, method = "bilinear")
  }
  precip_2018 - precip_normal
}

#' Mean 2018 precipitation anomaly per municipality (2018 minus long-term normal)
#'
#' @param shp sf of municipalities with `municipality`.
#' @param precip_2018 terra SpatRaster of 2018 annual precipitation.
#' @param precip_normal terra SpatRaster of the climatological annual normal.
#' @return tibble(municipality, precip_anomaly_2018)
#' @export
precip_anomaly <- function(shp, precip_2018, precip_normal) {
  anom <- precip_anomaly_raster(precip_2018, precip_normal)
  extract_raster_stat(shp, anom, "precip_anomaly_2018", fun = "mean")
}

#' Distance (km) from each municipality centroid to a border line
#'
#' Distances are computed in a metric CRS (EPSG:3116, MAGNA-SIRGAS / Colombia
#' Bogota zone). Island municipalities (San Andres 88001, Providencia 88564)
#' return large finite distances; callers may set them NA downstream.
#'
#' @param shp sf of municipalities with `municipality`.
#' @param border sf of the Colombia-Venezuela boundary (LINESTRING/MULTILINESTRING).
#' @return tibble(municipality, dist_venezuela_border_km)
#' @export
dist_to_border <- function(shp, border) {
  shp_m    <- sf::st_transform(shp, 3116)
  border_m <- sf::st_transform(border, 3116)
  ctr <- suppressWarnings(sf::st_centroid(shp_m))
  d   <- sf::st_distance(ctr, sf::st_union(border_m))
  tibble::tibble(
    municipality = shp$municipality,
    dist_venezuela_border_km = as.numeric(d) / 1000
  )
}

#' Build the full municipality-level geospatial covariate table
#'
#' Each of the 8 continuous (mean-aggregated) covariates is produced in two
#' aggregations: area-weighted (e.g. `ntl_mean`) and, when `weighted = TRUE`,
#' population-weighted using `layers$pop_2018` as the weight (suffix `_pw`, e.g.
#' `ntl_mean_pw`) — the value the average resident experiences. `dist_to_border`
#' (geometry) and the pop-derived columns are not weighted.
#'
#' @param shp sf of municipalities with integer `municipality`.
#' @param layers named list of single-layer terra SpatRasters:
#'   ntl, travel_time, dem, builtup, cropland, forest, precip_2018,
#'   precip_normal, pop_2018, pop_2015.
#' @param border sf of the Colombia-Venezuela boundary.
#' @param weighted logical; if TRUE (default) also emit population-weighted `_pw`
#'   columns for the 8 continuous covariates.
#' @return tibble keyed by `municipality` (ascending), one row per muni. Includes
#'   the flagged `pop_change_2015_2018` column (caller decides whether to use it).
#' @export
build_geo_aux <- function(shp, layers, border, weighted = TRUE) {
  reduce_join <- function(acc, x) merge(acc, x, by = "municipality", all.x = TRUE)

  # Continuous covariates share one value raster between the area-weighted and
  # population-weighted aggregations.
  value_rasters <- list(
    ntl_mean            = layers$ntl,
    travel_time_city    = layers$travel_time,
    terrain_ruggedness  = terra::terrain(layers$dem, v = "TRI"),
    built_up_share      = layers$builtup,
    cropland_share      = layers$cropland,
    forest_share        = layers$forest,
    precip_annual_2018  = layers$precip_2018,
    precip_anomaly_2018 = precip_anomaly_raster(layers$precip_2018,
                                                layers$precip_normal)
  )

  parts <- list(dist_to_border(shp, border))
  for (nm in names(value_rasters)) {
    parts[[length(parts) + 1L]] <-
      extract_raster_stat(shp, value_rasters[[nm]], nm, "mean")
    if (weighted) {
      parts[[length(parts) + 1L]] <-
        extract_weighted_mean(shp, value_rasters[[nm]], layers$pop_2018,
                              paste0(nm, "_pw"))
    }
  }

  # pop density (persons / km^2) and the flagged relative pop change
  area_km2 <- as.numeric(sf::st_area(shp)) / 1e6
  pop18 <- extract_raster_stat(shp, layers$pop_2018, "pop_2018_sum", "sum")
  pop15 <- extract_raster_stat(shp, layers$pop_2015, "pop_2015_sum", "sum")
  pop <- tibble::tibble(
    municipality        = shp$municipality,
    pop_density_2018    = pop18$pop_2018_sum / area_km2,
    pop_change_2015_2018 = (pop18$pop_2018_sum - pop15$pop_2015_sum) /
                            pmax(pop15$pop_2015_sum, 1)
  )
  parts <- c(parts, list(pop))

  out <- Reduce(reduce_join, parts)
  out <- out[order(out$municipality), , drop = FALSE]
  stopifnot(!is.unsorted(out$municipality), !anyDuplicated(out$municipality),
            nrow(out) == nrow(shp))
  tibble::as_tibble(out)
}

#' Population-weighted mean of every layer of a raster, in one pass
#'
#' The multi-band twin of extract_weighted_mean(), with the same steps in the
#' same order (reproject weights if the CRS differs, crop values to the
#' weights' extent, resample weights onto the value grid with "sum", NA
#' weights -> 0, exactextractr weighted_mean), so each column equals a
#' single-layer call; test-buildings.R pins that. The one difference is
#' speed: when the weights already share the value grid, the resample (an
#' identity) is skipped, and all layers share one extraction.
#'
#' @param rast multi-layer SpatRaster; layer names become column names
#' @return tibble(municipality, <layer names>)
#' @export
extract_weighted_mean_multi <- function(shp, rast, weights) {
  shp <- sf::st_transform(shp, terra::crs(rast))
  if (!terra::same.crs(rast, weights)) weights <- terra::project(weights, rast)
  if (isTRUE(terra::compareGeom(rast, weights, stopOnError = FALSE))) {
    w <- weights
  } else {
    rast <- terra::crop(rast, terra::ext(weights))
    w <- tryCatch(terra::resample(weights, rast, method = "sum"),
                  error = function(e) terra::resample(weights, rast, method = "bilinear"))
  }
  w <- terra::subst(w, NA, 0)
  vals <- exactextractr::exact_extract(rast, shp, "weighted_mean",
                                       weights = w, progress = FALSE)
  vals <- as.data.frame(vals)
  if (terra::nlyr(rast) == 1L) names(vals) <- names(rast)
  else names(vals) <- sub("^weighted_mean\\.", "", names(vals))
  out <- tibble::tibble(municipality = shp$municipality)
  for (nm in names(rast)) out[[nm]] <- as.numeric(vals[[nm]])
  out
}
