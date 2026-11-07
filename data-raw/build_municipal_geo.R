# data-raw/build_municipal_geo.R
# Run build_geo_aux() on the cached rasters -> municipal_geo.parquet (1,122 rows).
# Aggregate output only.
#
# Raw sources: public geo rasters (nighttime lights, travel time, elevation,
# landcover, CHIRPS precipitation, WorldPop) pre-fetched into
# SAECOL_DATA_DIR/geo_cache (see data-raw/README.md) plus the GADM country
# boundaries, downloaded on the fly by geodata::gadm() into the same cache
# directory.

suppressPackageStartupMessages({
  library(saecol2026); library(terra); library(sf); library(arrow)
})

geo_dir <- saecol2026::output_path("geo_cache")
rd <- function(f) terra::rast(file.path(geo_dir, f))
# CHIRPS encodes missing data as -9999 and terra may not pick it up from the
# GeoTIFF NoData tag; precipitation is non-negative, so mask negatives -> NA on
# the RAW rasters (the anomaly CAN be legitimately negative, so it must be built
# from clean inputs). Without this, NoData cells drag muni means negative.
rd_precip <- function(f) { r <- rd(f); r[r < 0] <- NA; r }

# Muni shapefile (same source + key as build_queens_W).
shp <- sf::st_as_sf(colmaps::municipios)
shp$municipality <- as.integer(as.character(shp$id))
shp <- shp[order(shp$municipality), c("municipality", "geometry")]

layers <- list(
  ntl           = rd("ntl_2018.tif"),
  travel_time   = rd("travel_time_2015.tif"),
  dem           = rd("elevation.tif"),
  builtup       = rd("landcover_built.tif"),
  cropland      = rd("landcover_cropland.tif"),
  forest        = rd("landcover_trees.tif"),
  precip_2018   = rd_precip("chirps_2018.tif"),
  precip_normal = rd_precip("precip_normal.tif"),
  pop_2018      = rd("worldpop_2018.tif"),
  pop_2015      = rd("worldpop_2015.tif")
)

# Colombia-Venezuela border: shared boundary of the two GADM-0 polygons.
co <- sf::st_as_sf(geodata::gadm("COL", level = 0, path = geo_dir))
ve <- sf::st_as_sf(geodata::gadm("VEN", level = 0, path = geo_dir))
border <- sf::st_intersection(sf::st_boundary(co), sf::st_boundary(ve))

geo <- build_geo_aux(shp, layers, border)

# Islands (San Andres 88001, Providencia 88564) keep their valid (large) centroid
# distance to the mainland border. NAs arise only where a source raster is NoData
# over a muni (e.g. accessibility ocean cells); they are counted below and
# median-imputed in the assembly step (data-raw/assemble_geo_X.R).

out_path <- file.path(saecol2026::output_path("aux_tables"), "municipal_geo.parquet")
arrow::write_parquet(geo, out_path)

# Aggregate diagnostics (no row-level data).
cat(sprintf("municipal_geo.parquet: %d rows x %d cols\n", nrow(geo), ncol(geo)))
cat("NA counts per column:\n"); print(colSums(is.na(geo)))
cat("\nSummary (curated headline columns):\n")
print(summary(geo[, setdiff(names(geo),
      c("municipality", "pop_change_2015_2018"))]))
