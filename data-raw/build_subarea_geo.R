#!/usr/bin/env Rscript
# data-raw/build_subarea_geo.R
# Geospatial covariates of the sub-area matrix, on WorldPop constrained R2024B.
#
#   subarea    the geo block at census-block / rural-section resolution, plus the
#              polygon area of each sub-area
#              -> aux_tables/subarea_geo_parts/dept_<dd>.parquet (checkpoints)
#              -> aux_tables/subarea_geo.parquet
#   municipal  the same layers over municipal polygons: the value a sub-area
#              takes when its block-level extraction is missing
#              -> aux_tables/municipal_geo_2018con.parquet
#
#   Rscript data-raw/build_subarea_geo.R <subarea|municipal>
#
# Both passes run build_geo_aux(), the recipe of the area-level geo block
# (data-raw/build_municipal_geo.R), on the same rasters; only the polygons and
# the population raster differ from the area-level block, which weights by the
# unconstrained WorldPop. The sub-area matrix keeps the published geo selection,
# kept_cols("geo", level = "subarea"), under the layer names without the census
# year (geo_subarea_names()).
#
# Inputs under SAECOL_DATA_DIR: aux_tables/subarea_pop.parquet
# (build_subarea_pop.R); geo_cache/ rasters ntl_2018, travel_time_2015,
# elevation, landcover_{built,cropland,trees}, chirps_2018, precip_normal,
# worldpop_r2024b_{2018,2015} (see data-raw/README.md); the GADM boundaries,
# fetched into geo_cache by geodata::gadm().
# Env: MGN_DIR, the unzipped DANE MGN 2018 GeoPackages (default: census_path());
#      DEPTS, a comma-separated department subset (checkpoints only, no
#      national assembly).
#
# Prints counts and coverage only.

suppressPackageStartupMessages({
  library(saecol2026); library(sf); library(terra); library(arrow); library(data.table)
})
args <- commandArgs(trailingOnly = TRUE)
mode <- if (length(args)) args[[1]] else stop("usage: build_subarea_geo.R <subarea|municipal>")

GEO <- output_path("geo_cache")
AUX <- output_path("aux_tables", create = TRUE)
kept_sub <- kept_cols("geo", level = "subarea")

rd <- function(f) terra::rast(file.path(GEO, f))
# CHIRPS stores no-data as -9999 without a NoData tag that terra reads, and
# precipitation is non-negative, so negative cells are set to NA.
rd_precip <- function(f) { r <- rd(f); r[r < 0] <- NA; r }

# Layer list in build_geo_aux()'s names; the population weight and the
# population for density are the constrained R2024B rasters.
geo_layers <- function(mask_precip) {
  pr <- if (mask_precip) rd_precip else rd
  list(ntl           = rd("ntl_2018.tif"),
       travel_time   = rd("travel_time_2015.tif"),
       dem           = rd("elevation.tif"),
       builtup       = rd("landcover_built.tif"),
       cropland      = rd("landcover_cropland.tif"),
       forest        = rd("landcover_trees.tif"),
       precip_2018   = pr("chirps_2018.tif"),
       precip_normal = pr("precip_normal.tif"),
       pop_2018      = rd("worldpop_r2024b_2018.tif"),
       pop_2015      = rd("worldpop_r2024b_2015.tif"))
}

# Colombia-Venezuela border: shared boundary of the two GADM-0 polygons.
venezuela_border <- function() {
  co <- sf::st_as_sf(geodata::gadm("COL", level = 0, path = GEO))
  ve <- sf::st_as_sf(geodata::gadm("VEN", level = 0, path = GEO))
  sf::st_intersection(sf::st_boundary(co), sf::st_boundary(ve))
}

if (mode == "subarea") {
  # MGN polygons are not all valid on the sphere, so s2 is off: st_area() then
  # gives ellipsoidal areas (lwgeom), which differ from spherical ones by a
  # near-constant factor. Distances are computed in a projected CRS and the
  # zonal statistics on the raster grid, so neither depends on the setting.
  sf::sf_use_s2(FALSE)
  MGN_DIR <- Sys.getenv("MGN_DIR", census_path())
  gp <- function(f) {
    p <- file.path(MGN_DIR, f)
    if (!file.exists(p)) stop("missing MGN GeoPackage: ", p, "\nset MGN_DIR")
    p
  }
  CKPT <- file.path(AUX, "subarea_geo_parts")
  dir.create(CKPT, showWarnings = FALSE, recursive = TRUE)

  pop <- as.data.table(read_parquet(file.path(AUX, "subarea_pop.parquet")))
  pop[, cls := substr(code, 6, 6)]
  cat("sub-areas to place:", format(nrow(pop), big.mark = ","), "\n")

  depts <- sort(unique(substr(pop$code, 1, 2)))
  SUBSET <- nzchar(Sys.getenv("DEPTS"))
  if (SUBSET) depts <- intersect(depts, strsplit(Sys.getenv("DEPTS"), ",")[[1]])
  todo <- depts[!file.exists(file.path(CKPT, sprintf("dept_%s.parquet", depts)))]
  if (length(todo)) {
    F_MZ <- gp("MGN_URB_MANZANA2018geopackage.gpkg")
    F_RS <- gp("MGN_RUR_SECCION2018geopackage.gpkg")
    F_US <- gp("MGN_URB_SECCION2018geopackage.gpkg")
    layers <- geo_layers(mask_precip = TRUE)
    border <- venezuela_border()
  }

  t_all <- Sys.time()
  for (dd in todo) {
    t0 <- Sys.time()
    g <- subarea_geometry(dd, pop, F_MZ, F_RS, F_US)
    if (is.null(g)) { cat("  dept", dd, ": no polygons\n"); next }
    g <- subarea_union_codes(g)
    nd <- attr(g, "n_duplicate")
    # build_geo_aux() keys on `municipality`; the sub-area code is passed under
    # that name and renamed back.
    g$municipality <- g$code
    geo <- as.data.table(build_geo_aux(g[, c("municipality", attr(g, "sf_column"))],
                                       layers, border))
    setnames(geo, geo_subarea_names(names(geo)))
    setnames(geo, "municipality", "code")
    # Polygon area of each sub-area, from the polygons the densities were
    # divided by, so pop_density * area_km2 is the WorldPop zonal sum.
    a <- data.table(code = g$code, area_km2 = as.numeric(sf::st_area(g)) / 1e6)
    geo <- merge(geo, a[, .(area_km2 = sum(area_km2)), by = code], by = "code", all.x = TRUE)
    write_parquet_safe(tibble::as_tibble(geo), file.path(CKPT, sprintf("dept_%s.parquet", dd)))
    cat(sprintf("  dept %s: %6d sub-areas (%d duplicate-code rows merged), %.1f min\n",
                dd, nrow(geo), nd,
                as.numeric(difftime(Sys.time(), t0, units = "mins"))))
    utils::flush.console()
  }
  if (SUBSET) {
    cat("\nDEPTS subset: checkpoints written, national assembly skipped.\n")
    quit(save = "no")
  }

  geo <- rbindlist(lapply(list.files(CKPT, "\\.parquet$", full.names = TRUE),
                          function(f) as.data.table(read_parquet(f))), fill = TRUE)
  stopifnot(!anyDuplicated(geo$code))
  miss <- setdiff(kept_sub, names(geo))
  if (length(miss)) stop("geo block is missing selected columns: ", paste(miss, collapse = ", "))

  out <- merge(pop[, .(code, municipality, stratum, n_households_census, n_persons_census)],
               geo[, c("code", kept_sub, "area_km2"), with = FALSE], by = "code", all.x = TRUE)
  cc <- stats::complete.cases(out[, ..kept_sub])
  cat("\n--- coverage ---\n")
  cat("sub-areas:", format(nrow(out), big.mark = ","),
      "| complete on all", length(kept_sub), "geo covariates:",
      format(sum(cc), big.mark = ","), sprintf("(%.2f%%)\n", 100 * mean(cc)))
  cat("households covered:", sprintf("%.2f%%\n",
      100 * sum(out$n_households_census[cc]) / sum(out$n_households_census)))
  cat("sub-areas with a polygon area:", format(sum(!is.na(out$area_km2)), big.mark = ","), "\n")
  p_out <- file.path(AUX, "subarea_geo.parquet")
  write_parquet_safe(tibble::as_tibble(out), p_out)
  cat("\nwrote:", p_out, "\n")
  cat("elapsed:", round(as.numeric(difftime(Sys.time(), t_all, units = "mins")), 1), "min\n")

} else if (mode == "municipal") {
  # Municipal polygons with s2 on (spherical areas), as for the area-level geo
  # block.
  sf::sf_use_s2(TRUE)
  shp <- sf::st_as_sf(colmaps::municipios)
  shp$municipality <- as.integer(as.character(shp$id))
  shp <- shp[order(shp$municipality), c("municipality", "geometry")]
  # Precipitation is read as stored, CHIRPS' -9999 cells unmasked, which gives
  # the fallback values of the published sub-area matrix. The population
  # weights are zero on all but a few coastal no-data cells, so only the
  # precipitation columns of the municipalities holding those cells take them
  # in; the block-level pass masks them.
  geo <- as.data.table(build_geo_aux(shp, geo_layers(mask_precip = FALSE),
                                     venezuela_border(), weighted = TRUE))
  setnames(geo, geo_subarea_names(names(geo)))
  geo <- geo[, c("municipality", kept_sub), with = FALSE]
  # A municipality left without a value (a layer with no data over it) takes the
  # median over municipalities, as in the area-level geo block.
  for (k in kept_sub) {
    n_na <- sum(is.na(geo[[k]]))
    if (n_na) {
      set(geo, which(is.na(geo[[k]])), k, stats::median(geo[[k]], na.rm = TRUE))
      cat(sprintf("median-imputed %d NA in %s\n", n_na, k))
    }
  }
  stopifnot(nrow(geo) == 1122L, !anyDuplicated(geo$municipality), !anyNA(geo))
  p_out <- file.path(AUX, "municipal_geo_2018con.parquet")
  write_parquet_safe(tibble::as_tibble(geo), p_out)
  cat(sprintf("wrote %s: %d municipalities x %d geo columns\n", p_out, nrow(geo), length(kept_sub)))

} else stop("unknown mode: ", mode)
