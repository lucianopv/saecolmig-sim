test_that("extract_raster_stat returns area-weighted mean per municipality", {
  shp <- make_geo_shp()
  r   <- make_geo_raster(base = 1)            # munis -> 1, 11, 21
  out <- extract_raster_stat(shp, r, "x", fun = "mean")

  expect_named(out, c("municipality", "x"))
  expect_equal(nrow(out), 3L)
  expect_equal(out$municipality, c(5001L, 5002L, 5003L))
  expect_equal(out$x, c(1, 11, 21), tolerance = 1e-6)
})

test_that("extract_raster_stat supports sum", {
  shp <- make_geo_shp()
  r   <- terra::rast(xmin = 0, xmax = 3, ymin = 0, ymax = 1,
                     resolution = 0.5, crs = "EPSG:4326")
  terra::values(r) <- 1                       # every cell = 1
  out <- extract_raster_stat(shp, r, "n", fun = "sum")
  # each 1x1 muni covers 2x2 = 4 half-degree cells -> sum 4
  expect_equal(out$n, c(4, 4, 4), tolerance = 1e-6)
})

test_that("extract_raster_stat sum returns NA when a polygon has only NoData", {
  shp <- make_geo_shp()
  r <- make_geo_raster(); terra::values(r) <- NA          # all NoData
  out <- extract_raster_stat(shp, r, "n", fun = "sum")
  expect_true(all(is.na(out$n)))
})

test_that("extract_weighted_mean equals the area mean under uniform weights", {
  shp <- make_geo_shp()
  r   <- make_geo_raster(base = 1)                 # constant per muni: 1, 11, 21
  w   <- make_geo_raster(); terra::values(w) <- 1  # uniform population
  out <- extract_weighted_mean(shp, r, w, "x_pw")
  expect_named(out, c("municipality", "x_pw"))
  expect_equal(out$x_pw, c(1, 11, 21), tolerance = 1e-6)
})

test_that("extract_weighted_mean shifts toward where population concentrates", {
  shp <- make_geo_shp()
  r <- make_geo_raster(); xs <- terra::xFromCell(r, seq_len(terra::ncell(r)))
  terra::values(r) <- ifelse((xs %% 1) >= 0.5, 10, 0)  # within-muni: left 0, right 10
  w <- make_geo_raster(); terra::values(w) <- ifelse((xs %% 1) >= 0.5, 1, 0)
  area <- extract_raster_stat(shp, r, "x", "mean")
  pw   <- extract_weighted_mean(shp, r, w, "x_pw")
  expect_equal(area$x,  c(5, 5, 5),    tolerance = 1e-6)   # half 0, half 10
  expect_equal(pw$x_pw, c(10, 10, 10), tolerance = 1e-6)   # weight only on the 10s
})

test_that("extract_weighted_mean aligns a non-conforming weights grid", {
  shp <- make_geo_shp()
  r   <- make_geo_raster(base = 1)             # constant per muni: 1, 11, 21
  # offset origin + finer res + wider extent -> non-nesting grid; would raise
  # "Incompatible extents" without the crop + resample alignment.
  w <- terra::rast(xmin = -0.07, xmax = 3.07, ymin = -0.07, ymax = 1.07,
                   resolution = 0.1, crs = "EPSG:4326")
  terra::values(w) <- 1
  out <- extract_weighted_mean(shp, r, w, "x_pw")
  # value is constant within each muni, so any weighting reproduces 1, 11, 21
  expect_equal(out$x_pw, c(1, 11, 21), tolerance = 1e-6)
})

test_that("extract_weighted_mean treats no-data weight cells as zero, not NA", {
  shp <- make_geo_shp()
  r   <- make_geo_raster(base = 1)             # constant per muni: 1, 11, 21
  w   <- make_geo_raster(); terra::values(w) <- 1   # 12x4 = 48 cells
  w[c(1, 2, 30)] <- NA                          # like WorldPop ocean / border no-data
  out <- extract_weighted_mean(shp, r, w, "x_pw")
  expect_false(anyNA(out$x_pw))                 # no NA-propagation from no-data
  expect_equal(out$x_pw, c(1, 11, 21), tolerance = 1e-6)
})

test_that("terrain_ruggedness is ~0 over a flat DEM and positive over a rough one", {
  shp  <- make_geo_shp()
  flat <- make_geo_raster(base = 100); terra::values(flat) <- 100
  out_flat <- terrain_ruggedness(shp, flat)
  expect_named(out_flat, c("municipality", "terrain_ruggedness"))
  expect_true(all(out_flat$terrain_ruggedness < 1e-6))

  rough <- flat
  set.seed(1)
  terra::values(rough) <- 100 + 500 * stats::runif(terra::ncell(rough))
  out_rough <- terrain_ruggedness(shp, rough)
  expect_true(all(out_rough$terrain_ruggedness > 0))
})

test_that("precip_anomaly is the muni-mean difference of the two rasters", {
  shp    <- make_geo_shp()
  p2018  <- make_geo_raster(base = 1200); terra::values(p2018) <- 1200
  normal <- make_geo_raster(base = 1000); terra::values(normal) <- 1000
  out <- precip_anomaly(shp, p2018, normal)
  expect_named(out, c("municipality", "precip_anomaly_2018"))
  expect_equal(out$precip_anomaly_2018, c(200, 200, 200), tolerance = 1e-6)
})

test_that("precip_anomaly resamples a mismatched-grid normal", {
  shp    <- make_geo_shp()
  p2018  <- terra::rast(xmin = 0, xmax = 3, ymin = 0, ymax = 1,
                        resolution = 0.25, crs = "EPSG:4326")
  terra::values(p2018) <- 1200
  normal <- terra::rast(xmin = 0, xmax = 3, ymin = 0, ymax = 1,
                        resolution = 0.5, crs = "EPSG:4326")   # coarser grid
  terra::values(normal) <- 1000
  out <- precip_anomaly(shp, p2018, normal)
  expect_equal(out$precip_anomaly_2018, c(200, 200, 200), tolerance = 1e-6)
})

test_that("precip_anomaly reprojects a normal in a different CRS", {
  shp   <- make_geo_shp()
  p2018 <- make_geo_raster(); terra::values(p2018) <- 1200
  normal <- terra::project(make_geo_raster(), "EPSG:3857")   # different CRS
  terra::values(normal) <- 1000
  out <- precip_anomaly(shp, p2018, normal)
  expect_equal(out$precip_anomaly_2018, c(200, 200, 200), tolerance = 5)
})

test_that("dist_to_border is ~0 for a muni on the border and grows with distance", {
  shp <- make_geo_shp()                       # munis at x in [0,1],[1,2],[2,3]
  # vertical border line at x = 3 (right edge, touching muni 5003)
  border <- sf::st_sf(geometry = sf::st_sfc(
    sf::st_linestring(rbind(c(3, -1), c(3, 2))), crs = 4326))
  out <- dist_to_border(shp, border)
  expect_named(out, c("municipality", "dist_venezuela_border_km"))
  # centroid distances: 5003 (x=2.5) closest, 5001 (x=0.5) farthest
  expect_lt(out$dist_venezuela_border_km[3], out$dist_venezuela_border_km[2])
  expect_lt(out$dist_venezuela_border_km[2], out$dist_venezuela_border_km[1])
  expect_true(all(out$dist_venezuela_border_km > 0))
})

test_that("build_geo_aux returns the full contract, key-ordered, no extra NA", {
  shp <- make_geo_shp()
  border <- sf::st_sf(geometry = sf::st_sfc(
    sf::st_linestring(rbind(c(3, -1), c(3, 2))), crs = 4326))

  mk <- function(v) { r <- make_geo_raster(); terra::values(r) <- v; r }
  layers <- list(
    ntl          = mk(5),
    travel_time  = mk(120),
    dem          = mk(800),
    builtup      = mk(0.1),
    cropland     = mk(0.3),
    forest       = mk(0.4),
    precip_2018  = mk(1200),
    precip_normal= mk(1000),
    pop_2018     = mk(50),
    pop_2015     = mk(40)
  )

  out <- build_geo_aux(shp, layers, border)

  expect_equal(out$municipality, c(5001L, 5002L, 5003L))  # ascending, guarded
  expect_true(all(c("ntl_mean", "travel_time_city", "terrain_ruggedness",
                    "dist_venezuela_border_km", "built_up_share",
                    "pop_density_2018", "cropland_share", "forest_share",
                    "precip_annual_2018", "precip_anomaly_2018",
                    "pop_change_2015_2018") %in% names(out)))
  expect_equal(nrow(out), 3L)
  expect_false(anyNA(out$ntl_mean))
  expect_equal(out$precip_anomaly_2018, c(200, 200, 200), tolerance = 1e-6)
  expect_true(all(out$pop_change_2015_2018 > 0))          # 50 vs 40 -> +0.25

  # population-weighted twins exist; for constant-per-muni inputs they match the
  # area-weighted values (every cell in a muni is equal, so weighting is a no-op)
  pw_cols <- c("ntl_mean_pw", "travel_time_city_pw", "terrain_ruggedness_pw",
               "built_up_share_pw", "cropland_share_pw", "forest_share_pw",
               "precip_annual_2018_pw", "precip_anomaly_2018_pw")
  expect_true(all(pw_cols %in% names(out)))
  expect_equal(out$ntl_mean_pw, out$ntl_mean, tolerance = 1e-6)
  expect_equal(out$precip_anomaly_2018_pw, c(200, 200, 200), tolerance = 1e-6)
  # weighted = FALSE drops the _pw columns
  out_unw <- build_geo_aux(shp, layers, border, weighted = FALSE)
  expect_false(any(pw_cols %in% names(out_unw)))
})
