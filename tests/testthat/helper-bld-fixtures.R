# Synthetic fixtures for the building-raster helpers. No network, no real data.

# A manifest list with the same shape as the real Open Buildings manifests
# (for example 89_EPSG_32618_2018_06_30.json).
make_bld_manifest <- function() {
  src <- function(uri, tx, ty) list(
    uris = list(uri),
    affineTransform = list(scaleX = 0.5, translateX = tx, scaleY = -0.5, translateY = ty),
    dimensions = list(width = 25000L, height = 25000L))
  list(name = "x", uriPrefix = "gs://open-buildings-temporal-data/v1/geotiffs/89",
       tilesets = list(list(id = "", sources = list(
         src("25c_2018_06_30/tile_A.tif", 500000, 1012500),
         src("25c_2018_06_30/tile_B.tif", 512500, 1012500)))))
}

# A 400 x 400 native tile (200 m square) in EPSG:32618, origin on the 100 m
# lattice, with a nodata block. Returns list(path, tile) where `tile` is a
# one-row index record pointing at it.
#   presence p = 0.5, height h = 4  in columns 1-200 (x < 500100), else 0
#   count      = 0.001 everywhere valid, except 0.003 in rows 101-200 x cols 201-400
#   nodata     = -99 in rows 1-100 x cols 301-400, all three bands (default)
#
# `nodata_blocks` is a list of list(rows=, cols=, bands=) entries, each setting
# -99 on just the named band(s) ("cnt"/"pres"/"hgt") over that row/col range
# — used to test partial (band-disagreeing) nodata masks. The default
# reproduces the single, all-three-bands block above exactly.
make_bld_tile <- function(dir = tempfile("bldtile_"),
                          nodata_blocks = list(list(rows = 1:100, cols = 301:400,
                                                     bands = c("cnt", "pres", "hgt")))) {
  dir.create(dir, showWarnings = FALSE)
  n <- 400L
  cnt <- matrix(0.001, n, n); pres <- matrix(0, n, n); hgt <- matrix(0, n, n)
  pres[, 1:200] <- 0.5; hgt[, 1:200] <- 4
  cnt[101:200, 201:400] <- 0.003
  mats <- list(cnt = cnt, pres = pres, hgt = hgt)
  for (blk in nodata_blocks) {
    for (m in blk$bands) mats[[m]][blk$rows, blk$cols] <- -99
  }
  cnt <- mats$cnt; pres <- mats$pres; hgt <- mats$hgt
  r <- terra::rast(nrows = n, ncols = n, nlyrs = 3, xmin = 500000, xmax = 500200,
                   ymin = 1000000, ymax = 1000200, crs = "EPSG:32618")
  terra::values(r) <- cbind(as.vector(t(cnt)), as.vector(t(hgt)), as.vector(t(pres)))
  names(r) <- c("building_fractional_count", "building_height", "building_presence")
  path <- file.path(dir, "tile.tif")
  terra::writeRaster(r, path, NAflag = -99, datatype = "FLT4S", overwrite = TRUE)
  tile <- data.table::data.table(uri = "local/tile.tif", epsg = 32618L, resx = 0.5,
                                 resy = 0.5, nx = n, ny = n, xmin = 500000,
                                 xmax = 500200, ymin = 1000000, ymax = 1000200,
                                 tile_id = "32618_tile")
  list(path = path, tile = tile, dir = dir)
}

# Support for the "impossible VIF" test in test-buildings.R: a VIF < 1 can
# only arise from an R^2 slightly ABOVE 1, which R's own summary.lm() cannot
# produce (r.squared = mss/(mss+rss), a ratio of two literal sums of squares,
# is bounded to [0,1] even under floating-point rounding). So the test fabricates
# the impossible value directly via a stubbed stats::lm() (see
# `local_mocked_bindings(lm = ..., .package = "stats")`) that tags one
# column's fit with this S3 class; this method intercepts `summary()` on it
# and reports r.squared = 1 + 1e-6. Registered explicitly into .GlobalEnv
# (not a plain top-level `<-`) because S3 dispatch from inside the package
# namespace falls back to .GlobalEnv specifically, not to whatever
# environment testthat happens to evaluate this helper file in.
assign("summary.bld_test_impossible_lm",
       function(object, ...) list(r.squared = 1 + 1e-6),
       envir = globalenv())
