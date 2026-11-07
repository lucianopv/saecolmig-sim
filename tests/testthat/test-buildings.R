test_that("bld_parse_manifest concatenates uriPrefix and uri with no separator", {
  d <- bld_parse_manifest(make_bld_manifest(), "89_EPSG_32618_2018_06_30.json")
  expect_equal(nrow(d), 2L)
  expect_identical(d$uri[1], "v1/geotiffs/8925c_2018_06_30/tile_A.tif")
  expect_identical(d$epsg, c(32618L, 32618L))
  expect_equal(d$xmax[1], 500000 + 25000 * 0.5)
  expect_equal(d$ymin[1], 1012500 - 25000 * 0.5)
  expect_equal(d$resy[1], 0.5)
})

test_that("bld_to_wgs84 adds a lon/lat box that contains the tile", {
  d <- bld_to_wgs84(bld_parse_manifest(make_bld_manifest(), "a_EPSG_32618_2018_06_30.json"))
  expect_true(all(c("lon_min", "lon_max", "lat_min", "lat_max") %in% names(d)))
  expect_true(all(d$lon_min < d$lon_max & d$lat_min < d$lat_max))
  expect_true(all(d$lon_min > -76 & d$lon_max < -74))   # zone 18 central band
})

test_that("bld_zones_for_bbox covers both hemispheres when the box crosses 0", {
  z <- bld_zones_for_bbox(c(-81.84, -4.23, -66.84, 15.91))
  expect_identical(z, c(32617L, 32618L, 32619L, 32717L, 32718L, 32719L))
  expect_identical(bld_zones_for_bbox(c(-75, 1, -74, 2)), 32618L)
})

test_that("bld_tile_index filters manifests by year and zone, via the injected fetcher", {
  listing <- list(items = list(
    list(name = "v1/manifests/89_EPSG_32618_2018_06_30.json"),
    list(name = "v1/manifests/89_EPSG_32618_2019_06_30.json"),
    list(name = "v1/manifests/10_EPSG_32650_2018_06_30.json")))
  seen <- character()
  fetch <- function(u) {
    seen <<- c(seen, u)
    if (grepl("/o\\?", u)) listing else make_bld_manifest()
  }
  d <- bld_tile_index(2018L, c(-75, 1, -74, 2), fetch = fetch)
  expect_equal(nrow(d), 2L)
  expect_equal(sum(grepl("manifests/89_EPSG_32618_2018", seen)), 1L)
  expect_false(any(grepl("2019|32650", seen[-1])))
})

test_that("bld_tile_id is zone-qualified and bld_tile_url is a vsicurl path", {
  expect_identical(bld_tile_id(32618L, "v1/geotiffs/8925c/tile_A.tif"), "32618_tile_A")
  expect_identical(bld_tile_url("v1/x.tif"),
    "/vsicurl/https://storage.googleapis.com/open-buildings-temporal-data/v1/x.tif")
})

test_that("bld_select_tiles keeps footprints that meet the buffered boundary", {
  # a 1-degree square in central Colombia (the buffer is applied in EPSG:3116,
  # which is only meaningful near Colombia)
  idx <- data.table::data.table(
    uri = c("a/t1.tif", "a/t2.tif", "a/t3.tif"), epsg = 32618L,
    lon_min = c(-74.5, -73.999, -70), lon_max = c(-74.1, -73.5, -69),
    lat_min = c(4.1, 4.1, 0), lat_max = c(4.4, 4.4, 1))
  sq <- sf::st_sf(geometry = sf::st_sfc(sf::st_polygon(list(rbind(
    c(-75, 4), c(-74, 4), c(-74, 5), c(-75, 5), c(-75, 4)))), crs = 4326))
  # t2 starts ~110 m east of the square: in with a 2 km buffer, out with 0
  expect_identical(bld_select_tiles(idx, sq, buffer_m = 2000)$tile_id,
                   c("32618_t1", "32618_t2"))
  expect_identical(bld_select_tiles(idx, sq, buffer_m = 0)$tile_id, "32618_t1")
})

test_that("bld_select_tiles refuses duplicate tile ids", {
  idx <- data.table::data.table(uri = c("a/t1.tif", "b/t1.tif"), epsg = 32618L,
                                lon_min = -74.5, lon_max = -74.4, lat_min = 4.1, lat_max = 4.2)
  sq <- sf::st_sf(geometry = sf::st_sfc(sf::st_polygon(list(rbind(
    c(-75, 4), c(-74, 4), c(-74, 5), c(-75, 5), c(-75, 4)))), crs = 4326))
  expect_error(bld_select_tiles(idx, sq), "duplicate tile_id")
})

test_that("bld_task_rows slices a static list into blocks", {
  expect_identical(bld_task_rows(60L, 1L, 25L), 1:25)
  expect_identical(bld_task_rows(60L, 3L, 25L), 51:60)
  expect_identical(bld_task_rows(60L, 4L, 25L), integer(0))
})

test_that("bld_tile_vrt exposes count, presence, p*h, p*h^2 with nodata propagated", {
  skip_if(!nzchar(Sys.which("gdalinfo")), "GDAL CLI not available")
  fx <- make_bld_tile()
  vrt <- bld_tile_vrt(fx$tile, tempfile(fileext = ".vrt"), src = fx$path)
  r <- terra::rast(vrt)
  expect_equal(terra::nlyr(r), 4L)
  v <- terra::values(r)                         # cell order: row-major
  cell <- function(row, col) (row - 1L) * 400L + col
  expect_equal(unname(v[cell(150, 50), ]), c(0.001, 0.5, 2, 8), tolerance = 1e-6)
  expect_equal(unname(v[cell(150, 350), ]), c(0.003, 0, 0, 0), tolerance = 1e-6)
  expect_true(all(is.na(v[cell(50, 350), ])))   # nodata block, all four bands
})

test_that("bld_warp_tile averages valid pixels only and writes atomically", {
  skip_if(!nzchar(Sys.which("gdalwarp")), "gdalwarp not available")
  fx <- make_bld_tile()
  out <- file.path(fx$dir, "out.tif")
  before_vrt <- list.files(tempdir(), "\\.vrt$")
  # same CRS, 100 m cells: an exact 200x200 block aggregation, 2 x 2 output
  st <- bld_warp_tile(fx$tile, out, t_srs = "EPSG:32618", res = 100, src = fx$path)
  expect_identical(st$status, "ok")
  expect_false(file.exists(paste0(out, ".part.tif")))
  expect_identical(list.files(tempdir(), "\\.vrt$"), before_vrt)  # no leftover temp VRT
  r <- terra::rast(out)
  expect_identical(names(r), bld_band_names)
  v <- terra::values(r)                         # TL, TR, BL, BR
  # TR: rows 1-200, cols 201-400; rows 1-100 x cols 301-400 are nodata (all
  # four bands), so for count 10,000 px at 0.001 and 20,000 px at 0.003
  # remain -> 0.0023333; presence/height are 0 throughout TR regardless
  # (right half of the tile), so p_h/p_h2 stay exactly 0, not NaN or leaked
  # garbage from the nodata sub-block.
  expect_equal(unname(v[2, ]), c((1e4 * 0.001 + 2e4 * 0.003) / 3e4, 0, 0, 0), tolerance = 1e-6)
  expect_equal(unname(v[1, ]), c(0.001, 0.5, 2, 8), tolerance = 1e-6)   # TL
  expect_equal(unname(v[4, ]), c(0.001, 0, 0, 0), tolerance = 1e-6)     # BR
  expect_true(all(v >= 0))                      # no -99 leak
  expect_false(anyNA(v))
})

test_that("bld_warp_tile masks each band by its own nodata sentinel, not all four at once", {
  skip_if(!nzchar(Sys.which("gdalwarp")), "gdalwarp not available")
  # Three non-overlapping, single-band nodata blocks inside the TL output
  # cell (rows 1-200, cols 1-200, the p > 0 / h > 0 region): height-only,
  # count-only, presence-only. If gdalwarp required every band to agree
  # before masking a pixel (UNIFIED_SRC_NODATA=YES, the default once
  # -srcnodata is passed), these partially-invalid pixels would be averaged
  # in as if valid, leaking -99/NaN-derived garbage into all four bands.
  fx <- make_bld_tile(nodata_blocks = list(
    list(rows = 1:50,    cols = 1:200, bands = "hgt"),
    list(rows = 51:100,  cols = 1:200, bands = "cnt"),
    list(rows = 101:150, cols = 1:200, bands = "pres")))
  out <- file.path(fx$dir, "partial.tif")
  st <- bld_warp_tile(fx$tile, out, t_srs = "EPSG:32618", res = 100, src = fx$path)
  expect_identical(st$status, "ok")
  r <- terra::rast(out)
  v <- terra::values(r)
  expect_false(anyNA(v))
  # Every band's own masking is independent, so every valid pixel in each
  # band is the same constant (0.001 / 0.5 / 2 / 8) and the block average
  # equals it exactly, regardless of which of the three blocks hit which band.
  expect_equal(unname(v[1, ]), c(0.001, 0.5, 2, 8), tolerance = 1e-6)
})

test_that("bld_warp_tile reports failure without leaving a file", {
  skip_if(!nzchar(Sys.which("gdalwarp")), "gdalwarp not available")
  fx <- make_bld_tile()
  out <- file.path(fx$dir, "bad.tif")
  st <- bld_warp_tile(fx$tile, out, t_srs = "EPSG:32618", res = 100,
                      src = file.path(fx$dir, "does_not_exist.tif"))
  expect_identical(st$status, "failed")
  expect_false(file.exists(out))
  expect_false(file.exists(paste0(out, ".part.tif")))
})

test_that("extract_weighted_mean_multi equals extract_weighted_mean band by band", {
  shp <- make_geo_shp()
  a <- make_geo_raster(base = 1); b <- make_geo_raster(base = 5) * 2
  st <- c(a, b); names(st) <- c("a", "b")
  w_same <- make_geo_raster(base = 0) + 3           # same grid
  w_diff <- terra::rast(xmin = 0, xmax = 3, ymin = 0, ymax = 1,
                        resolution = 0.125, crs = "EPSG:4326")
  terra::values(w_diff) <- seq_len(terra::ncell(w_diff))  # finer, non-uniform
  for (w in list(w_same, w_diff)) {
    m <- extract_weighted_mean_multi(shp, st, w)
    expect_named(m, c("municipality", "a", "b"))
    expect_equal(m$a, extract_weighted_mean(shp, a, w, "a")$a, tolerance = 1e-10)
    expect_equal(m$b, extract_weighted_mean(shp, b, w, "b")$b, tolerance = 1e-10)
  }
})

test_that("bld_derive_covariates recovers the presence-weighted mean and SD", {
  # two pixel groups inside one polygon: p = 0.2 at h = 3, p = 0.6 at h = 7,
  # equal pixel counts -> pw means are simple averages of the per-pixel products
  p <- c(0.2, 0.6); h <- c(3, 7)
  d <- bld_derive_covariates(pw_count = 0.001, pw_p = mean(p),
                             pw_ph = mean(p * h), pw_ph2 = mean(p * h^2))
  wm <- sum(p * h) / sum(p)
  wsd <- sqrt(sum(p * (h - wm)^2) / sum(p))
  expect_named(d, bld_cov_names)
  expect_equal(d$building_density_pw, 4000)          # 0.001 per 0.25 m2 -> per km2
  expect_equal(d$building_presence_pw, 0.4)
  expect_equal(d$building_height_pw, wm)
  expect_equal(d$building_height_sd_pw, wsd)
})

test_that("bld_derive_covariates: no buildings -> zeros; NA stays NA; rounding clamps", {
  d <- bld_derive_covariates(c(0, NA, 0.001), c(0, NA, 0.5), c(0, NA, 2),
                             c(0, NA, 8 - 1e-9))    # variance -4e-9 from rounding
  expect_equal(unlist(d[1, ]), setNames(c(0, 0, 0, 0), bld_cov_names))
  expect_true(all(is.na(unlist(d[2, ]))))
  expect_equal(d$building_height_sd_pw[3], 0)
})

test_that("bld_pw_covariates composes the two steps on a 4-band raster", {
  shp <- make_geo_shp()
  cnt <- make_geo_raster(base = 0) * 0 + 0.001
  pr  <- make_geo_raster(base = 0) * 0 + 0.5
  bld <- c(cnt, pr, pr * 4, pr * 16); names(bld) <- bld_band_names
  w   <- make_geo_raster(base = 0) + 1
  out <- bld_pw_covariates(shp, bld, w)
  expect_named(out, c("municipality", bld_cov_names))
  expect_equal(out$building_density_pw, rep(4000, 3), tolerance = 1e-6)
  expect_equal(out$building_height_pw, rep(4, 3), tolerance = 1e-9)
  expect_equal(out$building_height_sd_pw, rep(0, 3), tolerance = 1e-6)
})

test_that("bld_prune_vif drops only prunable columns, the worst first", {
  set.seed(1)
  n <- 400
  x1 <- rnorm(n); x2 <- rnorm(n)
  df <- data.frame(x1 = x1, x2 = x2, x3 = x1 + rnorm(n, sd = 0.01),  # base itself collinear
                   b1 = x1 + rnorm(n, sd = 0.02),                  # collinear with base
                   b2 = rnorm(n),
                   b3 = NA_real_)
  df$b3 <- 2 * df$b2 + rnorm(n, sd = 0.01)                          # collinear with b2
  res <- bld_prune_vif(df, base = c("x1", "x2", "x3"), prunable = c("b1", "b2", "b3"))
  expect_true(all(res$dropped$column %in% c("b1", "b2", "b3")))     # base never dropped
  expect_false("b1" %in% res$kept)
  expect_equal(sum(c("b2", "b3") %in% res$kept), 1L)
  expect_equal(res$dropped$step, seq_len(nrow(res$dropped)))
  expect_true(all(res$final_vif <= 10))
  expect_identical(names(res$final_vif), res$kept)
})

test_that("bld_prune_vif never drops the last prunable column (assemble_geo_X.R rule)", {
  set.seed(2)
  x1 <- rnorm(200)
  df <- data.frame(x1 = x1, b1 = x1 + rnorm(200, sd = 1e-3))
  res <- bld_prune_vif(df, base = "x1", prunable = "b1")
  expect_identical(res$kept, "b1")
  expect_gt(res$final_vif[["b1"]], 10)                              # recorded, not hidden
})

test_that("bld_prune_vif uses complete cases only", {
  df <- data.frame(x1 = c(1:9, NA), b1 = c(2, 1, 4, 3, 6, 5, 8, 7, 10, 9),
                   b2 = c(1, 3, 2, 5, 4, 7, 6, 9, 8, 10))
  expect_silent(res <- bld_prune_vif(df, base = "x1", prunable = c("b1", "b2")))
  expect_true(length(res$kept) >= 1L)
})

test_that("bld_prune_vif treats a non-finite VIF (e.g. a constant column) as infinite and drops it without erroring", {
  set.seed(3)
  n <- 100
  x1 <- rnorm(n)
  df <- data.frame(x1 = x1,
                   b1 = x1 + rnorm(n, sd = 0.05),  # normal, finite VIF
                   b2 = 0)                          # all-zero after complete cases ->
                                                     # exact 0/0 R^2 -> VIF non-finite
  expect_silent(res <- bld_prune_vif(df, base = "x1", prunable = c("b1", "b2")))
  expect_false("b2" %in% res$kept)                  # non-finite VIF is dropped first
  expect_identical(res$dropped$column[1], "b2")
  expect_true(length(res$kept) >= 1L)                # base untouched, a prunable survives
  expect_true("x1" %in% names(df))
})

test_that("bld_prune_vif treats a below-1 VIF (impossible; floating-point only) as infinite and drops it first", {
  # Reproducing R^2 > 1 from a real lm() fit is not reliable (see the helper
  # fixture comment), so stub stats::lm() to tag b2's fit with a class whose
  # summary() reports r.squared = 1 + 1e-6 -> VIF = 1/(1-1.000001) ~ -1e6, a
  # FINITE negative number. A guard on non-finite values alone would let it
  # through: smaller than b1's genuine ~1, b2 would never be the largest VIF,
  # so b1 would be dropped and the impossible b2 kept.
  set.seed(4)
  n <- 30
  df <- data.frame(x1 = rnorm(n), b1 = rnorm(n), b2 = rnorm(n))
  real_lm <- stats::lm
  testthat::local_mocked_bindings(
    lm = function(formula, data, ...) {
      fit <- real_lm(formula, data = data, ...)
      if (identical(all.vars(formula)[1], "b2"))
        class(fit) <- c("bld_test_impossible_lm", class(fit))
      fit
    },
    .package = "stats"
  )
  res <- bld_prune_vif(df, base = "x1", prunable = c("b1", "b2"))
  expect_identical(res$dropped$column[1], "b2")     # impossible VIF dropped FIRST
  expect_true(is.infinite(res$dropped$vif[1]))
  expect_true("b1" %in% res$kept)
})

test_that("bld_tile_src picks the local file only when a source dir is given", {
  t <- data.table::data.table(uri = "v1/geotiffs/89abc/tile_A.tif", epsg = 32618L,
                              tile_id = "32618_tile_A")
  expect_identical(bld_tile_src(t), bld_tile_url(t$uri))
  expect_identical(bld_tile_src(t, NULL), bld_tile_url(t$uri))
  expect_identical(bld_tile_src(t, ""), bld_tile_url(t$uri))
  expect_identical(bld_tile_src(t, "/data/src"), "/data/src/32618_tile_A.tif")
})

test_that("bld_warp_tile warps from a local source file", {
  skip_if(!nzchar(Sys.which("gdalwarp")), "gdalwarp not available")
  fx <- make_bld_tile()
  dir <- file.path(fx$dir, "src"); dir.create(dir, showWarnings = FALSE)
  file.copy(fx$path, file.path(dir, paste0(fx$tile$tile_id, ".tif")))
  out <- file.path(fx$dir, "local.tif")
  st <- bld_warp_tile(fx$tile, out, t_srs = "EPSG:32618", res = 100,
                      src = bld_tile_src(fx$tile, dir))
  expect_identical(st$status, "ok")
  v <- terra::values(terra::rast(out))
  expect_equal(unname(v[1, ]), c(0.001, 0.5, 2, 8), tolerance = 1e-6)
  expect_false(anyNA(v))
})

test_that("bld_pw_covariates on constant bands returns the covariates those bands imply", {
  shp <- make_geo_shp()
  cnt <- make_geo_raster(base = 0) * 0 + 0.001
  pr  <- make_geo_raster(base = 0) * 0 + 0.4
  # p_h and p_h2 consistent with a height of 5 m everywhere buildings stand
  bld <- c(cnt, pr, pr * 5, pr * 25); names(bld) <- bld_band_names
  w_flat <- make_geo_raster(base = 0) + 1
  pw <- bld_pw_covariates(shp, bld, w_flat)
  # density = 0.001 buildings per 0.25 m2 pixel x 4e6 = 4000 per km2; presence
  # 0.4; height 5 m; height SD 0. The density and the SD carry the floating-point
  # noise of exactextractr's weighted accumulation (the SD through the
  # E[X^2] - E[X]^2 cancellation near 0), so they get a looser tolerance than
  # presence and height.
  expect_equal(pw$building_density_pw, rep(4000, 3), tolerance = 1e-3)
  expect_equal(pw$building_presence_pw, rep(0.4, 3), tolerance = 1e-6)
  expect_equal(pw$building_height_pw, rep(5, 3), tolerance = 1e-6)
  expect_equal(pw$building_height_sd_pw, rep(0, 3), tolerance = 1e-3)
})
