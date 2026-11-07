# Google Open Buildings 2.5D Temporal -> a national raster on the WorldPop
# 3-arcsec grid, and the population-weighted building covariates built from it.
# The pipeline is data-raw/build_buildings_raster.R (raster),
# data-raw/build_area_geo_bld.R (municipal table, area-level matrix) and
# data-raw/build_subarea_geo_blk_bld.R (sub-area table, sub-area matrix).
#
# The tile-index code (manifest listing and parsing, UTM zone helpers) is
# vendored from GeoLink (github.com/SSA-Statistical-Team-Projects/GeoLink),
# branch `buildings-temporal`, commit a1efba53. Nothing else from that package
# is used.

.bld_bucket_api <- "https://storage.googleapis.com/storage/v1/b/open-buildings-temporal-data/o"
.bld_base       <- "https://storage.googleapis.com/open-buildings-temporal-data/"

#' Parse one Open Buildings manifest into tile footprints (native UTM)
#' @param m parsed manifest (jsonlite::fromJSON(simplifyVector = FALSE))
#' @param name manifest file name, carrying `EPSG_<code>_`
#' @export
bld_parse_manifest <- function(m, name) {
  epsg <- as.integer(sub(".*EPSG_(\\d+)_.*", "\\1", name))
  # uriPrefix and the per-source uri concatenate with NO separator; a
  # file.path() join 404s.
  pref <- sub("^gs://open-buildings-temporal-data/", "", m$uriPrefix)
  data.table::rbindlist(lapply(m$tilesets[[1]]$sources, function(s) {
    a <- s$affineTransform; d <- s$dimensions
    data.table::data.table(
      uri  = paste0(pref, s$uris[[1]]), epsg = epsg,
      resx = a$scaleX, resy = abs(a$scaleY),
      nx = as.integer(d$width), ny = as.integer(d$height),
      xmin = a$translateX, xmax = a$translateX + d$width * a$scaleX,
      ymax = a$translateY, ymin = a$translateY - d$height * abs(a$scaleY))
  }))
}

#' Add WGS84 bounding boxes to a tile table
#' @export
bld_to_wgs84 <- function(dt) {
  data.table::rbindlist(lapply(split(dt, dt$epsg), function(g) {
    pts <- sf::st_as_sf(
      data.frame(x = c(g$xmin, g$xmax, g$xmin, g$xmax),
                 y = c(g$ymin, g$ymin, g$ymax, g$ymax),
                 id = rep(seq_len(nrow(g)), 4)),
      coords = c("x", "y"), crs = g$epsg[1])
    ll <- sf::st_coordinates(sf::st_transform(pts, 4326))
    d  <- data.table::data.table(id = pts$id, X = ll[, 1], Y = ll[, 2])
    b  <- d[, list(lon_min = min(X), lon_max = max(X),
                   lat_min = min(Y), lat_max = max(Y)), by = "id"]
    data.table::setorderv(b, "id")
    cbind(g, b[, !"id", with = FALSE])
  }))
}

#' UTM zone EPSG codes a WGS84 bbox c(xmin, ymin, xmax, ymax) can fall in
#' @export
bld_zones_for_bbox <- function(bbox) {
  z <- floor((c(bbox[1], bbox[3]) + 180) / 6) + 1
  z <- seq(min(z), max(z))
  out <- integer(0)
  if (bbox[4] >= 0) out <- c(out, 32600L + z)
  if (bbox[2] <  0) out <- c(out, 32700L + z)
  sort(unique(as.integer(out)))
}

#' Tile index for one year, restricted to the zones a bbox needs
#'
#' @param fetch function(url) returning the parsed JSON; injectable for tests
#' @export
bld_tile_index <- function(year, bbox,
                           fetch = function(u) jsonlite::fromJSON(u, simplifyVector = FALSE)) {
  nms <- character(0); tok <- NULL
  repeat {
    u <- paste0(.bld_bucket_api, "?prefix=v1/manifests/&maxResults=1000",
                if (!is.null(tok)) paste0("&pageToken=", tok) else "")
    j <- fetch(u)
    nms <- c(nms, vapply(j$items, function(x) basename(x$name), character(1)))
    tok <- j$nextPageToken
    if (is.null(tok)) break
  }
  nms <- nms[grepl(paste0("_", year, "_06_30\\.json$"), nms)]
  z   <- bld_zones_for_bbox(bbox)
  nms <- nms[grepl(paste0("EPSG_(", paste(z, collapse = "|"), ")_"), nms)]
  if (!length(nms)) stop("no Open Buildings manifests for ", year, " in zones ",
                         paste(z, collapse = ","))
  idx <- data.table::rbindlist(lapply(nms, function(nm)
    bld_parse_manifest(fetch(paste0(.bld_base, "v1/manifests/", nm)), nm)))
  bld_to_wgs84(idx)
}

#' Zone-qualified tile identifier (tile basenames are only unique per zone)
#' @export
bld_tile_id <- function(epsg, uri) {
  sprintf("%d_%s", as.integer(epsg), sub("\\.tif$", "", basename(uri)))
}

#' GDAL path for a tile uri from the index
#' @export
bld_tile_url <- function(uri) paste0("/vsicurl/", .bld_base, uri)

#' Tiles whose WGS84 footprint meets a boundary buffered by `buffer_m`
#'
#' The buffer is applied in EPSG:3116 (MAGNA-SIRGAS / Colombia Bogota), a
#' metric CRS, so "2 km" means 2 km.
#' @param boundary sf polygon(s), any CRS
#' @export
bld_select_tiles <- function(idx, boundary, buffer_m = 2000) {
  idx <- data.table::copy(data.table::as.data.table(idx))
  idx[, tile_id := bld_tile_id(epsg, uri)]
  if (anyDuplicated(idx$tile_id)) stop("duplicate tile_id in the index")
  foot <- sf::st_as_sfc(lapply(seq_len(nrow(idx)), function(i) sf::st_polygon(list(rbind(
    c(idx$lon_min[i], idx$lat_min[i]), c(idx$lon_max[i], idx$lat_min[i]),
    c(idx$lon_max[i], idx$lat_max[i]), c(idx$lon_min[i], idx$lat_max[i]),
    c(idx$lon_min[i], idx$lat_min[i]))))), crs = 4326)
  b <- sf::st_union(sf::st_transform(boundary, 3116))
  if (buffer_m > 0) b <- sf::st_buffer(b, buffer_m)
  b <- sf::st_transform(b, 4326)
  old <- sf::sf_use_s2(FALSE); on.exit(sf::sf_use_s2(old), add = TRUE)
  hit <- lengths(sf::st_intersects(foot, b)) > 0
  idx[hit]
}

#' Row indices of one block of a static tile list
#' @export
bld_task_rows <- function(n, task_id, block) {
  from <- (as.integer(task_id) - 1L) * as.integer(block) + 1L
  if (from > n) return(integer(0))
  seq.int(from, min(from + as.integer(block) - 1L, n))
}

#' Output band order of every warped tile and of the national raster
#' @export
bld_band_names <- c("count", "presence", "p_h", "p_h2")

#' Set GDAL network/cache defaults, leaving any value already set alone
#' @export
bld_gdal_env <- function() {
  def <- c(GDAL_HTTP_MAX_RETRY = "5", GDAL_HTTP_RETRY_DELAY = "2",
           GDAL_DISABLE_READDIR_ON_OPEN = "EMPTY_DIR", GDAL_CACHEMAX = "512")
  unset <- def[!nzchar(Sys.getenv(names(def)))]
  if (length(unset)) do.call(Sys.setenv, as.list(unset))
  invisible(unset)
}

#' Path of a GDAL command-line tool, or an error naming it
#' @export
bld_gdal_bin <- function(name) {
  p <- Sys.which(name)
  if (!nzchar(p)) stop(name, " not found on PATH (the GDAL command-line tools are required)")
  unname(p)
}

#' Write the 4-band VRT for one tile
#'
#' Bands: 1 count (file band 1), 2 presence (file band 3), 3 p*h and 4 p*h^2,
#' computed by GDAL's `mul` pixel function directly on the native tile's
#' presence/height bands. File band order (1 count, 2 height, 3 presence)
#' was read from the tiles' band descriptions.
#'
#' `p_h`/`p_h2` declare `<NoDataValue>nan</NoDataValue>` (not -99, unlike
#' `count`/`presence`), while their input `ComplexSource`s still carry
#' `<NODATA>-99</NODATA>` as usual: GDAL fills a derived band's per-source
#' read buffer with the *band's own* declared NoDataValue before a skipped
#' -99 source pixel would otherwise leave the literal raw value in it, so
#' that pixel enters `mul()` as NaN. This is an implementation detail pinned
#' by tests/testthat/test-buildings.R; switching a source to a plain `SimpleSource` (i.e.
#' dropping its explicit `<NODATA>`) would break it on GDAL >= 3.12 (next
#' paragraph). Below GDAL 3.12, `mul` has no NoData awareness of its own —
#' it just multiplies whatever it is handed, so it silently ignores the
#' `propagateNoData` argument, and `NaN * anything` is NaN anyway, because
#' that is ordinary floating-point arithmetic, not a GDAL feature. Starting
#' with GDAL 3.12, `mul` gained a built-in NoData value of its own (read
#' from the derived band's declared NoDataValue, i.e. nan) and, by default,
#' *skips* any operand equal to it — treats it as a factor of 1 — instead of
#' propagating it, so a nodata height or presence pixel would silently
#' corrupt `p_h`/`p_h2` on 3.12+ without `propagateNoData="true"`. The
#' argument is therefore added unconditionally: a no-op below 3.12,
#' load-bearing at and above it. `count`/`presence` end up with the raw -99
#' sentinel at an invalid pixel, while `p_h`/`p_h2` end up with NaN there —
#' `bld_warp_tile` masks each band against its own sentinel (see its doc) so
#' a pixel invalid in only one band does not contaminate the other three.
#' This flat, single-VRT design is markedly cheaper per tile than routing
#' height/presence through a sibling `replace_nodata` VRT first.
#' @param src GDAL path of the native tile (defaults to its /vsicurl/ URL)
#' @export
bld_tile_vrt <- function(tile, path, src = bld_tile_url(tile$uri)) {
  # A relative `src` (e.g. a staged tile passed as a relative path) would
  # bake that relative path literally into the VRT XML, which then only
  # opens correctly from the one working directory it was written in. A
  # /vsicurl/ URL is left alone; anything else is resolved to absolute.
  if (!startsWith(src, "/vsicurl/")) src <- normalizePath(src, mustWork = FALSE)
  srcx <- function(b) sprintf(paste0(
    '    <ComplexSource><SourceFilename relativeToVRT="0">%s</SourceFilename>',
    '<SourceBand>%d</SourceBand>\n',
    '      <SrcRect xOff="0" yOff="0" xSize="%d" ySize="%d"/>',
    '<DstRect xOff="0" yOff="0" xSize="%d" ySize="%d"/><NODATA>-99</NODATA></ComplexSource>'),
    src, b, tile$nx, tile$ny, tile$nx, tile$ny)
  band <- function(i, desc, bands, fn = NULL, nodata = "-99") {
    head <- if (is.null(fn)) sprintf('  <VRTRasterBand dataType="Float32" band="%d">', i) else
      sprintf('  <VRTRasterBand dataType="Float32" band="%d" subClass="VRTDerivedRasterBand">', i)
    fnx <- if (is.null(fn)) "" else sprintf(paste0(
      '\n    <PixelFunctionType>%s</PixelFunctionType>',
      '<PixelFunctionArguments propagateNoData="true"/>'), fn)
    paste0(head,
           sprintf('\n    <Description>%s</Description><NoDataValue>%s</NoDataValue>', desc, nodata),
           fnx, "\n",
           paste(vapply(bands, srcx, ""), collapse = "\n"),
           "\n  </VRTRasterBand>")
  }
  xml <- paste0(
    sprintf('<VRTDataset rasterXSize="%d" rasterYSize="%d">\n', tile$nx, tile$ny),
    sprintf('  <SRS>EPSG:%d</SRS>\n', tile$epsg),
    sprintf('  <GeoTransform>%.6f, %.6f, 0, %.6f, 0, %.6f</GeoTransform>\n',
            tile$xmin, tile$resx, tile$ymax, -tile$resy),
    band(1L, "count", 1L), "\n",
    band(2L, "presence", 3L), "\n",
    band(3L, "p_h", c(3L, 2L), "mul", nodata = "nan"), "\n",
    band(4L, "p_h2", c(3L, 2L, 2L), "mul", nodata = "nan"), "\n</VRTDataset>\n")
  writeLines(xml, path)
  path
}

#' Warp one tile's 4-band VRT onto a regular grid, atomically
#'
#' Writes `<out_tif>.part.tif`, checks it opens with four bands, then renames.
#' A failure leaves no file behind and is returned, never raised, so a block of
#' tiles can finish its other tiles and exit non-zero at the end.
#'
#' `count`/`presence` carry the raw -99 sentinel; `p_h`/`p_h2` carry NaN (see
#' `bld_tile_vrt`), hence the per-band `-srcnodata` list below ("-99 -99 nan
#' nan", matching `bld_band_names`) rather than a single value. That alone is
#' not enough: passing `-srcnodata` makes gdalwarp default
#' `UNIFIED_SRC_NODATA` to `YES`, under which a pixel is masked only when
#' *every* band matches its nodata value at once, so a pixel invalid in just
#' one band (e.g. height only) would be averaged into all four bands as if
#' valid. `-wo UNIFIED_SRC_NODATA=NO` restores per-band masking, so `-r
#' average` excludes an invalid pixel from only the band(s) where it is
#' actually invalid.
#' @export
bld_warp_tile <- function(tile, out_tif, t_srs = "EPSG:4326", res = 1 / 1200,
                          src = bld_tile_url(tile$uri),
                          gdalwarp = bld_gdal_bin("gdalwarp")) {
  t0   <- Sys.time()
  vrt  <- bld_tile_vrt(tile, tempfile(fileext = ".vrt"), src = src)
  part <- paste0(out_tif, ".part.tif")
  on.exit(unlink(c(vrt, if (!file.exists(out_tif)) part)), add = TRUE)
  unlink(part)
  args <- c("-overwrite", "-q", "-t_srs", t_srs,
            "-tr", sprintf("%.17g", res), sprintf("%.17g", res), "-tap",
            "-r", "average", "-srcnodata", shQuote("-99 -99 nan nan"),
            "-dstnodata", "-99", "-wo", "UNIFIED_SRC_NODATA=NO", "-wm", "512",
            "-co", "COMPRESS=DEFLATE", "-co", "PREDICTOR=3", "-co", "TILED=YES",
            shQuote(vrt), shQuote(part))
  out <- suppressWarnings(system2(gdalwarp, args, stdout = TRUE, stderr = TRUE))
  st  <- attr(out, "status"); if (is.null(st)) st <- 0L
  ok  <- st == 0L && file.exists(part) &&
    isTRUE(tryCatch(terra::nlyr(terra::rast(part)) == 4L, error = function(e) FALSE))
  msg <- if (ok) "" else paste(utils::tail(out, 3), collapse = " | ")
  if (ok) {
    ok <- file.rename(part, out_tif)
    if (!ok) msg <- "rename failed"
  } else unlink(part)
  data.frame(tile_id = if (is.null(tile$tile_id)) NA_character_ else tile$tile_id,
             status = if (ok) "ok" else "failed",
             secs = as.numeric(difftime(Sys.time(), t0, units = "secs")),
             message = msg, stringsAsFactors = FALSE)
}

#' Names of the four building covariates, in output order
#' @export
bld_cov_names <- c("building_density_pw", "building_presence_pw",
                   "building_height_pw", "building_height_sd_pw")

#' Building covariates from population-weighted band means
#'
#' pw_ph / pw_p is the presence-weighted mean height and
#' pw_ph2 / pw_p - height^2 the presence-weighted variance, exactly, because
#' all four means share the same weights. Density converts the mean fractional
#' count per 0.25 m2 native pixel to buildings per km2 (x 4e6).
#' A polygon with pw_p == 0 has no buildings where people live: height and SD
#' are 0 (not NA) so no unit is dropped. NA inputs stay NA (zero weight or no
#' data); the caller backfills them. Rounding can make the variance slightly
#' negative; it is clamped at 0.
#' @export
bld_derive_covariates <- function(pw_count, pw_p, pw_ph, pw_ph2) {
  has <- !is.na(pw_p)
  pos <- has & pw_p > 0
  height <- ifelse(pos, pw_ph / pw_p, ifelse(has, 0, NA_real_))
  v      <- ifelse(pos, pw_ph2 / pw_p - height^2, ifelse(has, 0, NA_real_))
  out <- data.frame(pw_count * 4e6, pw_p, height, sqrt(pmax(v, 0)))
  names(out) <- bld_cov_names
  out
}

#' Population-weighted building covariates over polygons
#'
#' @param shp sf with an id column named `municipality`
#' @param bld 4-layer SpatRaster in `bld_band_names` order
#' @param weights WorldPop SpatRaster
#' @export
bld_pw_covariates <- function(shp, bld, weights) {
  stopifnot(terra::nlyr(bld) == 4L)
  names(bld) <- bld_band_names
  m <- extract_weighted_mean_multi(shp, bld, weights)
  d <- bld_derive_covariates(m$count, m$presence, m$p_h, m$p_h2)
  tibble::as_tibble(cbind(m["municipality"], d))
}

#' GDAL path of a tile: a local copy when one is staged, else the bucket
#'
#' Where there is no outbound internet access, the tiles can be staged locally
#' beforehand: `fetch` mode stages them as `<src_dir>/<tile_id>.tif`;
#' everything else keeps streaming from `/vsicurl/`.
#' @export
bld_tile_src <- function(tile, src_dir = NULL) {
  if (is.null(src_dir) || !nzchar(src_dir)) return(bld_tile_url(tile$uri))
  file.path(src_dir, paste0(tile$tile_id, ".tif"))
}

#' VIF-prune a block of columns against a fixed base
#'
#' vif_j = 1 / (1 - R^2_j), regressing prunable column j on the base plus the
#' other remaining prunable columns, on complete cases. While the largest VIF
#' exceeds `threshold`, drop that column; stop with one column left. This is
#' the rule of the geo and buildings selections (see [kept_cols()]). Base
#' columns are never dropped. Final VIFs of the kept
#' columns are recomputed and returned, including one above the threshold if
#' it is the survivor. A non-finite VIF (NaN/Inf, e.g. a 0/0 R^2 from a
#' constant or all-zero prunable column after complete-case filtering) is
#' treated as infinite so it is dropped first, rather than erroring in the
#' `max(v) <= threshold` comparison. A VIF BELOW 1 is likewise impossible
#' (VIF = 1/(1-R^2) with R^2 in \[0,1\] for a genuine fit): it can only arise
#' from floating-point error giving an R^2 slightly above 1, which yields a
#' finite NEGATIVE number. Left unguarded, that negative value looks smaller
#' than every honest positive VIF, so `which.max()` could never select it for
#' removal -- an impossible-but-small-looking column would then be treated as
#' the LEAST collinear one and kept, silently, in every downstream model
#' built from this selection. Both non-finite and below-1 values are forced
#' to infinite so they are dropped first, exactly like a genuine Inf.
#' @export
bld_prune_vif <- function(df, base, prunable, threshold = 10) {
  d <- df[stats::complete.cases(df[, c(base, prunable)]), c(base, prunable), drop = FALSE]
  vif_of <- function(target, preds) {
    f  <- stats::reformulate(sprintf("`%s`", preds), response = as.name(target))
    r2 <- suppressWarnings(summary(stats::lm(f, data = d))$r.squared)
    1 / (1 - r2)
  }
  vifs <- function(cols) {
    v <- vapply(cols, function(g) vif_of(g, c(base, setdiff(cols, g))), numeric(1))
    v[!is.finite(v) | v < 1] <- Inf
    v
  }
  cols <- prunable
  dropped <- data.frame(step = integer(0), column = character(0), vif = numeric(0))
  repeat {
    if (length(cols) <= 1L) break
    v <- vifs(cols)
    if (max(v) <= threshold) break
    k <- names(which.max(v))
    dropped <- rbind(dropped, data.frame(step = nrow(dropped) + 1L, column = k,
                                         vif = max(v)))
    cols <- setdiff(cols, k)
  }
  list(kept = cols, dropped = dropped, final_vif = vifs(cols))
}
