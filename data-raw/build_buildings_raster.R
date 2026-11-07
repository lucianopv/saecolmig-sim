#!/usr/bin/env Rscript
# data-raw/build_buildings_raster.R
# Google Open Buildings 2.5D Temporal -> a national 4-band raster (count,
# presence, p*h, p*h^2) on the constrained WorldPop 3-arcsec grid, written to
# geo_cache/buildings_<year>_3as.tif. The covariates are then taken from it
# by build_area_geo_bld.R (municipal) and build_subarea_geo_blk_bld.R
# (sub-area). Functions: R/buildings.R.
#
# Each tile (9,750 for Colombia in 2018) is warped on its own, so `run` takes
# one block of 25 tiles (390 blocks, about 1 h each, about 2.3 min a tile) and
# the blocks can run separately or in parallel. Finished tiles are skipped, so
# a rerun only does missing work. Needs the
# GDAL command-line tools (gdalwarp, gdalbuildvrt, gdal_translate) on PATH;
# the national raster was built with GDAL 3.11.
#
#   list [dept]            frozen tile list (Colombia, or one department)
#   run <task_id> <block>  warp one block of tiles
#   fetch <task_id> <block>  download those tiles to $BLD_WORK_DIR/src (resumable)
#   gate [n_check] [block]  completeness + exact-native cross-check (block is
#                          only used to phrase the rerun hint, default 25)
#   mosaic                 zone-overlap check, national COG on the WorldPop grid
#
# Env: BLD_YEAR (2018), BLD_WORK_DIR (~/.cache/saecol2026/buildings_<year>),
#      BLD_TAG ("" ; forced to _d<dept> by `list <dept>`).
#
# BLD_SRC_DIR: when set, `run` warps from <BLD_SRC_DIR>/<tile_id>.tif instead of
# the bucket (for machines with no internet access). `fetch` writes exactly
# that layout, so BLD_SRC_DIR=$BLD_WORK_DIR/src pairs with it.
#
# The target grid is the constrained WorldPop R2024B 2018 raster, expected at
# geo_cache/worldpop_r2024b_2018.tif:
#   https://data.worldpop.org/GIS/Population/Global_2015_2030/R2024B/2018/COL/v1/100m/constrained/col_pop_2018_CN_100m_R2024B_v1.tif
# Public raster data; prints counts and summary statistics only.

suppressPackageStartupMessages({
  library(saecol2026); library(data.table); library(sf); library(terra)
})
args <- commandArgs(trailingOnly = TRUE)
mode <- if (length(args)) args[[1]] else stop("usage: build_buildings_raster.R <list|run|fetch|gate|mosaic> ...")
YEAR <- as.integer(Sys.getenv("BLD_YEAR", "2018"))
WORK <- Sys.getenv("BLD_WORK_DIR",
                   file.path(Sys.getenv("HOME"), ".cache", "saecol2026",
                             paste0("buildings_", YEAR)))
TAG  <- Sys.getenv("BLD_TAG", "")
if (mode == "list" && length(args) >= 2L) TAG <- paste0("_d", args[[2]])
LIST <- file.path(WORK, sprintf("bld_tiles_%d%s.rds", YEAR, TAG))
TILES <- file.path(WORK, "tiles")
SRCD <- Sys.getenv("BLD_SRC_DIR", "")
dir.create(TILES, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(WORK, "logs"), showWarnings = FALSE)
bld_gdal_env()
sf::sf_use_s2(FALSE)
say <- function(...) { cat(format(Sys.time(), "%H:%M:%S"), " ", ..., "\n", sep = ""); flush.console() }

colombia <- function(dept = NULL) {
  if (is.null(dept)) {
    g <- geodata::gadm("COL", level = 0, path = output_path("geo_cache"))
    return(sf::st_as_sf(g))
  }
  m <- sf::st_as_sf(colmaps::municipios)
  m[substr(sprintf("%05d", as.integer(as.character(m$id))), 1, 2) == dept, ]
}

if (mode == "list") {
  dept <- if (length(args) >= 2L) args[[2]] else NULL
  b <- colombia(dept)
  idx <- bld_tile_index(YEAR, c(-81.84, -4.23, -66.84, 15.91))
  sel <- bld_select_tiles(idx, b, buffer_m = 2000)
  data.table::setorderv(sel, c("epsg", "tile_id"))
  if (is.null(dept) && abs(nrow(sel) - 9750L) > 0.01 * 9750L)
    stop("tile count ", nrow(sel), " differs from the measured 9,750 by more than 1%")
  saveRDS(sel, LIST)
  say(sprintf("%d tiles (%s) -> %s", nrow(sel),
              paste(names(table(sel$epsg)), table(sel$epsg), collapse = ", "), LIST))
  say(sprintf("blocks of 25 tiles: %d", ceiling(nrow(sel) / 25)))

} else if (mode == "run") {
  task  <- as.integer(args[[2]]); block <- as.integer(args[[3]])
  tl    <- readRDS(LIST)
  rows  <- bld_task_rows(nrow(tl), task, block)
  gw    <- bld_gdal_bin("gdalwarp")
  say("gdalwarp: ", paste(system2(gw, "--version", stdout = TRUE), collapse = " "))
  res   <- list()
  for (i in rows) {
    t   <- tl[i]
    out <- file.path(TILES, paste0(t$tile_id, ".tif"))
    done <- file.exists(out) &&
      isTRUE(tryCatch(terra::nlyr(terra::rast(out)) == 4L, error = function(e) FALSE))
    src <- bld_tile_src(t, SRCD)
    st <- if (done) {
      data.frame(tile_id = t$tile_id, status = "skipped", secs = 0,
                 message = "", stringsAsFactors = FALSE)
    } else if (nzchar(SRCD) && !file.exists(src)) {
      data.frame(tile_id = t$tile_id, status = "failed", secs = 0,
                 message = paste("staged tile missing:", src),
                 stringsAsFactors = FALSE)
    } else {
      bld_warp_tile(t, out, src = src, gdalwarp = gw)
    }
    res[[length(res) + 1L]] <- st
    say(sprintf("task %d  %-28s %-8s %6.1f s %s", task, st$tile_id, st$status,
                st$secs, st$message))
  }
  log <- if (length(res)) do.call(rbind, res) else data.frame()
  data.table::fwrite(log, file.path(WORK, "logs", sprintf("task_%04d%s.csv", task, TAG)))
  if (any(log$status == "failed")) quit(save = "no", status = 1L)

} else if (mode == "fetch") {
  task  <- as.integer(args[[2]]); block <- as.integer(args[[3]])
  dest  <- if (nzchar(SRCD)) SRCD else file.path(WORK, "src")
  dir.create(dest, showWarnings = FALSE, recursive = TRUE)
  tl    <- readRDS(LIST)
  rows  <- bld_task_rows(nrow(tl), task, block)
  curl  <- Sys.which("curl"); if (!nzchar(curl)) stop("curl not found on PATH")
  res   <- list()
  for (i in rows) {
    t <- tl[i]
    f <- file.path(dest, paste0(t$tile_id, ".tif"))
    # a tile is "done" when it opens with 3 bands at the index's dimensions --
    # a truncated download fails that, so a resumed fetch re-pulls it
    ok <- file.exists(f) && isTRUE(tryCatch({
      r <- terra::rast(f)
      terra::nlyr(r) == 3L && terra::ncol(r) == t$nx && terra::nrow(r) == t$ny
    }, error = function(e) FALSE))
    if (ok) { res[[length(res) + 1L]] <- data.frame(tile_id = t$tile_id,
              status = "skipped", secs = 0, bytes = file.size(f)); next }
    t0  <- Sys.time()
    tmp <- paste0(f, ".part")
    st  <- system2(curl, c("-sS", "-f", "--retry", "5", "--retry-delay", "2",
                           "-o", shQuote(tmp),
                           shQuote(sub("^/vsicurl/", "", bld_tile_url(t$uri)))))
    good <- st == 0L && file.exists(tmp) && isTRUE(tryCatch({
      r <- terra::rast(tmp)
      terra::nlyr(r) == 3L && terra::ncol(r) == t$nx && terra::nrow(r) == t$ny
    }, error = function(e) FALSE))
    if (good) file.rename(tmp, f) else unlink(tmp)
    res[[length(res) + 1L]] <- data.frame(tile_id = t$tile_id,
      status = if (good) "ok" else "failed",
      secs = as.numeric(difftime(Sys.time(), t0, units = "secs")),
      bytes = if (good) file.size(f) else 0)
  }
  log <- if (length(res)) do.call(rbind, res) else data.frame()
  data.table::fwrite(log, file.path(WORK, "logs", sprintf("fetch_%04d%s.csv", task, TAG)))
  say(sprintf("fetch task %d: ok %d | skipped %d | failed %d | %.1f GB",
              task, sum(log$status == "ok"), sum(log$status == "skipped"),
              sum(log$status == "failed"), sum(log$bytes) / 1e9))
  if (any(log$status == "failed")) quit(save = "no", status = 1L)

} else if (mode == "gate") {
  n_check <- if (length(args) >= 2L) as.integer(args[[2]]) else 50L
  block   <- if (length(args) >= 3L) as.integer(args[[3]]) else 25L
  tl  <- readRDS(LIST)
  f   <- file.path(TILES, paste0(tl$tile_id, ".tif"))
  ok  <- vapply(f, function(p) file.exists(p) &&
                  isTRUE(tryCatch(terra::nlyr(terra::rast(p)) == 4L,
                                  error = function(e) FALSE)), logical(1))
  miss <- tl$tile_id[!ok]
  miss_file <- file.path(WORK, sprintf("missing_tiles%s.txt", TAG))
  writeLines(miss, miss_file)
  say(sprintf("tiles present and readable: %d of %d (missing %d -> missing_tiles%s.txt)",
              sum(ok), nrow(tl), length(miss), TAG))
  if (length(miss)) {
    # Completeness is a gate criterion: an incomplete tile set must not pass.
    tasks <- unique(ceiling(match(miss, tl$tile_id) / block))
    say("blocks to rerun: ", paste(tasks, collapse = ","))
    stop(length(miss), " tile(s) missing -- see ", miss_file)
  }
  # -99 leak: with the nodata flag set, a leak reads as a negative value
  mins <- vapply(f[ok], function(p)
    min(terra::minmax(terra::rast(p), compute = TRUE)[1, ]), numeric(1))
  say(sprintf("minimum value over all tiles: %.4g (must be >= 0)", min(mins)))
  if (any(mins < 0)) stop("negative values in ", sum(mins < 0), " tiles: -99 leak")
  # exact-native cross-check on random populated cells
  set.seed(20260917)
  rows <- sample(which(ok), min(n_check, sum(ok)))
  rel <- rbindlist(lapply(rows, function(i) {
    r <- terra::rast(f[i]); names(r) <- bld_band_names
    cand <- which(terra::values(r[["presence"]], mat = FALSE) > 0)
    if (!length(cand)) return(NULL)
    cell <- cand[sample.int(length(cand), 1L)]
    tmpl <- terra::rast(r[[1]]); terra::values(tmpl) <- NA; tmpl[cell] <- 1
    poly <- sf::st_transform(sf::st_as_sf(terra::as.polygons(tmpl)), tl$epsg[i])
    nat  <- terra::rast(bld_tile_url(tl$uri[i]))
    names(nat) <- c("count", "h", "p")
    ex <- exactextractr::exact_extract(nat, poly, function(df) {
      cv <- df$coverage_fraction
      ok <- !is.na(df$count) & !is.na(df$p) & !is.na(df$h)
      data.frame(count = sum(df$count[ok] * cv[ok]) / sum(cv[ok]),
                 presence = sum(df$p[ok] * cv[ok]) / sum(cv[ok]),
                 p_h = sum(df$p[ok] * df$h[ok] * cv[ok]) / sum(cv[ok]),
                 p_h2 = sum(df$p[ok] * df$h[ok]^2 * cv[ok]) / sum(cv[ok]))
    }, summarize_df = TRUE, progress = FALSE)
    w  <- terra::values(r)[cell, ]
    exv <- unlist(ex[1, bld_band_names])
    dv  <- abs(exv - w)
    data.table(band = bld_band_names, d_abs = dv, d_rel = dv / pmax(abs(exv), 1e-12))
  }))
  if (nrow(rel)) {
    s <- rel[, .(n = .N, median_rel = median(d_rel), p95_rel = quantile(d_rel, 0.95),
                 max_rel = max(d_rel), median_abs = median(d_abs),
                 p95_abs = quantile(d_abs, 0.95), max_abs = max(d_abs)), by = band]
    print(as.data.frame(s), digits = 3)
    if (any(s$median_rel > 0.005)) stop("median relative difference above 0.5% on a band")
    # Absolute criterion: the median-relative check alone never looks at the
    # tail, so the largest absolute difference of each band must also stay
    # below its threshold (count 1e-5, presence 1e-3, p_h 1e-2, p_h2 1e-1).
    abs_thresh <- c(count = 1e-5, presence = 1e-3, p_h = 1e-2, p_h2 = 1e-1)
    bad <- s[s$max_abs > abs_thresh[as.character(s$band)], ]
    if (nrow(bad))
      stop("max absolute difference above threshold on band(s): ",
           paste(sprintf("%s (%.3g > %.3g)", bad$band, bad$max_abs,
                         abs_thresh[as.character(bad$band)]), collapse = ", "))
  }
  say("GATE PASSED")

} else if (mode == "mosaic") {
  tl  <- readRDS(LIST)
  f   <- file.path(TILES, paste0(tl$tile_id, ".tif"))
  if (!all(file.exists(f))) stop(sum(!file.exists(f)), " tiles missing: run `gate` first")
  gbv <- bld_gdal_bin("gdalbuildvrt"); gw <- bld_gdal_bin("gdalwarp")
  gt  <- bld_gdal_bin("gdal_translate")
  say("gdalwarp: ", paste(system2(gw, "--version", stdout = TRUE), collapse = " "))
  vrt_of <- function(files, path) {
    lst <- tempfile(fileext = ".txt"); writeLines(files, lst)
    st <- system2(gbv, c("-q", "-overwrite", "-input_file_list", shQuote(lst), shQuote(path)))
    if (st != 0) stop("gdalbuildvrt failed for ", path)
    path
  }
  zones <- split(f, tl$epsg)
  zv <- lapply(names(zones), function(z)
    vrt_of(zones[[z]], file.path(WORK, sprintf("zone_%s%s.vrt", z, TAG))))
  names(zv) <- names(zones)

  # Zone-overlap agreement. Both coordinates are drawn from `e`, the true
  # extent intersection of the two zone VRTs (the real 17/18 overlap is lon
  # [-77.9492, -77.8858], the 18/19 one [-71.9858, -71.8258]), not from a
  # nominal meridian, which would sample few or no overlap cells.
  pairs <- list(c("32617", "32618"), c("32618", "32619"), c("32718", "32719"),
                c("32618", "32718"), c("32619", "32719"))
  set.seed(1)
  # Gate on the MEDIAN ABSOLUTE difference per band, not relative -- these
  # four bands are already on natural, bounded scales. A zone pair with at
  # least 20 sampled cells fails when a band's median exceeds its threshold
  # (count 1e-4, presence 0.01, p_h 0.1, p_h2 0.1).
  overlap_thresh <- c(count = 1e-4, presence = 0.01, p_h = 0.1, p_h2 = 0.1)
  agree <- rbindlist(lapply(pairs, function(p) {
    if (!all(p %in% names(zv)))
      return(data.table(pair = paste(p, collapse = "/"), band = bld_band_names, n = 0L))
    a <- terra::rast(zv[[p[1]]]); b <- terra::rast(zv[[p[2]]])
    names(a) <- bld_band_names; names(b) <- bld_band_names
    e <- terra::intersect(terra::ext(a), terra::ext(b))
    if (is.null(e))
      return(data.table(pair = paste(p, collapse = "/"), band = bld_band_names, n = 0L))
    n  <- 2000
    xy <- cbind(runif(n, e$xmin, e$xmax), runif(n, e$ymin, e$ymax))
    va <- terra::extract(a, xy)[, bld_band_names]; vb <- terra::extract(b, xy)[, bld_band_names]
    both <- stats::complete.cases(va) & stats::complete.cases(vb) & va[, "presence"] > 0
    if (!any(both))
      return(data.table(pair = paste(p, collapse = "/"), band = bld_band_names, n = 0L))
    dd <- abs(as.matrix(va[both, , drop = FALSE]) - as.matrix(vb[both, , drop = FALSE]))
    data.table(pair = paste(p, collapse = "/"), band = bld_band_names, n = sum(both),
               median_abs = apply(dd, 2, median), p95_abs = apply(dd, 2, quantile, 0.95))
  }), fill = TRUE)
  say("zone-overlap agreement (cells with presence > 0, both coords sampled from the true overlap):")
  print(as.data.frame(agree), digits = 3)
  bad <- agree[!is.na(agree$median_abs) & agree$n >= 20 &
                 agree$median_abs > overlap_thresh[as.character(agree$band)], ]
  if (nrow(bad))
    stop("zone-overlap median absolute difference above threshold: ",
         paste(sprintf("%s %s=%.4g (> %.4g, n=%d)", bad$pair, bad$band, bad$median_abs,
                       overlap_thresh[as.character(bad$band)], bad$n), collapse = "; "))

  nat <- vrt_of(f, file.path(WORK, sprintf("national%s.vrt", TAG)))
  # Band descriptions. Per-tile outputs carry
  # Description = count/presence/p_h/p_h2, but gdalbuildvrt (GDAL 3.11.5)
  # drops them: the VRT's <VRTRasterBand> elements come back with no
  # <Description>. Patch them back into the VRT XML before
  # warping, so a reader of the final raster does not have to rely on band
  # order alone.
  vrt_xml <- readLines(nat)
  for (i in seq_along(bld_band_names)) {
    tag <- sprintf('<VRTRasterBand dataType="Float32" band="%d">', i)
    hit <- grep(tag, vrt_xml, fixed = TRUE)
    if (length(hit) == 1L)
      vrt_xml[hit] <- paste0(vrt_xml[hit],
        sprintf("\n    <Description>%s</Description>", bld_band_names[i]))
  }
  writeLines(vrt_xml, nat)

  wp_f <- file.path(output_path("geo_cache"), "worldpop_r2024b_2018.tif")
  if (!file.exists(wp_f)) stop("missing the target grid ", wp_f, " (URL in this script's header)")
  wp  <- terra::rast(wp_f)
  e   <- as.vector(terra::ext(wp))
  res <- 1 / 1200
  # Grid-lattice assertion. `compareGeom(r, wp)` on the FINAL raster below
  # always holds because
  # -te/-ts are literally read from `wp` -- that only proves gdalwarp obeyed
  # its own arguments, never that `-r near` was a valid, 1:1 resampling in
  # the first place. `-r near` is a correct cell-for-cell map only when the
  # SOURCE (nat) and TARGET (wp) grids sit on the same lattice, so check that
  # directly on the pre-warp national VRT before spending the warp on it.
  # unname() matters: terra::ext()$xmin returns a NAMED scalar ("xmin"), and
  # naming it again below (c(x = ...)) would silently combine into "x.xmin"
  # rather than "x", making on_lattice["x"] return NA on lookup even though
  # the underlying value (and therefore the stop-condition check, which does
  # not index by name) is fine.
  frac_cells <- function(v, r) { g <- (unname(v) / r) %% 1; ifelse(g > 0.5, g - 1, g) }
  nat_e <- terra::ext(terra::rast(nat))
  on_lattice  <- c(x = frac_cells(nat_e$xmin, res), y = frac_cells(nat_e$ymax, res))
  grid_offset <- c(x = frac_cells(nat_e$xmin - e[1], res), y = frac_cells(nat_e$ymax - e[4], res))
  say(sprintf(paste("grid lattice check: national VRT off-lattice by x=%.2e y=%.2e cells |",
                    "offset vs WorldPop grid x=%.2e y=%.2e cells"),
              on_lattice["x"], on_lattice["y"], grid_offset["x"], grid_offset["y"]))
  if (any(abs(on_lattice) > 1e-6))
    stop("national VRT origin is not on the 1/1200 deg lattice: ",
         paste(sprintf("%.3e", on_lattice), collapse = ", "), " cells off")
  if (any(abs(grid_offset) > 1e-3))
    stop("national VRT grid is offset from the WorldPop grid by more than 1e-3 cell: ",
         paste(sprintf("%.3e", grid_offset), collapse = ", "))

  out  <- file.path(output_path("geo_cache"), sprintf("buildings_%d_3as%s.tif", YEAR, TAG))
  part <- paste0(out, ".part.tif"); raw <- paste0(out, ".raw.tif")
  unlink(c(part, raw))
  # Two-step warp -> COG, not a direct `-of COG` warp: GDAL 3.11.5's COG
  # driver drops band <Description> metadata on CreateCopy from a warp, while
  # a plain GTiff warp keeps them and a subsequent gdal_translate to COG then
  # carries them through.
  st <- system2(gw, c("-q", "-overwrite", "-te", sprintf("%.12f", e[c(1, 3, 2, 4)]),
                      "-ts", terra::ncol(wp), terra::nrow(wp), "-r", "near",
                      "-dstnodata", "-99", "-co", "COMPRESS=DEFLATE",
                      "-co", "PREDICTOR=3", "-co", "TILED=YES", "-co", "BIGTIFF=YES",
                      shQuote(nat), shQuote(raw)))
  if (st != 0 || !file.exists(raw)) stop("national warp failed")
  st2 <- system2(gt, c("-q", "-of", "COG", "-co", "COMPRESS=DEFLATE",
                       "-co", "PREDICTOR=3", "-co", "BIGTIFF=YES",
                       shQuote(raw), shQuote(part)))
  unlink(raw)
  if (st2 != 0 || !file.exists(part)) stop("COG translation failed")
  file.rename(part, out)
  r <- terra::rast(out)
  if (!identical(names(r), bld_band_names))
    stop("national raster band names came back as ", paste(names(r), collapse = ", "),
         " instead of ", paste(bld_band_names, collapse = ", "),
         " -- band descriptions did not survive the VRT edit + warp/translate")
  stopifnot(terra::compareGeom(r, wp, stopOnError = FALSE))

  # Provenance: the exact 9,750 tile URIs exist only in
  # this work directory's frozen list. Re-deriving them from the live GCS
  # manifests can silently return a different set that the 1% count gate in
  # `list` mode would not notice, so freeze a copy next to the raster itself.
  list_copy <- file.path(output_path("geo_cache"), basename(LIST))
  file.copy(LIST, list_copy, overwrite = TRUE)
  say("froze tile list alongside the raster: ", list_copy)

  # population-weighted coverage: share of WorldPop people in cells with no data
  cov <- terra::global(c(wp, is.na(r[["presence"]]) * wp), "sum", na.rm = TRUE)
  say(sprintf("national raster: %s | people in NA building cells: %.4f%% of %s",
              out, 100 * cov$sum[2] / cov$sum[1], format(round(cov$sum[1]), big.mark = ",")))
  say("MOSAIC DONE")
} else stop("unknown mode: ", mode)
