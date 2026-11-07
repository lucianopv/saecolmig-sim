#!/usr/bin/env Rscript
# data-raw/build_subarea_geo_blk_bld.R
# The sub-area covariate matrix of configuration `geo_blk_bld`: per census block
# or rural section, the nine census shares, an urban flag, the geo block at
# sub-area resolution and the buildings block.
#
#   bld       population-weighted building covariates per sub-area
#             -> aux_tables/subarea_bld_parts_<year>/dept_<dd>.parquet (checkpoints)
#             -> aux_tables/subarea_bld.parquet
#   assemble  -> subarea/X_geo_blk_bld.parquet
#
#   Rscript data-raw/build_subarea_geo_blk_bld.R <bld|assemble>
#
# Inputs under SAECOL_DATA_DIR, from the other builders:
#   aux_tables/subarea_pop.parquet            build_subarea_pop.R
#   aux_tables/subarea_census.parquet         build_subarea_census.R
#   aux_tables/subarea_geo.parquet            build_subarea_geo.R subarea
#   aux_tables/municipal_geo_2018con.parquet  build_subarea_geo.R municipal
#   aux_tables/municipal_bld_2018con.parquet  build_area_geo_bld.R municipal
#   geo_cache/buildings_<year>_3as.tif        build_buildings_raster.R
#   geo_cache/worldpop_r2024b_2018.tif        constrained WorldPop R2024B
# Env: MGN_DIR (default: census_path()), DEPTS (department subset, checkpoints
#      only), BLD_YEAR (2018), BLD_TAG (suffix of the raster and outputs, as in
#      build_buildings_raster.R).
#
# Prints counts, quantiles and correlations only.

suppressPackageStartupMessages({
  library(saecol2026); library(sf); library(terra); library(arrow); library(data.table)
})
args <- commandArgs(trailingOnly = TRUE)
mode <- if (length(args)) args[[1]] else stop("usage: build_subarea_geo_blk_bld.R <bld|assemble>")
YEAR <- as.integer(Sys.getenv("BLD_YEAR", "2018"))
TAG  <- Sys.getenv("BLD_TAG", "")
GEO  <- output_path("geo_cache")
AUX  <- output_path("aux_tables", create = TRUE)
say  <- function(...) { cat(format(Sys.time(), "%H:%M:%S"), " ", ..., "\n", sep = ""); flush.console() }
rd_aux <- function(f) {
  p <- file.path(AUX, f)
  if (!file.exists(p)) stop("missing input: ", p)
  as.data.table(read_parquet(p))
}

if (mode == "bld") {
  # The recipe of the municipal building covariates (build_area_geo_bld.R):
  # a WorldPop-weighted zonal mean of the four building bands, then the
  # presence-weighted height and its SD, over the sub-area polygons. Weights are
  # the constrained R2024B raster, as for the sub-area geo block.
  sf::sf_use_s2(FALSE)
  MGN_DIR <- Sys.getenv("MGN_DIR", census_path())
  gp <- function(f) {
    p <- file.path(MGN_DIR, f)
    if (!file.exists(p)) stop("missing MGN GeoPackage: ", p, "\nset MGN_DIR")
    p
  }
  F_MZ <- gp("MGN_URB_MANZANA2018geopackage.gpkg")
  F_RS <- gp("MGN_RUR_SECCION2018geopackage.gpkg")
  F_US <- gp("MGN_URB_SECCION2018geopackage.gpkg")
  pop <- as.data.table(read_parquet(file.path(AUX, "subarea_pop.parquet"),
                                    col_select = c("code", "municipality")))
  pop[, cls := substr(code, 6, 6)]

  bp <- file.path(GEO, sprintf("buildings_%d_3as%s.tif", YEAR, TAG))
  if (!file.exists(bp)) stop("missing building raster: ", bp)
  CK <- file.path(AUX, sprintf("subarea_bld_parts_%d%s", YEAR, TAG))
  dir.create(CK, showWarnings = FALSE)
  # Every checkpoint of one directory must come from one building raster: its
  # identity is recorded on first use and checked on every later run.
  stamp_f   <- file.path(CK, "source.rds")
  stamp_now <- list(file = basename(bp), size = file.size(bp),
                    mtime = format(file.info(bp)$mtime, tz = "UTC"))
  if (file.exists(stamp_f)) {
    if (!identical(readRDS(stamp_f), stamp_now))
      stop("checkpoints in ", CK, " come from a different building raster; ",
           "delete the directory and run again")
  } else {
    saveRDS(stamp_now, stamp_f)
  }

  b <- terra::rast(bp); names(b) <- bld_band_names
  w <- terra::rast(file.path(GEO, "worldpop_r2024b_2018.tif"))
  depts <- sort(unique(substr(pop$code, 1, 2)))
  SUBSET <- nzchar(Sys.getenv("DEPTS"))
  if (SUBSET) depts <- intersect(depts, strsplit(Sys.getenv("DEPTS"), ",")[[1]])
  dup_total <- 0L
  for (dd in depts) {
    f_ck <- file.path(CK, sprintf("dept_%s.parquet", dd))
    if (file.exists(f_ck)) { say("dept ", dd, " cached"); next }
    t0 <- Sys.time()
    g <- subarea_geometry(dd, pop, F_MZ, F_RS, F_US)
    if (is.null(g)) { say("dept ", dd, ": no polygons"); next }
    g <- subarea_union_codes(g)
    nd <- attr(g, "n_duplicate")
    dup_total <- dup_total + nd
    g$municipality <- g$code     # bld_pw_covariates() keys on `municipality`
    # Building bands and weights share one grid, so both are cropped to the
    # department and the crops line up.
    e  <- terra::ext(as.numeric(sf::st_bbox(sf::st_transform(g, 4326)))[c(1, 3, 2, 4)])
    d  <- as.data.table(bld_pw_covariates(g[, "municipality"],
                                          terra::crop(b, e, snap = "out"),
                                          terra::crop(w, e, snap = "out")))
    setnames(d, "municipality", "code")
    write_parquet_safe(d, f_ck)
    say(sprintf("dept %s: %6d sub-areas (%d duplicate-code rows merged), %.1f min",
                dd, nrow(d), nd,
                as.numeric(difftime(Sys.time(), t0, units = "mins"))))
  }
  say(sprintf("duplicate-code rows merged: %d", dup_total))
  if (SUBSET) {
    say("DEPTS subset: checkpoints written, national table skipped.")
    quit(save = "no")
  }
  parts <- rbindlist(lapply(list.files(CK, "\\.parquet$", full.names = TRUE),
                            function(f) as.data.table(read_parquet(f))))
  stopifnot(!anyDuplicated(parts$code))
  X <- merge(pop[, .(code, municipality)], parts, by = "code", all.x = TRUE)

  # A sub-area without a block-level value (no polygon, or a band with no data
  # over it) takes its municipality's value on the same constrained weights.
  # The four bands go missing independently, so a row is flagged if any of them
  # was missing.
  mcon <- rd_aux(sprintf("municipal_bld_%dcon%s.parquet", YEAR, TAG))
  X[, bld_backfilled := !stats::complete.cases(X[, ..bld_cov_names])]
  for (k in bld_cov_names) {
    mv <- mcon[[k]][match(X$municipality, mcon$municipality)]
    i  <- which(is.na(X[[k]]) & !is.na(mv))
    if (length(i)) set(X, i, k, mv[i])
  }
  say(sprintf("sub-areas: %s | filled from municipality: %s (%.3f%%) | still NA: %d",
              format(nrow(X), big.mark = ","), format(sum(X$bld_backfilled), big.mark = ","),
              100 * mean(X$bld_backfilled), sum(!stats::complete.cases(X[, ..bld_cov_names]))))
  say(sprintf("sub-areas with no buildings where people live: %s (%.2f%%)",
              format(sum(X$building_presence_pw == 0, na.rm = TRUE), big.mark = ","),
              100 * mean(X$building_presence_pw == 0, na.rm = TRUE)))
  q <- rbindlist(lapply(bld_cov_names, function(k) data.table(
    covariate = k, n_na = sum(is.na(X[[k]])),
    p05 = stats::quantile(X[[k]], .05, na.rm = TRUE),
    median = stats::median(X[[k]], na.rm = TRUE),
    p95 = stats::quantile(X[[k]], .95, na.rm = TRUE))))
  print(as.data.frame(q), digits = 4)

  # Agreement with the municipal layer: sub-areas aggregated to municipality
  # with their WorldPop weight (pop_density * area_km2 of the sub-area geo
  # block), ranked against the municipal values.
  p_geo <- file.path(AUX, "subarea_geo.parquet")
  if (file.exists(p_geo)) {
    sg <- as.data.table(read_parquet(p_geo, col_select = c("code", "pop_density", "area_km2")))
    Y  <- merge(X, sg[, .(code, wp = pop_density * area_km2)], by = "code", all.x = TRUE)
    agg <- Y[!is.na(wp) & wp > 0, lapply(.SD, function(v) sum(v * wp) / sum(wp)),
             by = municipality, .SDcols = bld_cov_names]
    cmp <- merge(agg, mcon, by = "municipality", suffixes = c("_agg", "_muni"))
    say("sub-areas aggregated to municipality vs the municipal layer:")
    print(data.frame(covariate = bld_cov_names,
                     spearman = vapply(bld_cov_names, function(k)
                       stats::cor(cmp[[paste0(k, "_agg")]], cmp[[paste0(k, "_muni")]],
                                  method = "spearman", use = "complete.obs"), 0),
                     pct_munis_varying = vapply(bld_cov_names, function(k)
                       100 * Y[, .(v = uniqueN(get(k)) > 1L), by = municipality][, mean(v)], 0)),
          digits = 3, row.names = FALSE)
  }
  stopifnot(nrow(X) == nrow(pop), !anyNA(X[, ..bld_cov_names]))
  # Ranges: presence is a share, density and the height moments are not negative.
  stopifnot(all(X$building_presence_pw >= 0 & X$building_presence_pw <= 1),
            all(X$building_density_pw >= 0), all(X$building_height_pw >= 0),
            all(X$building_height_sd_pw >= 0))
  out <- file.path(AUX, paste0("subarea_bld", TAG, ".parquet"))
  write_parquet_safe(X, out)
  say("wrote ", out)

} else if (mode == "assemble") {
  keys     <- c("code", "municipality", "stratum", "n_households_census", "n_persons_census")
  kept_geo <- kept_cols("geo", level = "subarea")
  kept_bld <- kept_cols("bld", level = "subarea")

  cen <- rd_aux("subarea_census.parquet")
  geo <- rd_aux("subarea_geo.parquet")
  mg  <- rd_aux("municipal_geo_2018con.parquet")
  bld <- rd_aux("subarea_bld.parquet")
  shares <- grep("^share_", names(cen), value = TRUE)
  stopifnot(length(shares) == 9L, all(kept_geo %in% names(geo)), all(kept_geo %in% names(mg)),
            all(kept_bld %in% names(bld)), "area_km2" %in% names(geo))

  # The sources must describe the same sub-areas in the same order, or a
  # column-wise assembly would pair the wrong rows.
  cen[, municipality := as.integer(municipality)]
  geo[, municipality := as.integer(municipality)]
  bld[, municipality := as.integer(municipality)]
  mg[,  municipality := as.integer(municipality)]
  setkey(cen, code); setkey(geo, code); setkey(bld, code)
  stopifnot(identical(cen$code, geo$code), identical(cen$code, bld$code),
            identical(cen$municipality, bld$municipality))
  for (k in keys) stopifnot(identical(cen[[k]], geo[[k]]))

  # Geo: a sub-area whose block-level extraction is missing takes its
  # municipality's value of the same layer on the same constrained raster. The
  # gaps are counted and reported, so the share of sub-areas that carry a
  # municipal value is known.
  touched <- rowSums(is.na(geo[, ..kept_geo])) > 0
  fills <- integer(0)
  for (k in kept_geo) {
    mv <- mg[[k]][match(geo$municipality, mg$municipality)]
    i  <- which(is.na(geo[[k]]) & !is.na(mv))
    fills[k] <- length(i)
    if (length(i)) set(geo, i, k, mv[i])
  }
  cat("--- block-level geo gaps, filled from municipality ----------------\n")
  for (k in names(fills)[fills > 0]) {
    cat(sprintf("  %-26s %s of %s sub-areas (%.3f%%)\n", k, format(fills[[k]], big.mark = ","),
                format(nrow(geo), big.mark = ","), 100 * fills[[k]] / nrow(geo)))
  }
  cat(sprintf("  sub-areas touched: %s (%.3f%%), holding %.3f%% of households\n",
              format(sum(touched), big.mark = ","), 100 * mean(touched),
              100 * sum(geo$n_households_census[touched]) / sum(geo$n_households_census)))

  # WorldPop aggregation weight: density times polygon area recovers the zonal
  # population sum. Structural, never a covariate.
  geo[, pop_worldpop := pop_density * area_km2]

  X <- cbind(geo[, c(keys, "area_km2", "pop_worldpop"), with = FALSE],
             cen[, ..shares],
             data.table(is_urban = as.integer(cen$stratum == "urban")),
             geo[, ..kept_geo],
             bld[, ..kept_bld])
  xcols <- c(shares, "is_urban", kept_geo, kept_bld)
  stopifnot(!anyDuplicated(names(X)), identical(setdiff(names(X), c(keys, "area_km2", "pop_worldpop")), xcols))

  # Prediction universe: the sub-areas complete on every covariate.
  Xu <- subarea_universe(X)
  cat("\n--- coverage ------------------------------------------------------\n")
  cat(sprintf("sub-areas %s | complete on all %d covariates %s (%.2f%%)\n",
              format(nrow(X), big.mark = ","), length(xcols),
              format(nrow(Xu), big.mark = ","), 100 * nrow(Xu) / nrow(X)))
  cat(sprintf("households covered %.2f%% | municipalities %d of %d | urban / rural %s / %s\n",
              100 * sum(Xu$n_households_census) / sum(X$n_households_census),
              uniqueN(Xu$municipality), uniqueN(X$municipality),
              format(sum(Xu$stratum == "urban"), big.mark = ","),
              format(sum(Xu$stratum == "rural"), big.mark = ",")))
  cat(sprintf("sub-areas with a WorldPop weight: %.2f%%\n",
              100 * mean(!is.na(Xu$pop_worldpop) & Xu$pop_worldpop > 0)))
  # Municipal estimates aggregate over the predicted sub-areas, so every
  # municipality must keep at least one.
  stopifnot(uniqueN(Xu$municipality) == uniqueN(X$municipality))
  M <- cbind(1, as.matrix(Xu[, ..xcols]))
  cat(sprintf("design matrix: cond %.3e | rcond(X'X) %.3e | rank %d of %d\n",
              kappa(M, exact = TRUE), rcond(crossprod(M)), qr(M)$rank, ncol(M)))

  out <- file.path(output_path("subarea", create = TRUE), "X_geo_blk_bld.parquet")
  write_parquet_safe(tibble::as_tibble(Xu), out)
  cat("\nwrote:", out, "\n")
  cat("covariates:", length(xcols), "( census", length(shares) + 1L, "+ geo", length(kept_geo),
      "+ buildings", length(kept_bld), ")\n")

} else stop("unknown mode: ", mode)
