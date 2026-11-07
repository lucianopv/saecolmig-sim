# Synthetic fixtures for geo_aux tests. Three adjacent unit-square municipalities
# in a plain lon/lat grid, plus helpers to build matching SpatRasters with known
# values. No real data, no network.

make_geo_shp <- function() {
  # 3 municipalities: squares [0,1]x[0,1], [1,2]x[0,1], [2,3]x[0,1]
  poly <- function(xmin) {
    sf::st_polygon(list(rbind(
      c(xmin, 0), c(xmin + 1, 0), c(xmin + 1, 1), c(xmin, 1), c(xmin, 0))))
  }
  sf::st_sf(
    municipality = c(5001L, 5002L, 5003L),
    geometry = sf::st_sfc(poly(0), poly(1), poly(2), crs = 4326)
  )
}

# A constant-per-municipality raster: cell value = base + 10 * column-block index,
# so muni 5001 -> base, 5002 -> base+10, 5003 -> base+20 under area-weighted mean.
make_geo_raster <- function(base = 1, res = 0.25) {
  r <- terra::rast(xmin = 0, xmax = 3, ymin = 0, ymax = 1,
                   resolution = res, crs = "EPSG:4326")
  xs <- terra::xFromCell(r, seq_len(terra::ncell(r)))
  terra::values(r) <- base + 10 * floor(pmin(xs, 2.999))
  r
}
