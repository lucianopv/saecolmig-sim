# Published covariate selections, tuned XGB hyperparameters and the sub-area
# prediction universe.

# Column names build_geo_aux() produces, from the three-municipality fixture.
geo_aux_names <- function() {
  shp <- make_geo_shp()
  border <- sf::st_sf(geometry = sf::st_sfc(
    sf::st_linestring(rbind(c(3, -1), c(3, 2))), crs = 4326))
  mk <- function(v) { r <- make_geo_raster(); terra::values(r) <- v; r }
  layers <- list(ntl = mk(5), travel_time = mk(120), dem = mk(800),
                 builtup = mk(0.1), cropland = mk(0.3), forest = mk(0.4),
                 precip_2018 = mk(1200), precip_normal = mk(1000),
                 pop_2018 = mk(50), pop_2015 = mk(40))
  names(build_geo_aux(shp, layers, border))
}

test_that("the published selections name covariates the geo and building builders produce", {
  geo <- kept_cols("geo")
  bld <- kept_cols("bld")
  expect_gt(length(geo), 0L)
  expect_gt(length(kept_cols("bld")), 0L)
  expect_false(anyDuplicated(c(geo, bld)) > 0L)
  expect_true(all(bld %in% bld_cov_names))
  expect_true(all(geo %in% setdiff(geo_aux_names(), "municipality")))
})

test_that("kept_cols returns the published geo and building selections", {
  expect_identical(kept_cols("geo"), c(
    "ntl_mean_pw", "travel_time_city_pw", "terrain_ruggedness_pw",
    "dist_venezuela_border_km", "pop_density_2018", "cropland_share_pw",
    "forest_share_pw", "precip_annual_2018_pw", "precip_anomaly_2018_pw"))
  expect_identical(kept_cols("bld"), c(
    "building_density_pw", "building_height_pw", "building_height_sd_pw"))
  expect_identical(kept_cols(), kept_cols("geo"))
})

test_that("at sub-area level the geo selection carries the layer names without the year", {
  expect_identical(kept_cols("geo", level = "subarea"), c(
    "ntl_mean_pw", "travel_time_city_pw", "terrain_ruggedness_pw",
    "dist_venezuela_border_km", "pop_density", "cropland_share_pw",
    "forest_share_pw", "precip_annual_pw", "precip_anomaly_pw"))
  expect_identical(kept_cols("bld", level = "subarea"), kept_cols("bld"))
})

test_that("geo_subarea_names maps every year-labelled geo column and leaves the rest", {
  nm <- geo_aux_names()
  out <- geo_subarea_names(nm)
  expect_length(out, length(nm))
  expect_false(any(grepl("2018|2015", out)))
  expect_false(anyDuplicated(out) > 0L)
  expect_identical(out[!grepl("2018|2015", nm)], nm[!grepl("2018|2015", nm)])
  expect_identical(geo_subarea_names(c("pop_change_2015_2018", "precip_anomaly_2018")),
                   c("pop_change", "precip_anomaly"))
})

test_that("an unknown selection or level is refused, listing the choices", {
  # match.arg() lists the allowed values; the rest of its message is translated.
  expect_error(kept_cols("census"), "bld")
  expect_error(kept_cols("geo", level = "block"), "subarea")
})

test_that("XGB hyperparameters cover the six outcomes", {
  p <- xgb_sub_params("geo_blk_bld")
  expect_setequal(names(p), sae_outcomes$outcome)
  expect_true(all(vapply(p, function(x) all(c("eta", "max_depth", "nrounds") %in% names(x)), TRUE)))
})

test_that("XGB hyperparameters carry the tuned values with integer rounds and depth", {
  p <- xgb_sub_params()
  expect_identical(names(p), sae_outcomes$outcome)
  expect_identical(vapply(p, function(x) x$eta, numeric(1)), c(
    migrated_last_12m = 0.05, migrated_internal_12m = 0.05,
    migrated_intl_12m = 0.05, pers_migrated_12m = 0.10,
    pers_migrated_internal_12m = 0.10, pers_migrated_intl_12m = 0.02))
  expect_identical(vapply(p, function(x) x$max_depth, integer(1)), c(
    migrated_last_12m = 2L, migrated_internal_12m = 3L, migrated_intl_12m = 2L,
    pers_migrated_12m = 2L, pers_migrated_internal_12m = 2L,
    pers_migrated_intl_12m = 2L))
  for (x in p) {
    expect_identical(names(x), c("nrounds", "eta", "max_depth", "min_child_weight",
                                 "subsample", "colsample_bytree", "lambda", "alpha"))
    expect_identical(x$nrounds, 200L)
    expect_identical(x[c("min_child_weight", "subsample", "colsample_bytree",
                         "lambda", "alpha")],
                     list(min_child_weight = 1, subsample = 0.8,
                          colsample_bytree = 0.8, lambda = 1, alpha = 0))
  }
})

test_that("a configuration without a tuned table is refused", {
  expect_error(xgb_sub_params("geo"),
               "no tuned XGB hyperparameters ship for configuration 'geo'")
})

test_that("a hyperparameter table round-trips through its file", {
  p <- xgb_sub_params()
  f <- withr::local_tempfile(fileext = ".tsv")
  expect_identical(write_xgb_sub_params(p, f), f)
  expect_identical(xgb_sub_params(file = f), structure(p, source = f))
})

test_that("a hyperparameter table missing an outcome is refused, naming it", {
  p <- xgb_sub_params()
  f <- withr::local_tempfile(fileext = ".tsv")
  write_xgb_sub_params(p[-3], f)
  expect_error(xgb_sub_params(file = f), "migrated_intl_12m")
})

test_that("a hyperparameter table with a missing value is refused, naming the column", {
  p <- xgb_sub_params()
  p$migrated_last_12m$eta <- NA_real_
  f <- withr::local_tempfile(fileext = ".tsv")
  write_xgb_sub_params(p, f)
  expect_error(xgb_sub_params(file = f), "missing values in eta")
})

test_that("writing an empty or unnamed hyperparameter list is refused", {
  f <- withr::local_tempfile(fileext = ".tsv")
  expect_error(write_xgb_sub_params(list(), f), "non-empty list named by outcome")
  expect_error(write_xgb_sub_params(unname(xgb_sub_params()), f),
               "non-empty list named by outcome")
  expect_false(file.exists(f))
})

test_that("the sub-area universe is the complete cases of the covariate columns", {
  pop <- data.table::data.table(
    code = c("a", "b", "c", "d"), municipality = 1:4, stratum = "urban",
    n_households_census = 1:4, n_persons_census = 1:4,
    area_km2 = c(1, NA, 1, 1), pop_worldpop = c(1, NA, 1, NA),
    x1 = c(1, 2, NA, 4), x2 = c(1, 2, 3, NA))
  u <- subarea_universe(pop)
  # NA in a structural column (area, WorldPop weight) does not exclude a row.
  expect_identical(u$code, c("a", "b"))
  expect_true(data.table::is.data.table(u))
  expect_identical(names(u), names(pop))
  v <- subarea_universe(as.data.frame(pop))
  expect_false(data.table::is.data.table(v))
  expect_identical(v$code, c("a", "b"))
})

test_that("the covariate columns of a sub-area matrix are all but the structural ones", {
  pop <- data.frame(code = "a", x1 = 1, municipality = 1L, stratum = "urban",
                    n_households_census = 1L, n_persons_census = 1L,
                    area_km2 = 1, pop_worldpop = 1, x2 = 2)
  expect_identical(subarea_covariate_names(pop), c("x1", "x2"))
})

test_that("a sub-area matrix without covariate columns is refused", {
  pop <- data.frame(code = "a", municipality = 1L, stratum = "urban",
                    n_households_census = 1L, n_persons_census = 1L)
  expect_error(subarea_universe(pop), "covariate")
})
