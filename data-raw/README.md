# data-raw/

These scripts build the input files of the simulation from the raw CNPV 2018
census, the DANE MGN 2018 geographic framework and public rasters. Every
output is written under `SAECOL_DATA_DIR` through `saecol2026::output_path()`;
no script embeds a machine path.

## Environment

| variable | used by | meaning |
|---|---|---|
| `SAECOL_DATA_DIR` | all | data root (see `.Renviron.example`); raw sources and outputs live under it (layout below) |
| `MGN_DIR` | `build_subarea_geo.R subarea`, `build_subarea_geo_blk_bld.R bld` | directory of the unzipped MGN 2018 GeoPackages `MGN_URB_MANZANA2018geopackage.gpkg`, `MGN_RUR_SECCION2018geopackage.gpkg`, `MGN_URB_SECCION2018geopackage.gpkg` (default: `census_path()`) |
| `DEPTS` | the same two passes | comma-separated department codes: write those departments' checkpoints only, without the national table (smoke runs, or splitting the work); a run without `DEPTS` then assembles all checkpoints |
| `BLD_YEAR` | buildings scripts | Open Buildings epoch (2018) |
| `BLD_TAG` | buildings scripts | suffix of the raster and of every buildings output (e.g. `_d63` for a one-department run) |
| `BLD_WORK_DIR`, `BLD_SRC_DIR` | `build_buildings_raster.R` | tile work directory; local tile source for machines without internet |
| `CONFIG`, `XFILE`, `TUNE_REPS` | `tune_xgb_sub.R` | output name, sub-area matrix to tune on, number of tuning replicates (5) |

## Raw-source layout

Under `SAECOL_DATA_DIR`:

- `Census2018/` — the unzipped CNPV 2018 department files (read through
  `census_path()` by `build_frame.R`, `build_auxiliary.R`,
  `build_subarea_census.R`).
- `geo_cache/` — public rasters, one file each, clipped to Colombia or
  national (the sources are listed in the top-level `README.md`):

  | file | layer |
  |---|---|
  | `ntl_2018.tif` | night-time lights, NASA Black Marble VNP46A4 (annual), 2018, from NASA LAADS DAAC: `blackmarbler::bm_raster(roi_sf = <Colombia boundary>, product_id = "VNP46A4", date = 2018, bearer = <Earthdata token>)` (the use of `blackmarbler` in Suggests) |
  | `travel_time_2015.tif` | travel time to cities, Malaria Atlas Project (Weiss et al. 2018): `malariaAtlas::getRaster(dataset_id = "Accessibility__201501_Global_Travel_Time_to_Cities", shp = <Colombia boundary>)` (the use of `malariaAtlas` in Suggests) |
  | `elevation.tif` | SRTM elevation, `geodata::elevation_30s(country = "COL")` |
  | `landcover_built.tif`, `landcover_cropland.tif`, `landcover_trees.tif` | ESA WorldCover fractions, `geodata::landcover()` |
  | `chirps_2018.tif` | CHIRPS v2.0 annual precipitation, 2018 (`chirps-v2.0.2018.tif`) |
  | `precip_normal.tif` | WorldClim 2.1 precipitation, the 12 monthly layers of `geodata::worldclim_global(var = "prec", res = 2.5)` summed |
  | `worldpop_2015.tif`, `worldpop_2018.tif` | unconstrained WorldPop, `https://data.worldpop.org/GIS/Population/Global_2000_2020/<year>/COL/col_ppp_<year>.tif` |
  | `worldpop_r2024b_2015.tif`, `worldpop_r2024b_2018.tif` | constrained WorldPop R2024B, `https://data.worldpop.org/GIS/Population/Global_2015_2030/R2024B/<year>/COL/v1/100m/constrained/col_pop_<year>_CN_100m_R2024B_v1.tif` |
  | `buildings_2018_3as.tif` | written by `build_buildings_raster.R` (below) |

  GADM boundaries are fetched into the same directory by `geodata::gadm()`.
- The MGN 2018 GeoPackages, unzipped under `MGN_DIR`: DANE Geoportal,
  <https://geoportal.dane.gov.co>, Marco Geoestadístico Nacional 2018, block
  (manzana) and section (sección) levels (`MGN_URB_MANZANA2018geopackage.zip`, `MGN_RUR_SECCION2018geopackage.zip`,
  `MGN_URB_SECCION2018geopackage.zip`).

To keep raw sources elsewhere, symlink `Census2018/` and `geo_cache/` into
`SAECOL_DATA_DIR`.

## Build chain

Run in this order; each script reads outputs of earlier ones. Paths are
relative to `SAECOL_DATA_DIR`.

Frame, truth and strata:

| script | reads | writes |
|---|---|---|
| `build_frame.R` | `Census2018/` | `frame/national_household_frame.parquet` |
| `compute_truth.R`, `compute_truth_full.R` | the frame | `aux_tables/truth.parquet`, `aux_tables/truth_full.parquet` |
| `build_auxiliary.R` | `Census2018/` | `aux_tables/municipal_X.parquet` (census covariates) |
| `build_nn_strata.R` | `municipal_X.parquet` | `aux_tables/muni_stratum_nn_extended.parquet` |

Area level, configuration `geo_bld`:

| script | reads | writes |
|---|---|---|
| `build_municipal_geo.R` | `geo_cache/` | `aux_tables/municipal_geo.parquet` |
| `assemble_geo_X.R` | `municipal_X`, `municipal_geo` | `aux_tables/municipal_X_geo.parquet` |
| `build_buildings_raster.R` | Open Buildings tiles, `worldpop_r2024b_2018.tif` (target grid) | `geo_cache/buildings_2018_3as.tif` |
| `build_area_geo_bld.R municipal` | the buildings raster, `worldpop_2018.tif`, `worldpop_r2024b_2018.tif` | `aux_tables/municipal_bld_2018.parquet`, `aux_tables/municipal_bld_2018con.parquet` |
| `build_area_geo_bld.R assemble` | `municipal_X_geo`, `municipal_bld_2018` | `aux_tables/municipal_X_geo_bld.parquet` |

Sub-area level, configuration `geo_blk_bld`:

| script | reads | writes |
|---|---|---|
| `build_subarea_pop.R` | the frame | `aux_tables/subarea_pop.parquet` |
| `build_subarea_census.R` | `Census2018/`, `subarea_pop` | `aux_tables/subarea_census.parquet` |
| `build_subarea_geo.R subarea` | MGN, `geo_cache/`, `subarea_pop` | `aux_tables/subarea_geo.parquet` (checkpoints in `aux_tables/subarea_geo_parts/`) |
| `build_subarea_geo.R municipal` | `geo_cache/` | `aux_tables/municipal_geo_2018con.parquet` |
| `build_subarea_geo_blk_bld.R bld` | MGN, the buildings raster, `worldpop_r2024b_2018.tif`, `subarea_pop`, `municipal_bld_2018con` | `aux_tables/subarea_bld.parquet` (checkpoints in `aux_tables/subarea_bld_parts_2018/`) |
| `build_subarea_geo_blk_bld.R assemble` | `subarea_census`, `subarea_geo`, `municipal_geo_2018con`, `subarea_bld` | `subarea/X_geo_blk_bld.parquet` |
| `tune_xgb_sub.R` (optional) | `X_geo_blk_bld` (or `XFILE`), the frame, the strata | `aux_tables/xgb_sub_params_<CONFIG>.tsv` |

`W_queens_municipal` (the queens-contiguity weights of `fh_spatial`) has no
script here: `saecol2026::build_queens_W()` builds it on first use.

The raster passes are long: the national buildings raster is built in blocks
of tiles that can run separately (below); `build_subarea_geo.R subarea` and
`build_subarea_geo_blk_bld.R bld` read hundreds of thousands of polygons,
department by department, and resume from their checkpoints.
`tune_xgb_sub.R` runs `6 x TUNE_REPS` cross-validated grid searches of 32
configurations each.

Besides the packages in `DESCRIPTION` Imports, the builders use `colmaps`
(municipal polygons), `geodata` (GADM boundaries) and `lwgeom` (ellipsoidal
areas of the MGN polygons), all in Suggests. Fetching the rasters of
`geo_cache/` uses `geodata`, `blackmarbler` (night-time lights) and
`malariaAtlas` (travel time), also in Suggests; no script here downloads them.

## Published selections and hyperparameters

The geo and buildings columns of both configurations are fixed by
`inst/extdata/kept_geo_cols.txt` and `inst/extdata/kept_bld_cols.txt`
(`kept_cols()`), the selections of the published fits. Each came from a VIF
prune (threshold 10, `bld_prune_vif()`) against the other covariates of its
parent matrix; a prune against other covariates can keep other columns, so the
builders apply the lists. `assemble_geo_X.R` and `build_area_geo_bld.R vif`
print the prune on this repository's matrices beside the lists. At sub-area
level the geo columns carry the layer names without the census year
(`kept_cols("geo", level = "subarea")`).

`inst/extdata/xgb_sub_params_geo_blk_bld.tsv` (`xgb_sub_params()`) holds the
tuned XGB hyperparameters, one row per outcome, selected on the published
sub-area matrix. `tune_xgb_sub.R` repeats the selection on another matrix;
pass its table to the runner with `scripts/run_replicates.R --xgb-params=PATH`
(in R, `xgb_sub_params(file = ...)`).

## The sub-area matrix

`subarea/X_geo_blk_bld.parquet` has one row per census block (urban) or rural
section and these columns, in order:

- structural: `code`, `municipality`, `stratum`, `n_households_census`,
  `n_persons_census`, `area_km2`, `pop_worldpop` (never covariates;
  `pop_worldpop = pop_density * area_km2` is the WorldPop population of the
  sub-area);
- census: the nine `share_*` columns and `is_urban`;
- geo: `kept_cols("geo", level = "subarea")`, extracted over the sub-area
  polygons on constrained WorldPop R2024B weights; a sub-area whose extraction
  is missing takes its municipality's value of the same layer on the same
  raster (`municipal_geo_2018con`);
- buildings: `kept_cols("bld")`, extracted the same way, missing values taken
  from `municipal_bld_2018con`.

It holds the sub-areas complete on every covariate (`subarea_universe()`),
which is the prediction universe of the sub-area methods.

## Building covariates

Open Buildings 2.5D Temporal (epoch of 30 June 2018, public bucket
`open-buildings-temporal-data`, no credentials) becomes the buildings block in
three steps. Functions: `R/buildings.R`, plus `extract_weighted_mean_multi()`
in `R/geo_aux.R`.

```sh
# 1. National raster: 9,750 tiles warped onto the WorldPop 3-arcsec grid.
#    `run` warps one block of 25 tiles (about 1 h); 390 blocks in all, which
#    can run separately or in parallel. Finished tiles are skipped, so a rerun
#    only does missing work.
Rscript data-raw/build_buildings_raster.R list
Rscript data-raw/build_buildings_raster.R run <task_id> 25    # task_id 1..390
Rscript data-raw/build_buildings_raster.R gate                # completeness + native cross-check
Rscript data-raw/build_buildings_raster.R mosaic              # -> geo_cache/buildings_2018_3as.tif
# 2. Municipal tables (both WorldPop weightings) and the geo_bld matrix
Rscript data-raw/build_area_geo_bld.R municipal
Rscript data-raw/build_area_geo_bld.R assemble
# 3. Sub-area table (after build_subarea_geo.R), then the geo_blk_bld matrix
Rscript data-raw/build_subarea_geo_blk_bld.R bld
Rscript data-raw/build_subarea_geo_blk_bld.R assemble
```

- Step 1 needs the GDAL command-line tools on `PATH` (built with GDAL 3.11)
  and, as its target grid, `geo_cache/worldpop_r2024b_2018.tif` (URL in the
  script's header). A machine without internet access during the warp can stage
  the tiles first with the `fetch` mode and set `BLD_SRC_DIR`.
- Step 2 weights by population: `municipal_bld_2018` uses the unconstrained
  `worldpop_2018.tif`, as the area-level geo block does, and is the one the
  `geo_bld` matrix takes; `municipal_bld_2018con` uses the constrained raster
  and supplies the municipal value for sub-areas without one.

Without the raw census and raster sources, the runner can take covariate
matrices built elsewhere, in the layouts described in `DATA_DICTIONARY.md`.
