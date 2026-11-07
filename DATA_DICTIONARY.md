# Data dictionary

This document gives the layout of every file the replicate runner reads and
writes: path, grain (what one row is), key columns, column types and row
universe, and the script that writes each input. This repository holds code
only: no census, survey or covariate data ships with it. It ships the sampling
design tables and a few small text tables (column selections, XGB
hyperparameters), listed in section 1. Inputs are built from public sources by
the scripts in `data-raw/` (see `data-raw/README.md`) or supplied by the user in
the layouts below.

Every input path is relative to the data root, the folder named by the
environment variable `SAECOL_DATA_DIR` and resolved by `saecol2026::output_path()`.
Outputs go to a directory chosen at run time (`out_dir`).

Conventions used throughout:

- `municipality` (called `area_id` in estimate files) is the integer DIVIPOLA
  code: 2-digit department x 1000 + 3-digit municipality, so a code has four or
  five digits. A zero-padded text code such as `"05001"` is read as `5001`.
- The **universe** is the set of municipalities of the census household frame
  (1,122 for Colombia). Every estimate table covers the whole universe: a
  municipality drawn in a replicate is `in_sample`, any other is `out_of_sample`.
- Types are those of the parquet files the scripts write (`int32`, `double`,
  `string`); R reads them as integer, numeric and character.

## 1. Files at a glance

| File (under the data root) | One row is | Key | Rows | Written by |
|---|---|---|---|---|
| `frame/national_household_frame.parquet` | a census household | `hh_id` | one per household | `data-raw/build_frame.R` |
| `aux_tables/muni_stratum_nn_extended.parquet` | a municipality | `municipality` | universe | `data-raw/build_nn_strata.R` |
| `aux_tables/truth_full.parquet` | a municipality and an outcome | `area_id`, `outcome` | universe x 6 | `data-raw/compute_truth_full.R` |
| `aux_tables/municipal_X_<area_config>.parquet` | a municipality | `municipality` | universe | `data-raw/build_area_geo_bld.R assemble` (configuration `geo_bld`) |
| `subarea/X_<subarea_config>.parquet` | a census block or rural section | `code` | one per sub-area | `data-raw/build_subarea_geo_blk_bld.R assemble` (configuration `geo_blk_bld`) |
| `aux_tables/W_queens_municipal.parquet` (optional) | a municipality (a matrix row) | `municipality` | universe | `saecol2026::build_queens_W()` |

`sim_inputs()` (`R/sim_inputs.R`) reads these files together with the shipped
tables below, validates the matrices and returns one list for the runner.
`<area_config>` and `<subarea_config>` default to `geo_bld` and `geo_blk_bld`.

Files shipped with the package, in `inst/extdata/`:

| File | Content | Read by |
|---|---|---|
| `ecv2024_design.yml` | The ECV 2024 sampling design: municipalities (first stage) selected with probability proportional to size without replacement within strata, then systematic selection of segments of about ten dwellings, every dwelling in a selected segment enumerated; the calibration margins. | `parse_ecv_design()`, called by `sim_inputs()` and `data-raw/tune_xgb_sub.R` |
| `ecv2024_strata_realised.csv` | One row per stratum of the ECV 2024 (357): `stratum`, `n_psu` (municipalities sampled), `n_seg` (segments sampled), `n_dwell` (dwellings sampled), as realised in the survey. | `parse_ecv_design()`, through the `strata` entry of the YAML |
| `ecv2024_muni_to_stratum.csv` | The 636 municipalities the ECV 2024 sampled and their stratum: `municipality`, `stratum`. | `data-raw/build_nn_strata.R` |
| `kept_geo_cols.txt`, `kept_bld_cols.txt` | One column name per line: the geo and buildings columns that enter the covariate matrices (section 2.4). | `kept_cols()`, used by the builders |
| `xgb_sub_params_<subarea_config>.tsv` | Tuned XGB hyperparameters, one row per outcome (section 2.7). | `xgb_sub_params()`, called by `sim_inputs()` |

## 2. Inputs

### 2.1 `frame/national_household_frame.parquet` — household frame

One row per household enumerated in the CNPV 2018 census, before any sampling.
It carries the migration answers and no weights, strata or sampling-unit ids;
those are attached when a replicate is drawn.

| Column | Type | Definition |
|---|---|---|
| `dept` | int32 | Department code (2 digits). |
| `muni` | int32 | Municipality code within the department (3 digits). |
| `municipality` | int32 | DIVIPOLA code, `dept * 1000 + muni`. |
| `area` | string | `"urban"` (cabecera) or `"rural"` (all other), from the census area class. |
| `seccion` | string | First 14 characters of the census block code: department, municipality, class, rural sector, rural section and centro poblado. |
| `manzana` | string | Census block code (22 characters). |
| `hh_id` | string | Household id built from the census dwelling and household numbers. |
| `n_persons` | int32 | Household size; missing for a small share of households. |
| `migrated_last_12m`, `migrated_internal_12m`, `migrated_intl_12m` | int32 (0/1) | Household indicators: any member moved in the last 12 months (to any destination, within Colombia, abroad). They overlap; never sum them. |
| `n_pers_resp` | int32 | Household members with a recorded migration answer. |
| `n_mig_total`, `n_mig_internal`, `n_mig_intl` | int32 | Members who moved (to any destination, within Colombia, abroad). |

A person has moved when the census answer to where they lived a year earlier is
another municipality of the country (internal) or another country
(international); a person who was not yet born, lived in the same municipality
or gave no answer has not moved, and the unknown answers are not counted in
`n_pers_resp`. `n_persons` is the household size of the household file; the
other household columns are built by joining the person records to it, so the
four person counts (`n_pers_resp`, `n_mig_total`, `n_mig_internal`,
`n_mig_intl`) and the three indicators are missing (NA) for a household that
has no person record.

### 2.2 `aux_tables/muni_stratum_nn_extended.parquet` — stratum map

The ECV 2024 stratum of each municipality of the universe. Municipalities the
public ECV design assigns to a stratum keep it; each of the others takes the
stratum of its nearest neighbour in the same department on a 12-feature
standardised census profile (`extend_strata_nn()`).

| Column | Type | Definition |
|---|---|---|
| `municipality` | int32 | DIVIPOLA code. |
| `dept` | int32 | Department code. |
| `stratum` | int32 | ECV stratum code. |
| `source` | string | `"known"` (assigned by the design) or `"nn"` (taken from the nearest neighbour). |

Households of a municipality that is absent from this map cannot be drawn; the
municipality stays in the universe and is always out of sample.

### 2.3 `aux_tables/truth_full.parquet` — census truth

The census value of each outcome for each municipality, computed from the frame
(`compute_truth_all()`): the mean of the household indicator, or the ratio of
the migrant count to `n_pers_resp` for the person outcomes.

| Column | Type | Definition |
|---|---|---|
| `area_id` | int32 | DIVIPOLA code. |
| `theta_true` | double | The census value. |
| `n_total` | int32 | Census units behind it: households, or the sum of `n_pers_resp` for a person outcome. |
| `outcome` | string | One of the six outcomes (section 3). |

The runner refuses a truth file that lacks any (municipality, outcome)
combination of the universe.

### 2.4 `aux_tables/municipal_X_geo_bld.parquet` — area-level covariate matrix

One row per municipality of the universe: the key, census composition and
geospatial layers. Columns of the matrix the builders write:

| Columns | Type | Definition |
|---|---|---|
| `municipality` | int32 | DIVIPOLA code (the key). |
| `dept`, `muni` | int32 | The two parts of the code. Not covariates; dropped when the matrix is read. |
| `n_persons`, `n_households` | int32 | Census population and household count. Not covariates; dropped when the matrix is read (section 2.8). |
| `share_urban` | double | Share of persons in the urban area class (cabecera; `UA_CLASE == 1`). |
| `share_working_age` | double | Share of persons aged 15 to 64 (age-group code `P_EDADR` from 4 to 13). |
| `share_elderly` | double | Share of persons aged 65 and over (`P_EDADR` of 14 or more). |
| `share_secondary_completed` | double | Mean over **all** persons of the indicator "age 25 or over and education-level code at least 7": `mean(P_NIVEL_ANOSR >= 7 & P_EDADR >= 6)`. The denominator is every person, not only those 25 and over. |
| `share_illiterate` | double | Mean over **all** persons of the indicator "age 15 or over and cannot read or write": `mean(P_ALFABETA == 2 & P_EDADR >= 4)`. The denominator is every person, not only those 15 and over. |
| `share_economically_active` | double | Share of persons who worked in the reference week (`P_TRABAJO == 1`). |
| `share_afro`, `share_indigenous` | double | Share of persons who identify as Afro-Colombian (`PA1_GRP_ETNIC == 5`), as indigenous (`PA1_GRP_ETNIC == 1`). |
| `share_piped_water` | double | Share of households that cook with public-aqueduct water (`H_AGUA_COCIN == 1`). |
| `share_electricity` | double | Share of dwellings with electricity service (`VA_EE == 1`). |
| geo layers, `kept_cols("geo")` | double | Nine layers over the municipality polygons, in this order: `ntl_mean_pw`, `travel_time_city_pw`, `terrain_ruggedness_pw`, `dist_venezuela_border_km`, `pop_density_2018`, `cropland_share_pw`, `forest_share_pw`, `precip_annual_2018_pw`, `precip_anomaly_2018_pw`. |
| buildings layers, `kept_cols("bld")` | double | `building_density_pw`, `building_height_pw`, `building_height_sd_pw`. |

All `share_*` columns lie in [0, 1]. Each is the mean of an indicator over the
records of the area (persons, households or dwellings), leaving out a record on
which the indicator is unknown. The age-group code `P_EDADR` counts five-year
groups from 1 for 0 to 4 years, so code 4 starts at 15, code 6 at 25, code 13
ends at 64 and code 14 starts at 65. For the two conjunctions, a person is
counted 0 when either condition is false, even if the other answer is missing
(a person aged under 25 is a 0 whatever the education answer); the indicator is
unknown, and the person left out of the mean, only when a missing answer
prevents deciding it (for example a missing education answer on a person aged 25
or over). The geo and buildings layers are
extracted over the municipality polygons from public rasters:

- A column ending in `_pw` is a population-weighted mean: the mean of the
  raster cells weighted by the fraction of each cell inside the polygon times
  the WorldPop population of the cell, which is the value the average resident
  experiences. The other layers are not population-weighted
  (`pop_density_2018` is the WorldPop 2018 population over the polygon area;
  `dist_venezuela_border_km` is the distance from the municipality centroid to
  the Colombia-Venezuela border, in a metric projection).
- `ntl_mean_pw`: night-time light radiance (Black Marble annual product, 2018).
  `travel_time_city_pw`: minutes of motorised travel to the nearest city.
  `terrain_ruggedness_pw`: terrain ruggedness index of an elevation model.
  `cropland_share_pw`, `forest_share_pw`: fractions of cropland and tree cover.
  `precip_annual_2018_pw`: 2018 CHIRPS precipitation, mm per year.
  `precip_anomaly_2018_pw`: that precipitation minus the 1970-2000 WorldClim
  normal.
- `building_density_pw`: buildings per km2. `building_height_pw`,
  `building_height_sd_pw`: presence-weighted mean and standard deviation of
  building height in metres (0 where no building is present). They come from the
  Open Buildings 2.5D Temporal layer.
- A municipality whose raster is missing (for example accessibility over the
  islands) takes the median over municipalities, so the matrix is complete.

Which geo and buildings columns enter is fixed by `inst/extdata/kept_geo_cols.txt`
and `inst/extdata/kept_bld_cols.txt`. The layers and their sources are listed in
`data-raw/README.md`.

### 2.5 `subarea/X_geo_blk_bld.parquet` — sub-area covariate matrix

A sub-area is a census block (urban areas) or a rural section (the block code
cut to its first 11 characters: department, municipality, class, rural sector
and rural section). The sub-area methods model
sub-areas and aggregate to municipalities. One row per sub-area:

| Column | Type | Definition |
|---|---|---|
| `code` | string | Sub-area code, the key. Always text: a numeric code would lose leading zeros. 22 characters for a block, 11 for a rural section. |
| `municipality` | int32 | DIVIPOLA code of the municipality it belongs to. |
| `stratum` | string | Area class, `"urban"` or `"rural"`. (It is not the ECV stratum; it matches sample households to sub-areas.) |
| `n_households_census`, `n_persons_census` | int32 | Census households and persons. They weight the roll-up of sub-area predictions to a municipality (households for household outcomes, persons for person outcomes). |
| `area_km2` | double | Polygon area. |
| `pop_worldpop` | double | WorldPop population of the sub-area; the alternative roll-up weight that needs no census count. |
| `share_working_age`, `share_elderly`, `share_secondary_completed`, `share_illiterate`, `share_economically_active`, `share_afro`, `share_indigenous`, `share_piped_water`, `share_electricity` | double | Same definitions as in section 2.4, over the records of the sub-area; a share over no record is missing. |
| `is_urban` | int32 | 1 for a block, 0 for a rural section. |
| geo layers, `kept_cols("geo", level = "subarea")` | double | The area-level layers extracted over the sub-area polygons (population-weighted on constrained WorldPop weights), named without the census year: `ntl_mean_pw`, `travel_time_city_pw`, `terrain_ruggedness_pw`, `dist_venezuela_border_km`, `pop_density`, `cropland_share_pw`, `forest_share_pw`, `precip_annual_pw`, `precip_anomaly_pw`. A sub-area whose extraction is missing takes its municipality's value of the same layer. |
| buildings layers | double | `building_density_pw`, `building_height_pw`, `building_height_sd_pw`, extracted the same way. |

The first seven columns are structural; every other column is a covariate.
The **prediction universe** is the set of sub-areas complete on every
covariate (`subarea_universe()`); the others are left out, and the runner
reports how many. Every municipality of the universe must keep at least one
sub-area.

### 2.6 `aux_tables/W_queens_municipal.parquet` — spatial weights (optional)

Row-standardised queens-contiguity weights over the universe, used only by the
`fh_spatial` method; the runner loads them when the file exists. A
`municipality` column (int32) followed by one `double` column per
municipality, named by its code, in the same order as the rows. Cell (i, j) is
the weight of municipality j as a neighbour of i; each non-empty row sums to 1,
and an island with no neighbours has an all-zero row.

### 2.7 `inst/extdata/xgb_sub_params_<subarea_config>.tsv` — XGB hyperparameters

Shipped with the package. Tab-separated, one row per outcome. Columns: `outcome`
(string), then one numeric column per hyperparameter: `nrounds`, `eta`,
`max_depth`, `min_child_weight`, `subsample`, `colsample_bytree`, `lambda`,
`alpha`. Read with `xgb_sub_params()`; `sim_inputs()` takes the table of
`subarea_config`. `data-raw/tune_xgb_sub.R` repeats the selection on another
covariate matrix and writes a table in this layout to
`aux_tables/xgb_sub_params_<CONFIG>.tsv`, which is read with
`xgb_sub_params(file = ...)`.

### 2.8 Supplying covariates

Two ways, both checked by the rules below before any fit.

1. **Extra columns joined on a key** (`extra_area`, `extra_subarea`). A parquet
   or CSV file, or a data frame, keyed by `municipality` (area level) or by
   `code` (sub-area level; a column named `subarea` is read as `code`) with one
   or more covariate columns. They are appended to the built matrix.
2. **A complete matrix in place of the built one** (`area_matrix`,
   `subarea_matrix`), in the layouts of sections 2.4 and 2.5. Its covariate
   columns may be any numeric columns; the structural columns of section 2.5 are
   required at sub-area level.

CSV files are read with the key columns as text, so leading zeros survive.

Rules, applied to the built, supplied and extended matrices alike. An error
names the file (when the input came from a path), the problem and the count.

| Rule | Behaviour |
|---|---|
| `municipality` is a whole-number code; text such as `"05001"` is accepted | Coerced to integer. A value that is missing or not a whole number stops the run. |
| One row per municipality of the universe (area level) | A missing, repeated or out-of-universe municipality stops the run. |
| Covariates are numeric and finite (both levels) | A column that is not numeric, or that holds an infinite value, stops the run, named. |
| Covariates have no missing value (area level) | A column with a missing value stops the run, named with its count. (At sub-area level a sub-area with a missing covariate leaves the prediction universe, below.) |
| `dept`, `muni`, `n_persons`, `n_households` are not covariates | Dropped, with a message. They are identifiers or sizes; a size that spans orders of magnitude against shares near one makes the cross-product matrix of the robust Fay-Herriot fit ill-conditioned. |
| A covariate column without any observed value | Dropped, with a message. |
| `code` is character, unique and non-empty (sub-area level) | Otherwise the run stops. |
| `stratum` is `"urban"` or `"rural"`; the two census counts are present and complete | Otherwise the run stops. |
| A sub-area missing any covariate | Left out of the prediction universe; the count is reported. |
| Every universe municipality keeps a sub-area; no sub-area lies outside the universe | Otherwise the run stops. |
| The sub-area codes of the frame's households are found in the matrix | The frame codes that are missing, and the share of households in them, are reported (sub-areas incomplete on a covariate are expected to be missing). If they hold more than 5% of the frame's households, or the matrix matches none, the run stops with the count. |
| Truth: no missing or non-finite `theta_true`, no repeated (`area_id`, `outcome`) pair, a row for every municipality of the universe and every outcome | Otherwise the run stops with the count. |
| Extra covariates: no name already in the matrix, no repeated key, a row for every key of the matrix | A name that collides, a repeated key or an unmatched key stops the run. Rows for keys outside the matrix are ignored, with a message. |

## 3. Outcomes

Six estimands, registered in `saecol2026::sae_outcomes`. The `outcome` column of
every file uses these names.

| Level | Overall | Within Colombia | Abroad |
|---|---|---|---|
| Household indicator (0/1, any member; overlapping) | `migrated_last_12m` | `migrated_internal_12m` | `migrated_intl_12m` |
| Person ratio (migrants over `n_pers_resp`; additive) | `pers_migrated_12m` | `pers_migrated_internal_12m` | `pers_migrated_intl_12m` |

Household indicators are survey-weighted means; person ratios are ratio
estimators of the migrant count over `n_pers_resp`.

## 4. Outputs

### 4.1 `<out_dir>/direct/rep_NNNN.parquet` — direct estimates of a replicate

One row per municipality drawn in replicate `NNNN` and outcome, for all six
outcomes whichever outcomes are fitted, computed from the calibrated sample (seed
`100000 + NNNN`) by `direct_estimate_all()` and `smooth_variance_all()`; written
by `sim_run_replicate()` (`R/sim_run.R`).

| Column | Type | Definition |
|---|---|---|
| `area_id` | int32 | DIVIPOLA code. |
| `theta_direct` | double | Design-based estimate: weighted mean (household outcomes) or weighted ratio (person outcomes) under the sample design, with the calibrated weight. |
| `var_direct` | double | Taylor-linearised sampling variance of `theta_direct`. |
| `n_eff` | double | Effective sample size, `(sum w)^2 / sum w^2` over the contributing households. |
| `n_dwell` | int32 | Sampled households in the municipality. |
| `outcome` | string | One of the six outcomes. |
| `var_direct_smoothed` | double | Variance smoothed by a generalised variance function, fitted per outcome; the sampling variance the area-level models take as known. Municipalities with a zero raw variance get the predicted value. |

### 4.2 `<out_dir>/fits/<method>__<config>/rep_NNNN.parquet` — per-fit results

One file per method, covariate configuration and replicate, written by
`sim_run_replicate()` in the long layout below, rows ordered by outcome as in
`sae_outcomes`: one row per municipality of the universe and outcome. A failed
fit is recorded, not omitted: one row per municipality drawn in the replicate
(`area_set = in_sample`) with `fit_failed = TRUE`, the message in `error_msg`
and missing estimates. The columns and their order are `sim_fit_columns`.

Both kinds of file are written under a temporary name (`rep_NNNN.parquet.tmp.<pid>`)
and renamed when complete, so a file that exists was written whole; temporary
files of the same target older than an hour, left by an interrupted process,
are removed when the file is written (younger ones may belong to a process
that is still writing). A rerun keeps a direct file that reads and holds the six
outcomes, and keeps the outcomes already in a per-fit file that reads, has every
column and holds rows of its replicate; it fits only the outcomes that are
missing and rewrites the file with the kept and the new rows. A replicate whose
files are all complete is skipped without drawing its sample.

A failed fit is kept by a rerun. A failure that comes from the method (a
singular fit, a bootstrap that does not converge) recurs under the same seed;
one that comes from the environment (a process killed for memory, a file system
error) does not. `retry_failed = TRUE` (`scripts/run_replicates.R --retry-failed`)
fits again the requested outcomes whose fit failed and rewrites the file; failed
rows of outcomes not requested stay in it. Missing packages,
CmdStan not found and Stan models that do not compile stop the run before the
first draw, so they never reach these files. Each replicate reports, per method,
the number of requested outcomes whose fit failed; `sim_run()` returns them in
its `failed_<method>` columns and the script prints their totals.

| Column | Type | Definition |
|---|---|---|
| `area_id` | int32 | DIVIPOLA code. |
| `theta_hat` | double | Point estimate. For `bayes_fh_beta` out of sample it is tipsae's `theta_oos`, mu, not the in-sample theta = mu (1 - p0 - p1) + p1; with exact zeros mu exceeds the proportion where zeros are common, above all for the international outcomes. The evaluation scores `theta_hat`. |
| `mse` | double | The method's own reported uncertainty, as a variance: emdi MSE for `fh_arcsin`, posterior variance for `bayes_fh_beta`, povmap bootstrap variance for the sub-area methods. |
| `ci_lo`, `ci_hi` | double | Bounds of the 95% interval. |
| `area_set` | string | `in_sample` or `out_of_sample` for this replicate. |
| `fit_failed` | logical | TRUE when the fit raised an error or returned no usable estimates. |
| `error_msg` | string | The error message (first 200 characters): the message a fit function reports in its own diagnostics when it has one (an error inside povmap, for the sub-area methods), otherwise the error raised by the fit or by the contract check; missing on success. |
| `config` | string | Covariate configuration: `geo_bld` for the area-level methods, `geo_blk_bld` for the sub-area methods. |
| `method` | string | Method code (`fh_arcsin`, `bayes_fh_beta`, `xgb_sub`, ...). |
| `replicate` | integer | Replicate number. |
| `outcome` | string | One of the six outcomes. |
| `seed` | integer | The fit's seed: `sae_fit_seed(config, method, outcome, replicate)`, or the seed `settings` give the method instead. |
| `wall_secs` | double | Wall-clock seconds of the fit. |
| `weight` | string | Roll-up weight of a sub-area method (`census` or `worldpop`); missing for area-level methods. |
| `n_divergent` | double | Divergent transitions of the Bayesian fit; missing for other methods and for failed fits. |
| `max_rhat` | double | Largest rank-normalised split R-hat (`posterior::rhat`) over the variables the fit keeps. `bayes_fh_beta` under ExtBeta: every variable of the model and its generated quantities on CmdStan; on rstan `beta0`, `beta`, `v`, `sigma_v`, `mu`, `p0`, `p1`, `lambda_EB`, `theta`, `theta_oos`, `y_rep`, `log_lik` and `lp__`. `bayes_fh_heavy`: `sigma_v` and `beta`. Missing for other methods and for failed fits. |
| `sigma2_v` | double | Estimated variance of the area random effect: `sigma2_v` of a Fay-Herriot fit, the squared random-effect standard deviation of `megb_sub`; missing for other methods and for failed fits. |
| `povmap_sha` | string | Commit of the povmap installation, for the sub-area methods. |
| `mse_type` | string | Which bootstrap uncertainty `megb_sub` reports: `var` (variance of the centred bootstrap prediction errors) or `mse` (with the bootstrap bias term); missing for other methods. |
| `pkg_versions` | string | Versions of the packages behind the fit (`sim_pkg_versions()`), e.g. `emdi2=2.2.4@4a50df5;povmap=1.0.1@2b3959a;xgboost=3.1.2.1;tipsae=1.0.3;cmdstan=2.38.0`: `@` and the first seven characters of the commit for a package installed from a repository commit, `NA` for one that is not installed. |

### 4.3 External per-fit file

Estimates of a method fitted elsewhere (the repository does not refit it) enter
the evaluation in this layout, a parquet or CSV file or a data frame, one row per
replicate, outcome and municipality of the universe:

| Column | Type | Definition |
|---|---|---|
| `replicate` | integer | Replicate number. |
| `outcome` | string | One of the six outcomes. |
| `municipality` | integer | DIVIPOLA code. |
| `sampled` | logical | TRUE when the municipality was drawn in that replicate. |
| `estimate` | double | Point estimate. |
| `estimated_variance` | double | Reported variance; not negative. |
| `lower`, `upper` | double | Interval bounds; `lower` is not above `upper`. |
| `truth` | double | Census value; must equal `theta_true` of `truth_full.parquet`. |

The evaluation reads it with `sim_read_external()`, which maps it to the
per-fit shape of section 4.2 (`area_id` from `municipality`, `area_set` from
`sampled`, `theta_hat` from `estimate`, `mse` from `estimated_variance`,
`ci_lo`/`ci_hi` from `lower`/`upper`, `fit_failed = FALSE`,
`config = "external"`, and the method label given to it) and refuses missing values, an
unknown outcome, a `sampled` that is not logical (or 0/1), repeated
(`outcome`, `replicate`, `municipality`) rows, a negative variance, a lower
bound above its upper bound and unknown municipality codes. `sim_evaluate()`
also refuses a truth that differs from the truth file, a replicate range that
differs from the run's, and a drawn set that differs from the run's.

### 4.4 Evaluation tables

`sim_evaluate()` reads the direct and per-fit files of a run folder and,
optionally, an external per-fit file, joins the truth and returns aggregates
only; `scripts/evaluate.R` writes each table to `<out>/<table>.parquet`. The
measures are defined in `R/uncertainty_eval.R` and `R/uncertainty_display.R`.
Unless stated, a table is keyed by `method`, `config`, `outcome` and
`area_set`, and every table carries `column`, the label of the covariate set
(`geo` for the configurations here). `method` holds the display label of
`sim_method_labels` (`FH-arcsin`, `Bayes-Beta`, `MEGB`, `XGB`; a method without
one keeps its code) or the label given to the external file; `config` holds the
configuration the method was fitted on (`geo_bld`, `geo_blk_bld`, or `external`
for the external file). Before any table is computed, a fit (method, outcome
and replicate) is dropped whole when one of its rows has `fit_failed` or its
`max_rhat` exceeds 1.05, and rows with a non-finite estimate, MSE or bound
are dropped. A row of the estimate table is one fit of one municipality in one
replicate; `emp_mse` is the realised MSE of a municipality's estimates over
replicates, `rep_mse` the mean of the method's reported `mse`.

| Table | Grain | Content |
|---|---|---|
| `area_summary` | one municipality, per key | `n_reps`, `theta_true`; `emp_mse` (mean squared error against the truth), `emp_var` (variance of the estimates around their mean), `bias`, `bias_share` (`bias^2 / emp_mse`); `rep_mse`, `rep_rmse` (mean of the root reported MSE); `ratio` (`rep_mse / emp_mse`), `ratio_var` (`rep_mse / emp_var`), `rrmse_mse` (RMSE of the reported MSE around `emp_mse`, relative to it), `rb_rmse` (`rep_rmse / sqrt(emp_mse) - 1`); `coverage` and `width_mean` of the 95% interval, `width_rel` (`width_mean` over the width of a normal 95% interval at the realised RMSE, `2 z sqrt(emp_mse)` with `z = qnorm(0.975)`), `iscore` (interval score); `n_dwell_med`, the median dwellings drawn over the municipality's replicates in the set (0 out of sample). The ratios are missing where `emp_mse` (or `emp_var`) is 0. |
| `set_summary` | per key | Pooled over rows: `n_rows`, `coverage`, `miss_lo`, `miss_hi` (the truth is below, above the interval), `width_med`, `iscore`, `ratio_pooled` (`mean(mse) / mean(squared error)`), `rmse`, `rmse_reported`. Over municipalities with at least 20 replicates: `n_areas`, `n_ratio_undefined` (those whose `ratio` is missing because `emp_mse` is 0), quartiles of `ratio` (`ratio_med`, `ratio_q25`, `ratio_q75`) and shares below 0.5 and above 2, medians of `rrmse_mse`, `ratio_var`, `bias_share`, `width_rel`, `rb_rmse` (and its mean), the rank correlation `spearman_rmse` of `sqrt(rep_mse)` (the root of the mean reported MSE) with the realised RMSE `sqrt(emp_mse)`, and the shares of municipalities with coverage under 0.90 and under 0.80. |
| `import_check` | one row | `source` (file name of the external per-fit file, `<data frame>` for one given in memory), `rows` (rows of the external method that enter the tables) and `truth_max_abs_diff` (largest absolute difference between the external file's truth and the truth file), all three missing without an external file; over every evaluated method, `rows_kept` (rows that enter the tables), `fits_dropped` (fits dropped whole), of which `fits_dropped_failed` (a row with `fit_failed`) and `fits_dropped_rhat` (split R-hat above 1.05 and no failed row), and `rows_nonfinite` (rows dropped for a non-finite estimate, MSE or bound). |
| `fit_check` | per method, config and outcome | What happened to the fits (one per replicate) of each method and outcome: `fits_read`; `failed` (dropped whole for a row with `fit_failed`); `dropped_rhat` (dropped whole for a split R-hat above 1.05, no failed row); `dropped_nonfinite` (left without a row once rows with a non-finite estimate, MSE or bound are dropped); `kept` (fits with rows in the tables, `fits_read - failed - dropped_rhat - dropped_nonfinite`); `rows_nonfinite` (rows dropped for a non-finite value in the fits not dropped whole). |
| `by_bin` | per key and bin | `bin_var` names the binning: `size`, the dwellings drawn in the municipality and replicate, in eight quantile bins of the dwellings drawn over the rows of the direct files (fewer where quantiles coincide), each labelled in `bin` by its range, with the out-of-sample rows as a first bin of their own, `0 (out of sample)`; `prevalence`, the quintile of the municipality's true value within its outcome (`bin` 1 to at most 5: quintile limits that coincide give fewer bins). `bin_order` orders the bins of a `bin_var` (1 is the out-of-sample bin of `size`). Pooled over the rows of the bin: `n`, `coverage`, `width_med`, `ratio_pooled`, `iscore`, `miss_lo`, `miss_hi`, `rmse_empirical`, `rmse_estimated`, `truth_mean`; `x_med`, the median dwellings drawn of a `size` bin (missing for `prevalence`). |
| `cv_flags` | per key | Whether an estimate is publishable at a coefficient-of-variation threshold (0.20): `n` rows; `share_pub_reported` (share of rows whose reported CV, `sqrt(mse) / theta_hat`, is below the threshold) and `share_pub_true` (share whose municipality has a realised CV, `sqrt(emp_mse) / theta_true`, below it); `false_publish` (share of rows published on the reported CV whose realised CV is not below the threshold) and `missed_publish` (share of rows publishable on the realised CV that the reported CV does not publish). A non-positive estimate or truth is never publishable. |
| `mcse` | per key | Monte Carlo standard errors from a bootstrap over whole replicates (200 draws): `mcse_coverage`, `mcse_ratio_pooled`, `mcse_iscore`. |
| `reliability` | per key and one decile of the reported RMSE (`rbin`) | `n`, `rmse_reported`, `rmse_realised` within the decile. |
| `ratio_summary` | per method, outcome and area set | The ratio of estimated to empirical RMSE, `rep_rmse / sqrt(emp_mse)`, per municipality (those with `emp_mse` of 0 left out): quartiles `p25`, `med`, `p75` over municipalities, and `below`, the share of municipalities with a ratio under 1. |
| `tracking_summary` | per method, outcome and area set | `spearman`, the rank correlation of `rep_rmse` with `sqrt(emp_mse)` across municipalities, and `below`, the share with `rep_rmse` under `sqrt(emp_mse)`. |
| `size_medians` | per method, outcome, kind and sample-size group | In-sample municipalities grouped by their median dwellings drawn (groups of at least 25 municipalities, formed on the sizes of the first method in the order of `sim_method_labels`): `n_dwell` (median size) and `rmse` (median) for `kind` "Empirical RMSE" (`sqrt(emp_mse)`) and "Estimated RMSE" (`rep_rmse`). |
| `coverage_summary` | per method and outcome | In-sample municipalities: `overall` (coverage weighted by replicates), `share_ge_095` (share of municipalities with coverage of at least 0.95) and `share_lt_080` (share under 0.80). |

The last four tables are computed from `area_summary` for every outcome.

The evaluation stops, with a message that says which, when the external truth
differs from the truth file by more than 1e-9, the external replicates differ
from those requested, a method has no file for a requested replicate or no fit
of an evaluated outcome in a replicate, a fit that did not fail lacks one row
per municipality of the universe, a fit has no drawn municipality, the drawn
municipalities differ across methods, a method keeps no fit of an outcome
(naming how many failed, exceeded the R-hat threshold or had no finite row), or
a drawn municipality has no direct estimate. It warns, once, naming the methods
and their counts, when the methods keep different numbers of fits of an
outcome: they are then compared over different replicates.

### 4.5 `<out_dir>/run_manifest.json` — what the results of a folder are computed from

Written by the runner the first time it writes into `out_dir` and compared on
every later run into the same folder; a run whose manifest differs stops before
its first draw and names the differing fields, so results computed from
different inputs never mix in one folder.

| Field | Content | Compared |
|---|---|---|
| `area_config`, `subarea_config` | Configuration labels. | yes |
| `settings` | `sim_settings()` with the XGB hyperparameters: every argument the fits receive besides data and seed. | yes |
| `data.area_matrix_md5`, `data.subarea_matrix_md5` | MD5 of the area-level and sub-area matrices as the fits read them (after validation and extra covariates). | yes |
| `sources.<input>` | `path` as given, `size` in bytes and `md5` of the file, for `frame`, `stratum_map`, `area_matrix`, `subarea_matrix`, `extra_area`, `extra_subarea`, `spatial_weights` and `xgb_params` (`null` when not used; `"<data frame>"` for a matrix given in memory). | size and MD5 |
| `pkg_versions.<package>` | Versions as in the `pkg_versions` column. | yes, unless either is `NA` |
