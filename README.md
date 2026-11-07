# saecolmig-sim

Code for a design-based simulation that compares small area estimators of
migration for the 1,122 municipalities of Colombia. Samples are drawn from the
2018 census (CNPV 2018) under the design of the Encuesta Nacional de Calidad de
Vida (ECV) 2024; each sample is fitted with area-level and sub-area methods, and
the point estimates and the uncertainty each method reports are scored against
the census values.

The repository holds code only. It ships no census, survey or covariate data
and no results, and none are available from it on request. The only tables it
holds are the design tables and small text tables in `inst/extdata/` (the ECV
2024 strata and the sample realised in each, the covariate selections, the
tuned XGB hyperparameters). The inputs are built from public sources by the
scripts in `data-raw/`, or supplied by the user in the layouts of
`DATA_DICTIONARY.md`.

- `METHODS.md`: estimands, sampling design, direct estimator, covariates,
  methods and uncertainty measures.
- `DATA_DICTIONARY.md`: the layout of every file the runner reads and writes.
- `data-raw/README.md`: the build chain from raw sources to inputs.

## What the code does

1. **Replicates.** Replicate `r` draws a sample from the census household
   frame under the ECV 2024 design with seed `100000 + r`, calibrates its
   weights, computes design-based direct estimates of six migration outcomes
   for the drawn municipalities and smooths their sampling variances.
2. **Fits.** Each method is fitted to each outcome of each replicate with its
   own seed and predicts every municipality, drawn (in sample) or not (out of
   sample):
   - FH-arcsin: a Fay-Herriot model with the arcsine transformation
     (`fh_arcsin`, emdi2), on the municipal covariate matrix `geo_bld`;
   - Bayes-Beta: a Bayesian Fay-Herriot model with the extended Beta
     likelihood of tipsae (`bayes_fh_beta`, sampled by CmdStan), on
     `geo_bld`;
   - XGB: gradient-boosted predictions of census blocks and rural sections,
     aggregated to municipality (`xgb_sub`, povmap), on the sub-area matrix
     `geo_blk_bld`;
   - MEGB: mixed-effects gradient boosting at the same sub-area level
     (`megb_sub`, povmap). In the published comparison the MEGB estimates
     came from an external per-fit file (format below); this repository also
     ships its own MEGB wrapper.
3. **Evaluation.** Estimated against empirical RMSE, the relative bias of the
   estimated RMSE, interval coverage, width and interval score, the bias share
   of the empirical MSE, rank correlation, results by sample size, a check of
   publication at a coefficient of variation below 20%, and Monte Carlo
   standard errors, over replicates 1-100.

## Installation

R 4.4 or later. Clone the repository and restore the package library from the
lockfile, which pins every package the runner, the methods, the evaluation and
the tests load:

```r
install.packages("renv")
renv::restore()
```

The lockfile and `DESCRIPTION` pin emdi2 at commit `4a50df5`, povmap at commit
`2b3959a`, xgboost 3.1.2.1, tipsae 1.0.3 and cmdstanr at commit `da99e2b`. The
command-line scripts load the installed package, so install it into the
restored library with `R CMD INSTALL .` from the repository root (`make
install`); the tests run from the source tree.

Bayes-Beta samples with CmdStan. The code was verified with CmdStan 2.38.0:

```r
cmdstanr::install_cmdstan(version = "2.38.0")
```

Then copy `.Renviron.example` to `.Renviron` and set two variables:

- `SAECOL_DATA_DIR`: the folder that holds the inputs (see
  `DATA_DICTIONARY.md`).
- `CMDSTAN`: the CmdStan folder, for example `~/.cmdstan/cmdstan-2.38.0`. Set
  it explicitly. When it is unset, the runner takes the newest
  `cmdstan-<version>` folder under `~/.cmdstan` (cmdstanr's default install
  location) and stops if there is none; a CmdStan installed anywhere else is
  found only through `CMDSTAN`.

The building-covariate scripts also need the GDAL command-line tools on `PATH`
(the published building raster was built with GDAL 3.11).

## Data sources

All raw inputs are public. Place them under `SAECOL_DATA_DIR` as described in
`data-raw/README.md`, which gives the file names the scripts expect.

| Input | Source |
|---|---|
| CNPV 2018 (Censo Nacional de Población y Vivienda): person, household, dwelling and geographic-key (MGN) files, one ZIP file per department | DANE microdata catalog, <https://microdatos.dane.gov.co>. Unzip under `Census2018/`. |
| ECV 2024 | DANE microdata catalog 861, <https://microdatos.dane.gov.co/index.php/catalog/861>. Only its design is used, and that ships in `inst/extdata/` (the strata, the 636 sampled municipalities and the sample realised in each stratum), so the survey microdata are not needed to run the code. |
| Census blocks and sections: the `MGN_URB_MANZANA2018`, `MGN_RUR_SECCION2018` and `MGN_URB_SECCION2018` GeoPackages | DANE Geoportal, <https://geoportal.dane.gov.co>, Marco Geoestadístico Nacional (MGN) 2018, block (manzana) and section (sección) levels, downloaded as `MGN_URB_MANZANA2018geopackage.zip`, `MGN_RUR_SECCION2018geopackage.zip` and `MGN_URB_SECCION2018geopackage.zip`. Unzip into the folder named by `MGN_DIR`. |
| Municipal boundaries | the `colmaps` R package (`colmaps::municipios`), pinned in the lockfile. |
| Country boundaries (the Colombia-Venezuela border) | GADM, fetched by `geodata::gadm()`. |
| Night-time lights 2018 | NASA Black Marble VNP46A4 (annual), served by NASA LAADS DAAC, <https://ladsweb.modaps.eosdis.nasa.gov>. Fetch it with `blackmarbler::bm_raster(roi_sf = <Colombia boundary>, product_id = "VNP46A4", date = 2018, bearer = <token>)`, where the token is a NASA Earthdata bearer token; this is what `blackmarbler` in `DESCRIPTION` Suggests is for. Write the raster as `geo_cache/ntl_2018.tif`. |
| Travel time to cities 2015 | Malaria Atlas Project, Weiss et al. (2018): `malariaAtlas::getRaster(dataset_id = "Accessibility__201501_Global_Travel_Time_to_Cities", shp = <Colombia boundary>)`. |
| Elevation (for terrain ruggedness) | SRTM, `geodata::elevation_30s()`. |
| Built-up, cropland and tree cover fractions | ESA WorldCover, `geodata::landcover()`. |
| Precipitation 2018 | CHIRPS v2.0 annual, <https://data.chc.ucsb.edu/products/CHIRPS-2.0/global_annual/tifs/>. |
| Precipitation normal 1970-2000 | WorldClim 2.1, `geodata::worldclim_global(var = "prec", res = 2.5)`, summed over the months. |
| Population 2015 and 2018 | WorldPop, unconstrained Global 2000-2020 (`https://data.worldpop.org/GIS/Population/Global_2000_2020/<year>/COL/col_ppp_<year>.tif`) and constrained R2024B (`https://data.worldpop.org/GIS/Population/Global_2015_2030/R2024B/<year>/COL/v1/100m/constrained/col_pop_<year>_CN_100m_R2024B_v1.tif`). |
| Buildings 2018 | Google Open Buildings 2.5D Temporal, public bucket `open-buildings-temporal-data` (no credentials), read by `data-raw/build_buildings_raster.R`. |

## Covariates

The repository builds two covariate matrices from these sources: the
area-level matrix `geo_bld` (census composition, nine geospatial layers and
three building layers per municipality) and the sub-area matrix `geo_blk_bld`
(the same kinds of covariates per census block or rural section). The layer
selections and the build steps are in `METHODS.md` and `data-raw/README.md`.

The published results also used a block of municipal covariates from an
administrative registry that was available to the project under restricted
access. That block is not distributed, and the code does not build it: the
matrices built here are the published ones without it. Other covariates can be
used in two ways, both checked before any fit (rules in `DATA_DICTIONARY.md`,
section 2.8):

- add columns to the built matrices, joined by municipality code or sub-area
  code: `--extra-area=PATH` and `--extra-subarea=PATH`;
- supply a complete matrix in place of the built one, in the layouts of
  `DATA_DICTIONARY.md`: `--area-matrix=PATH` and `--subarea-matrix=PATH`.

Results with other covariates differ from the published ones. The XGB
hyperparameters shipped in `inst/extdata/xgb_sub_params_geo_blk_bld.tsv` were
tuned on the published sub-area matrix; the runner uses them for any sub-area
matrix unless another table is given with `--xgb-params=PATH`, and
`data-raw/tune_xgb_sub.R` repeats the tuning on another matrix.

The per-fit `config` column keeps the labels `geo_bld` and `geo_blk_bld` even
when a run uses `--area-matrix` or `--subarea-matrix`; `run_manifest.json`
records the MD5 of the matrices and of the other input files that distinguish
the run.

## What to expect with the same inputs

With the same inputs (frame, stratum map, truth and covariate matrices) and the
package versions of the lockfile:

- the samples and the truth are identical; the direct estimates and their
  smoothed variances are identical or agree up to floating-point rounding
  (about 1e-14 relative);
- XGB and MEGB give identical fits on the same platform; across the two x86_64
  Linux machines measured below they were bitwise identical;
- FH-arcsin gives the same fits up to floating-point rounding, which depends on
  the BLAS library and the compiler;
- Bayes-Beta gives identical fits only on the same toolchain (CmdStan version,
  C++ compiler and platform); on another toolchain its estimates differ within
  the Monte Carlo error of the sampler.

Measured agreement with the published runs: the fits below were repeated with
this repository's code on the published covariate matrices (given with
`--area-matrix` and `--subarea-matrix`), on x86_64 Linux (Fedora 43) with R
4.5.3, OpenBLAS, CmdStan 2.38.0 and g++ 15.3; the published runs were computed
on another machine and toolchain. The table gives the absolute difference
between the two, over every municipality of each fit, in and out of sample
(MEGB is this repository's `megb_sub`; the two outcomes are
`migrated_last_12m` and `migrated_intl_12m`).

| Method | Fits | Estimate: max | Estimate: median | MSE: max | MSE: median |
|---|---|---|---|---|---|
| FH-arcsin | replicates 1-3, six outcomes | 9.2e-11 | 1.2e-13 | 1.0e-12 | 1.2e-14 |
| XGB | replicates 1-3, six outcomes | 0 | 0 | 0 | 0 |
| MEGB | replicate 1, two outcomes | 0 | 0 | 0 | 0 |
| Bayes-Beta | replicates 1-2, two outcomes | 1.8e-3 | 1.4e-4 | 2.7e-4 | 3.8e-6 |

For Bayes-Beta, the difference in the estimate divided by the posterior
standard deviation has median 0.017 and 95th percentile 0.049 (largest 0.10):
in sample median 0.016 and 95th percentile 0.045, out of sample 0.018 and
0.053. These are a few hundredths of a posterior standard deviation,
consistent with the scale of the sampler's Monte Carlo error. The posterior
standard deviations themselves differ by a median of 2.2% (95th percentile
8.6%).

## Running replicates

```sh
Rscript scripts/run_replicates.R --reps=1-100 --out=results/geo
```

This runs the three refitted methods (`fh_arcsin`, `bayes_fh_beta`,
`xgb_sub`) on all six outcomes. `--methods` chooses others (for example
`--methods=fh_arcsin,bayes_fh_beta,xgb_sub,megb_sub` adds this repository's
MEGB), `--outcomes` a subset of outcomes, and `--workers=N` runs N replicates
at a time in forked processes (on Windows, where they are not available, it
falls back to one replicate at a time, with a message; a Bayes-Beta fit runs
its three chains in parallel, so N workers can use about 3N cores). The
covariate options are those of the previous sections; the header of
`scripts/run_replicates.R` lists every option.

- A rerun completes what is missing and skips what is complete, so a range of
  replicates can be split across jobs and an interrupted run restarted with
  the same command. `--retry-failed` fits again the outcomes whose fit failed.
- Use one `--out` folder per covariate set and settings. The folder keeps a run
  manifest (inputs, settings, package versions), and a run whose inputs,
  settings or package versions differ from it is refused.
- The package checks, the CmdStan check, the Stan model compilation and the
  manifest check run at the start of every run, even when every replicate is
  already complete.
- The BLAS and OpenMP libraries must run one thread per process; the script
  sets `OMP_NUM_THREADS=1` and `OPENBLAS_NUM_THREADS=1` before any package
  loads.

## Evaluation

```sh
Rscript scripts/evaluate.R --run=results/geo \
  --truth="$SAECOL_DATA_DIR/aux_tables/truth_full.parquet" \
  --reps=1-100 --out=results/geo/eval
```

`--external=PATH` adds the estimates of a method fitted elsewhere, labelled
`--external-label` (default `MEGB`; when the run itself holds `megb_sub`, give
the external file another label). `--methods` restricts the evaluation to some
of the run's methods, and `--column` sets the covariate-set label written to
every table (default `geo`). The published comparison is this evaluation over
replicates 1-100 of the three refitted methods together with the MEGB
estimates of an external per-fit file. The command writes one parquet file per
table to `--out`:

| Output | Content | Written by (computed in) |
|---|---|---|
| `area_summary` | per municipality: empirical MSE, bias and variance, estimated MSE and RMSE, their ratios, RB-RMSE, coverage, interval width and score, bias share | `scripts/evaluate.R` (`ur_area_summary()`) |
| `set_summary` | per method, outcome and area set: pooled coverage, width, interval score, pooled MSE ratio and RMSE, and the distribution of the per-municipality measures | `scripts/evaluate.R` (`ur_set_summary()`) |
| `import_check` | rows of the external file and its largest truth difference; fits and rows dropped | `scripts/evaluate.R` (`sim_evaluate()`) |
| `fit_check` | per method and outcome: fits read, failed, dropped for R-hat or for non-finite rows, kept | `scripts/evaluate.R` (`sim_evaluate()`) |
| `by_bin` | the pooled measures by sample-size bin and by prevalence quintile | `scripts/evaluate.R` (`ur_by_bin()`) |
| `cv_flags` | publication at a coefficient of variation below 20%: reported against realised | `scripts/evaluate.R` (`ur_cv_flags()`) |
| `mcse` | Monte Carlo standard errors of coverage, pooled MSE ratio and interval score | `scripts/evaluate.R` (`ur_mcse()`) |
| `reliability` | estimated against empirical RMSE by decile of the estimated RMSE | `scripts/evaluate.R` (`ur_reliability()`) |
| `ratio_summary` | quartiles across municipalities of estimated over empirical RMSE | `scripts/evaluate.R` (`ur_ratio_summary()`) |
| `tracking_summary` | rank correlation of estimated and empirical RMSE across municipalities | `scripts/evaluate.R` (`ur_tracking_summary()`) |
| `size_medians` | median empirical and estimated RMSE by sample-size group | `scripts/evaluate.R` (`ur_size_medians()`) |
| `coverage_summary` | in-sample coverage overall, and the shares of municipalities at or above 0.95 and below 0.80 | `scripts/evaluate.R` (`ur_coverage_summary()`) |

Every column is defined in `DATA_DICTIONARY.md`, section 4.4, and every measure
in `METHODS.md`, section 7.

## External per-fit file

Estimates of a method fitted outside the runner enter the evaluation as a
parquet or CSV file with one row per replicate, outcome and municipality of
the universe:

| Column | Type | Content |
|---|---|---|
| `replicate` | integer | replicate number |
| `outcome` | string | one of the six outcomes |
| `municipality` | integer | DIVIPOLA code |
| `sampled` | logical | TRUE when the municipality was drawn in the replicate |
| `estimate` | double | point estimate |
| `estimated_variance` | double | the method's reported variance, not negative |
| `lower`, `upper` | double | bounds of the 95% interval |
| `truth` | double | census value; must equal `theta_true` of the truth file |

The file is refused, with the reason, when a value is missing, an outcome or a
municipality is unknown, a row repeats, a variance is negative or a lower bound
lies above its upper bound; the evaluation stops when its truth, replicates or
drawn municipalities differ from the run's (`DATA_DICTIONARY.md`, section 4.3).

## Tests

```r
testthat::test_local()
```

Tests that read the census or raster sources skip when `SAECOL_DATA_DIR` is
unset. The Bayesian tests need CmdStan installed and `CMDSTAN` set.

## Citation and licence

Cite the software with the metadata in `CITATION.cff`.

The package's own code is licensed under the GNU General Public License,
version 2 or (at your option) any later version (`License: GPL (>= 2)`; the
GPL-2 text is in `LICENSE`). The Stan model files in `inst/stan/tipsae/` are
copied verbatim from tipsae 1.0.3 and are licensed under the GPL, version 3 or
later (their licence header is in `inst/stan/tipsae/include/license.stan`). The
distributed whole, the package together with those files, is therefore
conveyed under the GNU General Public License, version 3 or any later version.

The methods rest on the Fay-Herriot area-level model; the arcsine transformation
with the bias-corrected back-transformation of Hadam et al. (2020) (emdi2); the
extended Beta model of De Nicolò et al. (2024), as implemented in tipsae; the
robust Fay-Herriot estimator of Sinha and Rao and the synthetic-predictor MSE of
Rao and Molina (2015); gradient boosting and mixed-effects gradient boosting for
small areas in povmap, with the relative bias of the estimated RMSE as the MERF
and MEGB literature defines it (Krennmair and Schmid 2022; Messer and Schmid
2024); the generalised variance function of You and Hidiroglou (2023);
calibration estimators (Deville and Särndal 1992), here raking, that is
iterative proportional fitting of the design weights to the household and person
totals of each stratum; and the interval score of Gneiting and Raftery (2007).
`METHODS.md` cites each where it is used.

Deville, J.-C. and Särndal, C.-E. (1992). Calibration estimators in survey
sampling. *Journal of the American Statistical Association* 87(418), 376-382.

## Contributors

- Luciano Perfetti Villa (author, maintainer)
- David Newhouse (contributor)

Repository: <https://github.com/lucianopv/saecolmig-sim>
