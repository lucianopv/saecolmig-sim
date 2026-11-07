# Methods

What the simulation estimates and how. Column-level layouts are in
`DATA_DICTIONARY.md`; the ECV 2024 design facts the sampler uses are in
`inst/extdata/ecv2024_design.yml`. Function names refer to the package
`saecol2026` in `R/`.

## 1. Estimands

Six outcomes for each municipality, from the CNPV 2018 person question
`PA_VIVIA_1ANO` (where the person lived one year earlier): code 3 (another
municipality of Colombia) is an internal migrant, code 4 (another country) an
international migrant; codes 1 (not yet born) and 2 (same municipality) are
not migrants, and code 9 (no answer) is unknown.

| Level | Any destination | Within Colombia | Abroad |
|---|---|---|---|
| Household indicator: share of households with at least one member who moved; the three overlap, never sum them | `migrated_last_12m` | `migrated_internal_12m` | `migrated_intl_12m` |
| Person ratio: migrants over persons with a known answer (`n_pers_resp`); additive | `pers_migrated_12m` | `pers_migrated_internal_12m` | `pers_migrated_intl_12m` |

The registry is `sae_outcomes` (`R/zz_outcome_registry.R`). The **universe** is
the 1,122 municipalities of the census frame. The census value of an outcome
(the truth, `compute_truth_all()`) is the mean of the household indicator over
the households of the municipality with a known value, or the ratio of the
migrant count to `n_pers_resp` summed over its households; it is computed once,
over the whole frame, for every municipality of the universe.

## 2. Sampling design

Each replicate `r` is a sample drawn from the census household frame under the
design of the ECV 2024 (DANE, microdata catalog 861), with seed `100000 + r`
(`draw_replicate()`, `apply_calibration()`). Replicates 1-100 enter the
published comparison.

- **Strata.** The ECV 2024 stratifies municipalities within department and
  area (cabecera / resto). The 636 municipalities it sampled keep their
  stratum (`inst/extdata/ecv2024_muni_to_stratum.csv`); each other
  municipality takes the stratum of its nearest neighbour in the same
  department on twelve standardised census features (ten census shares, the
  household count and the log of the population; `extend_strata_nn()`). The
  ECV has 357 strata.
- **Stage 1: municipalities.** Within each stratum, the realised number of
  municipalities of the ECV 2024 (`inst/extdata/ecv2024_strata_realised.csv`)
  is drawn with probability proportional to the number of households, without
  replacement (`sampling::UPmaxentropy`, maximum-entropy sampling). A
  municipality whose inclusion probability reaches 1 is taken with certainty;
  a stratum that draws one municipality uses a single draw proportional to
  size. `pi_psu` is the first-stage inclusion probability.
- **Stage 2: segments.** The households of a drawn municipality are ordered by
  their census block code (which carries the sector and section at fixed
  positions: characters 7-11 for rural records, 15-20 for urban ones) and
  cut into consecutive segments of 10, the ECV's segments of about ten
  dwellings. A systematic sample (fractional interval: a random start in
  [0, N/n) and step N/n) selects n of the N segments, n being the stratum's
  realised segments over its realised municipalities, rounded (at least 1, at
  most N). `pi_seg = n / N`.
- **Stage 3: none.** Every household and every person of a selected segment
  is enumerated.
- **Weights.** The design weight `1 / (pi_psu * pi_seg)` is calibrated
  (Deville and Särndal 1992) by raking, that is iterative proportional fitting
  of the weights, to the census household count and person total of each
  stratum (`survey::calibrate(calfun = "raking")`). The ECV calibrates to
  sex x age cells and households per dwelling, which makes its factor constant
  within the household; the two stratum totals approximate it at
  total-persons granularity. Households with an unknown size (0.066% of the
  census) count zero persons in the person margin, on both sides.

The design, the realised allocation and the calibration margins are parsed
from `inst/extdata/ecv2024_design.yml` by `parse_ecv_design()`.

## 3. Direct estimator and variance smoothing

For each replicate and outcome, `direct_estimate_all()` estimates each drawn
municipality under the design `svydesign(ids = ~ssu_id, strata = ~stratum,
weights = ~weight_cal, nest = TRUE)` with `survey.lonely.psu = "average"`:
the design clusters on segments, so the setting applies to the strata in
which a single segment is drawn, which take the average variance contribution
of the other strata. The estimate is the weighted mean of the
household indicator (`svymean`) or the weighted ratio of the person counts
(`svyratio`), with its Taylor-linearised variance clustered on the segment.
`n_eff` is the Kish effective sample size `(sum w)^2 / sum w^2` and `n_dwell`
the number of sampled households.

The sampling variances are smoothed by a generalised variance function fitted
separately for each outcome (`smooth_variance_all(method = "ghs")`; pooling
outcomes of different prevalence would bias the fit), in the spirit of You and
Hidiroglou (2023): an OLS fit of `log(var_direct) ~ log(n_eff) + theta_direct
+ theta_direct^2` on the municipalities with a positive variance, whose
exponentiated prediction is `var_direct_smoothed` for every municipality with
an estimate, including those whose raw variance is zero (frequent for the
international outcomes). Where fewer than five municipalities of an outcome have
a positive variance, the raw variances are kept and the others take the smallest
positive one (at least 1e-8). The untransformed area-level models (`fh`,
`fh_spatial`, `fh_robust`, `bayes_fh_heavy`) take `var_direct_smoothed` as the
known sampling variance; `fh_arcsin` and `bayes_fh_beta` use it differently
(sections 5.1 and 5.2).

## 4. Covariates

**Area level, configuration `geo_bld`**
(`aux_tables/municipal_X_geo_bld.parquet`):

- ten census shares per municipality (`build_aux_municipal()`): urban
  residence, working age (15-64), elderly (65 and over), secondary education
  completed and aged 25 or over, illiterate and aged 15 or over, worked in
  the reference week, Afro-Colombian, indigenous (person shares), piped
  water (household share) and electricity (dwelling share);
- nine geospatial layers, population-weighted means over the municipal polygons
  with the unconstrained WorldPop 2018 weights unless stated: night-time lights,
  travel time to cities, terrain ruggedness, distance from the centroid to the
  Venezuelan border (not weighted), population density (WorldPop 2018 population
  over the polygon area), cropland and tree cover fractions, 2018 precipitation
  and its anomaly against the 1970-2000 normal (`build_geo_aux()`);
- three building layers, population-weighted from the Open Buildings 2.5D
  Temporal raster of 2018 on the WorldPop 3-arcsecond grid: buildings per km2
  and the presence-weighted mean and standard deviation of building height
  (`bld_pw_covariates()`).

A municipality without a value for a geospatial or building layer (a raster
with no data over it) takes the median over municipalities. The geo and
building columns are fixed by `inst/extdata/kept_geo_cols.txt` and
`kept_bld_cols.txt` (`kept_cols()`), the selections of a VIF prune with
threshold 10 (`bld_prune_vif()`) against the other covariates of the parent
matrix: the geo block without the built-up fraction and the building block
without building presence. A prune against other covariates can keep other
columns, so the builders apply the lists. The population and household counts
are not covariates.

**Sub-area level, configuration `geo_blk_bld`**
(`subarea/X_geo_blk_bld.parquet`). A sub-area is a census block in urban areas
and a rural section (the block code cut to its first 11 characters: department,
municipality, class, rural sector and rural section) elsewhere. Per sub-area:
the nine census shares other than urban residence (constant within a sub-area),
an urban flag, the same nine geospatial layers and three building layers
extracted over the sub-area polygons of the DANE MGN 2018, weighted by the
constrained WorldPop R2024B 2018 raster. A sub-area whose extraction is missing
takes its municipality's value of the same layer on the same constrained weights
(`municipal_geo_2018con`, `municipal_bld_2018con`). The prediction universe of
the sub-area methods is the set of sub-areas complete on every covariate
(`subarea_universe()`); the census household and person counts weight the
roll-up to municipality.

The municipal fallback of the sub-area geo block
(`data-raw/build_subarea_geo.R municipal`) reads the CHIRPS precipitation
rasters as stored, with spherical geometry on as for the area-level block: the
CHIRPS no-data cells (stored as -9999) are not masked there, while the
block-level pass and the area-level geo block mask them. The constrained
population weights are zero on all but a few coastal no-data cells, so only the
precipitation columns (`precip_annual_pw`, `precip_anomaly_pw`) of the
municipalities holding such cells take them in, and only for the sub-areas that
fall back to their municipality's value. This reproduces the published
sub-area matrix; masking those cells would move these values slightly.

The published matrices of both levels also held a block of municipal
covariates from an administrative registry under restricted access, which this
code does not build (see `README.md`).

## 5. Methods

Every method returns, for every municipality of the universe, a point
estimate, its reported uncertainty as a variance (`mse`) and a 95% interval,
in the shape checked by `sae_fit_contract()`; drawn municipalities are in
sample and the others out of sample. Each fit draws from its own seed,
`sae_fit_seed(config, method, outcome, replicate)` = 1e7 x configuration index +
1e6 x method index + 1e4 x outcome index + replicate (`R/fit_seed.R`).

### 5.1 FH-arcsin (`fh_arcsin`, area level)

The Fay-Herriot model on the arcsine scale, fitted by REML with `emdi2::fh`
(`transformation = "arcsin"`, `backtransformation = "bc"`): the direct
estimate is transformed with `asin(sqrt(.))`, its sampling variance on that
scale is `1 / (4 n_eff)` (`eff_smpsize = n_eff`, with `n_eff` floored at 1),
so the smoothed variance `var_direct_smoothed` plays no part, and the
predictions are back-transformed with the bias correction of Hadam et al.
(2020) that emdi2 implements. Before the transformation the direct estimates
are clipped to [1 / (2 n_eff), 1 - 1 / (2 n_eff)]. The MSE is a parametric
bootstrap with B = 500 and the fit's seed; where a bootstrap leaves a
municipality's MSE missing it takes the median MSE of the in-sample
municipalities. The interval is the estimate plus or minus 1.96 root MSE, cut to
[0, 1]. Out of sample emdi2 returns the back-transformed synthetic prediction
`sin^2(x'beta)` without the bias correction.

### 5.2 Bayes-Beta (`bayes_fh_beta`, area level)

The Bayesian Fay-Herriot model with the extended Beta likelihood of De Nicolò
et al. (2024), as implemented in tipsae 1.0.3. The direct estimate of area d
has an extended Beta likelihood with mean parameter
`mu_d = logit^-1(x_d' beta + v_d)`, a normal random effect `v_d`, a
dispersion smoothed by `tipsae::smoothing()`, and probabilities p0_d and p1_d
of the exact values 0 and 1 that depend on mu_d, a parameter lambda and the
area sample size m_d (`n_dwell`); its in-sample estimand is
`theta_d = mu_d (1 - p0_d - p1_d) + p1_d`. tipsae's
Stan model (vendored verbatim in `inst/stan/tipsae/`) is sampled by CmdStan on
the data tipsae assembles (`capture_tipsae_standata()`), with 3 chains of 2,000
iterations (half warm-up), `adapt_delta = 0.95` and `max_treedepth = 10`. The
direct estimates, exact zeros included, enter the likelihood as they are: the
extended Beta has a point mass at 0 and fits a zero directly. The dispersion is
smoothed by `tipsae::smoothing()` with `var_direct_smoothed` as its raw
variance, and an area with a zero estimate enters that regression at zero. The
chains start at zero on the unconstrained scale except the intercept, which
starts at the logit of the mean direct estimate (`extbeta_init_intercept()`),
because from tipsae's own start the zero mass of an area with an exact zero is
0.5^m, which is 0 in double precision from about 1,075 sampled households; the
rstan backend of `fit_bayes_fh_beta()` starts its chains at those values. The
estimate, its variance and the interval are the posterior mean, variance and
2.5% and 97.5% quantiles. Out of sample the estimate is tipsae's `theta_oos`,
which is mu at a new draw of the random effect, with no m_d, p0 or p1 (the limit
of theta for large m). mu is the mean of the Beta part of the mixture, which
leaves out the point masses at 0 and 1; it exceeds theta where the zero mass is
large, that is where zeros are common, above all for the international outcomes.
The evaluation scores `theta_hat`, which out of sample is mu. The evaluation
drops a fit whose `max_rhat` exceeds 1.05: the largest rank-normalised split
R-hat of the `posterior` package over the variables the backend keeps (every
variable of the model and its generated quantities on CmdStan; on rstan `beta0`,
`beta`, `v`, `sigma_v`, `mu`, `p0`, `p1`, `lambda_EB`, `theta`, `theta_oos`,
`y_rep`, `log_lik` and `lp__`).

### 5.3 XGB (`xgb_sub`, sub-area level)

Gradient-boosted regression of the sub-area outcome on the sub-area covariates
with `povmap::xgb()`, census blocks and rural sections as sub-domains of their
municipality. The sample enters as one row per sub-area realised in the
replicate: the design-weighted mean of the household indicator, or the
weighted ratio of the person counts, with the summed weight
(`collapse_sample_to_subarea()`). Predictions over the prediction universe are
aggregated to municipality with the census household counts (household
outcomes) or person counts (person outcomes). No transformation is applied (a
sub-area holds about one segment, so most sub-area outcomes are exactly 0 or
1); the estimates are not benchmarked; the uncertainty is povmap's residual
bootstrap with B = 200. The hyperparameters are tuned per outcome
(`data-raw/tune_xgb_sub.R`): a five-fold cross-validated grid over learning rate
(0.02, 0.05, 0.10, 0.20), depth (2, 3, 4, 6) and minimum child weight (1, 5)
with 200 rounds, subsample and column sample 0.8, lambda 1 and alpha 0, on each
of the first five replicates, each parameter then set to its most frequent
winner.

### 5.4 MEGB (`megb_sub`, sub-area level)

Mixed-effects gradient boosting (`povmap::megb()`) on the same inputs as XGB,
with a municipal random effect. The estimator is cross-fitted: ten folds hold
out whole municipalities, a sampled municipality is predicted by the fold model
that held it out and an unsampled one by the average of the fold models, with
early stopping. The uncertainty is povmap's parametric bootstrap with B = 100,
reported as the variance of the centred bootstrap prediction errors with a
basic bootstrap interval (`mse_type = "var"`, which leaves each municipality's
bootstrap bias out). These settings are passed explicitly, so a global option
cannot change the estimator. In the published comparison the MEGB estimates
came from an external per-fit file.

### 5.5 Other registered methods

`fh` (the Fay-Herriot model with `emdi2::fh`), `fh_spatial` (with a SAR random
effect on the queens-contiguity weights of `build_queens_W()`), `fh_robust`
(the robust Fay-Herriot estimator of Sinha and Rao, `saeRobust::rfh()` with
Huber constant 1.345; out of sample the synthetic prediction with the
leading-term MSE `sigma2_v + x' Var(beta) x` of Rao and Molina (2015, section
6.2.2)) and `bayes_fh_heavy` (a Gaussian Fay-Herriot model in Stan with a
half-Normal prior on the random-effect scale) satisfy the same contract and
can be run with `--methods`; they are not part of the published comparison.

## 6. Running and recording fits

`scripts/run_replicates.R` (`sim_run()`) draws each replicate, writes its
direct estimates and one file of per-fit rows per method, and records a failed
fit as rows with `fit_failed = TRUE` rather than stopping. The per-fit rows
carry the seed, the package versions and the povmap commit of each fit
(`DATA_DICTIONARY.md`, section 4.2).

## 7. Uncertainty measures

`scripts/evaluate.R` (`sim_evaluate()`) joins the per-fit rows of replicates
1-100 to the truth. A fit is dropped whole when it failed or its `max_rhat`
(section 5.2) exceeds 1.05, and rows with a non-finite estimate, MSE or bound
are dropped (`ur_drop_bad_fits()`); `fit_check` counts both. Every measure is
computed separately by method, outcome and area set (in sample or out of
sample). For a municipality with estimates `theta_hat_r`, reported variances
`mse_r` and intervals `[l_r, u_r]` over its `R` replicates in the set, and truth
`theta` (`ur_area_summary()`):

- empirical MSE `emp_mse = mean_r (theta_hat_r - theta)^2`, which splits into
  `bias^2 + emp_var`, with `bias = mean_r (theta_hat_r - theta)` and `emp_var`
  the variance of the estimates around their mean; the bias share is
  `bias^2 / emp_mse`;
- estimated MSE `rep_mse = mean_r mse_r` and estimated RMSE
  `rep_rmse = mean_r sqrt(mse_r)`; their ratios to the empirical values,
  `ratio = rep_mse / emp_mse` and `ratio_var = rep_mse / emp_var`;
- the relative bias of the estimated RMSE,
  `rb_rmse = rep_rmse / sqrt(emp_mse) - 1`, as the MERF and MEGB literature
  defines it (Krennmair and Schmid 2022; Messer and Schmid 2024);
- `rrmse_mse`, the RMSE of `mse_r` around `emp_mse`, relative to `emp_mse`;
- coverage, the share of intervals that hold the truth, and the mean interval
  width; `width_rel`, the mean width over that of a normal 95% interval at the
  empirical RMSE, `2 z sqrt(emp_mse)` with `z = qnorm(0.975)`;
- the interval score of Gneiting and Raftery (2007),
  `(u - l) + (2 / alpha) (l - theta)_+ + (2 / alpha) (theta - u)_+` with
  `alpha = 0.05`, averaged over replicates.

The ratios are missing where `emp_mse` (or `emp_var`) is zero. Over the
municipalities of a set (`ur_set_summary()`): coverage, width, interval score,
the pooled ratio `mean(mse) / mean(squared error)` and the pooled empirical
and estimated RMSE over all rows; and, over the municipalities with at least 20
replicates in the set, the quartiles of `ratio` and its shares below 0.5 and
above 2, the medians of the other per-municipality measures, the mean of
`rb_rmse`, the Spearman correlation of estimated and empirical RMSE, and the
shares of municipalities with coverage under 0.90 and 0.80.

- **By sample size and prevalence** (`ur_by_bin()`): the pooled measures within
  eight quantile bins of the households drawn per municipality and replicate,
  with the out-of-sample rows as a bin of their own, and within quintiles of
  the true value.
- **Publication at CV < 20%** (`ur_cv_flags()`): an estimate is published when
  its reported coefficient of variation `sqrt(mse) / theta_hat` is below 0.20
  and publishable when the realised one, `sqrt(emp_mse) / theta`, is; the table
  gives the share published on each, the share of published estimates that are
  not publishable and the share of publishable ones that are not published.
- **Reliability** (`ur_reliability()`): estimated and empirical RMSE within
  deciles of the estimated RMSE.
- **Monte Carlo standard errors** (`ur_mcse()`): of coverage, the pooled ratio
  and the interval score, from 200 bootstrap resamples of whole replicates
  shared by all methods.
- **Display summaries** (`R/uncertainty_display.R`), per outcome: the quartiles
  of `rep_rmse / sqrt(emp_mse)` across municipalities and the share below 1
  (`ur_ratio_summary()`); the rank correlation of `rep_rmse` with
  `sqrt(emp_mse)` (`ur_tracking_summary()`); the median empirical and estimated
  RMSE by groups of at least 25 in-sample municipalities with similar median
  sample sizes (`ur_size_medians()`); and in-sample coverage overall and the
  shares of municipalities at or above 0.95 and below 0.80
  (`ur_coverage_summary()`).

The evaluation requires every method to see the same drawn municipalities in
each replicate and outcome (`ur_assert_same_sets()`).

## References

Cited by author and year, as in the code, except the calibration and
back-transformation papers:

- De Nicolò et al. (2024): the extended Beta model of `bayes_fh_beta`.
- Deville, J.-C. and Särndal, C.-E. (1992). Calibration estimators in survey
  sampling. *Journal of the American Statistical Association* 87(418),
  376-382: the calibration of the replicate weights.
- Gneiting and Raftery (2007): the interval score.
- Hadam, S., Wuerz, N. and Kreutzmann, A.-K. (2020). Estimating regional
  unemployment with mobile network data for Functional Urban Areas in Germany.
  Refubium - Freie Universitaet Berlin Repository, 1-28: the bias-corrected
  back-transformation of the arcsine model.
- Krennmair and Schmid (2022); Messer and Schmid (2024): mixed-effects random
  forests and gradient boosting for small areas, and the relative bias of the
  estimated RMSE.
- Rao and Molina (2015, section 6.2.2): the synthetic-predictor MSE.
- Sinha and Rao: the robust Fay-Herriot estimator of `fh_robust`.
- Weiss et al. (2018): travel time to cities.
- You and Hidiroglou (2023): generalised variance functions.
