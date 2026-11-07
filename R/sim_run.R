#' Replicate runner
#'
#' One replicate `r` of the simulation is: draw a sample from the household
#' frame under the ECV design with seed `100000 + r` and calibrate its weights
#' ([sim_draw()]); compute the direct estimate of every outcome in every drawn
#' municipality and smooth the sampling variances with a generalised variance
#' function fitted per outcome; collapse the sample to sub-areas for the
#' sub-area methods; fit each requested method to each requested outcome with
#' the fit's own seed ([sim_fit_one()], seed from [sae_fit_seed()]); and write
#' the direct estimates and one file of per-fit rows per method
#' ([sim_run_replicate()]). [sim_run()] runs a range of replicates, in forked
#' worker processes if asked. File layouts are in `DATA_DICTIONARY.md`.
#'
#' A municipality drawn in the replicate is in sample; every other municipality
#' of the universe is out of sample and is predicted from its covariates alone.
#' Area-level methods are fitted on the area-level matrix (configuration
#' `inputs$area_config`), sub-area methods on the sub-area matrix
#' (`inputs$subarea_config`); the configuration names the output directory and
#' enters the seed.
#'
#' Every method except `fh_robust` draws its random numbers from its own seed,
#' so its fit gives the same result whichever process runs it and in whatever
#' order; the bootstrap MSE of `fh_robust` draws from the random-number stream
#' of the R process and so depends on what ran before it. The BLAS and OpenMP
#' libraries must run one thread per process
#' (`OMP_NUM_THREADS=1`, `OPENBLAS_NUM_THREADS=1`, set before any package
#' loads); several threads per fit oversubscribe the cores when fits run in
#' parallel. `scripts/run_replicates.R` sets them.
#'
#' @name sim_run_overview
NULL

#' Methods refitted for the published comparison
#'
#' Fay-Herriot with the arcsine transformation, Bayesian Fay-Herriot with the
#' extended Beta likelihood, and sub-area XGBoost.
#' @export
sim_v3_methods <- c("fh_arcsin", "bayes_fh_beta", "xgb_sub")

#' Columns of a per-fit results file
#'
#' The columns [sim_fit_one()] returns and `<out_dir>/fits/` files hold, in
#' this order; `DATA_DICTIONARY.md` defines each.
#' @export
sim_fit_columns <- c(
  "area_id", "theta_hat", "mse", "ci_lo", "ci_hi", "area_set", "fit_failed",
  "error_msg", "config", "method", "replicate", "outcome", "seed", "wall_secs",
  "weight", "n_divergent", "max_rhat", "sigma2_v", "povmap_sha", "mse_type",
  "pkg_versions")

# Columns of a direct-estimates file, in order.
.sim_direct_columns <- c("area_id", "theta_direct", "var_direct", "n_eff",
                         "n_dwell", "outcome", "var_direct_smoothed")

# Methods that sample with CmdStan whatever their settings.
.sim_cmdstan_methods <- "bayes_fh_heavy"

#' Fixed settings of each method
#'
#' The arguments every fit of a method receives besides its data and its seed,
#' in one place. Area-level methods take the arguments of [method_args()] with
#' these laid over them; sub-area methods take these directly.
#'
#' * `fh_arcsin`: `B = 500` bootstrap replicates for the MSE of the
#'   back-transformed estimate.
#' * `bayes_fh_beta`: the extended Beta likelihood, 3 chains of 2,000
#'   iterations (half warm-up) sampled by CmdStan. The direct estimates, exact
#'   zeros included, go to the likelihood unmodified, and out of sample the
#'   estimate is tipsae's `theta_oos`. The CmdStan sampler runs with
#'   `adapt_delta = 0.95` and `max_treedepth = 10`; every chain starts at zero
#'   on the unconstrained scale except the intercept, which starts at the logit
#'   of the mean direct estimate ([fit_bayes_fh_beta_cmdstan()],
#'   [extbeta_init_intercept()]).
#' * `xgb_sub`: `B = 200` residual-bootstrap replicates, sub-areas rolled up to
#'   municipality by census counts, and the tuned hyperparameters of each
#'   outcome.
#' * `megb_sub`: `B = 100` bootstrap replicates, census counts, and the
#'   variance of the centred bootstrap prediction errors as its uncertainty
#'   (`mse_type = "var"`), with povmap's default boosting parameters.
#'
#' Methods without an entry run with the defaults of [method_args()] and of
#' their fit function.
#'
#' @param xgb_params named list, outcome -> list of XGB hyperparameters, as
#'   [xgb_sub_params()] returns; `NULL` leaves `xgb_sub` without
#'   hyperparameters, and the runner then refuses it
#' @return named list, method -> list of arguments
#' @export
sim_settings <- function(xgb_params = NULL) {
  list(
    fh_arcsin = list(B = 500L),
    bayes_fh_beta = list(likelihood = "ExtBeta", chains = 3L, iter = 2000L,
                         backend = "cmdstan"),
    xgb_sub = list(B = 200L, weight = "census", params = xgb_params),
    megb_sub = list(B = 100L, weight = "census", mse_type = "var"))
}

# A replicate number as an integer from 1 to 9999 (the range of the seeds).
.sim_check_rep <- function(rep) {
  ok <- length(rep) == 1L && is.numeric(rep) && !is.na(rep) && rep == round(rep) &&
    rep >= 1 && rep <= 9999
  if (!ok) .sim_stop("a replicate must be a whole number from 1 to 9999")
  as.integer(rep)
}

# Covariate configuration a method is fitted on.
.sim_config <- function(method, inputs) {
  if (method %in% sae_methods_subarea) inputs$subarea_config else inputs$area_config
}

.sim_fit_path <- function(out_dir, method, config, rep) {
  file.path(out_dir, "fits", paste0(method, "__", config), sprintf("rep_%04d.parquet", rep))
}

.sim_direct_path <- function(out_dir, rep) {
  file.path(out_dir, "direct", sprintf("rep_%04d.parquet", rep))
}

.sim_xgb_message <- function(config, outcomes) {
  paste0("xgb_sub needs tuned hyperparameters for every outcome it fits (missing: ",
         paste(outcomes, collapse = ", "), "); none ship for sub-area configuration '",
         config, "'. Tune them with data-raw/tune_xgb_sub.R and pass the table: ",
         "scripts/run_replicates.R --xgb-params=PATH, or ",
         "inputs$tune <- xgb_sub_params(file = PATH) before calling the runner.")
}

# Packages each method needs at fit time. The FH family fits with emdi2 (the
# robust variant with saeRobust); Bayes-Beta builds its data with tipsae and
# rstan and, on the cmdstan backend, samples with cmdstanr; the sub-area methods
# fit with povmap on xgboost.
.sim_method_packages <- function(method, settings) {
  switch(method,
         fh = , fh_arcsin = , fh_spatial = "emdi2",
         fh_robust = "saeRobust",
         bayes_fh_beta = c("tipsae", "rstan",
                           if (identical(settings$bayes_fh_beta$backend, "cmdstan")) "cmdstanr"
                           else "posterior"),
         bayes_fh_heavy = "cmdstanr",
         xgb_sub = , megb_sub = c("povmap", "xgboost"),
         character(0))
}

# Whether a package can be loaded; loads its namespace, so forked workers start
# with it loaded.
.sim_has_pkg <- function(pkg) requireNamespace(pkg, quietly = TRUE)

# Whether any requested method samples with CmdStan.
.sim_uses_cmdstan <- function(methods, settings) {
  any(methods %in% .sim_cmdstan_methods) ||
    ("bayes_fh_beta" %in% methods && identical(settings$bayes_fh_beta$backend, "cmdstan"))
}

# Stop on a request the runner cannot carry out, before anything is drawn: an
# unknown method, outcome or configuration, a missing input, missing XGB
# hyperparameters, CmdStan not found, or a package a method needs that does not
# load. A fit that fails for one of these reasons would fail again on every
# replicate, so the run is refused instead of recording the failures.
.sim_check_request <- function(methods, outcomes, inputs, settings) {
  if (!length(methods)) .sim_stop("no method requested")
  if (!length(outcomes)) .sim_stop("no outcome requested")
  bad <- setdiff(methods, sae_methods_implemented)
  if (length(bad)) {
    .sim_stop("unknown method(s): ", paste(bad, collapse = ", "), "; implemented: ",
              paste(sae_methods_implemented, collapse = ", "))
  }
  bad <- setdiff(outcomes, sae_outcomes$outcome)
  if (length(bad)) {
    .sim_stop("unknown outcome(s): ", paste(bad, collapse = ", "), "; outcomes: ",
              paste(sae_outcomes$outcome, collapse = ", "))
  }
  if (anyDuplicated(methods)) .sim_stop("a method is requested more than once")
  if (anyDuplicated(outcomes)) .sim_stop("an outcome is requested more than once")
  configs <- vapply(methods, .sim_config, character(1), inputs = inputs)
  sae_fit_seed(configs, methods, outcomes[1], 1L)
  if (any(methods %in% sae_methods_subarea) && is.null(inputs$pop_sub)) {
    .sim_stop("the sub-area methods need the sub-area matrix (inputs$pop_sub)")
  }
  if ("fh_spatial" %in% methods && is.null(inputs$W)) {
    .sim_stop("fh_spatial needs the spatial weights aux_tables/W_queens_municipal.parquet")
  }
  if ("xgb_sub" %in% methods) {
    miss <- outcomes[vapply(outcomes, function(o) is.null(settings$xgb_sub$params[[o]]),
                            logical(1))]
    if (length(miss)) .sim_stop(.sim_xgb_message(inputs$subarea_config, miss))
  }
  # CMDSTAN before anything loads cmdstanr (the package check below loads it).
  if (.sim_uses_cmdstan(methods, settings) && !nzchar(Sys.getenv("CMDSTAN")) &&
      !"cmdstanr" %in% loadedNamespaces()) {
    sim_cmdstan_env(required = TRUE)
  }
  need <- lapply(methods, .sim_method_packages, settings = settings)
  names(need) <- methods
  pkgs <- unique(unlist(need))
  ok <- vapply(pkgs, .sim_has_pkg, logical(1))
  if (!all(ok)) {
    lost <- pkgs[!ok]
    who <- vapply(names(need)[vapply(need, function(p) any(p %in% lost), logical(1))],
                  function(m) paste0(m, " (", paste(intersect(need[[m]], lost), collapse = ", "), ")"),
                  character(1))
    .sim_stop("packages that the requested methods need do not load: ",
              paste(lost, collapse = ", "), "; needed by ", paste(who, collapse = ", "),
              ". Install them (renv::restore() installs the pinned versions).")
  }
  invisible(TRUE)
}

# Compile the Stan models of the requested methods once, in this process, so a
# compile error stops the run before anything is written and forked workers
# share the compiled model instead of compiling it at the same time.
.sim_prepare_models <- function(methods, settings) {
  build <- function(label, f) {
    tryCatch(f(), error = function(e) {
      .sim_stop("the Stan model of ", label, " could not be compiled: ", conditionMessage(e))
    })
  }
  if ("bayes_fh_beta" %in% methods && identical(settings$bayes_fh_beta$backend, "cmdstan")) {
    build("bayes_fh_beta", .tipsae_cmdstan_model)
  }
  if ("bayes_fh_heavy" %in% methods) build("bayes_fh_heavy", .bayes_fh_heavy_model)
  invisible(TRUE)
}

# ---- run manifest ----------------------------------------------------------

.sim_manifest_file <- "run_manifest.json"

# What the results of a run are computed from: the configuration labels, the
# settings of every method, the MD5 of the area-level and sub-area matrices the
# fits read, the source of each input (path as given, size, MD5), the source of
# the XGB hyperparameters and the package versions.
.sim_manifest <- function(inputs, settings) {
  pv <- strsplit(strsplit(sim_pkg_versions(), ";", fixed = TRUE)[[1]], "=", fixed = TRUE)
  pv <- stats::setNames(lapply(pv, `[`, 2L), vapply(pv, `[`, character(1), 1L))
  src <- inputs$sources
  tune_file <- attr(inputs$tune, "source")
  src["xgb_params"] <- list(
    if (is.null(inputs$tune)) NULL
    else if (!is.null(tune_file) && file.exists(tune_file)) .sim_source(tune_file)
    else list(path = "<list>", size = NA_real_, md5 = .sim_md5_tune(inputs$tune)))
  list(area_config = inputs$area_config,
       subarea_config = inputs$subarea_config,
       settings = settings,
       data = list(area_matrix_md5 = .sim_md5_object(inputs$X_area),
                   subarea_matrix_md5 = if (!is.null(inputs$pop_sub)) .sim_md5_object(inputs$pop_sub)),
       sources = src,
       pkg_versions = pv)
}

# MD5 of hyperparameters given as a list: the values only, without the path
# recorded in the attribute `source`, so the same values read from two files
# hash alike.
.sim_md5_tune <- function(tune) {
  x <- unclass(tune)
  attr(x, "source") <- NULL
  .sim_md5_object(x)
}

.sim_json <- function(x) {
  as.character(jsonlite::toJSON(x, auto_unbox = TRUE, digits = NA, null = "null",
                                na = "null", pretty = TRUE))
}

# The compared fields of a manifest parsed from JSON, as JSON text keyed by
# field: labels, settings, matrix contents, the size and MD5 of each source (the
# path is recorded but not compared) and each package version.
.sim_manifest_fields <- function(m) {
  js <- function(x) as.character(jsonlite::toJSON(x, auto_unbox = TRUE, digits = NA,
                                                  null = "null", na = "null"))
  out <- c(area_config = js(m$area_config), subarea_config = js(m$subarea_config),
           settings = js(m$settings))
  for (k in names(m$data)) out[[paste0("data.", k)]] <- js(m$data[[k]])
  for (k in names(m$sources)) {
    v <- m$sources[[k]]
    out[[paste0("sources.", k)]] <- js(if (is.null(v)) NULL else v[c("size", "md5")])
  }
  for (k in names(m$pkg_versions)) out[[paste0("pkg_versions.", k)]] <- js(m$pkg_versions[[k]])
  out
}

# Write the run manifest of `out_dir` on first use; on later use stop when the
# inputs, settings or package versions differ from those recorded, naming the
# differing fields. A package version recorded as NA (not installed, or CmdStan
# not located) is not compared.
.sim_check_manifest <- function(out_dir, inputs, settings) {
  path <- file.path(out_dir, .sim_manifest_file)
  txt <- .sim_json(.sim_manifest(inputs, settings))
  if (!file.exists(path)) {
    .sim_write_atomic(path, function(tmp) writeLines(txt, tmp))
    return(invisible(path))
  }
  old <- tryCatch(jsonlite::fromJSON(path, simplifyVector = FALSE), error = function(e) NULL)
  if (is.null(old)) .sim_stop(path, " cannot be read; write the results to another folder")
  a <- .sim_manifest_fields(old)
  b <- .sim_manifest_fields(jsonlite::fromJSON(txt, simplifyVector = FALSE))
  keys <- union(names(a), names(b))
  na_pkg <- function(v) is.na(v) || v %in% c("\"NA\"", "null")
  differ <- keys[vapply(keys, function(k) {
    x <- if (k %in% names(a)) a[[k]] else NA_character_
    y <- if (k %in% names(b)) b[[k]] else NA_character_
    if (startsWith(k, "pkg_versions.") && (na_pkg(x) || na_pkg(y))) return(FALSE)
    !identical(x, y)
  }, logical(1))]
  if (length(differ)) {
    .sim_stop(out_dir, " holds results computed from other inputs, settings or package ",
              "versions (", .sim_manifest_file, " differs in: ", paste(differ, collapse = ", "),
              "); write these results to another folder")
  }
  invisible(path)
}

# Every check and preparation that precedes the first draw of a run.
.sim_prepare <- function(methods, outcomes, inputs, settings, out_dir) {
  .sim_check_request(methods, outcomes, inputs, settings)
  .sim_prepare_models(methods, settings)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  .sim_check_manifest(out_dir, inputs, settings)
  invisible(TRUE)
}

#' Locate CmdStan through the CMDSTAN environment variable
#'
#' cmdstanr reads the CmdStan installation from `CMDSTAN`. When the variable is
#' unset at the time cmdstanr loads, cmdstanr first probes for an installation
#' under the Windows Subsystem for Linux by starting a process; on Linux that
#' start fails and the failed process handle can later end the whole R process
#' without a message. This function therefore sets `CMDSTAN`, before cmdstanr
#' loads, to the newest `cmdstan-<version>` folder in cmdstanr's default install
#' location `<home>/.cmdstan` when the variable is unset.
#'
#' @param home home folder whose `.cmdstan` folder is searched
#' @param required stop, with a message on how to set `CMDSTAN`, when no
#'   installation is found
#' @return the CmdStan folder in use, or `""` when none is found and `required`
#'   is FALSE
#' @export
sim_cmdstan_env <- function(home = Sys.getenv("HOME"), required = FALSE) {
  cur <- Sys.getenv("CMDSTAN")
  if (nzchar(cur)) return(cur)
  dirs <- if (nzchar(home)) Sys.glob(file.path(home, ".cmdstan", "cmdstan-*")) else character(0)
  dirs <- dirs[dir.exists(dirs)]
  if (length(dirs)) {
    ver <- suppressWarnings(numeric_version(sub("^cmdstan-", "", basename(dirs)), strict = FALSE))
    pick <- dirs[order(ver, na.last = FALSE)][length(dirs)]
    Sys.setenv(CMDSTAN = pick)
    return(pick)
  }
  if (required) {
    .sim_stop("CmdStan was not found. Set the environment variable CMDSTAN to a CmdStan ",
              "installation folder before cmdstanr loads (cmdstanr::install_cmdstan() ",
              "installs one under ~/.cmdstan/cmdstan-<version>), for example ",
              "CMDSTAN=~/.cmdstan/cmdstan-2.38.0 in .Renviron or the shell.")
  }
  ""
}

# CmdStan version from the makefile of the installation CMDSTAN names, else from
# cmdstanr when it is already loaded; NA when neither is available.
.sim_cmdstan_version <- function() {
  path <- Sys.getenv("CMDSTAN")
  if (nzchar(path)) {
    mk <- file.path(path, "makefile")
    if (file.exists(mk)) {
      hit <- grep("^CMDSTAN_VERSION[[:space:]]*:=", readLines(mk, warn = FALSE), value = TRUE)
      if (length(hit)) return(trimws(sub("^CMDSTAN_VERSION[[:space:]]*:=", "", hit[1])))
    }
  }
  if ("cmdstanr" %in% loadedNamespaces()) {
    v <- tryCatch(as.character(cmdstanr::cmdstan_version(error_on_NA = FALSE)),
                  error = function(e) NA_character_)
    if (length(v) == 1L && !is.na(v)) return(v)
  }
  NA_character_
}

#' Versions of the packages behind the fits
#'
#' One string, `name=version` pairs separated by `;`, for emdi2 (Fay-Herriot
#' fits), povmap and xgboost (sub-area methods), tipsae (the extended Beta
#' model) and CmdStan (its sampler). A package installed from a repository
#' commit carries the first seven characters of the commit after `@`, because
#' a version number does not change between commits. `NA` marks a package that
#' is not installed.
#'
#' @return character scalar, e.g.
#'   `"emdi2=2.2.4@4a50df5;povmap=1.0.1@2b3959a;xgboost=3.1.2.1;tipsae=1.0.3;cmdstan=2.38.0"`
#' @export
sim_pkg_versions <- function() {
  one <- function(p) {
    d <- suppressWarnings(utils::packageDescription(p))
    if (!inherits(d, "packageDescription")) return(paste0(p, "=NA"))
    v <- d$Version
    sha <- d$RemoteSha
    if (!is.null(sha) && grepl("^[0-9a-f]{7,40}$", sha)) v <- paste0(v, "@", substr(sha, 1L, 7L))
    paste0(p, "=", v)
  }
  pk <- vapply(c("emdi2", "povmap", "xgboost", "tipsae"), one, character(1))
  cs <- .sim_cmdstan_version()
  paste(c(pk, paste0("cmdstan=", if (is.na(cs)) "NA" else cs)), collapse = ";")
}

#' Draw one replicate
#'
#' Draws the sample of replicate `rep` with seed `100000 + rep` from the frame
#' under the ECV design, reusing the frame's segments ([draw_replicate()]), and
#' rakes its design weights to the household and person totals of each stratum
#' ([apply_calibration()]). The direct estimates of all six outcomes come from
#' [direct_estimate_all()] by municipality, with the sampling variances smoothed
#' by a generalised variance function fitted per outcome
#' ([smooth_variance_all()], `method = "ghs"`). The sub-area sample is the
#' calibrated sample collapsed to the sub-areas of `inputs$pop_sub`
#' ([collapse_sample_to_subarea()]) for each outcome; households in sub-areas
#' outside the prediction universe are left out.
#'
#' @param rep replicate number, 1-9999
#' @param inputs list from [sim_inputs()]
#' @return list: `sample` (the calibrated sample, with `weight_cal`), `direct`
#'   (tibble `area_id, theta_direct, var_direct, n_eff, n_dwell, outcome,
#'   var_direct_smoothed`, one row per drawn municipality and outcome) and
#'   `smp_sub` (data frame `code, y, wt, n_hh, municipality, outcome`, one row
#'   per realised sub-area and outcome; `NULL` without a sub-area matrix)
#' @export
sim_draw <- function(rep, inputs) {
  rep <- .sim_check_rep(rep)
  smp <- draw_replicate(inputs$frame, inputs$design, seed = 100000L + rep,
                        segments = inputs$segments)
  smp <- apply_calibration(smp, inputs$pop_totals)
  direct <- smooth_variance_all(direct_estimate_all(smp, area_var = "municipality"),
                                method = "ghs")
  smp_sub <- NULL
  if (!is.null(inputs$pop_sub)) {
    parts <- lapply(sae_outcomes$outcome, function(o) {
      sa <- collapse_sample_to_subarea(smp, inputs$pop_sub, o)
      data.table::set(sa, j = "outcome", value = o)
      sa
    })
    smp_sub <- as.data.frame(data.table::rbindlist(parts, fill = TRUE))
  }
  list(sample = smp, direct = direct, smp_sub = smp_sub)
}

# First metric of `names` present in a fit's diagnostics; NA when absent. Fit
# functions name some metrics differently, hence the alternatives.
.sim_diag_metric <- function(dg, names) {
  if (is.null(dg) || !"metric" %in% names(dg)) return(NA_real_)
  hit <- intersect(names, dg$metric)
  if (!length(hit)) return(NA_real_)
  as.numeric(dg$value[match(hit[1], dg$metric)])
}

# The error message a fit function reports in its own diagnostics (metric
# `<name>_error_message`, text in `character_value`); NA when there is none.
.sim_own_error <- function(dg) {
  if (is.null(dg) || !all(c("metric", "character_value") %in% names(dg))) return(NA_character_)
  i <- which(grepl("_error_message$", dg$metric) & !is.na(dg$character_value))
  if (length(i)) as.character(dg$character_value[i[1]]) else NA_character_
}

#' Fit one method to one outcome of one replicate
#'
#' The in-sample municipalities are those with a direct estimate of `outcome`
#' in `direct`; every other municipality of the universe is predicted out of
#' sample. An area-level method takes the direct estimates and the area-level
#' covariates of the drawn municipalities (in sample) and of the others (out of
#' sample), with the arguments of [method_args()] and `settings[[method]]` laid
#' over them. A sub-area method takes the direct estimates of `outcome` (which
#' identify the drawn municipalities), the sub-area sample of `outcome`, the
#' sub-area matrix, its settings and the tuned hyperparameters of `outcome`.
#' The seed is [sae_fit_seed()] of the configuration, method, outcome and
#' replicate, unless `settings` give the method another; the `seed` column
#' records the seed the fit used.
#'
#' The fit runs through [safe_method_fit()]: a fit that fails is recorded as
#' NA rows for the in-sample municipalities with `fit_failed = TRUE` and its
#' message in `error_msg`. When the fit function reports its own error message
#' in its diagnostics (a povmap error inside a sub-area method, for example),
#' that message is recorded.
#'
#' @param method a method of `sae_methods_implemented`
#' @param outcome an outcome of `sae_outcomes`
#' @param rep replicate number
#' @param direct direct estimates of the replicate, long by outcome
#'   ([sim_draw()]`$direct`)
#' @param smp_sub sub-area sample of the replicate, long by outcome
#'   ([sim_draw()]`$smp_sub`); unused by area-level methods
#' @param inputs list from [sim_inputs()]
#' @param settings per-method settings, [sim_settings()] with the tuned XGB
#'   hyperparameters of `inputs$tune`
#' @return tibble with the columns `sim_fit_columns`: one row per municipality
#'   of the universe (or per in-sample municipality for a failed fit). Columns
#'   that do not apply to the method are NA.
#' @export
sim_fit_one <- function(method, outcome, rep, direct, smp_sub, inputs,
                        settings = sim_settings(inputs$tune)) {
  if (length(method) != 1L || length(outcome) != 1L) {
    .sim_stop("sim_fit_one fits one method to one outcome")
  }
  if (!method %in% sae_methods_implemented) .sim_stop("unknown method: ", method)
  if (!outcome %in% sae_outcomes$outcome) .sim_stop("unknown outcome: ", outcome)
  rep <- .sim_check_rep(rep)
  subarea <- method %in% sae_methods_subarea
  config <- .sim_config(method, inputs)
  seed <- sae_fit_seed(config, method, outcome, rep)
  direct <- as.data.frame(direct)
  drawn_ids <- sort(unique(direct$area_id[direct$outcome == outcome]))
  if (!length(drawn_ids)) .sim_stop("no direct estimates for outcome ", outcome)
  s <- settings[[method]]

  if (subarea) {
    smp <- as.data.frame(smp_sub)
    smp <- smp[smp$outcome == outcome, , drop = FALSE]
    if (!nrow(smp)) {
      .sim_stop("the sub-area sample of replicate ", rep, " holds no rows for outcome ", outcome)
    }
    args <- list(direct  = direct[direct$outcome == outcome, , drop = FALSE],
                 smp_sub = smp,
                 pop_sub = inputs$pop_sub,
                 outcome = outcome,
                 seed    = seed)
    args <- c(args, s[setdiff(names(s), "params")])
    # Hyperparameters are tuned per outcome.
    if ("params" %in% names(s)) {
      p <- s$params[[outcome]]
      if (is.null(p)) .sim_stop(.sim_xgb_message(config, outcome))
      args$params <- p
    }
  } else {
    X <- inputs$X_area
    aux_in  <- X[X$municipality %in% drawn_ids, , drop = FALSE]
    aux_oos <- X[X$municipality %in% setdiff(inputs$universe, drawn_ids), , drop = FALSE]
    args <- method_args(method, direct = direct, aux_in = aux_in, aux_oos = aux_oos,
                        W_matrix_full = inputs$W, outcome = outcome, seed = seed)
    if (length(s)) args <- utils::modifyList(args, s)
  }

  t0 <- Sys.time()
  est <- safe_method_fit(method, args, area_ids = drawn_ids, diagnostics = TRUE)
  wall <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  dg <- attr(est, "diagnostics")
  n <- nrow(est)
  err <- as.character(est$error_msg)
  if (any(est$fit_failed)) {
    own <- .sim_own_error(dg)
    if (!is.na(own)) err <- rep(substr(own, 1L, 200L), n)
    dg <- NULL
  }
  megb <- identical(method, "megb_sub")
  chr <- function(on, value) if (on) as.character(value) else NA_character_

  tibble::tibble(
    area_id        = as.integer(est$area_id),
    theta_hat      = as.numeric(est$theta_hat),
    mse            = as.numeric(est$mse),
    ci_lo          = as.numeric(est$ci_lo),
    ci_hi          = as.numeric(est$ci_hi),
    area_set       = as.character(est$area_set),
    fit_failed     = as.logical(est$fit_failed),
    error_msg      = err,
    config         = config,
    method         = method,
    replicate      = rep,
    outcome        = outcome,
    seed           = as.integer(.null_or(args[["seed"]], seed)),
    wall_secs      = wall,
    weight         = chr(subarea, .null_or(args$weight, "census")),
    n_divergent    = .sim_diag_metric(dg, c("n_divergent", "n_divergences")),
    max_rhat       = .sim_diag_metric(dg, c("max_rhat", "rhat_max")),
    # The area random-effect variance: sigma2_v of the Fay-Herriot fits, the
    # squared random-effect standard deviation of MEGB.
    sigma2_v       = .sim_diag_metric(dg, c("sigma2_v", "ran_eff_var")),
    povmap_sha     = chr(subarea, .povmap_sha()),
    mse_type       = chr(megb, .null_or(args$mse_type, "var")),
    pkg_versions   = sim_pkg_versions())
}

# Age, in seconds, after which a temporary file of a target is taken to be left
# by an interrupted process. A process that is still writing keeps its file's
# modification time recent, and processes that share an output folder (jobs of a
# scheduler on several machines) cannot see each other's process ids.
.sim_tmp_stale_secs <- 3600

# Write `path` through `writer(tmp)` under a temporary name, then rename it, so
# an interrupted write never leaves a file that reads as complete. Temporary
# files of the same target older than `.sim_tmp_stale_secs` are removed; younger
# ones may belong to another process still writing and are left alone.
.sim_write_atomic <- function(path, writer) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  tmp <- paste0(path, ".tmp.", Sys.getpid())
  on.exit(if (file.exists(tmp)) unlink(tmp), add = TRUE)
  writer(tmp)
  if (!file.rename(tmp, path)) .sim_stop("could not write ", path)
  prefix <- paste0(basename(path), ".tmp.")
  cand <- list.files(dirname(path))
  cand <- cand[startsWith(cand, prefix) &
                 grepl("^[0-9]+$", substring(cand, nchar(prefix) + 1L))]
  if (length(cand)) {
    full <- file.path(dirname(path), cand)
    age <- as.numeric(difftime(Sys.time(), file.mtime(full), units = "secs"))
    stale <- full[!is.na(age) & age > .sim_tmp_stale_secs]
    if (length(stale)) unlink(stale)
  }
  invisible(path)
}

.sim_write_parquet <- function(df, path) {
  .sim_write_atomic(path, function(tmp) arrow::write_parquet(tibble::as_tibble(df), tmp))
}

# The rows of a results file; NULL when it is missing, cannot be read or lacks
# one of `cols`.
.sim_read_rows <- function(path, cols) {
  if (!file.exists(path)) return(NULL)
  d <- tryCatch(as.data.frame(arrow::read_parquet(path)), error = function(e) NULL)
  if (is.null(d) || !all(cols %in% names(d))) return(NULL)
  d[, cols, drop = FALSE]
}

# Outcomes of `outcomes` with a failed fit in the per-fit rows `d`.
.sim_failed_outcomes <- function(d, outcomes) {
  if (is.null(d) || !nrow(d)) return(character(0))
  intersect(outcomes, unique(d$outcome[d$fit_failed %in% TRUE]))
}

# "fh_arcsin 0, xgb_sub 2", methods in alphabetical order.
.sim_failed_text <- function(failed) {
  f <- failed[order(names(failed))]
  paste0(names(f), " ", f, collapse = ", ")
}

# One replicate after the run checks: keep what is complete, draw, fit what is
# missing, write. Returns list(status, failed), `failed` the number of
# requested outcomes with a failed fit per method.
.sim_replicate <- function(rep, inputs, methods, outcomes, out_dir, settings,
                           retry_failed = FALSE) {
  t0 <- Sys.time()
  dpath <- .sim_direct_path(out_dir, rep)
  d_old <- .sim_read_rows(dpath, .sim_direct_columns)
  direct_ok <- !is.null(d_old) && all(sae_outcomes$outcome %in% d_old$outcome)

  fpaths <- vapply(methods, function(m) .sim_fit_path(out_dir, m, .sim_config(m, inputs), rep),
                   character(1))
  kept <- lapply(fpaths, function(p) {
    d <- .sim_read_rows(p, sim_fit_columns)
    if (is.null(d) || !nrow(d) || !all(d$replicate %in% rep)) return(NULL)
    if (retry_failed) {
      # Only the requested outcomes are refitted; failed rows of other outcomes
      # stay in the file.
      d <- d[!d$outcome %in% .sim_failed_outcomes(d, outcomes), , drop = FALSE]
    }
    d
  })
  todo <- lapply(methods, function(m) setdiff(outcomes, unique(kept[[m]]$outcome)))
  names(todo) <- methods
  n_fit <- sum(lengths(todo))
  count_failed <- function(rows) {
    vapply(methods, function(m) length(.sim_failed_outcomes(rows[[m]], outcomes)), integer(1))
  }
  if (direct_ok && n_fit == 0L) {
    failed <- count_failed(kept)
    message(sprintf("replicate %d: complete, skipped; failed fits: %s", rep,
                    .sim_failed_text(failed)))
    return(invisible(list(status = "skipped", failed = failed)))
  }

  dr <- sim_draw(rep, inputs)
  if (!direct_ok) .sim_write_parquet(dr$direct[, .sim_direct_columns], dpath)
  final <- kept
  for (m in methods[lengths(todo) > 0L]) {
    new <- lapply(todo[[m]], function(o) {
      sim_fit_one(m, o, rep, dr$direct, dr$smp_sub, inputs, settings)
    })
    rows <- as.data.frame(data.table::rbindlist(c(list(kept[[m]]), new), use.names = TRUE))
    rows <- rows[order(match(rows$outcome, sae_outcomes$outcome), method = "radix"), , drop = FALSE]
    .sim_write_parquet(rows, fpaths[[m]])
    final[[m]] <- rows
  }
  failed <- count_failed(final)
  message(sprintf("replicate %d: %d fit(s) in %.1f s; failed fits: %s", rep, n_fit,
                  as.numeric(difftime(Sys.time(), t0, units = "secs")),
                  .sim_failed_text(failed)))
  invisible(list(status = "done", failed = failed))
}

#' Run one replicate and write its files
#'
#' Writes `<out_dir>/direct/rep_NNNN.parquet` (direct estimates of all six
#' outcomes) and, for each method, `<out_dir>/fits/<method>__<config>/rep_NNNN.parquet`
#' (the rows of [sim_fit_one()] for every requested outcome, ordered by
#' outcome as in `sae_outcomes`). Each file is written under a temporary name
#' and renamed when complete; temporary files of the same target older than an
#' hour, left by an interrupted process, are removed then.
#'
#' Before the first draw the request is checked: methods, outcomes and
#' configurations are known, the inputs and the XGB hyperparameters the
#' methods need are present, CmdStan is located, and the packages each method
#' needs load. The Stan models of the Bayesian methods are compiled once, and
#' the run manifest `<out_dir>/run_manifest.json` is written, or compared with
#' the one already there: results computed from other inputs, settings or
#' package versions are refused, naming what differs.
#'
#' The replicate can be restarted: a direct file that reads and holds every
#' outcome is kept, and a per-fit file that reads, has every column and holds
#' rows of this replicate keeps its outcomes. Only the outcomes missing from a
#' method's file are fitted, and the file is rewritten with the kept and the
#' new rows. A replicate whose files are all complete is skipped without a draw.
#'
#' A failed fit is kept unless `retry_failed` is TRUE. A failure that comes from
#' the method (a singular fit, a bootstrap that does not converge) recurs under
#' the same seed, so refitting it repeats the failure; one that comes from the
#' environment (a process killed for memory, a file system error) does not.
#' With `retry_failed = TRUE` the requested outcomes whose fit failed are fitted
#' again and the file is rewritten; failed rows of outcomes not requested stay.
#'
#' @param rep replicate number, 1-9999
#' @param inputs list from [sim_inputs()]
#' @param methods methods to fit, from `sae_methods_implemented`
#' @param outcomes outcomes to fit, from `sae_outcomes$outcome`
#' @param out_dir output directory
#' @param settings per-method settings, see [sim_settings()]
#' @param retry_failed fit again the requested outcomes whose fit failed
#' @return list, invisibly: `status` (`"done"` or `"skipped"`) and `failed`,
#'   the number of requested outcomes with a failed fit per method
#' @export
sim_run_replicate <- function(rep, inputs, methods = sim_v3_methods,
                              outcomes = sae_outcomes$outcome, out_dir,
                              settings = sim_settings(inputs$tune), retry_failed = FALSE) {
  rep <- .sim_check_rep(rep)
  .sim_prepare(methods, outcomes, inputs, settings, out_dir)
  .sim_replicate(rep, inputs, methods, outcomes, out_dir, settings, isTRUE(retry_failed))
}

#' Run a range of replicates
#'
#' Runs each replicate of `reps` as [sim_run_replicate()] does, after one round
#' of checks, model compilation and run manifest for the whole range. With
#' `workers > 1` the replicates run in that many forked worker processes at a
#' time (not available on Windows, where they run one after the other); each
#' Bayes-Beta fit also runs its chains in parallel, so a worker may use as many
#' cores as chains. A replicate that stops with an error does not stop the
#' others; after all have run, the errors are raised together. Rerunning the
#' same call completes what is missing and skips the rest.
#'
#' @param reps replicate numbers, 1-9999
#' @inheritParams sim_run_replicate
#' @param workers number of worker processes
#' @param verbose print one line per replicate; FALSE also silences the
#'   messages of the fit functions
#' @return tibble, invisibly: `replicate`, `status` (`"done"` or `"skipped"`),
#'   `message`, and one column `failed_<method>` per method with the number of
#'   requested outcomes whose fit failed; stops when a replicate failed
#' @export
sim_run <- function(reps, inputs, methods = sim_v3_methods, outcomes = sae_outcomes$outcome,
                    out_dir, workers = 1L, settings = sim_settings(inputs$tune),
                    verbose = TRUE, retry_failed = FALSE) {
  if (!length(reps)) .sim_stop("no replicate requested")
  reps <- unique(vapply(reps, .sim_check_rep, integer(1)))
  workers <- suppressWarnings(as.integer(workers))
  if (length(workers) != 1L || is.na(workers) || workers < 1L) {
    .sim_stop("workers must be a positive whole number")
  }
  .sim_prepare(methods, outcomes, inputs, settings, out_dir)
  if (workers > 1L && .Platform$OS.type == "windows") {
    message("forked workers are not available on Windows; replicates run one after the other")
    workers <- 1L
  }
  retry_failed <- isTRUE(retry_failed)

  one <- function(r) {
    run <- function() .sim_replicate(r, inputs, methods, outcomes, out_dir, settings, retry_failed)
    res <- tryCatch(if (verbose) run() else suppressMessages(run()),
                    error = function(e) e)
    if (inherits(res, "error")) {
      list(status = "error", message = conditionMessage(res), failed = NULL)
    } else {
      list(status = res$status, message = NA_character_, failed = res$failed)
    }
  }
  res <- if (workers > 1L) {
    parallel::mclapply(reps, one, mc.cores = min(workers, length(reps)),
                       mc.preschedule = FALSE)
  } else {
    lapply(reps, one)
  }
  # A worker that ended without a result (killed, out of memory) returns NULL
  # or a try-error.
  res <- lapply(res, function(x) {
    if (is.list(x) && !is.null(x$status)) return(x)
    list(status = "error",
         message = if (inherits(x, "try-error")) as.character(x) else
           "the worker process ended without a result",
         failed = NULL)
  })
  st <- tibble::tibble(replicate = reps,
                       status = vapply(res, `[[`, character(1), "status"),
                       message = vapply(res, `[[`, character(1), "message"))
  for (m in sort(methods)) {
    st[[paste0("failed_", m)]] <- vapply(res, function(x) {
      if (is.null(x$failed)) NA_integer_ else unname(x$failed[[m]])
    }, integer(1))
  }
  bad <- st[st$status == "error", , drop = FALSE]
  if (nrow(bad)) {
    .sim_stop(nrow(bad), " of ", nrow(st), " replicate(s) failed:\n",
              paste0("replicate ", bad$replicate, ": ", bad$message, collapse = "\n"))
  }
  invisible(st)
}

# Command-line options of a script: each element of `argv` is `--name=value`
# for a name in `values`, or `--name` for a name in `flags`. Returns a named
# list, name -> value (NA for a flag); `stop_usage` reports what is malformed.
.sim_parse_options <- function(argv, values, flags, stop_usage) {
  if (!length(argv)) stop_usage("no arguments")
  is_flag <- argv %in% paste0("--", flags)
  is_value <- grepl("^--[a-z-]+=.+$", argv)
  bad <- argv[!is_flag & !is_value]
  if (length(bad)) stop_usage("not of the form --name=value: ", paste(bad, collapse = " "))
  nm <- ifelse(is_flag, sub("^--", "", argv), sub("^--([a-z-]+)=.*$", "\\1", argv))
  val <- ifelse(is_flag, NA_character_, sub("^--[a-z-]+=", "", argv))
  unk <- setdiff(nm, c(values, flags))
  if (length(unk)) stop_usage("unknown option(s): ", paste0("--", unk, collapse = " "))
  valued_flag <- nm[is_value][nm[is_value] %in% flags]
  if (length(valued_flag)) stop_usage("--", valued_flag[1], " takes no value")
  if (anyDuplicated(nm)) stop_usage("an option is given more than once")
  stats::setNames(as.list(val), nm)
}

# Replicate numbers from a range, a list or both ("1-100", "1,4,7-9").
.sim_parse_reps <- function(x, stop_usage) {
  parts <- trimws(strsplit(x, ",", fixed = TRUE)[[1]])
  ok <- grepl("^[0-9]+(-[0-9]+)?$", parts)
  if (!all(ok)) stop_usage("--reps: cannot read ", paste(parts[!ok], collapse = ", "))
  bounds <- lapply(parts, function(p) suppressWarnings(as.numeric(strsplit(p, "-", fixed = TRUE)[[1]])))
  # The bounds are checked before a range is expanded, so a large range never
  # builds a long vector.
  if (any(vapply(bounds, function(b) any(b < 1 | b > 9999), logical(1)))) {
    stop_usage("--reps: replicates run from 1 to 9999")
  }
  unique(unlist(lapply(seq_along(parts), function(i) {
    b <- as.integer(bounds[[i]])
    if (length(b) == 2L && b[2] < b[1]) stop_usage("--reps: decreasing range ", parts[i])
    seq.int(b[1], b[length(b)])
  })))
}

# Usage line of scripts/run_replicates.R.
.sim_usage <- paste(
  "Rscript scripts/run_replicates.R --reps=1-100 --out=DIR [--methods=...] [--outcomes=...]",
  "[--workers=N] [--area-matrix=PATH] [--subarea-matrix=PATH] [--extra-area=PATH]",
  "[--extra-subarea=PATH] [--xgb-params=PATH] [--retry-failed]")

#' Read the arguments of the replicate script
#'
#' Parses the command line of `scripts/run_replicates.R`: options of the form
#' `--name=value`, and the flag `--retry-failed`. `--reps` (a range, a list or
#' both: `"1-100"`, `"1,4,7-9"`) and `--out` are required; `--methods` and
#' `--outcomes` are comma-separated lists checked against the registries; a
#' path option must name an existing file.
#'
#' @param argv character vector, as `commandArgs(trailingOnly = TRUE)`
#' @return list: `reps` (integer), `out`, `methods`, `outcomes`, `workers`
#'   (integer), `area_matrix`, `subarea_matrix`, `extra_area`, `extra_subarea`,
#'   `xgb_params` (paths or `NULL`) and `retry_failed` (logical)
#' @export
sim_parse_args <- function(argv) {
  stop_usage <- function(...) .sim_stop(..., "\nusage: ", .sim_usage)
  values <- c("reps", "out", "methods", "outcomes", "workers", "area-matrix",
              "subarea-matrix", "extra-area", "extra-subarea", "xgb-params")
  flags <- "retry-failed"
  opt <- .sim_parse_options(argv, values, flags, stop_usage)
  nm <- names(opt)
  for (k in c("reps", "out")) if (is.null(opt[[k]])) stop_usage("--", k, " is required")
  lst <- function(x) trimws(strsplit(x, ",", fixed = TRUE)[[1]])
  reps <- .sim_parse_reps(opt$reps, stop_usage)

  methods <- if (is.null(opt$methods)) sim_v3_methods else lst(opt$methods)
  bad <- setdiff(methods, sae_methods_implemented)
  if (length(bad)) stop_usage("--methods: unknown ", paste(bad, collapse = ", "))
  outcomes <- if (is.null(opt$outcomes)) sae_outcomes$outcome else lst(opt$outcomes)
  bad <- setdiff(outcomes, sae_outcomes$outcome)
  if (length(bad)) stop_usage("--outcomes: unknown ", paste(bad, collapse = ", "))

  workers <- 1L
  if (!is.null(opt$workers)) {
    if (!grepl("^[0-9]+$", opt$workers) || as.integer(opt$workers) < 1L) {
      stop_usage("--workers must be a positive whole number")
    }
    workers <- as.integer(opt$workers)
  }
  path <- function(k) {
    x <- opt[[k]]
    if (!is.null(x) && !file.exists(x)) stop_usage("--", k, ": file not found: ", x)
    x
  }
  list(reps = reps, out = opt$out, methods = methods, outcomes = outcomes, workers = workers,
       area_matrix = path("area-matrix"), subarea_matrix = path("subarea-matrix"),
       extra_area = path("extra-area"), extra_subarea = path("extra-subarea"),
       xgb_params = path("xgb-params"), retry_failed = "retry-failed" %in% nm)
}
