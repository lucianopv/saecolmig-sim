#' Evaluation of per-fit results
#'
#' [sim_evaluate()] reads the direct and per-fit files a replicate run wrote
#' (layouts in `DATA_DICTIONARY.md`, section 4), optionally the estimates of a
#' method fitted elsewhere ([sim_read_external()]), scores every estimate
#' against the census truth and returns the uncertainty measures of
#' `R/uncertainty_eval.R` and the display summaries of `R/uncertainty_display.R`
#' as aggregate tables. `scripts/evaluate.R` writes them to parquet files.
#'
#' @name sim_evaluate_overview
NULL

#' Display labels of the methods
#'
#' Method code -> the label the evaluation tables carry in `method`, in the
#' order the comparison shows the methods. A method without an entry keeps its
#' code as its label.
#' @export
sim_method_labels <- c(fh_arcsin = "FH-arcsin", bayes_fh_beta = "Bayes-Beta",
                       megb_sub = "MEGB", xgb_sub = "XGB")

# Columns of the external per-fit layout.
.sim_external_columns <- c("replicate", "outcome", "municipality", "sampled", "estimate",
                           "estimated_variance", "lower", "upper", "truth")

# Columns of a per-fit file the evaluation reads; `max_rhat` is read when present.
.sim_eval_columns <- c("area_id", "outcome", "config", "method", "replicate", "area_set",
                       "theta_hat", "mse", "ci_lo", "ci_hi", "fit_failed")
.sim_eval_optional <- "max_rhat"

# Columns of the per-fit layout that are text or integer; the others are double
# except `fit_failed`.
.sim_fit_character <- c("area_set", "error_msg", "config", "method", "outcome", "weight",
                        "povmap_sha", "mse_type", "pkg_versions")
.sim_fit_integer <- c("area_id", "replicate", "seed")

#' Read an external per-fit file
#'
#' Estimates of a method fitted outside the replicate runner enter the
#' evaluation in the layout of `DATA_DICTIONARY.md`, section 4.3: one row per
#' replicate, outcome and municipality, with `replicate`, `outcome`,
#' `municipality`, `sampled`, `estimate`, `estimated_variance`, `lower`,
#' `upper` and `truth`. They are mapped to the per-fit layout of the runner
#' (`sim_fit_columns`): `area_id` from `municipality` (a zero-padded text code
#' such as `"05001"` reads as 5001), `area_set` from `sampled`, `theta_hat`
#' from `estimate`, `mse` from `estimated_variance`, `ci_lo` and `ci_hi` from
#' `lower` and `upper`, and `fit_failed = FALSE`; columns the layout has and the
#' file does not are missing.
#'
#' The file is refused, with a message that says why, when a column is absent,
#' a value is missing, an outcome is not one of `sae_outcomes`, `sampled` is
#' not logical (or 0/1), a (`outcome`, `replicate`, `municipality`) row appears
#' more than once, a variance is negative, a lower bound lies above its upper
#' bound, or (with `universe`) a municipality code is not in the universe.
#'
#' @param x path to a parquet or CSV file, or a data frame
#' @param method the method name its rows carry in `method`
#' @param config the configuration name its rows carry in `config`
#' @param universe municipality codes of the universe; `NULL` skips the check of
#'   the codes
#' @return tibble with the columns `sim_fit_columns` followed by `truth`, the
#'   file's census value of each row
#' @export
sim_read_external <- function(x, method, config = "external", universe = NULL) {
  got <- .sim_read_input(x, "external per-fit file")
  where <- .sim_where("external per-fit file", got$source)
  e <- got$data
  absent <- setdiff(.sim_external_columns, names(e))
  if (length(absent)) .sim_stop(where, ": missing column(s): ", paste(absent, collapse = ", "))
  has_na <- Reduce(`|`, lapply(e[.sim_external_columns], is.na))
  if (any(has_na)) {
    n <- sum(has_na)
    .sim_stop(where, ": ", .sim_np(n, "row has", "rows have"), " a missing value")
  }
  area_id <- .sim_as_integer_key(e$municipality, "municipality", where)
  replicate <- .sim_as_integer_key(e$replicate, "replicate", where)
  unk <- setdiff(unique(as.character(e$outcome)), sae_outcomes$outcome)
  if (length(unk)) {
    .sim_stop(where, ": unknown outcome(s): ", paste(unk, collapse = ", "),
              "; outcomes are those of sae_outcomes")
  }
  sampled <- e$sampled
  if (is.numeric(sampled) && all(sampled %in% c(0, 1))) sampled <- sampled == 1
  if (!is.logical(sampled)) {
    .sim_stop(where, ": `sampled` must be TRUE or FALSE (or 1 or 0)")
  }
  num <- c("estimate", "estimated_variance", "lower", "upper", "truth")
  non_num <- num[!vapply(e[num], is.numeric, logical(1))]
  if (length(non_num)) {
    .sim_stop(where, ": column(s) not numeric: ", paste(non_num, collapse = ", "))
  }
  n_dup <- sum(duplicated(data.frame(o = e$outcome, r = replicate, a = area_id)))
  if (n_dup) {
    .sim_stop(where, ": ", .sim_np(n_dup, "duplicate row", "duplicate rows"),
              " of (outcome, replicate, municipality)")
  }
  n_neg <- sum(e$estimated_variance < 0)
  if (n_neg) .sim_stop(where, ": ", .sim_np(n_neg, "row has", "rows have"), " a negative variance")
  n_inv <- sum(e$lower > e$upper)
  if (n_inv) {
    .sim_stop(where, ": ", .sim_np(n_inv, "row has its", "rows have their"),
              " lower bound above the upper bound")
  }
  if (!is.null(universe)) {
    unknown <- sort(setdiff(unique(area_id), universe))
    if (length(unknown)) {
      .sim_stop(where, ": ", .sim_np(length(unknown), "municipality code is", "municipality codes are"),
                " not in the universe", .sim_first(unknown))
    }
  }
  n <- nrow(e)
  out <- lapply(sim_fit_columns, function(k) {
    if (k %in% .sim_fit_character) rep(NA_character_, n)
    else if (k %in% .sim_fit_integer) rep(NA_integer_, n)
    else if (k == "fit_failed") rep(FALSE, n)
    else rep(NA_real_, n)
  })
  names(out) <- sim_fit_columns
  out <- tibble::as_tibble(out)
  out$area_id <- area_id
  out$theta_hat <- as.numeric(e$estimate)
  out$mse <- as.numeric(e$estimated_variance)
  out$ci_lo <- as.numeric(e$lower)
  out$ci_hi <- as.numeric(e$upper)
  out$area_set <- ifelse(sampled, "in_sample", "out_of_sample")
  out$config <- rep(as.character(config), n)
  out$method <- rep(as.character(method), n)
  out$replicate <- replicate
  out$outcome <- as.character(e$outcome)
  out$truth <- as.numeric(e$truth)
  out
}

# "3" or "1-3" or "1-3, 7".
.sim_rep_text <- function(r) {
  r <- sort(unique(r))
  br <- cumsum(c(1L, diff(r) != 1L))
  paste(vapply(split(r, br), function(x) {
    if (length(x) == 1L) as.character(x) else paste0(x[1], "-", x[length(x)])
  }, character(1)), collapse = ", ")
}

# The census truth and its universe: a truth table without missing values or
# repeated area-outcome pairs, holding every outcome for every area it names.
.sim_eval_truth <- function(truth) {
  tr <- .sim_read_truth(truth, integer(0))
  universe <- sort(unique(tr$area_id))
  tr <- .sim_read_truth(tr, universe)
  list(truth = data.table::as.data.table(tr)[, c("area_id", "outcome", "theta_true"), with = FALSE],
       universe = universe)
}

# The run's methods: one row per `fits/<method>__<config>` folder, restricted to
# `methods` and in their order when given, otherwise in the order of
# sim_method_labels and then by name. A method fitted on two configurations in
# one folder is refused: every summary reads one configuration per method.
.sim_eval_methods <- function(run_dir, methods) {
  dirs <- list.dirs(file.path(run_dir, "fits"), full.names = FALSE, recursive = FALSE)
  dirs <- dirs[grepl("^.+__.+$", dirs)]
  found <- data.frame(dir = dirs, method = sub("__.*$", "", dirs),
                      config = sub("^.*?__", "", dirs, perl = TRUE), stringsAsFactors = FALSE)
  if (!is.null(methods)) {
    miss <- setdiff(methods, found$method)
    if (length(miss)) {
      .sim_stop(run_dir, ": no per-fit results for method ", paste(miss, collapse = ", "),
                " (a folder fits/<method>__<config>)")
    }
    found <- found[found$method %in% methods, , drop = FALSE]
    found <- found[order(match(found$method, methods)), , drop = FALSE]
  } else {
    found <- found[order(match(found$method, names(sim_method_labels)), found$method), , drop = FALSE]
  }
  twice <- unique(found$method[duplicated(found$method)])
  if (length(twice)) {
    .sim_stop(run_dir, ": method(s) ", paste(twice, collapse = ", "), " fitted on more than one ",
              "configuration; evaluate one configuration per folder")
  }
  found
}

# Per-fit rows of one method folder for `reps`, columns `.sim_eval_columns` (and
# `max_rhat` when present). Returns NULL with the replicates that have no file
# in attr "missing".
.sim_eval_read_fits <- function(run_dir, dir, method, config, reps) {
  paths <- file.path(run_dir, "fits", dir, sprintf("rep_%04d.parquet", reps))
  miss <- reps[!file.exists(paths)]
  if (length(miss)) return(structure(list(), missing = miss))
  rows <- lapply(seq_along(reps), function(i) {
    d <- data.table::as.data.table(arrow::read_parquet(paths[i]))
    absent <- setdiff(.sim_eval_columns, names(d))
    if (length(absent)) {
      .sim_stop(paths[i], ": missing column(s): ", paste(absent, collapse = ", "))
    }
    d <- d[, intersect(c(.sim_eval_columns, .sim_eval_optional), names(d)), with = FALSE]
    if (any(d$replicate != reps[i]) || any(d$method != method) || any(d$config != config)) {
      .sim_stop(paths[i], ": holds rows of another replicate, method or configuration")
    }
    d
  })
  data.table::rbindlist(rows, use.names = TRUE, fill = TRUE)
}

# Each fit of `d` (rows of one source) is present for every outcome and
# replicate, and every fit that did not fail has one row per municipality of the
# universe. A failed fit holds rows of the drawn municipalities only; it is
# dropped whole later.
.sim_eval_check_fits <- function(d, where, outcomes, reps, universe) {
  want <- data.table::CJ(outcome = outcomes, replicate = reps)
  have <- unique(d[, c("outcome", "replicate"), with = FALSE])
  miss <- want[!have, on = c("outcome", "replicate")]
  if (nrow(miss)) {
    .sim_stop(where, ": no fit for ", nrow(miss), " of the (outcome, replicate) combinations ",
              "evaluated (first: ", miss$outcome[1], ", replicate ", miss$replicate[1], ")")
  }
  n_u <- length(universe)
  per <- d[, .(n = .N, n_id = data.table::uniqueN(area_id), known = all(area_id %in% universe),
               failed = any(fit_failed %in% TRUE)), by = c("outcome", "replicate")]
  bad <- per[!failed & (n != n_u | n_id != n_u | !known)]
  if (nrow(bad)) {
    .sim_stop(where, ": ", nrow(bad), " fit(s) without one row per municipality of the ",
              "universe (first: ", bad$outcome[1], ", replicate ", bad$replicate[1], ", ",
              bad$n[1], " rows for ", n_u, " municipalities)")
  }
  invisible(TRUE)
}

# Every method sees the same drawn municipalities in each (outcome, replicate).
# The comparison of sets reads in-sample rows only, so a fit with no drawn
# municipality would leave no set to compare; it is refused first.
.sim_eval_check_drawn <- function(d) {
  fits <- unique(d[, c("method", "outcome", "replicate"), with = FALSE])
  drawn <- unique(d[area_set == "in_sample", c("method", "outcome", "replicate"), with = FALSE])
  none <- fits[!drawn, on = c("method", "outcome", "replicate")]
  if (nrow(none)) {
    .sim_stop(nrow(none), " fit(s) have no drawn municipality (first: ", none$method[1], ", ",
              none$outcome[1], ", replicate ", none$replicate[1], ")")
  }
  ok <- tryCatch(ur_assert_same_sets(d), error = function(e) e)
  if (!inherits(ok, "error")) return(invisible(TRUE))
  fp <- d[area_set == "in_sample", .(fp = paste(sort(unique(area_id)), collapse = ",")),
          by = c("method", "outcome", "replicate")]
  cells <- fp[, .(n = data.table::uniqueN(fp)), by = c("outcome", "replicate")][n > 1L]
  first <- fp[cells[1], on = c("outcome", "replicate")]
  groups <- vapply(split(first$method, first$fp), paste, character(1), collapse = ", ")
  .sim_stop("the drawn municipalities differ across methods in ", nrow(cells),
            " (outcome, replicate) combination(s) (first: ", cells$outcome[1], ", replicate ",
            cells$replicate[1], ": ", paste0("{", groups, "}", collapse = " vs "), ")")
}

# What happened to the fits of each method and outcome: fits read, dropped whole
# (failed, or split R-hat above the threshold), left without a finite row, and
# kept, with the rows dropped for a non-finite value. `d` holds the rows before
# ur_drop_bad_fits(), `cl` its result; rows follow the method order `ord` and
# the outcome order `outcomes`.
.sim_eval_fit_check <- function(d, cl, ord, outcomes) {
  gk <- c("method", "config", "outcome")
  fk <- c(gk, "replicate")
  read <- d[, .(fits_read = data.table::uniqueN(replicate), rows_read = .N), by = gk]
  drp <- cl$dropped[, .(failed = sum(reason == "failed"), dropped_rhat = sum(reason == "rhat")),
                    by = gk]
  in_drp <- d[cl$dropped, on = fk, nomatch = 0L][, .(rows_dropped = .N), by = gk]
  kept <- cl$d[, .(kept = data.table::uniqueN(replicate), rows_kept = .N), by = gk]
  fc <- Reduce(function(x, y) merge(x, y, by = gk, all.x = TRUE), list(read, drp, in_drp, kept))
  for (k in c("failed", "dropped_rhat", "rows_dropped", "kept", "rows_kept")) {
    data.table::set(fc, i = which(is.na(fc[[k]])), j = k, value = 0L)
  }
  fc[, `:=`(dropped_nonfinite = fits_read - failed - dropped_rhat - kept,
            rows_nonfinite = rows_read - rows_dropped - rows_kept)]
  data.table::setorderv(fc[, `:=`(o1 = match(method, ord), o2 = match(outcome, outcomes))],
                        c("o1", "o2"))
  fc[, c(gk, "fits_read", "failed", "dropped_rhat", "dropped_nonfinite", "kept",
         "rows_nonfinite"), with = FALSE]
}

# Stop when a method keeps no fit of an outcome, saying why; warn, once, when
# the methods keep different numbers of fits of an outcome, because they are
# then compared over different replicates.
.sim_eval_check_kept <- function(fc) {
  none <- fc[kept == 0L]
  if (nrow(none)) {
    why <- sprintf("%s keeps no fit of %s (of %d fits: %d failed, %d with a split R-hat above 1.05, %d without a finite row)",
                   none$method, none$outcome, none$fits_read, none$failed, none$dropped_rhat,
                   none$dropped_nonfinite)
    .sim_stop(paste(why, collapse = "; "))
  }
  uneven <- fc[, .(n = data.table::uniqueN(kept)), by = outcome][n > 1L]$outcome
  if (length(uneven)) {
    txt <- vapply(uneven, function(o) {
      x <- fc[outcome == o]
      paste0(o, " (", paste(x$method, x$kept, collapse = ", "), ")")
    }, character(1))
    warning("methods keep different numbers of fits of the same outcome and are compared ",
            "over different replicates: ", paste(txt, collapse = "; "), call. = FALSE)
  }
  invisible(TRUE)
}

# Dwellings drawn per municipality, outcome and replicate, from the direct files.
.sim_eval_read_direct <- function(run_dir, reps) {
  paths <- file.path(run_dir, "direct", sprintf("rep_%04d.parquet", reps))
  miss <- reps[!file.exists(paths)]
  if (length(miss)) {
    .sim_stop(file.path(run_dir, "direct"), ": no file for replicate ", .sim_rep_text(miss))
  }
  # Parquet files written with R attributes that arrow cannot restore warn
  # "Invalid metadata$r"; that warning carries no information about the values
  # and is muffled, any other surfaces.
  quiet <- function(expr) withCallingHandlers(expr, warning = function(w) {
    if (grepl("Invalid metadata\\$r", conditionMessage(w))) invokeRestart("muffleWarning")
  })
  data.table::rbindlist(lapply(seq_along(reps), function(i) {
    d <- quiet(data.table::as.data.table(
      arrow::read_parquet(paths[i], col_select = c("area_id", "outcome", "n_dwell"))))
    d[, replicate := reps[i]]
  }))
}

#' Evaluate the per-fit results of a replicate run
#'
#' Scores the estimates of each method against the census truth over the
#' replicates `reps`, in sample (municipalities drawn in a replicate) and out of
#' sample (the others), and returns the aggregate tables of
#' `DATA_DICTIONARY.md`, section 4.4.
#'
#' The steps:
#'
#' 1. Read the truth; its municipalities are the universe. Read the per-fit
#'    rows of each method of `run_dir/fits/` (all methods, or those named in
#'    `methods`) and, when given, the external per-fit file, whose rows take the
#'    method label `external_label` and the configuration `"external"`.
#' 2. Gates, each stopping with a message: the external truth equals the truth
#'    file within 1e-9; the external replicates are `reps`; every method has a
#'    fit of every outcome evaluated in every replicate, and every fit that did
#'    not fail has one row per municipality of the universe (a failed fit holds
#'    the rows of the drawn municipalities and is dropped in step 3); every fit
#'    has at least one drawn municipality, and every method sees the same drawn
#'    municipalities in each outcome and replicate ([ur_assert_same_sets()]).
#' 3. Drop failed and non-converged fits whole, and rows with a non-finite
#'    estimate, MSE or bound ([ur_drop_bad_fits()]); label the methods
#'    ([sim_method_labels]). The table `fit_check` counts, per method and
#'    outcome, the fits read, dropped and kept. The evaluation stops when a
#'    method keeps no fit of an outcome, saying why, and warns when the methods
#'    keep different numbers of fits of an outcome, since they are then
#'    compared over different replicates.
#' 4. Join the truth and the dwellings drawn (`n_dwell`) of each drawn
#'    municipality from `run_dir/direct/`; out of sample `n_dwell` is 0. Bin the
#'    rows by sample size (eight quantile bins of `n_dwell` over the direct rows
#'    of `reps`, the out-of-sample rows as their own first bin) and by the
#'    quintile of the true value within its outcome.
#' 5. Compute the row metrics and the tables: per municipality
#'    ([ur_area_summary()]), per set ([ur_set_summary()]), by bin
#'    ([ur_by_bin()]), the publication flags ([ur_cv_flags()]), the Monte Carlo
#'    standard errors ([ur_mcse()]) and the reliability curve
#'    ([ur_reliability()]); then, for every outcome, the display summaries
#'    ([ur_ratio_summary()], [ur_tracking_summary()], [ur_size_medians()],
#'    [ur_coverage_summary()]). The sample-size groups of [ur_size_medians()]
#'    take the sizes of the first method in the order of [sim_method_labels],
#'    so they do not depend on the order in which methods are named.
#'
#' Each method keeps the configuration it was fitted on in `config`; `column`
#' names the covariate set the configurations share.
#'
#' @param run_dir output directory of a replicate run
#' @param truth census truth, a path to `truth_full.parquet` (or CSV) or a data
#'   frame with `area_id`, `theta_true`, `n_total`, `outcome`
#' @param reps replicates to evaluate
#' @param methods method codes of the run to evaluate; `NULL` evaluates every
#'   method with results in `run_dir/fits/`
#' @param external `NULL`, or estimates of a method fitted elsewhere in the
#'   layout of [sim_read_external()]: a path or a data frame
#' @param external_label method label of the external rows; no method of the
#'   run may carry it
#' @param column label of the covariate set, written to every table
#' @return named list of data.tables: `area_summary`, `set_summary`,
#'   `import_check`, `fit_check`, `by_bin`, `cv_flags`, `mcse`, `reliability`,
#'   `ratio_summary`, `tracking_summary`, `size_medians`, `coverage_summary`
#' @export
sim_evaluate <- function(run_dir, truth, reps = 1:100, methods = NULL, external = NULL,
                         external_label = "MEGB", column = "geo") {
  if (!dir.exists(run_dir)) .sim_stop(run_dir, ": folder not found")
  if (!length(reps)) .sim_stop("no replicate requested")
  reps <- sort(unique(vapply(reps, .sim_check_rep, integer(1))))
  tt <- .sim_eval_truth(truth)
  TR <- tt$truth; universe <- tt$universe

  ## ---- per-fit rows of the run and of the external file --------------------
  runs <- .sim_eval_methods(run_dir, methods)
  if (!nrow(runs) && is.null(external)) .sim_stop(run_dir, ": no per-fit results in fits/")
  labels <- ifelse(runs$method %in% names(sim_method_labels),
                   unname(sim_method_labels[runs$method]), runs$method)
  if (anyDuplicated(labels)) {
    .sim_stop("methods ", paste(runs$method, collapse = ", "), " share a display label")
  }
  arms <- list(); where <- character(0)
  missing_files <- character(0)
  for (i in seq_len(nrow(runs))) {
    d <- .sim_eval_read_fits(run_dir, runs$dir[i], runs$method[i], runs$config[i], reps)
    miss <- attr(d, "missing")
    if (length(miss)) {
      missing_files <- c(missing_files, paste0("fits/", runs$dir[i], ": no file for replicate ",
                                               .sim_rep_text(miss)))
      next
    }
    arms[[labels[i]]] <- d
    where[[labels[i]]] <- paste0("fits/", runs$dir[i])
  }
  if (length(missing_files)) .sim_stop(run_dir, ": ", paste(missing_files, collapse = "; "))

  ext_source <- NA_character_; tdiff <- NA_real_
  if (!is.null(external)) {
    hit <- runs$method[labels == external_label | runs$method == external_label]
    if (length(hit)) {
      .sim_stop("the external label ", external_label, " is also run method ",
                paste(hit, collapse = ", "), "; leave it out with `methods` or choose another ",
                "external_label")
    }
    ext_source <- if (is.data.frame(external)) "<data frame>" else basename(external)
    e <- data.table::as.data.table(sim_read_external(external, external_label,
                                                      universe = universe))
    tdiff <- merge(e[, c("area_id", "outcome", "truth"), with = FALSE], TR,
                   by = c("area_id", "outcome"))[, suppressWarnings(max(abs(truth - theta_true)))]
    if (!is.finite(tdiff) || tdiff > 1e-9) {
      .sim_stop("the truth of the external per-fit file differs from the truth file by up to ",
                signif(tdiff, 3))
    }
    if (!setequal(unique(e$replicate), reps)) {
      .sim_stop("the external per-fit file holds replicates ", .sim_rep_text(e$replicate),
                "; the evaluated replicates are ", .sim_rep_text(reps))
    }
    keep <- intersect(c(.sim_eval_columns, .sim_eval_optional), names(e))
    arms[[external_label]] <- e[, keep, with = FALSE]
    where[[external_label]] <- "external per-fit file"
  }

  ## ---- gates on the raw rows -----------------------------------------------
  outcomes <- unique(unlist(lapply(arms, function(d) unique(d$outcome))))
  no_truth <- setdiff(outcomes, TR$outcome)
  if (length(no_truth)) .sim_stop("no truth for outcome(s) ", paste(no_truth, collapse = ", "))
  outcomes <- c(intersect(sae_outcomes$outcome, outcomes), setdiff(outcomes, sae_outcomes$outcome))
  for (m in names(arms)) .sim_eval_check_fits(arms[[m]], where[[m]], outcomes, reps, universe)
  # Rows carry the method label from here on; `ord` orders the labels as
  # sim_method_labels, then the others in the order read.
  ord <- c(intersect(unname(sim_method_labels), names(arms)),
           setdiff(names(arms), sim_method_labels))
  all_d <- data.table::rbindlist(arms, use.names = TRUE, fill = TRUE, idcol = ".label")
  all_d[, method := .label][, .label := NULL]
  rm(arms)
  .sim_eval_check_drawn(all_d)

  ## ---- drop bad fits, join truth, dwellings drawn and bins -------------------
  cl <- ur_drop_bad_fits(all_d)
  fit_check <- .sim_eval_fit_check(all_d, cl, ord, outcomes)
  rm(all_d)
  .sim_eval_check_kept(fit_check)
  d <- cl$d
  DW <- .sim_eval_read_direct(run_dir, reps)
  # One set of size breaks for every outcome: the dwellings drawn in a
  # municipality do not depend on the outcome, so the bins are the same
  # everywhere.
  size_br <- ur_size_breaks(DW$n_dwell, n_bins = 8L)
  PREV <- TR[, .(area_id, prev_bin = cut(theta_true, unique(stats::quantile(theta_true, 0:5 / 5)),
                                         include.lowest = TRUE, labels = FALSE)), by = outcome]
  d <- merge(d, TR, by = c("area_id", "outcome"))
  d <- merge(d, DW, by = c("area_id", "outcome", "replicate"), all.x = TRUE)
  d[area_set == "out_of_sample", n_dwell := 0L]
  if (d[area_set == "in_sample", anyNA(n_dwell)]) {
    .sim_stop(run_dir, ": in-sample rows without a direct estimate in direct/ (no n_dwell)")
  }
  d[, size_bin := ur_size_bin(n_dwell, area_set, size_br)]
  d <- merge(d, PREV, by = c("area_id", "outcome"))
  ur_add_row_metrics(d)

  ## ---- tables ----------------------------------------------------------------
  area <- ur_area_summary(d)
  sets <- ur_set_summary(d, area)
  # A size bin lies within one area set (the out-of-sample rows are a bin of
  # their own), so grouping by the set as well leaves the size bins unchanged.
  by_size <- ur_by_bin(d, "size_bin", by = UR_KEYS, x = "n_dwell")
  by_size[, `:=`(bin_var = "size", bin = as.character(size_bin),
                 bin_order = as.integer(size_bin), size_bin = NULL)]
  by_prev <- ur_by_bin(d, "prev_bin", by = UR_KEYS)
  by_prev[, `:=`(bin_var = "prevalence", bin = as.character(prev_bin),
                 bin_order = as.integer(prev_bin), x_med = NA_real_, prev_bin = NULL)]
  by_bin <- data.table::rbindlist(list(by_size, by_prev), use.names = TRUE)
  data.table::setcolorder(by_bin, c(UR_KEYS, "bin_var", "bin", "bin_order"))
  check <- data.table::data.table(
    source = ext_source,
    rows = if (is.null(external)) NA_integer_ else d[method == external_label, .N],
    truth_max_abs_diff = tdiff,
    rows_kept = nrow(d),
    fits_dropped = nrow(cl$dropped),
    fits_dropped_failed = sum(cl$dropped$reason == "failed"),
    fits_dropped_rhat = sum(cl$dropped$reason == "rhat"),
    rows_nonfinite = cl$rows_nonfinite)

  # The display summaries, for every outcome, with `method` ordered as
  # sim_method_labels so that the reference method of the size groups is fixed.
  a_disp <- data.table::copy(area)[, method := factor(method, levels = ord)]
  per_outcome <- function(f) {
    data.table::rbindlist(lapply(outcomes, function(o) {
      x <- f(a_disp, o)
      x[, `:=`(method = as.character(method), outcome = o)]
      data.table::setcolorder(x, c("method", "outcome"))
    }))
  }
  out <- list(
    area_summary = area,
    set_summary = sets,
    import_check = check,
    fit_check = fit_check,
    by_bin = by_bin,
    cv_flags = ur_cv_flags(d, area),
    mcse = ur_mcse(ur_rep_sums(d)),
    reliability = ur_reliability(d),
    ratio_summary = per_outcome(ur_ratio_summary),
    tracking_summary = per_outcome(ur_tracking_summary),
    size_medians = per_outcome(ur_size_medians),
    coverage_summary = per_outcome(ur_coverage_summary))
  lapply(out, function(x) data.table::set(x, j = "column", value = column)[])
}

# Usage line of scripts/evaluate.R.
.sim_eval_usage <- paste(
  "Rscript scripts/evaluate.R --run=DIR --truth=PATH --reps=1-100 --out=DIR",
  "[--external=PATH] [--external-label=MEGB] [--methods=...] [--column=geo]")

#' Read the arguments of the evaluation script
#'
#' Parses the command line of `scripts/evaluate.R`: options of the form
#' `--name=value`. `--run` (a folder), `--truth` (a file), `--reps` (a range, a
#' list or both: `"1-100"`, `"1,4,7-9"`) and `--out` are required;
#' `--external` names a file, `--external-label` its method label (default
#' `MEGB`; refused without `--external`), `--methods` a comma-separated list of
#' the run's methods and `--column` the covariate-set label (default `geo`).
#'
#' @param argv character vector, as `commandArgs(trailingOnly = TRUE)`
#' @return list: `run`, `truth`, `reps` (integer), `out`, `external` (path or
#'   `NULL`), `external_label`, `methods` (`NULL` for every method of the run)
#'   and `column`
#' @export
sim_parse_eval_args <- function(argv) {
  stop_usage <- function(...) .sim_stop(..., "\nusage: ", .sim_eval_usage)
  opt <- .sim_parse_options(argv, c("run", "truth", "reps", "out", "external", "external-label",
                                    "methods", "column"), character(0), stop_usage)
  # Options are read by exact name: `$` would match a prefix of a longer name
  # (`external` of `external-label`).
  opt_of <- function(k) opt[[k, exact = TRUE]]
  for (k in c("run", "truth", "reps", "out")) if (is.null(opt_of(k))) stop_usage("--", k, " is required")
  if (!dir.exists(opt_of("run"))) stop_usage("--run: folder not found: ", opt_of("run"))
  for (k in c("truth", "external")) {
    if (!is.null(opt_of(k)) && !file.exists(opt_of(k))) stop_usage("--", k, ": file not found: ", opt_of(k))
  }
  if (!is.null(opt_of("external-label")) && is.null(opt_of("external"))) {
    stop_usage("--external-label names the method of --external, which is not given")
  }
  lst <- function(x) trimws(strsplit(x, ",", fixed = TRUE)[[1]])
  list(run = opt_of("run"), truth = opt_of("truth"), reps = .sim_parse_reps(opt_of("reps"), stop_usage),
       out = opt_of("out"), external = opt_of("external"),
       external_label = if (is.null(opt_of("external-label"))) "MEGB" else opt_of("external-label"),
       methods = if (is.null(opt_of("methods"))) NULL else lst(opt_of("methods")),
       column = if (is.null(opt_of("column"))) "geo" else opt_of("column"))
}
