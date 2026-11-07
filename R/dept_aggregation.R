# R/dept_aggregation.R
#
# Department-level agreement: municipal model estimates rolled up to the
# department with census population weights, compared with the department's
# design-based direct estimate and with the census truth.

#' Population-weighted rollup of municipal model estimates to Departamento
#'
#' dept is derived from area_id via the DIVIPOLA convention (area_id =
#' dept * 1000 + muni, so dept = area_id \%/\% 1000L).
#'
#' @param area_id integer vector of municipality codes
#' @param theta_hat numeric vector of municipal point estimates
#' @param n_persons numeric vector of static census population per
#'   municipality (same length as area_id/theta_hat)
#' @return tibble(dept, theta_dept, n_persons_dept, n_munis)
#' @export
aggregate_theta_to_dept <- function(area_id, theta_hat, n_persons) {
  stopifnot(length(area_id) == length(theta_hat),
            length(area_id) == length(n_persons))
  dt <- data.table::data.table(
    dept      = area_id %/% 1000L,
    theta_hat = theta_hat,
    n_persons = n_persons
  )
  out <- dt[, .(
    theta_dept     = sum(theta_hat * n_persons) / sum(n_persons),
    n_persons_dept = sum(n_persons),
    n_munis        = .N
  ), by = dept]
  tibble::as_tibble(out)
}

#' Union in-sample and out-of-sample per-rep municipal estimates
#'
#' Combines results_long-shaped in-sample rows with the matching
#' out-of-sample synthetic-prediction rows for the same rep, tagging each
#' with area_set. Errors if a (area_id, outcome, replicate) triple appears
#' in both inputs (would double-count a municipality's population weight in
#' the downstream dept rollup).
#'
#' @param in_sample,oos tibbles/data.frames with columns area_id, outcome,
#'   replicate, theta_hat (mse optional, carried through if present in both)
#' @return tibble with area_id, outcome, replicate, theta_hat (+ mse if
#'   present), area_set ("in_sample"/"out_of_sample")
#' @export
combine_in_sample_oos <- function(in_sample, oos) {
  in_sample <- as.data.frame(in_sample)
  oos       <- as.data.frame(oos)
  keep <- intersect(colnames(in_sample), colnames(oos))
  keep <- intersect(keep, c("area_id", "outcome", "replicate", "theta_hat",
                            "mse"))
  stopifnot(all(c("area_id", "outcome", "replicate", "theta_hat") %in% keep))
  a <- in_sample[, keep]; a$area_set <- "in_sample"
  b <- oos[, keep];       b$area_set <- "out_of_sample"
  combined <- rbind(a, b)
  key_cols <- c("area_id", "outcome", "replicate")
  dup <- duplicated(combined[, key_cols]) |
         duplicated(combined[, key_cols], fromLast = TRUE)
  if (any(dup)) {
    stop("combine_in_sample_oos: ", sum(dup), " overlapping rows between ",
         "in_sample and oos on (area_id, outcome, replicate)")
  }
  tibble::as_tibble(combined)
}

#' Departamento-level agreement summary: direct vs. population-weighted model
#'
#' @param dept_direct tibble(dept, outcome, replicate, estimate) —
#'   dept-level direct estimates, one row per rep (direct_two_levels_all(),
#'   level == "dept").
#' @param dept_model tibble(dept, outcome, replicate, theta_dept) —
#'   population-weighted municipal rollup, one row per rep
#'   (aggregate_theta_to_dept() applied per rep/outcome).
#' @param dept_truth tibble(dept, outcome, theta_true) — fixed census truth
#'   (no replicate; compute_truth_two_levels_all(), level == "dept").
#' @return tibble(dept, outcome, bias_direct, rmse_direct, bias_model,
#'   rmse_model, median_gap, n_reps) where median_gap = median(theta_dept -
#'   estimate) across reps (model minus direct).
#' @export
summarise_dept_comparison <- function(dept_direct, dept_model, dept_truth) {
  dept_direct <- as.data.frame(dept_direct)
  dept_model  <- as.data.frame(dept_model)
  dept_truth  <- as.data.frame(dept_truth)
  keys <- c("dept", "outcome")
  dd <- merge(dept_direct, dept_truth[, c(keys, "theta_true")], by = keys)
  dm <- merge(dept_model,  dept_truth[, c(keys, "theta_true")], by = keys)
  dd <- data.table::as.data.table(dd)
  dm <- data.table::as.data.table(dm)

  direct_summary <- dd[, .(
    bias_direct = mean(estimate - theta_true),
    rmse_direct = sqrt(mean((estimate - theta_true)^2)),
    n_reps      = .N
  ), by = keys]

  model_summary <- dm[, .(
    bias_model = mean(theta_dept - theta_true),
    rmse_model = sqrt(mean((theta_dept - theta_true)^2))
  ), by = keys]

  gap <- merge(
    dd[, c(keys, "replicate", "estimate"), with = FALSE],
    dm[, c(keys, "replicate", "theta_dept"), with = FALSE],
    by = c(keys, "replicate")
  )
  gap_summary <- gap[, .(
    median_gap = stats::median(theta_dept - estimate)
  ), by = keys]

  out <- merge(direct_summary, model_summary, by = keys)
  out <- merge(out, gap_summary, by = keys)
  tibble::as_tibble(out)
}
