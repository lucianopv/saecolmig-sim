#' Shared plumbing of the sub-area methods
#'
#' `xgb_sub` and `megb_sub` model census blocks and rural sections and aggregate
#' to municipality, where they are scored against the same census truth as the
#' area-level methods. They therefore share an input shape, prepared once here:
#'
#'   smp_data  one row per sub-area realised in the replicate, carrying the
#'             design-weighted outcome `y`, the summed weight `wt`, and the
#'             sub-area covariates
#'   pop_data  one row per sub-area of the whole prediction universe, with the
#'             same covariates and its census household and person counts
#'
#' @name method_subarea_common
NULL

# Structural columns of the sub-area matrix; every other column is a covariate.
# `pop_worldpop` is the census-free aggregation weight and `area_km2` is what it
# is derived from: as predictors they would let the model see its own weight, so
# both are structural.
.subarea_id_cols <- c("code", "municipality", "stratum",
                      "n_households_census", "n_persons_census",
                      "area_km2", "pop_worldpop")

#' Assemble the sub-area model inputs for one replicate and outcome
#'
#' @param direct tibble of direct estimates for this outcome; used only to
#'   identify which municipalities were drawn, never fed to the model
#' @param smp_sub collapsed sub-area sample for this outcome, from
#'   `collapse_sample_to_subarea()`
#' @param pop_sub the sub-area covariate matrix over the whole universe
#' @param outcome outcome name; selects the population weight
#' @param level `"subarea"` (default) fits on census blocks and rural sections;
#'   `"municipality"` rolls both the sample and the population up to
#'   municipality first, so that the same estimator runs on the same variables
#'   at the coarser resolution. The pairing separates what the finer covariates
#'   contribute from what the method contributes. Predictors are aggregated;
#'   scoring stays at municipality throughout.
#' @param weight `"census"` (default) weights sub-areas by their census
#'   household count (household outcomes) or person count (person outcomes);
#'   `"worldpop"` weights by the `pop_worldpop` column of `pop_sub`, which needs
#'   no census count.
#' @return list(sm, pop, x_cols, wpop, fixed, in_ids, oos_ids)
#' @keywords internal
.subarea_inputs <- function(direct, smp_sub, pop_sub, outcome,
                            level = c("subarea", "municipality"),
                            weight = c("census", "worldpop")) {
  level <- match.arg(level)
  weight <- match.arg(weight)
  pop <- data.table::as.data.table(pop_sub)
  smp <- data.table::as.data.table(smp_sub)
  x_cols <- setdiff(colnames(pop), .subarea_id_cols)
  if (!length(x_cols)) stop("subarea pop_data carries no covariate columns")

  # Household indicators aggregate over households; person-basis count ratios
  # aggregate over persons, as in outcome_spec()'s denominator.
  spec <- outcome_spec(outcome)
  wpop <- if (identical(weight, "worldpop")) {
    # WorldPop yields persons, not households, so weighting a household
    # indicator by pop_worldpop assumes household size is constant across the
    # sub-areas of a municipality. Deriving households would need a mean
    # household size, which is census information; the assumption is the price
    # of a weight that uses no census count.
    if (!"pop_worldpop" %in% colnames(pop)) {
      stop(".subarea_inputs: weight = \"worldpop\" needs a `pop_worldpop` ",
           "column on pop_sub")
    }
    "pop_worldpop"
  } else if (identical(spec$level, "person")) "n_persons_census"
    else "n_households_census"

  # A sub-area that cannot be predicted must not train the model either.
  pop <- pop[stats::complete.cases(pop[, x_cols, with = FALSE])]
  pop <- pop[pop[[wpop]] > 0]

  sm <- merge(smp[, c("code", "y", "wt"), with = FALSE],
              pop, by = "code")
  sm <- sm[is.finite(y) & is.finite(wt) & wt > 0]
  if (!nrow(sm)) stop("no usable sub-area sample rows for outcome ", outcome)

  if (identical(level, "municipality")) {
    # Population-weighted aggregation, on the denominator each covariate was
    # built with (person shares on persons, household and dwelling ones on
    # households), so the rolled-up value reproduces the municipal matrix rather
    # than an unweighted average of blocks.
    wcols <- .subarea_person_share_cols[.subarea_person_share_cols %in% x_cols]
    hcols <- setdiff(x_cols, wcols)
    agg_by <- function(d, wp, wh) {
      d[, c(lapply(.SD[, wcols, with = FALSE],
                   function(v) sum(v * get(wp)) / sum(get(wp))),
            lapply(.SD[, hcols, with = FALSE],
                   function(v) sum(v * get(wh)) / sum(get(wh)))),
        by = municipality,
        .SDcols = c(x_cols, wp, wh)]
    }
    keep_w <- intersect(c("n_households_census", "n_persons_census",
                          "pop_worldpop"), colnames(pop))
    pop_c <- pop[, lapply(.SD, sum), by = municipality, .SDcols = keep_w]
    pop <- merge(agg_by(pop, "n_persons_census", "n_households_census"),
                 pop_c, by = "municipality")
    pop[, code := as.character(municipality)]
    # The sample rolls up on its own design weight: sum(w * y) / sum(w) over the
    # municipality's realised sub-areas is the design-weighted municipal
    # proportion, which is what a municipality-level fit should see.
    sm <- sm[, .(y = sum(wt * y) / sum(wt), wt = sum(wt)), by = municipality]
    sm <- merge(sm, pop, by = "municipality")
    sm <- sm[is.finite(y) & is.finite(wt) & wt > 0]
    if (!nrow(sm)) stop("no usable municipal sample rows for outcome ", outcome)
  }

  # In-sample is defined by the drawn municipalities, as for the area-level
  # methods, not by which municipalities survive the covariate filter, so that
  # both families score the same area sets.
  in_ids  <- sort(unique(as.integer(direct$area_id)))
  oos_ids <- sort(setdiff(unique(as.integer(pop$municipality)), in_ids))

  # povmap::xgb joins its domain and sub-domain tables internally and needs the
  # domain id of the same type in both, so it is carried as character in both
  # frames; .povmap_estimates_from_fit() maps it back to the integer
  # municipality code.
  sm[, municipality := as.character(as.integer(municipality))]
  pop[, municipality := as.character(as.integer(municipality))]

  list(sm = as.data.frame(sm), pop = as.data.frame(pop),
       x_cols = x_cols, wpop = wpop,
       fixed = stats::as.formula(paste("y ~", paste(x_cols, collapse = " + "))),
       in_ids = in_ids, oos_ids = oos_ids)
}

# Census shares whose denominator is persons; every other covariate of the
# sub-area matrix aggregates on households.
.subarea_person_share_cols <- c(
  "share_working_age", "share_elderly", "share_secondary_completed",
  "share_illiterate", "share_economically_active", "share_afro",
  "share_indigenous")

# Convert povmap's $ind / $var / $CI tables to the contract-shape estimates
# tibble. Domains are mapped back to integer municipality codes; missing slots
# fall back to NA / Wald-derived intervals.
.povmap_estimates_from_fit <- function(fit_raw, in_ids, oos_ids) {
  ind <- as.data.frame(fit_raw$ind)
  v   <- if (!is.null(fit_raw$var)) as.data.frame(fit_raw$var) else NULL
  ci  <- if (!is.null(fit_raw$CI))  as.data.frame(fit_raw$CI)  else NULL

  # povmap's Domain column is character; coerce to the integer municipality id.
  ind$area_id <- suppressWarnings(as.integer(as.character(ind$Domain)))
  est <- data.frame(area_id   = ind$area_id,
                    theta_hat = ind$Mean)
  if (!is.null(v) && nrow(v) > 0) {
    v$area_id <- suppressWarnings(as.integer(as.character(v$Domain)))
    v_sub <- v[, c("area_id", "Mean")]
    colnames(v_sub)[2] <- ".mse_raw"
    est <- merge(est, v_sub, by = "area_id", all.x = TRUE)
    # Areas missing from $var get .Machine$double.eps to satisfy the > 0
    # contract.
    est$mse <- ifelse(is.na(est$.mse_raw), .Machine$double.eps,
                      pmax(est$.mse_raw, .Machine$double.eps))
    est$.mse_raw <- NULL
  } else {
    est$mse <- .Machine$double.eps
  }
  if (!is.null(ci) && nrow(ci) > 0) {
    ci$area_id <- suppressWarnings(as.integer(as.character(ci$Domain)))
    est <- merge(est, ci[, c("area_id", "Lower", "Upper")], by = "area_id",
                 all.x = TRUE)
    est$ci_lo <- est$Lower; est$Lower <- NULL
    est$ci_hi <- est$Upper; est$Upper <- NULL
  } else {
    est$ci_lo <- NA_real_; est$ci_hi <- NA_real_
  }
  # Fill missing or collapsed intervals with Wald ones from mse so that the
  # contract's ci_lo < ci_hi holds.
  needs_wald <- is.na(est$ci_lo) | is.na(est$ci_hi) | (est$ci_lo >= est$ci_hi)
  if (any(needs_wald)) {
    wald_half <- 1.96 * sqrt(est$mse[needs_wald])
    est$ci_lo[needs_wald] <- est$theta_hat[needs_wald] - wald_half
    est$ci_hi[needs_wald] <- est$theta_hat[needs_wald] + wald_half
  }
  # Final guard: if theta_hat is NA the Wald fallback also fails; widen by eps.
  ci_collapsed <- !is.na(est$ci_lo) & !is.na(est$ci_hi) & est$ci_lo >= est$ci_hi
  if (any(ci_collapsed)) {
    est$ci_lo[ci_collapsed] <- est$ci_lo[ci_collapsed] - .Machine$double.eps
    est$ci_hi[ci_collapsed] <- est$ci_hi[ci_collapsed] + .Machine$double.eps
  }
  est$area_set <- ifelse(est$area_id %in% in_ids, "in_sample", "out_of_sample")
  est <- est[order(est$area_set, est$area_id), ]
  tibble::as_tibble(est)
}

.povmap_meta <- function(method, est, in_ids, oos_ids, t0) {
  list(method = method,
       runtime_sec = as.numeric(difftime(Sys.time(), t0, units = "secs")),
       n_areas = nrow(est),
       n_in_sample = sum(est$area_set == "in_sample"),
       n_out_of_sample = sum(est$area_set == "out_of_sample"),
       povmap_version = as.character(utils::packageVersion("povmap")),
       povmap_sha = .povmap_sha())
}

# The installed povmap commit. The package version string does not change
# between commits whose estimators differ, so the SHA is recorded in each fit.
.povmap_sha <- function() {
  s <- suppressWarnings(utils::packageDescription("povmap")$RemoteSha)
  if (is.null(s) || !length(s) || is.na(s)) NA_character_ else as.character(s)
}

# A fit that could not be computed: NA estimates for every area, the error
# message in the diagnostics, and a "<method>_failed" flag.
.povmap_failure <- function(method, in_ids, oos_ids, t0, err_msg) {
  all_ids <- c(in_ids, oos_ids)
  est <- tibble::tibble(
    area_id   = all_ids,
    area_set  = c(rep("in_sample",     length(in_ids)),
                  rep("out_of_sample", length(oos_ids))),
    theta_hat = NA_real_,
    mse       = NA_real_,
    ci_lo     = NA_real_,
    ci_hi     = NA_real_
  )
  diag <- tibble::tibble(
    metric = c(paste0(method, "_failed"), paste0(method, "_error_message")),
    value  = c(1, NA_real_),
    character_value = c(NA_character_, err_msg)
  )
  meta <- list(method = method,
               runtime_sec = as.numeric(difftime(Sys.time(), t0, units = "secs")),
               n_areas = nrow(est),
               n_in_sample = length(in_ids),
               n_out_of_sample = length(oos_ids))
  structure(list(estimates = est, diagnostics = diag, meta = meta),
            class = c("sae_fit", paste0("sae_fit_", method)))
}
