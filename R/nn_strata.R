#' Standardise feature columns to z-scores
#'
#' @param df data.frame containing an identifier column and feature columns.
#' @param id_col character. Name of the identifier column to preserve.
#' @param feature_cols character vector. Names of feature columns to z-score.
#' @return data.frame with the same row order, the identifier column kept,
#'   and each feature column transformed to mean 0, sd 1. Errors out if any
#'   feature column has zero variance.
#' @export
zscore_features <- function(df, id_col, feature_cols) {
  out <- df[, c(id_col, feature_cols), drop = FALSE]
  for (col in feature_cols) {
    x <- df[[col]]
    s <- sd(x, na.rm = TRUE)
    if (!is.finite(s) || s == 0) {
      stop(sprintf("Feature '%s' has zero variance — cannot z-score.", col))
    }
    out[[col]] <- (x - mean(x, na.rm = TRUE)) / s
  }
  out
}

#' Leave-one-out 1-nearest-neighbour stratum predictor
#'
#' For each row, holds it out, restricts candidates to other rows with the
#' same `cell_col` value, computes Euclidean distance over `feature_cols`,
#' and assigns the label of the nearest candidate.
#'
#' @param features data.frame with id, cell, label, and feature columns.
#' @param id_col character. Identifier column name (e.g. "municipality").
#' @param cell_col character. Cell constraint column (e.g. "dept").
#' @param label_col character. Label column to predict (e.g. "stratum").
#' @param feature_cols character vector. Feature columns; assumed already
#'   z-scored.
#' @return data.frame with one row per input, columns:
#'   `<id_col>`, `<cell_col>`, `true_<label_col>`, `predicted_<label_col>`,
#'   `candidate_count`, `distance_to_nn`, `tie_flag`.
#' @export
predict_loo_nn <- function(features, id_col, cell_col, label_col,
                           feature_cols) {
  n <- nrow(features)
  X <- as.matrix(features[, feature_cols, drop = FALSE])
  ids <- features[[id_col]]
  cells <- features[[cell_col]]
  labels <- features[[label_col]]

  predicted <- rep(NA, n)
  predicted <- if (is.integer(labels)) as.integer(predicted) else predicted
  cand_count <- integer(n)
  dist_to_nn <- rep(NA_real_, n)
  tie_flag <- logical(n)

  for (i in seq_len(n)) {
    same_cell <- which(cells == cells[i])
    candidates <- setdiff(same_cell, i)
    cand_count[i] <- length(candidates)
    if (length(candidates) == 0L) next
    diffs <- sweep(X[candidates, , drop = FALSE], 2, X[i, ], FUN = "-")
    d2 <- rowSums(diffs * diffs)
    min_d2 <- min(d2)
    dist_to_nn[i] <- sqrt(min_d2)
    tie_flag[i] <- sum(d2 == min_d2) > 1L
    nn_idx <- candidates[which.min(d2)]
    predicted[i] <- labels[nn_idx]
  }

  out <- data.frame(
    id = ids,
    cell = cells,
    true_label = labels,
    predicted_label = predicted,
    candidate_count = cand_count,
    distance_to_nn = dist_to_nn,
    tie_flag = tie_flag,
    stringsAsFactors = FALSE
  )
  names(out)[1] <- id_col
  names(out)[2] <- cell_col
  names(out)[3] <- paste0("true_", label_col)
  names(out)[4] <- paste0("predicted_", label_col)
  out
}

#' Compute LOO 1-NN classification metrics
#'
#' @param results data.frame from `predict_loo_nn()` with columns
#'   `cell`, `true_<label>`, `predicted_<label>`, `<candidate_count>`,
#'   plus an identifier.
#' @param label_col character. Stem of the label columns (e.g. "stratum"),
#'   used to construct `true_<label>` and `predicted_<label>` names.
#' @param cell_col character. Cell column name (used to compute the
#'   non-trivial subset).
#' @param candidate_count_col character. Column giving candidate-set size
#'   per row.
#' @return list with named scalars:
#'   * `n_total`: total input rows
#'   * `n_predictable`: rows with `candidate_count > 0`
#'   * `n_singleton_cell`: rows with `candidate_count == 0`
#'   * `exact_match_rate`: matches / n_predictable
#'   * `n_non_trivial`: rows whose cell contains >1 distinct true label
#'   * `non_trivial_exact_match_rate`: matches over the non-trivial subset
#'   * `kendall_tau`: stats::cor(method = "kendall") over predictable rows
#' @export
loo_nn_metrics <- function(results, label_col, cell_col,
                           candidate_count_col) {
  true_col  <- paste0("true_",      label_col)
  pred_col  <- paste0("predicted_", label_col)
  n_total <- nrow(results)
  predictable <- results[[candidate_count_col]] > 0L
  n_predictable <- sum(predictable)
  n_singleton <- n_total - n_predictable

  pred <- results[predictable, ]
  matches <- pred[[true_col]] == pred[[pred_col]]
  exact_match_rate <- mean(matches)

  # Non-trivial cells: cells with > 1 distinct true label.
  cell_counts <- tapply(pred[[true_col]], pred[[cell_col]],
                        function(x) length(unique(x)))
  non_trivial_cells <- names(cell_counts)[cell_counts > 1L]
  nt <- pred[as.character(pred[[cell_col]]) %in% non_trivial_cells, ]
  if (nrow(nt) > 0L) {
    non_trivial_rate <- mean(nt[[true_col]] == nt[[pred_col]])
  } else {
    non_trivial_rate <- NA_real_
  }

  tau <- stats::cor(pred[[true_col]], pred[[pred_col]],
                    method = "kendall")

  list(
    n_total = n_total,
    n_predictable = n_predictable,
    n_singleton_cell = n_singleton,
    exact_match_rate = exact_match_rate,
    n_non_trivial = nrow(nt),
    non_trivial_exact_match_rate = non_trivial_rate,
    kendall_tau = tau
  )
}

#' Build a 1,122-muni → stratum map by extending DANE's known mapping via 1-NN
#'
#' For each muni in `aux` not present in `known`, predict its stratum by
#' 1-nearest-neighbour Euclidean distance over `feature_cols`, restricted
#' to in-frame munis in the same `dept`. The features are z-scored
#' internally over the full `aux` set before distance computation. In-frame
#' munis keep their `known` stratum unchanged.
#'
#' @param aux data.frame of all munis with id, dept, and feature columns.
#' @param known data.frame with `<id_col>` and `<stratum_col>` for the
#'   in-frame munis whose stratum is fixed.
#' @param id_col character. Identifier column name.
#' @param dept_col character. Cell column name (within which the NN is
#'   computed).
#' @param stratum_col character. Stratum column name.
#' @param feature_cols character vector. Feature columns to z-score and
#'   use for distance.
#' @return data.frame with columns `<id_col>`, `<dept_col>`,
#'   `<stratum_col>`, `source` (`"known"` for in-frame, `"nn"` for OOS).
#'   One row per muni in `aux`.
#' @export
extend_strata_nn <- function(aux, known, id_col, dept_col, stratum_col,
                             feature_cols) {
  stopifnot(all(c(id_col, dept_col, feature_cols) %in% colnames(aux)))
  stopifnot(all(c(id_col, stratum_col) %in% colnames(known)))
  zaux <- zscore_features(aux, id_col = id_col,
                          feature_cols = feature_cols)
  zaux[[dept_col]] <- aux[[dept_col]]
  in_frame_ids <- known[[id_col]]
  in_frame_mask <- zaux[[id_col]] %in% in_frame_ids
  in_frame <- zaux[in_frame_mask, , drop = FALSE]
  oos      <- zaux[!in_frame_mask, , drop = FALSE]
  stratum_lookup <- setNames(known[[stratum_col]], as.character(known[[id_col]]))
  in_frame[[stratum_col]] <- stratum_lookup[as.character(in_frame[[id_col]])]
  in_frame$source <- "known"

  if (nrow(oos) == 0L) {
    out <- in_frame[, c(id_col, dept_col, stratum_col, "source"),
                    drop = FALSE]
    return(out)
  }

  oos[[stratum_col]] <- NA_integer_
  oos$source <- "nn"
  X_in  <- as.matrix(in_frame[, feature_cols, drop = FALSE])
  X_oos <- as.matrix(oos[, feature_cols, drop = FALSE])
  for (i in seq_len(nrow(oos))) {
    same_dept <- which(in_frame[[dept_col]] == oos[[dept_col]][i])
    if (length(same_dept) == 0L) {
      stop(sprintf(
        "Muni %s (dept %s) has no in-frame candidate in its dept.",
        as.character(oos[[id_col]][i]), as.character(oos[[dept_col]][i])))
    }
    diffs <- sweep(X_in[same_dept, , drop = FALSE], 2, X_oos[i, ], FUN = "-")
    d2 <- rowSums(diffs * diffs)
    nn_idx <- same_dept[which.min(d2)]
    oos[[stratum_col]][i] <- in_frame[[stratum_col]][nn_idx]
  }

  out_rows <- rbind(
    in_frame[, c(id_col, dept_col, stratum_col, "source"), drop = FALSE],
    oos[,      c(id_col, dept_col, stratum_col, "source"), drop = FALSE])
  out_rows[order(out_rows[[id_col]]), , drop = FALSE]
}

#' Simulate stage-1 PPS-WOR draws over a universe of strata
#'
#' For each stratum, draws `n_psu` PSUs without replacement with
#' probability proportional to `size`, repeated for `S` independent
#' replicates. Uses `sampling::UPmaxentropy`, as [draw_psus()] does; a
#' stratum whose allocation covers all its PSUs short-circuits to a
#' deterministic take-all.
#'
#' @param universe data.frame with one row per muni: id, stratum, size.
#' @param allocation data.frame with one row per stratum: stratum and
#'   `n_psu` (the stage-1 allocation per stratum).
#' @param id_col character. Identifier column name.
#' @param stratum_col character. Stratum column name.
#' @param size_col character. Size measure column name (e.g. n_dwellings).
#' @param n_psu_col character. Allocation column name in `allocation`.
#' @param S integer. Number of independent replicates.
#' @param base_seed integer. Per-rep seed = `base_seed + s - 1` for
#'   `s = 1, ..., S`.
#' @return data.frame with one row per muni: id, stratum, `n_drawn`
#'   (times drawn across S reps), `emp_pi` = n_drawn / S, `theoretical_pi`,
#'   `n_reps` = S.
#' @export
simulate_stage1_pps <- function(universe, allocation,
                                id_col, stratum_col, size_col,
                                n_psu_col, S, base_seed = 100000L) {
  alloc_lookup <- setNames(allocation[[n_psu_col]],
                            as.character(allocation[[stratum_col]]))
  n_drawn <- integer(nrow(universe))
  theoretical_pi <- numeric(nrow(universe))
  strata <- unique(universe[[stratum_col]])

  for (s in seq_len(S)) {
    set.seed(base_seed + s - 1L)
    for (h in strata) {
      idx_h <- which(universe[[stratum_col]] == h)
      sizes_h <- universe[[size_col]][idx_h]
      n_psu_h <- alloc_lookup[[as.character(h)]]
      if (is.null(n_psu_h) || is.na(n_psu_h)) {
        next  # no allocation for this stratum: skip
      }
      if (n_psu_h >= length(idx_h)) {
        n_drawn[idx_h] <- n_drawn[idx_h] + 1L
        next
      }
      pik <- sampling::inclusionprobabilities(sizes_h, n_psu_h)
      # UPmaxentropy crashes when any pik == 1 but n_psu < N (certainty units
      # arise when one PSU is very large). Handle by separating certainty units,
      # then drawing from the remainder.
      certain <- pik >= 1 - .Machine$double.eps * 10
      if (any(certain)) {
        n_drawn[idx_h[certain]] <- n_drawn[idx_h[certain]] + 1L
        n_remain <- n_psu_h - sum(certain)
        if (n_remain <= 0L) next
        idx_rem <- idx_h[!certain]
        sizes_rem <- sizes_h[!certain]
        pik_rem <- sampling::inclusionprobabilities(sizes_rem, n_remain)
        if (n_remain >= length(idx_rem)) {
          n_drawn[idx_rem] <- n_drawn[idx_rem] + 1L
        } else {
          drawn_rem <- sampling::UPmaxentropy(pik_rem)
          n_drawn[idx_rem] <- n_drawn[idx_rem] + as.integer(drawn_rem)
        }
      } else {
        drawn <- sampling::UPmaxentropy(pik)
        n_drawn[idx_h] <- n_drawn[idx_h] + as.integer(drawn)
      }
    }
  }

  for (h in strata) {
    idx_h <- which(universe[[stratum_col]] == h)
    sizes_h <- universe[[size_col]][idx_h]
    n_psu_h <- alloc_lookup[[as.character(h)]]
    if (is.null(n_psu_h) || is.na(n_psu_h)) {
      theoretical_pi[idx_h] <- NA_real_
      next
    }
    if (n_psu_h >= length(idx_h)) {
      theoretical_pi[idx_h] <- 1
    } else {
      theoretical_pi[idx_h] <-
        sampling::inclusionprobabilities(sizes_h, n_psu_h)
    }
  }

  out <- data.frame(
    id = universe[[id_col]],
    stratum = universe[[stratum_col]],
    n_drawn = n_drawn,
    emp_pi = n_drawn / S,
    theoretical_pi = theoretical_pi,
    n_reps = S
  )
  names(out)[1] <- id_col
  names(out)[2] <- stratum_col
  out
}
