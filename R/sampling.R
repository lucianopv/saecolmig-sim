#' Count dwellings per PSU within stratum
#'
#' Per the ECV 2024 design (DDI pp.2-3, 11), the PSU is the municipality;
#' the function is generic over the PSU column name, for other designs.
#'
#' @param frame household-level frame (must contain stratum and the PSU column)
#' @param psu_var character column name for the PSU (default "municipality")
#' @return tibble with one row per (stratum, psu) and column n_dwellings
#' @export
compute_psu_sizes <- function(frame, psu_var = "municipality") {
  if (!psu_var %in% colnames(frame)) {
    stop("PSU column '", psu_var, "' not found in frame.")
  }
  dt <- data.table::as.data.table(frame)
  out <- dt[, .(n_dwellings = .N),
            by = c("stratum", psu_var)]
  data.table::setnames(out, psu_var, "psu_id")
  tibble::as_tibble(out)
}

#' PPS-without-replacement draw of PSUs within one stratum
#'
#' Uses sampling::UPmaxentropy with inclusion probabilities proportional to
#' the dwellings count, with short-circuits for degenerate cases:
#'   - n_psu == nrow(sizes): take all PSUs (forced inclusion), pi_psu = 1.
#'   - n_psu == 1: a single multinomial PPS draw (no UPmaxentropy, which is
#'     numerically unstable for small samples).
#'
#' @param sizes tibble with columns (stratum, psu_id, n_dwellings) for one stratum
#' @param n_psu integer number of PSUs to select
#' @param seed integer
#' @return tibble of selected PSUs with first-stage inclusion probability pi_psu
#' @export
draw_psus <- function(sizes, n_psu, seed) {
  set.seed(seed)
  N <- nrow(sizes)
  if (n_psu > N) {
    stop("Requested n_psu (", n_psu, ") exceeds available PSUs (",
         N, ") in stratum ", sizes$stratum[1])
  }
  if (n_psu == N) {
    # Forced inclusion / take-all stratum
    out <- sizes
    out$pi_psu <- 1
    return(tibble::as_tibble(out))
  }
  if (n_psu == 1) {
    # Single PPS draw: probabilities proportional to size, sum to 1
    p <- sizes$n_dwellings / sum(sizes$n_dwellings)
    pick <- sample.int(N, size = 1, prob = p)
    out <- sizes[pick, , drop = FALSE]
    out$pi_psu <- p[pick]   # = inclusion prob for n_psu = 1
    return(tibble::as_tibble(out))
  }
  pik <- sampling::inclusionprobabilities(sizes$n_dwellings, n_psu)
  # Certainty units (pik == 1) crash sampling::UPmaxentropy. Force-include
  # them up front and call UPmaxentropy only on the sub-certainty remainder.
  cert <- pik >= 1 - 10 * .Machine$double.eps
  if (any(cert)) {
    n_cert <- sum(cert)
    if (n_cert == n_psu) {
      out <- sizes[cert, , drop = FALSE]
      out$pi_psu <- 1
      return(tibble::as_tibble(out))
    }
    sizes_sub <- sizes[!cert, , drop = FALSE]
    n_remain  <- n_psu - n_cert
    pik_sub   <- sampling::inclusionprobabilities(sizes_sub$n_dwellings,
                                                  n_remain)
    # Second-pass certainty check: if n_remain == nrow(sizes_sub), all
    # pik_sub are 1 and UPmaxentropy crashes with "subscript out of bounds".
    cert_sub <- pik_sub >= 1 - 10 * .Machine$double.eps
    if (all(cert_sub)) {
      drawn_sub <- sizes_sub
      drawn_sub$pi_psu <- 1
    } else {
      s_sub <- sampling::UPmaxentropy(pik_sub)
      drawn_sub <- sizes_sub[as.logical(s_sub), , drop = FALSE]
      drawn_sub$pi_psu <- pik_sub[as.logical(s_sub)]
    }
    cert_rows <- sizes[cert, , drop = FALSE]
    cert_rows$pi_psu <- 1
    return(tibble::as_tibble(rbind(cert_rows, drawn_sub)))
  }
  s   <- sampling::UPmaxentropy(pik)
  out <- sizes[as.logical(s), ]
  out$pi_psu <- pik[as.logical(s)]
  tibble::as_tibble(out)
}

#' Build ~10-dwelling segments within a single PSU
#'
#' Implements ECV's USM construction: order the dwellings within the PSU by
#' (sector, seccion, manzana) - approximated by ordering on the manzana code,
#' which already encodes the sector and section (characters 7-11 for rural
#' records, 15-20 for urban ones) - then chunk into
#' consecutive runs of `target_size` dwellings. The last segment may be
#' smaller than `target_size`.
#'
#' Returns the input frame augmented with a `segment_id` column unique within
#' the PSU.
#'
#' @param psu_frame household-level frame for one PSU (must contain manzana)
#' @param target_size integer; default 10 (per DDI)
#' @return psu_frame with added segment_id column
#' @export
build_segments <- function(psu_frame, target_size = 10L) {
  if (!"manzana" %in% colnames(psu_frame)) {
    stop("psu_frame must contain a 'manzana' column for ordering.")
  }
  ord <- order(as.character(psu_frame$manzana), psu_frame$hh_id)
  psu_frame <- psu_frame[ord, , drop = FALSE]
  n <- nrow(psu_frame)
  seg_idx <- ((seq_len(n) - 1L) %/% target_size) + 1L
  psu_frame$segment_id <- paste(psu_frame$municipality[1], seg_idx, sep = "_S")
  tibble::as_tibble(psu_frame)
}

#' Pre-build segments for every PSU in the frame
#'
#' Segments depend only on frame structure (deterministic given target_size).
#' Computing them once and reusing across replicates removes the largest
#' per-rep cost in draw_replicate.
#'
#' Only the columns the sampler downstream needs are kept (geo keys, hh_id,
#' segment_id, the outcomes, stratum, n_persons), which keeps the cache small
#' when it is copied to parallel worker processes.
#'
#' @param frame household-level frame (must contain psu_var, manzana, hh_id)
#' @param psu_var character; PSU column (default "municipality")
#' @param target_size integer; default 10
#' @return list keyed by character(psu_id); each element is the PSU's slim
#'   frame augmented with a segment_id column
#' @export
build_all_segments <- function(frame, psu_var = "municipality",
                               target_size = 10L) {
  keep_cols <- intersect(
    c("dept", "muni", "municipality", "area", "stratum",
      "seccion", "manzana", "hh_id", "n_persons",
      "migrated_last_12m", "migrated_internal_12m", "migrated_intl_12m",
      "n_pers_resp", "n_mig_total", "n_mig_internal", "n_mig_intl"),
    colnames(frame)
  )
  fr <- as.data.frame(frame)[, keep_cols, drop = FALSE]
  by_psu <- split(fr, fr[[psu_var]])
  out <- lapply(by_psu, build_segments, target_size = target_size)
  names(out) <- as.character(names(by_psu))
  out
}

#' Systematic selection of segments within a PSU
#'
#' Implements the Fractional Interval Method (DDI p.11) for systematic
#' selection of n_seg segments from N total. Each segment has the same
#' inclusion probability n_seg / N. The function returns the dwelling rows
#' for the selected segments, with the second-stage inclusion probability
#' pi_seg attached.
#'
#' @param psu_frame_segmented output of build_segments() (one PSU, with segment_id)
#' @param n_seg integer number of segments to select
#' @param seed integer
#' @return tibble of dwellings in the selected segments with pi_seg
#' @export
draw_segments <- function(psu_frame_segmented, n_seg, seed) {
  set.seed(seed)
  if (!"segment_id" %in% colnames(psu_frame_segmented)) {
    stop("psu_frame_segmented must contain a 'segment_id' column ",
         "(call build_segments first).")
  }
  segs <- unique(psu_frame_segmented$segment_id)
  N <- length(segs)
  if (n_seg > N) {
    stop("Requested n_seg (", n_seg, ") exceeds available segments (",
         N, ") in PSU ", psu_frame_segmented$municipality[1])
  }
  # Fractional Interval Method: random start in [0, N/n_seg), step N/n_seg.
  step    <- N / n_seg
  start   <- runif(1, min = 0, max = step)
  picks   <- floor(start + step * (seq_len(n_seg) - 1L)) + 1L
  picks   <- pmin(picks, N)               # guard against floating-point overshoot
  chosen  <- segs[picks]
  out     <- psu_frame_segmented[psu_frame_segmented$segment_id %in% chosen,
                                 , drop = FALSE]
  out$pi_seg <- n_seg / N
  tibble::as_tibble(out)
}

#' Draw one ECV-design replicate sample
#'
#' Implements the ECV 2024 two-stage design:
#'   stage 1: PPS-WOR selection of municipalities within stratum (size = dwellings)
#'   stage 2: systematic selection of ~10-dwelling segments within municipality
#'   (no within-segment selection: every dwelling/household in the segment is
#'   enumerated)
#'
#' Returns the drawn dwellings with both inclusion probabilities and the
#' design weight (pre-calibration). Calibration is applied separately by
#' apply_calibration().
#'
#' @param frame household-level national frame (output of build_national_frame),
#'   must contain `stratum`, `municipality`, `manzana`, `hh_id`,
#'   `migrated_last_12m`. The frame may additionally carry migration-type
#'   indicator columns (`migrated_internal_12m`, `migrated_intl_12m`) and
#'   person-level count columns (`n_pers_resp`, `n_mig_total`,
#'   `n_mig_internal`, `n_mig_intl`); these flow through the segments cache
#'   and the drawn sample unchanged.
#' @param design_spec object from parse_ecv_design()
#' @param seed integer; replicate seed (e.g. 100000 + i)
#' @param segments optional list (output of build_all_segments()); when supplied
#'   the per-PSU build_segments() call is skipped and pre-segmented frames are
#'   looked up by character(psu_id). Massive speed-up across replicates.
#' @return tibble of sampled households
#' @export
draw_replicate <- function(frame, design_spec, seed, segments = NULL) {
  set.seed(seed)
  required <- c("stratum", "municipality", "manzana", "hh_id")
  missing  <- setdiff(required, colnames(frame))
  if (length(missing) > 0) {
    stop("frame is missing required columns: ", paste(missing, collapse = ", "))
  }
  frame_dt <- data.table::as.data.table(frame)
  by_muni  <- if (is.null(segments)) split(frame_dt, frame_dt$municipality) else NULL

  sizes <- compute_psu_sizes(frame, psu_var = "municipality")
  per_stratum <- split(sizes, sizes$stratum)
  target_size       <- .null_or(design_spec$stage2$target_segment_size, 10L)
  global_n_per_psu  <- .null_or(design_spec$stage2$n_per_psu,           1L)
  has_per_stratum_n <- "n_seg_per_psu" %in% colnames(design_spec$strata)

  out_list <- list()
  for (stratum_name in names(per_stratum)) {
    row_idx <- which(design_spec$strata$stratum == stratum_name)
    if (length(row_idx) == 0) next
    n_psu_target  <- design_spec$strata$n_psu[row_idx]
    n_seg_per_psu <- if (has_per_stratum_n) {
      design_spec$strata$n_seg_per_psu[row_idx]
    } else {
      global_n_per_psu
    }

    psu_seed   <- seed + 1L
    drawn_psus <- draw_psus(per_stratum[[stratum_name]],
                            n_psu = n_psu_target, seed = psu_seed)

    for (i in seq_len(nrow(drawn_psus))) {
      psu_row    <- drawn_psus[i, ]
      psu_key    <- as.character(psu_row$psu_id)
      seg_seed   <- seed + 1000L + i

      if (!is.null(segments) && !is.null(segments[[psu_key]])) {
        segmented <- segments[[psu_key]]
      } else {
        segmented <- build_segments(by_muni[[psu_key]],
                                    target_size = target_size)
      }

      n_seg_eff  <- min(n_seg_per_psu, length(unique(segmented$segment_id)))
      drawn_segs <- draw_segments(segmented, n_seg = n_seg_eff, seed = seg_seed)
      drawn_segs$pi_psu <- psu_row$pi_psu
      drawn_segs$psu_id <- psu_row$psu_id
      drawn_segs$ssu_id <- drawn_segs$segment_id
      out_list[[length(out_list) + 1L]] <- drawn_segs
    }
  }

  out <- data.table::rbindlist(out_list, use.names = TRUE, fill = TRUE)
  out[, weight   := 1 / (pi_psu * pi_seg)]
  out[, rep_seed := seed]
  tibble::as_tibble(out)
}

#' Calibrate sample weights to stratum household (and optionally person) totals
#'
#' Uses raking (`survey::calibrate(..., calfun = "raking")`) to adjust design
#' weights so the calibrated weights reproduce population margins within each
#' stratum.
#'
#' **Single margin:** when `pop_totals` contains only `stratum` and `N`,
#' calibration uses `formula = ~stratum` so the calibrated weights sum to the
#' household count `N` within every stratum.
#'
#' **Dual margin:** when `pop_totals` also contains
#' `N_persons` *and* the sample carries an `n_persons` column, calibration
#' rakes to per-stratum person-total constraints —
#' `sum(weight_cal * n_persons) == N_persons` — while keeping the household
#' constraint. The approach mirrors the real ECV's calibration (person
#' sex×age cells + households-per-dwelling row producing a single factor
#' valid at both levels) at total-persons granularity without requiring the
#' full sex×age breakdown at simulation stage.
#'
#' **NA handling (dual-margin):** households with `NA` `n_persons` (e.g.
#' HA_TOT_PER missing in the census, ~0.066\% of households) contribute 0 to
#' the person-margin calibration variable. This mirrors the `na.rm = TRUE`
#' convention used when computing the frame person totals: those households'
#' persons are excluded from BOTH sides of the margin, so the constraint
#' stays coherent and `survey::calibrate` is not poisoned by NAs in the
#' model matrix. Internally a temporary `.n_persons_cal` column carries the
#' NA→0 recode; it is removed before the sample is returned.
#'
#' The `population` vector passed to `survey::calibrate` is constructed to
#' match the model-matrix column order exactly:
#'   `(Intercept)`, stratum dummies (all sorted levels except the reference),
#'   then `stratum<level>:.n_persons_cal` for every sorted level.
#' `pop_totals` is sorted by stratum internally before extracting offsets so
#' the reference level (first alphabetically) aligns with R's factor encoding.
#'
#' Note: the real ECV applies Deville-Sarndal linear-assisted calibration to
#' sex×age×dwelling-count cells with DANE population projections. The dual-
#' margin raking here is a simulation-stage approximation of it.
#'
#' @param sample tibble from draw_replicate (must include stratum, weight, psu_id;
#'   must include n_persons for dual-margin mode)
#' @param pop_totals data.frame with columns (stratum, N) for single-margin mode,
#'   or (stratum, N, N_persons) for dual-margin mode
#' @return sample with new column weight_cal
#' @export
apply_calibration <- function(sample, pop_totals) {
  # Restrict to strata present in the sample and sort for alignment
  pop_totals <- pop_totals[pop_totals$stratum %in% unique(sample$stratum), ,
                           drop = FALSE]
  pop_totals <- pop_totals[order(pop_totals$stratum), , drop = FALSE]

  has_pers <- "N_persons" %in% colnames(pop_totals) &&
              "n_persons" %in% colnames(sample)
  if (has_pers) {
    # NA household sizes contribute 0 to the person margin, mirroring the
    # na.rm = TRUE person totals (those persons are excluded on both sides).
    sample$.n_persons_cal <- ifelse(is.na(sample$n_persons), 0,
                                    sample$n_persons)
  }

  des <- survey::svydesign(ids = ~psu_id, strata = ~stratum,
                           weights = ~weight, nest = TRUE, data = sample)

  # Named household-count vector: stratum<level> for each sorted level
  hh_tot <- setNames(pop_totals$N, paste0("stratum", pop_totals$stratum))

  if (!has_pers) {
    # Single-margin: ~stratum model matrix has (Intercept) + k-1 dummies.
    # population = c(sum(N), N[-reference]) with names matching column names.
    totals <- c(`(Intercept)` = sum(pop_totals$N), hh_tot[-1])
    cal <- survey::calibrate(des, formula = ~stratum, population = totals,
                             calfun = "raking", maxit = 100)
  } else {
    # Dual-margin: ~stratum + stratum:.n_persons_cal model matrix columns are:
    #   (Intercept)
    #   stratum<level2>, ..., stratum<levelK>   [k-1 dummies, reference dropped]
    #   stratum<level1>:.n_persons_cal, ..., stratum<levelK>:.n_persons_cal
    # where levels are sorted alphabetically (R default factor encoding).
    pers_tot <- setNames(pop_totals$N_persons,
                         paste0("stratum", pop_totals$stratum,
                                ":.n_persons_cal"))
    totals <- c(`(Intercept)` = sum(pop_totals$N), hh_tot[-1], pers_tot)
    cal <- survey::calibrate(des,
                             formula    = ~stratum + stratum:.n_persons_cal,
                             population = totals,
                             calfun     = "raking", maxit = 100)
  }

  sample$weight_cal <- weights(cal)
  sample$.n_persons_cal <- NULL
  sample
}
