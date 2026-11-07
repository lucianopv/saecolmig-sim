#' Inputs of the replicate runner
#'
#' The runner draws samples from the household frame under the ECV design,
#' estimates the six outcomes for the municipalities of the frame (the
#' universe) and scores the estimates against the census truth. This file holds
#' the one entry point to everything it reads, [sim_inputs()], and the format
#' checks that run before any fit: a matrix with a missing municipality, a
#' missing covariate value or a zero-padded text code is refused (or coerced)
#' here, naming the file, the problem and the count, rather than failing inside
#' a model.
#'
#' Covariates reach the runner in two ways. Extra columns are joined to the
#' built matrix on its key (`extra_area`, `extra_subarea`); or a complete matrix
#' replaces the built one (`area_matrix`, `subarea_matrix`). Both are checked
#' with the same rules. The file layouts are in `DATA_DICTIONARY.md`.
#'
#' @name sim_inputs_overview
NULL

# Columns of the municipal matrix that are not covariates. `dept` and `muni` are
# the two parts of the municipality code. `n_persons` and `n_households` are
# sizes that span orders of magnitude against shares near one; a covariate of
# that scale squares into the cross-product matrix of the robust Fay-Herriot fit
# and leaves it numerically singular, so no area-level method takes them.
.sim_area_drop <- c("dept", "muni", "n_persons", "n_households")

# Sub-area class labels carried by the `stratum` column of the sub-area matrix;
# collapse_sample_to_subarea() matches sample rows to sub-areas through them.
.sim_subarea_classes <- c("urban", "rural")

# "1 value" / "3 values", and the matching verb.
.sim_np <- function(n, one, many) paste(n, if (n == 1L) one else many)
.sim_isare <- function(n) if (n == 1L) "is" else "are"

# First few offending codes, for the messages that name them.
.sim_first <- function(x, k = 5L) {
  paste0(" (first: ", paste(utils::head(x, k), collapse = ", "), ")")
}

# "area covariate matrix" or "area covariate matrix (path)".
.sim_where <- function(what, source = NULL) {
  if (is.null(source) || !nzchar(source)) what else sprintf("%s (%s)", what, source)
}

.sim_stop <- function(...) stop(..., call. = FALSE)

# Integer codes from integer, whole-number double, or digit-string input (so a
# zero-padded text code such as "05001" is read as 5001). A value that is NA or
# is not a whole number stops with its count.
.sim_as_integer_key <- function(x, name, where) {
  if (is.factor(x)) x <- as.character(x)
  if (is.character(x)) {
    x <- trimws(x)
    num <- suppressWarnings(as.numeric(ifelse(grepl("^[0-9]+$", x), x, NA_character_)))
  } else if (is.numeric(x)) {
    num <- as.numeric(x)
  } else {
    .sim_stop(where, ": `", name, "` must hold integer codes, not ", class(x)[1])
  }
  n_na <- sum(is.na(x))
  if (n_na) {
    .sim_stop(where, ": ", .sim_np(n_na, "value", "values"), " of `", name, "` ",
              .sim_isare(n_na), " missing")
  }
  bad <- is.na(num) | !is.finite(num) | num != round(num) |
    abs(num) > .Machine$integer.max
  if (any(bad)) {
    n_bad <- sum(bad)
    .sim_stop(where, ": ", .sim_np(n_bad, "value", "values"), " of `", name, "` ",
              .sim_isare(n_bad), " not a whole number")
  }
  as.integer(num)
}

# Stop unless `path` is a readable file; `hint` says how to produce it.
.sim_need_file <- function(path, what, hint = NULL) {
  if (!file.exists(path)) {
    .sim_stop(what, " (", path, "): file not found", if (!is.null(hint)) paste0("; ", hint))
  }
  invisible(path)
}

# A table from a data frame or from a parquet / CSV path. Key columns of a CSV
# are read as text so that leading zeros survive; they are coerced afterwards.
# Returns list(data, source), `source` being the path or NULL.
.sim_read_input <- function(x, what, hint = NULL) {
  if (is.data.frame(x)) return(list(data = as.data.frame(x), source = NULL))
  if (!is.character(x) || length(x) != 1L || is.na(x)) {
    .sim_stop(what, ": expected a file path or a data frame")
  }
  .sim_need_file(x, what, hint)
  d <- if (grepl("\\.parquet$", x, ignore.case = TRUE)) {
    as.data.frame(arrow::read_parquet(x))
  } else if (grepl("\\.csv$", x, ignore.case = TRUE)) {
    key <- intersect(c("municipality", "code", "subarea", "area_id"),
                     names(data.table::fread(x, nrows = 0L)))
    cc <- if (length(key)) list(character = key) else NULL
    as.data.frame(data.table::fread(x, colClasses = cc, data.table = FALSE))
  } else {
    .sim_stop(what, " (", x, "): only .parquet and .csv files are read")
  }
  list(data = d, source = x)
}

# The sub-area key is `code`; a table that names it `subarea` is renamed.
.sim_subarea_key <- function(d) {
  if (!"code" %in% names(d) && "subarea" %in% names(d)) names(d)[names(d) == "subarea"] <- "code"
  d
}

# Names of the columns of `d` (among `cols`) that have no observed value.
.sim_empty_cols <- function(d, cols) {
  cols[vapply(d[cols], function(v) all(is.na(v)), logical(1))]
}

# Stop on covariate columns that are not numeric, hold NA or are infinite.
.sim_check_covariates <- function(d, cols, where, allow_na = FALSE) {
  non_num <- cols[!vapply(d[cols], is.numeric, logical(1))]
  if (length(non_num)) {
    .sim_stop(where, ": ", .sim_np(length(non_num), "covariate column is",
                                   "covariate columns are"),
              " not numeric: ", paste(non_num, collapse = ", "))
  }
  inf <- vapply(d[cols], function(v) sum(is.infinite(v)), integer(1))
  if (any(inf > 0L)) {
    inf <- inf[inf > 0L]
    .sim_stop(where, ": covariate columns with infinite values: ",
              paste0(names(inf), " (", inf, ")", collapse = ", "))
  }
  if (!allow_na) {
    na <- vapply(d[cols], function(v) sum(is.na(v)), integer(1))
    if (any(na > 0L)) {
      na <- na[na > 0L]
      .sim_stop(where, ": covariate columns with missing values: ",
                paste0(names(na), " (", na, " ", ifelse(na == 1L, "row", "rows"), ")",
                       collapse = ", "))
    }
  }
  invisible(d)
}

#' Check an area-level covariate matrix
#'
#' A usable matrix has the integer key `municipality`, exactly one row for each
#' code of the universe, and numeric covariates without a missing or infinite
#' value. The key may arrive as a whole-number double or as a zero-padded text
#' code and is stored as integer. Identifier and size columns (`dept`, `muni`,
#' `n_persons`, `n_households`) are dropped: the two code parts are not
#' covariates, and sizes that span orders of magnitude against shares near one
#' make the cross-product matrix of the robust Fay-Herriot fit ill-conditioned.
#' A column with no observed value carries no information and is dropped too.
#' Row and column order are kept.
#'
#' @param X data frame
#' @param universe integer municipality codes the matrix must cover
#' @param source file the matrix was read from, named in the messages
#' @return `X` as a plain data frame with the checks applied
#' @export
sim_validate_area_matrix <- function(X, universe, source = NULL) {
  where <- .sim_where("area covariate matrix", source)
  if (!is.data.frame(X)) .sim_stop(where, ": expected a data frame")
  X <- as.data.frame(X)
  if (!"municipality" %in% names(X)) .sim_stop(where, ": no `municipality` column")
  X$municipality <- .sim_as_integer_key(X$municipality, "municipality", where)

  n_dup <- sum(duplicated(X$municipality))
  if (n_dup) {
    .sim_stop(where, ": ", .sim_np(n_dup, "municipality code appears",
                                   "municipality codes appear"), " more than once")
  }
  universe <- sort(unique(as.integer(universe)))
  absent <- setdiff(universe, X$municipality)
  if (length(absent)) {
    .sim_stop(where, ": ", .sim_np(length(absent), "municipality", "municipalities"),
              " of the universe ", .sim_isare(length(absent)), " missing",
              .sim_first(absent))
  }
  extra <- setdiff(X$municipality, universe)
  if (length(extra)) {
    .sim_stop(where, ": ", .sim_np(length(extra), "municipality", "municipalities"),
              " of the matrix ", .sim_isare(length(extra)), " not in the universe",
              .sim_first(extra))
  }

  drop_id <- intersect(.sim_area_drop, names(X))
  covs <- setdiff(names(X), c("municipality", drop_id))
  empty <- .sim_empty_cols(X, covs)
  if (length(drop_id)) {
    message(where, ": identifier and size columns are not covariates and are dropped: ",
            paste(drop_id, collapse = ", "))
  }
  if (length(empty)) {
    message(where, ": columns without an observed value are dropped: ",
            paste(empty, collapse = ", "))
  }
  covs <- setdiff(covs, empty)
  if (!length(covs)) {
    .sim_stop(where, ": no covariate column besides the key and the identifier columns")
  }
  .sim_check_covariates(X, covs, where)

  X <- X[, setdiff(names(X), c(drop_id, empty)), drop = FALSE]
  rownames(X) <- NULL
  X
}

#' Check a sub-area covariate matrix
#'
#' A sub-area is a census block (urban) or a rural section. A usable matrix has
#' the character key `code` (one row per sub-area; a numeric code would lose its
#' leading zeros, and a column named `subarea` is read as `code`), the integer
#' `municipality`, the class `stratum` (`"urban"` or `"rural"`), the census
#' household and person counts `n_households_census` and `n_persons_census` that
#' weight the roll-up to municipalities, and numeric covariates. The prediction
#' universe is the set of sub-areas complete on every covariate
#' ([subarea_universe()]); the others are left out and their number is reported.
#' Covariate columns with no observed value are dropped, as at area level.
#'
#' @param pop_sub data frame
#' @param source file the matrix was read from, named in the messages
#' @return the sub-area universe as a plain data frame
#' @export
sim_validate_subarea_matrix <- function(pop_sub, source = NULL) {
  where <- .sim_where("sub-area matrix", source)
  if (!is.data.frame(pop_sub)) .sim_stop(where, ": expected a data frame")
  P <- .sim_subarea_key(as.data.frame(pop_sub))
  need <- c("code", "municipality", "stratum", "n_households_census", "n_persons_census")
  absent <- setdiff(need, names(P))
  if (length(absent)) {
    .sim_stop(where, ": ", .sim_np(length(absent), "required column is", "required columns are"),
              " missing: ", paste(absent, collapse = ", "))
  }

  if (is.factor(P$code)) P$code <- as.character(P$code)
  if (!is.character(P$code)) {
    .sim_stop(where, ": `code` must be character (a numeric code loses its leading zeros), not ",
              class(P$code)[1])
  }
  n_blank <- sum(is.na(P$code) | !nzchar(P$code))
  if (n_blank) {
    .sim_stop(where, ": ", .sim_np(n_blank, "value", "values"), " of `code` ",
              .sim_isare(n_blank), " missing")
  }
  n_dup <- sum(duplicated(P$code))
  if (n_dup) {
    .sim_stop(where, ": ", .sim_np(n_dup, "sub-area code appears", "sub-area codes appear"),
              " more than once")
  }
  P$municipality <- .sim_as_integer_key(P$municipality, "municipality", where)
  n_cls <- sum(!P$stratum %in% .sim_subarea_classes)
  if (n_cls) {
    .sim_stop(where, ": ", .sim_np(n_cls, "value", "values"), " of `stratum` ",
              .sim_isare(n_cls), " not \"urban\" or \"rural\"")
  }
  sizes <- intersect(c("n_households_census", "n_persons_census", "area_km2", "pop_worldpop"),
                     names(P))
  .sim_check_covariates(P, sizes, where, allow_na = TRUE)
  for (cc in c("n_households_census", "n_persons_census")) {
    n_na <- sum(is.na(P[[cc]]))
    if (n_na) {
      .sim_stop(where, ": `", cc, "` has ", .sim_np(n_na, "missing value", "missing values"))
    }
  }

  covs <- setdiff(names(P), .subarea_id_cols)
  empty <- .sim_empty_cols(P, covs)
  if (length(empty)) {
    message(where, ": covariate columns without an observed value are dropped: ",
            paste(empty, collapse = ", "))
    P <- P[, setdiff(names(P), empty), drop = FALSE]
    covs <- setdiff(covs, empty)
  }
  if (!length(covs)) {
    .sim_stop(where, ": no covariate column besides the structural columns (",
              paste(.subarea_id_cols, collapse = ", "), ")")
  }
  .sim_check_covariates(P, covs, where, allow_na = TRUE)

  keep <- subarea_universe(P)
  n_out <- nrow(P) - nrow(keep)
  if (n_out) {
    message(where, ": ", .sim_np(n_out, "sub-area", "sub-areas"),
            " incomplete on a covariate ", if (n_out == 1L) "is" else "are",
            " left out of the prediction universe")
  }
  if (!nrow(keep)) .sim_stop(where, ": no sub-area is complete on every covariate")
  rownames(keep) <- NULL
  keep
}

#' Join extra covariates to a covariate matrix
#'
#' A left join on `key` that keeps the row order of `X`. It stops, instead of
#' suffixing, when an extra column has the name of a column already in `X`; when
#' a key of `X` has no row in `extra` (its covariate would be missing); and when
#' `extra` repeats a key. Rows of `extra` for keys outside `X` are ignored.
#' A text key is coerced to the type of the key of `X`.
#'
#' @param X data frame with the column `key`
#' @param extra data frame with the column `key` and at least one other column
#' @param key name of the join column
#' @param source file `extra` was read from, named in the messages
#' @return `X` with the extra columns appended
#' @export
sim_add_covariates <- function(X, extra, key, source = NULL) {
  where <- .sim_where("extra covariates", source)
  X <- as.data.frame(X); extra <- as.data.frame(extra)
  if (!key %in% names(X)) .sim_stop("the matrix has no `", key, "` column")
  if (!key %in% names(extra)) .sim_stop(where, ": no `", key, "` column")
  new <- setdiff(names(extra), key)
  if (!length(new)) .sim_stop(where, ": no covariate column besides `", key, "`")
  clash <- intersect(new, names(X))
  if (length(clash)) {
    .sim_stop(where, ": ", .sim_np(length(clash), "column name is", "column names are"),
              " already in the matrix: ", paste(clash, collapse = ", "))
  }
  extra[[key]] <- if (is.integer(X[[key]])) {
    .sim_as_integer_key(extra[[key]], key, where)
  } else {
    as.character(extra[[key]])
  }
  n_dup <- sum(duplicated(extra[[key]]))
  if (n_dup) {
    .sim_stop(where, ": ", .sim_np(n_dup, "key appears", "keys appear"), " more than once")
  }
  idx <- match(X[[key]], extra[[key]])
  n_none <- sum(is.na(idx))
  if (n_none) {
    .sim_stop(where, ": ", .sim_np(n_none, "row of the matrix has", "rows of the matrix have"),
              " no match in the extra covariates (key `", key, "`)")
  }
  n_unused <- sum(!extra[[key]] %in% X[[key]])
  if (n_unused) {
    message(where, ": ", .sim_np(n_unused, "row", "rows"), " with a key outside the matrix ",
            .sim_isare(n_unused), " ignored")
  }
  for (nm in new) X[[nm]] <- extra[[nm]][idx]
  X
}

# The census truth: one row per area and outcome for every area of the universe.
.sim_read_truth <- function(path, universe) {
  got <- .sim_read_input(path, "truth", hint = "write it with data-raw/compute_truth_full.R")
  where <- .sim_where("truth", got$source)
  tr <- got$data
  absent <- setdiff(c("area_id", "theta_true", "n_total", "outcome"), names(tr))
  if (length(absent)) {
    .sim_stop(where, ": missing column(s): ", paste(absent, collapse = ", "))
  }
  tr$area_id <- .sim_as_integer_key(tr$area_id, "area_id", where)
  n_na <- sum(!is.finite(tr$theta_true))
  if (n_na) {
    .sim_stop(where, ": ", .sim_np(n_na, "value", "values"), " of `theta_true` ",
              .sim_isare(n_na), " missing or not finite")
  }
  n_dup <- sum(duplicated(tr[, c("area_id", "outcome")]))
  if (n_dup) {
    .sim_stop(where, ": ", .sim_np(n_dup, "area-outcome pair appears", "area-outcome pairs appear"),
              " more than once")
  }
  want <- expand.grid(area_id = universe, outcome = sae_outcomes$outcome,
                      stringsAsFactors = FALSE)
  n_miss <- sum(!paste(want$area_id, want$outcome) %in% paste(tr$area_id, tr$outcome))
  if (n_miss) {
    .sim_stop(where, ": ", .sim_np(n_miss, "area-outcome combination", "area-outcome combinations"),
              " of the universe ", .sim_isare(n_miss), " missing")
  }
  tr
}

# Largest share of the frame's households that may sit in sub-areas the matrix
# does not hold. A sub-area incomplete on a covariate is left out of the
# prediction universe on purpose and the sample households in it are dropped
# from the sub-area fits, so a small gap is expected; a matrix keyed by another
# code system, or built for other sub-areas, matches few or none of the frame's.
.sim_max_unmatched_share <- 0.05

# The sub-area codes of the frame's households must be found in the matrix.
# Reports the gap; stops when the households in missing codes exceed the share
# above (all of them, when the matrix matches no code).
.sim_check_subarea_cover <- function(frame, pop_sub, where) {
  codes <- data.table::data.table(
    code = subarea_code(frame$area, frame$manzana, frame$seccion))[, .N, by = code]
  miss <- codes[is.na(code) | !code %in% pop_sub$code]
  if (!nrow(miss)) return(invisible(NULL))
  share <- sum(miss$N) / nrow(frame)
  what <- paste0(.sim_np(nrow(miss), "sub-area code", "sub-area codes"), " of the frame (holding ",
                 sprintf("%.2f%%", 100 * share), " of its households)")
  if (share > .sim_max_unmatched_share) {
    .sim_stop(where, ": ", what, " ", .sim_isare(nrow(miss)), " not in the matrix's prediction ",
              "universe, more than the ", sprintf("%g%%", 100 * .sim_max_unmatched_share),
              " of households allowed; check that the matrix is keyed by the frame's sub-area codes ",
              "(the 22-character block code in urban areas, the 11-character rural section code ",
              "elsewhere) and that its covariates are complete")
  }
  message(where, ": ", what, " ", .sim_isare(nrow(miss)), " not in the prediction universe; ",
          "sampled households in them are left out of the sub-area fits")
  invisible(NULL)
}

# Row-standardised queens-contiguity weights as `build_queens_W()` caches them:
# a `municipality` column, then one weight column per municipality in the same
# order.
.sim_read_W <- function(path, universe) {
  df <- as.data.frame(arrow::read_parquet(path))
  ids <- df$municipality
  w_cols <- setdiff(names(df), "municipality")
  if (!identical(w_cols, as.character(ids))) {
    .sim_stop("spatial weights (", path, "): the weight columns do not follow the row order")
  }
  n_miss <- length(setdiff(universe, ids))
  if (n_miss) {
    .sim_stop("spatial weights (", path, "): ", .sim_np(n_miss, "municipality", "municipalities"),
              " of the universe ", .sim_isare(n_miss), " missing")
  }
  W <- as.matrix(df[, w_cols, drop = FALSE])
  dimnames(W) <- list(as.character(ids), as.character(ids))
  W
}

# MD5 of an R object, from its serialisation.
.sim_md5_object <- function(x) {
  tmp <- tempfile()
  on.exit(unlink(tmp), add = TRUE)
  con <- file(tmp, "wb")
  serialize(x, con, version = 3L)
  close(con)
  unname(tools::md5sum(tmp))
}

# Where an input came from: the path as given, its size in bytes and the MD5 of
# its content; a data frame given in memory is recorded as "<data frame>" with
# the MD5 of its serialisation. NULL for an input that was not given.
.sim_source <- function(x) {
  if (is.null(x)) return(NULL)
  if (is.data.frame(x)) {
    return(list(path = "<data frame>", size = NA_real_,
                md5 = .sim_md5_object(as.data.frame(x))))
  }
  list(path = x, size = unname(as.numeric(file.size(x))), md5 = unname(tools::md5sum(x)))
}

# Tuned XGB hyperparameters of a sub-area configuration, or NULL when none ship
# for that label (a user matrix under its own label).
.sim_tune <- function(config) {
  f <- system.file("extdata", sprintf("xgb_sub_params_%s.tsv", config), package = "saecol2026")
  if (!nzchar(f)) {
    message("no tuned XGB hyperparameters ship for sub-area configuration '", config,
            "'; `tune` is NULL (tune with data-raw/tune_xgb_sub.R)")
    return(NULL)
  }
  xgb_sub_params(config)
}

#' The sampling inputs of the simulation
#'
#' Builds what [draw_replicate()] and [apply_calibration()] take: the household
#' frame with its ECV stratum, the design, the per-municipality segments and the
#' per-stratum household and person totals.
#'
#' Households of municipalities that have no stratum in the stratum map cannot be
#' drawn and are left out of the frame; the municipalities stay in `universe`
#' (every municipality of the frame as read), which is the set that is estimated
#' and scored. The stratum is stored as text. Segments of about ten dwellings
#' are formed once per municipality (the first-stage unit) from the frame, so
#' that a replicate only selects among them. `pop_totals` holds, for each
#' stratum sorted by label, the number of households `N` and the number of
#' persons `N_persons` (households without a recorded size count zero), the
#' margins the calibration rakes the design weights to.
#'
#' @param data_dir data root; reads `frame/national_household_frame.parquet` and
#'   `aux_tables/muni_stratum_nn_extended.parquet`
#' @return list: `frame` (data frame), `design` (from [parse_ecv_design()]),
#'   `segments`, `pop_totals` (data frame `stratum, N, N_persons`) and
#'   `universe` (sorted integer municipality codes)
#' @export
sim_sampling_inputs <- function(data_dir = output_path()) {
  fp <- file.path(data_dir, "frame", "national_household_frame.parquet")
  mp <- file.path(data_dir, "aux_tables", "muni_stratum_nn_extended.parquet")
  .sim_need_file(fp, "household frame", "build it with data-raw/build_frame.R")
  .sim_need_file(mp, "stratum map", "build it with data-raw/build_nn_strata.R")

  fr  <- data.table::as.data.table(arrow::read_parquet(fp))
  m2s <- data.table::as.data.table(arrow::read_parquet(mp))
  need_fr <- union(c("municipality", "area", "seccion", "manzana", "hh_id", "n_persons"),
                   stats::na.omit(c(sae_outcomes$indicator, sae_outcomes$numerator,
                                    sae_outcomes$denominator)))
  absent <- setdiff(need_fr, names(fr))
  if (length(absent)) {
    .sim_stop("household frame (", fp, "): missing column(s): ", paste(absent, collapse = ", "))
  }
  absent <- setdiff(c("municipality", "stratum"), names(m2s))
  if (length(absent)) {
    .sim_stop("stratum map (", mp, "): missing column(s): ", paste(absent, collapse = ", "))
  }
  fr[, municipality := .sim_as_integer_key(municipality, "municipality",
                                           paste0("household frame (", fp, ")"))]
  m2s[, municipality := .sim_as_integer_key(municipality, "municipality",
                                            paste0("stratum map (", mp, ")"))]
  if (anyDuplicated(m2s$municipality)) {
    .sim_stop("stratum map (", mp, "): ", .sim_np(sum(duplicated(m2s$municipality)),
                                                  "municipality appears", "municipalities appear"),
              " more than once")
  }

  universe <- sort(unique(fr$municipality))
  n_unstrat <- length(setdiff(universe, m2s$municipality))
  if (n_unstrat) {
    message("household frame: ", .sim_np(n_unstrat, "municipality has", "municipalities have"),
            " no stratum in the stratum map; ", if (n_unstrat == 1L) "its" else "their",
            " households are not drawn")
  }

  # Frame restricted to stratified municipalities, with the stratum attached.
  fr <- fr[municipality %in% m2s$municipality]
  fr <- merge(fr, m2s[, .(municipality, stratum)], by = "municipality", all.x = TRUE)
  fr[, stratum := as.character(stratum)]

  design <- parse_ecv_design()
  frame <- as.data.frame(fr)
  segments <- build_all_segments(frame, psu_var = "municipality",
                                 target_size = .null_or(design$stage2$target_segment_size, 10L))

  # Calibration margins: households and persons per stratum.
  pt <- fr[, .(N = .N, N_persons = sum(n_persons, na.rm = TRUE)), by = stratum][order(stratum)]
  pop_totals <- data.frame(stratum = pt$stratum, N = pt$N, N_persons = pt$N_persons,
                           stringsAsFactors = FALSE)
  list(frame = frame, design = design, segments = segments,
       pop_totals = pop_totals, universe = universe)
}

#' Read and check every input of the replicate runner
#'
#' One call returns what a replicate needs. The universe is the set of
#' municipalities of the household frame. The area-level matrix is
#' `aux_tables/municipal_X_<area_config>.parquet` and the sub-area matrix
#' `subarea/X_<subarea_config>.parquet`, unless `area_matrix` / `subarea_matrix`
#' name a file or give a data frame in their place; each is validated
#' ([sim_validate_area_matrix()], [sim_validate_subarea_matrix()]), then
#' `extra_area` / `extra_subarea` are joined on `municipality` / `code`
#' ([sim_add_covariates()]) and the result is validated again. Every municipality
#' of the universe must keep at least one sub-area, because the sub-area methods
#' estimate a municipality as a weighted mean over its sub-areas. The sub-area
#' codes of the frame's households must also be found in the matrix: households
#' in codes it lacks (sub-areas incomplete on a covariate) are reported, and more
#' than 5\% of the frame's households in missing codes stops the run. The truth
#' must have no missing or non-finite `theta_true`, no repeated (`area_id`,
#' `outcome`) pair and a row for every municipality of the universe and every
#' outcome.
#'
#' The spatial weights of `fh_spatial` are read from
#' `aux_tables/W_queens_municipal.parquet` when that file exists and are `NULL`
#' otherwise. `tune` holds the per-outcome XGB hyperparameters shipped for
#' `subarea_config` ([xgb_sub_params()]) and is `NULL`, with a message, when none
#' ship for that label.
#'
#' `sources` records where the frame, the stratum map, the two matrices, the
#' extra covariates and the spatial weights came from (path as given, size in
#' bytes, MD5 of the content; `NULL` for an input not given), so that a run can
#' state what its results were computed from (see [sim_run()]).
#'
#' @param data_dir data root (default `output_path()`, from `SAECOL_DATA_DIR`)
#' @param area_config label of the area-level covariate set; names the built file
#' @param subarea_config label of the sub-area covariate set; names the built file
#'   and selects the XGB hyperparameters
#' @param area_matrix path (parquet or CSV) or data frame used in place of the
#'   built area-level matrix
#' @param subarea_matrix path or data frame used in place of the built sub-area
#'   matrix
#' @param extra_area path or data frame of extra area-level covariates, keyed by
#'   `municipality`
#' @param extra_subarea path or data frame of extra sub-area covariates, keyed by
#'   `code` (or `subarea`)
#' @return list: `frame`, `design`, `segments`, `pop_totals` (see
#'   [sim_sampling_inputs()]), `truth` (`area_id, theta_true, n_total, outcome`),
#'   `X_area`, `universe`, `W`, `pop_sub`, `tune`, `area_config`,
#'   `subarea_config`, `sources`
#' @export
sim_inputs <- function(data_dir = output_path(), area_config = "geo_bld",
                       subarea_config = "geo_blk_bld", area_matrix = NULL,
                       subarea_matrix = NULL, extra_area = NULL, extra_subarea = NULL) {
  samp <- sim_sampling_inputs(data_dir)
  universe <- samp$universe
  truth <- .sim_read_truth(file.path(data_dir, "aux_tables", "truth_full.parquet"), universe)

  # Area-level matrix.
  a_in <- if (is.null(area_matrix)) {
    file.path(data_dir, "aux_tables", sprintf("municipal_X_%s.parquet", area_config))
  } else area_matrix
  a <- .sim_read_input(a_in, "area covariate matrix",
                       hint = if (is.null(area_matrix)) {
                         "build it with data-raw/build_area_geo_bld.R or pass `area_matrix`"
                       })
  X_area <- sim_validate_area_matrix(a$data, universe, a$source)
  if (!is.null(extra_area)) {
    e <- .sim_read_input(extra_area, "extra area covariates")
    X_area <- sim_add_covariates(X_area, e$data, "municipality", e$source)
    src <- c(a$source, e$source)
    X_area <- sim_validate_area_matrix(X_area, universe,
                                       if (length(src)) paste(src, collapse = " + "))
  }

  # Sub-area matrix and its prediction universe.
  s_in <- if (is.null(subarea_matrix)) {
    file.path(data_dir, "subarea", sprintf("X_%s.parquet", subarea_config))
  } else subarea_matrix
  s <- .sim_read_input(s_in, "sub-area matrix",
                       hint = if (is.null(subarea_matrix)) {
                         "build it with data-raw/build_subarea_geo_blk_bld.R or pass `subarea_matrix`"
                       })
  pop_sub <- sim_validate_subarea_matrix(s$data, s$source)
  s_src <- s$source
  if (!is.null(extra_subarea)) {
    e <- .sim_read_input(extra_subarea, "extra sub-area covariates")
    pop_sub <- sim_add_covariates(pop_sub, .sim_subarea_key(e$data), "code", e$source)
    s_src <- if (length(c(s$source, e$source))) paste(c(s$source, e$source), collapse = " + ")
    pop_sub <- sim_validate_subarea_matrix(pop_sub, s_src)
  }
  where <- .sim_where("sub-area matrix", s_src)
  outside <- setdiff(unique(pop_sub$municipality), universe)
  if (length(outside)) {
    .sim_stop(where, ": ", .sim_np(length(outside), "municipality", "municipalities"),
              " of the matrix ", .sim_isare(length(outside)), " not in the universe",
              .sim_first(outside))
  }
  lacking <- setdiff(universe, pop_sub$municipality)
  if (length(lacking)) {
    .sim_stop(where, ": ", .sim_np(length(lacking), "municipality of the universe has",
                                   "municipalities of the universe have"),
              " no sub-area complete on the covariates", .sim_first(lacking))
  }

  .sim_check_subarea_cover(samp$frame, pop_sub, where)

  w_path <- file.path(data_dir, "aux_tables", "W_queens_municipal.parquet")
  W <- if (file.exists(w_path)) .sim_read_W(w_path, universe) else NULL

  sources <- list(
    frame = .sim_source(file.path(data_dir, "frame", "national_household_frame.parquet")),
    stratum_map = .sim_source(file.path(data_dir, "aux_tables",
                                        "muni_stratum_nn_extended.parquet")),
    area_matrix = .sim_source(a_in),
    subarea_matrix = .sim_source(s_in),
    extra_area = .sim_source(extra_area),
    extra_subarea = .sim_source(extra_subarea),
    spatial_weights = if (file.exists(w_path)) .sim_source(w_path))

  list(frame = samp$frame, design = samp$design, segments = samp$segments,
       pop_totals = samp$pop_totals, truth = truth, X_area = X_area,
       universe = universe, W = W, pop_sub = pop_sub, tune = .sim_tune(subarea_config),
       area_config = area_config, subarea_config = subarea_config, sources = sources)
}
