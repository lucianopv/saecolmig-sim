#' Census composition covariates at sub-area resolution
#'
#' Computes the same nine census shares as [build_dept_aux()] (all but
#' `share_urban`, which is constant within a sub-area) but grouped by sub-area
#' (census block or rural section) instead of municipality. The variable list is
#' identical to the municipal matrix, so the
#' sub-area methods differ from the area-level methods in resolution and
#' method only, not in which covariates they were given.
#'
#' Records reach a sub-area through the MGN geographic file on the same
#' five-column key the household frame uses
#' (`U_DPTO`, `U_MPIO`, `UA_CLASE`, `COD_ENCUESTAS`, `U_VIVIENDA`), so a
#' household lands in the same block here and in `build_subarea_pop()`.
#'
#' @name subarea_census
NULL

# Person, household and dwelling columns needed for the nine shares.
.subarea_per_cols <- c("U_DPTO", "U_MPIO", "UA_CLASE", "COD_ENCUESTAS",
                       "U_VIVIENDA", "P_EDADR", "P_ALFABETA",
                       "P_NIVEL_ANOSR", "P_TRABAJO", "PA1_GRP_ETNIC")
.subarea_hog_cols <- c("U_DPTO", "U_MPIO", "UA_CLASE", "COD_ENCUESTAS",
                       "U_VIVIENDA", "H_AGUA_COCIN")
.subarea_viv_cols <- c("U_DPTO", "U_MPIO", "UA_CLASE", "COD_ENCUESTAS",
                       "U_VIVIENDA", "VA_EE")
.subarea_mgn_cols <- c("U_DPTO", "U_MPIO", "UA_CLASE", "COD_ENCUESTAS",
                       "U_VIVIENDA", "COD_DANE_ANM")

.mgn_key <- c("U_DPTO", "U_MPIO", "UA_CLASE", "COD_ENCUESTAS", "U_VIVIENDA")

# Attach COD_DANE_ANM and derive the sub-area code. Urban records carry the
# whole 22-character block code; rural records are held at rural section
# resolution,
# which is substr(., 1, 11) of the same string (the frame's `seccion` is
# substr(., 1, 14), so the two truncations agree).
.attach_subarea_code <- function(dt, mgn) {
  out <- merge(dt, mgn, by = .mgn_key, all.x = TRUE)
  out[, manzana := as.character(COD_DANE_ANM)]
  out[, code := data.table::fifelse(UA_CLASE == 1L, manzana,
                                    substr(manzana, 1L, 11L))]
  out[!is.na(code) & nzchar(code)]
}

#' Build sub-area census covariates for one department
#'
#' @param dept_code two-digit DANE department code
#' @param n_max optional row cap on the person file, for smoke tests
#' @return data.table, one row per sub-area
#' @export
build_dept_census_subarea <- function(dept_code, n_max = NA_integer_) {
  dept <- sprintf("%02d", as.integer(dept_code))
  dirs <- list.files(saecol2026::census_path(), pattern = paste0("^", dept),
                     full.names = TRUE)
  dirs <- dirs[file.info(dirs)$isdir]
  if (!length(dirs)) stop("Department folder not found for code ", dept_code)
  d <- dirs[1]

  f <- function(pat) {
    p <- list.files(d, pattern = pat, recursive = TRUE, full.names = TRUE)[1]
    if (is.na(p)) stop("Missing ", pat, " in ", d)
    p
  }
  mgn <- data.table::fread(f("^CNPV2018_MGN_.*\\.CSV$"),
                           select = .subarea_mgn_cols, showProgress = FALSE)
  per <- data.table::fread(f("^CNPV2018_5PER_.*\\.CSV$"),
                           select = .subarea_per_cols,
                           nrows = if (is.na(n_max)) Inf else n_max,
                           showProgress = FALSE)
  hog <- data.table::fread(f("^CNPV2018_2HOG_.*\\.CSV$"),
                           select = .subarea_hog_cols, showProgress = FALSE)
  viv <- data.table::fread(f("^CNPV2018_1VIV_.*\\.CSV$"),
                           select = .subarea_viv_cols, showProgress = FALSE)

  per <- .attach_subarea_code(per, mgn)
  hog <- .attach_subarea_code(hog, mgn)
  viv <- .attach_subarea_code(viv, mgn)

  # Same definitions as build_dept_aux(); only the grouping differs.
  per_agg <- per[, .(
    n_persons_sub             = .N,
    share_working_age         = mean(P_EDADR >= 4L & P_EDADR <= 13L, na.rm = TRUE),
    share_elderly             = mean(P_EDADR >= 14L, na.rm = TRUE),
    share_secondary_completed = mean(P_NIVEL_ANOSR >= 7L & P_EDADR >= 6L, na.rm = TRUE),
    share_illiterate          = mean(P_ALFABETA == 2L & P_EDADR >= 4L, na.rm = TRUE),
    share_economically_active = mean(P_TRABAJO == 1L, na.rm = TRUE),
    share_afro                = mean(PA1_GRP_ETNIC == 5L, na.rm = TRUE),
    share_indigenous          = mean(PA1_GRP_ETNIC == 1L, na.rm = TRUE)
  ), by = code]

  hog_agg <- hog[, .(n_households_sub = .N,
                     share_piped_water = mean(H_AGUA_COCIN == 1L, na.rm = TRUE)),
                 by = code]
  viv_agg <- viv[, .(share_electricity = mean(VA_EE == 1L, na.rm = TRUE)),
                 by = code]

  out <- merge(per_agg, hog_agg, by = "code", all = TRUE)
  out <- merge(out,     viv_agg, by = "code", all = TRUE)
  out[, dept := as.integer(dept)]

  # A share over zero records is NaN, which is not a share of zero. Normalise
  # it to NA so that complete-case filtering sees it as missing.
  for (cc in setdiff(colnames(out), c("code", "dept"))) {
    data.table::set(out, which(is.nan(out[[cc]])), cc, NA_real_)
  }
  out[]
}

#' Build sub-area census covariates for the whole country
#'
#' @param dept_codes department codes; defaults to every unzipped department
#' @param verbose print per-department progress
#' @return data.table, one row per sub-area nationally
#' @export
build_census_subarea <- function(dept_codes = NULL, verbose = TRUE) {
  if (is.null(dept_codes)) {
    e  <- list.files(saecol2026::census_path(), pattern = "^[0-9]{2}")
    e  <- e[file.info(file.path(saecol2026::census_path(), e))$isdir]
    dept_codes <- sort(unique(substr(e, 1, 2)))
  }
  parts <- vector("list", length(dept_codes))
  for (i in seq_along(dept_codes)) {
    t0 <- Sys.time()
    parts[[i]] <- build_dept_census_subarea(dept_codes[i])
    if (verbose) {
      cat(sprintf("  [%2d/%2d] dept %s: %7d sub-areas (%.1f min)\n",
                  i, length(dept_codes), dept_codes[i], nrow(parts[[i]]),
                  as.numeric(difftime(Sys.time(), t0, units = "mins"))))
      utils::flush.console()
    }
  }
  data.table::rbindlist(parts, fill = TRUE)
}
