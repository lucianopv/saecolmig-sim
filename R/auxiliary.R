#' Build per-municipality auxiliary variables from one CNPV 2018 department
#'
#' Reads PER + HOG + VIV CSVs with the columns needed for demographic,
#' education, labour, housing, and ethnicity aggregates, then aggregates
#' to municipality.
#'
#' Indicator definitions (DANE codings, CNPV 2018):
#'   share_urban             : UA_CLASE == 1 (cabecera)
#'   share_working_age       : P_EDADR in [4, 13]  (5-yr groups: 15-19, ..., 60-64)
#'   share_elderly           : P_EDADR >= 14       (65+)
#'   share_secondary_completed: P_NIVEL_ANOSR >= 7 and P_EDADR >= 6 (25+),
#'                             over all persons
#'   share_illiterate        : P_ALFABETA == 2 and P_EDADR >= 4 (15+),
#'                             over all persons
#'   share_economically_active: P_TRABAJO == 1
#'   share_afro              : PA1_GRP_ETNIC == 5
#'   share_indigenous        : PA1_GRP_ETNIC == 1
#'   share_piped_water       : H_AGUA_COCIN == 1 (acueducto publico)
#'   share_electricity       : VA_EE == 1 (vivienda con servicio de energia)
#'
#' @param dept_code character or integer (e.g. "05")
#' @param n_max integer or NA; if not NA, limits PER rows read (testing aid)
#' @return tibble with one row per municipality
#' @export
build_dept_aux <- function(dept_code, n_max = NA_integer_) {
  dept <- sprintf("%02d", as.integer(dept_code))
  candidate_dirs <- list.files(census_path(),
                               pattern = paste0("^", dept), full.names = TRUE)
  candidate_dirs <- candidate_dirs[file.info(candidate_dirs)$isdir]
  if (length(candidate_dirs) == 0) {
    stop("Department folder not found for code ", dept_code)
  }
  dept_dir <- candidate_dirs[1]

  per_csv <- list.files(dept_dir, pattern = "^CNPV2018_5PER_.*\\.CSV$",
                        recursive = TRUE, full.names = TRUE)[1]
  hog_csv <- list.files(dept_dir, pattern = "^CNPV2018_2HOG_.*\\.CSV$",
                        recursive = TRUE, full.names = TRUE)[1]
  viv_csv <- list.files(dept_dir, pattern = "^CNPV2018_1VIV_.*\\.CSV$",
                        recursive = TRUE, full.names = TRUE)[1]

  per_cols <- c("U_DPTO", "U_MPIO", "UA_CLASE", "P_EDADR", "P_ALFABETA",
                "P_NIVEL_ANOSR", "P_TRABAJO", "PA1_GRP_ETNIC")
  hog_cols <- c("U_DPTO", "U_MPIO", "UA_CLASE", "H_AGUA_COCIN")
  viv_cols <- c("U_DPTO", "U_MPIO", "UA_CLASE", "VA_EE")

  per <- data.table::fread(file = per_csv, select = per_cols,
                           nrows = if (is.na(n_max)) Inf else n_max,
                           showProgress = FALSE)
  hog <- data.table::fread(file = hog_csv, select = hog_cols,
                           showProgress = FALSE)
  viv <- data.table::fread(file = viv_csv, select = viv_cols,
                           showProgress = FALSE)

  per_agg <- per[, .(
    n_persons                 = .N,
    share_urban               = mean(UA_CLASE == 1L, na.rm = TRUE),
    share_working_age         = mean(P_EDADR >= 4L & P_EDADR <= 13L, na.rm = TRUE),
    share_elderly             = mean(P_EDADR >= 14L, na.rm = TRUE),
    share_secondary_completed = mean(P_NIVEL_ANOSR >= 7L &
                                       P_EDADR >= 6L, na.rm = TRUE),
    share_illiterate          = mean(P_ALFABETA == 2L &
                                       P_EDADR >= 4L, na.rm = TRUE),
    share_economically_active = mean(P_TRABAJO == 1L, na.rm = TRUE),
    share_afro                = mean(PA1_GRP_ETNIC == 5L, na.rm = TRUE),
    share_indigenous          = mean(PA1_GRP_ETNIC == 1L, na.rm = TRUE)
  ), by = .(U_DPTO, U_MPIO)]

  hog_agg <- hog[, .(
    n_households       = .N,
    share_piped_water  = mean(H_AGUA_COCIN == 1L, na.rm = TRUE)
  ), by = .(U_DPTO, U_MPIO)]

  viv_agg <- viv[, .(
    share_electricity = mean(VA_EE == 1L, na.rm = TRUE)
  ), by = .(U_DPTO, U_MPIO)]

  out <- merge(per_agg, hog_agg, by = c("U_DPTO", "U_MPIO"), all.x = TRUE)
  out <- merge(out,     viv_agg, by = c("U_DPTO", "U_MPIO"), all.x = TRUE)

  data.table::setnames(out, c("U_DPTO", "U_MPIO"), c("dept", "muni"))
  out[, municipality := dept * 1000L + muni]
  tibble::as_tibble(out)
}

#' Build a national municipal X matrix
#'
#' Iterates over every department directory under census_path() and stacks
#' the per-dept aggregates from build_dept_aux(). The matrix holds the census
#' covariates only; other layers are joined with [bind_aux_layer()].
#'
#' @param dept_codes character vector or NULL (all available)
#' @return tibble with one row per municipality: `municipality` and the census
#'   covariates
#' @export
build_aux_municipal <- function(dept_codes = NULL) {
  if (is.null(dept_codes)) {
    entries <- list.files(census_path(), pattern = "^[0-9]{2}")
    is_dir  <- file.info(file.path(census_path(), entries))$isdir
    entries <- entries[is_dir]
    dept_codes <- unique(substr(entries, 1, 2))
  }
  parts <- lapply(unique(dept_codes), build_dept_aux)
  do.call(rbind, parts)
}

#' Left-join an extra auxiliary layer onto the municipal X matrix
#'
#' Used to add remote-sensing or other layers to aux without modifying the
#' core builder. The layer must contain a `municipality` column.
#'
#' @param aux tibble from build_aux_municipal()
#' @param layer tibble keyed by `municipality`
#' @return aux with new layer columns appended
#' @export
bind_aux_layer <- function(aux, layer) {
  if (!"municipality" %in% colnames(layer)) {
    stop("layer must contain a 'municipality' column.")
  }
  merge(aux, layer, by = "municipality", all.x = TRUE)
}
