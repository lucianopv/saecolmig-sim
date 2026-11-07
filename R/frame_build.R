#' Derive migration flag and type columns from PA_VIVIA_1ANO
#'
#' DANE coding (CNPV 2018, person-level question
#'   "Hace un anio donde vivia esta persona?"):
#'   1 = aun no habia nacido (not yet born one year ago) - not a migrant
#'   2 = en este municipio (same municipality) - not a migrant
#'   3 = en otro municipio del pais (different Colombian municipality)
#'       - internal migrant
#'   4 = en otro pais (abroad) - international migrant
#'   9 = no informa - NA
#'
#' Persons with codes 3 or 4 are migrants (migrated_last_12m = 1);
#' codes 1 or 2 are non-migrants (0). Code 9 propagates to NA.
#' Any code outside 1, 2, 3, 4 and 9 raises an error.
#'
#' @param df data.frame/tibble with column PA_VIVIA_1ANO (character or integer)
#' @return df with three new integer columns:
#'   \describe{
#'     \item{migrated_last_12m}{1 if internal or international migrant, 0 if non-migrant, NA if unknown}
#'     \item{mig_internal}{1 if internal migrant (code 3), 0 if non-migrant or international, NA if unknown}
#'     \item{mig_intl}{1 if international migrant (code 4), 0 if non-migrant or internal, NA if unknown}
#'   }
#' @export
add_migrant_flag <- function(df) {
  v <- as.integer(df$PA_VIVIA_1ANO)
  bad <- !is.na(v) & !(v %in% c(1L, 2L, 3L, 4L, 9L))
  if (any(bad)) {
    stop("Unexpected PA_VIVIA_1ANO codes: ",
         paste(sort(unique(v[bad])), collapse = ", "))
  }
  known <- !is.na(v) & v %in% c(1L, 2L, 3L, 4L)
  df$migrated_last_12m <- ifelse(known, as.integer(v %in% c(3L, 4L)),
                                 NA_integer_)
  df$mig_internal      <- ifelse(known, as.integer(v == 3L), NA_integer_)
  df$mig_intl          <- ifelse(known, as.integer(v == 4L), NA_integer_)
  df
}

#' Aggregate person-level migrant flags to household counts + indicators
#'
#' Counts are over persons with known migration status; `n_mig_total =
#' n_mig_internal + n_mig_intl` by construction (person types are mutually
#' exclusive). Indicators use any-member semantics (a household with both
#' types counts in both `migrated_internal_12m` and `migrated_intl_12m`).
#' Households whose members are all NA get counts 0 / indicators 0, as
#' `any(..., na.rm = TRUE)` gives.
#'
#' @param per data.frame/data.table of persons after add_migrant_flag(),
#'   with geo keys U_DPTO, U_MPIO, UA_CLASE, COD_ENCUESTAS, U_VIVIENDA,
#'   P_NROHOG.
#' @return data.table, one row per household
#' @export
aggregate_persons_to_household <- function(per) {
  dt <- data.table::as.data.table(per)
  out <- dt[, .(
    n_pers_resp    = sum(!is.na(migrated_last_12m)),
    n_mig_internal = sum(mig_internal == 1L, na.rm = TRUE),
    n_mig_intl     = sum(mig_intl == 1L, na.rm = TRUE)
  ), by = .(U_DPTO, U_MPIO, UA_CLASE, COD_ENCUESTAS,
            U_VIVIENDA, P_NROHOG)]
  out[, n_mig_total           := n_mig_internal + n_mig_intl]
  out[, migrated_last_12m     := as.integer(n_mig_total > 0L)]
  out[, migrated_internal_12m := as.integer(n_mig_internal > 0L)]
  out[, migrated_intl_12m     := as.integer(n_mig_intl > 0L)]
  out[]
}

#' Read one department's CNPV 2018 person + household + geographic files
#' and return a household-level frame with the migrant flag and PSU/SSU keys
#'
#' Joins:
#'   PER (aggregated to household) -> HOG (one row per household) -> MGN (geo).
#'
#' Geographic derivation:
#'   seccion = first 14 chars of COD_DANE_ANM: department (1-2), municipality
#'   (3-5), class (6), rural sector (7-9), rural section (10-11) and centro
#'   poblado (12-14); the urban sector and section are characters 15-20.
#'   manzana = COD_DANE_ANM (the canonical unique manzana code).
#' For rural dispersed areas where COD_DANE_ANM does not encode a real manzana,
#' the same code is reused as both seccion and manzana grouping; the sampler
#' treats both as opaque identifiers, so this is consistent.
#'
#' @param dept_code character or integer dept code (e.g. "05")
#' @param n_max integer or NA; if not NA, limits PER rows read (testing aid)
#' @return tibble with one row per household
#' @export
read_dept_frame <- function(dept_code, n_max = NA_integer_) {
  dept <- sprintf("%02d", as.integer(dept_code))
  candidate_dirs <- list.files(census_path(),
                               pattern = paste0("^", dept), full.names = TRUE)
  candidate_dirs <- candidate_dirs[file.info(candidate_dirs)$isdir]
  if (length(candidate_dirs) == 0) {
    stop("Department folder not found for code ", dept_code,
         " under ", census_path(),
         ". Unzip the corresponding ZIP first.")
  }
  dept_dir <- candidate_dirs[1]

  per_csv <- list.files(dept_dir, pattern = "^CNPV2018_5PER_.*\\.CSV$",
                        recursive = TRUE, full.names = TRUE)[1]
  hog_csv <- list.files(dept_dir, pattern = "^CNPV2018_2HOG_.*\\.CSV$",
                        recursive = TRUE, full.names = TRUE)[1]
  mgn_csv <- list.files(dept_dir, pattern = "^CNPV2018_MGN_.*\\.CSV$",
                        recursive = TRUE, full.names = TRUE)[1]

  per_cols <- c("U_DPTO", "U_MPIO", "UA_CLASE", "COD_ENCUESTAS",
                "U_VIVIENDA", "P_NROHOG", "P_NRO_PER", "PA_VIVIA_1ANO")
  hog_cols <- c("U_DPTO", "U_MPIO", "UA_CLASE", "COD_ENCUESTAS",
                "U_VIVIENDA", "H_NROHOG", "HA_TOT_PER")
  mgn_cols <- c("U_DPTO", "U_MPIO", "UA_CLASE", "COD_ENCUESTAS",
                "U_VIVIENDA", "COD_DANE_ANM")

  per <- data.table::fread(per_csv, select = per_cols,
                           nrows = if (is.na(n_max)) Inf else n_max,
                           showProgress = FALSE)
  hog <- data.table::fread(hog_csv, select = hog_cols, showProgress = FALSE)
  mgn <- data.table::fread(mgn_csv, select = mgn_cols, showProgress = FALSE)

  per <- add_migrant_flag(per)
  per_hh <- aggregate_persons_to_household(per)

  fr <- merge(hog, per_hh,
              by.x = c("U_DPTO", "U_MPIO", "UA_CLASE", "COD_ENCUESTAS",
                       "U_VIVIENDA", "H_NROHOG"),
              by.y = c("U_DPTO", "U_MPIO", "UA_CLASE", "COD_ENCUESTAS",
                       "U_VIVIENDA", "P_NROHOG"),
              all.x = TRUE)

  fr <- merge(fr, mgn,
              by = c("U_DPTO", "U_MPIO", "UA_CLASE", "COD_ENCUESTAS",
                     "U_VIVIENDA"),
              all.x = TRUE)

  data.table::setnames(fr,
    c("U_DPTO", "U_MPIO", "UA_CLASE", "HA_TOT_PER"),
    c("dept", "muni", "area_code", "n_persons"))

  fr[, area := ifelse(area_code == 1L, "urban", "rural")]
  fr[, manzana := as.character(COD_DANE_ANM)]
  # seccion = first 14 chars of the block code: dept(2) + muni(3) + class(1) +
  # rural sector(3) + rural section(2) + centro poblado(3)
  fr[, seccion := substr(manzana, 1, 14)]
  # DIVIPOLA municipal code: 2-digit dept * 1000 + 3-digit muni (unique nationally)
  fr[, municipality := dept * 1000L + muni]
  fr[, hh_id := paste(COD_ENCUESTAS, U_VIVIENDA, H_NROHOG, sep = "_")]

  tibble::as_tibble(fr[, .(dept, muni, municipality, area,
                            seccion, manzana, hh_id,
                            n_persons,
                            migrated_last_12m, migrated_internal_12m,
                            migrated_intl_12m,
                            n_pers_resp, n_mig_total,
                            n_mig_internal, n_mig_intl)])
}

#' Build a national household-level frame by reading every department under census_path()
#'
#' @param dept_codes character vector or NULL (all available dept folders)
#' @return one tibble with all households
#' @export
build_national_frame <- function(dept_codes = NULL) {
  if (is.null(dept_codes)) {
    entries <- list.files(census_path(), pattern = "^[0-9]{2}")
    is_dir  <- file.info(file.path(census_path(), entries))$isdir
    entries <- entries[is_dir]
    dept_codes <- unique(substr(entries, 1, 2))
  }
  parts <- lapply(unique(dept_codes), read_dept_frame)
  do.call(rbind, parts)
}
