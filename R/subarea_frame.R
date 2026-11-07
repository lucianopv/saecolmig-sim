#' Census block and rural section geography for the sub-area methods
#'
#' The area-level FH family models municipalities. The povmap sub-area methods
#' model a finer unit, the census block or rural section, and aggregate up to
#' the 1,122 municipalities.
#'
#' A sub-area is
#'
#'   urban -> `manzana`, the 22-character `COD_DANE_ANM` census block
#'   rural -> `substr(seccion, 1, 11)`, the first 11 characters of the
#'            14-character `seccion` key: department, municipality, class,
#'            rural sector and rural section
#'
#' Rural blocks are not meaningful units in dispersed areas, which is why the
#' rural side is carried at rural section rather than block resolution. The
#' two code families cannot collide: urban codes are 22 characters and rural
#' codes 11.
#'
#' @name subarea_frame
NULL

#' Sub-area code for a set of frame or sample rows
#'
#' @param area character vector, "urban" or "rural"
#' @param manzana character vector, the 22-char COD_DANE_ANM block code
#' @param seccion character vector, the 14-char `seccion` key (the first 14
#'   characters of the block code)
#' @return character vector of sub-area codes, NA where `area` is neither
#'   "urban" nor "rural"
#' @export
subarea_code <- function(area, manzana, seccion) {
  stopifnot(length(area) == length(manzana), length(area) == length(seccion))
  out <- rep(NA_character_, length(area))
  u <- !is.na(area) & area == "urban"
  r <- !is.na(area) & area == "rural"
  out[u] <- as.character(manzana)[u]
  out[r] <- substr(as.character(seccion)[r], 1L, 11L)
  out
}

#' Build the sub-area census population frame
#'
#' Aggregates the national household frame to sub-area. This is the
#' `pop_data` side of the povmap call and the prediction universe of every
#' sub-area method: it must cover all sub-areas, including those in
#' municipalities a replicate does not draw.
#'
#' @param frame household-level national frame (data.frame/data.table) with
#'   `municipality`, `area`, `manzana`, `seccion`, `n_persons`
#' @return data.table, one row per sub-area, with `code`, `municipality`,
#'   `stratum` ("urban"/"rural"), `n_households_census` and
#'   `n_persons_census`, keyed and sorted by `code`
#' @export
build_subarea_pop <- function(frame) {
  fr <- data.table::as.data.table(frame)
  required <- c("municipality", "area", "manzana", "seccion", "n_persons")
  missing  <- setdiff(required, colnames(fr))
  if (length(missing)) {
    stop("build_subarea_pop: frame is missing ", paste(missing, collapse = ", "))
  }
  fr <- fr[, c(required), with = FALSE]
  fr[, code := subarea_code(area, manzana, seccion)]
  if (anyNA(fr$code)) {
    stop("build_subarea_pop: ", sum(is.na(fr$code)),
         " rows have an area that is neither 'urban' nor 'rural'")
  }
  pop <- fr[, .(
    municipality        = data.table::first(municipality),
    stratum             = data.table::first(area),
    n_households_census = .N,
    n_persons_census    = sum(n_persons, na.rm = TRUE)
  ), by = code]

  # A code carries its municipality and urban/rural class in its own prefix,
  # so neither may vary within a code; a violation would mis-assign households
  # to the wrong domain.
  chk <- fr[, .(n_muni = data.table::uniqueN(municipality),
                n_area = data.table::uniqueN(area)), by = code]
  if (any(chk$n_muni > 1L) || any(chk$n_area > 1L)) {
    stop("build_subarea_pop: ", sum(chk$n_muni > 1L | chk$n_area > 1L),
         " sub-area codes span more than one municipality or area class")
  }
  data.table::setkey(pop, code)
  pop[]
}

#' Collapse a replicate sample to sub-area aggregates
#'
#' Produces the `smp_data` side of the povmap call: one design-weighted row per
#' sub-area realised in the replicate.
#'
#' Sample rows are matched to sub-areas through the population frame rather
#' than recomputed: a row joins on `manzana` if that block is a known urban
#' sub-area, otherwise on the first 11 characters of `seccion` if that is a
#' known rural one. Rows matching neither are dropped, and their count is
#' returned in the result's attributes.
#'
#' Household-basis outcomes give the weighted share of households; person-basis
#' outcomes give the ratio of weighted person counts, matching `outcome_spec()`.
#'
#' @param sample_dt one replicate sample, with `manzana`, `seccion`,
#'   `weight_cal` and the outcome columns
#' @param pop sub-area population frame from `build_subarea_pop()`
#' @param outcome outcome name from the outcome registry
#' @return data.table with `code`, `municipality`, `y`, `wt`, `n_hh`;
#'   attribute `n_unmatched` records rows that joined to no sub-area
#' @export
collapse_sample_to_subarea <- function(sample_dt, pop, outcome) {
  s   <- data.table::as.data.table(sample_dt)
  pop <- data.table::as.data.table(pop)
  spec <- outcome_spec(outcome)

  ucodes <- pop[stratum == "urban", code]
  rcodes <- pop[stratum == "rural", code]
  s[, code := NA_character_]
  s[as.character(manzana) %in% ucodes, code := as.character(manzana)]
  s[is.na(code) & substr(as.character(seccion), 1L, 11L) %in% rcodes,
    code := substr(as.character(seccion), 1L, 11L)]
  n_unmatched <- sum(is.na(s$code))
  s <- s[!is.na(code)]

  if (identical(spec$level, "person")) {
    # Person-basis outcomes are additive count ratios: weight the numerator
    # and the denominator, never the per-household ratio.
    sp <- s[n_pers_resp > 0]
    sa <- sp[, .(y  = sum(weight_cal * get(spec$numerator)) /
                      sum(weight_cal * get(spec$denominator)),
                 wt = sum(weight_cal * get(spec$denominator)),
                 n_hh = .N), by = code]
  } else {
    sa <- s[, .(y  = sum(weight_cal * get(spec$indicator)) / sum(weight_cal),
                wt = sum(weight_cal),
                n_hh = .N), by = code]
  }
  sa <- sa[is.finite(y) & is.finite(wt) & wt > 0]
  sa <- merge(sa, pop[, .(code, municipality)], by = "code")
  data.table::setattr(sa, "n_unmatched", n_unmatched)
  sa[]
}
