#' Assemble sub-area polygons for one department
#'
#' One definition of the sub-area polygons, shared by every pass that measures
#' them. Sub-areas come from three MGN sources with different key columns, and
#' centros poblados need dissolving because several urban sections share one
#' 11-character code. Every pass reads its polygons here, so anything divided
#' by their areas (population density) is consistent with the extraction.
#'
#' s2 must be off before calling this: the repair runs under planar GEOS, which
#' tolerates defects that s2 rejects.
#'
#' @param dd two-character department code
#' @param pop sub-area population table with `code` and `cls`
#'   (`cls = substr(code, 6, 6)`: "1" urban manzana, "2" centro poblado,
#'   "3" rural section)
#' @param F_MZ,F_RS,F_US paths to the urban-manzana, rural-section and
#'   urban-section MGN geopackages
#' @return an sf data frame with `code` and `geom`, or NULL if the department
#'   has no matching polygons
#' @export
subarea_geometry <- function(dd, pop, F_MZ, F_RS, F_US) {
  lyr <- function(p) sf::st_layers(p)$name[1]
  want_u <- pop[cls == "1" & substr(code, 1, 2) == dd, code]
  want_3 <- pop[cls == "3" & substr(code, 1, 2) == dd, code]
  want_2 <- pop[cls == "2" & substr(code, 1, 2) == dd, code]
  parts <- list()

  if (length(want_u)) {
    g <- sf::st_read(F_MZ, quiet = TRUE, query = sprintf(
      'SELECT MANZ_CCNCT AS code, geom FROM "%s" WHERE DPTO_CCDGO = \'%s\'',
      lyr(F_MZ), dd))
    parts[[length(parts) + 1L]] <- g[g$code %in% want_u, ]
  }
  if (length(want_3)) {
    g <- sf::st_read(F_RS, quiet = TRUE, query = sprintf(
      'SELECT SECR_CCNCT AS code, geom FROM "%s" WHERE DPTO_CCDGO = \'%s\'',
      lyr(F_RS), dd))
    parts[[length(parts) + 1L]] <- g[g$code %in% want_3, ]
  }
  if (length(want_2)) {
    # Centros poblados: several urban sections can share one 11-character code,
    # so they must be dissolved before extraction or the same sub-area would be
    # measured more than once.
    g <- sf::st_read(F_US, quiet = TRUE, query = sprintf(
      'SELECT SECU_CCNCT, geom FROM "%s" WHERE DPTO_CCDGO = \'%s\'',
      lyr(F_US), dd))
    if (nrow(g)) {
      g$code <- substr(g$SECU_CCNCT, 1, 11)
      g <- g[g$code %in% want_2, c("code", attr(g, "sf_column"))]
      if (nrow(g)) g <- stats::aggregate(g, by = list(code = g$code),
                                         FUN = function(z) z[1], do_union = TRUE)
      g <- g[, c("code", attr(g, "sf_column"))]
    }
    parts[[length(parts) + 1L]] <- g
  }
  parts <- Filter(function(z) !is.null(z) && nrow(z) > 0L, parts)
  if (!length(parts)) return(NULL)
  # st_read yields "geom"; aggregate() renames the geometry column, and rbind()
  # of sf objects requires all parts to agree on it, so normalise it.
  parts <- lapply(parts, function(z) {
    z <- z[, "code"]
    sf::st_geometry(z) <- "geom"
    z
  })
  out <- do.call(rbind, parts)
  out <- out[!duplicated(out$code), ]

  # MGN polygons carry minor topology defects (duplicate vertices within a
  # ring) which planar GEOS tolerates and s2 rejects ("Loop 0 is not valid").
  # Repair under GEOS, where st_make_valid is robust, then hand valid geometry
  # back to s2 for the area computation. The repair runs on every polygon:
  # geometry that GEOS reports as valid can still fail the s2 conversion, and
  # st_make_valid leaves valid input unchanged.
  out <- sf::st_make_valid(out)
  # make_valid can turn a polygon into a collection (polygon + stray line);
  # keep only the polygonal part, which is what st_area and the extractors use.
  out <- suppressWarnings(sf::st_collection_extract(out, "POLYGON"))
  out <- out[!sf::st_is_empty(out), ]
  out
}

#' Merge the polygons of a sub-area code into one
#'
#' The validity repair in [subarea_geometry()] can split one sub-area into
#' several polygons, which then arrive as several rows with the same code.
#' Measuring each row separately would count that sub-area more than once (and
#' give it several values), so the rows of a code are unioned into one polygon
#' before any extraction. Input with one row per code is returned unchanged.
#'
#' @param g sf data frame with a `code` column
#' @return `g` with one row per code; attribute `n_duplicate` holds the number
#'   of rows merged away
#' @export
subarea_union_codes <- function(g) {
  nd <- sum(duplicated(g$code))
  if (nd) {
    g <- stats::aggregate(g, by = list(code_u = g$code), FUN = function(z) z[1],
                          do_union = TRUE)
    g <- g[, "code"]
  }
  attr(g, "n_duplicate") <- nd
  g
}
