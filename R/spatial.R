#' Build a same-department adjacency matrix for a set of municipalities
#'
#' A coarse alternative to the queens-contiguity matrix of [build_queens_W()]:
#' it needs no polygons, only the municipality codes. The runner does not use
#' it.
#'
#' Each municipality (DIVIPOLA dept * 1000 + muni) is connected with equal
#' weight to every other municipality sharing its department. Singletons (the
#' only muni in their dept among `area_ids`) get a zero row. The result is
#' row-normalised so each non-empty row sums to 1.
#'
#' @param area_ids integer vector of DIVIPOLA municipal codes
#' @return a sparse Matrix (Matrix::Matrix) of dimension n x n
#' @keywords internal
build_dept_adjacency <- function(area_ids) {
  if (!requireNamespace("Matrix", quietly = TRUE)) {
    stop("Package 'Matrix' is required for build_dept_adjacency().")
  }
  n     <- length(area_ids)
  depts <- area_ids %/% 1000L
  ids_c <- as.character(area_ids)

  W <- Matrix::Matrix(0, nrow = n, ncol = n, sparse = TRUE,
                      dimnames = list(ids_c, ids_c))
  for (d in unique(depts)) {
    members <- which(depts == d)
    if (length(members) > 1L) {
      m <- length(members)
      block <- matrix(1 / (m - 1L), nrow = m, ncol = m)
      diag(block) <- 0
      W[members, members] <- block
    }
  }
  W
}

#' Build a row-standardised queens-contiguity W matrix for Colombian municipalities
#'
#' Reads `colmaps::municipios` (a `SpatialPolygonsDataFrame` with 1,122
#' municipalities and DIVIPOLA codes in `id`), converts to `sf`, computes
#' queens contiguity via `spdep::poly2nb(queen = TRUE)`, row-standardises
#' via `spdep::nb2mat(style = "W", zero.policy = TRUE)`, and caches at
#' `aux_tables/W_queens_municipal.parquet` under `SAECOL_DATA_DIR`.
#'
#' Isolated municipalities (San Andres y Providencia, dept 88) get an all-zero
#' row; `fh_spatial` predicts them synthetically from `X * beta`.
#'
#' @param force_recompute if TRUE, ignore the on-disk cache and rebuild
#' @return a 1,122 x 1,122 numeric matrix with rownames/colnames =
#'   character DIVIPOLA codes (sorted ascending)
#' @export
build_queens_W <- function(force_recompute = FALSE) {
  out_path <- file.path(output_path("aux_tables", create = TRUE),
                       "W_queens_municipal.parquet")

  if (!force_recompute && file.exists(out_path)) {
    df <- read_parquet_safe(out_path)
    ids <- df$municipality
    weight_cols <- setdiff(colnames(df), "municipality")
    # Guard the cache round-trip: the parquet column order must equal the
    # ascending DIVIPOLA codes of the `municipality` column; a cache written in
    # another order stops here instead of mis-indexing the matrix.
    stopifnot(identical(weight_cols, as.character(ids)))
    W <- as.matrix(df[, weight_cols, drop = FALSE])
    rownames(W) <- as.character(ids)
    colnames(W) <- as.character(ids)
    return(W)
  }

  if (!requireNamespace("colmaps", quietly = TRUE)) {
    stop("Package 'colmaps' is required for build_queens_W(). ",
         "Install via: remotes::install_github('nebulae-co/colmaps')")
  }
  if (!requireNamespace("sf", quietly = TRUE)) {
    stop("Package 'sf' is required for build_queens_W().")
  }
  if (!requireNamespace("spdep", quietly = TRUE)) {
    stop("Package 'spdep' is required for build_queens_W().")
  }

  m <- sf::st_as_sf(colmaps::municipios)
  m$id <- as.integer(as.character(m$id))
  m <- m[order(m$id), ]

  nb <- spdep::poly2nb(m, queen = TRUE)
  W  <- spdep::nb2mat(nb, style = "W", zero.policy = TRUE)

  rownames(W) <- as.character(m$id)
  colnames(W) <- as.character(m$id)

  # Diagnostics (no row-level data; just summary stats).
  rs <- rowSums(W)
  n_isolated <- sum(abs(rs - 0) < 1e-8)
  cat(sprintf("build_queens_W: %d munis, %d isolated, ",
              nrow(W), n_isolated))
  card <- sapply(nb, length)
  cat(sprintf("median n_neighbours = %.0f, IQR = [%.0f, %.0f]\n",
              stats::median(card), stats::quantile(card, 0.25),
              stats::quantile(card, 0.75)))

  # Cache as parquet (one column per neighbour weight + an id column).
  df <- as.data.frame(W)
  df <- cbind(municipality = as.integer(rownames(W)), df)
  write_parquet_safe(df, out_path)

  W
}
