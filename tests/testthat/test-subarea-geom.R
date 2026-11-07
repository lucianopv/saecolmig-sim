# subarea_geometry (R/subarea_geom.R) on small synthetic geopackages.

.sq <- function(x0, y0, w = 1) {
  sf::st_polygon(list(rbind(c(x0, y0), c(x0 + w, y0), c(x0 + w, y0 + w),
                            c(x0, y0 + w), c(x0, y0))))
}

# Three MGN sources for department "05", plus one block of department "08".
.fixture_mgn <- function(dir) {
  a <- paste0("050011", strrep("0", 15), "1")      # urban block, 22 characters
  b <- paste0("050011", strrep("0", 15), "2")      # urban block, bow-tie polygon
  r <- "05001300001"                               # rural section, 11 characters
  cp <- "05001200002"                              # centro poblado, 11 characters
  other <- paste0("080011", strrep("0", 15), "1")  # block of another department
  bowtie <- sf::st_polygon(list(rbind(c(0, 5), c(2, 7), c(2, 5), c(0, 7), c(0, 5))))
  mk <- function(geoms, ...) sf::st_sf(..., geom = sf::st_sfc(geoms, crs = 3116))
  f_mz <- file.path(dir, "mz.gpkg"); f_rs <- file.path(dir, "rs.gpkg")
  f_us <- file.path(dir, "us.gpkg")
  sf::st_write(mk(list(.sq(0, 0), bowtie, .sq(10, 10)),
                  MANZ_CCNCT = c(a, b, other), DPTO_CCDGO = c("05", "05", "08")),
               f_mz, layer = "mz", quiet = TRUE)
  sf::st_write(mk(list(.sq(3, 0)), SECR_CCNCT = r,
                  DPTO_CCDGO = "05"),
               f_rs, layer = "rs", quiet = TRUE)
  sf::st_write(mk(list(.sq(6, 0), .sq(7, 0)),
                  SECU_CCNCT = paste0(cp, c("001", "002")),
                  DPTO_CCDGO = c("05", "05")),
               f_us, layer = "us", quiet = TRUE)
  pop <- data.table::data.table(
    code = c(a, b, r, cp, other, "05001300099"),
    cls  = c("1", "1", "3", "2", "1", "3"))
  list(pop = pop, F_MZ = f_mz, F_RS = f_rs, F_US = f_us,
       codes = c(a = a, b = b, r = r, cp = cp))
}

test_that("subarea_geometry returns one repaired polygon per known sub-area", {
  skip_if_not_installed("sf")
  s2 <- suppressMessages(sf::sf_use_s2(FALSE))
  withr::defer(suppressMessages(sf::sf_use_s2(s2)))
  m <- .fixture_mgn(withr::local_tempdir())
  # the input holds a self-intersecting polygon, which the output must not
  expect_false(all(sf::st_is_valid(sf::st_read(m$F_MZ, quiet = TRUE))))
  g <- subarea_geometry("05", m$pop, m$F_MZ, m$F_RS, m$F_US)
  expect_s3_class(g, "sf")
  expect_identical(names(g), c("code", "geom"))
  # four sub-areas of department 05; the other department's block is left out,
  # and a rural code with no polygon is simply absent
  expect_setequal(g$code, unname(m$codes))
  expect_identical(anyDuplicated(g$code), 0L)
  expect_true(all(sf::st_is_valid(g)))
  expect_true(all(sf::st_geometry_type(g) %in% c("POLYGON", "MULTIPOLYGON")))
  expect_false(any(sf::st_is_empty(g)))
})

test_that("subarea_geometry dissolves a centro poblado split over several sections", {
  skip_if_not_installed("sf")
  s2 <- suppressMessages(sf::sf_use_s2(FALSE))
  withr::defer(suppressMessages(sf::sf_use_s2(s2)))
  m <- .fixture_mgn(withr::local_tempdir())
  g <- subarea_geometry("05", m$pop, m$F_MZ, m$F_RS, m$F_US)
  cp <- g[g$code == m$codes[["cp"]], ]
  expect_identical(nrow(cp), 1L)
  # two unit squares side by side make one polygon of area 2
  expect_equal(as.numeric(sf::st_area(cp)), 2, tolerance = 1e-9)
})

test_that("subarea_geometry returns NULL for a department without polygons", {
  skip_if_not_installed("sf")
  s2 <- suppressMessages(sf::sf_use_s2(FALSE))
  withr::defer(suppressMessages(sf::sf_use_s2(s2)))
  m <- .fixture_mgn(withr::local_tempdir())
  expect_null(subarea_geometry("11", m$pop, m$F_MZ, m$F_RS, m$F_US))
})

test_that("subarea_union_codes merges the rows of a code into one polygon", {
  skip_if_not_installed("sf")
  g <- sf::st_sf(code = c("a", "b", "a"),
                 geom = sf::st_sfc(list(.sq(0, 0), .sq(5, 0), .sq(2, 0)), crs = 3116))
  u <- subarea_union_codes(g)
  expect_identical(attr(u, "n_duplicate"), 1L)
  expect_identical(sort(u$code), c("a", "b"))
  area <- stats::setNames(as.numeric(sf::st_area(u)), u$code)
  expect_equal(unname(area[c("a", "b")]), c(2, 1))
})

test_that("subarea_union_codes leaves one-row-per-code polygons unchanged", {
  skip_if_not_installed("sf")
  g <- sf::st_sf(code = c("b", "a"),
                 geom = sf::st_sfc(list(.sq(5, 0), .sq(0, 0)), crs = 3116))
  u <- subarea_union_codes(g)
  expect_identical(attr(u, "n_duplicate"), 0L)
  attr(u, "n_duplicate") <- NULL
  expect_identical(u, g)
})
