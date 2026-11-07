.fixture_subarea <- function() {
  # Two municipalities. Muni 1 has 2 urban blocks + 1 rural section, muni 2 has
  # 1 urban block + 1 rural section. The `seccion` key is 14 chars, block codes 22,
  # both carrying the municipality in their prefix, as in the real frame.
  sec <- function(m, s) sprintf("%05d1%04d%04d", m, s, s)      # 14 chars
  blk <- function(m, s, b) paste0(sec(m, s), sprintf("%08d", b))  # 22 chars
  rows <- rbind(
    data.frame(municipality = 1L, area = "urban", seccion = sec(1, 1),
               manzana = blk(1, 1, 1), n_persons = c(3, 4)),
    data.frame(municipality = 1L, area = "urban", seccion = sec(1, 2),
               manzana = blk(1, 2, 2), n_persons = 5),
    data.frame(municipality = 1L, area = "rural", seccion = sec(1, 3),
               manzana = blk(1, 3, 3), n_persons = c(2, 2, 6)),
    data.frame(municipality = 2L, area = "urban", seccion = sec(2, 1),
               manzana = blk(2, 1, 1), n_persons = 7),
    data.frame(municipality = 2L, area = "rural", seccion = sec(2, 2),
               manzana = blk(2, 2, 2), n_persons = c(1, 9))
  )
  data.table::as.data.table(rows)
}

test_that("subarea_code uses the block when urban and the truncated section when rural", {
  a <- c("urban", "rural")
  m <- c(strrep("9", 22), strrep("8", 22))
  s <- c(strrep("7", 14), "12345678901234")
  out <- subarea_code(a, m, s)
  expect_identical(out[1], strrep("9", 22))          # urban -> whole manzana
  expect_identical(out[2], "12345678901")            # rural -> first 11 chars
  expect_identical(nchar(out), c(22L, 11L))
})

test_that("subarea_code returns NA for an unrecognised area class", {
  expect_true(is.na(subarea_code("suburban", strrep("9", 22), strrep("7", 14))))
})

test_that("urban and rural codes cannot collide: they differ in width", {
  f <- .fixture_subarea()
  cd <- subarea_code(f$area, f$manzana, f$seccion)
  expect_length(intersect(cd[f$area == "urban"], cd[f$area == "rural"]), 0L)
})

test_that("build_subarea_pop aggregates households and persons per sub-area", {
  pop <- build_subarea_pop(.fixture_subarea())
  expect_equal(nrow(pop), 5L)                        # 3 urban blocks + 2 rural sections
  expect_setequal(colnames(pop),
                  c("code", "municipality", "stratum",
                    "n_households_census", "n_persons_census"))
  expect_equal(sum(pop$n_households_census), nrow(.fixture_subarea()))
  expect_equal(sum(pop$n_persons_census), sum(.fixture_subarea()$n_persons))
  expect_setequal(pop$stratum, c("urban", "rural"))
  # every sub-area sits in exactly one municipality
  expect_equal(nrow(unique(pop[, c("code", "municipality")])), nrow(pop))
})

test_that("build_subarea_pop errors on a missing column rather than part-building", {
  f <- .fixture_subarea()
  f[, manzana := NULL]
  expect_error(build_subarea_pop(f), "missing")
})

test_that("build_subarea_pop rejects a code spanning two municipalities", {
  f <- .fixture_subarea()
  # Force a collision: give a muni-2 row a muni-1 block code.
  f[municipality == 2L & area == "urban", manzana := f[municipality == 1L][1]$manzana]
  expect_error(build_subarea_pop(f), "more than one municipality")
})

test_that("collapse_sample_to_subarea weights a household outcome by weight_cal", {
  pop <- build_subarea_pop(.fixture_subarea())
  s   <- .fixture_subarea()
  s[, `:=`(weight_cal = c(2, 2, 1, 3, 3, 3, 4, 5, 5),
           migrated_last_12m = c(1, 0, 1, 0, 1, 1, 0, 1, 0),
           n_pers_resp = n_persons, n_mig_total = 0L)]
  sa <- collapse_sample_to_subarea(s, pop, "migrated_last_12m")

  expect_true(all(c("code", "municipality", "y", "wt", "n_hh") %in% colnames(sa)))
  expect_equal(sum(sa$n_hh), nrow(s))
  expect_true(all(sa$y >= 0 & sa$y <= 1))
  # first urban block: weights 2,2 and outcomes 1,0 -> 2/4
  first_blk <- s[1]$manzana
  expect_equal(sa[code == first_blk]$y, 0.5)
  expect_equal(sa[code == first_blk]$wt, 4)
})

test_that("collapse_sample_to_subarea uses the count ratio for person outcomes", {
  pop <- build_subarea_pop(.fixture_subarea())
  s   <- .fixture_subarea()
  s[, `:=`(weight_cal = 2, n_pers_resp = n_persons,
           n_mig_intl = ifelse(n_persons > 4, 1L, 0L))]
  sa <- collapse_sample_to_subarea(s, pop, "pers_migrated_intl_12m")
  # sum(w*num)/sum(w*den) on the first block: persons 3,4 -> no migrants
  expect_equal(sa[code == s[1]$manzana]$y, 0)
  # the denominator is weighted persons, not weighted households
  expect_equal(sa[code == s[1]$manzana]$wt, 2 * (3 + 4))
})

test_that("collapse_sample_to_subarea drops unmatched rows and records how many", {
  pop <- build_subarea_pop(.fixture_subarea())
  s   <- .fixture_subarea()
  s[, `:=`(weight_cal = 1, migrated_last_12m = 1)]
  # A block and section belonging to no known sub-area.
  extra <- s[1]
  extra[, `:=`(manzana = strrep("0", 22), seccion = strrep("0", 14))]
  s2 <- rbind(s, extra)
  sa <- collapse_sample_to_subarea(s2, pop, "migrated_last_12m")
  expect_equal(attr(sa, "n_unmatched"), 1L)
  expect_equal(sum(sa$n_hh), nrow(s))
})

test_that("collapse_sample_to_subarea leaves the caller's sample untouched", {
  pop <- build_subarea_pop(.fixture_subarea())
  s   <- .fixture_subarea()
  s[, `:=`(weight_cal = 1, migrated_last_12m = 1)]
  before <- data.table::copy(s)
  collapse_sample_to_subarea(s, pop, "migrated_last_12m")
  expect_identical(s, before)
})

test_that("collapse_sample_to_subarea joins a rural row through its truncated section", {
  pop <- build_subarea_pop(.fixture_subarea())
  s   <- .fixture_subarea()
  s[, `:=`(weight_cal = 1, migrated_last_12m = c(0, 0, 0, 1, 0, 0, 0, 0, 1))]
  sa <- collapse_sample_to_subarea(s, pop, "migrated_last_12m")
  rural_code <- substr(s[municipality == 1L & area == "rural"]$seccion[1], 1L, 11L)
  expect_identical(sa[code == rural_code]$n_hh, 3L)
  expect_equal(sa[code == rural_code]$y, 1 / 3)
  expect_identical(sa[code == rural_code]$municipality, 1L)
})
