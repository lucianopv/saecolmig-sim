test_that("add_migrant_flag derives flag + type columns from PA_VIVIA_1ANO", {
  raw <- tibble::tibble(
    hh_id         = c("a", "b", "c", "d", "e", "f"),
    muni          = 100L,
    PA_VIVIA_1ANO = c("1", "2", "3", "4", "9", NA)
  )
  out <- add_migrant_flag(raw)
  # 1 = not yet born (non-migrant); 2 = same muni (non-migrant)
  # 3 = other Colombian muni (internal); 4 = abroad (international)
  # 9 = no informa (NA); NA propagates
  expect_equal(out$migrated_last_12m,
               c(0L, 0L, 1L, 1L, NA_integer_, NA_integer_))
  expect_equal(out$mig_internal,
               c(0L, 0L, 1L, 0L, NA_integer_, NA_integer_))
  expect_equal(out$mig_intl,
               c(0L, 0L, 0L, 1L, NA_integer_, NA_integer_))
})

test_that("add_migrant_flag stops on unexpected codes", {
  raw <- tibble::tibble(PA_VIVIA_1ANO = c("1", "7"))
  expect_error(add_migrant_flag(raw), "Unexpected PA_VIVIA_1ANO")
})

test_that("read_dept_frame returns the household-level schema", {
  skip_if_no_data_dir()
  skip_if_not(file.exists(file.path(census_path(), "05Antioquia",
                                    "CNPV2018_MGN_A2_05.CSV")),
              "Antioquia not unzipped")
  fr <- read_dept_frame("05", n_max = 2000L)
  expected_cols <- c("dept", "muni", "municipality", "area", "seccion",
                     "manzana", "hh_id", "n_persons",
                     "migrated_last_12m", "migrated_internal_12m",
                     "migrated_intl_12m", "n_pers_resp", "n_mig_total",
                     "n_mig_internal", "n_mig_intl")
  expect_true(all(expected_cols %in% colnames(fr)))
  expect_true(all(fr$dept == 5L))
  expect_true(all(fr$area %in% c("urban", "rural")))
  # seccion is the 14-char prefix of manzana
  expect_true(all(nchar(fr$seccion[!is.na(fr$seccion)]) == 14))
  # municipality follows DIVIPOLA: dept(2) * 1000 + muni(3)
  expect_true(all(fr$municipality == fr$dept * 1000L + fr$muni))
})

test_that("aggregate_persons_to_household produces counts and indicators", {
  per <- tibble::tibble(
    U_DPTO = 5L, U_MPIO = 1L, UA_CLASE = 1L,
    COD_ENCUESTAS = c(1L, 1L, 1L, 2L, 2L, 3L),
    U_VIVIENDA    = 1L,
    P_NROHOG      = 1L,
    PA_VIVIA_1ANO = c("2", "3", "4",   # hh1: both types
                      "2", "9",        # hh2: no migrants, one NA
                      "4")             # hh3: intl only
  )
  per <- add_migrant_flag(per)
  hh  <- aggregate_persons_to_household(per)
  hh  <- hh[order(hh$COD_ENCUESTAS), ]
  expect_equal(hh$n_pers_resp,    c(3L, 1L, 1L))
  expect_equal(hh$n_mig_internal, c(1L, 0L, 0L))
  expect_equal(hh$n_mig_intl,     c(1L, 0L, 1L))
  expect_equal(hh$n_mig_total,    c(2L, 0L, 1L))
  expect_equal(hh$migrated_last_12m,     c(1L, 0L, 1L))
  expect_equal(hh$migrated_internal_12m, c(1L, 0L, 0L))
  expect_equal(hh$migrated_intl_12m,     c(1L, 0L, 1L))
  # overlap: hh1 counts in BOTH internal and intl indicators
  expect_true(hh$migrated_internal_12m[1] + hh$migrated_intl_12m[1] >
              hh$migrated_last_12m[1])
})
