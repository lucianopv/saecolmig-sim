# Census covariates at sub-area resolution (R/subarea_census.R), on a small
# synthetic CNPV layout whose shares can be counted by hand.

test_that("build_dept_census_subarea reproduces the hand-counted shares per sub-area", {
  root <- withr::local_tempdir()
  codes <- .write_cnpv_fixture(root)
  withr::local_envvar(SAECOL_DATA_DIR = root)
  out <- build_dept_census_subarea("05")
  expect_s3_class(out, "data.table")
  expect_identical(nrow(out), 3L)
  expect_setequal(out$code, c(codes$u1, codes$u2, codes$r))
  expect_true(all(nchar(out$code[out$code != codes$r]) == 22L))
  expect_identical(nchar(codes$r), 11L)
  expect_true(all(out$dept == 5L))
  g <- function(cd, col) out[[col]][out$code == cd]

  # block U1: four persons, two households, two dwellings
  expect_identical(g(codes$u1, "n_persons_sub"), 4L)
  expect_identical(g(codes$u1, "n_households_sub"), 2L)
  expect_equal(g(codes$u1, "share_working_age"), 2 / 4)      # ages 4 and 6
  expect_equal(g(codes$u1, "share_elderly"), 1 / 4)          # age 14
  expect_equal(g(codes$u1, "share_secondary_completed"), 2 / 4)
  expect_equal(g(codes$u1, "share_illiterate"), 2 / 4)       # ages 4 and 6, alfa 2
  expect_equal(g(codes$u1, "share_economically_active"), 2 / 4)
  expect_equal(g(codes$u1, "share_afro"), 1 / 4)
  expect_equal(g(codes$u1, "share_indigenous"), 1 / 4)
  expect_equal(g(codes$u1, "share_piped_water"), 1 / 2)
  expect_equal(g(codes$u1, "share_electricity"), 1)

  # block U2: two persons, one household; the work variable is entirely missing
  expect_identical(g(codes$u2, "n_persons_sub"), 2L)
  expect_equal(g(codes$u2, "share_working_age"), 1 / 2)
  expect_equal(g(codes$u2, "share_elderly"), 1 / 2)
  expect_equal(g(codes$u2, "share_secondary_completed"), 1 / 2)
  expect_equal(g(codes$u2, "share_illiterate"), 0)
  expect_equal(g(codes$u2, "share_piped_water"), 1)
  expect_equal(g(codes$u2, "share_electricity"), 0)

  # rural section: two dwellings with different block codes share the section code
  expect_identical(g(codes$r, "n_persons_sub"), 3L)
  expect_identical(g(codes$r, "n_households_sub"), 2L)
  expect_equal(g(codes$r, "share_working_age"), 2 / 3)
  expect_equal(g(codes$r, "share_elderly"), 0)
  expect_equal(g(codes$r, "share_secondary_completed"), 1 / 3)
  expect_equal(g(codes$r, "share_illiterate"), 1 / 3)
  expect_equal(g(codes$r, "share_indigenous"), 2 / 3)
  expect_equal(g(codes$r, "share_piped_water"), 1 / 2)
  expect_equal(g(codes$r, "share_electricity"), 1)
})

test_that("a share over no non-missing value is NA, not NaN", {
  root <- withr::local_tempdir()
  codes <- .write_cnpv_fixture(root)
  withr::local_envvar(SAECOL_DATA_DIR = root)
  out <- build_dept_census_subarea("05")
  x <- out$share_economically_active[out$code == codes$u2]
  expect_true(is.na(x))
  expect_false(is.nan(x))
  expect_false(any(vapply(out, function(v) is.numeric(v) && any(is.nan(v)), TRUE)))
})

test_that("n_max reads only the first person records", {
  root <- withr::local_tempdir()
  .write_cnpv_fixture(root)
  withr::local_envvar(SAECOL_DATA_DIR = root)
  out <- build_dept_census_subarea("05", n_max = 4L)
  expect_identical(sum(out$n_persons_sub, na.rm = TRUE), 4L)
})

test_that("build_dept_census_subarea names a missing department", {
  root <- withr::local_tempdir()
  .write_cnpv_fixture(root)
  withr::local_envvar(SAECOL_DATA_DIR = root)
  expect_error(build_dept_census_subarea("99"), "Department folder not found")
})

test_that("build_census_subarea stacks the requested departments", {
  root <- withr::local_tempdir()
  .write_cnpv_fixture(root)
  withr::local_envvar(SAECOL_DATA_DIR = root)
  all_dept <- build_census_subarea(verbose = FALSE)
  one <- build_census_subarea("05", verbose = FALSE)
  expect_identical(nrow(all_dept), 3L)
  expect_identical(all_dept, one)
})
