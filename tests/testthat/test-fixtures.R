test_that("make_synthetic_frame returns the expected shape", {
  fr <- make_synthetic_frame()
  expect_equal(nrow(fr), 240)
  expect_setequal(unique(fr$stratum), c("99_urban", "99_rural"))
  # 4 municipalities per stratum -> 8 unique
  expect_equal(length(unique(fr$municipality)), 8)
  # 6 manzanas per municipality -> 8 * 6 = 48
  expect_equal(length(unique(fr$manzana)), 48)
  expect_equal(length(unique(fr$hh_id)), 240)
  # binary indicators are 0/1
  expect_true(all(fr$migrated_last_12m     %in% c(0L, 1L)))
  expect_true(all(fr$migrated_internal_12m %in% c(0L, 1L)))
  expect_true(all(fr$migrated_intl_12m     %in% c(0L, 1L)))
  # counts are non-negative integers
  expect_true(all(fr$n_persons     >= 1L))
  expect_true(all(fr$n_pers_resp   >= 1L))
  expect_true(all(fr$n_mig_total   >= 0L))
  expect_true(all(fr$n_mig_internal >= 0L))
  expect_true(all(fr$n_mig_intl    >= 0L))
  # consistency: total = internal + intl
  expect_equal(fr$n_mig_total, fr$n_mig_internal + fr$n_mig_intl)
  # indicator consistency: migrated_last_12m == 1 iff n_mig_total > 0
  expect_equal(fr$migrated_last_12m, as.integer(fr$n_mig_total > 0L))
  expect_equal(fr$migrated_internal_12m, as.integer(fr$n_mig_internal > 0L))
  expect_equal(fr$migrated_intl_12m,     as.integer(fr$n_mig_intl > 0L))
  # n_pers_resp == n_persons (full response in synthetic fixture)
  expect_equal(fr$n_pers_resp, fr$n_persons)
})
