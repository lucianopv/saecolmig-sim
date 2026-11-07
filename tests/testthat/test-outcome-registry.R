test_that("sae_outcomes registers 6 outcomes with consistent columns", {
  expect_equal(nrow(sae_outcomes), 6L)
  expect_setequal(sae_outcomes$level, c("household", "person"))
  hh <- sae_outcomes[sae_outcomes$level == "household", ]
  pp <- sae_outcomes[sae_outcomes$level == "person", ]
  expect_true(all(!is.na(hh$indicator)))
  expect_true(all(!is.na(pp$numerator)) && all(!is.na(pp$denominator)))
  expect_true(all(pp$denominator == "n_pers_resp"))
})

test_that("outcome_spec returns one row and errors on unknown", {
  o <- outcome_spec("pers_migrated_intl_12m")
  expect_equal(o$level, "person")
  expect_equal(o$numerator, "n_mig_intl")
  expect_error(outcome_spec("nope"), "Unknown outcome")
})
