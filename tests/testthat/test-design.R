test_that("parse_ecv_design loads the shipped YAML and returns an ecv_design_spec object", {
  design <- parse_ecv_design()
  expect_s3_class(design, "ecv_design_spec")
  expect_true(all(c("stage1", "stage2", "stage3", "strata", "calibration")
                  %in% names(design)))
  # Per ECV 2024 DDI: PSU = municipality, SSU = ~10-dwelling segment.
  expect_equal(design$stage1$unit, "municipality")
  expect_equal(design$stage2$unit, "segment")
  expect_true(is.data.frame(design$strata))
  expect_true(all(c("stratum", "n_psu") %in% colnames(design$strata)))
})

test_that("parse_ecv_design errors on missing path", {
  expect_error(parse_ecv_design("/nonexistent/path.yml"),
               "Design YAML not found")
})
