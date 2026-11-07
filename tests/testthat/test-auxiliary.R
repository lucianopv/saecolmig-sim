test_that("build_dept_aux returns one row per municipality with core columns", {
  skip_if_no_data_dir()
  skip_if_not(file.exists(file.path(census_path(), "05Antioquia",
                                    "CNPV2018_5PER_A2_05.CSV")),
              "Antioquia not unzipped")
  aux <- build_dept_aux("05", n_max = 5000L)
  required_cols <- c("dept", "muni", "municipality",
                     "n_persons", "n_households",
                     "share_urban", "share_working_age", "share_elderly",
                     "share_secondary_completed", "share_illiterate",
                     "share_economically_active",
                     "share_piped_water", "share_electricity",
                     "share_afro", "share_indigenous")
  expect_true(all(required_cols %in% colnames(aux)),
              info = paste("missing:",
                           paste(setdiff(required_cols, colnames(aux)),
                                 collapse = ", ")))
  expect_equal(nrow(aux), length(unique(aux$municipality)))
  share_cols <- grep("^share_", colnames(aux), value = TRUE)
  for (col in share_cols) {
    v <- aux[[col]]
    expect_true(all(v >= 0 & v <= 1, na.rm = TRUE), info = col)
  }
})

test_that("bind_aux_layer left-joins on municipality and preserves rows", {
  base  <- tibble::tibble(municipality = c(5001L, 5002L, 5003L),
                          share_urban  = c(0.9, 0.5, 0.1))
  layer <- tibble::tibble(municipality = c(5001L, 5003L, 5099L),
                          ndvi         = c(0.4, 0.7, 0.2))
  out <- bind_aux_layer(base, layer)
  expect_equal(nrow(out), 3)
  expect_true("ndvi" %in% colnames(out))
  expect_equal(out$ndvi[out$municipality == 5002L], NA_real_)
  expect_equal(out$ndvi[out$municipality == 5001L], 0.4)
})

test_that("bind_aux_layer errors when layer lacks municipality column", {
  base  <- tibble::tibble(municipality = 5001L, share_urban = 0.9)
  layer <- tibble::tibble(muni_code = 5001L, ndvi = 0.4)
  expect_error(bind_aux_layer(base, layer),
               "must contain a 'municipality' column")
})

test_that("the municipal covariate builder takes the department codes and nothing else", {
  expect_identical(names(formals(build_aux_municipal)), "dept_codes")
})
