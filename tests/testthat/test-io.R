test_that("output_path errors clearly when SAECOL_DATA_DIR is unset", {
  old <- Sys.getenv("SAECOL_DATA_DIR", unset = NA)
  withr::defer({
    if (is.na(old)) Sys.unsetenv("SAECOL_DATA_DIR") else Sys.setenv(SAECOL_DATA_DIR = old)
  })
  Sys.unsetenv("SAECOL_DATA_DIR")
  expect_error(output_path("x"), "SAECOL_DATA_DIR")
})

test_that("output_path resolves nested folders under SAECOL_DATA_DIR and creates them on demand", {
  old <- Sys.getenv("SAECOL_DATA_DIR", unset = NA)
  withr::defer({
    if (is.na(old)) Sys.unsetenv("SAECOL_DATA_DIR") else Sys.setenv(SAECOL_DATA_DIR = old)
  })
  tmp_root <- withr::local_tempdir()
  Sys.setenv(SAECOL_DATA_DIR = tmp_root)

  p <- output_path("samples", create = TRUE)
  expect_true(grepl("samples", p, fixed = TRUE))
  expect_true(startsWith(p, tmp_root))
  expect_true(dir.exists(p))
})

test_that("census_path and ecv_path resolve under SAECOL_DATA_DIR", {
  old <- Sys.getenv("SAECOL_DATA_DIR", unset = NA)
  withr::defer({
    if (is.na(old)) Sys.unsetenv("SAECOL_DATA_DIR") else Sys.setenv(SAECOL_DATA_DIR = old)
  })
  tmp_root <- withr::local_tempdir()
  Sys.setenv(SAECOL_DATA_DIR = tmp_root)

  expect_true(grepl("Census2018", census_path(), fixed = TRUE))
  expect_true(grepl("ECV2024", ecv_path(), fixed = TRUE))
})

test_that("parquet round-trip preserves a small tibble", {
  tmp <- withr::local_tempfile(fileext = ".parquet")
  df  <- tibble::tibble(a = 1:3, b = letters[1:3])
  write_parquet_safe(df, tmp)
  out <- read_parquet_safe(tmp)
  expect_equal(as.data.frame(out), as.data.frame(df))
})

test_that("no path component output_path() produces is a Windows reserved device name", {
  old <- Sys.getenv("SAECOL_DATA_DIR", unset = NA)
  withr::defer({
    if (is.na(old)) Sys.unsetenv("SAECOL_DATA_DIR") else Sys.setenv(SAECOL_DATA_DIR = old)
  })
  tmp_root <- withr::local_tempdir()
  Sys.setenv(SAECOL_DATA_DIR = tmp_root)

  # Win32 refuses these as file OR directory names, on any filesystem, and a
  # zip holding such an entry will not extract there either. The check is on
  # the component, case-insensitively, with any extension stripped.
  reserved <- c("con", "prn", "aux", "nul",
                paste0("com", 1:9), paste0("lpt", 1:9))
  is_reserved <- function(path) {
    parts <- strsplit(path, "[/\\\\]")[[1]]
    parts <- sub("\\..*$", "", parts)
    any(tolower(parts) %in% reserved)
  }

  # Every subfolder this package writes to or reads from.
  folders <- c("aux_tables", "frame", "direct", "samples", "diagnostics",
               "estimates", "figures", "reports", "geo_cache", "geo",
               "Census2018", "ECV2024", "ecv2024")
  for (f in folders) {
    expect_false(is_reserved(f), label = paste0("folder name '", f, "'"))
  }
  # And the resolved paths themselves, so a reserved name in any of them fails
  # here.
  expect_false(is_reserved(output_path("aux_tables")))
  expect_false(is_reserved(census_path()))
  expect_false(is_reserved(ecv_path()))
})

test_that("output_path maps an \"aux\" component onto aux_tables on every platform", {
  old <- Sys.getenv("SAECOL_DATA_DIR", unset = NA)
  withr::defer({
    if (is.na(old)) Sys.unsetenv("SAECOL_DATA_DIR") else Sys.setenv(SAECOL_DATA_DIR = old)
  })
  tmp_root <- withr::local_tempdir()
  Sys.setenv(SAECOL_DATA_DIR = tmp_root)

  # The mapping is platform-independent, so a path names the same directory
  # wherever it is built.
  expect_identical(output_path("aux"), output_path("aux_tables"))
  expect_false(grepl("(^|/)aux$", output_path("aux")))
  # Nested components are remapped too, not just the first.
  expect_identical(output_path("aux", "truth.parquet"),
                   output_path("aux_tables", "truth.parquet"))
})
