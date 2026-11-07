# tipsae compiles its Stan models when it is installed from source, and its
# configure step loads rstantools (>= 2.1.1), which tipsae lists only under
# Suggests. A lockfile built from the Imports, Depends and LinkingTo fields
# leaves rstantools out, and a restore that builds tipsae from source then stops
# at tipsae; the lockfile therefore records rstantools explicitly.
# renv does not follow Suggests, so no recorded requirement puts rstantools
# before tipsae: a restore installs it first only through renv's install order
# (renv 1.1.4, from an empty cache).
test_that("the lockfile records rstantools beside tipsae", {
  lock <- test_path("..", "..", "renv.lock")
  skip_if_not(file.exists(lock), "renv.lock is not in the source tree")
  pk <- jsonlite::fromJSON(lock, simplifyVector = FALSE)$Packages
  expect_true("tipsae" %in% names(pk))
  expect_true("rstantools" %in% names(pk))
  v <- if (is.null(pk$rstantools$Version)) "0" else pk$rstantools$Version
  expect_true(utils::compareVersion(v, "2.1.1") >= 0)
})
