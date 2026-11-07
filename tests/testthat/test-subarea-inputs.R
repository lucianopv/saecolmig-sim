# Inputs shared by the sub-area methods (R/method_subarea_common.R).

test_that(".subarea_inputs picks the population weight from the outcome level", {
  f <- .fixture_sub()
  hh <- saecol2026:::.subarea_inputs(f$direct, f$smp, f$pop, "migrated_last_12m")
  pp <- saecol2026:::.subarea_inputs(f$direct, f$smp, f$pop, "pers_migrated_12m")
  expect_identical(hh$wpop, "n_households_census")
  expect_identical(pp$wpop, "n_persons_census")
})

test_that(".subarea_inputs treats every non-identifier column as a covariate", {
  f <- .fixture_sub()
  inp <- saecol2026:::.subarea_inputs(f$direct, f$smp, f$pop, "migrated_last_12m")
  expect_setequal(inp$x_cols, c("z1", "z2", "z3"))
  expect_false(any(c("code", "municipality", "stratum",
                     "n_households_census", "n_persons_census") %in% inp$x_cols))
  expect_identical(deparse(inp$fixed), "y ~ z1 + z2 + z3")
})

test_that("the aggregation weight and its source are never covariates", {
  f <- .fixture_sub()
  f$pop[, `:=`(area_km2 = 2, pop_worldpop = n_persons_census * 1.1)]
  inp <- saecol2026:::.subarea_inputs(f$direct, f$smp, f$pop, "migrated_last_12m")
  expect_setequal(inp$x_cols, c("z1", "z2", "z3"))
})

test_that(".subarea_inputs defines in-sample by the drawn municipalities", {
  f <- .fixture_sub()
  inp <- saecol2026:::.subarea_inputs(f$direct, f$smp, f$pop, "migrated_last_12m")
  expect_identical(inp$in_ids, sort(f$drawn))
  expect_identical(inp$oos_ids, sort(setdiff(seq_len(f$n_muni), f$drawn)))
  expect_length(intersect(inp$in_ids, inp$oos_ids), 0L)
  expect_setequal(c(inp$in_ids, inp$oos_ids), seq_len(f$n_muni))
})

test_that(".subarea_inputs carries the domain id as character for povmap's join", {
  f <- .fixture_sub()
  inp <- saecol2026:::.subarea_inputs(f$direct, f$smp, f$pop, "migrated_last_12m")
  # povmap joins domains to sub-domains and errors on a character/integer mismatch
  expect_type(inp$sm$municipality, "character")
  expect_type(inp$pop$municipality, "character")
  expect_identical(sort(unique(inp$pop$municipality)),
                   sort(as.character(seq_len(f$n_muni))))
})

test_that(".subarea_inputs drops sub-areas that cannot be predicted", {
  f <- .fixture_sub()
  f$pop[1, z1 := NA_real_]
  inp <- saecol2026:::.subarea_inputs(f$direct, f$smp, f$pop, "migrated_last_12m")
  expect_false(f$pop$code[1] %in% inp$pop$code)
  expect_false(f$pop$code[1] %in% inp$sm$code)
})

test_that(".subarea_inputs drops sample rows without a positive finite weight or outcome", {
  f <- .fixture_sub()
  f$smp[1, wt := 0]
  f$smp[2, y := NA_real_]
  inp <- saecol2026:::.subarea_inputs(f$direct, f$smp, f$pop, "migrated_last_12m")
  expect_false(any(f$smp$code[1:2] %in% inp$sm$code))
  expect_identical(nrow(inp$sm), nrow(f$smp) - 2L)
})

test_that(".subarea_inputs errors rather than fitting on nothing", {
  f <- .fixture_sub()
  f$smp[, wt := 0]
  expect_error(
    saecol2026:::.subarea_inputs(f$direct, f$smp, f$pop, "migrated_last_12m"),
    "no usable sub-area sample rows")
})

test_that(".subarea_inputs refuses a population without covariate columns", {
  f <- .fixture_sub()
  bare <- f$pop[, .(code, municipality, stratum,
                    n_households_census, n_persons_census)]
  expect_error(saecol2026:::.subarea_inputs(f$direct, f$smp, bare,
                                            "migrated_last_12m"),
               "no covariate columns")
})

test_that("municipality level rolls the sample and covariates up by hand-computed weights", {
  h <- .fixture_sub_hand()
  inp <- saecol2026:::.subarea_inputs(h$direct, h$smp, h$pop, "migrated_last_12m",
                                      level = "municipality")
  sm <- inp$sm[order(inp$sm$municipality), ]
  pp <- inp$pop[order(inp$pop$municipality), ]
  expect_identical(sm$municipality, c("1", "2"))
  # design-weighted outcome: (4 * 0.2 + 12 * 0.6) / 16 = 0.5; municipality 2 has one row
  expect_equal(sm$y, c(0.5, 0.1))
  expect_equal(sm$wt, c(16, 10))
  # person shares weight by persons: (0.5 * 40 + 0.3 * 20) / 60 and (0.2 * 20 + 0.6 * 60) / 80
  expect_equal(pp$share_working_age, c(26 / 60, 0.5, (0.1 * 10 + 0.4 * 50) / 60))
  # household quantities weight by households: (1 * 10 + 3 * 30) / 40 and (5 * 20 + 7 * 20) / 40
  expect_equal(pp$z1, c(2.5, 6, (2 * 10 + 4 * 30) / 40))
  expect_equal(pp$n_households_census, c(40, 40, 40))
  expect_equal(pp$n_persons_census, c(60, 80, 60))
  # each municipality is its own single cell
  expect_identical(pp$code, c("1", "2", "3"))
  expect_identical(inp$in_ids, c(1L, 2L))
  expect_identical(inp$oos_ids, 3L)
  # the sub-area level keeps one row per sub-area
  sub <- saecol2026:::.subarea_inputs(h$direct, h$smp, h$pop, "migrated_last_12m")
  expect_identical(nrow(sub$sm), 3L)
  expect_identical(nrow(sub$pop), 6L)
})

test_that("the census-free weight needs its column and replaces the census count", {
  f <- .fixture_sub()
  expect_error(
    saecol2026:::.subarea_inputs(f$direct, f$smp, f$pop, "migrated_last_12m",
                                 weight = "worldpop"),
    "pop_worldpop")
  f$pop[, pop_worldpop := n_persons_census * 1.1]
  inp <- saecol2026:::.subarea_inputs(f$direct, f$smp, f$pop, "migrated_last_12m",
                                      weight = "worldpop")
  expect_identical(inp$wpop, "pop_worldpop")
  expect_true("pop_worldpop" %in% colnames(inp$pop))
})
