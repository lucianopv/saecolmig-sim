# ---- area-level matrix --------------------------------------------------------

test_that("area matrix with a missing municipality is refused", {
  X <- fixture_area_matrix(); X <- X[-1, ]
  expect_error(sim_validate_area_matrix(X, fixture_universe()), "1 municipality.*missing")
})

test_that("zero-padded character municipality codes are coerced", {
  X <- fixture_area_matrix(); X$municipality <- sprintf("%05d", X$municipality)
  expect_type(sim_validate_area_matrix(X, fixture_universe())$municipality, "integer")
})

test_that("NA covariates are refused with the column named", {
  X <- fixture_area_matrix(); X[[2]][3] <- NA
  expect_error(sim_validate_area_matrix(X, fixture_universe()), names(X)[2])
})

test_that("extra covariates that collide with existing names are refused", {
  X <- fixture_area_matrix(); extra <- X[, c("municipality", names(X)[2])]
  expect_error(sim_add_covariates(X, extra, "municipality"), names(X)[2])
})

test_that("identifier columns are dropped from the covariates", {
  X <- fixture_area_matrix(); X$n_persons <- 1e6
  v <- suppressMessages(sim_validate_area_matrix(X, fixture_universe()))
  expect_false("n_persons" %in% names(v))
})

test_that("identifier, size and empty columns leave only the key and the covariates", {
  X <- fixture_area_matrix()
  X$dept <- X$municipality %/% 1000L; X$muni <- X$municipality %% 1000L
  X$n_persons <- 1e6; X$n_households <- 3e5
  X$empty <- NA_real_
  msgs <- testthat::capture_messages(
    v <- sim_validate_area_matrix(X, fixture_universe()))
  expect_identical(names(v), names(fixture_area_matrix()))
  expect_match(paste(msgs, collapse = " "), "n_persons")
})

test_that("rows keep the order of the supplied matrix and the result is a data frame", {
  X <- fixture_area_matrix(); X <- X[rev(seq_len(nrow(X))), ]
  v <- sim_validate_area_matrix(tibble::as_tibble(X), fixture_universe())
  expect_s3_class(v, "data.frame"); expect_false(inherits(v, "tbl_df"))
  expect_identical(v$municipality, rev(fixture_universe()))
})

test_that("the refusal names the file the matrix came from", {
  X <- fixture_area_matrix(); X <- X[-(1:2), ]
  expect_error(sim_validate_area_matrix(X, fixture_universe(), source = "area_mat.parquet"),
               "area_mat\\.parquet.*2 municipalities")
})

test_that("duplicated and out-of-universe municipalities are counted", {
  X <- fixture_area_matrix()
  d <- rbind(X, X[1:2, ])
  expect_error(sim_validate_area_matrix(d, fixture_universe()), "2 .*more than once")
  more <- X[1:3, ]; more$municipality <- 9001L:9003L
  e <- rbind(X, more)
  expect_error(sim_validate_area_matrix(e, fixture_universe()), "3 municipalities .*not in the universe")
})

test_that("municipality codes that are not whole numbers are refused", {
  X <- fixture_area_matrix(); X$municipality <- as.numeric(X$municipality); X$municipality[2] <- 5002.5
  expect_error(sim_validate_area_matrix(X, fixture_universe()), "1 .*not .*whole")
  Y <- fixture_area_matrix(); Y$municipality <- as.character(Y$municipality); Y$municipality[4] <- "05x04"
  expect_error(sim_validate_area_matrix(Y, fixture_universe()), "1 .*not .*whole")
  Z <- fixture_area_matrix(); Z$municipality[5] <- NA
  expect_error(sim_validate_area_matrix(Z, fixture_universe()), "1 .*missing .*municipality|municipality.*missing")
})

test_that("whole-number double codes are stored as integer", {
  X <- fixture_area_matrix(); X$municipality <- as.numeric(X$municipality)
  expect_type(sim_validate_area_matrix(X, fixture_universe())$municipality, "integer")
})

test_that("non-numeric and non-finite covariates are refused with the column named", {
  X <- fixture_area_matrix(); X$flag <- "a"
  expect_error(sim_validate_area_matrix(X, fixture_universe()), "not numeric.*flag")
  Y <- fixture_area_matrix(); Y[[3]][2] <- Inf
  expect_error(sim_validate_area_matrix(Y, fixture_universe()), paste0("infinite.*", names(Y)[3]))
})

test_that("a matrix with no covariate left is refused", {
  X <- fixture_area_matrix()[, "municipality", drop = FALSE]; X$n_persons <- 1
  expect_error(suppressMessages(sim_validate_area_matrix(X, fixture_universe())),
               "no covariate")
})

test_that("a matrix without the municipality column is refused", {
  X <- fixture_area_matrix(); names(X)[1] <- "id"
  expect_error(sim_validate_area_matrix(X, fixture_universe()), "municipality")
})

# ---- adding covariates --------------------------------------------------------

test_that("extra covariates land on their own municipality, whatever the row order", {
  X <- fixture_area_matrix()
  extra <- data.frame(municipality = rev(X$municipality), w1 = seq_len(nrow(X)))
  out <- sim_add_covariates(X, extra, "municipality")
  expect_identical(out$municipality, X$municipality)
  expect_identical(out$w1, rev(seq_len(nrow(X))))
  expect_identical(names(out), c(names(X), "w1"))
})

test_that("extra covariates with a zero-padded character key join an integer key", {
  X <- fixture_area_matrix()
  extra <- data.frame(municipality = sprintf("%05d", X$municipality), w1 = 1:12)
  expect_identical(sim_add_covariates(X, extra, "municipality")$w1, 1:12)
})

test_that("extra covariates missing a key of the matrix are refused with the count", {
  X <- fixture_area_matrix()
  extra <- data.frame(municipality = X$municipality[-(1:3)], w1 = 1)
  expect_error(sim_add_covariates(X, extra, "municipality", source = "extra.csv"),
               "extra\\.csv.*3 .*no match")
})

test_that("extra covariates with a repeated key are refused", {
  X <- fixture_area_matrix()
  extra <- data.frame(municipality = c(X$municipality, X$municipality[1]), w1 = 1)
  expect_error(sim_add_covariates(X, extra, "municipality"), "1 .*more than once")
})

test_that("extra rows for keys outside the matrix are ignored", {
  X <- fixture_area_matrix()
  extra <- data.frame(municipality = c(X$municipality, 9999L), w1 = 1)
  out <- suppressMessages(sim_add_covariates(X, extra, "municipality"))
  expect_identical(out$municipality, X$municipality)
})

test_that("extra covariates need the key column and at least one covariate", {
  X <- fixture_area_matrix()
  expect_error(sim_add_covariates(X, data.frame(id = 1, w1 = 1), "municipality"), "municipality")
  expect_error(sim_add_covariates(X, X[, "municipality", drop = FALSE], "municipality"), "no covariate")
})

# ---- sub-area matrix ----------------------------------------------------------

test_that("a complete sub-area matrix passes unchanged", {
  P <- fixture_subarea_matrix()
  v <- sim_validate_subarea_matrix(P)
  expect_s3_class(v, "data.frame")
  expect_identical(v, P)
  expect_type(v$code, "character"); expect_type(v$municipality, "integer")
})

test_that("sub-areas with a missing covariate leave the universe", {
  P <- fixture_subarea_matrix(); P$z1[c(2, 7)] <- NA
  msgs <- testthat::capture_messages(v <- sim_validate_subarea_matrix(P))
  expect_identical(nrow(v), nrow(P) - 2L)
  expect_false(any(P$code[c(2, 7)] %in% v$code))
  expect_match(paste(msgs, collapse = " "), "2 sub-areas")
})

test_that("sub-area keys must be character codes and unique", {
  P <- fixture_subarea_matrix(); P$code <- seq_len(nrow(P))
  expect_error(sim_validate_subarea_matrix(P), "code.*character")
  Q <- fixture_subarea_matrix(); Q$code[2] <- Q$code[1]
  expect_error(sim_validate_subarea_matrix(Q), "1 .*more than once")
})

test_that("the sub-area key may be named subarea", {
  P <- fixture_subarea_matrix(); names(P)[names(P) == "code"] <- "subarea"
  v <- sim_validate_subarea_matrix(P)
  expect_true("code" %in% names(v)); expect_false("subarea" %in% names(v))
})

test_that("sub-area matrices lacking a structural column name it", {
  P <- fixture_subarea_matrix(); P$n_persons_census <- NULL
  expect_error(sim_validate_subarea_matrix(P, source = "blocks.parquet"),
               "blocks\\.parquet.*n_persons_census")
  Q <- fixture_subarea_matrix(); Q$municipality <- NULL
  expect_error(sim_validate_subarea_matrix(Q), "municipality")
})

test_that("sub-area class, counts and covariates are checked", {
  P <- fixture_subarea_matrix(); P$stratum[1:2] <- "periurban"
  expect_error(sim_validate_subarea_matrix(P), "2 .*urban.*rural")
  Q <- fixture_subarea_matrix(); Q$n_households_census[4] <- NA
  expect_error(sim_validate_subarea_matrix(Q), "n_households_census.*1")
  R <- fixture_subarea_matrix(); R$z3 <- "a"
  expect_error(sim_validate_subarea_matrix(R), "not numeric.*z3")
})

test_that("sub-area covariates without any observed value are dropped", {
  P <- fixture_subarea_matrix(); P$empty <- NA_real_
  msgs <- testthat::capture_messages(v <- sim_validate_subarea_matrix(P))
  expect_false("empty" %in% names(v)); expect_identical(nrow(v), nrow(P))
  expect_match(paste(msgs, collapse = " "), "empty")
})

test_that("a sub-area matrix with no covariate is refused", {
  P <- fixture_subarea_matrix()[, c("code", "municipality", "stratum",
                                    "n_households_census", "n_persons_census")]
  expect_error(sim_validate_subarea_matrix(P), "covariate")
})

# ---- the inputs of the runner -------------------------------------------------

test_that("sim_inputs returns the documented elements", {
  inp <- sim_inputs(data_dir = fixture_data_dir())
  expect_identical(names(inp), c("frame", "design", "segments", "pop_totals", "truth",
                                 "X_area", "universe", "W", "pop_sub", "tune",
                                 "area_config", "subarea_config", "sources"))
  expect_identical(inp$universe, fixture_universe())
  expect_identical(inp$area_config, "geo_bld")
  expect_identical(inp$subarea_config, "geo_blk_bld")
  expect_identical(names(inp$tune), sae_outcomes$outcome)
  expect_null(inp$W)
  expect_s3_class(inp$design, "ecv_design_spec")
  expect_s3_class(inp$frame, "data.frame")
  expect_true("stratum" %in% names(inp$frame))
  expect_type(inp$frame$stratum, "character")
  expect_setequal(names(inp$truth), c("area_id", "theta_true", "n_total", "outcome"))
  expect_identical(nrow(inp$truth), length(inp$universe) * nrow(sae_outcomes))
})

test_that("the area matrix is the file's covariates over the universe", {
  inp <- sim_inputs(data_dir = fixture_data_dir())
  expect_identical(inp$X_area, fixture_area_matrix())
  expect_identical(inp$pop_sub, fixture_subarea_matrix())
})

test_that("segments hold every household once under their municipality", {
  inp <- sim_inputs(data_dir = fixture_data_dir())
  expect_identical(names(inp$segments), as.character(fixture_universe()))
  n_seg_rows <- sum(vapply(inp$segments, nrow, integer(1)))
  expect_identical(n_seg_rows, nrow(inp$frame))
  expect_true(all(vapply(inp$segments, function(s) "segment_id" %in% names(s), logical(1))))
})

test_that("pop_totals are the households and persons of each stratum", {
  inp <- sim_inputs(data_dir = fixture_data_dir())
  pt <- inp$pop_totals
  expect_identical(names(pt), c("stratum", "N", "N_persons"))
  expect_identical(pt$stratum, sort(unique(inp$frame$stratum)))
  expect_identical(sum(pt$N), nrow(inp$frame))
  expect_equal(sum(pt$N_persons), sum(inp$frame$n_persons, na.rm = TRUE))
  by_s <- split(inp$frame, inp$frame$stratum)
  expect_equal(pt$N, as.integer(vapply(by_s[pt$stratum], nrow, integer(1))))
})

test_that("a sample drawn on the inputs calibrates to the stratum totals", {
  inp <- sim_inputs(data_dir = fixture_data_dir())
  s <- apply_calibration(
    draw_replicate(inp$frame, inp$design, seed = 100001L, segments = inp$segments),
    inp$pop_totals)
  tot <- tapply(s$weight_cal, s$stratum, sum)
  expect_equal(as.numeric(tot), inp$pop_totals$N[match(names(tot), inp$pop_totals$stratum)],
               tolerance = 1e-6)
})

test_that("segments do not change the sample that is drawn", {
  inp <- sim_inputs(data_dir = fixture_data_dir())
  a <- draw_replicate(inp$frame, inp$design, seed = 100002L, segments = inp$segments)
  b <- draw_replicate(inp$frame, inp$design, seed = 100002L)
  expect_identical(sort(a$hh_id), sort(b$hh_id))
})

test_that("municipalities without a stratum stay in the universe but are not sampled", {
  root <- fixture_data_dir()
  m <- fixture_strata_map()[-1, ]
  write_parquet_safe(m, file.path(root, "aux_tables", "muni_stratum_nn_extended.parquet"))
  msgs <- testthat::capture_messages(inp <- sim_inputs(data_dir = root))
  expect_false(fixture_universe()[1] %in% inp$frame$municipality)
  expect_identical(inp$universe, fixture_universe())
  expect_match(paste(msgs, collapse = " "), "1 municipality")
})

test_that("a supplied area matrix replaces the built one and is checked", {
  root <- fixture_data_dir()
  X <- fixture_area_matrix(); X$extra_col <- seq_len(nrow(X))
  inp <- sim_inputs(data_dir = root, area_matrix = X)
  expect_true("extra_col" %in% names(inp$X_area))
  f <- file.path(withr::local_tempdir(), "mine.parquet")
  write_parquet_safe(X[-1, ], f)
  expect_error(sim_inputs(data_dir = root, area_matrix = f),
               "mine\\.parquet.*1 municipality")
})

test_that("a supplied area matrix does not need the built file", {
  root <- fixture_data_dir()
  unlink(file.path(root, "aux_tables", "municipal_X_geo_bld.parquet"))
  expect_error(sim_inputs(data_dir = root), "municipal_X_geo_bld\\.parquet")
  inp <- sim_inputs(data_dir = root, area_matrix = fixture_area_matrix())
  expect_identical(inp$universe, fixture_universe())
})

test_that("a missing input file is named", {
  root <- fixture_data_dir()
  unlink(file.path(root, "aux_tables", "truth_full.parquet"))
  expect_error(sim_inputs(data_dir = root), "truth_full\\.parquet")
  root2 <- fixture_data_dir()
  expect_error(sim_inputs(data_dir = root2, area_matrix = file.path(root2, "nope.parquet")),
               "nope\\.parquet")
})

test_that("a CSV matrix is read with its codes as text", {
  root <- fixture_data_dir()
  f <- file.path(withr::local_tempdir(), "sub.csv")
  utils::write.csv(fixture_subarea_matrix(), f, row.names = FALSE)
  inp <- sim_inputs(data_dir = root, subarea_matrix = f)
  expect_identical(inp$pop_sub$code, fixture_subarea_matrix()$code)
  g <- file.path(dirname(f), "area.csv")
  X <- fixture_area_matrix(); X$municipality <- sprintf("%05d", X$municipality)
  utils::write.csv(X, g, row.names = FALSE)
  expect_equal(sim_inputs(data_dir = root, area_matrix = g)$X_area, fixture_area_matrix())
})

test_that("extra covariates join the area and the sub-area matrices", {
  root <- fixture_data_dir()
  ea <- data.frame(municipality = fixture_universe(), own_a = seq_len(12))
  sa <- fixture_subarea_matrix()
  es <- data.frame(code = sa$code, own_s = seq_len(nrow(sa)))
  f <- file.path(withr::local_tempdir(), "extra_area.parquet")
  write_parquet_safe(ea, f)
  inp <- sim_inputs(data_dir = root, extra_area = f, extra_subarea = es)
  expect_identical(inp$X_area$own_a, seq_len(12))
  expect_identical(inp$pop_sub$own_s, seq_len(nrow(sa)))
  expect_error(sim_inputs(data_dir = root, extra_area = ea[-1, ]), "1 .*no match")
  clash <- data.frame(code = sa$code, z1 = 1)
  expect_error(sim_inputs(data_dir = root, extra_subarea = clash), "z1")
})

test_that("extra sub-area covariates may be keyed by subarea and cover only the universe", {
  root <- fixture_data_dir()
  P <- fixture_subarea_matrix(); P$z1[1] <- NA
  inc <- P[!is.na(P$z1), ]
  es <- data.frame(subarea = inc$code, own_s = 1)
  suppressMessages(inp <- sim_inputs(data_dir = root, subarea_matrix = P, extra_subarea = es))
  expect_identical(nrow(inp$pop_sub), nrow(P) - 1L)
})

test_that("a sub-area matrix that loses a municipality is refused", {
  root <- fixture_data_dir()
  P <- fixture_subarea_matrix(); P <- P[P$municipality != 5001L, ]
  expect_error(sim_inputs(data_dir = root, subarea_matrix = P), "1 municipality.*no sub-area")
  Q <- fixture_subarea_matrix(); Q$z1[Q$municipality == 8002L] <- NA
  expect_error(suppressMessages(sim_inputs(data_dir = root, subarea_matrix = Q)),
               "1 municipality.*no sub-area")
  R <- fixture_subarea_matrix(); R$municipality[1] <- 9999L
  expect_error(sim_inputs(data_dir = root, subarea_matrix = R), "not in the universe")
})

test_that("truth that lacks an area or an outcome is refused", {
  root <- fixture_data_dir()
  t <- fixture_truth_full(); t <- t[!(t$area_id == 5002L & t$outcome == sae_outcomes$outcome[2]), ]
  write_parquet_safe(t, file.path(root, "aux_tables", "truth_full.parquet"))
  expect_error(sim_inputs(data_dir = root), "truth_full\\.parquet.*1 area")
  u <- fixture_truth_full(); u$n_total <- NULL
  write_parquet_safe(u, file.path(root, "aux_tables", "truth_full.parquet"))
  expect_error(sim_inputs(data_dir = root), "n_total")
})

test_that("truth with a missing value or a repeated area-outcome pair is refused", {
  root <- fixture_data_dir()
  t <- fixture_truth_full(); t$theta_true[4] <- NA
  write_parquet_safe(t, file.path(root, "aux_tables", "truth_full.parquet"))
  expect_error(sim_inputs(data_dir = root),
               "truth_full\\.parquet.*1 value of `theta_true` is missing or not finite")
  t$theta_true[4] <- Inf
  write_parquet_safe(t, file.path(root, "aux_tables", "truth_full.parquet"))
  expect_error(sim_inputs(data_dir = root),
               "truth_full\\.parquet.*1 value of `theta_true` is missing or not finite")
  u <- fixture_truth_full(); u <- rbind(u, u[1:2, ])
  write_parquet_safe(u, file.path(root, "aux_tables", "truth_full.parquet"))
  expect_error(sim_inputs(data_dir = root), "truth_full\\.parquet.*2 area-outcome pairs appear more than once")
})

test_that("a sub-area matrix that matches none of the frame's sub-areas is refused", {
  root <- fixture_data_dir()
  P <- fixture_subarea_matrix(); P$code <- paste0("x", P$code)
  expect_error(sim_inputs(data_dir = root, subarea_matrix = P),
               paste0("36 sub-area codes of the frame \\(holding 100.00% of its households\\) ",
                      "are not in the matrix's prediction universe, more than the 5% of ",
                      "households allowed"))
})

test_that("a sub-area matrix that leaves out many of the frame's households is refused", {
  root <- fixture_data_dir()
  P <- fixture_subarea_matrix()
  first_block <- !duplicated(P$municipality)
  expect_error(sim_inputs(data_dir = root, subarea_matrix = P[!first_block, ]),
               paste0("12 sub-area codes of the frame \\(holding 25.00% of its households\\) ",
                      "are not in the matrix's prediction universe"))
})

test_that("a small gap in the sub-area codes of the frame is reported, not refused", {
  root <- fixture_data_dir()
  P <- fixture_subarea_matrix(); gone <- P$code[1]
  msgs <- testthat::capture_messages(inp <- sim_inputs(data_dir = root, subarea_matrix = P[-1, ]))
  expect_false(gone %in% inp$pop_sub$code)
  expect_match(paste(msgs, collapse = " "),
               "1 sub-area code of the frame \\(holding 2.08% of its households\\) is not in the prediction universe")
})

test_that("queens weights are read when cached and left out otherwise", {
  root <- fixture_data_dir()
  ids <- fixture_universe()
  W <- matrix(0, 12, 12, dimnames = list(NULL, as.character(ids)))
  W[cbind(1:11, 2:12)] <- 1
  write_parquet_safe(cbind(data.frame(municipality = ids), as.data.frame(W)),
                     file.path(root, "aux_tables", "W_queens_municipal.parquet"))
  inp <- sim_inputs(data_dir = root)
  expect_identical(dim(inp$W), c(12L, 12L))
  expect_identical(rownames(inp$W), as.character(ids))
  expect_identical(colnames(inp$W), as.character(ids))
})

test_that("a configuration without shipped hyperparameters leaves tune empty", {
  root <- fixture_data_dir()
  msgs <- testthat::capture_messages(
    inp <- sim_inputs(data_dir = root, subarea_config = "own_cfg",
                      subarea_matrix = fixture_subarea_matrix()))
  expect_null(inp$tune)
  expect_match(paste(msgs, collapse = " "), "own_cfg")
  expect_identical(inp$subarea_config, "own_cfg")
})

test_that("the data directory defaults to SAECOL_DATA_DIR", {
  root <- fixture_data_dir()
  withr::local_envvar(SAECOL_DATA_DIR = root)
  expect_identical(sim_inputs()$universe, fixture_universe())
  withr::local_envvar(SAECOL_DATA_DIR = NA)
  expect_error(sim_inputs(), "SAECOL_DATA_DIR")
})
