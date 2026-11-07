test_that("build_queens_W returns a 1,122 x 1,122 matrix with the right structure", {
  skip_if_no_data_dir()
  skip_if_not_installed("colmaps")
  skip_if_not_installed("spdep")
  skip_if_not_installed("sf")
  W <- build_queens_W(force_recompute = TRUE)
  expect_true(is.matrix(W))
  expect_equal(nrow(W), 1122L)
  expect_equal(ncol(W), 1122L)
  expect_equal(rownames(W), colnames(W))
  expect_equal(diag(W), rep(0, 1122L), ignore_attr = TRUE)
})

test_that("build_queens_W rows sum to 0 (isolated) or 1 (row-standardised)", {
  skip_if_no_data_dir()
  skip_if_not_installed("colmaps")
  W <- build_queens_W(force_recompute = TRUE)
  rs <- rowSums(W)
  expect_true(all(abs(rs - 0) < 1e-8 | abs(rs - 1) < 1e-8))
  n_isolated <- sum(abs(rs - 0) < 1e-8)
  expect_lte(n_isolated, 5L)
})

test_that("build_queens_W indexes match aux$municipality DIVIPOLA convention", {
  skip_if_no_data_dir()
  skip_if_not_installed("colmaps")
  W <- build_queens_W(force_recompute = TRUE)
  ids <- as.integer(rownames(W))
  expect_equal(ids, sort(ids))
  expect_gte(ids[1], 5001L)
  expect_lte(tail(ids, 1L), 99999L)
})

# build_dept_adjacency(): every municipality is a neighbour, with equal
# weight, of every other municipality of its department.
test_that("build_dept_adjacency preserves dept-cluster invariants", {
  skip_if_not_installed("Matrix")
  ids <- c(5001L, 5002L, 5003L, 11001L, 13001L, 13002L)
  W <- build_dept_adjacency(ids)
  expect_true(inherits(W, "Matrix") || inherits(W, "dgCMatrix"))
  expect_equal(dim(W), c(length(ids), length(ids)))
  expect_equal(rownames(W), as.character(ids))

  row_5001 <- W[as.character(5001L), ]
  expect_equal(sum(row_5001), 1)
  expect_equal(unname(row_5001[as.character(5002L)]), 0.5)
  expect_equal(unname(row_5001[as.character(5003L)]), 0.5)

  # Singleton dept (11001) gets a zero row.
  expect_equal(sum(W[as.character(11001L), ]), 0)
})
