test_that("matrix scaling and its inverse fully numeric", {
  nrow <- 1000
  ncol <- 10
  mtx <- matrix(runif(nrow*ncol, -5, 5), nrow, ncol)
  wt <- runif(nrow, 0, 5000)
  mtx_scaled <- scale_mtx(mtx, wt)
  mtx_inv_scaled <- scale_mtx_inv(mtx_scaled,
                                  center = attr(mtx_scaled, "scaled:center"),
                                  scale = attr(mtx_scaled, "scaled:scale"))
  expect_equal(mtx, mtx_inv_scaled, check.attributes = FALSE)
})

test_that("matrix scaling and its inverse only dummies", {
  nrow <- 1000
  ncol <- 10
  mtx <- matrix(sample(c(0,1), nrow*ncol, replace = T), nrow, ncol)
  wt <- runif(nrow, 0, 5000)
  mtx_scaled <- scale_mtx(mtx, wt)
  mtx_inv_scaled <- scale_mtx_inv(mtx_scaled,
                                  center = attr(mtx_scaled, "scaled:center"),
                                  scale = attr(mtx_scaled, "scaled:scale"))
  expect_equal(mtx, mtx_inv_scaled, check.attributes = FALSE)
  expect_true(sum(abs(mtx - mtx_inv_scaled) > 0) == 0)
})

test_that("matrix scaling and its inverse mixed numeric and dummies", {
  nrow <- 1000
  ncol <- 10
  mtx <- matrix(runif(nrow*ncol, -5, 5), nrow, ncol)
  mtx <- cbind(mtx, matrix(sample(c(0,1), nrow*ncol, replace = T), nrow, ncol))
  wt <- runif(nrow, 0, 5000)
  mtx_scaled <- scale_mtx(mtx, wt)
  mtx_inv_scaled <- scale_mtx_inv(mtx_scaled,
                                  center = attr(mtx_scaled, "scaled:center"),
                                  scale = attr(mtx_scaled, "scaled:scale"))
  expect_equal(mtx, mtx_inv_scaled, check.attributes = FALSE)
  expect_true(sum(abs(mtx[,11:20] - mtx_inv_scaled[,11:20]) > 0) == 0)
})
