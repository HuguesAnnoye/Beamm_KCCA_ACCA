test_that("compatibility matrix works 1", {
  df.don <- data.frame(x1 = factor(c("a", "b", "c"), levels = c("a", "b", "c")),
                       x2 = factor(c("c", "b", "a"), levels = c("a", "b", "c")))
  df.rec <- data.frame(x1 = factor(c("a", "b", "c"), levels = c("a", "b", "c")),
                       x2 = factor(c("c", "b", "a"), levels = c("a", "b", "c")))
  comp.mtx <- compute_mtx_compatibilities(df.rec, df.don, c("x1", "x2"), index = FALSE)
  comp.idx <- compute_mtx_compatibilities(df.rec, df.don, c("x1", "x2"), index = TRUE)
  expect_equal(comp.mtx, diag(rep(1,3)))
  expect_equal(comp.idx, cbind(1:3,1:3), check.attributes = FALSE)
})

test_that("compatibility matrix works 2", {
  df.don <- data.frame(x1 = factor(c("a"), levels = c("a", "b", "c")),
                       x2 = factor(c("c"), levels = c("a", "b", "c")))
  df.rec <- data.frame(x1 = factor(c("a", "b", "c"), levels = c("a", "b", "c")),
                       x2 = factor(c("a", "b", "a"), levels = c("a", "b", "c")))
  comp.mtx <- compute_mtx_compatibilities(df.rec, df.don, c("x1", "x2"), index = FALSE)
  comp.idx <- compute_mtx_compatibilities(df.rec, df.don, c("x1", "x2"), index = TRUE)
  expect_equal(comp.mtx, rbind(rep(1,3)))
  expect_equal(comp.idx, cbind(rep(1,3),1:3), check.attributes = FALSE)
})

test_that("compatibility matrix works 3", {
  df.don <- data.frame(x1 = factor(c("a", "b", "c"), levels = c("a", "b", "c")),
                       x2 = factor(c("c", "b", "a"), levels = c("a", "b", "c")))
  df.rec <- data.frame(x1 = factor(c("a"), levels = c("a", "b", "c")),
                       x2 = factor(c("c"), levels = c("a", "b", "c")))
  comp.mtx <- compute_mtx_compatibilities(df.rec, df.don, c("x1", "x2"), index = FALSE)
  comp.idx <- compute_mtx_compatibilities(df.rec, df.don, c("x1", "x2"), index = TRUE)
  expect_equal(comp.mtx, cbind(c(1,0,0)))
  expect_equal(comp.idx, rbind(c(1,1)), check.attributes = FALSE)
})

test_that("ompute_mtx_n_diff_CV_rcpp works 1", {
  df.don <- data.frame(x1 = factor(c("a"), levels = c("a", "b", "c")),
                       x2 = factor(c("c"), levels = c("a", "b", "c")))
  df.rec <- data.frame(x1 = factor(c("a", "b", "c"), levels = c("a", "b", "c")),
                       x2 = factor(c("a", "b", "a"), levels = c("a", "b", "c")))
  mtx.don <- df2mtx(df.don)
  mtx.rec <- df2mtx(df.rec)
  n_diff_CV <- compute_mtx_n_diff_CV_rcpp(mtx.don, mtx.rec)/2
  expect_equal(n_diff_CV, rbind(c(1,2,2)))
})

test_that("ompute_mtx_n_diff_CV_rcpp works 2", {
  df.don <- data.frame(x1 = factor(c("a", "b", "c"), levels = c("a", "b", "c")),
                       x2 = factor(c("c", "b", "a"), levels = c("a", "b", "c")))
  df.rec <- data.frame(x1 = factor(c("a"), levels = c("a", "b", "c")),
                       x2 = factor(c("c"), levels = c("a", "b", "c")))
  mtx.don <- df2mtx(df.don)
  mtx.rec <- df2mtx(df.rec)
  n_diff_CV <- compute_mtx_n_diff_CV_rcpp(mtx.don, mtx.rec)/2
  expect_equal(n_diff_CV, cbind(c(0,2,2)))
})

test_that("ompute_mtx_n_diff_CV_rcpp works 2", {
  df.don <- data.frame(x1 = factor(c("a", "b", "c"), levels = c("a", "b", "c")),
                       x2 = factor(c("c", "b", "a"), levels = c("a", "b", "c")))
  df.rec <- data.frame(x1 = factor(c("a", "b", "c"), levels = c("a", "b", "c")),
                       x2 = factor(c("c", "b", "a"), levels = c("a", "b", "c")))
  mtx.don <- df2mtx(df.don)
  mtx.rec <- df2mtx(df.rec)
  n_diff_CV <- compute_mtx_n_diff_CV_rcpp(mtx.don, mtx.rec)/2
  test <- matrix(rep(2,9), 3, 3)
  diag(test) <- 0
  expect_equal(n_diff_CV, test, check.attributes = FALSE)
})
