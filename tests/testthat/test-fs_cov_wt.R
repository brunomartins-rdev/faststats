library(testthat)
library(faststats)

set.seed(123)

test_that("fs_cov_wt matches stats::cov.wt for unweighted unbiased", {
  x <- matrix(rnorm(60), nrow = 20, ncol = 3)

  r_base <- stats::cov.wt(x)
  r_cpp <- fs_cov_wt(x) # use wrapper

  expect_equal(r_cpp$cov, r_base$cov, tolerance = 1e-12)
  expect_equal(as.numeric(r_cpp$center), r_base$center, tolerance = 1e-12)
  expect_equal(r_cpp$n.obs, r_base$n.obs)
})

test_that("fs_cov_wt matches stats::cov.wt for weighted ML with correlation", {
  x <- matrix(rnorm(60), nrow = 20, ncol = 3)
  wt <- runif(nrow(x))

  r_base <- stats::cov.wt(x, wt = wt, method = "ML", cor = TRUE)
  r_cpp <- fs_cov_wt(x, wt = wt, method = "ML", cor = TRUE)

  expect_equal(r_cpp$cov, r_base$cov, tolerance = 1e-12)
  expect_equal(r_cpp$cor, r_base$cor, tolerance = 1e-12)
  expect_equal(as.numeric(r_cpp$center), r_base$center, tolerance = 1e-12)
  expect_equal(r_cpp$n.obs, r_base$n.obs)
  expect_equal(r_cpp$wt, r_base$wt, tolerance = 1e-12)
})

test_that("fs_cov_wt works with user-specified center", {
  x <- matrix(rnorm(60), nrow = 20, ncol = 3)
  center <- colMeans(x)

  r_cpp <- fs_cov_wt(x, center = center)

  # Covariance computed should be same as subtracting provided center manually
  X <- sweep(x, 2, center, "-")
  cov_manual <- crossprod(X) / (nrow(x) - 1)

  expect_equal(r_cpp$cov, cov_manual, tolerance = 1e-12)
})
