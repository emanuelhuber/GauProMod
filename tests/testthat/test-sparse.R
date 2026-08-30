test_that("crossDist_sparse matches dense crossDist when rmax = Inf (regression: index-scrambling bug)", {
  set.seed(1)
  xr <- runif(8, -5, 5)
  yr <- runif(6, -5, 5)
  Dsp <- GauProMod:::crossDist_sparse(matrix(xr, ncol = 1), matrix(yr, ncol = 1), rmax = Inf)
  Dd  <- crossDist(xr, yr, M = NULL, use_symmetry = FALSE)
  expect_equal(as.matrix(Dsp), unname(Dd), tolerance = 1e-10, ignore_attr = TRUE)
})

test_that("crossDist_sparse respects the cutoff and stores the correct entries at the correct positions", {
  set.seed(2)
  xr <- runif(8, -5, 5)
  yr <- runif(6, -5, 5)
  Dd <- crossDist(xr, yr, M = NULL, use_symmetry = FALSE)
  cutoff <- 3
  Dsp <- GauProMod:::crossDist_sparse(matrix(xr, ncol = 1), matrix(yr, ncol = 1), rmax = cutoff)
  Dsp_dense <- as.matrix(Dsp)
  kept <- Dd <= cutoff
  expect_true(all(abs(Dsp_dense[kept] - Dd[kept]) < 1e-10))
  expect_true(all(Dsp_dense[!kept] == 0))
})

test_that("covm(sparse = TRUE) works without an explicit library(Matrix) call", {
  x <- seq(0, 5, length.out = 6)
  Ksparse <- covm(x, x, list(kernel = "matern_3_2", l = 0.5, h = 1),
                   sparse = TRUE, cutoff = 0.3)
  expect_s4_class(Ksparse, "dgCMatrix")
})

test_that("covm(sparse = TRUE) matches the dense result at kept entries and is zero beyond cutoff", {
  set.seed(3)
  n <- 10
  xx <- sort(runif(n, -5, 5))
  covM <- list(kernel = "matern", l = 1.3, h = 1.7, v = 1.5)
  cutoff <- 2
  K_dense  <- covm(xx, xx, covM, use_symmetry = TRUE)
  K_sparse <- covm(xx, xx, covM, sparse = TRUE, cutoff = cutoff, use_symmetry = TRUE)
  K_sparse_dense <- as.matrix(K_sparse)
  D_full <- as.matrix(crossDist(xx, xx, M = NULL, use_symmetry = TRUE))
  kept <- D_full <= cutoff
  expect_true(all(abs(K_sparse_dense[kept] - as.matrix(K_dense)[kept]) < 1e-10))
  expect_true(all(K_sparse_dense[!kept] == 0))
})

test_that("covm(sparse = TRUE) errors clearly for feature/Gram kernels and for d != 0", {
  x <- c(-1, 0, 1)
  expect_error(covm(x, x, list(kernel = "linear", h = 1, c = 0), sparse = TRUE))
  expect_error(covm(x, x, list(kernel = "gaussian", l = 1, h = 1), sparse = TRUE, d = 1))
})
