test_that("linear kernel: square case matches h^2*x*y + c^2 (regression: W-dimension bug)", {
  x <- c(-4, -3, -1, 0, 4)
  covModel <- list(kernel = "linear", h = 1.5, c = 0)
  K <- covm(x, x, covModel, use_symmetry = TRUE)
  manual <- outer(x, x, function(a, b) 1.5^2 * (a * b) + 0^2)
  expect_equal(unname(as.matrix(K)), manual, tolerance = 1e-10)
})

test_that("linear kernel: non-square case works (regression: W-dimension bug)", {
  x <- c(-4, -3, -1, 0, 4)
  z <- seq(-10, 10, length.out = 50)
  covModel <- list(kernel = "linear", h = 1.5, c = 0.5)
  K <- covm(x, z, covModel, use_symmetry = FALSE)
  expect_equal(dim(K), c(length(x), length(z)))
  manual <- outer(x, z, function(a, b) 1.5^2 * ((a - 0.5) * (b - 0.5)) + 0.5^2)
  expect_equal(unname(as.matrix(K)), manual, tolerance = 1e-10)
})

test_that("linear kernel ignores 'b' (parameter was removed as unused/dead)", {
  x <- c(-2, -1, 0, 1, 2)
  covModel_with_b    <- list(kernel = "linear", h = 1.2, c = 0.3, b = 999)
  covModel_without_b <- list(kernel = "linear", h = 1.2, c = 0.3)
  K1 <- covm(x, x, covModel_with_b, use_symmetry = TRUE)
  K2 <- covm(x, x, covModel_without_b, use_symmetry = TRUE)
  expect_equal(as.matrix(K1), as.matrix(K2))
})

test_that("polynomial kernel matches h^2*(x*y+c)^degree (regression: missing Y arg + wrong dispatch)", {
  x <- c(-4, -3, -1, 0, 4)
  for (degree in c(2, 3)) {
    covModel <- list(kernel = "polynomial", degree = degree, h = 1, c = 1)
    K <- covm(x, x, covModel, use_symmetry = TRUE)
    manual <- outer(x, x, function(a, b) (a * b + 1)^degree)
    expect_equal(unname(as.matrix(K)), manual, tolerance = 1e-10)
  }
})

test_that("polynomial kernel: non-square case works", {
  x <- c(-4, -3, -1, 0, 4)
  z <- seq(-10, 10, length.out = 50)
  covModel <- list(kernel = "polynomial", degree = 2, h = 1, c = 1)
  K <- covm(x, z, covModel, use_symmetry = FALSE)
  expect_equal(dim(K), c(length(x), length(z)))
  manual <- outer(x, z, function(a, b) (a * b + 1)^2)
  expect_equal(unname(as.matrix(K)), manual, tolerance = 1e-10)
})

test_that("dead kernel_wrapper() was removed and is not reachable", {
  expect_false(exists("kernel_wrapper", where = asNamespace("GauProMod"), inherits = FALSE))
})
