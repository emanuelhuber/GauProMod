test_that("gaussian kernel: square and non-square cases", {
  x <- c(-4, -3, -1, 0, 4)
  z <- seq(-10, 10, length.out = 30)
  covModel <- list(kernel = "gaussian", l = 1.3, h = 1.7)
  K1 <- covm(x, x, covModel, use_symmetry = TRUE)
  K2 <- covm(x, z, covModel, use_symmetry = FALSE)
  expect_equal(dim(K1), c(5, 5))
  expect_equal(dim(K2), c(5, 30))
  D <- as.matrix(dist(x))
  manual <- 1.7^2 * exp(-0.5 * (D / 1.3)^2)
  expect_equal(unname(as.matrix(K1)), unname(manual), tolerance = 1e-10)
})

test_that("matern kernel-name aliases resolve to matern with the documented fixed v", {
  x <- sort(runif(10, -5, 5))
  aliases <- list(exponential = 0.5, matern_1_2 = 0.5,
                   matern_3_2 = 1.5, matern_5_2 = 2.5)
  for (nm in names(aliases)) {
    K_alias  <- covm(x, x, list(kernel = nm, l = 1, h = 1), use_symmetry = TRUE)
    K_manual <- covm(x, x, list(kernel = "matern", l = 1, h = 1, v = aliases[[nm]]),
                      use_symmetry = TRUE)
    expect_equal(as.matrix(K_alias), as.matrix(K_manual), tolerance = 1e-10,
                 label = paste("alias", nm))
  }
})

# not yet implemented in C++
# test_that("power_exp: d=0 matches h^2*exp(-(r/l)^v) for several v", {
#   x <- c(-3, -1.2, 0, 0.5, 2.8)
#   D <- as.matrix(dist(x))
#   l <- 1.7; h <- 1.3
#   for (v in c(2, 1.5, 1, 0.6)) {
#     K <- covm(x, x, list(kernel = "power_exp", l = l, h = h, v = v), use_symmetry = TRUE)
#     manual <- h^2 * exp(-(D / l)^v)
#     expect_equal(unname(as.matrix(K)), unname(manual), tolerance = 1e-10,
#                  label = paste("v =", v))
#   }
# })
# 
# test_that("power_exp: d=1 matches numerical d/dx (same convention as 'gaussian')", {
#   l <- 1.7; h <- 1.3; xx <- 0.3; yy <- 1.7; eps <- 1e-6
#   for (v in c(2, 1.5, 1, 0.6)) {
#     para <- list(kernel = "power_exp", l = l, h = h, v = v)
#     Kd1 <- covm(xx, yy, para, d = 1, use_symmetry = FALSE)
#     Kp  <- covm(xx + eps, yy, para, d = 0, use_symmetry = FALSE)
#     Km  <- covm(xx - eps, yy, para, d = 0, use_symmetry = FALSE)
#     num <- (Kp - Km) / (2 * eps)
#     expect_equal(as.numeric(Kd1), as.numeric(num), tolerance = 1e-5,
#                  label = paste("v =", v))
#   }
# })
# 
# test_that("power_exp: d=2 matches numerical second derivative away from r=0", {
#   l <- 1.7; h <- 1.3; xx <- 0.3; yy <- 1.7; eps <- 1e-4
#   for (v in c(2, 1.5, 1, 0.6)) {
#     para <- list(kernel = "power_exp", l = l, h = h, v = v)
#     Kd2 <- covm(xx, yy, para, d = 2, use_symmetry = FALSE)
#     r0 <- abs(xx - yy)
#     f <- function(r) h^2 * exp(-(r / l)^v)
#     num <- (f(r0 + eps) - 2 * f(r0) + f(r0 - eps)) / eps^2
#     # this codebase's d=2 convention returns -f''(r) (matches kGaussian)
#     expect_equal(as.numeric(Kd2), -num, tolerance = 1e-4, label = paste("v =", v))
#   }
# })
# 
# test_that("power_exp: d=2 at r=0 is finite only for v in {1, 2}, +/-Inf otherwise", {
#   para <- function(v) list(kernel = "power_exp", l = 1.7, h = 1.3, v = v)
#   expect_true(is.finite(covm(0, 0, para(2), d = 2, use_symmetry = TRUE)))
#   expect_true(is.finite(covm(0, 0, para(1), d = 2, use_symmetry = TRUE)))
#   expect_true(is.infinite(covm(0, 0, para(1.5), d = 2, use_symmetry = TRUE)))
#   expect_true(is.infinite(covm(0, 0, para(0.6), d = 2, use_symmetry = TRUE)))
# })
# 
# test_that("power_exp: d=1 at r=0 is exactly 0 (matches w=0 there)", {
#   for (v in c(2, 1.5, 1, 0.6)) {
#     para <- list(kernel = "power_exp", l = 1.7, h = 1.3, v = v)
#     expect_equal(as.numeric(covm(0, 0, para, d = 1, use_symmetry = TRUE)), 0)
#   }
# })
# 
# test_that("power_exp: v outside (0, 2] and missing v both error clearly", {
#   x <- c(-1, 0, 1)
#   expect_error(covm(x, x, list(kernel = "power_exp", l = 1, h = 1, v = 2.5),
#                      use_symmetry = TRUE))
#   expect_error(covm(x, x, list(kernel = "power_exp", l = 1, h = 1),
#                      use_symmetry = TRUE))
# })
# 
# test_that("power_exp works end-to-end in gpCond() and gpFit()", {
#   obs <- list(x = c(-4, -3, -1, 0, 4), y = c(-2, 0, 1, 2, 0))
#   targ <- list(x = seq(-10, 10, length.out = 20))
#   A <- gpCond(obs, targ, list(list(kernel = "power_exp", l = 2, h = 1.5, v = 1.3)),
#               sigma = 0.2)
#   expect_length(A$mean, 20)
#   expect_true(all(is.finite(A$mean)))
# })
