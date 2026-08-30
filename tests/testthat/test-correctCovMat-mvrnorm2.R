make_bad_sigma <- function() {
  set.seed(42)
  p <- 6
  A <- matrix(rnorm(p * p), p, p)
  Sigma_pd <- A %*% t(A) + diag(p)
  scales <- c(1, 4, 9, 0.25, 16, 2)
  Sigma_pd <- diag(sqrt(scales)) %*% cov2cor(Sigma_pd) %*% diag(sqrt(scales))
  Sigma_bad <- Sigma_pd
  Sigma_bad[1, 2] <- Sigma_bad[2, 1] <- Sigma_bad[1, 2] * 3
  Sigma_bad
}

test_that("correctCovMat returns a lower-triangular L with L %*% t(L) approx Sigma", {
  Sigma_bad <- make_bad_sigma()
  eig <- eigen(Sigma_bad, symmetric = TRUE, only.values = TRUE)$values
  expect_lt(min(eig), 0)  # confirm the input really isn't PD

  L <- GauProMod:::correctCovMat(Sigma_bad)
  expect_true(all(L[upper.tri(L)] == 0))
  Sigma_repaired <- L %*% t(L)
  expect_true(min(eigen(Sigma_repaired, symmetric = TRUE, only.values = TRUE)$values) >= -1e-12)
})

test_that("correctCovMat preserves the original variances (does not silently rescale to a correlation matrix)", {
  Sigma_bad <- make_bad_sigma()
  L <- GauProMod:::correctCovMat(Sigma_bad)
  Sigma_repaired <- L %*% t(L)
  expect_equal(diag(Sigma_repaired), diag(Sigma_bad), tolerance = 1e-8)
})

test_that("correctCovMat's maxit is actually enforced (regression: unbounded while loop)", {
  Sigma_bad <- make_bad_sigma()
  expect_error(GauProMod:::correctCovMat(Sigma_bad, maxit = 0))
})

test_that("correctCovMat rejects a non-positive diagonal", {
  bad <- matrix(c(-1, 0, 0, 1), 2, 2)
  expect_error(GauProMod:::correctCovMat(bad))
})

test_that("mvrnorm2 works for n > 1 (regression: non-conformable arguments crash)", {
  set.seed(1)
  p <- 5
  A <- matrix(rnorm(p * p), p, p)
  Sigma <- A %*% t(A) + diag(p)
  mu <- c(1, -2, 0.5, 3, -1)

  samp <- mvrnorm2(2000, mu, Sigma)
  expect_equal(dim(samp), c(p, 2000))
  expect_equal(rowMeans(samp), mu, tolerance = 0.1)
  expect_equal(cov(t(samp)), Sigma, tolerance = 0.15)
})

test_that("mvrnorm2 still works for n = 1", {
  set.seed(1)
  p <- 4
  A <- matrix(rnorm(p * p), p, p)
  Sigma <- A %*% t(A) + diag(p)
  mu <- rep(0, p)
  samp <- mvrnorm2(1, mu, Sigma)
  expect_equal(dim(samp), c(p, 1))
})

test_that("mvrnorm2's correctCovMat fallback produces samples with the right covariance (regression: wrong triangular convention)", {
  set.seed(1)
  Sigma_bad <- make_bad_sigma()
  mu <- rep(0, nrow(Sigma_bad))
  n <- 20000
  samp <- mvrnorm2(n, mu, Sigma_bad)
  L <- GauProMod:::correctCovMat(Sigma_bad)
  Sigma_repaired <- L %*% t(L)
  expect_equal(cov(t(samp)), Sigma_repaired, tolerance = 0.15)
})

test_that("gpSim() output has the right dimensions for n > 1 (regression: non-conformable arguments crash)", {
  obs <- list(x = c(-4, -3, -1, 0, 4), y = c(-2, 0, 1, 2, 0))
  targ <- list(x = seq(-10, 10, length.out = 15))
  A <- gpCond(obs, targ, list(list(kernel = "gaussian", l = 2, h = 1.5)), sigma = 0.2)

  n <- 300
  sims <- gpSim(A, n = n)
  expect_equal(dim(sims), c(length(A$xstar), 1 + n))

  sim_vals <- sims[, -1]
  expect_equal(rowMeans(sim_vals), A$mean, tolerance = 0.2)
})

test_that("gpSim(L = ...) branch produces samples with the right mean", {
  obs <- list(x = c(-4, -3, -1, 0, 4), y = c(-2, 0, 1, 2, 0))
  targ <- list(x = seq(-10, 10, length.out = 15))
  A <- gpCond(obs, targ, list(list(kernel = "gaussian", l = 2, h = 1.5)), sigma = 0.2)

  L <- t(chol(A$cov))  # lower-triangular, L %*% t(L) = cov
  n <- 300
  sims <- gpSim(A, L = L, n = n)
  expect_equal(dim(sims), c(length(A$xstar), 1 + n))
  sim_vals <- sims[, -1]
  expect_equal(rowMeans(sim_vals), A$mean, tolerance = 0.2)
})
