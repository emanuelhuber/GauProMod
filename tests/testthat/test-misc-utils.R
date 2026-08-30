test_that("crossDist works with the documented default use_symmetry=FALSE (regression: self-referential default)", {
  x <- c(-1, 0, 1)
  expect_no_error(D <- crossDist(x, x))
  expect_equal(as.matrix(D), as.matrix(dist(x)), tolerance = 1e-10, ignore_attr = TRUE)
})

test_that("invm matches solve() (regression: fed lower-triangular factor to chol2inv, which expects upper)", {
  set.seed(1)
  p <- 5
  A <- matrix(rnorm(p * p), p, p)
  Sigma <- A %*% t(A) + diag(p)
  expect_equal(invm(Sigma), solve(Sigma), tolerance = 1e-8)
})

test_that("setPosTime's obs/xstar use unique locations/times, not pre-expanded (regression: double-expansion crash in gpCond space-time)", {
  xy <- matrix(c(0, 0, 1, 1), ncol = 2, byrow = TRUE)
  tt <- c(1, 2, 5)  # irregular gaps
  val <- 1:6
  xystar <- matrix(c(0.5, 0.5), ncol = 2)
  res <- setPosTime(xy, tt, val, xystar)
  # x/t must stay UNIQUE (unexpanded) -- gpCond() does its own space-time
  # Kronecker expansion internally
  expect_equal(nrow(res$obs$x), nrow(xy))
  expect_equal(res$obs$t, tt)
  expect_equal(nrow(res$xstar$x), nrow(xystar))
  expect_equal(res$xstar$t, tt)
  expect_length(res$obs$y, nrow(xy) * length(tt))
})

test_that("setPosTime output plugs directly into gpCond() for a space-time model, with correctly-sized output (regression: double-expansion crash)", {
  xy <- matrix(c(0, 0, 1, 1), ncol = 2, byrow = TRUE)
  tt <- c(1, 2, 5)
  val <- c(1, 2, 1.5, 3, 4, 3.5)
  xystar <- matrix(c(0.5, 0.5), ncol = 2)
  res <- setPosTime(xy, tt, val, xystar)

  targ <- list(x = res$xstar$x, t = res$xstar$t)
  covModels <- list(list(kernel = "gaussian", l = 1, h = 1),
                     list(kernel = "gaussian", l = 2, h = 1))
  A <- gpCond(res$obs, targ, covModels, sigma = 0.1, sigmat = 0.1)
  # 1 target location x 3 target times = 3 predictions (previously this
  # came out as 3*3 = 9 due to gpCond() expanding the already-expanded
  # setPosTime() output a second time)
  expect_length(A$mean, nrow(xystar) * length(tt))
  expect_true(all(is.finite(A$mean)))
})

test_that("matGrid/vecGrid produce the documented layout", {
  g <- matGrid(1:3, 10:11)
  expect_equal(g$X, matrix(rep(1:3, each = 2), nrow = 2))
  expect_equal(g$Y, matrix(rep(10:11, times = 3), nrow = 2))

  v <- vecGrid(1:3, 10:11)
  expect_equal(dim(v), c(6, 2))
  expect_equal(v[1, ], c(1, 10))
})

test_that("cholfac returns the lower factor: L %*% t(L) == x", {
  set.seed(1)
  p <- 4
  A <- matrix(rnorm(p * p), p, p)
  Sigma <- A %*% t(A) + diag(p)
  L <- cholfac(Sigma)
  expect_true(all(L[upper.tri(L)] == 0))
  expect_equal(L %*% t(L), Sigma, tolerance = 1e-8)
})
