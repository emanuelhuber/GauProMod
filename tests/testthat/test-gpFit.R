obs <- list(x = c(-4, -3, -1, 0, 4), y = c(-2, 0, 1, 2, 0))
targ <- list(x = seq(-10, 10, length.out = 50))

test_that("gpLogLik matches gpCond()$logLik exactly (zero-mean case)", {
  covModel <- list(kernel = "gaussian", l = 1.3, h = 1.7)
  A <- gpCond(obs, targ, list(covModel), sigma = 0.2, onlyMean = FALSE)
  ll_lean <- gpLogLik(obs, list(covModel), sigma = 0.2)
  expect_equal(ll_lean, A$logLik, tolerance = 1e-10)
})

test_that("gpLogLik matches gpCond()$logLik exactly (monomial mean function, op=1)", {
  covModel <- list(kernel = "gaussian", l = 1.3, h = 1.7)
  A <- gpCond(obs, targ, list(covModel), sigma = 0.2, op = 1, onlyMean = FALSE)
  ll_lean <- gpLogLik(obs, list(covModel), sigma = 0.2, op = 1)
  expect_equal(ll_lean, A$logLik, tolerance = 1e-10)
})

test_that("gpFit() converges and improves on the starting log-likelihood", {
  set_params <- function(theta) {
    list(covModels = list(list(kernel = "gaussian",
                                l = exp(theta[["log_l"]]),
                                h = exp(theta[["log_h"]]))),
         sigma = exp(theta[["log_sigma"]]))
  }
  theta0 <- c(log_l = log(1), log_h = log(1), log_sigma = log(0.5))
  ll0 <- -gpNegLogLik(theta0, obs, set_params)

  fit <- gpFit(obs, theta0, set_params, method = "Nelder-Mead")
  expect_equal(fit$convergence, 0)
  expect_gte(fit$logLik, ll0)
})

test_that("gpFit() gradient is near zero at a well-conditioned optimum", {
  set.seed(1)
  n <- 25
  x <- sort(runif(n, -10, 10))
  y <- sin(x) + 0.3 * x + rnorm(n, sd = 0.15)
  obs2 <- list(x = x, y = y)

  set_params <- function(theta) {
    list(covModels = list(list(kernel = "gaussian",
                                l = exp(theta[["log_l"]]),
                                h = exp(theta[["log_h"]]))),
         sigma = exp(theta[["log_sigma"]]))
  }
  theta0 <- c(log_l = log(2), log_h = log(1), log_sigma = log(0.3))
  fit <- gpFit(obs2, theta0, set_params, method = "L-BFGS-B")
  expect_equal(fit$convergence, 0)

  g <- GauProMod:::.centralDiffGrad(gpNegLogLik, fit$par,
                                     obs = obs2, set_params = set_params,
                                     op = 0, bc = NULL, penalty = 1e10)
  expect_lt(sqrt(sum(g^2)), 1e-3)
})

test_that("gpFit() can hold a parameter fixed (e.g. matern v)", {
  set.seed(1)
  n <- 25
  x <- sort(runif(n, -10, 10))
  y <- sin(x) + 0.3 * x + rnorm(n, sd = 0.15)
  obs2 <- list(x = x, y = y)

  set_params <- function(theta) {
    list(covModels = list(list(kernel = "matern",
                                l = exp(theta[["log_l"]]),
                                h = exp(theta[["log_h"]]),
                                v = 1.5)),                # fixed, not in theta
         sigma = exp(theta[["log_sigma"]]))
  }
  theta0 <- c(log_l = log(2), log_h = log(1), log_sigma = log(0.3))
  fit <- gpFit(obs2, theta0, set_params, method = "L-BFGS-B")
  expect_equal(fit$fitted$covModels[[1]]$v, 1.5)
  expect_setequal(names(fit$par), c("log_l", "log_h", "log_sigma"))
})

test_that("gpFit() validates the 'gr' argument", {
  set_params <- function(theta) {
    list(covModels = list(list(kernel = "gaussian",
                                l = exp(theta[["log_l"]]),
                                h = exp(theta[["log_h"]]))),
         sigma = exp(theta[["log_sigma"]]))
  }
  theta0 <- c(log_l = log(1), log_h = log(1), log_sigma = log(0.5))
  expect_error(gpFit(obs, theta0, set_params, gr = "not-a-valid-option"))
})

test_that("gpNegLogLik returns the penalty (not an error) for a non-PD trial point", {
  # duplicate observation locations with zero noise -> singular Kxx
  obs_dup <- list(x = c(0, 0, 1, 2), y = c(1, 1, 2, 3))
  set_params <- function(theta) {
    list(covModels = list(list(kernel = "gaussian", l = theta[["l"]], h = theta[["h"]])),
         sigma = 0)
  }
  val <- gpNegLogLik(c(l = 1, h = 1), obs_dup, set_params, penalty = 12345)
  expect_equal(val, 12345)
})
