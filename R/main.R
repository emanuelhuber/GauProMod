######## GAUSSIAN PROCESS FUNCTIONS #########

#' @useDynLib GauProMod
#' @importFrom Rcpp sourceCpp
NULL

# observations:
# obs = list("x" = position,
# 			 "t" = time (optional),
# 			 "y" = observed value)
# taget
# targ = list("x" = position,
# 			 "t" = time (optional))
# bc <- list(x = cbind(c( 0.5, 2.5,  5, 7.5, 9.5),
#                      c(0.5,  0.5, 0.5,  0.5, 0.5)),
#            v = cbind(1*c( 1,   1,  1,   1,  1),
#                      0*c(1,  1, 1,  1, 1)),
#            y =  0* c( 1,   1,  1,   1,  1),
#             sigma = 0)
# monomial functions 1D: op = 1 (linear), op = 2 (quadratic), op = 3 (cube)
# monomial functions 2D:
# op = 2 -> linear mean function
# op = 5 -> quadratic mean function
#
# RETURNS A LIST containing:
# - "mean" > mean function evaluated at xstar
# - "cov" > covariance function evaluated at xstar
# - "xstar"
#' Conditional Gaussian Process simulation
#'
#' Compute the posterior (conditional) mean and covariance of a Gaussian
#' Process at target locations \code{targ}, given noisy observations
#' \code{obs}, one or two covariance models \code{covModels} (spatial
#' only, or spatial x temporal via a Kronecker product), an optional
#' polynomial mean function, and optional derivative/boundary
#' constraints.
#'
#' @param obs list of observations. For a spatial-only GP
#'   (\code{length(covModels) == 1}): \code{x} (locations -- a numeric
#'   vector for 1D, or a matrix with one row per location for 2D/3D) and
#'   \code{y} (observed values, same length/number of rows as \code{x}).
#'   For a space-time GP (\code{length(covModels) == 2}): \code{x} lists
#'   only the \code{nxy} \emph{unique} spatial locations and \code{t}
#'   only the \code{nt} \emph{unique} times (the full space-time
#'   covariance is built internally as a Kronecker product of the
#'   spatial and temporal covariances, not from a pre-expanded
#'   one-row-per-(location,time) input); \code{y} then has the full
#'   \code{nxy * nt} observed values, ordered as all \code{nt} times for
#'   the first location, then all \code{nt} times for the second
#'   location, and so on. See \code{\link{setPosTime}()} for a helper
#'   that builds this space-time layout from "every location observed
#'   at every time" data.
#' @param targ list of target locations to predict at, in the same
#'   convention as \code{obs}: \code{x} (unique locations), and, for a
#'   space-time model, optionally \code{t} (unique target times,
#'   defaults to \code{obs$t}, i.e. predict at the same times as
#'   observed).
#' @param covModels a list of one covariance model (spatial-only GP) or
#'   two covariance models (space-time GP: \code{covModels[[1]]} for
#'   space, \code{covModels[[2]]} for time; the full covariance is their
#'   Kronecker product). Each covariance model is a list as documented
#'   in \code{\link{covm}()}'s \code{covModel} argument.
#' @param sigma observation noise standard deviation: either a single
#'   value (applied to every observation) or a vector matching
#'   \code{length(obs$y)}.
#' @param op order of an optional polynomial mean function fitted
#'   alongside the GP (rather than assuming a zero mean): \code{0}
#'   (default) for no mean function; for 1D \code{obs$x}, \code{1}/
#'   \code{2}/\code{3} for a linear/quadratic/cubic polynomial; for 2D
#'   \code{obs$x}, \code{2} for a linear (\eqn{1, x_1, x_2}) or \code{5}
#'   for a quadratic (\eqn{1, x_1, x_2, x_1^2, x_2^2, x_1 x_2}) surface.
#' @param bc optional derivative/boundary constraints: a list with
#'   \code{x} (locations of the constraints), \code{v} (unit direction
#'   vector(s) the derivative is taken along, for 2D+ locations),
#'   \code{y} (the constrained derivative values, e.g. \code{0} for a
#'   zero-gradient/no-flux boundary), and \code{sigma} (noise on the
#'   constraint). Adds these as extra "observations" of the GP's
#'   derivative, jointly conditioned on with \code{obs}.
#' @param sigmat for a space-time model (\code{length(covModels) == 2}
#'   only), an additional noise standard deviation applied along the
#'   time dimension (added to the diagonal of the temporal noise term).
#' @param onlyMean if \code{TRUE}, skip computing the posterior
#'   covariance and log marginal likelihood and return only the
#'   posterior mean (faster; useful when only point predictions are
#'   needed, e.g. inside a loop). Default \code{FALSE}.
#'
#' @return a list with:
#'   \item{mean}{posterior mean at \code{targ} (length
#'     \code{nrow(targ$x)}, or \code{nrow(targ$x) * length(targ$t)} for
#'     a space-time model).}
#'   \item{cov}{posterior covariance matrix at \code{targ} (omitted if
#'     \code{onlyMean = TRUE}).}
#'   \item{logLik}{log marginal likelihood of \code{obs} under
#'     \code{covModels} (omitted if \code{onlyMean = TRUE}; see
#'     \code{\link{gpLogLik}()} for a version that computes this alone,
#'     without the cost of also predicting at \code{targ} -- the right
#'     choice when fitting hyperparameters, e.g. with
#'     \code{\link{gpFit}()}).}
#'   \item{xstar}{the target locations predictions correspond to, in
#'     the same order as \code{mean}/\code{cov} (for a space-time model,
#'     \code{targ$x} repeated for each time, column-bound with the
#'     matching \code{targ$t}).}
#'
#' @seealso [covm()], [gpLogLik()], [gpFit()], [gpSim()], [setPosTime()]
#' @name gpCond
#' @export
gpCond <- function(obs, targ, covModels, sigma=0, op = 0 , bc = NULL,
                   sigmat = 0, onlyMean = FALSE){
  
  Kxx       <- covm( obs$x,  obs$x, covModels[[1]], use_symmetry = TRUE)
  Kstar     <- covm( obs$x, targ$x, covModels[[1]])
  Kstarstar <- covm(targ$x, targ$x, covModels[[1]], use_symmetry = TRUE)
  if(length(sigma) == 1){
    sigma     <- rep(sigma, ncol(Kxx))
  }else if(length(sigma) != ncol(Kxx)){
    stop("length of sigma must be equal to 1 or to the number of observations")
  }
  y <- obs$y
  
  # if there are derivative
  if(!is.null(bc)){
    # if length(dim (obs$x)) == 1 or   length(dim (obs$y)) == 1
    # dx = 1 (bc$v is Null and will not be used)
    Kdxx  <- covm(obs$x, bc$x, covModels[[1]] , d = 1, dx = bc$v)
    Kdxdx <- covm( bc$x, bc$x, covModels[[1]] , d = 2, dx = bc$v, use_symmetry = TRUE)
    sigma <- c(sigma, rep(bc$sigma, ncol(Kdxx)))
    # Kxx is symmetric, but the derivative blocks are NOT symmetric (K_dxdx != K_dxdx^T)
    # K_obs,dx = Kdxx (needs w = sign(x-x'))
    # K_dx,obs = -Kdxx^T (needs w = sign(x'-x))
    Kxx  <- rbind(cbind(   Kxx,  Kdxx ),
                  cbind(t(-Kdxx), Kdxdx))
    Kdxstar <- covm(targ$x, bc$x, covModels[[1]] , d = 1, dx=bc$v)
    # K_obs,star is Kstar. We need K_dx,star: K_dx,star = -t(K_star,dx)
    Kstar   <- rbind(Kstar, t(-Kdxstar))
    y <- c(y, bc$y)
  }
  # if space-time or space-space GP
  if(length(covModels) == 2){
    if(is.null(targ$t)){
      targ$t <- obs$t
    }
    nt <- length(obs$t)
    xstar <- targ$x[rep(seq_len(nrow(targ$x)),each=nt),]
    nxy <- nrow(obs$x)
    x <- obs$x[rep(seq_len(nxy), each=nt),]
    #     y <- obs$y
    AA <- cbind(xstar,targ$t)
    Ktt <- covm(obs$t, obs$t, covModels[[2]], use_symmetry = TRUE)
    Knoise <- diag(sigma^2)
    Ktnoise <- diag(rep(sigmat^2, length(obs$t)))
    Kxx <-  (Kxx + Knoise) %x% (Ktt + Ktnoise)
    # structure of Kxx: Kxx[1,1]*Knoise[1,1] Kxx[1,1]*Knoise[1,2] ...
    #                   Kxx[1,1]*Knoise[2,1] Kxx[1,1]*Knoise[2,2] ...
    #                   ...
    Kstarstar <- Kstarstar %x%
      covm(targ$t, targ$t, covModels[[2]], use_symmetry = TRUE)
    Kstar <- Kstar %x% covm(obs$t, targ$t, covModels[[2]])
    if(!is.null(bc)){
      # see the matching fix/comment in gpLogLik(): each = nt, not nt
      y <- c(obs$y, rep(bc$y, each = nt))
    }
  }else{
    Kxx <-  Kxx  +  diag(sigma^2)
    x <- obs$x
    xstar <- targ$x
    AA <- targ$x
  }
  # if monomial functions
  if(op > 0){   # monomial mean functions to be estimated
    H <- Hmat(x,op)
    if(!is.null(bc)){
      nbc <- nrow(bc$x)
      if(length(covModels) == 2){
        Hdx <- Hmat(bc$x[rep(seq_len(nbc), each=nt),], op, dx = bc$v)
        H <- cbind(H, Hdx)
      }else{
        Hdx <- Hmat(bc$x, op, dx = bc$v)
        H <- cbind(H, (Hdx))
      }
    }
    Hstar <- Hmat(xstar,op)
    A2 <- GPpredmean_rcpp(Kxx, Kstar, Kstarstar, y, H, Hstar, only_mean = onlyMean)
    # A2[["logLik"]] <- A2$logLik
  }else{
    A2 <- GPpred_rcpp(Kxx, Kstar, Kstarstar, y, only_mean = onlyMean)
    # logLik <- A2$logLik# 1 - sum(log(A2$logLik2)) - nrow(Kxx)/2 * log(2*pi)
    # A2[["logLik"]] <- A2$logLik
  }
  A2[["xstar"]] <- AA
  return(A2)
}



#' Log marginal likelihood of a Gaussian Process model
#'
#' Builds the same observation covariance matrix \code{Kxx} that
#' \code{\link{gpCond}} would (including noise, derivative/boundary
#' constraints \code{bc}, the space-time Kronecker structure, and the
#' monomial mean basis \code{H} when \code{op > 0}), but -- unlike
#' \code{gpCond} -- it never builds \code{Kstar} or \code{Kstarstar}, i.e.
#' it does no work that depends on \code{targ} at all. That makes it the
#' right building block for hyperparameter optimization: each evaluation
#' only pays for what the marginal likelihood actually needs.
#'
#' @param obs list with observation locations \code{x} (and optionally
#'   \code{t}) and observed values \code{y}, as used by \code{\link{gpCond}}.
#' @param covModels list of one (spatial) or two (spatial + temporal)
#'   covariance model specifications, as used by \code{\link{gpCond}}.
#' @param sigma observation noise standard deviation(s); see
#'   \code{\link{gpCond}}.
#' @param op monomial mean function order (0 = no mean function); see
#'   \code{\link{gpCond}}.
#' @param bc optional derivative/boundary constraints; see
#'   \code{\link{gpCond}}.
#' @param sigmat optional time-noise standard deviation for space-time
#'   models; see \code{\link{gpCond}}.
#'
#' @return the scalar log marginal likelihood (\code{-Inf} if the
#'   resulting covariance matrix is not positive definite for the given
#'   hyperparameters -- e.g. a length scale that is too small).
#'
#' @seealso [gpCond()], [gpNegLogLik()], [gpFit()]
#' @name gpLogLik
#' @export
gpLogLik <- function(obs, covModels, sigma = 0, op = 0, bc = NULL, sigmat = 0){
  
  Kxx <- covm(obs$x, obs$x, covModels[[1]], use_symmetry = TRUE)
  if(length(sigma) == 1){
    sigma <- rep(sigma, ncol(Kxx))
  }else if(length(sigma) != ncol(Kxx)){
    stop("length of sigma must be equal to 1 or to the number of observations")
  }
  y <- obs$y
  x <- obs$x
  
  # if there are derivative/boundary constraints
  if(!is.null(bc)){
    Kdxx  <- covm(obs$x, bc$x, covModels[[1]], d = 1, dx = bc$v)
    Kdxdx <- covm( bc$x, bc$x, covModels[[1]], d = 2, dx = bc$v, use_symmetry = TRUE)
    sigma <- c(sigma, rep(bc$sigma, ncol(Kdxx)))
    Kxx  <- rbind(cbind(   Kxx,  Kdxx ),
                  cbind(t(-Kdxx), Kdxdx))
    y <- c(y, bc$y)
  }
  
  # if space-time or space-space GP
  if(length(covModels) == 2){
    nt  <- length(obs$t)
    nxy <- nrow(obs$x)
    x   <- obs$x[rep(seq_len(nxy), each = nt), ]
    Ktt     <- covm(obs$t, obs$t, covModels[[2]], use_symmetry = TRUE)
    Knoise  <- diag(sigma^2)
    Ktnoise <- diag(rep(sigmat^2, length(obs$t)))
    Kxx <- (Kxx + Knoise) %x% (Ktt + Ktnoise)
    if(!is.null(bc)){
      # each = nt (not the bare nt used previously): the Kronecker structure
      # of Kxx blocks each bc row over its nt time repeats (bc row 1 for all
      # nt times, then bc row 2 for all nt times, ...), matching the bc$x
      # expansion used for Hdx below. rep(bc$y, nt) instead cycled the whole
      # vector nt times, misaligning y against Kxx's rows whenever nbc > 1.
      y <- c(obs$y, rep(bc$y, each = nt))
    }
  }else{
    Kxx <- Kxx + diag(sigma^2)
  }
  
  # if monomial mean function
  if(op > 0){
    H <- Hmat(x, op)
    if(!is.null(bc)){
      nbc <- nrow(bc$x)
      if(length(covModels) == 2){
        Hdx <- Hmat(bc$x[rep(seq_len(nbc), each = nt), ], op, dx = bc$v)
      }else{
        Hdx <- Hmat(bc$x, op, dx = bc$v)
      }
      H <- cbind(H, Hdx)
    }
    ll <- gpLogLikMean_rcpp(Kxx, y, H)
  }else{
    ll <- gpLogLik_rcpp(Kxx, y)
  }
  return(ll)
}

#' Negative log marginal likelihood, ready for optim()
#'
#' Thin wrapper around \code{\link{gpLogLik}} with the signature
#' \code{optim()} expects: a numeric parameter vector \code{theta} in,
#' a single scalar out. How \code{theta} maps onto \code{covModels} /
#' \code{sigma} / \code{sigmat} is entirely up to \code{set_params()},
#' so any parameterization works -- e.g. optimizing on a log scale to
#' keep length scales and standard deviations positive.
#'
#' Cholesky failures (non-positive-definite covariance matrices, which
#' happen routinely for bad trial hyperparameters) are caught and turned
#' into \code{penalty} instead of an error, so \code{optim()} can step
#' through them without the whole fit aborting.
#'
#' @param theta numeric parameter vector, on whatever scale
#'   \code{set_params} expects.
#' @param obs observations; see \code{\link{gpLogLik}}.
#' @param set_params function(theta) -> list with elements
#'   \code{covModels}, \code{sigma} and (optionally) \code{sigmat}, i.e.
#'   the arguments \code{\link{gpLogLik}} needs.
#' @param op,bc passed through to \code{\link{gpLogLik}}.
#' @param penalty value returned when the trial \code{theta} yields a
#'   non-positive-definite covariance matrix (default \code{1e10}, i.e.
#'   a very bad but finite objective for the optimizer to move away from).
#'
#' @return scalar negative log marginal likelihood.
#' @seealso [gpLogLik()], [gpFit()]
#' @name gpNegLogLik
#' @export
gpNegLogLik <- function(theta, obs, set_params, op = 0, bc = NULL, penalty = 1e10){
  pars <- set_params(theta)
  sigmat <- if(!is.null(pars$sigmat)) pars$sigmat else 0
  ll <- try(gpLogLik(obs, pars$covModels, sigma = pars$sigma, op = op,
                     bc = bc, sigmat = sigmat),
            silent = TRUE)
  if(inherits(ll, "try-error") || !is.finite(ll)){
    return(penalty)
  }
  return(-ll)
}


#' Fit GP hyperparameters by maximum marginal likelihood
#'
#' Convenience wrapper that runs \code{\link[stats]{optim}} on
#' \code{\link{gpNegLogLik}}, then reports the fitted \code{covModels} /
#' \code{sigma} alongside the usual \code{optim()} output.
#'
#' By default a gradient is supplied to \code{optim()} (see \code{gr}
#' below) rather than leaving \code{gr = NULL}: gradient-based methods
#' like \code{"L-BFGS-B"} still work without one, but \code{optim()}'s
#' fallback is a forward difference with a single fixed absolute step
#' (ndeps, default \code{1e-3}) applied identically to every
#' parameter. That is a coarse approximation once parameters differ in
#' scale (e.g. a log-length-scale near 0 vs. a log-noise-variance near
#' -5), and it is often the reason a fit reports slow or unreliable
#' convergence as the number of hyperparameters grows.
#'
#' @param obs observations; see \code{\link{gpLogLik}}.
#' @param theta0 named numeric vector of starting values, on the scale
#'   \code{set_params} expects.
#' @param set_params function(theta) -> list(covModels=, sigma=,
#'   sigmat=); see \code{\link{gpNegLogLik}}. A typical choice optimizes
#'   on the log scale, e.g. for a Gaussian kernel:
#'   \preformatted{
#'   set_params <- function(theta){
#'     list(covModels = list(list(kernel = "gaussian",
#'                                 l = exp(theta[["log_l"]]),
#'                                 h = exp(theta[["log_h"]]))),
#'          sigma = exp(theta[["log_sigma"]]))
#'   }
#'   }
#' @param op,bc passed through to \code{\link{gpLogLik}}.
#' @param method optimization method passed to \code{optim()} (default
#'   \code{"L-BFGS-B"}).
#' @param gr gradient of the objective, passed to \code{optim()}.
#'   \itemize{
#'     \item \code{"central"} (default): a generic central-difference
#'       gradient of \code{\link{gpNegLogLik}}, with a per-parameter step
#'       scaled to that parameter's magnitude (\code{grad_rel_step},
#'       \code{grad_abs_step} below). Works for any kernel/\code{set_params}
#'       with no kernel-specific derivative code, and -- because
#'       \code{\link{gpLogLik}} is cheap -- costs little relative to the
#'       gain in convergence reliability. This is a numerical, not
#'       analytic, gradient: for kernels/methods where an exact analytic
#'       gradient is worth deriving, supply it via the next option instead.
#'     \item a function \code{function(theta, ...)} returning the gradient
#'       directly: use this to supply an exact analytic gradient (e.g.
#'       derived from \eqn{\partial K/\partial \theta} for your kernel).
#'     \item \code{NULL}: let \code{optim()} use its own default.
#'   }
#' @param grad_rel_step,grad_abs_step step-size controls used when
#'   \code{gr = "central"}: the step for parameter \code{i} is
#'   \code{max(abs(theta[i]) * grad_rel_step, grad_abs_step)}.
#' @param ... further arguments passed to \code{optim()} (e.g.
#'   \code{lower}, \code{upper}, \code{control}).
#'
#' @return the list returned by \code{optim()}, with two extra elements:
#'   \code{fitted} (the result of \code{set_params(optim_par)}, i.e. the
#'   fitted \code{covModels}/\code{sigma}/\code{sigmat}) and \code{logLik}
#'   (the maximized log marginal likelihood, \code{-value}).
#'
#' @seealso [gpLogLik()], [gpNegLogLik()], [gpCond()]
#' @name gpFit
#' @export
gpFit <- function(obs, theta0, set_params, op = 0, bc = NULL,
                  method = "L-BFGS-B", gr = "central",
                  grad_rel_step = 1e-4, grad_abs_step = 1e-6,
                  ..., penalty = 1e10){
  
  gr_fun <- NULL
  if(is.function(gr)){
    gr_fun <- gr
  }else if(identical(gr, "central")){
    gr_fun <- function(theta, obs, set_params, op = 0, bc = NULL, penalty = 1e10){
      .centralDiffGrad(gpNegLogLik, theta,
                       obs = obs, set_params = set_params, op = op, bc = bc,
                       penalty = penalty,
                       rel_step = grad_rel_step, abs_step = grad_abs_step)
    }
  }else if(!is.null(gr)){
    stop('gr must be "central", a gradient function, or NULL')
  }
  
  fit <- optim(par = theta0, fn = gpNegLogLik, gr = gr_fun,
               obs = obs, set_params = set_params, op = op, bc = bc,
               penalty = penalty, method = method, ...)
  fit$fitted <- set_params(fit$par)
  fit$logLik <- -fit$value
  return(fit)
}

#' Central-difference gradient of a scalar function
#'
#' Generic numerical-gradient helper used by \code{\link{gpFit}} when
#' \code{gr = "central"}. Central differences have \eqn{O(h^2)} error
#' (vs. \eqn{O(h)} for the forward differences \code{optim()} falls back
#' to), and the step for each parameter is scaled to that parameter's own
#' magnitude rather than using one fixed absolute step for every
#' parameter -- both matter once hyperparameters span different scales
#' (e.g. a log-length-scale near 0 next to a log-noise-variance near -5).
#'
#' @param fn function(theta, ...) -> scalar.
#' @param theta numeric vector at which to evaluate the gradient.
#' @param ... further arguments passed on to \code{fn}.
#' @param rel_step,abs_step per-parameter step is
#'   \code{max(abs(theta[i]) * rel_step, abs_step)}.
#' @return numeric gradient vector, same length/names as \code{theta}.
#' @keywords internal
.centralDiffGrad <- function(fn, theta, ..., rel_step = 1e-4, abs_step = 1e-6){
  h <- pmax(abs(theta) * rel_step, abs_step)
  g <- numeric(length(theta))
  for(i in seq_along(theta)){
    thp <- theta; thp[i] <- thp[i] + h[i]
    thm <- theta; thm[i] <- thm[i] - h[i]
    g[i] <- (fn(thp, ...) - fn(thm, ...)) / (2 * h[i])
  }
  names(g) <- names(theta)
  return(g)
}


#' Lower Cholesky factor
#'
#' Compute the lower Cholesky factor \code{L} of a symmetric
#' positive-definite matrix \code{x}, such that \code{x = L \%*\% t(L)}.
#' (Note: this is the transpose of base R's \code{\link[base]{chol}()},
#' which returns the upper factor \code{U} with \code{x = t(U) \%*\% U}.)
#'
#' @param x a symmetric positive-definite matrix.
#' @return the lower Cholesky factor \code{L} (\code{L \%*\% t(L) = x}).
#' @name cholfac
#' @export
cholfac <- function(x){
  #   return(cholnew_rcpp(x))
  storage.mode(x) <- "numeric"
  return(cholfac_rcpp(x))
}

#
#
# # Cholesky update
# cholUp <- function(L, x) {
#   p <- nrow(R)
#   stopifnot(is.matrix(R) && p==ncol(R))
#   stopifnot(is.numeric(x) && length(x)==p)
#   L <- .Fortran(dchud, R, p, p, x, 0, 0, 0, 0, 0, numeric(p), numeric(p))
#   return(L[[1]])
# }
#
# # Computes Q such that Q^T Q = R^T R - x x^T.
#
# # Cholesky downdate
# CholDo <- function(R, x) {
#   p <- as.integer(nrow(R))
#   z <- as.integer(0)
#   R <- as.matrix(R)
#   x <- as.numeric(x)
#   stopifnot(p==ncol(R) && p==length(x))
#   L <- .Fortran(dchdd, R, p, p, x, z, z, z, z, z, numeric(p), numeric(p),
#                 integer(1))
#   info <- L[[12]]
#   if (info==-1)
#     stop("downdating produced a non-positive-definite matrix")
#   return(L[[1]])
# }




# L = Cholesky factor (lower matrix)
# A = list with mean and covariance
#' Simulate a Gaussian Process
#'
#' Draw random realizations from a (conditional) Gaussian Process, given
#' its mean and covariance as computed by \code{\link{gpCond}()}.
#'
#' @param A a list with elements \code{mean} (length \code{p}),
#'   \code{cov} (\code{p x p} covariance matrix) and \code{xstar}
#'   (target locations) -- typically the output of \code{\link{gpCond}()}
#'   called with \code{onlyMean = FALSE}.
#' @param L optional pre-computed lower Cholesky factor of \code{A$cov}
#'   (\code{L \%*\% t(L) = A$cov}, as returned by \code{\link{cholfac}()}).
#'   If supplied, this is used directly instead of resampling it from
#'   \code{A$cov} via \code{\link{mvrnorm2}()} -- useful when drawing
#'   from the same covariance repeatedly, to avoid recomputing the
#'   factorization each time.
#' @param n number of realizations to draw (default 1).
#' @return a matrix with \code{n + 1} columns: the first column is
#'   \code{A$xstar} (the locations), followed by \code{n} columns of
#'   simulated values at those locations (one realization per column).
#' @seealso [gpCond()], [mvrnorm2()]
#' @name gpSim
#' @export
gpSim <- function(A, L = NULL, n = 1){
  if(is.null(L)){
    ystar <-  try(mvrnorm2(n, A$mean, A$cov),silent=TRUE)
    if(inherits(ystar, "try-error")){
      cat("Error with Cholesky decomposition...\n")
      # MASS::mvrnorm() returns samples as rows (an n x p matrix, or --
      # when n == 1 -- a bare length-p vector due to R's automatic
      # dimension dropping). Reshape to the p x n ("locations x
      # samples") convention used everywhere else here / in mvrnorm2().
      raw <- MASS::mvrnorm(n, A$mean, A$cov)
      p <- length(A$mean)
      ystar <- if(is.matrix(raw)) t(raw) else matrix(raw, nrow = p, ncol = n)
    }
  }else{
    p <- length(A$mean)
    std <- matrix(rnorm(p * n), nrow = p, ncol = n)
    ystar <- A$mean + (L %*% std)
  }
  return(cbind(A$xstar,ystar))
}



#' Reshape space-time observations and targets for gpCond()
#'
#' Helper to build the \code{obs}/\code{targ} lists that
#' \code{\link{gpCond}()} expects for a space-time model
#' (\code{length(covModels) == 2}), from a "wide" input where every
#' spatial location in \code{xy} was observed at every time in \code{tt}.
#'
#' @param xy matrix of spatial observation locations, one row per
#'   location (\code{nxy} rows).
#' @param tt numeric vector of the \code{nt} observation times (may be
#'   irregularly spaced; \code{NA}s are dropped when computing gaps).
#' @param val observed values, length \code{nxy * nt}, ordered as all
#'   \code{nt} times for the first location in \code{xy}, then all
#'   \code{nt} times for the second location, and so on (i.e.
#'   \code{val <- as.vector(t(observation_matrix))} for an
#'   \code{nxy x nt} matrix of observations).
#' @param xystar matrix of target spatial locations to predict at (same
#'   number of columns as \code{xy}), predicted at the same \code{nt}
#'   time points as the observations.
#' @return a list with elements \code{obs} and \code{xstar}, each a list
#' Reshape space-time observations and targets for gpCond()
#'
#' Helper to build the \code{obs}/\code{targ} lists that
#' \code{\link{gpCond}()} expects for a space-time model
#' (\code{length(covModels) == 2}), from a "wide" input where every
#' spatial location in \code{xy} was observed at every time in \code{tt}.
#'
#' IMPORTANT: for a space-time model, \code{gpCond()} expects \code{x}
#' and \code{t} to each list only the \emph{unique} spatial locations /
#' time points (\code{nxy} rows / \code{nt} values respectively) -- it
#' builds the full space-time covariance internally as a Kronecker
#' product of the purely-spatial and purely-temporal covariances, rather
#' than expecting a pre-expanded one-row-per-(location,time) input.
#' Only \code{y} (the observed values) needs the full \code{nxy * nt}
#' length. Passing pre-expanded locations/times (as an earlier version
#' of this function did) causes \code{gpCond()} to expand them a second
#' time internally, silently producing wrong (oversized) covariance
#' matrices -- this can be a numerically severe bug: an all-C++ Cholesky/
#' Kronecker computation on a wrongly-shaped matrix can even crash R
#' (segfault) rather than raising a clean R-level error.
#'
#' @param xy matrix of spatial observation locations, one row per
#'   location (\code{nxy} rows).
#' @param tt numeric vector of the \code{nt} observation times (may be
#'   irregularly spaced; \code{NA}s are dropped when computing gaps).
#' @param val observed values, length \code{nxy * nt}, ordered as all
#'   \code{nt} times for the first location in \code{xy}, then all
#'   \code{nt} times for the second location, and so on (i.e.
#'   \code{val <- as.vector(t(observation_matrix))} for an
#'   \code{nxy x nt} matrix of observations).
#' @param xystar matrix of target spatial locations to predict at (same
#'   number of columns as \code{xy}), predicted at the same \code{nt}
#'   time points as the observations.
#' @return a list with elements \code{obs} (with \code{x}: the
#'   \code{nxy} unique locations, \code{t}: the \code{nt} unique times,
#'   and \code{y}: the \code{nxy * nt} observed values) and \code{xstar}
#'   (with \code{x}: the unique target locations and \code{t}: the same
#'   \code{nt} time points) -- pass these directly as \code{gpCond()}'s
#'   \code{obs} argument and as \code{targ = list(x = xstar$x, t = xstar$t)}.
#' @examples
#' xy <- matrix(c(0, 0, 1, 1), ncol = 2, byrow = TRUE)  # 2 locations
#' tt <- c(1, 2, 5)                                     # 3 times
#' val <- c(1, 2, 1.5, 3, 4, 3.5)                        # 2 * 3 = 6 values
#' xystar <- matrix(c(0.5, 0.5), ncol = 2)               # 1 target location
#' rt <- setPosTime(xy, tt, val, xystar)
#' covModels <- list(list(kernel = "gaussian", l = 1, h = 1),
#'                    list(kernel = "gaussian", l = 2, h = 1))
#' A <- gpCond(rt$obs, list(x = rt$xstar$x, t = rt$xstar$t), covModels,
#'             sigma = 0.1, sigmat = 0.1)
#' @name setPosTime
#' @export
setPosTime <-function(xy, tt, val, xystar){
  tsteps <- c(1, 1+cumsum(diff(tt[!is.na(tt)])))
  nxy <- nrow(xy)
  nt <- length(tsteps)
  if(nt != length(tt)){
    stop("problem with time\n")
  }
  # observations: x and t list only the UNIQUE locations/times (gpCond()
  # does its own space-time expansion internally via a Kronecker
  # product); only y needs the full nxy * nt length.
  obs <- list()
  obs$x <- xy
  obs$t <- tsteps
  obs$y <- val
  # xstar (points at which we want to predict, at the SAME time points
  # as the observations): likewise unexpanded.
  xstar <- list()
  xstar$x <- xystar
  xstar$t <- tsteps
  
  return(list("obs"=obs, "xstar"=xstar))
}

#' Multi-variate Gaussian simulation
#'
#' A more robust alternative to \code{MASS::mvrnorm()}: attempts a
#' Cholesky factorization of \code{Sigma} directly, and falls back to
#' \code{\link{correctCovMat}} (nearest-PD repair) if \code{Sigma} isn't
#' exactly positive definite (as can happen with covariance matrices
#' built from numerically near-singular kernels).
#'
#' @param n number of samples to draw.
#' @param mu mean vector, length \code{p}.
#' @param Sigma covariance matrix, \code{p x p}.
#' @return a \code{p x n} matrix: one length-\code{p} sample per column.
#'   (Note: this is the transpose of \code{MASS::mvrnorm()}'s
#'   convention, which returns \code{n x p} with one sample per row.)
#' @name mvrnorm2
#' @export
mvrnorm2 <- function(n, mu, Sigma){
  p <- length(mu)
  # 	cholStatus <- try(SChol <- chol(Sigma),silent=TRUE)
  cholStatus <- try(SChol <- cholfac_rcpp(Sigma),silent=TRUE)
  cholError <- inherits(cholStatus, "try-error")
  if(cholError){
    SChol <- correctCovMat(Sigma)
  }
  # one length-p N(0,I) draw per column, so a single matrix
  # multiplication transforms all n samples at once
  std <- matrix(rnorm(p * n), nrow = p, ncol = n)
  realz <- SChol %*% std
  return(mu + realz)
}


# 	http://comisef.wikidot.com/tutorial:repairingcorrelation
# http://www.r-bloggers.com/fixing-non-positive-definite-correlation-
# matrices-using-r-2/
#  Rebonato and Jackel, "The most general methodology for creating a valid
# correlation matrix for risk management and option pricing purposes",
# Journal of Risk, Vol 2, No 2, 2000
#
#' Repair a near positive-definite covariance matrix
#'
#' Nearest-correlation-matrix repair (Rebonato & Jackel: clip negative
#' eigenvalues to zero, reconstruct, iterate) applied to \code{Sigma},
#' returning the lower Cholesky factor \code{L} such that
#' \code{L \%*\% t(L) approx Sigma} -- the same convention used by
#' \code{\link{cholfac}}/\code{cholfac_rcpp}, so the output of this
#' function is a drop-in fallback wherever those are used (as it is in
#' \code{\link{mvrnorm2}}).
#'
#' The Rebonato-Jackel algorithm operates on *correlation* matrices. To
#' avoid silently discarding \code{Sigma}'s actual variances, this
#' function rescales \code{Sigma} to a correlation matrix first, repairs
#' that, and rescales the repaired factor back to \code{Sigma}'s original
#' scale before returning it.
#'
#' @param Sigma covariance matrix to repair (must have a strictly
#'   positive diagonal).
#' @param maxit maximum number of eigenvalue-clipping iterations before
#'   giving up (default 100); previously this loop had no cap and could
#'   spin forever on a pathological input.
#' @param jitter diagonal jitter tried as a last resort if \code{maxit}
#'   iterations of eigenvalue clipping still don't yield a matrix that
#'   \code{chol()} accepts (can happen from floating-point error even
#'   after clipping).
#' @return the lower Cholesky factor \code{L} (\code{L \%*\% t(L) approx
#'   Sigma}).
#' @name correctCovMat
#' @export
correctCovMat <- function(Sigma, maxit = 100, jitter = 1e-10){
  p <- nrow(Sigma)
  if(any(diag(Sigma) <= 0) || anyNA(diag(Sigma))){
    stop("Sigma must have a strictly positive diagonal (variances).")
  }
  d <- sqrt(diag(Sigma))
  # work on the correlation matrix so the repair doesn't discard
  # Sigma's original variances
  R <- Sigma / outer(d, d)
  
  iter <- 0
  cholError <- TRUE
  newR <- R
  while(cholError && iter < maxit){
    iter <- iter + 1
    # compute eigenvectors/-values
    E <- eigen(newR, symmetric = TRUE)
    # replace negative eigenvalues by zero
    E$values <- pmax(E$values, 0)
    # reconstruct correlation matrix
    newR <- E$vectors %*% diag(E$values, nrow = p) %*% t(E$vectors)
    newR <- newR / sqrt(diag(newR) %*% t(diag(newR)))
    cholStatus <- try(u <- chol(newR), silent = TRUE)
    cholError <- inherits(cholStatus, "try-error")
  }
  
  if(cholError){
    # eigenvalue clipping alone didn't converge to something chol()
    # accepts within maxit iterations; try a small diagonal jitter
    # before giving up.
    cholStatus <- try(u <- chol(newR + diag(jitter, p)), silent = TRUE)
    cholError <- inherits(cholStatus, "try-error")
    if(cholError){
      stop("correctCovMat: failed to repair Sigma into a positive ",
           "definite matrix after ", maxit,
           " eigenvalue-clipping iterations plus jitter.")
    }
  }
  
  # u is upper triangular with t(u) %*% u = R (base R chol() convention).
  # Rescale back to Sigma's original variances and convert to the
  # lower-triangular L / L %*% t(L) = Sigma convention used elsewhere
  # in the package: if V = u %*% diag(d) (still upper triangular, since
  # right-multiplying by a diagonal matrix only rescales columns), then
  # t(V) %*% V = diag(d) %*% t(u) %*% u %*% diag(d) = diag(d) %*% R %*% diag(d) = Sigma,
  # so L = t(V) = diag(d) %*% t(u) satisfies L %*% t(L) = Sigma.
  L <- diag(d, nrow = p) %*% t(u)
  return(L)
}


# # for point aligned in a grid (constant dx and dy)
# # The FFT method improves further over CHOL:
# mvrnormGrid <- function(n, mu, Sigma){
# 	M <- length(mu)
# 	N <- n
# 	d <- dim(SIGMA.Y.ext.row1)
# 	dp <- prod(d)
# 	sdp <- sqrt(dp)
# 	prefix <- sqrt(Re(fft(SIGMA.Y.ext.row1, TRUE)))
# 	t2 <- Sys.time()
# 	std <- rnorm(dp)
# 	realz <- prefix * (fft(matrix(std, d[1], d[2]))/sdp)
# 	realz <- as.vector(Re(fft(realz, TRUE)/sdp)[1:M, 1:N])
# 	realz[!inside.owin(x = cent[, 1], y = cent[, 2], w = W)] <- NA
# 	realization.fft.29 <- matrix(realz, M, N, byrow = TRUE)
# }

## for GP with polynomial basis functions 1D
# where the H matrix collects the h(x) vectors for all training
# h(x) are a set of fixed basis functions
# return a matrix of nrow = op + 1 and ncol= number of observation
# dx for derivative!!
Hmat <- function(x,op, dx = NULL){
  # 1D
  if(is.null(dim(x))){
    if(!(op %in% c(1,2,3))){
      stop("Polynomial order should be 1, 2 or 3\n")
    }
    if(is.null(dx)){
      return(t(sapply(0:op,function(a,x) x^a,x)))
    }else{
      HH <- t(sapply(0:op,function(a,x){
        ifelse(x == 0, 0, a*x^(a-1))
      },x))
      if(length(x) == 1){
        dim(HH) <- c(op+1,1)
      }
      return(HH)
    }
    # 2D
  }else{
    if(!(op %in% c(2,5))){
      stop("Polynomial order should be 2 or 5\n")
    }
    if(is.null(dx)){
      HH <- matrix(1,ncol=nrow(x), nrow=op+1)
      if(op ==2 || op == 5){
        HH[2,] <- x[,1]
        HH[3,] <- x[,2]
      }
      if(op == 5){
        HH[4,] <- x[,1]^2
        HH[5,] <- x[,2]^2
        HH[6,] <- x[,1]*x[,2]
      }
    }else{
      HH <- matrix(0,ncol=nrow(x), nrow=op+1)
      dxn <- sqrt(apply(dx^2,1,sum))
      if(op ==2 || op == 5){
        HH[2,] <- 1*dx[,1]/dxn
        HH[3,] <- 1*dx[,2]/dxn
      }
      if(op == 5){
        HH[4,] <- 2*x[,1]*dx[,1]/dxn
        HH[5,] <- 2*x[,2]*dx[,2]/dxn
        HH[6,] <- dx[,1]*x[,2]/dxn + dx[,2]*x[,1]/dxn
      }
    }
    return(HH)
  }
}

# derivative of Hmat
dHmat <- function(x,op, dx=c(1,1)){
  # 1D
  if(is.null(dim(x))){
    if(!(op %in% c(1,2,3))){
      stop("Polynomial order should be 1, 2 or 3\n")
    }
    HH <- t(sapply(0:op,function(a,x){
      ifelse(x == 0, 0, a*x^(a-1))
    },x))
    if(length(x) == 1){
      dim(HH) <- c(op+1,1)
      return(HH)
    }else{
      return(HH)
    }
    # 2D
  }else{
    if(!(op %in% c(2,5))){
      stop("Polynomial order should be 2 or 5\n")
    }
    HH <- matrix(0,ncol=nrow(x), nrow=op+1)
    dxn <- sqrt(apply(dx^2,1,sum))
    if(op ==2 || op == 5){
      HH[2,] <- 1*dx[,1]/dxn
      HH[3,] <- 1*dx[,2]/dxn
    }
    if(op == 5){
      HH[4,] <- 2*x[,1]*dx[,1]/dxn
      HH[5,] <- 2*x[,2]*dx[,2]/dxn
      HH[6,] <- dx[,1]*x[,2]/dxn + dx[,2]*x[,1]/dxn
    }
    return(HH)
  }
}


# x <- obs$x
# t(sapply(0:op,function(a,x) x^a,x))


## MISC
# return two matrix, one for x, one for y
#' Create a 2D coordinate grid
#'
#' Build the two coordinate matrices (\code{X}, \code{Y}) of the
#' rectangular grid formed by every combination of \code{x} and
#' \code{y}, analogous to \code{base::expand.grid()} but returning
#' matrices instead of a data frame (useful e.g. for \code{image()}/
#' \code{contour()}-style plotting).
#'
#' @param x numeric vector of x-coordinates (columns of the grid).
#' @param y numeric vector of y-coordinates (rows of the grid).
#' @return a list with elements \code{X} and \code{Y}, each an
#'   \code{length(y) x length(x)} matrix: \code{X[i,j] = x[j]} and
#'   \code{Y[i,j] = y[i]}.
#' @examples
#' matGrid(1:3, 1:2)
#' @name matGrid
#' @export
matGrid <- function(x,y){
  if (!is.numeric(x) || !is.numeric(y))
    stop("Arguments 'x' and 'y' must be numeric vectors.")
  x <- c(x)
  y <- c(y)
  n <- length(x)
  m <- length(y)
  X <- matrix(rep(x, each = m), nrow = m, ncol = n)
  Y <- matrix(rep(y, times = n), nrow = m, ncol = n)
  return(list(X = X, Y = Y))
}

#' Create a 2D coordinate grid as a two-column matrix
#'
#' Like \code{\link{matGrid}()}, but returns the grid of every
#' combination of \code{x} and \code{y} flattened into a single
#' two-column matrix -- the format expected as \code{x}/\code{obs$x}/
#' \code{targ$x} by \code{\link{covm}()}/\code{\link{gpCond}()} for 2D
#' problems.
#'
#' @param x numeric vector of x-coordinates.
#' @param y numeric vector of y-coordinates.
#' @return a \code{(length(x)*length(y)) x 2} matrix, one row per grid
#'   point, with columns \code{x} and \code{y}.
#' @examples
#' vecGrid(1:3, 1:2)
#' @name vecGrid
#' @export
vecGrid <- function(x,y){
  XY <- matGrid(x,y)
  A <- matrix(nrow=length(XY$X),ncol=2)
  A[,1] <- as.vector(XY$X)
  A[,2] <- as.vector(XY$Y)
  return(A)
}



# Cholesky decomposition provides an effcient and numerically stable method
# for solving equations of the form AX = Y when A is a symmetric,
# positive-definite matrix. The modified Cholesky decomposition is even
# better, because it avoids taking scalar square roots. It is the
# recommended method for forming the term (HPH’ + R)^-1 H
# in the conventional Kalman filter without explicitly inverting a matrix.
# That is, if one decomposes HPH’ + R as UDU’ , then
# (UDU’)(HPH’ + R)^-1 H = H.
# It then suffices to solve
# UDU’X = H
# for X.
# Inverse of a Positive Semi-definite matrix
# Inverse of a Positive Semi-definite matrix
#' Inverse of a symmetric positive-definite matrix
#'
#' Inverts \code{x} via its Cholesky factorization
#' (\code{\link{cholfac}()}, faster and more numerically stable than a
#' general-purpose inverse for a symmetric positive-definite matrix),
#' falling back to a direct matrix inversion if the Cholesky
#' factorization fails (e.g. \code{x} is not positive definite).
#'
#' @param x a symmetric positive-definite matrix.
#' @return the inverse of \code{x}.
#' @name invm
#' @export
invm <- function(x){
  cholx <- try(cholfac(x),silent=TRUE)
  if(inherits(cholx, "try-error")){
    cat("Error with the Cholesky decomposition\n")
    return(rcppeigen_invert_matrix(x))
  }else{
    # cholfac() returns the LOWER factor L (L %*% t(L) = x), but
    # chol2inv() expects the UPPER factor U (t(U) %*% U = x, base R's
    # chol() convention) -- feeding it L instead silently produced a
    # wrong inverse. t(L) is the matching upper factor.
    return(chol2inv(t(cholx)))
  }
}