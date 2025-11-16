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
      y <- c(obs$y, rep(bc$y, nt))
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



#' Return the lower Cholesky factor
#'
#' Return the lower Cholesky factor L such that X = L t(L)
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
#' @name gpSim
#' @export
gpSim <- function(A, L = NULL, n = 1){
  if(is.null(L)){
    ystar <-  try(mvrnorm2(n, A$mean, A$cov),silent=TRUE)
    if(class(ystar) == "try-error"){
      cat("Error with Cholesky decomposition...\n")
      ystar <-  MASS::mvrnorm(n, A$mean, A$cov)
    }
  }else{
    p <- length(A$mean)
    std <- rnorm(p * n)
    realz <- crossprod(t(L),  std)
    ystar <-  A$mean + matrix(realz, nrow=p, ncol=n, byrow = TRUE)
  }
  return(cbind(A$xstar,ystar))
}



# format correctly the position matrix (xy) and
# the time (in fact time-step), as well as the
# observed value and the target position matrix (xystar)
#' Reshape target
#'
#' @name setPosTime
#' @export
setPosTime <-function(xy, tt, val, xystar){
	tsteps <- c(1, 1+cumsum(diff(tt[!is.na(tt)])))
	nxy <- nrow(xy)
	nt <- length(tsteps)
	if(nt != length(tt)){
	  stop("problem with time\n")
	}
	# observations
	obs <- list()
	obs$x <- xy[rep(seq_len(nxy), each=nt),]
	obs$t <- rep(tsteps,nxy)
	obs$y <- val
	# xstar (points at which we want to predict)
	xstar <- list()
	xstar$x <- xystar[rep(seq_len(nrow(xystar)),each=nt),]
	xstar$t <- rep(seq_len(nt),nrow(xystar))

	return(list("obs"=obs, "xstar"=xstar))
}

# Tilman M. Davies and David J. Bryant (2013)
# On Circulant Embedding for Gaussian Random
# Fields in R. Journal of Statistical Software 55(9).
#' Multi-variate Gaussian simulation
#'
#' A more robute alternative to the \code{mvrnorm} function.
#' @name mvrnorm2
#' @export
mvrnorm2 <- function(n, mu, Sigma){
	p <- length(mu)
# 	cholStatus <- try(SChol <- chol(Sigma),silent=TRUE)
	cholStatus <- try(SChol <- cholfac_rcpp(Sigma),silent=TRUE)
	cholError <- ifelse(class(cholStatus) == "try-error", TRUE, FALSE)
	if(cholError){
	  cat("error\n")
		SChol <- correctCovMat(Sigma)
	}
	std <- rnorm(p * n)
	realz <- crossprod(t(SChol),  std)
	return(mu + matrix(realz, nrow=p, ncol=n, byrow = TRUE))
}


# 	http://comisef.wikidot.com/tutorial:repairingcorrelation
# http://www.r-bloggers.com/fixing-non-positive-definite-correlation-
# matrices-using-r-2/
#  Rebonato and Jackel, “The most general methodology for creating a valid
# correlation matrix for risk management and option pricing purposes”,
# Journal of Risk, Vol 2, No 2, 2000
correctCovMat <- function(Sigma){
	iter <- 0
	cholError <- TRUE
	newSigma <- Sigma
	while (cholError) {
		cat("*")
		# compute eigenvectors/-values
		E   <- eigen(newSigma, symmetric = TRUE)
		# replace negative eigenvalues by zero
		E$values   <- pmax(E$values,0)
		# reconstruct correlation matrix
		newSigma  <- E$vectors %*% diag(E$values) %*% t(E$vectors)
		newSigma <- newSigma/sqrt(diag(newSigma) %*% t(diag(newSigma)))
		cholStatus <- try(u <- chol(newSigma), silent = TRUE)
		cholError <- ifelse(class(cholStatus) == "try-error", TRUE, FALSE)
	}
	return(u)
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
#' Create grid
#'
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

#' Create vecgrid
#'
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
#' Inverse matrix
#'
#' This function first try the Cholesky decomposition
#' @name invm
#' @export
invm <- function(x){
	cholx <- try(cholfac(x),silent=TRUE)
	if(class(cholx) == "try-error"){
	  cat("Error with the Cholesky decomposition\n")
	  return(rcppeigen_invert_matrix(x))
	}else{
	  return(chol2inv(cholx))
	}
}






