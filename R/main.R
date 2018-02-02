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
                    sigmat = 0){

  Kxx       <- covm( obs$x,  obs$x, covModels[[1]])
  Kstar     <- covm( obs$x, targ$x, covModels[[1]])
  Kstarstar <- covm(targ$x, targ$x, covModels[[1]])
  if(length(sigma) == 1){
    sigma     <- rep(sigma, ncol(Kxx))
  }else if(length(sigma) != ncol(Kxx)){
    stop("length of sigma must be equal to 1 or to the number of observations")
  }
  y         <- obs$y

  # if there are derivative
  if(!is.null(bc)){
    Kdxx  <- covm(obs$x, bc$x, covModels[[1]] , d = 1, dx = bc$v)
    Kdxdx <- covm( bc$x, bc$x, covModels[[1]] , d = 2, dx = bc$v)
    sigma <- c(sigma, rep(bc$sigma, ncol(Kdxx)))
    Kxx  <- rbind(cbind(   Kxx,  Kdxx ),
                  cbind(t(-Kdxx), Kdxdx))
    Kdxstar <- covm(targ$x, bc$x, covModels[[1]] , d = 1, dx=bc$v)
    Kstar   <- rbind(Kstar, t(Kdxstar))
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
    Ktt <- covm(obs$t, obs$t, covModels[[2]])
    Knoise <- diag(sigma^2)
    Ktnoise <- diag(rep(sigmat^2, length(obs$t)))
    Kxx <-  (Kxx + Knoise) %x% (Ktt + Ktnoise)
    # structure of Kxx: Kxx[1,1]*Knoise[1,1] Kxx[1,1]*Knoise[1,2] ...
    #                   Kxx[1,1]*Knoise[2,1] Kxx[1,1]*Knoise[2,2] ...
    #                   ...
    Kstarstar <- Kstarstar %x%
                 covm(targ$t, targ$t, covModels[[2]])
#     Kstar <- covm(obs$x, targ$x, covModels[[1]]) %x%
    Kstar <- Kstar %x% covm(obs$t, targ$t, covModels[[2]])
    if(!is.null(bc)){
#       nbc <- nrow(bc$x)
#       y <- numeric(length = nt*(nxy + nbc))
#       yo <- rep((1:nxy),nt) +
#               rep((nbc+nxy)*(0:(nt-1)),each=nxy)
#       ybc <- rep((nxy+1):(nbc+nxy),nt) +
#               rep((nbc+nxy)*(0:(nt-1)),each=nbc)
#       y[yo] <- obs$y
#       y[ybc] <- rep(bc$y, nt)

      y <- c(obs$y, rep(bc$y, nt))
    }
  }else{
    Kxx <-  Kxx  +  diag(sigma^2)
#     Kstarstar <- covm(targ$x, targ$x, covModels[[1]])
#     Kstar <- covm(obs$x, targ$x, covModels[[1]])
#     y <- obs$y
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
#         H0 <- Hmat(obs$x,op)
#         Hdx <- Hmat(bc$x, op, dx = bc$v)
#         H0 <- cbind(H0, Hdx)
#         H <- matrix(rep((H0), nt), nrow = nrow(H0), byrow = FALSE)
#         dim(H) = 3 x 30
#         dim(Hdx) = 3 x 5
        Hdx <- Hmat(bc$x[rep(seq_len(nbc), each=nt),], op, dx = bc$v)
        H <- cbind(H, Hdx)
      }else{
        Hdx <- Hmat(bc$x, op, dx = bc$v)
        H <- cbind(H, (Hdx))
      }
    }
    Hstar <- Hmat(xstar,op)
    A2 <- GPpredmean_rcpp(Kxx, Kstar, Kstarstar, y, H, Hstar)
    m <- qr(H)$rank
    logLik <- - A2$logLik1 - sum(log(A2$logLik2)) - sum(log(A2$logLik3)) -
                (nrow(Kxx)- m)/2 * log(2*pi)
    A2["logLik1"] <- NULL
    A2["logLik2"] <- NULL
    A2["logLik3"] <- NULL
    A2[["logLik"]] <- logLik
  }else{
    A2 <- GPpred_rcpp(Kxx, Kstar, Kstarstar, y)
    logLik <- - A2$logLik1 - sum(log(A2$logLik2)) - nrow(Kxx)/2 * log(2*pi)
    A2["logLik1"] <- NULL
    A2["logLik2"] <- NULL
    A2[["logLik"]] <- logLik
  }
  A2[["xstar"]] <- AA
  return(A2)
}


gpCondOld <- function(obs, targ, covModels, sigma=0, op = 0 , bc = NULL){

  if(length(covModels) == 2){
    if(is.null(targ$t)){
	    targ$t <- obs$t
	  }
	  nt <- length(obs$t)
	  xstar <- targ$x[rep(seq_len(nrow(targ$x)),each=nt),]
    nxy <- nrow(obs$x)
	  x <- obs$x[rep(seq_len(nxy), each=nt),]
	  y <- obs$y
	  AA <- cbind(xstar,targ$t)
	  Kxx <-  covm(obs$x, obs$x, covModels[[1]]) %x%
            covm(obs$t, obs$t, covModels[[2]]) +
            diag(sigma^2, length(obs$y))
	  Kstarstar <- covm(targ$x, targ$x, covModels[[1]]) %x%
                 covm(targ$t, targ$t, covModels[[2]])
	  Kstar <- covm(obs$x, targ$x, covModels[[1]]) %x%
             covm(obs$t, targ$t, covModels[[2]])
  }else{
	  Kxx <-  covm(obs$x, obs$x,covModels[[1]])  +  diag(sigma^2, length(obs$y))
	  Kstarstar <- covm(targ$x, targ$x,covModels[[1]])
	  Kstar <- covm(obs$x, targ$x,covModels[[1]])
	  x <- obs$x
	  y <- obs$y
	  xstar <- targ$x
	  AA <- targ$x
  }
  if(op > 0){
	  H <- Hmat(x,op)
	  Hstar <- Hmat(xstar,op)
	  A2 <- GPpredmean_rcpp(Kxx,Kstar,Kstarstar,y, H, Hstar)
  }else{
	  A2 <- GPpred_rcpp(Kxx, Kstar, Kstarstar, y)
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
	cholx <- try(chol(x),silent=TRUE)
	if(class(cholx) == "try-error"){
	  cat("Error with the Cholesky decomposition\n")
	  return(rcppeigen_invert_matrix(x))
	}else{
	  return(chol2inv(cholx))
	}
}




##--- COVARIANCE MATRIX
# x & y = points where to compute the covariance
# covModel = covariance model
# d = derivative? 0 = no derivative, 1 = first derivative, 2 = second deriv.
# dx = ?
#' Covariance matrix
#'
#' Create a covariance matrix according to the kernel parametrisation
#' @name covm
#' @export
covm <- function(x, y, covModel, d = 0, dx = 1, ...){
#   outer(x,y, covModel$kernel,covModel)
  if(length(covModel$type) == 1 && length(covModel$kernel) == 0){
    covModel[["kernel"]] <- covModel$type
    warning("In covModel, rename 'type' into 'kernel'.\n")
  }
  if(is.null(dim(x))){
      #XY <- outer(x, y, function(x, y){ sqrt((x - y)^2)})
    M <- NULL
  }else{
    # M = identity matrix
    M <- diag(rep(1L, ncol(x)))
    # scaling
    if(!is.null(covModel$scale)){
      if(length(covModel$scale) == ncol(x)){
        M <- M %*% diag(covModel$scale)
      }else{
        stop(paste0("'covModel$scale' must have length identical ",
                    "to the number of position coordinates!\n"))
      }
    }
    #CHECK scaling/rotation where when?
    # rotation
    if(!is.null(covModel$rot)){
      if(ncol(x) == 2){
        if(length(covModel$rot) == 1){
          mrot <- covModel$rot
          M <- M %*% matrix(c(cos(mrot), - sin(mrot),
                            sin(mrot), cos(mrot)),
                            ncol = 2, nrow = 2, byrow = TRUE)
        }else{
          stop(paste0("'covModel$rot' must have length one!\n"))
        }
      }else if(ncol(x) == 3){
        if(length(covModel$rot) == 2){
          a1 <- covModel$rot[1]
          a2 <- covModel$rot[2]
          M <- M %*% matrix(c(cos(a1)*cos(a2), -sin(a1), -cos(a1)*sin(a2),
                              sin(a1)*cos(a2),  cos(a1), -sin(a1)*sin(a2),
                              sin(a2), 0, cos(a2)),
                            ncol = 2, nrow = 2, byrow = TRUE)
        }else{
          stop(paste0("'covModel$rot' must have length two!\n"))
        }
      }else{
        stop(paste0("'covModel$rot' must have length one!\n"))
      }
    }
  }
  if(covModel$kernel == "linear"){
    if(!is.null(dim(x)) && dim(x)[2] > 1 && !is.null(M)){
      L <- cholfac(M)
      x <- x %*% (L)
      y <- y %*% (L)
    }
    # kernelName <- paste0("k", toupper(substr(covModel$kernel,0,1)),
    #                      substr(covModel$kernel,2,nchar(covModel$kernel)) )
    kernelName <- .kernelName(covModel$kernel)
    KK <- do.call(kernelName, list(x, y, covModel, d = d, w = 1, ...))
    return(KK)
  }else{
    XY <- crossDist(x, y, M)
    #if(dim(x)[2] == 2){
    #    XY <- dist2(x, y, M)
    #}else{
    #  XY <- distn(x,y, M)
    #}
    if(d == 1){
      if(is.null(dim(x)) && is.null(dim(y))){
        w0 <- sign(outer(x, y, "-"))
  #         w1 <- outer(rep(1,length(x)), dx, function(x,y){sign2(y)})
  #         w <- (w0*w1)
        w <- (w0)
      }else if(dim(x) > 1 && dim(y) > 1){
        r1 <- outer(x[,1], y[,1], "-")
        r2 <- outer(x[,2], y[,2], "-")
        rn <- sqrt(r1^2 + r2^2)
        v1 <-  outer(x[,1], dx[,1], function(x,y) y)
        v2 <-  outer(x[,1], dx[,2], function(x,y) y)
        vn <- sqrt(v1^2 + v2^2)
        w <- (r1*v1 + r2*v2)/(rn*vn)
        w[(rn*vn) == 0] <- 0
      }
    }else if(d == 2){
      if(is.null(dim(x)) && is.null(dim(y))){
        w <- 1
      }else if(dim(x) > 1 && dim(y) > 1){
        r1 <- outer(x[,1], y[,1], "-")
        r2 <- outer(x[,2], y[,2], "-")
        rn <- sqrt(r1^2 + r2^2)
        v1 <-  outer(dx[,1], dx[,1], function(x,y) x)
        v2 <-  outer(dx[,2], dx[,2], function(x,y) x)
        vn <- sqrt(v1^2 + v2^2)
        u1 <-  outer(dx[,1], dx[,1], function(x,y) y)
        u2 <-  outer(dx[,2], dx[,2], function(x,y) y)
        un <- sqrt(u1^2 + u2^2)
        w1 <- (r1*v1 + r2*v2)/(rn*vn)
        w1[(rn*vn) == 0] <- 0
        w2 <- (r1*u1 + r2*u2)/(rn*un)
        w2[(rn*un) == 0] <- 0
        w <- w1 * w2
      }
    }else{
      w = 1
    }
    # kernelName <- paste0("k", toupper(substr(covModel$kernel,0,1)),
    #                      substr(covModel$kernel,2,nchar(covModel$kernel)) )
    kernelName <- .kernelName(covModel$kernel)
    KK <- do.call(kernelName, list(XY, covModel, d = d, w = w, ...))
    return(KK)
  }
}



sign2 <- function(x){
  ifelse(x == 0, 1, sign(x))
}


dcovm <- function(x, y, covModel, ...){
  if(length(covModel$type) == 1 && length(covModel$kernel) == 0){
    covModel[["kernel"]] <- covModel$type
    warning("In covModel, rename 'type' into 'kernel'.\n")
  }
  XY <- crossDist(x, y)
	# if(is.null(dim(x))){
	#	XY <- outer(x, y, "-")
	# }else if(dim(x)[2] == 2){
	#	XY <- dist2(x,y)
	# }else{
	#	XY <- distn(x,y)
	# }
  # kernelName <- paste0("k", toupper(substr(covModel$kernel,0,1)),
  #                     substr(covModel$kernel,2,nchar(covModel$kernel)) )
  kernelName <- .kernelName(covModel$kernel)
  KK <- do.call(kernelName, list(XY, covModel,  ...))
  return(KK)
}


.kernelName <- function(kname){
  return(paste0("k", toupper(substr(kname,0,1)),
                substr(kname,2,nchar(kname)) ))
}

#' Cross-distance between two matrix
#' 
#' Compute the distance between every rows of two matrix. The returned distance
#' has for dimension: nrow(X) x ncol(Y).
#' If M is the identity matrix (by default), the distance is isotropic, if not
#' the distance is anisotropic.
#' @param X a matrix
#' @param Y a matrix with same number of columns as X
#' @param M a positive semidefinite matrix (nrow(M) = ncol(M) = ncol(X))
#' @name covm
#' @export
crossDist <- function(X, Y, M = NULL){
  if(!identical(ncol(X), ncol(Y))){
    stop("X and Y must have identical dimensions!\n")
  }
  if(is.null(dim(X))){
    # return( outer(X, Y, "-") )
    return( outer(X, Y, function(X, Y){ sqrt((X - Y)^2)}))
  }else if(dim(X)[2] == 2){
    return( dist2(X, Y, M)) 
  }else{
    return( distn(X, Y, M) )
  }
}

# distance for ncol(X) > 2
distn <- function(X, Y, M){
  if(!is.null(M)){
    L <- cholfac(M)
    X <- X %*% (L)
    Y <- Y %*% (L)  
  }
	return( apply(outer(X,t(Y),"-"),c(1,4),
	           function(x)sqrt(sum(diag(x*x)))))
}

dist2 <- function(X,Y, M){
  if(!is.null(M)){
    L <- cholfac(M)
    X <- X %*% (L)
    Y <- Y %*% (L)
  }
	nx <- nrow(X)
	ny <- nrow(Y)
	matx1 <- matrix(rep(X[, 1], ny), nx ,  ny)
	maty1 <- matrix(rep(Y[, 1], nx), nx ,  ny, byrow = TRUE)
	matx2 <- matrix(rep(X[, 2], ny), nx ,  ny)
	maty2 <- matrix(rep(Y[, 2], nx), nx ,  ny, byrow = TRUE)
	D <- sqrt((matx1 - maty1)^2 + (matx2 - maty2)^2)
	return(D)
}

dist2oldschool <- function(X,Y,M){
  if(!is.null(dim(X)) && dim(X)[2] > 1 && !is.null(M)){
    L <- cholfac(M)
    X <- X %*% (L)
    Y <- Y %*% (L)
  }
  nx <- nrow(X)
  ny <- nrow(Y)
  Dref <- matrix(nrow = nx, ncol=ny)
  for(i in 1:nx){
    for(j in 1:ny){
      U <-  (X[i,,drop=FALSE] - Y[j,,drop=FALSE]) 
      Dref[i,j] <- U %*% t(U)
    }
  }
  return(Dref)
}


#' @export
se <- function(x, y, covModel, d = 0, w = 1){
  warning("Deprecated function! Use 'kGaussian' instead!\n")
  kGaussian(r, covModel, d, w)
}

#' @export
kSe <- function(r, covModel, d = 0, w = 1){
  warning("Deprecated function! Use 'kGaussian' instead!\n")
  kGaussian(r = r, para = covModel, d = d, w = w)
}

#' Kernels (covariance functions) for Gaussian process
#'
#' Squared Exponential Covariance Function (or radial basis or Gaussian)
#' over-smoothness, infinitely differentiable at h=0
#' @name kernels
#' @rdname kernels
#' @export
kGaussian <- function(r, para, d = 0, w = 1){
	l <- para$l
	h <- para$h
	u <- -0.5 * (r / l)^2
    K <-  exp(u)
    K[abs(u) < .Machine$double.eps^0.5] <- 1
    # same results as in Solak et al.,
    # Derivative observations in Gaussian Process Models of Dynamic Systems
#     if(d == 1){
#       K <-   (r/l) * K
#     }else if(d == 2){
#       K <- (1 - r^2/l)/l * K
#     }
#     Pardo-Igúzquiza, E., and M. Chica-Olmo (2004),
#     Estimation of gradients from sparse data by universal kriging,
#     Water Resour. Res., 40, W12418, doi:10.1029/2004WR003081.
    if(d == 1){
      K <-   (w*r/l^2) * K
#       K <-   (r/l^2) * K
    }else if(d == 2){
      K <- (1 - w*(r/l)^2)*(1/l^2) * K
#       K <- (1 -(r/l)^2)*(1/l^2) * K
    }
    K <- K*h^2
	return(K)
}

#' @export
linear <- function(x, y, covModel, d = 0, w = 1){
  warning("Deprecated function! Use 'kLinear' instead!\n")
  kLinear(x = x, y = y, para = covModel, d = d, w = w)
}

#' @name kLinear
#' @rdname kernels
#' @export
kLinear <- function(x, y, para, d = 0, w = 1){
  b <- para$b
  h <- para$h
  cc <- para$c
  x <- x - cc
  y <- y - cc
  if(is.null(dim(x))){
    XY <- outer(x, y, "*")
  }else{
    XY <- x %*% t(y)
  }
  K <- XY * h^2 + b^2
  return(K)
}


# Generalized Cauchy model
# gCauchy <- function(){
#   C = ss.*(1 + ha.^2).^(-m.alpha);
# }

#' @export
matern <- function(r, covModel, d = 0, w = 1){
  warning("Deprecated function! Use 'kMatern' instead!\n")
  kMatern(r = r, para = covModel, d = d, w = w)
}

# MATERN
# adjustable smoothness via parameter v
# second-order smooth fiel: v >= 2
# C. E. Rasmussen & C. K. I. Williams, Gaussian Processes for Machine Learning,
# the MIT Press, 2006, ISBN 026218253X. c 2006
# Massachusetts Institute of Technology. www.GaussianProcess.org/gpml
# Matern covariance matrix


#' @name kMatern
#' @rdname kernels
#' @export
kMatern <- function(r, para, d = 0, w = 1){
  l <- para$l
  v <- para$v
  h <- para$h
  u <-  r /l
  if(d == 0){
    if( v == 1/2 ){
      K <- exp(- r/l )
    }else if( v == 3/2){
      u <- sqrt(3) * r / l
      K <- (1 + u) * exp(- u)
    }else if(v==5/2){
      u <- sqrt(5) * r / l
      K <- (1 + u + (5 * r^2) / (3 * l^2)) * exp(- u)
    }else{
      gamV <- gamma(v)
      besK <- besselK(u, v)
      K <- ( 2^(1-v) / gamV )* u^v  * besK
    }
    K[ abs(u) < .Machine$double.eps ^ 0.5 ] <- 1
    K <- K*h^2
  }else if(d == 1){
    common <- h^2 * ( 2^(1-v) / gamma(v) ) * u^v
    K <- w * (1/l) *  common * besselK(u, v -1)
    K[u < .Machine$double.eps ^ 0.5 ] <- 0
#     common = ss.*2^(1 - m.nu)./gamma(m.nu).*hag.^m.nu;
#     hpg = hproj(lag > eps);
#     C(lag > eps) = hpg.*common./a.*besselk(m.nu-1,hag);
#     C(lag <= eps) = 0.0;  % specify manually
  }else if(d == 2){
#     cc = ss.*2^(1 - m.nu)./gamma(m.nu).*hag.^(m.nu - 2);
    k1 <- h^2 * ( 2^(1-v) / gamma(v) ) * u^(v - 2)
    k2 <- ( u * besselK(u, v -1)  -  u^2 * w * besselK(u, v -2))
    K <-  (1/l)^2 *k1 * k2
    k0 <- h^2 * gamma(v)/(2*l^2 * (v-1)*gamma(v))
    K[u < .Machine$double.eps ^ 0.5 ] <-  k0
  }
  return(K)
#   ifelse(u>0,( 2^(1-v) / gamV )* u^v * gamV * besK,1)
}



