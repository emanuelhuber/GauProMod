
##--- COVARIANCE MATRIX
# x & y = points where to compute the covariance
# covModel = covariance model
# d = derivative? 0 = no derivative, 1 = first derivative, 2 = second deriv.
# dx = for 2D data: he gradient derivative, i.e., a unit vector normal to the 
# no-flow boundary (or tangent to a constant-head boundary)
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
    # M <- diag(rep(1L, ncol(x)))
    M <- diag(ncol(x))
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
    kernelName <- .kernelName(covModel$kernel)
    KK <- do.call(kernelName, list(x, y, covModel, d = d, w = 1, ...))
    return(KK)
  }else{
    XY <- crossDist(x, y, M)
    
    # DERIVATIVE WEIGHTS
    if(d == 1){
      if(is.null(dim(x)) && is.null(dim(y))){
        # 1D Case: w must be the derivative of the distance |x-y| w.r.t the second variable (y).
        # d/dy |x-y| = -sign(x-y). 
        # Your kGaussian uses w*r/l^2, so we need w = d/dr * d/dy(r) * ...
        # d/dy(r) = -sign(x-y).
        # For K_obs,dx (Kdxx), we need d/dy. The matrix is outer(x,y,"-").
        # The weight should be -sign(outer(x, y, "-")).
        w <- -sign(outer(x, y, "-"))
        # w <- (w0)
      }else if(length(dim(x)) > 1 && length(dim(y)) > 1){
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
      }else if(length(dim(x)) > 1 && length(dim(y)) > 1){
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


#' Return covariance as a function of distance
#'
#' @param r vector of distance
#' @param covModel Covariance mdoel
#' @name covfx
#' @export
#' @examples
#' covModel <- list(kernel="matern",
#'                  l = 5,     # correlation length
#'                  v = 1,     # smoothness
#'                  h = 2.45   # std. deviation
#' )
#' r <- seq(0, 20, by = 0.1)
#' myCov <- covfx(r = r, covModel = covModel)
#' plot(r, myCov, type = "l", ylim = c(0, max(myCov)),
#'      ylab = "covariance", xlab = "distance", xaxs = "i", yaxs = "i")
covfx <- function(r, covModel){
  kernelName <- .kernelName(covModel$kernel)
  do.call(kernelName, list(r, covModel))
}


sign2 <- function(x){
  ifelse(x == 0, 1, sign(x))
}


# dcovm <- function(x, y, covModel, ...){
#   if(length(covModel$type) == 1 && length(covModel$kernel) == 0){
#     covModel[["kernel"]] <- covModel$type
#     warning("In covModel, rename 'type' into 'kernel'.\n")
#   }
#   XY <- crossDist(x, y)
#   kernelName <- .kernelName(covModel$kernel)
#   KK <- do.call(kernelName, list(XY, covModel,  ...))
#   return(KK)
# }


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
#' @param X a matrix or vector
#' @param Y a matrix or vector with same number of columns as X
#' @param M a positive semidefinite matrix (nrow(M) = ncol(M) = ncol(X))
#' @name crossDist
#' @export
crossDist <- function(X, Y, M = NULL){
  # if(!identical(ncol(X), ncol(Y))){
  #   stop("X and Y must have identical dimensions!\n")
  # }
  # # case 1D
  # if(is.null(dim(X))){
  #   # return( outer(X, Y, "-") )
  #   return( outer(X, Y, function(X, Y){ sqrt((X - Y)^2)}))
  # }else{
  #   # case 2D
  #   if(!is.null(M)){
  #     L <- cholfac(M)
  #     X <- X %*% L
  #     Y <- Y %*% L
  #   }
  #   Xn <- rowSums(X^2)
  #   Yn <- rowSums(Y^2)
  #   D2 <- outer(Xn, Yn, "+") - 2 * tcrossprod(X, Y)
  #   D2[D2 < 0] <- 0  # numerical stability
  #   sqrt(D2)
  # }
  # Ensure X and Y are matrices for RcppEigen
  X_mat <- as.matrix(X)
  Y_mat <- as.matrix(Y)
  
  # Original check
  if(!identical(ncol(X_mat), ncol(Y_mat))){
    stop("X and Y must have identical number of columns!")
  }
  
  # Call the Rcpp version
  # The Rcpp function handles both 1D and 2D cases efficiently
  return(crossDist_rcpp(X_mat, Y_mat, M))
}
# crossDist <- function(X, Y, M = NULL){
#   if(!identical(ncol(X), ncol(Y))){
#     stop("X and Y must have identical dimensions!\n")
#   }
#   if(is.null(dim(X))){
#     # return( outer(X, Y, "-") )
#     return( outer(X, Y, function(X, Y){ sqrt((X - Y)^2)}))
#   }else if(dim(X)[2] == 2){
#     return( dist2(X, Y, M))
#   }else{
#     return( distn(X, Y, M) )
#   }
# }

# # distance for ncol(X) > 2
# distn <- function(X, Y, M){
#   if(!is.null(M)){
#     L <- cholfac(M)
#     X <- X %*% (L)
#     Y <- Y %*% (L)
#   }
#   return( apply(outer(X,t(Y),"-"),c(1,4),
#                 function(x)sqrt(sum(diag(x*x)))))
# }
#
# dist2 <- function(X,Y, M){
#   if(!is.null(M)){
#     L <- cholfac(M)
#     X <- X %*% (L)
#     Y <- Y %*% (L)
#   }
#   nx <- nrow(X)
#   ny <- nrow(Y)
#   matx1 <- matrix(rep(X[, 1], ny), nx ,  ny)
#   maty1 <- matrix(rep(Y[, 1], nx), nx ,  ny, byrow = TRUE)
#   matx2 <- matrix(rep(X[, 2], ny), nx ,  ny)
#   maty2 <- matrix(rep(Y[, 2], nx), nx ,  ny, byrow = TRUE)
#   D <- sqrt((matx1 - maty1)^2 + (matx2 - maty2)^2)
#   return(D)
# }

# dist2oldschool <- function(X,Y,M){
#   if(!is.null(dim(X)) && dim(X)[2] > 1 && !is.null(M)){
#     L <- cholfac(M)
#     X <- X %*% (L)
#     Y <- Y %*% (L)
#   }
#   nx <- nrow(X)
#   ny <- nrow(Y)
#   Dref <- matrix(nrow = nx, ncol=ny)
#   for(i in 1:nx){
#     for(j in 1:ny){
#       U <-  (X[i,,drop=FALSE] - Y[j,,drop=FALSE]) 
#       Dref[i,j] <- U %*% t(U)
#     }
#   }
#   return(Dref)
# }


#' #' @export
#' se <- function(x, y, covModel, d = 0, w = 1){
#'   warning("Deprecated function! Use 'kGaussian' instead!\n")
#'   kGaussian(r, covModel, d, w)
#' }

#' #' @export
#' kSe <- function(r, covModel, d = 0, w = 1){
#'   warning("Deprecated function! Use 'kGaussian' instead!\n")
#'   kGaussian(r = r, para = covModel, d = d, w = w)
#' }

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
  #     Pardo-Iguzquiza, E., and M. Chica-Olmo (2004),
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

#' #' @export
#' linear <- function(x, y, covModel, d = 0, w = 1){
#'   warning("Deprecated function! Use 'kLinear' instead!\n")
#'   kLinear(x = x, y = y, para = covModel, d = d, w = w)
#' }

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

#' #' @export
#' matern <- function(r, covModel, d = 0, w = 1){
#'   warning("Deprecated function! Use 'kMatern' instead!\n")
#'   kMatern(r = r, para = covModel, d = d, w = w)
#' }

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
