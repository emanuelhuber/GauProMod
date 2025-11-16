

#' Covariance matrix
#'
#' Compute a covariance (kernel) matrix — or its derivatives — between two sets of
#' locations using a specified kernel model.  
#' The wrapper handles **anisotropy**, **rotation**, and **derivative weights**,
#' and calls a C++ backend for efficient evaluation.
#'
#' @param x Numeric vector or matrix of locations.  
#'   - If a vector, it is treated as 1D positions.  
#'   - If a matrix, it must have one location per row and the same number of columns as `y`.
#'
#' @param y Numeric vector or matrix of locations.  
#'   If omitted, defaults to `y = x`.
#'
#' @param covModel A named list describing the covariance model.  
#'   **Required fields:**
#'   - `kernel`: character, one of `"gaussian"`, `"exponential"`, `"matern_3_2"`,
#'     `"matern_5_2"`, `"power_exp"`, `"linear"`.  
#'     (A legacy `covModel$type` is also accepted and will be renamed.)
#'   - `l`: numeric, the length scale (positive).  
#'   - `h`: numeric, the marginal standard deviation (≥ 0).  
#'
#'   **Optional fields:**
#'   - `v`: smoothness exponent (for `"power_exp"` or `"matern"`).  
#'   - `b`: bias for the `"linear"` kernel.  
#'   - `scale`: numeric vector of length `ncol(x)` for per-dimension scaling.  
#'   - `rot`: rotation angle in radians (2D only).  
#'
#' @param d Integer derivative order:  
#'   - `0`: covariance (default)  
#'   - `1`: first radial derivative ∂K/∂r  
#'   - `2`: second radial derivative ∂²K/∂r²  
#'   
#'   When `d > 0`, the function multiplies the radial derivative by a
#'   weight matrix `w` so that the result corresponds to  
#'   ∂K/∂y = (dK/dr) × (dr/dy).
#'
#' @param dx For directional derivatives (2D only), either a length-2 numeric
#'   unit vector or an `nrow(x) × 2` matrix of unit direction vectors.  
#'   Ignored when `d = 0`.  
#'   For 1D data, `dx` is not used — the wrapper computes the correct sign-based
#'   weight automatically.
#'
#' @param sparse Logical; if `TRUE`, return a sparse matrix (`Matrix::dgCMatrix`)
#'   using a distance cutoff. Default is `FALSE`.
#'
#' @param cutoff Numeric; distance threshold used when `sparse = TRUE`.  
#'   Entries corresponding to distances greater than `cutoff` are set to zero.  
#'   When `cutoff <= 0`, sparsity is not applied.
#'
#' @param ... Additional arguments passed to kernel functions.
#'
#' @details
#' The function constructs an anisotropy matrix `M` (from `covModel$scale` and `covModel$rot`)
#' and computes pairwise distances:
#'
#' - For non-linear stationary kernels (e.g., `"gaussian"`, `"exponential"`, `"matern_3_2"`,
#'   `"matern_5_2"`, `"power_exp"`), it computes  
#'   *r* = ‖*x₁* − *y₁*‖ (possibly anisotropic) and passes the distance matrix to C++.
#'
#' - For the `"linear"` kernel, the dot-product matrix `x %*% t(y)` is used instead of distances.
#'
#' Derivative orders `d = 1` or `d = 2` compute radial derivatives multiplied by
#' geometric weights `w` so that the output corresponds to directional partial
#' derivatives with respect to `y`.
#'
#' Parameter ordering in the backend:  
#' `params = c(l, h, ...)`, with additional entries if needed:
#'
#' - For `"power_exp"`: append `v`.  
#' - For `"linear"`: append `b` (default 0 if missing).
#'
#' @return
#' - A dense numeric matrix (`matrix`) when `sparse = FALSE`.  
#' - A sparse matrix (`Matrix::dgCMatrix`) when `sparse = TRUE`.  
#'
#' If `d > 0`, the entries represent directional derivatives of the covariance.
#'
#' @seealso
#' - [crossDist()] for distance computation  
#' - [covfx()] for covariance as a function of distance  
#' - The underlying C++ wrapper `.covm_rcpp_wrapper`
#'
#' @examples
#' # --- 1D Gaussian covariance
#' x <- seq(0, 5, length.out = 6)
#' covModel <- list(kernel = "gaussian", l = 1, h = 2)
#' K <- covm(x, x, covModel)
#' dim(K)
#'
#' # --- 2D Matern (3/2)
#' pts <- expand.grid(seq(0, 1, by = 0.5), seq(0, 1, by = 0.5))
#' covModel2 <- list(kernel = "matern_3_2", l = 0.5, h = 1)
#' K2 <- covm(as.matrix(pts), as.matrix(pts), covModel2)
#'
#' # --- Linear kernel (dot product)
#' covLin <- list(kernel = "linear", l = 1, h = 1.5, b = 0.1)
#' Klin <- covm(pts, pts, covLin)
#'
#' # --- First derivative in 1D
#' Kd1 <- covm(x, x, covModel, d = 1)
#'
#' # --- Sparse computation with cutoff
#' Ksparse <- covm(x, x, covModel2, sparse = TRUE, cutoff = 0.3)
#'
#' @export
covm <- function(x, y, covModel, d = 0, dx = 1, use_symmetry = FALSE, ...){
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
    KK <- do.call(kernelName, list(x, y, covModel, d = d, w = 1, use_symmetry = use_symmetry, ...))
    return(KK)
  }else{
    XY <- crossDist(x, y, M, use_symmetry = use_symmetry)

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
    kernelName <- .kernelName(covModel$kernel)
    KK <- do.call(kernelName, list(XY, covModel, d = d, w = w, use_symmetry = use_symmetry, ...))
    return(KK)
  }
}

    # # Compute derivative weights w only if d > 0
    # w_arg <- NULL
    # if (d > 0) {
    #   if (ncol(x_mat) == 1) {
    #     # 1D: outer differences sign
    #     w_arg <- -sign(outer(as.vector(x_mat), as.vector(y_mat), "-"))
    #   } else {
    #     # 2D case: follow your original computation
    #     if (d == 1) {
    #       r1 <- outer(x_mat[,1], y_mat[,1], "-")
    #       r2 <- outer(x_mat[,2], y_mat[,2], "-")
    #       rn <- sqrt(r1^2 + r2^2)
    #       # dx can be single unit vector or matrix of vectors per row
    #       if (is.null(dim(dx))) {
    #         v1 <- outer(rep(1, nrow(x_mat)), rep(dx[1], nrow(y_mat)), function(a,b) b) # simple
    #         v2 <- outer(rep(1, nrow(x_mat)), rep(dx[2], nrow(y_mat)), function(a,b) b)
    #       } else {
    #         # dx provided per location; ensure dims match
    #         dx_mat <- as.matrix(dx)
    #         if (nrow(dx_mat) != nrow(x_mat)) stop("dx must have same number of rows as x when provided per point.")
    #         v1 <- outer(x_mat[,1], dx_mat[,1], function(a,b) b)
    #         v2 <- outer(x_mat[,2], dx_mat[,2], function(a,b) b)
    #       }
    #       vn <- sqrt(v1^2 + v2^2)
    #       w_tmp <- (r1 * v1 + r2 * v2) / (rn * vn)
    #       w_tmp[(rn * vn) == 0] <- 0
    #       w_arg <- w_tmp
    #     } else if (d == 2) {
    #       # your original code: product of two directional cosines
    #       r1 <- outer(x_mat[,1], y_mat[,1], "-")
    #       r2 <- outer(x_mat[,2], y_mat[,2], "-")
    #       rn <- sqrt(r1^2 + r2^2)
    #       # we interpret dx as direction vectors; if single vector, broadcast
    #       if (is.null(dim(dx))) {
    #         v1 <- outer(rep(1, nrow(x_mat)), rep(dx[1], nrow(y_mat)), function(a,b) a) # note swap but consistent with original
    #         v2 <- outer(rep(1, nrow(x_mat)), rep(dx[2], nrow(y_mat)), function(a,b) a)
    #         u1 <- v1; u2 <- v2
    #       } else {
    #         dx_mat <- as.matrix(dx)
    #         if (nrow(dx_mat) != nrow(x_mat)) stop("dx must have same number of rows as x when provided per point.")
    #         v1 <- outer(dx_mat[,1], dx_mat[,1], function(a,b) a) # matches your original style
    #         v2 <- outer(dx_mat[,2], dx_mat[,2], function(a,b) a)
    #         u1 <- outer(dx_mat[,1], dx_mat[,1], function(a,b) b)
    #         u2 <- outer(dx_mat[,2], dx_mat[,2], function(a,b) b)
    #       }
    #       vn <- sqrt(v1^2 + v2^2)
    #       un <- sqrt(u1^2 + u2^2)
    #       w1 <- (r1 * v1 + r2 * v2) / (rn * vn); w1[(rn * vn) == 0] <- 0
    #       w2 <- (r1 * u1 + r2 * u2) / (rn * un); w2[(rn * un) == 0] <- 0
    #       w_arg <- w1 * w2
    #     } else {
    #       w_arg <- matrix(1, nrow = nrow(x_mat), ncol = nrow(y_mat)) # fallback
    #     }
    #   }
    # }
    

# Return covariance as a function of distance
#'
# @param r vector of distance
# @param covModel Covariance mdoel
# @name covfx
# @export
# @examples
# covModel <- list(kernel="matern",
#                  l = 5,     # correlation length
#                  v = 1,     # smoothness
#                  h = 2.45   # std. deviation
# )
# r <- seq(0, 20, by = 0.1)
# myCov <- covfx(r = r, covModel = covModel)
# plot(r, myCov, type = "l", ylim = c(0, max(myCov)),
#      ylab = "covariance", xlab = "distance", xaxs = "i", yaxs = "i")
covfx <- function(r, covModel){
  kernelName <- .kernelName(covModel$kernel)
  do.call(kernelName, list(r, covModel))
}




sign2 <- function(x){
  ifelse(x == 0, 1, sign(x))
}


.kernelName <- function(kname){
  return(paste0("k", toupper(substr(kname,0,1)),
                substr(kname,2,nchar(kname)) ))
}

# Generic R wrapper for all kernels
kernel_wrapper <- function(X, Y, para, d = 0, w = 1, kernel_type, use_symmetry = FALSE){
  # Extract common parameters, defaulting to 0 if missing
  l      <- if(!is.null(para$l)) para$l else 0
  h      <- if(!is.null(para$h)) para$h else 0
  v      <- if(!is.null(para$v)) para$v else 0
  degree <- if(!is.null(para$degree)) para$degree else 0
  c      <- if(!is.null(para$c)) para$c else 0
  
  # Expand scalar w to full weight matrix
  W_mat <- if(is.matrix(w)) w else matrix(w, nrow=nrow(X), ncol=nrow(Y))
  
  # Shift X and Y if needed (for linear or polynomial)
  if(kernel_type %in% c("kLinear", "kPolynomial")){
    X <- X - c
    Y <- Y - c
  }
  
  # Call the corresponding Rcpp kernel function
  K <- switch(kernel_type,
              kGaussian   = kGaussian_rcpp(X, Y, l, h, v, degree, c, d, W_mat, use_symmetry),
              kLinear     = kLinear_rcpp(X, Y, l, h, v, degree, c, d, W_mat, use_symmetry),
              kPolynomial = kPolynomial_rcpp(X, Y, l, h, v, degree, c, d, W_mat, use_symmetry),
              kMatern     = kMatern_rcpp(X, Y, l, h, v, degree, c, d, W_mat, use_symmetry),
              kCauchy     = kCauchy_rcpp(X, Y, l, h, v, degree, c, d, W_mat, use_symmetry),
              kTriangular = kTriangular_rcpp(X, Y, l, h, v, degree, c, d, W_mat, use_symmetry),
              kSpherical  = kSpherical_rcpp(X, Y, l, h, v, degree, c, d, W_mat, use_symmetry),
              stop("Unknown kernel type"))
  
  return(K)
}

# General helper to ensure W is a matrix of correct dimensions
make_W <- function(X, w) {
  if (is.matrix(w)) {
    return(w)
  } else {
    return(matrix(w, nrow = nrow(X), ncol = ncol(X)))
  }
}

# ---------------- Distance-based kernels ----------------
kGaussian <- function(r, para, d = 0, w = 1, use_symmetry = FALSE){
  l <- para$l
  h <- para$h
  v <- 0       # not used for Gaussian
  degree <- 0
  c <- 0
  W_mat <- make_W(r, w)
  K <- kernel_dispatch_auto_rcpp(r, r, l, h, v, degree, c, d, W_mat, "gaussian", use_symmetry)
  return(K)
}

kMatern <- function(r, para, d = 0, w = 1, use_symmetry = FALSE){
  l <- para$l
  h <- para$h
  v <- para$v
  degree <- 0
  c <- 0
  W_mat <- make_W(r, w)
  K <- kernel_dispatch_auto_rcpp(r, r, l, h, v, degree, c, d, W_mat, "matern", use_symmetry)
  return(K)
}

kCauchy <- function(r, para, d = 0, w = 1, use_symmetry = FALSE){
  l <- para$l
  h <- para$h
  v <- para$v
  degree <- 0
  c <- 0
  W_mat <- make_W(r, w)
  K <- kernel_dispatch_auto_rcpp(r, r, l, h, v, degree, c, d, W_mat, "cauchy", use_symmetry)
  return(K)
}

kTriangular <- function(r, para, d = 0, w = 1, use_symmetry = FALSE){
  l <- para$l
  h <- para$h
  v <- 0
  degree <- 0
  c <- 0
  W_mat <- make_W(r, w)
  K <- kernel_dispatch_auto_rcpp(r, r, l, h, v, degree, c, d, W_mat, "triangular", use_symmetry)
  return(K)
}

kSpherical <- function(r, para, d = 0, w = 1, use_symmetry = FALSE){
  l <- para$l
  h <- para$h
  v <- 0
  degree <- 0
  c <- 0
  W_mat <- make_W(r, w)
  K <- kernel_dispatch_auto_rcpp(r, r, l, h, v, degree, c, d, W_mat, "spherical", use_symmetry)
  return(K)
}

# ---------------- Feature/Gram kernels ----------------
kLinear <- function(X, Y, para, d = 0, w = 1, use_symmetry = FALSE){
  b <- para$b    # scale factor
  h <- para$h
  c <- para$c    # bias
  v <- 0
  degree <- 1
  W_mat <- make_W(X,w)
  # Subtract bias from features if desired
  K <- kernel_dispatch_auto_rcpp(X - c, Y - c, b, h, v, degree, c, d, W_mat, "linear", use_symmetry)
  return(K)
}

kPolynomial <- function(X, para, d = 0, w = 1, use_symmetry = FALSE){
  degree <- para$degree
  h <- para$h
  c <- para$c
  l <- 1    # length scale unused for polynomial
  v <- 0
  W_mat <- make_W(X,  w)
  K <- kernel_dispatch_auto_rcpp(X - c, Y - c, l, h, v, degree, c, d, W_mat, "polynomial", use_symmetry)
  return(K)
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
crossDist <- function(X, Y, M = NULL, use_symmetry = use_symmetry){
  # Ensure X and Y are matrices for RcppEigen
  # X_mat <- as.matrix(X)
  # Y_mat <- as.matrix(Y)
  
  # Original check
  # if(!identical(ncol(X_mat), ncol(Y_mat))){
  #   stop("X and Y must have identical number of columns!")
  # }
  # Call the Rcpp version
  # The Rcpp function handles both 1D and 2D cases efficiently
  return(crossDist_rcpp(as.matrix(X), as.matrix(Y), M, use_symmetry = use_symmetry))
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
