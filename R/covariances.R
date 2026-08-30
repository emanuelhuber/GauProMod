

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
#'   - `scale`: numeric vector of length `ncol(x)` for per-dimension scaling.  
#'   - `rot`: rotation angle in radians (2D only).  
#'
#' @param d Integer derivative order:  
#'   - `0`: covariance (default)  
#'   - `1`: first derivative with respect to `x`, ∂K/∂x  
#'   - `2`: the negative of the second radial derivative, -∂²K/∂r²  
#'         (this is what is needed for the diagonal/self-covariance of a
#'         derivative-observation process, e.g. via `bc` in
#'         [gpCond()]; it does not depend on `w`/`dx`, unlike `d = 1`)
#'   
#'   For `d = 1`, the result is `w * (dK/dr) * (dr/dx)` for a weight
#'   matrix `w` built internally from `x`, `y` (and, for 2D+, `dx`) --
#'   see the `dx` argument below. Support for `d`/derivatives varies by
#'   kernel: not all kernels implement `d = 1`/`d = 2` (e.g. `"linear"`/
#'   `"polynomial"` only support `d = 0`), and `sparse = TRUE` only
#'   supports `d = 0`.
#'
#' @param dx For directional derivatives (2D only), either a length-2 numeric
#'   unit vector or an `nrow(x) × 2` matrix of unit direction vectors.  
#'   Ignored when `d = 0`.  
#'   For 1D data, `dx` is not used — the wrapper computes the correct sign-based
#'   weight automatically.
#'
#' @param sparse Logical; if `TRUE`, return a sparse matrix (`Matrix::dgCMatrix`)
#'   using a distance cutoff. Only supported for distance-based kernels
#'   (`"gaussian"`, `"matern"` and its `"matern_3_2"`/`"matern_5_2"`/
#'   `"exponential"` aliases, `"cauchy"`, `"triangular"`, `"spherical"`)
#'   and `d = 0`; not supported for `"linear"`/`"polynomial"` or for
#'   derivatives (`d = 1` or `d = 2`). Default is `FALSE`.
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
#' - For the `"linear"`/`"polynomial"` (Gram/feature) kernels, the
#'   dot-product matrix `x %*% t(y)` is used instead of distances.
#'
#' Derivative orders `d = 1` or `d = 2` compute derivatives multiplied by
#' geometric weights `w` (for `d = 1`) so that the output corresponds to
#' directional partial derivatives with respect to `x` -- see the `d`
#' argument above for the exact convention.
#'
#' Parameter ordering in the backend:  
#' `params = c(l, h, ...)`, with additional entries if needed:
#'
#' - For `"power_exp"`: append `v`.  
#'
#' @return
#' - A dense numeric matrix (`matrix`) when `sparse = FALSE`.  
#' - A sparse matrix (`Matrix::dgCMatrix`) when `sparse = TRUE`.  
#'
#' If `d > 0`, the entries represent directional derivatives of the covariance.
#'
#' @seealso
#' - [crossDist()] for distance computation  
#' - [gpCond()], which uses `covm()` internally to build the observation,
#'   cross-, and target covariance matrices  
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
#' covLin <- list(kernel = "linear", h = 1.5, c = 0.1)
#' Klin <- covm(pts, pts, covLin)
#'
#' # --- First derivative in 1D
#' Kd1 <- covm(x, x, covModel, d = 1)
#'
#' # --- Sparse computation with cutoff
#' Ksparse <- covm(x, x, covModel2, sparse = TRUE, cutoff = 0.3)
#'
#' @export
covm <- function(x, y, covModel, d = 0, dx = 1, use_symmetry = FALSE,
                  sparse = FALSE, cutoff = 0, ...){
  #   outer(x,y, covModel$kernel,covModel)
  if(length(covModel$type) == 1 && length(covModel$kernel) == 0){
    covModel[["kernel"]] <- covModel$type
    warning("In covModel, rename 'type' into 'kernel'.\n")
  }
  covModel <- .resolveKernelAlias(covModel)

  # Coerce matrix-like inputs (e.g. a data.frame from expand.grid()) to
  # real numeric matrices, so downstream matrix algebra (e.g. x %*% L
  # for scaling/rotation) doesn't fail on non-matrix input. Plain
  # vectors (dim(x) == NULL) are left untouched.
  if(!is.null(dim(x))) x <- as.matrix(x)
  if(!is.null(dim(y))) y <- as.matrix(y)

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

  # ---- Sparse path: distance-based kernels only, d = 0 only ----
  if(isTRUE(sparse)){
    if(covModel$kernel %in% c("linear", "polynomial")){
      stop("'sparse = TRUE' is only supported for distance-based kernels ",
           "(e.g. \"gaussian\", \"matern\", \"cauchy\", \"triangular\", ",
           "\"spherical\"), not for the feature/Gram kernels ",
           "\"linear\"/\"polynomial\".")
    }
    if(d != 0){
      stop("'sparse = TRUE' is currently only supported for d = 0 ",
           "(the covariance itself, not its derivatives).")
    }
    if(!requireNamespace("Matrix", quietly = TRUE)){
      stop("'sparse = TRUE' requires the 'Matrix' package.")
    }
    x_mat <- if(is.null(dim(x))) matrix(x, ncol = 1) else as.matrix(x)
    y_mat <- if(is.null(dim(y))) matrix(y, ncol = 1) else as.matrix(y)
    rmax  <- if(is.null(cutoff) || cutoff <= 0) Inf else cutoff
    Xsp   <- crossDist_sparse(x_mat, y_mat, rmax = rmax, M = M)
    Wsp   <- Xsp
    Wsp@x[] <- 1   # all-ones weights, same sparsity pattern as Xsp

    l <- if(!is.null(covModel$l)) covModel$l else 1
    h <- if(!is.null(covModel$h)) covModel$h else 1
    v <- if(!is.null(covModel$v)) covModel$v else 0
    KK <- kernel_dispatch_auto_rcpp(Xsp, Xsp, l, h, v, 0, 0, d, Wsp,
                                     covModel$kernel, use_symmetry)
    return(KK)
  }

  # Gram/feature kernels: operate on the raw feature vectors X, Y
  # directly (inner products), not on pairwise Euclidean distances.
  if(covModel$kernel %in% c("linear", "polynomial")){
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

# Documented kernel-name aliases (see covm()'s @param covModel) that
# resolve to an existing kernel implementation with a fixed smoothness
# 'v' rather than being separate functions in their own right.
# Previously "matern_3_2", "matern_5_2" and "exponential" were advertised
# in covm()'s documentation but .kernelName() would map them to
# nonexistent functions (e.g. "matern_3_2" -> "kMatern_3_2").
.resolveKernelAlias <- function(covModel){
  alias_v <- list(exponential = 0.5, matern_1_2 = 0.5,
                   matern_3_2  = 1.5, matern_5_2  = 2.5)
  k <- covModel$kernel
  if(!is.null(k) && k %in% names(alias_v)){
    covModel$kernel <- "matern"
    covModel$v <- alias_v[[k]]
  }
  return(covModel)
}

# General helper to ensure W is a matrix of correct dimensions
make_W <- function(X, w, Y = NULL) {
  if (is.matrix(w)) {
    return(w)
  } 
  # else {
  #   return(matrix(w, nrow = nrow(X), ncol = ncol(X)))
  # }
  # If Y is omitted, X is assumed to already represent
  # the pairwise matrix dimensions.
  if (is.null(Y)) {
    nx <- nrow(X)
    ny <- ncol(X)
  } else {
    nx <- if (is.null(dim(X))) length(X) else nrow(X)
    ny <- if (is.null(dim(Y))) length(Y) else nrow(Y)
  }
  
  matrix(w, nrow = nx, ncol = ny)
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

# @details The `"power_exp"` kernel: \eqn{K(r) = h^2 \exp(-(r/l)^v)},
#   \eqn{0 < v \le 2}. Requires `covModel$l`, `covModel$h` and
#   `covModel$v`. Note this is the standard power-exponential
#   parameterization (no factor of 0.5 in the exponent), so `v = 2` is
#   qualitatively similar to, but not numerically identical to, the
#   `"gaussian"` kernel with the same `l`. Derivatives (`d = 1`, `d = 2`)
#   are supported, but for `v < 1` the first derivative and for any
#   `v` other than `1` or `2` the second derivative are not
#   mathematically finite exactly at distance 0 (a genuine property of
#   this kernel family, not a bug) -- see `kPowerExp()`/`kPowerExp_rcpp_fast`
#   for details. `sparse = TRUE` is not currently supported for this
#   kernel.
kPower_exp <- function(r, para, d = 0, w = 1, use_symmetry = FALSE){
  l <- para$l
  h <- para$h
  v <- para$v
  if(is.null(v)){
    stop("covModel$v is required for the \"power_exp\" kernel (0 < v <= 2).")
  }
  degree <- 0
  c <- 0
  W_mat <- make_W(r, w)
  K <- kernel_dispatch_auto_rcpp(r, r, l, h, v, degree, c, d, W_mat, "power_exp", use_symmetry)
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
  h <- para$h
  # c <- para$c    # bias
  v <- 0
  degree <- 1
  if(is.null(dim(X))) dim(X) <- c(length(X), 1)
  if(is.null(dim(Y))) dim(Y) <- c(length(Y), 1)
  W_mat <- make_W(X, w, Y)
  # Subtract bias from features if desired
  # (the 'l' slot of kernel_dispatch_auto_rcpp is unused for the linear
  # kernel, so we pass 0 as a placeholder)
  K <- kernel_dispatch_auto_rcpp(X - para$c, Y - para$c, 
                                 0, h, v, degree, para$c, d, W_mat, 
                                 "linear", use_symmetry)
  return(K)
}

kPolynomial <- function(X, Y, para, d = 0, w = 1, use_symmetry = FALSE){
  degree <- para$degree
  h <- para$h
  if(is.null(dim(X))) dim(X) <- c(length(X), 1)
  if(is.null(dim(Y))) dim(Y) <- c(length(Y), 1)
  W_mat <- make_W(X, w, Y)
  # NOTE: unlike kLinear, X/Y are NOT shifted by para$c here: the C++
  # polynomial kernel is h^2 * (X %*% t(Y) + c)^degree, i.e. c enters
  # additively *inside* the dot product, not as a shift of the features.
  # ('l'/'v' slots of kernel_dispatch_auto_rcpp are unused for polynomial,
  # so 0 is passed as a placeholder for each)
  K <- kernel_dispatch_auto_rcpp(X, Y, 0, h, 0, degree, para$c, d, W_mat,
                                 "polynomial", use_symmetry)
  return(K)
}


#' Cross-distance between two matrices
#'
#' Compute the pairwise (Euclidean, or Mahalanobis if `M` is supplied)
#' distance between every row of `X` and every row of `Y`.
#' The returned distance matrix has dimension `nrow(X) x nrow(Y)`.
#' If `M` is `NULL` (the default), the distance is isotropic Euclidean;
#' otherwise the distance is anisotropic (Mahalanobis-type, using `M` as
#' the metric).
#' @param X a matrix or vector.
#' @param Y a matrix or vector with the same number of columns as `X`.
#' @param M optional positive semidefinite matrix (`nrow(M) = ncol(M) =
#'   ncol(X)`) defining an anisotropic metric. `NULL` (default) for plain
#'   isotropic Euclidean distance.
#' @param use_symmetry logical; if `TRUE`, only the upper triangle is
#'   computed and mirrored (faster, and exact only when `X` and `Y`
#'   represent the same set of points). Default `FALSE`.
#' @return a numeric matrix of dimension `nrow(X) x nrow(Y)`.
#' @examples
#' crossDist(c(-1, 0, 1), c(-1, 0, 1), use_symmetry = TRUE)
#' @name crossDist
#' @export
crossDist <- function(X, Y, M = NULL, use_symmetry = FALSE){
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
