
if(!require("devtools")) install.packages("devtools")
devtools::install_github("emanuelhuber/GauProMod")



library(GauProMod)
library(plot3D)
library(RColorBrewer)

M <- matrix(c(1L, 0L, 0L, 1L), nrow = 2)
cholfac(M)

M1 <- M
storage.mode(M1) <- "numeric"


obs <- list(x = cbind(c(2.17, 7.92, 8.98, 7.77, 2.79, 5.36, 4.27, 3.07, 6.31),
                      c(1.33, 7.24, 4.26, 2.67, 6.17, 8.04, 3.18, 5.63, 8.33)),
            y = c(2.60, 1.48, 1.36, 8.61, 1.00, 1.58, 8.42, 8.39, 1.50))

targ <- list(x = cbind(c(2.17, 7.92, 8.98, 7.77, 2.79, 5.36, 4.27, 3.07, 6.31, 
                         3.74, 5.93, 7.19, 6.61, 5.54, 2.27, 1.61, 4.02, 1.06),
                       c(1.33, 7.24, 4.26, 2.67, 6.17, 8.04, 3.18, 5.63, 8.33,
                         6.34, 3.68, 6.82, 1.79, 8.60, 7.73, 5.35, 2.45, 4.92))
)


  # Matern kernel
  covModel <- list(kernel="matern",
                   l = 5,     # correlation length
                   v = 1,     # smoothness
                   h = 2.45   # std. deviation
  )

# 2D linear mean function
op <- 2

# Gaussian likelihood
sigma <- 0.2


GP <- gpCond(obs = obs, targ = targ, covModels=list(pos=covModel), 
             sigma = sigma, op = op)

Kxx       <- covm( obs$x,  obs$x, covModels[[1]])
Kstar     <- covm( obs$x, targ$x, covModels[[1]])
Kstarstar <- covm(targ$x, targ$x, covModels[[1]])

# Error in cholfac_rcpp(x) : Wrong R type for mapped matrix

x <- obs$x
y <- obs$x
covModel <- covModels[[1]]
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
    M <- as.numeric(M)
    XY <- crossDist(x, y, M)
    
    X <- x
    Y <- y
    M
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
    
    
    
    
    
    #if(dim(x)[2] == 2){
    #    XY <- dist2(x, y, M)
    #}else{
    #  XY <- distn(x,y, M)
    #}
    if(d == 1){
      if(is.null(dim(x)) && is.null(dim(y))){
        w0 <- sign(outer(x, y, "-"))
        # w1 <- outer(rep(1,length(x)), dx, function(x,y){sign2(y)})
        # w <- (w0*w1)
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



#observations
obs <- list(x=c(-4, -3, -1, 0, 4),
            y=c(-2,  0,  1, 2, 0))
# targets
targ <- list("x"=seq(-10,10,len=200))

GP <- gpCond(obs = obs, targ = targ, covModels=list(pos=covModel), 
             sigma = sigma, op = op)
names(GP)






names(GP)
# GP$mean   = mean value at location xstar
# GP$cov    = covariance matrix of the conditioned GP
# GP$logLik = log-likelihood of the conditioned GP
# GP$xstar  = x-coordinates at which the GP is simulated

# mean
Ymean <- GP$mean
# standard deviation
YSD <- sqrt(diag(GP$cov))
ysdplus <- Ymean - 1.95* YSD
ysdminus <- Ymean + 1.95* YSD