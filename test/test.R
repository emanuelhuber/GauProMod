
if(!require("devtools")) install.packages("devtools")
devtools::install_github("emanuelhuber/GauProMod")

devtools::install_local("/media/huber/Elements/UNIBAS/software/codeR/package_GauProMod/GauProMod")



library(GauProMod)
library(plot3D)
library(RColorBrewer)

#------------------------------------------------------------------------------#


#observations
obs <- list(x=c(-4, -3, -1, 0, 4),
            y=c(-2,  0,  1, 2, 0))
# targets
targ <- list("x"=seq(-10,10,len=200))

covModel <- list(kernel="gaussian",
                 l = 0.7,   # correlation length
                 h = 1.)  # std. deviation

# squared exponential kernel (Gaussian)
covModel <- list(kernel="gaussian",
                 l = 0.5,   # correlation length
                 h = 0.25)  # std. deviation

A1 <- covm( obs$x, targ$x, covModel)
A2 <- covm( targ$x, obs$x, covModel)
A1 - A2

# linear mean function
op <- 0

sigma <- 0.2



GP <- gpCond(obs = obs, targ = targ, covModels=list(pos=covModel), 
             sigma = sigma, op = op)
names(GP)
# GP$mean   = mean value at location xstar
# GP$cov    = covariance matrix of the conditioned GP
# GP$logLik = log-likelihood of the conditioned GP
# GP$xstar  = x-coordinates at which the GP is simulated

#--- plot mean +/- sd
xp <-(GP$mean + sqrt(diag(GP$cov)))  # mean + sd
xm <-(GP$mean - sqrt(diag(GP$cov)))  # mean - sd

# initialise the plot
plot(cbind(obs$x, obs$y), type="p", xlab="x", ylab="y", 
     xlim = range(c(obs$x, targ$x)), ylim = range(c(xp, xm, obs$y)),
     pch = 20, col = "black", main = GP$logLik) 
lines(GP$xstar, GP$mean,col="red")  # mean
lines(GP$xstar, xm,lty=3)            # + sd
lines(GP$xstar, xp,lty=3)            # - sd
legend("topleft", legend = c("obs", "mean", "sd"), lty = c(NA, 1, 3),
       pch = c(20, NA, NA), col=c("black", "red", "black"), bty="n")


# random conditional simulations
# cholesky factorisation
L <- cholfac(GP$cov)
# random simulation
ystar <- gpSim(GP , L = L)

plot(rbind(cbind(obs$x, obs$y), ystar), type="n", xlab="x", ylab="y") 
lines(ystar, col = "blue")
points(cbind(obs$x, obs$y), col = "black", pch = 20)
legend("topleft", legend = c("obs", "GP sim"), lty = c(NA, 1),
       pch = c(20, NA), col=c("black", "blue"), bty="n")


# 1D GP with derivative

obs <- list(x=c(-4, -3, -1, 0, 4),
            y=c(-2,  0,  1, 2, 0))
# targets
targ <- list("x"=seq(-10,10,len=200))

covModel <- list(kernel = "matern",
                 l = 0.25,
                 v = 3.5,
                 h = 0.55)

bc <- list(x = c(-4.5, -2, 0,  3, 4.5),
           y = c(   0,  1, 0, -1,   0),
           sigma = 0)

sigma <- 0.1

GP <- gpCond(obs = obs, targ = targ, covModels=list(pos=covModel), 
             sigma = sigma, op = 0)
GP2 <- gpCond(obs = obs, targ = targ, covModels=list(pos=covModel), 
              sigma = sigma, op = 0, bc = bc)

#--- plot mean +/- sd
xp <-(GP$mean + sqrt(diag(GP$cov)))  # mean + sd
xm <-(GP$mean - sqrt(diag(GP$cov)))  # mean - sd
xp2 <-(GP2$mean + sqrt(diag(GP2$cov)))  # mean + sd
xm2 <-(GP2$mean - sqrt(diag(GP2$cov)))  # mean - sd


plot(cbind(obs$x, obs$y), type="p", xlab="x", ylab="y", 
     xlim = range(c(obs$x, targ$x)), ylim = range(c(xp, xm, obs$y)),
     pch = 20, col = "black", main = "without derivatives") 
lines(GP$xstar, GP$mean,col="red")  # mean
lines(GP$xstar, xm,lty=3)            # + sd
lines(GP$xstar, xp,lty=3)            # - sd
legend("topleft", legend = c("obs", "mean", "sd"), lty = c(NA, 1, 3),
       pch = c(20, NA, NA), col=c("black", "red", "black"), bty="n")

plot(cbind(obs$x, obs$y), type="p", xlab="x", ylab="y", 
     xlim = range(c(obs$x, targ$x)), ylim = range(c(xp, xm, obs$y)),
     pch = 20, col = "black", main = "with derivatives") 
lines(GP2$xstar, GP2$mean, col = "red")  # mean
lines(GP2$xstar, xm2,lty=3)            # + sd
lines(GP2$xstar, xp2,lty=3)            # - sd
legend("topleft", legend = c("obs", "mean", "sd"), lty = c(NA, 1, 3),
       pch = c(20, NA, NA), col=c("black", "red", "black"), bty="n")
y0 <- GP2$mean[sapply(bc$x, function(x, a) which.min(abs(x - a)), GP2$xstar)]
arrows(x0 = bc$x - 1/2, y0 = y0 - bc$y/2, 
       x1 = bc$x + 1/2, y1 = y0 + bc$y/2,
       length = 0.15, col = "dodgerblue2", lwd = 2)


covModel <- list(kernel = "matern",
                 l = 0.25,
                 v = 3.5,
                 h = 0.55)

# squared exponential kernel (Gaussian)
covModel <- list(kernel="gaussian",
                 l = 0.25,   # correlation length
                 h = 0.5)  # std. deviation

# - `x` is the locations where we set the derivative of the Gaussian field, 
# - `v` the gradient derivative, i.e., a unit vector normal to the no-flow boundary (or tangent to a constant-head boundary)
# - `y` is the value of the gradient (in this case `0` meaning that the gradient is flat)
# - `sigma` the standard deviation reprensenting the uncertainty on the value of the gradient (i.e., `y`).
# 

bc <- list(x = c(-4.5, -2, 0,  3, 4.5),
           y = c(   0,  1, 0, -1,   0),
           sigma = 0)

sigma <- 0.1

GP <- gpCond(obs = obs, targ = targ, covModels=list(pos=covModel), 
             sigma = sigma, op = 0)
GP2 <- gpCond(obs = obs, targ = targ, covModels=list(pos=covModel), 
              sigma = sigma, op = 0, bc = bc)


#--- plot mean +/- sd
xp <-(GP$mean + sqrt(diag(GP$cov)))  # mean + sd
xm <-(GP$mean - sqrt(diag(GP$cov)))  # mean - sd
xp2 <-(GP2$mean + sqrt(diag(GP2$cov)))  # mean + sd
xm2 <-(GP2$mean - sqrt(diag(GP2$cov)))  # mean - sd

plot(cbind(obs$x, obs$y), type="p", xlab="x", ylab="y", 
     xlim = range(c(obs$x, targ$x)), ylim = range(c(xp, xm, obs$y)),
     pch = 20, col = "black", main = "without derivatives") 
lines(GP$xstar, GP$mean,col="red")  # mean
lines(GP$xstar, xm,lty=3)            # + sd
lines(GP$xstar, xp,lty=3)            # - sd
legend("topleft", legend = c("obs", "mean", "sd"), lty = c(NA, 1, 3),
       pch = c(20, NA, NA), col=c("black", "red", "black"), bty="n")


plot(cbind(obs$x, obs$y), type="p", xlab="x", ylab="y", 
     xlim = range(c(obs$x, targ$x)), ylim = range(c(xp, xm, obs$y)),
     pch = 20, col = "black", main = "with derivatives") 
lines(GP2$xstar, GP2$mean, col = "red")  # mean
lines(GP2$xstar, xm2,lty=3)            # + sd
lines(GP2$xstar, xp2,lty=3)            # - sd
legend("topleft", legend = c("obs", "mean", "sd"), lty = c(NA, 1, 3),
       pch = c(20, NA, NA), col=c("black", "red", "black"), bty="n")
y0 <- GP2$mean[sapply(bc$x, function(x, a) which.min(abs(x - a)), GP2$xstar)]
arrows(x0 = bc$x - 1/2, y0 = y0 - bc$y/2, 
       x1 = bc$x + 1/2, y1 = y0 + bc$y/2,
       length = 0.15, col = "dodgerblue2", lwd = 2)





#------------------------------------------------------------------------------#

M <- matrix(c(1L, 0L, 0L, 1L), nrow = 2)
cholfac(M)

M1 <- M
storage.mode(M1) <- "numeric"
cholfac(M1)



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



par(mfrow = c(1,1))
ylim <- range(Ymean, obs$y)
plot3D::scatter3D(x = targ$x[,1], y = targ$x[,2], z = Ymean, clim = ylim, 
                  pch = 20)
plot3D::arrows3D(x0 = targ$x[,1], y0 = targ$x[,2], z0 = ysdminus, 
                 x1 = targ$x[,1], y1 = targ$x[,2], z1 = ysdplus,
                 length=0.05, angle=90, code=3, add = TRUE, col = "black")
# large dots = observations
plot3D::scatter3D(x = obs$x[,1], y = obs$x[,2], z = obs$y, add = TRUE, 
                  pch = 20, cex = 3, clim = ylim)





par(mfrow = c(1, 2))
ylim <- range(ysdplus, ysdminus, obs$y)

plot(targ$x[,1], Ymean, type = "p", ylim = ylim, pch = 3, cex = 0.5)
arrows(targ$x[,1], ysdminus, targ$x[,1], ysdplus, length=0.05, angle=90, code=3)
points(obs$x[,1], obs$y, col = "dodgerblue", pch = 20)

plot(targ$x[,2], Ymean, type = "p", ylim = ylim, pch = 3, cex = 0.5)
arrows(targ$x[,2], ysdminus, targ$x[,2], ysdplus, length=0.05, angle=90, code=3)
points(obs$x[,2], obs$y, col = "dodgerblue", pch = 20)


L <- cholfac(GP$cov)
ystar <- gpSim(GP , L = L)

par(mfrow = c(1,1))
ylim <- range(ystar[,3], obs$y)
plot3D::scatter3D(x = targ$x[,1], y = targ$x[,2], z = ystar[,3], clim = ylim, 
                  pch = 18)
# dots = observations
plot3D::scatter3D(x = obs$x[,1], y = obs$x[,2], z = obs$y, add = TRUE, 
                  pch = 20, cex = 3, clim = ylim)



#observations
obs <- list(x = cbind(c(2.17, 7.92, 8.98, 7.77, 2.79, 5.36, 4.27, 3.07, 6.31, 
                        3.74, 5.93, 7.19, 6.61, 5.54, 2.27, 1.61, 4.02, 1.06),
                      c(1.33, 7.24, 4.26, 2.67, 6.17, 8.04, 3.18, 5.63, 8.33,
                        6.34, 3.68, 6.82, 1.79, 8.60, 7.73, 5.35, 2.45, 4.92)),
            y = c(2.60, 1.48, 1.36, 8.61, 1.00, 1.58, 8.42, 8.39, 1.50, 
                  9.05, 1.14, 1.49, 9.19, 1.32, 1.03, 6.41, 6.16, 5.42))




# targets (=2D mesh)
vx <- seq(0, 10, by = 0.5)
vy <- seq(0, 10, by = 0.5)
targ <- list(x = vecGrid(vx, vy))


# linear kernel
covModel <- list(kernel="linear",
                 b = 1,         # slope
                 h = 0.5,       # std. deviation
                 c = 1)         # constant

# Matern kernel
covModel <- list(kernel="matern",
                 l = 1,       # correlation length
                 v = 2.5,     # smoothness
                 h = 2.45)    # std. deviation

# squared exponential kernel (Gaussian)
covModel <- list(kernel="gaussian",
                 l = 0.5,   # correlation length
                 h = 0.25)  # std. deviation


# 2D quadratic mean function
op <- 5
# zero-mean function (no trend)
op <- 0 
# 2D linear mean function
op <- 2

sigma <- 0.2


GP <- gpCond(obs = obs, targ = targ, covModels=list(pos=covModel), 
             sigma = sigma, op = op)
names(GP)
# GP$mean   = mean value at location xstar
# GP$cov    = covariance matrix of the conditioned GP
# GP$logLik = log-likelihood of the conditioned GP
#


# mean
Ymean <- matrix(GP$mean, nrow = length(vx), ncol = length(vy), byrow = TRUE)
# standard deviation
YSD <- matrix(sqrt(diag(GP$cov)), nrow = length(vx), ncol = length(vy), 
              byrow = TRUE)

par(mfrow = c(2,2))
plot3D::image2D(x = vx, y = vy, z = Ymean, asp=1)
points(obs$x, col="white",pch=3)
title(main = "mean")

plot3D::contour2D(x = vx, y = vy, Ymean, asp=1)
points(obs$x, col="black",pch=3)
rect(vx[1], vy[1], vx[length(vx)], vy[length(vy)])
title(main = "mean")

plot3D::image2D(x = vx, y = vy, z = YSD, asp=1)
points(obs$x, col="white",pch=3)
title(main = "standard deviation")

plot3D::contour2D(x = vx, y = vy, YSD, asp=1)
points(obs$x, col="black",pch=3)
rect(vx[1], vy[1], vx[length(vx)], vy[length(vy)])
title(main = "standard deviation")

# anisotropy

covModelAni <- list(kernel="matern",
                    l = 1,       # correlation length
                    v = 2.5,     # smoothness
                    h = 2.45,
                    scale = c(1, 0.25))    # std. deviation
# 2D linear mean function
op <- 2

GP <- gpCond(obs = obs, targ = targ, covModels=list(pos=covModelAni), 
             sigma = sigma, op = op)
names(GP)
# GP$mean   = mean value at location xstar
# GP$cov    = covariance matrix of the conditioned GP
# GP$logLik = log-likelihood of the conditioned GP
# GP$xstar  = x-coordinates at which the GP is simulated

# mean
YmeanAni <- matrix(GP$mean, nrow = length(vx), ncol = length(vy), byrow = TRUE)
# standard deviation
YSDAni <- matrix(sqrt(diag(GP$cov)), nrow = length(vx), ncol = length(vy), 
                 byrow = TRUE)

par(mfrow = c(2,2))
plot3D::image2D(x = vx, y = vy, z = Ymean, asp = 1)
points(obs$x, col="white",pch=3)
title(main = "isotropic GP: mean ")

plot3D::contour2D(x = vx, y = vy, Ymean, asp = 1)
points(obs$x, col="black",pch=3)
rect(vx[1], vy[1], vx[length(vx)], vy[length(vy)])
title(main = "isotropic GP: mean")

plot3D::image2D(x = vx, y = vy, z = YmeanAni, asp = 1)
points(obs$x, col="white",pch=3)
title(main = "anisotropic GP: mean ")

plot3D::contour2D(x = vx, y = vy, YmeanAni, asp = 1)
points(obs$x, col="black",pch=3)
rect(vx[1], vy[1], vx[length(vx)], vy[length(vy)])
title(main = "anisotropic GP: mean")

#---------- anisotropy + scaling

covModelAni2 <- list(kernel="matern",
                     l = 1,       # correlation length
                     v = 2.5,     # smoothness
                     h = 2.45,
                     scale = c(1, 0.25),
                     rot = c(1.0))    # std. deviation
# 2D linear mean function
op <- 2

GP <- gpCond(obs = obs, targ = targ, covModels=list(pos=covModelAni2), 
             sigma = sigma, op = op)
names(GP)
# GP$mean   = mean value at location xstar
# GP$cov    = covariance matrix of the conditioned GP
# GP$logLik = log-likelihood of the conditioned GP
# GP$xstar  = x-coordinates at which the GP is simulated
# mean
YmeanAni2 <- matrix(GP$mean, nrow = length(vx), ncol = length(vy), byrow = TRUE)
# standard deviation
YSDAni2 <- matrix(sqrt(diag(GP$cov)), nrow = length(vx), ncol = length(vy), 
                  byrow = TRUE)

par(mfrow = c(2,2))
plot3D::image2D(x = vx, y = vy, z = YmeanAni, asp = 1)
points(obs$x, col="white",pch=3)
title(main = "anisotropic GP (scale): mean ")

plot3D::contour2D(x = vx, y = vy, YmeanAni, asp = 1)
points(obs$x, col="black",pch=3)
rect(vx[1], vy[1], vx[length(vx)], vy[length(vy)])
title(main = "anisotropic GP (scale): mean")

plot3D::image2D(x = vx, y = vy, z = YmeanAni2, asp = 1)
points(obs$x, col="white",pch=3)
title(main = "anisotropic GP (scale + rotation): mean ")

plot3D::contour2D(x = vx, y = vy, YmeanAni2, asp = 1)
points(obs$x, col="black",pch=3)
rect(vx[1], vy[1], vx[length(vx)], vy[length(vy)])
title(main = "anisotropic GP (scale + rotation): mean")


# ----------------------------

obs <- list(x = cbind(c(2.17, 7.92, 8.98, 7.77, 2.79, 5.36, 4.27, 3.07, 6.31),
                      c(1.33, 7.24, 4.26, 2.67, 6.17, 8.04, 3.18, 5.63, 8.33)),
            y = c(2.60, 1.48, 1.36, 8.61, 1.00, 1.58, 8.42, 8.39, 1.50))



# Matern kernel

vx <- seq(0, 10, by = 0.25)
vy <- seq(0, 10, by = 0.25)
targ <- list(x = vecGrid(vx, vy))


covModel <- list(kernel="matern",
                 l = 1,
                 v = 2.5,
                 h = 3)

op <- 5    
sigma <- 0.05


GP <- gpCond(obs = obs, targ = targ, covModels=list(pos=covModel), 
             sigma = sigma, op = op)
# mean
Ymean <- matrix(GP$mean, nrow = length(vx), ncol = length(vy), byrow = TRUE)


# - `x` is the locations where we set the derivative of the Gaussian field, 
# - `v` the gradient derivative, i.e., a unit vector normal to the no-flow boundary (or tangent to a constant-head boundary)
# - `y` is the value of the gradient (in this case `0` meaning that the gradient is flat)
# - `sigma` the standard deviation reprensenting the uncertainty on the value of the gradient (i.e., `y`).
# 

#  no-flow boundary top and bottom
bc <- list(x = cbind(c( rep(seq(0.5, 9.5, by = 2), 2)),
                     c(rep(0, 5), rep(10, 5))),
           v = cbind( rep(0, 10),
                      rep(1, 10) ),
           y =  rep(0, 10),
           sigma = 0) 



GPbc <- gpCond(obs = obs, targ = targ, covModels=list(pos=covModel), 
               sigma = sigma, op = op, bc = bc)
# mean
Ymeanbc <- matrix(GPbc$mean, nrow = length(vx), ncol = length(vy), byrow = TRUE)

par(mfrow = c(1,2))


plot3D::contour2D(x = vx, y = vy, Ymean, asp=1, nlevels = 20, col = "black", lwd = 2, 
                  xaxs = "i", yaxs = "i", labcex = 1)
points(obs$x, pch=20, col = "red", cex = 2)
#rect(vx[1], vy[1], vx[length(vx)], vy[length(vy)])
title(main = "mean GP")



plot3D::contour2D(x = vx, y = vy, Ymeanbc, asp=1, nlevels = 20, col = "black", lwd = 2, 
                  xaxs = "i", yaxs = "i", labcex = 1)
points(obs$x, pch=20, col = "red", cex = 2)
#rect(vx[1], vy[1], vx[length(vx)], vy[length(vy)])
title(main = "mean GP with no-flow\n boundaries")

points(bc$x, pch = 4)
arrows(bc$x[,1] - bc$v[,2]/2, bc$x[,2] - bc$v[,1]/2, 
       bc$x[,1] + bc$v[,2]/2, bc$x[,2] + bc$v[,1]/2,
       length = 0.15, col = "dodgerblue2", lwd = 2)



#------------------------------------ LARGER SIZE ----------------------------#

#observations
obs <- list(x = cbind(c(2.17, 7.92, 8.98, 7.77, 2.79, 5.36, 4.27, 3.07, 6.31, 
                        3.74, 5.93, 7.19, 6.61, 5.54, 2.27, 1.61, 4.02, 1.06),
                      c(1.33, 7.24, 4.26, 2.67, 6.17, 8.04, 3.18, 5.63, 8.33,
                        6.34, 3.68, 6.82, 1.79, 8.60, 7.73, 5.35, 2.45, 4.92)),
            y = c(2.60, 1.48, 1.36, 8.61, 1.00, 1.58, 8.42, 8.39, 1.50, 
                  9.05, 1.14, 1.49, 9.19, 1.32, 1.03, 6.41, 6.16, 5.42))




# targets (=2D mesh)
vx <- seq(0, 10, by = 0.1)
vy <- seq(0, 10, by = 0.1)
targ <- list(x = vecGrid(vx, vy))


# linear kernel
covModel <- list(kernel="linear",
                 b = 1,         # slope
                 h = 0.5,       # std. deviation
                 c = 1)         # constant

# Matern kernel
covModel <- list(kernel="matern",
                 l = 1,       # correlation length
                 v = 2.5,     # smoothness
                 h = 2.45)    # std. deviation

# squared exponential kernel (Gaussian)
covModel <- list(kernel="gaussian",
                 l = 0.5,   # correlation length
                 h = 0.25)  # std. deviation


# 2D quadratic mean function
op <- 5
# zero-mean function (no trend)
op <- 0 
# 2D linear mean function
op <- 2

sigma <- 0.2

starttime <- Sys.time()
GP <- gpCond(obs = obs, targ = targ, covModels=list(pos=covModel), 
             sigma = sigma, op = op)
endtime <- Sys.time()

endtime - starttime
# Time difference of 19.06291 secs

GP$logLik
# -712.986  <- version gemini
# -434.9469 <- version gemini modified by me (but not sure if it is correct...)



names(GP)
# GP$mean   = mean value at location xstar
# GP$cov    = covariance matrix of the conditioned GP
# GP$logLik = log-likelihood of the conditioned GP
#


# mean
Ymean <- matrix(GP$mean, nrow = length(vx), ncol = length(vy), byrow = TRUE)
# standard deviation
YSD <- matrix(sqrt(diag(GP$cov)), nrow = length(vx), ncol = length(vy), 
              byrow = TRUE)

par(mfrow = c(2,2))
plot3D::image2D(x = vx, y = vy, z = Ymean, asp=1)
points(obs$x, col="white",pch=3)
title(main = "mean")

plot3D::contour2D(x = vx, y = vy, Ymean, asp=1)
points(obs$x, col="black",pch=3)
rect(vx[1], vy[1], vx[length(vx)], vy[length(vy)])
title(main = "mean")

plot3D::image2D(x = vx, y = vy, z = YSD, asp=1)
points(obs$x, col="white",pch=3)
title(main = "standard deviation")

plot3D::contour2D(x = vx, y = vy, YSD, asp=1)
points(obs$x, col="black",pch=3)
rect(vx[1], vy[1], vx[length(vx)], vy[length(vy)])
title(main = "standard deviation")

