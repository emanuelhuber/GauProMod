# GauProMod
R functions for Gaussian process (GP) modelling. The core functions are coded 
in C++ and based on the EIGEN library (through RcppEigen)

## Notes
Currently implemented:
- [x] GP conditioned to data
- [x] Space-time Gaussian process
- [x] GP with monomial mean functions
- [x] GP conditioned to derivative observations
- [ ] add function for hyper-parameter estimation
- [ ] add other covariance models

This is an ongoing project.
If you have any questions, don't hesitate to contact me:

emanuel.huber@stanford.edu

Thank you!

## How to install/load

```r
library(devtools)
devtools::install_github("emanuelhuber/GauProMod")
library(GauProMod)
```

## Short tutorial
### Load libraries
```r
library(GauProMod)
library(plot3D)
library(RColorBrewer)
```

### 1D Gaussian Process Modelling

#### Observations and  target
The observations are defined by a list with `x` the positions of the 
observations and `y` the observed values. The targets are the positions `x` 
where to simulate the Gaussian Process.

```r
#observations
obs <- list(x=c(-4, -3, -1, 0, 4),
            y=c(-2,  0,  1, 2, 0))
# targets
targ <- list("x"=seq(-10,10,len=200))
```
#### Covariance function, mean function and likelihood
To build the covariance functions, the following kernels are available
```r
# linear kernel
covModel <- list(type="linear",
                 b = 1,         # slope
                 h = 1.5,       # std. deviation
                 c = 0)         # constant
                 
# Matern kernel
covModel <- list(type="matern",
                 l = 1,       # correlation length
                 v = 2.5,     # smoothness
                 h = 2.45)    # std. deviation

# squared exponential kernel (Gaussian)
covModel <- list(type="se",
                 l = 0.5,   # correlation length
                 h = 0.25)  # std. deviation
```

The following mean function (or basis functions) are available 
(see Rasmussen and Williams (2006), Section 2.7): 

```r
# quadratic mean function
op <- 2
# linear mean function
op <- 1
# zero-mean function (no trend)
op <- 0 
```

Because nothing is perfectly observed, it makes sense to account for uncertainty
in the observation. Gaussian likelihood, defined by the standard deviation 
`sigma` (erreur uncertainty) is the only form of likelihood 
currently implemented in GauProMod

```r
# standard deviation measurement error
# Gaussian likelihood
sigma <- 0.2
```


#### Conditional Gaussian Process modelling

```r
para <- gpCond(obs = obs, targ = targ, covModels=list(pos=covModel), 
               sigma = sigma, op = op)
names(para)
# para$mean   = mean value at location xstar
# para$cov    = covariance matrix of the conditioned GP
# para$logLik = log-likelihood of the conditioned GP
# para$xstar  = x-coordinates at which the GP is simulated
```

Plot the mean function plus/minus the standard deviation

```r
#--- plot mean +/- sd
xp <-(para$mean + sqrt(diag(para$cov)))  # mean + sd
xm <-(para$mean - sqrt(diag(para$cov)))  # mean - sd

# initialise the plot
plot(cbind(obs$x, obs$y), type="p", xlab="x", ylab="y", 
     xlim = range(c(obs$x, targ$x)), ylim = range(c(xp, xm, obs$y)),
     pch = 20, col = "black") 
lines(para$xstar, para$mean,col="red")  # mean
lines(para$xstar, xm,lty=3)            # + sd
lines(para$xstar, xp,lty=3)            # - sd
legend("topleft", legend = c("obs", "mean", "sd"), lty = c(NA, 3, 1),
       pch = c(20, NA, NA), col=c("black", "red", "black"), bty="n")
```
       
Random conditional simulation

```r
# cholesky factorisation
L <- cholfac(para$cov)
# random simulation
ystar <- gpSim(para , L = L)
```

You can also directly use `ystar <- gpSim(para)` without the argument `L` (the
Cholesky factor) but each time you will call `gpSim(para)`, `gpSim` will 
compute again internally the Cholesky factor. So, if you plan to run many 
unconditional simulations, it is faster to first compute the Cholesky factor
and then run several time `gpSim` with the argument `L`.

Plot the random simulation:
```r
plot(rbind(cbind(obs$x, obs$y), ystar), type="n", xlab="x", ylab="y") 
lines(ystar, col = "blue")
points(cbind(obs$x, obs$y), col = "black", pch = 20)
legend("topleft", legend = c("obs", "GP sim"), lty = c(NA, 1),
       pch = c(20, NA), col=c("black", "blue"), bty="n")
```

### 2D Gaussian Process Modelling
To understand everything, please read the previous section ("1D Gaussian
Process Modelling").

#### Observations and  target
Here, the element `x` of the observation list is a matrix corresponding to
the coordinates of the observations points (East/North coordinates or 
x/y coordinates).



```r
#observations
obs <- list(x = cbind(c(2.17, 7.92, 8.98, 7.77, 2.79, 5.36, 4.27, 3.07, 6.31, 
                       3.74, 5.93, 7.19, 6.61, 5.54, 2.27, 1.61, 4.02, 1.06),
                     c(1.33, 7.24, 4.26, 2.67, 6.17, 8.04, 3.18, 5.63, 8.33,
                       6.34, 3.68, 6.82, 1.79, 8.60, 7.73, 5.35, 2.45, 4.92)),
            y = c(2.60, 1.48, 1.36, 8.61, 1.00, 1.58, 8.42, 8.39, 1.50, 
                  9.05, 1.14, 1.49, 9.19, 1.32, 1.03, 6.41, 6.16, 5.42))
```

The target is defined by a regular grid defined by two orthogonal vectors.
The function `vecGrid`returns a two-columns matrix corresponding to the
coordinates of each element of the grid.

```r
# targets
vx <- seq(0, 10, by = 0.5)
vy <- seq(0, 10, by = 0.5)
targ <- list(x = vecGrid(vx, vy))
```

#### Covariance function, mean function and likelihood
To build the covariance functions, the same kernels as in the previously 
defined are available:
```r
# linear kernel
covModel <- list(type="linear",
                 b = 1,         # slope
                 h = 1.5,       # std. deviation
                 c = 0)         # constant
                 
# Matern kernel
covModel <- list(type="matern",
                 l = 1,       # correlation length
                 v = 2.5,     # smoothness
                 h = 2.45)    # std. deviation

# squared exponential kernel (Gaussian)
covModel <- list(type="se",
                 l = 0.5,   # correlation length
                 h = 0.25)  # std. deviation
```

Note that the 2D mean functions (or basis functions) are differently defined: 

```r
# 2D quadratic mean function
op <- 5
# 2D linear mean function
op <- 2
# zero-mean function (no trend)
op <- 0 
```

Standard deviation (measurement error):
```r
# Gaussian likelihood
sigma <- 0.2
```

#### Conditional Gaussian Process modelling

```r
para <- gpCond(obs = obs, targ = targ, covModels=list(pos=covModel), 
               sigma = sigma, op = op)
names(para)
# para$mean   = mean value at location xstar
# para$cov    = covariance matrix of the conditioned GP
# para$logLik = log-likelihood of the conditioned GP
# para$xstar  = x-coordinates at which the GP is simulated
```

Plot the mean function 
```r
Ymean <- matrix(para$mean, nrow = length(vx), ncol = length(vy), byrow = TRUE)

plot3D::image2D(x = vx, y = vy, z = Ymean)
points(obs$x, col="white",pch=3)


plot3D::contour2D(x = vx, y = vy, Ymean)
points(obs$x, col="black",pch=3)
rect(vx[1], vy[1], vx[length(vx)], vy[length(vy)])
```

Plot the standard deviation

```r
YSD <- matrix(sqrt(diag(para$cov)), nrow = length(vx), ncol = length(vy), 
              byrow = TRUE)
              
plot3D::image2D(x = vx, y = vy, z = YSD)
points(obs$x, col="white",pch=3)


plot3D::contour2D(x = vx, y = vy, YSD)
points(obs$x, col="black",pch=3)
rect(vx[1], vy[1], vx[length(vx)], vy[length(vy)])
```
       
Random conditional simulation

```r
L <- cholfac(para$cov)
ystar <- gpSim(para , L = L)

Ysim <- matrix(ystar[,3], nrow = length(vx), ncol = length(vy), byrow = TRUE)

plot3D::image2D(x = vx, y = vy, z = Ysim)
points(obs$x, col="white",pch=3)


plot3D::contour2D(x = vx, y = vy, Ysim)
points(obs$x, col="black",pch=3)
rect(vx[1], vy[1], vx[length(vx)], vy[length(vy)])
```


## References
Rasmussen C.E. and Williams C.K.I. (2006), Gaussian Processes for Machine 
Learning, the MIT Press, ISBN 026218253X.
www.GaussianProcess.org/gpml
