---
layout: default
title: Home
date: 2018-03-07
---

# Gaussian Process Modelling

<!--
<p class="message">
  Gaussian Process Modelling
</p>
-->

R functions for Gaussian process (GP) modelling. The core functions are coded 
in C++ and based on the EIGEN library (through RcppEigen)

# Notes
Current features:
- [x] Posterior Gaussian Process with Gaussian likelihood (Gaussian process
      conditioned to noise-free and noisy observations)
- [x] Space-time Gaussian process 
- [x] Gaussian Process with monomial mean functions with vague Gaussian prior
      on the coefficient parameters.
- [x] Gaussian Process conditioned to derivative observations
- [x] Covariance function: Matern, Gaussian, linear
- [x] Anisotropic covariance functions (scale and rotation)
- [x] Log marginal likelihood of the Gaussian process
- [x] Cross-matrix distance (distance between every rows of each matrix):
      `crossdist(x,y,M)` (with `M` a positive semidefinite matrix for
      anisotropic distances)
- [ ] maximum likelihood hyper-parameter estimation
- [ ] McMC hyper-parameter sampling
- [ ] spatially varying covariance function
- [ ] Gaussian Process approximations (to deal with larger data set)
- [ ] add other covariance models

This is an ongoing project. If you have any questions, requirements, suggestions, don't hesitate to contact me: emanuel.huber@alumni.ethz.ch

Thank you!

# How to install/load

## 1. Install R

Download R from the [R Cran website](http://cran.r-project.org) and install it.
2. Optionally install a R-editor:
  * [Notepad++](https://notepad-plus-plus.org/) combined with [NppToR](https://sourceforge.net/projects/npptor/)
  * [Rstudio](https://www.rstudio.com/)
  * [RKward](https://rkward.kde.org/)
3. If necessary, learn some R basics:
  * [An interactive introduction to R](http://tryr.codeschool.com)
  * [Short R introduction](http://cran.r-project.org/doc/contrib/Torfs+Brauer-Short-R-Intro.pdf) 
  * [R course](http://www.rochester.edu/college/psc/thestarlab/help/rcourse/R-Course.pdf)


## 2. Install a R-editor [optionally]

Possible R-editor choices:
* [Rstudio](https://www.rstudio.com/) (recommanded if you start)
* [Notepad++](https://notepad-plus-plus.org/) combined with [NppToR](https://sourceforge.net/projects/npptor/)
* [RKward](https://rkward.kde.org/)


## 3. Install & load `GauProMod`

```r
if(!require("devtools")) install.packages("devtools")
devtools::install_github("emanuelhuber/GauProMod")
library(GauProMod)
```



<!--

2. [Learn some R basics](02_rbasics)
3. [Learn to use RStudio](03_rstudio)
-->

<!--
$$\forall x \in R$$
-->