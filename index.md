---
layout: default
title: Home
date: 2026-10-30
---

# Gaussian Process Modelling

<!--
<p class="message">
  Gaussian Process Modelling
</p>
-->

R functions for Gaussian process (GP) modelling. The core functions are coded 
in C++ and based on the EIGEN library (through RcppEigen)

## Feature Implementation Checklist

### Core Capabilities & Completed Features

- [x] **Posterior Gaussian Process with Gaussian Likelihood**
  - Gaussian process conditioned on both noise-free and noisy observations (`gpCond`).
- [x] **Space-Time Gaussian Process Support**
  - Supported via `gpCond` and coordinate mapping utilities such as `setPosTime`.
- [x] **Gaussian Process with Monomial Mean Functions**
  - Supports vague Gaussian priors on coefficient parameters (`gpLogLikMean_rcpp`).
- [x] **Gaussian Process Conditioned on Derivative Observations**
  - Supports derivative data conditioning for specialized spatial-temporal models.
- [x] **Anisotropic Covariance Functions**
  - Generalized scale parameters and custom rotation matrices $M$.
- [x] **Log Marginal Likelihood Evaluation**
  - Fast likelihood evaluation via `gpLogLik`, `gpLogLik_rcpp`, and `gpLogLikMean_rcpp`.
- [x] **Cross-Matrix Distance Calculations**
  - Efficient distance computations: `crossDist(x, y, M)`, `crossDist_rcpp`, and `crossDist_sparse`.
- [x] **Flexible Kernel Selection**
  - Matérn 3/2 & 5/2
  - Gaussian / Squared Exponential (`sqex`)
  - Exponential
  - Linear
- [x] **Maximum Likelihood Hyper-Parameter Estimation (`gpFit`)**
  - Optimization via `gpNegLogLik` and numerical gradients (`.centralDiffGrad`).
- [x] **Gaussian Process Simulation & Trajectory Generation (`gpSim`)**
  - Simulates Gaussian process sample paths across evaluation grids.
- [x] **Robust Numerical Stabilization**
  - Automatic jitter addition (`cholfac`) and positive-definiteness adjustments (`correctCovMat`).
- [x] **High-Performance C++ Backend (Rcpp Integration)**
  - C++ acceleration for prediction, likelihood evaluations, distance matrices, and Cholesky decompositions (`src/`).
- [x] **Efficient Multivariate Normal Sampling (`mvrnorm2`)**
  - Fast sampling routine for high-dimensional normal distributions.
- [x] **Evaluation Grid Utilities**
  - Domain grid generation tools (`matGrid`, `vecGrid`).


### Planned & Future Enhancements

- [ ] **MCMC Hyper-Parameter Sampling**: Bayesian inference via Markov Chain Monte Carlo for full posterior uncertainty over hyper-parameters.
- [ ] **Spatially Varying (Non-Stationary) Covariance Functions**: Support for spatially dynamic length-scales and non-stationary kernels.
- [ ] **Large-Scale Gaussian Process Approximations**: Scalable GP methods for large datasets (e.g., Sparse GPs, Inducing Point Methods, FITC, VFE).
- [ ] **Moving / Rolling Covariance Functions for Time Series**: Online/sequential covariance updates with exponential forgetting mechanisms for dynamic time-series models.

## Core Technical Implementation Details

| Feature / Utility | Module / Rcpp Source | Description |
| :--- | :--- | :--- |
| **Model Fitting** | `gpFit`, `gpNegLogLik` | Fits covariance parameters and mean model coefficients using numerical gradient optimization. |
| **Covariance Kernels** | `covm.cpp`, `covm` | Computes covariance matrices for Matérn (3/2, 5/2), Gaussian, Exponential, and Linear kernels. |
| **Likelihood Calculation** | `gpLogLik.cpp`, `gpLogLik` | Evaluates the log marginal likelihood of conditioned Gaussian processes. |
| **Simulation** | `gpSim` | Generates sample trajectories/realizations across user-defined spatial-temporal grids. |
| **Numerical Stability** | `cholfac`, `correctCovMat` | Handles ill-conditioned covariance matrices by adding diagonal jitter or correcting non-positive eigenvalues. |
| **Grid Generation** | `matGrid`, `vecGrid` | Constructs evaluation coordinate matrices for predictions and spatial-temporal visualization. |


This is an ongoing project. If you have any questions, requirements, suggestions, 
don't hesitate to contact me (in english, french or german):
<br/><a href="mailto:emanuel.huber@alumni.ethz.ch">emanuel.huber@alumni.ethz.ch</a></p>


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
