// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::depends(RcppThread)]]
// [[Rcpp::plugins(openmp)]]
#include <RcppEigen.h>
#include <RcppThread.h>
#include <cmath>   // For M_PI
#include <limits>  // For -Inf

// ---------------------------------------------------------------------
// Lean log marginal likelihood evaluation for hyperparameter optimization.
//
// GPpred_rcpp() / GPpredmean_rcpp() compute the log marginal likelihood as
// a *byproduct* of prediction, so they require Kstar (obs x targ) and
// Kstarstar (targ x targ) even when only the likelihood value is wanted.
// During an optimization loop (e.g. inside optim()) those target-related
// covariance matrices are pure overhead: they can dominate the cost per
// iteration and are thrown away immediately. The functions below take
// only K (= Kxx, the observation covariance) and y (and H for the
// monomial-mean case), and return the scalar log-likelihood directly.
//
// Both return -Inf (rather than throwing) when the Cholesky decomposition
// fails, so a caller like optim() just sees a very poor objective value
// for that trial point instead of crashing the whole optimization run.
// ---------------------------------------------------------------------

//' Log marginal likelihood (zero/constant mean GP)
//'
//' Equivalent to the \code{logLik} value returned by \code{GPpred_rcpp},
//' but computed from \code{K} and \code{y} alone (no \code{Kstar} /
//' \code{Kstarstar} needed). Intended for use inside a hyperparameter
//' optimizer, where only the scalar likelihood matters and target
//' covariances would be wasted work.
//'
//' @param K observation covariance matrix (n x n), including the noise
//'   variance added to the diagonal (e.g. \code{Kxx + diag(sigma^2)}).
//' @param y observed values (length n).
//' @return the log marginal likelihood, or \code{-Inf} if \code{K} is not
//'   positive definite for the current hyperparameters.
// [[Rcpp::export]]
double gpLogLik_rcpp(const Eigen::Map<Eigen::MatrixXd>& K,
                      const Eigen::Map<Eigen::VectorXd>& y) {

  int n = K.rows();

  Eigen::LLT<Eigen::MatrixXd> lltOfK(K);
  if (lltOfK.info() != Eigen::Success) {
    // Not positive definite for these hyperparameters: signal a very
    // poor fit rather than stopping the optimizer.
    return -std::numeric_limits<double>::infinity();
  }
  Eigen::MatrixXd L = lltOfK.matrixL();

  // a = L^-1 y
  Eigen::VectorXd a = L.triangularView<Eigen::Lower>().solve(y);

  // logLik = -0.5 * ||a||^2 - sum(log(L_ii)) - (n/2) * log(2*pi)
  double logLik_T1 = -0.5 * a.squaredNorm();
  double logLik_T2 = -L.diagonal().array().log().sum();
  double logLik_T3 = -0.5 * n * std::log(2.0 * M_PI);

  return logLik_T1 + logLik_T2 + logLik_T3;
}

//' Log marginal likelihood (GP with monomial/basis mean function)
//'
//' Equivalent to the \code{logLik} value returned by \code{GPpredmean_rcpp}
//' (the REML-style marginal likelihood that integrates out the basis-
//' function coefficients), computed from \code{K}, \code{y} and \code{H}
//' alone -- no \code{Kstar} / \code{Kstarstar} / \code{Hstar} needed.
//'
//' @param K observation covariance matrix (m x m), noise already added.
//' @param y observed values (length m).
//' @param H basis-function design matrix (k x m), as built by \code{Hmat()}.
//' @return the log marginal likelihood, or \code{-Inf} if \code{K} or the
//'   basis-function system is not positive definite for the current
//'   hyperparameters.
// [[Rcpp::export]]
double gpLogLikMean_rcpp(const Eigen::Map<Eigen::MatrixXd>& K,
                          const Eigen::Map<Eigen::VectorXd>& y,
                          const Eigen::Map<Eigen::MatrixXd>& H) {

  int m = K.rows();
  int k = H.rows();

  Eigen::LLT<Eigen::MatrixXd> lltOfK(K);
  if (lltOfK.info() != Eigen::Success) {
    return -std::numeric_limits<double>::infinity();
  }
  Eigen::MatrixXd L = lltOfK.matrixL();

  // d = L^-1 H^T  (m x k)
  Eigen::MatrixXd d = L.triangularView<Eigen::Lower>().solve(H.adjoint());

  // HKHT = H K^-1 H^T = d^T d  (k x k)
  Eigen::MatrixXd HKHT(k, k);
  HKHT.setZero().selfadjointView<Eigen::Lower>().rankUpdate(d.adjoint());

  Eigen::LLT<Eigen::MatrixXd> lltOfHKHT(HKHT);
  if (lltOfHKHT.info() != Eigen::Success) {
    return -std::numeric_limits<double>::infinity();
  }
  Eigen::MatrixXd L2 = lltOfHKHT.matrixL();

  // a = L^-1 y
  Eigen::VectorXd a = L.triangularView<Eigen::Lower>().solve(y);

  // B_HK = H K^-1 y = d^T a  (k x 1)
  Eigen::VectorXd B_HK = d.adjoint() * a;

  // LB = L2^-1 * B_HK
  Eigen::VectorXd LB = L2.triangularView<Eigen::Lower>().solve(B_HK);

  double logLik_T1 = -0.5 * (a.squaredNorm() - LB.squaredNorm());
  double logLik_T2 = -L.diagonal().array().log().sum();
  double logLik_T3 = -L2.diagonal().array().log().sum();
  double logLik_T4 = -0.5 * (m - k) * std::log(2.0 * M_PI);

  return logLik_T1 + logLik_T2 + logLik_T3 + logLik_T4;
}
