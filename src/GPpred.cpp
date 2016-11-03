// To run:  library(Rcpp); sourceCpp("backsolve_eigen_test.cpp")
// Backsolve triangular lower
// [[Rcpp::depends(RcppEigen)]]
//#include <Rcpp.h>
#include <RcppEigen.h>

// #define twoM_PI       6.2831853071795864769252867665590057683943388015061   /* pi */

using namespace Eigen;

// [[Rcpp::export]]

Rcpp::List GPpred_rcpp(const Eigen::Map<Eigen::MatrixXd>& K, 
                       const Eigen::Map<Eigen::MatrixXd>& Kstar, 
                       const Eigen::Map<Eigen::MatrixXd>& Kstarstar, 
                       const Eigen::Map<Eigen::MatrixXd>& y) { 
      
  int n = K.rows();
  int kk = Kstarstar.cols();

  Eigen::MatrixXd L(n, n);
  Eigen::VectorXd M(kk);	// Mean
  Eigen::MatrixXd C(kk, kk);	// covariance
  Eigen::MatrixXd bt(kk, kk);	
  Eigen::MatrixXd btb(kk, kk);	
  Eigen::VectorXd a(kk);
  Eigen::VectorXd logLik1(1);
  Eigen::VectorXd logLik2(n);

  // chol(K)	
  L =  K.adjoint().llt().matrixL();	
  bt = (L.triangularView<Lower>().solve(Kstar)).adjoint();
  a = L.triangularView<Lower>().solve(y);
  // Mean
  M = bt * a;
  btb = MatrixXd(kk,kk).setZero().selfadjointView<Lower>().rankUpdate(bt);
  // covariance
  C = Kstarstar - btb;
  // log-likelihood
  logLik1 = 0.5 * a.adjoint() * a;
  logLik2 = L.diagonal();
  return Rcpp::List::create(Rcpp::Named("mean") = M,
  					  		          Rcpp::Named("cov") = C,
						                Rcpp::Named("logLik1") = logLik1,
						                Rcpp::Named("logLik2") = logLik2);
}
