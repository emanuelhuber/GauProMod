// To run:  library(Rcpp); sourceCpp("backsolve_eigen_test.cpp")
// Backsolve triangular lower
// [[Rcpp::depends(RcppEigen)]]
//#include <Rcpp.h>
#include <RcppEigen.h>

using namespace Eigen;

// [[Rcpp::export]]

Rcpp::List GPpredmean_rcpp(const Eigen::Map<Eigen::MatrixXd> K, 
                           const Eigen::Map<Eigen::MatrixXd> Kstar, 
                           const Eigen::Map<Eigen::MatrixXd> Kstarstar, 
                           const Eigen::Map<Eigen::VectorXd> y, 
                           const Eigen::Map<Eigen::MatrixXd> H, 
                           const Eigen::Map<Eigen::MatrixXd> Hstar) { 
      
  int m = K.rows();			// m observations
  int n = Kstarstar.cols();	// n targets
  int k = H.rows();			// k basis function (for mean function)
//  int nn = Hstar.cols();

  Eigen::MatrixXd L(m, m);		// chol(K)
  Eigen::MatrixXd d(m, k);
  Eigen::MatrixXd d2(k, m);
  Eigen::MatrixXd HKHT(k, k);	// H K^-1 H^T
  Eigen::MatrixXd L2(k, k);	// chol(HKHT) = chol( H K^-1 H^T )
  Eigen::VectorXd a(m);
  Eigen::VectorXd w(k);
  Eigen::MatrixXd b(m, n);
  Eigen::VectorXd LB(k);
  Eigen::VectorXd B(k);
  Eigen::MatrixXd R(k, n);
  Eigen::MatrixXd e(k, n);
  Eigen::MatrixXd btb(n, n);
  Eigen::MatrixXd ete(n, n);	
  Eigen::VectorXd M(n);		// Mean
  Eigen::MatrixXd C(n, n);	// covariance
  Eigen::VectorXd logLik1(1);
  Eigen::VectorXd logLik2(m);
  Eigen::VectorXd logLik3(k);
	
  
  // cholesky factor L such that K = LL^T	
  L =  K.adjoint().llt().matrixL();	
  // d=L^-1 H^T (forward substitiution)
  d = L.triangularView<Lower>().solve(H.adjoint());
  //dT = L.triangularView<Upper>().solve(H);
  // HKHT = H K^-1 H^T = d^T d
  HKHT = MatrixXd(k,k).setZero().selfadjointView<Lower>().rankUpdate(d.adjoint());
  // cholesky factor L2 such that HKHT = L2 L2^T
  L2 =  HKHT.adjoint().llt().matrixL();
  // 
  a = L.triangularView<Lower>().solve(y);
  //
  b = L.triangularView<Lower>().solve(Kstar);
	
  // Compute B = (H K⁻1 H^T)^-1 (H K^-1 y)
  LB = L2.triangularView<Lower>().solve(d.adjoint() * a);
  B = L2.transpose().triangularView<Upper>().solve(LB);	
  
  // Compute R
  R = Hstar - d.adjoint() * b;
	    
  //
  e = L2.triangularView<Lower>().solve(R);	

  // mean
  M = b.adjoint() * a + R.adjoint() * B;
	
  // Covariance
  btb = MatrixXd(n,n).setZero().selfadjointView<Lower>().rankUpdate(b.adjoint());
  ete = MatrixXd(n,n).setZero().selfadjointView<Lower>().rankUpdate(e.adjoint());
  C = Kstarstar - btb + ete;
	
  // log-likelihood
  // second term in the log-likelihood 
  d2 =  L2.triangularView<Lower>().solve(H);
  w =	d2 * L.adjoint().triangularView<Lower>().solve(a);
  logLik1 = 0.5 * a.adjoint()*a;
  logLik1 += 0.5 * w.adjoint()*w;
  logLik2 = L.diagonal();
  logLik3 = L2.diagonal();


  return Rcpp::List::create(Rcpp::Named("mean") = M,
  					        Rcpp::Named("cov") = C,
  					        Rcpp::Named("logLik1") = logLik1,
  					        Rcpp::Named("logLik2") = logLik2,
  					        Rcpp::Named("logLik3") = logLik3);
}

// http://lists.r-forge.r-project.org/pipermail/rcpp-devel/2014-June/007781.html
// Statistics & Software Consulting
// GKX Group, GKX Associates Inc.
// tel: 1-877-GKX-GROUP
// email: ggrothendieck at gmail.com


