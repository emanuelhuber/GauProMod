// To run:  library(Rcpp); sourceCpp("backsolve_eigen_test.cpp")
// Backsolve triangular lower
// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>

using namespace Eigen;

// [[Rcpp::export]]
Eigen::MatrixXd cholupdateL_rcpp(const Eigen::Map<Eigen::MatrixXd> L, 
                                 const Eigen::Map<Eigen::MatrixXd> V12, 
                                 const Eigen::Map<Eigen::MatrixXd> V22) { 
      
  int k = L.rows();
  int k2 = V22.rows();
  Eigen::MatrixXd S(k, k2);
  Eigen::MatrixXd U(k2, k2);
  Eigen::MatrixXd M(k2, k2);
  Eigen::MatrixXd Lup(k+k2, k+k2);
  
  Lup.setZero();
  S = L.triangularView<Lower>().solve(V12);
  M = V22 -  S.adjoint() * S ;
  U =  M.adjoint().llt().matrixL();
  Lup.topLeftCorner(k,k) = L;
  Lup.bottomLeftCorner(k2,k) = S.adjoint();
  Lup.bottomRightCorner(k2,k2) = U;
  return Lup;
}



