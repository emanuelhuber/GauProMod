// To run:  library(Rcpp); sourceCpp("backsolve_eigen_test.cpp")
// Backsolve triangular lower
// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::depends(RcppThread)]] // Include RcppThread for robust parallelization in R
// [[Rcpp::plugins(openmp)]]      // Enable OpenMP for Eigen's internal routines
#include <RcppEigen.h>
#include <RcppThread.h>  

using namespace Eigen;

// [[Rcpp::export]]
Eigen::MatrixXd cholupdateL_rcpp(const Eigen::Map<Eigen::MatrixXd>& L, 
                                 const Eigen::Map<Eigen::MatrixXd>& V12, 
                                 const Eigen::Map<Eigen::MatrixXd>& V22) { 
      
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


// [[Rcpp::export]]
Eigen::MatrixXd cholfac_rcpp(const Eigen::Map<Eigen::MatrixXd>& A){
  // Eigen::MatrixXd L = A.llt().matrixL();
  // return L;
  // A.llt().matrixL() computes the lower triangular factor L such that A = L*L.transpose()
  Eigen::LLT<Eigen::MatrixXd> lltOfA(A);
  if(lltOfA.info() != Success) {
    Rcpp::stop("Cholesky decomposition failed. Matrix M must be positive semi-definite.");
  }
  return lltOfA.matrixL();
}


//' Cross-distance between two matrices (RcppEigen version)
 //'
 //' Compute the Mahalanobis or Euclidean distance between every row of two matrices.
 //' @param X a matrix (or vector, handled as a matrix)
 //' @param Y a matrix (or vector) with the same number of columns as X
 //' @param M a positive semi-definite matrix for Mahalanobis distance, or NULL for Euclidean.
 //' @return A distance matrix of dimension nrow(X) x nrow(Y).
 // [[Rcpp::export]]
 Eigen::MatrixXd crossDist_rcpp(const Eigen::Map<Eigen::MatrixXd>& X,
                                const Eigen::Map<Eigen::MatrixXd>& Y,
                                const Rcpp::Nullable<Eigen::Map<Eigen::MatrixXd>>& M = R_NilValue) {
   
   // Check dimensions
   if (X.cols() != Y.cols()) {
     Rcpp::stop("X and Y must have the identical number of columns!");
   }
   
   // --- Case 1D (handled internally by Eigen if X is a vector/1-column matrix) ---
   // If the data is truly 1D (1 column), the R implementation uses outer(X, Y, "-").
   // The 2D (matrix) approach handles the 1D case efficiently as well:
   // D2 = (Xn + Yn.t()) - 2 * (X * Y.t())
   
   Eigen::MatrixXd X_mat = X;
   Eigen::MatrixXd Y_mat = Y;
   
   // --- Mahalanobis Transformation (If M is provided) ---
   if (M.isNotNull()) {
     const Eigen::Map<Eigen::MatrixXd>& M_mat = Rcpp::as<Eigen::Map<Eigen::MatrixXd>>(M);
     
     if (M_mat.rows() != X.cols() || M_mat.cols() != X.cols()) {
       Rcpp::stop("Matrix M must be square and have dimensions equal to the number of columns in X/Y.");
     }
     
     // Compute Cholesky factor L (M = L * L.transpose())
     Eigen::MatrixXd L = cholfac_rcpp(M_mat);
     
     // Apply transformation: Z = X %*% L (Mahalanobis distance)
     X_mat = X * L;
     Y_mat = Y * L;
   }
   
   // --- Compute Squared Euclidean Distance (General Case) ---
   
   // Xn = rowSums(X^2)
   Eigen::VectorXd Xn = X_mat.rowwise().squaredNorm();
   // Yn = rowSums(Y^2)
   Eigen::VectorXd Yn = Y_mat.rowwise().squaredNorm();
   
   // D2 = outer(Xn, Yn, "+") - 2 * tcrossprod(X, Y)
   // In Eigen: D2 = Xn * Eigen::RowVectorXd::Ones(Y_mat.rows()) + Eigen::VectorXd::Ones(X_mat.rows()) * Yn.transpose()
   // A much cleaner way using broadcasting:
   Eigen::MatrixXd D2 = Xn.replicate(1, Y_mat.rows()) + Yn.transpose().replicate(X_mat.rows(), 1);
   
   // D2 -= 2 * (X_mat * Y_mat.transpose())
   D2 -= 2.0 * (X_mat * Y_mat.transpose());
   
   // Numerical stability: D2[D2 < 0] <- 0
   for (int i = 0; i < D2.size(); ++i) {
     if (D2(i) < 0.0) {
       D2(i) = 0.0;
     }
   }
   
   // sqrt(D2)
   return D2.array().sqrt();
 }
