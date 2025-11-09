
// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::depends(RcppThread)]] // Include RcppThread for robust parallelization in R
// [[Rcpp::plugins(openmp)]]      // Enable OpenMP for Eigen's internal routines
#include <RcppEigen.h>
#include <RcppThread.h>          
#include <cmath> // For M_PI

// [[Rcpp::export]]
Rcpp::List GPpredmean_rcpp(const Eigen::Map<Eigen::MatrixXd>& K,
                           const Eigen::Map<Eigen::MatrixXd>& Kstar,
                           const Eigen::Map<Eigen::MatrixXd>& Kstarstar,
                           const Eigen::Map<Eigen::VectorXd>& y, // Correctly mapped as VectorXd
                           const Eigen::Map<Eigen::MatrixXd>& H,    // k x m
                           const Eigen::Map<Eigen::MatrixXd>& Hstar, // k x n
                           bool only_mean = false)
{
  
  // Dimensions
  int m = K.rows();             // m observations
  int n = Kstarstar.cols();     // n targets
  int k = H.rows();             // k basis functions (Note: H.cols() is m)
  
  // --- 1. Factorization K = L L^T (Highly Parallelized) ---
  // Use LLT object for efficiency and stability
  Eigen::LLT<Eigen::MatrixXd> lltOfK(K);
  if (lltOfK.info() != Eigen::Success) {
    Rcpp::stop("Cholesky decomposition of K failed. K must be positive definite.");
  }
  Eigen::MatrixXd L = lltOfK.matrixL();
  
  // --- 2. Compute Intermediate Components ---
  
  // d = L^-1 H^T. (d is m x k)
  // Solves L*d = H^T. Parallelized.
  Eigen::MatrixXd d = L.triangularView<Eigen::Lower>().solve(H.adjoint());
  
  // HKHT = H K^-1 H^T = d^T * d. (HKHT is k x k)
  // Matrix-Matrix multiplication, highly parallelized.
  Eigen::MatrixXd HKHT(k, k);
  HKHT.setZero().selfadjointView<Eigen::Lower>().rankUpdate(d.adjoint());
  
  // --- 3. Factorization HKHT = L2 L2^T (Parallelized) ---
  Eigen::LLT<Eigen::MatrixXd> lltOfHKHT(HKHT);
  if (lltOfHKHT.info() != Eigen::Success) {
    Rcpp::stop("Cholesky decomposition of HKHT failed. Matrix is singular.");
  }
  Eigen::MatrixXd L2 = lltOfHKHT.matrixL();
  
  // --- 4. Mean Parameter Estimation (B_hat) ---
  
  // a = L^-1 y. (a is m x 1)
  Eigen::VectorXd a = L.triangularView<Eigen::Lower>().solve(y); 
  
  // b = L^-1 Kstar. (b is m x n)
  // Multiple triangular solves, parallelized.
  Eigen::MatrixXd b = L.triangularView<Eigen::Lower>().solve(Kstar);
  
  // H K^-1 y = d^T * a. (k x 1)
  Eigen::VectorXd B_HK = d.adjoint() * a; 
  
  // LB = L2^-1 * (H K^-1 y). (LB is k x 1)
  // Solves L2 * LB = B_HK.
  Eigen::VectorXd LB = L2.triangularView<Eigen::Lower>().solve(B_HK); 
  
  // B = L2^-T * LB. (B is k x 1, which is beta_hat)
  Eigen::VectorXd B = L2.transpose().triangularView<Eigen::Upper>().solve(LB);
  
  // --- 5. Prediction Components ---
  
  // R = Hstar - d^T * b. (R is k x n)
  // This is Hstar - H K^-1 Kstar
  Eigen::MatrixXd R = Hstar - d.adjoint() * b;
  
  // Mean M = b^T * a + R^T * B. (M is n x 1)
  Eigen::VectorXd M = b.adjoint() * a + R.adjoint() * B;
  
  if (only_mean) {
    return Rcpp::List::create(Rcpp::Named("mean") = M);
  }
  
  // e = L2^-1 R. (e is k x n)
  // Solves L2 * e = R. Parallelized.
  Eigen::MatrixXd e = L2.triangularView<Eigen::Lower>().solve(R);
  
  // Covariance C = Kstarstar - b^T * b + e^T * e. (C is n x n)
  
  // b^T * b (btb)
  Eigen::MatrixXd btb(n, n);
  btb.setZero().selfadjointView<Eigen::Lower>().rankUpdate(b.adjoint());
  
  // e^T * e (ete)
  Eigen::MatrixXd ete(n, n);
  ete.setZero().selfadjointView<Eigen::Lower>().rankUpdate(e.adjoint());
  
  Eigen::MatrixXd C = Kstarstar - btb + ete;
  
  // --- 6. Corrected Log-Likelihood ---
  
  // Quadratic Term: -1/2 * ( ||a||^2 - ||LB||^2 ).
  // This correctly computes the quadratic term using the estimated mean.
  double logLik_T1 = -0.5 * (a.squaredNorm() - LB.squaredNorm());
  
  // Eigen::MatrixXd d2 =  L2.triangularView<Eigen::Lower>().solve(H);
  // Eigen::VectorXd w =	d2 * L.adjoint().triangularView<Eigen::Lower>().solve(a);
  // double logLik_T1 = -0.5 * (a.squaredNorm() - w.squaredNorm());
  // logLik1 -= 0.5 * w.adjoint()*w;
  
  // Determinant Term 1: -1/2 * log(|K|) = -sum(log(L_ii))
  double logLik_T2 = -L.diagonal().array().log().sum();
  
  // Determinant Term 2: -1/2 * log(|H K^-1 H^T|) = -sum(log(L2_ii))
  double logLik_T3 = -L2.diagonal().array().log().sum();
  
  // Constant Term: - (m-k)/2 * log(2*pi)
  double logLik_T4 = -0.5 * (m - k) * std::log(2.0 * M_PI);
  
  double logLik = logLik_T1 + logLik_T2 + logLik_T3 + logLik_T4;
  
  // Return results
  return Rcpp::List::create(Rcpp::Named("mean") = M,
                            Rcpp::Named("cov") = C,
                            Rcpp::Named("logLik") = logLik);
}
