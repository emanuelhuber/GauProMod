// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::depends(RcppThread)]] // Include RcppThread for parallel context
// [[Rcpp::plugins(openmp)]]      // Enable OpenMP compiler directives
#include <RcppEigen.h>
#include <RcppThread.h>          // Use RcppThread instead of just Rcpp for multi-threading
#include <cmath> // For M_PI

// K: K_observed_observed (K_oo), K*: K_observed_unobserved (K_ou), K**: K_unobserved_unobserved (K_uu)
// y: observed values (y_o)
// [[Rcpp::export]]
Rcpp::List GPpred_rcpp(const Eigen::Map<Eigen::MatrixXd>& K,
                       const Eigen::Map<Eigen::MatrixXd>& Kstar,
                       const Eigen::Map<Eigen::MatrixXd>& Kstarstar,
                       const Eigen::Map<Eigen::VectorXd>& y) {
  
  // IMPORTANT: Temporarily disable Eigen's multi-threading if running on a single core is required by Rcpp/R
  // However, for pure performance, we let Eigen manage the threads.
  // Eigen::setNbThreads(RcppThread::get=num_threads()); // Optional, to link to RcppThread settings
  
  int n = K.rows(); // Number of observed points
  // int k_star = Kstarstar.rows(); // Number of prediction points (not directly needed)
  
  // 1. Cholesky Decomposition of K: K = L L^T
  // This step is **automatically parallelized** by Eigen (via OpenMP/MKL if configured)
  Eigen::LLT<Eigen::MatrixXd> lltOfK(K);
  
  // Check if decomposition succeeded (for robustness, though less critical for performance)
  if (lltOfK.info() != Eigen::Success) {
    Rcpp::stop("Cholesky decomposition failed. Matrix K is not positive definite.");
  }
  
  Eigen::MatrixXd L = lltOfK.matrixL(); 
  
  // 2. Compute components using L.triangularView<Lower>().solve(X)
  
  // a = L^-1 * y (vector solve)
  // The triangular solve is efficient, but parallel speed-up is modest.
  Eigen::VectorXd a = L.triangularView<Eigen::Lower>().solve(y); 
  
  // V = L^-1 * Kstar (multiple triangular solves, parallelization is effective here)
  Eigen::MatrixXd V = L.triangularView<Eigen::Lower>().solve(Kstar);
  Eigen::MatrixXd beta_t = V.transpose(); 
  
  // --- Prediction Components ---
  
  // Mean M = beta_t * a 
  // Matrix-vector multiplication, automatically parallelized.
  Eigen::VectorXd M = beta_t * a;
  
  // VTV = V^T * V 
  // Matrix-Matrix multiplication, the most parallelized step in prediction.
  Eigen::MatrixXd VTV = V.transpose() * V; 
  
  // Covariance C = Kstarstar - VTV
  Eigen::MatrixXd C = Kstarstar - VTV;
  
  
  // --- Log-Likelihood Component ---
  // The following steps involve vector operations (e.g., squaredNorm, sum) 
  // which are also automatically parallelized by Eigen's array operations.
  
  // log_likelihood = -0.5 * ||a||^2 - sum(log(L_ii)) - (n/2) * log(2*pi)
  
  // Term 1: -0.5 * a^T * a
  double logLik_T1 = -0.5 * a.squaredNorm(); 
  // Term 2: - log(|L|) = - sum(log(L_ii))
  double logLik_T2 = -L.diagonal().array().log().sum();
  // Term 3: Constant term
  double logLik_T3 = -0.5 * n * std::log(2.0 * M_PI);
  
  double logLik = logLik_T1 + logLik_T2 + logLik_T3;
  
  // Return results
  return Rcpp::List::create(Rcpp::Named("mean") = M,
                            Rcpp::Named("cov") = C,
                            Rcpp::Named("logLik") = logLik);
}