// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::plugins(openmp)]]
#include <RcppEigen.h>
#include <RcppThread.h>
#include <RcppBessel.h>
#include <cmath>
#include <limits>
#include <string>
#include <stdexcept>
#include <Rmath.h> // R::bessel_k, R::gammafn
#include <Eigen/Core>
#include <omp.h>

using namespace Rcpp;
using namespace Eigen;

// --------------------- Helper Functions ---------------------
inline double sqr(double x) { return x * x; }
const double TOL = std::sqrt(std::numeric_limits<double>::epsilon());

// ---------------------------------------------------------------------
// RAII helper: temporarily set Eigen to single-threaded while OpenMP
// parallel regions run, then restore previous thread count.
// This avoids nested parallelism (Eigen internal threads + OpenMP).
// ---------------------------------------------------------------------
struct EigenThreadGuard {
  int prev_threads;
  EigenThreadGuard() {
    prev_threads = Eigen::nbThreads();
    if (prev_threads != 1) Eigen::setNbThreads(1);
  }
  ~EigenThreadGuard() {
    if (Eigen::nbThreads() != prev_threads)
      Eigen::setNbThreads(prev_threads);
  }
};

// Kernel type enum for faster dispatch
enum class KernelType { GAUSSIAN, MATERN, CAUCHY, TRIANGULAR, SPHERICAL, LINEAR, POLYNOMIAL };

inline KernelType kernel_string_to_enum(const std::string& kernel) {
  if (kernel == "gaussian") return KernelType::GAUSSIAN;
  if (kernel == "matern") return KernelType::MATERN;
  if (kernel == "cauchy") return KernelType::CAUCHY;
  if (kernel == "triangular") return KernelType::TRIANGULAR;
  if (kernel == "spherical") return KernelType::SPHERICAL;
  if (kernel == "linear") return KernelType::LINEAR;
  if (kernel == "polynomial") return KernelType::POLYNOMIAL;
  Rcpp::stop("Unknown kernel type: " + kernel);
  return KernelType::GAUSSIAN; // unreachable
}

// --------------------- Matern Scalar Kernel ---------------------
double kMatern_scalar(double r, double l, double h, double v, int d, double w) {
  if (l <= 0 || h < 0) return 0.0;
  
  if (r < TOL * l) {
    if (d==0) return h*h;
    if (d==1) return 0.0;
    if (d==2) return (v <= 1.0) ? R_PosInf : (h*h/(l*l))*(1.0/(2.0*(v-1.0)));
  }
  
  double u = r / l;
  double h_sq = h*h;
  double l_sq = l*l;
  const double gamV = R::gammafn(v);
  const double two_pow_1_minus_v = std::pow(2.0, 1.0 - v); // micro-opt: precompute
  
  if(d==0) {
    if (std::abs(v-0.5)<1e-9) return h_sq * std::exp(-u);
    if (std::abs(v-1.5)<1e-9) return h_sq * (1+std::sqrt(3.0)*u) * std::exp(-std::sqrt(3.0)*u);
    if (std::abs(v-2.5)<1e-9) return h_sq * (1 + std::sqrt(5.0)*u + 5.0*r*r/(3*l_sq)) * std::exp(-std::sqrt(5.0)*u);
    double besK = R::bessel_k(u,v,1);
    return h_sq * two_pow_1_minus_v/gamV * std::pow(u,v) * besK;
  }
  else if(d==1) {
    double besK = R::bessel_k(u,v-1.0,1);
    return w*h_sq*two_pow_1_minus_v/gamV*std::pow(u,v)*besK / l;
  }
  else if(d==2) {
    double besK1 = R::bessel_k(u,v-1.0,1);
    double besK2 = R::bessel_k(u,v-2.0,1);
    double k1 = h_sq * two_pow_1_minus_v/gamV*std::pow(u,v-2.0);
    return k1 * (u*besK1 - u*u*w*besK2) / (l*l);
  }
  Rcpp::stop("Unsupported derivative order for Matern kernel.");
  return 0.0;
}

// --------------------- Dense Kernels (vectorized with tolerance fix) ---------------------
Eigen::MatrixXd kGaussian_rcpp_fast(const Eigen::MatrixXd &R, double l, double h, int d,
                                    const Eigen::MatrixXd &W, bool use_symmetry) {
  const double h_sq = h*h;
  const double inv_l_sq = 1.0 / (l * l);
  
  if(!use_symmetry) {
 
    Eigen::ArrayXXd r_arr = R.array();
    Eigen::ArrayXXd w_arr = W.array();
    // Compute u = -0.5 * (r / l)^2
    Eigen::ArrayXXd r_div_l = r_arr / l;
    Eigen::ArrayXXd u = -0.5 * r_div_l.square();
    Eigen::ArrayXXd K = (r_arr < TOL).select(1, u.exp());
    
    if(d==1){
      K = (w_arr * r_arr * inv_l_sq) * K;
    } else if(d==2){
      K = (1.0 - r_div_l.square()) * inv_l_sq * K;
    } 
    
    K = K * h_sq;
    return K.matrix();
  }
  
  const int n = R.rows();
  if(n != R.cols()) Rcpp::stop("Symmetric kernel requires square R");
  
  Eigen::MatrixXd K = Eigen::MatrixXd::Zero(n, n);
  // --- Compute only upper triangle using triangularView ---
  for(int i=0; i<n; ++i) {
    int len = n - i;
    Eigen::ArrayXd R_row = R.block(i, i, 1, len).transpose().array();
    Eigen::ArrayXd W_row = W.block(i, i, 1, len).transpose().array();
    Eigen::ArrayXd K_row = (R_row < TOL).select(1, (-0.5 * R_row.square() * inv_l_sq).exp());

    if(d==1) {
      // K_row = h_sq * W_row * R_row * inv_l_sq * (-0.5 * R_row.square() * inv_l_sq).exp();
      K_row = W_row * R_row * inv_l_sq * K_row;
    }
    else if(d==2) {
      // K_row = h_sq * (1.0 - R_row.square() * inv_l_sq) * (-0.5 * R_row.square() * inv_l_sq).exp();
      K_row = (1.0 - R_row.square() * inv_l_sq) * inv_l_sq * K_row;
    }
    else if(d>2){
      Rcpp::stop("Unsupported derivative order");
    }
    
    K_row = K_row * h_sq;
    // Fill upper triangle
    K.block(i,i,1,len) = K_row.transpose();
    // Mirror to lower triangle
    if(len>1) K.block(i+1,i,len-1,1) = K_row.tail(len-1);
  }
  return K;
  // if(use_symmetry) K = 0.5*(K + K.transpose());
}

Eigen::MatrixXd kCauchy_rcpp_fast(const Eigen::MatrixXd &R, double l, double h, double nu, 
                                  int d,
                                  const Eigen::MatrixXd &W, bool use_symmetry) {
  const double h_sq = h*h;
  const double l_sq = l*l;
  Eigen::ArrayXXd r_arr = R.array();
  Eigen::ArrayXXd w_arr = W.array();
  Eigen::ArrayXXd base = 1.0 + r_arr.square()/l_sq;
  Eigen::ArrayXXd K;
  
  if(d==0){
    // K = base.pow(-nu);
    K = (r_arr < TOL).select(1, base.pow(-nu));
  } 
  else if(d==1) {
    // K = w_arr*2.0*nu*h_sq*r_arr/l_sq*base.pow(-nu-1.0);
    K = (r_arr < TOL).select(1, base.pow(-nu-1.0));
    K = w_arr*2.0*nu*r_arr/l_sq* K;
  }
  else if(d==2){
    Eigen::ArrayXXd K2;
    K = (r_arr < TOL).select(1, base.pow(-nu-2.0));
    K2 = (r_arr < TOL).select(1, base.pow(-nu-1.0));
    // K = w_arr*2.0*nu*h_sq/l_sq*((nu+1.0)*2.0*r_arr.square()/l_sq*base.pow(-nu-2.0) - base.pow(-nu-1.0));
    K = w_arr*2.0*nu/l_sq*((nu+1.0)*2.0*r_arr.square()/l_sq*K - K2);
  } 
  
  K = K * h_sq;
  // if(use_symmetry) K = 0.5*(K + K.transpose());
  return K.matrix();
}

Eigen::MatrixXd kTriangular_rcpp_fast(const Eigen::MatrixXd &R,double l,double h,int d,
                                      const Eigen::MatrixXd &W,bool use_symmetry) {
  const double h_sq=h*h;
  Eigen::ArrayXXd r_arr = R.array();
  Eigen::ArrayXXd w_arr = W.array();
  Eigen::ArrayXXd mask = (r_arr < l).cast<double>();
  Eigen::ArrayXXd K;
  
  if(d==0) K = mask*(h_sq*(1.0 - r_arr/l));
  else if(d==1) K = mask*(w_arr*h_sq/l);
  else K = Eigen::ArrayXXd::Zero(R.rows(), R.cols());
  
  if(d==0) K = (r_arr < TOL).select(h_sq, K);
  else if(d>0) K = (r_arr < TOL).select(0.0, K);
  
  // if(use_symmetry) K = 0.5*(K + K.transpose());
  return K.matrix();
}

Eigen::MatrixXd kSpherical_rcpp_fast(const Eigen::MatrixXd &R,double l,double h,int d,
                                     const Eigen::MatrixXd &W,bool use_symmetry) {
  const double h_sq=h*h;
  Eigen::ArrayXXd r_arr = R.array();
  Eigen::ArrayXXd w_arr = W.array();
  Eigen::ArrayXXd x = r_arr/l;
  Eigen::ArrayXXd mask = (r_arr < l).cast<double>();
  Eigen::ArrayXXd K;
  
  if(d==0) K = mask*(h_sq*(1.0 - 1.5*x + 0.5*x.pow(3)));
  else if(d==1) K = mask*(w_arr*h_sq*1.5/l*(1.0 - x.square()));
  else K = mask*(-h_sq*3.0/(l*l)*x);
  
  if(d==0) K = (r_arr < TOL).select(h_sq, K);
  else K = (r_arr < TOL).select(0.0, K);
  
  // if(use_symmetry) K = 0.5*(K + K.transpose());
  return K.matrix();
}

// --------------------- Gram/Dense Linear & Polynomial ---------------------
Eigen::MatrixXd kLinear_rcpp_fast(const Eigen::MatrixXd &X,const Eigen::MatrixXd &Y,
                                  double h,double c,int d,
                                  const Eigen::MatrixXd &W,bool use_symmetry) {
  Eigen::ArrayXXd base = (X*Y.transpose()).array();
  Eigen::ArrayXXd K;
  
  if(d==0) K = h*h*base + c*c;
  else if(d==1) K = W.array()*h*h;
  else K = Eigen::ArrayXXd::Zero(X.rows(), Y.rows());
  
  // if(use_symmetry) K = 0.5*(K + K.transpose());
  return K.matrix();
}

Eigen::MatrixXd kPolynomial_rcpp_fast(const Eigen::MatrixXd &X,
                                      const Eigen::MatrixXd &Y,
                                      double h,
                                      int degree,
                                      double c,
                                      int d,
                                      const Eigen::MatrixXd &W,
                                      bool use_symmetry) {
  Eigen::ArrayXXd base = (X * Y.transpose()).array() + c;
  Eigen::ArrayXXd K(X.rows(), Y.rows());
  
  if(d == 0) {
    K = h*h * base.pow(degree);
  }
  else if(d == 1) {
    K = W.array() * h*h * degree * base.pow(degree - 1);
  }
  else if(d == 2) {
    // Use select instead of ?: for Eigen compatibility
    Eigen::ArrayXXd K_zero = Eigen::ArrayXXd::Zero(X.rows(), Y.rows());
    if(degree < 2) {
      K = K_zero;
    } else {
      K = W.array() * h*h * degree * (degree - 1) * base.pow(degree - 2);
    }
  }
  else {
    Rcpp::stop("Unsupported derivative order 'd'.");
  }
  // if(use_symmetry) K = 0.5 * (K + K.transpose());
  return K.matrix();
}

// --------------------- Sparse Kernels with tolerance check ---------------------
template <typename KernelFun>
Eigen::SparseMatrix<double> sparse_kernel_generic(const Eigen::SparseMatrix<double>& R,
                                                  const Eigen::SparseMatrix<double>& W,
                                                  KernelFun fun) {
  std::vector<Eigen::Triplet<double>> triplets;
  triplets.reserve(R.nonZeros());
  
  for(int k=0;k<R.outerSize();++k){
    for(Eigen::SparseMatrix<double>::InnerIterator it(R,k); it; ++it){
      double r = it.value();
      double w = W.coeff(it.row(), it.col());
      if(r < TOL) r = 0.0;
      double val = fun(r,w);
      if(val != 0.0) triplets.emplace_back(it.row(), it.col(), val);
    }
  }
  Eigen::SparseMatrix<double> K(R.rows(),R.cols());
  K.setFromTriplets(triplets.begin(),triplets.end());
  return K;
}

// Computes the Matern kernel for a sparse distance matrix.
// 
// @param R Sparse distance matrix (n x m), non-negative entries.
// @param l Length scale (>0).
// @param h Marginal standard deviation (>0).
// @param v Smoothness parameter (nu > 0).
// @param d Derivative order: 0, 1, 2.
// @param W Sparse or dense weight matrix, same dimensions as R. If dense, it is converted internally.
// @param use_symmetry If true, enforces symmetry (only works for square R/W).
// @return SparseMatrix<double> containing the kernel values (or derivatives).
Eigen::SparseMatrix<double> kMatern_sparse(const Eigen::SparseMatrix<double> &R,
                                           double l,
                                           double h,
                                           double v,
                                           int d,
                                           const Eigen::SparseMatrix<double> &W,
                                           bool use_symmetry = false) {
  
  const int n_rows = R.rows();
  const int n_cols = R.cols();
  
  if (W.rows() != n_rows || W.cols() != n_cols)
    Rcpp::stop("Weight matrix W must have same dimensions as sparse R.");
  
  if (use_symmetry && (n_rows != n_cols))
    Rcpp::stop("use_symmetry requires square matrices.");
  
  std::vector<Triplet<double>> triplets;
  triplets.reserve(R.nonZeros());

  // Ensure Eigen uses a single thread inside our OpenMP parallel region
  EigenThreadGuard eigen_guard;
  // Optional OpenMP parallelization
  #pragma omp parallel
  {
    std::vector<Triplet<double>> local_triplets;
    // Heuristic reserve per-thread to reduce reallocations
    // local_triplets.reserve(std::max(16, R.nonZeros() / (std::max(1, omp_get_max_threads()))));
    local_triplets.reserve(std::max(Eigen::Index(16), R.nonZeros() / (std::max(Eigen::Index(1), Eigen::Index(omp_get_max_threads())))));
    
    #pragma omp for nowait
    for (int col = 0; col < R.outerSize(); ++col) {
      for (SparseMatrix<double>::InnerIterator it(R, col); it; ++it) {
        const int i = it.row();
        const int j = it.col();
        const double r = it.value();
        if (r < 0.0) Rcpp::stop("Distance matrix must be non-negative.");

        double w_ij = W.coeff(i, j);
        double Kval = kMatern_scalar(r, l, h, v, d, w_ij);

        if (Kval != 0.0) local_triplets.emplace_back(i, j, Kval);
      }
    }

    #pragma omp critical
    triplets.insert(triplets.end(), local_triplets.begin(), local_triplets.end());
  }

  SparseMatrix<double> Ksp(n_rows, n_cols);
  if (!triplets.empty()) Ksp.setFromTriplets(triplets.begin(), triplets.end());

  if (use_symmetry) {
    // enforce symmetry by copying upper/lower triangle
    SparseMatrix<double> Ksym = Ksp;
    for (int k = 0; k < Ksp.outerSize(); ++k) {
      for (SparseMatrix<double>::InnerIterator it(Ksp, k); it; ++it) {
        int i = it.row();
        int j = it.col();
        if (i != j) Ksym.coeffRef(j, i) = it.value();
      }
    }
    return Ksym;
  }
  return Ksp;
}

// --------------------- Sparse Gaussian ---------------------
SparseMatrix<double> kGaussian_sparse(const SparseMatrix<double>& R, double l, double h, int d,
                                      const SparseMatrix<double>& W) {
  double h_sq=h*h, inv_l_sq=1.0/(l*l);
  return sparse_kernel_generic(R,W,[=](double r,double w){
    if(d==0) return h_sq*std::exp(-0.5*r*r*inv_l_sq);
    if(d==1) return w*h_sq*r*inv_l_sq*std::exp(-0.5*r*r*inv_l_sq);
    if(d==2) return w*h_sq*(1.0-r*r*inv_l_sq)*inv_l_sq*std::exp(-0.5*r*r*inv_l_sq);
    return 0.0;
  });
}

// --------------------- Sparse Cauchy ---------------------
SparseMatrix<double> kCauchy_sparse(const SparseMatrix<double>& R,double l,double h,double nu,int d,
                                    const SparseMatrix<double>& W){
  double h_sq=h*h,l_sq=l*l;
  return sparse_kernel_generic(R,W,[=](double r,double w){
    double u=r*r/l_sq, base=1.0+u;
    if(d==0) return h_sq*std::pow(base,-nu);
    if(d==1) return w*2.0*nu*h_sq*r/l_sq*std::pow(base,-nu-1.0);
    if(d==2){ double C=2.0*nu*h_sq/l_sq; return -( -C*std::pow(base,-nu-1.0) + C*(2*r*r/l_sq)*(nu+1.0)*std::pow(base,-nu-2.0)); }
    return 0.0;
  });
}

// --------------------- Sparse Triangular ---------------------
SparseMatrix<double> kTriangular_sparse(const SparseMatrix<double>& R,double l,double h,int d,
                                        const SparseMatrix<double>& W){
  double h_sq=h*h;
  return sparse_kernel_generic(R,W,[=](double r,double w){
    if(r>=l) return 0.0;
    if(d==0) return h_sq*(1.0 - r/l);
    if(d==1) return w*h_sq/l;
    if(d==2) return 0.0;
    return 0.0;
  });
}

// --------------------- Sparse Spherical ---------------------
SparseMatrix<double> kSpherical_sparse(const SparseMatrix<double>& R,double l,double h,int d,
                                       const SparseMatrix<double>& W){
  double h_sq=h*h;
  return sparse_kernel_generic(R,W,[=](double r,double w){
    if(r>=l) return 0.0;
    double x=r/l;
    if(d==0) return h_sq*(1.0-1.5*x+0.5*x*x*x);
    if(d==1) return w*h_sq*1.5/l*(1.0-x*x);
    if(d==2) return -h_sq*3.0/(l*l)*x;
    return 0.0;
  });
}

// --------------------- Sparse Linear ---------------------
SparseMatrix<double> kLinear_sparse(const SparseMatrix<double>& X,const SparseMatrix<double>& Y,
                                    double h,double c,int d,const SparseMatrix<double>& W,bool use_symmetry){
  SparseMatrix<double> K(X.rows(),Y.cols());
  std::vector<Triplet<double>> triplets;
  triplets.reserve(W.nonZeros());
  for(int k=0;k<W.outerSize();++k){
    for(SparseMatrix<double>::InnerIterator it(W,k);it; ++it){
      int i=it.row(), j=it.col(); double w=it.value();
      double val=0.0;
      if(d==0) val=h*h*(X.row(i)*Y.row(j).transpose()).sum() + c*c;
      else if(d==1) val=w*h*h;
      else if(d==2) val=0.0;
      triplets.emplace_back(i,j,val);
    }
  }
  K.setFromTriplets(triplets.begin(),triplets.end());
  if(use_symmetry) K=(K+SparseMatrix<double>(K.transpose()))*0.5;
  return K;
}

// --------------------- Sparse Polynomial ---------------------
SparseMatrix<double> kPolynomial_sparse(const SparseMatrix<double>& X,const SparseMatrix<double>& Y,
                                        double h,int degree,double c,int d,const SparseMatrix<double>& W,bool use_symmetry){
  SparseMatrix<double> K(X.rows(),Y.cols());
  std::vector<Triplet<double>> triplets;
  triplets.reserve(W.nonZeros());
  for(int k=0;k<W.outerSize();++k){
    for(SparseMatrix<double>::InnerIterator it(W,k);it;++it){
      int i=it.row(), j=it.col(); double w=it.value();
      double base=(X.row(i)*Y.row(j).transpose()).sum()+c;
      double val=0.0;
      if(d==0) val=h*h*std::pow(base,degree);
      else if(d==1) val=w*h*h*degree*std::pow(base,degree-1);
      else if(d==2) val=(degree<2)?0.0:w*h*h*degree*(degree-1)*std::pow(base,degree-2);
      triplets.emplace_back(i,j,val);
    }
  }
  K.setFromTriplets(triplets.begin(),triplets.end());
  if(use_symmetry) K=(K+SparseMatrix<double>(K.transpose()))*0.5;
  return K;
}

// -----------------------------------------------------------------------------
// --- Threaded general path for non-vectorized kernels (keeps Matern etc.) ---
// -----------------------------------------------------------------------------
Eigen::MatrixXd kernel_matrix_threaded_general(
    const MatrixXd &R, double l, double h, double v, int d,
    const MatrixXd &W, KernelType kernel_type, bool use_symmetry) {
  
  const int n_rows = R.rows();
  const int n_cols = R.cols();
  
  if (W.rows() != n_rows || W.cols() != n_cols)
    Rcpp::stop("R and W must have identical dimensions.");
  
  if (use_symmetry && (n_rows != n_cols))
    Rcpp::stop("use_symmetry==true requires square R/W matrices.");
  
  Eigen::MatrixXd K = Eigen::MatrixXd::Zero(n_rows, n_cols);
  
  // Ensure Eigen uses a single thread inside OpenMP-parallel region
  EigenThreadGuard eigen_guard;
  #pragma omp parallel for schedule(static)
  for (int i = 0; i < n_rows; ++i) {
    int jstart = 0;
    if (use_symmetry) jstart = i;
    for (int j = jstart; j < n_cols; ++j) {
      double r = R(i, j);
      double w = W(i, j);
      double Kval = 0.0;
      
      switch(kernel_type) {
        case KernelType::MATERN:
              Kval = kMatern_scalar(r, l, h, v, d, w);
              break;
          default:
              #pragma omp critical
              { Rcpp::stop("Unsupported kernel in threaded general path."); }
      }
      K(i, j) = Kval;
    }
  }
  
  if (use_symmetry) {
    for (int i = 0; i < n_rows; ++i) {
      for (int j = i + 1; j < n_cols; ++j) {
        K(j, i) = K(i, j);
      }
    }
  }
  return K;
}

// -----------------------------------------------------------------------------
// --- Kernel matrix wrapper: dispatch to vectorized fast paths where possible
// -----------------------------------------------------------------------------

//' Kernel Dispatcher: Gaussian, Matern, Cauchy, Linear, Polynomial, Spherical
//'
//' Computes a kernel matrix for various kernels with optional derivatives.
//'
//' @param X Numeric matrix (n x p) or distance matrix (n x m) depending on kernel.
//' @param Y Numeric matrix (m x p) or distance matrix (n x m) depending on kernel.
//' @param l Length scale (used for Gaussian, Matern, Cauchy, etc.).
//' @param h Marginal standard deviation / scaling.
//' @param v Parameter for Matern (nu) or Cauchy (nu) kernels.
//' @param degree Degree of polynomial kernel (for polynomial kernel).
//' @param c Bias term (for linear/polynomial kernels).
//' @param d Derivative order: 0 = kernel, 1 = first derivative, 2 = second derivative.
//' @param W Weight matrix (n x m), used to scale derivatives.
//' @param kernel Kernel type: "gaussian", "matern", "cauchy", "linear", "polynomial", "spherical".
//' @param use_symmetry Logical, if TRUE enforces symmetry (only valid for square X/Y).
//' @return Kernel matrix (n x m) or derivatives as matrix.
Eigen::MatrixXd kernel_matrix_rcpp(
   const Eigen::MatrixXd  &X,
   const Eigen::MatrixXd  &Y,
   double l, 
   double h, 
   double v, 
   int degree,
   double c, 
   int d, 
   const Eigen::MatrixXd &W, 
   const std::string &kernel, 
   bool use_symmetry) {
 
 if(X.rows() != W.rows() || X.cols() != W.cols())
   Rcpp::stop("Weight matrix W must match output dimensions.");
 
 if(use_symmetry && (X.rows() != Y.rows()))
   Rcpp::stop("use_symmetry requires X and Y to have same number of rows.");

 // Convert kernel string to enum for faster dispatch
  KernelType ktype = kernel_string_to_enum(kernel);
  switch(ktype) {
      case KernelType::GAUSSIAN: return kGaussian_rcpp_fast(X, l, h, d, W, use_symmetry);
      case KernelType::CAUCHY: return kCauchy_rcpp_fast(X, l, h, v, d, W, use_symmetry);
      case KernelType::TRIANGULAR: return kTriangular_rcpp_fast(X, l, h, d, W, use_symmetry);
      case KernelType::SPHERICAL: return kSpherical_rcpp_fast(X, l, h, d, W, use_symmetry);
      case KernelType::LINEAR: return kLinear_rcpp_fast(X, Y, h, c, d, W, use_symmetry);
      case KernelType::POLYNOMIAL: return kPolynomial_rcpp_fast(X, Y, h, degree, c, d, W, use_symmetry);
      case KernelType::MATERN: {
        // KernelType ktype = kernel_string_to_enum(kernel); // inside braces
        return kernel_matrix_threaded_general(X, l, h, v, d, W, ktype, use_symmetry);
      }
      default:
        Rcpp::stop("Unsupported kernel in threaded general path.");
      }
}

// Forward declaration
SEXP kernel_dispatch_auto_rcpp(SEXP X_s,
                               SEXP Y_s,
                               double l,
                               double h,
                               double v,
                               int degree,
                               double c,
                               int d,
                               SEXP W_s,
                               std::string kernel,
                               bool use_symmetry);

//' Kernel Dispatcher (Dense + Sparse Support)
//'
//' Computes kernel matrices for a variety of kernels, supporting both dense feature matrices
//' and sparse distance matrices (dgCMatrix). Automatically dispatches to optimized implementations
//' depending on the kernel type and input format.
//'
//' Supported kernels:
//' * Distance-based kernels (sparse/dense): `"gaussian"`, `"matern"`, `"cauchy"`, `"triangular"`, `"spherical"`
//' * Feature/Gram kernels (dense only): `"linear"`, `"polynomial"`
//'
//' Derivatives are supported via the `d` argument:
//' * `d = 0`: Kernel values
//' * `d = 1`: First derivative
//' * `d = 2`: Second derivative
//'
//' Symmetry enforcement (`use_symmetry = TRUE`) works for square distance matrices only.
//'
//' @param X Dense numeric matrix or sparse dgCMatrix (distance matrix for distance kernels)
//' @param Y Dense numeric matrix or sparse dgCMatrix (same dimensions as X)
//' @param l Length scale (>0) for distance kernels
//' @param h Scale factor / marginal standard deviation
//' @param v Smoothness parameter (`nu`) for Matern or shape for Cauchy
//' @param degree Polynomial degree (integer >=1) for polynomial kernel
//' @param c Bias term for linear/polynomial kernels
//' @param d Derivative order: 0, 1, or 2
//' @param W Dense or sparse weight matrix, same dimensions as X/Y, scales derivatives elementwise
//' @param kernel Kernel type (string): `"gaussian"`, `"matern"`, `"cauchy"`, `"triangular"`, `"spherical"`, `"linear"`, `"polynomial"`
//' @param use_symmetry Logical; if TRUE, enforces symmetry (only valid for square X/Y distance matrices)
//' @return Kernel matrix (dense `MatrixXd` if inputs are dense, sparse `dgCMatrix` if inputs are sparse)
//' @examples
//' # Dense Gaussian kernel
//' X <- matrix(rnorm(20), 5, 4)
//' W <- matrix(1, nrow(X), nrow(X))
//' k <- kernel_dispatch_auto_rcpp(X, X, l=1, h=1, v=0, degree=0, c=0, d=0, W, "gaussian", TRUE)
//'
//' # Sparse distance-based Matern kernel
//' library(Matrix)
//' R <- as(Matrix(dist(matrix(rnorm(25),5,5))), "dgCMatrix")
//' Wsp <- Matrix(1,5,5,sparse=TRUE)
//' Ksp <- kernel_dispatch_auto_rcpp(R, R, l=1, h=1, v=1.5, degree=0, c=0, d=0, Wsp, "matern", TRUE)
// [[Rcpp::export]]
SEXP kernel_dispatch_auto_rcpp(SEXP X_s,
                              SEXP Y_s,
                              double l,
                              double h,
                              double v,
                              int degree,
                              double c,
                              int d,
                              SEXP W_s,
                              std::string kernel,
                              bool use_symmetry = false) {
 
 if (h < 0) Rcpp::stop("Scale 'h' must be >= 0.");
 if (d < 0 || d > 2) Rcpp::stop("Derivative order 'd' must be 0, 1 or 2.");
 
 auto is_dgCMatrix = [](SEXP s)->bool {
   return Rf_isS4(s) && Rcpp::RObject(s).inherits("dgCMatrix");
 };
 
 bool is_distance_kernel = (kernel == "gaussian" || kernel == "matern" ||
                            kernel == "cauchy" || kernel == "triangular" ||
                            kernel == "spherical");
 
 // Feature kernels cannot use sparse input
 if (!is_distance_kernel && (is_dgCMatrix(X_s) || is_dgCMatrix(Y_s))) {
   Rcpp::stop("Feature/gram kernels (linear/polynomial) require dense feature matrices X and Y.");
 }
 
 // ---------------- Dense path ----------------
 if (!is_dgCMatrix(X_s) && !is_dgCMatrix(Y_s)) {
   Eigen::MatrixXd X = Rcpp::as<Eigen::MatrixXd>(X_s);
   Eigen::MatrixXd Y = Rcpp::as<Eigen::MatrixXd>(Y_s);
   Eigen::MatrixXd W = Rcpp::as<Eigen::MatrixXd>(W_s);
   
   Eigen::MatrixXd K = kernel_matrix_rcpp(X, Y, l, h, v, degree, c, d, W, kernel, use_symmetry);
   return Rcpp::wrap(K);
 }
 
 // ---------------- Sparse path ----------------
 if (!is_distance_kernel) {
   Rcpp::stop("Sparse inputs supported only for distance-based kernels.");
 }
 
 Eigen::SparseMatrix<double> Xsp = Rcpp::as<Eigen::SparseMatrix<double>>(X_s);
 Eigen::SparseMatrix<double> Ysp = Rcpp::as<Eigen::SparseMatrix<double>>(Y_s);
 
 if (Xsp.rows() != Ysp.rows() || Xsp.cols() != Ysp.cols()) {
   Rcpp::stop("Sparse X and Y must have identical dimensions.");
 }
 
 const int n_rows = Xsp.rows();
 const int n_cols = Xsp.cols();
 bool W_is_sparse = is_dgCMatrix(W_s);
 Eigen::SparseMatrix<double> Wsp;
 Eigen::MatrixXd Wdense;
 if (W_is_sparse) {
   Wsp = Rcpp::as<Eigen::SparseMatrix<double>>(W_s);
   if (Wsp.rows() != n_rows || Wsp.cols() != n_cols)
     Rcpp::stop("Sparse W must match dimensions of X/Y.");
 } else {
   Wdense = Rcpp::as<Eigen::MatrixXd>(W_s);
   if (Wdense.rows() != n_rows || Wdense.cols() != n_cols)
     Rcpp::stop("Weight matrix W must match dimensions of X/Y.");
 }
 
 // Convert dense W to sparse if needed
 Eigen::SparseMatrix<double> Wsp_converted = W_is_sparse ? Wsp : Wdense.sparseView();
 
 // Dispatch kernels
 if (kernel == "gaussian") {
   Eigen::SparseMatrix<double> Ksp = kGaussian_sparse(Xsp, l, h, d, Wsp_converted);
   return Rcpp::wrap(Ksp);
 } else if (kernel == "cauchy") {
   Eigen::SparseMatrix<double> Ksp = kCauchy_sparse(Xsp, l, h, v, d, Wsp_converted);
   return Rcpp::wrap(Ksp);
 } else if (kernel == "triangular") {
   Eigen::SparseMatrix<double> Ksp = kTriangular_sparse(Xsp, l, h, d, Wsp_converted);
   return Rcpp::wrap(Ksp);
 } else if (kernel == "spherical") {
   Eigen::SparseMatrix<double> Ksp = kSpherical_sparse(Xsp, l, h, d, Wsp_converted);
   return Rcpp::wrap(Ksp);
 } else if (kernel == "matern") {
   Eigen::SparseMatrix<double> Ksp = kMatern_sparse(Xsp, l, h, v, d, Wsp_converted, use_symmetry);
   return Rcpp::wrap(Ksp);
 } else {
   Rcpp::stop("Unsupported kernel for sparse dispatch: " + kernel);
 }
 return R_NilValue;
}

 
