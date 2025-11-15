// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::depends(RcppThread)]]
// [[Rcpp::plugins(openmp)]]
#include <RcppEigen.h>
#include <RcppThread.h>
#include <RcppBessel.h>
#include <cmath>
#include <limits>
#include <string>
#include <stdexcept>
#include <Rmath.h> // R::bessel_k, R::gammafn

// Inline helpers
inline double sqr(double x) { return x * x; }

// -----------------------------------------------------------------------------
// --- Existing Matern scalar (unchanged) -------------------------------------
// -----------------------------------------------------------------------------
double kMatern_scalar(double r, double l, double h, double v, int d, double w) {
  if (l <= 0 || h < 0) return 0.0;
  
  if (r < std::sqrt(std::numeric_limits<double>::epsilon()) * l) {
    if (d == 0) return h * h;
    if (d == 1) return 0.0;
    if (d == 2) {
      if (v <= 1.0) Rcpp::stop("Second derivative limit is infinite for Matern nu <= 1.");
      double k0 = (h * h / (l * l)) * (1.0 / (2.0 * (v - 1.0)));
      return k0;
    }
  }
  
  double u_t = r / l;
  double gamV = R::gammafn(v);
  double h_sq = h * h;
  double l_sq = l * l;
  
  if (d == 0) {
    if (std::abs(v - 0.5) < 1e-9) {
      return h_sq * std::exp(-u_t);
    } else if (std::abs(v - 1.5) < 1e-9) {
      double u_sqrt3 = std::sqrt(3.0) * u_t;
      return h_sq * (1.0 + u_sqrt3) * std::exp(-u_sqrt3);
    } else if (std::abs(v - 2.5) < 1e-9) {
      double u_sqrt5 = std::sqrt(5.0) * u_t;
      double r_sq = r * r;
      return h_sq * (1.0 + u_sqrt5 + (5.0 * r_sq) / (3.0 * l_sq)) * std::exp(-u_sqrt5);
    }
  }
  
  double common_factor_base = std::pow(2.0, 1.0 - v) / gamV;
  
  if (d == 0) {
    double besK_v = R::bessel_k(u_t, v, 1);
    double K = common_factor_base * std::pow(u_t, v) * besK_v;
    return K * h_sq;
    
  } else if (d == 1) {
    double besK_v_minus_1 = R::bessel_k(u_t, v - 1.0, 1);
    double common = h * h * common_factor_base * std::pow(u_t, v);
    double K = w * (1.0 / l) * common * besK_v_minus_1;
    return K;
    
  } else if (d == 2) {
    double k1_base = common_factor_base * std::pow(u_t, v - 2.0);
    double besK_v_minus_1 = R::bessel_k(u_t, v - 1.0, 1);
    double besK_v_minus_2 = R::bessel_k(u_t, v - 2.0, 1);
    
    double k1 = h * h * k1_base;
    double u_t_sq = u_t * u_t;
    double k2 = (u_t * besK_v_minus_1) - (u_t_sq * w * besK_v_minus_2);
    double K = (1.0 / (l * l)) * k1 * k2;
    
    return K;
  }
  
  Rcpp::stop("Unsupported derivative order 'd' for Matern kernel.");
  return 0.0;
}

// -----------------------------------------------------------------------------
// --- Vectorized Cauchy, Linear (distance), Spherical kernels -----------------
// -----------------------------------------------------------------------------

// Vectorized Cauchy kernel (works for rectangular R/W)
Eigen::MatrixXd kCauchy_rcpp_fast(const Eigen::Map<Eigen::MatrixXd> &R,
                                  double l, double h, double nu, int d,
                                  const Eigen::Map<Eigen::MatrixXd> &W,
                                  bool use_symmetry) {
  const int n_rows = R.rows();
  const int n_cols = R.cols();
  if (W.rows() != n_rows || W.cols() != n_cols) Rcpp::stop("R and W must have identical dimensions.");
  if (use_symmetry && (n_rows != n_cols)) Rcpp::stop("use_symmetry==true requires square R/W matrices.");
  if (l <= 0 || h < 0) return Eigen::MatrixXd::Zero(n_rows, n_cols);
  
  const double h_sq = h * h;
  const double l_sq = l * l;
  
  Eigen::ArrayXXd r_arr = R.array();
  Eigen::ArrayXXd w_arr = W.array();
  Eigen::ArrayXXd u = (r_arr.square()) / l_sq;         // (r/l)^2
  Eigen::ArrayXXd base = (1.0 + u);                   // 1 + (r/l)^2
  
  if (d == 0) {
    Eigen::ArrayXXd arr = h_sq * base.pow(-nu);
    Eigen::MatrixXd K = arr.matrix();
    if (use_symmetry) K = (K + K.transpose()) * 0.5;
    return K;
  } else if (d == 1) {
    // return w * (-dK/dr) = w * h^2 * (2*nu*r / l^2) * base^{-nu-1}
    Eigen::ArrayXXd arr = w_arr * (h_sq * (2.0 * nu) * r_arr / l_sq) * base.pow(-nu - 1.0);
    Eigen::MatrixXd K = arr.matrix();
    if (use_symmetry) K = (K + K.transpose()) * 0.5;
    return K;
  } else if (d == 2) {
    // d2K/dr2 derived earlier: compute elementwise then return -d2
    // C = 2 nu h^2 / l^2
    double C = (2.0 * nu * h_sq) / l_sq;
    Eigen::ArrayXXd term1 = -C * base.pow(-nu - 1.0);
    Eigen::ArrayXXd term2 = C * ((2.0 * r_arr.square()) / l_sq) * (nu + 1.0) * base.pow(-nu - 2.0);
    Eigen::ArrayXXd d2 = term1 + term2;
    Eigen::ArrayXXd arr = -d2; // return -d2K/dr2
    Eigen::MatrixXd K = arr.matrix();
    if (use_symmetry) K = (K + K.transpose()) * 0.5;
    return K;
  }
  
  Rcpp::stop("Unsupported derivative order 'd' for Cauchy vectorized kernel.");
  return Eigen::MatrixXd::Zero(n_rows, n_cols);
}

// Vectorized distance-based triangular kernel
Eigen::MatrixXd kTriangular_rcpp_fast(const Eigen::Map<Eigen::MatrixXd> &R,
                                          double l, double h, int d,
                                          const Eigen::Map<Eigen::MatrixXd> &W,
                                          bool use_symmetry) {
  const int n_rows = R.rows();
  const int n_cols = R.cols();
  if (W.rows() != n_rows || W.cols() != n_cols) Rcpp::stop("R and W must have identical dimensions.");
  if (use_symmetry && (n_rows != n_cols)) Rcpp::stop("use_symmetry==true requires square R/W matrices.");
  if (l <= 0 || h < 0) return Eigen::MatrixXd::Zero(n_rows, n_cols);
  
  const double h_sq = h * h;
  Eigen::ArrayXXd r_arr = R.array();
  Eigen::ArrayXXd w_arr = W.array();
  
  // mask for inside support r < l
  Eigen::ArrayXXd mask = (r_arr < l).cast<double>();
  
  if (d == 0) {
    Eigen::ArrayXXd arr = h_sq * (1.0 - r_arr / l);
    arr = mask * arr; // zero outside support
    Eigen::MatrixXd K = arr.matrix();
    if (use_symmetry) K = (K + K.transpose()) * 0.5;
    return K;
  } else if (d == 1) {
    // -dK/dr = h^2 / l, so return w * (h^2 / l) inside support
    Eigen::ArrayXXd arr = w_arr * (h_sq / l);
    arr = mask * arr;
    Eigen::MatrixXd K = arr.matrix();
    if (use_symmetry) K = (K + K.transpose()) * 0.5;
    return K;
  } else if (d == 2) {
    // -d2K/dr2 = 0 everywhere (piecewise linear)
    Eigen::ArrayXXd arr = Eigen::ArrayXXd::Zero(n_rows, n_cols);
    Eigen::MatrixXd K = arr.matrix();
    return K;
  }
  
  Rcpp::stop("Unsupported derivative order 'd' for Triangular distance kernel.");
  return Eigen::MatrixXd::Zero(n_rows, n_cols);
}

// Vectorized spherical kernel
Eigen::MatrixXd kSpherical_rcpp_fast(const Eigen::Map<Eigen::MatrixXd> &R,
                                     double l, double h, int d,
                                     const Eigen::Map<Eigen::MatrixXd> &W,
                                     bool use_symmetry) {
  const int n_rows = R.rows();
  const int n_cols = R.cols();
  if (W.rows() != n_rows || W.cols() != n_cols) Rcpp::stop("R and W must have identical dimensions.");
  if (use_symmetry && (n_rows != n_cols)) Rcpp::stop("use_symmetry==true requires square R/W matrices.");
  if (l <= 0 || h < 0) return Eigen::MatrixXd::Zero(n_rows, n_cols);
  
  const double h_sq = h * h;
  Eigen::ArrayXXd r_arr = R.array();
  Eigen::ArrayXXd w_arr = W.array();
  Eigen::ArrayXXd x = r_arr / l;
  Eigen::ArrayXXd mask = (r_arr < l).cast<double>();
  
  if (d == 0) {
    Eigen::ArrayXXd arr = h_sq * (1.0 - 1.5 * x + 0.5 * x.pow(3));
    arr = mask * arr;
    Eigen::MatrixXd K = arr.matrix();
    if (use_symmetry) K = (K + K.transpose()) * 0.5;
    return K;
  } else if (d == 1) {
    // -dK/dr = h^2*(1.5/l)*(1 - x^2)
    Eigen::ArrayXXd arr = w_arr * (h_sq * (1.5 / l) * (1.0 - x.square()));
    arr = mask * arr;
    Eigen::MatrixXd K = arr.matrix();
    if (use_symmetry) K = (K + K.transpose()) * 0.5;
    return K;
  } else if (d == 2) {
    // -d2K/dr2 = - h^2 * (3/(l*l)) * x
    Eigen::ArrayXXd arr = - (h_sq * (3.0 / (l * l)) * x);
    arr = mask * arr;
    Eigen::MatrixXd K = arr.matrix();
    if (use_symmetry) K = (K + K.transpose()) * 0.5;
    return K;
  }
  
  Rcpp::stop("Unsupported derivative order 'd' for spherical kernel.");
  return Eigen::MatrixXd::Zero(n_rows, n_cols);
}

// -----------------------------------------------------------------------------
// --- Existing Gaussian fast path (unchanged) --------------------------------
// -----------------------------------------------------------------------------
Eigen::MatrixXd kGaussian_rcpp_fast(const Eigen::Map<Eigen::MatrixXd> &R,
                                    double l, double h, int d,
                                    const Eigen::Map<Eigen::MatrixXd> &W,
                                    bool use_symmetry) {
  const int n_rows = R.rows();
  const int n_cols = R.cols();
  
  if (W.rows() != n_rows || W.cols() != n_cols) {
    Rcpp::stop("R and W must have identical dimensions.");
  }
  if (use_symmetry && (n_rows != n_cols)) {
    Rcpp::stop("use_symmetry==true requires square R/W matrices.");
  }
  
  const double h_sq = h * h;
  if (l <= 0 || h < 0) {
    return Eigen::MatrixXd::Zero(n_rows, n_cols);
  }
  
  const double inv_l = 1.0 / l;
  const double inv_l_sq = inv_l * inv_l;
  
  if (d == 0) {
    Eigen::ArrayXXd r_div_l = R.array() * inv_l;
    Eigen::ArrayXXd u = -0.5 * r_div_l.square();
    Eigen::ArrayXXd arr = u.exp() * h_sq;
    Eigen::MatrixXd K = arr.matrix();
    if (use_symmetry) K = (K + K.transpose()) * 0.5;
    return K;
  }
  
  Eigen::ArrayXXd r_arr = R.array();
  Eigen::ArrayXXd w_arr = W.array();
  Eigen::ArrayXXd r_div_l = r_arr * inv_l;
  Eigen::ArrayXXd base = ( (-0.5) * r_div_l.square() ).exp();
  
  if (d == 1) {
    Eigen::ArrayXXd Karr = (w_arr * r_arr * inv_l_sq) * base * h_sq;
    Eigen::MatrixXd K = Karr.matrix();
    if (use_symmetry) K = (K + K.transpose()) * 0.5;
    return K;
  }
  
  if (d == 2) {
    Eigen::ArrayXXd Karr = ( (1.0 - w_arr * r_div_l.square()) * inv_l_sq ) * base * h_sq;
    Eigen::MatrixXd K = Karr.matrix();
    if (use_symmetry) K = (K + K.transpose()) * 0.5;
    return K;
  }
  
  Rcpp::stop("Unsupported derivative order 'd' in Gaussian fast path.");
  return Eigen::MatrixXd::Zero(n_rows, n_cols);
}

//' Polynomial (Gram) Kernel with Derivatives
 //'
 //' Computes the polynomial kernel K(X, Y) = h^2 * (X %*% t(Y) + c)^degree
 //' and optionally supports first and second derivatives with respect to input entries.
 //'
 //' @param X Numeric matrix (n x p).
 //' @param Y Numeric matrix (m x p).
 //' @param h Scaling factor for the kernel.
 //' @param degree Integer degree of the polynomial (>=1).
 //' @param c Bias term added inside the polynomial.
 //' @param d Derivative order: 0 = kernel value, 1 = first derivative, 2 = second derivative.
 //' @param W Weight matrix (n x m), used to scale derivatives elementwise.
 //' @param use_symmetry Logical, if TRUE enforces symmetry (only valid if X and Y have same number of rows).
 //' @return Numeric matrix (n x m) containing kernel values or derivatives.
 Eigen::MatrixXd kPolynomial_rcpp_fast(const Eigen::Map<Eigen::MatrixXd> X,
                                  const Eigen::Map<Eigen::MatrixXd> Y,
                                  double h, int degree, double c,
                                  int d,
                                  const Eigen::Map<Eigen::MatrixXd> W,
                                  bool use_symmetry) {
   if(degree < 1) Rcpp::stop("Polynomial degree must be >= 1.");
   
   int n = X.rows();
   int m = Y.rows();
   
   if (W.rows() != n || W.cols() != m)
     Rcpp::stop("Weight matrix W must have same dimensions as output kernel.");
   
   Eigen::MatrixXd K(n, m);
   
   // Compute Gram matrix with bias
   Eigen::ArrayXXd base = (X * Y.transpose()).array() + c;
   
   if(d == 0){
     K = (h * h) * base.pow(degree).matrix();
     if(use_symmetry){
       if(n != m) Rcpp::stop("use_symmetry requires X and Y to have same number of rows");
       K = (K + K.transpose()) * 0.5;
     }
   } else if(d == 1){
     // First derivative: h^2 * degree * (base)^(degree-1) * W
     Eigen::ArrayXXd arr = W.array() * (h * h * degree) * base.pow(degree - 1);
     K = arr.matrix();
     if(use_symmetry){
       if(n != m) Rcpp::stop("use_symmetry requires X and Y to have same number of rows");
       K = (K + K.transpose()) * 0.5;
     }
   } else if(d == 2){
     // Second derivative: h^2 * degree * (degree-1) * (base)^(degree-2) * W
     if(degree < 2){
       // degree 1 polynomial => second derivative zero
       K = Eigen::MatrixXd::Zero(n, m);
     } else {
       Eigen::ArrayXXd arr = W.array() * (h * h * degree * (degree - 1)) * base.pow(degree - 2);
       K = arr.matrix();
     }
     if(use_symmetry){
       if(n != m) Rcpp::stop("use_symmetry requires X and Y to have same number of rows");
       K = (K + K.transpose()) * 0.5;
     }
   } else {
     Rcpp::stop("Unsupported derivative order 'd'. Must be 0, 1, or 2.");
   }
   
   return K;
 }

//' Linear (Gram) Kernel with Derivatives
 //'
 //' Computes the linear kernel K(X, Y) = h^2 * (X %*% t(Y)) + b^2
 //' and optionally supports first and second derivatives with respect to input entries.
 //'
 //' @param X Numeric matrix (n x p).
 //' @param Y Numeric matrix (m x p).
 //' @param h Scaling factor for the inner product term.
 //' @param b Bias term added to the kernel.
 //' @param d Derivative order: 0 = kernel value, 1 = first derivative, 2 = second derivative.
 //' @param w Weight matrix (n x m), used to scale derivatives elementwise.
 //' @param use_symmetry Logical, if TRUE enforces symmetry (only valid if X and Y have same number of rows).
 //' @return Numeric matrix (n x m) containing kernel values or derivatives.
 Eigen::MatrixXd kLinear_rcpp_fast(const Eigen::Map<Eigen::MatrixXd> X,
                              const Eigen::Map<Eigen::MatrixXd> Y,
                              double h, 
                              double b,
                              int d,
                              const Eigen::Map<Eigen::MatrixXd> W,
                              bool use_symmetry) {
   
   int n = X.rows();
   int m = Y.rows();
   
   if (W.rows() != n || W.cols() != m)
     Rcpp::stop("Weight matrix W must have same dimensions as output kernel.");
   
   Eigen::MatrixXd K(n, m);
   
   if(d == 0){
     // Kernel value: h^2 * X Y^T + b^2
     K.noalias() = h * h * (X * Y.transpose());
     K.array() += b * b;
     // Symmetry enforcement
     if(use_symmetry){
       if(n != m) Rcpp::stop("use_symmetry requires X and Y to have same number of rows");
       K = (K + K.transpose()) * 0.5;
     }
   } else if(d == 1){
     // First derivative w.r.t each entry is constant h^2
     // Multiply elementwise by W
     K = W.array() * (h * h);
     if(use_symmetry){
       if(n != m) Rcpp::stop("use_symmetry requires X and Y to have same number of rows");
       K = (K + K.transpose()) * 0.5;
     }
   } else if(d == 2){
     // Second derivative is zero
     K = Eigen::MatrixXd::Zero(n, m);
   } else {
     Rcpp::stop("Unsupported derivative order 'd'. Must be 0, 1, or 2.");
   }
   
   return K;
 }

// -----------------------------------------------------------------------------
// --- Threaded general path for non-vectorized kernels (keeps Matern etc.) ---
// -----------------------------------------------------------------------------
Eigen::MatrixXd kernel_matrix_threaded_general(
    const Eigen::Map<Eigen::MatrixXd> &R, double l, double h, double v, int d,
    const Eigen::Map<Eigen::MatrixXd> &W, const std::string &kernel, bool use_symmetry) {
  
  const int n_rows = R.rows();
  const int n_cols = R.cols();
  
  if (W.rows() != n_rows || W.cols() != n_cols)
    Rcpp::stop("R and W must have identical dimensions.");
  
  if (use_symmetry && (n_rows != n_cols))
    Rcpp::stop("use_symmetry==true requires square R/W matrices.");
  
  Eigen::MatrixXd K = Eigen::MatrixXd::Zero(n_rows, n_cols);
  
  const double inv_l = 1.0 / l;
  const double inv_l_sq = inv_l * inv_l;
  const double h_sq = h * h;
  
  #pragma omp parallel for schedule(static)
  for (int i = 0; i < n_rows; ++i) {
    int jstart = 0;
    if (use_symmetry) jstart = i;
    for (int j = jstart; j < n_cols; ++j) {
      double r = R(i, j);
      double w = W(i, j);
      double Kval = 0.0;
      
      if (kernel == "matern") {
        Kval = kMatern_scalar(r, l, h, v, d, w);
      } else {
        // safety fallback (should not be hit for Gaussian/Cauchy/Triangular/Spherical)
        #pragma omp critical
        {
          Rcpp::stop("Unsupported kernel in threaded general path.");
        }
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
    const Eigen::Map<Eigen::MatrixXd> X,
    const Eigen::Map<Eigen::MatrixXd> Y,
    double l, 
    double h, 
    double v, 
    int degree,
    double c, 
    int d, 
    const Eigen::Map<Eigen::MatrixXd> &W, 
    const std::string &kernel, 
    bool use_symmetry) {
  
  if(X.rows() != W.rows() || X.cols() != W.cols())
    Rcpp::stop("Weight matrix W must match output dimensions.");
  
  if(use_symmetry && (X.rows() != Y.rows()))
    Rcpp::stop("use_symmetry requires X and Y to have same number of rows.");

  
  // Fast vectorized paths
  if (kernel == "gaussian") {
    return kGaussian_rcpp_fast(X, l, h, d, W, use_symmetry);
  } else if (kernel == "cauchy") {
    return kCauchy_rcpp_fast(X, l, h, v, d, W, use_symmetry);
  } else if (kernel == "triangular") {
    return kTriangular_rcpp_fast(X, l, h, d, W, use_symmetry);
  } else if (kernel == "spherical") {
    return kSpherical_rcpp_fast(X, l, h, d, W, use_symmetry);
  }else if(kernel == "linear"){
    return kLinear_rcpp_fast(X, Y, h, c, d, W, use_symmetry);
  } else if(kernel == "polynomial"){
    return kPolynomial_rcpp_fast(X, Y, h, degree, c, d, W, use_symmetry);
  } else if(kernel == "matern"){
    // Fallback threaded path (Matern etc.)
    return kernel_matrix_threaded_general(X, l, h, v, d, W, kernel, use_symmetry);
  }
}



// -----------------------------------------------------------------------------
// --- Exported wrappers -------------------------------------------------------
// -----------------------------------------------------------------------------

// ------------------- Gaussian Kernel -------------------
//' Gaussian Kernel Matrix
//'
//' Computes the Gaussian (squared exponential) kernel matrix for a given distance matrix.
//'
//' @param R Distance matrix (n x m). Must be non-negative.
//' @param l Length scale parameter (>0).
//' @param h Marginal standard deviation (scale factor).
//' @param d Derivative order: 0 = kernel, 1 = first derivative, 2 = second derivative.
//' @param W Weight matrix (n x m) to scale derivatives elementwise.
//' @param use_symmetry Logical; if TRUE, enforces symmetry (only valid if R is square).
//' @return Kernel matrix (n x m) or its derivative.
//' @examples
//' R <- as.matrix(dist(matrix(rnorm(10*2), 10, 2)))
//' W <- matrix(1, nrow(R), ncol(R))
//' kGaussian_rcpp(R, l=1, h=1, d=0, W, use_symmetry=TRUE)
// [[Rcpp::export]]
Eigen::MatrixXd kGaussian_rcpp(const Eigen::Map<Eigen::MatrixXd> R,
                              double l,
                              double h,
                              int d,
                              const Eigen::Map<Eigen::MatrixXd> W,
                              bool use_symmetry) {
 return kernel_matrix_rcpp(R, R, l, h, 0.0, 0, 0.0, d, W, "gaussian", use_symmetry);
}

// ------------------- Matern Kernel -------------------
//' Matern Kernel Matrix
//'
//' Computes the Matern kernel matrix for a given distance matrix.
//'
//' @param R Distance matrix (n x m).
//' @param l Length scale parameter (>0).
//' @param h Marginal standard deviation (scale factor).
//' @param v Smoothness parameter (nu > 0).
//' @param d Derivative order: 0 = kernel, 1 = first derivative, 2 = second derivative.
//' @param W Weight matrix (n x m) to scale derivatives elementwise.
//' @param use_symmetry Logical; if TRUE, enforces symmetry (only valid if R is square).
//' @return Kernel matrix (n x m) or its derivative.
// [[Rcpp::export]]
Eigen::MatrixXd kMatern_rcpp(const Eigen::Map<Eigen::MatrixXd> R,
                            double l,
                            double h,
                            double v,
                            int d,
                            const Eigen::Map<Eigen::MatrixXd> W,
                            bool use_symmetry) {
 return kernel_matrix_rcpp(R, R, l, h, v, 0, 0.0, d, W, "matern", use_symmetry);
}

// ------------------- Triangular Kernel -------------------
//' Triangular Kernel Matrix
//'
//' Computes the triangular kernel matrix (distance-based).
//'
//' @param R Distance matrix (n x m).
//' @param l Range/length scale (>0).
//' @param h Scale factor.
//' @param d Derivative order: 0 = kernel, 1 = first derivative, 2 = second derivative.
//' @param W Weight matrix (n x m) to scale derivatives.
//' @param use_symmetry Logical; if TRUE, enforces symmetry.
//' @return Kernel matrix (n x m) or its derivative.
// [[Rcpp::export]]
Eigen::MatrixXd kTriangular_rcpp(const Eigen::Map<Eigen::MatrixXd> R,
                                double l,
                                double h,
                                int d,
                                const Eigen::Map<Eigen::MatrixXd> W,
                                bool use_symmetry) {
 return kernel_matrix_rcpp(R, R, l, h, 0.0, 0, 0.0, d, W, "triangular", use_symmetry);
}

// ------------------- Spherical Kernel -------------------
//' Spherical Kernel Matrix
//'
//' Computes the spherical kernel matrix (distance-based).
//'
//' @param R Distance matrix (n x m).
//' @param l Range/length scale (>0).
//' @param h Scale factor.
//' @param d Derivative order: 0 = kernel, 1 = first derivative, 2 = second derivative.
//' @param W Weight matrix (n x m) to scale derivatives.
//' @param use_symmetry Logical; if TRUE, enforces symmetry.
//' @return Kernel matrix (n x m) or its derivative.
// [[Rcpp::export]]
Eigen::MatrixXd kSpherical_rcpp(const Eigen::Map<Eigen::MatrixXd> R,
                               double l,
                               double h,
                               int d,
                               const Eigen::Map<Eigen::MatrixXd> W,
                               bool use_symmetry) {
 return kernel_matrix_rcpp(R, R, l, h, 0.0, 0, 0.0, d, W, "spherical", use_symmetry);
}

// ------------------- Cauchy Kernel -------------------
//' Cauchy Kernel Matrix
//'
//' Computes the Cauchy kernel matrix (distance-based).
//'
//' @param R Distance matrix (n x m).
//' @param l Length scale (>0).
//' @param h Scale factor.
//' @param v Shape parameter (nu > 0).
//' @param d Derivative order: 0 = kernel, 1 = first derivative, 2 = second derivative.
//' @param W Weight matrix (n x m) to scale derivatives.
//' @param use_symmetry Logical; if TRUE, enforces symmetry.
//' @return Kernel matrix (n x m) or its derivative.
// [[Rcpp::export]]
Eigen::MatrixXd kCauchy_rcpp(const Eigen::Map<Eigen::MatrixXd> R,
                            double l,
                            double h,
                            double v,
                            int d,
                            const Eigen::Map<Eigen::MatrixXd> W,
                            bool use_symmetry) {
 return kernel_matrix_rcpp(R, R, l, h, v, 0, 0.0, d, W, "cauchy", use_symmetry);
}

// ------------------- Linear Gram Kernel -------------------
//' Linear (Gram) Kernel Matrix
//'
//' Computes the linear kernel matrix: K = h^2 * X %*% t(Y) + c^2.
//'
//' @param X Feature matrix (n x p).
//' @param Y Feature matrix (m x p).
//' @param h Scale factor for the inner product.
//' @param c Bias term added to the kernel.
//' @param d Derivative order: 0 = kernel, 1 = first derivative, 2 = second derivative.
//' @param W Weight matrix (n x m) to scale derivatives.
//' @param use_symmetry Logical; if TRUE, enforces symmetry (requires X and Y have same rows).
//' @return Kernel matrix (n x m) or its derivative.
// [[Rcpp::export]]
Eigen::MatrixXd kLinear_rcpp(const Eigen::Map<Eigen::MatrixXd> X,
                            const Eigen::Map<Eigen::MatrixXd> Y,
                            double h,
                            double c,
                            int d,
                            const Eigen::Map<Eigen::MatrixXd> W,
                            bool use_symmetry) {
 return kernel_matrix_rcpp(X, Y, 0.0, h, 0.0, 0, c, d, W, "linear", use_symmetry);
}

// ------------------- Polynomial Gram Kernel -------------------
//' Polynomial (Gram) Kernel Matrix
//'
//' Computes the polynomial kernel: K = h^2 * (X %*% t(Y) + c)^degree.
//'
//' @param X Feature matrix (n x p).
//' @param Y Feature matrix (m x p).
//' @param h Scale factor.
//' @param degree Polynomial degree (integer >=1).
//' @param c Bias term added inside the polynomial.
//' @param d Derivative order: 0 = kernel, 1 = first derivative, 2 = second derivative.
//' @param W Weight matrix (n x m) to scale derivatives.
//' @param use_symmetry Logical; if TRUE, enforces symmetry (requires X and Y have same rows).
//' @return Kernel matrix (n x m) or its derivative.
// [[Rcpp::export]]
Eigen::MatrixXd kPolynomial_rcpp(const Eigen::Map<Eigen::MatrixXd> X,
                                const Eigen::Map<Eigen::MatrixXd> Y,
                                double h,
                                int degree,
                                double c,
                                int d,
                                const Eigen::Map<Eigen::MatrixXd> W,
                                bool use_symmetry) {
 return kernel_matrix_rcpp(X, Y, 0.0, h, 0.0, degree, c, d, W, "polynomial", use_symmetry);
}
