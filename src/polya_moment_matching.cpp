// -----------------------------------------------------------------------------
// Random draws from a Polya proposal.
// Note: The parameters are set via moment matching
// -----------------------------------------------------------------------------

#include <RcppEigen.h>

// [[Rcpp::depends(RcppEigen)]]

#include <iostream>
#include <Eigen/Dense>
#include <Eigen/IterativeLinearSolvers>
#include <cmath>
#include "polya_moment_matching.h"


/**
 * This function optimizes $f(x) = \left(2 \frac{\psi(x) - \psi(c-x)}{2x-c}    \
 * -\lambda\right)^2$, where $\psi()$ denotes the digamma function. The        \
 * optimization is performed using Newton's method. Such an optimum is needed  \
 * to find the parameters of the Polya proposal via moment matching as in      \
 * https://github.com/changwoo-lee/logisticbeta-reproduce/                     \
 *                                                                             \
 * @param c double.                                                            \
 * @param l double.                                                            \
 * @param sum_e1T float.                                                       \
 *                                                                             \ 
 * @return Eigen::VectorXf: Minimizer of the function f(x).                    \
 */
// [[Rcpp::export]]
double optimize_moment_matching(double c, double l) {
  
  // Optimization settings -----------------------------------------------------
  const double tol = 1e-08; 
  const double eps = 1e-05;
  const int maxit = 100;
  
  // Midpoint check ------------------------------------------------------------
  double q_mid = 2.0 * R::trigamma(c / 2.0);
  if (l <= q_mid + tol) return c / 2.0;
  
  // Initialization of Newton's routine ----------------------------------------
  double x = c / 4.0;
  
  // Main Newton's loop --------------------------------------------------------
  for (int i = 0; i < maxit; i++) {
    // Constants
    double d = 2.0 * x - c;
    double psi_x = R::digamma(x);
    double psi_cx = R::digamma(c - x);
    double tri_x = R::trigamma(x);
    double tri_cx = R::trigamma(c - x);
    double delta = psi_x - psi_cx;
    // Check for convergence
    double root; // Root of the function
    if (std::abs(d) < 1e-8) {
      root = 2.0 * R::trigamma(c / 2.0) - l;
    } else {
      root = 2.0 * delta / d - l;
    }
    if (std::abs(root) < tol) break;
    // Newton's update
    double h = 2.0 * delta - l * d; 
    double num = d * h;
    double den = 2 * (d * (tri_x + tri_cx) - 2 * delta);
    if (!R_finite(den) || std::abs(den) < 1e-14) break;
    double x_new = x - num / den;
    // Domain enforcement: keep t ∈ (0 + eps, c/2 - eps)
    if (!R_finite(x_new) || x_new <= 0 + eps) {
      x_new = x / 2.0;   // safe fallback
    } else if (x_new >= c / 2 - eps) {
      x_new = 0.5 * (x + c / 2 - eps);
    }
    // Step-size stopping rule
    if (std::abs(x_new - x) < 1e-12) {
      x = x_new;
      break;
    }
    x = x_new;
  }
  
  return x;
}


/**
 * This function obtains random draws from a Polya propsal, where the          \
 * parameters of the Polya distribution are set via moment                     \
 * matching as in https://github.com/changwoo-lee/logisticbeta-reproduce/      \
 *                                                                             \
 * @param a_vec Eigen::VectorXd. Vector of "a" values                          \
 * @param b_vec Eigen::VectorXd. Vector of "b" values                          \
 * @param lambda_vec Eigen::VectorXd. Vector of "lambda" values                \
 *                                                                             \ 
 * @return Eigen::MatrixXd: A matrix containing the Polya draws (first row),   \
 * the new "a" parameters (second row) and the new "b" parameters (third row). \
 */
// [[Rcpp::export]]
Eigen::MatrixXd rand_polya_proposal(const Eigen::VectorXd & a_vec,
                                    const Eigen::VectorXd & b_vec,
                                    const Eigen::VectorXd & lambda_vec) {
  
  // Pre-compute constants and prepare the returns  ----------------------------
  const int T_tot = a_vec.size();
  const int trunc = 200;
  Rcpp::IntegerVector k = Rcpp::seq(0, trunc);
  Rcpp::NumericVector k_num(k.begin(), k.end());
  Eigen::MatrixXd out(3, T_tot);
  
  // Iterate over time points --------------------------------------------------
  for (int t = 0; t < T_tot; ++t) {
    
    // Find the parameters of the Polya proposal via moment matching
    double c = a_vec[t] + b_vec[t];
    double a_prime = optimize_moment_matching(c, lambda_vec[t]);
    double b_prime = c - a_prime;
    
    // Sample from the Polya distribution
    Rcpp::NumericVector denom = (a_prime + k_num) * (b_prime + k_num) / 2.0;
    Rcpp::NumericVector v = Rcpp::rexp(trunc + 1, 1) / denom;
    out(0, t) = Rcpp::sum(v); // Polya random draw 
    out(1, t) = a_prime; // Proposal a parameter 
    out(2, t) = b_prime; // Proposal b parameter
  }
  
  return out;
}
