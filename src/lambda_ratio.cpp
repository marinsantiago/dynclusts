// -----------------------------------------------------------------------------
// "Ratio log-L()" function for lambda
// -----------------------------------------------------------------------------

#include <RcppEigen.h>

// [[Rcpp::depends(RcppEigen)]]

#include <iostream>
#include <Eigen/Dense>
#include <Eigen/IterativeLinearSolvers>
#include <cmath>
#include "lambda_ratio.h"


/**
 * Function to evaluate "L" in the acceptance ratio for lambda, using the      \
 * Metropolis-Hastings algorithm proposed by:                                  \                    
 * Lee, Zito, Sang & Dunson, D. B. (2025+) <https://doi.org/10.1214/25-BA1541> \
 *                                                                             \
 * @param center_lambda Eigen::VectorXd: Center of the distribution.           \
 * @param cov_lambda Eigen::MatrixXd: Covariance of the distribution.          \
 *                                                                             \
 * @return double: "L" based on the input lambda.                              \
 */
// [[Rcpp::export]]
double lambda_ratio_log_l(const Eigen::VectorXd & center_lambda,
                          const Eigen::MatrixXd & cov_lambda) {
  
  // Pre-compute constants and prepare the returns -----------------------------
  const int d = center_lambda.size();
  const double lnSqrt2Pi = 0.5 * std::log(2 * M_PI);
  double out;
  
  // Cholesky decomposition ----------------------------------------------------
  typedef Eigen::LLT<Eigen::MatrixXd> Chol;
  Chol chol(cov_lambda);
  
  // Log-determinant -----------------------------------------------------------
  const Chol::Traits::MatrixL & L = chol.matrixL();
  double log_det_L = L.toDenseMatrix().diagonal().array().log().sum();
  
  // Compute quadratic form ----------------------------------------------------
  double quadform = (L.solve(center_lambda)).squaredNorm();
  
  // Evaluate log-density ------------------------------------------------------
   out = -d * lnSqrt2Pi - log_det_L - 0.5 * quadform;
  
  return out; 
}
