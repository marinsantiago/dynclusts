// -----------------------------------------------------------------------------
// Multivariate normal sampler
// -----------------------------------------------------------------------------

#include <RcppEigen.h>

// [[Rcpp::depends(RcppEigen)]]

#include <iostream>
#include <Eigen/Dense>
#include <Eigen/IterativeLinearSolvers>
#include <cmath>
#include "mvn_sampler.h"
#include "utils.h"


/**
 * Function to generate random draws from a multivariate normal distribution,  \
 * using the algorithm proposed by                                             \
 * Rue (2001) <https://doi.org/10.1111/1467-9868.00288>, which is based on a   \
 * Cholesky factorization of the precision matrix.                             \
 *                                                                             \
 * Important: For computational efficiency, this function uses single          \ 
 * precision (floats)! If needed, the function "rand_mvnorm_db()" below        \
 * uses double precision.                                                      \
 *                                                                             \
 * @param precision_matrix Eigen::MatrixXf: Precision matrix.                  \
 * @param parameter_vector Eigen::VectorXf: Parameter vector such that         \
 *  precision_matrix^(-1) * parameter_vector = mean of the MVN distribution.   \
 *                                                                             \
 * @return Eigen::VectorXf: Random draws from the MVN distribution.            \
 */
// [[Rcpp::export]]
Eigen::VectorXf rand_mvnorm(const Eigen::MatrixXf & precision_matrix,
                            const Eigen::VectorXf & parameter_vector) {
  
  // Pre-compute constants and prepare the returns -----------------------------
  const int mvnorm_dim = parameter_vector.size();
  Eigen::VectorXf out;
  
  // Cholesky factorization ----------------------------------------------------
  Eigen::LLT<Eigen::MatrixXf> llt(precision_matrix);
  // Solve the inner and outer systems using the Cholesky factorization
  Eigen::VectorXf mu_tilde = llt.solve(parameter_vector);
  
  // Sample z ~ MVN(0, I_p) ----------------------------------------------------
  Rcpp::NumericVector z_Rcpp = Rcpp::rnorm(mvnorm_dim, 0.0, 1.0);
  Eigen::VectorXf z = NumVec_to_EigenVec_float(z_Rcpp);
  
  // Transform z into a draw from our MVN distribution -------------------------
  Eigen::MatrixXf Lt = llt.matrixU();
  out = mu_tilde + Lt.triangularView<Eigen::Upper>().solve(z);
  
  return out;
}


/**
 * Function to generate random draws from a multivariate normal distribution,  \
 * using the algorithm proposed by                                             \
 * Rue (2001) <https://doi.org/10.1111/1467-9868.00288>, which is based on a   \
 * Cholesky factorization of the precision matrix.                             \
 *                                                                             \
 * Important: This function uses double precision, as such is more precise,    \
 * but computationally more expensive.                                         \
 *                                                                             \
 * @param precision_matrix Eigen::MatrixXd: Precision matrix.                  \
 * @param parameter_vector Eigen::VectorXd: Parameter vector such that         \
 *  precision_matrix^(-1) * parameter_vector = mean of the MVN distribution.   \
 *                                                                             \
 * @return Eigen::VectorXd: Random draws from the MVN distribution.            \
 */
// [[Rcpp::export]]
Eigen::VectorXd rand_mvnorm_db(const Eigen::MatrixXd & precision_matrix,
                               const Eigen::VectorXd & parameter_vector) {
  
  // Pre-compute constants and prepare the returns -----------------------------
  const int mvnorm_dim = parameter_vector.size();
  Eigen::VectorXd out;
  
  // Cholesky factorization ----------------------------------------------------
  Eigen::LLT<Eigen::MatrixXd> llt(precision_matrix);
  // Solve the inner and outer systems using the Cholesky factorization
  Eigen::VectorXd mu_tilde = llt.solve(parameter_vector);

  // Sample z ~ MVN(0, I_p) ----------------------------------------------------
  Rcpp::NumericVector z_Rcpp = Rcpp::rnorm(mvnorm_dim, 0.0, 1.0);
  Eigen::VectorXd z = NumVec_to_EigenVec(z_Rcpp);
  
  // Transform z into a draw from our MVN distribution -------------------------
  Eigen::MatrixXd Lt = llt.matrixU();
  out = mu_tilde + Lt.triangularView<Eigen::Upper>().solve(z);
  
  return out;
}
