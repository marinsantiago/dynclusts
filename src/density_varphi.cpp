// -----------------------------------------------------------------------------
// Log-full-conditional density of var-phi and its gradient
// -----------------------------------------------------------------------------

#include <RcppEigen.h>

// [[Rcpp::depends(RcppEigen)]]

#include <iostream>
#include <Eigen/Dense>
#include <Eigen/IterativeLinearSolvers>
#include <cmath>
#include "density_varphi.h"


/**
 * This function evaluates the log-full-conditional density of varphi and its  \
 * first derivative, i.e., the gradient. See the Appendix of the manuscript    \
 * for additional details.                                                     \
 *                                                                             \
 * Important: For computational efficiency, this function uses single          \ 
 * precision (floats)!                                                         \
 *                                                                             \
 * @param varphi float: Value of varphi to evaluate the gradient at.           \
 * @param logdet_K_check float: Log-determinant of K.check.                    \                                   \
 * @param K_check_inv Eigen::MatrixXf: Inverse matrix of K.check.              \
 * @param K_check Eigen::MatrixXf: Current matrix K.check.                     \
 * @param D Eigen::MatrixXf: Matrix of squared distances.                      \
 * @param gamm Eigen::VectorXf: Vector of location-specific random effects.    \ 
 * @param tau2 float: Current value of tau2.                                   \
 * @param a_phi float: Value of a_phi.                                         \
 * @param b_phi float: Value of b_phi.                                         \  
 *                                                                             \
 * @return Eigen::VectorXf: A vector containing the log-full conditional       \
 * density of varphi and its gradient.                                         \
 */
// [[Rcpp::export]]
Eigen::VectorXf logdens_varphi(const float varphi, 
                               const float logdet_K_check, 
                               const Eigen::MatrixXf K_check_inv,
                               const Eigen::MatrixXf K_check,
                               const Eigen::MatrixXf D, 
                               const Eigen::VectorXf gamm, 
                               const float tau2, 
                               const float a_phi, 
                               const float b_phi) {
  
  // Pre-compute constants and prepare the returns  ----------------------------
  const Eigen::RowVectorXf gamma_t = gamm.transpose();
  const float inv_varphi3 = 1.0f / (varphi * varphi * varphi);
  Eigen::VectorXf out(2); // To store log-density and its gradient
  
  // Pre-compute (K_check_inv %*% gamma) / tau2 --------------------------------
  const Eigen::VectorXf K_gamma_tau2 = (K_check_inv * gamm) / tau2;
  
  // Compute quadratic form ----------------------------------------------------
  const float quad_form = gamma_t * K_gamma_tau2;
  
  // Compute log-density -------------------------------------------------------
  const float logdens = (a_phi - 1.0f) * std::logf(varphi) - b_phi * varphi - 
    (logdet_K_check + quad_form) / 2.0f;
  
  // Compute Hadamard product --------------------------------------------------
  const Eigen::MatrixXf H = (K_check.array() * D.array()) * inv_varphi3;
  
  // Compute Lambda.check^(-1) * Hadamard --------------------------------------
  const Eigen::MatrixXf M = K_check_inv * H;
  
  // Compute h1 and h2 ---------------------------------------------------------
  const float h1 = M.trace(); 
  const float h2 = gamma_t * M * K_gamma_tau2;
  
  // Compute the gradient of the log-density -----------------------------------
  const float grad = -b_phi + (a_phi - 1.0f) / varphi - (h1 - h2) / 2.0f;

  // Prepare the returns -------------------------------------------------------
  out[0] = logdens;
  out[1] = grad;
  
  return out;
}
