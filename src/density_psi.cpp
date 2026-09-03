// -----------------------------------------------------------------------------
// Log-full-conditional density of psi and its gradient
// -----------------------------------------------------------------------------

#include <RcppEigen.h>

// [[Rcpp::depends(RcppEigen)]]

#include <iostream>
#include <Eigen/Dense>
#include <Eigen/IterativeLinearSolvers>
#include <cmath>
#include "density_psi.h"


/**
 * This function evaluates the log-full-conditional density of psi and its     \
 * first derivative, i.e., the gradient. See the Appendix of the manuscript    \
 * for additional details.                                                     \
 *                                                                             \
 * Important: For computational efficiency, this function uses single          \ 
 * precision (floats)!                                                         \
 *                                                                             \
 * @param psi float.                                                           \
 * @param H1T1 float.                                                          \                                \
 * @param sum_e1T float.                                                       \
 * @param sum_e_mid float.                                                     \
 * @param sum_e_dif float.                                                     \
 *                                                                             \ 
 * @return Eigen::VectorXf: A vector containing the log-full conditional       \
 * density of psi and its gradient.                                            \
 */
// [[Rcpp::export]]
Eigen::VectorXf logdens_psi(const float psi, 
                            const float H1T1, 
                            const float sum_e1T,
                            const float sum_e_mid,
                            const float sum_e_dif) {
  
  // Pre-compute constants and prepare the returns  ----------------------------
  Eigen::VectorXf out(2); // To store log-density and its gradient
  const float psi2 = std::pow(psi, 2.0f);
  const float psi2_diff = 1.0f - psi2;
  const float psi2_diff_sq = std::pow(psi2_diff, 2.0f);
  const float psi2_add = 1.0f + psi2;
  
  // Compute log-density -------------------------------------------------------
  float logdens = -H1T1 * std::log(psi2_diff) / 2.0f;
  logdens -= (sum_e1T + psi2_add * sum_e_mid) / (2.0f * psi2_diff);
  logdens += 2.0f * psi * sum_e_dif / (2.0f * psi2_diff);
  
  // Compute gradient ----------------------------------------------------------
  float grad = psi * H1T1 / psi2_diff;
  grad -= psi * (sum_e1T + 2.0f * sum_e_mid) / psi2_diff_sq;
  grad += psi2_add * sum_e_dif / psi2_diff_sq;
  
  // Prepare the returns -------------------------------------------------------
  out[0] = logdens;
  out[1] = grad;
  
  return out;
}
