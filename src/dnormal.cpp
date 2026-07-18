// -----------------------------------------------------------------------------
// Density function of the normal distribution
// -----------------------------------------------------------------------------

#include <RcppEigen.h>

// [[Rcpp::depends(RcppEigen)]]

#include <iostream>
#include <Eigen/Dense>
#include <cmath>
#include "dnormal.h"


/**
 * This function computes the element-wise probability density function (PDF) of
 * a normal distribution.
 *
 * @param y Eigen::ArrayXXd: Observed values as a matrix.
 * @param mu Eigen::ArrayXXd: Means of the distribution as a matrix.
 * @param sigma2 Eigen::ArrayXXd: Variance of the distribution as a scalar.
 *
 * @return Eigen::ArrayXXd: Evaluated density at each y.
 */
// [[Rcpp::export]]
Eigen::ArrayXXd dnormal(const Eigen::ArrayXXd & y,
                        const Eigen::ArrayXXd & mu,
                        const double sigma2) {
  
  static const double log2pi = std::log(2.0 * M_PI);
  Eigen::ArrayXXd exponent = -0.5 * ((y - mu).square() / sigma2);
  double denom = std::log(std::sqrt(sigma2)) + 0.5 * log2pi;
  return (exponent - denom).exp(); // Returns density
}
