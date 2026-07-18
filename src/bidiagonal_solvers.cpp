// -----------------------------------------------------------------------------
// Bi-diagonal solvers
// -----------------------------------------------------------------------------

#include <Rcpp.h>
#include "bidiagonal_solvers.h"

using namespace Rcpp;


/**
 * Solves a linear system Ax = b where A is a lower bi-diagonal matrix via     \
 * forward substitution. The matrix A is represented only by its diagonal and  \
 * sub-diagonal vectors.                                                       \
 *                                                                             \ 
 * @param a Rcpp::NumericVector Diagonal elements of A.                        \
 * @param l Rcpp::NumericVector Sub-diagonal elements of A.                    \
 * @param b Rcpp::NumericVector Right-hand side vector.                        \
 *                                                                             \ 
 * @return Rcpp::NumericVector Solution to the system, x.                      \
 *                                                                             \ 
 * @note Assumes that all diagonal elements in 'a' are non-zero.               \
 */
// [[Rcpp::export]]
Rcpp::NumericVector solve_lower_bidiagonal(const Rcpp::NumericVector & a,
                                           const Rcpp::NumericVector & l,
                                           const Rcpp::NumericVector & b) {
  
  // Pre-compute constants and prepare the returns -----------------------------
  int n = b.size();
  Rcpp::NumericVector x(n);
  
  // Forward substitution ------------------------------------------------------
  x[0] = b[0] / a[0];
  for (int i = 1; i < n; ++i) {
    x[i] = (b[i] - l[i - 1] * x[i - 1]) / a[i];
  }
  
  return x;
}


/**
 * Solves a linear system Ax = b where A is a upper bi-diagonal matrix via     \
 * backward substitution. The matrix A is represented only by its diagonal and \
 * super-diagonal vectors.                                                     \
 *                                                                             \ 
 * @param a Rcpp::NumericVector Diagonal elements of A.                        \
 * @param u Rcpp::NumericVector Super-diagonal elements of A.                  \
 * @param b Rcpp::NumericVector Right-hand side vector.                        \
 *                                                                             \ 
 * @return Rcpp::NumericVector Solution to the system, x.                      \
 *                                                                             \ 
 * @note Assumes that all diagonal elements in 'a' are non-zero.               \
 */
// [[Rcpp::export]]
Rcpp::NumericVector solve_upper_bidiagonal(const Rcpp::NumericVector & a, 
                                           const Rcpp::NumericVector & u, 
                                           const Rcpp::NumericVector & b) {
  
  // Pre-compute constants and prepare the returns -----------------------------
  int n = b.size();
  Rcpp::NumericVector x(n);
  
  // Backward substitution -----------------------------------------------------
  x[n - 1] = b[n - 1] / a[n - 1];
  for (int i = n - 2; i >= 0; --i) {
    x[i] = (b[i] - u[i] * x[i + 1]) / a[i];
  }
  
  return x;
}
