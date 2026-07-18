// -----------------------------------------------------------------------------
// Utils and linear algebra helpers
// -----------------------------------------------------------------------------

#include <RcppEigen.h>

// [[Rcpp::depends(RcppEigen)]]

#include <iostream>
#include <Eigen/Dense>
#include <cmath>
#include <vector>
#include "utils.h"


/**
 * Counts the number of unique elements in each row of a numeric matrix.       \
 *                                                                             \
 * @param mat Rcpp::NumericMatrix: Matrix where each row will be analyzed.     \
 *                                                                             \
 * @return Rcpp::IntegerVector: Count of unique elements per row.              \
 */
// [[Rcpp::export]]
Rcpp::IntegerVector row_unique_counts(Rcpp::NumericMatrix mat) {
  
  // Prepare the returns
  int nrow = mat.nrow();
  int ncol = mat.ncol();
  Rcpp::IntegerVector result(nrow);
  
  for (int i = 0; i < nrow; ++i) {
    std::unordered_set<double> unique_vals;
    for (int j = 0; j < ncol; ++j) {
      unique_vals.insert(mat(i, j));
    }
    result[i] = unique_vals.size();
  }
  
  return result;
}


/**
 * Function to convert an Rcpp::NumericVector into an Eigen::VectorXd.         \ 
 *                                                                             \
 * @param x Rcpp::NumericVector: A vector to convert into Eigen::VectorXd      \
 *                                                                             \
 * @return Eigen::VectorXd: Converted vector.                                  \
 */
// [[Rcpp::export]]
Eigen::VectorXd NumVec_to_EigenVec(Rcpp::NumericVector & x) {
  Eigen::Map<Eigen::VectorXd> out(Rcpp::as<Eigen::Map<Eigen::VectorXd>>(x));
  return out;
}


/**
 * Function to convert an Rcpp::NumericVector into an Eigen::VectorXf (float)! \ 
 *                                                                             \
 * @param x Rcpp::NumericVector: A vector to convert into Eigen::VectorXf      \
 *                                                                             \
 * @return Eigen::VectorXf: Converted vector (as float)!                       \
 */
// [[Rcpp::export]]
Eigen::VectorXf NumVec_to_EigenVec_float(const Rcpp::NumericVector & x) {
  Eigen::VectorXf out(x.size());
  for (int i = 0; i < x.size(); ++i) {
    out[i] = static_cast<float>(x[i]);
  }
  return out;
}


/**
 * Computes the Cholesky factorization of a symmetric positive-definite        \ 
 * square matrix.                                                              \ 
 *                                                                             \
 * @param M Eigen::MatrixXd: Symmetric positive-definite square matrix.        \
 *                                                                             \
 * @return Eigen::MatrixXd: Cholesky factor of M.                              \
 */
// [[Rcpp::export]]
Eigen::MatrixXf chol_float(const Eigen::MatrixXf & M) {
  // Perform Cholesky decomposition
  Eigen::LLT<Eigen::MatrixXf> llt(M);
  // Check if decomposition was successful
  if (llt.info() != Eigen::Success) {
    Rcpp::stop("Matrix is not positive definite.");
  }
  // Return the Cholesky factor of M
  return llt.matrixU();
}

