// -----------------------------------------------------------------------------
// Random integer sampling
// -----------------------------------------------------------------------------

#include <RcppEigen.h>

// [[Rcpp::depends(RcppEigen)]]

#include <iostream>
#include <Eigen/Dense>
#include <Eigen/IterativeLinearSolvers>
#include <cmath>
#include <vector>
#include <algorithm>
#include <random>
#include "int_sampling.h"
using namespace Rcpp;


/**
 * Random sampling with replacement.                                           \
 *                                                                             \
 * @param K int: A positive number, the number of items to choose from.        \
 * @param num_samples int: A non-negative integer giving the number of items   \
 *  to choose.                                                                 \
 * @param probs NumericVector: a vector of probability weights for obtaining   \
 *  the elements of the vector being sampled.                                  \
 *                                                                             \
 * @return IntegerVector: An integer vector of length `num_samples` with       \
 *  elements from 1:K                                                          \
 */
// [[Rcpp::export]]
Rcpp::IntegerVector int_sampling(const int & K,
                                 const int & num_samples,
                                 const Rcpp::NumericVector & probs) {
  
  // Prepare the returns -------------------------------------------------------
  Rcpp::IntegerVector out(num_samples);
  
  // Pre-compute constants -----------------------------------------------------
  Rcpp::NumericVector cum_probs = Rcpp::cumsum(probs);
  Rcpp::NumericVector random_values = Rcpp::runif(num_samples);
  auto cum_probs_begin = cum_probs.begin();
  auto cum_probs_end = cum_probs.end();
  
  for (int i = 0; i < num_samples; ++i) {
    out[i] = std::lower_bound(cum_probs_begin,
                              cum_probs_end,
                              random_values[i]) - cum_probs_begin + 1;
  }
  
  return out;
}

/**
 * Samples one index from each row of a probability matrix.                    \
 *                                                                             \
 * @param probs Eigen::MatrixXd: Matrix of probabilities (or weights).         \
 *                                                                             \
 * @return IntegerVector of length n_rows, where each element is the sampled   \
 *   index from the corresponding row.                                         \
 *                                                                             \
 */
// [[Rcpp::export]]
Rcpp::IntegerVector int_sampling_rows(const Rcpp::NumericMatrix & probs_mat) {
  
  // Prepare the returns -------------------------------------------------------
  int n_rows = probs_mat.rows();
  int n_cols = probs_mat.cols();
  Rcpp::IntegerVector samples(n_rows);
  
  // Iterate over the rows of probs_mat ----------------------------------------
  for (int i = 0; i < n_rows; ++i) {
    Rcpp::NumericVector row = probs_mat.row(i);
    row = row / sum(row);
    Rcpp::IntegerVector int_sample = int_sampling(n_cols, 1, row);
    samples[i] = int_sample[0];
  }
  
  return samples;
}
