#ifndef INT_SAMPLING_H
#define INT_SAMPLING_H

#include <RcppEigen.h>

Rcpp::IntegerVector int_sampling(const int & K,
                                 const int & num_samples,
                                 const Rcpp::NumericVector & probs);

Rcpp::IntegerVector int_sampling_rows(const Rcpp::NumericMatrix & probs_mat);

#endif  // INT_SAMPLING_H
