#ifndef UTILS_H
#define UTILS_H

#include <RcppEigen.h>

Rcpp::IntegerVector row_unique_counts(Rcpp::NumericMatrix mat);

Eigen::VectorXd NumVec_to_EigenVec(Rcpp::NumericVector & x);

Eigen::VectorXf NumVec_to_EigenVec_float(const Rcpp::NumericVector & x);

Eigen::MatrixXf chol_float(const Eigen::MatrixXf & M);

#endif  // UTILS_H
