#ifndef BIDIAGONAL_SOLVERS_H
#define BIDIAGONAL_SOLVERS_H

#include <RcppEigen.h>

Rcpp::NumericVector solve_lower_bidiagonal(const Rcpp::NumericVector & a,
                                           const Rcpp::NumericVector & l,
                                           const Rcpp::NumericVector & b);


Rcpp::NumericVector solve_upper_bidiagonal(const Rcpp::NumericVector & a, 
                                           const Rcpp::NumericVector & u, 
                                           const Rcpp::NumericVector & b);

#endif  // BIDIAGONAL_SOLVERS_H
