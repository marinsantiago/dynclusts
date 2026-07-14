#ifndef LAMBDA_RATIO_H
#define LAMBDA_RATIO_H

#include <RcppEigen.h>

double lambda_ratio_log_l(const Eigen::VectorXd & center_lambda,
                          const Eigen::MatrixXd & cov_lambda);

#endif  // LAMBDA_RATIO_H
