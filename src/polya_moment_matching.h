#ifndef POLYA_MOMENT_MATCHING_H
#define POLYA_MOMENT_MATCHING_H

#include <RcppEigen.h>

double optimize_moment_matching(double c, double l);

Eigen::MatrixXd rand_polya_proposal(const Eigen::VectorXd & a_vec,
                                    const Eigen::VectorXd & b_vec,
                                    const Eigen::VectorXd & lambda_vec);

#endif  // POLYA_MOMENT_MATCHING_H
