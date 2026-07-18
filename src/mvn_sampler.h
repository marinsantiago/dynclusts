#ifndef MVN_SAMPLER_H
#define MVN_SAMPLER_H

#include <RcppEigen.h>

Eigen::VectorXf rand_mvnorm(const Eigen::MatrixXf & precision_matrix,
                            const Eigen::VectorXf & parameter_vector);

Eigen::VectorXd rand_mvnorm_db(const Eigen::MatrixXd & precision_matrix,
                               const Eigen::VectorXd & parameter_vector);

#endif  // MVN_SAMPLER_H
