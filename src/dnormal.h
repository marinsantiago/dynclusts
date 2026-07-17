#ifndef DNORMAL_H
#define DNORMAL_H

#include <RcppEigen.h>

Eigen::ArrayXXd dnormal_disjoint(const Eigen::ArrayXXd & y,
                                 const Eigen::ArrayXXd & mu,
                                 const Eigen::ArrayXXd & sigma2);

Eigen::ArrayXXd dnormal_joint(const Eigen::ArrayXXd & y,
                              const Eigen::ArrayXXd & mu,
                              const double sigma2);

#endif  // DNORMAL_H
