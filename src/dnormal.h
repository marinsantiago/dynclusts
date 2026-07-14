#ifndef DNORMAL_H
#define DNORMAL_H

#include <RcppEigen.h>

Eigen::ArrayXXd dnormal(const Eigen::ArrayXXd & y,
                        const Eigen::ArrayXXd & mu,
                        const double sigma2);

#endif  // DNORMAL_H
