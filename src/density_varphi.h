#ifndef DENSITY_VARPHI_H
#define DENSITY_VARPHI_H

#include <RcppEigen.h>

Eigen::VectorXf logdens_varphi(const float varphi, 
                               const float logdet_K_check, 
                               const Eigen::MatrixXf K_check_inv,
                               const Eigen::MatrixXf K_check,
                               const Eigen::MatrixXf D, 
                               const Eigen::VectorXf gamm, 
                               const float tau2, 
                               const float a_phi, 
                               const float b_phi);

#endif  // DENSITY_VARPHI_H
