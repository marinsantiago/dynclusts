#ifndef DENSITY_PSI_H
#define DENSITY_PSI_H

#include <RcppEigen.h>

Eigen::VectorXf logdens_psi(const float psi, 
                            const float H1T1, 
                            const float sum_e1T,
                            const float sum_e_mid,
                            const float sum_e_dif);

#endif  // DENSITY_PSI_H
