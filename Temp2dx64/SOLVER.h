#ifndef TEMP2DX64_SOLVER_H
#define TEMP2DX64_SOLVER_H

#include "GlobData.h"
#include "FeSM_heat.h"
#include <Eigen/Sparse>
#include <Eigen/SparseLU>

extern GlobData data;
typedef Eigen::Triplet<double> T;

void SOLVER();

#endif //TEMP2DX64_SOLVER_H
