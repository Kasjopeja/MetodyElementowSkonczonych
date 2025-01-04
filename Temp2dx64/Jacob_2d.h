#ifndef TEMP2DX64_JACOB_2D_H
#define TEMP2DX64_JACOB_2D_H

#include <vector>
#include <iostream>
#include <cmath>
#include "GlobData.h"

extern GlobData data;

void Jacob_2d(std::vector<std::vector<double>>& J_,
              std::vector<std::vector<double>>& J_inv,
              int P, int N_p, int NEL,
              const std::vector<std::vector<double>>& N1,
              const std::vector<std::vector<double>>& N2,
              double& DetJ);
void Inv_MAT(int N, const std::vector<std::vector<double>>& Mat, std::vector<std::vector<double>>& Inv);

#endif //TEMP2DX64_JACOB_2D_H
