#ifndef TEMP2DX64_JACOB_2D_H
#define TEMP2DX64_JACOB_2D_H

#include <vector>
#include <iostream>
#include <cmath>
#include "GlobData.h"

extern GlobData data;

void Jacob_2d(double J[2][2],
              double J_inv[2][2],
              int p,
              int N_p,
              int NBN,
              const std::vector<std::vector<double>> &N1,
              const std::vector<std::vector<double>> &N2,
              const std::vector<double> &X,
              const std::vector<double> &Y,
              double &DetJ);
void Inv_MAT(int N, const std::vector<std::vector<double>>& Mat, std::vector<std::vector<double>>& Inv);

#endif //TEMP2DX64_JACOB_2D_H
