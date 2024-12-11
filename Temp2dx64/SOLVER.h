#ifndef TEMP2DX64_SOLVER_H
#define TEMP2DX64_SOLVER_H

#include "GlobData.h"

extern GlobData data;

void FeSM_heat(int);
void DLSAQS(int, std::vector<std::vector<double>>&, int, int, std::vector<double>&, std::vector<double>&);
void SOLVER();

#endif //TEMP2DX64_SOLVER_H
