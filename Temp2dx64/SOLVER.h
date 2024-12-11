#ifndef TEMP2DX64_SOLVER_H
#define TEMP2DX64_SOLVER_H

#include "GlobData.h"
#include "FeSM_heat.h"

extern GlobData data;

void DLSAQS(int, std::vector<std::vector<double>>&, int, int, std::vector<double>&, std::vector<double>&);
void SOLVER();

#endif //TEMP2DX64_SOLVER_H
