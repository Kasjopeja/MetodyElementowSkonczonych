#include <iostream>
#include <cmath>

#include "GlobData.h"
#include "IniEL4.h"
#include "InpData.h"
#include "GenGrid2d.h"
#include "SaveGridToVTK.h"
#include "SOLVER.h"

// dane globalne
GlobData data;

int main() {
    // Deklaracje zmiennych
    double Asr, dTauMax, TauP;
    int n, Ntau, i, iErr;

    // Inicjalizacja
    IniEL4();
    InpData();
    GenGrid2d();

    SetControlPoints();
    ALLOCATE_Matrix();
    WriteControlPointsBegin();

    GlobData test = data; // DEBBUG

    // Obliczenia
    Asr = data.mK / (data.mC * data.mR);
    data.mdTime = pow(data.mB0 / (1e3 * data.mNhB), 2) / (0.5 * Asr);
    WriteControlPoints();

    Ntau = static_cast<int>(data.mTime / data.mdTime);
    data.mTau = 0.0;

    for (n = 1; n <= Ntau; ++n) {
        data.mTau += data.mdTime;

        SOLVER();

        WriteControlPoints();
    }

    // Zakończenie programu
    SaveGridToVTK("C:\\Users\\Lenovo\\CLionProjects\\Metody Elementow Skonczonych\\grid_output.vtk");
    return 0;
}
