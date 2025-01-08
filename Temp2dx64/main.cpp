#include <iostream>
#include <cmath>

#include "GlobData.h"
#include "IniEL4.h"
#include "InpData.h"
#include "GenGrid2d.h"
#include "SaveGridToVTK.h"
#include "SOLVER.h"
#include <ctime>

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
    SaveResultToVTK("grid_step_0.vtk");
    //SaveGridToVTK("C:\\Users\\Lenovo\\CLionProjects\\Metody Elementow Skonczonych\\grid_output.vtk");

    SetControlPoints();
    ALLOCATE_Matrix();
    WriteControlPointsBegin();

    GlobData test = data; // DEBBUG

    // Obliczenia
//    Asr = data.mK / (data.mC * data.mR);
//    data.mdTime = pow(data.mB0 / (1e3 * data.mNhB), 2) / (0.5 * Asr);
    WriteControlPoints();


    Ntau = static_cast<int>(data.mTime / data.mdTime);
    data.mTau = 0.0;

    clock_t start = clock();
    for (n = 1; n <= Ntau; ++n) {
        data.mTau += data.mdTime;

        SOLVER();
        std::string filename = "grid_step_" + std::to_string(n) + ".vtk";
        SaveResultToVTK(filename);

        WriteControlPoints();
    }
    std::cout << "Czas wykonywania: " << clock() - start  << std:: endl;

    // Zakończenie programu
    SaveGridToVTK("C:\\Users\\Lenovo\\CLionProjects\\Metody Elementow Skonczonych\\grid_output.vtk");
    return 0;
}
