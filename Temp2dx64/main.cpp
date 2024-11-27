#include <iostream>
#include <cmath>

#include "GlobData.h"
#include "IniEL4.h"
#include "InpData.h"
#include "GenGrid2d.h"

// dane globalne
GlobData data;

// Tymczasowe deklaracjie funkcji
void GenGrid2d(double mH0, double mB0, int mNhH, int mNhB, Gr2d &mGr) {};
void SetControlPoints() {};
void WriteControlPointsBegin() {};
void WriteControlPoints() {};
void SOLVER() {};

int main() {
    // Deklaracje zmiennych
    double Asr, dTauMax, TauP;
    int n, Ntau, i, iErr;

    // Inicjalizacja
    IniEL4();
    InpData();
    GenGrid2d(data.mH0, data.mB0, data.mNhH, data.mNhB, data.mGr);
    SetControlPoints();
    ALLOCATE_Matrix();
    WriteControlPointsBegin();

    // Obliczenia
    Asr = data.mK / (data.mC * data.mR);
    double mdTime = pow(data.mB0 / (1e3 * data.mNhB), 2) / (0.5 * Asr);
    WriteControlPoints();

    Ntau = static_cast<int>(data.mTime / mdTime) + 1;
    mdTime = data.mTime / static_cast<double>(Ntau);
    double mTau = 0.0;

    for (n = 1; n <= Ntau; ++n) {
        mTau += mdTime;

        SOLVER();

        WriteControlPoints();
    }

    // Zakończenie programu
    std::cout << "Closing files and deallocating memory." << std::endl;
    return 0;
}
