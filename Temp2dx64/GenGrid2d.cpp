#include "GenGrid2d.h"

// Funkcja generacji siatki MES
void GenGrid2d() {
    data.mGr.nh = data.mNhH * data.mNhB;
    data.mGr.ne = (data.mNhH - 1) * (data.mNhB - 1);
    data.mGr.nbn = 4;
    data.mGr.ncn = 4;
    data.mGr.nhPov = data.mNhB;

    data.mGr.ND.resize(data.mGr.nh);
    data.mGr.EL.resize(data.mGr.ne);

    double dx = data.mB0 / (data.mNhB - 1);
    double dy = data.mH0 / (data.mNhH - 1);

    // Generacja węzłów
    int nodeIndex = 0;
    for (int i = 0; i < data.mNhB; ++i) {
        double x = i * dx;
        for (int j = 0; j < data.mNhH; ++j) {
            double y = j * dy;
            data.mGr.ND[nodeIndex++] = {x, y, data.mTbegin, 0.0, 0}; // Współrzędne, temperatura, status
        }
    }

    // Generacja elementów
    int elementIndex = 0;
    for (int i = 0; i < data.mNhB - 1; ++i) {
        for (int j = 0; j < data.mNhH - 1; ++j) {
            int i1 = i * data.mNhH + j;
            int i2 = (i + 1) * data.mNhH + j;
            int i3 = (i + 1) * data.mNhH + j + 1;
            int i4 = i * data.mNhH + j + 1;

            data.mGr.EL[elementIndex++] = {{i1, i2, i3, i4}, 0, {}}; // Indeksy węzłów, liczba powierzchni, lista powierzchni
        }
    }

    // Aktualizacja statusów węzłów
    for (auto &node : data.mGr.ND) {
        double x = node.x;
        double y = node.y;

        if (x >= data.mB0 - 1e-5 || x <= 1e-5 || y >= data.mH0 - 1e-5 || y <= 1e-5) {
            node.status = 1; // Węzeł brzegowy
        }
    }

    // Aktualizacja powierzchni w elementach
    for (auto &element : data.mGr.EL) {
        int St[4] = {0};
        for (int i = 0; i < 4; ++i) {
            St[i] = data.mGr.ND[element.nop[i]].status;
        }

        std::vector<bool> St_OK(4);
        St_OK[1] = (St[0] >= 1) && (St[1] >= 1);
        St_OK[2] = (St[1] >= 1) && (St[2] >= 1);
        St_OK[3] = (St[2] >= 1) && (St[3] >= 1);
        St_OK[0] = (St[3] >= 1) && (St[0] >= 1);

        for (int i = 0; i < 4; ++i) {
            if (St_OK[i]) {
                element.Npov++;
                element.aPov[element.Npov - 1] = i + 1; // Dodanie lokalnego numeru powierzchni
            }
        }
    }
}

// Funkcja zapisująca punkty kontrolne
void WriteControlPoints() {
    std::ofstream files[2] = {
            std::ofstream("OutDataT.txt", std::ios::app),
            std::ofstream("OutDataCR.txt", std::ios::app)
    };

    for (std::ofstream& file : files) {
        if (!file.is_open()) {
            std::cerr << "Nie można otworzyć plików do zapisu punktów kontrolnych!" << std::endl;
            return;
        }
    }

    for (int f = 0; f < 2; ++f) {
        files[f] << std::fixed << std::setprecision(4) << data.mTau << " ";
        for (int i = 0; i < 9; ++i) {
            int pointIndex = data.mContrPoints[i];
            files[f] << std::fixed << std::setprecision(1)
                     << (f == 0 ? data.mGr.ND[pointIndex].t : data.mGr.ND[pointIndex].CR)
                     << " ";
        }
        files[f] << "\n";
    }
}


// Funkcja ustawiająca punkty kontrolne
void SetControlPoints() {
    data.mcpX[0] = 0.0;
    data.mcpX[1] = data.mB0 / 2.0;
    data.mcpX[2] = data.mB0;
    data.mcpX[3] = 0.0;
    data.mcpX[4] = data.mB0 / 2.0;
    data.mcpX[5] = data.mB0;
    data.mcpX[6] = 0.0;
    data.mcpX[7] = data.mB0 / 2.0;
    data.mcpX[8] = data.mB0;

    data.mcpY[0] = 0.0;
    data.mcpY[1] = 0.0;
    data.mcpY[2] = 0.0;
    data.mcpY[3] = data.mH0 / 2.0;
    data.mcpY[4] = data.mH0 / 2.0;
    data.mcpY[5] = data.mH0 / 2.0;
    data.mcpY[6] = data.mH0;
    data.mcpY[7] = data.mH0;
    data.mcpY[8] = data.mH0;

    for (int j = 0; j < 9; ++j) {
        double Rmin = 1e10;
        for (int i = 0; i < data.mGr.nh; ++i) {
            double Rr = std::sqrt(
                    std::pow(data.mcpX[j] - data.mGr.ND[i].x, 2) +
                    std::pow(data.mcpY[j] - data.mGr.ND[i].y, 2));
            if (Rr <= Rmin) {
                data.mContrPoints[j] = i;
                Rmin = Rr;
            }
        }
    }
}

void WriteControlPointsBegin() {
    std::ofstream files[2] = {
            std::ofstream("OutDataT.txt"),
            std::ofstream("OutDataCR.txt")
    };

    for (std::ofstream& file : files) {
        if (!file.is_open()) {
            std::cerr << "Failed to open output file!" << std::endl;
            return;
        }

        file << " **********************************************************\n";
        file << " * PROGRAM PLATE2d                                        *\n";
        file << " ********************************************************** \n";
        file << "\n\n";
        file << " ********** Coordinates of the control points ************\n";
        for (int i = 0; i < 9; ++i) {
            file << " No=" << std::setw(4) << data.mContrPoints[i]
                 << " X=" << std::fixed << std::setprecision(2) << std::setw(8) << data.mcpX[i]
                 << " Y=" << std::fixed << std::setprecision(2) << std::setw(8) << data.mcpY[i] << "\n";
        }

        file << " ********************************************************** \n";
    }
}