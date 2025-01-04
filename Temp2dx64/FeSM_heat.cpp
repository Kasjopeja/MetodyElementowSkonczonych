#include "FeSM_heat.h"

void FeSM_heat(int NEL) {
    // Reset element matrix and vector
    for (int i = 0; i < 4; ++i) {
        for (int j = 0; j < 4; ++j) {
            data.est[i][j] = 0.0;
        }
    }
    for (int i = 0; i < 4; ++i) {
        data.r[i] = 0.0;
    }

    PRE_heat_mat(NEL);
    PRE_heat_pov_mat(NEL);
}


void PRE_heat_mat(int NEL) {
    std::vector<std::vector<double>> J_(2, std::vector<double>(2, 0.0));
    std::vector<std::vector<double>> J_inv(2, std::vector<double>(2, 0.0));

    double DetJ, Ni, Nn, Hin, Cin;
    double T0p;
    std::vector<double> Ndx(4), Ndy(4), Temp_0(4);

    for (int i = 0; i < data.mEL4.nbn; ++i) {
        int Id = std::abs(data.mGr.EL[NEL - 1].nop[i]);
        Temp_0[i] = data.mGr.ND[Id - 1].t;
    }

    for (int P = 0; P < data.mEL4.N_p; ++P) {
        Jacob_2d(J_, J_inv, P, data.mEL4.N_p, NEL, data.mEL4.N1, data.mEL4.N2, DetJ);
        T0p = 0.0;

        for (int i = 0; i < data.mEL4.nbn; ++i) {
            Ndx[i] = data.mEL4.N1[i][P] * J_inv[0][0] + data.mEL4.N2[i][P] * J_inv[0][1];
            Ndy[i] = data.mEL4.N1[i][P] * J_inv[1][0] + data.mEL4.N2[i][P] * J_inv[1][1];
            Ni = data.mEL4.Nf[i][P];
            T0p += Temp_0[i] * Ni;
        }

        DetJ = std::abs(DetJ) * data.mEL4.W[P];
        for (int n = 0; n < data.mEL4.nbn; ++n) {
            for (int i = 0; i < data.mEL4.nbn; ++i) {
                Ni = data.mEL4.Nf[i][P];
                Nn = data.mEL4.Nf[n][P];
                Hin = data.mK * (Ndx[n] * Ndx[i] + Ndy[n] * Ndy[i]) * DetJ;
                Cin = data.mC * data.mR * Nn * Ni * DetJ;
                data.est[n][i] += Hin + Cin / data.mdTime;
                data.r[n] += (Cin / data.mdTime) * T0p;
            }
        }
    }
}

void PRE_heat_pov_mat(int NEL) {
    double DetJ, Ni, Nn, Pn;
    std::vector<double> X(4), Y(4);

    for (int i = 0; i < 4; ++i) {
        int Id = std::abs(data.mGr.EL[NEL - 1].nop[i]);
        X[i] = data.mGr.ND[Id - 1].x;
        Y[i] = data.mGr.ND[Id - 1].y;
    }

    for (int iPov = 0; iPov < data.mGr.EL[NEL - 1].Npov; ++iPov) {
        int id = data.mGr.EL[NEL - 1].aPov[iPov];
        switch (id) {
            case 1:
                DetJ = std::sqrt(std::pow(X[3] - X[0], 2) + std::pow(Y[3] - Y[0], 2));
                break;
            case 2:
                DetJ = std::sqrt(std::pow(X[0] - X[1], 2) + std::pow(Y[0] - Y[1], 2));
                break;
            case 3:
                DetJ = std::sqrt(std::pow(X[1] - X[2], 2) + std::pow(Y[1] - Y[2], 2));
                break;
            case 4:
                DetJ = std::sqrt(std::pow(X[2] - X[3], 2) + std::pow(Y[2] - Y[3], 2));
                break;
        }

        for (int P = 0; P < 2; ++P) {
            for (int n = 0; n < 4; ++n) {
                for (int i = 0; i < 4; ++i) {
                    Ni = data.mEL4.Sf[id - 1].Nf[i][P];
                    Nn = data.mEL4.Sf[id - 1].Nf[n][P];
                    data.est[n][i] += data.mAlfa * Nn * Ni * DetJ;
                }
                Pn = data.mAlfa * data.mT_otoczenia * Nn * DetJ;
                data.r[n] += Pn;
            }
        }
    }
}
