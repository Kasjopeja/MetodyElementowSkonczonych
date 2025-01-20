#include "GlobData.h"
extern GlobData data;

void IniEL4()
{
    double Alfa = 1.0 / sqrt(3.0);
    data.mEL4.nbn = 4;   // Liczba węzłów na element
    data.mEL4.N_p = 4;   // Liczba punktów całkowania

    data.mEL4.N1.resize(data.mEL4.nbn, std::vector<double>(data.mEL4.N_p));
    data.mEL4.N2.resize(data.mEL4.nbn, std::vector<double>(data.mEL4.N_p));
    data.mEL4.Nf.resize(data.mEL4.nbn, std::vector<double>(data.mEL4.N_p));
    data.mEL4.P.resize(data.mEL4.N_p);
    data.mEL4.W.resize(data.mEL4.N_p, 1.0);
    data.mEL4.L.resize(data.mEL4.nbn);

    // Lokalne współrzędne punktów całkowania
    data.mEL4.P[0].E = -Alfa; data.mEL4.P[0].N = -Alfa;
    data.mEL4.P[1].E = Alfa; data.mEL4.P[1].N = -Alfa;
    data.mEL4.P[2].E = Alfa; data.mEL4.P[2].N = Alfa;
    data.mEL4.P[3].E = -Alfa; data.mEL4.P[3].N = Alfa;

    // Set local coordinates for element nodes
    data.mEL4.L[0].E = -1; data.mEL4.L[0].N = -1;
    data.mEL4.L[1].E = 1; data.mEL4.L[1].N = -1;
    data.mEL4.L[2].E = 1; data.mEL4.L[2].N = 1;
    data.mEL4.L[3].E = -1; data.mEL4.L[3].N = 1;

    // Funkcje kształtu i ich pochodne po współrzędnych lokalnych
    for (int iP = 0; iP < data.mEL4.N_p; iP++) {
        double L1 = data.mEL4.P[iP].E;
        double L2 = data.mEL4.P[iP].N;

        data.mEL4.Nf[0][iP] = 0.25 * (1 - L1) * (1 - L2);
        data.mEL4.N1[0][iP] = -0.25 * (1 - L2);
        data.mEL4.N2[0][iP] = -0.25 * (1 - L1);

        data.mEL4.Nf[1][iP] = 0.25 * (1 + L1) * (1 - L2);
        data.mEL4.N1[1][iP] = 0.25 * (1 - L2);
        data.mEL4.N2[1][iP] = -0.25 * (1 + L1);

        data.mEL4.Nf[2][iP] = 0.25 * (1 + L1) * (1 + L2);
        data.mEL4.N1[2][iP] = 0.25 * (1 + L2);
        data.mEL4.N2[2][iP] = 0.25 * (1 + L1);

        data.mEL4.Nf[3][iP] = 0.25 * (1 - L1) * (1 + L2);
        data.mEL4.N1[3][iP] = -0.25 * (1 + L2);
        data.mEL4.N2[3][iP] = 0.25 * (1 - L1);
    }

    for (int i = 0; i < 4; i++) {
        data.mEL4.Sf[i].N_p = 2; // Dwa punkty całkowania na brzegach
        data.mEL4.Sf[i].P.resize(2);
        data.mEL4.Sf[i].W.resize(2, 1.0);
        data.mEL4.Sf[i].N1.resize(data.mEL4.nbn, std::vector<double>(2));
        data.mEL4.Sf[i].N2.resize(data.mEL4.nbn, std::vector<double>(2));
        data.mEL4.Sf[i].Nf.resize(data.mEL4.nbn, std::vector<double>(2));
        data.mEL4.Sf[i].UZEL.resize(2);
    }

    // Lokalne indexy węzłów
    data.mEL4.Sf[0].UZEL[0] = 4; data.mEL4.Sf[0].UZEL[1] = 1;
    data.mEL4.Sf[1].UZEL[0] = 1; data.mEL4.Sf[1].UZEL[1] = 2;
    data.mEL4.Sf[2].UZEL[0] = 2; data.mEL4.Sf[2].UZEL[1] = 3;
    data.mEL4.Sf[3].UZEL[0] = 3; data.mEL4.Sf[3].UZEL[1] = 4;

    // Współrzędne punktów całkowania dla brzegów
    data.mEL4.Sf[0].P[0].E = -1; data.mEL4.Sf[0].P[0].N = Alfa;
    data.mEL4.Sf[0].P[1].E = -1; data.mEL4.Sf[0].P[1].N = -Alfa;

    data.mEL4.Sf[1].P[0].E = -Alfa; data.mEL4.Sf[1].P[0].N = -1;
    data.mEL4.Sf[1].P[1].E = Alfa; data.mEL4.Sf[1].P[1].N = -1;

    data.mEL4.Sf[2].P[0].E = 1; data.mEL4.Sf[2].P[0].N = -Alfa;
    data.mEL4.Sf[2].P[1].E = 1; data.mEL4.Sf[2].P[1].N = Alfa;

    data.mEL4.Sf[3].P[0].E = Alfa; data.mEL4.Sf[3].P[0].N = 1;
    data.mEL4.Sf[3].P[1].E = -Alfa; data.mEL4.Sf[3].P[1].N = 1;

    // Funkcje kształtu dla brzegów
    for (int i = 0; i < 4; i++) {
        for (int iP = 0; iP < 2; iP++) {
            double e = data.mEL4.Sf[i].P[iP].E;
            double n = data.mEL4.Sf[i].P[iP].N;

            data.mEL4.Sf[i].Nf[0][iP] = 0.25 * (1 - e) * (1 - n);
            data.mEL4.Sf[i].Nf[1][iP] = 0.25 * (1 + e) * (1 - n);
            data.mEL4.Sf[i].Nf[2][iP] = 0.25 * (1 + e) * (1 + n);
            data.mEL4.Sf[i].Nf[3][iP] = 0.25 * (1 - e) * (1 + n);

            data.mEL4.Sf[i].N1[0][iP] = -0.25 * (1 - n);
            data.mEL4.Sf[i].N1[1][iP] = 0.25 * (1 - n);
            data.mEL4.Sf[i].N1[2][iP] = 0.25 * (1 + n);
            data.mEL4.Sf[i].N1[3][iP] = -0.25 * (1 + n);

            data.mEL4.Sf[i].N2[0][iP] = -0.25 * (1 - e);
            data.mEL4.Sf[i].N2[1][iP] = -0.25 * (1 + e);
            data.mEL4.Sf[i].N2[2][iP] = 0.25 * (1 + e);
            data.mEL4.Sf[i].N2[3][iP] = 0.25 * (1 - e);
        }
    }
}