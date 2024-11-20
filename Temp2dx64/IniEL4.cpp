#include "GlobData.cpp"
extern GlobData data;

void IniEL4()
{
    double Alfa = 1.0 / 2.0;
    data.mEL4.nbn = 4;   // Number of boundary nodes
    data.mEL4.N_p = 4;   // Number of integration points

    // Allocate memory for element arrays
    data.mEL4.N1.resize(data.mEL4.nbn, std::vector<double>(data.mEL4.N_p));
    data.mEL4.N2.resize(data.mEL4.nbn, std::vector<double>(data.mEL4.N_p));
    data.mEL4.Nf.resize(data.mEL4.nbn, std::vector<double>(data.mEL4.N_p));
    data.mEL4.P.resize(data.mEL4.N_p);
    data.mEL4.W.resize(data.mEL4.N_p, 1.0);
    data.mEL4.L.resize(data.mEL4.nbn);

    // Set coordinates for integration points
    data.mEL4.P[0] = {-Alfa, -Alfa}; // Point 1
    data.mEL4.P[1] = { Alfa, -Alfa}; // Point 2
    data.mEL4.P[2] = { Alfa,  Alfa}; // Point 3
    data.mEL4.P[3] = {-Alfa,  Alfa}; // Point 4

    // Set local coordinates for element nodes
    data.mEL4.L[0] = {-1, -1}; // Node 1
    data.mEL4.L[1] = { 1, -1}; // Node 2
    data.mEL4.L[2] = { 1,  1}; // Node 3
    data.mEL4.L[3] = {-1,  1}; // Node 4

    // Calculate shape functions and their derivatives
    for (int iP = 0; iP < data.mEL4.N_p; ++iP) {
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

    // Initialize data for boundary surfaces
    for (int i = 0; i < 4; ++i) {
        data.mEL4.Sf[i].N_p = 2; // Two integration points per surface
        data.mEL4.Sf[i].P.resize(2);
        data.mEL4.Sf[i].W.resize(2, 1.0);
        data.mEL4.Sf[i].N1.resize(data.mEL4.nbn, std::vector<double>(2));
        data.mEL4.Sf[i].N2.resize(data.mEL4.nbn, std::vector<double>(2));
        data.mEL4.Sf[i].Nf.resize(data.mEL4.nbn, std::vector<double>(2));
        data.mEL4.Sf[i].UZEL.resize(2);
    }

    // Set local node numbers for each boundary surface
    data.mEL4.Sf[0].UZEL = {4, 1}; // Surface 1
    data.mEL4.Sf[1].UZEL = {1, 2}; // Surface 2
    data.mEL4.Sf[2].UZEL = {2, 3}; // Surface 3
    data.mEL4.Sf[3].UZEL = {3, 4}; // Surface 4

    // Set integration points for boundary surfaces
    data.mEL4.Sf[0].P[0] = {-1,  Alfa}; data.mEL4.Sf[0].P[1] = {-1, -Alfa};
    data.mEL4.Sf[1].P[0] = {-Alfa, -1}; data.mEL4.Sf[1].P[1] = { Alfa, -1};
    data.mEL4.Sf[2].P[0] = { 1, -Alfa}; data.mEL4.Sf[2].P[1] = { 1,  Alfa};
    data.mEL4.Sf[3].P[0] = { Alfa,  1}; data.mEL4.Sf[3].P[1] = {-Alfa,  1};

    // Compute shape functions for boundary surfaces
    for (int i = 0; i < 4; ++i) {
        for (int iP = 0; iP < 2; ++iP) {
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