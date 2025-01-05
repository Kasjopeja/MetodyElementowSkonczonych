#include "FeSM_heat.h"


// Główna funkcja odpowiednik Fortranowej: SUBROUTINE FeSM_heat(NEL)
void FeSM_heat(int NEL)
{
    // W Fortranie:
    //   est = 0
    //   r = 0
    //   CALL PRE_heat_mat(NEL)
    //   CALL PRE_heat_pov_mat(NEL)

    // Wyzerowanie tablic data.est i data.r
    // (po 4 węzły, ale w docelowym programie może być inna liczba)
    for(int i = 0; i < 4; i++)
    {
        data.r[i] = 0.0;
        for(int j = 0; j < 4; j++)
        {
            data.est[i][j] = 0.0;
        }
    }
    // Obliczenie wkładu "objętościowego" (macierzy i wektora)
    PRE_heat_mat(NEL);

    // Obliczenie wkładu warunków brzegowych
    PRE_heat_pov_mat(NEL);

}

// -----------------------------------------------------------------------------
// Funkcja pomocnicza: oblicza część wkładu do macierzy i wektora r
// (ten fragment w Fortranie: PRE_heat_mat(NEL)).
//
// Tłumaczenie z Fortrana:
//     - wczytanie X(i), Y(i), Temp_0(i) z węzłów
//     - pętla po punktach całkowania P
//     - obliczenia pochodnych N1, N2 (Ndx, Ndy)
//     - składanie do macierzy 'est' i wektora 'r'
void PRE_heat_mat(int NEL)
{
    // Najpierw wyzerujemy / zdefiniujemy pomocnicze tablice
    // (w Fortranie: REAL*8, DIMENSION(4):: X, Y, Temp_0, Ndx, Ndy)
    double Temp_0[4];
    std::vector<double> Y(data.mEL4.nbn);
    std::vector<double> X(data.mEL4.nbn);
    double Ndx[4], Ndy[4];

    // Odczyt współrzędnych i temperatur węzłów
    // (Uwaga: w Fortranie i=1..nbn, w C++ i=0..nbn-1)
    for(int i = 0; i < data.mEL4.nbn; ++i)
    {
        // Fortran: Id = ABS( mGr%EL(NEL)%NOP(I) )
        // U nas:   int Id = std::abs(data.mGr.EL[NEL].nop[i]) - 1;
        // Jeśli nop[] przechowuje 1-based index węzła.
        int Id = std::abs(data.mGr.EL[NEL].nop[i]);

        X[i]      = data.mGr.ND[Id].x;
        Y[i]      = data.mGr.ND[Id].y;
        Temp_0[i] = data.mGr.ND[Id].t;
    }

    // Pętla po punktach całkowania P=1..N_p (tu N_p = data.mEL4.N_p).
    for(int p = 0; p < data.mEL4.N_p; ++p)
    {
        // Obliczenie Jacobiego i jego odwrotności
        double J[2][2], J_inv[2][2];
        double DetJ = 0.0;

        Jacob_2d(J, J_inv,
                 p,
                 data.mEL4.N_p,
                 data.mEL4.nbn,
                 data.mEL4.N1,
                 data.mEL4.N2,
                 X, Y,
                 DetJ);

        // T0p (średnia temperatura w punkcie całkowania)
        double T0p = 0.0;

        // Wyliczenie Ndx, Ndy, zsumowanie T0p
        for(int i = 0; i < data.mEL4.nbn; ++i)
        {
            // Fortran: Ndx(i) = N1(i,p)*J_inv(1,1) + N2(i,p)*J_inv(1,2)
            // w C++:
            double Ni1 = data.mEL4.N1[i][p];
            double Ni2 = data.mEL4.N2[i][p];

            Ndx[i] = Ni1 * J_inv[0][0] + Ni2 * J_inv[0][1];
            Ndy[i] = Ni1 * J_inv[1][0] + Ni2 * J_inv[1][1];

            double NiShape = data.mEL4.Nf[i][p];
            T0p += Temp_0[i] * NiShape;
        }

        // Modyfikacja DetJ: w Fortranie -> detJ = abs(detJ)*mEL4%W(p)
        DetJ = std::abs(DetJ) * data.mEL4.W[p];

        // Teraz składanie do macierzy est(n,i) i wektora r(n).
        for(int n = 0; n < data.mEL4.nbn; ++n)
        {
            for(int i = 0; i < data.mEL4.nbn; ++i)
            {
                double Ni = data.mEL4.Nf[i][p];
                double Nn = data.mEL4.Nf[n][p];

                double Hin = data.mK * (Ndx[n]*Ndx[i] + Ndy[n]*Ndy[i]) * DetJ;
                double Cin = data.mC * data.mR * Nn * Ni * DetJ;

                data.est[n][i] += Hin + Cin / data.mdTime;
                data.r[n]      += (Cin / data.mdTime) * T0p;

            }

            // Przykład ewentualnego dodatkowego źródła Q:
            //  Q = 0.8 * T0p * Hp  ...
            //  data.r[n] += (Nn * Q * DetJ);
        }

    }
}

// -----------------------------------------------------------------------------
// Funkcja pomocnicza: warunki brzegowe (fortran: PRE_heat_pov_mat(NEL))
void PRE_heat_pov_mat(int NEL)
{
    double X[4], Y[4];
    for (int I = 0; I < 4; I++) {
        int Id = std::abs(data.mGr.EL[NEL].nop[I]);
        X[I] = data.mGr.ND[Id].x;
        Y[I] = data.mGr.ND[Id].y;
    }

    int nPov = data.mGr.EL[NEL].Npov;
    for (int iPov = 0; iPov < nPov; iPov++) {
        int id = data.mGr.EL[NEL].aPov[iPov];
        double DetJ = 0.0;

        switch (id) {
            case 1:
                DetJ = std::sqrt((X[3] - X[0]) * (X[3] - X[0]) + (Y[3] - Y[0]) * (Y[3] - Y[0]));
                break;
            case 2:
                DetJ = std::sqrt((X[0] - X[1]) * (X[0] - X[1]) + (Y[0] - Y[1]) * (Y[0] - Y[1]));
                break;
            case 3:
                DetJ = std::sqrt((X[1] - X[2]) * (X[1] - X[2]) + (Y[1] - Y[2]) * (Y[1] - Y[2]));
                break;
            case 4:
                DetJ = std::sqrt((X[2] - X[3]) * (X[2] - X[3]) + (Y[2] - Y[3]) * (Y[2] - Y[3]));
                break;
            default:
                break;
        }

        std::cout << "Powierzchnia: " << id << ", DetJ: " << DetJ << std::endl;

        for (int P = 0; P < 2; ++P) {
            for (int n = 0; n < 4; ++n) {
                for (int i = 0; i < 4; ++i) {
                    double Ni = data.mEL4.Sf[id-1].Nf[i][P];
                    double Nn = data.mEL4.Sf[id-1].Nf[n][P];
                    data.est[n][i] += data.mAlfa * Nn * Ni * DetJ;
                }
                double Pn = data.mAlfa * data.mT_otoczenia * data.mEL4.Sf[id-1].Nf[n][P] * DetJ;
                data.r[n] += Pn;
                std::cout << "Pn = " << Pn << " dla powierzchni " << id << std::endl;
            }
        }
    }
}


