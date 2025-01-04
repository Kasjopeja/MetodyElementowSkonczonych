#include "Jacob_2d.h"

// Funkcja odwracania macierzy 2x2 metodą eliminacji Gaussa
// (bez zmian względem poprzedniej wersji, ale można skrócić
//  do 2x2, jeśli kto woli — tutaj zachowujemy podejście NxN).
void Inv_MAT(int N, const std::vector<std::vector<double>>& Mat, std::vector<std::vector<double>>& Inv)
{
    std::vector<std::vector<double>> A(N, std::vector<double>(2*N, 0.0));

    // Inicjalizacja macierzy powiększonej
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            A[i][j] = Mat[i][j];
        }
        A[i][N + i] = 1.0;
    }

    // Eliminacja Gaussa
    for (int PivRow = 0; PivRow < N; ++PivRow) {
        double PivElt = A[PivRow][PivRow];

        if (std::abs(PivElt) < 1e-12) {
            int K = PivRow + 1;
            while (K < N && std::abs(A[K][PivRow]) < 1e-12) {
                ++K;
            }
            if (K == N) {
                std::cerr << "Error: Singular matrix, cannot compute inverse.\n";
                return;
            }
            std::swap(A[PivRow], A[K]);
            PivElt = A[PivRow][PivRow];
        }

        // Normalizacja wiersza głównego
        for (int j = 0; j < 2*N; ++j) {
            A[PivRow][j] /= PivElt;
        }

        // Eliminacja w pozostałych wierszach
        for (int TarRow = 0; TarRow < N; ++TarRow) {
            if (TarRow != PivRow) {
                double TarElt = A[TarRow][PivRow];
                for (int j = 0; j < 2*N; ++j) {
                    A[TarRow][j] -= A[PivRow][j] * TarElt;
                }
            }
        }
    }

    // Odczytanie odwrotności z macierzy powiększonej
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            Inv[i][j] = A[i][N + j];
        }
    }
}

// -----------------------------------------
// Jacob_2d: oblicza macierz Jakobiego (J_) i jej wyznacznik (DetJ)
//           w punkcie całkowania P dla elementu NEL.
//           Następnie wywołuje Inv_MAT(2, J_, J_inv) by uzyskać J_inv.

void Jacob_2d(std::vector<std::vector<double>>& J_,
              std::vector<std::vector<double>>& J_inv,
              int P, int N_p, int NEL,
              const std::vector<std::vector<double>>& N1,
              const std::vector<std::vector<double>>& N2,
              double& DetJ)
{

    // 2) W pętli sumujemy wkład dN/dξ * X(i) itd.
    //    UWAGA: Fortran indeksuje od 1, C++ od 0.
    //    Fortran: do i=1, NBN
    //    C++    : i=0; i< NBN
    for (int i = 0; i < data.mEL4.nbn; ++i) {
        // Indeks węzła globalnego:
        int globalNodeId = std::abs(data.mGr.EL[NEL - 1].nop[i]) - 1;

        // Współrzędne węzła i-tego
        double x = data.mGr.ND[globalNodeId].x;
        double y = data.mGr.ND[globalNodeId].y;

        // Pochodne funkcji kształtu w punkcie P
        double dN_dXi  = data.mEL4.N1[i][P];  // ~ dNi/dξ
        double dN_dEta = data.mEL4.N2[i][P]; // ~ dNi/dη

        // Odpowiedniki fortranowych:
        //   J_(1,1) += N1(i,P)*X(i)
        //   J_(1,2) += N1(i,P)*Y(i)
        //   J_(2,1) += N2(i,P)*X(i)
        //   J_(2,2) += N2(i,P)*Y(i)
        J_[0][0] += dN_dXi  * x;  // ∂x/∂ξ
        J_[0][1] += dN_dXi  * y;  // ∂y/∂ξ
        J_[1][0] += dN_dEta * x;  // ∂x/∂η
        J_[1][1] += dN_dEta * y;  // ∂y/∂η
    }

    // 3) Wyznacznik Jakobiego:
    //    DetJ = J_(1,1)*J_(2,2) - J_(1,2)*J_(2,1)
    //    Tu 0-based:
    DetJ = J_[0][0] * J_[1][1] - J_[0][1] * J_[1][0];

    // 4) Odwrócenie macierzy (matematycznie: J_inv = J_^{-1})
    Inv_MAT(2, J_, J_inv);

    // 5) Jeśli chcemy wartości J_, J_inv, DetJ zapisać
    //    do "globalnej" pamięci (np. używane w innym miejscu),
    //    to dopiszmy je w polach data:
    data.est[0][0] = J_[0][0];
    data.est[0][1] = J_[0][1];
    data.est[1][0] = J_[1][0];
    data.est[1][1] = J_[1][1];

    // Choć "est" było w oryginale macierzą do PDE (sztywność/pojemność),
    // tutaj używamy tylko tymczasowo, by zobrazować, że wartości się gdzieś
    // zapisują. W realnym kodzie lepiej wstawić do dedykowanych zmiennych,
    // np. data.J_ i data.J_inv, data.DetJ.
    data.r[0] = DetJ;
}
