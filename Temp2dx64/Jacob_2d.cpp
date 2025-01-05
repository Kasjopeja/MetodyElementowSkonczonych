#include "Jacob_2d.h"

void Inv_MAT(int N,
             const std::vector<std::vector<double>> &Mat,
             std::vector<std::vector<double>> &Inv)
{
    // ----------------------------------------------------------------------
    // Oryginał fortranowy tworzył tablicę A o wymiarze N × (2N)
    // i do jej lewej połowy kopiował Mat, a do prawej połowy — macierz jednostkową.
    // Potem wykonywał eliminację Gaussa-Jordana.
    // ----------------------------------------------------------------------

    // Sprawdzamy, czy rozmiar Mat i Inv to NxN
    // (dla bezpieczeństwa).
    if ((int)Mat.size() != N || (int)Mat[0].size() != N) {
        std::cerr << "Inv_MAT: wymiar macierzy Mat jest nieprawidłowy!\n";
        return;
    }
    if ((int)Inv.size() != N || (int)Inv[0].size() != N) {
        std::cerr << "Inv_MAT: wymiar macierzy Inv jest nieprawidłowy!\n";
        return;
    }

    // Tworzymy "A" o wymiarze N × (2N). W C++: wektor wektorów.
    std::vector<std::vector<double>> A(N, std::vector<double>(2*N, 0.0));

    // Inicjalizacja:
    //   A[ i ][ j ] = Mat[i][j] dla j w [0..N-1],
    //   A[ i ][ N+i ] = 1.0 (macierz jednostkowa w prawej części).
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            A[i][j] = Mat[i][j];
        }
        A[i][N + i] = 1.0; // jednostkowa w kolumnie N+i
    }

    // Gauss-Jordan
    for (int pivRow = 0; pivRow < N; pivRow++)
    {
        double pivElt = A[pivRow][pivRow];
        if (pivElt == 0.0) {
            // Znajdź wiersz K > pivRow, w którym A[K][pivRow] != 0
            int K = pivRow + 1;
            while (K < N && A[K][pivRow] == 0.0) {
                K++;
            }
            if (K == N) {
                // Nie znaleziono niezerowego pivota
                std::cerr << "Inv_MAT: Couldn't find a non-zero pivot. Solution is rubbish.\n";
                return;
            } else {
                // Zamieniamy wiersze pivRow i K
                std::swap(A[pivRow], A[K]);
                // Ustawiamy pivot
                pivElt = A[pivRow][pivRow];
            }
        }

        // Dzielenie całego wiersza przez pivElt
        for (int col = 0; col < 2*N; col++) {
            A[pivRow][col] /= pivElt;
        }

        // Zerowanie elementów w tej samej kolumnie w innych wierszach
        for (int tarRow = 0; tarRow < N; tarRow++) {
            if (tarRow != pivRow) {
                double tarElt = A[tarRow][pivRow];
                for (int col = 0; col < 2*N; col++) {
                    A[tarRow][col] -= A[pivRow][col] * tarElt;
                }
            }
        }
    }

    // Po zakończeniu – prawa połowa A to macierz odwrotna
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            Inv[i][j] = A[i][N + j];
        }
    }
}

void Jacob_2d(double J[2][2],
              double J_inv[2][2],
              int p,
              int N_p,
              int NBN,
              const std::vector<std::vector<double>> &N1,
              const std::vector<std::vector<double>> &N2,
              const std::vector<double> &X,
              const std::vector<double> &Y,
              double &DetJ)
{
    for (int row = 0; row < 2; row++) {
        for (int col = 0; col < 2; col++) {
            J[row][col] = 0.0;
            J_inv[row][col] = 0.0;
        }
    }

    for (int i = 0; i < NBN; i++) {
        J[0][0] += N1[i][p] * X[i];
        J[0][1] += N1[i][p] * Y[i];
        J[1][0] += N2[i][p] * X[i];
        J[1][1] += N2[i][p] * Y[i];
    }

    DetJ = J[0][0] * J[1][1] - J[0][1] * J[1][0];
    if (DetJ <= 0.0) {
        std::cerr << "Błąd: Ujemny lub zerowy DetJ = " << DetJ << std::endl;
    }

    std::vector<std::vector<double>> Mat(2, std::vector<double>(2, 0.0));
    Mat[0][0] = J[0][0];
    Mat[0][1] = J[0][1];
    Mat[1][0] = J[1][0];
    Mat[1][1] = J[1][1];

    std::vector<std::vector<double>> Inv(2, std::vector<double>(2, 0.0));

    Inv_MAT(2, Mat, Inv);

    J_inv[0][0] = Inv[0][0];
    J_inv[0][1] = Inv[0][1];
    J_inv[1][0] = Inv[1][0];
    J_inv[1][1] = Inv[1][1];
}
