#include "SOLVER.h"

typedef Eigen::Triplet<double> T;

void SOLVER() {
    // Lista tripletów do budowy macierzy rzadkiej
    std::vector<T> tripletList;
    tripletList.reserve(data.mGr.ne * 16); // Szacunkowa liczba tripletów

    // Inicjalizacja wektora obciążeń i wektora niewiadomych
    data.mB.setZero();
    data.mX.setZero();

    // Loop over elements
    for (int NEL = 1; NEL <= data.mGr.ne; ++NEL) {
        std::vector<int> nk(data.mEL4.nbn);

        // Przypisywanie numerów węzłów
        for (int i = 0; i < data.mEL4.nbn; ++i) {
            nk[i] = data.mGr.EL[NEL - 1].nop[i];
        }

         FeSM_heat(NEL);

        for (int i = 0; i < data.mEL4.nbn; ++i) {
            int ii = nk[i]; // Indeks wiersza
            for (int j = 0; j < data.mEL4.nbn; ++j) {
                int jj = nk[j]; // Indeks kolumny

                // Zakładam, że macierz jest symetryczna
                if (jj >= ii && (data.mLDA + ii - jj) <= data.mLDA) {
                    tripletList.emplace_back(ii, jj, data.est[i][j]);
                }
            }
            data.mB(ii) += data.r[i];
        }
    }

    // Wypełnianie macierzy rzadkiej za pomocą tripletów
    data.mA.setFromTriplets(tripletList.begin(), tripletList.end());
    data.mA.makeCompressed(); // Optymalizacja przechowywania

    // Rozwiązywanie układu równań za pomocą SimplicialLLT (dla symetrycznych, dodatnio określonych macierzy)
    Eigen::SimplicialLLT<Eigen::SparseMatrix<double>> solver;
    solver.compute(data.mA);

    if (solver.info() != Eigen::Success) {
        std::cerr << "Błąd podczas dekompozycji macierzy A." << std::endl;
        return;
    }

    // Rozwiązanie układu równań A x = B
    data.mX = solver.solve(data.mB);

    if (solver.info() != Eigen::Success) {
        std::cerr << "Błąd podczas rozwiązywania układu równań." << std::endl;
        return;
    }

    // Aktualizacja węzłów na podstawie rozwiązania
    for (int i = 0; i < data.mGr.nh; ++i) {
        data.mGr.ND[i].CR = (data.mGr.ND[i].t - data.mX(i)) / data.mdTime;
        data.mGr.ND[i].t = data.mX(i);
    }
}


