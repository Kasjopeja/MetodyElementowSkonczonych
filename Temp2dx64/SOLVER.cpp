#include "SOLVER.h"

typedef Eigen::Triplet<double> T;

void SOLVER() {
    // 1. Zerowanie macierzy i wektorów globalnych
    data.mA.setZero();
    data.mB.setZero();
    data.mX.setZero();
    std::cout << "Macierze i wektory globalne wyzerowane." << std::endl;

    // 2. Pętla po elementach – składanie (assembly) do macierzy globalnej
    for (int NEL = 0; NEL < data.mGr.ne; ++NEL) {
        std::array<int, 4> nk{}; // Lokalne indeksy węzłów elementu

        for (int i = 0; i < data.mEL4.nbn; i++) {
            nk[i] = data.mGr.EL[NEL].nop[i];
        }

        FeSM_heat(NEL); // Obliczanie macierzy lokalnych est(i,j) i wektora r(i)

        // Debug: wypisz macierz lokalną i wektor r dla elementu
        std::cout << "Element " << NEL << ": macierz lokalna i wektor r:" << std::endl;
        for (int i = 0; i < data.mEL4.nbn; ++i) {
            for (int j = 0; j < data.mEL4.nbn; ++j) {
                std::cout << "est[" << i << "][" << j << "] = " << data.est[i][j] << " ";
            }
            std::cout << "r[" << i << "] = " << data.r[i] << std::endl;
        }

        // Składanie do macierzy globalnej
        for (int i = 0; i < data.mEL4.nbn; i++) {
            int ii = nk[i];
            for (int j = 0; j < data.mEL4.nbn; j++) {
                int jj = nk[j];
                int iB = data.mLDA + ii - jj;

                if (jj >= ii && iB <= data.mLDA && iB >= 0) {
                    data.mA.coeffRef(ii, jj) += data.est[i][j];
                }
            }
            data.mB(ii) += data.r[i];
        }
    }

    // Debug: wypisz globalną macierz A i wektor B
    std::cout << "Globalna macierz A:" << std::endl;
    for (int k = 0; k < data.mA.outerSize(); ++k) {
        for (Eigen::SparseMatrix<double>::InnerIterator it(data.mA, k); it; ++it) {
            std::cout << "A(" << it.row() << ", " << it.col() << ") = " << it.value() << std::endl;
        }
    }
    std::cout << "Globalny wektor B:" << std::endl;
    for (int i = 0; i < data.mB.size(); ++i) {
        std::cout << "B[" << i << "] = " << data.mB[i] << std::endl;
    }

    // 3. Rozwiązywanie układu równań
    data.mA.makeCompressed();
    Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
    solver.analyzePattern(data.mA);
    solver.factorize(data.mA);
    if (solver.info() != Eigen::Success) {
        std::cerr << "SparseLU: factorization failed!" << std::endl;
        return;
    }

    data.mX = solver.solve(data.mB);
    if (solver.info() != Eigen::Success) {
        std::cerr << "SparseLU: solving failed!" << std::endl;
        return;
    }

    // Debug: wypisz rozwiązanie
    std::cout << "Rozwiązanie (wektor X):" << std::endl;
    for (int i = 0; i < data.mX.size(); ++i) {
        std::cout << "X[" << i << "] = " << data.mX[i] << std::endl;
    }

    // 4. Aktualizacja temperatur i pochodnej CR w każdym węźle
    for (int i = 0; i < data.mGr.nh; i++) {
        data.mGr.ND[i].CR = (data.mGr.ND[i].t - data.mX(i)) / data.mdTime;
        data.mGr.ND[i].t = data.mX(i);
    }

    // Debug: wypisz nowe wartości temperatur w węzłach
    std::cout << "Zaktualizowane temperatury i CR w węzłach:" << std::endl;
    for (int i = 0; i < data.mGr.nh; ++i) {
        std::cout << "Node " << i << ": T = " << data.mGr.ND[i].t << ", CR = " << data.mGr.ND[i].CR << std::endl;
    }
}



