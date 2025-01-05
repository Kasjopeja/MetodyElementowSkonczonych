#include "SOLVER.h"
#define DEBUG1

typedef Eigen::Triplet<double> T;

void SOLVER() {
    // 1. Zerowanie macierzy i wektorów globalnych
    data.mA.setZero();
    data.mB.setZero();
    data.mX.setZero();
    //std::cout << "Macierze i wektory globalne wyzerowane." << std::endl;

    // 2. Pętla po elementach – składanie (assembly) do macierzy globalnej
    for (int NEL = 0; NEL < data.mGr.ne; ++NEL) {
        int nk[4]; // Lokalne indeksy węzłów elementu

        for (int i = 0; i < data.mEL4.nbn; i++) {
            nk[i] = data.mGr.EL[NEL].nop[i];
            //std::cout << nk[i] << " " <<std::endl;
        }

        //std::cout << "---------" <<std::endl;

        FeSM_heat(NEL); // Obliczanie macierzy lokalnych est(i,j) i wektora r(i)

        // Debug: wypisz macierz lokalną i wektor r dla elementu

#ifdef DEBUG
        std::cout << "Element " << NEL << ": macierz lokalna i wektor r:" << std::endl;
        for (int i = 0; i < data.mEL4.nbn; ++i) {
            for (int j = 0; j < data.mEL4.nbn; ++j) {
                std::cout << "est[" << i << "][" << j << "] = " << data.est[i][j] << " ";
            }
            std::cout << "r[" << i << "] = " << data.r[i] << std::endl;
        }
#endif
        // Składanie do macierzy globalnej
        for (int i = 0; i < data.mEL4.nbn; ++i) {
            int ii = nk[i];
            for (int j = 0; j < data.mEL4.nbn; ++j) {
                int jj = nk[j];

                    data.mA.coeffRef(ii, jj) += data.est[i][j];
                    //std::cout << data.mA.coeffRef(ii, jj) << std::endl;
            }
            data.mB(ii) += data.r[i];
        }
    }
#ifdef DEBUG
    // Debug: wypisz globalną macierz A i wektor B
    std::cout << "Globalna macierz A (wypisano wszystkie wartości, w tym zera):" << std::endl;
    for (int i = 0; i < data.mA.rows(); ++i) {
        for (int j = 0; j < data.mA.cols(); ++j) {
            double value = data.mA.coeff(i, j); // Uzyskanie wartości na pozycji (i, j)
            std::cout << "A(" << i << ", " << j << ") = " << value << std::endl;
        }
    }

    std::cout << "Globalny wektor B:" << std::endl;
    for (int i = 0; i < data.mB.size(); ++i) {
        std::cout << "B[" << i << "] = " << data.mB[i] << std::endl;
    }
#endif

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

#ifdef DEBUG
    // Debug: wypisz rozwiązanie
    std::cout << "Rozwiązanie (wektor X):" << std::endl;
    for (int i = 0; i < data.mX.size(); ++i) {
        std::cout << "X[" << i << "] = " << data.mX[i] << std::endl;
    }

#endif

    // 4. Aktualizacja temperatur i pochodnej CR w każdym węźle
    for (int i = 0; i < data.mGr.nh; i++) {
        data.mGr.ND[i].CR = (data.mGr.ND[i].t - data.mX(i)) / data.mdTime;
        data.mGr.ND[i].t = data.mX(i);
    }

#ifdef DEBUG
    // Debug: wypisz nowe wartości temperatur w węzłach
    std::cout << "Zaktualizowane temperatury i CR w węzłach:" << std::endl;
    for (int i = 0; i < data.mGr.nh; ++i) {
        std::cout << "Node " << i << ": T = " << data.mGr.ND[i].t << ", CR = " << data.mGr.ND[i].CR << std::endl;
    }

#endif
}



