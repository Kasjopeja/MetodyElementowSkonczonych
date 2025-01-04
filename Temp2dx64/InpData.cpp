#include <iostream>
#include <fstream>
#include <vector>
#include "InpData.h"

extern GlobData data;

void InpData() {
    std::ifstream infile("C:\\Users\\Lenovo\\CLionProjects\\Metody Elementow Skonczonych\\Temp2dx64\\dane.txt");

    if (!infile) {
        std::cerr << "Blad: Nie udalo sie otworzyc pliku 'dane.txt'" << std::endl;
        return;
    }
    infile >> data.mTbegin >> data.mTime >> data.mdTime >> data.mT_otoczenia >> data.mAlfa
           >> data.mH0 >> data.mB0 >> data.mNhH >> data.mNhB >> data.mC >> data.mK >> data.mR;
    infile.close();

    std::cout << "Dane wczytane poprawnie z pliku 'dane.txt'" << std::endl;
}

void ALLOCATE_Matrix() {
    int mLDA = 0;
    int NeMaxB = 0;

    // Iterowanie przez elementy w siatce
    for (size_t NEL = 0; NEL < data.mGr.EL.size(); ++NEL) {
        const auto &element = data.mGr.EL[NEL];
        std::vector<int> nk(data.mEL4.nbn);

        // Przypisywanie numerów węzłów
        for (int i = 0; i < data.mEL4.nbn; ++i) {
            nk[i] = element.nop[i];
        }

        // Obliczanie szerokości pasma macierzy
        for (int i = 0; i < data.mEL4.nbn; ++i) {
            int ii = nk[i];
            for (int j = 0; j < data.mEL4.nbn; ++j) {
                int jj = nk[j];
                int jB = jj - ii + 1;
                if (jB >= mLDA) {
                    mLDA = jB;
                    NeMaxB = NEL;
                }
            }
        }
    }

    data.mLDA = mLDA;

    try {
        // Alokacja macierzy rzadkiej za pomocą Eigen
        data.mA.resize(data.mGr.nh, data.mGr.nh); // Rozmiar: nh x nh
        data.mA.reserve(Eigen::VectorXi::Constant(data.mGr.nh, data.mLDA)); // Rezerwacja miejsca

        // Inicjalizacja wektorów
        data.mB = Eigen::VectorXd::Zero(data.mGr.nh);
        data.mX = Eigen::VectorXd::Zero(data.mGr.nh);
    } catch (const std::bad_alloc &e) {
        std::cerr << "Błąd alokacji pamięci: " << e.what() << std::endl;
    }
}